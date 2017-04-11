!==========================================================================================!
! Forestry.f90. These subroutines calculate the area to be taken from each patch to meet   !
!               the biomass demand, in case this site is to be harvest using biomass       !
!               targets.                                                                   !
!                                                                                          !
!               The site is harvested up to a target biomass.  Patches are harvested first !
!               from those above the harvest_age, with equal rates.  If the biomass target !
!               is not met, the remaining biomass is harvested from the patches below the  !
!               minimum age, starting with the oldest.  Harvest targets are taken from the !
!               land use transition matrix dataset, which must be formatted as George      !
!               Hurtt's GLM output.                                                        !
!                                                                                          !
! NOTICE:  These subroutines have not been thoroughly tested. Please                       !
!          report problems to David Medvigy, medvigy@post.harvard.edu.                     !
!==========================================================================================!
module forestry



   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the disturbance rates associated with logging, when biomass !
   ! demands are provided instead of actual disturbance rates.                             !
   !---------------------------------------------------------------------------------------!
   subroutine find_lambda_harvest(cpoly,isi,onsp,harvestable_agb,lambda_harvest)
      use ed_state_vars        , only : polygontype                & ! structure
                                      , sitetype                   & ! structure
                                      , patchtype                  & ! structure
                                      , allocate_sitetype          & ! subroutine
                                      , deallocate_sitetype        & ! subroutine
                                      , copy_sitetype              ! ! subroutine
      use disturb_coms         , only : ianth_disturb              & ! intent(in)
                                      , lutime                     & ! intent(in)
                                      , min_patch_area             & ! intent(in)
                                      , plantation_rotation        & ! intent(in)
                                      , min_harvest_biomass        & ! intent(in)
                                      , mature_harvest_age         & ! intent(in)
                                      , min_oldgrowth              ! ! intent(in)
      use fuse_fiss_utils      , only : terminate_patches          ! ! subroutine
      use ed_max_dims          , only : n_pft                      & ! intent(in)
                                      , n_dbh                      ! ! intent(in)
      use grid_coms            , only : nzg                        & ! intent(in)
                                      , nzs                        ! ! intent(in)
      use detailed_coms        , only : idetailed                  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)             , target        :: cpoly
      integer                       , intent(in)    :: isi
      integer                       , intent(in)    :: onsp
      real, dimension(onsp)         , intent(inout) :: harvestable_agb
      real, dimension(onsp)         , intent(inout) :: lambda_harvest
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)                , pointer       :: csite
      type(patchtype)               , pointer       :: cpatch
      integer                                       :: ipa
      integer                                       :: ico
      integer                                       :: ipft
      integer                                       :: ilu
      logical                                       :: is_oldgrowth
      logical                                       :: is_mature
      logical                                       :: print_detailed
      real                                          :: primary_harvest_target
      real                                          :: secondary_harvest_target
      real                                          :: site_harvest_target
      real                                          :: site_harvestable_agb
      real                                          :: area_mature_primary
      real                                          :: area_mature_secondary
      real                                          :: area_mature_plantation
      real                                          :: hvagb_mature_primary
      real                                          :: hvagb_mature_secondary
      real                                          :: hvagb_mature_plantation
      real                                          :: harvest_deficit
      !------------------------------------------------------------------------------------!


      !----- Nothing to do here in case anthropogenic disturbance is turned off. ----------!
      if (ianth_disturb == 0) return
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Find out whether to print detailed information on screen.                    !
      !------------------------------------------------------------------------------------!
      print_detailed = btest(idetailed,6)
      !------------------------------------------------------------------------------------!


      !----- Link to the current site. ----------------------------------------------------!
      csite => cpoly%site(isi)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Set biomass targets based on current rates and unapplied harvest from         !
      ! previous years (memory).  These are in kgC/m2.  In case DBH-based logging is       !
      ! applied, then harvest target is set to zero.                                       !
      !------------------------------------------------------------------------------------!
      primary_harvest_target   = cpoly%primary_harvest_target  (isi)                       &
                               + cpoly%primary_harvest_memory  (isi)
      secondary_harvest_target = cpoly%secondary_harvest_target(isi)                       &
                               + cpoly%primary_harvest_memory  (isi)
      site_harvest_target      = primary_harvest_target + secondary_harvest_target
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find total harvestable biomass density in kgC/m2.  The harvestable biomass is !
      ! not necessarily the same as the total site biomass, because it is possible that    !
      ! the demand is for trees above a minimum size (which is very common practice in the !
      ! tropics).                                                                          !
      !------------------------------------------------------------------------------------!
      site_harvestable_agb = 0.
      hpat_loop: do ipa=1,onsp
         cpatch => csite%patch(ipa)
         ilu = csite%dist_type(ipa)

         !---------------------------------------------------------------------------------!
         !      Save some flags here so we can account for the different types of patches. !
         !---------------------------------------------------------------------------------!
         ilu           = csite%dist_type(ipa)
         is_oldgrowth  = csite%age(ipa) >= min_oldgrowth(ilu)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Go over each cohort, seek harvestable biomass.  We don't bother looking at  !
         ! croplands and pastures.                                                         !
         !---------------------------------------------------------------------------------!
         select case(ilu)
         case (2:7)
            hcoh_loop: do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)
               if (cpatch%dbh(ico) >= cpoly%mindbh_harvest(ipft,isi)) then
                  !----- Cohort is harvestable. -------------------------------------------!
                  harvestable_agb(ipa) = harvestable_agb(ipa)                              &
                                       + cpatch%nplant(ico) * cpatch%btimber(ico)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do hcoh_loop
         end select
         !---------------------------------------------------------------------------------!

         !----- Update site biomass only when the patch has sufficient biomass. -----------!
         if (harvestable_agb(ipa) >= min_harvest_biomass) then
            site_harvestable_agb = site_harvestable_agb                                    &
                                 + harvestable_agb(ipa) * csite%area(ipa)
         end if
         !---------------------------------------------------------------------------------!
      end do hpat_loop
      !------------------------------------------------------------------------------------!



      !====================================================================================!
      !====================================================================================!
      !     Find out whether any harvest can occur at this site.  For harvest to occur,    !
      ! the site must meet two criteria:                                                   !
      ! a. Harvestable biomass must be greater than min_harvest_biomass;                   !
      ! b. Target biomass must be greater than min_harvest_biomass.                        !
      !                                                                                    !
      !    In case one or both criteria are not met, we add the harvest target to memory.  !
      !------------------------------------------------------------------------------------!
      if ( site_harvestable_agb <  min_harvest_biomass                 .or.                &
           site_harvest_target  <= site_harvestable_agb * min_patch_area      ) then
         cpoly%primary_harvest_target  (isi) = 0.0
         cpoly%secondary_harvest_target(isi) = 0.0
         cpoly%primary_harvest_memory  (isi) = primary_harvest_target
         cpoly%secondary_harvest_memory(isi) = secondary_harvest_target
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !       Print the inventory and the target.                                       !
         !---------------------------------------------------------------------------------!
         if (print_detailed) then
            write (unit=*,fmt='(a)'      )     ' '
            write (unit=*,fmt='(a)'      )     '------------------------------------------'
            write (unit=*,fmt='(a)'      )     ' FORESTRY.  NO HARVEST THIS YEAR...'
            write (unit=*,fmt='(a)'      )     ' '
            write (unit=*,fmt='(a,1x,i5)')     ' ISI                   = ',isi
            write (unit=*,fmt='(a,1x,es12.5)') ' PRIMARY TARGET AGB    = '                 &
                                                            , primary_harvest_target
            write (unit=*,fmt='(a,1x,es12.5)') ' SECONDARY TARGET AGB  = '                 &
                                                            , secondary_harvest_target
            write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL TARGET AGB      = '                 &
                                                            , site_harvest_target
            write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL HARVESTABLE AGB = '                 &
                                                            , site_harvestable_agb
            write (unit=*,fmt='(a)'      )     '------------------------------------------'
            write (unit=*,fmt='(a)'      )     ' '
         end if
         !---------------------------------------------------------------------------------!
         return
      end if
      !------------------------------------------------------------------------------------!



      !------ Compute current stocks of timber in mature forests. -------------------------!
      call inventory_mature_forests(cpoly,isi,onsp,harvestable_agb                         &
                                   ,area_mature_primary   ,hvagb_mature_primary            &
                                   ,area_mature_secondary ,hvagb_mature_secondary          &
                                   ,area_mature_plantation,hvagb_mature_plantation)
      !------------------------------------------------------------------------------------!



      !------ Compute the mature-forest harvest rates. ------------------------------------!
      call mature_forest_harvest(cpoly,isi,onsp                                            &
                                ,hvagb_mature_primary,hvagb_mature_secondary               &
                                ,hvagb_mature_plantation,primary_harvest_target            &
                                ,secondary_harvest_target,harvestable_agb,lambda_harvest   &
                                ,harvest_deficit)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Compute the disturbance rate applied to young patches to help meet the        !
      ! biomass demands.                                                                   !
      !------------------------------------------------------------------------------------!
      call young_forest_harvest(cpoly,isi,onsp,harvestable_agb,lambda_harvest              &
                               ,harvest_deficit)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Print the inventory and the target.                                          !
      !------------------------------------------------------------------------------------!
      if (print_detailed) then
         write (unit=*,fmt='(a)'      )     ' '
         write (unit=*,fmt='(a)'      )     '---------------------------------------------'
         write (unit=*,fmt='(a)'      )     ' FORESTRY.  HARVEST RATES'
         write (unit=*,fmt='(a)'      )     ' '
         write (unit=*,fmt='(a,1x,i5)')     ' ISI                    = ',isi
         write (unit=*,fmt='(a,1x,es12.5)') ' PRIMARY TARGET AGB     = '                   &
                                                         , primary_harvest_target
         write (unit=*,fmt='(a,1x,es12.5)') ' SECONDARY TARGET AGB   = '                   &
                                                         , secondary_harvest_target
         write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL TARGET AGB       = '                   &
                                                         , site_harvest_target
         write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL HARVESTABLE AGB  = '                   &
                                                         , site_harvestable_agb
         write (unit=*,fmt='(a,1x,es12.5)') ' HV AGB (PRIMARY)       = '                   &
                                                         , hvagb_mature_primary
         write (unit=*,fmt='(a,1x,es12.5)') ' HV AGB (SECONDARY)     = '                   &
                                                         , hvagb_mature_secondary
         write (unit=*,fmt='(a,1x,es12.5)') ' HV AGB (PLANTATION)    = '                   &
                                                         , hvagb_mature_plantation
         write (unit=*,fmt='(a,1x,es12.5)') ' HV AREA (PRIMARY)      = '                   &
                                                         , area_mature_primary
         write (unit=*,fmt='(a,1x,es12.5)') ' HV AREA (SECONDARY)    = '                   &
                                                         , area_mature_secondary
         write (unit=*,fmt='(a,1x,es12.5)') ' HV AREA (PLANTATION)   = '                   &
                                                         , area_mature_plantation
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a,1x,es12.5)') ' HARVEST DEFICIT = ', harvest_deficit
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a)'          ) '---------------------------------------------'
         write (unit=*,fmt='(7(a,1x))'    ) '  IPA','   LU','         AGE','        AREA'  &
                                            ,'      LAMBDA','      MATURE','  OLD_GROWTH'
         write (unit=*,fmt='(a)'      )     '---------------------------------------------'
         do ipa=1,onsp

            !------------------------------------------------------------------------------!
            !      Save some flags here so we can account for the different types of       !
            ! patches.                                                                     !
            !------------------------------------------------------------------------------!
            ilu           = csite%dist_type(ipa)
            select case (ilu)
            case (2)
               is_mature    = csite%age(ipa) >= plantation_rotation
               is_oldgrowth = .false.
            case (3:7)
               is_mature    = csite%age(ipa) >= mature_harvest_age
               is_oldgrowth = csite%age(ipa) >= min_oldgrowth(ilu)
            case default
               is_mature    = .false.
               is_oldgrowth = .false.
            end select
            !------------------------------------------------------------------------------!
            write (unit=*,fmt='(2(i5,1x),3(f12.7,1x),2(11x,l1,1x))')                       &
               ipa,csite%dist_type(ipa),csite%age(ipa),csite%area(ipa)                     &
                  ,lambda_harvest(ipa),is_mature,is_oldgrowth
         end do
         write (unit=*,fmt='(a)'      )     '---------------------------------------------'
         write (unit=*,fmt='(a)'      )     ' '
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Reset the primary forest memory, and assign any remaining deficit to the      !
      ! secondary forest memory.                                                           !
      !------------------------------------------------------------------------------------!
      cpoly%primary_harvest_memory  (isi) = 0.0
      cpoly%secondary_harvest_memory(isi) = harvest_deficit
      !------------------------------------------------------------------------------------!

      return
   end subroutine find_lambda_harvest
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This sub-routine calculates the area and total biomass associated with primary     !
   ! forests, secondary forests, and forest plantations (both mature and young).           !
   !---------------------------------------------------------------------------------------!
   subroutine inventory_mature_forests(cpoly,isi,onsp,harvestable_agb                      &
                                      ,area_mature_primary   , hvagb_mature_primary        &
                                      ,area_mature_secondary , hvagb_mature_secondary      &
                                      ,area_mature_plantation, hvagb_mature_plantation)
      use ed_state_vars , only : polygontype         & ! structure
                               , sitetype            ! ! structure
      use disturb_coms  , only : plantation_rotation & ! intent(in)
                               , mature_harvest_age  & ! intent(in)
                               , min_harvest_biomass & ! intent(in)
                               , min_oldgrowth       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)    , target      :: cpoly
      integer              , intent(in)  :: isi
      integer              , intent(in)  :: onsp
      real, dimension(onsp), intent(in)  :: harvestable_agb
      real                 , intent(out) :: area_mature_primary
      real                 , intent(out) :: hvagb_mature_primary
      real                 , intent(out) :: area_mature_secondary
      real                 , intent(out) :: hvagb_mature_secondary
      real                 , intent(out) :: area_mature_plantation
      real                 , intent(out) :: hvagb_mature_plantation
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer     :: csite
      integer                            :: ipa
      integer                            :: ilu
      logical                            :: is_oldgrowth
      logical                            :: is_mature
      logical                            :: is_rotation
      !------------------------------------------------------------------------------------!



      !----- Initialize inventory. --------------------------------------------------------!
      area_mature_primary     = 0.0
      area_mature_secondary   = 0.0
      area_mature_plantation  = 0.0
      hvagb_mature_primary    = 0.0
      hvagb_mature_secondary  = 0.0
      hvagb_mature_plantation = 0.0
      !------------------------------------------------------------------------------------!

      csite => cpoly%site(isi)
      patchloop: do ipa=1,csite%npatches

         !----- Skip the patch if the biomass is low. -------------------------------------!
         if (harvestable_agb(ipa) < min_harvest_biomass) cycle patchloop
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Save some flags here so we can account for the different types of patches. !
         !---------------------------------------------------------------------------------!
         ilu           = csite%dist_type(ipa)
         is_oldgrowth  = csite%age(ipa) >= min_oldgrowth(ilu)
         is_mature     = csite%age(ipa) >= mature_harvest_age
         is_rotation   = csite%age(ipa) >= plantation_rotation
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Run some checks to determine the type of forest and whether it has reached !
         ! maturity.  TO DO: check whether George Hurtt's data consider burnt patches      !
         ! primary or secondary.                                                           !
         !---------------------------------------------------------------------------------!
         select case (csite%dist_type(ipa))
         case (2)
            !---- Forest plantation. ------------------------------------------------------!
            if (is_rotation) then
               area_mature_plantation  = area_mature_plantation + csite%area(ipa)
               hvagb_mature_plantation = hvagb_mature_plantation                           &
                                       + harvestable_agb(ipa)   * csite%area(ipa)
            end if
            !------------------------------------------------------------------------------!

         case (3)
            !---- Treefall always goes to "primary" forest. -------------------------------!
            if (is_mature) then
               area_mature_primary     = area_mature_primary    + csite%area(ipa)
               hvagb_mature_primary    = hvagb_mature_primary                              &
                                       + harvestable_agb(ipa)   * csite%area(ipa)
            end if
            !------------------------------------------------------------------------------!

         case (4:7)
            !------------------------------------------------------------------------------!
            !     Other disturbances, assume secondary in case the patch is not considered !
            ! old-growth, and primary otherwise.                                           !
            !------------------------------------------------------------------------------!
            if (is_mature .and. is_oldgrowth) then
               area_mature_primary     = area_mature_primary    + csite%area(ipa)
               hvagb_mature_primary    = hvagb_mature_primary                              &
                                       + harvestable_agb(ipa)   * csite%area(ipa)
            else if (is_mature) then
               area_mature_secondary   = area_mature_secondary  + csite%area(ipa)
               hvagb_mature_secondary  = hvagb_mature_secondary                            &
                                       + harvestable_agb(ipa)   * csite%area(ipa)
            end if
            !------------------------------------------------------------------------------!
         case default
            !----- Probably pasture/croplands, skip the patch. ----------------------------!
            continue
            !------------------------------------------------------------------------------!
         end select 
         !---------------------------------------------------------------------------------!
      end do patchloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine inventory_mature_forests
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine adjusts the mature forest disturbance rates associated with      !
   ! harvest.  In case the target biomass for primary forests is not met, it first tries   !
   ! to take the remaining biomass from old secondary forests before going to immature     !
   ! forests.                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine mature_forest_harvest(cpoly,isi,onsp,hvagb_mature_primary                    &
                                   ,hvagb_mature_secondary,hvagb_mature_plantation         &
                                   ,primary_harvest_target,secondary_harvest_target        &
                                   ,harvestable_agb,lambda_harvest,harvest_deficit)
      use ed_state_vars , only : polygontype         & ! structure
                               , sitetype            ! ! structure
      use disturb_coms  , only : plantation_rotation & ! intent(in)
                               , mature_harvest_age  & ! intent(in)
                               , min_harvest_biomass & ! intent(in)
                               , min_oldgrowth       ! ! intent(in)
      use consts_coms   , only : lnexp_max
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)    , target        :: cpoly
      integer              , intent(in)    :: isi
      integer              , intent(in)    :: onsp
      real                 , intent(in)    :: hvagb_mature_primary
      real                 , intent(in)    :: hvagb_mature_secondary
      real                 , intent(in)    :: hvagb_mature_plantation
      real                 , intent(in)    :: primary_harvest_target
      real                 , intent(inout) :: secondary_harvest_target
      real, dimension(onsp), intent(in)    :: harvestable_agb
      real, dimension(onsp), intent(inout) :: lambda_harvest
      real                 , intent(out)   :: harvest_deficit
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer       :: csite
      integer                              :: ipa
      integer                              :: ilu
      logical                              :: is_oldgrowth
      logical                              :: is_mature
      logical                              :: is_rotation
      real                                 :: lambda_mature_primary
      real                                 :: lambda_mature_plantation
      real                                 :: lambda_mature_secondary
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find harvesting rate in mature primary vegetation.  In case there is not enough !
      ! biomass harvest, harvest all primary vegetation then add the unmet biomass to the  !
      ! target for secondary vegetation.                                                   !
      !------------------------------------------------------------------------------------!
      if (hvagb_mature_primary > primary_harvest_target) then
         lambda_mature_primary = primary_harvest_target / hvagb_mature_primary
      else
         lambda_mature_primary    = lnexp_max
         harvest_deficit          = primary_harvest_target   - hvagb_mature_primary
         secondary_harvest_target = secondary_harvest_target + harvest_deficit
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find harvesting rate for mature plantations and mature secondary forests.     !
      ! First try to remove all biomass from plantations.  In case there isn't sufficient  !
      ! biomass, harvest all plantations then remove the unmet biomass from mature         !
      ! secondary vegetation.  In case there isn't enough biomass, leave the remaining     !
      ! target in harvest_deficit, which will be used to determine harvesting from young   !
      !  forests.                                                                          !
      !------------------------------------------------------------------------------------!
      if (hvagb_mature_plantation > secondary_harvest_target) then
         lambda_mature_plantation = secondary_harvest_target / hvagb_mature_plantation
         lambda_mature_secondary  = 0.0
         harvest_deficit          = 0.0
      else
         lambda_mature_plantation = lnexp_max
         harvest_deficit          = secondary_harvest_target - hvagb_mature_plantation

         if(hvagb_mature_secondary > harvest_deficit) then
            lambda_mature_secondary = harvest_deficit / hvagb_mature_secondary
            harvest_deficit         = 0.0
         else
            lambda_mature_secondary = lnexp_max
            harvest_deficit         = harvest_deficit - hvagb_mature_secondary
         end if
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over patches.                                                            ! 
      !------------------------------------------------------------------------------------!
      csite => cpoly%site(isi)
      patch_loop: do ipa=1,onsp
     
         !----- Skip patch in case the biomass is less than the minimum for harvesting. ---!
         if (harvestable_agb(ipa) < min_harvest_biomass) cycle patch_loop
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Save some flags here so we can account for the different types of patches. !
         !---------------------------------------------------------------------------------!
         ilu           = csite%dist_type(ipa)
         is_oldgrowth  = csite%age(ipa) >= min_oldgrowth(ilu)
         is_mature     = csite%age(ipa) >= mature_harvest_age
         is_rotation   = csite%age(ipa) >= plantation_rotation
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Find out whether to harvest this patch.                                      !
         !---------------------------------------------------------------------------------!
         select case (csite%dist_type(ipa))
         case (2)
            !----- Forest plantation. -----------------------------------------------------!
            if ( is_rotation ) then
               lambda_harvest(ipa) = lambda_mature_plantation
            end if
            !------------------------------------------------------------------------------!

         case (3)
            !----- Primary forest. --------------------------------------------------------!
            if ( is_mature ) then
               lambda_harvest(ipa) = lambda_mature_primary
            end if
            !------------------------------------------------------------------------------!

         case (4:7)
            !------------------------------------------------------------------------------!
            !      Other types of forest.  We decide whether they are old_growth           !
            ! ("primary") or secondary depending on their age.                             !
            !------------------------------------------------------------------------------!
            if ( is_oldgrowth .and. is_mature ) then
               lambda_harvest(ipa) = lambda_mature_primary
            elseif ( is_mature ) then
               lambda_harvest(ipa) = lambda_mature_secondary
            end if
            !------------------------------------------------------------------------------!
         case default
            !----- Agriculture.  Do not log. ----------------------------------------------!
            continue
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end do patch_loop
      !------------------------------------------------------------------------------------!

      return
   end subroutine mature_forest_harvest
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine finds the area of mature patches that has to be harvested to meet  !
   ! the demand for biomass.                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine young_forest_harvest(cpoly,isi,onsp,harvestable_agb,lambda_harvest           &
                                  ,harvest_deficit)
      use ed_state_vars     , only : polygontype          & ! structure
                                   , sitetype             & ! structure
                                   , patchtype            ! ! structure
      use disturb_coms      , only : plantation_rotation  & ! intent(in)
                                   , mature_harvest_age   & ! intent(in)
                                   , min_harvest_biomass  & ! intent(in)
                                   , min_oldgrowth        ! ! intent(in)
      use consts_coms       , only : lnexp_max            ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)                  , target        :: cpoly
      integer                            , intent(in)    :: isi
      integer                            , intent(in)    :: onsp
      real, dimension(onsp)              , intent(in)    :: harvestable_agb
      real, dimension(onsp)              , intent(inout) :: lambda_harvest
      real                               , intent(inout) :: harvest_deficit
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)                     , pointer       :: csite
      integer                                            :: ipa
      integer                                            :: ilu
      logical                                            :: is_harvestable
      logical                                            :: is_plantation
      logical                                            :: is_secondary
      logical                                            :: is_primary
      logical                                            :: is_young
      !------------------------------------------------------------------------------------!



      !----- Select patch. ----------------------------------------------------------------!
      csite => cpoly%site(isi)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      First loop, try to obtain the biomass from young forest plantations.          !
      !------------------------------------------------------------------------------------!
      patch_loop_fopl: do ipa=1,onsp
         !----- Check whether we can harvest this patch. ----------------------------------!
         is_harvestable = harvestable_agb(ipa) >= min_harvest_biomass
         is_plantation  = csite%dist_type(ipa) == 2
         is_young       = csite%age(ipa)       <  plantation_rotation
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Harvest this patch if it qualifies.                                        !
         !---------------------------------------------------------------------------------!
         if (  is_harvestable  .and. is_plantation .and. is_young) then
            !----- Immature patch is harvestable.  Check how much to harvest. -------------!
            if( (csite%area(ipa) * harvestable_agb(ipa)) > harvest_deficit) then
               !---------------------------------------------------------------------------!
               !      Biomass target has been met, partially harvest the patch, then quit  !
               ! the sub-routine.                                                          !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = log(1./ (1. - harvest_deficit / harvestable_agb(ipa)))
               harvest_deficit     = 0.0
               return
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !      Biomass target has not been met, harvest the entire patch, and keep  !
               ! searching for biomass.                                                    !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = lnexp_max
               harvest_deficit     = harvest_deficit                                       &
                                   - csite%area(ipa) * harvestable_agb(ipa)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do patch_loop_fopl
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Second loop, try to obtain the biomass from secondary forests.                !
      !------------------------------------------------------------------------------------!
      patch_loop_2ary: do ipa=1,onsp
         !----- Check whether we can harvest this patch. ----------------------------------!
         ilu            = csite%dist_type(ipa)
         is_harvestable = harvestable_agb(ipa) >= min_harvest_biomass
         select case (ilu)
         case (4:7)
            is_secondary = csite%age(ipa) <  min_oldgrowth(ilu)
         case default
            is_secondary = .false.
         end select
         is_young = csite%age(ipa) <  mature_harvest_age
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Harvest this patch if it is has enough biomass, is a secondary forest, and !
         ! is young.                                                                       !
         !---------------------------------------------------------------------------------!
         if (  is_harvestable  .and. is_secondary .and. is_young) then
            !----- Immature patch is harvestable.  Check how much to harvest. -------------!
            if( (csite%area(ipa) * harvestable_agb(ipa)) > harvest_deficit) then
               !---------------------------------------------------------------------------!
               !      Biomass target has been met, partially harvest the patch, then quit  !
               ! the sub-routine.                                                          !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = log(1./ (1. - harvest_deficit / harvestable_agb(ipa)))
               harvest_deficit     = 0.0
               return
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !      Biomass target has not been met, harvest the entire patch, and keep  !
               ! searching for biomass.                                                    !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = lnexp_max
               harvest_deficit     = harvest_deficit                                       &
                                   - csite%area(ipa) * harvestable_agb(ipa)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do patch_loop_2ary
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Final loop, try to obtain the biomass from primary forests.                   !
      !------------------------------------------------------------------------------------!
      patch_loop_1ary: do ipa=1,onsp
         !----- Check whether we can harvest this patch. ----------------------------------!
         ilu            = csite%dist_type(ipa)
         is_harvestable = harvestable_agb(ipa) >= min_harvest_biomass
         select case (ilu)
         case (3)
            is_primary = .true.
         case (4:7)
            is_primary = csite%age(ipa) >= min_oldgrowth(ilu)
         case default
            is_primary = .false.
         end select
         is_young = csite%age(ipa) <  mature_harvest_age
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Harvest this patch if it is has enough biomass, is a primary forest, and   !
         ! is young.                                                                       !
         !---------------------------------------------------------------------------------!
         if (is_harvestable .and. is_primary .and. is_young) then
            !----- Immature patch is harvestable.  Check how much to harvest. -------------!
            if( (csite%area(ipa) * harvestable_agb(ipa)) > harvest_deficit) then
               !---------------------------------------------------------------------------!
               !      Biomass target has been met, partially harvest the patch, then quit  !
               ! the sub-routine.                                                          !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = log(1./ (1. - harvest_deficit / harvestable_agb(ipa)))
               harvest_deficit     = 0.0
               return
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !      Biomass target has not been met, harvest the entire patch, and keep  !
               ! searching for biomass.                                                    !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = lnexp_max
               harvest_deficit     = harvest_deficit                                       &
                                   - csite%area(ipa) * harvestable_agb(ipa)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do patch_loop_1ary
      !------------------------------------------------------------------------------------!

      return
   end subroutine young_forest_harvest
   !=======================================================================================!
   !=======================================================================================!
end module forestry
!==========================================================================================!
!==========================================================================================!
