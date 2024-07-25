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
!==========================================================================================!
module forestry



   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the disturbance rates associated with logging, when biomass !
   ! demands are provided instead of actual disturbance rates.                             !
   !---------------------------------------------------------------------------------------!
   subroutine find_lambda_harvest(cpoly,isi,onsp,lambda_harv_og_max,lambda_harvest)
      use ed_state_vars        , only : polygontype                & ! structure
                                      , sitetype                   & ! structure
                                      , patchtype                  & ! structure
                                      , allocate_sitetype          & ! subroutine
                                      , deallocate_sitetype        & ! subroutine
                                      , copy_sitetype              ! ! subroutine
      use disturb_coms         , only : ianth_disturb              & ! intent(in)
                                      , min_patch_area             & ! intent(in)
                                      , plantation_rotation        & ! intent(in)
                                      , min_harvest_biomass        & ! intent(in)
                                      , mature_harvest_age         & ! intent(in)
                                      , min_oldgrowth              ! ! intent(in)
      use fuse_fiss_utils      , only : terminate_patches          ! ! subroutine
      use ed_max_dims          , only : n_pft                      & ! intent(in)
                                      , n_dbh                      ! ! intent(in)
      use detailed_coms        , only : idetailed                  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)             , target        :: cpoly
      integer                       , intent(in)    :: isi
      integer                       , intent(in)    :: onsp
      real                          , intent(inout) :: lambda_harv_og_max
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
      real, dimension(onsp)                         :: pat_hvmax_btimber
      real, dimension(onsp)                         :: pat_hvpot_btimber
      real                                          :: primary_harvest_target
      real                                          :: secondary_harvest_target
      real                                          :: site_harvest_target
      real                                          :: site_hvmax_btimber
      real                                          :: site_hvpot_btimber
      real                                          :: area_mature_primary
      real                                          :: area_mature_secondary
      real                                          :: area_mature_plantation
      real                                          :: hvmax_mature_primary
      real                                          :: hvmax_mature_secondary
      real                                          :: hvmax_mature_plantation
      real                                          :: hvpot_mature_primary
      real                                          :: hvpot_mature_secondary
      real                                          :: hvpot_mature_plantation
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


      !----- Initialise book keeping variables. -------------------------------------------!
      pat_hvmax_btimber(:) = 0.0
      pat_hvpot_btimber(:) = 0.0
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
      ! tropics).  Also, we distinguish between the maximum harvestable biomass and the    !
      ! potential, because it may be required to leave behind a fraction of commercial     !
      ! trees that meet the logging size.                                                  !
      !------------------------------------------------------------------------------------!
      site_hvmax_btimber = 0.
      site_hvpot_btimber = 0.
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
         case (1,8)
            !---- Pasture or cropland. Do nothing. ----------------------------------------!
            continue
            !------------------------------------------------------------------------------!
         case (2:7)
            hcoh_loop: do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)
               if (cpatch%dbh(ico) >= cpoly%mindbh_harvest(ipft,isi)) then
                  !----- Cohort is harvestable. -------------------------------------------!
                  pat_hvmax_btimber(ipa) = pat_hvmax_btimber(ipa)                          &
                                         + cpoly%prob_harvest(ipft,isi)                    &
                                         * cpatch%nplant(ico) * cpatch%btimber(ico)
                  pat_hvpot_btimber(ipa) = pat_hvpot_btimber(ipa)                          &
                                         + cpatch%nplant(ico) * cpatch%btimber(ico)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do hcoh_loop
         end select
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Update site harvestable biomass only when the patch has sufficient         !
         ! harvestable biomass (after accounting for prob_harvest).                        !
         !---------------------------------------------------------------------------------!
         if (pat_hvmax_btimber(ipa) >= min_harvest_biomass) then
            site_hvmax_btimber = site_hvmax_btimber                                        &
                               + pat_hvmax_btimber(ipa) * csite%area(ipa)
            site_hvpot_btimber = site_hvpot_btimber                                        &
                               + pat_hvpot_btimber(ipa) * csite%area(ipa)
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
      if ( site_hvmax_btimber  <  min_harvest_biomass                 .or.                 &
           site_harvest_target <= site_hvmax_btimber * min_patch_area      ) then
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
            write (unit=*,fmt='(a,1x,i5)')     ' ISI                           = ',isi
            write (unit=*,fmt='(a,1x,es12.5)') ' PRIMARY TARGET BIOMASS        = '         &
                                                            , primary_harvest_target
            write (unit=*,fmt='(a,1x,es12.5)') ' SECONDARY TARGET BIOMASS      = '         &
                                                            , secondary_harvest_target
            write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL TARGET BIOMASS          = '         &
                                                            , site_harvest_target
            write (unit=*,fmt='(a,1x,es12.5)') ' MAXIMUM HARVESTABLE BIOMASS   = '         &
                                                            , site_hvmax_btimber
            write (unit=*,fmt='(a,1x,es12.5)') ' POTENTIAL HARVESTABLE BIOMASS = '         &
                                                            , site_hvpot_btimber
            write (unit=*,fmt='(a)'      )     '------------------------------------------'
            write (unit=*,fmt='(a)'      )     ' '
         end if
         !---------------------------------------------------------------------------------!


         return
      end if
      !------------------------------------------------------------------------------------!



      !------ Compute current stocks of timber in mature forests. -------------------------!
      call inventory_mature_forests(cpoly,isi,onsp,pat_hvmax_btimber,pat_hvpot_btimber     &
                                   ,area_mature_primary    ,hvmax_mature_primary           &
                                   ,hvpot_mature_primary   ,area_mature_secondary          &
                                   ,hvmax_mature_secondary ,hvpot_mature_secondary         &
                                   ,area_mature_plantation ,hvmax_mature_plantation        &
                                   ,hvpot_mature_plantation)
      !------------------------------------------------------------------------------------!



      !------ Compute the mature-forest harvest rates. ------------------------------------!
      call mature_forest_harvest(cpoly,isi,onsp                                            &
                                ,hvmax_mature_primary,hvpot_mature_primary                 &
                                ,hvmax_mature_secondary,hvpot_mature_secondary             &
                                ,hvmax_mature_plantation,hvpot_mature_plantation           &
                                ,primary_harvest_target,secondary_harvest_target           &
                                ,pat_hvmax_btimber,lambda_harv_og_max,lambda_harvest       &
                                ,harvest_deficit)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Compute the disturbance rate applied to young patches to help meet the        !
      ! biomass demands.                                                                   !
      !------------------------------------------------------------------------------------!
      call young_forest_harvest(cpoly,isi,onsp,pat_hvmax_btimber,pat_hvpot_btimber         &
                               ,lambda_harv_og_max,lambda_harvest,harvest_deficit)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Print the inventory and the target.                                          !
      !------------------------------------------------------------------------------------!
      if (print_detailed) then
         write (unit=*,fmt='(a)'      )     ' '
         write (unit=*,fmt='(a)'      )     '---------------------------------------------'
         write (unit=*,fmt='(a)'      )     ' FORESTRY.  HARVEST RATES'
         write (unit=*,fmt='(a)'      )     ' '
         write (unit=*,fmt='(a,1x,i5)')     ' ISI                           = ',isi
         write (unit=*,fmt='(a,1x,es12.5)') ' PRIMARY TARGET BIOMASS        = '            &
                                                         , primary_harvest_target
         write (unit=*,fmt='(a,1x,es12.5)') ' SECONDARY TARGET BIOMASS      = '            &
                                                         , secondary_harvest_target
         write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL TARGET BIOMASS          = '            &
                                                         , site_harvest_target
         write (unit=*,fmt='(a,1x,es12.5)') ' MAXIMUM HARVESTABLE BIOMASS   = '            &
                                                         , site_hvmax_btimber
         write (unit=*,fmt='(a,1x,es12.5)') ' POTENTIAL HARVESTABLE BIOMASS = '            &
                                                         , site_hvpot_btimber
         write (unit=*,fmt='(a,1x,es12.5)') ' HVMAX BIOMASS (PRIMARY)       = '            &
                                                         , hvmax_mature_primary
         write (unit=*,fmt='(a,1x,es12.5)') ' HVMAX BIOMASS (SECONDARY)     = '            &
                                                         , hvmax_mature_secondary
         write (unit=*,fmt='(a,1x,es12.5)') ' HVMAX BIOMASS (PLANTATION)    = '            &
                                                         , hvmax_mature_plantation
         write (unit=*,fmt='(a,1x,es12.5)') ' HVPOT BIOMASS (PRIMARY)       = '            &
                                                         , hvpot_mature_primary
         write (unit=*,fmt='(a,1x,es12.5)') ' HVPOT BIOMASS (SECONDARY)     = '            &
                                                         , hvpot_mature_secondary
         write (unit=*,fmt='(a,1x,es12.5)') ' HVPOT BIOMASS (PLANTATION)    = '            &
                                                         , hvpot_mature_plantation
         write (unit=*,fmt='(a,1x,es12.5)') ' HARVEST AREA (PRIMARY)        = '            &
                                                         , area_mature_primary
         write (unit=*,fmt='(a,1x,es12.5)') ' HARVEST AREA (SECONDARY)      = '            &
                                                         , area_mature_secondary
         write (unit=*,fmt='(a,1x,es12.5)') ' HARVEST AREA (PLANTATION)     = '            &
                                                         , area_mature_plantation
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a,1x,es12.5)') ' HARVEST DEFICIT               = '            &
                                                         , harvest_deficit
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a,1x,es12.5)') ' MAXIMUM LAMBDA (OLD-GROWTH)   = '            &
                                                         , lambda_harv_og_max
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a)'          ) '---------------------------------------------'
         write (unit=*,fmt='(9(a,1x))'    ) '  IPA','   LU','         AGE','        AREA'  &
                                            ,'   PAT_HVMAX','   PAT_HVPOT','      LAMBDA'  &
                                            ,'      MATURE','  OLD_GROWTH'
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
            write (unit=*,fmt='(2(i5,1x),5(f12.7,1x),2(11x,l1,1x))')                       &
               ipa,csite%dist_type(ipa),csite%age(ipa),csite%area(ipa)                     &
                  ,pat_hvmax_btimber(ipa),pat_hvpot_btimber(ipa),lambda_harvest(ipa)       &
                  ,is_mature,is_oldgrowth
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
   subroutine inventory_mature_forests( cpoly,isi,onsp,pat_hvmax_btimber,pat_hvpot_btimber &
                                      , area_mature_primary   , hvmax_mature_primary       &
                                      , hvpot_mature_primary  , area_mature_secondary      &
                                      , hvmax_mature_secondary, hvpot_mature_secondary     &
                                      , area_mature_plantation, hvmax_mature_plantation    &
                                      , hvpot_mature_plantation)
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
      real, dimension(onsp), intent(in)  :: pat_hvmax_btimber
      real, dimension(onsp), intent(in)  :: pat_hvpot_btimber
      real                 , intent(out) :: area_mature_primary
      real                 , intent(out) :: hvmax_mature_primary
      real                 , intent(out) :: hvpot_mature_primary
      real                 , intent(out) :: area_mature_secondary
      real                 , intent(out) :: hvmax_mature_secondary
      real                 , intent(out) :: hvpot_mature_secondary
      real                 , intent(out) :: area_mature_plantation
      real                 , intent(out) :: hvmax_mature_plantation
      real                 , intent(out) :: hvpot_mature_plantation
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer     :: csite
      integer                            :: ipa
      integer                            :: ilu
      logical                            :: is_oldgrowth
      logical                            :: is_mature
      logical                            :: is_rotation
      !------------------------------------------------------------------------------------!



      !----- Initialise inventory. --------------------------------------------------------!
      area_mature_primary     = 0.0
      area_mature_secondary   = 0.0
      area_mature_plantation  = 0.0
      hvmax_mature_primary    = 0.0
      hvpot_mature_primary    = 0.0
      hvmax_mature_secondary  = 0.0
      hvpot_mature_secondary  = 0.0
      hvmax_mature_plantation = 0.0
      hvpot_mature_plantation = 0.0
      !------------------------------------------------------------------------------------!

      csite => cpoly%site(isi)
      patchloop: do ipa=1,csite%npatches

         !----- Skip the patch if the biomass is low. -------------------------------------!
         if (pat_hvmax_btimber(ipa) < min_harvest_biomass) cycle patchloop
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
               hvmax_mature_plantation = hvmax_mature_plantation                           &
                                       + pat_hvmax_btimber(ipa)   * csite%area(ipa)
               hvpot_mature_plantation = hvpot_mature_plantation                           &
                                       + pat_hvpot_btimber(ipa)   * csite%area(ipa)
            end if
            !------------------------------------------------------------------------------!

         case (3)
            !---- Treefall always goes to "primary" forest. -------------------------------!
            if (is_mature) then
               area_mature_primary     = area_mature_primary    + csite%area(ipa)
               hvmax_mature_primary    = hvmax_mature_primary                              &
                                       + pat_hvmax_btimber(ipa)   * csite%area(ipa)
               hvpot_mature_primary    = hvpot_mature_primary                              &
                                       + pat_hvpot_btimber(ipa)   * csite%area(ipa)
            end if
            !------------------------------------------------------------------------------!

         case (4:7)
            !------------------------------------------------------------------------------!
            !     Other disturbances, assume secondary in case the patch is not considered !
            ! old-growth, and primary otherwise.                                           !
            !------------------------------------------------------------------------------!
            if (is_mature .and. is_oldgrowth) then
               area_mature_primary     = area_mature_primary    + csite%area(ipa)
               hvmax_mature_primary    = hvmax_mature_primary                              &
                                       + pat_hvmax_btimber(ipa)   * csite%area(ipa)
               hvpot_mature_primary    = hvpot_mature_primary                              &
                                       + pat_hvpot_btimber(ipa)   * csite%area(ipa)
            else if (is_mature) then
               area_mature_secondary   = area_mature_secondary  + csite%area(ipa)
               hvmax_mature_secondary  = hvmax_mature_secondary                            &
                                       + pat_hvmax_btimber(ipa)   * csite%area(ipa)
               hvpot_mature_secondary  = hvpot_mature_secondary                            &
                                       + pat_hvpot_btimber(ipa)   * csite%area(ipa)
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
   subroutine mature_forest_harvest(cpoly,isi,onsp                                         &
                                   ,hvmax_mature_primary,hvpot_mature_primary              &
                                   ,hvmax_mature_secondary,hvpot_mature_secondary          &
                                   ,hvmax_mature_plantation,hvpot_mature_plantation        &
                                   ,primary_harvest_target,secondary_harvest_target        &
                                   ,pat_hvmax_btimber,lambda_harv_og_max,lambda_harvest    &
                                   ,harvest_deficit)
      use ed_state_vars , only : polygontype         & ! structure
                               , sitetype            ! ! structure
      use disturb_coms  , only : plantation_rotation & ! intent(in)
                               , mature_harvest_age  & ! intent(in)
                               , min_harvest_biomass & ! intent(in)
                               , min_oldgrowth       ! ! intent(in)
      use consts_coms   , only : lnexp_max           & ! intent(in)
                               , almost_one          & ! intent(in)
                               , tiny_num            ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)    , target        :: cpoly
      integer              , intent(in)    :: isi
      integer              , intent(in)    :: onsp
      real                 , intent(in)    :: hvmax_mature_primary
      real                 , intent(in)    :: hvpot_mature_primary
      real                 , intent(in)    :: hvmax_mature_secondary
      real                 , intent(in)    :: hvpot_mature_secondary
      real                 , intent(in)    :: hvmax_mature_plantation
      real                 , intent(in)    :: hvpot_mature_plantation
      real                 , intent(in)    :: primary_harvest_target
      real                 , intent(inout) :: secondary_harvest_target
      real, dimension(onsp), intent(in)    :: pat_hvmax_btimber
      real                 , intent(in)    :: lambda_harv_og_max
      real, dimension(onsp), intent(inout) :: lambda_harvest
      real                 , intent(out)   :: harvest_deficit
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer       :: csite
      integer                              :: ipa
      integer                              :: ilu
      logical                              :: is_oldgrowth
      logical                              :: is_mature
      logical                              :: is_rotation
      real                                 :: f_harvest
      real                                 :: f_harv_prim_max
      real                                 :: harvest_actual
      real                                 :: lambda_mature_primary
      real                                 :: lambda_mature_plantation
      real                                 :: lambda_mature_secondary
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Make sure the disturbance rate for primary forests doesn't exceed the maximum  !
      ! disturbance taken from the land use instructions for this year.                    !
      !------------------------------------------------------------------------------------!
      f_harv_prim_max = 1. - exp (-lambda_harv_og_max)
      !------------------------------------------------------------------------------------!


      !---- Find harvesting rate in mature primary vegetation. ----------------------------!
      if ( hvpot_mature_primary > tiny_num ) then
         !---------------------------------------------------------------------------------!
         !    Find the total biomass to be harvested from mature primary vegetation.  The  !
         ! actual harvest cannot exceed the target, cannot exceed the actual biomass, and  !
         ! must be bounded by the disturbance rate target.                                 !
         !---------------------------------------------------------------------------------!
         harvest_actual = min( f_harv_prim_max * hvpot_mature_primary                      &
                             , hvmax_mature_primary, primary_harvest_target )
         f_harvest      = harvest_actual / hvpot_mature_primary
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Check fraction to be harvest, to avoid floating point exceptions.            !
         !---------------------------------------------------------------------------------!
         if (f_harvest > almost_one) then
            !---- Apply sweeping logging disturbance across primary forests. --------------!
            f_harvest                = 1.0
            harvest_actual           = hvpot_mature_primary
            lambda_mature_primary    = lnexp_max
            harvest_deficit          = primary_harvest_target   - harvest_actual
            secondary_harvest_target = secondary_harvest_target + harvest_deficit
            !------------------------------------------------------------------------------!
         else
            !---- Find the disturbance rate from the fraction of biomass to be harvested. -!
            lambda_mature_primary    = log(1./(1.-f_harvest))
            harvest_deficit          = primary_harvest_target   - harvest_actual
            secondary_harvest_target = secondary_harvest_target + harvest_deficit
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     No mature primary forest to harvest. Do not log from any primary forest and !
         ! try to harvest secondary forest instead.                                        !
         !---------------------------------------------------------------------------------!
         f_harvest                = 0.0
         harvest_actual           = 0.0
         lambda_mature_primary    = 0.0
         secondary_harvest_target = secondary_harvest_target + primary_harvest_target
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find harvest disturbance rates for mature forest plantations.  We always give !
      ! preference for logging plantations over second-growth forests.                     !
      !------------------------------------------------------------------------------------!
      if ( hvpot_mature_plantation > tiny_num ) then
         !---------------------------------------------------------------------------------!
         !    Find the total biomass to be harvested from mature forest plantation.  The   !
         ! actual harvest cannot exceed the target or the actual biomass.                  !
         !---------------------------------------------------------------------------------!
         harvest_actual = min( secondary_harvest_target, hvmax_mature_plantation           &
                             , hvpot_mature_plantation                           )
         f_harvest      = harvest_actual / hvpot_mature_plantation
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Check fraction to be harvest, to avoid floating point exceptions.            !
         !---------------------------------------------------------------------------------!
         if (f_harvest > almost_one) then
            !---- Apply sweeping logging disturbance across primary forests. --------------!
            f_harvest                = 1.0
            harvest_actual           = hvpot_mature_plantation
            lambda_mature_plantation = lnexp_max
            harvest_deficit          = max(0.,secondary_harvest_target - harvest_actual)
            !------------------------------------------------------------------------------!
         else
            !---- Find the disturbance rate from the fraction of biomass to be harvested. -!
            lambda_mature_plantation = log(1./(1.-f_harvest))
            harvest_deficit          = max(0.,secondary_harvest_target - harvest_actual)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     No forest plantation to harvest. Do not log from any plantation and try to  !
         ! harvest secondary forest instead.                                               !
         !---------------------------------------------------------------------------------!
         f_harvest                = 0.0
         harvest_actual           = 0.0
         lambda_mature_plantation = 0.0
         harvest_deficit          = secondary_harvest_target
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find disturbance rates for secondary forests other than forest plantation.    !
      !------------------------------------------------------------------------------------!
      if ( hvpot_mature_secondary > tiny_num ) then
         !---------------------------------------------------------------------------------!
         !    Find the total biomass to be harvested from mature second-growth forest.     !
         ! The actual harvest cannot exceed the target or the actual biomass.              !
         !---------------------------------------------------------------------------------!
         harvest_actual = min( harvest_deficit, hvmax_mature_secondary                     &
                             , hvpot_mature_secondary                  )
         f_harvest      = harvest_actual / hvpot_mature_secondary
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Check fraction to be harvest, to avoid floating point exceptions.            !
         !---------------------------------------------------------------------------------!
         if (f_harvest > almost_one) then
            !---- Apply sweeping logging disturbance across primary forests. --------------!
            f_harvest               = 1.0
            harvest_actual          = hvpot_mature_secondary
            lambda_mature_secondary = lnexp_max
            harvest_deficit         = max(0.,harvest_deficit - harvest_actual)
            !------------------------------------------------------------------------------!
         else
            !---- Find the disturbance rate from the fraction of biomass to be harvested. -!
            lambda_mature_secondary = log(1./(1.-f_harvest))
            harvest_deficit         = max(0.,harvest_deficit - harvest_actual)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     No mature second-growth forest to harvest. Do not log from any plantation   !
         ! and try to harvest secondary forest instead.                                    !
         !---------------------------------------------------------------------------------!
         f_harvest               = 0.0
         harvest_actual          = 0.0
         lambda_mature_secondary = 0.0
         ! harvest_deficit       = harvest_deficit
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over patches.                                                            ! 
      !------------------------------------------------------------------------------------!
      csite => cpoly%site(isi)
      patch_loop: do ipa=1,onsp
     
         !----- Skip patch in case the biomass is less than the minimum for harvesting. ---!
         if (pat_hvmax_btimber(ipa) < min_harvest_biomass) cycle patch_loop
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
   subroutine young_forest_harvest(cpoly,isi,onsp,pat_hvmax_btimber,pat_hvpot_btimber      &
                                  ,lambda_harv_og_max,lambda_harvest,harvest_deficit)
      use ed_state_vars     , only : polygontype          & ! structure
                                   , sitetype             & ! structure
                                   , patchtype            ! ! structure
      use disturb_coms      , only : plantation_rotation  & ! intent(in)
                                   , mature_harvest_age   & ! intent(in)
                                   , min_harvest_biomass  & ! intent(in)
                                   , min_oldgrowth        ! ! intent(in)
      use consts_coms       , only : lnexp_max            & ! intent(in)
                                   , almost_one           & ! intent(in)
                                   , tiny_num             ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)                  , target        :: cpoly
      integer                            , intent(in)    :: isi
      integer                            , intent(in)    :: onsp
      real, dimension(onsp)              , intent(in)    :: pat_hvmax_btimber
      real, dimension(onsp)              , intent(in)    :: pat_hvpot_btimber
      real, dimension(onsp)              , intent(inout) :: lambda_harvest
      real                               , intent(in)    :: lambda_harv_og_max
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
      real                                               :: f_harvest
      real                                               :: f_harv_prim_max
      real                                               :: harvest_actual
      real                                               :: site_hvmax_btimber
      real                                               :: site_hvpot_btimber
      !------------------------------------------------------------------------------------!



      !----- Select patch. ----------------------------------------------------------------!
      csite => cpoly%site(isi)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Make sure the disturbance rate for primary forests doesn't exceed the maximum  !
      ! disturbance taken from the land use instructions for this year.                    !
      !------------------------------------------------------------------------------------!
      f_harv_prim_max = 1. - exp (-lambda_harv_og_max)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      First loop, try to obtain the biomass from young forest plantations.          !
      !------------------------------------------------------------------------------------!
      patch_loop_fypl: do ipa=1,onsp
         !----- Check whether we can harvest this patch. ----------------------------------!
         is_harvestable = pat_hvmax_btimber(ipa) >= min_harvest_biomass
         is_plantation  = csite%dist_type(ipa)   == 2
         is_young       = csite%age(ipa)         <  plantation_rotation
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Harvest this patch if it qualifies.                                        !
         !---------------------------------------------------------------------------------!
         if (  is_harvestable  .and. is_plantation .and. is_young) then
            !----- Find the site-level harvestable biomass. -------------------------------!
            site_hvmax_btimber = pat_hvmax_btimber(ipa) * csite%area(ipa)
            site_hvpot_btimber = pat_hvpot_btimber(ipa) * csite%area(ipa)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check if there is enough biomass to harvest in this patch.               !
            !------------------------------------------------------------------------------!
            if ( site_hvpot_btimber > tiny_num ) then
               !---------------------------------------------------------------------------!
               !     Check how much biomass can be extracted from this patch, up to the    !
               ! target biomass.                                                           !
               !---------------------------------------------------------------------------!
               harvest_actual = min(harvest_deficit,site_hvmax_btimber,site_hvpot_btimber)
               f_harvest      = harvest_actual / site_hvpot_btimber
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Check fraction to be harvest, to avoid floating point exceptions.      !
               !---------------------------------------------------------------------------!
               if (f_harvest > almost_one) then
                  !---- Apply sweeping logging disturbance across primary forests. --------!
                  f_harvest            = 1.0
                  harvest_actual       = site_hvpot_btimber
                  lambda_harvest(ipa)  = lnexp_max
                  harvest_deficit      = max(0.,harvest_deficit - harvest_actual)
                  !------------------------------------------------------------------------!
               else
                  !---- Find the disturbance rate. ----------------------------------------!
                  lambda_harvest(ipa)  = log(1./(1.-f_harvest))
                  harvest_deficit      = max(0.,harvest_deficit - harvest_actual)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

               !---- If we successfully met harvest demands, quit routine. ----------------!
               if ( harvest_deficit == 0. ) return
               !---------------------------------------------------------------------------!
            else
               !----- Patch does not have enough biomass. ---------------------------------!
               f_harvest         = 0.0
               harvest_actual    = 0.0
               ! harvest_deficit = harvest_deficit
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do patch_loop_fypl
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Second loop, try to obtain the biomass from secondary forests.                !
      !------------------------------------------------------------------------------------!
      patch_loop_2ary: do ipa=1,onsp
         !----- Check whether we can harvest this patch. ----------------------------------!
         ilu            = csite%dist_type(ipa)
         is_harvestable = pat_hvmax_btimber(ipa) >= min_harvest_biomass
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
            !----- Find the site-level harvestable biomass. -------------------------------!
            site_hvmax_btimber = pat_hvmax_btimber(ipa) * csite%area(ipa)
            site_hvpot_btimber = pat_hvpot_btimber(ipa) * csite%area(ipa)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check if there is enough biomass to harvest in this patch.               !
            !------------------------------------------------------------------------------!
            if ( site_hvpot_btimber > tiny_num ) then
               !---------------------------------------------------------------------------!
               !     Check how much biomass can be extracted from this patch, up to the    !
               ! target biomass.                                                           !
               !---------------------------------------------------------------------------!
               harvest_actual = min(harvest_deficit,site_hvmax_btimber,site_hvpot_btimber)
               f_harvest      = harvest_actual / site_hvpot_btimber
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Check fraction to be harvest, to avoid floating point exceptions.      !
               !---------------------------------------------------------------------------!
               if (f_harvest > almost_one) then
                  !---- Apply sweeping logging disturbance across primary forests. --------!
                  f_harvest            = 1.0
                  harvest_actual       = site_hvpot_btimber
                  lambda_harvest(ipa)  = lnexp_max
                  harvest_deficit      = max(0.,harvest_deficit - harvest_actual)
                  !------------------------------------------------------------------------!
               else
                  !---- Find the disturbance rate. ----------------------------------------!
                  lambda_harvest(ipa)  = log(1./(1.-f_harvest))
                  harvest_deficit      = max(0.,harvest_deficit - harvest_actual)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

               !---- If we successfully met harvest demands, quit routine. ----------------!
               if ( harvest_deficit == 0. ) return
               !---------------------------------------------------------------------------!
            else
               !----- Patch does not have enough biomass. ---------------------------------!
               f_harvest         = 0.0
               harvest_actual    = 0.0
               ! harvest_deficit = harvest_deficit
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
         is_harvestable = pat_hvmax_btimber(ipa) >= min_harvest_biomass
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
            !----- Find the site-level harvestable biomass. -------------------------------!
            site_hvmax_btimber = pat_hvmax_btimber(ipa) * csite%area(ipa)
            site_hvpot_btimber = pat_hvpot_btimber(ipa) * csite%area(ipa)
            !------------------------------------------------------------------------------!





            !------------------------------------------------------------------------------!
            !     Check if there is enough biomass to harvest in this patch.               !
            !------------------------------------------------------------------------------!
            if ( site_hvpot_btimber > tiny_num ) then
               !---------------------------------------------------------------------------!
               !     Check how much biomass can be extracted from this patch, up to the    !
               ! target biomass.                                                           !
               !---------------------------------------------------------------------------!
               harvest_actual = min( harvest_deficit, site_hvmax_btimber                   &
                                   , f_harv_prim_max * site_hvpot_btimber )
               f_harvest      = harvest_actual / site_hvpot_btimber
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Check fraction to be harvest, to avoid floating point exceptions.      !
               !---------------------------------------------------------------------------!
               if (f_harvest > almost_one) then
                  !---- Apply sweeping logging disturbance across primary forests. --------!
                  f_harvest            = 1.0
                  harvest_actual       = site_hvpot_btimber
                  lambda_harvest(ipa)  = lnexp_max
                  harvest_deficit      = max(0.,harvest_deficit - harvest_actual)
                  !------------------------------------------------------------------------!
               else
                  !---- Find the disturbance rate. ----------------------------------------!
                  lambda_harvest(ipa)  = log(1./(1.-f_harvest))
                  harvest_deficit      = max(0.,harvest_deficit - harvest_actual)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

               !---- If we successfully met harvest demands, quit routine. ----------------!
               if ( harvest_deficit == 0. ) return
               !---------------------------------------------------------------------------!
            else
               !----- Patch does not have enough biomass. ---------------------------------!
               f_harvest         = 0.0
               harvest_actual    = 0.0
               ! harvest_deficit = harvest_deficit
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!




            !----- Immature patch is harvestable.  Check how much to harvest. -------------!
            if (almost_one * site_hvmax_btimber > harvest_deficit) then
               !---------------------------------------------------------------------------!
               !      Biomass target has been met, harvest the patch, then quit the sub-   !
               ! -routine.                                                                 !
               !---------------------------------------------------------------------------!
               f_harvest           = harvest_deficit / site_hvpot_btimber
               lambda_harvest(ipa) = log(1./ (1. - f_harvest))
               harvest_deficit     = 0.0
               return
               !---------------------------------------------------------------------------!
            else if (almost_one * site_hvpot_btimber > site_hvmax_btimber) then
               !---------------------------------------------------------------------------!
               !      Biomass target has not been met, harvest the entire patch, and keep  !
               ! searching for biomass.                                                    !
               !---------------------------------------------------------------------------!
               f_harvest           = site_hvmax_btimber / site_hvpot_btimber
               lambda_harvest(ipa) = log(1./ (1. - f_harvest))
               harvest_deficit     = harvest_deficit - site_hvmax_btimber
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !      Biomass target has not been met, harvest the entire patch, and keep  !
               ! searching for biomass.                                                    !
               !---------------------------------------------------------------------------!
               lambda_harvest(ipa) = lnexp_max
               harvest_deficit     = harvest_deficit - site_hvmax_btimber
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
