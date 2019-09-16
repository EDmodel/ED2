module forestry
  contains

!==========================================================================================!
! Forestry.f90. These subroutines calculate the area to be taken from each patch to meet   !
!               the biomass demand, in case this site is to be harvest using biomass       !
!               targets.                                                                   !                                                                                          !
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
subroutine find_harvest_area(cpoly,isi,onsp,harvestable_agb,pot_area_harv)
   use ed_state_vars        , only : polygontype                & ! structure
                                   , sitetype                   & ! structure
                                   , patchtype                  & ! structure
                                   , allocate_sitetype          & ! subroutine
                                   , deallocate_sitetype        & ! subroutine
                                   , copy_sitetype              ! ! subroutine
   use disturb_coms         , only : ianth_disturb              & ! intent(in)
                                   , lutime                     & ! intent(in)
                                   , min_patch_area             & ! intent(in)
                                   , min_harvest_biomass        & ! intent(in)
                                   , plantation_rotation        & ! intent(in)
                                   , mature_harvest_age         ! ! intent(in)
!                                    , min_oldgrowth              ! ! intent(in) ! New; need to add
   use fuse_fiss_utils      , only : terminate_patches          ! ! subroutine
   use ed_max_dims          , only : n_pft                      & ! intent(in)
                                   , n_dbh                      ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)             , target        :: cpoly
   integer                       , intent(in)    :: isi
   integer                       , intent(in)    :: onsp
   real, dimension(onsp)         , intent(inout) :: harvestable_agb
   real, dimension(onsp)         , intent(inout) :: pot_area_harv
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)                , pointer       :: csite
   type(patchtype)               , pointer       :: cpatch
   integer                                       :: ipa
   integer                                       :: ico
   integer                                       :: ipft
   integer                                       :: ilu
   real, dimension(n_pft)                        :: mindbh_harvest
   real                                          :: primary_harvest_target
   real                                          :: secondary_harvest_target
   real                                          :: site_harvest_target
   real                                          :: site_harvestable_agb
!       real                                          :: site_hvmax_btimber
!       real                                          :: site_hvpot_btimber   
   real                                          :: area_mature_primary
   real                                          :: hvagb_mature_primary
   real                                          :: area_mature_secondary
   real                                          :: hvagb_mature_secondary
   real                                          :: area_mature_plantation
   real                                          :: hvagb_mature_plantation
   real                                          :: lambda_mature_primary
   real                                          :: lambda_mature_secondary
   real                                          :: lambda_mature_plantation
   real                                          :: harvest_deficit
!    real, dimension(onsp)                         :: pat_hvmax_btimber
!    real, dimension(onsp)                         :: pat_hvpot_btimber

   !----- Local constants. ----------------------------------------------------------------!
   logical                       , parameter     :: print_detail = .false.
   !---------------------------------------------------------------------------------------!

   select case (ianth_disturb)
   !----- Nothing to do because not basing harvest off of biomass -------------------------!
   case (0) ! Nothing to do because anthropogenic disturbance is turned off
	   return
   case (2) ! Harvest is based on size + area

	   !----- Link to the current site. -------------------------------------------------------!
	   csite => cpoly%site(isi)
	   !---------------------------------------------------------------------------------------!

	   cpoly%primary_harvest_memory  (isi) = 0.0
	   cpoly%secondary_harvest_memory(isi) = 0.0
	   
	   lambda_mature_plantation = cpoly%disturbance_rates(2,2,isi)
	   lambda_mature_primary = cpoly%disturbance_rates(6,6,isi)

	   !---------------------------------------------------------------------------------------!
       !      Loop over patches.                                                               ! 
       !---------------------------------------------------------------------------------------!
       patch_loop: do ipa=1,onsp
  
		  !------------------------------------------------------------------------------------!
		  !    Find out whether to harvest this patch.                                         !
		  !------------------------------------------------------------------------------------!
		  select case (csite%dist_type(ipa))
		  case (2)
			 !----- Forest plantation. --------------------------------------------------------!
			 if ( csite%age(ipa) > plantation_rotation ) then
				pot_area_harv(ipa) = csite%area(ipa) * lambda_mature_plantation
			 end if
			 !---------------------------------------------------------------------------------!
		  case (6)
			 !----- Primary/Secondary forest. -------------------------------------------------!
			 if ( csite%age(ipa) > mature_harvest_age  ) then
				pot_area_harv(ipa) = csite%area(ipa) * cpoly%disturbance_rates(6,6,isi)
			 end if
		  case default
			 !----- Agriculture.  Do not log. -------------------------------------------------!
			 continue
			 !---------------------------------------------------------------------------------!
		  end select
		  !------------------------------------------------------------------------------------!
	   end do patch_loop
	   !---------------------------------------------------------------------------------------!

	  write (unit=*,fmt='(a)'      )     ' '
	  write (unit=*,fmt='(a)'      )     '------------------------------------------------'
	  write (unit=*,fmt='(a)'      )     ' FORESTRY.  HARVEST RATES'
	  write (unit=*,fmt='(a)'      )     ' '
	  write (unit=*,fmt='(a,1x,i5)')     ' ISI                    = ',isi
	  write (unit=*,fmt='(a,1x,es12.5)') ' HV LAMBDA (PRIMARY)    = '                      &
													  , lambda_mature_primary
	  write (unit=*,fmt='(a,1x,es12.5)') ' HV LAMBDA (PLANTATION) = '                      &
													  , lambda_mature_plantation
	  write (unit=*,fmt='(a)'          ) ' '
	  write (unit=*,fmt='(a)'          ) '------------------------------------------------'
	  write (unit=*,fmt='(5(a,1x))'    ) '  IPA','   LU','         AGE','        AREA'     &
														,'     HV_AREA'
	  write (unit=*,fmt='(a)'      )     '------------------------------------------------'
	  do ipa=1,onsp
		 write (unit=*,fmt='(2(i5,1x),3(f12.7,1x))') ipa,csite%dist_type(ipa)              &
												   ,csite%age(ipa),csite%area(ipa)         &
												   ,pot_area_harv(ipa)
	  end do
	  write (unit=*,fmt='(a)'      )     '------------------------------------------------'
	  write (unit=*,fmt='(a)'      )     ' '

   !---------------------------------------------------------------------------------------!

   case (1)
	   !----- Link to the current site. -------------------------------------------------------!
	   csite => cpoly%site(isi)
	   !---------------------------------------------------------------------------------------!



	   !---------------------------------------------------------------------------------------!
	   !      Set biomass targets based on current rates and unapplied harvest from previous   !
	   ! years (memory).  These are in kgC/m2.  In case DBH-based logging is applied, then     !
	   ! harvest target is set to zero.                                                        !
	   !---------------------------------------------------------------------------------------!
	   primary_harvest_target   = cpoly%primary_harvest_target  (isi)                          &
								+ cpoly%primary_harvest_memory  (isi)
	   secondary_harvest_target = cpoly%secondary_harvest_target(isi)                          &
								+ cpoly%primary_harvest_memory  (isi)
	   site_harvest_target      = primary_harvest_target + secondary_harvest_target
	   !---------------------------------------------------------------------------------------!


	   !---------------------------------------------------------------------------------------!
	   !      Find total harvestable biomass density in kgC/m2.  The harvestable biomass is    !
	   ! not necessarily the same as the total site biomass, because it is possible that the   !
	   ! demand is for trees above a minimum size (which is very common practice in the        !
	   ! tropics).                                                                             !
	   !---------------------------------------------------------------------------------------!
	   site_harvestable_agb = 0.
	   hpat_loop: do ipa=1,onsp
		  cpatch => csite%patch(ipa)
		  ilu = csite%dist_type(ipa)

		  !----- Select the minimum DBH depending on the forest category. ---------------------!
		  select case(ilu)
		  case (1)
			 !----- Agriculture.  No timber harvesting here. ----------------------------------!
			 mindbh_harvest(:) = huge(1.)
			 !---------------------------------------------------------------------------------!
		  case (2:6)
			 !----- Harvesting of some sort.  ALl use same variable now -----------------------!
			 mindbh_harvest(:) = cpoly%mindbh_harvest(:,isi)
			 !---------------------------------------------------------------------------------!
		  end select
		  !------------------------------------------------------------------------------------!

		  !------------------------------------------------------------------------------------!
		  !     Go over each cohort, seek harvestable biomass.                                 !
		  !------------------------------------------------------------------------------------!
		  hcoh_loop: do ico=1,cpatch%ncohorts
			 ipft = cpatch%pft(ico)
			 if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
				!----- Cohort is harvestable. -------------------------------------------------!
				harvestable_agb(ipa) = harvestable_agb(ipa)                                    &
									 + cpatch%nplant(ico) * cpatch%agb(ico)
				!------------------------------------------------------------------------------!
			 end if
			 !---------------------------------------------------------------------------------!
		  end do hcoh_loop
		  !------------------------------------------------------------------------------------!

		  !----- Update site biomass only when the patch has sufficient biomass. --------------!
		  if (harvestable_agb(ipa) >= min_harvest_biomass) then
			 site_harvestable_agb = site_harvestable_agb                                       &
								  + harvestable_agb(ipa) * csite%area(ipa)
		  end if
		  !------------------------------------------------------------------------------------!
	   end do hpat_loop
	   !---------------------------------------------------------------------------------------!



	   !=======================================================================================!
	   !=======================================================================================!
	   !     Find out whether any harvest can occur at this site.  For harvest to occur, the   !
	   ! site must meet two criteria:                                                          !
	   ! a. Harvestable biomass must be greater than min_harvest_biomass;                      !
	   ! b. Target biomass must be greater than min_harvest_biomass.                           !
	   !                                                                                       !
	   !    In case one or both criteria are not met, we add the harvest target to memory.     !
	   !---------------------------------------------------------------------------------------!
	   if ( site_harvestable_agb <  min_harvest_biomass                 .or.                   &
			site_harvest_target  <= site_harvestable_agb * min_patch_area      ) then
		  cpoly%primary_harvest_target  (isi) = 0.0
		  cpoly%secondary_harvest_target(isi) = 0.0
		  cpoly%primary_harvest_memory  (isi) = primary_harvest_target
		  cpoly%secondary_harvest_memory(isi) = secondary_harvest_target
		  !------------------------------------------------------------------------------------!


		  !------------------------------------------------------------------------------------!
		  !       Print the inventory and the target.                                          !
		  !------------------------------------------------------------------------------------!
		  if (print_detail) then
			 write (unit=*,fmt='(a)'      )     ' '
			 write (unit=*,fmt='(a)'      )     '---------------------------------------------'
			 write (unit=*,fmt='(a)'      )     ' FORESTRY.  HARVEST WON''T OCCUR THIS YEAR'
			 write (unit=*,fmt='(a)'      )     ' '
			 write (unit=*,fmt='(a,1x,i5)')     ' ISI                   = ',isi
			 write (unit=*,fmt='(a,1x,es12.5)') ' PRIMARY TARGET AGB    = '                    &
															 , primary_harvest_target
			 write (unit=*,fmt='(a,1x,es12.5)') ' SECONDARY TARGET AGB  = '                    &
															 , secondary_harvest_target
			 write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL TARGET AGB      = '                    &
															 , site_harvest_target
			 write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL HARVESTABLE AGB = '                    &
															 , site_harvestable_agb
			 write (unit=*,fmt='(a)'      )     '---------------------------------------------'
			 write (unit=*,fmt='(a)'      )     ' '
		  end if
		  !------------------------------------------------------------------------------------!
		  return
	   end if
	   !---------------------------------------------------------------------------------------!



	   !------ Compute current stocks of agb in mature forests. -------------------------------!
	   call inventory_mat_forests(cpoly,isi,onsp,harvestable_agb                               &
								 ,area_mature_primary   ,hvagb_mature_primary                  &
								 ,area_mature_secondary ,hvagb_mature_secondary                &
								 ,area_mature_plantation,hvagb_mature_plantation)
	   !---------------------------------------------------------------------------------------!



	   !------ Compute the mature-forest harvest rates. ---------------------------------------!
	   call mat_forest_harv_rates(hvagb_mature_primary,hvagb_mature_secondary                  &
								 ,hvagb_mature_plantation,primary_harvest_target               &
								 ,secondary_harvest_target,lambda_mature_primary               &
								 ,lambda_mature_secondary,lambda_mature_plantation             &
								 ,harvest_deficit)                                    
	   !---------------------------------------------------------------------------------------!



	   !---------------------------------------------------------------------------------------!
	   !      Compute the area lost by mature patches through harvest.                         !
	   !---------------------------------------------------------------------------------------!
	   call area_harvest_mature(cpoly,isi,onsp,harvestable_agb,pot_area_harv                   &
							   ,lambda_mature_primary,lambda_mature_secondary                  &
							   ,lambda_mature_plantation)
	   !---------------------------------------------------------------------------------------!



	   !---------------------------------------------------------------------------------------!
	   !      Compute the area lost by immature patches to meet the biomass demand.            !
	   !---------------------------------------------------------------------------------------!
	   call area_harvest_immature(cpoly,isi,onsp,harvestable_agb,pot_area_harv,harvest_deficit)
	   !---------------------------------------------------------------------------------------!



	   !---------------------------------------------------------------------------------------!
	   !       Print the inventory and the target.                                             !
	   !---------------------------------------------------------------------------------------!
	   if (print_detail) then
		  write (unit=*,fmt='(a)'      )     ' '
		  write (unit=*,fmt='(a)'      )     '------------------------------------------------'
		  write (unit=*,fmt='(a)'      )     ' FORESTRY.  HARVEST RATES'
		  write (unit=*,fmt='(a)'      )     ' '
		  write (unit=*,fmt='(a,1x,i5)')     ' ISI                    = ',isi
		  write (unit=*,fmt='(a,1x,es12.5)') ' PRIMARY TARGET AGB     = '                      &
														  , primary_harvest_target
		  write (unit=*,fmt='(a,1x,es12.5)') ' SECONDARY TARGET AGB   = '                      &
														  , secondary_harvest_target
		  write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL TARGET AGB       = '                      &
														  , site_harvest_target
		  write (unit=*,fmt='(a,1x,es12.5)') ' TOTAL HARVESTABLE AGB  = '                      &
														  , site_harvestable_agb
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV AGB (PRIMARY)       = '                      &
														  , hvagb_mature_primary
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV AGB (SECONDARY)     = '                      &
														  , hvagb_mature_secondary
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV AGB (PLANTATION)    = '                      &
														  , hvagb_mature_plantation
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV AREA (PRIMARY)      = '                      &
														  , area_mature_primary
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV AREA (SECONDARY)    = '                      &
														  , area_mature_secondary
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV AREA (PLANTATION)   = '                      &
														  , area_mature_plantation
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV LAMBDA (PRIMARY)    = '                      &
														  , lambda_mature_primary
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV LAMBDA (SECONDARY)  = '                      &
														  , lambda_mature_secondary
		  write (unit=*,fmt='(a,1x,es12.5)') ' HV LAMBDA (PLANTATION) = '                      &
														  , lambda_mature_plantation
		  write (unit=*,fmt='(a)'          ) ' '
		  write (unit=*,fmt='(a,1x,es12.5)') ' HARVEST DEFICIT = ', harvest_deficit
		  write (unit=*,fmt='(a)'          ) ' '
		  write (unit=*,fmt='(a)'          ) '------------------------------------------------'
		  write (unit=*,fmt='(5(a,1x))'    ) '  IPA','   LU','         AGE','        AREA'     &
															,'     HV_AREA'
		  write (unit=*,fmt='(a)'      )     '------------------------------------------------'
		  do ipa=1,onsp
			 write (unit=*,fmt='(2(i5,1x),3(f12.7,1x))') ipa,csite%dist_type(ipa)              &
													   ,csite%age(ipa),csite%area(ipa)         &
													   ,pot_area_harv(ipa)
		  end do
		  write (unit=*,fmt='(a)'      )     '------------------------------------------------'
		  write (unit=*,fmt='(a)'      )     ' '
	   end if
	   !---------------------------------------------------------------------------------------!


	   !---------------------------------------------------------------------------------------!
	   !      Reset the primary forest memory, and assign any remaining deficit to the         !
	   ! secondary forest memory.                                                              !
	   !---------------------------------------------------------------------------------------!
	   cpoly%primary_harvest_memory  (isi) = 0.0
	   cpoly%secondary_harvest_memory(isi) = harvest_deficit
	   !---------------------------------------------------------------------------------------!
   end select
   return
end subroutine find_harvest_area
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This sub-routine calculates the area and total biomass associated with mature primary !
! forests, mature secondary forests, and mature forest plantations.                        !
!------------------------------------------------------------------------------------------!
subroutine inventory_mat_forests(cpoly,isi,onsp,harvestable_agb                            &
                                ,area_mature_primary   , hvagb_mature_primary              &
                                ,area_mature_secondary , hvagb_mature_secondary            &
                                ,area_mature_plantation, hvagb_mature_plantation )
   use ed_state_vars , only : polygontype         & ! structure
                            , sitetype            ! ! structure
   use disturb_coms  , only : plantation_rotation & ! intent(in)
                            , mature_harvest_age  & ! intent(in)
                            , min_harvest_biomass ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
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
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)       , pointer     :: csite
   integer                            :: ipa
   !---------------------------------------------------------------------------------------!



   !----- Initialize inventory. -----------------------------------------------------------!
   area_mature_primary     = 0.0
   area_mature_secondary   = 0.0
   area_mature_plantation  = 0.0
   hvagb_mature_primary    = 0.0
   hvagb_mature_secondary  = 0.0
   hvagb_mature_plantation = 0.0
   !---------------------------------------------------------------------------------------!

   csite => cpoly%site(isi)
   patchloop: do ipa=1,csite%npatches

      !----- Skip the patch if the biomass is low. ----------------------------------------!
      if (harvestable_agb(ipa) < min_harvest_biomass) cycle patchloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Run some checks to determine the type of forest and whether it has reached    !
      ! maturity.  TO DO: check whether George Hurtt's data consider burnt patches primary !
      ! or secondary.                                                                      !
      !------------------------------------------------------------------------------------!
      select case (csite%dist_type(ipa))
      case (2)
         !---- Forest plantation. ---------------------------------------------------------!
         if (csite%age(ipa) > plantation_rotation) then
            area_mature_plantation   = area_mature_plantation + csite%area(ipa)
            hvagb_mature_plantation  = hvagb_mature_plantation                             &
                                     + harvestable_agb(ipa)   * csite%area(ipa)
         end if
         !---------------------------------------------------------------------------------!

      case (3)
         !---- Primary forest. ------------------------------------------------------------!
         if (csite%age(ipa) > mature_harvest_age) then
            area_mature_primary      = area_mature_primary + csite%area(ipa)
            hvagb_mature_primary     = hvagb_mature_primary                                &
                                     + harvestable_agb(ipa)   * csite%area(ipa)
         end if
         !---------------------------------------------------------------------------------!

      case (4,5,6)
         !---- Secondary forest. ----------------------------------------------------------!
         if (csite%age(ipa) > mature_harvest_age) then
            area_mature_secondary   = area_mature_secondary   + csite%area(ipa)
            hvagb_mature_secondary  = hvagb_mature_secondary                               &
                                    + harvestable_agb(ipa)    * csite%area(ipa)
         end if
         !---------------------------------------------------------------------------------!
      case default
         !----- Probably pasture/croplands, skip the patch. -------------------------------!
         continue
         !---------------------------------------------------------------------------------!
      end select 
      !------------------------------------------------------------------------------------!
   end do patchloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine inventory_mat_forests
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine adjusts the mature forest disturbance rates associated with         !
! harvest.  In case the target biomass for primary forests is not met, it first tries to   !
! take the remaining biomass from old secondary forests before going to immature forests.  !
!------------------------------------------------------------------------------------------!
subroutine mat_forest_harv_rates(hvagb_mature_primary,hvagb_mature_secondary               &
                                ,hvagb_mature_plantation,primary_harvest_target            &
                                ,secondary_harvest_target,lambda_mature_primary            &
                                ,lambda_mature_secondary,lambda_mature_plantation          &
                                ,harvest_deficit)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)    :: hvagb_mature_primary
   real, intent(in)    :: hvagb_mature_secondary
   real, intent(in)    :: hvagb_mature_plantation
   real, intent(in)    :: primary_harvest_target
   real, intent(inout) :: secondary_harvest_target
   real, intent(out)   :: lambda_mature_primary
   real, intent(out)   :: lambda_mature_plantation
   real, intent(out)   :: lambda_mature_secondary
   real, intent(out)   :: harvest_deficit
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Find harvesting rate in mature primary vegetation.  In case there is not enough    !
   ! biomass harvest, harvest all primary vegetation then add the unmet biomass to the     !
   ! target for secondary vegetation.                                                      !
   !---------------------------------------------------------------------------------------!
   if (hvagb_mature_primary > primary_harvest_target) then
      lambda_mature_primary = primary_harvest_target / hvagb_mature_primary
   else
      lambda_mature_primary    = 1.0
      harvest_deficit          = primary_harvest_target   - hvagb_mature_primary
      secondary_harvest_target = secondary_harvest_target + harvest_deficit
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find harvesting rate for mature plantations and mature secondary forests.        !
   ! First try to remove all biomass from plantations.  In case there isn't sufficient     !
   ! biomass, harvest all plantations then remove the unmet biomass from mature secondary  !
   ! vegetation.  In case there isn't enough biomass, leave the remaining target in        !
   ! harvest_deficit, which will be used to determine harvesting from immature forests.    !
   !---------------------------------------------------------------------------------------!
   if (hvagb_mature_plantation > secondary_harvest_target) then
      lambda_mature_plantation = secondary_harvest_target / hvagb_mature_plantation
      lambda_mature_secondary  = 0.0
      harvest_deficit          = 0.0
   else
      lambda_mature_plantation = 1.0
      harvest_deficit          = secondary_harvest_target - hvagb_mature_plantation

      if(hvagb_mature_secondary > harvest_deficit) then
         lambda_mature_secondary = harvest_deficit / hvagb_mature_secondary
         harvest_deficit         = 0.0
      else
         lambda_mature_secondary = 1.0
         harvest_deficit         = harvest_deficit - hvagb_mature_secondary
      end if
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine mat_forest_harv_rates
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the area of mature patches that is going to be harvested.       !
!------------------------------------------------------------------------------------------!
subroutine area_harvest_mature(cpoly,isi,onsp,harvestable_agb,pot_area_harv                &
                              ,lambda_mature_primary,lambda_mature_secondary               &
                              ,lambda_mature_plantation)
   use ed_state_vars    , only : polygontype          & ! structure
                               , sitetype             ! ! structure
   use disturb_coms     , only : mature_harvest_age   & ! intent(in)
                               , plantation_rotation  & ! intent(in)
                               , min_harvest_biomass  ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)                  , target        :: cpoly
   integer                            , intent(in)    :: isi
   integer                            , intent(in)    :: onsp
   real, dimension(onsp)              , intent(inout) :: harvestable_agb
   real, dimension(onsp)              , intent(inout) :: pot_area_harv
   real                               , intent(in)    :: lambda_mature_plantation
   real                               , intent(in)    :: lambda_mature_secondary
   real                               , intent(in)    :: lambda_mature_primary
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)                     , pointer       :: csite
   integer                                            :: ipa
   !---------------------------------------------------------------------------------------!



   !----- Select patch. -------------------------------------------------------------------!
   csite => cpoly%site(isi)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Loop over patches.                                                               ! 
   !---------------------------------------------------------------------------------------!
   patch_loop: do ipa=1,onsp
  
      !----- Skip patch in case the biomass is less than the minimum for harvesting. ------!
      if (harvestable_agb(ipa) < min_harvest_biomass) cycle patch_loop
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find out whether to harvest this patch.                                         !
      !------------------------------------------------------------------------------------!
      select case (csite%dist_type(ipa))
      case (2)
         !----- Forest plantation. --------------------------------------------------------!
         if ( csite%age(ipa) > plantation_rotation ) then
            pot_area_harv(ipa) = csite%area(ipa) * lambda_mature_plantation
         end if
         !---------------------------------------------------------------------------------!

      case (3)
         !----- Primary forest. -----------------------------------------------------------!
         if ( csite%age(ipa) > mature_harvest_age  ) then
            pot_area_harv(ipa) = csite%area(ipa) * lambda_mature_primary
         end if
         !---------------------------------------------------------------------------------!

      case (4,5,6)
         !----- Secondary forest. ---------------------------------------------------------!
         if ( csite%age(ipa) > mature_harvest_age  ) then
            pot_area_harv(ipa) = csite%area(ipa) * lambda_mature_secondary
         end if
         !---------------------------------------------------------------------------------!
      case default
         !----- Agriculture.  Do not log. -------------------------------------------------!
         continue
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
   end do patch_loop
   !---------------------------------------------------------------------------------------!


   return
end subroutine area_harvest_mature
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the area of mature patches that has to be harvested to meet the !
! demand for biomass.                                                                      !
!------------------------------------------------------------------------------------------!
subroutine area_harvest_immature(cpoly,isi,onsp,harvestable_agb,pot_area_harv              &
                                ,harvest_deficit)
   use ed_state_vars     , only : polygontype          & ! structure
                                , sitetype             & ! structure
                                , patchtype            ! ! structure
   use disturb_coms      , only : plantation_rotation  & ! intent(in)
                                , mature_harvest_age   & ! intent(in)
                                , min_harvest_biomass  ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)                  , target        :: cpoly
   integer                            , intent(in)    :: isi
   integer                            , intent(in)    :: onsp
   real, dimension(onsp)              , intent(in)    :: harvestable_agb
   real, dimension(onsp)              , intent(inout) :: pot_area_harv
   real                               , intent(inout) :: harvest_deficit
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)                     , pointer       :: csite
   integer                                            :: ipa
   !---------------------------------------------------------------------------------------!



   !----- Select patch. -------------------------------------------------------------------!
   csite => cpoly%site(isi)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      First loop, try to obtain the biomass from forest plantations.                   !
   !---------------------------------------------------------------------------------------!
   patch_loop_fopl: do ipa=1,onsp
      !------------------------------------------------------------------------------------!
      !      Harvest this patch if it is has enough biomass, is a plantation, and is       !
      ! immature.                                                                          !
      !------------------------------------------------------------------------------------!
      if (  harvestable_agb(ipa) >= min_harvest_biomass  .and.                             &
            csite%dist_type(ipa) == 2                    .and.                             &
            csite%age(ipa)       <  plantation_rotation        ) then
         !----- Immature patch is harvestable.  Check how much to harvest. ----------------!
         if( (csite%area(ipa) * harvestable_agb(ipa)) > harvest_deficit) then
            !------------------------------------------------------------------------------!
            !      Biomass target has been met, partially harvest the patch, then quit the !
            ! sub-routine.                                                                 !
            !------------------------------------------------------------------------------!
            pot_area_harv(ipa) = harvest_deficit / harvestable_agb(ipa)
            harvest_deficit    = 0.0
            return
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !      Biomass target has not been met, harvest the entire patch, and keep     !
            ! searching for biomass.                                                       !
            !------------------------------------------------------------------------------!
            pot_area_harv(ipa) = csite%area(ipa)
            harvest_deficit    = harvest_deficit - csite%area(ipa) * harvestable_agb(ipa)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do patch_loop_fopl
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Second loop, try to obtain the biomass from secondary forests.                   !
   !---------------------------------------------------------------------------------------!
   patch_loop_2ary: do ipa=1,onsp
      !------------------------------------------------------------------------------------!
      !      Harvest this patch if it is has enough biomass, is a secondary forest, and is !
      ! immature.                                                                          !
      !------------------------------------------------------------------------------------!
      if (  harvestable_agb(ipa) >= min_harvest_biomass  .and.                             &
            csite%dist_type(ipa) >= 4                    .and.                             &
            csite%dist_type(ipa) <= 6                    .and.                             &
            csite%age(ipa)       <  mature_harvest_age         ) then
         !----- Immature patch is harvestable.  Check how much to harvest. ----------------!
         if( (csite%area(ipa) * harvestable_agb(ipa)) > harvest_deficit) then
            !------------------------------------------------------------------------------!
            !      Biomass target has been met, partially harvest the patch, then quit the !
            ! sub-routine.                                                                 !
            !------------------------------------------------------------------------------!
            pot_area_harv(ipa) = harvest_deficit / harvestable_agb(ipa)
            harvest_deficit    = 0.0
            return
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !      Biomass target has not been met, harvest the entire patch, and keep     !
            ! searching for biomass.                                                       !
            !------------------------------------------------------------------------------!
            pot_area_harv(ipa) = csite%area(ipa)
            harvest_deficit    = harvest_deficit - csite%area(ipa) * harvestable_agb(ipa)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do patch_loop_2ary
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Final loop, try to obtain the biomass from primary forests.                      !
   !---------------------------------------------------------------------------------------!
   patch_loop_1ary: do ipa=1,onsp
      !------------------------------------------------------------------------------------!
      !      Harvest this patch if it is has enough biomass, is a primary forest, and is   !
      ! immature.                                                                          !
      !------------------------------------------------------------------------------------!
      if (  harvestable_agb(ipa) >= min_harvest_biomass  .and.                             &
            csite%dist_type(ipa) == 3                    .and.                             &
            csite%age(ipa)       <  mature_harvest_age         ) then
         !----- Immature patch is harvestable.  Check how much to harvest. ----------------!
         if( (csite%area(ipa) * harvestable_agb(ipa)) > harvest_deficit) then
            !------------------------------------------------------------------------------!
            !      Biomass target has been met, partially harvest the patch, then quit the !
            ! sub-routine.                                                                 !
            !------------------------------------------------------------------------------!
            pot_area_harv(ipa) = harvest_deficit / harvestable_agb(ipa)
            harvest_deficit    = 0.0
            return
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !      Biomass target has not been met, harvest the entire patch, and keep     !
            ! searching for biomass.                                                       !
            !------------------------------------------------------------------------------!
            pot_area_harv(ipa) = csite%area(ipa)
            harvest_deficit    = harvest_deficit - csite%area(ipa) * harvestable_agb(ipa)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do patch_loop_1ary
   !---------------------------------------------------------------------------------------!

   return
end subroutine area_harvest_immature
!==========================================================================================!
!==========================================================================================!

end module forestry