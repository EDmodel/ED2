!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main driver for the longer-term vegetation dynamics.  This    !
! has become a file by itself to reduce the number of sub-routines that are doubled        !
! between ED-2.1 stand alone and the coupled model.                                        !
!------------------------------------------------------------------------------------------!
subroutine vegetation_dynamics(new_month,new_year)
   use grid_coms        , only : ngrids
   use ed_misc_coms     , only : current_time           & ! intent(in)
                               , dtlsm                  & ! intent(in)
                               , frqsum                 & ! intent(in)
                               , ied_init_mode          ! ! intent(in)
   use disturb_coms     , only : include_fire           ! ! intent(in)
   use disturbance_utils, only : apply_disturbances     & ! subroutine
                               , site_disturbance_rates ! ! subroutine
   use fuse_fiss_utils  , only : fuse_patches           ! ! subroutine
   use ed_state_vars    , only : edgrid_g               & ! intent(inout)
                               , edtype                 ! ! variable type
   use growth_balive    , only : dbalive_dt             & ! subroutine
                               , dbalive_dt_eq_0        ! ! subroutine
   use consts_coms      , only : day_sec                & ! intent(in)
                               , yr_day                 ! ! intent(in)
   use mem_polygons     , only : maxpatch               ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical     , intent(in)   :: new_month
   logical     , intent(in)   :: new_year
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype), pointer      :: cgrid
   real                       :: tfact1
   real                       :: tfact2
   integer                    :: doy
   integer                    :: ip
   integer                    :: isite
   integer                    :: ifm
   !----- External functions. -------------------------------------------------------------!
   integer     , external     :: julday
   !---------------------------------------------------------------------------------------!

   !----- Find the day of year. -----------------------------------------------------------!
   doy = julday(current_time%month, current_time%date, current_time%year)
  
   !----- Time factor for normalizing daily variables updated on the DTLSM step. ----------!
   tfact1 = dtlsm / day_sec
   !----- Time factor for averaging dailies. ----------------------------------------------!
   tfact2 = 1.0 / yr_day

   !----- Apply events. -------------------------------------------------------------------!
   call prescribed_event(current_time%year,doy)

  
   !---------------------------------------------------------------------------------------!
   !   Loop over all domains.                                                              !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids

      cgrid => edgrid_g(ifm) 

      !------------------------------------------------------------------------------------!
      !     The following block corresponds to the daily time-step.                        !
      !------------------------------------------------------------------------------------!
      !----- Standardise the fast-scale uptake and respiration, for growth rates. ---------!
      call normalize_ed_daily_vars(cgrid, tfact1)
      !----- Update phenology and growth of live tissues. ---------------------------------!
      select case (ied_init_mode)
      case (-8)
         !----- Special case, in which we don't solve the actual vegetation dynamics. -----!
         call phenology_driver_eq_0(cgrid,doy,current_time%month, tfact1)
         call dbalive_dt_eq_0(cgrid,tfact2)
      case default
         call phenology_driver(cgrid,doy,current_time%month, tfact1)
         call dbalive_dt(cgrid,tfact2)
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following block corresponds to the monthly time-step:                      !
      !------------------------------------------------------------------------------------!
      if (new_month) then

         !----- Update the mean workload counter. -----------------------------------------!
         call update_workload(cgrid)

         !----- Update the growth of the structural biomass. ------------------------------!
         call structural_growth(cgrid, current_time%month)

         !----- Solve the reproduction rates. ---------------------------------------------!
         call reproduction(cgrid,current_time%month)

         !----- Update the fire disturbance rates. ----------------------------------------!
         if (include_fire /= 0) then
            call fire_frequency(current_time%month,cgrid)
         end if

         !----- Update the disturbance rates. ---------------------------------------------!
         call site_disturbance_rates(current_time%month, current_time%year, cgrid)

         !----- This is actually the yearly time-step, apply the disturbances. ------------!
         if (new_year) then
            call apply_disturbances(cgrid)
         end if
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     This should be done every day, but after the longer-scale steps.  We update    !
      ! the carbon and nitrogen pools, and re-set the daily variables.                     !
      !------------------------------------------------------------------------------------!
      call update_C_and_N_pools(cgrid)
      call zero_ed_daily_vars(cgrid)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Fuse patches last, after all updates have been applied.  This reduces the     !
      ! number of patch variables that actually need to be fused.                          !
      !------------------------------------------------------------------------------------!
      if(new_year) then
         if (maxpatch >= 0) call fuse_patches(cgrid,ifm)
      end if
      !------------------------------------------------------------------------------------!



      !----- Recalculate the AGB and basal area at the polygon level. ---------------------!
      call update_polygon_derived_props(cgrid)
      call print_C_and_N_budgets(cgrid)
      !------------------------------------------------------------------------------------!
   end do

   return
end subroutine vegetation_dynamics
!==========================================================================================!
!==========================================================================================!
