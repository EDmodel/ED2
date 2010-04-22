!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the timestep in ED.  Formerly called from the main      !
! BRAMS driver, it is now called from BRAMS timestep driver, so we ensure synchronicity,   !
! particularly when multiple grids exist.                                                  !
!------------------------------------------------------------------------------------------!
subroutine ed_timestep()
  
   use grid_dims   , only : maxgrds      ! ! intent(in)
   use mem_grid    , only : ngrid        & ! intent(in)
                          , time         & ! intent(in)
                          , dtlt         ! ! intent
   use ed_node_coms, only : mynum        ! ! intent(in)
   use ed_misc_coms, only : dtlsm        & ! intent(in)
                          , current_time ! ! intent(inout)
   use mem_edcp    , only : edtime1      & ! intent(out)
                          , edtime2      ! ! intent(out)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)                            :: thistime
   real                                    :: wtime_start
   real                                    :: cputime1
   real                                    :: cputime2
   real                                    :: wtime1
   real                                    :: wtime2
   integer                                 :: thismonth
   integer                                 :: thisyear
   integer                                 :: thisdate
   !----- Local constants. ----------------------------------------------------------------!
   logical                    , parameter  :: print_banner = .false.
   !----- Locally saved variable to control when ED should be called. ---------------------!
   logical, dimension(maxgrds), save       :: first_time = .true.
   !----- External function. --------------------------------------------------------------!
   real                       , external   :: walltime
   !---------------------------------------------------------------------------------------!


  
   !---------------------------------------------------------------------------------------!
   !     Now the solve the water fluxes.  This is called every time step.                  !
   !---------------------------------------------------------------------------------------!
   call simple_lake_model(time,dtlt)
     
   !----- Now we check whether this is the time to call ED. -------------------------------!
   if ( mod(time+dble(dtlt),dble(dtlsm)) < dble(dtlt) .or. first_time(ngrid) ) then

      thisyear  = current_time%year
      thismonth = current_time%month
      thisdate  = current_time%date
      thistime  = current_time%time

      if (print_banner) then
         wtime_start=walltime(0.)

         !----- Finding the CPU and model timing information. -----------------------------!
         call timing(1,cputime1)
         wtime1 = walltime(wtime_start)
      end if


      !------------------------------------------------------------------------------------!
      !     Transfer the fluxes from the terrestrial and water models to the previous time !
      ! arrays.                                                                            !
      !------------------------------------------------------------------------------------!
      call copy_fluxes_future_2_past(ngrid)
      call copy_atm2lsm(ngrid,.false.)
      
      !----- Call the actual model driver. ------------------------------------------------!
      call ed_coup_model(ngrid)

      edtime1  = time
      edtime2  = time + dble(dtlsm)

      !----- Copy the fluxes from ED to LEAF. ---------------------------------------------!
      call copy_fluxes_lsm2atm(ngrid)
      
      if (print_banner) then
         wtime2 = walltime(wtime_start)
         call timing(2,cputime2)

         if (mynum == 1) then
            write (unit=*,fmt='(a,i4,2(a,i2.2),a,i4.4,1x,f6.0,a,2(1x,f7.3,a))')            &
                 ' ED2 LSM Timestep; Grid ',ngrid,'; Sim time  ',thismonth, '-',thisdate   &
                ,'-',thisyear,thistime,'s; Wall',wtime2-wtime1,'s; CPU',cputime2-cputime1  &
                ,'s'
         end if
      end if


      
      !------------------------------------------------------------------------------------!
      !      If this is the first time the routine is called, then the ed_fluxp_g arrays   !
      ! (the flux arrays from previous time) are zero.  So just this once, copy the future !
      ! arrays after they have been populated with real values, into the past arrays.  the !
      ! transfer function needs both in order to do a temporal interpolation.              !
      !------------------------------------------------------------------------------------!
      if (first_time(ngrid)) call copy_fluxes_future_2_past(ngrid)

      first_time(ngrid)=.false.
   end if
  

   !---------------------------------------------------------------------------------------!
   !      Update the leaf_g surface flux arrays of tstar,rstar and ustar.  This gets       !
   ! called every step because it is a time interpolation.                                 !
   !---------------------------------------------------------------------------------------!
   call transfer_ed2leaf(ngrid,time)



   return
end subroutine ed_timestep
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the ED model at one time step.  Notice that differently !
! from the offline model, this is called for each grid.  Also the time change is done out- !
! side this subroutine to avoid confusion when we run nested grid simulations.             !
!------------------------------------------------------------------------------------------!
subroutine ed_coup_model(ifm)
   use ed_max_dims  , only : maxgrds            ! ! intent(in)
   use ed_misc_coms , only : integration_scheme & ! intent(in)
                           , simtime            & ! variable type
                           , current_time       & ! intent(inout)
                           , frqfast            & ! intent(in)
                           , frqstate           & ! intent(in)
                           , out_time_fast      & ! intent(in)
                           , dtlsm              & ! intent(in)
                           , ifoutput           & ! intent(in)
                           , isoutput           & ! intent(in)
                           , idoutput           & ! intent(in)
                           , imoutput           & ! intent(in)
                           , itoutput           & ! intent(in)
                           , iyoutput           & ! intent(in)
                           , frqsum             & ! intent(inout)
                           , unitfast           & ! intent(in)
                           , unitstate          & ! intent(in)
                           , imontha            & ! intent(in)
                           , iyeara             & ! intent(in)
                           , outstate           & ! intent(in)
                           , outfast            & ! intent(in)
                           , nrec_fast          & ! intent(in)
                           , nrec_state         & ! intent(in)
                           , outputmonth        ! ! intent(in)
   use grid_coms    , only : ngrids             & ! intent(in)
                           , istp               & ! intent(in)
                           , time               & ! intent(inout)
                           , timmax             ! ! intent(in)
   use ed_state_vars, only : edgrid_g           & ! intent(inout)
                           , edtype             & ! variable type
                           , patchtype          & ! variable type
                           , filltab_alltypes   ! ! subroutine
   use rk4_driver   , only : rk4_timestep       ! ! subroutine
   use ed_node_coms , only : mynum              & ! intent(in)
                           , nnodetot           ! ! intent(in)
   use mem_polygons , only : maxpatch           & ! intent(in)
                           , maxcohort          ! ! intent(in)
   use consts_coms  , only : day_sec            ! ! intent(in)
   use io_params    , only : ioutput            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)                   :: ifm
   !----- Local variables. ----------------------------------------------------------------!
   logical                               :: analysis_time
   logical                               :: new_day
   logical                               :: new_month
   logical                               :: new_year
   logical                               :: the_end
   logical                               :: history_time
   logical                               :: annual_time
   logical                               :: mont_analy_time
   logical                               :: dail_analy_time
   logical                               :: reset_time
   integer                               :: ndays
   integer                               :: jfm
   !----- External functions. -------------------------------------------------------------!
   integer                    , external :: num_days
   !----- Locally saved variables. --------------------------------------------------------!
   logical                    , save     :: first_time = .true.
   logical                    , save     :: writing_dail
   logical                    , save     :: writing_mont
   logical                    , save     :: writing_year
   logical, dimension(maxgrds), save     :: calledgrid
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Saving some of the output control variables. The test can be done only once.      !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      writing_dail      = idoutput > 0
      writing_mont      = imoutput > 0
      writing_year      = iyoutput > 0
      calledgrid(:)     = .false.
      first_time        = .false.
   end if

   !---------------------------------------------------------------------------------------!
   !     Flagging that this grid has been called.                                          !
   !---------------------------------------------------------------------------------------!
   calledgrid(ifm) = .true.

   istp                = istp + 1
   out_time_fast       = current_time
   out_time_fast%month = -1


   !----- Updating the cohort status (stable to be solved or unstable). -------------------!
   call flag_stable_cohorts(edgrid_g(ifm))

   !----- Radiation scheme. ---------------------------------------------------------------!
   call radiate_driver(edgrid_g(ifm))
   
   !----- Solve the enthalpy, water, and carbon budgets. ----------------------------------!
   select case (integration_scheme)
   case (0)
      call euler_timestep(edgrid_g(ifm))
   case (1,2)
      call rk4_timestep(edgrid_g(ifm),ifm)
   end select

   !---------------------------------------------------------------------------------------!
   !     Update the daily averages if daily or monthly analysis are needed.                !
   !---------------------------------------------------------------------------------------!
   if (writing_dail .or. writing_mont) then
      call integrate_ed_daily_output_state(edgrid_g(ifm))
   end if


   !---------------------------------------------------------------------------------------!
   !     The remainder of this subroutine is called only once, and it is done after all    !
   ! grids have been called.                                                               !
   !---------------------------------------------------------------------------------------!
   if (all(calledgrid(1:ngrids))) then
   
      !----- Reset all grids to false for next time step. ---------------------------------!
      calledgrid(1:ngrids) = .false.
   
      !----- Updating the time. -----------------------------------------------------------!
      time = time + dble(dtlsm)
      call update_model_time_dm(current_time, dtlsm)

      !------------------------------------------------------------------------------------!
      !     Checking whether now is any of those special times...                          !
      !------------------------------------------------------------------------------------!
      new_day         = current_time%time  <  dtlsm
      new_month       = current_time%date  == 1  .and. new_day
      new_year        = current_time%month == 1  .and. new_month
      mont_analy_time = new_month .and. writing_mont
      annual_time     = new_month .and. writing_year .and.                                 &
                        current_time%month == outputmonth
      dail_analy_time = new_day   .and. writing_dail
      reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
      the_end         = mod(time,timmax)       < dble(dtlsm)

      !----- Checking whether this is time to write fast analysis output or not. ----------!
      select case (unitfast)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         analysis_time   = mod(current_time%time, frqfast) < dtlsm .and.                   &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. ioutput /= 0)
      case (2)   !----- Months, analysis time is at the new month -------------------------!
         analysis_time   = new_month .and.                                                 &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. ioutput /= 0) .and.      &
                           mod(real(12+current_time%month-imontha),frqfast) == 0.
      case (3) !----- Year, analysis time is at the same month as initial time ------------!
         analysis_time   = new_month  .and.                                                &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. ioutput /= 0) .and.      &
                           current_time%month == imontha .and.                             &
                           mod(real(current_time%year-iyeara),frqfast) == 0.
      end select

      !----- Checking whether this is time to write restart output or not. ----------------!
      select case(unitstate)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         history_time   = mod(current_time%time, frqstate) < dtlsm .and.                   &
                          (isoutput /= 0 .or. ioutput /= 0)
      case (2)   !----- Months, history time is at the new month --------------------------!
         history_time   = new_month .and. isoutput /= 0 .and.                              &
                          mod(real(12+current_time%month-imontha),frqstate) == 0.
      case (3) !----- Year, history time is at the same month as initial time -------------!
         history_time   = new_month .and. isoutput /= 0 .and.                              &
                          current_time%month == imontha .and.                              &
                          mod(real(current_time%year-iyeara),frqstate) == 0.
      end select

      !------------------------------------------------------------------------------------!
      !    Updating nrec_fast and nrec_state if it is a new month and outfast/outstate are !
      ! monthly and frqfast/frqstate are daily or by seconds.                              !
      !------------------------------------------------------------------------------------!
      if (new_month) then
         ndays = num_days(current_time%month,current_time%year)
         if (outfast  == -2.) nrec_fast  = ndays*ceiling(day_sec/frqfast)
         if (outstate == -2.) nrec_state = ndays*ceiling(day_sec/frqstate)
      end if


      !----- Call the output subroutine. --------------------------------------------------!
      call ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,annual_time     &
                    ,writing_dail,writing_mont,history_time,the_end)

      !------------------------------------------------------------------------------------!
      !     If this is the time to write the output, then send the data back to the LEAF   !
      ! arrays.                                                                            !
      !------------------------------------------------------------------------------------!
      if (analysis_time) then
         do jfm = 1, ngrids
            call copy_avgvars_to_leaf(jfm)
         end do
      end if

      !------------------------------------------------------------------------------------!
      ! Reset time happens every frqsum. This is to avoid variables to build up when       !
      ! history and analysis are off.  Put outside ed_output so we have a chance to copy   !
      ! some of these to BRAMS structures.                                                 !
      !------------------------------------------------------------------------------------!
      if (reset_time) then
         do jfm=1,ngrids
            call reset_averaged_vars(edgrid_g(jfm))
         end do
      end if

      !------------------------------------------------------------------------------------!
      !     Checking whether this is the beginning of a new simulated day.  Longer-scale   !
      ! processes, those which are updated once a day, once a month, or once a year, are   !
      ! updated here.  Since the subroutines here have internal grid loops, we only call   !
      ! this part of the code when all grids have reached this point.                      !
      !------------------------------------------------------------------------------------!
      if (new_day) then
         
         !----- Do phenology, growth, mortality, recruitment, disturbance. ----------------!
         call vegetation_dynamics(new_month,new_year)
         
         !---------------------------------------------------------------------------------!
         !    First day of a month.  On the monthly timestep we have performed various     !
         ! fusion, fission, and extinction calls.  Therefore the var-table's pointer       !
         ! vectors must be updated, and the global definitions of the total numbers must   !
         ! be exported to every node.                                                      !
         !---------------------------------------------------------------------------------!
         if (new_month .and. (maxpatch >= 0 .or. maxcohort >= 0)) then
            call filltab_alltypes

            !---- Also, we must re-allocate integration buffer if RK4 is being used. ------!
            if (integration_scheme == 1 .or. integration_scheme == 2) then
               call initialize_rk4patches(0)
            end if
         end if
      end if

      if (new_day .and. new_month) then
         do jfm=1,ngrids
            call updateHydroParms(edgrid_g(jfm))
         end do
      end if
   

      !----- Update Lateral hydrology. ----------------------------------------------------!
      call calcHydroSubsurface()
      call calcHydroSurface()
      call writeHydro()

   end if

   return
end subroutine ed_coup_model
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_model_time_dm(ctime,dtlong)

   use ed_misc_coms, only : simtime
   use  consts_coms, only : day_sec
   implicit none

   type(simtime) :: ctime
   real, intent(in) :: dtlong
   logical, external :: isleap
   integer, dimension(12) :: daymax
  
   daymax=(/31,28,31,30,31,30,31,31,30,31,30,31/)


   ctime%time = ctime%time + dtlong
  
   if (ctime%time >= day_sec)then
      ctime%time = ctime%time - day_sec
      ctime%date = ctime%date + 1

      ! Before checking, adjust for leap year
      if (isleap(ctime%year)) daymax(2) = 29
    
      if (ctime%date > daymax(ctime%month)) then
         ctime%date  = 1
         ctime%month = ctime%month + 1
      
         if(ctime%month == 13)then
            ctime%month = 1
            ctime%year = ctime%year + 1
         endif
      endif
      
   elseif(ctime%time < 0.0)then
      ctime%time = ctime%time + day_sec
      ctime%date = ctime%date - 1

      if(ctime%date == 0)then
         ctime%month = ctime%month - 1
         
         if(ctime%month == 0)then
            ctime%month = 12
            ctime%year = ctime%year - 1
            ctime%date = daymax(12)

         else
            if (isleap(ctime%year)) daymax(2) = 29
            ctime%date = daymax(ctime%month)
         end if
      end if
   end if

   return
end subroutine update_model_time_dm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the longer-term vegetation dynamics.                      !
!------------------------------------------------------------------------------------------!
subroutine vegetation_dynamics(new_month,new_year)
   use grid_coms        , only : ngrids
   use ed_misc_coms     , only : current_time           & ! intent(in)
                               , dtlsm                  & ! intent(in)
                               , frqsum                 ! ! intent(in)
   use disturb_coms     , only : include_fire           ! ! intent(in)
   use disturbance_utils, only : apply_disturbances     & ! subroutine
                               , site_disturbance_rates ! ! subroutine
   use fuse_fiss_utils  , only : fuse_patches           ! ! subroutine
   use ed_state_vars    , only : edgrid_g               & ! intent(inout)
                               , edtype                 ! ! variable type
   use growth_balive    , only : dbalive_dt             ! ! subroutine
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

  
   do ifm=1,ngrids

      cgrid => edgrid_g(ifm) 
      call normalize_ed_daily_vars(cgrid, tfact1)
      call phenology_driver(cgrid,doy,current_time%month, tfact1)
      call dbalive_dt(cgrid,tfact2)

      if(new_month)then

         call update_workload(cgrid)

         call structural_growth(cgrid, current_time%month)
         call reproduction(cgrid,current_time%month)

         if(include_fire /= 0) then
            call fire_frequency(current_time%month,cgrid)
         end if

         call site_disturbance_rates(current_time%month, current_time%year, cgrid)

         if(new_year) then
            !write (unit=*,fmt='(a)') '### Apply_disturbances...'
            call apply_disturbances(cgrid)
         end if
         
      end if

      call update_C_and_N_pools(cgrid)
      call zero_ed_daily_vars(cgrid)

      !------------------------------------------------------------------------------------!
      !      Fuse patches last, after all updates have been applied.  This reduces the     !
      ! number of patch variables that actually need to be fused.                          !
      !------------------------------------------------------------------------------------!
      if(new_year) then
         if (maxpatch >= 0) call fuse_patches(cgrid,ifm)
      end if

      !----- Recalculate the AGB and basal area at the polygon level. ---------------------!
      call update_polygon_derived_props(cgrid)
      call print_C_and_N_budgets(cgrid)

   end do

   return
end subroutine vegetation_dynamics
!==========================================================================================!
!==========================================================================================!
