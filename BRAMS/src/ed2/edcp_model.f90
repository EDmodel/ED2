!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the timestep in ED.  Formerly called from the main      !
! BRAMS driver, it is now called from BRAMS timestep driver, so we ensure synchronicity,   !
! particularly when multiple grids exist.                                                  !
!------------------------------------------------------------------------------------------!
subroutine ed_timestep()
  
   use grid_dims    , only : maxgrds             ! ! intent(in)
   use mem_grid     , only : ngrid               & ! intent(in)
                           , time                & ! intent(in)
                           , dtlt                ! ! intent
   use ed_node_coms , only : mynum               ! ! intent(in)
   use ed_misc_coms , only : dtlsm               & ! intent(in)
                           , current_time        ! ! intent(inout)
   use mem_edcp     , only : edtime1             & ! intent(out)
                           , edtime2             ! ! intent(out)
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
   integer                                 :: thishour
   integer                                 :: thismin
   integer                                 :: thissec
   !----- Local constants. ----------------------------------------------------------------!
   logical                    , parameter  :: print_banner = .false.
   !----- Locally saved variable to control when ED should be called. ---------------------!
   logical, dimension(maxgrds), save       :: first_time = .true.
   !----- External function. --------------------------------------------------------------!
   real                       , external   :: walltime
   !---------------------------------------------------------------------------------------!
     
   !----- Now we check whether this is the time to call ED. -------------------------------!
   if ( mod(time+dble(dtlt),dble(dtlsm)) < dble(dtlt) .or. first_time(ngrid) ) then

      thisyear  = current_time%year
      thismonth = current_time%month
      thisdate  = current_time%date
      thistime  = current_time%time
      thishour  = current_time%hour
      thismin   = current_time%min
      thissec   = current_time%sec

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
      

      edtime1  = time
      edtime2  = time + dble(dtlsm)

      !----- Call the actual model driver, for water, and for land. -----------------------!
      call simple_lake_model()
      call ed_coup_model(ngrid)

      !----- Copy the fluxes from ED to LEAF. ---------------------------------------------!
      call copy_fluxes_lsm2atm(ngrid)
      
      if (print_banner) then
         wtime2 = walltime(wtime_start)
         call timing(2,cputime2)

         if (mynum == 1) then
            write (unit=*,fmt='(a,i4,2(a,i2.2),a,i4.4,1x,3(i2.2,a),1x,a,2(1x,f7.3,a))')    &
                 ' ED2 LSM Timestep; Grid ',ngrid                                          &
                ,'; Sim time  ',thismonth, '-',thisdate,'-',thisyear                       &
                ,thishour,':',thismin,':',thissec,'UTC;'                                   &
                ,'Wall',wtime2-wtime1,'s; CPU',cputime2-cputime1,'s'
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
   !      Update the leaf_g surface flux arrays of tstar,rstar and ustar.  This is called  !
   ! every step because it is a time interpolation.                                        !
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
   use ed_max_dims  , only : maxgrds               ! ! intent(in)
   use ed_misc_coms , only : ivegt_dynamics        & ! intent(in)
                           , integration_scheme    & ! intent(in)
                           , simtime               & ! variable type
                           , current_time          & ! intent(inout)
                           , frqfast               & ! intent(in)
                           , frqstate              & ! intent(in)
                           , out_time_fast         & ! intent(in)
                           , dtlsm                 & ! intent(in)
                           , ifoutput              & ! intent(in)
                           , isoutput              & ! intent(in)
                           , idoutput              & ! intent(in)
                           , imoutput              & ! intent(in)
                           , iqoutput              & ! intent(in)
                           , itoutput              & ! intent(in)
                           , iyoutput              & ! intent(in)
                           , writing_dail          & ! intent(in)
                           , writing_mont          & ! intent(in)
                           , writing_dcyc          & ! intent(in)
                           , writing_year          & ! intent(in)
                           , writing_eorq          & ! intent(in)
                           , writing_long          & ! intent(in)
                           , frqsum                & ! intent(inout)
                           , unitfast              & ! intent(in)
                           , unitstate             & ! intent(in)
                           , imontha               & ! intent(in)
                           , iyeara                & ! intent(in)
                           , outstate              & ! intent(in)
                           , outfast               & ! intent(in)
                           , nrec_fast             & ! intent(in)
                           , nrec_state            & ! intent(in)
                           , outputmonth           ! ! intent(in)
   use grid_coms    , only : ngrids                & ! intent(in)
                           , istp                  & ! intent(in)
                           , time                  & ! intent(inout)
                           , timmax                ! ! intent(in)
   use ed_state_vars, only : edgrid_g              & ! intent(inout)
                           , edtype                & ! variable type
                           , patchtype             & ! variable type
                           , filltab_alltypes      & ! subroutine
                           , filltables            ! ! intent(in)
   use rk4_driver   , only : rk4_timestep          ! ! subroutine
   use rk4_coms     , only : record_err            & ! intent(out)
                           , print_detailed        & ! intent(out)
                           , print_thbnd           ! ! intent(out)
   use ed_node_coms , only : mynum                 & ! intent(in)
                           , nnodetot              ! ! intent(in)
   use mem_polygons , only : maxpatch              & ! intent(in)
                           , maxcohort             ! ! intent(in)
   use consts_coms  , only : day_sec               ! ! intent(in)
   use io_params    , only : ioutput               ! ! intent(in)
   use average_utils, only : update_ed_yearly_vars ! ! sub-routine
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
   logical                               :: dcycle_time
   logical                               :: annual_time
   logical                               :: mont_analy_time
   logical                               :: dail_analy_time
   logical                               :: dcyc_analy_time
   logical                               :: reset_time
   integer                               :: ndays
   integer                               :: jfm
   !----- External functions. -------------------------------------------------------------!
   integer                    , external :: num_days
   !----- Locally saved variables. --------------------------------------------------------!
   logical                    , save     :: first_time = .true.
   logical, dimension(maxgrds), save     :: calledgrid
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Save some of the output control variables, and settings that should never happen  !
   ! in a coupled model.  The test can be done only once.                                  !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      filltables        = .false.
      record_err        = .false.
      print_detailed    = .false.
      print_thbnd       = .false.
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


   !---------------------------------------------------------------------------------------!
   !     At this point, all meteorologic driver data for the land surface model has been   !
   ! updated for the current timestep.  Perform the time space average for the output      !
   ! diagnostic.                                                                           !
   !---------------------------------------------------------------------------------------!
   call integrate_ed_fmean_met_vars(edgrid_g(ifm))
   !---------------------------------------------------------------------------------------!


   !----- Solve the enthalpy, water, and carbon budgets. ----------------------------------!
   select case (integration_scheme)
   case (0)
      call euler_timestep(edgrid_g(ifm))
   case (1)
      call rk4_timestep(edgrid_g(ifm),ifm)
   case (2)
      call heun_timestep(edgrid_g(ifm))
   case (3)
      call hybrid_timestep(edgrid_g(ifm))
   end select


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
      !     Check whether now is any of those special times...                             !
      !------------------------------------------------------------------------------------!
      new_day         = current_time%time  <  dtlsm
      new_month       = current_time%date  == 1  .and. new_day
      new_year        = current_time%month == 1  .and. new_month
      mont_analy_time = new_month .and. writing_mont
      dcyc_analy_time = new_month .and. writing_dcyc
      annual_time     = new_month .and. writing_year .and.                                 &
                        current_time%month == outputmonth
      dail_analy_time = new_day   .and. writing_dail
      reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
      the_end         = mod(time,timmax)       < dble(dtlsm)

      !----- Check whether this is time to write fast analysis output or not. -------------!
      select case (unitfast)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         analysis_time   = mod(current_time%time, frqfast) < dtlsm .and.                   &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. ioutput /= 0)
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      case (2)   !----- Months, analysis time is at the new month -------------------------!
         analysis_time   = new_month .and.                                                 &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. ioutput /= 0) .and.      &
                           mod(real(12+current_time%month-imontha),frqfast) == 0.
         dcycle_time     = .false.
      case (3) !----- Year, analysis time is at the same month as initial time ------------!
         analysis_time   = new_month  .and.                                                &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. ioutput /= 0) .and.      &
                           current_time%month == imontha .and.                             &
                           mod(real(current_time%year-iyeara),frqfast) == 0.
         dcycle_time     = .false.
      end select

      !----- Check whether this is time to write restart output or not. -------------------!
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
      !    Update nrec_fast and nrec_state if it is a new month and outfast/outstate are   !
      ! monthly and frqfast/frqstate are daily or by seconds.                              !
      !------------------------------------------------------------------------------------!
      if (new_month) then
         ndays = num_days(current_time%month,current_time%year)
         if (outfast  == -2.) nrec_fast  = ndays*ceiling(day_sec/frqfast)
         if (outstate == -2.) nrec_state = ndays*ceiling(day_sec/frqstate)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the model output driver.                                                  !
      !------------------------------------------------------------------------------------!
      call ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,dcyc_analy_time &
                    ,annual_time,history_time,dcycle_time,the_end)
      !------------------------------------------------------------------------------------!



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



      !------------------------------------------------------------------------------------!
      !      Reset time happens every frqsum. This is to avoid variables to build up when  !
      ! history and analysis are off.  Put outside ed_output so we have a chance to copy   !
      ! some of these to BRAMS structures.                                                 !
      !------------------------------------------------------------------------------------!
      if (reset_time) then
         do jfm=1,ngrids
            call zero_ed_fmean_vars(edgrid_g(jfm))
         end do
      end if

      !------------------------------------------------------------------------------------!
      !     Check whether this is the beginning of a new simulated day.  Longer-scale      !
      ! processes, those which are updated once a day, once a month, or once a year, are   !
      ! updated here.  Since the subroutines here have internal grid loops, we only call   !
      ! this part of the code when all grids have reached this point.                      !
      !------------------------------------------------------------------------------------!
      if (new_day) then

         !---------------------------------------------------------------------------------!
         !     Compute phenology, growth, mortality, recruitment, disturbance, and check   !
         ! whether we will apply them to the ecosystem or not.                             !
         !---------------------------------------------------------------------------------!
         select case (ivegt_dynamics)
         case (0)
            !------------------------------------------------------------------------------!
            !     Dummy vegetation dynamics, we compute the tendencies but we don't really !
            ! apply to the vegetation, so they will remain constant throughout the entire  !
            ! simulation.                                                                  !
            !------------------------------------------------------------------------------!
            call vegetation_dynamics_eq_0(new_month,new_year)
            !------------------------------------------------------------------------------!

         case (1)
            !------------------------------------------------------------------------------!
            !     Actual vegetation dynamics, we compute the tendencies and apply to the   !
            ! vegetation.                                                                  !
            !------------------------------------------------------------------------------!
            call vegetation_dynamics(new_month,new_year)
            !------------------------------------------------------------------------------!

         end select
         !---------------------------------------------------------------------------------!


         
         !---------------------------------------------------------------------------------!
         !    First day of a month.  On the monthly timestep we have performed various     !
         ! fusion, fission, and extinction calls.  Therefore the var-table's pointer       !
         ! vectors must be updated, and the global definitions of the total numbers must   !
         ! be exported to every node.                                                      !
         !---------------------------------------------------------------------------------!
         if (new_month) then
            !----- Flag to run the subroutine filltab_alltypes. ---------------------------!
            if (maxcohort >= 0 .or. maxpatch >= 0) filltables=.true.

            !---- Also, we must re-allocate the cohort-level integration buffer. ----------!
            call initialize_rk4patches(.false.)
         end if
      end if

      if (new_day .and. new_month) then
         do jfm=1,ngrids
            call updateHydroParms(edgrid_g(jfm))
         end do
      end if
   

      if (analysis_time .and. new_month .and. new_day .and. current_time%month == 6) then
         do jfm = 1,ngrids
            call update_ed_yearly_vars(edgrid_g(jfm))
         end do
      end if

      !----- Update lateral hydrology. ----------------------------------------------------!
      call calcHydroSubsurface()
      call calcHydroSurface()
      call writeHydro()

   end if

   return
end subroutine ed_coup_model
!==========================================================================================!
!==========================================================================================!
