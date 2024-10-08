!==========================================================================================!
!==========================================================================================!
!    This file contains the main driver for ED2.                                           !
!------------------------------------------------------------------------------------------!
! SUBROUTINE: ED_MODEL
!
!> \brief   Begins, updates, and outputs results from ecosystem simulation.
!> \details Coordinates meteorological forcing data, time-stepping, radiation, vegetation
!>          dynamics, hydraulogy, and the HDF5 output subsystem.
!> \author  Translated from ED1 by David Medvigy, Ryan Knox and Marcos Longo
!------------------------------------------------------------------------------------------!
subroutine ed_model()
   use ed_misc_coms        , only : simtime                     & ! structure
                                  , fmtrest                     & ! intent(in)
                                  , ivegt_dynamics              & ! intent(in)
                                  , integration_scheme          & ! intent(in)
                                  , current_time                & ! intent(in)
                                  , frqfast                     & ! intent(in)
                                  , frqstate                    & ! intent(in)
                                  , out_time_fast               & ! intent(in)
                                  , dtlsm                       & ! intent(in)
                                  , ifoutput                    & ! intent(in)
                                  , isoutput                    & ! intent(in)
                                  , iqoutput                    & ! intent(in)
                                  , itoutput                    & ! intent(in)
                                  , iooutput                    & ! intent(in)
                                  , restore_file                & ! intent(in)
                                  , frqsum                      & ! intent(in)
                                  , unitfast                    & ! intent(in)
                                  , unitstate                   & ! intent(in)
                                  , imontha                     & ! intent(in)
                                  , iyeara                      & ! intent(in)
                                  , outstate                    & ! intent(in)
                                  , outfast                     & ! intent(in)
                                  , nrec_fast                   & ! intent(in)
                                  , nrec_state                  & ! intent(in)
                                  , runtype                     & ! intent(in)
                                  , month_yrstep                & ! intent(in)
                                  , writing_dail                & ! intent(in)
                                  , writing_mont                & ! intent(in)
                                  , writing_dcyc                & ! intent(in)
                                  , writing_eorq                & ! intent(in)
                                  , writing_long                & ! intent(in)
                                  , writing_year                ! ! intent(in)
   use ed_init             , only : remove_obstime              & ! sub-routine
                                  , is_obstime                  ! ! sub-routine
   use grid_coms           , only : ngrids                      & ! intent(in)
                                  , istp                        & ! intent(in)
                                  , time                        & ! intent(in)
                                  , timmax                      ! ! intent(in)
   use ed_state_vars       , only : edgrid_g                    & ! intent(in)
                                  , edtype                      & ! intent(in)
                                  , patchtype                   & ! intent(in)
                                  , filltab_alltypes            & ! intent(in)
                                  , filltables                  ! ! intent(in)
   use rk4_driver          , only : rk4_timestep                ! ! sub-routine
   use rk4_coms            , only : integ_err                   & ! intent(in)
                                  , integ_lab                   & ! intent(in)
                                  , record_err                  & ! intent(inout)
                                  , print_detailed              & ! intent(inout)
                                  , nerr                        & ! intent(in)
                                  , errmax_fout                 & ! intent(in)
                                  , sanity_fout                 & ! intent(in)
                                  , alloc_integ_err             & ! subroutine
                                  , assign_err_label            & ! subroutine
                                  , reset_integ_err             ! ! subroutine
   use ed_node_coms        , only : mynum                       & ! intent(in)
                                  , nnodetot                    ! ! intent(in)
   use mem_polygons        , only : n_ed_region                 & ! intent(in)
                                  , n_poi                       ! ! intent(in)
   use consts_coms         , only : day_sec                     ! ! intent(in)
   use average_utils       , only : update_ed_yearly_vars       & ! sub-routine
                                  , zero_ed_dmean_vars          & ! sub-routine
                                  , zero_ed_mmean_vars          & ! sub-routine
                                  , zero_ed_qmean_vars          & ! sub-routine
                                  , zero_ed_fmean_vars          & ! sub-routine
                                  , integrate_ed_fmean_met_vars & ! sub-routine
                                  , zero_ed_yearly_vars         ! ! sub-routine
   use edio                , only : ed_output                   ! ! sub-routine
   use ed_met_driver       , only : read_met_drivers            & ! sub-routine
                                  , update_met_drivers          ! ! sub-routine
   use euler_driver        , only : euler_timestep              ! ! sub-routine
   use heun_driver         , only : heun_timestep               ! ! sub-routine
   use hybrid_driver       , only : hybrid_timestep             ! ! sub-routine
   use lsm_hyd             , only : updateHydroParms            & ! sub-routine
                                  , calcHydroSubsurface         & ! sub-routine
                                  , calcHydroSurface            & ! sub-routine
                                  , writeHydro                  ! ! sub-routine
   use radiate_driver      , only : canopy_radiation            ! ! sub-routine
   use rk4_integ_utils     , only : initialize_rk4patches       & ! sub-routine
                                  , initialize_misc_stepvars    ! ! sub-routine
   use stable_cohorts      , only : flag_stable_cohorts         ! ! sub-routine
   use update_derived_utils, only : update_model_time_dm        ! ! sub-routine
   use budget_utils        , only : ed_init_budget              ! ! sub-routine
   use vegetation_dynamics , only : veg_dynamics_driver         ! ! sub-routine
   use ed_type_init        , only : ed_init_viable              ! ! sub-routine
   use soil_respiration    , only : zero_litter_inputs          ! ! sub-routine
#if defined(RAMS_MPI)
   use mpi
#endif
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   type(simtime)      :: daybefore
   character(len=28)  :: fmthead
   character(len=32)  :: fmtcntr
   integer            :: ifm
   integer            :: nn
   integer            :: ndays
   integer            :: dbndays
   integer            :: obstime_idx
   logical            :: last_step
   logical            :: analysis_time
   logical            :: observation_time
   logical            :: new_day
   logical            :: new_month
   logical            :: new_year
   logical            :: history_time
   logical            :: dcycle_time
   logical            :: annual_time
   logical            :: mont_analy_time
   logical            :: dail_analy_time
   logical            :: dcyc_analy_time
   logical            :: reset_time
   logical            :: past_one_day
   logical            :: past_one_month
   logical            :: printbanner
   logical            :: veget_dyn_on
   real               :: wtime_start
   real               :: t1
   real               :: wtime1
   real               :: wtime2
   real               :: t2
   real               :: wtime_tot
   real               :: dbndaysi
   real               :: gr_tfact0
   !----- Local variables (MPI only). -----------------------------------------------------!
#if defined(RAMS_MPI)
   integer            :: ierr
#endif
   !----- Local constants. ----------------------------------------------------------------!
   logical         , parameter :: whos_slow=.false. ! This will print out node numbers
                                                    !    during synchronization, so you
                                                    !    can find out which node is the
                                                    !    slow one
   !----- External functions. -------------------------------------------------------------!
   real    , external :: walltime ! Wall time
   integer , external :: num_days ! Number of days in the current month
   !---------------------------------------------------------------------------------------!

   past_one_day   = .false.
   past_one_month = .false.
   filltables     = .false.

   !----- Print the hour banner only for regional runs that aren't massively parallel. ----!
   printbanner = n_ed_region > 0 .and. edgrid_g(1)%npolygons > 50 .and. mynum == 1

   !----- Run with vegetation dynamics turned on?  ----------------------------------------!
   veget_dyn_on = ivegt_dynamics == 1

   wtime_start=walltime(0.)
   istp = 0

   !---------------------------------------------------------------------------------------!
   !     If we are going to record the integrator errors, here is the time to open it for  !
   ! the first time and write the header.  But just before we do it, we check whether this !
   ! is a single POI run, the only case where we will allow this recording.                !
   !---------------------------------------------------------------------------------------!
   record_err     = record_err     .and. n_ed_region == 0 .and. n_poi == 1
   print_detailed = print_detailed .and. n_ed_region == 0 .and. n_poi == 1
   if(record_err) then
      !----- Initialise the error structures. ---------------------------------------------!
      call alloc_integ_err()
      call reset_integ_err()
      call assign_err_label()

      !----- Define the formats for both the header and the actual output. ----------------!
      write(fmthead,fmt='(a,i3.3,a)')  '(a4,1x,2(a3,1x),',nerr,'(a13,1x))'
      write(fmtcntr,fmt='(a,i3.3,a)')  '(i4.4,1x,2(i3.2,1x),',nerr,'(i13,1x))'

      open  (unit=77,file=trim(errmax_fout),form='formatted',status='replace')
      write (unit=77,fmt=fmthead) 'YEAR','MON','DAY',(integ_lab(nn),nn=1,nerr)
      close (unit=77,status='keep')

      open  (unit=78,file=trim(sanity_fout),form='formatted',status='replace')
      write (unit=78,fmt=fmthead) 'YEAR','MON','DAY',(integ_lab(nn),nn=1,nerr)
      close (unit=78,status='keep')
   end if

   out_time_fast     = current_time
   out_time_fast%month = -1

   !---------------------------------------------------------------------------------------!
   !      If this is not a history restart, then zero out the long term diagnostics.       !
   !---------------------------------------------------------------------------------------!
   select case (trim(runtype))
   case ('HISTORY')
      continue
   case default

      do ifm=1,ngrids
         if (writing_long) call zero_ed_dmean_vars(edgrid_g(ifm))
         if (writing_eorq) call zero_ed_mmean_vars(edgrid_g(ifm))
         if (writing_dcyc) call zero_ed_qmean_vars(edgrid_g(ifm))
      end do

      !----- Long-term dynamics structure. ------------------------------------------------!
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      The fast analysis is always reset, including history runs.                       !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids
      call zero_ed_fmean_vars(edgrid_g(ifm))
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Allocate memory to the integration patch, Euler now utilises the RK4 buffers too.  !
   !---------------------------------------------------------------------------------------!
   call initialize_rk4patches(.true.)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Initialize some stepping variables.                                                !
   !---------------------------------------------------------------------------------------!
   call initialize_misc_stepvars()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Here we must initialise or reset a group of variables.                           !
   ! 1.  Variable is_viable must be set to .true..  This variable is not saved in          !
   !     ed_init_history, and the default is .false., which would eliminate all cohorts.   !
   ! 2.  In the case of a initial simulation, we must reset all budget fluxes and set all  !
   !     budget stocks.  This should not be done in HISTORY initialisation, all variables  !
   !     should be read from history.                                                      !
   ! 3.  Litter inputs must be reset in the HISTORY initialisation, in case the history    !
   !     file is at midnight UTC (daily time step).  These variables are normally reset    !
   !     after writing the output so they are meaningful in the output.  Because history   !
   !     files are written before the inputs are reset, the inputs would be double-counted !
   !     in the second day of simulation.                                                  !
   !---------------------------------------------------------------------------------------!
   select case (trim(runtype))
   case ('INITIAL')
      do ifm=1,ngrids
         call ed_init_budget(edgrid_g(ifm),.true.)
         call ed_init_viable(edgrid_g(ifm))
      end do
   case ('HISTORY')
      new_day         = current_time%time < dtlsm
      do ifm=1,ngrids
         call flag_stable_cohorts(edgrid_g(ifm),.true.)
         call ed_init_viable(edgrid_g(ifm))      
         if (new_day) then
            call zero_litter_inputs(edgrid_g(ifm))
         end if
      end do
   end select
   !---------------------------------------------------------------------------------------!


   if (ifoutput /= 0) call h5_output('INST')

   if (isoutput /= 0) then
      select case (trim(runtype))
      case ('INITIAL')
         call h5_output('HIST')
      case ('HISTORY')
         call h5_output('CONT')
      end select
   end if

   if (writing_year ) call h5_output('YEAR')

   !----- Start the timesteps. ------------------------------------------------------------!
   if (mynum == 1) write(unit=*,fmt='(a)') ' === Time integration starts (model) ==='


   timestep: do while (time < timmax)

      istp = istp + 1

      !------------------------------------------------------------------------------------!
      !   CPU timing information & model timing information.                               !
      !------------------------------------------------------------------------------------!
      call timing(1,t1)
      wtime1=walltime(wtime_start)
      !------------------------------------------------------------------------------------!

      if (current_time%time < dtlsm .and. mynum == 1) then
           write (unit=*,fmt='(a,3x,2(i2.2,a),i4.4,a,3(i2.2,a))')                          &
              ' - Simulating:',current_time%month,'/',current_time%date,'/'                &
                              ,current_time%year,' ',current_time%hour,':'                 &
                              ,current_time%min,':',current_time%sec,' UTC'
      end if


      !----- Define which cohorts are to be solved prognostically. ------------------------!
      do ifm=1,ngrids
         call flag_stable_cohorts(edgrid_g(ifm),.false.)
      end do
      !------------------------------------------------------------------------------------!


      !----- Solve the radiation profile. -------------------------------------------------!
      do ifm=1,ngrids
          call canopy_radiation(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     At this point, all meteorologic driver data for the land surface model has     !
      ! been updated for the current timestep.  Perform the time average for the output    !
      ! diagnostic.                                                                        !
      !------------------------------------------------------------------------------------!
      do ifm=1,ngrids
         call integrate_ed_fmean_met_vars(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!


      !----- Solve the photosynthesis and biophysics. -------------------------------------!
      select case (integration_scheme)
      case (0)
         do ifm=1,ngrids
            call euler_timestep(edgrid_g(ifm))
         end do
      case (1)
         do ifm=1,ngrids
            call rk4_timestep(edgrid_g(ifm))
         end do
      case (2)
         do ifm=1,ngrids
            call heun_timestep(edgrid_g(ifm))
         end do
      case (3)
         do ifm=1,ngrids
            call hybrid_timestep(edgrid_g(ifm))
         end do
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the model time.                                                         !
      !------------------------------------------------------------------------------------!
      time=time+dble(dtlsm)
      call update_model_time_dm(current_time, dtlsm)
      !------------------------------------------------------------------------------------!


      !----- Check whether it is some special time... -------------------------------------!
      new_day         = current_time%time < dtlsm
      if (.not. past_one_day .and. new_day) past_one_day=.true.

      new_month       = current_time%date == 1  .and. new_day
      if (.not. past_one_month .and. new_month) past_one_month=.true.

      new_year        = current_time%month == month_yrstep .and. new_month
      mont_analy_time = new_month .and. writing_mont
      dail_analy_time = new_day   .and. writing_dail
      dcyc_analy_time = new_month .and. writing_dcyc
      reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
      annual_time     = new_year .and. writing_year
      last_step       = time >= timmax
      !------------------------------------------------------------------------------------!



      !----- Check whether this is time to write fast analysis output or not. -------------!
      select case (unitfast)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         analysis_time   = mod(current_time%time, frqfast) < dtlsm .and.                   &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. iooutput /= 0)
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      case (2)   !----- Months, analysis time is at the new month -------------------------!
         analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /=0) .and.         &
                           mod(real(12+current_time%month-imontha),frqfast) == 0.
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      case (3) !----- Year, analysis time is at the same month as initial time ------------!
         analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /= 0) .and.        &
                           current_time%month == imontha .and.                             &
                           mod(real(current_time%year-iyeara),frqfast) == 0.
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      end select
      !------------------------------------------------------------------------------------!



      !----- Check whether it is an observation time --------------------------------------!
      if (iooutput == 0 .or. unitfast /= 0) then
         !------ Observation_time is not used when unitfast /= 0 or iooutput is 0. --------!
         observation_time = .false. 
         !---------------------------------------------------------------------------------!
      else
         !------ check whether it is the observation time. --------------------------------!
         call is_obstime(current_time%year,current_time%month,current_time%date            &
                        ,current_time%time,observation_time,obstime_idx)
         !---------------------------------------------------------------------------------!

         !------ Get rid of the obstime record if observation_time is true. ---------------!
         if (observation_time) then
            call remove_obstime(obstime_idx)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !----- Check whether this is time to write restart output or not. -------------------!
      select case(unitstate)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         history_time   = mod(current_time%time, frqstate) < dtlsm .and. isoutput /= 0
      case (2)   !----- Months, history time is at the new month --------------------------!
         history_time   = new_month .and. isoutput /= 0 .and.                              &
                          mod(real(12+current_time%month-imontha),frqstate) == 0.
      case (3) !----- Year, history time is at the same month as initial time -------------!
         history_time   = new_month .and. isoutput /= 0 .and.                              &
                          current_time%month == imontha .and.                              &
                          mod(real(current_time%year-iyeara),frqstate) == 0.
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      If this is the last step, write the history even if it is not the a typical   !
      ! time step to write history.                                                        !
      !------------------------------------------------------------------------------------!
      history_time      = history_time .or. last_step
      !------------------------------------------------------------------------------------!



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



      !----- Check if this is the beginning of a new simulated day. -----------------------!
      if (new_day) then


         if (record_err) then

            open (unit=77,file=trim(errmax_fout),form='formatted',access='append'          &
                 ,status='old')
            write (unit=77,fmt=fmtcntr) current_time%year,current_time%month               &
                                       ,current_time%date,(integ_err(nn,1),nn=1,nerr)
            close(unit=77,status='keep')

            open (unit=78,file=trim(sanity_fout),form='formatted',access='append'          &
                 ,status='old')
            write (unit=78,fmt=fmtcntr) current_time%year,current_time%month               &
                                       ,current_time%date,(integ_err(nn,2),nn=1,nerr)
            close(unit=78,status='keep')

            call reset_integ_err()
         end if


         !----- Find the number of days in this month and the previous month. -------------!
         call yesterday_info(current_time,daybefore,dbndays,dbndaysi)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     This cap limits the growth rate depending on the day of the month.  This is !
         ! to account for the different time scales between heartwood growth (monthly) and !
         ! growth of the other tissues (daily).  This factor ensures that growth           !
         ! respiration is evenly distributed during the month, as opposed to have a spike  !
         ! on the second day of every month, when live tissues are growing to catch up the !
         ! allometry after heartwood biomass had increased.  This factor grows as the time !
         ! step approaches the end of the month, but the biomass increment will be the     !
         ! same every day and trees will be back on allometry be the time of the following !
         ! month in case storage is not limiting.
         !---------------------------------------------------------------------------------!
         gr_tfact0 = 1.0 / (dbndays - daybefore%date + 1)
         !------------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Compute phenology, growth, mortality, recruitment, disturbance, and check   !
         ! whether we will apply them to the ecosystem or not.                             !
         !---------------------------------------------------------------------------------!
         call veg_dynamics_driver(new_month,new_year,gr_tfact0,veget_dyn_on)
         !---------------------------------------------------------------------------------!

         !----- First day of a month. -----------------------------------------------------!
         if (new_month) then

            !------------------------------------------------------------------------------!
            !      On the monthly timestep we have performed various fusion/fission calls. !
            ! Therefore the var-table's pointer vectors must be updated, and the global    !
            ! definitions of the total numbers must be exported to all nodes.              !
            !      Also, if we do not need to fill the tables until we do I/O, so instead  !
            ! of running this routine every time the demographics change, we set this flag !
            ! and run the routine when the next IO occurs.                                 !
            !------------------------------------------------------------------------------!
            if (nnodetot > 1) then
               if (mynum == 1) write(unit=*,fmt='(a)')                                     &
                                               '-- Monthly node synchronization - waiting'
               if (whos_slow ) then
                  write(unit=*,fmt='(a,1x,i5,1x,a,1x,f7.1)') 'Node',mynum                  &
                                                            ,'time', walltime(wtime_start)
               end if
#if defined(RAMS_MPI)
               call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
               if (mynum == 1) write(unit=*,fmt='(a)') '-- Synchronized.'
            end if

            filltables=.true.   ! call filltab_alltypes

            !----- Read new met driver files only if this is the first timestep. ----------!
            call read_met_drivers()
            !----- Re-allocate integration buffer. ----------------------------------------!
            call initialize_rk4patches(.false.)
         end if
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Update the yearly variables.                                                  !
      !------------------------------------------------------------------------------------!
      if (analysis_time .and. new_year .and. new_day) then
         do ifm = 1,ngrids
            call update_ed_yearly_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the model output driver.                                                  !
      !------------------------------------------------------------------------------------!
      call ed_output(observation_time,analysis_time,new_day,new_year,dail_analy_time       &
                    ,mont_analy_time,dcyc_analy_time,annual_time,history_time,dcycle_time)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Write a file with the current history time.                                    !
      !------------------------------------------------------------------------------------!
      if (history_time .and. mynum == 1) then
         open (unit=18,file=trim(restore_file),form='formatted',status='replace'           &
              ,action='write')
         write(unit=18,fmt=fmtrest) current_time%year,current_time%month,current_time%date &
                                   ,current_time%hour,current_time%min
         close(unit=18,status='keep')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Reset time happens every frqsum.  This is to avoid variables to build up when !
      ! history and analysis are off.  This should be done outside ed_output so I have a   !
      ! chance to copy some of these to BRAMS structures.                                  !
      !------------------------------------------------------------------------------------!
      if (reset_time) then
         do ifm=1,ngrids
            call zero_ed_fmean_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Reset inputs to soil carbon.                                                  !
      !------------------------------------------------------------------------------------!
      if (new_day) then
         do ifm=1,ngrids
            call zero_litter_inputs(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update the meteorological driver, and the hydrology parameters.               !
      !------------------------------------------------------------------------------------!
      do ifm=1,ngrids
         call update_met_drivers(edgrid_g(ifm))
      end do
      if (new_day .and. new_month) then
         do ifm = 1,ngrids
            call updateHydroParms(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update the yearly variables.                                                  !
      !------------------------------------------------------------------------------------!
      !if (analysis_time .and. new_month .and. new_day .and. current_time%month == 6) then
      !   do ifm = 1,ngrids
      !      call update_ed_yearly_vars(edgrid_g(ifm))
      !   end do
      !end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update lateral hydrology.                                                     !
      !------------------------------------------------------------------------------------!
      call calcHydroSubsurface()
      call calcHydroSurface()
      call writeHydro()
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update wall time.                                                             !
      !------------------------------------------------------------------------------------!
      wtime2=walltime(wtime_start)
      call timing(2,t2)
      if (printbanner) then
         write (unit=*,fmt='(a,i10,a,i2.2,a,i2.2,a,i4.4,a,f6.0,2(a,f7.3),a)')              &
             ' Timestep ',istp,'; Sim time  '                                              &
            ,current_time%month,'-',current_time%date,'-',current_time%year,' '            &
            ,current_time%time,'s; Wall',wtime2-wtime1,'s; CPU',t2-t1,'s'
      end if
   end do timestep
   !---------------------------------------------------------------------------------------!

   wtime_tot=walltime(wtime_start)
   write(unit=*,fmt='(a,1x,f10.1,1x,a)') ' === Time integration ends; Total elapsed time=' &
                                        ,wtime_tot," ==="
   return
end subroutine ed_model
!==========================================================================================!
!==========================================================================================!
