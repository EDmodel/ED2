!==========================================================================================!
!==========================================================================================!
!    This file contains the main driver for ED2.                                           !
!------------------------------------------------------------------------------------------!
subroutine ed_model()
  
   use ed_misc_coms  , only : integration_scheme  & ! intent(in)
                            , current_time        & ! intent(in)
                            , frqfast             & ! intent(in)
                            , frqstate            & ! intent(in)
                            , out_time_fast       & ! intent(in)
                            , dtlsm               & ! intent(in)
                            , ifoutput            & ! intent(in)
                            , isoutput            & ! intent(in)
                            , idoutput            & ! intent(in)
                            , imoutput            & ! intent(in)
                            , iqoutput            & ! intent(in)
                            , iyoutput            & ! intent(in)
                            , itoutput            & ! intent(in)
                            , frqsum              & ! intent(in)
                            , unitfast            & ! intent(in)
                            , unitstate           & ! intent(in)
                            , imontha             & ! intent(in)
                            , iyeara              & ! intent(in)
                            , outstate            & ! intent(in)
                            , outfast             & ! intent(in)
                            , nrec_fast           & ! intent(in)
                            , nrec_state          & ! intent(in)
                            , ffilout             & ! intent(in)
                            , runtype
   use ed_misc_coms  , only : outputMonth         & ! intent(in)
                            , fast_diagnostics    ! ! intent(in)
   use grid_coms     , only : ngrids              & ! intent(in)
                            , istp                & ! intent(in)
                            , time                & ! intent(in)
                            , timmax              & ! intent(in)
                            , nnxp                & ! intent(in)
                            , nnyp                & ! intent(in)
                            , nzs                 & ! intent(in)
                            , nzg                 ! ! intent(in)
   use ed_state_vars , only : edgrid_g            & ! intent(in)
                            , edtype              & ! intent(in)
                            , patchtype           & ! intent(in)
                            , filltab_alltypes    & ! intent(in)
                            , filltables          ! ! intent(in)
   use rk4_driver    , only : rk4_timestep        ! ! intent(in)
   use rk4_coms      , only : checkbudget         & ! intent(in)
                            , integ_err           & ! intent(in)
                            , integ_lab           & ! intent(in)
                            , record_err          & ! intent(inout)
                            , print_detailed      & ! intent(inout)
                            , nerr                & ! intent(in)
                            , errmax_fout         & ! intent(in)
                            , sanity_fout         & ! intent(in)
                            , alloc_integ_err     & ! subroutine
                            , assign_err_label    & ! subroutine
                            , reset_integ_err     ! ! subroutine
   use ed_node_coms  , only : mynum               & ! intent(in)
                            , nnodetot            ! ! intent(in)
   use disturb_coms  , only : include_fire        ! ! intent(in)
   use mem_polygons  , only : n_ed_region         & ! intent(in)
                            , n_poi               & ! intent(in)
                            , maxpatch            & ! intent(in)
                            , maxcohort           ! ! intent(in)
   use consts_coms   , only : day_sec             ! ! intent(in)
   implicit none
   !----- Common blocks. ------------------------------------------------------------------!
   include 'mpif.h'
   !----- Local variables. ----------------------------------------------------------------!
   character(len=28)  :: fmthead
   character(len=32)  :: fmtcntr
   integer            :: ifm
   integer            :: i
   integer            :: ierr
   integer            :: nn
   integer            :: ndays
   logical            :: analysis_time
   logical            :: new_day
   logical            :: new_month
   logical            :: new_year
   logical            :: the_end
   logical            :: writing_dail
   logical            :: writing_mont
   logical            :: writing_dcyc
   logical            :: writing_year
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
   real               :: wtime_start
   real               :: t1
   real               :: wtime1
   real               :: wtime2
   real               :: t2
   real               :: wtime_tot
   !----- Local constants. ----------------------------------------------------------------!
   logical         , parameter :: whos_slow=.false. ! This will print out node numbers 
                                                    !    during synchronization, so you 
                                                    !    can find out which node is the 
                                                    !    slow one
   !----- External functions. -------------------------------------------------------------!
   integer , external :: julday   ! Get the elapsed # of days since ED time origin.
   real    , external :: walltime ! Wall time
   integer , external :: num_days ! Number of days in the current month
   !---------------------------------------------------------------------------------------!

   past_one_day   = .false.
   past_one_month = .false.
   filltables     = .false.
   
   !----- Print the hour banner only for regional runs that aren't massively parallel. ----!
   printbanner = n_ed_region > 0 .and. edgrid_g(1)%npolygons > 50 .and. mynum == 1


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

   writing_dail      = idoutput > 0
   writing_mont      = imoutput > 0
   writing_dcyc      = iqoutput > 0
   writing_year      = iyoutput > 0
   out_time_fast     = current_time
   out_time_fast%month = -1

   !---------------------------------------------------------------------------------------!
   !     Checking if the user has indicated a need for any of the fast flux diagnostic     !
   ! variables, these are used in conditions of ifoutput,idoutput and imoutput conditions. !
   ! If they are not >0, then set the logical, fast_diagnostics to false.                  !
   !---------------------------------------------------------------------------------------!
   fast_diagnostics = checkbudget   .or. ifoutput /= 0 .or. idoutput /= 0 .or.             &
                      imoutput /= 0 .or. iqoutput /= 0 .or. itoutput /= 0

   !---------------------------------------------------------------------------------------!
   !      If this is not a history restart, then zero out the long term diagnostics.       !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      
      if (writing_mont .or. writing_dcyc) then
         do ifm=1,ngrids
            call zero_ed_monthly_output_vars(edgrid_g(ifm))
            call zero_ed_daily_output_vars(edgrid_g(ifm))
         end do
      elseif (writing_dail) then
         do ifm=1,ngrids
            call zero_ed_daily_output_vars(edgrid_g(ifm))
         end do
      end if

      !----- Output initial state. --------------------------------------------------------!
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
      

   endif
   
   !---------------------------------------------------------------------------------------!
   !    Allocate memory to the integration patch, Euler now utilises the RK4 buffers too.  !
   !---------------------------------------------------------------------------------------!
   call initialize_rk4patches(.true.)
   !---------------------------------------------------------------------------------------!

   do ifm=1,ngrids
      call reset_averaged_vars(edgrid_g(ifm))
   end do
   
   
   if (ifoutput /= 0) call h5_output('INST')
   if (isoutput /= 0) call h5_output('HIST')
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
           write (unit=*,fmt='(a,3x,2(i2.2,a),i4.4,a,f8.2)')                               &
              ' - Simulating:',current_time%month,'/',current_time%date,'/'                &
                              ,current_time%year,' ',current_time%time
      end if


      !----- Define which cohorts are to be solved prognostically. ------------------------!
      do ifm=1,ngrids
         call flag_stable_cohorts(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!


      !----- Solve the radiation profile. -------------------------------------------------!
      do ifm=1,ngrids
          call radiate_driver(edgrid_g(ifm))
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
            call rk4_timestep(edgrid_g(ifm),ifm)
         end do
      case (2)
         do ifm=1,ngrids
            call heun_timestep(edgrid_g(ifm))
         end do
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the daily averages if daily or monthly analysis are needed.             !
      !------------------------------------------------------------------------------------!
      if (writing_dail .or. writing_mont .or. writing_dcyc) then
         do ifm=1,ngrids
            call integrate_ed_daily_output_state(edgrid_g(ifm))
         end do
      end if
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

      new_year        = current_time%month == 1 .and. new_month
      mont_analy_time = new_month .and. writing_mont
      dail_analy_time = new_day   .and. writing_dail
      dcyc_analy_time = new_month .and. writing_dcyc
      reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
      the_end         = mod(time,timmax) < dble(dtlsm)
      annual_time     = new_month .and. writing_year .and.                                 &
                        current_time%month == outputMonth

      !----- Check whether this is time to write fast analysis output or not. -------------!
      select case (unitfast)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         analysis_time   = mod(current_time%time, frqfast) < dtlsm .and.                   &
                           (ifoutput /= 0 .or. itoutput /=0)
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
      !    Update nrec_fast and nrec_state if it is a new month and outfast/outstate are   !
      ! monthly and frqfast/frqstate are daily or by seconds.                              !
      !------------------------------------------------------------------------------------!
      if (new_month) then
         ndays=num_days(current_time%month,current_time%year)
         if (outfast  == -2.) nrec_fast  = ndays*ceiling(day_sec/frqfast)
         if (outstate == -2.) nrec_state = ndays*ceiling(day_sec/frqstate)
      end if

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

         !----- Do phenology, growth, mortality, recruitment, disturbance. ----------------!
         call vegetation_dynamics(new_month,new_year)

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
               call MPI_Barrier(MPI_COMM_WORLD,ierr)
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
      !     Call the model output driver.                                                  !
      !------------------------------------------------------------------------------------!
      call ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,dcyc_analy_time &
                    ,annual_time,writing_dail,writing_mont,writing_dcyc,history_time       &
                    ,dcycle_time,the_end)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Reset time happens every frqsum.  This is to avoid variables to build up when !
      ! history and analysis are off.  This should be done outside ed_output so I have a   !
      ! chance to copy some of these to BRAMS structures.                                  !
      !------------------------------------------------------------------------------------!
      if (reset_time) then
         do ifm=1,ngrids
            call reset_averaged_vars(edgrid_g(ifm))
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
      if (analysis_time .and. new_month .and. new_day .and. current_time%month == 6) then
         do ifm = 1,ngrids
            call update_ed_yearly_vars(edgrid_g(ifm))
         end do
      end if
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
