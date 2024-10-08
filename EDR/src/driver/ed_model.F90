!==========================================================================================!
!==========================================================================================!
!    This file contains the main driver for ED2.                                           !
!------------------------------------------------------------------------------------------!
subroutine ed_model()
   use ed_met_driver
   use heun_driver
   use euler_driver
   use update_derived_props_module
   use ism_hyd
   use stable_cohorts
   use rk4_integ_utils
   use radiate_driver_module
   use hybrid_driver
  
   use ed_misc_coms  , only : ivegt_dynamics              & ! intent(in)
                            , integration_scheme          & ! intent(in)
                            , current_time                & ! intent(in)
                            , frqfast                     & ! intent(in)
                            , frqstate                    & ! intent(in)
                            , out_time_fast               & ! intent(in)
                            , dtlsm                       & ! intent(in)
                            , ifoutput                    & ! intent(in)
                            , isoutput                    & ! intent(in)
                            , idoutput                    & ! intent(in)
                            , imoutput                    & ! intent(in)
                            , iqoutput                    & ! intent(in)
                            , iyoutput                    & ! intent(in)
                            , itoutput                    & ! intent(in)
                            , frqsum                      & ! intent(in)
                            , unitfast                    & ! intent(in)
                            , unitstate                   & ! intent(in)
                            , imontha                     & ! intent(in)
                            , iyeara                      & ! intent(in)
                            , outstate                    & ! intent(in)
                            , outfast                     & ! intent(in)
                            , nrec_fast                   & ! intent(in)
                            , nrec_state                  & ! intent(in)
                            , ffilout                     & ! intent(in)
                            , runtype                     ! ! intent(in)
   use ed_misc_coms  , only : outputMonth                 & ! intent(in)
                            , fast_diagnostics            & ! intent(in)
                            , writing_dail                & ! intent(in)
                            , writing_mont                & ! intent(in)
                            , writing_dcyc                & ! intent(in)
                            , writing_eorq                & ! intent(in)
                            , writing_long                & ! intent(in)
                            , writing_year                ! ! intent(in)
   use grid_coms     , only : ngrids                      & ! intent(in)
                            , istp                        & ! intent(in)
                            , time                        & ! intent(in)
                            , timmax                      & ! intent(in)
                            , nnxp                        & ! intent(in)
                            , nnyp                        & ! intent(in)
                            , nzs                         & ! intent(in)
                            , nzg                         ! ! intent(in)
   use ed_state_vars , only : edgrid_g                    & ! intent(in)
                            , edtype                      & ! intent(in)
                            , patchtype                   & ! intent(in)
                            , filltab_alltypes            & ! intent(in)
                            , filltables                  ! ! intent(in)
   use rk4_driver    , only : rk4_timestep                ! ! intent(in)
   use rk4_coms      , only : checkbudget                 & ! intent(in)
                            , integ_err                   & ! intent(in)
                            , integ_lab                   & ! intent(in)
                            , record_err                  & ! intent(inout)
                            , print_detailed              & ! intent(inout)
                            , nerr                        & ! intent(in)
                            , errmax_fout                 & ! intent(in)
                            , sanity_fout                 & ! intent(in)
                            , alloc_integ_err             & ! subroutine
                            , assign_err_label            & ! subroutine
                            , reset_integ_err             ! ! subroutine
   use ed_node_coms  , only : mynum                       & ! intent(in)
                            , nnodetot                    ! ! intent(in)
   use mem_polygons  , only : n_ed_region                 & ! intent(in)
                            , n_poi                       & ! intent(in)
                            , maxpatch                    & ! intent(in)
                            , maxcohort                   ! ! intent(in)
   use consts_coms   , only : day_sec                     ! ! intent(in)
   use average_utils , only : update_ed_yearly_vars       & ! sub-routine
                            , zero_ed_dmean_vars          & ! sub-routine
                            , zero_ed_mmean_vars          & ! sub-routine
                            , zero_ed_qmean_vars          & ! sub-routine
                            , zero_ed_fmean_vars          & ! sub-routine
                            , integrate_ed_fmean_met_vars & ! sub-routine
                            , zero_ed_yearly_vars         ! ! sub-routine
#if defined(RAMS_MPI)
   use mpi
#endif
   implicit none
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

   out_time_fast     = current_time
   out_time_fast%month = -1

   !---------------------------------------------------------------------------------------!
   !      If this is not a history restart, then zero out the long term diagnostics.       !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      
      do ifm=1,ngrids
         if (writing_long) call zero_ed_dmean_vars(edgrid_g(ifm))
         if (writing_eorq) call zero_ed_mmean_vars(edgrid_g(ifm))
         if (writing_dcyc) call zero_ed_qmean_vars(edgrid_g(ifm))
      end do

      !----- Output initial state. --------------------------------------------------------!
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
   end if
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
         call flag_stable_cohorts(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!


      !----- Solve the radiation profile. -------------------------------------------------!
      do ifm=1,ngrids
          call radiate_driver(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!

goto 100

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

   wtime_tot=walltime(wtime_start)
   write(unit=*,fmt='(a,1x,f10.1,1x,a)') ' === Time integration ends; Total elapsed time=' &
                                        ,wtime_tot," ===" 
   return

100 continue
end subroutine ed_model
!==========================================================================================!
!==========================================================================================!
