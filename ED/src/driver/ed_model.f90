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
                            , integ_err           & ! intent(in)
                            , record_err          & ! intent(in)
                            , err_label           & ! intent(in)
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
   use rk4_coms      , only : checkbudget         ! ! intent(in)
   use ed_node_coms  , only : mynum               & ! intent(in)
                            , nnodetot            ! ! intent(in)
   use disturb_coms  , only : include_fire        ! ! intent(in)
   use mem_sites     , only : n_ed_region         & ! intent(in)
                            , maxpatch            & ! intent(in)
                            , maxcohort           ! ! intent(in)
   use consts_coms   , only : day_sec             ! ! intent(in)
   implicit none
   !----- Common blocks. ------------------------------------------------------------------!
   include 'mpif.h'
   !----- Local variables. ----------------------------------------------------------------!
   character(len=10)  :: c0
   character(len=512) :: integ_fname
   integer            :: ifm,i
   integer            :: ierr
   integer            :: ndays
   logical            :: analysis_time, new_day, new_month, new_year, the_end
   logical            :: writing_dail,writing_mont,writing_year,history_time, annual_time
   logical            :: mont_analy_time,dail_analy_time,reset_time
   logical            :: past_one_day,past_one_month
   logical            :: printbanner
   real               :: wtime_start,t1,wtime1,wtime2,t2,wtime_tot
   !----- Local constants. ----------------------------------------------------------------!
   character(len=*), parameter :: h="**(model)**"   ! Header.
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
   if(record_err) then
      integ_err = 0_8
      integ_fname = "integrator.log"
      open  (unit=77,file=trim(integ_fname),form="formatted",status="replace")     
      write (unit=77,fmt='(a)') "num  name  ERMAX  IFLAG"
      close (unit=77,status='keep')
   end if
   writing_dail      = idoutput > 0
   writing_mont      = imoutput > 0
   writing_year      = iyoutput > 0
   out_time_fast     = current_time
   out_time_fast%month = -1

   !---------------------------------------------------------------------------------------!
   !     Checking if the user has indicated a need for any of the fast flux diagnostic
   ! variables, these are used in conditions of ifoutput,idoutput and imoutput conditions.
   ! If they are not >0, then set the logical, fast_diagnostics to false.
   !---------------------------------------------------------------------------------------!
   fast_diagnostics = checkbudget .or. ifoutput /= 0 .or. idoutput /= 0 .or. imoutput /= 0 .or. itoutput /= 0

   !------------------------------------------------------------------------!
   ! If this is not a history restart - then zero out the
   ! long term diagnostics
   !------------------------------------------------------------------------!
   if (trim(runtype) .ne. 'HISTORY') then
      
      if (writing_mont) then
         do ifm=1,ngrids
            call zero_ed_monthly_output_vars(edgrid_g(ifm))
            call zero_ed_daily_output_vars(edgrid_g(ifm))
         end do
      elseif (writing_dail) then
         do ifm=1,ngrids
            call zero_ed_daily_output_vars(edgrid_g(ifm))
         end do
      end if

      !!Output Initial State
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
      

   endif
   
   !    Allocate memory to the integration patch

   if (integration_scheme == 1 .or. integration_scheme == 2) then
      call initialize_rk4patches(1)
   end if

   do ifm=1,ngrids
      call reset_averaged_vars(edgrid_g(ifm))
   end do
   
   
   if (writing_year) call h5_output('YEAR')

   !         Start the timesteps

   if ( mynum == 1) then
      write(*,"(/,a,/)") " === Time integration starts (model) ==="
   endif


   timestep: do while (time < timmax)
      
      istp = istp + 1
      !     begtime=time
      
      !   CPU timing information & model timing information
      !   ===================================================

      call timing(1,t1)

      wtime1=walltime(wtime_start)

      if(current_time%time < dtlsm .and. mynum == 1)  &
           write(*,'(a,3x,2(i2.2,a),i4.4,a,f8.2)')'Simulating:',current_time%month,'/',  &
           current_time%date,'/',current_time%year,' ',current_time%time

      do ifm=1,ngrids
         call flag_stable_cohorts(edgrid_g(ifm))
      end do

      do ifm=1,ngrids
          call radiate_driver(edgrid_g(ifm))
      end do

      ! THEN, DO THE PHOTOSYNTHESIS AND BIOPHYSICS.
      select case (integration_scheme)
      case (0)
         do ifm=1,ngrids
            call euler_timestep(edgrid_g(ifm))
         end do
      case (1,2)
         do ifm=1,ngrids
            call rk4_timestep(edgrid_g(ifm),ifm)
         end do
      end select
      
      !-------------------------------------------------------------------!
      ! Update the daily averages if daily or monthly analysis are needed !
      !-------------------------------------------------------------------!
      if (writing_dail .or. writing_mont) then
         do ifm=1,ngrids
            call integrate_ed_daily_output_state(edgrid_g(ifm))
         end do
      end if
      
      time=time+dble(dtlsm)
      call update_model_time_dm(current_time, dtlsm)

      !----- Checking whether it is some special time... -----------------------------------!
      new_day         = current_time%time < dtlsm
      if (.not. past_one_day .and. new_day) past_one_day=.true.
      
      new_month       = current_time%date == 1  .and. new_day
      if (.not. past_one_month .and. new_month) past_one_month=.true.

      new_year        = current_time%month == 1 .and. new_month
      mont_analy_time = new_month .and. writing_mont
      dail_analy_time = new_day   .and. writing_dail
      reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
      the_end         = mod(time,timmax) < dble(dtlsm)
      annual_time     = new_month .and. writing_year .and. current_time%month == outputMonth

      !----- Checking whether this is time to write fast analysis output or not. -----------!
      select case (unitfast)
      case (0,1) !----- Now both are in seconds --------------------------------------------!
         analysis_time   = mod(current_time%time, frqfast) < dtlsm .and. (ifoutput /= 0 .or. itoutput /=0)
      case (2)   !----- Months, analysis time is at the new month --------------------------!
         analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /=0) .and.&
                           mod(real(12+current_time%month-imontha),frqfast) == 0.
      case (3) !----- Year, analysis time is at the same month as initial time -------------!
         analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /= 0) .and.         &
                           current_time%month == imontha .and.                              &
                           mod(real(current_time%year-iyeara),frqfast) == 0.
      end select

      !----- Checking whether this is time to write restart output or not. -----------------!
      select case(unitstate)
      case (0,1) !----- Now both are in seconds --------------------------------------------!
         history_time   = mod(current_time%time, frqstate) < dtlsm .and. isoutput /= 0
      case (2)   !----- Months, history time is at the new month ---------------------------!
         history_time   = new_month .and. isoutput /= 0 .and.                               &
                          mod(real(12+current_time%month-imontha),frqstate) == 0.
      case (3) !----- Year, history time is at the same month as initial time --------------!
         history_time   = new_month .and. isoutput /= 0 .and.                               &
                          current_time%month == imontha .and.                               &
                          mod(real(current_time%year-iyeara),frqstate) == 0.
      end select

      !-------------------------------------------------------------------------------------!
      !    Updating nrec_fast and nrec_state if it is a new month and outfast/outstate are  !
      ! monthly and frqfast/frqstate are daily or by seconds.                               !
      !-------------------------------------------------------------------------------------!
      if (new_month) then
         ndays=num_days(current_time%month,current_time%year)
         if (outfast  == -2.) nrec_fast  = ndays*ceiling(day_sec/frqfast)
         if (outstate == -2.) nrec_state = ndays*ceiling(day_sec/frqstate)
      end if

      ! Check if this is the beginning of a new simulated day.
      if(new_day)then
         if(record_err) then
            open(unit=77,file=trim(integ_fname),form="formatted",access="append",status="old")
            do i = 1,46
               if(sum(integ_err(i,1:2)) .gt. 0_8)then                 
                  write(unit=77,fmt='(2(i4,1x),a,2(1x,i7))') mynum,i,trim(err_label(i)),integ_err(i,1:2)
!                  print*,i,trim(err_label(i)),integ_err(i,1:2)
               endif
            enddo
            close(unit=77,status='keep')
            integ_err = 0_8
         endif

         ! Do phenology, growth, mortality, recruitment, disturbance.
         call vegetation_dynamics(new_month,new_year)

         ! First day of a month.
         if(new_month)then

            ! On the monthly timestep we have performed various
            ! fusion/fission calls. Therefore the var-table's pointer
            ! vectors must be updated, and the global definitions
            ! of the total numbers must be exported to all nodes
            
            ! If maxpatch and maxcohort are both negative, the number of patches and 
            ! cohorts remain the same throughout the run, no need to call it.

            ! Also, if we do not need to fill the tables until we do I/O, so instead of
            ! running this routine every time the demographics change, we set this flag
            ! and run the routine when the next IO occurs.

            if (nnodetot>1) then
               if ( mynum == 1) write(*,"(a)")'-- Monthly node synchronization - waiting'
               if ( whos_slow ) write(*,"(a,i5,a,f7.1)")'Node',mynum,' time', walltime(wtime_start)
               call MPI_Barrier(MPI_COMM_WORLD,ierr)
               if ( mynum == 1) write(*,"(a)")'-- Synchronized'
            endif

            if (maxcohort >= 0 .or. maxpatch >= 0) filltables=.true.   ! call filltab_alltypes

            ! Read new met driver files only if this is the first timestep 
            call read_met_drivers()
            
            ! Re-allocate integration buffer
            if (integration_scheme == 1 .or. integration_scheme == 2) then
               call initialize_rk4patches(0)
            end if
         endif
         
      endif

      !   Call the model output driver 
      !   ====================================================
      call ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,annual_time &
                    ,writing_dail,writing_mont,history_time,the_end)

      ! Reset time happens every frqsum. This is to avoid variables to build up when
      ! history and analysis are off. Put outside ed_output so I have a chance to copy
      ! some of these to BRAMS structures.
      if(reset_time) then    
         do ifm=1,ngrids
            call reset_averaged_vars(edgrid_g(ifm))
         end do
      end if

      do ifm=1,ngrids
         call update_met_drivers(edgrid_g(ifm))
      end do
      if(new_day .and. new_month)then
         ! Loop all grids
         do ifm = 1,ngrids
            call updateHydroParms(edgrid_g(ifm))
         end do
      endif
      

      if(analysis_time)then
         if(new_month .and. new_day)then
            if(current_time%month == 6)then
               do ifm = 1,ngrids
                     
                  call update_ed_yearly_vars(edgrid_g(ifm))
                  enddo
                  !call zero_ed_yearly_vars(polygon_list_g(1)%first_polygon)
            endif
         endif
      endif

      !!Update Lateral Hydrology
      call calcHydroSubsurface()
      call calcHydroSurface()
      call writeHydro()


      wtime2=walltime(wtime_start)
      call TIMING(2,T2)
      
      if(printbanner) then
         write(*,"(a,i10,a,i2.2,a,i2.2,a,i4.4,a,f6.0,2(a,f7.3),a)") &
              ' Timestep ',istp,&
              '; Sim time  ',current_time%month, &
              '-',current_time%date,           &
              '-',current_time%year,           &
              ' ',current_time%time,           &         
              's; Wall',wtime2-wtime1,&
              's; CPU',t2-t1,&
              's'
      endif



   end do timestep
   

   wtime_tot=walltime(wtime_start)
   write(c0,"(f10.1)") wtime_tot
   write(*,"(/,a,/)") " === Time integration ends; Total elapsed time="//&
        &trim(adjustl(c0))//" ===" 
   

   return
end subroutine ed_model
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_model_time_dm(ctime,dtlong)

   use ed_misc_coms, only: simtime
   use consts_coms, only : day_sec
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
subroutine vegetation_dynamics(new_month,new_year)

  ! NB:  (1) Each subroutine has its own loops over polygons and sites.  Once
  ! we switch to arrays, this will improve vectorization.
  !      (2) Do not change the order of the subroutine calls below.  The 
  !          calculation of budgets depends on the order.

  use ed_node_coms,only:mynum,nnodetot
  use grid_coms, only: ngrids
  use ed_misc_coms, only: current_time, dtlsm,frqsum,ied_init_mode
  use disturb_coms, only: include_fire
  use disturbance_utils, only: apply_disturbances, site_disturbance_rates
  use fuse_fiss_utils, only : fuse_patches
  use ed_state_vars,only : edgrid_g,edtype
  use growth_balive,only : dbalive_dt, dbalive_dt_eq_0
  use consts_coms, only : day_sec,yr_day
  use mem_sites, only: maxpatch

  implicit none
  include 'mpif.h'
  logical, intent(in)   :: new_month,new_year
  integer               :: doy
  integer, external     :: julday
  real                  :: tfact1,tfact2
  integer               :: ifm
  type(edtype), pointer :: cgrid
  integer               :: ierr
  logical, save         :: first_time = .true.

  ! find the day of year
  doy = julday(current_time%month, current_time%date, current_time%year)
  
  ! Time factor for normalizing daily variables updated on the DTLSM step.
  tfact1 = dtlsm / day_sec

  ! Time factor for averaging dailies 
  tfact2 = 1.0 / yr_day

  !! Apply Events
  call prescribed_event(current_time%year,doy)

  
  do ifm=1,ngrids

     cgrid => edgrid_g(ifm) 
!     write (unit=*,fmt='(a)') '~~~ Normalize_ed_daily_vars...'
     call normalize_ed_daily_vars(cgrid, tfact1)
     
!     write (unit=*,fmt='(a)') '~~~ Phenology_driver...'
     if (ied_init_mode == -8) then
        call phenology_driver_eq_0(cgrid,doy,current_time%month, tfact1)
     else
        call phenology_driver(cgrid,doy,current_time%month, tfact1)
     end if
     
!     write (unit=*,fmt='(a)') '~~~ Dbalive_dt...'
     if (ied_init_mode == -8) then
        call dbalive_dt_eq_0(cgrid,tfact2)
     else
        call dbalive_dt(cgrid,tfact2)
     end if

     if (new_month) then

!        write (unit=*,fmt='(a)') '^^^ Update_workload...'
        call update_workload(cgrid)

!        write (unit=*,fmt='(a)') '^^^ Structural_growth...'
        if (ied_init_mode == -8) then
           call structural_growth_eq_0(cgrid, current_time%month)
        else
           call structural_growth(cgrid, current_time%month)
        end if


!        write (unit=*,fmt='(a)') '^^^ Reproduction...'
        call reproduction(cgrid,current_time%month)

        if(include_fire /= 0) then
!           write (unit=*,fmt='(a)') '^^^ Fire_frequency...'
           call fire_frequency(current_time%month,cgrid)
        end if

!        write (unit=*,fmt='(a)') '^^^ Site_disturbance_rates...'
        call site_disturbance_rates(current_time%month,   &
             current_time%year, cgrid)

        if(new_year) then

!           write (unit=*,fmt='(a)') '### Apply_disturbances...'
           call apply_disturbances(cgrid)
        end if
        
     end if

     !     write (unit=*,fmt='(a)') '~~~ Update_C_and_N_pools...'
     call update_C_and_N_pools(cgrid)
     
     
     !  write (unit=*,fmt='(a)') '~~~ Zero_ed_daily_vars...'
     call zero_ed_daily_vars(cgrid)
     

     ! Fuse patches last, after all updates have been applied.  This reduces
     ! the number of patch variables that actually need to be fused.  
     if(new_year) then
!        write (unit=*,fmt='(a)') '### Fuse_patchesar...'
        if (maxpatch >= 0) call fuse_patches(cgrid,ifm)
        first_time =.false.
     end if

     ! Recalculate the agb and basal area at the polygon level
!     write (unit=*,fmt='(a)') '~~~ Update_polygon_derived_props...'
     call update_polygon_derived_props(cgrid)

!     write (unit=*,fmt='(a)') '~~~ Print_C_and_N_budgets...'
     call print_C_and_N_budgets(cgrid)
  end do
  return
end subroutine vegetation_dynamics
!==========================================================================================!
!==========================================================================================!
