!=================================================================
!
!  Main offline model driver for ED2
!
!
!=================================================================
subroutine ed_model()
  
  use misc_coms, only: integration_scheme, current_time, frqfast, frqstate    &
                      , out_time_fast, dtlsm, ifoutput, isoutput, idoutput    &
                      , imoutput, iyoutput,frqsum,unitfast,unitstate, imontha &
                      , iyeara, outstate,outfast, nrec_fast, nrec_state
  use ed_misc_coms, only: outputMonth

  use grid_coms, only : &
       ngrids,          &
       istp,            &
       time,            &
       timmax,          &
       nnxp,            &
       nnyp,            &
       nzs,             &
       nzg
  
  use ed_state_vars,only: edgrid_g, &
       integration_buff_g,          &
       edtype,                      &
       patchtype,                   &
       filltab_alltypes
  use rk4_driver_ar,only: rk4_timestep_ar
  use ed_node_coms,only:mynum,nnodetot
  use disturb_coms, only: include_fire
  use mem_sites, only : n_ed_region
  use consts_coms, only: day_sec

  implicit none

  include 'mpif.h'

  integer :: npass,nndtflg,icm,ifm,nfeed,i
  real :: wtime_start,t1,wtime1,wtime2,t2,wtime_tot
  integer :: doy,ierr

  integer,external :: julday
  
  real, external :: walltime
  character(len=10) :: c0, c1, c2
  character(len=*), parameter :: h="**(model)**"
  character(len=512) :: header_ext
  real           :: timefac_sst
  real           :: dtwater
  real, external :: update_sst_factor
  real :: ccont,ctemp,stemp,swat,lai

  real :: tfact1
  integer :: ipa,ico
  logical :: analysis_time, new_day, new_month, new_year, the_end
  logical :: writing_dail,writing_mont,writing_year,history_time, annual_time
  logical :: mont_analy_time,dail_analy_time,reset_time
  logical :: past_one_day,past_one_month
  logical :: printbanner
  integer :: ndays
  integer, external :: num_days
  
  past_one_day   = .false.
  past_one_month = .false.
  
  ! This should keep both SOI and regional runners happy...
  printbanner = n_ed_region > 0 .and. mynum == 1
  
  wtime_start=walltime(0.)
  istp = 0
  
  writing_dail      = idoutput > 0
  writing_mont      = imoutput > 0
  writing_year      = iyoutput > 0
  out_time_fast     = current_time
  out_time_fast%month = -1

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

  !    Allocate memory to the integration patch

  if(integration_scheme == 1) call initialize_rk4patches_ar(1)

  do ifm=1,ngrids
     call reset_averaged_vars(edgrid_g(ifm))
  end do
  
  !!Output Initial State
  do ifm=1,ngrids
     call update_ed_yearly_vars_ar(edgrid_g(ifm))
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


     if (nnodetot>1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

     call timing(1,t1)

     wtime1=walltime(wtime_start)

     if(current_time%time < dtlsm .and. mynum == 1)  &
          write(*,'(a,3x,2(i2.2,a),i4.4)')'Simulating:',current_time%month,'/',  &
          current_time%date,'/',current_time%year

     do ifm=1,ngrids
         call radiate_driver_ar(edgrid_g(ifm))
     end do

     ! THEN, DO THE PHOTOSYNTHESIS AND BIOPHYSICS.
     if(integration_scheme == 0)then
        do ifm=1,ngrids
           call euler_timestep_ar(edgrid_g(ifm))
        end do
     elseif(integration_scheme == 1)then
        do ifm=1,ngrids
           call rk4_timestep_ar(edgrid_g(ifm),ifm,integration_buff_g)
        end do
     endif
     
     !-------------------------------------------------------------------!
     ! Update the daily averages if daily or monthly analysis are needed !
     !-------------------------------------------------------------------!
     if (writing_dail .or. writing_mont) then
        do ifm=1,ngrids
           call integrate_ed_daily_output_state(edgrid_g(ifm))
        end do
     end if
     
     time=time+dtlsm
     call update_model_time_dm(current_time, dtlsm)

     !----- Checking whether it is some special time... -----------------------------------!
     new_day         = current_time%time < dtlsm
     if (.not. past_one_day .and. new_day) past_one_day=.true.
     
     new_month       = current_time%date == 1  .and. new_day
     if (.not. past_one_month .and. new_month) past_one_month=.true.

     new_year        = current_time%month == 1 .and. new_month
     mont_analy_time = new_month .and. writing_mont
     dail_analy_time = new_day   .and. writing_dail
     reset_time      = mod(time,dble(frqsum)) < dtlsm
     the_end         = mod(time,timmax) < dtlsm
     annual_time     = new_month .and. writing_year .and. current_time%month == outputMonth

     !----- Checking whether this is time to write fast analysis output or not. -----------!
     select case (unitfast)
     case (0,1) !----- Now both are in seconds --------------------------------------------!
        analysis_time   = mod(current_time%time, frqfast) < dtlsm .and. ifoutput /= 0
     case (2)   !----- Months, analysis time is at the new month --------------------------!
        analysis_time   = new_month .and. ifoutput /= 0 .and.                              &
                          mod(real(12+current_time%month-imontha),frqfast) == 0.
     case (3) !----- Year, analysis time is at the same month as initial time -------------!
        analysis_time   = new_month .and. ifoutput /= 0 .and.                              &
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
        if (outfast  == -2.) nrec_fast  = ndays*day_sec/frqfast
        if (outstate == -2.) nrec_state = ndays*day_sec/frqstate
     end if

     !   Call the model output driver 
     !   ====================================================
     call ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,annual_time &
                   ,writing_dail,writing_mont,history_time,reset_time,the_end)

     ! Check if this is the beginning of a new simulated day.
     if(new_day)then

        ! Do phenology, growth, mortality, recruitment, disturbance.
        call vegetation_dynamics(new_month,new_year)

        ! First day of a month.
        if(new_month)then
           
           ! On the monthly timestep we have performed various
           ! fusion/fission calls. Therefore the var-table's pointer
           ! vectors must be updated, and the global definitions
           ! of the total numbers must be exported to all nodes
           
           call filltab_alltypes

           ! Read new met driver files only if this is the first timestep 
           call read_met_drivers_array()
           
           ! Re-allocate integration buffer
           if(integration_scheme == 1) call initialize_rk4patches_ar(0)
        endif
        
     endif

     do ifm=1,ngrids
        call update_met_drivers_array(edgrid_g(ifm))
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
                 call update_ed_yearly_vars_ar(edgrid_g(ifm))
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

   use misc_coms, only: simtime
   use consts_coms, only : day_sec
   implicit none

   type(simtime) :: ctime
   real, intent(in) :: dtlong
   logical, external :: isleap
   real, dimension(12) :: daymax
  
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

  use grid_coms, only: ngrids
  use misc_coms, only: current_time, dtlsm,frqsum
  use disturb_coms, only: include_fire
  use disturbance_utils_ar, only: apply_disturbances_ar, site_disturbance_rates_ar
  use fuse_fiss_utils_ar, only : fuse_patches_ar
  use ed_state_vars,only : edgrid_g,filltab_alltypes,edtype
  use growth_balive_ar,only : dbalive_dt_ar
  use consts_coms, only : day_sec,yr_day
  implicit none

  logical, intent(in)   :: new_month,new_year
  integer               :: doy
  integer, external     :: julday
  real                  :: tfact1,tfact2
  integer               :: ip
  integer               :: isite
  integer               :: ifm
  type(edtype), pointer :: cgrid

  ! find the day of year
  doy = julday(current_time%month, current_time%date, current_time%year)
  
  ! Time factor for normalizing daily variables updated on the DTLSM step.
  tfact1 = dtlsm / day_sec

  ! Time factor for averaging dailies 
  tfact2 = 1.0 / yr_day

  
  do ifm=1,ngrids

     cgrid => edgrid_g(ifm) 
!     write (unit=*,fmt='(a)') '~~~ Normalize_ed_daily_vars...'
     call normalize_ed_daily_vars(cgrid, tfact1)
     
!     write (unit=*,fmt='(a)') '~~~ Phenology_driver_ar...'
     call phenology_driver_ar(cgrid,doy,current_time%month, tfact1)
     
!     write (unit=*,fmt='(a)') '~~~ Dbalive_dt_ar...'
     call dbalive_dt_ar(cgrid,tfact2)
     
     
     if(new_month)then

!        write (unit=*,fmt='(a)') '^^^ Structural_growth_ar...'
        call structural_growth_ar(cgrid, current_time%month)


!        write (unit=*,fmt='(a)') '^^^ Reproduction_ar...'
        call reproduction_ar(cgrid,current_time%month)

        ! NB: FIRE CURRENTLY OCCURS AT THE SITE LEVEL. MIKE: MAYBE
        ! YOU HAVE SOME IDEAS HERE?

        if(include_fire == 1) then
!           write (unit=*,fmt='(a)') '^^^ Fire_frequency_ar...'
           call fire_frequency_ar(current_time%month,cgrid)
        end if
!        write (unit=*,fmt='(a)') '^^^ Site_disturbance_rates_ar...'
        call site_disturbance_rates_ar(current_time%month,   &
             current_time%year, cgrid)

        if(new_year) then
!           write (unit=*,fmt='(a)') '### Apply_disturbances_ar...'
           call apply_disturbances_ar(cgrid)
        end if
        
     end if

!     write (unit=*,fmt='(a)') '~~~ Update_C_and_N_pools_ar...'
     call update_C_and_N_pools_ar(cgrid)

!     write (unit=*,fmt='(a)') '~~~ Zero_ed_daily_vars...'
     call zero_ed_daily_vars(cgrid)

     ! Fuse patches last, after all updates have been applied.  This reduces
     ! the number of patch variables that actually need to be fused.  
     if(new_year) then
!        write (unit=*,fmt='(a)') '### Fuse_patchesar...'
        call fuse_patches_ar(cgrid)
     end if

     ! Recalculate the agb and basal area at the polygon level
!     write (unit=*,fmt='(a)') '~~~ Update_polygon_derived_props_ar...'
     call update_polygon_derived_props_ar(cgrid)

!     write (unit=*,fmt='(a)') '~~~ Print_C_and_N_budgets...'
     call print_C_and_N_budgets(cgrid)
  end do
  return
end subroutine vegetation_dynamics
!==========================================================================================!
!==========================================================================================!
