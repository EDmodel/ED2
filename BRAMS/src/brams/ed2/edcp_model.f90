subroutine ed_timestep(timel,dtlong)
  
  use mem_grid,only:ngrids,jdim
  use mem_leaf,only:isfcl
  use ed_node_coms,only:nnodetot,mynum
  use misc_coms,only:dtlsm,current_time
  use mem_edcp,only:edtime1,edtime2

  implicit none
  include 'mpif.h'
  integer :: ifm
  real(kind=8) :: timel
  real :: dtlong
  real :: wtime_start,t1,wtime1,wtime2,t2,wtime_tot
  real, external :: walltime
  integer :: ierr
  logical,save :: first = .true.


  ! Activate the ED2 model if time has reached the LSM timestep

  if (isfcl.ne.5)return
    
  if ( mod(timel+dble(dtlong),dble(dtlsm)) < dble(dtlong) .or. first ) then

     wtime_start=walltime(0.)
     !   CPU timing information & model timing information
     !   ===================================================

     call timing(1,t1)
     wtime1=walltime(wtime_start)

     ! Transfer the fluxes from the terrestrial and ocean models
     ! to the previous time arrays
     do ifm = 1,ngrids
        
        call copy_fluxes_future_2_past(ifm)
        
        call newgrid(ifm)
     
        call copy_atm2lsm(ifm,.false.)
     enddo
     
     call ed_coup_model()

     edtime1  = timel
     edtime2  = timel+dble(dtlsm)
     
     
     do ifm = 1,ngrids
        
        call newgrid(ifm)
     
        call simple_lake_model(ifm,dtlong)
     
        call copy_fluxes_lsm2atm(ifm)
     enddo
     
     wtime2=walltime(wtime_start)
     call TIMING(2,T2)
     
     if(mynum.eq.1) then
        write(*,"(a,a,i2.2,a,i2.2,a,i4.4,a,f6.0,2(a,f7.3),a)") &
             ' ED2 LSM Timestep ',&
             '; Sim time  ',current_time%month, &
             '-',current_time%date,           &
             '-',current_time%year,           &
             ' ',current_time%time,           &         
             's; Wall',wtime2-wtime1,&
             's; CPU',t2-t1,&
             's'
     endif

     
     ! If this is the first time the routine is called, then
     ! the ed_fluxp_g arrays (the flux arrays from previous time)
     ! are zero.  So just this once, copy the future arrays
     ! after they have been populated with real values, into
     ! the past arrays.  the transfer function needs both
     ! in order to do a temporal interpolation
     ! ---------------------------------------------------------

     if(first)then
        do ifm=1,ngrids
           call copy_fluxes_future_2_past(ifm)
        enddo
     endif


     first=.false.
  endif
  

  ! Update the leaf_g surface flux arrays of tstar,rstar and ustar
  ! This gets called every step because it is a time interpolation
  ! --------------------------------------------------------------

  do ifm = 1,ngrids
     call newgrid(ifm)
     call transfer_ed2leaf(ifm,timel)
  enddo

  ! Call a barrier
!  call MPI_Barrier(MPI_COMM_WORLD,ierr)



  return
end subroutine ed_timestep

!==================================================================

subroutine ed_coup_model()
  
  use misc_coms, only: integration_scheme, current_time, frqfast, frqstate     &
                      , out_time_fast, dtlsm, ifoutput, isoutput, idoutput     &
                      , imoutput, iyoutput, frqsum,unitfast,unitstate, imontha &
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
  
  integer :: doy,ierr

  integer,external :: julday
  
  
  character(len=10) :: c0, c1, c2
  character(len=*), parameter :: h="**(model)**"
  character(len=512) :: header_ext
  real,parameter :: frqcflx = 3600
  real           :: timefac_sst
  real           :: dtwater
  real, external :: update_sst_factor
  real :: ccont,ctemp,stemp,swat,lai

  real :: tfact1
  integer :: ipa,ico
  logical :: analysis_time, new_day, new_month, new_year, the_end
  logical :: writing_dail,writing_mont,writing_year,history_time,annual_time
  logical :: mont_analy_time,dail_analy_time,reset_time
  logical :: printbanner
  integer :: ndays
  integer, external :: num_days

  logical,save :: past_one_day   = .false.
  logical,save :: past_one_month = .false.
  
  istp = istp + 1
  
  writing_dail      = idoutput > 0
  writing_mont      = imoutput > 0
  writing_year      = iyoutput > 0
  out_time_fast     = current_time
  out_time_fast%month = -1


  !         Start the timesteps

  do ifm=1,ngrids
     call radiate_driver_ar(edgrid_g(ifm))
  end do
  
  if(integration_scheme == 0)then
     do ifm=1,ngrids
        call euler_timestep_ar(edgrid_g(ifm))
     end do
  elseif(integration_scheme == 1)then
     do ifm=1,ngrids
        call rk4_timestep_ar(edgrid_g(ifm),integration_buff_g)
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
  
  ! Checking whether it is some special time...
  new_day         = current_time%time < dtlsm
  if (.not. past_one_day .and. new_day) past_one_day=.true.
  
  new_month       = current_time%date == 1  .and. new_day
  if (.not. past_one_month .and. new_month) past_one_month=.true.
  
  new_year        = current_time%month == 1 .and. new_month
  mont_analy_time = new_month .and. writing_mont
  annual_time     = new_month .and. writing_year .and. current_time%month == outputMonth
  dail_analy_time = new_day   .and. writing_dail
  reset_time      = mod(time,dble(frqsum)) < dtlsm
  the_end         = mod(time,timmax) < dtlsm

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
!        call read_met_drivers_array()
        
        ! Re-allocate integration buffer
        if(integration_scheme == 1) call initialize_rk4patches_ar(0)
     endif
     
  endif
  
!  do ifm=1,ngrids
!     call update_met_drivers_array(edgrid_g(ifm))
!  end do

  if(new_day .and. new_month)then
     ! Loop all grids
     do ifm = 1,ngrids
        call updateHydroParms(edgrid_g(ifm))
     end do
  endif
  
  
  if(analysis_time)then
     if(new_month .and. new_day)then
        if(current_time%month == 6)then
           !              call update_ed_yearly_vars(polygon_list_g(1)%first_polygon)
           !              call zero_ed_yearly_vars(polygon_list_g(1)%first_polygon)
        endif
     endif
  endif

  !!Update Lateral Hydrology
  call calcHydroSubsurface()
  call calcHydroSurface()
  call writeHydro()
  

  return
end subroutine ed_coup_model
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

subroutine simple_lake_model(ifm,dtlong)

  use misc_coms,only:dtlsm

  use node_mod,only:ja,jz,ia,iz

  use consts_coms,only:stefan,cpi,vonk,cp,grav,p00,rocp,cpor
  
  use mem_edcp,only:wgrids_g,wgridf_g

  !------- Transfer these arrays to the polygons ----!
  use mem_leaf,   only: leaf_g
  use mem_basic,  only: basic_g
  use mem_radiate,only: radiate_g
  use mem_cuparm, only: cuparm_g
  use mem_micro,  only: micro_g
  use mem_grid,   only: zt,grid_g,dzt,zm,if_adap,jdim
  use therm_lib,  only: rslif
  !--------------------------------------------------!

  implicit none

  real, intent(in) :: dtlong
  real :: dtwb
  real :: cosz
  real :: exner
  real :: prss
  real :: ustar,tstar,rstar,thetacan,water_rsat,zts,water_rough
  real :: vels_pat,dtllohcc,dtllowcc
  real :: b,csm,csh,d,a2,c1,ri,fm,fh,c2,cm,ch,c3
  real :: ustaro,d_vel,d_veln,vel_new,delz
  integer :: ifix,ifixu
  real :: zoverl,bot,top,last
  real :: dtll_factor,dtll,hcapcan,wcapcan
  real :: z0fac_water,pis,rdi,idt
  integer :: niter_leaf,niter_can,nsubsteps,n

  real,parameter :: dtwbmax = 90.0  ! Water body timesteps can go no longer than 90 seconds
  real,external :: vertical_vel_flux
  real, parameter :: ubmin = .25   &    ! should use ubmin=1.0 for convec case
       ,ustmin = .1     !                 ubmin=0.1 for stable case

  real,parameter :: emiss_w = 0.97    ! emissivity of water (super gross approximation!)

  integer :: m1,m2,m3,i,j,m1max,ifm
  integer :: k1w,k2w,k3w,k2u,k2u_1,k2v,k2v_1
  real :: up_mean,vp_mean,pi0_mean,dn0_mean
  real :: rv_mean,theta_mean
  real :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2
  real :: canopy_water_vapor
  real :: canopy_tempk
  real :: sflux_u,sflux_v,sflux_w,sflux_t,sflux_r

  ! Note: Water albedo from Atwater and Bell (1981), excerpted from the LEAF3 scheme

  ! Using dtlong unless it is a very coarse timestep
  dtwb=min(dtlong,dtwbmax)
  
  ! The target step size for the water body suface flux is about 30 seconds
  
  nsubsteps = ceiling(dtlsm/dtwb)

  hcapcan = 2.0e4
  wcapcan = 2.0e1
  
  z0fac_water = .016 / grav

  ! Transfer atmospheric information to the ed lsm and fill those
  ! values.  Requires some mean calculations.  Also perform a 
  ! sanity check on the pressure; exit if it is unacceptable.
  !------------------------------------------------------------

  ! Use sib2 routine raddrv to calculate diffuse shortwave radiation
  !-----------------------------------------------------------------

  ! Prepare the atm data fields
  
  
  do j=ja,jz
     do i=ia,iz
        if (if_adap == 1) then
           
           ! Shaved Eta coordinate system
           !--------------------------------------
           
           k2w = nint(grid_g(ifm)%flpw(i,j))
           k1w = k2w - 1
           k3w = k2w + 1
           
           k2u   = nint(grid_g(ifm)%flpu(i,j))
           k2u_1 = nint(grid_g(ifm)%flpu(i-1,j))
           
           k2v   = nint(grid_g(ifm)%flpv(i,j))
           k2v_1 = nint(grid_g(ifm)%flpv(i,j-jdim))
           
           topma_t = .25 * (grid_g(ifm)%topma(i,j) + grid_g(ifm)%topma(i-1,j)  &
                + grid_g(ifm)%topma(i,j-jdim) + grid_g(ifm)%topma(i-1,j-jdim))
           
           ! weights for lowest predicted points, relative to points above them
           
           wtw = (zm(k2w) - topma_t) * dzt(k2w)
           wtu1 = grid_g(ifm)%aru(k2u_1,i-1,j)   / grid_g(ifm)%aru(k2u_1+1,i-1,j)
           wtu2 = grid_g(ifm)%aru(k2u,i,j)       / grid_g(ifm)%aru(k2u+1,i,j)
           wtv1 = grid_g(ifm)%arv(k2v_1,i,j-jdim) / grid_g(ifm)%arv(k2v_1+1,i,j-jdim)
           wtv2 = grid_g(ifm)%arv(k2v,i,j)       / grid_g(ifm)%arv(k2v+1,i,j)
           
           theta_mean   =  wtw * basic_g(ifm)%theta(k2w,i,j) + (1. - wtw)  * basic_g(ifm)%theta(k3w,i,j)
           
           rv_mean      =  wtw * basic_g(ifm)%rv(k2w,i,j)    + (1. - wtw)  * basic_g(ifm)%rv(k3w,i,j)
           
           up_mean      = (wtu1        * basic_g(ifm)%up(k2u_1,i-1,j)    &
                +  (1. - wtu1) * basic_g(ifm)%up(k2u_1+1,i-1,j)  &
                +  wtu2        * basic_g(ifm)%up(k2u,i,j)        &
                +  (1. - wtu2) * basic_g(ifm)%up(k2u+1,i,j)) * .5
        
           vp_mean      = (wtv1        * basic_g(ifm)%vp(k2v_1,i,j-jdim)    &
                +  (1. - wtv1) * basic_g(ifm)%vp(k2v_1+1,i,j-jdim)  &
                +  wtv2        * basic_g(ifm)%vp(k2v,i,j)          &
                +  (1. - wtv2) * basic_g(ifm)%vp(k2v+1,i,j)) * .5
           
           if (wtw >= .5) then
              pi0_mean   = ((wtw - .5) * (basic_g(ifm)%pp(k1w,i,j) + basic_g(ifm)%pi0(k1w,i,j))  &
                   + (1.5 - wtw) * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j)))
              dn0_mean   = (wtw - .5)  * basic_g(ifm)%dn0(k1w,i,j)  &
                   + (1.5 - wtw) * basic_g(ifm)%dn0(k2w,i,j)
           else
              pi0_mean  = ((wtw + .5) * (basic_g(ifm)%pp(k2w,i,j) + basic_g(ifm)%pi0(k2w,i,j))  &
                   + (.5 - wtw) * (basic_g(ifm)%pp(k3w,i,j) + basic_g(ifm)%pi0(k3w,i,j)))
              dn0_mean  = (wtw + .5) * basic_g(ifm)%dn0(k2w,i,j)  &
                   + (.5 - wtw) * basic_g(ifm)%dn0(k3w,i,j)
           endif
           
           
        else
           
           ! Terrain following coordinate system
           !--------------------------------------       
           
           theta_mean   = basic_g(ifm)%theta(2,i,j)
           rv_mean      = basic_g(ifm)%rv(2,i,j)
           
           up_mean       = (basic_g(ifm)%up(2,i,j) + basic_g(ifm)%up(2,i-1,j)) * 0.5
           vp_mean       = (basic_g(ifm)%vp(2,i,j) + basic_g(ifm)% vp(2,i,j-jdim)) * 0.5
           pi0_mean      = (basic_g(ifm)%pp(1,i,j) + basic_g(ifm)%pp(2,i,j) & 
                +basic_g(ifm)%pi0(1,i,j)+basic_g(ifm)%pi0(2,i,j)) * 0.5
           
           if(pi0_mean.ne.pi0_mean)then
              print*,'bad pi0 mean'
              print*,i,j,basic_g(ifm)%pp(1,i,j),basic_g(ifm)%pp(2,i,j), &
                   basic_g(ifm)%pi0(1,i,j),basic_g(ifm)%pi0(2,i,j)
              stop
           endif
           
           dn0_mean     = (basic_g(ifm)%dn0(1,i,j) + basic_g(ifm)%dn0(2,i,j)) * 0.5
           
        endif

        pi0_mean = pi0_mean  ! Convert to pascals

        cosz  = radiate_g(ifm)%cosz(i,j)
        exner = cp * (pi0_mean*100.0 / p00)**rocp

        ustar              = wgridf_g(ifm)%ustar(i,j)
        canopy_tempk       = wgrids_g(ifm)%canopy_tempk(i,j)
        canopy_water_vapor = wgrids_g(ifm)%canopy_water_vapor(i,j)

        sflux_u = 0.0
        sflux_v = 0.0
        sflux_w = 0.0
        sflux_r = 0.0
        sflux_t = 0.0

        idt = 0
        do n = 1,nsubsteps
           
           dtllohcc = min(dtlsm-idt,dtwb)/hcapcan
           dtllowcc = min(dtlsm-idt,dtwb)/wcapcan

           idt = idt+dtwb

           pis =  pi0_mean * cpi
           
           prss = pis ** cpor * p00

           water_rsat  = rslif(prss, leaf_g(ifm)%seatp(i,j))
           
           water_rough = max(z0fac_water * ustar ** 2,.0001)
           
           thetacan = canopy_tempk / pis
           
           if(thetacan.lt.theta_mean) then
              vels_pat = max(0.1,sqrt(up_mean**2 + vp_mean**2))
           else
              vels_pat = max(1.0,sqrt(up_mean**2 + vp_mean**2))
           endif
           
           zts = zt(2) + grid_g(ifm)%rtgt(i,j)
           b = 5.
           csm = 7.5
           csh = 5.
           d = 5.
           
           ! a2 is the drag coefficient in neutral conditions, here same for h/m
           ! ri is the bulk richardson numer, eq. 3.45 in Garratt
           
           a2 = (vonk / log(zts / water_rough)) ** 2
           c1 = a2 * vels_pat
           ri = grav * zts * (theta_mean - thetacan)  &
                / (.5 * (theta_mean + thetacan) * vels_pat * vels_pat)
           
           if (theta_mean - thetacan > 0.) then   ! STABLE CASE
              
              fm = 1. / (1. + (2 * b * ri / sqrt(1 + d * ri)))
              fh = 1. / (1. + (3 * b * ri * sqrt(1 + d * ri)))
              
           else                            ! UNSTABLE CASE
              
              c2 = b * a2 * sqrt(zts / water_rough * (abs(ri)))
              cm = csm * c2
              ch = csh * c2
              fm = (1. - 2 * b * ri / (1. + 2 * cm))
              fh = (1. - 3 * b * ri / (1. + 3 * ch))
              
           endif
           
           ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
           c3 = c1 * fh / ustar
           rstar = c3 * (rv_mean - canopy_water_vapor)
           tstar = c3 * (theta_mean - thetacan)
           
           ! Calculate the water surface fluxes using the friction velocities
           ! ---------------------------------------------------------------
           ! Limit ustar so that the flux cannot take more than 1/2 velocity in a timestep
           
           ifixu=0
           ustaro=ustar
           delz = 2.*zts
           d_vel =  - ustar * ustar *dtlsm / delz
           vel_new = vels_pat + d_vel
           if(vel_new < .5 * vels_pat) then
              ifixu=1
              d_veln = .5 * vels_pat
              ustar=sqrt(d_veln*delz/dtlsm)
           endif

           ! Calculate the heat,moisture and momentum fluxes
           ! -----------------------------------------------

           sflux_u = sflux_u - ustar*ustar*up_mean/vels_pat
           sflux_v = sflux_v - ustar*ustar*vp_mean/vels_pat
           sflux_t = sflux_t - ustar*tstar
           sflux_r = sflux_r - ustar*rstar
           sflux_w = sflux_w + vertical_vel_flux(grav * zts * cpi * exner / (theta_mean*pis)    &
                ,tstar,ustar)
           
           ! Update the sea surface air temperature and water vapor mixing ratio
           ! -------------------------------------------------------------------
           
           rdi = .2 * ustar
           
           canopy_tempk  = canopy_tempk        &
                + dtllohcc * dn0_mean * cp                          &
                * ( (leaf_g(ifm)%seatp(i,j) -  canopy_tempk) * rdi    &
                + ustar * tstar * pis)
           
           bot = ( water_rsat - canopy_water_vapor) * rdi
           top = ustar * rstar
           last = canopy_water_vapor
           
           canopy_water_vapor = canopy_water_vapor &
                + dtllowcc * dn0_mean * (( water_rsat-canopy_water_vapor) * rdi  &
                + ustar * rstar )
           
           if(canopy_water_vapor .lt. 0.001 ) then
              print*,"WATER VAPOR IN CAS ABOVE WATER-BODIES IS SCREWY"
              print*,i,j,n
              print*,pi0_mean, leaf_g(ifm)%seatp(i,j)
              print*,canopy_water_vapor,water_rsat,rdi,dtllowcc,dn0_mean,ustar,rstar
              print*,rv_mean,c3,c1,fh,bot,top,last
              print*,grid_g(ifm)%glat,grid_g(ifm)%glon
              stop
           endif
     
        enddo

        ! Transfer model scalars back to global arrays
        ! --------------------------------------------

        wgridf_g(ifm)%ustar(i,j) = ustar
        wgridf_g(ifm)%rstar(i,j) = rstar
        wgridf_g(ifm)%tstar(i,j) = tstar

        wgridf_g(ifm)%sflux_u(i,j) = dn0_mean*sflux_u/real(nsubsteps)
        wgridf_g(ifm)%sflux_v(i,j) = dn0_mean*sflux_v/real(nsubsteps)
        wgridf_g(ifm)%sflux_w(i,j) = dn0_mean*sflux_w/real(nsubsteps)
        wgridf_g(ifm)%sflux_t(i,j) = dn0_mean*sflux_t/real(nsubsteps)
        wgridf_g(ifm)%sflux_r(i,j) = dn0_mean*sflux_r/real(nsubsteps)

        wgridf_g(ifm)%albedt(i,j) = min(max(-.0139 + .0467 * tan(acos(cosz)),.03),.999)
        wgridf_g(ifm)%rlongup(i,j) = emiss_w * stefan * leaf_g(ifm)%seatp(i,j)**4
        wgrids_g(ifm)%canopy_tempk(i,j) = canopy_tempk
        wgrids_g(ifm)%canopy_water_vapor(i,j) = canopy_water_vapor
        

     enddo
  enddo

!  print*,"LAKE MODEL",ja,ja,ia,iz
!  print*,"ustar:", wgridf_g(ifm)%ustar(2,2)," tstar: " &
!       , wgridf_g(ifm)%tstar(2,2)," rstar:", wgridf_g(ifm)%rstar(2,2)
!  print*,"CANOPY TEMP:",wgrids_g(ifm)%canopy_tempk(2,2)
!  print*,"SEA TEMP:",leaf_g(ifm)%seatp(2,2)
!  print*,"WATER VAPOR ",wgrids_g(ifm)%canopy_water_vapor(ia:iz,ja:jz)
!  print*,"U MOMENTUM FLUX",wgridf_g(ifm)%sflux_u(2,2)
!  print*,"HEAT FLUX:",wgridf_g(ifm)%sflux_t(2,2)
!  print*,"MOISTURE FLUX:",wgridf_g(ifm)%sflux_r(2,2)

  return
end subroutine simple_lake_model
