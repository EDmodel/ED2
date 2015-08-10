!==============================================================================
! Driver and core for numerically integrating the canopy via BDF2
! BDF2 is a trapezoidal split step method (backward and forward)
!
! To solve for the next step, we must solve a set of linear equations
! using data on the current and previous step.
!
! U_{n+1} = [4U_{n}+U_{n-1}][3I-2dtA]^{-1}
! 
! The crux of this method is solving the for the instantaneous derivatives
! dU/dt = A*U_{n+1} + B          (where B is a forcing like SW radiation
!                                 energy flux from the ABL, or energy flux
!                                 from the ground.  The mass and energy
!                                 flux from the ground is solved elsewhere
!
! One difficulty in creating the A matrix, is that some of these processes
! are nearly binary...is there ponded water or not at the next step.
! Therefore, it may be necessary to to a sanity check on the existence of
! surface water and test to see if it was correct at the end of the step.
!
!==============================================================================

subroutine bdf2_solver(cpatch,yprev,ycurr,ynext,dydt,dtf,dtb)

  use grid_coms,only         : nzg,nzs
  use rk4_coms,only          : effarea_evap,effarea_heat,effarea_transp, &
       rk4site,rk4patchtype,rk4aux,bdf2patchtype,checkbudget,ibranch_thermo
  use ed_misc_coms,only    : fast_diagnostics
  use ed_state_vars,only   : patchtype
  use therm_lib8,only      : reducedpress8,idealdenssh8
  use ed_therm_lib,only    : ed_grndvap8
  use consts_coms,only     : cpdry8,p00i8,rdry8, &
       rocv8,cpocv8,cliq8,cice8,rocp8,cpocv,epim1,   &
       alli8,t3ple8,alvi38,wdns8,cph2o8,tsupercool_liq8,pi18
  use therm_lib8,only      : uint2tl8,uextcm2tl8,tq2enthalpy8, tl2uint8, cmtl2uext8
  use soil_coms, only      : soil8,dslz8
  !$ use omp_lib
  implicit none
  
  ! define the previous,current and next patch states
  integer                                    :: nstate
  real(kind=8),pointer,dimension(:)          :: Y
  real(kind=8),pointer,dimension(:)          :: Yf
  real(kind=8),pointer,dimension(:,:)        :: A
  real(kind=8),pointer,dimension(:)          :: B
  type(patchtype),target                     :: cpatch
  type(bdf2patchtype), target                :: yprev
  type(rk4patchtype), target                 :: ycurr
  type(rk4patchtype), target                 :: ynext
  type(rk4patchtype), target                 :: dydt
  real(kind=8),intent(in)                    :: dtf        ! dt for n -> n+1
  real(kind=8),intent(in)                    :: dtb        ! dt for n -> n-1

  integer :: id_tcan,id_tveg,ico
  integer :: k,ksn
  integer :: i,j
  real(kind=8)                  :: shctop
  real(kind=8)                  :: rhoc  ! canopy air density [kg/m3]
  real(kind=8)                  :: dc    ! canopy depth       [m]
  real(kind=8)                  :: qc    ! canopy spec. hum   [kg/kg]
  real(kind=8)                  :: Tg    ! soil temperature   [K]
  real(kind=8)                  :: gg    ! ground-can conductivity [m/s]
  real(kind=8)                  :: mgc   ! mass flux ground->can   [kg/s/m2]
  real(kind=8)                  :: ga    ! atm-canopy conductivity [m/s]
  real(kind=8)                  :: Ta    ! atm temperature         [K]
  real(kind=8)                  :: mac   ! mass flux atm->can      [kg/s/m2]
  real(kind=8)                  :: dqcdt ! time partial of qc      [kg/kg/s]
  real(kind=8)                  :: xc    ! lumped term (canopy)    [-]
  real(kind=8)                  :: href  ! reference spec. enth.   [J/kg]
  real(kind=8)                  :: gv    ! veg-canopy conductivity [m/s]
  real(kind=8),dimension(300)   :: gv_hfa   ! effective area of heat times g
  real(kind=8),dimension(300)   :: mlc   ! mass flux leaf->can      [kg/m2/s]
  real(kind=8),dimension(300)   :: mwc   ! mass flux wood->can
  real(kind=8),dimension(300)   :: mtr   ! mass flux transp veg->can [kg/m2/s]
  real(kind=8),dimension(300)   :: fliq  ! liquid fraction leaf surf [kg/kg]
  real(kind=8),dimension(300)   :: qv    ! water mass on vegetation  [kg/m2]
  real(kind=8),dimension(300)   :: dqvdt ! time partial of qv  [kg/m2/s]
  real(kind=8),dimension(300)   :: xv    ! lumped term (vegetation)
  real(kind=8),dimension(300)   :: hflx  ! enth flux from bunch of stuff
  logical,dimension(300)        :: resolve ! self-explanatory

  real(kind=8)  :: veg_temp,leaf_temp,wood_temp
  real(kind=8)  :: eflxac
  real(kind=8)  :: qwflxac
  real(kind=8)  :: hflxac
  real(kind=8)  :: qwflxlc,qwflxlc_tot
  real(kind=8)  :: qwflxwc,qwflxwc_tot
  real(kind=8)  :: qtransp,qtransp_tot
  real(kind=8)  :: hflxlc,hflxlc_tot
  real(kind=8)  :: hflxwc,hflxwc_tot
  integer       :: ibuff

  ibuff = 1
  !$ ibuff = OMP_get_thread_num()+1

  !==========================================================================!
  !                                                                          !
  ! Part 1: Build the forecasted transition matrix A                         !
  !                                                                          !
  ! Assume that Y takes the following form                                   !
  ! Y = [ T_can T_veg1 T_veg2 ... T_vegN ]                                   !
  !                                                                          !
  !==========================================================================!

  !==========================================================================!
  ! Row and Column indices for the transition matrix operators               !
  !==========================================================================!
  id_tcan   = 1

  !==========================================================================!
  ! Heat flux to the atmosphere... few ways to go about this                 !
  ! 1) use the heat flux at step time n                                      !
  ! 2) use atmospheric temp at step n, but incorporate the canopy temp at n+1!
  !    to calculate flux, uses conductivity at n                             !
  ! 3) use forward time step of atmospheric temp as well (a little overboard)!
  ! 4) also would use the forward euler canopy and leaf temps to recalc      !
  !    aerodynamic conductances                                              !
  ! **** Using method 2 right now                                            !
  !==========================================================================!

  !--------------------------------------------------------------------------!
  ! Ground Temp is estimated via an explicit method, it represents           !
  ! the effective ground surface temperature, be it standing water           !
  ! or moist/dry soil                                                        !
  !--------------------------------------------------------------------------!
  
!!  ksn=0
!!  do k=1,ycurr%nlev_sfcwater
!!     if(ycurr%sfcwater_mass(k)>1.d-5)ksn=k
!!  end do
  
  shctop = soil8(rk4site%ntext_soil(nzg))%slcpd
  call uextcm2tl8(ynext%soil_energy(nzg),ynext%soil_water(nzg)*wdns8,shctop &
       ,ynext%soil_tempk(nzg),ynext%soil_fracliq(nzg))

  ksn=ycurr%nlev_sfcwater
  if (ksn>0) then
     call uint2tl8(ycurr%sfcwater_energy(ksn)/ycurr%sfcwater_mass(ksn),              &
          ycurr%sfcwater_tempk(ksn),ycurr%sfcwater_fracliq(ksn))
  endif


  k=max(1,ksn)

  call ed_grndvap8(ksn,ynext%soil_water(nzg),ynext%soil_tempk(nzg)        &
       ,ynext%soil_fracliq(nzg),ycurr%sfcwater_tempk(k)                   &
       ,ycurr%sfcwater_fracliq(k),ycurr%snowfac,ycurr%can_prss,ynext%can_shv &
       ,ynext%ground_shv,ynext%ground_ssh,ynext%ground_temp               &
       ,ynext%ground_fliq,ynext%ggsoil)
  

  !===========================================================================!
  ! Determine the leaf-branch thermodynamics scheme to be used.
  ! Set some temporary variables.
  ! Calculate the number of elements.
  !===========================================================================!


  resolve(:)=.false.
  id_tveg = 1
  select case(ibranch_thermo)
  case(0,2)
     do ico=1,cpatch%ncohorts
        if (ycurr%leaf_resolvable(ico)) then
           gv_hfa(ico) = ycurr%leaf_gbh(ico)/(ycurr%can_rhos*cpdry8) &
                *effarea_heat*ycurr%lai(ico)
           mlc(ico)  = ycurr%wflxlc(ico)
           mtr(ico)  = ycurr%wflxtr(ico)
           fliq(ico) = ycurr%leaf_fliq(ico)
           qv(ico)   = ynext%leaf_water(ico)
           dqvdt(ico)= dydt%leaf_water(ico)
           xv(ico)   = ycurr%leaf_fliq(ico)*cliq8*qv(ico) + &
                       (1.d0-fliq(ico))*cice8*qv(ico)     + &
                       ycurr%leaf_hcap(ico)
           hflx(ico) = ycurr%hflx_lrsti(ico)
           resolve(ico) = .true.
           id_tveg   = id_tveg + 1
        end if
     end do
     
     if(ibranch_thermo.eq.2)then
        do ico=1,cpatch%ncohorts
           if (ycurr%wood_resolvable(ico)) then
              id_tveg   = id_tveg + 1
           end if
        end do
     end if

  case(1)
     do ico=1,cpatch%ncohorts
        gv_hfa(ico) = 0.0
        mlc(ico)    = 0.0
        mtr(ico)    = 0.0
        fliq(ico)   = -1.0
        qv(ico)     = 0.0
        dqvdt(ico)  = 0.0
        hflx(ico)   = 0.0
        if (ycurr%leaf_resolvable(ico) .or. ycurr%wood_resolvable(ico))then
           id_tveg = id_tveg+1
           resolve(ico) = .true.
        end if

        if (ycurr%leaf_resolvable(ico)) then
           gv_hfa(ico) = ycurr%leaf_gbh(ico)/(ycurr%can_rhos*cpdry8) &
                *effarea_heat*ycurr%lai(ico)
           mlc(ico)  = ycurr%wflxlc(ico)
           mtr(ico)  = ycurr%wflxtr(ico)
           fliq(ico) = ycurr%leaf_fliq(ico)
           qv(ico)   = ynext%leaf_water(ico)
           dqvdt(ico)= dydt%leaf_water(ico)
           hflx(ico) = ycurr%hflx_lrsti(ico)
        end if
        if (ycurr%wood_resolvable(ico)) then
           gv_hfa(ico) = gv_hfa(ico)+(ycurr%wood_gbh(ico)/(ycurr%can_rhos*cpdry8) &
                *effarea_heat*ycurr%wai(ico))
           mlc(ico)  = mlc(ico)+ycurr%wflxwc(ico)
           if(fliq(ico)>-0.99)then
              fliq(ico) = (fliq(ico)+ycurr%wood_fliq(ico))/2.0
           else
              fliq(ico) = ycurr%wood_fliq(ico)
           end if
           
           qv(ico)   = qv(ico)+ynext%wood_water(ico)
           dqvdt(ico)= dqvdt(ico)+dydt%wood_water(ico)
           hflx(ico) = hflx(ico)+ycurr%hflx_wrsti(ico)
        end if

        ! We use the wood and the leaf heat capacities regardless 
        ! of their ability to be solved.
        ! ---------------------------------------------------------

        xv(ico)   = ycurr%leaf_fliq(ico)*cliq8*ynext%leaf_water(ico)        + &
                    (1.d0-ycurr%leaf_fliq(ico))*cice8*ynext%leaf_water(ico) + &
                    ycurr%leaf_hcap(ico)                                    + &
                    ycurr%wood_fliq(ico)*cliq8*ynext%wood_water(ico)        + &
                    (1.d0-ycurr%wood_fliq(ico))*cice8*ynext%wood_water(ico) + &
                    ycurr%wood_hcap(ico)

     end do

  end select

  nstate = id_tveg

  !------ Allocate the matrices used for the linear operations --------------!
  
  allocate(Y(nstate))
  allocate(Yf(nstate))
  allocate(A(nstate,nstate))
  allocate(B(nstate))

  !----- Initialize the matrices used for the linear operations -------------!

  A = 0.d0
  B = 0.d0
  Y = 0.d0
  Yf= 0.d0


  !----- Set temporary canopy level variables used for readability ----------!
  
  rhoc  = ycurr%can_rhos
  dc    = ycurr%can_depth
  qc    = ynext%can_shv
  Tg    = ycurr%ground_temp
  gg    = ycurr%ggnet
  mgc   = ycurr%wflxsc + ycurr%wflxgc
  ga    = ycurr%ggbare
  Ta    = rk4site%atm_theta*rk4site%atm_exner/cpdry8
  mac   = ycurr%wflxac
  dqcdt = dydt%can_shv
  xc    = rhoc*dc*((1.d0-qc)*cpdry8+qc*cph2o8)
  href  = t3ple8*cice8+alvi38-cph2o8*t3ple8

  ! USES NEW TG
!!  B(id_tcan) = (1.d0/xc)*       &
!!       (gg*rhoc*cpdry8*Tg       &
!!       + mgc*(href+cph2o8*Tg)    &
!!       + ga*rhoc*cpdry8*Ta      &
!!       + mac*href               &
!!       + mac*cph2o8*0.5d0*Ta     &
!!       + sum(ycurr%wflxlc)*href &
!!       + sum(ycurr%wflxtr)*href &
!!       + sum(ycurr%wflxwc)*href &
!!       - dc*rhoc*href*dqcdt)

  ! USES ycurr HFLXSC+HFLXGC
  B(id_tcan) = (1.d0/xc)*       &
       ( ycurr%hflxsc           & 
       + ycurr%hflxgc           &
       + ycurr%qwflxsc          &
       + ycurr%qwflxgc          &
  !    + mgc*(href+cph2o8*Tg)   &
       + ga*rhoc*cpdry8*Ta      &
       + mac*href               &
       + mac*cph2o8*0.5d0*Ta     &
       + sum(ycurr%wflxlc)*href &
       + sum(ycurr%wflxtr)*href &
       + sum(ycurr%wflxwc)*href &
       - dc*rhoc*href*dqcdt)



  ! USES NEW TG
!!  A(id_tcan,id_tcan) = (1.d0/xc)*           &
!!       (-gg*rhoc*cpdry8                     &
!!       -ga*rhoc*cpdry8                      &
!!       +0.5d0*mac*cph2o8)                   &
!!       -(dqcdt*dc*rhoc/(xc**2.d0))*           &
!!       ((1.d0-qc)*cpdry8 + qc*cph2o8)*(cph2o8-cpdry8)
  
  ! USES ycurr HFLXGC
  A(id_tcan,id_tcan) = (1.d0/xc)*           &
       (-ga*rhoc*cpdry8                      &
       +0.5d0*mac*cph2o8)                   &
       -(dqcdt*dc*rhoc/(xc**2.d0))*           &
       ((1.d0-qc)*cpdry8 + qc*cph2o8)*(cph2o8-cpdry8)
  



  Y(id_tcan) = (3.d0+(dtf/dtb))*ycurr%can_temp -        &
       (dtf/dtb)*yprev%can_temp + 2.d0*B(id_tcan)*dtf



  !===========================================================================!
  ! Derivatives associated with vegetation <-> canopy-air fluxes
  !===========================================================================!
  
  id_tveg = id_tcan

  do ico=1,cpatch%ncohorts

     if (resolve(ico))then
        
        id_tveg = id_tveg+1

        ! dTc/dt ~ Tc)
        A(id_tcan,id_tcan) = A(id_tcan,id_tcan) - gv_hfa(ico)*rhoc*cpdry8/xc
        
        ! A(dTc/dt ~ Tv)
        A(id_tcan,id_tveg) = A(id_tcan,id_tveg)  &
             + (1.d0/xc)*(gv_hfa(ico)*rhoc*cpdry8 + mlc(ico)*cph2o8 + mtr(ico)*cph2o8)
        
        ! B(dTv/dt)
        B(id_tveg) = (1.d0/xv(ico))*           &
             ( hflx(ico)                       &
             - href*mlc(ico)                   &
             - href*mtr(ico)                   &
             + (fliq(ico)*cliq8*t3ple8 - fliq(ico)*cice8*t3ple8 - fliq(ico)*alli8)*dqvdt(ico))
        

        ! A(dTv/dt ~ Tc)
        A(id_tveg,id_tcan) = A(id_tveg,id_tcan) + &
             (gv_hfa(ico)*rhoc*cpdry8)/xv(ico)


        ! A(dTv/dt ~ Tv)
        A(id_tveg,id_tveg) = A(id_tveg,id_tveg) + &
             (1.d0/xv(ico))*( (fliq(ico)*cice8 - fliq(ico)*cliq8 - cice8)*dqvdt(ico) &
             - gv_hfa(ico)*rhoc*cpdry8 - (mlc(ico)+mtr(ico))*cph2o8)


        Y(id_tveg) = (3.d0+dtf/dtb)*ycurr%leaf_temp(ico) - &
             (dtf/dtb)*yprev%leaf_temp(ico) + &
             2.d0*B(id_tveg)*dtf


     end if
        
  end do

  ! ============================================================!
  ! If this is a thermo case where wood is separate, do another
  ! ============================================================!
  if(ibranch_thermo.eq.2)then

     do ico=1,cpatch%ncohorts
        if (ycurr%wood_resolvable(ico)) then
           gv_hfa(ico) = ycurr%wood_gbh(ico)/(ycurr%can_rhos*cpdry8) &
                *effarea_heat*ycurr%wai(ico)
           mlc(ico)  = ycurr%wflxwc(ico)
           mtr(ico)  = 0.0
           fliq(ico) = ycurr%wood_fliq(ico)
           qv(ico)   = ynext%wood_water(ico)
           dqvdt(ico)= dydt%wood_water(ico)
           xv(ico)   = ycurr%wood_fliq(ico)*cliq8*qv(ico) + &
                (1.d0-fliq(ico))*cice8*qv(ico) +          &
                ycurr%wood_hcap(ico)
           hflx(ico) = ycurr%hflx_wrsti(ico)
        end if
     end do
     
     do ico=1,cpatch%ncohorts
        if (ycurr%wood_resolvable(ico)) then

           id_tveg = id_tveg+1

           ! dTc/dt ~ Tc)
           A(id_tcan,id_tcan) = A(id_tcan,id_tcan) - gv_hfa(ico)*rhoc*cpdry8/xc
        
           ! A(dTc/dt ~ Tv)
           A(id_tcan,id_tveg) = A(id_tcan,id_tveg)  &
                + (1.d0/xc)*(gv_hfa(ico)*rhoc*cpdry8 + mlc(ico)*cph2o8 + mtr(ico)*cph2o8)
           
           ! B(dTv/dt)
           B(id_tveg) = (1.d0/xv(ico))*           &
                ( hflx(ico)                       &
                - href*mlc(ico)                   &
                - href*mtr(ico)                   &
                + (fliq(ico)*cliq8*t3ple8 - fliq(ico)*cice8*t3ple8 &
                - fliq(ico)*alli8)*dqvdt(ico))
        
           
           ! A(dTv/dt ~ Tc)
           A(id_tveg,id_tcan) = A(id_tveg,id_tcan) + &
                (gv_hfa(ico)*rhoc*cpdry8)/xv(ico)
           
           
           ! A(dTv/dt ~ Tv)
           A(id_tveg,id_tveg) = A(id_tveg,id_tveg) + &
                (1.d0/xv(ico))*( (fliq(ico)*cice8 - fliq(ico)*cliq8 - cice8)*dqvdt(ico) &
                - gv_hfa(ico)*rhoc*cpdry8 - (mlc(ico)+mtr(ico))*cph2o8)

           Y(id_tveg) = (3.d0+dtf/dtb)*ycurr%wood_temp(ico) - &
                (dtf/dtb)*yprev%wood_temp(ico) + &
                2.d0*B(id_tveg)*dtf
           
        end if
     end do
  end if


  ! Create the matrix  [3I-2Adt]
  ! And create the upper portion of the equation [4U_{n}+U_{n-1}+2Bdt]

  do i=1,nstate
     do j=1,nstate
        A(i,j)=-A(i,j)*2.d0*dtf
     end do
     A(i,i) = 3.d0+A(i,i)
  end do
  

  call selective_gaussian_2body(A,Yf,Y,nstate)

  Y=Yf

  ! ------------------------------------------------------------------------!
  ! Set the leaf and canopy temepratures in the memory buffer (rk4type)     !
  ! Update vars that otherwise would not have been diagnostic if we were    !
  ! using a different scheme.  For instance, rk4 and forward euler use      !
  ! can_enthalpy
  !-------------------------------------------------------------------------!
  
  ! Send the canopy and leaf temperatures to the forward 
  ! Be sure to update leaf_energy and canopy ln-theta


  ynext%can_temp = Y(1)
  ynext%can_enthalpy = (1.d0-qc)*cpdry8*Y(1) + qc*(href + cph2o8*Y(1))

  ! ------------------------------------------------------------------------!
  ! Note: Significant assumption being made here.  The partial derivative   !
  ! of the leaf and wood temperature assumed that phase was constant.       !
  ! This is true unless of course the water starts oscillating around the   !
  ! freezing point.  This solver will assume the veg system will heat and   !
  ! cool without phase change, which of course is not true, but nothing is  !
  ! perfect, so deal with it.  AFter ynext%leaf_Temp crosses t3ple8, it will!
  ! change the liquid fraction to either 0 or 1 depending on the direction  !
  ! of change.                                                              !
  ! One possible alternative, is to determine the change in energy assuming !
  ! no phase change, and then apply that energy to the phase change at the  !
  ! triple point.  If the change in energy is less than phase change, set   !
  ! the ynext temp to t3ple and linearly scale the liquid fraction.         !
  ! ------------------------------------------------------------------------!

  qwflxlc_tot = 0.d0
  qwflxwc_tot = 0.d0
  qtransp_tot = 0.d0
  hflxlc_tot  = 0.d0
  hflxwc_tot  = 0.d0
  

  id_tveg=1
  
  select case(ibranch_thermo)
  case(0,2)
     do ico=1,cpatch%ncohorts
        if (resolve(ico)) then
           id_tveg=id_tveg+1
           ynext%leaf_temp(ico) = Y(id_tveg)
           
           if(ynext%leaf_temp(ico) < t3ple8) then
              ynext%leaf_fliq(ico) = 0.d0
              ynext%leaf_energy(ico) = &
                   ynext%leaf_water(ico)*cice8*ynext%leaf_temp(ico) +&
                   ynext%leaf_temp(ico)*ycurr%leaf_hcap(ico)
           else
              ynext%leaf_fliq(ico) = 1.d0
              ynext%leaf_energy(ico) = ynext%leaf_temp(ico)*              &
                   (ycurr%leaf_hcap(ico)+ynext%leaf_water(ico)*cliq8) - &
                   ynext%leaf_water(ico)*cliq8*tsupercool_liq8
           end if
           
           ! Back calculate the latent and sensible heat fluxes of leaves
           ! ========================================================================
           
           ! First calculate the effective qflxlc
           qwflxlc = ycurr%wflxlc(ico)*tq2enthalpy8(ycurr%leaf_temp(ico),1.d0,.true.)        
           
           ! Then effective transpiraiton
           qtransp = ycurr%wflxtr(ico)*tq2enthalpy8(ycurr%leaf_temp(ico),1.d0,.true.)
           
           ! Use the resulting change in leaf energy to back-caculate what heat
           ! flux would had been
           hflxlc  = ycurr%hflx_lrsti(ico) - qwflxlc - qtransp - &
                (ynext%leaf_energy(ico)-ycurr%leaf_energy(ico))/dtf
           
           qwflxlc_tot = qwflxlc_tot + qwflxlc
           qtransp_tot = qtransp_tot + qtransp
           hflxlc_tot  = hflxlc_tot  + hflxlc
           
        end if
     end do
     
     if(ibranch_thermo.eq.2) then
        do ico=1,cpatch%ncohorts
           if (ycurr%wood_resolvable(ico)) then
              id_tveg=id_tveg+1
              ynext%wood_temp(ico) = Y(id_tveg)
              if(ynext%wood_temp(ico) < t3ple8) then
                 ynext%wood_fliq(ico) = 0.d0
                 ynext%wood_energy(ico) = &
                      ynext%wood_water(ico)*cice8*ynext%wood_temp(ico) +&
                      ynext%wood_temp(ico)*ycurr%wood_hcap(ico)
              else
                 ynext%wood_fliq(ico) = 1.d0
                 
                 ynext%wood_energy(ico) = ynext%wood_temp(ico)*              &
                      (ycurr%wood_hcap(ico)+ynext%wood_water(ico)*cliq8) - &
                      ynext%wood_water(ico)*cliq8*tsupercool_liq8
              end if
              
              ! Back calculate the latent and sensible heat fluxes of wood
              ! ========================================================================
              qwflxwc = ycurr%wflxwc(ico)*tq2enthalpy8(ycurr%wood_temp(ico),1.d0,.true.)
              
              hflxwc  = ycurr%hflx_wrsti(ico) - qwflxwc - &
                   (ynext%wood_energy(ico)-ycurr%wood_energy(ico))/dtf
              
              qwflxwc_tot = qwflxwc_tot + qwflxwc
              hflxwc_tot  = hflxwc_tot  + hflxwc
           end if  !(if(wood_resolvable))
        end do
     end if
     
  case(1)  !select(ibranch_thermo)
     
!     print*,""

     do ico=1,cpatch%ncohorts
        if (resolve(ico)) then
           
           id_tveg=id_tveg+1

           ! Update with the new temperature
           ynext%leaf_temp(ico) = Y(id_tveg)
           ynext%wood_temp(ico) = Y(id_tveg)
           
           if(ynext%leaf_temp(ico) < t3ple8) then
              ynext%leaf_fliq(ico) = 0.d0
              ynext%wood_fliq(ico) = 0.d0
           else
              ynext%leaf_fliq(ico) = 1.d0
              ynext%wood_fliq(ico) = 1.d0
           end if

           ! Back calculate the latent and sensible heat fluxes of leaves
           ! ========================================================================

           ynext%veg_energy(ico)  = 0.0
           ynext%wood_energy(ico) = 0.0
           ynext%leaf_energy(ico) = 0.0

           if (ycurr%leaf_resolvable(ico) .or. ycurr%wood_resolvable(ico)) then

              ynext%leaf_energy(ico) = cmtl2uext8(     &
                ycurr%leaf_hcap (ico),              &
                ynext%leaf_water(ico),              &
                ynext%leaf_temp (ico),              &
                ynext%leaf_fliq (ico) )
              
              ynext%veg_energy(ico)=ynext%leaf_energy(ico)

              qwflxlc = ycurr%wflxlc(ico) * &
                   tq2enthalpy8(ycurr%leaf_temp(ico),1.d0,.true.)        
              
              hflxlc  = ycurr%hflx_lrsti(ico) - qwflxlc - &
                   (ynext%leaf_energy(ico)-ycurr%leaf_energy(ico))/dtf
              
              qwflxlc_tot = qwflxlc_tot + qwflxlc
              hflxlc_tot  = hflxlc_tot  + hflxlc

              ynext%wood_energy(ico) = cmtl2uext8(     &
                   ycurr%wood_hcap (ico),              &
                   ynext%wood_water(ico),              &
                   ynext%wood_temp (ico),              &
                   ynext%wood_fliq (ico) )   

              ynext%veg_energy(ico) = ynext%veg_energy(ico)+ &
                   ynext%wood_energy(ico)
              
              qwflxwc = ycurr%wflxwc(ico) * &
                   tq2enthalpy8(ycurr%wood_temp(ico),1.d0,.true.)        
              
              hflxwc  = ycurr%hflx_wrsti(ico) - qwflxwc - &
                   (ynext%wood_energy(ico)-ycurr%wood_energy(ico))/dtf
              
              qwflxwc_tot = qwflxwc_tot + qwflxwc
              hflxwc_tot  = hflxwc_tot  + hflxwc
           end if

           ! How does the temperature of the vegetation derived from
           ! its energy compare?


           call uextcm2tl8(ynext%veg_energy(ico),ynext%veg_water(ico)  &
                ,ycurr%veg_hcap(ico),veg_temp,fliq(ico))

           if(abs(veg_temp-ynext%leaf_temp(ico))>1d-10)then
              print*,"ISSUE IN ENERGY CONSERVATION,BDF2"
              stop
           end if
              


        end if
     end do
  end select


  ! Update eulerian based budget fluxes
  ! ============================================================
  

  eflxac = rk4aux(ibuff)%hcapcan*(ynext%can_enthalpy-ycurr%can_enthalpy)/dtf  - &
       (dydt%avg_sensible_gc + ycurr%qwflxsc + ycurr%qwflxgc    + &
       hflxlc_tot + qwflxlc_tot + qtransp_tot + hflxwc_tot + qwflxwc_tot)

!  print*,hflxlc_tot,qwflxlc_tot,qtransp_tot,hflxwc_tot,qwflxwc_tot
  
  qwflxac = ycurr%wflxac * tq2enthalpy8(0.5*(ycurr%can_temp+Ta),1.d0,.true.)
  
  hflxac  = eflxac-qwflxac

  if(checkbudget)then 
     
     ! Remove the previous integration
     ynext%ebudget_loss2atm   = ynext%ebudget_loss2atm            - &
          dydt%ebudget_loss2atm*dtf
     
     ! Add the new integration
     dydt%ebudget_loss2atm    = -eflxac
     ynext%ebudget_loss2atm   = ynext%ebudget_loss2atm - eflxac*dtf
     
  end if
  
  if (fast_diagnostics .or. checkbudget ) then     

     ! Update the sensible heat flux diagnostic
     ynext%avg_sensible_ac = ynext%avg_sensible_ac - dydt%avg_sensible_ac*dtf
     ynext%avg_sensible_ac = ynext%avg_sensible_ac + (hflxac)*dtf

  end if


  ! Make corrections to sensible heat flux and tstar
  ! ======================================================
  
  ! Remove the previous increment
  ynext%tpwp      = ynext%tpwp-dydt%tpwp*dtf
  ynext%avg_tstar = ynext%avg_tstar-dydt%tstar*dtf
  

  ! Make the current increment
  ynext%avg_tstar = ynext%tstar + &
       dtf*(hflxac/(rhoc*ycurr%ustar*ycurr%can_exner))
   

  ynext%tpwp = ynext%tpwp -(hflxac/(rhoc*ycurr%can_exner))*dtf
  

  ! Free memory
  deallocate(Y,Yf,A,B)


  return
end subroutine bdf2_solver

!================================================================

subroutine selective_gaussian_2body(ad,yd,xd,np)

  implicit none
  
  integer, intent(in)                       :: np
  real(kind=8), dimension(np),intent(out)   :: yd
  real(kind=8), dimension(np,np),intent(in) :: ad
  real(kind=8), dimension(np),intent(in)    :: xd
  real(kind=8)                              :: ydl
  integer                                   :: i


  
  if (np>1) then

     yd(1) = xd(1)
     ydl   = ad(1,1)

     do i=2,np
        yd(1) = yd(1)-ad(1,i)*xd(i)/ad(i,i)
        ydl = ydl-ad(1,i)*ad(i,1)/ad(i,i)
     end do

     ! Solved the first term by gaussian substitution
     yd(1) = yd(1)/ydl

     ! Now do the remaining terms

     do i=2,np
        yd(i) = (xd(i)-ad(i,1)*yd(1))/ad(i,i)
     end do

  else

     ! Trivial solution, no cohorts

     yd(1) = xd(1)/ad(1,1)

  end if



  return
end subroutine selective_gaussian_2body

!================================================================


subroutine ludcmp_dble(ad,n,np,indx,d)

  implicit none
  
  real(kind=8), parameter :: tiny_offset = 1.0d-32
  integer, intent(in) :: n
  integer, intent(in) :: np
  real, intent(out) :: d
!  real, (kind=8),dimension(np,np), intent(inout) :: a
  real(kind=8), dimension(np,np),intent(inout) :: ad
  integer, dimension(n), intent(out) :: indx
  real(kind=8), dimension(np) :: vv
  integer :: i
  integer :: j
  integer :: k
  integer :: imax
  real(kind=8) :: sum
  real(kind=8) :: aamax
  real(kind=8) :: dum

!  ad = dble(a)

  d = 1.0

  do i = 1, n
     aamax = 0.0d0
     do j = 1, n
        if(abs(ad(i,j)) > aamax)aamax = abs(ad(i,j))
     enddo
     if(aamax == 0.0d0)then
        print*,'singular matrix in ludcmp'
        do j=1,n
           print*,i,j,abs(ad(i,j)),aamax
        enddo
        stop
     endif
     vv(i) = 1.0d0 / aamax
  enddo

  do j = 1, n
     if(j.gt.1)then
        do i = 1, j - 1
           sum = ad(i,j)
           if(i > 1)then
              do k = 1, i - 1
                 sum = sum - ad(i,k) * ad(k,j)
              enddo
              ad(i,j) = sum
           endif
        enddo
     endif
     aamax = 0.0d0
     do i=j,n
        sum = ad(i,j)
        if (j > 1)then
           do k = 1, j - 1
              sum = sum - ad(i,k) * ad(k,j)
           enddo
           ad(i,j) = sum
        endif
        dum = vv(i) * abs(sum)
        if(dum >= aamax)then
           imax = i
           aamax = dum
        endif
     enddo
     if(j /= imax)then
        do k = 1, n
           dum = ad(imax,k)
           ad(imax,k) = ad(j,k)
           ad(j,k) = dum
        enddo
        d = -d
        vv(imax) = vv(j)
     endif
     indx(j) = imax
     if(j /= n)then
        if(ad(j,j) == 0.0d0) ad(j,j) = tiny_offset
        dum = 1.0d0 / ad(j,j)
        do i = j + 1, n
           ad(i,j) = ad(i,j) * dum
        enddo
     endif
  enddo
  if(ad(n,n) == 0.0d0)ad(n,n) = tiny_offset

!  a = real(ad)

  return
end subroutine ludcmp_dble

!=============================================================

subroutine lubksb_dble(ad,n,np,indx,bd)
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: np
  integer, dimension(n), intent(in) :: indx
!  real(, dimension(n), intent(inout) :: b
  real(kind=8), dimension(n),intent(inout) :: bd
!  real(kind=8), dimension(np,np), intent(in) :: a
  real(kind=8), dimension(np,np),intent(in) :: ad
  integer :: ii
  integer :: i
  integer :: ll
  real(kind=8) :: sum
  integer :: j

!  ad = dble(a)
!  bd = dble(b)

  ii = 0
  do i=1,n
     ll = indx(i)
     sum = bd(ll)
     bd(ll) = bd(i)
     if(ii /= 0)then
        do j=ii,i-1
           sum = sum - ad(i,j) * bd(j)
        enddo
     elseif(sum /= 0.0d0)then
        ii = i
     endif
     bd(i) = sum
  enddo
  do i=n,1,-1
     sum = bd(i)
     if ( i < n )then
        do j=i+1,n
           sum = sum - ad(i,j) * bd(j)
        enddo
     endif
     bd(i) = sum / ad(i,i)
  enddo

!  b = real(bd)

  return
end subroutine lubksb_dble



! Updated 10/24/2001.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE MIGS (A,N,X,INDX)
!
! Subroutine to invert matrix A(N,N) with the inverse stored
! in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  REAL, INTENT (INOUT), DIMENSION (N,N):: A
  REAL, INTENT (OUT), DIMENSION (N,N):: X
  REAL, DIMENSION (N,N) :: B
!
  DO I = 1, N
    DO J = 1, N
      B(I,J) = 0.0
    END DO
  END DO
  DO I = 1, N
    B(I,I) = 1.0
  END DO
!
  CALL ELGS (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      DO K = 1, N
        B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
      END DO
    END DO
  END DO
!
  DO I = 1, N
    X(N,I) = B(INDX(N),I)/A(INDX(N),N)
    DO J = N-1, 1, -1
      X(J,I) = B(INDX(J),I)
      DO K = J+1, N
        X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
      END DO
      X(J,I) =  X(J,I)/A(INDX(J),J)
    END DO
  END DO
END SUBROUTINE MIGS
!
SUBROUTINE ELGS (A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K,ITMP
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  REAL :: C1,PI,PI1,PJ
  REAL, INTENT (INOUT), DIMENSION (N,N) :: A
  REAL, DIMENSION (N) :: C
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Find the rescaling factors, one from each row
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      C1 = AMAX1(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! Search the pivoting (largest) element from each column
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
END SUBROUTINE ELGS





