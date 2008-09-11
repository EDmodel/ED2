!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine timestep()

  use mem_basic, only: &
       basic_g  ! intent(inout)

  use node_mod, only: &
       mzp, mxp, myp,  & ! intent(in)
       ia, iz, ja, jz, & ! intent(in)
       i0, j0,         & ! intent(in)
       izu, jzv,       & ! intent(in)
       mynum,          & ! intent(in)
       ibcon,          & ! intent(in)
       ipara             ! intent(in)

  !use mem_radiate, only: ! Not used

  use mem_cuparm, only: &
       nnqparm,         & ! intent(in)
       if_cuinv           ! intent(in)

  use mem_varinit, only: &
       nud_type ! intent(in)

  use mem_turb, only: &
       if_urban_canopy, & ! intent(in)
       ihorgrad           ! intent(in)

  use mem_oda,   only: &
       if_oda ! intent(in)

  use micphys,   only: &
       level ! intent(in)

  use mem_grid, only: &
       ngrid,      & ! intent(in)
       time,       & ! intent(in)
       dtlong,     & ! intent(in)
       dtlongn,    & ! intent(in)
       iyeara,     & ! intent(in)
       imontha,    & ! intent(in)
       idatea,     & ! intent(in)
       grid_g,     & ! intent(inout)
       nxtnest,    & ! intent(in)
       if_adap,    & ! intent(in)
       dtlt,       & ! intent(in)
       istp,       & ! intent(in)
       jdim,       & ! intent(in)
       nzp,        & ! intent(in)
       f_thermo_e, & ! intent(in)
       f_thermo_w, & ! intent(in)
       f_thermo_s, & ! intent(in)
       f_thermo_n    ! intent(in)

  use shcu_vars_const, only: & ! For Shallow Cumulus Paramet.
       nnshcu                ! ! intent(in)

  use mem_scalar, only: & ! For SiB
       scalar_g ! intent(in)

  use mem_leaf, only: & ! For SiB
       isfcl ! intent(in)

  ! CATT
  use catt_start, only: &
       catt ! INTENT(IN)

  ! For CATT
  !!use burns, only : queimadas
  use emission_source_map, only: &
       burns ! Subroutine

 ! TEB_SPM
  use teb_spm_start, only: &
       TEB_SPM ! INTENT(IN)

  ! For TEB_SPM
  use mem_emiss, only: &
       ichemi,         & ! INTENT(IN)
       isource           ! INTENT(IN)

  ! ALF
  ! Necessary in new advection scheme
  use advect_kit, only : &
       calc_advec          ! Subroutine

  ! For specific optimization depending the type of machine
!  use machine_arq, only: &
!       machine ! INTENT(IN)

  ! [MLO - For Exner function and massflux computation:
  use mem_mass, only : iexev, imassflx
  ! MLO]

  implicit none

  real :: t1,w1
  real, external :: cputime

  integer, save :: ncall=0

  logical, parameter :: acct = .false. ! To Activate acctimes
  
  logical, parameter :: banneron=.false. ! To print position banner


  

  if (acct) call acctimes('init',0,' ',t1,w1)

  !        +-------------------------------------------------------------+
  !        |   Timestep driver for the hybrid non-hydrostatic time-split |
  !        |      model.                                                 |
  !        +-------------------------------------------------------------+


  ! MLO - Initialize microphysics tables if the level is 3. This should probably 
  !       be placed outside the subroutine timestep for efficiency, but I will 
  !       leave it here.
  if (level == 3) call micro_1st()

  !  Zero out all tendency arrays.   
  !--------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling tend0...'
  call TEND0()          
  if (acct) call acctimes('accu',1,'TEND0',t1,w1)

  if (CATT==1) then
     !#########################################################################
     !srf-fev-2003-large and subgrid scale forcing for shallow and deep cumulus
     !if( NNQPARM(ngrid).eq.2 .and. NNSHCU(ngrid)==1 ) call prepare_lsf(1)
     !#########################################################################

     t1=cputime(w1)
     !--------------------------------srf
     !!call queimadas(mzp,mxp,myp,ia,iz,ja,jz)
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling burns...'
     call burns(ngrid, mzp, mxp, myp, ia, iz, ja, jz, &
       scalar_g, time, iyeara, imontha, idatea)
     !--------------------------------srf
     if (acct) call acctimes('accu',2,'BURNS',t1,w1)
  endif

  !  Thermodynamic diagnosis   
  !--------------------------------
  t1=cputime(w1)
  if (level  /=  3) then
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling thermo(supsat)...'
     call THERMO(mzp,mxp,myp,ia,iz,ja,jz,'SUPSAT') 
  endif
  if (acct) call acctimes('accu',3,'THERMO',t1,w1)

  !--------Medvigy's mass conservation fix ---------------
  if (iexev == 2) then
     t1=cputime(w1)
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling exevolve(adv)...'
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV')
     if (acct) call acctimes('accu',4,'EXEVOLVE_ADV',t1,w1)
  end if
  !------------------------------------------------------- MLO - SRF]
  
  !  Radiation parameterization
  !--------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling radiate...'
  call RADIATE(mzp,mxp,myp,ia,iz,ja,jz,mynum) 
  if (acct) call acctimes('accu',5,'RADIATE',t1,w1)

  !  Surface layer, soil and veggie model
  !----------------------------------------
  t1=cputime(w1)
  if (isfcl <= 2) then
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling sfclyr...'
     call SFCLYR(mzp,mxp,myp,ia,iz,ja,jz,ibcon)
  elseif (isfcl == 3) then
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling sfclyr_sib...'
     call sfclyr_sib(mzp,mxp,myp,ia,iz,ja,jz,ibcon)
  endif
  

  if (acct) call acctimes('accu',6,'SFCLYR',t1,w1)

  ! CO2  bio source
  if (ISFCL == 3) then
     t1=cputime(w1)
     !! New version of SiB-BRAMS
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling CO2_biosource...'
     call co2_biosource(mzp,mxp,myp,ia,iz,ja,jz,ngrid,  &
          scalar_g(1,ngrid)%sclt(1),basic_g(ngrid)%dn0,grid_g(ngrid)%rtgt)
     if (acct) call acctimes('accu',7,'co2_biosource',t1,w1)
  endif
  !----------------------------------------

  if (CATT==1) then
     !#########################################################################
     !   print*,'- dry deposition'
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling drydep_driver...'
     call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz)     
     !#########################################################################
     !#########################################################################
     !srf-fev-2003-large and subgrid scale forcing for shallow and deep cumulus
     !if( NNQPARM(ngrid).eq.2 .and. NNSHCU(ngrid)==1 ) call prepare_lsf(2)
     !#########################################################################
  endif

  !  Send boundaries to adjoining nodes
  !-------------------------------------------
  t1=cputime(w1)
  if (ipara  ==  1) then
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling node_sendlbc...'
     call node_sendlbc()
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling node_sendcyclic...'
     if (ngrid  ==  1) call node_sendcyclic(1)
  endif
  if (acct) call acctimes('accu',8,'SendLBC',t1,w1)

  !  Coriolis terms
  !  ----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling corlos...'
  call CORLOS(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv) 
  if (acct) call acctimes('accu',9,'CORLOS',t1,w1)

  !  Velocity advection
  !----------------------------------------
  t1=cputime(w1)
  ! Use Optmized advection only in SX-6, for the moment
!  if (machine==0) then
     ! If Generic IA32 use old Advction Scheme
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling advectc(v)...'
     call ADVECTc('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!  elseif (machine==1) then
     ! Using optmized advection scheme only in SX-6
!     call calc_advec('V',ngrid,mzp,mxp,myp)
!  endif
  if (acct) call acctimes('accu',10,'ADVECTv',t1,w1)

  !  Cumulus parameterization
  !----------------------------------------
  !t1=cputime(w1)
  !if(NNQPARM(ngrid) == 1 .or. IF_CUINV == 1) call cuparm()      
  !if (acct) call acctimes('accu',7,'CUPARM',t1,w1)

  !  Urban canopy parameterization
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling urban_canopy...'
  if( IF_URBAN_CANOPY == 1) call urban_canopy()      
  if (acct) call acctimes('accu',11,'URBAN_CANOPY',t1,w1)

  !  Analysis nudging and boundary condition
  !------------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling datassim...'
  if(NUD_TYPE > 0) call DATASSIM()  
  if (acct) call acctimes('accu',12,'DATASSIM',t1,w1)

  !  Observation data assimilation 
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling oda_nudge...'
  if(IF_ODA == 1) call oda_nudge()  
  if (acct) call acctimes('accu',13,'DATASSIM',t1,w1)

  !  Nested grid boundaries
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling nstbdriv...'
  if(nxtnest(ngrid).ge.1) call nstbdriv()  
  if (acct) call acctimes('accu',14,'NSTBDRIV',t1,w1)

  !  Rayleigh friction for theta
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling rayft...'
  call RAYFT()           
  if (acct) call acctimes('accu',15,'RAYFT',t1,w1)

  !  Get the overlap region between parallel nodes
  !---------------------------------------------------
  if(ipara == 1) then      
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling node_getlbc...'
     call node_getlbc()  
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling node_getcyclic...'
     if (ngrid  ==  1) call node_getcyclic(1)
  endif
  if (acct) call acctimes('accu',16,'GETlbc',t1,w1)

  ! Exner function correction
  !----------------------------------------
  if (iexev == 2) then
     t1=cputime(w1)
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling thermo(micro)...'
     call thermo(mzp,mxp,myp,1,mxp,1,myp,'MICRO') 
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling exevolve(tha)...'
     call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA')
     if (acct) call acctimes('accu',17,'EXEVOLVE_THA',t1,w1)
  end if


  !  Sub-grid diffusion terms
  !----------------------------------------
  t1=cputime(w1)
  if ((if_adap==0).and.(ihorgrad==2)) then
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling diffuse_brams31...'
     call diffuse_brams31() !call optimized subroutine
  else
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling diffuse...'
     call diffuse()
  endif
  if (acct) call acctimes('accu',18,'DIFFUSE',t1,w1)

  !if (CATT==1) then
     !#########################################################################
     !srf-fev-2003-large and subgrid scale forcing for shallow and deep cumulus
     !srf - neste ponto LSF = radiation + pbl(vertical) diffusion,
     !srf - este eh o forcing para o shallow
     !if( NNQPARM(ngrid).eq.2 .and. NNSHCU(ngrid)==1 ) call prepare_lsf(3)
     !#########################################################################
  !endif

  !  Velocity advection
  !----------------------------------------
  t1=cputime(w1)
  ! Use Optmized advection only in SX-6, for the moment
 ! if (machine==0) then
     ! If Generic IA32 use old Advction Scheme
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling advectc(t)...'
     call ADVECTc('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
 ! elseif (machine==1) then
     ! Using optmized advection scheme only in SX-6
 !    if (ngrid<=2) then
 !       call calc_advec('T',ngrid,mzp,mxp,myp)
 !    else
 !       call ADVECTc('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
 !    endif
 ! endif
  if (acct) call acctimes('accu',19,'ADVECTs',t1,w1)

  if (imassflx == 1) then
    t1=cputime(w1)
    if (banneron) write(unit=*,fmt='(a)') '     [-] Calling prep_advflx_to_mass...'
    call prep_advflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,ngrid)
    if (acct) call acctimes('accu',20,'ADVEC_TO_MASS',t1,w1)
  end if

  !--------------------------------------------------------------------------------!
  !MLO - Cumulus parameterization.                                                 !
  !--------------------------------------------------------------------------------!
  t1=cputime(w1)
  if(NNQPARM(ngrid) == 1 .or. IF_CUINV == 1) call cuparm()      
  if(NNSHCU(ngrid) == 1) call SHCUPA()
  if (CATT==1) then
     if(nnshcu(ngrid)==2) call cuparm_grell_catt(2) 
     if(NNQPARM(ngrid)==2) CALL cuparm_grell_catt(1) 
  else if (CATT==0) then
     if(nnshcu(ngrid) == 2)  call cuparm_grell(2)
     if(nnqparm(ngrid) == 2) call cuparm_grell(1)
  end if
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling cuparm...'
  if (acct) call acctimes('accu',21,'CUPARM',t1,w1)
  !--------------------------------------------------------------------------------!

  if (TEB_SPM==1) then
     ! Update urban emissions
     !----------------------------------------
     if(isource==1)then
        t1=cputime(w1)
        if (banneron) write(unit=*,fmt='(a)') '     [-] Calling sources_teb...'
        call sources_teb(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt)
        if (acct) CALL acctimes('accu',22,'EMISS',t1,w1)
     endif
     !  Update chemistry
     !----------------------------------------
     if(ichemi==1)then
        t1=cputime(w1)
        if (banneron) write(unit=*,fmt='(a)') '     [-] Calling ozone...'
        call ozone(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt)
        if (acct) CALL acctimes('accu',23,'OZONE',t1,w1)
     endif
  endif

  !  Update scalars
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling predtr...'
  call PREDTR()          
  if (acct) call acctimes('accu',24,'PREDTR',t1,w1)

  !  Moisture variables positive definite
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling negadj1...'
  call negadj1(mzp,mxp,myp) 
  if (acct) call acctimes('accu',25,'NEGADJ1',t1,w1)

  !  Microphysics
  !----------------------------------------
  t1=cputime(w1)
  if (level==3) then
!     if (machine==1 .and. TEB_SPM==0) then
!        ! Optimized version only for SX-6
!        call micro_opt()
!     else
        ! Original Version used in a Generic IA32 machine
        if (banneron) write(unit=*,fmt='(a)') '     [-] Calling micro...'
        call micro()
!     endif
  endif
  if (acct) call acctimes('accu',26,'MICRO',t1,w1)

  !  Thermodynamic diagnosis
  !----------------------------------------
  t1=cputime(w1)
  if (level /= 3) then
     if (banneron) write(unit=*,fmt='(a)') '     [-] Calling thermo(micro)...'
     call THERMO(mzp,mxp,myp,1,mxp,1,myp,'MICRO') 
  endif
  if (acct) call acctimes('accu',27,'THERMO',t1,w1)


  !----- Medvigy's mass conservation fix ------
  ! Right before calling "TRSETS", add a function call:
  if (iexev == 2) then
      t1=cputime(w1)
      if (banneron) write(unit=*,fmt='(a)') '     [-] Calling exevolve(ths)...'
      call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS')
      if (acct) call acctimes('accu',28,'EXEVOLVE_THS',t1,w1)
  end if
  !-----------------------------------------------------------------------
  

  !  Apply scalar b.c.'s
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling trsets...'
  call TRSETS()          
  if (acct) call acctimes('accu',29,'TRSETS',t1,w1)

  !  Lateral velocity boundaries - radiative
  !-------------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling latbnd...'
  call LATBND()
  if (acct) call acctimes('accu',30,'LATBND',t1,w1)

  !  First stage Asselin filter
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling hadvance(1)...'
  call HADVANCE(1)     
  if (acct) call acctimes('accu',31,'HADVANCE',t1,w1)

  !  Buoyancy term for w equation
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling buoyancy...'
  call BUOYANCY()
  if (acct) call acctimes('accu',32,'BUOYANCY',t1,w1)

  !  Acoustic small timesteps
  !----------------------------------------
  t1=cputime(w1)
  !call ACOUSTIC()
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling acoustic_new...'
  call ACOUSTIC_new()
  if (acct) call acctimes('accu',33,'ACOUSTIC',t1,w1)

  !  Last stage of Asselin filter
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling hadvance(2)...'
  call HADVANCE(2)

  !  Velocity/pressure boundary conditions
  !----------------------------------------
  t1=cputime(w1)
  if (banneron) write(unit=*,fmt='(a)') '     [-] Calling vpsets...'
  call VPSETS()          
  if (acct) call acctimes('accu',34,'HADVANCE',t1,w1)

  if (acct) then
     if(mod(istp,4) == 0) then
        call acctimes('prin',-1,' ',t1,w1)
     endif
  endif

  ! Call THERMO on the boundaries
  ![MLO - The logical flags are vectors so I'sending the information about the grid only...
  call thermo_boundary_driver((time+dtlongn(ngrid)), dtlongn(ngrid), &
       f_thermo_e(ngrid), f_thermo_w(ngrid), f_thermo_s(ngrid), f_thermo_n(ngrid), &
       nzp, mxp, myp, jdim)


!!$  call mass_flux(nzp,nxp,nyp,mzp,mxp,myp,a(iup),a(ivp),a(iwp) &
!!$       ,a(idn0),a(irtgu),a(irtgv),a(idyu),a(idxv),a(ipp),a(ipi0))

  return
end subroutine timestep

!*************************************************************************

subroutine acctimes(action,num,string,t1,w1)

  use mem_all
  use node_mod

  implicit none
  integer :: num
  real :: t1,w1
  character(len=*) :: string,action

  include 'interface.h'

  integer, parameter :: num_times=100
  real, save :: rtimes(num_times),wtimes(num_times)
  character(len=8),save :: crtimes(num_times)

  real, external :: cputime,walltime,valugp
  real :: sumtime,pcpu,fsecs,cpuinc,ww
  integer :: i,j,npts,ip,jp,kp

  if(action(1:4)=='init') then
     crtimes(1:num_times)=' '
     rtimes(1:num_times)=0.
     wtimes(1:num_times)=0.
     basic_g(ngrid)%cputime(1:mxp,1:myp)=0.
  elseif(action(1:4)=='prin') then
     if(mynum/=2.and.mynum/=0) return
     sumtime=0.
     do i=1,num_times
        sumtime=sumtime+rtimes(i)
     enddo

     write(6,*) '======= total CPU =====',sumtime,'=========='
     do i=1,34
        write(6,'(a10,i4,''-'',a12,f10.3,f7.2,2f9.3)') &
             'Timings-',mynum,crtimes(i)  &
             ,rtimes(i),rtimes(i)/sumtime*100.,wtimes(i)  &
             ,wtimes(i)-rtimes(i)
     enddo
     sumtime=0.
     do j=1,myp
        do i=1,mxp
           sumtime=sumtime+basic_g(ngrid)%cputime(i,j)
        enddo
     enddo
     write(6,'(a,2i5,f10.5)') 'Total CPU secs -ngrid,mynum,secs:'  &
          ,ngrid,mynum,sumtime

  elseif(action(1:4)=='accu'.or.action(1:4)=='null') then
     ! Accumulate full times into tables
     pcpu=(cputime(ww)-t1)
     rtimes(num)=rtimes(num) + pcpu
     crtimes(num)=string

     fsecs=72559200.
!     wtimes(num)=wtimes(num)+(walltime(fsecs)-w1)
     wtimes(num)=wtimes(num)+(ww-w1)

     if(action(1:4)=='accu') then
        ! Divide cpu time equally amoung columns and accumulate in
        ! basic_g(ngrid)%cputime(1,1)
        npts=(iz-ia+1)*(jz-ja+1)
        cpuinc=pcpu/float(npts)
        do j=ja,jz
           do i=ia,iz
              basic_g(ngrid)%cputime(i,j)=basic_g(ngrid)%cputime(i,j)+cpuinc
           enddo
        enddo
     endif
  endif

  !    only here for debugging purposes
  ip=2
  jp=2
  kp=2
  !   if( (mynum == 1.or.mynum == 0).and.ngrid == 2) then
  ! if( (mynum == 1.or.mynum == 0)) then
  !    print '(a10,2i3,f9.1,12e14.7)',string,2,mynum,time  &
       !       ,radiate_g(ngrid)%rshort(ip,jp) &
  !       ,radiate_g(ngrid)%rlongup(ip,jp) 
  !       ,basic_g(ngrid)%up(kp,ip,jp)  &
       !       ,valugp(mzp,mxp,myp,kp,ip,jp,tend%ut(1))  &
  !       ,basic_g(ngrid)%vp(kp,ip,jp)  &
       !       ,valugp(mzp,mxp,myp,kp,ip,jp,tend%vt(1))  &
  !       ,basic_g(ngrid)%wp(kp,ip,jp)  &
       !       ,valugp(mzp,mxp,myp,kp,ip,jp,tend%wt(1)) &
  !       ,basic_g(ngrid)%pp(kp,ip,jp)  &
       !       ,valugp(mzp,mxp,myp,kp,ip,jp,tend%pt(1))  &
  !       ,basic_g(ngrid)%pi0(kp,ip,jp)  
  !       ,basic_g(ngrid)%theta(kp,ip,jp) &
       !       ,basic_g(ngrid)%thp(kp,ip,jp) 
  !       ,grid_g(ngrid)%flpw(ip,jp)  &
       !              ,basic_g(2)%rtp(kp,ip,jp)  &
  !              ,basic_g(2)%rv(kp,ip,jp)  &
  !              ,valugp(mzp,mxp,myp,kp,ip,jp,tend%rtt(1)) 

  !  do jp=22,1,-1
  !     print'(i3,20g12.4)',jp,(basic_g(ngrid)%pp(kp,ip,jp)),ip=1,8)
  !  enddo

  !  endif
  return
end subroutine acctimes

!*************************************************************************

subroutine mass_flux(n1,n2,n3,m1,m2,m3,up,vp,wp  &
     ,dn0,rtgu,rtgv,dyu,dxv,pp,pi0)

  use mem_grid
  use rconstants

  implicit none
  integer :: n1,n2,n3,m1,m2,m3
  real :: up(m1,m2,m3),vp(m1,m2,m3),wp(m1,m2,m3)  &
       ,dn0(n1,n2,n3),rtgu(n2,n3),dyu(n2,n3),dxv(n2,n3)  &
       ,rtgv(n2,n3),pp(m1,m2,m3),pi0(n1,n2,n3)

  real, save :: aintmass=0.

  integer :: i,j,k
  real :: wmass,emass,smass,nmass,prtot,tmass,ppp,area

  !cc      if (mod(time,300.).gt..1) return

  !  west/east bound
  wmass=0.
  emass=0.
  do j=2,nyp-1
     do k=2,nzp-1
        i=1
        wmass=wmass +  &
             up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i+1,j))*.5
        i=nxp-1
        emass=emass -  &
             up(k,i,j)*rtgu(i,j)/(dyu(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i+1,j))*.5
     enddo
  enddo

  !  north/south bound
  smass=0.
  nmass=0.
  do i=2,nxp-1
     do k=2,nzp-1
        j=1
        smass=smass +  &
             vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i,j+1))*.5
        j=nyp-1
        nmass=nmass -  &
             vp(k,i,j)*rtgv(i,j)/(dxv(i,j)*dzt(k))  &
             *(dn0(k,i,j)+dn0(k,i,j+1))*.5
     enddo
  enddo

  k=2
  prtot=0.
  do j=2,nyp-1
     do i=2,nxp-1
        ppp= ( (pp(k,i,j)+pi0(k,i,j))/cp )**cpor*p00
        prtot=prtot+ppp/(dyu(i,j)*dxv(i,j))
     enddo
  enddo


  tmass=wmass+emass+smass+nmass
  aintmass=aintmass+tmass*dtlong
  area=(nxp-2)*deltax*(nyp-2)*deltay


  print*,'==============================='
  print*,' Mass flux - W, E, S, N'
  print*,  wmass,emass,smass,nmass
  print*, 'total (kg/(m2 s):',tmass/area
  print*, 'total (kg/m2):',aintmass/area
  print*, 'total pr change (pa):',aintmass/area*9.8
  print*, 'computed mean press:',prtot/area
  print*,'==============================='

  return
end subroutine mass_flux
