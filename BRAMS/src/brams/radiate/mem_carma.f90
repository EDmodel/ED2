!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_carma

  use grid_dims

  implicit none

  type carma_v

     real, pointer, dimension(:,:,:) :: &
          aot

  end type carma_v

  type (carma_v), allocatable :: carma(:), carma_m(:)

  !Used in radcomp_carma (radcom) - radriv
  real :: tdec
  real :: sdec
  real :: cdec
  real :: declin
  real :: rvr(nzpmax)
  real :: rtr(nzpmax)
  real :: dn0r(nzpmax)
  real :: pird(nzpmax)
  real :: prd(nzpmax)
  real :: temprd(nzpmax+1)  
  real :: fthrl(nzpmax)
  real :: dzmr(nzpmax)
  real :: dztr(nzpmax)
  real :: fthrs(nzpmax) 
  
  !Used for vetorization
  integer :: p_isize
  integer :: p_jsize

  real, allocatable, dimension(:)       :: p_surf 
  real, allocatable, dimension(:)       :: p_top  
  real, allocatable, dimension(:)       :: t_surf 
  real, allocatable, dimension(:)       :: tabove_aerad 
  real, allocatable, dimension(:)       :: solnet
  real, allocatable, dimension(:)       :: xirdown
    
  real, allocatable, dimension(:,:)     :: p      
  real, allocatable, dimension(:,:)     :: t      
  real, allocatable, dimension(:,:)     :: t_aerad 
  real, allocatable, dimension(:,:)     :: p_aerad 
  real, allocatable, dimension(:,:)     :: qv_aerad
!kmlnew
  real, allocatable, dimension(:,:)     :: LWL_aerad
  real, allocatable, dimension(:,:)     :: IWL_aerad
  real, allocatable, dimension(:,:)     :: LWP_aerad
  real, allocatable, dimension(:,:)     :: IWP_aerad
  real, allocatable, dimension(:)       :: xland_aerad
!srf  REAL, ALLOCATABLE, DIMENSION(:,:)     :: RAIN_aerad
  real, allocatable, dimension(:)       :: RAIN_aerad
!kmlnew   
  real, allocatable, dimension(:,:)     :: press   
  real, allocatable, dimension(:,:)     :: dpg     
  real, allocatable, dimension(:,:)     :: tt      
  real, allocatable, dimension(:,:)     :: rhoa  
  real, allocatable, dimension(:,:)     :: rsfx   
  real, allocatable, dimension(:,:)     :: heats_aerad 
  real, allocatable, dimension(:,:)     :: heati_aerad 
  real, allocatable, dimension(:,:)     :: rdh2o   
  real, allocatable, dimension(:,:)     :: ptempg  
  real, allocatable, dimension(:,:)     :: ptempt  
  real, allocatable, dimension(:,:)     :: u1i

  real, allocatable, dimension(:,:,:)   :: gc     
  real, allocatable, dimension(:,:,:)   :: pah2o  
  real, allocatable, dimension(:,:,:)   :: paco2  
  real, allocatable, dimension(:,:,:)   :: pao2   
  real, allocatable, dimension(:,:,:)   :: pao3   
  real, allocatable, dimension(:,:,:)   :: taugas 
  real, allocatable, dimension(:,:,:)   :: paray  
  real, allocatable, dimension(:,:,:)   :: tauaer 
  real, allocatable, dimension(:,:,:)   :: taul   
  real, allocatable, dimension(:,:,:)   :: taucld 
  real, allocatable, dimension(:,:,:)   :: wcld   
  real, allocatable, dimension(:,:,:)   :: gcld   
  real, allocatable, dimension(:,:,:)   :: w0   
  real, allocatable, dimension(:,:,:)   :: g0
!kmlnew
  real, allocatable, dimension(:,:,:)   :: taucldlw 
  real, allocatable, dimension(:,:,:)   :: wolc   
  real, allocatable, dimension(:,:,:)   :: gl 
  real, allocatable, dimension(:,:,:)   :: taucldice 
  real, allocatable, dimension(:,:,:)   :: woice   
  real, allocatable, dimension(:,:,:)   :: gice 
    
!kmlnew    
  real, allocatable, dimension(:,:,:)   :: opd   
  real, allocatable, dimension(:,:,:)   :: uopd   
  real, allocatable, dimension(:,:,:)   :: ptemp  
  real, allocatable, dimension(:,:,:)   :: slope  
  real, allocatable, dimension(:,:,:)   :: b1
  real, allocatable, dimension(:,:,:)   :: b2
  real, allocatable, dimension(:,:,:)   :: b3
  real, allocatable, dimension(:,:,:)   :: ak
  real, allocatable, dimension(:,:,:)   :: gami
  real, allocatable, dimension(:,:,:)   :: ee1
  real, allocatable, dimension(:,:,:)   :: el1
  real, allocatable, dimension(:,:,:)   :: em1
  real, allocatable, dimension(:,:,:)   :: el2
  real, allocatable, dimension(:,:,:)   :: em2
  real, allocatable, dimension(:,:,:)   :: af
  real, allocatable, dimension(:,:,:)   :: bf
  real, allocatable, dimension(:,:,:)   :: ef
  real, allocatable, dimension(:,:,:)   :: cp
  real, allocatable, dimension(:,:,:)   :: cpb
  real, allocatable, dimension(:,:,:)   :: cmb
  real, allocatable, dimension(:,:,:)   :: ck1
  real, allocatable, dimension(:,:,:)   :: ck2
  real, allocatable, dimension(:,:,:)   :: fnet
  real, allocatable, dimension(:,:,:)   :: tmi
  real, allocatable, dimension(:,:,:)   :: tmid
  real, allocatable, dimension(:,:,:)   :: tmiu
  real, allocatable, dimension(:,:,:)   :: direc
  real, allocatable, dimension(:,:,:)   :: directu
    
  real, allocatable, dimension(:,:,:,:) :: pc     
  real, allocatable, dimension(:,:,:,:) :: pc_aerad 
  real, allocatable, dimension(:,:,:,:) :: caer    
  real, allocatable, dimension(:,:,:,:) :: y3    

  integer,allocatable,dimension(:) :: isl_aerad 
  integer,allocatable,dimension(:) :: lla
  integer,allocatable,dimension(:) :: lls
  
  real,parameter :: emisir_aerad = 1.0
  real,parameter :: h2ocol_aerad = 0.01

  integer :: ir_aerad

contains

  subroutine alloc_carma(car,ng,nmxp,nmyp,nw)

    implicit none

    type (carma_v), intent(inout) :: car(:)
    integer, intent(in) :: ng,nw,nmxp,nmyp

    allocate (car(ng)%aot(nmxp,nmyp,nw))

    return
  end subroutine alloc_carma

  !---------------------------------------------------------------

  subroutine nullify_carma(car,ng)

    implicit none

    type (carma_v), intent(inout) :: car(:)
    integer :: ng

    if (associated(car(ng)%aot ))  nullify (car(ng)%aot )


    return
  end subroutine nullify_carma
  !---------------------------------------------------------------

  subroutine dealloc_carma(car,ng)

    implicit none

    type (carma_v), intent(inout) :: car(:)
    integer :: ng

    if (associated(car(ng)%aot ))  deallocate (car(ng)%aot )

    return
  end subroutine dealloc_carma

  !-----------------------------------------------------------------
  subroutine zero_carma(car,ng)

    implicit none

    type (carma_v), intent(inout) :: car(:)
    integer,intent(in) :: ng

    car(ng)%aot(:,:,:)=0.

  end subroutine zero_carma

  !---------------------------------------------------------------

  subroutine filltab_carma(cv,cvm,ng,imean,n1,n2,n3)

    implicit none

    integer, intent(in) :: ng,n1,n2,n3,imean
    integer :: npts
    character (len=7) :: sname
    type (carma_v) :: cv,cvm

    if (associated(cv%aot)) then
       npts=n1*n2*n3
       write(sname,'(a4)') 'AOT'
       call vtables2 (cv%aot(1,1,1),cvm%aot(1,1,1),ng, &
            npts, imean,sname//' :7:hist:anal:mpti:mpt3')
       ! Not necessary MPT1 - Comunication NODE to NODE on DTLONG
       ! Radiation is a Column oriented process
    endif

  end subroutine filltab_carma

  !---------------------------------------------------------------

  subroutine init_carma(ia,iz,ja,jz,m1,m2,m3)
    use mem_globrad, only: ntotal,nlayer,ndbl,ngauss,nrad
    use mem_aerad, only: nz_rad,ngas,nelem,ngroup,nbin
  
    integer,intent(IN) :: ia,iz,ja,jz,m1,m2,m3
    integer :: iend
  
    iend=(iz-ia+1)*(jz-ja+1)
    !2d
    allocate(p_surf(iend)) 
    allocate(p_top(iend))  
    allocate(t_surf(iend)) 
    allocate(tabove_aerad(iend))  
    allocate(solnet(iend))
    allocate(xirdown(iend))
    !3d  
    allocate(p(iend,m1))      
    allocate(t(iend,m1))      
    allocate(t_aerad(iend,m1)) 
    allocate(p_aerad(iend,m1)) 
    allocate(qv_aerad(iend,m1)) 
    !kmlnew
    allocate(LWL_aerad(iend,m1)) 
    allocate(IWL_aerad(iend,m1)) 
    allocate(LWP_aerad(iend,m1)) 
    allocate(IWP_aerad(iend,m1)) 
    allocate(xland_aerad(iend)) 
    !srf  ALLOCATE(RAIN_aerad(iend,m1))
    allocate(RAIN_aerad(iend))
    !kmlnew

    allocate(press(iend,nlayer))   
    allocate(dpg(iend,nlayer))    
    allocate(tt(iend,nlayer))     
    allocate(rdh2o(iend,nlayer))     
    allocate(rhoa(iend,m1))   
    allocate(rsfx(iend,ntotal))  
    allocate(heats_aerad(iend,nz_rad)) 
    allocate(heati_aerad(iend,nz_rad)) 
    allocate(ptempg(iend,ntotal))  
    allocate(ptempt(iend,ntotal)) 
    allocate(u1i(iend,ntotal)) 
    !4d
    allocate(gc(iend,m1,ngas))     
    allocate(pah2o(iend,ntotal,nlayer)) 
    allocate(paco2(iend,ntotal,nlayer)) 
    allocate(pao2(iend,ntotal,nlayer))  
    allocate(pao3(iend,ntotal,nlayer))  
    allocate(taugas(iend,ntotal,nlayer))
    allocate(paray(iend,ntotal,nlayer)) 
    allocate(tauaer(iend,ntotal,nlayer)) 
    allocate(taul(iend,ntotal,nlayer))  
    allocate(taucld(iend,ntotal,nlayer)) 
    allocate(wcld(iend,ntotal,nlayer))   
    allocate(gcld(iend,ntotal,nlayer)) 
    !kmlnew
    allocate(taucldlw (iend,ntotal,nlayer)) 
    allocate(wolc     (iend,ntotal,nlayer))   
    allocate(gl       (iend,ntotal,nlayer)) 
    allocate(taucldice(iend,ntotal,nlayer)) 
    allocate(woice    (iend,ntotal,nlayer))   
    allocate(gice     (iend,ntotal,nlayer)) 
    !kmlnew
    allocate(w0(iend,ntotal,nlayer))   
    allocate(g0(iend,ntotal,nlayer))   
    allocate(opd(iend,ntotal,nlayer))  
    allocate(uopd(iend,ntotal,nlayer)) 
    allocate(ptemp(iend,ntotal,nlayer))
    allocate(slope(iend,ntotal,nlayer))  
    allocate(b1(iend,ntotal,nlayer))
    allocate(b2(iend,ntotal,nlayer))
    allocate(b3(iend,ntotal,nlayer))
    allocate(ak(iend,ntotal,nlayer))
    allocate(gami(iend,ntotal,nlayer))
    allocate(ee1(iend,ntotal,nlayer))
    allocate(el1(iend,ntotal,nlayer))
    allocate(em1(iend,ntotal,nlayer))
    allocate(el2(iend,ntotal,nlayer))
    allocate(em2(iend,ntotal,nlayer))
    allocate(af(iend,ntotal,ndbl))
    allocate(bf(iend,ntotal,ndbl))
    allocate(ef(iend,ntotal,ndbl))
    allocate(cp(iend,ntotal,nlayer))
    allocate(cpb(iend,ntotal,nlayer))
    allocate(cmb(iend,ntotal,nlayer))
    allocate(ck1(iend,ntotal,nlayer))
    allocate(ck2(iend,ntotal,nlayer))
    allocate(fnet(iend,ntotal,nlayer))
    allocate(tmi(iend,ntotal,nlayer ))
    allocate(tmid(iend,ntotal,nlayer))
    allocate(tmiu(iend,ntotal,nlayer))
    allocate(direc(iend,ntotal,nlayer))
    allocate(directu(iend,ntotal,nlayer))
    !5d  
    allocate(pc(iend,m1,nbin,nelem))  
    allocate(pc_aerad(iend,nz_rad,nbin,ngroup)) 
    allocate(caer(iend,nlayer,nrad,ngroup))    
    allocate(y3(iend,ntotal,ngauss,nlayer))    
    !2d int
    allocate(isl_aerad(iend)) 
    allocate(lla(iend))
    allocate(lls(iend))

    !Zero all variables
    p_surf=0.0 
    p_top=0.0  
    t_surf=0.0 
    tabove_aerad=0.0  
    lla=0
    lls=0
    solnet=0.0
    xirdown=0.0
    !3d  
    p=0.0      
    t=0.0      
    t_aerad=0.0 
    p_aerad=0.0 
    qv_aerad=0.0 
    !kmlnew
    LWL_aerad=0.0 
    IWL_aerad=0.0 
    LWP_aerad=0.0 
    IWP_aerad=0.0
    xland_aerad=0.0
    RAIN_aerad=0.0 
    !kmlnew
    press=0.0   
    dpg=0.0    
    tt=0.0     
    rdh2o=0.0     
    rhoa=0.0   
    rsfx=0.0  
    heats_aerad=0.0 
    heati_aerad=0.0 
    ptempg=0.0 
    ptempt=0.0 
    u1i=0.0
    !4d
    gc=0.0     
    pah2o=0.0 
    paco2=0.0 
    pao2=0.0  
    pao3=0.0  
    taugas=0.0
    paray=0.0 
    tauaer=0.0 
    taul=0.0  
    taucld=0.0
    !kmlnew
    taucldlw=0.0
    wolc=0.0 
    gl=0.0 
    taucldice=0.0
    woice=0.0	
    gice=0.0 
    !kmlnew  
    wcld=0.0   
    gcld=0.0   
    w0=0.0   
    g0=0.0   
    opd=0.0  
    uopd=0.0 
    ptemp=0.0
    slope=0.0  
    b1=0.0
    b2=0.0
    b3=0.0
    ak=0.0
    gami=0.0
    ee1=0.0
    el1=0.0
    em1=0.0
    el2=0.0
    em2=0.0
    af=0.0
    bf=0.0
    ef=0.0
    cp=0.0
    cpb=0.0
    cmb=0.0
    ck1=0.0
    ck2=0.0
    fnet=0.0
    tmi=0.0
    tmid=0.0
    tmiu=0.0
    direc=0.0
    directu=0.0
    !5d  
    pc=0.0  
    pc_aerad=0.0 
    caer=0.0    
    y3=0.0    
    !2d int
    isl_aerad=0    
  end subroutine init_carma

  !---------------------------------------------------------------

  subroutine end_carma()

    !2d
    deallocate(p_surf) 
    deallocate(p_top)  
    deallocate(t_surf) 
    deallocate(tabove_aerad)  
    deallocate(lla)
    deallocate(lls)
    deallocate(solnet)
    deallocate(xirdown)
    !3d  
    deallocate(p)      
    deallocate(t)      
    deallocate(t_aerad) 
    deallocate(p_aerad) 
    deallocate(qv_aerad) 
    !kmlnew
    deallocate(LWL_aerad) 
    deallocate(IWL_aerad) 
    deallocate(LWP_aerad) 
    deallocate(IWP_aerad) 
    deallocate(xland_aerad)
    deallocate(RAIN_aerad)
    !kmlnew
    deallocate(press)   
    deallocate(dpg)    
    deallocate(tt)     
    deallocate(rdh2o)     
    deallocate(rhoa)   
    deallocate(rsfx)  
    deallocate(heats_aerad) 
    deallocate(heati_aerad) 
    deallocate(ptempg)  
    deallocate(ptempt) 
    deallocate(u1i) 
    !4d
    deallocate(gc)     
    deallocate(pah2o) 
    deallocate(paco2) 
    deallocate(pao2)  
    deallocate(pao3)  
    deallocate(taugas)
    deallocate(paray) 
    deallocate(tauaer) 
    deallocate(taul)  
    deallocate(taucld) 
    deallocate(wcld)   
    deallocate(gcld)
    !kmlnew
    deallocate(taucldlw) 
    deallocate(wolc)   
    deallocate(gl) 
    deallocate(taucldice) 
    deallocate(woice)   
    deallocate(gice) 
    !kmlnew

    deallocate(w0)   
    deallocate(g0)   
    deallocate(opd)  
    deallocate(uopd) 
    deallocate(ptemp)
    deallocate(slope)  
    deallocate(b1)
    deallocate(b2)
    deallocate(b3)
    deallocate(ak)
    deallocate(gami)
    deallocate(ee1)
    deallocate(el1)
    deallocate(em1)
    deallocate(el2)
    deallocate(em2)
    deallocate(af)
    deallocate(bf)
    deallocate(ef)
    deallocate(cp)
    deallocate(cpb)
    deallocate(cmb)
    deallocate(ck1)
    deallocate(ck2)
    deallocate(fnet)
    deallocate(tmi)
    deallocate(tmid)
    deallocate(tmiu)
    deallocate(direc)
    deallocate(directu)
    !5d  
    deallocate(pc)  
    deallocate(pc_aerad) 
    deallocate(caer)    
    deallocate(y3)    
    !2d int
    deallocate(isl_aerad) 
  end subroutine end_carma

end module mem_carma
