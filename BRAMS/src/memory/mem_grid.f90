!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_grid

  use grid_dims

  type grid_vars

     ! Variables to be dimensioned by (nxp,nyp)

     real, pointer :: topt(:,:)
     real, pointer :: topu(:,:)
     real, pointer :: topv(:,:)
     real, pointer :: topm(:,:)
     real, pointer :: topma(:,:)
     real, pointer :: topta(:,:)
     real, pointer :: rtgt(:,:)
     real, pointer :: rtgu(:,:)
     real, pointer :: rtgv(:,:)
     real, pointer :: rtgm(:,:)
     real, pointer :: f13t(:,:)
     real, pointer :: f13u(:,:)
     real, pointer :: f13v(:,:)
     real, pointer :: f13m(:,:)
     real, pointer :: f23t(:,:)
     real, pointer :: f23u(:,:)
     real, pointer :: f23v(:,:)
     real, pointer :: f23m(:,:)
     real, pointer :: dxt(:,:)
     real, pointer :: dxu(:,:)
     real, pointer :: dxv(:,:)
     real, pointer :: dxm(:,:)
     real, pointer :: dyt(:,:)
     real, pointer :: dyu(:,:)
     real, pointer :: dyv(:,:)
     real, pointer :: dym(:,:)
     real, pointer :: fmapt(:,:)
     real, pointer :: fmapu(:,:)
     real, pointer :: fmapv(:,:)
     real, pointer :: fmapm(:,:)
     real, pointer :: fmapti(:,:)
     real, pointer :: fmapui(:,:)
     real, pointer :: fmapvi(:,:)
     real, pointer :: fmapmi(:,:)
     real, pointer :: glat(:,:)
     real, pointer :: glon(:,:)
     real, pointer :: topzo(:,:)

     !  Variables for the ADAP coordinate

     real, pointer :: aru(:,:,:)
     real, pointer :: arv(:,:,:)
     real, pointer :: arw(:,:,:)
     real, pointer :: volu(:,:,:)
     real, pointer :: volv(:,:,:)
     real, pointer :: volw(:,:,:)
     real, pointer :: volt(:,:,:)
     ! MLO - Transformed in real since lpu/lpv/lpw go through vtables2, which expect
     !       real numbers 
     real, pointer :: flpu(:,:)
     real, pointer :: flpv(:,:)
     real, pointer :: flpw(:,:)
     
     !MLO - Mynum, just for debugging purposes
     real, pointer :: fmynum(:,:)
  end type grid_vars


  type (grid_vars), allocatable :: grid_g(:)
  type (grid_vars), allocatable :: gridm_g(:)

  character(len=64) :: expnme            ! experiment name
  integer :: ngrids                      ! how many grids
  integer :: ngridsh
  integer :: nxtnest(maxgrds)            ! next coarser grid number (0 if grid is not nested)
  real(kind=8) :: time, timmax 
  real :: dtlongn(maxgrds)               ! delta t long

  integer, target :: nnxp(maxgrds)       ! grid cells at x direction
  integer :: nnx(maxgrds)                ! nnxp - 1
  integer :: nnx1(maxgrds)               ! nnxp - 2
  integer :: nnx2(maxgrds)               ! nnxp - 3
  integer :: nstratx(maxgrds)            ! nest ratio for next coarser grid
  real :: deltaxn(maxgrds)               ! delta x
  real :: xtn(nxpmax,maxgrds)            ! x coordinate of cell center on polar stereographic projection
  real :: xmn(nxpmax,maxgrds)            ! x coordinate of higher cell boundary on polar stereographic projection
  integer :: ninest(maxgrds)             ! index on next coarser grid where this grid starts (lower southwest corner)
  !                                      ! ds to interpolate x direction from coarser to finner grids: 
  integer :: ipm(nxpmax,maxgrds)         ! next coarser grid cell index (icoarser) that contains this finer grid cell
  real :: ei1(nxpmax,maxgrds)            ! for icoarser-1 on 3 points interpolation
  real :: ei2(nxpmax,maxgrds)            ! for icoarser   on 3 points interpolation
  real :: ei3(nxpmax,maxgrds)            ! for icoarser+1 on 3 points interpolation
  real :: ei4(nxpmax,maxgrds)            ! for icoarser-2 on 4 points interpolation
  real :: ei5(nxpmax,maxgrds)            ! for icoarser-1 on 4 points interpolation
  real :: ei6(nxpmax,maxgrds)            ! for icoarser   on 4 points interpolation
  real :: ei7(nxpmax,maxgrds)            ! for icoarser+1 on 4 points interpolation

  integer, target :: nnyp(maxgrds)       ! grid cells at y direction
  integer :: nny(maxgrds)                ! nnyp - 1
  integer :: nny1(maxgrds)               ! nnyp - 2
  integer :: nny2(maxgrds)               ! nnyp - 3
  integer :: nstraty(maxgrds)            ! nest ratio for next coarser grid
  real :: deltayn(maxgrds)               ! delta y
  real :: ytn(nypmax,maxgrds)            ! y coordinate of cell center on polar stereographic projection
  real :: ymn(nypmax,maxgrds)            ! y coordinate of higher cell boundary on polar stereographic projection
  integer :: njnest(maxgrds)             ! index on next coarser grid where this grid starts (lower southwest corner)
  !                                      ! ds to interpolate y direction from coarser to finner grids: 
  integer :: jpm(nypmax,maxgrds)         ! next coarser grid cell index (jcoarser) that contains this finer grid cell
  real :: ej1(nypmax,maxgrds)            ! for jcoarser-1 on 3 points interpolation
  real :: ej2(nypmax,maxgrds)            ! for jcoarser   on 3 points interpolation
  real :: ej3(nypmax,maxgrds)            ! for jcoarser+1 on 3 points interpolation
  real :: ej4(nypmax,maxgrds)            ! for jcoarser-2 on 4 points interpolation
  real :: ej5(nypmax,maxgrds)            ! for jcoarser-1 on 4 points interpolation
  real :: ej6(nypmax,maxgrds)            ! for jcoarser   on 4 points interpolation
  real :: ej7(nypmax,maxgrds)            ! for jcoarser+1 on 4 points interpolation

  integer, target :: nnzp(maxgrds)       ! grid points z direction
  integer :: nnz(maxgrds)                ! nnzp - 1
  integer :: nnz1(maxgrds)               ! nnzp - 2
  real :: deltazn(maxgrds)               ! delta z
  real :: ztn(nzpmax,maxgrds)            ! z coordinate of interval center
  real :: zmn(nzpmax,maxgrds)            ! z coordinate of grid point
  integer :: nknest(maxgrds)             ! index on next coarser grid where this grid starts (lower level)
  !                                      ! ds to interpolate z direction (kcoarser) from coarser to finner grids: 
  integer :: kpm(nzpmax,maxgrds)         ! next coarser grid cell index that contains this finer grid cell
  real :: ek1(nzpmax,maxgrds)            ! for kcoarser-1 on 3 points interpolation
  real :: ek2(nzpmax,maxgrds)            ! for kcoarser   on 3 points interpolation
  real :: ek3(nzpmax,maxgrds)            ! for kcoarser+1 on 3 points interpolation
  real :: ek4(nzpmax,maxgrds)            ! for kcoarser-2 on 4 points interpolation
  real :: ek5(nzpmax,maxgrds)            ! for kcoarser-1 on 4 points interpolation
  real :: ek6(nzpmax,maxgrds)            ! for kcoarser   on 4 points interpolation
  real :: ek7(nzpmax,maxgrds)            ! for kcoarser+1 on 4 points interpolation

  integer :: nnxyp(maxgrds)              ! nnxp*nnyp (grid points at each vertical)
  integer :: nnxyzp(maxgrds)             ! nnxp*nnyp*nnzp (grid points at the air)
  integer :: nnxysp(maxgrds)             ! nnxp*nnyp*(nzg+nzs+3)*npatch (grid points beneath ground)

  real :: platn(maxgrds)                 ! pole latitude (degrees)
  real :: plonn(maxgrds)                 ! pole longitude (degrees)

  real :: centlat(maxgrds)               ! grid center latitude (degrees)
  real :: centlon(maxgrds)               ! grid center longitude (degrees)

  real :: htn(nzpmax,maxgrds)
  real :: ht2n(nzpmax,maxgrds)
  real :: ht4n(nzpmax,maxgrds)
  real :: hwn(nzpmax,maxgrds)
  real :: hw2n(nzpmax,maxgrds)
  real :: hw4n(nzpmax,maxgrds)
  real :: dztn(nzpmax,maxgrds)
  real :: dzmn(nzpmax,maxgrds)
  real :: dzt2n(nzpmax,maxgrds)
  real :: dzm2n(nzpmax,maxgrds)
  integer :: nxp
  integer :: nx
  integer :: nx1
  integer :: nx2
  integer :: nyp
  integer :: ny
  integer :: ny1
  integer :: ny2
  integer :: nzp
  integer :: nzpp
  integer :: nz
  integer :: nz1
  integer :: nxyzp
  integer :: nxyp
  integer :: nxysp
  integer :: nscl
  integer :: nsttop
  integer :: nstbot
  integer :: ndtrat

  integer :: jdim                        ! all horizontal grids are 1D (jdim=0) or 2D (jdim=1)
  real :: deltax
  real :: deltay
  real :: deltaz
  real :: ht(nzpmax)
  real :: ht2(nzpmax)
  real :: ht4(nzpmax)
  real :: hw(nzpmax)
  real :: hw2(nzpmax)
  real :: hw4(nzpmax)
  real :: zt(nzpmax)
  real :: zm(nzpmax)
  real :: dzt(nzpmax)
  real :: dzm(nzpmax)
  real :: dzt2(nzpmax)
  real :: dzm2(nzpmax)
  real :: xt(nxpmax)
  real :: xm(nxpmax)
  real :: yt(nypmax)
  real :: ym(nypmax)


  integer :: ngrid                       ! current grid
  integer :: nzg                         ! soil layers
  integer :: nzs                         ! snow layers
  integer :: npatch                      ! surface patches per grid cell
  integer :: if_adap
  integer :: itopo
  integer :: ihtran
  integer :: ngridc
  integer :: ngrido
  integer :: iscr1
  integer :: iscr2
  integer :: memsize
  integer :: iounit
  integer :: maxpro
  integer :: memscr
  integer :: memind
  integer :: iogrid
  integer :: maxpts
  integer :: maxnzp
  integer :: maxnxp
  integer :: maxnyp
  integer :: i2dvar
  integer :: memgrd(maxgrds)
  real :: ztop
  real :: dzrat
  real :: dzmax
  real :: eps
  integer :: impl
  integer :: ideltat

  integer :: iyeara
  integer :: imontha
  integer :: idatea
  integer :: ihoura
  integer :: itimea

  integer :: iyearz
  integer :: imonthz
  integer :: idatez
  integer :: ihourz
  integer :: itimez

  integer :: nacoust
  integer :: initial
  integer :: iflag
  integer :: nnacoust(maxgrds)
  real :: dimove(maxgrds)
  real :: djmove(maxgrds)
  real :: gridu(maxgrds)
  real :: gridv(maxgrds)
  real :: zz(nzpmax)
  real :: dtlong
  real :: sspct
  real :: polelat
  real :: polelon
  real :: cflxy(maxgrds)
  real :: cflz(maxgrds)
  character(len=16) :: runtype
  character(len=1)  :: timeunit
  integer :: isstp
  integer :: istp
  real    :: dts
  real    :: dtlt
  real    :: dtlv
  integer :: nestz1
  integer :: nestz2
  integer :: nndtrat(maxgrds)            ! delta t ratio (coarser/nested), indexed by nested)
  integer :: ngbegun(maxgrds)
  integer :: nnsttop(maxgrds)
  integer :: nnstbot(maxgrds)
  integer :: nstratz1(nzpmax)
  integer :: nstratz2(nzpmax)


  integer, parameter :: maxsched=2000    ! maximum number of nested timesteps for a dtlong time advance
  integer, parameter :: maxschent=5      ! number of events to be recorded at each nested timestep
  integer :: nsubs                       ! actual number of nested timesteps for a dtlong time advance
  integer :: isched(maxsched,maxschent)  ! nested timestep events (see modsched for detailed description)


  integer :: nrzflg
  integer :: nrz(nzpmax,maxgrds)
  real :: fbcf(nzpmax,maxgrds,4)
  integer :: iadvl
  integer :: iadvf
  integer :: lsflg
  integer :: ibnd
  integer :: jbnd
  integer :: icorflg
  integer :: nfpt
  integer :: naddsc
  integer :: iversion
  real :: distim
  real :: cphas

  ! Global simulation parameters

  integer :: nhemgrd2                    ! second hemispheric grid (0 if not global simulation)
  integer :: nhemt
  integer :: nhemu
  integer :: nhemv
  integer :: ihem1tt(4,maxhp)
  integer :: jhem1tt(4,maxhp)
  integer :: ihem1uu(4,maxhp)
  integer :: jhem1uu(4,maxhp)
  integer :: ihem1uv(4,maxhp)
  integer :: jhem1uv(4,maxhp)
  integer :: ihem1vu(4,maxhp)
  integer :: jhem1vu(4,maxhp)
  integer :: ihem1vv(4,maxhp)
  integer :: jhem1vv(4,maxhp)
  integer :: ihem2tt(maxhp)
  integer :: jhem2tt(maxhp)
  integer :: ihem2uu(maxhp)
  integer :: jhem2uu(maxhp)
  integer :: ihem2uv(maxhp)
  integer :: jhem2uv(maxhp)
  integer :: ihem2vu(maxhp)
  integer :: jhem2vu(maxhp)
  integer :: ihem2vv(maxhp)
  integer :: jhem2vv(maxhp)
  real :: hlatt(maxhp)
  real :: hlatu(maxhp)
  real :: hlatv(maxhp)
  real :: hlont(maxhp)
  real :: hlonu(maxhp)
  real :: hlonv(maxhp)
  real :: whem1tt(4,maxhp)
  real :: whem1uu(4,maxhp)
  real :: whem1uv(4,maxhp)
  real :: whem1vu(4,maxhp)
  real :: whem1vv(4,maxhp)

  ! ALF
  ! Flags to set the thermo call on the horizontal boundaries
  logical :: f_thermo_e(maxgrds)
  logical :: f_thermo_w(maxgrds)
  logical :: f_thermo_n(maxgrds)
  logical :: f_thermo_s(maxgrds)

contains





  subroutine alloc_grid(grid,n1,n2,n3,ng,if_adap)
    implicit none
    type (grid_vars) :: grid
    integer, intent(in) :: n1,n2,n3,ng,if_adap

    ! Allocate arrays based on options (if necessary)

    allocate (grid%topt(n2,n3))
    allocate (grid%topu(n2,n3))
    allocate (grid%topv(n2,n3))
    allocate (grid%topm(n2,n3))
    allocate (grid%topma(n2,n3))
    allocate (grid%topta(n2,n3))
    allocate (grid%rtgt(n2,n3))
    allocate (grid%rtgu(n2,n3))
    allocate (grid%rtgv(n2,n3))
    allocate (grid%rtgm(n2,n3))
    allocate (grid%f13t(n2,n3))
    allocate (grid%f13u(n2,n3))
    allocate (grid%f13v(n2,n3))
    allocate (grid%f13m(n2,n3))
    allocate (grid%f23t(n2,n3))
    allocate (grid%f23u(n2,n3))
    allocate (grid%f23v(n2,n3))
    allocate (grid%f23m(n2,n3))
    allocate (grid%dxt(n2,n3))
    allocate (grid%dxu(n2,n3))
    allocate (grid%dxv(n2,n3))
    allocate (grid%dxm(n2,n3))
    allocate (grid%dyt(n2,n3))
    allocate (grid%dyu(n2,n3))
    allocate (grid%dyv(n2,n3))
    allocate (grid%dym(n2,n3))
    allocate (grid%fmapt(n2,n3))
    allocate (grid%fmapu(n2,n3))
    allocate (grid%fmapv(n2,n3))
    allocate (grid%fmapm(n2,n3))
    allocate (grid%fmapti(n2,n3))
    allocate (grid%fmapui(n2,n3))
    allocate (grid%fmapvi(n2,n3))
    allocate (grid%fmapmi(n2,n3))
    allocate (grid%glat(n2,n3))
    allocate (grid%glon(n2,n3))
    allocate (grid%topzo(n2,n3))
    if (if_adap == 1) then
       allocate (grid%aru(n1,n2,n3)) 
       allocate (grid%arv(n1,n2,n3)) 
       allocate (grid%arw(n1,n2,n3))
       allocate (grid%volu(n1,n2,n3))
       allocate (grid%volv(n1,n2,n3))
       allocate (grid%volw(n1,n2,n3))
       allocate (grid%volt(n1,n2,n3))
    endif
    allocate (grid%flpu(n2,n3)) 
    allocate (grid%flpv(n2,n3)) 
    allocate (grid%flpw(n2,n3))
    allocate (grid%fmynum(n2,n3))
  end subroutine alloc_grid





  subroutine nullify_grid(grid)
    implicit none
    type (grid_vars) :: grid
    nullify (grid%topt)  ;  nullify (grid%topu)  ;  nullify (grid%topv)
    nullify (grid%topm)  ;  nullify (grid%topma) ;  nullify (grid%topta)   
    nullify (grid%rtgt)  ;  nullify (grid%rtgu)
    nullify (grid%rtgv)  ;  nullify (grid%rtgm)  ;  nullify (grid%f13t)
    nullify (grid%f13u)  ;  nullify (grid%f13v)  ;  nullify (grid%f13m)
    nullify (grid%f23t)  ;  nullify (grid%f23u)  ;  nullify (grid%f23v)
    nullify (grid%f23m)  ;  nullify (grid%dxt)   ;  nullify (grid%dxu)
    nullify (grid%dxv)   ;  nullify (grid%dxm)   ;  nullify (grid%dyt)
    nullify (grid%dyu)   ;  nullify (grid%dyv)   ;  nullify (grid%dym)
    nullify (grid%fmapt) ;  nullify (grid%fmapu) ;  nullify (grid%fmapv)
    nullify (grid%fmapm) ;  nullify (grid%fmapti);  nullify (grid%fmapui)
    nullify (grid%fmapvi);  nullify (grid%fmapmi);  nullify (grid%glat)
    nullify (grid%glon)  ;  nullify (grid%topzo)
    nullify (grid%aru)   ;  nullify (grid%arv)   ;  nullify (grid%arw)
    nullify (grid%volu)  ;  nullify (grid%volv)  ;  nullify (grid%volw)
    nullify (grid%volt)  
    nullify (grid%flpu)  ;  nullify (grid%flpv)  ;  nullify (grid%flpw)
    nullify (grid%fmynum)
  end subroutine nullify_grid




  subroutine dealloc_grid(grid)
    implicit none
    type (grid_vars) :: grid
    if (associated(grid%topt)  )    deallocate (grid%topt)
    if (associated(grid%topu)  )    deallocate (grid%topu)
    if (associated(grid%topv)  )    deallocate (grid%topv)
    if (associated(grid%topm)  )    deallocate (grid%topm)
    if (associated(grid%topma) )    deallocate (grid%topma)
    if (associated(grid%topta) )    deallocate (grid%topta)
    if (associated(grid%rtgt)  )    deallocate (grid%rtgt)
    if (associated(grid%rtgu)  )    deallocate (grid%rtgu)
    if (associated(grid%rtgv)  )    deallocate (grid%rtgv)
    if (associated(grid%rtgm)  )    deallocate (grid%rtgm)
    if (associated(grid%f13t)  )    deallocate (grid%f13t)
    if (associated(grid%f13u)  )    deallocate (grid%f13u)
    if (associated(grid%f13v)  )    deallocate (grid%f13v)
    if (associated(grid%f13m)  )    deallocate (grid%f13m)
    if (associated(grid%f23t)  )    deallocate (grid%f23t)
    if (associated(grid%f23u)  )    deallocate (grid%f23u)
    if (associated(grid%f23v)  )    deallocate (grid%f23v)
    if (associated(grid%f23m)  )    deallocate (grid%f23m)
    if (associated(grid%dxt)   )    deallocate (grid%dxt)
    if (associated(grid%dxu)   )    deallocate (grid%dxu)
    if (associated(grid%dxv)   )    deallocate (grid%dxv)
    if (associated(grid%dxm)   )    deallocate (grid%dxm)
    if (associated(grid%dyt)   )    deallocate (grid%dyt)
    if (associated(grid%dyu)   )    deallocate (grid%dyu)
    if (associated(grid%dyv)   )    deallocate (grid%dyv)
    if (associated(grid%dym)   )    deallocate (grid%dym)
    if (associated(grid%fmapt) )    deallocate (grid%fmapt)
    if (associated(grid%fmapu) )    deallocate (grid%fmapu)
    if (associated(grid%fmapv) )    deallocate (grid%fmapv)
    if (associated(grid%fmapm) )    deallocate (grid%fmapm)
    if (associated(grid%fmapti))    deallocate (grid%fmapti)
    if (associated(grid%fmapui))    deallocate (grid%fmapui)
    if (associated(grid%fmapvi))    deallocate (grid%fmapvi)
    if (associated(grid%fmapmi))    deallocate (grid%fmapmi)
    if (associated(grid%glat)  )    deallocate (grid%glat)
    if (associated(grid%glon)  )    deallocate (grid%glon)
    if (associated(grid%topzo) )    deallocate (grid%topzo)
    if (associated(grid%aru)   )    deallocate (grid%aru) 
    if (associated(grid%arv)   )    deallocate (grid%arv) 
    if (associated(grid%arw)   )    deallocate (grid%arw)
    if (associated(grid%volu)  )    deallocate (grid%volu)
    if (associated(grid%volv)  )    deallocate (grid%volv)
    if (associated(grid%volw)  )    deallocate (grid%volw)
    if (associated(grid%volt)  )    deallocate (grid%volt)
    if (associated(grid%flpu)  )    deallocate (grid%flpu) 
    if (associated(grid%flpv)  )    deallocate (grid%flpv) 
    if (associated(grid%flpw)  )    deallocate (grid%flpw)
    if (associated(grid%fmynum))    deallocate (grid%fmynum)
  end subroutine dealloc_grid




  subroutine filltab_grid(grid,gridm,imean,n1,n2,n3,ng)
    use var_tables
    implicit none
    type (grid_vars) :: grid,gridm
    integer, intent(in) :: imean,n1,n2,n3,ng
    integer :: npts
    real, pointer :: var,varm

    ! Fill pointers to arrays into variable tables

    npts=n2*n3
    if (associated(grid%topt)) &
         call vtables2 (grid%topt,gridm%topt,ng,npts,imean,  &
         'TOPT :2:hist:anal:mpti')      
    if (associated(grid%topu)) &
         call vtables2 (grid%topu,gridm%topu,ng, npts, imean,  &
         'TOPU :2:mpti')      
    if (associated(grid%topv)) &
         call vtables2 (grid%topv,gridm%topv,ng, npts, imean,  &
         'TOPV :2:mpti')      
    if (associated(grid%topm)) &
         call vtables2 (grid%topm,gridm%topm,ng, npts, imean,  &
         'TOPM :2:mpti')      
    if (associated(grid%topma)) &
         call vtables2 (grid%topma,gridm%topma,ng, npts, imean,  &
         'TOPMA :2:hist:anal:mpti')      
    if (associated(grid%topta)) &
         call vtables2 (grid%topta,gridm%topta,ng, npts, imean,  &
         'TOPTA :2:hist:anal:mpti')      
    if (associated(grid%rtgt)) &
         call vtables2 (grid%rtgt,gridm%rtgt,ng, npts, imean,  &
         'RTGT :2:mpti')      
    if (associated(grid%rtgu)) &
         call vtables2 (grid%rtgu,gridm%rtgu,ng, npts, imean,  &
         'RTGU :2:mpti')      
    if (associated(grid%rtgv)) &
         call vtables2 (grid%rtgv,gridm%rtgv,ng, npts, imean,  &
         'RTGV :2:mpti')      
    if (associated(grid%rtgm)) &
         call vtables2 (grid%rtgm,gridm%rtgm,ng, npts, imean,  &
         'RTGM :2:mpti')      
    if (associated(grid%f13t)) &
         call vtables2 (grid%f13t,gridm%f13t,ng, npts, imean,  &
         'F13T :2:mpti')      
    if (associated(grid%f13u)) &
         call vtables2 (grid%f13u,gridm%f13u,ng, npts, imean,  &
         'F13U :2:mpti')      
    if (associated(grid%f13v)) &
         call vtables2 (grid%f13v,gridm%f13v,ng, npts, imean,  &
         'F13V :2:mpti')      
    if (associated(grid%f13m)) &
         call vtables2 (grid%f13m,gridm%f13m,ng, npts, imean,  &
         'F13M :2:mpti')      
    if (associated(grid%f23t)) &
         call vtables2 (grid%f23t,gridm%f23t,ng, npts, imean,  &
         'F23T :2:mpti')      
    if (associated(grid%f23u)) &
         call vtables2 (grid%f23u,gridm%f23u,ng, npts, imean,  &
         'F23U :2:mpti')      
    if (associated(grid%f23v)) &
         call vtables2 (grid%f23v,gridm%f23v,ng, npts, imean,  &
         'F23V :2:mpti')      
    if (associated(grid%f23m)) &
         call vtables2 (grid%f23m,gridm%f23m,ng, npts, imean,  &
         'F23M :2:mpti')      
    if (associated(grid%dxt)) &
         call vtables2 (grid%dxt,gridm%dxt,ng, npts, imean,  &
         'DXT :2:mpti')      
    if (associated(grid%dxu)) &
         call vtables2 (grid%dxu,gridm%dxu,ng, npts, imean,  &
         'DXU :2:mpti')      
    if (associated(grid%dxv)) &
         call vtables2 (grid%dxv,gridm%dxv,ng, npts, imean,  &
         'DXV :2:mpti')      
    if (associated(grid%dxm)) &
         call vtables2 (grid%dxm,gridm%dxm,ng, npts, imean,  &
         'DXM :2:mpti')      
    if (associated(grid%dyt)) &
         call vtables2 (grid%dyt,gridm%dyt,ng, npts, imean,  &
         'DYT :2:mpti')      
    if (associated(grid%dyu)) &
         call vtables2 (grid%dyu,gridm%dyu,ng, npts, imean,  &
         'DYU :2:mpti')      
    if (associated(grid%dyv)) &
         call vtables2 (grid%dyv,gridm%dyv,ng, npts, imean,  &
         'DYV :2:mpti')      
    if (associated(grid%dym)) &
         call vtables2 (grid%dym,gridm%dym,ng, npts, imean,  &
         'DYM :2:mpti')      
    if (associated(grid%fmapt)) &
         call vtables2 (grid%fmapt,gridm%fmapt,ng, npts, imean,  &
         'FMAPT :2:mpti')      
    if (associated(grid%fmapu)) &
         call vtables2 (grid%fmapu,gridm%fmapu,ng, npts, imean,  &
         'FMAPU :2:mpti')      
    if (associated(grid%fmapv)) &
         call vtables2 (grid%fmapv,gridm%fmapv,ng, npts, imean,  &
         'FMAPV :2:mpti')      
    if (associated(grid%fmapm)) &
         call vtables2 (grid%fmapm,gridm%fmapm,ng, npts, imean,  &
         'FMAPM :2:mpti')      
    if (associated(grid%fmapti)) &
         call vtables2 (grid%fmapti,gridm%fmapti,ng, npts, imean,  &
         'FMAPTI :2:mpti')      
    if (associated(grid%fmapui)) &
         call vtables2 (grid%fmapui,gridm%fmapui,ng, npts, imean,  &
         'FMAPUI :2:mpti')      
    if (associated(grid%fmapvi)) &
         call vtables2 (grid%fmapvi,gridm%fmapvi,ng, npts, imean,  &
         'FMAPVI :2:mpti')      
    if (associated(grid%fmapmi)) &
         call vtables2 (grid%fmapmi,gridm%fmapmi,ng, npts, imean,  &
         'FMAPMI :2:mpti')      
    if (associated(grid%glat)) &
         call vtables2 (grid%glat,gridm%glat,ng, npts, imean,  &
         'GLAT :2:mpti:anal')      
    if (associated(grid%glon)) &
         call vtables2 (grid%glon,gridm%glon,ng, npts, imean,  &
         'GLON :2:mpti:anal')      
    if (associated(grid%topzo)) &
         call vtables2 (grid%topzo,gridm%topzo,ng, npts, imean,  &
         'TOPZO :2:mpti')      

    npts=n2*n3
    if (associated(grid%flpu)) &
         call vtables2 (grid%flpu,gridm%flpu,ng,npts,imean,  &
         'LPU :2:mpti')      
    if (associated(grid%flpv)) &
         call vtables2 (grid%flpv,gridm%flpv,ng,npts,imean,  &
         'LPV :2:mpti')      
    if (associated(grid%flpw)) &
         call vtables2 (grid%flpw,gridm%flpw,ng,npts,imean,  &
         'LPW :2:mpti')

    if (associated(grid%fmynum)) &
         call vtables2 (grid%fmynum,gridm%fmynum,ng,npts,imean,  &
         'MYNUM :2:mpti:anal')      

    npts=n1*n2*n3
    if (associated(grid%aru)) &
         call vtables2 (grid%aru,gridm%aru,ng,npts,imean,  &
         'ARU :3:mpti')      
    if (associated(grid%arv)) &
         call vtables2 (grid%arv,gridm%arv,ng,npts,imean,  &
         'ARV :3:mpti')      
    if (associated(grid%arw)) &
         call vtables2 (grid%arw,gridm%arw,ng,npts,imean,  &
         'ARW :3:mpti')      

    if (associated(grid%volu)) &
         call vtables2 (grid%volu,gridm%volu,ng,npts,imean,  &
         'VOLU :3:mpti')
    if (associated(grid%volv)) &
         call vtables2 (grid%volv,gridm%volv,ng,npts,imean,  &
         'VOLV :3:mpti')
    if (associated(grid%volw)) &
         call vtables2 (grid%volw,gridm%volw,ng,npts,imean,  &
         'VOLW :3:mpti')
    if (associated(grid%volt)) &
         call vtables2 (grid%volt,gridm%volt,ng,npts,imean,  &
         'VOLT :3:anal:mpti')
  end subroutine filltab_grid






  subroutine dump_mem_grid()

    ! dump_mem_grid: dumps mem_grid sizes and mappings for all grids

    implicit none
    integer :: ind, indGrid, j
    character(len=10) :: c0, c1
    character(len=*), parameter :: h="**(dump_mem_grid)**"

    do indGrid = 1, ngrids
       if (nxtnest(indGrid) == 0) then
          write(c0,"(i10)") indGrid
          write(*,"(a)") h//" Dumping outer grid number "//trim(adjustl(c0))
          write(c0,"(i10)") nnxp(indGrid)-2
          write(c1,"(f10.1)") deltaxn(indGrid)
          write(*,"(a)") h//" x axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" x cells: (index, higher cell boundary coordinate)"
          do ind = 1, nnxp(indGrid), 5
             do j = ind, min(ind+4,nnxp(indGrid))
                write(*,'(" (",i3,",",f10.1,");")',ADVANCE="NO") &
                     j,xmn(j,indGrid)
             end do
             write(*,"(1x)")
          end do
          write(c0,"(i10)") nnyp(indGrid)-2
          write(c1,"(f10.1)") deltayn(indGrid)
          write(*,"(a)") h//" y axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" y cells: (index, higher cell boundary coordinate)"
          do ind = 1, nnyp(indGrid), 5
             do j = ind, min(ind+4,nnyp(indGrid))
                write(*,'(" (",i3,",",f10.1,");")',ADVANCE="NO") &
                     j,ymn(j,indGrid)
             end do
             write(*,"(1x)")
          end do
       else
          write(*,"(a)") h
          write(c0,"(i10)") indGrid
          write(c1,"(i10)") nxtnest(indGrid)
          write(*,"(a)") h//" Dumping inner grid number "//trim(adjustl(c0))//&
               ", nested into grid number "//trim(adjustl(c1))
          write(c0,"(i10)") nnxp(indGrid)-2
          write(c1,"(f10.1)") deltaxn(indGrid)
          write(*,"(a)") h//" x axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" x cells: (index, higher cell boundary coordinate, coarser grid cell index)"
          do ind = 1, nnxp(indGrid), 5
             do j = ind, min(ind+4,nnxp(indGrid))
                write(*,'(" (",i3,",",f10.1,",",i3,");")',ADVANCE="NO") &
                     j,xmn(j,indGrid),ipm(j,indGrid)
             end do
             write(*,"(1x)")
          end do
          write(c0,"(i10)") nnyp(indGrid)-2
          write(c1,"(f10.1)") deltayn(indGrid)
          write(*,"(a)") h//" y axis has "//trim(adjustl(c0))//&
               " inner cells and 2 boundary cells of length "//trim(adjustl(c1))
          write(*,"(a)") h//" y cells: (index, higher cell boundary coordinate, coarser grid cell index)"
          do ind = 1, nnyp(indGrid), 5
             do j = ind, min(ind+4,nnyp(indGrid))
                write(*,'(" (",i3,",",f10.1,",",i3,");")',ADVANCE="NO") &
                     j,ymn(j,indGrid),jpm(j,indGrid)
             end do
             write(*,"(1x)")
          end do
       end if
    end do
  end subroutine dump_mem_grid
end module mem_grid
