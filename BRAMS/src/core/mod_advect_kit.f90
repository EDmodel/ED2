module advect_kit

  use mod_GhostBlockPartition, only : &
       GhostBlockPartition              ! intent(in) - Type
  use mod_GhostBlock, only : &
       GhostBlockPointer                ! intent(in) - Type

  implicit none

  private
  public :: advect_first_alloc
  public :: advect_dealloc
  public :: prepare_inv
  public :: calc_advec

  type GhostBlockPointer_ngrids
     type(GhostBlockPartition), pointer :: gbp
     type(GhostBlockPointer), pointer :: agb(:)
  end type GhostBlockPointer_ngrids
  type(GhostBlockPointer_ngrids), allocatable :: Magb(:)

  integer, allocatable :: tamBlk(:)
  integer :: nproc
  
contains

  subroutine advect_first_alloc(ngrids, mzp, mxp, myp)

    ! Modules for Block Partition
    use mod_GhostBlockPartition, only : &
         CreateGhostBlockPartition        ! Function
    use mod_GhostBlock, only : &
         CreateGhostBlock                 ! Subroutine

    ! Modules for data extration
    use mem_grid,    only: &
         if_adap              ! intent(in)

    use var_tables,  only: &
         num_scalar           ! intent(in)

    use node_mod,    only: &
         ibcon, ia, iz, izu, ja, jz, jzv ! intent(inout)

    implicit none
    
    ! Arguments:
    integer,          intent(in   ) :: ngrids
    integer,          intent(in   ) :: mzp(ngrids)
    integer,          intent(in   ) :: mxp(ngrids)
    integer,          intent(in   ) :: myp(ngrids)
    
    ! Local variables:
    integer :: iBlk, local_nBlk !,tamBlk(ngrids)
    integer :: ng
    integer :: nproc
    integer :: totalVector, vectorPerThread, left, factor
    
!!    !OMP
!!$    INTEGER :: NTHREADS, TID
!!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
    
    if (if_adap ==0) then
       ! Sigma-Z coordenate

       allocate(tamBlk(ngrids))
       
!!       !$OMP PARALLEL PRIVATE(NTHREADS, TID)
       
!!$     TID = OMP_GET_THREAD_NUM()
!!$     PRINT *, 'Hello World from thread = ', TID
!!$     !Only master thread does this
!!$     IF (TID .EQ. 0) THEN 
!!$        NTHREADS = OMP_GET_NUM_THREADS()
!!$        PRINT *, 'Number of threads = ', NTHREADS
!!$     END IF
!!$
!!$     nproc = NTHREADS
       nproc = 1
       
       ! Allocating Magb
       allocate(Magb(ngrids))

       do ng = 1, ngrids
          
          ! Defining block lenght
          ! original:
          !tamBlk(ng) = (mxp(ng)*(myp(ng)-2)-2)/nproc
          ! Adjusting the block lenght as close to 2000 as possible
          totalVector = mxp(ng)*(myp(ng)-2)-2
          vectorPerThread = totalVector/nproc
          left = totalVector - vectorPerThread*nproc
          if (left /= 0) then
             vectorPerThread = vectorPerThread + 1
          end if
          factor = 1
          do
             if (vectorPerThread <= 2000) then
                exit
             else
                factor = factor + 1
                vectorPerThread = totalVector/(factor*nproc)
                left = totalVector - vectorPerThread*factor*nproc
                if (left /= 0) then
                   vectorPerThread = vectorPerThread + 1
                end if
             end if
          end do
          tamBlk(ng) = vectorPerThread


          Magb(ng)%gbp => &
               CreateGhostBlockPartition(mxp(ng), myp(ng), tamBlk(ng))
          
          local_nBlk = Magb(ng)%gbp%nBlk

          allocate(Magb(ng)%agb(local_nBlk))

          call newgrid(ng)
          call node_index()
          
          do iBlk = 1, local_nBlk
        
             Magb(ng)%agb(iBlk)%p  =>  &
                  CreateGhostBlock(Magb(ng)%gbp, mzp(ng), num_scalar(ng), &
                  iBlk, ibcon, ia, iz, izu, ja, jz, jzv)
        
          end do
       end do
    endif

  end subroutine advect_first_alloc

  ! ------------------------------------------------------------------------

  subroutine advect_dealloc(ngrids)

  use mod_GhostBlockPartition, only : &
       DestroyGhostBlockPartition       ! Subroutine
  use mod_GhostBlock, only : &
       DestroyGhostBlock       ! Subroutine

    implicit none
    ! Arguments:
    integer,          intent(in   ) :: ngrids
    ! Local Variables:
    integer :: ng, iBlk

    ! desaloca
    do ng=1, ngrids
       do iBlk = 1, Magb(ng)%gbp%nBlk
          call DestroyGhostBlock(Magb(ng)%agb(iBlk)%p)
       end do
       deallocate(Magb(ng)%agb)
       call DestroyGhostBlockPartition(Magb(ng)%gbp)
    enddo
    !deallocate (gbp)
    deallocate (Magb)
     
  end subroutine advect_dealloc

  ! ------------------------------------------------------------------------

  subroutine prepare_inv(ngrids)

    use mod_GhostBlock, only : &
         PutTimeInvGhostBlock    ! Subroutine

        use mem_grid,    only: &
         if_adap,          &  ! intent(in)
         grid_g,           &  ! %dxt,%dxu,%dxv,%dyt,%dyu,%rtgt,%rtgu,%rtgv,
         ! %f13t,%f23t,%fmapt,%fmapu,%fmapv,%fmapui,%fmapvi,
         ! %rtgt,%rtgu,%rtgv,%fmapt,%fmapui,%fmapvi,%f13t,
         ! %f23t,%dxu,%dyv,%dxt,%dyt, all intent(in)
         dtlt,             &  ! intent(in)
         jdim,             &  ! intent(in)
         time,             &
         zt,               &  ! intent(in)
         zm,               &  ! intent(in)
         dzm,              &  ! intent(in)
         hw4,              &  ! intent(in)
         dzt                  ! intent(in)

    use mem_basic,   only: &
         basic_g              ! %dn0,%dn0u,%dn0v, all intent(in)
   
    implicit none
    
    ! Arguments:
    integer,          intent(in   ) :: ngrids

    ! Local variables:
    integer :: ng, iBlk, local_nBlk

    if (if_adap ==0) then
       ! Sigma-Z coordenate
       do ng = 1, ngrids
          local_nBlk = Magb(ng)%gbp%nBlk
          
          do iBlk = 1, local_nBlk

             call newgrid(ng)
        
             call PutTimeInvGhostBlock(grid_g(ng)%dxt, grid_g(ng)%dxu, &
                  grid_g(ng)%dxv, grid_g(ng)%dyt, grid_g(ng)%dyu, &
                  grid_g(ng)%dyv, grid_g(ng)%f13t, grid_g(ng)%f23t, hw4, &
                  grid_g(ng)%fmapt, grid_g(ng)%fmapu, grid_g(ng)%fmapv, &
                  grid_g(ng)%fmapui, grid_g(ng)%fmapvi, &
                  grid_g(ng)%rtgt, grid_g(ng)%rtgu, grid_g(ng)%rtgv, &
                  dzt, dzm, zt, zm, &
                  basic_g(ng)%dn0, basic_g(ng)%dn0u, basic_g(ng)%dn0v, &
                  Magb(ng)%agb(iBlk)%p)
        
          end do
       enddo
    endif

  end subroutine prepare_inv

  ! ------------------------------------------------------------------------

  subroutine calc_advec(varn, ngrid, mzp, mxp, myp)

    ! Modules for Block Partition
    use mod_GhostBlock, only : &
         GhostBlock ,          & ! intent(in) - Type
         PutTimeVarGhostBlock, & ! Subroutine
         PutTimeVarVGhostBlock, & ! Subroutine
         PutTimeVarSGhostBlock, & ! Subroutine
         GetTimeVarGhostBlock_oper, &  ! Subroutine
         GetTimeVarVGhostBlock_oper, &  ! Subroutine
         GetTimeVarSGhostBlock_oper   ! Subroutine

    ! Modules for data extration
    use mem_tend,    only: &
         tend_g              ! %ut (inout); %vt (inout); %wt (inout)
    use var_tables,  only: &
         num_scalar,       & ! intent(in)
         scalar_tab          ! %var_p (in); %var_t (inout)
    use mem_grid,    only: &
         if_adap,          &  ! intent(in)
         dtlt,             &  ! intent(in)
         itopo                ! intent(in)
    use mem_basic,   only: &
         basic_g              ! %uc,%vc,%wc,%up,%vp,%wp,%dn0,%dn0u,%dn0v, all intent(in)

    use node_mod, only : &
         ia, iz, ja, jz, izu, jzv ! intent(inout)

    implicit none
    ! Arguments:
    character(len=*), intent(in   ) :: varn
    integer,          intent(in   ) :: ngrid
    integer,          intent(in   ) :: mzp
    integer,          intent(in   ) :: mxp
    integer,          intent(in   ) :: myp

    ! Local variables:
    integer :: local_nBlk, iBlk, n, ijk, i, j, k
    type(GhostBlock), pointer :: p
    ! Local array for Scalars
    real    :: scalar_p(mzp, mxp, myp, num_scalar(ngrid))
    real    :: scalar_t(mzp, mxp, myp, num_scalar(ngrid))
    real    :: local_ut(mzp, mxp, myp)
    real    :: local_vt(mzp, mxp, myp)
    real    :: local_wt(mzp, mxp, myp)

    if (if_adap ==0) then
       ! SigmaZ coordenate: new optmized process

       call newgrid(ngrid)

       !!call node_index()

       if (varn .eq. 'T' .or. varn .eq. 'ALL') then
          ! Initiating scalar's local arrays

          !!scalar_p = 0.
          !!scalar_t = 0.

          do n = 1, num_scalar(ngrid)
             ijk = 0
             do j = 1, myp
                do i = 1, mxp
                   do k = 1, mzp
                      ijk = ijk + 1
                      scalar_p(k, i, j, n) = scalar_tab(n, ngrid)%a_var_p(ijk)
                      scalar_t(k, i, j, n) = scalar_tab(n, ngrid)%a_var_t(ijk)
                   enddo
                enddo
             enddo
          enddo

       endif
       
       if (varn .eq. 'V' .or. varn .eq. 'ALL') then
          ! Initialising tendencies local arrays

          ijk = 0
          do j = 1, myp
             do i = 1, mxp
                do k = 1, mzp
                   ijk = ijk + 1
                   local_ut(k, i, j) = tend_g(ngrid)%ut(k,i,j)
                   local_vt(k, i, j) = tend_g(ngrid)%vt(k,i,j)
                   local_wt(k, i, j) = tend_g(ngrid)%wt(k,i,j)
                enddo
             enddo
          enddo
          !

       endif
       
       local_nBlk = Magb(ngrid)%gbp%nBlk
          
       do iBlk = 1, local_nBlk
        
          if (varn .eq. 'V' .or. varn .eq. 'ALL') then
             ! Prepare to Advect  U  V, and W

             call PutTimeVarVGhostBlock(basic_g(ngrid)%up, basic_g(ngrid)%vp, &
                  basic_g(ngrid)%wp, basic_g(ngrid)%uc, basic_g(ngrid)%vc, &
                  basic_g(ngrid)%wc, local_ut, local_vt, local_wt, &
                  Magb(ngrid)%agb(iBlk)%p)

          endif

          if (varn .eq. 'T' .or. varn .eq. 'ALL') then
             ! Prepare to Advect  scalars

             call PutTimeVarSGhostBlock(basic_g(ngrid)%up, basic_g(ngrid)%vp, &
                  basic_g(ngrid)%wp, basic_g(ngrid)%uc, basic_g(ngrid)%vc, &
                  basic_g(ngrid)%wc, scalar_p, scalar_t, &
                  Magb(ngrid)%agb(iBlk)%p)

          endif


       enddo

!!       !$OMP PARALLEL PRIVATE(NTHREADS, TID)

!!       !$OMP DO PRIVATE(iBlk, p)
       do iBlk = 1, local_nBlk
        
          p => Magb(ngrid)%agb(iBlk)%p
        
          if (varn .eq. 'V' .or. varn .eq. 'ALL') then
             ! Advect  U  V, and W

             call advect_vec(p%part%iMax, p%ibDim,   &  
                  p%ibStart, p%ibEnd, p%ibEndU, p%ibEndV,     &
                  p%kMax, 2, p%kMax-1, p%kMax-2,          &
                  p%uc,p%vc,p%wc,p%ut,p%vt,p%wt,p%dn0,p%dn0u,p%dn0v,    &
                  p%dxt,p%dxu,p%dxv,p%dyt,p%dyu,p%dyv,            &
                  p%rtgt,p%rtgu,p%rtgv,p%f13t,p%f23t,           &
                  p%fmapt,p%fmapu,p%fmapv,p%fmapui,p%fmapvi,    &
                  itopo, p%hw4, p%dzt, p%dzm,               &
                  p%realData, p%noWorkU, p%noWorkV, &
                  ! DEBUG-ALF
                  p%iPerIb, p%jPerIb)
                  !
          endif
        
          if (varn .eq. 'T' .or. varn .eq. 'ALL') then
             ! Advect  scalars

             call advect_sca(p%ibDim, p%kMax, p%part%iMax, &
                  p%ibStart, p%ibEnd, 2, p%kMax-1, p%nScalars, &
                  dtlt, p%rtgt, p%rtgu, p%rtgv, p%dxu, p%dyv, p%dzm, &
                  p%zt, p%zm, p%hw4, &
                  p%dn0, p%dn0u, p%dn0v, p%f13t, p%f23t, p%dxt, p%dyt, p%dzt, &
                  p%fmapt, p%fmapui, p%fmapvi, &
                  p%up, p%vp, p%wp, p%uc, p%vc, p%wc, p%RealData, p%noGhostI, &
                  p%noBottonGhostI, p%scalarp, p%scalart, &
                  ! DEBUG-ALF
                  p%iPerIb, p%jPerIb)
                  !
          endif

       enddo
!!       !$OMP END DO

!!       !$OMP END PARALLEL

       do iBlk = 1, local_nBlk

          if (varn .eq. 'V' .or. varn .eq. 'ALL') then
             ! Get data from Advect  U  V, and W

             call GetTimeVarVGhostBlock_oper(local_ut, local_vt, local_wt, &
               Magb(ngrid)%agb(iBlk)%p)

          endif
          if (varn .eq. 'T' .or. varn .eq. 'ALL') then
             ! Get data from Advect  scalars

             call GetTimeVarSGhostBlock_oper(scalar_t, Magb(ngrid)%agb(iBlk)%p)

          endif

       enddo

       ! Updating tendencies arrays with local values

       if (varn .eq. 'V' .or. varn .eq. 'ALL') then
          ! Get data from Advect  U  V, and W

          ijk = 0
          do j = 1, myp
             do i = 1, mxp
                do k = 1, mzp
                   ijk = ijk + 1
                   ! ALF
                   if (i>=ia .and. i<=izu .and. j>=ja .and. j<=jz) &
                   tend_g(ngrid)%ut(k,i,j) = local_ut(k, i, j)
                   ! ALF
                   if (i>=ia .and. i<=iz .and. j>=ja .and. j<=jzv) &
                   tend_g(ngrid)%vt(k,i,j) = local_vt(k, i, j)
                   ! ALF
                   if (i>=ia .and. i<=iz .and. j>=ja .and. j<=jz) &
                   tend_g(ngrid)%wt(k,i,j) = local_wt(k, i, j)
                enddo
             enddo
          enddo

       endif

       if (varn .eq. 'T' .or. varn .eq. 'ALL') then
          ! Get data from Advect  scalars

          do n = 1, num_scalar(ngrid)
             ijk = 0
             do j = 1, myp
                do i = 1, mxp
                   do k = 1, mzp
                      ijk = ijk + 1
                      ! ALF-TESTE
                      if (i>=ia .and. i<=iz .and. j>=ja .and. j<=jz) &
                           scalar_tab(n, ngrid)%a_var_t(ijk) = &
                           scalar_t(k, i, j, n)
                   enddo
                enddo
             enddo
          enddo

       endif

    else
       
       ! Shaved ETA coordenate: call old subrotuine
       call ADVECTc(varn,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,1)

    endif

  end subroutine calc_advec

end module advect_kit
