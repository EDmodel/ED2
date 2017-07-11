module mod_GhostBlock
  use mod_GhostBlockPartition, only: GhostBlockPartition
  implicit none

  private
  public :: GhostBlock
  public :: GhostBlockPointer
  public :: CreateGhostBlock
  public :: DumpGhostBlock
  public :: DestroyGhostBlock
  public :: PutTimeInvGhostBlock
  public :: PutTimeVarGhostBlock
  public :: PutTimeVarVGhostBlock
  public :: PutTimeVarSGhostBlock
  public :: GetTimeVarGhostBlock
  public :: GetTimeVarGhostBlock_oper
  public :: GetTimeVarVGhostBlock_oper
  public :: GetTimeVarSGhostBlock_oper

  type GhostBlock
     type(GhostBlockPartition), pointer :: part
     integer          :: blkId
     integer          :: kMax
     integer          :: ibDim       ! dimensioned
     integer          :: ibStart     ! first index without ghost zone
     integer          :: ibEnd       ! last index without ghost zone
     integer          :: ibMax       ! last index with ghost zone
     integer          :: ibEndU      ! last index without ghost zone that contributes for u
     integer          :: ibEndV      ! last index without ghost zone that contributes for v
     integer          :: nScalars    ! number of scalars
     integer, pointer :: iPerIb(:)   ! map ib into i
     integer, pointer :: jPerIb(:)   ! map ib into j
     real,    pointer :: realData(:) ! 1 if outside ghost zone, 0 ow
     real,    pointer :: noWorkU(:)  ! 1 if computes u, 0 ow
     real,    pointer :: noWorkV(:)  ! 1 if computes v, 0 ow
     real,    pointer :: noGhostI(:) !
     real,    pointer :: noBottonGhostI(:) !
     real,    pointer :: dxt(:)
     real,    pointer :: dxu(:)
     real,    pointer :: dxv(:)
     real,    pointer :: dyt(:)
     real,    pointer :: dyu(:)
     real,    pointer :: dyv(:)
     real,    pointer :: f13t(:)
     real,    pointer :: f23t(:)
     real,    pointer :: hw4(:)
     real,    pointer :: fmapt(:)
     real,    pointer :: fmapu(:)
     real,    pointer :: fmapv(:)
     real,    pointer :: fmapui(:)
     real,    pointer :: fmapvi(:)
     real,    pointer :: rtgt(:)
     real,    pointer :: rtgu(:)
     real,    pointer :: rtgv(:)
     real,    pointer :: dzt(:)
     real,    pointer :: dzm(:)
     real,    pointer :: zt(:)
     real,    pointer :: zm(:)
     real,    pointer :: up(:,:)
     real,    pointer :: vp(:,:)
     real,    pointer :: wp(:,:)
     real,    pointer :: uc(:,:)
     real,    pointer :: vc(:,:)
     real,    pointer :: wc(:,:)
     real,    pointer :: ut(:,:)
     real,    pointer :: vt(:,:)
     real,    pointer :: wt(:,:)
     real,    pointer :: dn0(:,:)
     real,    pointer :: dn0u(:,:)
     real,    pointer :: dn0v(:,:)
     real,    pointer :: scalarp(:,:,:)
     real,    pointer :: scalart(:,:,:)
  end type GhostBlock


  type GhostBlockPointer
     type(GhostBlock), pointer :: p
  end type GhostBlockPointer

  integer, parameter :: CacheLineRealData=16
contains





  function CreateGhostBlock(part, kMax, nScalars, blkNbr, ibcon, &
       ia, iz, izu, ja, jz, jzv) result(gb)

    ! Arguments:
    type(GhostBlockPartition), pointer :: part
    integer, intent(in) :: kMax
    integer, intent(in) :: nScalars
    integer, intent(in) :: blkNbr
    integer, intent(in) :: ibcon
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: izu
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: jzv
    ! Local Data:
    type(GhostBlock), pointer :: gb
    character(len=*), parameter :: h="**(CreateGhostBlock)**"
    integer :: rem, ib, ij, i, j

    ! consistency

    if (.not. associated(part)) then
       write(*,"(a)") h//" GhostBlockPartition not associated"
       stop
    else if (  &
         .not. associated(part%blockSize) .or. &
         .not. associated(part%ijStart)   .or. &
         .not. associated(part%ijEnd)     ) then
       write(*,"(a)") h//" GhostBlockPartition not fully created"
       stop
    else
       if (    &
            (part%nBlk > size(part%blockSize)) .or. &
            (part%nBlk > size(part%ijStart))   .or. &
            (part%nBlk > size(part%ijEnd))    ) then
          write(*,"(a)") h//" inconsistent GhostBlockPartition "
          stop
       end if
    end if
    if (blkNbr <= 0 .or. blkNbr > part%nBlk) then
       write(*,"(a,i10)") h//" inconsistent blkNbr=",blkNbr
       stop
    else if (kMax <=0) then
       write(*,"(a,i10)") h//" inconsistent kMax=",kMax
       stop
    else if (nScalars <=0) then
       write(*,"(a,i10)") h//" inconsistent nScalars=",nScalars
       stop
    end if

    ! fill type variable

    allocate(gb)

    gb%part => part
    gb%blkId   = blkNbr
    gb%kMax    = kMax
    gb%ibStart = gb%part%iMax + 2
    gb%ibEnd   = gb%ibStart + gb%part%ijEnd(blkNbr) - gb%part%ijStart(blkNbr)
    gb%ibMax   = gb%ibEnd + gb%part%iMax + 1
    gb%nScalars= nScalars

    rem = mod(gb%ibMax,CacheLineRealData)
    if (rem == 0) then
       gb%ibDim = gb%ibMax
    else
       gb%ibDim = gb%ibMax + CacheLineRealData - rem
    end if

    ! Teste para SX-6
    ! gb%ibDim igual a numero impar
!!$    if (mod(gb%ibMax, 2) == 0) then
!!$       gb%ibDim = gb%ibMax + 1
!!$    else
!!$       gb%ibDim = gb%ibMax
!!$    endif


    allocate(gb%iPerIb(gb%ibDim))
    allocate(gb%jPerIb(gb%ibDim))
    allocate(gb%RealData(gb%ibDim))
    allocate(gb%noWorkU(gb%ibDim))
    allocate(gb%noWorkV(gb%ibDim))
    allocate(gb%noGhostI(gb%ibDim))
    allocate(gb%noBottonGhostI(gb%ibDim))

    ij = gb%part%ijStart(blkNbr)-gb%part%iMax-1
    do ib = 1, gb%ibMax
       gb%iPerIb(ib) = mod(ij,gb%part%iMax) + 1
       gb%jPerIb(ib) = ij/gb%part%iMax + 1
       ij = ij + 1
    end do

    gb%realData = 0.0
    gb%noWorkU  = 0.0
    gb%noWorkV  = 0.0
    gb%noGhostI = 0.0
    gb%noBottonGhostI = 0.0
    gb%ibEndU = 0
    gb%ibEndV = 0
    do ib = 1, gb%ibMax
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       if ( i>=2 .and. i<=gb%part%iMax-1 ) then
          gb%noGhostI(ib) = 1.0
       end if
       if ( i<=gb%part%iMax-1 ) then
          gb%noBottonGhostI(ib) = 1.0
       end if
    end do

    do ib = gb%ibStart, gb%ibEnd
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       if ( i>=2 .and. i<=gb%part%iMax-1 .and. &
            j>=2 .and. j<=gb%part%jMax-1 ) then
          gb%realData(ib) = 1.0
       end if
       if ( i>=ia .and. i<=izu .and. &
            j>=ja .and. j<=jz ) then
          gb%noWorkU(ib) = 1.0
          gb%ibEndU = ib
       end if
       if ( i>=ia .and. i<=iz .and. &
            j>=ja .and. j<=jzv ) then
          gb%noWorkV(ib) = 1.0
          gb%ibEndV = ib
       end if
    end do

    allocate(gb%dxt(gb%ibDim))
    allocate(gb%dxu(gb%ibDim))
    allocate(gb%dxv(gb%ibDim))
    allocate(gb%dyt(gb%ibDim))
    allocate(gb%dyu(gb%ibDim))
    allocate(gb%dyv(gb%ibDim))
    allocate(gb%f13t(gb%ibDim))
    allocate(gb%f23t(gb%ibDim))
    allocate(gb%hw4(gb%kMax))
    allocate(gb%fmapt(gb%ibDim))
    allocate(gb%fmapu(gb%ibDim))
    allocate(gb%fmapv(gb%ibDim))
    allocate(gb%fmapui(gb%ibDim))
    allocate(gb%fmapvi(gb%ibDim))
    allocate(gb%rtgt(gb%ibDim))
    allocate(gb%rtgu(gb%ibDim))
    allocate(gb%rtgv(gb%ibDim))
    allocate(gb%dzt(gb%kMax))
    allocate(gb%dzm(gb%kMax))
    allocate(gb%zt(gb%kMax))
    allocate(gb%zm(gb%kMax))
    allocate(gb%up(gb%ibDim,gb%kMax))
    allocate(gb%vp(gb%ibDim,gb%kMax))
    allocate(gb%wp(gb%ibDim,gb%kMax))
    allocate(gb%uc(gb%ibDim,gb%kMax))
    allocate(gb%vc(gb%ibDim,gb%kMax))
    allocate(gb%wc(gb%ibDim,gb%kMax))
    allocate(gb%ut(gb%ibDim,gb%kMax))
    allocate(gb%vt(gb%ibDim,gb%kMax))
    allocate(gb%wt(gb%ibDim,gb%kMax))
    allocate(gb%dn0(gb%ibDim,gb%kMax))
    allocate(gb%dn0u(gb%ibDim,gb%kMax))
    allocate(gb%dn0v(gb%ibDim,gb%kMax))
    allocate(gb%scalarp(gb%ibDim,gb%kMax,gb%nScalars))
    allocate(gb%scalart(gb%ibDim,gb%kMax,gb%nScalars))
  end function CreateGhostBlock






  subroutine DumpGhostBlock(p)
    type(GhostBlock), pointer :: p

    integer :: nRealData
    character(len=*), parameter :: h="**(DumpGhostBlock)**"
    character(len=20) :: c0, c1, c2, c3, c4
    integer :: i, j, ib

    if (associated(p)) then
       if ( associated(p%iPerIb)   .and. &
            associated(p%jPerIb)   .and. &
            associated(p%realData) .and. &
            associated(p%noWorkU)  .and. &
            associated(p%noWorkV)  ) then
          write(c0,"(i20)") p%ibDim
          write(c1,"(i20)") p%ibStart
          write(c2,"(i20)") p%ibEnd
          write(c3,"(i20)") p%ibMax
          write(c4,"(i20)") p%blkId
          write(*,"(a)") h//" Real data block # "//trim(adjustl(c4))//&
               &" starts at "//trim(adjustl(c1))//&
               &"; ends at "//trim(adjustl(c2))//"; ghost zone ends at "//&
               &trim(adjustl(c3))//"; dimensioned "//trim(adjustl(c0))
          nRealData = count(p%realData(1:p%ibMax) == 1.0)
          write(c0,"(i20)") (100*nRealData)/p%ibMax
          write(c1,"(i20)") nRealData
          write(c2,"(i20)") p%ibMax-nRealData
          write(*,"(a)") h//" "//trim(adjustl(c1))//" real data points ("//& 
               trim(adjustl(c0))//"%) and "//&
               trim(adjustl(c2))//" ghost zone points"
          write(*,"(a)") h//" ib, i, j, RealData, noWorkU, noWorkV, noGhostI, noBottonGhostI"
          do ib = 1, p%ibMax
             write(*,"(i8,2i4,5f3.0)") ib, p%iPerIb(ib), p%jPerIb(ib), p%RealData(ib), &
                  p%noWorkU(ib), p%noWorkV(ib), p%noGhostI(ib), p%noBottonGhostI(ib)
          end do
       else
          write(*,"(a)") h//" GhostBlock not fully created"
       end if
    else
       write(*,"(a)") h//" GhostBlock not associated"
    end if
  end subroutine DumpGhostBlock





  subroutine DestroyGhostBlock (gb)
    type(GhostBlock), pointer :: gb

    character(len=*), parameter :: h="**(DestroyGhostBlock)**"
    integer :: rem, ib, ij, i, j

    if (.not. associated(gb)) then
       write(*,"(a)") h//" GhostBlock not associated"
    else
       nullify(gb%part)              ! memory leak; part should be deallocated explicitelly elsewhere
       if (associated(gb%iPerIb)) deallocate(gb%iPerIb)
       if (associated(gb%jPerIb)) deallocate(gb%jPerIb)
       if (associated(gb%realData)) deallocate(gb%realData)
       if (associated(gb%noWorkU)) deallocate(gb%noWorkU)
       if (associated(gb%noWorkV)) deallocate(gb%noWorkV)
       if (associated(gb%noGhostI)) deallocate(gb%noGhostI)
       if (associated(gb%noBottonGhostI)) deallocate(gb%noBottonGhostI)
       if (associated(gb%dxt)) deallocate(gb%dxt)
       if (associated(gb%dxu)) deallocate(gb%dxu)
       if (associated(gb%dxv)) deallocate(gb%dxv)
       if (associated(gb%dyt)) deallocate(gb%dyt)
       if (associated(gb%dyu)) deallocate(gb%dyu)
       if (associated(gb%dyv)) deallocate(gb%dyv)
       if (associated(gb%f13t)) deallocate(gb%f13t)
       if (associated(gb%f23t)) deallocate(gb%f23t)
       if (associated(gb%hw4)) deallocate(gb%hw4)
       if (associated(gb%fmapt)) deallocate(gb%fmapt)
       if (associated(gb%fmapu)) deallocate(gb%fmapu)
       if (associated(gb%fmapv)) deallocate(gb%fmapv)
       if (associated(gb%fmapui)) deallocate(gb%fmapui)
       if (associated(gb%fmapvi)) deallocate(gb%fmapvi)
       if (associated(gb%rtgt)) deallocate(gb%rtgt)
       if (associated(gb%rtgu)) deallocate(gb%rtgu)
       if (associated(gb%rtgv)) deallocate(gb%rtgv)
       if (associated(gb%dzt)) deallocate(gb%dzt)
       if (associated(gb%dzm)) deallocate(gb%dzm)
       if (associated(gb%zt)) deallocate(gb%zt)
       if (associated(gb%zm)) deallocate(gb%zm)
       if (associated(gb%up)) deallocate(gb%up)
       if (associated(gb%vp)) deallocate(gb%vp)
       if (associated(gb%wp)) deallocate(gb%wp)
       if (associated(gb%uc)) deallocate(gb%uc)
       if (associated(gb%vc)) deallocate(gb%vc)
       if (associated(gb%wc)) deallocate(gb%wc)
       if (associated(gb%ut)) deallocate(gb%ut)
       if (associated(gb%vt)) deallocate(gb%vt)
       if (associated(gb%wt)) deallocate(gb%wt)
       if (associated(gb%dn0)) deallocate(gb%dn0)
       if (associated(gb%dn0u)) deallocate(gb%dn0u)
       if (associated(gb%dn0v)) deallocate(gb%dn0v)
       if (associated(gb%scalarp)) deallocate(gb%scalarp)
       if (associated(gb%scalart)) deallocate(gb%scalart)
       deallocate(gb)
       nullify(gb)
    end if
  end subroutine DestroyGhostBlock



  subroutine PutTimeInvGhostBlock(dxt, dxu, dxv, dyt, dyu, dyv, &
       f13t, f23t, hw4, fmapt, fmapu, fmapv, fmapui, fmapvi, &
       rtgt, rtgu, rtgv, dzt, dzm, zt, zm, dn0, dn0u, dn0v, gb)

    ! Arguments:
    real, intent(in) :: dxt(:,:)
    real, intent(in) :: dxu(:,:)
    real, intent(in) :: dxv(:,:)
    real, intent(in) :: dyt(:,:)
    real, intent(in) :: dyu(:,:)
    real, intent(in) :: dyv(:,:)
    real, intent(in) :: f13t(:,:)
    real, intent(in) :: f23t(:,:)
    real, intent(in) :: fmapt(:,:)
    real, intent(in) :: fmapu(:,:)
    real, intent(in) :: fmapv(:,:)
    real, intent(in) :: fmapui(:,:)
    real, intent(in) :: fmapvi(:,:)
    real, intent(in) :: rtgt(:,:)
    real, intent(in) :: rtgu(:,:)
    real, intent(in) :: rtgv(:,:)
    real, intent(in) :: hw4(:)
    real, intent(in) :: dzt(:)
    real, intent(in) :: dzm(:)
    real, intent(in) :: zt(:)
    real, intent(in) :: zm(:)
    real, intent(in) :: dn0(:,:,:)
    real, intent(in) :: dn0u(:,:,:)
    real, intent(in) :: dn0v(:,:,:)
    type(GhostBlock), pointer :: gb

    ! Local Variables:
    integer :: k, ib, i, j

    do ib = 1, gb%ibMax
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)

       gb%dxt(ib) = dxt(i,j)
       gb%dxu(ib) = dxu(i,j)
       gb%dxv(ib) = dxv(i,j)
       gb%dyt(ib) = dyt(i,j)
       gb%dyu(ib) = dyu(i,j)
       gb%dyv(ib) = dyv(i,j)
       gb%f13t(ib) = f13t(i,j)
       gb%f23t(ib) = f23t(i,j)
       gb%fmapt(ib) = fmapt(i,j)
       gb%fmapu(ib) = fmapu(i,j)
       gb%fmapv(ib) = fmapv(i,j)
       gb%fmapui(ib) = fmapui(i,j)
       gb%fmapvi(ib) = fmapvi(i,j)
       gb%rtgt(ib) = rtgt(i,j)
       gb%rtgu(ib) = rtgu(i,j)
       gb%rtgv(ib) = rtgv(i,j)

       do k = 1, gb%kMax
          gb%dn0(ib,k) = dn0(k,i,j)
          gb%dn0u(ib,k) = dn0u(k,i,j)
          gb%dn0v(ib,k) = dn0v(k,i,j)
       enddo
    end do
    do k = 1, gb%kMax
       gb%hw4(k) = hw4(k)
       gb%dzt(k) = dzt(k)
       gb%dzm(k) = dzm(k)
       gb%zt(k) = zt(k)
       gb%zm(k) = zm(k)
    end do

  end subroutine PutTimeInvGhostBlock






  subroutine PutTimeVarGhostBlock(up, vp, wp, uc, vc, wc, &
       ut, vt, wt, scalarp, scalart, gb)
    real,    intent(in) :: up(:,:,:)
    real,    intent(in) :: vp(:,:,:)
    real,    intent(in) :: wp(:,:,:)
    real,    intent(in) :: uc(:,:,:)
    real,    intent(in) :: vc(:,:,:)
    real,    intent(in) :: wc(:,:,:)
    real,    intent(in) :: ut(:,:,:)
    real,    intent(in) :: vt(:,:,:)
    real,    intent(in) :: wt(:,:,:)
    real,    intent(in) :: scalarp(:,:,:,:)
    real,    intent(in) :: scalart(:,:,:,:)
    type(GhostBlock), pointer :: gb

    integer :: ib, i, j, k, iSca, nScalars

    nScalars=size(scalarp,4)

    do ib = 1, gb%ibMax
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       do k = 1, gb%kMax
          gb%up(ib,k) = up(k,i,j)
          gb%vp(ib,k) = vp(k,i,j)
          gb%wp(ib,k) = wp(k,i,j)
          gb%uc(ib,k) = uc(k,i,j)
          gb%vc(ib,k) = vc(k,i,j)
          gb%wc(ib,k) = wc(k,i,j)
          gb%ut(ib,k) = ut(k,i,j)
          gb%vt(ib,k) = vt(k,i,j)
          gb%wt(ib,k) = wt(k,i,j)
       end do
    end do
    do iSca = 1, nScalars
       do ib = 1, gb%ibMax
          i = gb%iPerIb(ib)
          j = gb%jPerIb(ib)
          do k = 1, gb%kMax
             gb%scalarp(ib,k,iSca) = scalarp(k,i,j,iSca)
             gb%scalart(ib,k,iSca) = scalart(k,i,j,iSca)
          end do
       end do
    end do
  end subroutine PutTimeVarGhostBlock



  subroutine PutTimeVarVGhostBlock(up, vp, wp, uc, vc, wc, &
       ut, vt, wt, gb)

    ! Arguments:
    real,    intent(in) :: up(:,:,:)
    real,    intent(in) :: vp(:,:,:)
    real,    intent(in) :: wp(:,:,:)
    real,    intent(in) :: uc(:,:,:)
    real,    intent(in) :: vc(:,:,:)
    real,    intent(in) :: wc(:,:,:)
    real,    intent(in) :: ut(:,:,:)
    real,    intent(in) :: vt(:,:,:)
    real,    intent(in) :: wt(:,:,:)
    type(GhostBlock), pointer :: gb

    ! Local Variables:
    integer :: ib, i, j, k

    do ib = 1, gb%ibMax
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       do k = 1, gb%kMax
          gb%up(ib,k) = up(k,i,j)
          gb%vp(ib,k) = vp(k,i,j)
          gb%wp(ib,k) = wp(k,i,j)
          gb%uc(ib,k) = uc(k,i,j)
          gb%vc(ib,k) = vc(k,i,j)
          gb%wc(ib,k) = wc(k,i,j)
          gb%ut(ib,k) = ut(k,i,j)
          gb%vt(ib,k) = vt(k,i,j)
          gb%wt(ib,k) = wt(k,i,j)
       end do
    end do

  end subroutine PutTimeVarVGhostBlock




  subroutine PutTimeVarSGhostBlock(up, vp, wp, uc, vc, wc, &
       scalarp, scalart, gb)
    real,    intent(in) :: up(:,:,:)
    real,    intent(in) :: vp(:,:,:)
    real,    intent(in) :: wp(:,:,:)
    real,    intent(in) :: uc(:,:,:)
    real,    intent(in) :: vc(:,:,:)
    real,    intent(in) :: wc(:,:,:)
    real,    intent(in) :: scalarp(:,:,:,:)
    real,    intent(in) :: scalart(:,:,:,:)
    type(GhostBlock), pointer :: gb

    integer :: ib, i, j, k, iSca, nScalars

    do ib = 1, gb%ibMax
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       do k = 1, gb%kMax
          gb%up(ib,k) = up(k,i,j)
          gb%vp(ib,k) = vp(k,i,j)
          gb%wp(ib,k) = wp(k,i,j)
          gb%uc(ib,k) = uc(k,i,j)
          gb%vc(ib,k) = vc(k,i,j)
          gb%wc(ib,k) = wc(k,i,j)
       enddo
    enddo

    nScalars=size(scalarp,4)

    do iSca = 1, nScalars
       do ib = 1, gb%ibMax
          i = gb%iPerIb(ib)
          j = gb%jPerIb(ib)
          do k = 1, gb%kMax
             gb%scalarp(ib,k,iSca) = scalarp(k,i,j,iSca)
             gb%scalart(ib,k,iSca) = scalart(k,i,j,iSca)
          end do
       end do
    end do
  end subroutine PutTimeVarSGhostBlock





  subroutine GetTimeVarGhostBlock(up, vp, wp, uc, vc, wc, &
       ut, vt, wt, dn0, dn0u, dn0v, scalarp, scalart, gb)
    real,    intent(inout) :: up(:,:,:)
    real,    intent(inout) :: vp(:,:,:)
    real,    intent(inout) :: wp(:,:,:)
    real,    intent(inout) :: uc(:,:,:)
    real,    intent(inout) :: vc(:,:,:)
    real,    intent(inout) :: wc(:,:,:)
    real,    intent(inout) :: ut(:,:,:)
    real,    intent(inout) :: vt(:,:,:)
    real,    intent(inout) :: wt(:,:,:)
    real,    intent(inout) :: dn0(:,:,:)
    real,    intent(inout) :: dn0u(:,:,:)
    real,    intent(inout) :: dn0v(:,:,:)
    real,    intent(inout) :: scalarp(:,:,:,:)
    real,    intent(inout) :: scalart(:,:,:,:)
    type(GhostBlock), pointer :: gb

    integer :: ib, i, j, k, iSca, nScalars

    nScalars=size(scalarp,4)

    do ib = gb%ibStart, gb%ibEnd
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       do k = 1, gb%kMax
          up(k,i,j) = gb%up(ib,k)
          vp(k,i,j) = gb%vp(ib,k)
          wp(k,i,j) = gb%wp(ib,k)
          uc(k,i,j) = gb%uc(ib,k)
          vc(k,i,j) = gb%vc(ib,k)
          wc(k,i,j) = gb%wc(ib,k)
          ut(k,i,j) = gb%ut(ib,k)
          vt(k,i,j) = gb%vt(ib,k)
          wt(k,i,j) = gb%wt(ib,k)
          dn0(k,i,j) = gb%dn0(ib,k)
          dn0u(k,i,j) = gb%dn0u(ib,k)
          dn0v(k,i,j) = gb%dn0v(ib,k)
       end do
    end do
    do iSca = 1, nScalars
       do ib = gb%ibStart, gb%ibEnd
          i = gb%iPerIb(ib)
          j = gb%jPerIb(ib)
          do k = 1, gb%kMax
             scalarp(k,i,j,iSca) = gb%scalarp(ib,k,iSca) 
             scalart(k,i,j,iSca) = gb%scalart(ib,k,iSca) 
          end do
       end do
    end do
  end subroutine GetTimeVarGhostBlock




  subroutine GetTimeVarGhostBlock_oper(up, vp, wp, uc, vc, wc, &
       ut, vt, wt, dn0, dn0u, dn0v, scalarp, scalart, gb)
    real,    intent(in   ) :: up(:,:,:)
    real,    intent(in   ) :: vp(:,:,:)
    real,    intent(in   ) :: wp(:,:,:)
    real,    intent(in   ) :: uc(:,:,:)
    real,    intent(in   ) :: vc(:,:,:)
    real,    intent(in   ) :: wc(:,:,:)
    real,    intent(inout) :: ut(:,:,:)
    real,    intent(inout) :: vt(:,:,:)
    real,    intent(inout) :: wt(:,:,:)
    real,    intent(in   ) :: dn0(:,:,:)
    real,    intent(in   ) :: dn0u(:,:,:)
    real,    intent(in   ) :: dn0v(:,:,:)
    real,    intent(in   ) :: scalarp(:,:,:,:)
    real,    intent(inout) :: scalart(:,:,:,:)
    type(GhostBlock), pointer :: gb

    integer :: ib, i, j, k, iSca, nScalars

    nScalars=size(scalarp,4)

    do ib = gb%ibStart, gb%ibEnd
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       do k = 1, gb%kMax
          ut(k,i,j) = gb%ut(ib,k)
          vt(k,i,j) = gb%vt(ib,k)
          wt(k,i,j) = gb%wt(ib,k)
       end do
    end do
    do iSca = 1, nScalars
       do ib = gb%ibStart, gb%ibEnd
          i = gb%iPerIb(ib)
          j = gb%jPerIb(ib)
          do k = 1, gb%kMax
             scalart(k,i,j,iSca) = gb%scalart(ib,k,iSca) 
          end do
       end do
    end do
  end subroutine GetTimeVarGhostBlock_oper



  subroutine GetTimeVarVGhostBlock_oper(ut, vt, wt, gb)
    real,    intent(inout) :: ut(:,:,:)
    real,    intent(inout) :: vt(:,:,:)
    real,    intent(inout) :: wt(:,:,:)
    type(GhostBlock), pointer :: gb

    integer :: ib, i, j, k

    do ib = gb%ibStart, gb%ibEnd
       i = gb%iPerIb(ib)
       j = gb%jPerIb(ib)
       do k = 1, gb%kMax
          ut(k,i,j) = gb%ut(ib,k)
          vt(k,i,j) = gb%vt(ib,k)
          wt(k,i,j) = gb%wt(ib,k)
       end do
    end do
  end subroutine GetTimeVarVGhostBlock_oper



  subroutine GetTimeVarSGhostBlock_oper(scalart, gb)
    real,    intent(inout) :: scalart(:,:,:,:)
    type(GhostBlock), pointer :: gb

    integer :: ib, i, j, k, iSca, nScalars

    nScalars=size(scalart,4)

    do iSca = 1, nScalars
       do ib = gb%ibStart, gb%ibEnd
          i = gb%iPerIb(ib)
          j = gb%jPerIb(ib)
          do k = 1, gb%kMax
             scalart(k,i,j,iSca) = gb%scalart(ib,k,iSca) 
          end do
       end do
    end do
  end subroutine GetTimeVarSGhostBlock_oper


end module mod_GhostBlock
