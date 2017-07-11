module mod_GhostBlockPartition
  implicit none

  type GhostBlockPartition
     integer          :: iMax
     integer          :: jMax
     integer          :: nBlk
     integer, pointer :: blockSize(:)
     integer, pointer :: ijStart(:)
     integer, pointer :: ijEnd(:)
  end type GhostBlockPartition
contains



  ! CreateGhostBlockPartition:
  !    given ij plane dimensions (including ghost zone), partition
  !    ij plane into contiguous ij segments of length at most sizeBlk
  !    returns pointer variable of type GhostBlockPartition



  function CreateGhostBlockPartition (iMax, jMax, sizeBlk) result(res)
    integer, intent(in) :: iMax
    integer, intent(in) :: jMax
    integer, intent(in) :: sizeBlk
    type(GhostBlockPartition), pointer :: res

    character(len=*), parameter :: h="**(CreateGhostBlockPartition)**"
    integer :: innerSize
    integer :: nBlkTent
    integer :: rem
    integer :: ij, ijRem
    integer :: ib, jb

    ! input data consistency

    if (iMax <=0) then
       write(*,"(a,i10)") h//" iMax not positive, =",iMax
       stop
    else if (jMax <=0) then
       write(*,"(a,i10)") h//" jMax not positive, =",jMax
       stop
    else if (sizeBlk <=0) then
       write(*,"(a,i10)") h//" sizeBlk not positive, =",sizeBlk
       stop
    end if

    ! innerSize: # horizontal points from first real data to last
    !            real data counting internal ghost zone

    innerSize = (jMax-2)*iMax - 2 ! inner cols with ghost - both extremes

    ! # blocks to cover inner ghost zone

    nBlkTent = innerSize/sizeBlk
    rem  = mod(innerSize,sizeBlk)
    if (rem /= 0) then
       nBlkTent = nBlkTent + 1
    end if

    allocate(res)
    allocate(res%blockSize(nBlkTent))
    allocate(res%ijStart(nBlkTent))
    allocate(res%ijEnd(nBlkTent))

    ! store input data

    res%iMax = iMax
    res%jMax = jMax

    ! blockSize(jb): # of innerSize points per block

    do jb = 1, nBlkTent-1
       res%blockSize(jb) = sizeBlk
    end do

    if (rem == 0) then
       res%blockSize(nBlkTent) = sizeBlk
    else
       res%blockSize(nBlkTent) = rem
    end if

    ! visit all blocks, defining starting and ending real data points
    ! considering ghost zone; avoid starting and ending points 
    ! in internal ghost zone. Points are counted in the entire 
    ! ij plane.
    ! ij plane elements counted from 0 in array element order.

    ij = iMax + 1
    jb = 1
    do 

       ! map starting point; advance if ghost zone

       ijRem = mod(ij,iMax)
       if (ijRem == 0) then
          res%ijStart(jb) = ij + 1
       else if (ijRem == iMax-1) then
          res%ijStart(jb) = ij + 2
       else
          res%ijStart(jb) = ij
       end if

       ! map ending point; garantee inbounds; backup if ghost zone

       ij = min(res%ijStart(jb) + res%blockSize(jb) - 1, innerSize+iMax)
       ijRem = mod(ij,iMax)
       if (ijRem == 0) then
          res%ijEnd(jb) = ij - 2
       else if (ijRem == iMax-1) then
          res%ijEnd(jb) = ij - 1
       else
          res%ijEnd(jb) = ij
       end if

       ij = ij + 1

       ! if non-empty block, recompute block size and 
       !    go to next block;

       if (res%ijEnd(jb) >= res%ijStart(jb)) then
          res%blockSize(jb) = res%ijEnd(jb)-res%ijStart(jb)+1
          jb = jb + 1
       end if

       ! exit if exausted inner points

       if (ij > innerSize + iMax) then
          exit
       end if
    end do

    ! total number of blocks

    res%nBlk = jb - 1
  end function CreateGhostBlockPartition



  ! DestroyGhostBlockPartition
  !    destroy pointer variable of type GhostBlockPartition



  subroutine DestroyGhostBlockPartition(p)
    type(GhostBlockPartition), pointer :: p

    if (associated(p)) then
       if (associated(p%blockSize)) then
          deallocate(p%blockSize)
       end if
       if (associated(p%ijStart)) then
          deallocate(p%ijStart)
       end if
       if (associated(p%ijEnd)) then
          deallocate(p%ijEnd)
       end if
       deallocate(p)
    end if
    nullify(p)
  end subroutine DestroyGhostBlockPartition



  ! DumpGhostBlockPartition
  !    dumps pointer variable of type GhostBlockPartition



  subroutine DumpGhostBlockPartition(p)
    type(GhostBlockPartition), pointer :: p

    integer :: i
    integer :: sizeBlockSize, sizeIJStart, sizeIJEnd, nBlk
    character(len=*), parameter :: h="**(DumpGhostBlockPartition)**"
    character(len=20) :: c0, c1, c2

    if (associated(p)) then
       if ( .not. associated(p%blockSize) .or. &
            .not. associated(p%ijStart)   .or. &
            .not. associated(p%ijEnd)    ) then
          write(*,"(a)") h//" GhostBlockPartition not fully created"
       else
          nBlk = p%nBlk
          sizeBlockSize = size(p%blockSize)
          sizeIJStart   = size(p%ijStart)
          sizeIJEnd     = size(p%ijEnd)
          if ( (sizeBlockSize < nBlk) .or. &
               (sizeIJStart   < nBlk) .or. &
               (sizeIJEnd     < nBlk)   ) then
             write(*,"(a)") h//" inconsistent GhostBlockPartition "
          else
             write(c0,"(i20)") p%nBlk
             write(c1,"(i20)") p%iMax
             write(c2,"(i20)") p%jMax
             write(*,"(a)") h//" with "//trim(adjustl(c0))//" blocks; iMax="//&
                  trim(adjustl(c1))//"; jMax="//trim(adjustl(c2))
             write(*,"(a)") h//" block   blockSize     ijStart       ijEnd"
             do i = 1, nBlk
                write(*,"(a,i6,3i12)") h, i, p%blockSize(i), p%ijStart(i), p%ijEnd(i)
             end do
          end if
       end if
    else
       write(*,"(a)") h//" GhostBlockPartition not associated"
    end if
  end subroutine DumpGhostBlockPartition
end module mod_GhostBlockPartition
