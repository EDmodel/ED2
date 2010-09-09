module extras

  ! Used in CATT

  implicit none

  integer, parameter :: na_extra2d=4, na_extra3d=6

  type ext2d
     real,pointer,dimension(:,:) :: d2
  end type ext2d
  type ext3d
     real,pointer,dimension(:,:,:) :: d3
  end type ext3d

  type(ext2d),allocatable :: extra2d(:,:),extra2dm(:,:)
  ! extrad3d(indice,ngrid)
  type(ext3d),allocatable :: extra3d(:,:),extra3dm(:,:)

contains

  subroutine alloc_extra2d(scal,m1,m2,na2d,ngrid)

    implicit none

    type (ext2d),intent(INOUT) :: scal(:,:)
    integer,intent(IN) :: m1,m2 !Dimension of arrays
    integer,intent(IN) :: na2d ! number of 2d extras arrays without ngrid
    integer,intent(IN) :: ngrid
    integer :: j

    do j=1,na2d
       allocate(scal(j,ngrid)%d2(m1,m2))
    end do

  end subroutine alloc_extra2d

  !---------------------------------------------------------------

  subroutine alloc_extra3d(scal,m1,m2,m3,na3d,ngrid)

    implicit none

    type (ext3d),intent(INOUT) :: scal(:,:)
    integer,intent(IN) :: m1,m2,m3 !Dimension of arrays
    integer,intent(IN) :: na3d ! number of 2d extras arrays without ngrid
    integer,intent(IN) :: ngrid
    integer :: j

    do j=1,na3d
       allocate(scal(j,ngrid)%d3(m1,m2,m3))
    end do
  end subroutine alloc_extra3d

  !---------------------------------------------------------------

  subroutine dealloc_extra2d(scal,na2d,ngrid)

    implicit none

    type (ext2d),intent(INOUT) :: scal(:,:)
    integer,intent(in) :: ngrid,na2d
    integer :: nsc,nd

    !  Deallocate arrays

    do nd=1,na2d
       do nsc=1,ngrid
          if (associated(scal(nd,nsc)%d2))   deallocate (scal(nd,nsc)%d2)
       enddo
    end do

  end subroutine dealloc_extra2d

  !---------------------------------------------------------------

  subroutine dealloc_extra3d(scal,na3d,ngrid)

    implicit none

    type (ext3d),intent(INOUT) :: scal(:,:)
    integer,intent(in) :: ngrid,na3d
    integer :: nsc,nd

    !  Deallocate arrays

    do nd=1,na3d
       do nsc=1,ngrid
          if (associated(scal(nd,nsc)%d3))   deallocate (scal(nd,nsc)%d3)
       enddo
    end do

  end subroutine dealloc_extra3d

  !---------------------------------------------------------------

  subroutine nullify_extra2d(scal,na2d,ngrid)

    implicit none

    type (ext2d),intent(INOUT) :: scal(:,:)
    integer,intent(in) :: na2d,ngrid
    integer :: nsc,nd

    !  Deallocate arrays

    do nd=1,na2d
       do nsc=1,ngrid
          if (associated(scal(nd,ngrid)%d2))   nullify (scal(nd,ngrid)%d2)
       enddo
    enddo

    return
  end subroutine nullify_extra2d

  !---------------------------------------------------------------


  subroutine nullify_extra3d(scal,na3d,ngrid)

    implicit none

    type (ext3d),intent(INOUT) :: scal(:,:)
    integer,intent(in) :: na3d,ngrid
    integer :: nsc,nd

    !  Deallocate arrays

    do nd=1,na3d
       do nsc=1,ngrid
          if (associated(scal(nd,ngrid)%d3))   nullify (scal(nd,ngrid)%d3)
       enddo
    enddo

    return
  end subroutine nullify_extra3d

  !---------------------------------------------------------------
  subroutine filltab_extra2d(scal2,scalm2,imean,n1,n2,ng,na)

    use var_tables

    implicit none

    type (ext2d) :: scal2,scalm2
    integer, intent(in) :: imean,n1,n2,ng,na

    integer :: npts
    character (len=7) :: sname

    ! Fill pointers to arrays into variable tables

    if (associated(scal2%d2)) then
       npts=n1*n2
       write(sname,'(a2,i3.3)') 'd2', na
       call vtables2 (scal2%d2(1,1),scalm2%d2(1,1)  &
            ,ng, npts, imean,  &
            trim(sname)//' :2:hist:anal:mpti:mpt3') ! Default - Column oriented Proc.

    endif

  end subroutine filltab_extra2d

  !---------------------------------------------------------------
  subroutine filltab_extra3d(scal3,scalm3,imean,n1,n2,n3,ng,na)

    use var_tables

    implicit none

    type (ext3d) :: scal3,scalm3
    integer, intent(in) :: imean,n1,n2,n3,ng,na

    integer :: npts
    character (len=7) :: sname

    ! Fill pointers to arrays into variable tables

    if (associated(scal3%d3)) then
       npts=n1*n2*n3
       write(sname,'(a2,i3.3)') 'd3', na
       call vtables2 (scal3%d3(1,1,1),scalm3%d3(1,1,1)  &
            ,ng, npts, imean,  &
            sname//' :3:hist:anal:mpti:mpt3') ! Default - Column oriented Proc.

    endif

    return
  end subroutine filltab_extra3d

  !-----------------------------------------------------------------
  subroutine zero_extra3d(scal,na3d,ngrid)

    implicit none

    type (ext3d),intent(INOUT) :: scal(:,:)
    integer,intent(in) :: na3d,ngrid
    integer :: nsc

    !  Deallocate arrays

    do nsc=1,na3d
       scal(nsc,ngrid)%d3(:,:,:)=0.
    enddo

    return
  end subroutine zero_extra3d

  !-----------------------------------------------------------------
  subroutine zero_extra2d(scal,na2d,ngrid)

    implicit none

    type (ext2d),intent(INOUT) :: scal(:,:)
    integer,intent(in) :: na2d,ngrid
    integer :: nsc

    !  Deallocate arrays

    do nsc=1,na2d
       scal(nsc,ngrid)%d2(:,:)=0.
    enddo

    return
  end subroutine zero_extra2d


end module extras
