!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_scalar

  ! Added scalar variables and tendencies

  type scalar_vars
     real, pointer, dimension(:,:,:) :: sclp
     real, pointer, dimension(:,:) :: drydep
     real, pointer, dimension(:) :: sclt
     ! For CATT
     real, pointer, dimension(:,:) :: wetdep
     real, pointer, dimension(:,:,:) :: srcsc
  end type scalar_vars

  !    scal_p allocated by (maxsclr,ngrids)
  type (scalar_vars), allocatable :: scalar_g(:,:),scalarm_g(:,:)

  ! For CATT
  integer :: RECYCLE_TRACERS

contains
  !---------------------------------------------------------------

  subroutine alloc_scalar(scal,n1,n2,n3,naddsc)

    implicit none

    integer,intent(in) :: naddsc
    type (scalar_vars),dimension(naddsc) :: scal
    integer,intent(in) :: n1,n2,n3

    integer :: nsc

    ! print *,'Size of scal=' ,size(scal,1)
    ! Allocate arrays based on options (if necessary).
    do nsc=1,naddsc
       !print*,'escalar=',nsc,naddsc,n1,n2,n3
       allocate (scal(nsc)%sclp(n1,n2,n3))
       allocate (scal(nsc)%drydep(n2,n3))
       allocate (scal(nsc)%wetdep(n2,n3))
       allocate (scal(nsc)%srcsc(n1,n2,n3))
    enddo

    return
  end subroutine alloc_scalar

  !---------------------------------------------------------------

  subroutine dealloc_scalar(scal,naddsc)

    implicit none

    type (scalar_vars) :: scal(naddsc)
    integer :: naddsc
    integer :: nsc

    !  Deallocate arrays

    do nsc=1,naddsc
       if (associated(scal(nsc)%sclp))   deallocate (scal(nsc)%sclp)
       if (associated(scal(nsc)%drydep)) deallocate (scal(nsc)%drydep)
       ! For CATT
       if (associated(scal(nsc)%wetdep)) deallocate (scal(nsc)%wetdep)
       if (associated(scal(nsc)%srcsc)) deallocate (scal(nsc)%srcsc)
    enddo

    return
  end subroutine dealloc_scalar

  !---------------------------------------------------------------

  subroutine nullify_scalar(scal,naddsc)

    implicit none

    type (scalar_vars) :: scal(naddsc)

    integer :: naddsc
    integer :: nsc

    !  Deallocate arrays

    do nsc=1,naddsc
       if (associated(scal(nsc)%sclp))   nullify (scal(nsc)%sclp)
       if (associated(scal(nsc)%drydep)) nullify (scal(nsc)%drydep)
       ! For CATT
       if (associated(scal(nsc)%wetdep)) nullify (scal(nsc)%wetdep)
       if (associated(scal(nsc)%srcsc)) nullify (scal(nsc)%srcsc)
    enddo

    return
  end subroutine nullify_scalar

  !---------------------------------------------------------------

  subroutine filltab_scalar(scal,scalm,imean,n1,n2,n3,ng,na)

    use var_tables

    implicit none

    type (scalar_vars) :: scal,scalm
    integer, intent(in) :: imean,n1,n2,n3,ng,na

    integer :: npts  !,nsc,nptg
    character (len=15) :: sname

    ! ALF
    character(len=8) :: str_recycle

    ! ALF
    str_recycle = ''
    if (RECYCLE_TRACERS == 1) then
       str_recycle = ':recycle'
    endif

    ! Fill pointers to arrays into variable tables

    if (associated(scal%sclp)) then
       npts=n1*n2*n3

       write(sname,'(a4,i3.3)') 'SCLP',na
       call vtables2 (scal%sclp(1,1,1),scalm%sclp(1,1,1)  &
            ,ng, npts, imean,  &
            trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1'//trim(str_recycle))

       npts=n2*n3

       write(sname,'(a4,i3.3)') 'SCDD',na
       call vtables2 (scal%drydep(1,1),scalm%drydep(1,1)  &
            ,ng, npts, imean,  &
            trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1')

       write(sname,'(a6,i3.3)') 'wetdep',na
       call vtables2 (scal%wetdep(1,1),scalm%wetdep(1,1)  &
            ,ng, npts, imean,  &
            trim(sname)//' :2:hist:anal:mpti:mpt3:mpt1')

       npts=n1*n2*n3

       write(sname,'(a5,i3.3)') 'scrsc',na
       call vtables2 (scal%srcsc(1,1,1),scalm%srcsc(1,1,1)  &
            ,ng, npts, imean,  &
            trim(sname)//' :3:hist:anal:mpti:mpt3:mpt1')
    endif

  end subroutine filltab_scalar

end module mem_scalar
