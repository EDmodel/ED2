module ed_work_vars

  type work_vars

     ! Variables to be dimensioned by (nxp,nyp)

     real,    pointer :: glon(:,:)
     real,    pointer :: glat(:,:)
     real,    pointer :: work(:,:)
     logical, pointer :: land(:,:)
     real,    pointer :: landfrac(:,:)
     integer, pointer :: ntext(:,:)
     integer, pointer :: xatm(:,:)
     integer, pointer :: yatm(:,:)


  end type work_vars


  type (work_vars), allocatable :: work_e(:)

contains
!==========================================================================================!
!==========================================================================================!
  subroutine ed_alloc_work(worke,n2,n3)
    implicit none
    type (work_vars) :: worke
    integer, intent(in) :: n2,n3

    ! Allocate arrays based on options (if necessary)
    allocate (worke%glon(n2,n3))
    allocate (worke%glat(n2,n3))
    allocate (worke%xatm(n2,n3))
    allocate (worke%yatm(n2,n3))
    allocate (worke%work(n2,n3))
    allocate (worke%land(n2,n3))
    allocate (worke%landfrac(n2,n3))
    allocate (worke%ntext(n2,n3))
    return
  end subroutine ed_alloc_work
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
  subroutine ed_nullify_work(worke)
    implicit none
    type (work_vars) :: worke

    if (associated(worke%glon) ) nullify (worke%glon)
    if (associated(worke%glat) ) nullify (worke%glat)
    if (associated(worke%work) ) nullify (worke%work)
    if (associated(worke%land) ) nullify (worke%land)
    if (associated(worke%landfrac) ) nullify (worke%landfrac)
    if (associated(worke%ntext)) nullify (worke%ntext)
    if (associated(worke%xatm) ) nullify (worke%xatm)
    if (associated(worke%yatm) ) nullify (worke%yatm)

    return
  end subroutine ed_nullify_work
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
  subroutine ed_dealloc_work(worke)
    implicit none
    type (work_vars) :: worke

    if (associated(worke%glon)  )    deallocate (worke%glon)
    if (associated(worke%glat)  )    deallocate (worke%glat)
    if (associated(worke%work)  )    deallocate (worke%work)
    if (associated(worke%land)  )    deallocate (worke%land)
    if (associated(worke%landfrac)  )    deallocate (worke%landfrac)
    if (associated(worke%ntext) )    deallocate (worke%ntext)
    if (associated(worke%xatm)  )    deallocate (worke%xatm)
    if (associated(worke%yatm)  )    deallocate (worke%yatm)
    return
  end subroutine ed_dealloc_work

end module ed_work_vars
