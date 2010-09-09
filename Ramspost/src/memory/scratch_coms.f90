!==========================================================================================!
!==========================================================================================!
!    This module contains some variables that temporarily host variables during the post-  !
! processing.                                                                              !
!------------------------------------------------------------------------------------------!
module scratch_coms

   type scratch_vars
      real, dimension(:), pointer :: c
      real, dimension(:), pointer :: d
      real, dimension(:), pointer :: e
      real, dimension(:), pointer :: f
      real, dimension(:), pointer :: g
      real, dimension(:), pointer :: h
      real, dimension(:), pointer :: i
      real, dimension(:), pointer :: j
      real, dimension(:), pointer :: k
      real, dimension(:), pointer :: l
      real, dimension(:), pointer :: m
      real, dimension(:), pointer :: n
      real, dimension(:), pointer :: o
      real, dimension(:), pointer :: p
      real, dimension(:), pointer :: q
      real, dimension(:), pointer :: r
      real, dimension(:), pointer :: s
      real, dimension(:), pointer :: t
      real, dimension(:), pointer :: u
      real, dimension(:), pointer :: v
      real, dimension(:), pointer :: w
      real, dimension(:), pointer :: x
      real, dimension(:), pointer :: y
      real, dimension(:), pointer :: z
   end type scratch_vars

   !----- This is the scratch structure. --------------------------------------------------!
   type(scratch_vars) :: scr

   
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   subroutine alloc_scratch(scrvar,ndim)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(scratch_vars), intent(inout) :: scrvar
      integer           , intent(in)    :: ndim
      !------------------------------------------------------------------------------------!

      call dealloc_scratch(scrvar)
      call nullify_scratch(scrvar)

      allocate(scrvar%c(ndim))
      allocate(scrvar%d(ndim))
      allocate(scrvar%e(ndim))
      allocate(scrvar%f(ndim))
      allocate(scrvar%g(ndim))
      allocate(scrvar%h(ndim))
      allocate(scrvar%i(ndim))
      allocate(scrvar%j(ndim))
      allocate(scrvar%k(ndim))
      allocate(scrvar%l(ndim))
      allocate(scrvar%m(ndim))
      allocate(scrvar%n(ndim))
      allocate(scrvar%o(ndim))
      allocate(scrvar%p(ndim))
      allocate(scrvar%q(ndim))
      allocate(scrvar%r(ndim))
      allocate(scrvar%s(ndim))
      allocate(scrvar%t(ndim))
      allocate(scrvar%u(ndim))
      allocate(scrvar%v(ndim))
      allocate(scrvar%w(ndim))
      allocate(scrvar%x(ndim))
      allocate(scrvar%y(ndim))
      allocate(scrvar%z(ndim))

      return
   end subroutine alloc_scratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine nullify_scratch(scrvar)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(scratch_vars), intent(inout) :: scrvar
      !------------------------------------------------------------------------------------!

      if(associated(scrvar%c)) nullify(scrvar%c)
      if(associated(scrvar%d)) nullify(scrvar%d)
      if(associated(scrvar%e)) nullify(scrvar%e)
      if(associated(scrvar%f)) nullify(scrvar%f)
      if(associated(scrvar%g)) nullify(scrvar%g)
      if(associated(scrvar%h)) nullify(scrvar%h)
      if(associated(scrvar%i)) nullify(scrvar%i)
      if(associated(scrvar%j)) nullify(scrvar%j)
      if(associated(scrvar%k)) nullify(scrvar%k)
      if(associated(scrvar%l)) nullify(scrvar%l)
      if(associated(scrvar%m)) nullify(scrvar%m)
      if(associated(scrvar%n)) nullify(scrvar%n)
      if(associated(scrvar%o)) nullify(scrvar%o)
      if(associated(scrvar%p)) nullify(scrvar%p)
      if(associated(scrvar%q)) nullify(scrvar%q)
      if(associated(scrvar%r)) nullify(scrvar%r)
      if(associated(scrvar%s)) nullify(scrvar%s)
      if(associated(scrvar%t)) nullify(scrvar%t)
      if(associated(scrvar%u)) nullify(scrvar%u)
      if(associated(scrvar%v)) nullify(scrvar%v)
      if(associated(scrvar%w)) nullify(scrvar%w)
      if(associated(scrvar%x)) nullify(scrvar%x)
      if(associated(scrvar%y)) nullify(scrvar%y)
      if(associated(scrvar%z)) nullify(scrvar%z)
   end subroutine nullify_scratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine dealloc_scratch(scrvar)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(scratch_vars), intent(inout) :: scrvar
      !------------------------------------------------------------------------------------!

      if(associated(scrvar%c)) deallocate(scrvar%c)
      if(associated(scrvar%d)) deallocate(scrvar%d)
      if(associated(scrvar%e)) deallocate(scrvar%e)
      if(associated(scrvar%f)) deallocate(scrvar%f)
      if(associated(scrvar%g)) deallocate(scrvar%g)
      if(associated(scrvar%h)) deallocate(scrvar%h)
      if(associated(scrvar%i)) deallocate(scrvar%i)
      if(associated(scrvar%j)) deallocate(scrvar%j)
      if(associated(scrvar%k)) deallocate(scrvar%k)
      if(associated(scrvar%l)) deallocate(scrvar%l)
      if(associated(scrvar%m)) deallocate(scrvar%m)
      if(associated(scrvar%n)) deallocate(scrvar%n)
      if(associated(scrvar%o)) deallocate(scrvar%o)
      if(associated(scrvar%p)) deallocate(scrvar%p)
      if(associated(scrvar%q)) deallocate(scrvar%q)
      if(associated(scrvar%r)) deallocate(scrvar%r)
      if(associated(scrvar%s)) deallocate(scrvar%s)
      if(associated(scrvar%t)) deallocate(scrvar%t)
      if(associated(scrvar%u)) deallocate(scrvar%u)
      if(associated(scrvar%v)) deallocate(scrvar%v)
      if(associated(scrvar%w)) deallocate(scrvar%w)
      if(associated(scrvar%x)) deallocate(scrvar%x)
      if(associated(scrvar%y)) deallocate(scrvar%y)
      if(associated(scrvar%z)) deallocate(scrvar%z)
   end subroutine dealloc_scratch
   !=======================================================================================!
   !=======================================================================================!

end module scratch_coms
!==========================================================================================!
!==========================================================================================!
