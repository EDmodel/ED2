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

      nullify(scrvar%c)
      nullify(scrvar%d)
      nullify(scrvar%e)
      nullify(scrvar%f)
      nullify(scrvar%g)
      nullify(scrvar%h)
      nullify(scrvar%i)
      nullify(scrvar%j)
      nullify(scrvar%k)
      nullify(scrvar%l)
      nullify(scrvar%m)
      nullify(scrvar%n)
      nullify(scrvar%o)
      nullify(scrvar%p)
      nullify(scrvar%q)
      nullify(scrvar%r)
      nullify(scrvar%s)
      nullify(scrvar%t)
      nullify(scrvar%u)
      nullify(scrvar%v)
      nullify(scrvar%w)
      nullify(scrvar%x)
      nullify(scrvar%y)
      nullify(scrvar%z)
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






   !=======================================================================================!
   !=======================================================================================!
   subroutine zero_scratch(scrvar)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(scratch_vars), intent(inout) :: scrvar
      !------------------------------------------------------------------------------------!

      if(associated(scrvar%c)) scrvar%c(:) = 0.0
      if(associated(scrvar%d)) scrvar%d(:) = 0.0
      if(associated(scrvar%e)) scrvar%e(:) = 0.0
      if(associated(scrvar%f)) scrvar%f(:) = 0.0
      if(associated(scrvar%g)) scrvar%g(:) = 0.0
      if(associated(scrvar%h)) scrvar%h(:) = 0.0
      if(associated(scrvar%i)) scrvar%i(:) = 0.0
      if(associated(scrvar%j)) scrvar%j(:) = 0.0
      if(associated(scrvar%k)) scrvar%k(:) = 0.0
      if(associated(scrvar%l)) scrvar%l(:) = 0.0
      if(associated(scrvar%m)) scrvar%m(:) = 0.0
      if(associated(scrvar%n)) scrvar%n(:) = 0.0
      if(associated(scrvar%o)) scrvar%o(:) = 0.0
      if(associated(scrvar%p)) scrvar%p(:) = 0.0
      if(associated(scrvar%q)) scrvar%q(:) = 0.0
      if(associated(scrvar%r)) scrvar%r(:) = 0.0
      if(associated(scrvar%s)) scrvar%s(:) = 0.0
      if(associated(scrvar%t)) scrvar%t(:) = 0.0
      if(associated(scrvar%u)) scrvar%u(:) = 0.0
      if(associated(scrvar%v)) scrvar%v(:) = 0.0
      if(associated(scrvar%w)) scrvar%w(:) = 0.0
      if(associated(scrvar%x)) scrvar%x(:) = 0.0
      if(associated(scrvar%y)) scrvar%y(:) = 0.0
      if(associated(scrvar%z)) scrvar%z(:) = 0.0
   end subroutine zero_scratch
   !=======================================================================================!
   !=======================================================================================!

end module scratch_coms
!==========================================================================================!
!==========================================================================================!
