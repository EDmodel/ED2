!==========================================================================================!
!==========================================================================================!
!      This module contains some functions to generate random numbers.                     !
!------------------------------------------------------------------------------------------!
module random_utils
   implicit none

   !=======================================================================================!
   !=======================================================================================!



   contains


   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine initialises the random seed from the system clock at every run.   !
   ! You must call this at least once during the main program execution if you don't want  !
   ! results to look the same.                                                             !
   !                                                                                       !
   !     This sub-routine has been borrowed from:                                          !
   !     http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html                       !
   !---------------------------------------------------------------------------------------!
   subroutine init_random_seed()
      !------------------------------------------------------------------------------------!
      !      GETPID must be explicitly loaded from standard intel modules in case you are  !
      !  compiling with Intel.  In case you get an error, you may want to check your       !
      !  include.mk file, look for CMACH, and replace by one of the options below, or      !
      !  create a unique machine name (likely to require changes in rsys.F90 too).  In     !
      !  case you create a new name, please be creative so it doesn't conflict with other  !
      !  user's choice.                                                                    !
      !------------------------------------------------------------------------------------!
#if defined(ODYSSEY) || defined(SUNHPC) || defined(PC_INTEL)
      use ifport, only : getpid
#endif
      !------------------------------------------------------------------------------------!
      implicit none
      !----- Local variables. -------------------------------------------------------------!
      integer, dimension(:), allocatable :: seed
      integer                            :: i
      integer                            :: n
      integer                            :: ierr
      integer, dimension(8)              :: dt
      integer                            :: pid
      integer(kind=8)                    :: when
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Allocate seed.                                                                 !
      !------------------------------------------------------------------------------------!
      call random_seed(size=n)
      allocate(seed(n))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     First try if the OS provides a random number generator.                        !
      !------------------------------------------------------------------------------------!
      open (unit=72,file='/dev/urandom',access='stream',form='unformatted',action='read'   &
           ,status='old',iostat=ierr)
      select case (ierr)
      case (0)
         !----- System provides RN generator, use it. -------------------------------------!
         read(unit=72) seed
         close(unit=72,status='keep')
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !     System doesn't have random number generator.  Fallback to XOR:ing the       !
         ! current time and pid.  The PID is useful in case one launches multiple          !
         ! instances of the same program in parallel.                                      !
         !---------------------------------------------------------------------------------!
         call system_clock(when)
         if (when == 0) then
            call date_and_time(values=dt)
            when = ( dt(1) - 1970 ) * 365_8 * 24_8 * 60 * 60 * 1000                        &
                 +   dt(2)          * 31_8  * 24_8 * 60 * 60 * 1000                        &
                 +   dt(3)                  * 24_8 * 60 * 60 * 1000                        &
                 +   dt(5)                         * 60 * 60 * 1000                        &
                 +   dt(6)                              * 60 * 1000                        &
                 +   dt(7)                                   * 1000                        &
                 +   dt(8)
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Get Process ID.                                                            !
         !---------------------------------------------------------------------------------!
         pid  = getpid()
         when = ieor(when,int(pid,kind(when)))
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Set seeds.                                                                  !
         !---------------------------------------------------------------------------------!
         do i=1,n
            seed(i) = lcg(when)
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      call random_seed(put=seed)

      deallocate(seed)

      return
   end subroutine init_random_seed
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This simple PRNG might not be good enough for real work, but is sufficient for   !
   ! seeding a better PRNG.                                                                !
   !                                                                                       !
   !      This function has been borrowed from                                             !
   ! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html                           !
   !---------------------------------------------------------------------------------------!
   integer function lcg(s)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer(kind=8), intent(in) :: s
      !----- Local variables. -------------------------------------------------------------!
      integer(kind=8)             :: sloc
      !------------------------------------------------------------------------------------!

      if (s == 0) then
         sloc = 104729_8
      else
         sloc = mod(s, 4294967296_8)
      end if
      sloc = mod(sloc * 279470273_8, 4294967291_8)
      lcg  = int(mod(sloc, int(huge(0), 8)), kind(0))

      return
   end function lcg
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine samples data with or without replacement (real numbers).          !
   !                                                                                       !
   !  Input/Output variables:                                                              !
   ! nxi         (input)           -- Input vector length                                  !
   ! xi          (input)           -- Input vector, which will be used for sampling.       !
   ! nxo         (input)           -- Size of the output vector                            !
   ! xo          (output)          -- Output vector, with samples                          !
   ! replacement (input)           -- Sampling with replacement?                           !
   ! mask        (input, optional) -- Mask to skip some elements.                          !
   !---------------------------------------------------------------------------------------!
   subroutine fsample(nxi,xi,nxo,xo,replacement,mask)
      implicit none

      !----- Required arguments. ----------------------------------------------------------!
      integer                     , intent(in)           :: nxi
      real(kind=4), dimension(nxi), intent(in)           :: xi
      integer                     , intent(in)           :: nxo
      real(kind=4), dimension(nxo), intent(out)          :: xo
      logical                     , intent(in)           :: replacement
      logical     , dimension(nxi), intent(in), optional :: mask
      !----- Local variables. -------------------------------------------------------------!
      integer                                 :: o      ! Counter
      integer                                 :: s      ! Counter
      integer                                 :: nxs    ! # of points still selectable
      real(kind=4)                            :: rndnow ! Current random number
      real(kind=4), dimension(:), allocatable :: xs     ! Input data (still selectable)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise xs.  In case mask is absent, xs=xi and nxs=nxi, otherwise it will   !
      ! be the unmasked entries.  xs and nxs will remain the same in the case of sampling  !
      ! with replacement.  If this is sampling without replacement, we will reorganise the !
      ! vector and reduce nxs so selected items cannot be selected more than once.         !
      !------------------------------------------------------------------------------------!
      if (present(mask)) then
         nxs = count(mask)
         allocate(xs(nxs))
         xs = pack(xi,mask)
      else
         nxs = nxi
         allocate(xs(nxs))
         xs(:) = xi(:)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Sanity check: nxo cannot be greater than nxi if it is sampling without        !
      ! replacement                                                                        !
      !------------------------------------------------------------------------------------!
      if ((.not. replacement) .and. nxi < nxo) then
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a)')       ' Invalid settings (sampling without replacement).'
         write(unit=*,fmt='(a)')       ' NXO must be <= valid input data.'
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a,1x,i6)') ' NXI           = ',nxi
         write(unit=*,fmt='(a,1x,i6)') ' # of valid XI = ',nxs
         write(unit=*,fmt='(a,1x,i6)') ' NXO           = ',nxo
         write(unit=*,fmt='(a,1x,l1)') ' REPLACEMENT   = ',replacement
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         call fatal_error('Invalid sampling size','fsample','random_utils.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Populate the output vector.                                                    !
      !------------------------------------------------------------------------------------!
      oloop: do o=1,nxo
         !----- Randomly select a point (bounded between 1 and nxs). ----------------------!
         call random_number(rndnow)
         s     = 1 + mod(floor(real(nxs)*rndnow),nxs)
         xo(o) = xs(s)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    In case this is a sampling without replacement, copy the last selectable     !
         ! point to the point we have just selected.  Then we shrink the maximum number of !
         ! selectable points.                                                              !
         !---------------------------------------------------------------------------------!
         if (.not. replacement) then
            xs(  s) = xs(nxs)
            nxs     = nxs - 1
         end if
         !---------------------------------------------------------------------------------!
      end do oloop
      !------------------------------------------------------------------------------------!


      !----- Free memory. -----------------------------------------------------------------!
      deallocate(xs)
      !------------------------------------------------------------------------------------!

      return
   end subroutine fsample
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine samples data with or without replacement.                         !
   !                                                                                       !
   !  Input/Output variables:                                                              !
   ! nxi         (input)           -- Input vector length                                  !
   ! xi          (input)           -- Input vector, which will be used for sampling.       !
   ! nxo         (input)           -- Size of the output vector                            !
   ! xo          (output)          -- Output vector, with samples                          !
   ! replacement (input)           -- Sampling with replacement?                           !
   ! mask        (input, optional) -- Mask to skip some elements.                          !
   !---------------------------------------------------------------------------------------!
   subroutine isample(nxi,xi,nxo,xo,replacement,mask)
      implicit none

      !----- Required arguments. ----------------------------------------------------------!
      integer                     , intent(in)           :: nxi
      integer     , dimension(nxi), intent(in)           :: xi
      integer                     , intent(in)           :: nxo
      integer     , dimension(nxo), intent(out)          :: xo
      logical                     , intent(in)           :: replacement
      logical     , dimension(nxi), intent(in), optional :: mask
      !----- Local variables. -------------------------------------------------------------!
      integer                                 :: o      ! Counter
      integer                                 :: s      ! Counter
      integer                                 :: nxs    ! # of points still selectable
      real(kind=4)                            :: rndnow ! Current random number
      integer     , dimension(:), allocatable :: xs     ! Input data (still selectable)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise xs.  In case mask is absent, xs=xi and nxs=nxi, otherwise it will   !
      ! be the unmasked entries.  xs and nxs will remain the same in the case of sampling  !
      ! with replacement.  If this is sampling without replacement, we will reorganise the !
      ! vector and reduce nxs so selected items cannot be selected more than once.         !
      !------------------------------------------------------------------------------------!
      if (present(mask)) then
         nxs = count(mask)
         allocate(xs(nxs))
         xs = pack(xi,mask)
      else
         nxs = nxi
         allocate(xs(nxs))
         xs(:) = xi(:)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Sanity check: nxo cannot be greater than nxi if it is sampling without        !
      ! replacement                                                                        !
      !------------------------------------------------------------------------------------!
      if ((.not. replacement) .and. nxi < nxo) then
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a)')       ' Invalid settings (sampling without replacement).'
         write(unit=*,fmt='(a)')       ' NXO must be <= valid input data.'
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a,1x,i6)') ' NXI           = ',nxi
         write(unit=*,fmt='(a,1x,i6)') ' # of valid XI = ',nxs
         write(unit=*,fmt='(a,1x,i6)') ' NXO           = ',nxo
         write(unit=*,fmt='(a,1x,l1)') ' REPLACEMENT   = ',replacement
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         call fatal_error('Invalid sampling size','isample','random_utils.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Populate the output vector.                                                    !
      !------------------------------------------------------------------------------------!
      oloop: do o=1,nxo
         !----- Randomly select a point (bounded between 1 and nxs). ----------------------!
         call random_number(rndnow)
         s     = 1 + mod(floor(real(nxs)*rndnow),nxs)
         xo(o) = xs(s)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    In case this is a sampling without replacement, copy the last selectable     !
         ! point to the point we have just selected.  Then we shrink the maximum number of !
         ! selectable points.                                                              !
         !---------------------------------------------------------------------------------!
         if (.not. replacement) then
            xs(  s) = xs(nxs)
            nxs     = nxs - 1
         end if
         !---------------------------------------------------------------------------------!
      end do oloop
      !------------------------------------------------------------------------------------!


      !----- Free memory. -----------------------------------------------------------------!
      deallocate(xs)
      !------------------------------------------------------------------------------------!

      return
   end subroutine isample
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function randomly selects one item from a vector.                            !
   !                                                                                       !
   !  Input/Output variables:                                                              !
   ! nxi         (input)           -- Input vector length                                  !
   ! xi          (input)           -- Input vector, which will be used for sampling.       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function fpickone(nxi,xi)
      implicit none

      !----- Required arguments. ----------------------------------------------------------!
      integer                     , intent(in) :: nxi
      real(kind=4), dimension(nxi), intent(in) :: xi
      !----- Local variables. -------------------------------------------------------------!
      integer                                  :: i      ! Counter
      real(kind=4)                             :: rndnow ! Current random number
      !------------------------------------------------------------------------------------!



      !----- Randomly select a point (bounded between 1 and nxs). -------------------------!
      call random_number(rndnow)
      i        = 1 + mod(floor(real(nxi)*rndnow),nxi)
      fpickone = xi(i)
      !------------------------------------------------------------------------------------!

      return
   end function fpickone
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function randomly selects one item from a vector.                            !
   !                                                                                       !
   !  Input/Output variables:                                                              !
   ! nxi         (input)           -- Input vector length                                  !
   ! xi          (input)           -- Input vector, which will be used for sampling.       !
   !---------------------------------------------------------------------------------------!
   integer function ipickone(nxi,xi)
      implicit none

      !----- Required arguments. ----------------------------------------------------------!
      integer                     , intent(in) :: nxi
      integer     , dimension(nxi), intent(in) :: xi
      !----- Local variables. -------------------------------------------------------------!
      integer                                  :: i      ! Counter
      real(kind=4)                             :: rndnow ! Current random number
      !------------------------------------------------------------------------------------!



      !----- Randomly select a point (bounded between 1 and nxs). -------------------------!
      call random_number(rndnow)
      i        = 1 + mod(floor(real(nxi)*rndnow),nxi)
      ipickone = xi(i)
      !------------------------------------------------------------------------------------!

      return
   end function ipickone
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine generates a vector of random numbers bounded by xmin and xmax,    !
   ! assuming uniform distribution.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine runif(nx,xlwr,xupr,x)
      implicit none

      !----- Required arguments. ----------------------------------------------------------!
      integer                     , intent(in)  :: nx     ! Size of the x vector
      real(kind=4)                , intent(in)  :: xlwr   ! Lower bound
      real(kind=4)                , intent(in)  :: xupr   ! Upper bound
      real(kind=4), dimension(nx) , intent(out) :: x      ! Vector with the random numbers
      !----- Local variables. -------------------------------------------------------------!
      integer                                   :: o      ! Counter
      real(kind=4)                              :: rndnow ! Current random number
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Populate the output vector.                                                    !
      !------------------------------------------------------------------------------------!
      oloop: do o=1,nx
         !----- Call random_number (0 <= rndnow < 1). -------------------------------------!
         call random_number(rndnow)
         !---------------------------------------------------------------------------------!


         !------ Scale the random numbers. ------------------------------------------------!
         x(o) = xlwr + (xupr - xlwr) * rndnow
         !---------------------------------------------------------------------------------!
      end do oloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine runif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine generates a single random number bounded by xmin and xmax,        !
   ! assuming uniform distribution.                                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function runif_sca(xlwr,xupr)
      implicit none

      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4)                , intent(in)  :: xlwr   ! Lower bound
      real(kind=4)                , intent(in)  :: xupr   ! Upper bound
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                              :: rndnow ! Current random number
      !------------------------------------------------------------------------------------!



      !----- Call random_number (0 <= rndnow < 1). ----------------------------------------!
      call random_number(rndnow)
      !------------------------------------------------------------------------------------!


      !------ Scale the random number. ----------------------------------------------------!
      runif_sca = xlwr + (xupr - xlwr) * rndnow
      !------------------------------------------------------------------------------------!

      return
   end function runif_sca
   !=======================================================================================!
   !=======================================================================================!
end module random_utils
!==========================================================================================!
!==========================================================================================!
