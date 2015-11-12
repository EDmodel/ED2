!==========================================================================================!
!==========================================================================================!
!     This subroutine initialises the random seed from the system clock at every run.      !
! You must call this at least once during the main program execution if you don't want     !
! results to look the same.                                                                !
!                                                                                          !
!     This sub-routine has been borrowed from:                                             !
!     http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html                          !
!------------------------------------------------------------------------------------------!
subroutine init_random_seed()
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(:), allocatable :: seed
   integer                            :: i
   integer                            :: n
   integer                            :: ierr
   integer, dimension(8)              :: dt
   integer                            :: pid
   integer(kind=8)                    :: when
   !----- External functions. -------------------------------------------------------------!
   integer              , external    :: lcg
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Allocate seed.                                                                    !
   !---------------------------------------------------------------------------------------!
   call random_seed(size=n)
   allocate(seed(n))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     First try if the OS provides a random number generator.                           !
   !---------------------------------------------------------------------------------------!
   open (unit=72,file='/dev/urandom',access='stream',form='unformatted',action='read'      &
        ,status='old',iostat=ierr)
   select case (ierr)
   case (0)
      !----- System provides RN generator, use it. ----------------------------------------!
      read(unit=72) seed
      close(unit=72,status='keep')
      !------------------------------------------------------------------------------------!
   case default
      !------------------------------------------------------------------------------------!
      !     System doesn't have random number generator.  Fallback to XOR:ing the          !
      ! current time and pid.  The PID is useful in case one launches multiple             !
      ! instances of the same program in parallel.                                         !
      !------------------------------------------------------------------------------------!
      call system_clock(when)
      if (when == 0) then
         call date_and_time(values=dt)
         when = ( dt(1) - 1970 ) * 365_8 * 24_8 * 60 * 60 * 1000                           &
              +   dt(2)          * 31_8  * 24_8 * 60 * 60 * 1000                           &
              +   dt(3)                  * 24_8 * 60 * 60 * 1000                           &
              +   dt(5)                         * 60 * 60 * 1000                           &
              +   dt(6)                              * 60 * 1000                           &
              +   dt(7)                                   * 1000                           &
              +   dt(8)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Get Process ID.                                                               !
      !------------------------------------------------------------------------------------!
      pid  = getpid()
      when = ieor(when,int(pid,kind(when)))
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Set seeds.                                                                     !
      !------------------------------------------------------------------------------------!
      do i=1,n
         seed(i) = lcg(when)
      end do
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!

   call random_seed(put=seed)

   return
end subroutine init_random_seed
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This simple PRNG might not be good enough for real work, but is sufficient for      !
! seeding a better PRNG.                                                                   !
!                                                                                          !
!      This function has been borrowed from                                                !
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html                              !
!------------------------------------------------------------------------------------------!
integer function lcg(s)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer(kind=8), intent(in) :: s
   !----- Local variables. ----------------------------------------------------------------!
   integer(kind=8)             :: sloc
   !---------------------------------------------------------------------------------------!

   if (s == 0) then
      sloc = 104729_8
   else
      sloc = mod(s, 4294967296_8)
   end if
   sloc = mod(sloc * 279470273_8, 4294967291_8)
   lcg  = int(mod(sloc, int(huge(0), 8)), kind(0))

   return
end function lcg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine samples data with or without replacement (real numbers).             !
!------------------------------------------------------------------------------------------!
subroutine fsample(nxi,xi,nxo,xo,replacement)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nxi         ! Size of the x vector
   real(kind=4), dimension(nxi), intent(in)  :: xi          ! Vector to be sampled
   integer                     , intent(in)  :: nxo         ! Size of the output vector
   real(kind=4), dimension(nxo), intent(out) :: xo          ! Vector with the samples
   logical                     , intent(in)  :: replacement ! Sample with replacement?
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: i      ! Counter
   integer                                   :: o      ! Counter
   integer                                   :: s      ! Counter
   integer                                   :: nxs    ! # of points still selectable
   real(kind=4)                              :: rndnow ! Current random number
   real(kind=4), dimension(nxi)              :: xs     ! Input data (still selectable)
   real(kind=4)                              :: xph    ! Place holder
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Sanity check: nxo cannot be greater than nxi if it is sampling without           !
   ! replacement                                                                           !
   !---------------------------------------------------------------------------------------!
   if ((.not. replacement) .and. nxi < nxo) then
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      write(unit=*,fmt='(a)')       ' Invalid settings.  NXO cannot exceed NXI for '
      write(unit=*,fmt='(a)')       '    sampling without replacement.' 
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      write(unit=*,fmt='(a,1x,i6)') ' NXI         = ',nxi
      write(unit=*,fmt='(a,1x,i6)') ' NXO         = ',nxo
      write(unit=*,fmt='(a,1x,l1)') ' REPLACEMENT = ',replacement
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      call fatal_error('Invalid sampling size','fsample','random_utils.f90')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initially xs=xi and nxs=nxi.  This will remain true in the case of sampling      !
   ! with replacement.  If this is sampling without replacement, we will reorganise the    !
   ! vector and reduce nxs so selected items cannot be selected more than once.            !
   !---------------------------------------------------------------------------------------!
   nxs   = nxi
   xs(:) = xi(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Randomly select a point (bounded between 1 and nxs). -------------------------!
      call random_number(rndnow)
      s     = 1 + mod(floor(real(nxs)*rndnow),nxs)
      xo(o) = xs(s)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case this is a sampling without replacement, swap data between the index     !
      ! we have just selected and the index of the last selectable point.  Then we         !
      ! shrink the maximum number of selectable points so this point cannot be chosen      !
      ! again.                                                                             !
      !------------------------------------------------------------------------------------!
      if (.not. replacement) then
         xph     = xs(nxs)
         xs(nxs) = xs(  s)
         xs(  s) = xph
         nxs     = nxs - 1
      end if
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine fsample
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine samples data with or without replacement.                            !
!------------------------------------------------------------------------------------------!
subroutine isample(nxi,xi,nxo,xo,repl)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nxi         ! Size of the x vector
   integer     , dimension(nxi), intent(in)  :: xi          ! Vector to be sampled
   integer                     , intent(in)  :: nxo         ! Size of the output vector
   integer     , dimension(nxo), intent(out) :: xo          ! Vector with the samples
   logical                     , intent(in)  :: replacement ! Sample with replacement?
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: i      ! Counter
   integer                                   :: o      ! Counter
   integer                                   :: s      ! Counter
   integer                                   :: nxs    ! # of still selectable points
   real(kind=4)                              :: rndnow ! Current random number
   integer     , dimension(nxi)              :: xs     ! Input data (still selectable)
   integer                                   :: xph    ! Place holder
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Sanity check: nxo cannot be greater than nxi if it is sampling without           !
   ! replacement                                                                           !
   !---------------------------------------------------------------------------------------!
   if ((.not. replacement) .and. nxi < nxo) then
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      write(unit=*,fmt='(a)')       ' Invalid settings.  NXO cannot exceed NXI for '
      write(unit=*,fmt='(a)')       '    sampling without replacement.' 
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      write(unit=*,fmt='(a,1x,i6)') ' NXI         = ',nxi
      write(unit=*,fmt='(a,1x,i6)') ' NXO         = ',nxo
      write(unit=*,fmt='(a,1x,l1)') ' REPLACEMENT = ',replacement
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      call fatal_error('Invalid sampling size','isample','random_utils.f90')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initially xs=xi and nxs=nxi.  This will remain true in the case of sampling      !
   ! with replacement.  If this is sampling without replacement, we will reorganise the    !
   ! vector and reduce nxs so selected items cannot be selected more than once.            !
   !---------------------------------------------------------------------------------------!
   nxs   = nxi
   xs(:) = xi(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Randomly select a point (bounded between 1 and nxs). -------------------------!
      call random_number(rndnow)
      s     = 1 + mod(floor(real(nxs)*rndnow),nxs)
      xo(o) = xs(s)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case this is a sampling without replacement, swap data between the index     !
      ! we have just selected and the index of the last selectable point.  Then we         !
      ! shrink the maximum number of selectable points so this point cannot be chosen      !
      ! again.                                                                             !
      !------------------------------------------------------------------------------------!
      if (.not. replacement) then
         xph     = xs(nxs)
         xs(nxs) = xs(  s)
         xs(  s) = xph
         nxs     = nxs - 1
      end if
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine isample
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine samples data with or without replacement (double precision           !
! numbers).                                                                                !
!------------------------------------------------------------------------------------------!
subroutine fsample8(nxi,xi,nxo,xo,replacement)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nxi         ! Size of the x vector
   real(kind=8), dimension(nxi), intent(in)  :: xi          ! Vector to be sampled
   integer                     , intent(in)  :: nxo         ! Size of the output vector
   real(kind=8), dimension(nxo), intent(out) :: xo          ! Vector with the samples
   logical                     , intent(in)  :: replacement ! Sample with replacement?
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: i      ! Counter
   integer                                   :: o      ! Counter
   integer                                   :: s      ! Counter
   integer                                   :: nxs    ! # of points still selectable
   real(kind=4)                              :: rndnow ! Current random number
   real(kind=8), dimension(nxi)              :: xs     ! Input data (still selectable)
   real(kind=8)                              :: xph    ! Place holder
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Sanity check: nxo cannot be greater than nxi if it is sampling without           !
   ! replacement                                                                           !
   !---------------------------------------------------------------------------------------!
   if ((.not. replacement) .and. nxi < nxo) then
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      write(unit=*,fmt='(a)')       ' Invalid settings.  NXO cannot exceed NXI for '
      write(unit=*,fmt='(a)')       '    sampling without replacement.' 
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      write(unit=*,fmt='(a,1x,i6)') ' NXI         = ',nxi
      write(unit=*,fmt='(a,1x,i6)') ' NXO         = ',nxo
      write(unit=*,fmt='(a,1x,l1)') ' REPLACEMENT = ',replacement
      write(unit=*,fmt='(a)')       '--------------------------------------------------'
      call fatal_error('Invalid sampling size','fsample8','random_utils.f90')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initially xs=xi and nxs=nxi.  This will remain true in the case of sampling      !
   ! with replacement.  If this is sampling without replacement, we will reorganise the    !
   ! vector and reduce nxs so selected items cannot be selected more than once.            !
   !---------------------------------------------------------------------------------------!
   nxs   = nxi
   xs(:) = xi(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Randomly select a point (bounded between 1 and nxs). -------------------------!
      call random_number(rndnow)
      s     = 1 + mod(floor(real(nxs)*rndnow),nxs)
      xo(o) = xs(s)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case this is a sampling without replacement, swap data between the index     !
      ! we have just selected and the index of the last selectable point.  Then we         !
      ! shrink the maximum number of selectable points so this point cannot be chosen      !
      ! again.                                                                             !
      !------------------------------------------------------------------------------------!
      if (.not. replacement) then
         xph     = xs(nxs)
         xs(nxs) = xs(  s)
         xs(  s) = xph
         nxs     = nxs - 1
      end if
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine fsample8
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine generates a vector of random numbers bound by xmin and xmax,         !
! assuming uniform distribution.                                                           !
!------------------------------------------------------------------------------------------!
subroutine runif(nx,xlwr,xupr,x)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nx     ! Size of the x vector
   real(kind=4)                , intent(in)  :: xlwr   ! Lower bound
   real(kind=4)                , intent(in)  :: xupr   ! Upper bound
   real(kind=4), dimension(nx) , intent(out) :: x      ! Vector with the random numbers
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: o      ! Counter
   real(kind=4)                              :: rndnow ! Current random number
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Call random_number (bounded between 0-1). ------------------------------------!
      call random_number(rndnow)
      !------------------------------------------------------------------------------------!


      !------ Scale the random numbers. ---------------------------------------------------!
      x(o) = xlwr + (xupr - xlwr) * rndnow
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine runif
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine generates a vector of random numbers bound by xmin and xmax,         !
! assuming uniform distribution (double precision).                                        !
!------------------------------------------------------------------------------------------!
subroutine runif8(nx,xlwr,xupr,x)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nx     ! Size of the x vector
   real(kind=8)                , intent(in)  :: xlwr   ! Lower bound
   real(kind=8)                , intent(in)  :: xupr   ! Upper bound
   real(kind=8), dimension(nx) , intent(out) :: x      ! Vector with the random numbers
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: o      ! Counter
   real(kind=4)                              :: rndnow ! Current random number
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Call random_number (bounded between 0-1). ------------------------------------!
      call random_number(rndnow)
      !------------------------------------------------------------------------------------!


      !------ Scale the random numbers. ---------------------------------------------------!
      x(o) = xlwr + (xupr - xlwr) * dble(rndnow)
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine runif8
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine generates a vector of random numbers following the normal            !
! distribution.                                                                            !
!------------------------------------------------------------------------------------------!
subroutine rnorm(nx,xmean,xsdev,x)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nx         ! Size of the x vector
   real(kind=4)                , intent(in)  :: xmean      ! Mean
   real(kind=4)                , intent(in)  :: xsdev      ! Standard deviation
   real(kind=4), dimension(nx) , intent(out) :: x          ! Vector with the random numbers
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: o          ! Counter
   real(kind=4)                              :: rndnow     ! Current random number
   real(kind=4)                              :: xnorm      ! Normalised number
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)                , external    :: cdf2normal ! Quantile finder
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Call random_number (bounded between 0-1). ------------------------------------!
      call random_number(rndnow)
      xnorm = cdf2normal(rndnow)
      !------------------------------------------------------------------------------------!


      !------ Scale the random numbers. ---------------------------------------------------!
      x(o) = xmean + xnorm * xsdev
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine rnorm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine generates a vector of random numbers following the normal            !
! distribution.  This may improve later, right now we use the single precision number.     !
!------------------------------------------------------------------------------------------!
subroutine rnorm8(nx,xmean,xsdev,x)
   implicit none

   !----- Required arguments. -------------------------------------------------------------!
   integer                     , intent(in)  :: nx         ! Size of the x vector
   real(kind=8)                , intent(in)  :: xmean      ! Mean
   real(kind=8)                , intent(in)  :: xsdev      ! Standard deviation
   real(kind=8), dimension(nx) , intent(out) :: x          ! Vector with the random numbers
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: o          ! Counter
   real(kind=4)                              :: rndnow     ! Current random number
   real(kind=4)                              :: xnorm      ! Normalised number
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)                , external    :: cdf2normal ! Quantile finder
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Populate the output vector.                                                       !
   !---------------------------------------------------------------------------------------!
   oloop: do o=1,nxo
      !----- Call random_number (bounded between 0-1). ------------------------------------!
      call random_number(rndnow)
      xnorm = cdf2normal(rndnow)
      !------------------------------------------------------------------------------------!


      !------ Scale the random numbers. ---------------------------------------------------!
      x(o) = xmean + dble(xnorm) * xsdev
      !------------------------------------------------------------------------------------!
   end do oloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine rnorm8
!==========================================================================================!
!==========================================================================================!
