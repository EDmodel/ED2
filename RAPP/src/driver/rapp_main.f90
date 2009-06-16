!==========================================================================================!
!==========================================================================================!
!  RAPP - The ED Re-Analysis Pre-Processor.                                                !
!                                                                                          !
!     This program converts the reanalysis data to the standard ED meteorological driver   !
! format.  It currently supports only NCEP reanalysis , but it can be extended to other    !
! meteorological data sets.                                                                ! 
!                                                                                          !
!------------------------------------------------------------------------------------------!
program rapp_main
   use mod_maxdims, only : &
           maxstr          ! ! Maximum string length
   implicit none
   
   !---------------------------------------------------------------------------------------!
   !   Variable declaration section:                                                       !
   !---------------------------------------------------------------------------------------!
   integer                              :: numarg               ! # of arguments
   integer                              :: arg                  ! Argument counter
   character(len=maxstr), dimension(2)  :: arguments            ! List of arguments
   character(len=maxstr)                :: rapp_in              ! Namelist name
   character(len=12)                    :: c0                   ! Just to print the banner
   !
   real                                 :: w1,w2,wtime_start    ! For time calculation
   real, external                       :: walltime             ! To compute runtime
   !---------------------------------------------------------------------------------------!
   

   !---------------------------------------------------------------------------------------!
   ! 1.  Reading arguments:                                                                !
   !---------------------------------------------------------------------------------------!
   numarg = iargc()

   !---------------------------------------------------------------------------------------!
   ! 2.  Checking arguments                                                                !
   !---------------------------------------------------------------------------------------!
   select case (numarg) 
   case (0) ! No arguments provided, use default namelist name
      rapp_in='RAPP_IN'
      write (unit=*,fmt='(102a)') ('-',arg=1,102)
      write (unit=*,fmt='(a,1x,a,1x,a)')                                                     &
          ' No arguments, using the default namelist :',trim(rapp_in),'...'
      write (unit=*,fmt='(102a)') ('-',arg=1,102)
      write (unit=*,fmt=*)
   case (2) ! Correct # of arguments, checking whether they make sense or not.
      do arg=1,numarg
         call ugetarg(arg,arguments(arg))
      end do
      if (trim(arguments(1)) == '-f') then
         rapp_in=trim(arguments(2))
         write (unit=*,fmt='(102a)') ('-',arg=1,102)
         write (unit=*,fmt='(a,1x,a,1x,a)')                                                     &
             ' Using the user-defined namelist :',trim(rapp_in),'...'
         write (unit=*,fmt='(102a)') ('-',arg=1,102)
         write (unit=*,fmt=*)
      else
         call fatal_error('Invalid syntax! Usage: ./(rapp_executable) [-f (rapp_nlist)]' &
                         ,'rapp_main','rapp_main.f90')
      end if
   case default
      call fatal_error('Invalid syntax! Usage: ./(rapp_executable) [-f (rapp_nlist)]' &
                      ,'rapp_main','rapp_main.f90')
   
   end select

   !----- Initializing time ---------------------------------------------------------------! 
   write (unit=*,fmt='(a)') 'RAPP execution begins...'
   wtime_start=walltime(0.)
   w1=walltime(wtime_start)

   !---------------------------------------------------------------------------------------!
   ! 3.  Loading the namelist and storing the information at the appropriate modules.      !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)') 'Loading namelist...'
   call load_namelist(rapp_in)
   
   
   !---------------------------------------------------------------------------------------!
   ! 4.  Retrieving the files available, in case we are using RAMS/BRAMS.                  !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)') 'Calling the main driver...'
   call rapp_driver()

   !----- Getting the final time and printing the final time banner -----------------------!
   w2=walltime(wtime_start)
   write(c0,"(f12.2)") w2
   write(unit=*,fmt='(/,a,/)') ' === File conversion ends; Total elapsed time='//&
                               trim(adjustl(c0))//' ==='
  

   write (unit=*,fmt=*)

   stop '****** RAPP execution ends ******'
end program rapp_main
!==========================================================================================!
!==========================================================================================!
