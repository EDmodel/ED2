subroutine fatal_error(reason,subr,file)
!------------------------------------------------------------------------------------------!
!   Subroutine that outputs error messages and halts the execution.                        !
!------------------------------------------------------------------------------------------!
   implicit none
   character(len=*), intent(in) :: reason
   character(len=*), intent(in) :: subr,file
  
   write(unit=*,fmt='(a)') ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(a)') '---------------------------------------------------------------'
   write(unit=*,fmt='(a)') '                    !!! FATAL ERROR !!!                        '
   write(unit=*,fmt='(a)') '---------------------------------------------------------------'
   write(unit=*,fmt='(a,1x,a)')    '    ---> File:       ',trim(file)
   write(unit=*,fmt='(a,1x,a)')    '    ---> Subroutine: ',trim(subr)
   write (unit=*,fmt='(a,1x,a)')   '    ---> Reason:     ',trim(reason)
   write(unit=*,fmt='(a)') '---------------------------------------------------------------'
   write(unit=*,fmt='(a)') ' MEVI execution halts (see previous error message)...'
   write(unit=*,fmt='(a)') '---------------------------------------------------------------'
   stop 'fatal_error'
end subroutine fatal_error
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine opspec_fatal(reason,opssub)
!------------------------------------------------------------------------------------------!
!   This is just a first warning to be given in the standard output. Since the namelist    !
! may have more than one error, I list all the problems, then the model will stop.         !
!------------------------------------------------------------------------------------------!
   implicit none
   character(len=*), intent(in) :: reason,opssub

   write (unit=*,fmt='(a)')       ' '
   write (unit=*,fmt='(a)')       '----------------------------------------------------------------------------'
   write (unit=*,fmt='(3(a,1x))') '>>>> ',trim(opssub),' error! in your namelist!'
   write (unit=*,fmt='(a,1x,a)')  '    ---> Reason:     ',trim(reason)
   write (unit=*,fmt='(a)')       '----------------------------------------------------------------------------'
   write (unit=*,fmt='(a)')       ' '
   return
end subroutine opspec_fatal
!==========================================================================================!
!==========================================================================================!
