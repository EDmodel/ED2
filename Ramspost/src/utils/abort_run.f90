!==========================================================================================!
!==========================================================================================!
subroutine abort_run(reason,subr,file)
!------------------------------------------------------------------------------------------!
!   Subroutine based on RAMS, just to output error messages and halts the execution        !
! properly. You should always use this one, since it checks whether the run is running in  !
! parallel. If so, it will use MPI_Abort rather than stop, so it will exit rather than     !
! being frozen.                                                                            !
!------------------------------------------------------------------------------------------!
   implicit none
   character(len=*), intent(in) :: reason,subr,file
  
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') '                     !!! FATAL ERROR !!!                          '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a,1x,a)')    '    ---> File:       ',trim(file)
   write(unit=*,fmt='(a,1x,a)')    '    ---> Subroutine: ',trim(subr)
   write(unit=*,fmt='(a,1x,a)')    '    ---> Reason:     ',trim(reason)
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') ' BRAMS execution halts (see previous error message)...'
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   stop 'abort_run'
end subroutine abort_run
!==========================================================================================!
!==========================================================================================!
