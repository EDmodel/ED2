subroutine fatal_error(reason,subr,file)
!------------------------------------------------------------------------------------------!
!   Subroutine based on RAMS, just to output error messages and halts the execution        !
! properly. You should alway use this one, since it checks whether the run is running in   !
! parallel. If so, it will use MPI_Abort rather than stop, so it will exit rather than     !
!being frozen.                                                                             !
!------------------------------------------------------------------------------------------!
   use ed_node_coms   , only : nnodetot       & ! intent(in)
                             , mynum          ! ! intent(in)
   implicit none
   character(len=*), intent(in) :: reason
   character(len=*), intent(in) :: subr,file

   include 'mpif.h'
  
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') '                     !!! FATAL ERROR !!!                          '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   if (nnodetot > 1 .and. mynum /= nnodetot) then
      write(unit=*,fmt='(a,1x,i5,a)') ' On node: ',mynum,':'
   elseif (nnodetot > 1) then
      write(unit=*,fmt='(a)')         ' On the master node:'
   end if
   write (unit=*,fmt='(a,1x,a)')    '    ---> File:       ',trim(file)
   write (unit=*,fmt='(a,1x,a)')    '    ---> Subroutine: ',trim(subr)
   write (unit=*,fmt='(a,1x,a)')                      '    ---> Reason:     ',trim(reason)
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') ' ED execution halts (see previous error message)...'
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'

   !---------------------------------------------------------------------------------------!
   !     Remind the user of deprecated ED2IN choices...                                    !
   !---------------------------------------------------------------------------------------!
   if (nnodetot > 1) call MPI_Abort(MPI_COMM_WORLD, 1)
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

   include 'mpif.h'
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






!==========================================================================================!
!==========================================================================================!
subroutine warning(reason,subr,file)
!------------------------------------------------------------------------------------------!
!  Warning message, does not exit                                                          !
!------------------------------------------------------------------------------------------!

   use ed_node_coms, only: nnodetot,mynum
   implicit none
   character(len=*), intent(in) :: reason
   character(len=*), intent(in) :: subr,file

   include 'mpif.h'
  
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') '                       !!! WARNING !!!                            '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   if (nnodetot > 1 .and. mynum /= nnodetot) then
      write(unit=*,fmt='(a,1x,i5,a)') ' On node: ',mynum,':'
   elseif (nnodetot > 1) then
      write(unit=*,fmt='(a)')         ' On the master node:'
   end if
   ! Although it is optional, it should always be present 
   write(unit=*,fmt='(a,1x,a)')    '    ---> File:       ',trim(file)
   write(unit=*,fmt='(a,1x,a)')    '    ---> Subroutine: ',trim(subr)
   write (unit=*,fmt='(a,1x,a)')   '    ---> Reason:     ',trim(reason)
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
 end subroutine warning
!==========================================================================================!
!==========================================================================================!
