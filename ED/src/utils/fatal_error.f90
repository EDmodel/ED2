subroutine fatal_error(reason,subr,file)
!------------------------------------------------------------------------------------------!
!   Subroutine based on RAMS, just to output error messages and halts the execution        !
! properly. You should alway use this one, since it checks whether the run is running in   !
! parallel. If so, it will use MPI_Abort rather than stop, so it will exit rather than     !
!being frozen.                                                                             !
!------------------------------------------------------------------------------------------!
   use ed_node_coms   , only : nnodetot       & ! intent(in)
                             , mynum          ! ! intent(in)
   use canopy_air_coms, only : icanturb       ! ! intent(in)
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
   if (icanturb == -1) then
      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(a)') '------------------------------------------------------------'
      write(unit=*,fmt='(a)') '     I TOLD YOU NOT TO RUN WITH THE old ED-2.0 canopy       '
      write(unit=*,fmt='(a)') ' turbulence structure.  It''s always like that, I warn,     '
      write(unit=*,fmt='(a)') ' I try to convince it is a very bad idea to run me with     '
      write(unit=*,fmt='(a)') ' deprecated options, for what?  In the end, they never      '
      write(unit=*,fmt='(a)') ' listen to me, instead they ignore my advices and keep      '
      write(unit=*,fmt='(a)') ' insisting on these bad options.  Oh well, but at           '
      write(unit=*,fmt='(a)') ' least this time I can say that I told you, I told you,     '
      write(unit=*,fmt='(a)') ' I told you, I told you, I told you...  But hey, look,      '
      write(unit=*,fmt='(a)') ' don''t take this message personally, I''m still a nice     '
      write(unit=*,fmt='(a)') ' model and I promise I will try my best to please you with  '
      write(unit=*,fmt='(a)') ' good and exciting results should you use another ICANTURB  '
      write(unit=*,fmt='(a)') ' like 0 or 2 (1 is not that good either).                   '
      write(unit=*,fmt='(a)') '------------------------------------------------------------'
   end if
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
   character(len=*), intent(in), optional   :: subr,file

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
   if (present(file)) write(unit=*,fmt='(a,1x,a)')    '    ---> File:       ',trim(file)
   if (present(subr)) write(unit=*,fmt='(a,1x,a)')    '    ---> Subroutine: ',trim(subr)
   write (unit=*,fmt='(a,1x,a)')                      '    ---> Reason:     ',trim(reason)
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
 end subroutine warning
!==========================================================================================!
!==========================================================================================!
