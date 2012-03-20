subroutine abort_run(reason,subr,file)
!------------------------------------------------------------------------------------------!
!   Subroutine based on RAMS, just to output error messages and halts the execution        !
! properly. You should always use this one, since it checks whether the run is running in  !
! parallel. If so, it will use MPI_Abort rather than stop, so it will exit rather than     !
! being frozen.                                                                            !
!------------------------------------------------------------------------------------------!
   use node_mod, only: nmachs,mynum
   implicit none
   character(len=*), intent(in) :: reason,subr,file

   include 'mpif.h'
  
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') '                     !!! FATAL ERROR !!!                          '
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   if (nmachs > 0 .and. mynum > 0) then
      write(unit=*,fmt='(a,1x,i5,a)') ' On node: ',mynum,':'
   elseif (nmachs > 0) then
      write(unit=*,fmt='(a)')         ' On the master node:'
   end if
   write(unit=*,fmt='(a,1x,a)')    '    ---> File:       ',trim(file)
   write(unit=*,fmt='(a,1x,a)')    '    ---> Subroutine: ',trim(subr)
   write(unit=*,fmt='(a,1x,a)')    '    ---> Reason:     ',trim(reason)
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   write(unit=*,fmt='(a)') ' BRAMS execution halts (see previous error message)...'
   write(unit=*,fmt='(a)') '------------------------------------------------------------------'
   if (nmachs > 0) call MPI_Abort(MPI_COMM_WORLD, 1)
   stop 'abort_run'
end subroutine abort_run
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine opspec_mess(reason,opssub)
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
end subroutine opspec_mess
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine brams_fail_whale()

   implicit none

   write(unit=*,fmt='(a)') ''
   write(unit=*,fmt='(a)') ''
   write(unit=*,fmt='(a)') ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '    ___ _  _ ____    ____ ____ _ _       _ _ _ _  _ ____ _    ____         '
   write(unit=*,fmt='(a)') '     |  |__| |___    |___ |__| | |       | | | |__| |__| |    |___         '
   write(unit=*,fmt='(a)') '     |  |  | |___    |    |  | | |___    |_|_| |  | |  | |___ |___         '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '       _  _ ____ ____    ____ ____ ____ ____ _  _ ____ ___                 '
   write(unit=*,fmt='(a)') '       |__| |__| [__     |    |__/ |__| [__  |__| |___ |  \                '
   write(unit=*,fmt='(a)') '       |  | |  | ___]    |___ |  \ |  | ___] |  | |___ |__/                '
   write(unit=*,fmt='(a)') '                                                                           ' 
   write(unit=*,fmt='(a)') '       _ _  _ ___ ____    _   _ ____ _  _ ____    ____ _ _  _              '
   write(unit=*,fmt='(a)') '       | |\ |  |  |  |     \_/  |  | |  | |__/    [__  | |\/|              '
   write(unit=*,fmt='(a)') '       | | \|  |  |__|      |   |__| |__| |  \    ___] | |  |              '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '                                             .+shhhhhhhhhhyso/-`           '
   write(unit=*,fmt='(a)') '             `.-::///+oooooooooooo+/:.`     -hhhhhhhhhhhhhhhhhhhs+.        '
   write(unit=*,fmt='(a)') '        -/oyhhhhhhhhhhhhhhhhhhhhhhhhhhhyo:` -hhhhhhhhhhhhhhhhhhhhhhy/      '
   write(unit=*,fmt='(a)') '     -ohhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhs:-oyhhhhhhhhhhhhhhhhhhhhhy-    '
   write(unit=*,fmt='(a)') '   -yhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhy/ `-+yhhhhhhhhhhhhhhhhhhh+   '
   write(unit=*,fmt='(a)') '  +hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhy:   :yhhhhhhhhhhhhhhhhhho  '
   write(unit=*,fmt='(a)') ' ohhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhho`  `yhhhhhhhhhhhhhhhhhh+ '
   write(unit=*,fmt='(a)') '/hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhy-  .hhhhhhhhhhyshhhhhhh-'
   write(unit=*,fmt='(a)') 'yhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhs` +hhhhhhhh:  .hhhhhhs'
   write(unit=*,fmt='(a)') 'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh/`yhhhhhy.    shhhhhh'
   write(unit=*,fmt='(a)') 'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh+/+o+:     :hhhhhhh'
   write(unit=*,fmt='(a)') 'hhhhhhhhhhhhhhh::yhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhso////+ohhhhhhhhs'
   write(unit=*,fmt='(a)') 'yhhhhhhhhhhhhhh+/yhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh.'
   write(unit=*,fmt='(a)') '/hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh: '
   write(unit=*,fmt='(a)') ' ohhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhs.  '
   write(unit=*,fmt='(a)') '  /hhhhhhhhhhhoohhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhho-    '
   write(unit=*,fmt='(a)') '   `:/+++/////shhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhsosyyyys+:.       '
   write(unit=*,fmt='(a)') '    -hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh::+/              '
   write(unit=*,fmt='(a)') '     :/++++yhhhhhhhhhhs++oyhhhhhhhhhhs+++ohhhhhhs+/ohhhhhhhy:              '
   write(unit=*,fmt='(a)') '           /shhhhhs+-      `:+oso+/.       .::-     `:///-.                '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') '                                                                           '
   write(unit=*,fmt='(a)') ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'

   return
end subroutine brams_fail_whale
!==========================================================================================!
!==========================================================================================!
