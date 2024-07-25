!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

program main

  ! Main:
  !   defines execution strategy (sequential or parallel)
  !   enrolls processes at defined execution
  !   parses command line arguments, looking for namelist file
  !   dispatches processes according to execution strategy,
  !   invoking master/slave processes or full model process
  !   destroy processes
  !
#if defined(RAMS_MPI)
  use mpi
#endif

  implicit none

  ! command line arguments

  integer, parameter :: MAX_INPUT_ARGS=63  ! maximum # input arguments (includes MPI own args)
  integer :: numarg  ! actual # input args
  integer :: iargc   ! function to return numarg
  integer, parameter :: MAX_INPUT_ARG_LENGTH=256 ! maximum length of each input argument
  character(len=MAX_INPUT_ARG_LENGTH)   :: cargs(0:MAX_INPUT_ARGS)  ! args 
  character(len=2*MAX_INPUT_ARG_LENGTH) :: cargx ! scratch

  ! input namelist file name

  character(len=MAX_INPUT_ARG_LENGTH) :: name_name = 'RAMSIN'

  ! process rank and size (default single process running full model)

  integer :: machsize=0
  integer :: machnum=0

  ! ipara: execution strategy; default single process
  !    0 iff single process (no MPI run or MPI run with a single process)
  !    1 iff master-slave processes (MPI run with more than one process)
  integer :: ipara=0

  ! icall: this process function on execution strategy; default full model
  !    0 iff full model (no MPI run) or master on MPI run
  !    1 iff slave on MPI run
  integer :: icall=0

  ! summary of execution strategy and process function: 
  !           ipara=0        ipara=1
  ! icall=0   full model     master on master/slave run
  ! icall=1   impossible     slave  on master/slave run

  ! number of slave processes (master only!)

  integer :: nslaves

  ! scratch

  integer :: n
  character(len=*), parameter :: h="**(main)**"
  character(len=8) :: c0, c1

  ! For MPI interface 
  integer :: ierr
#if defined(RAMS_MPI)
  include 'interface.h'
#endif
  ! Get input arguments (required by C interface of MPI_Init)

  numarg=iargc()
  if (numarg > MAX_INPUT_ARGS) then
     write(c0,"(i8)") numarg
     write(c1,"(i8)") MAX_INPUT_ARGS
     write(*,"(a)") h//"ERROR**: input argument list length ("//&
          &trim(adjustl(c0))//")exceeds MAX_INPUT_ARGS ("//&
          &trim(adjustl(c1))//")"
     stop
  end if
  do n=0,numarg
     call ugetarg(n,cargx)
     if (len_trim(cargx) > MAX_INPUT_ARG_LENGTH) then
        write(c0,"(i8)") len_trim(cargx)
        write(c1,"(i8)") MAX_INPUT_ARG_LENGTH
        write(*,"(a)") h//"ERROR**: input argument data length ("//&
             &trim(adjustl(c0))//")exceeds MAX_INPUT_ARG_LENGTH ("//&
             &trim(adjustl(c1))//")"
        stop
     end if
     cargs(n)=trim(cargx)//char(0)
  enddo

  ! find out if sequential or MPI run; if MPI run, enroll this process.
  ! if sequential execution, machnum and machsize return untouched (both zero);
  ! if MPI execution, machnum returns process rank and machsize process size;

  numarg=numarg+1
#if defined(RAMS_MPI)
  call MPI_Init(ierr)                                              
  call MPI_Comm_rank(MPI_COMM_WORLD,machnum,ierr)                  
  call MPI_Comm_size(MPI_COMM_WORLD,machsize,ierr)                 
#else
   machnum=0
   machsize=1
#endif
  write (*,'(a)')       '«»«» Parallel info: «»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»'
  write (*,'(a,1x,i6)') '«» » Machnum  =',machnum
  write (*,'(a,1x,i6)') '«» » Machsize =',machsize
  write (*,'(a)')       '«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»«»'
  ! if MPI run, define master or slave process
  ! if sequential run, keep default (sigle process does full model)

  if (machsize > 1) then
     ipara = 1
     if (machnum /= 0) then
        icall=1
     end if
  endif

  ! master process gets number of slaves and sets process id

  if (icall == 0) then
     nslaves=machsize-1
  end if

  ! master process parse command line arguments looking for "-f <namelist filename>"

  if (icall == 0) then
     do n = 1, numarg
        if (cargs(n)(1:2) == '-f') then
           name_name = cargs(n+1)(1:len_trim(cargs(n+1))-1)
        end if
     end do
  end if

  ! dispatch processes

  if (icall == 0) then
     call rams_master (ipara, nslaves, machnum, name_name)
  else
     call rams_node()
  endif

  ! finishes execution

  if (ipara == 1) then
#if defined(RAMS_MPI)
     call MPI_Finalize(ierr)
#endif
  end if
  if (icall == 0) then
     write(*,"(a)") ' ****** BRAMS execution ends ******'
  end if
  stop
end program main
