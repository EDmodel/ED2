
program main

  ! Main:
  !   defines execution strategy (sequential or parallel)
  !   enrolls processes at defined execution
  !   parses command line arguments, looking for namelist file
  !   dispatches processes according to execution strategy,
  !   invoking master/slave processes or full model process
  !   destroy processes
  !

  implicit none

  ! command line arguments

  integer, parameter :: MAX_INPUT_ARGS=63  ! maximum # input arguments (includes MPI own args)
  integer :: numarg  ! actual # input args
  integer :: iargc   ! function to return numarg
  integer, parameter :: MAX_INPUT_ARG_LENGTH=256 ! maximum length of each input argument
  character(len=MAX_INPUT_ARG_LENGTH)   :: cargs(0:MAX_INPUT_ARGS)  ! args 
  character(len=2*MAX_INPUT_ARG_LENGTH) :: cargx ! scratch

  ! input namelist file name

  character(len=MAX_INPUT_ARG_LENGTH) :: name_name = 'ED2IN'

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

  ! isingle: force non mpi run
  !    0 iff normal run with potential MPI
  !    1 iff no MPI init is called
  integer :: isingle=0

  ! scratch

  integer :: n
  character(len=*), parameter :: h="**(main)**"
  character(len=8) :: c0, c1

  ! For MPI interface 
  integer :: ierr
  include 'mpif.h'

  ! Get input arguments (required by C interface of MPI_Init)

  numarg=iargc()
  if (numarg > MAX_INPUT_ARGS) then
     write(c0,"(i8)") numarg
     write(c1,"(i8)") MAX_INPUT_ARGS
     write(*,"(a)") h//"ERROR**: input argument list length ("//&
          trim(adjustl(c0))//")exceeds MAX_INPUT_ARGS ("//&
          trim(adjustl(c1))//")"
     stop
  end if
  do n=0,numarg
     call ugetarg(n,cargx)
     if (len_trim(cargx) > MAX_INPUT_ARG_LENGTH) then
        write(c0,"(i8)") len_trim(cargx)
        write(c1,"(i8)") MAX_INPUT_ARG_LENGTH
        write(*,"(a)") h//"ERROR**: input argument data length ("//&
             trim(adjustl(c0))//")exceeds MAX_INPUT_ARG_LENGTH ("//&
             trim(adjustl(c1))//")"
        stop
     end if
     cargs(n)=trim(cargx)//char(0)
  enddo

  ! find out if sequential or MPI run; if MPI run, enroll this process.
  ! if sequential execution, machnum and machsize return untouched (both zero);
  ! if MPI execution, machnum returns process rank and machsize process size;

  numarg=numarg+1
  do n = 1, numarg
      if (cargs(n)(1:2) == '-s') then
         isingle=1
      end if
  end do

  if (isingle == 0) then
     call MPI_Init(ierr)                                              
     call MPI_Comm_rank(MPI_COMM_WORLD,machnum,ierr)                  
     call MPI_Comm_size(MPI_COMM_WORLD,machsize,ierr)                 
  else
     machnum=0
     machsize=1
  end if
  write (*,'(a)')       '+--- Parallel info: -------------------------------------+'
  write (*,'(a,1x,i6)') '+  - Machnum  =',machnum
  write (*,'(a,1x,i6)') '+  - Machsize =',machsize
  write (*,'(a)')       '+--------------------------------------------------------+'
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

  ! Read the namelist and initialize the variables in the nodes if needed 
  if (icall == 0) then
     call ed_1st_master(ipara,machsize,nslaves,machnum,name_name)
  else
     call ed_1st_node(1)
  endif

  ! Calling Main driver: it allocates the structures, initializes the variables,
  ! calls the timestep driver, deals with I/O, cooks, does the laundry and irons.
  ! The stand-alone driver tells the master node that it actually has to get a job
  ! and do something with its life.  So the driver is passed a zero here, which
  ! tells the node with mynum = nnodetot-1 that the master is next in line for
  ! sequencing.
  call ed_driver()

  ! finishes execution

  if (ipara == 1) then
     call MPI_Finalize(ierr)
  end if
  if (icall == 0) then
     write(*,"(a)") ' ------ ED execution ends ------'
  end if
  stop
end program main
