!==========================================================================================!
!==========================================================================================!
!                    ****** Ecosystem Demography Model -- ED-2.2 ******                    !
!------------------------------------------------------------------------------------------!
!                                                                                          !
! Main program:                                                                            !
!   - Defines execution strategy (sequential or parallel)                                  !
!   - Enrolls processes at defined execution                                               !
!   - Parses command line arguments, looking for namelist file                             !
!   - Dispatches processes according to execution strategy,                                !
!     invoking master/slave processes or full model process                                !
!   - Destroy processes                                                                    !
!                                                                                          !
!------------------------------------------------------------------------------------------!
program main
   !$ use omp_lib
   implicit none

   !---------------------------------------------------------------------------------------!
   !      Local constants.                                                                 !
   !---------------------------------------------------------------------------------------!
   !----- Maximum number of input arguments, including MPI own arguments. -----------------!
   integer, parameter                    :: max_input_args       = 63
   !----- Maximum length of each input argument. ------------------------------------------!
   integer, parameter                    :: max_input_arg_length = 256
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: numarg                  ! actual # input args
   character(len=max_input_arg_length)   :: cargs(0:max_input_args) ! args 
   character(len=2*max_input_arg_length) :: cargx                   ! scratch
   character(len=max_input_arg_length)   :: name_name
   integer                               :: machsize
   integer                               :: machnum
   integer                               :: ipara
   integer                               :: icall
   integer                               :: nslaves
   integer                               :: isingle
   integer                               :: n
   integer                               :: ierr
   !------ Intrinsic function to return number of arguments (numarg). ---------------------!
   integer                               :: iargc
   !------ MPI interface. -----------------------------------------------------------------!
#if defined(RAMS_MPI)
   !----- OMP information. ----------------------------------------------------------------!
   integer                               :: max_threads      !<= omp_get_max_threads()
   integer                               :: num_procs        !<= omp_get_num_procs()
   integer                               :: thread
   integer                               :: cpu
   integer, dimension(64)                :: thread_use
   integer, dimension(64)                :: cpu_use
   integer, external                     :: findmycpu
   !---------------------------------------------------------------------------------------!
   include 'mpif.h'
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Default settings, which may change depending on the arguments.                   !
   !---------------------------------------------------------------------------------------!
      !----- Namelist. --------------------------------------------------------------------!
      name_name = 'ED2IN'
      !----- Process rank and size (default is single process running full model). --------!
      machsize = 0
      machnum  = 0
      !------------------------------------------------------------------------------------!
      ! IPARA: execution strategy; default single process                                  !
      !    0 iff single process (no MPI run or MPI run with a single process)              !
      !    1 iff master-slave processes (MPI run with more than one process)               !
      !------------------------------------------------------------------------------------!
      ipara    = 0
      !------------------------------------------------------------------------------------!
      ! ICALL: this process function on execution strategy; default full model
      !    0 iff full model (no MPI run) or master on MPI run
      !    1 iff slave on MPI run
      !------------------------------------------------------------------------------------!
      icall    = 0

      !------------------------------------------------------------------------------------!
      ! ISINGLE: force non mpi run                                                         !
      !    0 iff normal run with potential MPI                                             !
      !    1 iff no MPI init is called                                                     !
      !                                                                                    !
      ! MLO: Not sure if this will bring conflicts.  I implemented the non-MPI option as   !
      !      a preprocessor feature, and this skips the MPI commands everywhere in the     !
      !      code.  I tried to reconcile this to the best of my knowledge, it would be     !
      !      good if someone who knows MPI better than me could check.                     !
      !------------------------------------------------------------------------------------!
      isingle  = 0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Summary of execution strategy and process function:                                !
      !           ipara=0        ipara=1                                                   !
      ! icall=0   full model     master on master/slave run                                !
      ! icall=1   impossible     slave  on master/slave run                                !
      !------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Get input arguments.                                                              !
   !---------------------------------------------------------------------------------------!
   numarg=iargc()
   if (numarg > max_input_args) then
      write(unit=*,fmt='(a)')       '-----------------------------------------------------'
      write(unit=*,fmt='(a,1x,i6)') ' NUMARG         = ',numarg
      write(unit=*,fmt='(a,1x,i6)') ' MAX_INPUT_ARGS = ',max_input_args
      write(unit=*,fmt='(a)')       '-----------------------------------------------------'
      call fatal_error('Input argument list length (NUMARG) exceeds MAX_INPUT_ARGS!'       &
                      ,'main','edmain.F90')
   end if
   do n=0,numarg
      call ugetarg(n,cargx)
      if (len_trim(cargx) > max_input_arg_length) then
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a,1x,i6)') ' ARGUMENT NUMBER      = ',n
         write(unit=*,fmt='(a,1x,a)' ) ' ARGUMENT             = ',cargx
         write(unit=*,fmt='(a,1x,i6)') ' LEN_TRIM(ARGUMENT)   = ',len_trim(cargx)
         write(unit=*,fmt='(a,1x,i6)') ' MAX_INPUT_ARG_LENGTH = ',max_input_arg_length
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         call fatal_error('Input argument length exceeds MAX_INPUT_ARG_LENGTH!'            &
                         ,'main','edmain.F90')
      end if
      cargs(n)=trim(cargx)//char(0)
   end do
   numarg=numarg+1
   !---------------------------------------------------------------------------------------!



#if defined(RAMS_MPI)
   !---------------------------------------------------------------------------------------!
   !      Find out if sequential or MPI run; if MPI run, enroll this process.  If          !
   ! sequential execution, machnum and machsize return untouched (both zero); if MPI       !
   ! execution, machnum returns process rank and machsize process size.                    !
   !---------------------------------------------------------------------------------------!
   do n = 1, numarg
       if (cargs(n)(1:2) == '-s') then
          isingle=1
       end if
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Decide whether to call MPI.                                                       !
   !---------------------------------------------------------------------------------------!
   select case (isingle)
   case (0)
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,machnum,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,machsize,ierr)
   case default
      machnum  = 0
      machsize = 1
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Check OMP thread and processor use and availability.                                  !
   !                                                                                       !
   ! Note: One could use omp_get_num_threads() in loop, but that would depend on how many  !
   ! threads were open at the time of its call.                                            !
   !---------------------------------------------------------------------------------------!
   max_threads   = 1
   num_procs     = 1
   thread        = 1
   cpu           = 1
   thread_use(:) = 0
   cpu_use(:)    = 0

   !$ max_threads = omp_get_max_threads()
   !$ num_procs   = omp_get_num_procs()

   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(thread,cpu)
   do n = 1,max_threads
     !$ thread = omp_get_thread_num() + 1
     !$ cpu    = findmycpu() + 1

     thread_use(thread) = 1
     cpu_use(cpu)       = 1
   end do
   !$OMP END PARALLEL DO
   !---------------------------------------------------------------------------------------!

   write (*,'(a)')       '+---------------- MPI parallel info: --------------------+'
   write (*,'(a,1x,i6)') '+  - Machnum  =',machnum
   write (*,'(a,1x,i6)') '+  - Machsize =',machsize
   write (*,'(a)')       '+---------------- OMP parallel info: --------------------+'
   write (*,'(a,1x,i6)') '+  - thread  use: ', sum(thread_use)
   write (*,'(a,1x,i6)') '+  - threads max: ', max_threads
   write (*,'(a,1x,i6)') '+  - cpu     use: ', sum(cpu_use)
   write (*,'(a,1x,i6)') '+  - cpus    max: ', num_procs
   write (*,'(a)')       '+  Note: Max vals are for node, not sockets.'
   write (*,'(a)')       '+--------------------------------------------------------+'
#else


   !---------------------------------------------------------------------------------------!
   !   Set dummy values for all OMP variables.                                             !
   !---------------------------------------------------------------------------------------!
   isingle       = 1
   machnum       = 0
   machsize      = 1
   max_threads   = 1
   num_procs     = 1
   thread        = 1
   cpu           = 1
   thread_use(:) = 0
   cpu_use(:)    = 0
   do n = 1,max_threads
     thread_use(thread) = 1
     cpu_use(cpu)       = 1
   end do
   !---------------------------------------------------------------------------------------!

#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    If this is MPI run, define master or slave process, otherwise, keep default        !
   ! (single process does full model).                                                     !
   !---------------------------------------------------------------------------------------!
   if (machsize > 1) then
      ipara = 1
      if (machnum /= 0) then
         icall=1
      end if
   end if
   !---------------------------------------------------------------------------------------!


   !----- Master process gets number of slaves and sets process ID. -----------------------!
   if (icall == 0) then
      nslaves=machsize-1
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Master process parse command line arguments looking for "-f <namelist filename>". !
   !---------------------------------------------------------------------------------------!
   if (icall == 0) then
      do n = 1, numarg
         if (cargs(n)(1:2) == '-f') then
            name_name = cargs(n+1)(1:len_trim(cargs(n+1))-1)
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- Read the namelist and initialize the variables in the nodes if needed. ----------!
   if (icall == 0) then
      call ed_1st_master(ipara,machsize,nslaves,machnum,name_name)
   else
      call ed_1st_node(1)
   endif
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Call the main driver: it allocates the structures, initializes the variables,     !
   ! calls the timestep driver, deals with I/O, cooks, does the laundry and irons.  The    !
   ! stand-alone driver tells the master node that it actually has to get a job and do     !
   ! something with its life.  So the driver is passed a zero here, which tells the node   !
   ! with mynum = nnodetot-1 that the master is next in line for sequencing.               !
   !---------------------------------------------------------------------------------------!
   call ed_driver()
   !---------------------------------------------------------------------------------------!



   !----- Finishes execution. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (ipara == 1) then
      call MPI_Finalize(ierr)
   end if
#endif
   if (icall == 0) then
      write(unit=*,fmt='(a)') ' ------ ED-2.2 execution ends ------'
   end if
   !---------------------------------------------------------------------------------------!

   stop
end program main
!==========================================================================================!
!==========================================================================================!
