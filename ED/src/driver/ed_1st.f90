!==========================================================================================!
!==========================================================================================!
!    First subroutine. This subroutine is responsible for setting up the run so each node  !
! can work independently (or almost). This is the only subroutine that is called           !
! differently whether the node is a master or a slave (truly speaking, the model almost    !
! doesn't have master and slaves, except for reading the namelist met file header and to   !
! assign the polygons to each node. So ed_1st is the only subroutine that is actually      !
! called differently whether it is a "head" or "slave".                                    !
!                                                                                          !
! 1. ed_1st_master: run only by the master, it reads the namelist, the met drivers, and    !
!    splits the polygons accross the nodes, keeping some to itself. Then it sends the      !
!    namelist/met driver information to the slave nodes.                                   !
! 2. ed_1st_node: run only by the slaves, it gets the information from the master so the   !
!    modules can be properly initialized.                                                  !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine ed_1st_master (ipara, nnodestotal,nslaves, headnode_num, name_name)

   use ed_para_coms , only : iparallel          & ! intent(inout)
                           , machsize           ! ! intent(in)
   use ed_misc_coms , only : runtype            & ! intent(inout)
                           , iyeara             & ! intent(inout)
                           , imontha            & ! intent(inout)
                           , idatea             & ! intent(inout)
                           , itimea               ! intent(inout)
   use ed_state_vars, only : allocate_edglobals & ! subroutine
                           , filltab_alltypes   ! ! subroutine

   implicit none

   !----- Pre-compiled variables from MPI. ------------------------------------------------!
   include 'mpif.h'
   !----- Arguments. ----------------------------------------------------------------------!
   integer         , intent(in) :: ipara        ! 0 if sequential run; 1 if parallel run
   integer         , intent(in) :: nnodestotal  ! total number of nodes on any run
   integer         , intent(in) :: nslaves      ! number of slaves on a parallel run
   integer         , intent(in) :: headnode_num ! this process rank on a parallel run
   character(len=*), intent(in) :: name_name    ! namelist file name
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: nndtflg
   integer                      :: ierr
   real                         :: w1
   real                         :: wtime_start 
   !----- Local parameters, this sub-routine shan't ever be called by coupled runs. -------!
   logical         , parameter  :: masterworks = .true. ! Master should solve polygons
   logical         , parameter  :: standalone  = .true. ! ED will be run on its own.
   !----- External functions. -------------------------------------------------------------!
   real            , external   :: walltime             ! wall time
   !---------------------------------------------------------------------------------------!
   



   !----- Initialise wall time so we know how long it takes to setup the run. -------------!
   wtime_start = walltime(0.)
   w1          = walltime(wtime_start)
   
   !----- Save parallel flag (0 - serial, 1 - parallel). ----------------------------------!
   iparallel=ipara
   
   !----- Setup number of machines. -------------------------------------------------------!
   machsize = nnodestotal

   !----- Read the namelist file. ---------------------------------------------------------!
   write (unit=*,fmt='(a)') 'Reading namelist information'
   call read_nl(trim(name_name))

   write (unit=*,fmt='(a)') 'Copying namelist'
   call copy_nl('ALL_CASES')

   if (runtype == 'HISTORY') then
     call copy_nl('HISTORY')
   else
     call copy_nl('NOT_HISTORY')
   end if

   !---------------------------------------------------------------------------------------!
   !    Now that the namelist is loaded, I check whether all settings provided by the user !
   ! make sense or not.  This is done by the master only...                                !
   !---------------------------------------------------------------------------------------!
   call ed_opspec_grid()
   if (iparallel == 1) call ed_opspec_par()
   call ed_opspec_times()
   call ed_opspec_misc()

   !----- Read the met_driver namelist. ---------------------------------------------------!
   call read_met_driver_head()
   
   !----- Read the soil depth database. ---------------------------------------------------!
   call read_soil_depth()

   !----- Print the banner. ---------------------------------------------------------------!
   write(unit=*,fmt='(a)') '+------------------------------------------------------------+'
   write(unit=*,fmt='(a)') '|           Ecosystem Demography Model, version 2.2           '
   write(unit=*,fmt='(a)') '+------------------------------------------------------------+'
   write(unit=*,fmt='(a)') '|  Input namelist filename is '//trim(name_name)
   write(unit=*,fmt='(a)') '|'
   if (ipara == 0) then
      write(unit=*,fmt='(a)') '|  Single process execution on '//trim(runtype)//' run.'
   else
      write(unit=*,fmt='(a,1x,i4,3(1x,a))') '|  Parallel execution with master and'        &
                                           ,nslaves,'slave processes on',trim(runtype)     &
                                           ,'run.'
      if (iparallel == 0) then
         write(unit=*,fmt='(3(a,1x))')      '|  Converted into a sequential run for'       &
                                           ,trim(runtype),'run.'
      end if
   end if
   write(unit=*,fmt='(a)') '+------------------------------------------------------------+'
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now we are going to send the information to the nodes.  Even if the run is serial, !
   ! we still need to go through some of this routines to fill ed_node_coms and work load  !
   ! variables, because the master also integrates some polygons.  Since the integration   !
   ! subroutines do not care about whether the it is a head or a slave node (and they      !
   ! shouldn't), we must fill ed_node_coms in the head node as well.  For simplicity, the  !
   ! master is the last slave, so mynum(head)=nslaves+1.                                   !
   !---------------------------------------------------------------------------------------!
   if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

   call ed_masterput_processid(nslaves,headnode_num,masterworks,iparallel)
   call ed_masterput_nl(iparallel)

   call ed_masterput_met_header(iparallel)

   !---------------------------------------------------------------------------------------!
   !     The following subroutine does node decomposition, but it also does the initial    !
   ! read-in of the land-sea mask and soil textural class.                                 !
   !---------------------------------------------------------------------------------------!
   call ed_node_decomp(1,standalone,masterworks)
   if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

   call ed_masterput_poly_dims(iparallel,masterworks)

   if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

   call ed_masterput_worklist_info(iparallel)

   if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

   return
end subroutine ed_1st_master
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is called by all nodes other than the master, and sets up the         !
! polygons and parameters for the nodes.  This sub-routine won't be called if this is a    !
! serial run.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine ed_1st_node(init)
   use ed_mem_alloc, only : ed_memory_allocation ! ! subroutine
   implicit none
   !----- Pre-compiled variables from MPI. ------------------------------------------------!
   include 'mpif.h'
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in) :: init
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !----- Make sure the node is synchronised with all fellows. ----------------------------!
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   !----- Get the node information and the namelist. --------------------------------------!
   call ed_nodeget_processid(1)
   call ed_nodeget_nl

   !----- Get the meteorological driver header. -------------------------------------------!
   call ed_nodeget_met_header()
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   !----- Get the polygon dimensions and allocate structures. -----------------------------!
   call ed_nodeget_poly_dims()
   call ed_memory_allocation(2)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   !----- Get the work load structures. ---------------------------------------------------!
   call ed_nodeget_worklist_info()
   call MPI_Barrier(MPI_COMM_WORLD,ierr)

   return
end subroutine ed_1st_node
!==========================================================================================!
!==========================================================================================!
