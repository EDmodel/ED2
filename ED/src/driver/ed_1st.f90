!------------------------------------------------------------------------------------------!
!                                                                                          !
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

  use ed_para_coms, only:   &
              iparallel,    &
              machsize
  
  use misc_coms, only:  &
       runtype,         &
       iyeara,          &
       imontha,         &
       idatea,          &
       itimea
  
  use ed_state_vars, only: &
       allocate_edglobals, &    ! implicitly interfaced subroutine
       filltab_alltypes

  use mem_sites, only: n_soi

  implicit none

  !   Input arguments

  integer, intent(in) :: ipara              ! 0 iff sequential run; 1 iff parallel run
  integer, intent(in) :: nnodestotal        ! total number of nodes on any run
  integer, intent(in) :: nslaves            ! number of slaves on a parallel run
  integer, intent(in) :: headnode_num       ! this process rank on a parallel run
  character(len=*), intent(in) :: name_name ! namelist file name
  logical, parameter  :: masterworks=.true.,standalone=.true.
  
  !   CPU Timing variables

  real :: w1,w2,w3,wtime_start  ! wall time
  real, external :: walltime    ! wall time
  real :: t1,t2                 ! cpu time

  character(len=12) :: c1
  integer :: i,j,ifm,nndtflg
  integer :: ierr
  
  integer :: id1,id2,jd1,jd2
  character(len=24) :: fmts1

  !   MPI header
  
  include 'mpif.h'

  wtime_start=walltime(0.)
  w1=walltime(wtime_start)
  
  !    Setup iparallel

  iparallel=ipara
  
  !   Setup number of machines
  machsize = nnodestotal

  !    Read the namelist file
  write (unit=*,fmt='(a)') 'Reading namelist information'
  call read_nl(trim(name_name))
  
  write (unit=*,fmt='(a)') 'Copying namelist'
  call copy_nl('ALL_CASES')

  if (runtype =='HISTORY') then
    call copy_nl('HISTORY')
  else
    call copy_nl('NOT_HISTORY')
  end if

  ! Now that the namelist is loaded, I check if for consistency. This is done by the master only...
  call ed_opspec_grid()
  if (iparallel == 1) call ed_opspec_par()
  call ed_opspec_times()
  call ed_opspec_misc()



  ! Read the met_driver namelist 
  call read_met_driver_head()
  
  ! Read the soil depth database
  call read_soil_depth()

  !    Print the banner.

  write(unit=*,fmt='(a)') '+-------------------------------------------------------------------------+'
  write(unit=*,fmt='(a)') '|                    Ecosystem Demography Model 2                          '
  write(unit=*,fmt='(a)') '+-------------------------------------------------------------------------+'
  write(unit=*,fmt='(a)') '|  Input namelist filename is '//trim(name_name)
  write(unit=*,fmt='(a)') '|'
  if (ipara == 0) then
     write(unit=*,fmt='(a)') '|  Single process execution on '//trim(runtype)//' run.'
  else
     write(unit=*,fmt='(a,1x,i4,3(1x,a))') '|  Parallel execution with master and' &
                                     ,nslaves,'slave processes on',trim(runtype),'run.'
     if (iparallel == 0) then
        write(unit=*,fmt='(3(a,1x))') '|  Converted into a sequential run for' &
                                ,trim(runtype),'run.'
     end if
  end if
  write(unit=*,fmt='(a)') '+-------------------------------------------------------------------------+'


!------------------------------------------------------------------------------------------!
!    Now I am going to send the information to the nodes. Even if the run is serial, I     !
! still need to go through some of this routines to fill some ed_node_coms structures,     !
! because the master is also a slave now, in the sense that it also integrates some        !
! polygons. Since the integration subroutines do not care about whether the it is a head   !
! or a slave node (and they shouldn't), we must fill ed_node_coms in the head node as      !
! well. For simplicity, the master is the last slave, so mynum(head)=nslaves+1.            !
!------------------------------------------------------------------------------------------!
  if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call ed_masterput_processid(nslaves,headnode_num,masterworks,iparallel)
  call ed_masterput_nl(iparallel)
  
  call ed_masterput_met_header(iparallel)

  call ed_node_decomp(1,standalone,masterworks)
  if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call ed_masterput_grid_dimens(iparallel)
  call ed_masterput_poly_dims(iparallel)
  call ed_dump_Domain_Decomposition(masterworks)
  
  if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call ed_masterput_gridded_info(iparallel)

  if (iparallel == 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  

  return
end subroutine ed_1st_master
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
subroutine ed_1st_node(init)


  implicit none

  integer, intent(in) :: init

  include 'mpif.h'

  integer :: ierr

  !          get all initialization info from the master process
  !          -------------------------------------------------------

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call ed_nodeget_processid(1)
  call ed_nodeget_nl
  call ed_nodeget_met_header()
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call ed_nodeget_grid_dimens()
  call ed_nodeget_poly_dims()
  call ed_mem_alloc(2)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  call ed_nodeget_gridded_info()
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  return

end subroutine ed_1st_node
