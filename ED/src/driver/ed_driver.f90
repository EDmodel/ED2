!------------------------------------------------------------------------!
!                                                                        !
! Main subroutine for driving the initialization process for             !
! the Ecosystem Demography Model 2.  All compute nodes, including        !
! the node formerly known as master.                                     !
!                                                                        !
!------------------------------------------------------------------------!

subroutine ed_driver
  
  use grid_coms, only: &
       ngrids,          &
       time,            &
       timmax
  
  use ed_state_vars, only: &
       allocate_edglobals, &    ! implicitly interfaced subroutine
       filltab_alltypes, &
       edgrid_g
  
  use misc_coms, only: &
       iyeara,          &
       imontha,         &
       idatea,          &
       itimea,          &
       runtype
       
  use soil_coms, only: alloc_soilgrid
  use ed_node_coms , only: mynum,nnodetot,sendnum,recvnum

  implicit none

  real :: w1,w2,w3,wtime_start  ! wall time
  real, external :: walltime    ! wall time
  real :: t1,t2                 ! cpu time
  character(len=12) :: c0
  character(len=12) :: c1
  integer :: i,j,ifm,nndtflg,ipy

  integer :: id1,id2,jd1,jd2
  character(len=24) :: fmts1

  !   MPI header
  include 'mpif.h'

  integer                                   :: ierr
  integer :: igr
  integer :: ping 

  ping = 741776

  wtime_start=walltime(0.)
  w1=walltime(wtime_start)

  !---------------------------------------------------------------------------!
  ! STEP 1: Set the ED model parameters                                       !
  !---------------------------------------------------------------------------!
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ed_Ecosystem_Params...'
  call load_ed_ecosystem_params()
  
  !---------------------------------------------------------------------------!
  ! STEP 2: Overwrite the parameters in case a XML file is provided           !
  !---------------------------------------------------------------------------!

  ! THIS IS SHOULD ONLY BE TRUE FOR A STAND-ALONE RUN
  if (mynum == nnodetot-1) sendnum = 0

  if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,602,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
  
  call overwrite_with_xml_config(mynum)
  
  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,602,MPI_COMM_WORLD,ierr)
  if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  !---------------------------------------------------------------------------!
  ! STEP 3: Allocate soil grid arrays                                         !
  !---------------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Alloc_Soilgrid...'
  call alloc_soilgrid()
  
  !---------------------------------------------------------------------------------!
  ! STEP 4: Set some polygon-level basic information, such as lon/lat/soil texture  !
  !---------------------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Set_Polygon_Coordinates_Ar...'
  call set_polygon_coordinates_ar()
  
  !---------------------------------------------------------------------------!
  ! STEP 5: Initialize inherent soil and vegetation properties.               !
  !---------------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Sfcdata_ED...'
  call sfcdata_ed()
  
  
  if (trim(runtype) == 'HISTORY') then
       
 
     !-----------------------------------------------------------------------!
     ! STEP 6A: Initialize the model state as a replicate image of a previous
     !          state.
     !-----------------------------------------------------------------------!

     if (mynum == nnodetot-1) sendnum = 0
     
     if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,606,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  
     if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Init_Full_History_Restart...'
     call init_full_history_restart()
     
     if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,606,MPI_COMM_WORLD,ierr)
     if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     
  else
     
     !------------------------------------------------------------------------!
     ! STEP 6B: Initialize state properties of polygons/sites/patches/cohorts !
     !------------------------------------------------------------------------!
     
     if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ecosystem_State...'
     call load_ecosystem_state()

  end if

  !--------------------------------------------------------------------------------!
  ! STEP 7: Initialize hydrology related variables                                 !
  !--------------------------------------------------------------------------------!
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] initHydrology...'
  call initHydrology()


  !-----------------------------------------------------------------------!
  ! STEP 8: Initialize meteorological drivers                             !
  !-----------------------------------------------------------------------!
  
  if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  if (mynum == nnodetot-1) sendnum = 0

  if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,608,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Init_Met_Drivers_Array...'
  call init_met_drivers_array
  
  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Read_Met_Drivers_Init_Array...'
  call read_met_drivers_init_array


  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,608,MPI_COMM_WORLD,ierr)
  if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
  !-----------------------------------------------------------------------!
  ! STEP 9. Initialize ed fields that depend on the atmosphere
  !-----------------------------------------------------------------------!
  
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Ed_Init_Atm_Ar...'
  call ed_init_atm_ar
  
  !-----------------------------------------------------------------------!
  ! STEP 10. Initialized some derived variables. This must be done        !
  !          outside init_full_history_restart because it depends on some !
  !          meteorological variables that are initialized at step 9.     !
  !-----------------------------------------------------------------------!
  if (trim(runtype) == 'HISTORY') then
     do ifm=1,ngrids
        call update_derived_props(edgrid_g(ifm))
     end do
  end if
  
  !-----------------------------------------------------------------------!
  ! STEP 11. Fill the variable data-tables with all of the state          !
  !          data.  Also calculate the indexing of the vectors            !
  !          to allow for segmented I/O of hyperslabs and referencing     !
  !          of high level hierarchical data types with their parent      !
  !          types.                                                       !
  !-----------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Filltab_Alltypes...'
  call filltab_alltypes

  !-----------------------------------------------------------------------!
  ! STEP 12. Checking how the output was configure and determining the    !
  !          averaging frequency.                                         !
  !-----------------------------------------------------------------------!
  if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Finding frqsum...'
  call find_frqsum()
  
  
  !-----------------------------------------------------------------------!
  ! STEP 13. Getting the CPU time and printing the banner                 !
  !-----------------------------------------------------------------------!
  call timing(1,t1)
  w2=walltime(wtime_start)
  write(c0,'(f12.2)') t1
  write(c1,'(f12.2)') w2-w1
  write(*,'(/,a,/)') ' === Finish initialization; CPU(sec)='//&
       trim(adjustl(c0))//'; Wall(sec)='//trim(adjustl(c1))//&
       '; Time integration starts (ed_master) ===' 
  
  
  !-----------------------------------------------------------------------!
  ! STEP 14. Running the model or skipping if it is a zero time run       !
  !-----------------------------------------------------------------------!
  if (time < timmax  ) then
     call ed_model()
  end if

  return
end subroutine ed_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine find_frqsum()
   use misc_coms, only:  &
        isoutput,        &
        ifoutput,         &
        imoutput,        &
        idoutput,        &
        frqstate,         &
        frqfast,          &
        frqsum
   use consts_coms, only: day_sec

   implicit none 
   ! Determining which frequency I should use to normalize variables.
   if (ifoutput == 0 .and. isoutput == 0 .and. idoutput == 0 .and. imoutput == 0) then
         write(unit=*,fmt='(a)') '---------------------------------------------------------'
         write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write(unit=*,fmt='(a)') '---------------------------------------------------------'
         write(unit=*,fmt='(a)') ' You are running a simulation that will have no output...'
         frqsum=day_sec ! Just for avoiding build up, 
   elseif (isoutput == 0 .and. ifoutput == 0) then
       frqsum=day_sec
   elseif (isoutput > 0 .and. ifoutput == 0) then
       frqsum=min(frqstate,day_sec)
   elseif (ifoutput > 0 .and. isoutput == 0) then
       frqsum=min(frqfast,day_sec)
   else ! If both are on. If both are off the simulation won't reach this point.
       frqsum=min(min(frqstate,frqfast),day_sec)
   end if
   return
end subroutine find_frqsum
!==========================================================================================!
!==========================================================================================!
