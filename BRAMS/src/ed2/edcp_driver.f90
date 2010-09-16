!==========================================================================================!
!     Main subroutine for driving the initialization process for the Ecosystem Demography  !
! Model 2, when run in coupled mode.                                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_coup_driver()
  
  use grid_coms     , only : ngrids              & ! intent(in)
                           , time                & ! intent(in)
                           , timmax              ! ! intent(in)
  use mem_edcp      , only : wgrid_g             & ! intent(inout)
                           , ed_fluxp_g          & ! intent(inout)
                           , ed_fluxf_g          & ! intent(inout)
                           , ed_precip_g         ! ! intent(inout)
  use ed_state_vars , only : allocate_edglobals  & ! subroutine
                           , filltab_alltypes    & ! subroutine
                           , edgrid_g            ! ! subroutine

  use ed_misc_coms  , only : fast_diagnostics    & ! intent(in)
                           , iyeara              & ! intent(in)
                           , imontha             & ! intent(in)
                           , idatea              & ! intent(in)
                           , itimea              & ! intent(in)
                           , runtype             & ! intent(in)
                           , ifoutput            & ! intent(in)
                           , idoutput            & ! intent(in)
                           , imoutput            & ! intent(in)
                           , integration_scheme  ! ! intent(in)
  use ed_work_vars  , only : ed_dealloc_work     & ! subroutine
                           , work_e              ! ! intent(inout)
  use soil_coms     , only : alloc_soilgrid      ! ! subroutine
  use ed_node_coms  , only : mynum               & ! intent(in)
                           , nnodetot            & ! intent(in)
                           , sendnum             & ! intent(in)
                           , recvnum             & ! intent(in)
                           , mmxp                & ! intent(in)
                           , mmyp                ! ! intent(in)
  use io_params     , only : ioutput
  use rk4_coms      , only : checkbudget         ! ! intent(in)
  implicit none
  !----- Local variables. -----------------------------------------------------------------!
  character(len=12)           :: c0
  character(len=12)           :: c1
  integer                     :: i
  integer                     :: j
  integer                     :: ifm
  integer                     :: nndtflg
  integer                     :: ipy
  integer                     :: id1
  integer                     :: id2
  integer                     :: jd1
  integer                     :: jd2
  integer                     :: ierr
  integer                     :: igr
  integer                     :: ping 
  real                        :: wtime1
  real                        :: wtime2
  real                        :: wtime_start     ! wall time
  real                        :: cputime1
  !----- External function. ---------------------------------------------------------------!
  real             , external :: walltime    ! wall time
  !----- MPI header. ----------------------------------------------------------------------!
  include 'mpif.h'
  !----------------------------------------------------------------------------------------!
  
  
  ping        = 741776
  wtime_start = walltime(0.)
  wtime1      = walltime(wtime_start)

  !---------------------------------------------------------------------------!
  ! STEP 1: Set the special diagnostic parameter                              !
  !                                                                           !
  !     Checking if the user has indicated a need for any of the fast flux    !
  ! diagnostic variables, these are used in conditions of ifoutput,idoutput   !
  ! and imoutput conditions.                                                  !
  ! If they are not >0, then set the logical, fast_diagnostics to false.      !
  !---------------------------------------------------------------------------!
  fast_diagnostics = checkbudget .or. ifoutput /= 0 .or. idoutput /= 0 .or.   &
                     imoutput /= 0 .or. ioutput /= 0

  !---------------------------------------------------------------------------!
  ! STEP 2: Set the ED model parameters                                       !
  !---------------------------------------------------------------------------!
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ed_Ecosystem_Params...'
  call load_ed_ecosystem_params()
  
  !---------------------------------------------------------------------------!
  ! STEP 3: Overwrite the parameters in case a XML file is provided           !
  !---------------------------------------------------------------------------!

  if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,91,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
  
  call overwrite_with_xml_config(mynum)
  
  if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,91,MPI_COMM_WORLD,ierr)

  
  !---------------------------------------------------------------------------!
  ! STEP 4: Allocate soil grid arrays                                         !
  !---------------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Alloc_Soilgrid...'
  call alloc_soilgrid()
  
  !---------------------------------------------------------------------------------!
  ! STEP 5: Set some polygon-level basic information, such as lon/lat/soil texture  !
  !---------------------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Set_Polygon_Coordinates...'
  call set_polygon_coordinates()
  
  !---------------------------------------------------------------------------!
  ! STEP 6: Initialize inherent soil and vegetation properties.               !
  !---------------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Sfcdata_ED...'
  call sfcdata_ed()
  
  if (trim(runtype) == 'HISTORY') then
       
     !-----------------------------------------------------------------------!
     ! STEP 7A: Initialize the model state as a replicate image of a previous
     !          state.
     !-----------------------------------------------------------------------!

     if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,90,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  
     if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Init_Full_History_Restart...'
     call init_full_history_restart()
     
     if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,90,MPI_COMM_WORLD,ierr)
!     if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
     
  else
     
     !------------------------------------------------------------------------!
     ! STEP 7B: Initialize state properties of polygons/sites/patches/cohorts !
     !------------------------------------------------------------------------!
     
     if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Load_Ecosystem_State...'
     call load_ecosystem_state()

  end if

  !--------------------------------------------------------------------------------!
  ! STEP 8: Initialize hydrology related variables                                 !
  !--------------------------------------------------------------------------------!

  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Initializing Hydrology...'
  call initHydrology()

  !-----------------------------------------------------------------------!
  ! STEP 9: Inform edtypes which atmospheric cell to look at
  !          
  !-----------------------------------------------------------------------!
  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Setting the correspondent atmospheric cells...'
  do ifm=1,ngrids
     call set_edtype_atm(ifm)
  end do


  !-----------------------------------------------------------------------!
  ! STEP 10: Initialize the flux arrays that pass to the atmosphere
  !-----------------------------------------------------------------------!
  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Allocating Transfer Arrays...'
  allocate(wgrid_g(ngrids))
  allocate(ed_fluxp_g(ngrids))
  allocate(ed_fluxf_g(ngrids))
  allocate(ed_precip_g(ngrids))
  do ifm=1,ngrids
     call newgrid(ifm)
     call initialize_ed2leaf(ifm,mmxp(ifm),mmyp(ifm))
  enddo

  !-----------------------------------------------------------------------!
  ! STEP 11: Initialize meteorology                                       !
  !-----------------------------------------------------------------------!
  
  if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Initializing meteorology...'
  do ifm = 1,ngrids
     call newgrid(ifm)
     call copy_atm2lsm(ifm,.true.)
  enddo

     
  !-----------------------------------------------------------------------!
  ! STEP 12. Initialize ed fields that depend on the atmosphere
  !-----------------------------------------------------------------------!
  
  if (trim(runtype) /= 'HISTORY') then
     if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] ed_init_coup_atm...'
     call ed_init_coup_atm
  endif

  !-----------------------------------------------------------------------!
  ! STEP 13. Initialize upwelling long wave and albedo
  !            from sst and air temperature.
  !-----------------------------------------------------------------------!
  
  call ed_init_radiation



  !-----------------------------------------------------------------------!
  ! STEP 14. Initialized some derived variables. This must be done        !
  !          outside init_full_history_restart because it depends on some !
  !          meteorological variables that are initialized at step 9.     !
  !-----------------------------------------------------------------------!
  if (trim(runtype) == 'HISTORY') then
     do ifm=1,ngrids
        call update_derived_props(edgrid_g(ifm))
     end do
  end if
  
  !-----------------------------------------------------------------------!
  ! STEP 15. Fill the variable data-tables with all of the state          !
  !          data.  Also calculate the indexing of the vectors            !
  !          to allow for segmented I/O of hyperslabs and referencing     !
  !          of high level hierarchical data types with their parent      !
  !          types.                                                       !
  !-----------------------------------------------------------------------!
  
  if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Filltab_Alltypes...'
  call filltab_alltypes

  

  !-----------------------------------------------------------------------!
  ! STEP 16. Checking how the output was configure and determining the    !
  !          averaging frequency.                                         !
  !-----------------------------------------------------------------------!
  if (mynum == 1) write(unit=*,fmt='(a)') ' [+] Finding frqsum...'
  call find_frqsum()

  !-----------------------------------------------------------------------!
  ! STEP 17. Zero and initialize diagnostic output variables              !
  !-----------------------------------------------------------------------!

  if (imoutput > 0) then
     do ifm=1,ngrids
        call zero_ed_monthly_output_vars(edgrid_g(ifm))
        call zero_ed_daily_output_vars(edgrid_g(ifm))
     end do
  elseif (idoutput > 0) then
     do ifm=1,ngrids
        call zero_ed_daily_output_vars(edgrid_g(ifm))
     end do
  end if

  !-----------------------------------------------------------------------!
  ! STEP 18: Allocate memory to the integration patch
  !-----------------------------------------------------------------------!
  call initialize_rk4patches(.true.)

  do ifm=1,ngrids
     call reset_averaged_vars(edgrid_g(ifm))
  end do

  
  !-----------------------------------------------------------------------!
  ! STEP 19: Deallocate the work arrays
  !-----------------------------------------------------------------------!
  do ifm=1,ngrids
     call ed_dealloc_work(work_e(ifm))
  enddo

  !-----------------------------------------------------------------------!
  ! STEP 20. Getting the CPU time and printing the banner                 !
  !-----------------------------------------------------------------------!
  if (mynum == nnodetot) then
     call timing(1,cputime1)
     wtime2=walltime(wtime_start)
     write(c0,'(f12.2)') cputime1
     write(c1,'(f12.2)') wtime2-wtime1
     write(*,'(/,a,/)') ' === Finish initialization of the ED2 LSM; CPU(sec)='//&
          trim(adjustl(c0))//'; Wall(sec)='//trim(adjustl(c1))//&
          '; Time integration starts (ed_master) ===' 
  endif

  return
end subroutine ed_coup_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine find_frqsum()
   use ed_misc_coms, only:  &
        unitfast,        &
        unitstate,       &
        isoutput,        &
        itoutput,        &
        ifoutput,        &
        imoutput,        &
        idoutput,        &
        frqstate,        &
        frqfast,         &
        frqsum
   use consts_coms, only: day_sec
   use io_params  , only: ioutput
   implicit none 
   !---------------------------------------------------------------------------------------!
   ! Determining which frequency I should use to normalize variables. FRQSUM should never  !
   ! exceed 1 day.                                                                         !
   !---------------------------------------------------------------------------------------!
   if (ifoutput == 0 .and. isoutput == 0 .and. idoutput == 0 .and. imoutput == 0 .and.     &
       itoutput == 0 .and. ioutput  == 0 ) then
      write(unit=*,fmt='(a)') '---------------------------------------------------------'
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '---------------------------------------------------------'
      write(unit=*,fmt='(a)') ' You are running a simulation that will have no output...'
      frqsum=day_sec ! This avoids the number to get incredibly large.

   !---------------------------------------------------------------------------------------!
   !     Either no instantaneous output was requested, or the user is outputting it at     !
   ! monthly or yearly scale, force it to be one day.                                      !
   !---------------------------------------------------------------------------------------!
   elseif ((isoutput == 0 .and. ifoutput == 0 .and.                                        &
            itoutput == 0 .and. ioutput == 0 ) .or.                                        &
           (ifoutput == 0 .and. itoutput == 0 .and.                                        &
            isoutput  > 0 .and. unitstate > 1) .or.                                        &
           (isoutput == 0 .and. (ifoutput  > 0 .or. itoutput > 0) .and.                    &
            unitfast  > 1) .or.                                                            &
           ((ifoutput  > 0 .or. itoutput > 0) .and. isoutput  > 0 .and.                    &
             unitstate > 1 .and. unitfast > 1)                                             &
          ) then
      frqsum=day_sec
   !---------------------------------------------------------------------------------------!
   !    BRAMS output is on, and because this is a coupled run, the units are in seconds,   !
   ! test which frqsum to use (the smallest between frqfast and frqstate.).                !
   !---------------------------------------------------------------------------------------!
   elseif (ioutput > 0) then
      frqsum=min(min(frqstate,frqfast),day_sec)
   !---------------------------------------------------------------------------------------!
   !    Only restarts, and the unit is in seconds, test which frqsum to use.               !
   !---------------------------------------------------------------------------------------!
   elseif (ifoutput == 0 .and. itoutput == 0 .and. isoutput > 0) then
      frqsum=min(frqstate,day_sec)
   !---------------------------------------------------------------------------------------!
   !    Only fast analysis, and the unit is in seconds, test which frqsum to use.          !
   !---------------------------------------------------------------------------------------!
   elseif (isoutput == 0 .and. (ifoutput > 0 .or. itoutput > 0)) then
      frqsum=min(frqfast,day_sec)
   !---------------------------------------------------------------------------------------!
   !    Both are on and both outputs are in seconds or day scales. Choose the minimum      !
   ! between them and one day.                                                             !
   !---------------------------------------------------------------------------------------!
   elseif (unitfast < 2 .and. unitstate < 2) then 
      frqsum=min(min(frqstate,frqfast),day_sec)
   !---------------------------------------------------------------------------------------!
   !    Both are on but unitstate is in month or years. Choose the minimum between frqfast !
   ! and one day.                                                                          !
   !---------------------------------------------------------------------------------------!
   elseif (unitfast < 2) then 
      frqsum=min(frqfast,day_sec)
   !---------------------------------------------------------------------------------------!
   !    Both are on but unitfast is in month or years. Choose the minimum between frqstate !
   ! and one day.                                                                          !
   !---------------------------------------------------------------------------------------!
   else
      frqsum=min(frqstate,day_sec)
   end if

   return
end subroutine find_frqsum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine set_edtype_atm(ifm)

   use ed_work_vars,only:work_e
   use ed_state_vars,only:edtype,edgrid_g
   use ed_node_coms,only:mmxp,mmyp,iwest,jsouth
   
   implicit none

   integer :: ipy,i,j,ifm
   type(edtype),pointer :: cgrid

   cgrid=>edgrid_g(ifm)

   ipy = 0
   do j = 1,mmyp(ifm)
      do i=1,mmxp(ifm)
         if (work_e(ifm)%land(i,j)) then
            ipy = ipy + 1

            ! THE FOLLOWING VECTORS SHOULD STAY LOCAL! THEY ARE NEEDED FOR
            ! MATCHING POLYGONS WITH THE ATMOSPHERIC GRIDS
            cgrid%ilon(ipy) = work_e(ifm)%xatm(i,j)
            cgrid%ilat(ipy) = work_e(ifm)%yatm(i,j)


            ! THE FOLLOWING VECTORS SHOULD HAVE GLOBAL INDEXING FOR DIAGNOSTICS
            ! This was also set when the ed set polygon routine was celled
            ! but unlike the offline version xatm in the work array
            ! is locally indexed in the coupled version, so we need to add
            ! back in the tile offsets.  This will fix that and overwrite

            cgrid%yatm(ipy) = work_e(ifm)%yatm(i,j) + jsouth(ifm) - 1
            cgrid%xatm(ipy) = work_e(ifm)%xatm(i,j) + iwest(ifm) - 1

         endif
      end do
   end do

   return
end subroutine set_edtype_atm
!==========================================================================================!
!==========================================================================================!
