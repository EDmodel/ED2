!==========================================================================================!
!     Main subroutine for driving the initialization process for the Ecosystem Demography  !
! Model 2, when run in coupled mode.                                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_coup_driver()
   use grid_coms            , only : ngrids                & ! intent(in)
                                   , time                  & ! intent(in)
                                   , timmax                ! ! intent(in)
   use ed_state_vars        , only : allocate_edglobals    & ! subroutine
                                   , filltab_alltypes      & ! subroutine
                                   , edgrid_g              ! ! subroutine
   use ed_misc_coms         , only : fast_diagnostics      & ! intent(in)
                                   , iyeara                & ! intent(in)
                                   , imontha               & ! intent(in)
                                   , idatea                & ! intent(in)
                                   , itimea                & ! intent(in)
                                   , runtype               & ! intent(in)
                                   , ifoutput              & ! intent(in)
                                   , idoutput              & ! intent(in)
                                   , imoutput              & ! intent(in)
                                   , iqoutput              & ! intent(in)
                                   , isoutput              & ! intent(in)
                                   , iyoutput              & ! intent(in)
                                   , writing_long          & ! intent(in)
                                   , writing_eorq          & ! intent(in)
                                   , writing_dcyc          & ! intent(in)
                                   , runtype               ! ! intent(in)
   use ed_work_vars         , only : ed_dealloc_work       & ! subroutine
                                   , work_e                ! ! intent(inout)
   use soil_coms            , only : alloc_soilgrid        ! ! subroutine
   use ed_node_coms         , only : mynum                 & ! intent(in)
                                   , nnodetot              & ! intent(in)
                                   , sendnum               & ! intent(in)
                                   , recvnum               ! ! intent(in)
   use io_params            , only : ioutput               ! ! intent(in)
   use rk4_coms             , only : checkbudget           ! ! intent(in)
   use phenology_aux        , only : first_phenology       ! ! subroutine
   use average_utils        , only : update_ed_yearly_vars & ! sub-routine
                                   , zero_ed_fmean_vars    & ! sub-routine
                                   , zero_ed_dmean_vars    & ! sub-routine
                                   , zero_ed_qmean_vars    & ! sub-routine
                                   , zero_ed_mmean_vars    ! ! sub-routine
   use hrzshade_utils       , only : init_cci_variables    ! ! subroutine
   use canopy_radiation_coms, only : ihrzrad               ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
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
   !----- External function. --------------------------------------------------------------!
   real             , external :: walltime    ! wall time
   !----- MPI header. ---------------------------------------------------------------------!
   include 'mpif.h'
   !---------------------------------------------------------------------------------------!
   
   
   ping        = 741776
   wtime_start = walltime(0.)
   wtime1      = walltime(wtime_start)




   !---------------------------------------------------------------------------------------!
   !      Set most ED model parameters that do not come from the namelist (RAMSIN).        !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ed_Ecosystem_Params...'
   call load_ed_ecosystem_params()
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     If we are running EDBRAMS then we must make sure that fast diagnostic averages    !
   ! are found when ioutput is requested.                                                  !
   !---------------------------------------------------------------------------------------!
   fast_diagnostics = fast_diagnostics .or. ioutput  /= 0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Overwrite the parameters in case a XML file is provided.                          !
   !---------------------------------------------------------------------------------------!
   if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,91,MPI_COMM_WORLD              &
                                ,MPI_STATUS_IGNORE,ierr)
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
   call overwrite_with_xml_config(mynum)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,91,MPI_COMM_WORLD,ierr)

   !---------------------------------------------------------------------------------------!
   !     Allocate soil grid arrays.                                                        !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Alloc_Soilgrid...'
   call alloc_soilgrid()

   !---------------------------------------------------------------------------------------!
   !     Set some polygon-level basic information, such as lon/lat/soil texture.           !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) then
      write (unit=*,fmt='(a)') ' [+] Set_Polygon_Coordinates (coupled)...'
   end if
   call set_polygon_coordinates_edcp()

   !---------------------------------------------------------------------------------------!
   !      Initialize inherent soil and vegetation properties.                              !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Sfcdata_ED...'
   call sfcdata_ed()

   if (trim(runtype) == 'HISTORY') then

      !------------------------------------------------------------------------------------!
      !      Initialize the model state as a replicate image of a previous state.          !
      !------------------------------------------------------------------------------------!
      if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,90,MPI_COMM_WORLD           &
                                   ,MPI_STATUS_IGNORE,ierr)

      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Init_Full_History_Restart...'
      call init_full_history_restart()

      if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,90,MPI_COMM_WORLD    &
                                          ,ierr)
   else

      !------------------------------------------------------------------------------------!
      !      Initialize state properties of polygons/sites/patches/cohorts.                !
      !------------------------------------------------------------------------------------!
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ecosystem_State...'
      call load_ecosystem_state()
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      In case this simulation will use horizontal shading, initialise the landscape    !
   ! arrays.                                                                               !
   !---------------------------------------------------------------------------------------!
   select case (ihrzrad)
   case (1)
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Init_cci_variables...'
      call init_cci_variables()
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialize hydrology related variables.                                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Initializing Hydrology...'
   call initHydrology()


   !---------------------------------------------------------------------------------------!
   !      Initialize the flux arrays that pass to the atmosphere.                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Initialise flux arrays...'
   do ifm=1,ngrids
      call newgrid(ifm)
      call initialize_ed2leaf(ifm)
   end do

   !---------------------------------------------------------------------------------------!
   !      Initialize meteorology.                                                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Initializing meteorology...'
   do ifm = 1,ngrids
      call newgrid(ifm)
      call copy_atm2lsm(ifm,.true.)
   end do

   !---------------------------------------------------------------------------------------!
   !      Initialize ed fields that depend on the atmosphere.                              !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      if (mynum == nnodetot) then
         write (unit=*,fmt='(a)') ' [+] Initialise atmospheric fields...'
      end if
      call ed_init_atm()
   end if

   !---------------------------------------------------------------------------------------!
   !      Initialize upwelling long wave and albedo from sst and air temperature.          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) then
      write (unit=*,fmt='(a)') ' [+] Initialise radiation...'
   end if
   call ed_init_radiation()


   !---------------------------------------------------------------------------------------!
   !      Initialized some derived variables.  This must be done outside                   !
   ! init_full_history_restart because it depends on some meteorological variables that    !
   ! are initialized in ed_init_atm.                                                       !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) then
      write (unit=*,fmt='(a)') ' [+] Initialise derived properties...'
   end if
   do ifm=1,ngrids
      call update_derived_props(edgrid_g(ifm))
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise drought phenology.  This should be done after the soil moisture has   !
   ! been set up.                                                                          !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Initialise phenology...'
      do ifm=1,ngrids
         call first_phenology(edgrid_g(ifm))
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Fill the variable data-tables with all of the state data.  In addition, find     !
   ! the indexing of the vectors to allow for segmented I/O of hyperslabs and referencing  !
   ! of high level hierarchical data types with their parent types.                        !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Filltab_Alltypes...'
   call filltab_alltypes()
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Check how the output was configure and determining the averaging frequency.      !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Finding frqsum...'
   call find_frqsum()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Reset the diagnostic output variables, unless this is a history run, in which    !
   ! case the variables may be partially integrated.                                       !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Reset long-term means...'
      do ifm=1,ngrids
         if (writing_long) call zero_ed_dmean_vars(edgrid_g(ifm))
         if (writing_eorq) call zero_ed_mmean_vars(edgrid_g(ifm))
         if (writing_dcyc) call zero_ed_qmean_vars(edgrid_g(ifm))
      end do

      !----- Output Initial State. --------------------------------------------------------!
      if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Update annual means...'
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
   end if

   !---------------------------------------------------------------------------------------!
   !      Allocate memory to the integration patch.                                        !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Initialise RK4 patches...'
   call initialize_rk4patches(.true.)

   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Reset averaged variables...'
   do ifm=1,ngrids
      call zero_ed_fmean_vars(edgrid_g(ifm))
   end do


   !---------------------------------------------------------------------------------------!
   !      Output initial state.                                                            !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Output initial state...'
   if (ifoutput  /= 0) call h5_output('INST')
   if (isoutput  /= 0) then
      select case (trim(runtype))
      case ('INITIAL')
         call h5_output('HIST')
      case ('HISTORY')
         call h5_output('CONT')
      end select
   end if
   if (iyoutput /= 0) call h5_output('YEAR')

   !---------------------------------------------------------------------------------------!
   !      Deallocate the work arrays.                                                      !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Deallocate work arrays...'
   do ifm=1,ngrids
      call ed_dealloc_work(work_e(ifm))
   end do

   !---------------------------------------------------------------------------------------!
   !      Get the CPU time and print the banner.                                           !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Get CPU time...'
   if (mynum == nnodetot) then
      call timing(1,cputime1)
      wtime2=walltime(wtime_start)
      write(c0,fmt='(f12.2)') cputime1
      write(c1,fmt='(f12.2)') wtime2-wtime1
      write(unit=*,fmt='(/,a,/)') ' === Finish initialization of the ED2 LSM; CPU(sec)='// &
                                  trim(adjustl(c0))//'; Wall(sec)='//trim(adjustl(c1))//   &
                                  '; Time integration starts (ed_master) ===' 
   end if

   return
end subroutine ed_coup_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine finds which frequency the model should use to normalise averaged    !
! variables.  FRQSUM should never exceed one day to avoid build up and overflows.          !
!------------------------------------------------------------------------------------------!
subroutine find_frqsum()
   use ed_misc_coms, only : unitfast  & ! intent(in)
                          , unitstate & ! intent(in)
                          , isoutput  & ! intent(in)
                          , itoutput  & ! intent(in)
                          , ifoutput  & ! intent(in)
                          , imoutput  & ! intent(in)
                          , idoutput  & ! intent(in)
                          , iqoutput  & ! intent(in)
                          , frqstate  & ! intent(in)
                          , frqfast   & ! intent(in)
                          , frqsum    ! ! intent(out)
   use consts_coms , only : day_sec   ! ! intent(in)
   use io_params   , only : ioutput   ! ! intent(in)
   implicit none


   if (ifoutput == 0 .and. isoutput == 0 .and. idoutput == 0 .and. imoutput == 0 .and.     &
       iqoutput == 0 .and. itoutput == 0 .and. ioutput  == 0 ) then
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
   !    BRAMS output is on.  Units are in seconds, test which frqsum to use (the smallest  !
   ! between frqfast and frqstate.).                                                       !
   !---------------------------------------------------------------------------------------!
   elseif (ioutput > 0) then
      frqsum=min(min(frqstate,frqfast),day_sec)

   !---------------------------------------------------------------------------------------!
   !    Mean diurnal cycle is on.  Frqfast will be in seconds, so it is likely to be the   !
   ! smallest.  The only exception is if frqstate is more frequent thant frqfast, so we    !
   ! just need to check that too.                                                          !
   !---------------------------------------------------------------------------------------!
   elseif (iqoutput > 0) then
      if (unitstate == 0) then
         frqsum = min(min(frqstate,frqfast),day_sec)
      else
         frqsum = min(frqfast,day_sec)
      end if

   !---------------------------------------------------------------------------------------!
   !     Either no instantaneous output was requested, or the user is outputting it at     !
   ! monthly or yearly scale, force it to be one day.                                      !
   !---------------------------------------------------------------------------------------!
   elseif ((isoutput == 0 .and. ifoutput == 0 .and.                                        &
            itoutput == 0 .and. ioutput == 0 ) .or.                                        &
           (ifoutput == 0 .and. itoutput == 0 .and.                                        &
            isoutput  > 0 .and. unitstate > 1) .or.                                        &
           (isoutput == 0 .and.                                                            &
            (ifoutput  > 0 .or. itoutput > 0) .and. unitfast  > 1) .or.                    &
           ((ifoutput  > 0 .or. itoutput > 0) .and.                                        &
            isoutput  > 0 .and. unitstate > 1 .and. unitfast > 1)                          &
          ) then
      frqsum=day_sec

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
!     This sub-routine will assign the longitude, latitude, and the match between the      !
! polygon and the corresponding grid point in BRAMS.                                       !
!------------------------------------------------------------------------------------------!
subroutine set_polygon_coordinates_edcp()
   use grid_coms    , only : ngrids   & ! intent(in)
                           , nzg      ! ! intent(in)
   use ed_work_vars , only : work_v   & ! intent(in)
                           , work_e   ! ! intent(in)
   use ed_state_vars, only : edtype   & ! variable type
                           , edgrid_g & ! intent(inout)
                           , gdpy     ! ! intent(in)
   use ed_node_coms , only : iskip    & ! intent(in)
                           , jskip    & ! intent(in)
                           , mynum    ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype), pointer :: cgrid
   integer               :: ifm
   integer               :: ipy
   integer               :: npoly
   !---------------------------------------------------------------------------------------!

   gridloop: do ifm = 1, ngrids
      cgrid=>edgrid_g(ifm)

      npoly = gdpy(mynum,ifm)
      polyloop: do ipy = 1,npoly

         !----- The polygon co-ordinates. -------------------------------------------------!
         cgrid%lon(ipy)              = work_v(ifm)%glon(ipy)
         cgrid%lat(ipy)              = work_v(ifm)%glat(ipy)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Ilon and ilat are the pointer to make the match between the ED2 polygon    !
         ! and the corresponding atmospheric grid cell.  They are local to the nodes,      !
         ! since the parallel code keeps the polygon and the corresponding grid cell in    !
         ! the same node.                                                                  !
         !---------------------------------------------------------------------------------!
         cgrid%ilon(ipy) = work_v(ifm)%xid(ipy)
         cgrid%ilat(ipy) = work_v(ifm)%yid(ipy)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      xatm and yatm have the global indexing for diagnostics purposes.           !
         !---------------------------------------------------------------------------------!
         cgrid%xatm(ipy) = work_v(ifm)%xid(ipy) + iskip(ifm)
         cgrid%yatm(ipy) = work_v(ifm)%yid(ipy) + jskip(ifm)
         !---------------------------------------------------------------------------------!
      end do polyloop
   end do gridloop

   return
end subroutine set_polygon_coordinates_edcp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine allocates and initialises the coupled ED-BRAMS structures.  This is !
! now called from "rams_mem_alloc" sub-routine, because we will load a few variables from  !
! the history file in case this is a history run.                                          !
!------------------------------------------------------------------------------------------!
subroutine alloc_edcp_driver(ngrds,nmxp,nmyp)
   use mem_edcp , only : ed_fluxf_g        & ! intent(inout)
                       , ed_fluxp_g        & ! intent(inout)
                       , ed_fluxpm_g       & ! intent(inout)
                       , ed_precip_g       & ! intent(inout)
                       , ed_precipm_g      & ! intent(inout)
                       , alloc_edprecip    & ! sub-routine
                       , zero_edprecip     & ! sub-routine
                       , filltab_ed_precip & ! sub-routine
                       , alloc_edflux      & ! sub-routine
                       , zero_edflux       & ! sub-routine
                       , filltab_edflux    ! ! sub-routine
   use mem_leaf , only : dtleaf            ! ! intent(in)
   use leaf_coms, only : ustmin            ! ! intent(in)
   use mem_grid , only : npatch            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in) :: ngrds
   integer, dimension(ngrds), intent(in) :: nmxp
   integer, dimension(ngrds), intent(in) :: nmyp
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: ng
   !---------------------------------------------------------------------------------------!


   !----- Allocate the actual structures. -------------------------------------------------!
   allocate(ed_fluxp_g  (ngrds))
   allocate(ed_fluxpm_g (ngrds))
   allocate(ed_fluxf_g  (ngrds))
   allocate(ed_precip_g (ngrds))
   allocate(ed_precipm_g(ngrds))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Allocate the flux arrays from terrestrial and water bodies.  These arrays are    !
   ! the same size as the leaf arrays.                                                     !
   !---------------------------------------------------------------------------------------!
   do ng = 1, ngrds
      call alloc_edflux  (ed_fluxf_g  (ng), nmxp(ng), nmyp(ng), npatch)
      call alloc_edflux  (ed_fluxp_g  (ng), nmxp(ng), nmyp(ng), npatch)
      call alloc_edflux  (ed_fluxpm_g (ng),        1,        1,      1)
      call alloc_edprecip(ed_precip_g (ng), nmxp(ng), nmyp(ng))
      call alloc_edprecip(ed_precipm_g(ng),        1,        1)
      !----- Assign zeroes to the newly allocated matrices. -------------------------------!
      call zero_edflux  (ed_fluxf_g  (ng), ustmin, dtleaf)
      call zero_edflux  (ed_fluxp_g  (ng), ustmin, dtleaf)
      call zero_edflux  (ed_fluxpm_g (ng), ustmin, dtleaf)
      call zero_edprecip(ed_precip_g (ng))
      call zero_edprecip(ed_precipm_g(ng))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Fill the variable tables, and do not bother making average arrays.  We don't  !
      ! need it.                                                                           !
      !------------------------------------------------------------------------------------!
      call filltab_edflux   (ed_fluxp_g (ng),ed_fluxpm_g (ng),0,nmxp(ng),nmyp(ng),npatch,ng)
      call filltab_ed_precip(ed_precip_g(ng),ed_precipm_g(ng),0,nmxp(ng),nmyp(ng),ng)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine alloc_edcp_driver
!==========================================================================================!
!==========================================================================================!
