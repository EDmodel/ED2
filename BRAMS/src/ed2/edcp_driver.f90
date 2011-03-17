!==========================================================================================!
!     Main subroutine for driving the initialization process for the Ecosystem Demography  !
! Model 2, when run in coupled mode.                                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_coup_driver()
   use grid_coms     , only : ngrids              & ! intent(in)
                            , time                & ! intent(in)
                            , timmax              ! ! intent(in)
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
                            , iqoutput            & ! intent(in)
                            , iyoutput            & ! intent(in)
                            , integration_scheme  & ! intent(in)
                            , runtype             ! ! intent(in)
   use ed_work_vars  , only : ed_dealloc_work     & ! subroutine
                            , work_e              ! ! intent(inout)
   use soil_coms     , only : alloc_soilgrid      ! ! subroutine
   use ed_node_coms  , only : mynum               & ! intent(in)
                            , nnodetot            & ! intent(in)
                            , sendnum             & ! intent(in)
                            , recvnum             ! ! intent(in)
   use io_params     , only : ioutput             ! ! intent(in)
   use rk4_coms      , only : checkbudget         ! ! intent(in)
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
   ! STEP 1: Set the special diagnostic parameter                                          !
   !                                                                                       !
   !     Checking if the user has indicated a need for any of the fast flux diagnostic     !
   ! variables, these are used in conditions of ifoutput,idoutput and imoutput conditions. !
   ! If they are not >0, then set the logical, fast_diagnostics to false.                  !
   !---------------------------------------------------------------------------------------!
   fast_diagnostics = checkbudget   .or. ifoutput /= 0 .or. idoutput /= 0 .or.             &
                      iqoutput /= 0 .or. imoutput /= 0 .or. ioutput  /= 0 .or.             &
                      iyoutput /= 0

   !---------------------------------------------------------------------------------------!
   ! STEP 2: Set the ED model parameters                                                   !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ed_Ecosystem_Params...'
   call load_ed_ecosystem_params()
   
   !---------------------------------------------------------------------------------------!
   ! STEP 3: Overwrite the parameters in case a XML file is provided.                      !
   !---------------------------------------------------------------------------------------!
   if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,91,MPI_COMM_WORLD              &
                                ,MPI_STATUS_IGNORE,ierr)
   if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
   call overwrite_with_xml_config(mynum)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,91,MPI_COMM_WORLD,ierr)

   !---------------------------------------------------------------------------------------!
   ! STEP 4: Allocate soil grid arrays.                                                    !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Alloc_Soilgrid...'
   call alloc_soilgrid()

   !---------------------------------------------------------------------------------------!
   ! STEP 5: Set some polygon-level basic information, such as lon/lat/soil texture.       !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) then
      write (unit=*,fmt='(a)') ' [+] Set_Polygon_Coordinates (coupled)...'
   end if
   call set_polygon_coordinates_edcp()

   !---------------------------------------------------------------------------------------!
   ! STEP 6: Initialize inherent soil and vegetation properties.                           !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Sfcdata_ED...'
   call sfcdata_ed()

   if (trim(runtype) == 'HISTORY') then

      !------------------------------------------------------------------------------------!
      ! STEP 7A: Initialize the model state as a replicate image of a previous state.      !
      !------------------------------------------------------------------------------------!
      if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,90,MPI_COMM_WORLD           &
                                   ,MPI_STATUS_IGNORE,ierr)

      if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Init_Full_History_Restart...'
      call init_full_history_restart()

      if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,90,MPI_COMM_WORLD    &
                                          ,ierr)
   else

      !------------------------------------------------------------------------------------!
      ! STEP 7B: Initialize state properties of polygons/sites/patches/cohorts.            !
      !------------------------------------------------------------------------------------!
      if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Load_Ecosystem_State...'
      call load_ecosystem_state()
   end if

   !---------------------------------------------------------------------------------------!
   ! STEP 8: Initialize hydrology related variables.                                       !
   !---------------------------------------------------------------------------------------!
   if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Initializing Hydrology...'
   call initHydrology()


   !---------------------------------------------------------------------------------------!
   ! STEP 9: Initialize the flux arrays that pass to the atmosphere.                       !
   !---------------------------------------------------------------------------------------!
   if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Initialise flux arrays...'
   do ifm=1,ngrids
      call newgrid(ifm)
      call initialize_ed2leaf(ifm)
   end do

   !---------------------------------------------------------------------------------------!
   ! STEP 10: Initialize meteorology.                                                      !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Initializing meteorology...'
   do ifm = 1,ngrids
      call newgrid(ifm)
      call copy_atm2lsm(ifm,.true.)
   end do

   !---------------------------------------------------------------------------------------!
   ! STEP 11. Initialize ed fields that depend on the atmosphere.                          !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] ed_init_atm...'
      call ed_init_atm()
   end if

   !---------------------------------------------------------------------------------------!
   ! STEP 12. Initialize upwelling long wave and albedo from sst and air temperature.      !
   !---------------------------------------------------------------------------------------!
   call ed_init_radiation()


   !---------------------------------------------------------------------------------------!
   ! STEP 13. Initialized some derived variables.  This must be done outside               !
   !          init_full_history_restart because it depends on some meteorological          !
   !          variables that are initialized at step 9.                                    !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) == 'HISTORY') then
      do ifm=1,ngrids
         call update_derived_props(edgrid_g(ifm))
      end do
   end if

   !---------------------------------------------------------------------------------------!
   ! STEP 14. Fill the variable data-tables with all of the state data.  Also calculate    !
   !          the indexing of the vectors to allow for segmented I/O of hyperslabs and     !
   !          referencing of high level hierarchical data types with their parent types.   !
   !---------------------------------------------------------------------------------------!
   if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Filltab_Alltypes...'
   call filltab_alltypes()


   !---------------------------------------------------------------------------------------!
   ! STEP 15. Checking how the output was configure and determining the averaging          !
   !          frequency.                                                                   !
   !---------------------------------------------------------------------------------------!
   if (mynum == 1) write(unit=*,fmt='(a)') ' [+] Finding frqsum...'
   call find_frqsum()

   !---------------------------------------------------------------------------------------!
   ! STEP 16. Zero and initialize diagnostic output variables                              !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY') then
      if (imoutput > 0 .or. iqoutput > 0) then
         do ifm=1,ngrids
            call zero_ed_monthly_output_vars(edgrid_g(ifm))
            call zero_ed_daily_output_vars(edgrid_g(ifm))
         end do
      elseif (idoutput > 0) then
         do ifm=1,ngrids
            call zero_ed_daily_output_vars(edgrid_g(ifm))
         end do
      end if

      !----- Output Initial State. --------------------------------------------------------!
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
   end if

   !---------------------------------------------------------------------------------------!
   ! STEP 17: Allocate memory to the integration patch.                                    !
   !---------------------------------------------------------------------------------------!
   call initialize_rk4patches(.true.)

   do ifm=1,ngrids
      call reset_averaged_vars(edgrid_g(ifm))
   end do


   !---------------------------------------------------------------------------------------!
   ! STEP 18:  Output initial state.                                                       !
   !---------------------------------------------------------------------------------------!
   if (ioutput  > 0) call h5_output('HIST')
   if (ioutput  > 0) call h5_output('INST')
   if (iyoutput > 0) call h5_output('YEAR')

   !---------------------------------------------------------------------------------------!
   ! STEP 19: Deallocate the work arrays.                                                  !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids
      call ed_dealloc_work(work_e(ifm))
   end do

   !---------------------------------------------------------------------------------------!
   ! STEP 20. Get the CPU time and print the banner.                                       !
   !---------------------------------------------------------------------------------------!
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
!     This sub-routine will assign the longitude, latitude, and soil class for all         !
! non-empty grid points.  The difference between this and the original ED-2 stand alone    !
! sub-routine is that here we will also assign the match between the polygon and the       !
! corresponding grid point in BRAMS.                                                       !
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

         !----- Soil texture class. -------------------------------------------------------!
         cgrid%ntext_soil(1:nzg,ipy) = work_v(ifm)%ntext(ipy)


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
   use mem_edcp     , only : ed_fluxf_g        & ! intent(inout)
                           , ed_fluxp_g        & ! intent(inout)
                           , wgrid_g           & ! intent(inout)
                           , ed_precip_g       & ! intent(inout)
                           , ed_precipm_g      & ! intent(inout)
                           , alloc_edprecip    & ! sub-routine
                           , zero_edprecip     & ! sub-routine
                           , filltab_ed_precip & ! sub-routine
                           , alloc_edflux      & ! sub-routine
                           , zero_edflux       ! ! sub-routine

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in) :: ngrds
   integer, dimension(ngrds), intent(in) :: nmxp
   integer, dimension(ngrds), intent(in) :: nmyp
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: ng
   !---------------------------------------------------------------------------------------!


   !----- Allocate the actual structures. -------------------------------------------------!
   allocate(wgrid_g     (ngrds))
   allocate(ed_fluxp_g  (ngrds))
   allocate(ed_fluxf_g  (ngrds))
   allocate(ed_precip_g (ngrds))
   allocate(ed_precipm_g(ngrds))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Allocate the flux arrays from terrestrial and water bodies.  These arrays are    !
   ! the same size as the leaf arrays.                                                     !
   !---------------------------------------------------------------------------------------!
   do ng = 1, ngrds
      call alloc_edflux  (ed_fluxf_g  (ng), nmxp(ng), nmyp(ng))
      call alloc_edflux  (ed_fluxp_g  (ng), nmxp(ng), nmyp(ng))
      call alloc_edflux  (wgrid_g     (ng), nmxp(ng), nmyp(ng))
      call alloc_edprecip(ed_precip_g (ng), nmxp(ng), nmyp(ng))
      call alloc_edprecip(ed_precipm_g(ng),        1,        1)
      !----- Assign zeroes to the newly allocated matrices. -------------------------------!
      call zero_edflux  (ed_fluxf_g  (ng))
      call zero_edflux  (ed_fluxp_g  (ng))
      call zero_edflux  (wgrid_g     (ng))
      call zero_edprecip(ed_precip_g (ng))
      call zero_edprecip(ed_precipm_g(ng))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Fill the variable tables, and do not bother making average arrays.  We don't  !
      ! need it.                                                                           !
      !------------------------------------------------------------------------------------!
      call filltab_ed_precip(ed_precip_g(ng),ed_precipm_g(ng),0,nmxp(ng),nmyp(ng),ng)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine alloc_edcp_driver
!==========================================================================================!
!==========================================================================================!
