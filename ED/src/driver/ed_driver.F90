
!==========================================================================================!
!==========================================================================================!
!     Main subroutine that initialises the several structures for the Ecosystem Demography !
! Model 2.  In this version 2.1 all nodes solve some polygons, including the node formerly !
! known as master.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine ed_driver()
   use update_derived_props_module
   use ism_hyd
   use ed_met_driver
   use ed_init_full_history
   use ed_init
   use grid_coms         , only : ngrids              & ! intent(in)
                                , time                & ! intent(inout)
                                , timmax              ! ! intent(inout)
   use ed_state_vars     , only : allocate_edglobals  & ! sub-routine
                                , filltab_alltypes    & ! sub-routine
                                , edgrid_g            ! ! intent(inout)
   use ed_misc_coms      , only : runtype             ! ! intent(in)
   use soil_coms         , only : alloc_soilgrid      ! ! sub-routine
   use ed_node_coms      , only : mynum               & ! intent(in)
                                , nnodetot            & ! intent(in)
                                , sendnum             & ! intent(inout)
                                , recvnum             ! ! intent(in)
   use detailed_coms     , only : idetailed           & ! intent(in)
                                , patch_keep          ! ! intent(in)
   use phenology_aux     , only : first_phenology     ! ! subroutine
   implicit none
   !----- Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h' ! MPI commons
#endif
   !----- Local variables. ----------------------------------------------------------------!
   character(len=12)           :: c0
   character(len=12)           :: c1
   integer                     :: ifm
   integer                     :: ping
   real                        :: t1
   real                        :: w1
   real                        :: w2
   real                        :: wtime_start
   logical                     :: patch_detailed
   !----- Local variable (MPI only). ------------------------------------------------------!
#if defined(RAMS_MPI)
   integer                     :: ierr
#endif
   !----- External functions. -------------------------------------------------------------!
   real             , external :: walltime    ! wall time
   !---------------------------------------------------------------------------------------!

   ping = 741776

   !---------------------------------------------------------------------------------------!
   !      Set the initial time.                                                            !
   !---------------------------------------------------------------------------------------!
   wtime_start = walltime(0.)
   w1          = walltime(wtime_start)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Set most ED model parameters that do not come from the namelist (ED2IN).         !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ed_Ecosystem_Params...'
   call load_ed_ecosystem_params()
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Overwrite the parameters in case a XML file is provided                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot-1) sendnum = 0

#if defined(RAMS_MPI)
   if (mynum /= 1) then
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,80,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
   else
      write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
   end if
#else
   write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
#endif

   call overwrite_with_xml_config(mynum)

#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,80,MPI_COMM_WORLD,ierr)
   if (nnodetot /= 1 )    call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Allocate soil grid arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Alloc_Soilgrid...'
   call alloc_soilgrid()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Set some polygon-level basic information, such as lon/lat/soil texture.          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Set_Polygon_Coordinates...'
   call set_polygon_coordinates()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialize inherent soil and vegetation properties.                              !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Sfcdata_ED...'
   call sfcdata_ed()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   if (trim(runtype) == 'HISTORY' ) then
      !------------------------------------------------------------------------------------!
      !      Initialize the model state as a replicate image of a previous  state.         !
      !------------------------------------------------------------------------------------!
      if (mynum == nnodetot-1) sendnum = 0

#if defined(RAMS_MPI)
      if (mynum /= 1) then
         call MPI_RECV(ping,1,MPI_INTEGER,recvnum,81,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      else
         write (unit=*,fmt='(a)') ' [+] Init_Full_History_Restart...'
      end if
#else
      write (unit=*,fmt='(a)') ' [+] Init_Full_History_Restart...'
#endif
      call init_full_history_restart()

#if defined(RAMS_MPI)
      if (mynum < nnodetot ) then
         call MPI_Send(ping,1,MPI_INTEGER,sendnum,81,MPI_COMM_WORLD,ierr)
      end if

      if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      !------------------------------------------------------------------------------------!
   else

      !------------------------------------------------------------------------------------!
      !      Initialize state properties of polygons/sites/patches/cohorts.                !
      !------------------------------------------------------------------------------------!
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ecosystem_State...'
      call load_ecosystem_state()
      !------------------------------------------------------------------------------------!
   end if

   !---------------------------------------------------------------------------------------!
   !      In case the runs is going to produce detailed output, we eliminate all patches   !
   ! but the one to be analysed in detail.  Special cases:                                 !
   !  0 -- Keep all patches.                                                               !
   ! -1 -- Keep the one with the highest LAI                                               !
   ! -2 -- Keep the one with the lowest LAI                                                !
   !---------------------------------------------------------------------------------------!
   patch_detailed = ibclr(idetailed,5) > 0
   if (patch_detailed) then
      call exterminate_patches_except(patch_keep)
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Initialize meteorological drivers.                                               !
   !---------------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   if (mynum == nnodetot-1) sendnum = 0

#if defined(RAMS_MPI)
   if (mynum /= 1) then
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,82,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
   else
      write (unit=*,fmt='(a)') ' [+] Init_Met_Drivers...'
   end if
#else
   write (unit=*,fmt='(a)') ' [+] Init_Met_Drivers...'
#endif

   call init_met_drivers()
   if (mynum == 1) write (unit=*,fmt='(a)') ' [+] Read_Met_Drivers_Init...'
   call read_met_drivers_init


#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,82,MPI_COMM_WORLD,ierr)
   if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise the site-level meteorological forcing.                                !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Update_met_drivers...'
   do ifm=1,ngrids
      call update_met_drivers(edgrid_g(ifm))
   end do


   !---------------------------------------------------------------------------------------!
   !      Initialize ed fields that depend on the atmosphere.                              !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Ed_Init_Atm...'
   call ed_init_atm()
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Initialize hydrology related variables.                                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] initHydrology...'
   call initHydrology()
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Initialise some derived variables.  This must be done outside                    !
   ! init_full_history_restart because it depends on some meteorological variables that    !
   ! were not initialised until the sub-routine ed_init_atm was called.                    !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids
      call update_derived_props(edgrid_g(ifm))
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise drought phenology.  This should be done after the soil moisture has   !
   ! been set up.                                                                          !
   !---------------------------------------------------------------------------------------!
   if (runtype /= 'HISTORY') then
      do ifm=1,ngrids
         call first_phenology(edgrid_g(ifm))
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Fill the variable data-tables with all of the state data.  Also calculate the    !
   ! indexing of the vectors to allow for segmented I/O of hyperslabs and referencing of   !
   ! high level hierarchical data types with their parent types.                           !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Filltab_Alltypes...'
   call filltab_alltypes
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Check how the output was configured and determine the averaging frequency.       !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Finding frqsum...'
   call find_frqsum()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Get the CPU time and print the banner.                                           !
   !---------------------------------------------------------------------------------------!
   call timing(1,t1)
   w2 = walltime(wtime_start)
   if (mynum == nnodetot) then
      write(c0,'(f12.2)') t1
      write(c1,'(f12.2)') w2-w1
      write(unit=*,fmt='(/,a,/)') ' === Finish initialization; CPU(sec)='//                &
                                  trim(adjustl(c0))//'; Wall(sec)='//trim(adjustl(c1))//   &
                                  '; Time integration starts (ed_master) ==='
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! STEP 14. Run the model or skip if it is a zero time run.                              !
   !---------------------------------------------------------------------------------------!
   if (time < timmax) then
      call ed_model()
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine finds which frequency the model should use to normalise averaged    !
! variables.  FRQSUM should never exceed one day to avoid build up and overflows.          !
!------------------------------------------------------------------------------------------!
subroutine find_frqsum()
   use ed_misc_coms, only : unitfast   & ! intent(in)
                          , unitstate  & ! intent(in)
                          , isoutput   & ! intent(in)
                          , ifoutput   & ! intent(in)
                          , itoutput   & ! intent(in)
                          , imoutput   & ! intent(in)
                          , idoutput   & ! intent(in)
                          , iqoutput   & ! intent(in)
                          , frqstate   & ! intent(in)
                          , frqfast    & ! intent(in)
                          , frqsum     ! ! intent(out)
   use consts_coms, only: day_sec

   implicit none 
   if (ifoutput == 0 .and. isoutput == 0 .and. idoutput == 0 .and. imoutput == 0 .and.     &
       iqoutput == 0 .and. itoutput == 0 ) then
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
   elseif ((isoutput == 0  .and. (ifoutput == 0 .and. itoutput == 0)) .or.                 &
           ((ifoutput == 0 .and. itoutput == 0) .and.                                      &
             isoutput  > 0 .and. unitstate > 1) .or.                                       &
           (isoutput == 0 .and.                                                            &
            (ifoutput > 0 .or. itoutput > 0) .and. unitfast  > 1) .or.                     &
           ((ifoutput  > 0 .or. itoutput > 0) .and.                                        &
             isoutput  > 0 .and. unitstate > 1 .and. unitfast > 1)                         &
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
!    This sub-routine eliminates all patches except the one you want to save...  This      !
! shouldn't be used unless you are debugging the code.                                     !
!------------------------------------------------------------------------------------------!
subroutine exterminate_patches_except(keeppa)
   use ed_state_vars  , only : edgrid_g           & ! structure
                             , edtype             & ! structure
                             , polygontype        & ! structure
                             , sitetype           & ! structure
                             , patchtype          ! ! structure
   use grid_coms      , only : ngrids             ! ! intent(in)
   use fuse_fiss_utils, only : terminate_patches  ! ! sub-routine
   !----- Arguments -----------------------------------------------------------------------!
   integer                        , intent(in)  :: keeppa
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)                   , pointer     :: cgrid
   type(polygontype)              , pointer     :: cpoly
   type(sitetype)                 , pointer     :: csite
   type(patchtype)                , pointer     :: cpatch
   integer                                      :: ifm
   integer                                      :: ipy
   integer                                      :: isi
   integer                                      :: ipa
   integer                                      :: keepact
   real             , dimension(:), allocatable :: csite_lai
   !---------------------------------------------------------------------------------------!


   gridloop: do ifm=1,ngrids
      cgrid => edgrid_g(ifm)

      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            select case(keeppa)
            case (0)
               return
            case (-2)
               !----- Keep the one with the lowest LAI. -----------------------------------!
               allocate(csite_lai(csite%npatches))
               csite_lai(:) = 0.0
               keepm2loop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)
                  if (cpatch%ncohorts > 0) csite_lai(ipa) = sum(cpatch%lai)
               end do keepm2loop
               keepact = minloc(csite_lai,dim=1)
               deallocate(csite_lai)
               !---------------------------------------------------------------------------!
            case (-1)
               !----- Keep the one with the lowest LAI. -----------------------------------!
               allocate(csite_lai(csite%npatches))
               csite_lai(:) = 0.0
               keepm1loop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)
                  if (cpatch%ncohorts > 0) csite_lai(ipa) = sum(cpatch%lai)
               end do keepm1loop
               keepact = maxloc(csite_lai,dim=1)
               deallocate(csite_lai)
               !---------------------------------------------------------------------------!
            case default
               !----- Keep a fixed patch number. ------------------------------------------!
               keepact = keeppa
               
               if (keepact > csite%npatches) then
                  write(unit=*,fmt='(a)')       '-----------------------------------------'
                  write(unit=*,fmt='(a,1x,i6)') ' - IPY      = ',ipy
                  write(unit=*,fmt='(a,1x,i6)') ' - ISI      = ',isi
                  write(unit=*,fmt='(a,1x,i6)') ' - NPATCHES = ',csite%npatches
                  write(unit=*,fmt='(a,1x,i6)') ' - KEEPPA   = ',keeppa
                  write(unit=*,fmt='(a)')       '-----------------------------------------'
                  call fail_whale ()
                  call fatal_error('KEEPPA can''t be greater than NPATCHES'                &
                                  ,'exterminate_patches_except','ed_driver.f90')
               end if
            end select

            patchloop: do ipa=1,csite%npatches
               if (ipa == keepact) then
                  csite%area(ipa) = 1.0
               else
                  csite%area(ipa) = 0.0
               end if
            end do patchloop

            call terminate_patches(csite)

         end do siteloop
      end do polyloop
   end do gridloop

   return
end subroutine exterminate_patches_except
!==========================================================================================!
!==========================================================================================!
