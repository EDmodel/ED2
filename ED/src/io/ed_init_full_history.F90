!==========================================================================================!
!==========================================================================================!
!     This subroutine will set a full history start for the simulation.  In a full history !
! start (runtype = 'HISTORY'), we assumes that you are continuing with the exact same      !
! configuration as a given model simulation that wrote the file in which you are using.    !
! In this case, each node starts with a list of polygons, and searches the HDF history     !
! file for these polygons.  It expects to find at least one polygon in the file within 250 !
! meters of proximity (although it should be exactly at the same place).  If it does not   !
! find this it stops.                                                                      !
!     Assuming that it finds the polygons, the subroutine then fills each of these         !
! polygons with data, and traverses the data hierarchical tree that roots from each        !
! polygon, initializing the model to the exact same state as the end of the previous run.  !
! A search based restart does not expect to find exact matches, and may not use all of the !
! polygons in the file.  Nonetheless, it will traverse the tree from these polygons and    !
! populate the model states with what is found in the files tree.                          !
!------------------------------------------------------------------------------------------!
subroutine init_full_history_restart()
   use ed_max_dims      , only : n_pft                 & ! intent(in)
                               , str_len               ! ! intent(in)
   use ed_misc_coms     , only : sfilin                & ! intent(in)
                               , current_time          & ! intent(in)
                               , max_poihist_dist      ! ! intent(in)
   use ed_state_vars    , only : polygontype           & ! structure
                               , sitetype              & ! structure
                               , patchtype             & ! structure
                               , edtype                & ! structure
                               , edgrid_g              & ! structure
                               , allocate_sitetype     & ! subroutine
                               , allocate_patchtype    & ! subroutine
                               , allocate_polygontype  ! ! subroutine
   use soil_coms        , only : alloc_soilgrid        ! ! subroutine
   use grid_coms        , only : ngrids                ! ! intent(in)
   use phenology_startup, only : phenology_init        ! ! subroutine
   use ed_node_coms     , only : mynum                 & ! intent(in)
                               , nmachs                & ! intent(in)
                               , nnodetot              & ! intent(in)
                               , mchnum                & ! intent(in)
                               , machs                 ! ! intent(in)
   use hdf5
   use hdf5_coms        , only : file_id               & ! intent(inout)
                               , dset_id               & ! intent(inout)
                               , dspace_id             & ! intent(inout)
                               , plist_id              & ! intent(inout)
                               , globdims              & ! intent(inout)
                               , chnkdims              & ! intent(inout)
                               , chnkoffs              ! ! intent(inout)
   implicit none
   !------ Local variables. ---------------------------------------------------------------!
   type(edtype)                        , pointer     :: cgrid
   type(polygontype)                   , pointer     :: cpoly
   type(sitetype)                      , pointer     :: csite
   type(patchtype)                     , pointer     :: cpatch
   character(len=3)                                  :: cgr
   character(len=str_len)                            :: hnamel
   integer               , dimension(:), allocatable :: pysi_n
   integer               , dimension(:), allocatable :: pysi_id
   integer               , dimension(:), allocatable :: sipa_n
   integer               , dimension(:), allocatable :: sipa_id
   integer               , dimension(:), allocatable :: paco_n
   integer               , dimension(:), allocatable :: paco_id
   integer                                           :: ngr
   integer                                           :: ifpy
   integer                                           :: ipy
   integer                                           :: isi
   integer                                           :: ipa
   integer                                           :: py_index
   integer                                           :: si_index
   integer                                           :: pa_index
   integer                                           :: hdferr
   logical               , dimension(:), allocatable :: is_burnt
   logical                                           :: exists
   real                  , dimension(:), allocatable :: file_lats
   real                  , dimension(:), allocatable :: file_lons
   real                                              :: mindist
   real                                              :: polydist
   real(kind=8)                                      :: dbletime
   !------ Local constants. ---------------------------------------------------------------!
   character(len=1)                    , parameter   :: vnam = 'S'
   !------ External functions. ------------------------------------------------------------!
   real                                , external    :: dist_gc
   !---------------------------------------------------------------------------------------!



   write (unit=*,fmt='(a)') '-----------------------------------------------------'
   write (unit=*,fmt='(a)') '  Loading Full State (HISTORY)'


   !----- Open the HDF environment. -------------------------------------------------------!
   call h5open_f(hdferr)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Turn off automatic error printing.  This is done because there may be datasets    !
   ! that are not in the file, but it is OK.  If data that should be there are missing,    !
   ! ED2 error reporting will detect it.  If something that can't be found is missing, the !
   ! following call can be bypassed.  Note that automatic error reporting is turned back   !
   ! on in the end.                                                                        !
   !---------------------------------------------------------------------------------------!
   call h5eset_auto_f(0,hdferr)


   gridloop: do ngr=1,ngrids
      cgrid => edgrid_g(ngr)

      !------------------------------------------------------------------------------------!
      !     Make the history file name.                                                    !
      !------------------------------------------------------------------------------------!
      write(cgr,'(a1,i2.2)') 'g',ngr

      dbletime = dble(current_time%time)

      call makefnam(hnamel,sfilin(1),dbletime,current_time%year,current_time%month         &
                   ,current_time%date,0,vnam,cgr,'h5 ')
      inquire(file=trim(hnamel),exist=exists)

      if (.not.exists) then
         !----- History couldn't be found.  Stop the run. ---------------------------------!
         call fatal_error ('File '//trim(hnamel)//' not found.'                            &
                          ,'init_full_history_restart','ed_init_full_history.F90')
      else
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a,1x,i8)') 'Error opening HDF5 file - error - ',hdferr
            write(unit=*,fmt='(a,1x,a)' ) '- Filename: ',trim(hnamel)
            call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)            &
                            ,'init_full_history_restart','ed_init_full_history.F90')
         end if
      end if


      !------------------------------------------------------------------------------------!
      !     Retrieve global vector sizes.                                                  !
      !------------------------------------------------------------------------------------!
      globdims    = 0_8
      chnkdims    = 0_8
      chnkoffs    = 0_8
      globdims(1) = 1_8
      
      call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%npolygons_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%nsites_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%npatches_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%ncohorts_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Retrieve the mapping of the data tree.                                         !
      !------------------------------------------------------------------------------------!
      globdims    = 0_8
      globdims(1) = int(cgrid%npolygons_global,8)
      
      allocate(pysi_n(cgrid%npolygons_global))
      allocate(pysi_id(cgrid%npolygons_global))
      
      call h5dopen_f(file_id,'PYSI_N', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,pysi_n,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'PYSI_ID', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,pysi_id,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      globdims(1) = int(cgrid%nsites_global,8)
      
      allocate(sipa_n(cgrid%nsites_global))
      allocate(sipa_id(cgrid%nsites_global))
      
      call h5dopen_f(file_id,'SIPA_N', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,sipa_n,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'SIPA_ID', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,sipa_id,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      globdims(1) = int(cgrid%npatches_global,8)
      
      allocate(paco_n(cgrid%npatches_global))
      allocate(paco_id(cgrid%npatches_global))
      
      call h5dopen_f(file_id,'PACO_N', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,paco_n,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'PACO_ID', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,paco_id,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      !------------------------------------------------------------------------------------!
      
      
      !------------------------------------------------------------------------------------!
      !      Retrieve the polygon coordinates data.                                        !
      !------------------------------------------------------------------------------------!
      globdims(1) = int(cgrid%npolygons_global,8)
      allocate(file_lats(cgrid%npolygons_global))
      allocate(file_lons(cgrid%npolygons_global))
      
      call h5dopen_f(file_id,'LATITUDE', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_REAL,file_lats,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)

      call h5dopen_f(file_id,'LONGITUDE', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_REAL,file_lons,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)



      !------------------------------------------------------------------------------------!
      !     Loop the polygons in the model state and match them with those in the file.    !
      ! After the match, we must walk through the data from that polygon and initialize.   !
      ! We check the distance between the expected coordinates and the retrieved ones, and !
      ! they ought to be less than 250 metres apart, otherwise we can't use the polygon.   !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         py_index   = 0
         mindist    = huge(1.)
         do ifpy = 1,cgrid%npolygons_global
            polydist = dist_gc(file_lons(ifpy),cgrid%lon(ipy)                              &
                              ,file_lats(ifpy),cgrid%lat(ipy))
            if (polydist < mindist) then
               mindist    = polydist
               py_index   = ifpy
            end if
         end do

         !---------------------------------------------------------------------------------!
         !     Check whether the closest polygon is close.                                 !
         !---------------------------------------------------------------------------------!
         if (mindist > max_poihist_dist) then
            write (unit=*,fmt='(a)'          ) '------------------------------------------'
            write (unit=*,fmt='(a)'          ) ' None of the polygons in the history file'
            write (unit=*,fmt='(a)'          ) '    is enough close!  The model will stop!'
            write (unit=*,fmt='(a)'          ) ' ED Polygon:'
            write (unit=*,fmt='(a,1x,i8)'    ) ' - Polygon   :',ipy
            write (unit=*,fmt='(a,1x,es12.5)') ' - Longitude :',cgrid%lon(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') ' - Latitude  :',cgrid%lat(ipy)
            write (unit=*,fmt='(a)'          ) ' History file''s closest polygon:'
            write (unit=*,fmt='(a,1x,i8)'    ) ' - Polygon   :',py_index
            write (unit=*,fmt='(a,1x,es12.5)') ' - Longitude :',file_lons(py_index)
            write (unit=*,fmt='(a,1x,es12.5)') ' - Latitude  :',file_lats(py_index)
            write (unit=*,fmt='(a)'          ) '------------------------------------------'

            call fatal_error('Mismatch between polygon and dataset'                        &
                            ,'init_full_history_restart','ed_init_full_history.F90')
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Get all necessary polygon variables associated with this index for the     !
         ! current polygon.                                                                !
         !---------------------------------------------------------------------------------!
         call fill_history_grid(cgrid,ipy,py_index)

         if (pysi_n(py_index) > 0) then
            !----- Allocate the polygontype structure (site level). -----------------------!
            call allocate_polygontype(cpoly,pysi_n(py_index))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Get all necessary site variables associated with this index for the      !
            ! current polygon.                                                             !
            !------------------------------------------------------------------------------!
            allocate (is_burnt(pysi_n(py_index)))
            is_burnt(:) = .false.
            call fill_history_polygon(cpoly,pysi_id(py_index),cgrid%nsites_global          &
                                     ,pysi_n(py_index),is_burnt)
            
            siteloop: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)
               
               !------ Calculate the index of this site's data in the HDF. ----------------!
               si_index = pysi_id(py_index) + isi - 1

               if (sipa_n(si_index) > 0) then
                  !----- Allocate the sitetype structure (patch level). -------------------!
                  call allocate_sitetype(csite,sipa_n(si_index))

                  !------------------------------------------------------------------------!
                  !     Get all necessary site variables associated with this index for    !
                  ! the current site.                                                      !
                  !------------------------------------------------------------------------!
                  call fill_history_site(csite,sipa_id(si_index),cgrid%npatches_global     &
                                        ,is_burnt(isi))

                  patchloop: do ipa = 1,csite%npatches
                     cpatch => csite%patch(ipa)
                     pa_index = sipa_id(si_index) + ipa - 1

                     if (paco_n(pa_index) > 0) then
                        !----- Allocate the patchtype structure (cohort level). -----------!
                        call allocate_patchtype(cpatch,paco_n(pa_index))
                        
                        !------------------------------------------------------------------!
                        !     Get all necessary site variables associated with this index  !
                        ! for the current patch.                                           !
                        !------------------------------------------------------------------!
                        call fill_history_patch(cpatch,paco_id(pa_index)                   &
                                               ,cgrid%ncohorts_global)
                        !------------------------------------------------------------------!
                     else
                        cpatch%ncohorts = 0
                     endif
                  end do patchloop
               else
                  write (unit=*,fmt='(a)'          ) '------------------------------------'
                  write (unit=*,fmt='(a)'          ) ' Found a site with no patches.'
                  write (unit=*,fmt='(a)'          ) ' This is not allowed.'
                  write (unit=*,fmt='(a)'          ) ' ED Polygon and Site:'
                  write (unit=*,fmt='(a,1x,i8)'    ) ' - Polygon   :',ipy
                  write (unit=*,fmt='(a,1x,i8)'    ) ' - Site      :',isi
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Longitude :',cgrid%lon(ipy)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Latitude  :',cgrid%lat(ipy)
                  write (unit=*,fmt='(a)'          ) '------------------------------------'
                  call fatal_error('Attempted to load an empty site.'                      &
                                  ,'init_full_history_restart','ed_init_full_history.F90')
               end if

            end do siteloop
            deallocate (is_burnt)
            
         else
            write (unit=*,fmt='(a)'          ) '------------------------------------'
            write (unit=*,fmt='(a)'          ) ' Found a polygon with no sites.'
            write (unit=*,fmt='(a)'          ) ' This is not allowed.'
            write (unit=*,fmt='(a)'          ) ' ED Polygon:'
            write (unit=*,fmt='(a,1x,i8)'    ) ' - Polygon   :',ipy
            write (unit=*,fmt='(a,1x,es12.5)') ' - Longitude :',cgrid%lon(ipy)
            write (unit=*,fmt='(a,1x,es12.5)') ' - Latitude  :',cgrid%lat(ipy)
            write (unit=*,fmt='(a)'          ) '------------------------------------'
            call fatal_error('Attempted to load an empty polygon.'                         &
                            ,'init_full_history_restart','ed_init_full_history.F90')
         end if
      end do polyloop


      call h5fclose_f(file_id, hdferr)
      if (hdferr /= 0) then
          print*,hdferr
          call fatal_error('Could not close the HDF file'                                  &
                          ,'init_full_history_restart','ed_init_full_history.F90')
          
      end if

      deallocate(file_lats)
      deallocate(file_lons)
      deallocate(paco_n)
      deallocate(paco_id)
      deallocate(sipa_n)
      deallocate(sipa_id)
      deallocate(pysi_n)
      deallocate(pysi_id)

   end do gridloop

   !---------------------------------------------------------------------------------------!
   !     Turn automatic error reporting back on.  This is probably unnecessary, because    !
   ! the environment is about to be flushed.                                               !
   !---------------------------------------------------------------------------------------!
   call h5eset_auto_f(1,hdferr)
   !---------------------------------------------------------------------------------------!



   !----- Close the HDF environment. ------------------------------------------------------!
   call h5close_f(hdferr)
   !---------------------------------------------------------------------------------------!



   !----- Load the anthropogenic disturbance (or set them all to zero). -------------------!
   write(unit=*,fmt='(a,i2.2)') ' Checking anthropogenic disturbance.  Node: ',mynum
   call landuse_init()
   !---------------------------------------------------------------------------------------!


   !----- Load phenology in case it is prescribed (or set them with defaults). ------------!
   write(unit=*,fmt='(a,i2.2)') ' Checking prescribed phenology.  Node: ',mynum
   call phenology_init()

   return
end subroutine init_full_history_restart
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !---------------------------------------------------------------------------------------!

   if (cgrid%npolygons == 0) return


   !---------------------------------------------------------------------------------------!
   !      Split the routine into smaller routines to avoid segmentation violation.         !
   !---------------------------------------------------------------------------------------!
   call fill_history_grid_p11     (cgrid,ipy,py_index)
   call fill_history_grid_p11dmean(cgrid,ipy,py_index)
   call fill_history_grid_p11mmean(cgrid,ipy,py_index)
   call fill_history_grid_p12     (cgrid,ipy,py_index)
   call fill_history_grid_m11     (cgrid,ipy,py_index)
   call fill_history_grid_p19     (cgrid,ipy,py_index)
   call fill_history_grid_m12     (cgrid,ipy,py_index)
   call fill_history_grid_p146    (cgrid,ipy,py_index)
   !---------------------------------------------------------------------------------------!

   return
end subroutine fill_history_grid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_p11(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables.                                                      !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(cgrid%npolygons_global,8)
   chnkdims(1) = 1_8
   chnkoffs(1) = int(py_index - 1,8)
   memdims (1) = 1_8
   memoffs (1) = 0_8
   memsize (1) = 1_8
   !---------------------------------------------------------------------------------------!

   call hdf_getslab_i(cgrid%load_adjacency          (ipy:ipy)                              &
                     ,'LOAD_ADJACENCY '           ,dsetrank,iparallel,.true. ,foundvar)

   call hdf_getslab_r(cgrid%wbar                    (ipy:ipy)                              &
                     ,'WBAR '                     ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%Te                      (ipy:ipy)                              &
                     ,'TE '                       ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%zbar                    (ipy:ipy)                              &
                     ,'ZBAR '                     ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%sheat                   (ipy:ipy)                              &
                     ,'SHEAT '                    ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%baseflow                (ipy:ipy)                              &
                     ,'BASEFLOW '                 ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%runoff                  (ipy:ipy)                              &
                     ,'RUNOFF '                   ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%qrunoff                 (ipy:ipy)                              &
                     ,'QRUNOFF '                  ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%swliq                   (ipy:ipy)                              &
                     ,'SWLIQ '                    ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%total_agb               (ipy:ipy)                              &
                     ,'TOTAL_AGB '                ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_basal_area        (ipy:ipy)                              &
                     ,'TOTAL_BASAL_AREA '         ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_agb_growth        (ipy:ipy)                              &
                     ,'TOTAL_AGB_GROWTH '         ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_agb_mort          (ipy:ipy)                              &
                     ,'TOTAL_AGB_MORT '           ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_agb_recruit       (ipy:ipy)                              &
                     ,'TOTAL_AGB_RECRUIT '        ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_basal_area_growth (ipy:ipy)                              &
                     ,'TOTAL_BASAL_AREA_GROWTH '  ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_basal_area_mort   (ipy:ipy)                              &
                     ,'TOTAL_BASAL_AREA_MORT '    ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%total_basal_area_recruit(ipy:ipy)                              &
                     ,'TOTAL_BASAL_AREA_RECRUIT ' ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%cosz                    (ipy:ipy)                              &
                     ,'COSZ '                     ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%cbudget_initialstorage  (ipy:ipy)                              &
                     ,'CBUDGET_INITIALSTORAGE '   ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%cbudget_nep             (ipy:ipy)                              &
                     ,'CBUDGET_NEP '              ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%nbudget_initialstorage  (ipy:ipy)                              &
                     ,'NBUDGET_INITIALSTORAGE '   ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%Cleaf_litter_flux       (ipy:ipy)                              &
                     ,'CLEAF_LITTER_FLUX '        ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%Croot_litter_flux       (ipy:ipy)                              &
                     ,'CROOT_LITTER_FLUX '        ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%Nleaf_litter_flux       (ipy:ipy)                              &
                     ,'NLEAF_LITTER_FLUX '        ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%Nroot_litter_flux       (ipy:ipy)                              &
                     ,'NROOT_LITTER_FLUX '        ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%Nbiomass_uptake         (ipy:ipy)                              &
                     ,'NBIOMASS_UPTAKE '          ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%Ngross_min              (ipy:ipy)                              &
                     ,'NGROSS_MIN '               ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%Nnet_min                (ipy:ipy)                              &
                     ,'NNET_MIN '                 ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%fast_soil_c             (ipy:ipy)                              &
                     ,'FAST_SOIL_C_PY '           ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%slow_soil_c             (ipy:ipy)                              &
                     ,'SLOW_SOIL_C_PY '           ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%struct_soil_c           (ipy:ipy)                              &
                     ,'STRUCT_SOIL_C_PY '         ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%struct_soil_l           (ipy:ipy)                              &
                     ,'STRUCT_SOIL_L_PY '         ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%cwd_c                   (ipy:ipy)                              &
                     ,'CWD_C_PY '                 ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%fast_soil_n             (ipy:ipy)                              &
                     ,'FAST_SOIL_N_PY '           ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%mineral_soil_n          (ipy:ipy)                              &
                     ,'MINERAL_SOIL_N_PY '        ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cgrid%cwd_n                   (ipy:ipy)                              &
                     ,'CWD_N_PY '                 ,dsetrank,iparallel,.false.,foundvar)
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_p11
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_p11dmean(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables.                                                      !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(cgrid%npolygons_global,8)
   chnkdims(1) = 1_8
   chnkoffs(1) = int(py_index - 1,8)
   memdims (1) = 1_8
   memoffs (1) = 0_8
   memsize (1) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !       Try to load daily and monthly means if the user is going to write them.         !
   !---------------------------------------------------------------------------------------!
   !------ Daily means. -------------------------------------------------------------------!
   if (writing_long) then
      call hdf_getslab_r(cgrid%dmean_nppleaf        (ipy:ipy)                              &
                        ,'DMEAN_NPPLEAF_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nppfroot       (ipy:ipy)                              &
                        ,'DMEAN_NPPFROOT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nppsapwood     (ipy:ipy)                              &
                        ,'DMEAN_NPPSAPWOOD_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nppcroot       (ipy:ipy)                              &
                        ,'DMEAN_NPPCROOT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nppseeds       (ipy:ipy)                              &
                        ,'DMEAN_NPPSEEDS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nppwood        (ipy:ipy)                              &
                        ,'DMEAN_NPPWOOD_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nppdaily       (ipy:ipy)                              &
                        ,'DMEAN_NPPDAILY_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_A_decomp       (ipy:ipy)                              &
                        ,'DMEAN_A_DECOMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_Af_decomp      (ipy:ipy)                              &
                        ,'DMEAN_AF_DECOMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_co2_residual   (ipy:ipy)                              &
                        ,'DMEAN_CO2_RESIDUAL_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_energy_residual(ipy:ipy)                              &
                        ,'DMEAN_ENERGY_RESIDUAL_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_water_residual (ipy:ipy)                              &
                        ,'DMEAN_WATER_RESIDUAL_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_gpp            (ipy:ipy)                              &
                        ,'DMEAN_GPP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_npp            (ipy:ipy)                              &
                        ,'DMEAN_NPP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_resp      (ipy:ipy)                              &
                        ,'DMEAN_LEAF_RESP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_root_resp      (ipy:ipy)                              &
                        ,'DMEAN_ROOT_RESP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_growth_resp(ipy:ipy)                             &
                        ,'DMEAN_LEAF_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_root_growth_resp(ipy:ipy)                             &
                        ,'DMEAN_ROOT_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sapa_growth_resp(ipy:ipy)                             &
                        ,'DMEAN_SAPA_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sapb_growth_resp(ipy:ipy)                             &
                        ,'DMEAN_SAPB_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_storage_resp(ipy:ipy)                            &
                        ,'DMEAN_LEAF_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_root_storage_resp(ipy:ipy)                            &
                        ,'DMEAN_ROOT_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sapa_storage_resp(ipy:ipy)                            &
                        ,'DMEAN_SAPA_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sapb_storage_resp(ipy:ipy)                            &
                        ,'DMEAN_SAPB_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_plresp         (ipy:ipy)                              &
                        ,'DMEAN_PLRESP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_energy    (ipy:ipy)                              &
                        ,'DMEAN_LEAF_ENERGY_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_water     (ipy:ipy)                              &
                        ,'DMEAN_LEAF_WATER_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_hcap      (ipy:ipy)                              &
                        ,'DMEAN_LEAF_HCAP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_vpdef     (ipy:ipy)                              &
                        ,'DMEAN_LEAF_VPDEF_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_temp      (ipy:ipy)                              &
                        ,'DMEAN_LEAF_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_fliq      (ipy:ipy)                              &
                        ,'DMEAN_LEAF_FLIQ_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_gsw       (ipy:ipy)                              &
                        ,'DMEAN_LEAF_GSW_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_leaf_gbw       (ipy:ipy)                              &
                        ,'DMEAN_LEAF_GBW_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wood_energy    (ipy:ipy)                              &
                        ,'DMEAN_WOOD_ENERGY_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wood_water     (ipy:ipy)                              &
                        ,'DMEAN_WOOD_WATER_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wood_hcap      (ipy:ipy)                              &
                        ,'DMEAN_WOOD_HCAP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wood_temp      (ipy:ipy)                              &
                        ,'DMEAN_WOOD_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wood_fliq      (ipy:ipy)                              &
                        ,'DMEAN_WOOD_FLIQ_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wood_gbw       (ipy:ipy)                              &
                        ,'DMEAN_WOOD_GBW_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_fs_open        (ipy:ipy)                              &
                        ,'DMEAN_FS_OPEN_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_fsw            (ipy:ipy)                              &
                        ,'DMEAN_FSW_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_fsn            (ipy:ipy)                              &
                        ,'DMEAN_FSN_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_a_open         (ipy:ipy)                              &
                        ,'DMEAN_A_OPEN_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_a_closed       (ipy:ipy)                              &
                        ,'DMEAN_A_CLOSED_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_a_net          (ipy:ipy)                              &
                        ,'DMEAN_A_NET_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_a_light        (ipy:ipy)                              &
                        ,'DMEAN_A_LIGHT_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_a_rubp         (ipy:ipy)                              &
                        ,'DMEAN_A_RUBP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_a_co2          (ipy:ipy)                              &
                        ,'DMEAN_A_CO2_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_psi_open       (ipy:ipy)                              &
                        ,'DMEAN_PSI_OPEN_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_psi_closed     (ipy:ipy)                              &
                        ,'DMEAN_PSI_CLOSED_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_water_supply   (ipy:ipy)                              &
                        ,'DMEAN_WATER_SUPPLY_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_par_l          (ipy:ipy)                              &
                        ,'DMEAN_PAR_L_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_par_l_beam     (ipy:ipy)                              &
                        ,'DMEAN_PAR_L_BEAM_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_par_l_diff     (ipy:ipy)                              &
                        ,'DMEAN_PAR_L_DIFF_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rshort_l       (ipy:ipy)                              &
                        ,'DMEAN_RSHORT_L_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rlong_l        (ipy:ipy)                              &
                        ,'DMEAN_RLONG_L_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sensible_lc    (ipy:ipy)                              &
                        ,'DMEAN_SENSIBLE_LC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_vapor_lc       (ipy:ipy)                              &
                        ,'DMEAN_VAPOR_LC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_transp         (ipy:ipy)                              &
                        ,'DMEAN_TRANSP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_intercepted_al (ipy:ipy)                              &
                        ,'DMEAN_INTERCEPTED_AL_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wshed_lg       (ipy:ipy)                              &
                        ,'DMEAN_WSHED_LG_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rshort_w       (ipy:ipy)                              &
                        ,'DMEAN_RSHORT_W_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rlong_w        (ipy:ipy)                              &
                        ,'DMEAN_RLONG_W_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sensible_wc    (ipy:ipy)                              &
                        ,'DMEAN_SENSIBLE_WC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_vapor_wc       (ipy:ipy)                              &
                        ,'DMEAN_VAPOR_WC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_intercepted_aw (ipy:ipy)                              &
                        ,'DMEAN_INTERCEPTED_AW_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_wshed_wg       (ipy:ipy)                              &
                        ,'DMEAN_WSHED_WG_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rh             (ipy:ipy)                              &
                        ,'DMEAN_RH_PY               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_cwd_rh         (ipy:ipy)                              &
                        ,'DMEAN_CWD_RH_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nep            (ipy:ipy)                              &
                        ,'DMEAN_NEP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rk4step        (ipy:ipy)                              &
                        ,'DMEAN_RK4STEP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_available_water(ipy:ipy)                              &
                        ,'DMEAN_AVAILABLE_WATER_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_theiv      (ipy:ipy)                              &
                        ,'DMEAN_CAN_THEIV_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_theta      (ipy:ipy)                              &
                        ,'DMEAN_CAN_THETA_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_vpdef      (ipy:ipy)                              &
                        ,'DMEAN_CAN_VPDEF_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_temp       (ipy:ipy)                              &
                        ,'DMEAN_CAN_TEMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_shv        (ipy:ipy)                              &
                        ,'DMEAN_CAN_SHV_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_co2        (ipy:ipy)                              &
                        ,'DMEAN_CAN_CO2_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_rhos       (ipy:ipy)                              &
                        ,'DMEAN_CAN_RHOS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_prss       (ipy:ipy)                              &
                        ,'DMEAN_CAN_PRSS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_gnd_temp       (ipy:ipy)                              &
                        ,'DMEAN_GND_TEMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_gnd_shv        (ipy:ipy)                              &
                        ,'DMEAN_GND_SHV_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_can_ggnd       (ipy:ipy)                              &
                        ,'DMEAN_CAN_GGND_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sfcw_depth     (ipy:ipy)                              &
                        ,'DMEAN_SFCW_DEPTH_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sfcw_energy    (ipy:ipy)                              &
                        ,'DMEAN_SFCW_ENERGY_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sfcw_mass      (ipy:ipy)                              &
                        ,'DMEAN_SFCW_MASS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sfcw_temp      (ipy:ipy)                              &
                        ,'DMEAN_SFCW_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sfcw_fliq      (ipy:ipy)                              &
                        ,'DMEAN_SFCW_FLIQ_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rshort_gnd     (ipy:ipy)                              &
                        ,'DMEAN_RSHORT_GND_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_par_gnd        (ipy:ipy)                              &
                        ,'DMEAN_PAR_GND_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rlong_gnd      (ipy:ipy)                              &
                        ,'DMEAN_RLONG_GND_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rlongup        (ipy:ipy)                              &
                        ,'DMEAN_RLONGUP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_parup          (ipy:ipy)                              &
                        ,'DMEAN_PARUP_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_nirup          (ipy:ipy)                              &
                        ,'DMEAN_NIRUP_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rshortup       (ipy:ipy)                              &
                        ,'DMEAN_RSHORTUP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rnet           (ipy:ipy)                              &
                        ,'DMEAN_RNET_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_albedo         (ipy:ipy)                              &
                        ,'DMEAN_ALBEDO_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_albedo_par     (ipy:ipy)                              &
                        ,'DMEAN_ALBEDO_PAR_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_albedo_nir     (ipy:ipy)                              &
                        ,'DMEAN_ALBEDO_NIR_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_rlong_albedo   (ipy:ipy)                              &
                        ,'DMEAN_RLONG_ALBEDO_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_ustar          (ipy:ipy)                              &
                        ,'DMEAN_USTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_tstar          (ipy:ipy)                              &
                        ,'DMEAN_TSTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_qstar          (ipy:ipy)                              &
                        ,'DMEAN_QSTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_cstar          (ipy:ipy)                              &
                        ,'DMEAN_CSTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_carbon_ac      (ipy:ipy)                              &
                        ,'DMEAN_CARBON_AC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_carbon_st      (ipy:ipy)                              &
                        ,'DMEAN_CARBON_ST_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_vapor_gc       (ipy:ipy)                              &
                        ,'DMEAN_VAPOR_GC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_vapor_ac       (ipy:ipy)                              &
                        ,'DMEAN_VAPOR_AC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_throughfall    (ipy:ipy)                              &
                        ,'DMEAN_THROUGHFALL_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_runoff         (ipy:ipy)                              &
                        ,'DMEAN_RUNOFF_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_drainage       (ipy:ipy)                              &
                        ,'DMEAN_DRAINAGE_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sensible_gc    (ipy:ipy)                              &
                        ,'DMEAN_SENSIBLE_GC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sensible_ac    (ipy:ipy)                              &
                        ,'DMEAN_SENSIBLE_AC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_qthroughfall   (ipy:ipy)                              &
                        ,'DMEAN_QTHROUGHFALL_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_qrunoff        (ipy:ipy)                              &
                        ,'DMEAN_QRUNOFF_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_qdrainage      (ipy:ipy)                              &
                        ,'DMEAN_QDRAINAGE_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_theiv      (ipy:ipy)                              &
                        ,'DMEAN_ATM_THEIV_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_theta      (ipy:ipy)                              &
                        ,'DMEAN_ATM_THETA_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_temp       (ipy:ipy)                              &
                        ,'DMEAN_ATM_TEMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_vpdef      (ipy:ipy)                              &
                        ,'DMEAN_ATM_VPDEF_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_shv        (ipy:ipy)                              &
                        ,'DMEAN_ATM_SHV_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_rshort     (ipy:ipy)                              &
                        ,'DMEAN_ATM_RSHORT_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_rshort_diff(ipy:ipy)                              &
                        ,'DMEAN_ATM_RSHORT_DIFF_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_par        (ipy:ipy)                              &
                        ,'DMEAN_ATM_PAR_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_par_diff   (ipy:ipy)                              &
                        ,'DMEAN_ATM_PAR_DIFF_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_rlong      (ipy:ipy)                              &
                        ,'DMEAN_ATM_RLONG_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_vels       (ipy:ipy)                              &
                        ,'DMEAN_ATM_VELS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_rhos       (ipy:ipy)                              &
                        ,'DMEAN_ATM_RHOS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_prss       (ipy:ipy)                              &
                        ,'DMEAN_ATM_PRSS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_atm_co2        (ipy:ipy)                              &
                        ,'DMEAN_ATM_CO2_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_pcpg           (ipy:ipy)                              &
                        ,'DMEAN_PCPG_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_qpcpg          (ipy:ipy)                              &
                        ,'DMEAN_QPCPG_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_dpcpg          (ipy:ipy)                              &
                        ,'DMEAN_DPCPG_PY            ',dsetrank,iparallel,.false.,foundvar)
   end if
   return
end subroutine fill_history_grid_p11dmean
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_p11mmean(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables.                                                      !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(cgrid%npolygons_global,8)
   chnkdims(1) = 1_8
   chnkoffs(1) = int(py_index - 1,8)
   memdims (1) = 1_8
   memoffs (1) = 0_8
   memsize (1) = 1_8
   !---------------------------------------------------------------------------------------!



   !------ Monthly means. -----------------------------------------------------------------!
   if (writing_eorq) then
      call hdf_getslab_r(cgrid%mmean_fast_soil_c    (ipy:ipy)                              &
                        ,'MMEAN_FAST_SOIL_C_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_slow_soil_c    (ipy:ipy)                              &
                        ,'MMEAN_SLOW_SOIL_C_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_struct_soil_c  (ipy:ipy)                              &
                        ,'MMEAN_STRUCT_SOIL_C_PY    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_struct_soil_l  (ipy:ipy)                              &
                        ,'MMEAN_STRUCT_SOIL_L_PY    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_fast_soil_n    (ipy:ipy)                              &
                        ,'MMEAN_FAST_SOIL_N_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_cwd_c          (ipy:ipy)                              &
                        ,'MMEAN_CWD_C_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_mineral_soil_n (ipy:ipy)                              &
                        ,'MMEAN_MINERAL_SOIL_N_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_cwd_n          (ipy:ipy)                              &
                        ,'MMEAN_CWD_N_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_gpp            (ipy:ipy)                              &
                        ,'MMEAN_GPP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_npp            (ipy:ipy)                              &
                        ,'MMEAN_NPP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_resp      (ipy:ipy)                              &
                        ,'MMEAN_LEAF_RESP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_root_resp      (ipy:ipy)                              &
                        ,'MMEAN_ROOT_RESP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_growth_resp(ipy:ipy)                             &
                        ,'MMEAN_LEAF_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_root_growth_resp(ipy:ipy)                             &
                        ,'MMEAN_ROOT_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sapa_growth_resp(ipy:ipy)                             &
                        ,'MMEAN_SAPA_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sapb_growth_resp(ipy:ipy)                             &
                        ,'MMEAN_SAPB_GROWTH_RESP_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_storage_resp(ipy:ipy)                            &
                        ,'MMEAN_LEAF_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_root_storage_resp(ipy:ipy)                            &
                        ,'MMEAN_ROOT_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sapa_storage_resp(ipy:ipy)                            &
                        ,'MMEAN_SAPA_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sapb_storage_resp(ipy:ipy)                            &
                        ,'MMEAN_SAPB_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_plresp         (ipy:ipy)                              &
                        ,'MMEAN_PLRESP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_energy    (ipy:ipy)                              &
                        ,'MMEAN_LEAF_ENERGY_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_water     (ipy:ipy)                              &
                        ,'MMEAN_LEAF_WATER_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_hcap      (ipy:ipy)                              &
                        ,'MMEAN_LEAF_HCAP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_vpdef     (ipy:ipy)                              &
                        ,'MMEAN_LEAF_VPDEF_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_temp      (ipy:ipy)                              &
                        ,'MMEAN_LEAF_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_fliq      (ipy:ipy)                              &
                        ,'MMEAN_LEAF_FLIQ_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_gsw       (ipy:ipy)                              &
                        ,'MMEAN_LEAF_GSW_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_gbw       (ipy:ipy)                              &
                        ,'MMEAN_LEAF_GBW_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wood_energy    (ipy:ipy)                              &
                        ,'MMEAN_WOOD_ENERGY_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wood_water     (ipy:ipy)                              &
                        ,'MMEAN_WOOD_WATER_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wood_hcap      (ipy:ipy)                              &
                        ,'MMEAN_WOOD_HCAP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wood_temp      (ipy:ipy)                              &
                        ,'MMEAN_WOOD_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wood_fliq      (ipy:ipy)                              &
                        ,'MMEAN_WOOD_FLIQ_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wood_gbw       (ipy:ipy)                              &
                        ,'MMEAN_WOOD_GBW_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_fs_open        (ipy:ipy)                              &
                        ,'MMEAN_FS_OPEN_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_fsw            (ipy:ipy)                              &
                        ,'MMEAN_FSW_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_fsn            (ipy:ipy)                              &
                        ,'MMEAN_FSN_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_a_open         (ipy:ipy)                              &
                        ,'MMEAN_A_OPEN_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_a_closed       (ipy:ipy)                              &
                        ,'MMEAN_A_CLOSED_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_a_net          (ipy:ipy)                              &
                        ,'MMEAN_A_NET_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_a_light        (ipy:ipy)                              &
                        ,'MMEAN_A_LIGHT_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_a_rubp         (ipy:ipy)                              &
                        ,'MMEAN_A_RUBP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_a_co2          (ipy:ipy)                              &
                        ,'MMEAN_A_CO2_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_psi_open       (ipy:ipy)                              &
                        ,'MMEAN_PSI_OPEN_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_psi_closed     (ipy:ipy)                              &
                        ,'MMEAN_PSI_CLOSED_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_water_supply   (ipy:ipy)                              &
                        ,'MMEAN_WATER_SUPPLY_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_par_l          (ipy:ipy)                              &
                        ,'MMEAN_PAR_L_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_par_l_beam     (ipy:ipy)                              &
                        ,'MMEAN_PAR_L_BEAM_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_par_l_diff     (ipy:ipy)                              &
                        ,'MMEAN_PAR_L_DIFF_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rshort_l       (ipy:ipy)                              &
                        ,'MMEAN_RSHORT_L_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rlong_l        (ipy:ipy)                              &
                        ,'MMEAN_RLONG_L_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sensible_lc    (ipy:ipy)                              &
                        ,'MMEAN_SENSIBLE_LC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_vapor_lc       (ipy:ipy)                              &
                        ,'MMEAN_VAPOR_LC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_transp         (ipy:ipy)                              &
                        ,'MMEAN_TRANSP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_intercepted_al (ipy:ipy)                              &
                        ,'MMEAN_INTERCEPTED_AL_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wshed_lg       (ipy:ipy)                              &
                        ,'MMEAN_WSHED_LG_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rshort_w       (ipy:ipy)                              &
                        ,'MMEAN_RSHORT_W_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rlong_w        (ipy:ipy)                              &
                        ,'MMEAN_RLONG_W_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sensible_wc    (ipy:ipy)                              &
                        ,'MMEAN_SENSIBLE_WC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_vapor_wc       (ipy:ipy)                              &
                        ,'MMEAN_VAPOR_WC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_intercepted_aw (ipy:ipy)                              &
                        ,'MMEAN_INTERCEPTED_AW_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_wshed_wg       (ipy:ipy)                              &
                        ,'MMEAN_WSHED_WG_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppleaf        (ipy:ipy)                              &
                        ,'MMEAN_NPPLEAF_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppfroot       (ipy:ipy)                              &
                        ,'MMEAN_NPPFROOT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppsapwood     (ipy:ipy)                              &
                        ,'MMEAN_NPPSAPWOOD_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppcroot       (ipy:ipy)                              &
                        ,'MMEAN_NPPCROOT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppseeds       (ipy:ipy)                              &
                        ,'MMEAN_NPPSEEDS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppwood        (ipy:ipy)                              &
                        ,'MMEAN_NPPWOOD_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppdaily       (ipy:ipy)                              &
                        ,'MMEAN_NPPDAILY_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rh             (ipy:ipy)                              &
                        ,'MMEAN_RH_PY               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_cwd_rh         (ipy:ipy)                              &
                        ,'MMEAN_CWD_RH_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nep            (ipy:ipy)                              &
                        ,'MMEAN_NEP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rk4step        (ipy:ipy)                              &
                        ,'MMEAN_RK4STEP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_available_water(ipy:ipy)                              &
                        ,'MMEAN_AVAILABLE_WATER_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_theiv      (ipy:ipy)                              &
                        ,'MMEAN_CAN_THEIV_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_theta      (ipy:ipy)                              &
                        ,'MMEAN_CAN_THETA_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_vpdef      (ipy:ipy)                              &
                        ,'MMEAN_CAN_VPDEF_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_temp       (ipy:ipy)                              &
                        ,'MMEAN_CAN_TEMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_shv        (ipy:ipy)                              &
                        ,'MMEAN_CAN_SHV_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_co2        (ipy:ipy)                              &
                        ,'MMEAN_CAN_CO2_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_rhos       (ipy:ipy)                              &
                        ,'MMEAN_CAN_RHOS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_prss       (ipy:ipy)                              &
                        ,'MMEAN_CAN_PRSS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_gnd_temp       (ipy:ipy)                              &
                        ,'MMEAN_GND_TEMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_gnd_shv        (ipy:ipy)                              &
                        ,'MMEAN_GND_SHV_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_can_ggnd       (ipy:ipy)                              &
                        ,'MMEAN_CAN_GGND_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sfcw_depth     (ipy:ipy)                              &
                        ,'MMEAN_SFCW_DEPTH_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sfcw_energy    (ipy:ipy)                              &
                        ,'MMEAN_SFCW_ENERGY_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sfcw_mass      (ipy:ipy)                              &
                        ,'MMEAN_SFCW_MASS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sfcw_temp      (ipy:ipy)                              &
                        ,'MMEAN_SFCW_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sfcw_fliq      (ipy:ipy)                              &
                        ,'MMEAN_SFCW_FLIQ_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rshort_gnd     (ipy:ipy)                              &
                        ,'MMEAN_RSHORT_GND_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_par_gnd        (ipy:ipy)                              &
                        ,'MMEAN_PAR_GND_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rlong_gnd      (ipy:ipy)                              &
                        ,'MMEAN_RLONG_GND_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rlongup        (ipy:ipy)                              &
                        ,'MMEAN_RLONGUP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_parup          (ipy:ipy)                              &
                        ,'MMEAN_PARUP_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nirup          (ipy:ipy)                              &
                        ,'MMEAN_NIRUP_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rshortup       (ipy:ipy)                              &
                        ,'MMEAN_RSHORTUP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rnet           (ipy:ipy)                              &
                        ,'MMEAN_RNET_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_albedo         (ipy:ipy)                              &
                        ,'MMEAN_ALBEDO_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_albedo_par     (ipy:ipy)                              &
                        ,'MMEAN_ALBEDO_PAR_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_albedo_nir     (ipy:ipy)                              &
                        ,'MMEAN_ALBEDO_NIR_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_rlong_albedo   (ipy:ipy)                              &
                        ,'MMEAN_RLONG_ALBEDO_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_ustar          (ipy:ipy)                              &
                        ,'MMEAN_USTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_tstar          (ipy:ipy)                              &
                        ,'MMEAN_TSTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_qstar          (ipy:ipy)                              &
                        ,'MMEAN_QSTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_cstar          (ipy:ipy)                              &
                        ,'MMEAN_CSTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_carbon_ac      (ipy:ipy)                              &
                        ,'MMEAN_CARBON_AC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_carbon_st      (ipy:ipy)                              &
                        ,'MMEAN_CARBON_ST_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_vapor_gc       (ipy:ipy)                              &
                        ,'MMEAN_VAPOR_GC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_vapor_ac       (ipy:ipy)                              &
                        ,'MMEAN_VAPOR_AC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_throughfall    (ipy:ipy)                              &
                        ,'MMEAN_THROUGHFALL_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_runoff         (ipy:ipy)                              &
                        ,'MMEAN_RUNOFF_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_drainage       (ipy:ipy)                              &
                        ,'MMEAN_DRAINAGE_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sensible_gc    (ipy:ipy)                              &
                        ,'MMEAN_SENSIBLE_GC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sensible_ac    (ipy:ipy)                              &
                        ,'MMEAN_SENSIBLE_AC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_qthroughfall   (ipy:ipy)                              &
                        ,'MMEAN_QTHROUGHFALL_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_qrunoff        (ipy:ipy)                              &
                        ,'MMEAN_QRUNOFF_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_qdrainage      (ipy:ipy)                              &
                        ,'MMEAN_QDRAINAGE_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppleaf        (ipy:ipy)                              &
                        ,'MMEAN_NPPLEAF_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppfroot       (ipy:ipy)                              &
                        ,'MMEAN_NPPFROOT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppsapwood     (ipy:ipy)                              &
                        ,'MMEAN_NPPSAPWOOD_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppcroot       (ipy:ipy)                              &
                        ,'MMEAN_NPPCROOT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppseeds       (ipy:ipy)                              &
                        ,'MMEAN_NPPSEEDS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppwood        (ipy:ipy)                              &
                        ,'MMEAN_NPPWOOD_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_nppdaily       (ipy:ipy)                              &
                        ,'MMEAN_NPPDAILY_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_A_decomp       (ipy:ipy)                              &
                        ,'MMEAN_A_DECOMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_Af_decomp      (ipy:ipy)                              &
                        ,'MMEAN_AF_DECOMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_co2_residual   (ipy:ipy)                              &
                        ,'MMEAN_CO2_RESIDUAL_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_energy_residual(ipy:ipy)                              &
                        ,'MMEAN_ENERGY_RESIDUAL_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_water_residual (ipy:ipy)                              &
                        ,'MMEAN_WATER_RESIDUAL_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_theiv      (ipy:ipy)                              &
                        ,'MMEAN_ATM_THEIV_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_theta      (ipy:ipy)                              &
                        ,'MMEAN_ATM_THETA_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_temp       (ipy:ipy)                              &
                        ,'MMEAN_ATM_TEMP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_vpdef      (ipy:ipy)                              &
                        ,'MMEAN_ATM_VPDEF_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_shv        (ipy:ipy)                              &
                        ,'MMEAN_ATM_SHV_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_rshort     (ipy:ipy)                              &
                        ,'MMEAN_ATM_RSHORT_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_rshort_diff(ipy:ipy)                              &
                        ,'MMEAN_ATM_RSHORT_DIFF_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_par        (ipy:ipy)                              &
                        ,'MMEAN_ATM_PAR_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_par_diff   (ipy:ipy)                              &
                        ,'MMEAN_ATM_PAR_DIFF_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_rlong      (ipy:ipy)                              &
                        ,'MMEAN_ATM_RLONG_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_vels       (ipy:ipy)                              &
                        ,'MMEAN_ATM_VELS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_rhos       (ipy:ipy)                              &
                        ,'MMEAN_ATM_RHOS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_prss       (ipy:ipy)                              &
                        ,'MMEAN_ATM_PRSS_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_atm_co2        (ipy:ipy)                              &
                        ,'MMEAN_ATM_CO2_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_pcpg           (ipy:ipy)                              &
                        ,'MMEAN_PCPG_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_qpcpg          (ipy:ipy)                              &
                        ,'MMEAN_QPCPG_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_dpcpg          (ipy:ipy)                              &
                        ,'MMEAN_DPCPG_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_gpp            (ipy:ipy)                              &
                        ,'MMSQU_GPP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_npp            (ipy:ipy)                              &
                        ,'MMSQU_NPP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_plresp         (ipy:ipy)                              &
                        ,'MMSQU_PLRESP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_sensible_lc    (ipy:ipy)                              &
                        ,'MMSQU_SENSIBLE_LC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_vapor_lc       (ipy:ipy)                              &
                        ,'MMSQU_VAPOR_LC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_transp         (ipy:ipy)                              &
                        ,'MMSQU_TRANSP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_sensible_wc    (ipy:ipy)                              &
                        ,'MMSQU_SENSIBLE_WC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_vapor_wc       (ipy:ipy)                              &
                        ,'MMSQU_VAPOR_WC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_rh             (ipy:ipy)                              &
                        ,'MMSQU_RH_PY               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_cwd_rh         (ipy:ipy)                              &
                        ,'MMSQU_CWD_RH_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_nep            (ipy:ipy)                              &
                        ,'MMSQU_NEP_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_rlongup        (ipy:ipy)                              &
                        ,'MMSQU_RLONGUP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_parup          (ipy:ipy)                              &
                        ,'MMSQU_PARUP_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_nirup          (ipy:ipy)                              &
                        ,'MMSQU_NIRUP_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_rshortup       (ipy:ipy)                              &
                        ,'MMSQU_RSHORTUP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_rnet           (ipy:ipy)                              &
                        ,'MMSQU_RNET_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_albedo         (ipy:ipy)                              &
                        ,'MMSQU_ALBEDO_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_ustar          (ipy:ipy)                              &
                        ,'MMSQU_USTAR_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_carbon_ac      (ipy:ipy)                              &
                        ,'MMSQU_CARBON_AC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_carbon_st      (ipy:ipy)                              &
                        ,'MMSQU_CARBON_ST_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_vapor_gc       (ipy:ipy)                              &
                        ,'MMSQU_VAPOR_GC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_vapor_ac       (ipy:ipy)                              &
                        ,'MMSQU_VAPOR_AC_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_sensible_gc    (ipy:ipy)                              &
                        ,'MMSQU_SENSIBLE_GC_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmsqu_sensible_ac    (ipy:ipy)                              &
                        ,'MMSQU_SENSIBLE_AC_PY      ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_p11mmean
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_p12(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (nzg; npolygons).                                     !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims (1) = int(nzg,8)
   memsize (1) = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8

   globdims(2) = int(cgrid%npolygons_global,8)
   chnkdims(2) = 1_8
   chnkoffs(2) = int(py_index - 1,8)
   memdims (2) = 1_8
   memsize (2) = 1_8
   memoffs (2) = 0_8

   if (writing_long) then
      call hdf_getslab_r(cgrid%dmean_soil_energy(:,ipy)                                    &
                        ,'DMEAN_SOIL_ENERGY_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_soil_mstpot(:,ipy)                                    &
                        ,'DMEAN_SOIL_MSTPOT_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_soil_water (:,ipy)                                    &
                        ,'DMEAN_SOIL_WATER_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_soil_temp  (:,ipy)                                    &
                        ,'DMEAN_SOIL_TEMP_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_soil_fliq  (:,ipy)                                    &
                        ,'DMEAN_SOIL_FLIQ_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_smoist_gg  (:,ipy)                                    &
                        ,'DMEAN_SMOIST_GG_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_transloss  (:,ipy)                                    &
                        ,'DMEAN_TRANSLOSS_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%dmean_sensible_gg(:,ipy)                                    &
                        ,'DMEAN_SENSIBLE_GG_PY ',dsetrank,iparallel,.false.,foundvar)
   end if
   if (writing_eorq) then
      call hdf_getslab_r(cgrid%mmean_soil_energy(:,ipy)                                    &
                        ,'MMEAN_SOIL_ENERGY_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_soil_mstpot(:,ipy)                                    &
                        ,'MMEAN_SOIL_MSTPOT_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_soil_water (:,ipy)                                    &
                        ,'MMEAN_SOIL_WATER_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_soil_temp  (:,ipy)                                    &
                        ,'MMEAN_SOIL_TEMP_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_soil_fliq  (:,ipy)                                    &
                        ,'MMEAN_SOIL_FLIQ_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_smoist_gg  (:,ipy)                                    &
                        ,'MMEAN_SMOIST_GG_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_transloss  (:,ipy)                                    &
                        ,'MMEAN_TRANSLOSS_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_sensible_gg(:,ipy)                                    &
                        ,'MMEAN_SENSIBLE_GG_PY ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_p12
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_m11(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (ndcycle; npolygons).                                 !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(ndcycle,8)
   chnkdims(1) = int(ndcycle,8)
   memdims (1) = int(ndcycle,8)
   memsize (1) = int(ndcycle,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8

   globdims(2) = int(cgrid%npolygons_global,8)
   chnkdims(2) = 1_8
   chnkoffs(2) = int(py_index - 1,8)
   memdims (2) = 1_8
   memsize (2) = 1_8
   memoffs (2) = 0_8

   if (writing_dcyc) then
      call hdf_getslab_r(cgrid%qmean_gpp            (:,ipy)                                &
                        ,'QMEAN_GPP_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_npp            (:,ipy)                                &
                        ,'QMEAN_NPP_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_resp      (:,ipy)                                &
                        ,'QMEAN_LEAF_RESP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_root_resp      (:,ipy)                                &
                        ,'QMEAN_ROOT_RESP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_growth_resp(:,ipy)                               &
                        ,'QMEAN_LEAF_GROWTH_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_root_growth_resp(:,ipy)                               &
                        ,'QMEAN_ROOT_GROWTH_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sapa_growth_resp(:,ipy)                               &
                        ,'QMEAN_SAPA_GROWTH_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sapb_growth_resp(:,ipy)                               &
                        ,'QMEAN_SAPB_GROWTH_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_storage_resp(:,ipy)                              &
                        ,'QMEAN_LEAF_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_root_storage_resp(:,ipy)                              &
                        ,'QMEAN_ROOT_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sapa_storage_resp(:,ipy)                              &
                        ,'QMEAN_SAPA_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sapb_storage_resp(:,ipy)                              &
                        ,'QMEAN_SAPB_STORAGE_RESP_PY',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_plresp         (:,ipy)                                &
                        ,'QMEAN_PLRESP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_energy    (:,ipy)                                &
                        ,'QMEAN_LEAF_ENERGY_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_water     (:,ipy)                                &
                        ,'QMEAN_LEAF_WATER_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_hcap      (:,ipy)                                &
                        ,'QMEAN_LEAF_HCAP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_vpdef     (:,ipy)                                &
                        ,'QMEAN_LEAF_VPDEF_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_temp      (:,ipy)                                &
                        ,'QMEAN_LEAF_TEMP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_fliq      (:,ipy)                                &
                        ,'QMEAN_LEAF_FLIQ_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_gsw       (:,ipy)                                &
                        ,'QMEAN_LEAF_GSW_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_leaf_gbw       (:,ipy)                                &
                        ,'QMEAN_LEAF_GBW_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wood_energy    (:,ipy)                                &
                        ,'QMEAN_WOOD_ENERGY_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wood_water     (:,ipy)                                &
                        ,'QMEAN_WOOD_WATER_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wood_hcap      (:,ipy)                                &
                        ,'QMEAN_WOOD_HCAP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wood_temp      (:,ipy)                                &
                        ,'QMEAN_WOOD_TEMP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wood_fliq      (:,ipy)                                &
                        ,'QMEAN_WOOD_FLIQ_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wood_gbw       (:,ipy)                                &
                        ,'QMEAN_WOOD_GBW_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_fs_open        (:,ipy)                                &
                        ,'QMEAN_FS_OPEN_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_fsw            (:,ipy)                                &
                        ,'QMEAN_FSW_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_fsn            (:,ipy)                                &
                        ,'QMEAN_FSN_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_a_open         (:,ipy)                                &
                        ,'QMEAN_A_OPEN_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_a_closed       (:,ipy)                                &
                        ,'QMEAN_A_CLOSED_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_a_net          (:,ipy)                                &
                        ,'QMEAN_A_NET_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_a_light        (:,ipy)                                &
                        ,'QMEAN_A_LIGHT_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_a_rubp         (:,ipy)                                &
                        ,'QMEAN_A_RUBP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_a_co2          (:,ipy)                                &
                        ,'QMEAN_A_CO2_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_psi_open       (:,ipy)                                &
                        ,'QMEAN_PSI_OPEN_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_psi_closed     (:,ipy)                                &
                        ,'QMEAN_PSI_CLOSED_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_water_supply   (:,ipy)                                &
                        ,'QMEAN_WATER_SUPPLY_PY    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_par_l          (:,ipy)                                &
                        ,'QMEAN_PAR_L_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_par_l_beam     (:,ipy)                                &
                        ,'QMEAN_PAR_L_BEAM_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_par_l_diff     (:,ipy)                                &
                        ,'QMEAN_PAR_L_DIFF_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rshort_l       (:,ipy)                                &
                        ,'QMEAN_RSHORT_L_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rlong_l        (:,ipy)                                &
                        ,'QMEAN_RLONG_L_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sensible_lc    (:,ipy)                                &
                        ,'QMEAN_SENSIBLE_LC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_vapor_lc       (:,ipy)                                &
                        ,'QMEAN_VAPOR_LC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_transp         (:,ipy)                                &
                        ,'QMEAN_TRANSP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_intercepted_al (:,ipy)                                &
                        ,'QMEAN_INTERCEPTED_AL_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wshed_lg       (:,ipy)                                &
                        ,'QMEAN_WSHED_LG_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rshort_w       (:,ipy)                                &
                        ,'QMEAN_RSHORT_W_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rlong_w        (:,ipy)                                &
                        ,'QMEAN_RLONG_W_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sensible_wc    (:,ipy)                                &
                        ,'QMEAN_SENSIBLE_WC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_vapor_wc       (:,ipy)                                &
                        ,'QMEAN_VAPOR_WC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_intercepted_aw (:,ipy)                                &
                        ,'QMEAN_INTERCEPTED_AW_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_wshed_wg       (:,ipy)                                &
                        ,'QMEAN_WSHED_WG_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rh             (:,ipy)                                &
                        ,'QMEAN_RH_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_cwd_rh         (:,ipy)                                &
                        ,'QMEAN_CWD_RH_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_nep            (:,ipy)                                &
                        ,'QMEAN_NEP_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rk4step        (:,ipy)                                &
                        ,'QMEAN_RK4STEP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_available_water(:,ipy)                                &
                        ,'QMEAN_AVAILABLE_WATER_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_theiv      (:,ipy)                                &
                        ,'QMEAN_CAN_THEIV_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_theta      (:,ipy)                                &
                        ,'QMEAN_CAN_THETA_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_vpdef      (:,ipy)                                &
                        ,'QMEAN_CAN_VPDEF_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_temp       (:,ipy)                                &
                        ,'QMEAN_CAN_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_shv        (:,ipy)                                &
                        ,'QMEAN_CAN_SHV_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_co2        (:,ipy)                                &
                        ,'QMEAN_CAN_CO2_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_rhos       (:,ipy)                                &
                        ,'QMEAN_CAN_RHOS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_prss       (:,ipy)                                &
                        ,'QMEAN_CAN_PRSS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_gnd_temp       (:,ipy)                                &
                        ,'QMEAN_GND_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_gnd_shv        (:,ipy)                                &
                        ,'QMEAN_GND_SHV_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_can_ggnd       (:,ipy)                                &
                        ,'QMEAN_CAN_GGND_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sfcw_depth     (:,ipy)                                &
                        ,'QMEAN_SFCW_DEPTH_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sfcw_energy    (:,ipy)                                &
                        ,'QMEAN_SFCW_ENERGY_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sfcw_mass      (:,ipy)                                &
                        ,'QMEAN_SFCW_MASS_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sfcw_temp      (:,ipy)                                &
                        ,'QMEAN_SFCW_TEMP_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sfcw_fliq      (:,ipy)                                &
                        ,'QMEAN_SFCW_FLIQ_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rshort_gnd     (:,ipy)                                &
                        ,'QMEAN_RSHORT_GND_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_par_gnd        (:,ipy)                                &
                        ,'QMEAN_PAR_GND_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rlong_gnd      (:,ipy)                                &
                        ,'QMEAN_RLONG_GND_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rlongup        (:,ipy)                                &
                        ,'QMEAN_RLONGUP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_parup          (:,ipy)                                &
                        ,'QMEAN_PARUP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_nirup          (:,ipy)                                &
                        ,'QMEAN_NIRUP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rshortup       (:,ipy)                                &
                        ,'QMEAN_RSHORTUP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rnet           (:,ipy)                                &
                        ,'QMEAN_RNET_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_albedo         (:,ipy)                                &
                        ,'QMEAN_ALBEDO_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_albedo_par     (:,ipy)                                &
                        ,'QMEAN_ALBEDO_PAR_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_albedo_nir     (:,ipy)                                &
                        ,'QMEAN_ALBEDO_NIR_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_rlong_albedo   (:,ipy)                                &
                        ,'QMEAN_RLONG_ALBEDO_PY    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_ustar          (:,ipy)                                &
                        ,'QMEAN_USTAR_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_tstar          (:,ipy)                                &
                        ,'QMEAN_TSTAR_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_qstar          (:,ipy)                                &
                        ,'QMEAN_QSTAR_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_cstar          (:,ipy)                                &
                        ,'QMEAN_CSTAR_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_carbon_ac      (:,ipy)                                &
                        ,'QMEAN_CARBON_AC_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_carbon_st      (:,ipy)                                &
                        ,'QMEAN_CARBON_ST_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_vapor_gc       (:,ipy)                                &
                        ,'QMEAN_VAPOR_GC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_vapor_ac       (:,ipy)                                &
                        ,'QMEAN_VAPOR_AC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_throughfall    (:,ipy)                                &
                        ,'QMEAN_THROUGHFALL_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_runoff         (:,ipy)                                &
                        ,'QMEAN_RUNOFF_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_drainage       (:,ipy)                                &
                        ,'QMEAN_DRAINAGE_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sensible_gc    (:,ipy)                                &
                        ,'QMEAN_SENSIBLE_GC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sensible_ac    (:,ipy)                                &
                        ,'QMEAN_SENSIBLE_AC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_qthroughfall   (:,ipy)                                &
                        ,'QMEAN_QTHROUGHFALL_PY    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_qrunoff        (:,ipy)                                &
                        ,'QMEAN_QRUNOFF_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_qdrainage      (:,ipy)                                &
                        ,'QMEAN_QDRAINAGE_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_theiv      (:,ipy)                                &
                        ,'QMEAN_ATM_THEIV_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_theta      (:,ipy)                                &
                        ,'QMEAN_ATM_THETA_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_temp       (:,ipy)                                &
                        ,'QMEAN_ATM_TEMP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_vpdef      (:,ipy)                                &
                        ,'QMEAN_ATM_VPDEF_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_shv        (:,ipy)                                &
                        ,'QMEAN_ATM_SHV_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_rshort     (:,ipy)                                &
                        ,'QMEAN_ATM_RSHORT_PY      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_rshort_diff(:,ipy)                                &
                        ,'QMEAN_ATM_RSHORT_DIFF_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_par        (:,ipy)                                &
                        ,'QMEAN_ATM_PAR_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_par_diff   (:,ipy)                                &
                        ,'QMEAN_ATM_PAR_DIFF_PY    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_rlong      (:,ipy)                                &
                        ,'QMEAN_ATM_RLONG_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_vels       (:,ipy)                                &
                        ,'QMEAN_ATM_VELS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_rhos       (:,ipy)                                &
                        ,'QMEAN_ATM_RHOS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_prss       (:,ipy)                                &
                        ,'QMEAN_ATM_PRSS_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_atm_co2        (:,ipy)                                &
                        ,'QMEAN_ATM_CO2_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_pcpg           (:,ipy)                                &
                        ,'QMEAN_PCPG_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_qpcpg          (:,ipy)                                &
                        ,'QMEAN_QPCPG_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_dpcpg          (:,ipy)                                &
                        ,'QMEAN_DPCPG_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_gpp            (:,ipy)                                &
                        ,'QMSQU_GPP_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_npp            (:,ipy)                                &
                        ,'QMSQU_NPP_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_plresp         (:,ipy)                                &
                        ,'QMSQU_PLRESP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_sensible_lc    (:,ipy)                                &
                        ,'QMSQU_SENSIBLE_LC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_vapor_lc       (:,ipy)                                &
                        ,'QMSQU_VAPOR_LC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_transp         (:,ipy)                                &
                        ,'QMSQU_TRANSP_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_sensible_wc    (:,ipy)                                &
                        ,'QMSQU_SENSIBLE_WC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_vapor_wc       (:,ipy)                                &
                        ,'QMSQU_VAPOR_WC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_rh             (:,ipy)                                &
                        ,'QMSQU_RH_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_cwd_rh         (:,ipy)                                &
                        ,'QMSQU_CWD_RH_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_nep            (:,ipy)                                &
                        ,'QMSQU_NEP_PY             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_rlongup        (:,ipy)                                &
                        ,'QMSQU_RLONGUP_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_parup          (:,ipy)                                &
                        ,'QMSQU_PARUP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_nirup          (:,ipy)                                &
                        ,'QMSQU_NIRUP_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_rshortup       (:,ipy)                                &
                        ,'QMSQU_RSHORTUP_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_rnet           (:,ipy)                                &
                        ,'QMSQU_RNET_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_albedo         (:,ipy)                                &
                        ,'QMSQU_ALBEDO_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_ustar          (:,ipy)                                &
                        ,'QMSQU_USTAR_PY           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_carbon_ac      (:,ipy)                                &
                        ,'QMSQU_CARBON_AC_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_carbon_st      (:,ipy)                                &
                        ,'QMSQU_CARBON_ST_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_vapor_gc       (:,ipy)                                &
                        ,'QMSQU_VAPOR_GC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_vapor_ac       (:,ipy)                                &
                        ,'QMSQU_VAPOR_AC_PY        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_sensible_gc    (:,ipy)                                &
                        ,'QMSQU_SENSIBLE_GC_PY     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmsqu_sensible_ac    (:,ipy)                                &
                        ,'QMSQU_SENSIBLE_AC_PY     ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_m11
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_p19(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!







   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (13 months; npolygons).                               !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(13,8)
   chnkdims(1) = int(13,8)
   memdims (1) = int(13,8)
   memsize (1) = int(13,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(cgrid%npolygons_global,8)
   chnkdims(2) = 1_8
   chnkoffs(2) = int(py_index - 1,8)
   memdims (2) = 1_8
   memsize (2) = 1_8
   memoffs (2) = 0_8

   call hdf_getslab_r(cgrid%workload(:,ipy)                                                &
                     ,'WORKLOAD ',dsetrank,iparallel,.true.,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_p19
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_m12(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (nzg;ndcycle;npolygons).                              !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 3
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims (1) = int(nzg,8)
   memsize (1) = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8

   globdims(2) = int(ndcycle,8)
   chnkdims(2) = int(ndcycle,8)
   memdims (2) = int(ndcycle,8)
   memsize (2) = int(ndcycle,8)
   chnkoffs(2) = 0_8
   memoffs (2) = 0_8

   globdims(3)  = int(cgrid%npolygons_global,8)
   chnkdims(3)  = 1_8
   chnkoffs(3)  = int(py_index - 1,8)
   memdims (3)  = 1_8
   memsize (3)  = 1_8
   memoffs (3)  = 0_8

   if (writing_dcyc) then
      call hdf_getslab_r(cgrid%qmean_soil_energy(:,:,ipy)                                  &
                        ,'QMEAN_SOIL_ENERGY_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_soil_mstpot(:,:,ipy)                                  &
                        ,'QMEAN_SOIL_MSTPOT_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_soil_water (:,:,ipy)                                  &
                        ,'QMEAN_SOIL_WATER_PY  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_soil_temp  (:,:,ipy)                                  &
                        ,'QMEAN_SOIL_TEMP_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_soil_fliq  (:,:,ipy)                                  &
                        ,'QMEAN_SOIL_FLIQ_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_smoist_gg  (:,:,ipy)                                  &
                        ,'QMEAN_SMOIST_GG_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_transloss  (:,:,ipy)                                  &
                        ,'QMEAN_TRANSLOSS_PY   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%qmean_sensible_gg(:,:,ipy)                                  &
                        ,'QMEAN_SENSIBLE_GG_PY ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_m12
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid_p146(cgrid,ipy,py_index)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , n_age         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!


   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)   , target      :: cgrid
   integer        , intent(in)  :: ipy
   integer        , intent(in)  :: py_index
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: iparallel
   integer                      :: dsetrank
   logical                      :: foundvar
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (nzg;ndcycle;npolygons).                              !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 3
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims (1) = int(n_pft,8)
   memsize (1) = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8

   globdims(2) = int(n_dbh,8)
   chnkdims(2) = int(n_dbh,8)
   memdims (2) = int(n_dbh,8)
   memsize (2) = int(n_dbh,8)
   chnkoffs(2) = 0_8
   memoffs (2) = 0_8

   globdims(3)  = int(cgrid%npolygons_global,8)
   chnkdims(3)  = 1_8
   chnkoffs(3)  = int(py_index - 1,8)
   memdims (3)  = 1_8
   memsize (3)  = 1_8
   memoffs (3)  = 0_8

   call hdf_getslab_r(cgrid%nplant                   (:,:,ipy)                             &
                     ,'NPLANT_PY                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%agb                      (:,:,ipy)                             &
                     ,'AGB_PY                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%lai                      (:,:,ipy)                             &
                     ,'LAI_PY                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%wai                      (:,:,ipy)                             &
                     ,'WAI_PY                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%basal_area               (:,:,ipy)                             &
                     ,'BASAL_AREA_PY                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bdead                    (:,:,ipy)                             &
                     ,'BDEAD_PY                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%balive                   (:,:,ipy)                             &
                     ,'BALIVE_PY                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bleaf                    (:,:,ipy)                             &
                     ,'BLEAF_PY                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%broot                    (:,:,ipy)                             &
                     ,'BROOT_PY                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bsapwooda                (:,:,ipy)                             &
                     ,'BSAPWOODA_PY                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bsapwoodb                (:,:,ipy)                             &
                     ,'BSAPWOODB_PY                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bseeds                   (:,:,ipy)                             &
                     ,'BSEEDS_PY                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bstorage                 (:,:,ipy)                             &
                     ,'BSTORAGE_PY                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bdead_n                  (:,:,ipy)                             &
                     ,'BDEAD_N_PY                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%balive_n                 (:,:,ipy)                             &
                     ,'BALIVE_N_PY                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bleaf_n                  (:,:,ipy)                             &
                     ,'BLEAF_N_PY                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%broot_n                  (:,:,ipy)                             &
                     ,'BROOT_N_PY                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bsapwooda_n              (:,:,ipy)                             &
                     ,'BSAPWOODA_N_PY               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bsapwoodb_n              (:,:,ipy)                             &
                     ,'BSAPWOODB_N_PY               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bseeds_n                 (:,:,ipy)                             &
                     ,'BSEEDS_N_PY                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%bstorage_n               (:,:,ipy)                             &
                     ,'BSTORAGE_N_PY                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%leaf_maintenance         (:,:,ipy)                             &
                     ,'LEAF_MAINTENANCE_PY          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%root_maintenance         (:,:,ipy)                             &
                     ,'ROOT_MAINTENANCE_PY          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cgrid%leaf_drop                (:,:,ipy)                             &
                     ,'LEAF_DROP_PY                 ',dsetrank,iparallel,.true. ,foundvar)
   if (writing_eorq) then
      call hdf_getslab_r(cgrid%mmean_lai             (:,:,ipy)                             &
                        ,'MMEAN_LAI_PY              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_bleaf           (:,:,ipy)                             &
                        ,'MMEAN_BLEAF_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_broot           (:,:,ipy)                             &
                        ,'MMEAN_BROOT_PY            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_bstorage        (:,:,ipy)                             &
                        ,'MMEAN_BSTORAGE_PY         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_bleaf_n         (:,:,ipy)                             &
                        ,'MMEAN_BLEAF_N_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_broot_n         (:,:,ipy)                             &
                        ,'MMEAN_BROOT_N_PY          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_bstorage_n      (:,:,ipy)                             &
                        ,'MMEAN_BSTORAGE_N_PY       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_maintenance(:,:,ipy)                             &
                        ,'MMEAN_LEAF_MAINTENANCE_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_root_maintenance(:,:,ipy)                             &
                        ,'MMEAN_ROOT_MAINTENANCE_PY ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cgrid%mmean_leaf_drop       (:,:,ipy)                             &
                        ,'MMEAN_LEAF_DROP_PY        ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_grid_p146
!==========================================================================================!
!==========================================================================================!







!==========================================================================================!
!==========================================================================================!
!      This sub-routine loads all site-level variables from the history file.              !
!------------------------------------------------------------------------------------------!
subroutine fill_history_polygon(cpoly,pysi_index,nsites_global,nsites_now,is_burnt)
   use ed_state_vars, only : polygontype   ! ! structure
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_max_dims  , only : n_pft         & ! intent(in)
                           , n_dbh         & ! intent(in)
                           , max_site      & ! intent(in)
                           , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms    , only : file_id       & ! intent(inout)
                           , dset_id       & ! intent(inout)
                           , dspace_id     & ! intent(inout)
                           , plist_id      & ! intent(inout)
                           , globdims      & ! intent(inout)
                           , chnkdims      & ! intent(inout)
                           , chnkoffs      & ! intent(inout)
                           , cnt           & ! intent(inout)
                           , stride        & ! intent(inout)
                           , memdims       & ! intent(inout)
                           , memoffs       & ! intent(inout)
                           , memsize       & ! intent(inout)
                           , datatype_id   ! ! intent(inout)
   use ed_misc_coms , only : ndcycle       & ! intent(in)
                           , writing_long  & ! intent(in)
                           , writing_eorq  & ! intent(in)
                           , writing_dcyc  ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Arguments. ----------------------------------------------------------------------!
   type(polygontype)             , target        :: cpoly
   integer                       , intent(in)    :: pysi_index
   integer                       , intent(in)    :: nsites_global
   integer                       , intent(in)    :: nsites_now
   logical, dimension(nsites_now), intent(inout) :: is_burnt
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: iparallel
   integer                                       :: dsetrank
   integer                                       :: isi
   logical                                       :: foundvar
   integer                                       :: pidx
   integer                                       :: sidx
   integer, dimension(:)         , allocatable   :: nat_dist_type
   real, dimension(:,:,:)        , allocatable   :: tmp_dist_memory
   real, dimension(:,:,:)        , allocatable   :: tmp_dist_rates
   logical                                       :: old_histo
   !----- Local constants. ----------------------------------------------------------------!
   integer               , parameter   :: n_lu_old = 3
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables - Integers.                                           !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(nsites_global ,8)
   chnkdims(1) = int(cpoly%nsites  ,8)
   chnkoffs(1) = int(pysi_index - 1,8)
   memdims (1) = int(cpoly%nsites  ,8)
   memsize (1) = int(cpoly%nsites  ,8)
   memoffs (1) = 0_8
   !call hdf_getslab_i(cpoly%patch_count                                                    &
   !                  ,'PATCH_COUNT             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%sitenum                                                        &
                     ,'SITENUM                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%num_landuse_years                                              &
                     ,'NUM_LANDUSE_YEARS       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%hydro_next                                                     &
                     ,'HYDRO_NEXT              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%hydro_prev                                                     &
                     ,'HYDRO_PREV              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%plantation                                                     &
                     ,'PLANTATION_SI           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%agri_stocking_pft                                              &
                     ,'AGRI_STOCKING_PFT       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpoly%plantation_stocking_pft                                        &
                     ,'PLANTATION_STOCKING_PFT ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !       The _SI extension has been dropped as LSL and soil colour are no longer polygon !
   !  variables.  In case the history is from an old file, they may still exist, so we try !
   !  them first and if it doesn't work we try the one without extension.                  !
   !---------------------------------------------------------------------------------------!
   call hdf_getslab_i(cpoly%lsl                                                            &
                     ,'LSL_SI                  ',dsetrank,iparallel,.false.,foundvar)
   if (.not.foundvar) then
      call hdf_getslab_i(cpoly%lsl                                                         &
                        ,'LSL                  ',dsetrank,iparallel,.true. ,foundvar)
   end if
   call hdf_getslab_i(cpoly%ncol_soil                                                      &
                     ,'NCOL_SOIL_SI            ',dsetrank,iparallel,.false.,foundvar)
   if (.not.foundvar) then
      call hdf_getslab_i(cpoly%ncol_soil                                                   &
                        ,'NCOL_SOIL            ',dsetrank,iparallel,.true. ,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !     Test whether this is a new or an old history.                                     !
   ! MLO. Variable NAT_DIST_TYPE was deleted because now burnt patches are distinguished   !
   !      from tree fall patches.   Thus, if NAT_DIST_TYPE exists in the file, then it     !
   !      must be and old history file.                                                    !
   !---------------------------------------------------------------------------------------!
   allocate(nat_dist_type(cpoly%nsites))
   call hdf_getslab_i(nat_dist_type                                                        &
                     ,'NAT_DIST_TYPE           ',dsetrank,iparallel,.false.,foundvar)
   old_histo = foundvar
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables - Real.                                               !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(nsites_global ,8)
   chnkdims(1) = int(cpoly%nsites  ,8)
   chnkoffs(1) = int(pysi_index - 1,8)
   memdims (1) = int(cpoly%nsites  ,8)
   memsize (1) = int(cpoly%nsites  ,8)
   memoffs (1) = 0_8
   call hdf_getslab_r(cpoly%area                                                           &
                      ,'AREA_SI                    ' ,dsetrank,iparallel,.true. ,foundvar)
   ! call hdf_getslab_r(cpoly%patch_area                                                     &
   !                    ,'PATCH_AREA                 ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%elevation                                                      &
                      ,'ELEVATION                  ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%slope                                                          &
                      ,'SLOPE                      ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%aspect                                                         &
                      ,'ASPECT                     ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%TCI                                                            &
                      ,'TCI                        ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%pptweight                                                      &
                      ,'PPTWEIGHT                  ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%moist_W                                                        &
                      ,'MOIST_W                    ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%moist_f                                                        &
                      ,'MOIST_F                    ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%moist_tau                                                      &
                      ,'MOIST_TAU                  ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%moist_zi                                                       &
                      ,'MOIST_ZI                   ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%baseflow                                                       &
                      ,'BASEFLOW_SI                ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%min_monthly_temp                                               &
                      ,'MIN_MONTHLY_TEMP           ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%agri_stocking_density                                          &
                      ,'AGRI_STOCKING_DENSITY      ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%plantation_stocking_density                                    &
                      ,'PLANTATION_STOCKING_DENSITY' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%primary_harvest_target                                         &
                      ,'PRIMARY_HARVEST_TARGET     ' ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cpoly%secondary_harvest_target                                       &
                      ,'SECONDARY_HARVEST_TARGET   ' ,dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(cpoly%primary_harvest_memory                                         &
                      ,'PRIMARY_HARVEST_MEMORY     ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%secondary_harvest_memory                                       &
                      ,'SECONDARY_HARVEST_MEMORY   ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%ignition_rate                                                  &
                      ,'IGNITION_RATE              ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%rad_avg                                                        &
                      ,'RAD_AVG                    ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%daylight                                                       &
                      ,'DAYLIGHT                   ' ,dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%cosaoi                                                         &
                      ,'COSAOI                     ' ,dsetrank,iparallel,.true. ,foundvar)
   !------ Daily means. -------------------------------------------------------------------!
   if (writing_long) then
      call hdf_getslab_r(cpoly%dmean_atm_theiv                                             &
                        ,'DMEAN_ATM_THEIV_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_theta                                             &
                        ,'DMEAN_ATM_THETA_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_temp                                              &
                        ,'DMEAN_ATM_TEMP_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_vpdef                                             &
                        ,'DMEAN_ATM_VPDEF_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_shv                                               &
                        ,'DMEAN_ATM_SHV_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_rshort                                            &
                        ,'DMEAN_ATM_RSHORT_SI       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_rshort_diff                                       &
                        ,'DMEAN_ATM_RSHORT_DIFF_SI  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_par                                               &
                        ,'DMEAN_ATM_PAR_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_par_diff                                          &
                        ,'DMEAN_ATM_PAR_DIFF_SI     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_rlong                                             &
                        ,'DMEAN_ATM_RLONG_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_vels                                              &
                        ,'DMEAN_ATM_VELS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_rhos                                              &
                        ,'DMEAN_ATM_RHOS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_prss                                              &
                        ,'DMEAN_ATM_PRSS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_atm_co2                                               &
                        ,'DMEAN_ATM_CO2_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_pcpg                                                  &
                        ,'DMEAN_PCPG_SI             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_qpcpg                                                 &
                        ,'DMEAN_QPCPG_SI            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%dmean_dpcpg                                                 &
                        ,'DMEAN_DPCPG_SI            ',dsetrank,iparallel,.false.,foundvar)
   end if
   !------ Monthly means. -----------------------------------------------------------------!
   if (writing_eorq) then
      call hdf_getslab_r(cpoly%mmean_atm_theiv                                             &
                        ,'MMEAN_ATM_THEIV_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_theta                                             &
                        ,'MMEAN_ATM_THETA_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_temp                                              &
                        ,'MMEAN_ATM_TEMP_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_vpdef                                             &
                        ,'MMEAN_ATM_VPDEF_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_shv                                               &
                        ,'MMEAN_ATM_SHV_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_rshort                                            &
                        ,'MMEAN_ATM_RSHORT_SI       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_rshort_diff                                       &
                        ,'MMEAN_ATM_RSHORT_DIFF_SI  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_par                                               &
                        ,'MMEAN_ATM_PAR_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_par_diff                                          &
                        ,'MMEAN_ATM_PAR_DIFF_SI     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_rlong                                             &
                        ,'MMEAN_ATM_RLONG_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_vels                                              &
                        ,'MMEAN_ATM_VELS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_rhos                                              &
                        ,'MMEAN_ATM_RHOS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_prss                                              &
                        ,'MMEAN_ATM_PRSS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_atm_co2                                               &
                        ,'MMEAN_ATM_CO2_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_pcpg                                                  &
                        ,'MMEAN_PCPG_SI             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_qpcpg                                                 &
                        ,'MMEAN_QPCPG_SI            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%mmean_dpcpg                                                 &
                        ,'MMEAN_DPCPG_SI            ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!







   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (ndcycle;nsites).                                     !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(ndcycle,8)
   chnkdims(1) = int(ndcycle,8)
   memdims (1) = int(ndcycle,8)
   memsize (1) = int(ndcycle,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(nsites_global ,8)
   chnkdims(2) = int(cpoly%nsites  ,8)
   chnkoffs(2) = int(pysi_index - 1,8)
   memdims (2) = int(cpoly%nsites  ,8)
   memsize (2) = int(cpoly%nsites  ,8)
   memoffs (2) = 0_8
   if (writing_dcyc) then
      call hdf_getslab_r(cpoly%qmean_atm_theiv                                             &
                        ,'QMEAN_ATM_THEIV_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_theta                                             &
                        ,'QMEAN_ATM_THETA_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_temp                                              &
                        ,'QMEAN_ATM_TEMP_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_vpdef                                             &
                        ,'QMEAN_ATM_VPDEF_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_shv                                               &
                        ,'QMEAN_ATM_SHV_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_rshort                                            &
                        ,'QMEAN_ATM_RSHORT_SI       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_rshort_diff                                       &
                        ,'QMEAN_ATM_RSHORT_DIFF_SI  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_par                                               &
                        ,'QMEAN_ATM_PAR_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_par_diff                                          &
                        ,'QMEAN_ATM_PAR_DIFF_SI     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_rlong                                             &
                        ,'QMEAN_ATM_RLONG_SI        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_vels                                              &
                        ,'QMEAN_ATM_VELS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_rhos                                              &
                        ,'QMEAN_ATM_RHOS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_prss                                              &
                        ,'QMEAN_ATM_PRSS_SI         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_atm_co2                                               &
                        ,'QMEAN_ATM_CO2_SI          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_pcpg                                                  &
                        ,'QMEAN_PCPG_SI             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_qpcpg                                                 &
                        ,'QMEAN_QPCPG_SI            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpoly%qmean_dpcpg                                                 &
                        ,'QMEAN_DPCPG_SI            ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (n_pft;nsites).                                       !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims (1) = int(n_pft,8)
   memsize (1) = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(nsites_global,8)
   chnkdims(2) = int(cpoly%nsites,8)
   chnkoffs(2) = int(pysi_index - 1,8)
   memdims (2) = int(cpoly%nsites,8)
   memsize (2) = int(cpoly%nsites,8)
   memoffs (2) = 0_8

   call hdf_getslab_r(cpoly%green_leaf_factor                                              &
                     ,'GREEN_LEAF_FACTOR  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%leaf_aging_factor                                              &
                     ,'LEAF_AGING_FACTOR  ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (nzg;nsites).                                         !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims (1) = int(nzg,8)
   memsize (1) = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(nsites_global ,8)
   chnkdims(2) = int(cpoly%nsites  ,8)
   chnkoffs(2) = int(pysi_index - 1,8)
   memdims (2) = int(cpoly%nsites  ,8)
   memsize (2) = int(cpoly%nsites  ,8)
   memoffs (2) = 0_8

   !---------------------------------------------------------------------------------------!
   !       The _SI extension has been dropped as soil texture is no longer a polygon       !
   !  variable.  In case the history is from an old file, they may still exist, so we try  !
   !  it first and if it doesn't work we try the one without extension.                    !
   !---------------------------------------------------------------------------------------!
   call hdf_getslab_i(cpoly%ntext_soil                                                     &
                     ,'NTEXT_SOIL_SI ',dsetrank,iparallel,.false.,foundvar)
   if (.not. foundvar) then
      call hdf_getslab_i(cpoly%ntext_soil                                                  &
                        ,'NTEXT_SOIL ',dsetrank,iparallel,.true.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (n_pft;nsites).                                       !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = 12_8
   chnkdims(1) = 12_8
   memdims (1) = 12_8
   memsize (1) = 12_8
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(nsites_global ,8)
   chnkdims(2) = int(cpoly%nsites  ,8)
   chnkoffs(2) = int(pysi_index - 1,8)
   memdims (2) = int(cpoly%nsites  ,8)
   memsize (2) = int(cpoly%nsites  ,8)
   memoffs (2) = 0_8
   call hdf_getslab_r(cpoly%lambda_fire                                                    &
                     ,'LAMBDA_FIRE ',dsetrank,iparallel,.true.,foundvar)
   call hdf_getslab_r(cpoly%avg_monthly_pcpg                                               &
                     ,'AVG_MONTHLY_PCPG ',dsetrank,iparallel,.true.,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (n_dist_types;n_dist_types;nsites).  Because the      !
   ! number of land use types has changed, we must check whether this history file is old  !
   ! or new.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (old_histo) then
      !------------------------------------------------------------------------------------!
      !     Old history file, read the information to a temporary array.                   !
      !------------------------------------------------------------------------------------!
      dsetrank       = 3
      globdims(1:2)  = int(n_lu_old,8)
      chnkdims(1:2)  = int(n_lu_old,8)
      memdims (1:2)  = int(n_lu_old,8)
      memsize (1:2)  = int(n_lu_old,8)
      chnkoffs(1:2)  = 0_8
      memoffs (1:2)  = 0_8
      globdims(3)    = int(nsites_global ,8)
      chnkdims(3)    = int(cpoly%nsites  ,8)
      chnkoffs(3)    = int(pysi_index - 1,8)
      memdims(3)     = int(cpoly%nsites  ,8)
      memsize(3)     = int(cpoly%nsites  ,8)
      memoffs(3)     = 0_8
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Allocate temporary variable, we must translate the transition rates to the     !
      ! new matrix.                                                                        !
      !------------------------------------------------------------------------------------!
      allocate(tmp_dist_memory(n_lu_old,n_lu_old,cpoly%nsites))
      allocate(tmp_dist_rates (n_lu_old,n_lu_old,cpoly%nsites))
      cpoly%disturbance_memory(:,:,:) = 0.0
      cpoly%disturbance_rates (:,:,:) = 0.0

      call hdf_getslab_r(tmp_dist_memory                                                   &
                        ,'DISTURBANCE_MEMORY ',dsetrank,iparallel,.true.,foundvar)
      !------------------------------------------------------------------------------------!
      !       The _SI extension has been dropped as disturbance rate is no longer a        !
      ! polygon variable.  In case the history is from an old file, they may still exist,  !
      ! so we try it first and if it doesn't work we try the one without extension.        !
      !------------------------------------------------------------------------------------!
      call hdf_getslab_r(tmp_dist_rates                                                    &
                        ,'DISTURBANCE_RATES_SI ',dsetrank,iparallel,.false.,foundvar)
      if (.not. foundvar) then
         call hdf_getslab_r(tmp_dist_rates                                                 &
                           ,'DISTURBANCE_RATES ',dsetrank,iparallel,.true.,foundvar)
      end if
      !------------------------------------------------------------------------------------!


      do isi=1,cpoly%nsites
         !------ Save information on whether the site has been burnt. ---------------------!
         is_burnt(isi) = nat_dist_type(isi) == 1
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Find the indices for disturbance transitions based on the last fire regime  !
         ! and whether it is a forest plantation of an explored stand.  The translation    !
         ! won't be perfect, because we cannot retrieve all information needed.            !
         !---------------------------------------------------------------------------------!
         !----- "Primary forest index", tree fall or burnt. -------------------------------!
         select case (nat_dist_type(isi))
         case (0)
            pidx = 3
         case (1)
            pidx = 4
         end select
         !----- "Secondary forest index", logged or plantation. ---------------------------!
         select case (cpoly%plantation(isi))
         case (0)
            sidx = 6
         case (1)
            sidx = 2
         end select
         !---------------------------------------------------------------------------------!


         !----- Translate the memory transition matrix. -----------------------------------!
         cpoly%disturbance_memory(   1,   1,isi) = tmp_dist_memory(1,1,isi)
         cpoly%disturbance_memory(   1,sidx,isi) = tmp_dist_memory(1,2,isi)
         cpoly%disturbance_memory(   1,pidx,isi) = tmp_dist_memory(1,3,isi)
         cpoly%disturbance_memory(   5,   1,isi) = tmp_dist_memory(2,1,isi)
         cpoly%disturbance_memory(sidx,sidx,isi) = tmp_dist_memory(2,2,isi)
         cpoly%disturbance_memory(sidx,pidx,isi) = tmp_dist_memory(2,3,isi)
         cpoly%disturbance_memory(pidx,sidx,isi) = tmp_dist_memory(3,2,isi)
         cpoly%disturbance_memory(pidx,pidx,isi) = tmp_dist_memory(3,3,isi)
         !---------------------------------------------------------------------------------!




         !----- Translate the current transition matrix. ----------------------------------!
         cpoly%disturbance_rates (   1,   1,isi) = tmp_dist_rates (1,1,isi)
         cpoly%disturbance_rates (   1,sidx,isi) = tmp_dist_rates (1,2,isi)
         cpoly%disturbance_rates (   1,pidx,isi) = tmp_dist_rates (1,3,isi)
         cpoly%disturbance_rates (   5,   1,isi) = tmp_dist_rates (2,1,isi)
         cpoly%disturbance_rates (sidx,sidx,isi) = tmp_dist_rates (2,2,isi)
         cpoly%disturbance_rates (sidx,pidx,isi) = tmp_dist_rates (2,3,isi)
         cpoly%disturbance_rates (pidx,sidx,isi) = tmp_dist_rates (3,2,isi)
         cpoly%disturbance_rates (pidx,pidx,isi) = tmp_dist_rates (3,3,isi)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !----- Free memory. -----------------------------------------------------------------!
      deallocate(tmp_dist_memory)
      deallocate(tmp_dist_rates )
      !------------------------------------------------------------------------------------!

   else
      !------------------------------------------------------------------------------------!
      !     Compatible history file, just read the information.                            !
      !------------------------------------------------------------------------------------!
      dsetrank       = 3
      globdims(1:2)  = int(n_dist_types,8)
      chnkdims(1:2)  = int(n_dist_types,8)
      memdims (1:2)  = int(n_dist_types,8)
      memsize (1:2)  = int(n_dist_types,8)
      chnkoffs(1:2)  = 0_8
      memoffs (1:2)  = 0_8
      globdims(3)    = int(nsites_global ,8)
      chnkdims(3)    = int(cpoly%nsites  ,8)
      chnkoffs(3)    = int(pysi_index - 1,8)
      memdims(3)     = int(cpoly%nsites  ,8)
      memsize(3)     = int(cpoly%nsites  ,8)
      memoffs(3)     = 0_8
      !------------------------------------------------------------------------------------!

      call hdf_getslab_r(cpoly%disturbance_memory                                          &
                        ,'DISTURBANCE_MEMORY ',dsetrank,iparallel,.true.,foundvar)
      !------------------------------------------------------------------------------------!
      !       The _SI extension has been dropped as disturbance rate is no longer a        !
      ! polygon variable.  In case the history is from an old file, they may still exist,  !
      ! so we try it first and if it doesn't work we try the one without extension.        !
      !------------------------------------------------------------------------------------!
      call hdf_getslab_r(cpoly%disturbance_rates                                           &
                        ,'DISTURBANCE_RATES_SI ',dsetrank,iparallel,.false.,foundvar)
      if (.not. foundvar) then
         call hdf_getslab_r(cpoly%disturbance_rates                                        &
                           ,'DISTURBANCE_RATES ',dsetrank,iparallel,.true.,foundvar)
      end if
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (n_pft;n_dbh_types;nsites).                           !
   !---------------------------------------------------------------------------------------!

   dsetrank    = 3
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims (1) = int(n_pft,8)
   memsize (1) = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(n_dbh,8)
   chnkdims(2) = int(n_dbh,8)
   memdims (2) = int(n_dbh,8)
   memsize (2) = int(n_dbh,8)
   chnkoffs(2) = 0_8
   memoffs (2) = 0_8
   globdims(3) = int(nsites_global ,8)
   chnkdims(3) = int(cpoly%nsites  ,8)
   chnkoffs(3) = int(pysi_index - 1,8)
   memdims (3) = int(cpoly%nsites  ,8)
   memsize (3) = int(cpoly%nsites  ,8)
   memoffs (3) = 0_8

   call hdf_getslab_r(cpoly%basal_area                                                     &
                     ,'BASAL_AREA_SI     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%agb                                                            &
                     ,'AGB_SI            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%basal_area_growth                                              &
                     ,'BASAL_AREA_GROWTH ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%agb_growth                                                     &
                     ,'AGB_GROWTH        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%basal_area_mort                                                &
                     ,'BASAL_AREA_MORT   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%basal_area_cut                                                 &
                     ,'BASAL_AREA_CUT    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%agb_mort                                                       &
                     ,'AGB_MORT          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpoly%agb_cut                                                        &
                     ,'AGB_CUT           ',dsetrank,iparallel,.true. ,foundvar)


   !----- Free memory. --------------------------------------------------------------------!
   deallocate(nat_dist_type)
   !---------------------------------------------------------------------------------------!
   return
end subroutine fill_history_polygon
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_site(csite,sipa_index,npatches_global,is_burnt)
   use ed_state_vars      , only : sitetype      ! ! structure
   use grid_coms          , only : nzg           & ! intent(in)
                                 , nzs           ! ! intent(in)
   use ed_max_dims        , only : n_pft         & ! intent(in)
                                 , n_dbh         & ! intent(in)
                                 , max_site      & ! intent(in)
                                 , n_dist_types  ! ! intent(in)
   use hdf5
   use hdf5_coms          , only : file_id       & ! intent(inout)
                                 , dset_id       & ! intent(inout)
                                 , dspace_id     & ! intent(inout)
                                 , plist_id      & ! intent(inout)
                                 , globdims      & ! intent(inout)
                                 , chnkdims      & ! intent(inout)
                                 , chnkoffs      & ! intent(inout)
                                 , cnt           & ! intent(inout)
                                 , stride        & ! intent(inout)
                                 , memdims       & ! intent(inout)
                                 , memoffs       & ! intent(inout)
                                 , memsize       & ! intent(inout)
                                 , datatype_id   & ! intent(inout)
                                 , setsize       ! ! intent(inout)
   use ed_misc_coms       , only : ndcycle       & ! intent(in)
                                 , writing_long  & ! intent(in)
                                 , writing_eorq  & ! intent(in)
                                 , writing_dcyc  ! ! intent(in)
   use fusion_fission_coms, only : ff_nhgt       ! ! intent(in)
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)             , target       :: csite
   integer                    , intent(in)   :: sipa_index
   integer                    , intent(in)   :: npatches_global
   logical                    , intent(in)   :: is_burnt
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: iparallel
   integer                                   :: dsetrank
   logical                                   :: foundvar
   integer                                   :: ipa
   ! real(kind=8), dimension(:,:), allocatable :: buff
   integer     , dimension(:)  , allocatable :: plantation
   !---------------------------------------------------------------------------------------!



   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables - Integers.                                           !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(npatches_global,8)
   chnkdims(1) = int(csite%npatches,8)
   chnkoffs(1) = int(sipa_index - 1,8)
   memdims (1)  = int(csite%npatches,8)
   memsize (1)  = int(csite%npatches,8)
   memoffs (1)  = 0_8
   call hdf_getslab_i(csite%dist_type                                                      &
                     ,'DIST_TYPE                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(csite%nlev_sfcwater                                                  &
                     ,'NLEV_SFCWATER               ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Check whether this is an old history file (back when there were only 3 types).   !
   ! If this is the case, we must translate the indices to the new disturbance table.      !
   ! This is not going to be perfect because the number of disturbance types increased,    !
   ! but this should be rarely needed, as it is not a good idea to switch versions then    !
   ! continue to run the model using binary-history.                                       !
   !---------------------------------------------------------------------------------------!
   allocate(plantation(csite%npatches))
   call hdf_getslab_i(plantation                                                           &
                     ,'PLANTATION                  ',dsetrank,iparallel,.false. ,foundvar)

   !----- Old history files have variable plantation.  If so, check disturbance types. ----!
   if (foundvar) then
      !----- Go through all patches and decide whether to change the disturbance type. ----!
      do ipa=1,csite%npatches
         select case(csite%dist_type(ipa))
         case (2)
            if (plantation(ipa) == 0) csite%dist_type(ipa) = 6
         case (3)
            if (is_burnt            ) csite%dist_type(ipa) = 4
         end select
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   deallocate(plantation)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables - Real.                                               !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 1
   globdims(1) = int(npatches_global,8)
   chnkdims(1) = int(csite%npatches,8)
   chnkoffs(1) = int(sipa_index - 1,8)
   memdims (1)  = int(csite%npatches,8)
   memsize (1)  = int(csite%npatches,8)
   memoffs (1)  = 0_8
   call hdf_getslab_r(csite%area                                                           &
                     ,'AREA                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%age                                                            &
                     ,'AGE                         ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%fast_soil_C                                                    &
                     ,'FAST_SOIL_C                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%slow_soil_C                                                    &
                     ,'SLOW_SOIL_C                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%structural_soil_C                                              &
                     ,'STRUCTURAL_SOIL_C           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%structural_soil_L                                              &
                     ,'STRUCTURAL_SOIL_L           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%mineralized_soil_N                                             &
                     ,'MINERALIZED_SOIL_N          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%fast_soil_N                                                    &
                     ,'FAST_SOIL_N                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%sum_dgd                                                        &
                     ,'SUM_DGD                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%sum_chd                                                        &
                     ,'SUM_CHD                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_theiv                                                      &
                     ,'CAN_THEIV                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_vpdef                                                      &
                     ,'CAN_VPDEF                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_temp                                                       &
                     ,'CAN_TEMP                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_temp_pv                                                    &
                     ,'CAN_TEMP_PV                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_shv                                                        &
                     ,'CAN_SHV                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_co2                                                        &
                     ,'CAN_CO2                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_rhos                                                       &
                     ,'CAN_RHOS                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_prss                                                       &
                     ,'CAN_PRSS                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_theta                                                      &
                     ,'CAN_THETA                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%can_depth                                                      &
                     ,'CAN_DEPTH                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%opencan_frac                                                   &
                     ,'OPENCAN_FRAC                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ggbare                                                         &
                     ,'GGBARE                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ggveg                                                          &
                     ,'GGVEG                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ggnet                                                          &
                     ,'GGNET                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ggsoil                                                         &
                     ,'GGSOIL                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ground_shv                                                     &
                     ,'GROUND_SHV                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ground_ssh                                                     &
                     ,'GROUND_SSH                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ground_temp                                                    &
                     ,'GROUND_TEMP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ground_fliq                                                    &
                     ,'GROUND_FLIQ                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rough                                                          &
                     ,'ROUGH                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_l_max                                                      &
                     ,'PAR_L_MAX                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_l_beam_max                                                 &
                     ,'PAR_L_BEAM_MAX              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_l_diffuse_max                                              &
                     ,'PAR_L_DIFFUSE_MAX           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%avg_daily_temp                                                 &
                     ,'AVG_DAILY_TEMP              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%avg_monthly_gndwater                                           &
                     ,'AVG_MONTHLY_GNDWATER        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%avg_monthly_waterdef                                           &
                     ,'AVG_MONTHLY_WATERDEF        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_loss2atm                                               &
                     ,'WBUDGET_LOSS2ATM            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_denseffect                                             &
                     ,'WBUDGET_DENSEFFECT          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_precipgain                                             &
                     ,'WBUDGET_PRECIPGAIN          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_loss2runoff                                            &
                     ,'WBUDGET_LOSS2RUNOFF         ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_loss2drainage                                          &
                     ,'WBUDGET_LOSS2DRAINAGE       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_initialstorage                                         &
                     ,'WBUDGET_INITIALSTORAGE      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wbudget_residual                                               &
                     ,'WBUDGET_RESIDUAL            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_loss2atm                                               &
                     ,'EBUDGET_LOSS2ATM            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_denseffect                                             &
                     ,'EBUDGET_DENSEFFECT          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_prsseffect                                             &
                     ,'EBUDGET_PRSSEFFECT          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_loss2runoff                                            &
                     ,'EBUDGET_LOSS2RUNOFF         ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_loss2drainage                                          &
                     ,'EBUDGET_LOSS2DRAINAGE       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_netrad                                                 &
                     ,'EBUDGET_NETRAD              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_precipgain                                             &
                     ,'EBUDGET_PRECIPGAIN          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_initialstorage                                         &
                     ,'EBUDGET_INITIALSTORAGE      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ebudget_residual                                               &
                     ,'EBUDGET_RESIDUAL            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_initialstorage                                       &
                     ,'CO2BUDGET_INITIALSTORAGE    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_residual                                             &
                     ,'CO2BUDGET_RESIDUAL          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_loss2atm                                             &
                     ,'CO2BUDGET_LOSS2ATM          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_denseffect                                           &
                     ,'CO2BUDGET_DENSEFFECT        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_gpp                                                  &
                     ,'CO2BUDGET_GPP               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_plresp                                               &
                     ,'CO2BUDGET_PLRESP            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%co2budget_rh                                                   &
                     ,'CO2BUDGET_RH                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%today_A_decomp                                                 &
                     ,'TODAY_A_DECOMP              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%today_Af_decomp                                                &
                     ,'TODAY_AF_DECOMP             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%veg_rough                                                      &
                     ,'VEG_ROUGH                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%veg_height                                                     &
                     ,'VEG_HEIGHT                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%veg_displace                                                   &
                     ,'VEG_DISPLACE                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%fsc_in                                                         &
                     ,'FSC_IN                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ssc_in                                                         &
                     ,'SSC_IN                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ssl_in                                                         &
                     ,'SSL_IN                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%fsn_in                                                         &
                     ,'FSN_IN                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%total_plant_nitrogen_uptake                                    &
                     ,'TOTAL_PLANT_NITROGEN_UPTAKE ',dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(csite%mineralized_N_loss                                             &
                     ,'MINERALIZED_N_LOSS          ',dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(csite%mineralized_N_input                                            &
                     ,'MINERALIZED_N_INPUT         ',dsetrank,iparallel,.false.,foundvar)
   call hdf_getslab_r(csite%rshort_g                                                       &
                     ,'RSHORT_G                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rshort_g_beam                                                  &
                     ,'RSHORT_G_BEAM               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rshort_g_diffuse                                               &
                     ,'RSHORT_G_DIFFUSE            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_g                                                          &
                     ,'PAR_G                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_g_beam                                                     &
                     ,'PAR_G_BEAM                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_g_diffuse                                                  &
                     ,'PAR_G_DIFFUSE               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_b                                                          &
                     ,'PAR_B                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_b_beam                                                     &
                     ,'PAR_B_BEAM                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_b_diffuse                                                  &
                     ,'PAR_B_DIFFUSE               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%nir_b                                                          &
                     ,'NIR_B                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%nir_b_beam                                                     &
                     ,'NIR_B_BEAM                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%nir_b_diffuse                                                  &
                     ,'NIR_B_DIFFUSE               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rlong_g                                                        &
                     ,'RLONG_G                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rlong_s                                                        &
                     ,'RLONG_S                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%albedo                                                         &
                     ,'ALBEDO                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%albedo_par                                                     &
                     ,'ALBEDO_PAR                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%albedo_nir                                                     &
                     ,'ALBEDO_NIR                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rlong_albedo                                                   &
                     ,'RLONG_ALBEDO                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rnet                                                           &
                     ,'RNET                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rlongup                                                        &
                     ,'RLONGUP                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%parup                                                          &
                     ,'PARUP                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%nirup                                                          &
                     ,'NIRUP                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rshortup                                                       &
                     ,'RSHORTUP                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%total_sfcw_depth                                               &
                     ,'TOTAL_SFCW_DEPTH            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%snowfac                                                        &
                     ,'SNOWFAC                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%A_decomp                                                       &
                     ,'A_DECOMP                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%f_decomp                                                       &
                     ,'F_DECOMP                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rh                                                             &
                     ,'RH                          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%cwd_rh                                                         &
                     ,'CWD_RH                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%plant_ag_biomass                                               &
                     ,'PLANT_AG_BIOMASS            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ustar                                                          &
                     ,'USTAR                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%tstar                                                          &
                     ,'TSTAR                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%qstar                                                          &
                     ,'QSTAR                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%cstar                                                          &
                     ,'CSTAR                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%zeta                                                           &
                     ,'ZETA                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%ribulk                                                         &
                     ,'RIBULK                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%upwp                                                           &
                     ,'UPWP                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%qpwp                                                           &
                     ,'QPWP                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%cpwp                                                           &
                     ,'CPWP                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%tpwp                                                           &
                     ,'TPWP                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%wpwp                                                           &
                     ,'WPWP                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%htry                                                           &
                     ,'HTRY                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%hprev                                                          &
                     ,'HPREV                       ',dsetrank,iparallel,.true. ,foundvar)
   !----- Daily means. --------------------------------------------------------------------!
   if (writing_long) then
      call hdf_getslab_r(csite%dmean_A_decomp                                              &
                        ,'DMEAN_A_DECOMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_Af_decomp                                             &
                        ,'DMEAN_AF_DECOMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_co2_residual                                          &
                        ,'DMEAN_CO2_RESIDUAL_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_energy_residual                                       &
                        ,'DMEAN_ENERGY_RESIDUAL_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_water_residual                                        &
                        ,'DMEAN_WATER_RESIDUAL_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rh                                                    &
                        ,'DMEAN_RH_PA               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_cwd_rh                                                &
                        ,'DMEAN_CWD_RH_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_nep                                                   &
                        ,'DMEAN_NEP_PA              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rk4step                                               &
                        ,'DMEAN_RK4STEP_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_available_water                                       &
                        ,'DMEAN_AVAILABLE_WATER_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_theiv                                             &
                        ,'DMEAN_CAN_THEIV_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_theta                                             &
                        ,'DMEAN_CAN_THETA_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_vpdef                                             &
                        ,'DMEAN_CAN_VPDEF_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_temp                                              &
                        ,'DMEAN_CAN_TEMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_shv                                               &
                        ,'DMEAN_CAN_SHV_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_co2                                               &
                        ,'DMEAN_CAN_CO2_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_rhos                                              &
                        ,'DMEAN_CAN_RHOS_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_prss                                              &
                        ,'DMEAN_CAN_PRSS_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_gnd_temp                                              &
                        ,'DMEAN_GND_TEMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_gnd_shv                                               &
                        ,'DMEAN_GND_SHV_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_can_ggnd                                              &
                        ,'DMEAN_CAN_GGND_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sfcw_depth                                            &
                        ,'DMEAN_SFCW_DEPTH_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sfcw_energy                                           &
                        ,'DMEAN_SFCW_ENERGY_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sfcw_mass                                             &
                        ,'DMEAN_SFCW_MASS_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sfcw_temp                                             &
                        ,'DMEAN_SFCW_TEMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sfcw_fliq                                             &
                        ,'DMEAN_SFCW_FLIQ_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rshort_gnd                                            &
                        ,'DMEAN_RSHORT_GND_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_par_gnd                                               &
                        ,'DMEAN_PAR_GND_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rlong_gnd                                             &
                        ,'DMEAN_RLONG_GND_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rlongup                                               &
                        ,'DMEAN_RLONGUP_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_parup                                                 &
                        ,'DMEAN_PARUP_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_nirup                                                 &
                        ,'DMEAN_NIRUP_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rshortup                                              &
                        ,'DMEAN_RSHORTUP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rnet                                                  &
                        ,'DMEAN_RNET_PA             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_albedo                                                &
                        ,'DMEAN_ALBEDO_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_albedo_par                                            &
                        ,'DMEAN_ALBEDO_PAR_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_albedo_nir                                            &
                        ,'DMEAN_ALBEDO_NIR_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_rlong_albedo                                          &
                        ,'DMEAN_RLONG_ALBEDO_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_ustar                                                 &
                        ,'DMEAN_USTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_tstar                                                 &
                        ,'DMEAN_TSTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_qstar                                                 &
                        ,'DMEAN_QSTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_cstar                                                 &
                        ,'DMEAN_CSTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_carbon_ac                                             &
                        ,'DMEAN_CARBON_AC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_carbon_st                                             &
                        ,'DMEAN_CARBON_ST_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_vapor_gc                                              &
                        ,'DMEAN_VAPOR_GC_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_vapor_ac                                              &
                        ,'DMEAN_VAPOR_AC_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_throughfall                                           &
                        ,'DMEAN_THROUGHFALL_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_runoff                                                &
                        ,'DMEAN_RUNOFF_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_drainage                                              &
                        ,'DMEAN_DRAINAGE_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sensible_gc                                           &
                        ,'DMEAN_SENSIBLE_GC_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sensible_ac                                           &
                        ,'DMEAN_SENSIBLE_AC_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_qthroughfall                                          &
                        ,'DMEAN_QTHROUGHFALL_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_qrunoff                                               &
                        ,'DMEAN_QRUNOFF_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_qdrainage                                             &
                        ,'DMEAN_QDRAINAGE_PA        ',dsetrank,iparallel,.false.,foundvar)
   end if
   !----- Monthly means. ------------------------------------------------------------------!
   if (writing_eorq) then
      call hdf_getslab_r(csite%mmean_fast_soil_c                                           &
                        ,'MMEAN_FAST_SOIL_C_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_slow_soil_c                                           &
                        ,'MMEAN_SLOW_SOIL_C_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_struct_soil_c                                         &
                        ,'MMEAN_STRUCT_SOIL_C_PA    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_struct_soil_l                                         &
                        ,'MMEAN_STRUCT_SOIL_L_PA    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_fast_soil_n                                           &
                        ,'MMEAN_FAST_SOIL_N_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_mineral_soil_n                                        &
                        ,'MMEAN_MINERAL_SOIL_N_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_co2_residual                                          &
                        ,'MMEAN_CO2_RESIDUAL_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_energy_residual                                       &
                        ,'MMEAN_ENERGY_RESIDUAL_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_water_residual                                        &
                        ,'MMEAN_WATER_RESIDUAL_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rh                                                    &
                        ,'MMEAN_RH_PA               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_cwd_rh                                                &
                        ,'MMEAN_CWD_RH_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_nep                                                   &
                        ,'MMEAN_NEP_PA              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_A_decomp                                              &
                        ,'MMEAN_A_DECOMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_Af_decomp                                             &
                        ,'MMEAN_AF_DECOMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rk4step                                               &
                        ,'MMEAN_RK4STEP_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_available_water                                       &
                        ,'MMEAN_AVAILABLE_WATER_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_theiv                                             &
                        ,'MMEAN_CAN_THEIV_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_theta                                             &
                        ,'MMEAN_CAN_THETA_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_vpdef                                             &
                        ,'MMEAN_CAN_VPDEF_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_temp                                              &
                        ,'MMEAN_CAN_TEMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_shv                                               &
                        ,'MMEAN_CAN_SHV_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_co2                                               &
                        ,'MMEAN_CAN_CO2_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_rhos                                              &
                        ,'MMEAN_CAN_RHOS_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_prss                                              &
                        ,'MMEAN_CAN_PRSS_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_gnd_temp                                              &
                        ,'MMEAN_GND_TEMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_gnd_shv                                               &
                        ,'MMEAN_GND_SHV_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_can_ggnd                                              &
                        ,'MMEAN_CAN_GGND_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sfcw_depth                                            &
                        ,'MMEAN_SFCW_DEPTH_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sfcw_energy                                           &
                        ,'MMEAN_SFCW_ENERGY_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sfcw_mass                                             &
                        ,'MMEAN_SFCW_MASS_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sfcw_temp                                             &
                        ,'MMEAN_SFCW_TEMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sfcw_fliq                                             &
                        ,'MMEAN_SFCW_FLIQ_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rshort_gnd                                            &
                        ,'MMEAN_RSHORT_GND_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_par_gnd                                               &
                        ,'MMEAN_PAR_GND_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rlong_gnd                                             &
                        ,'MMEAN_RLONG_GND_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rlongup                                               &
                        ,'MMEAN_RLONGUP_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_parup                                                 &
                        ,'MMEAN_PARUP_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_nirup                                                 &
                        ,'MMEAN_NIRUP_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rshortup                                              &
                        ,'MMEAN_RSHORTUP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rnet                                                  &
                        ,'MMEAN_RNET_PA             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_albedo                                                &
                        ,'MMEAN_ALBEDO_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_albedo_par                                            &
                        ,'MMEAN_ALBEDO_PAR_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_albedo_nir                                            &
                        ,'MMEAN_ALBEDO_NIR_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_rlong_albedo                                          &
                        ,'MMEAN_RLONG_ALBEDO_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_ustar                                                 &
                        ,'MMEAN_USTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_tstar                                                 &
                        ,'MMEAN_TSTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_qstar                                                 &
                        ,'MMEAN_QSTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_cstar                                                 &
                        ,'MMEAN_CSTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_carbon_ac                                             &
                        ,'MMEAN_CARBON_AC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_carbon_st                                             &
                        ,'MMEAN_CARBON_ST_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_vapor_gc                                              &
                        ,'MMEAN_VAPOR_GC_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_vapor_ac                                              &
                        ,'MMEAN_VAPOR_AC_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_throughfall                                           &
                        ,'MMEAN_THROUGHFALL_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_runoff                                                &
                        ,'MMEAN_RUNOFF_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_drainage                                              &
                        ,'MMEAN_DRAINAGE_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sensible_gc                                           &
                        ,'MMEAN_SENSIBLE_GC_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sensible_ac                                           &
                        ,'MMEAN_SENSIBLE_AC_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_qthroughfall                                          &
                        ,'MMEAN_QTHROUGHFALL_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_qrunoff                                               &
                        ,'MMEAN_QRUNOFF_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_qdrainage                                             &
                        ,'MMEAN_QDRAINAGE_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_A_decomp                                              &
                        ,'MMEAN_A_DECOMP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_Af_decomp                                             &
                        ,'MMEAN_AF_DECOMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_co2_residual                                          &
                        ,'MMEAN_CO2_RESIDUAL_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_energy_residual                                       &
                        ,'MMEAN_ENERGY_RESIDUAL_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_water_residual                                        &
                        ,'MMEAN_WATER_RESIDUAL_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_rh                                                    &
                        ,'MMSQU_RH_PA               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_cwd_rh                                                &
                        ,'MMSQU_CWD_RH_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_nep                                                   &
                        ,'MMSQU_NEP_PA              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_rlongup                                               &
                        ,'MMSQU_RLONGUP_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_parup                                                 &
                        ,'MMSQU_PARUP_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_nirup                                                 &
                        ,'MMSQU_NIRUP_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_rshortup                                              &
                        ,'MMSQU_RSHORTUP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_rnet                                                  &
                        ,'MMSQU_RNET_PA             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_albedo                                                &
                        ,'MMSQU_ALBEDO_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_ustar                                                 &
                        ,'MMSQU_USTAR_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_carbon_ac                                             &
                        ,'MMSQU_CARBON_AC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_carbon_st                                             &
                        ,'MMSQU_CARBON_ST_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_vapor_gc                                              &
                        ,'MMSQU_VAPOR_GC_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_vapor_ac                                              &
                        ,'MMSQU_VAPOR_AC_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_sensible_gc                                           &
                        ,'MMSQU_SENSIBLE_GC_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmsqu_sensible_ac                                           &
                        ,'MMSQU_SENSIBLE_AC_PA      ',dsetrank,iparallel,.false.,foundvar)
   end if                                                                           
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (ndcycle,npatches).                                   !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(ndcycle,8)
   chnkdims(1) = int(ndcycle,8)
   memdims (1) = int(ndcycle,8)
   memsize (1) = int(ndcycle,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches ,8)
   chnkoffs(2) = int(sipa_index - 1 ,8)
   memdims (2) = int(csite%npatches ,8)
   memsize (2) = int(csite%npatches ,8)
   memoffs (2) = 0_8

   !----- Mean diel variables. ------------------------------------------------------------!
   if (writing_dcyc) then
      call hdf_getslab_r(csite%qmean_rh                                                    &
                        ,'QMEAN_RH_PA              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_cwd_rh                                                &
                        ,'QMEAN_CWD_RH_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_nep                                                   &
                        ,'QMEAN_NEP_PA             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rk4step                                               &
                        ,'QMEAN_RK4STEP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_available_water                                       &
                        ,'QMEAN_AVAILABLE_WATER_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_theiv                                             &
                        ,'QMEAN_CAN_THEIV_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_theta                                             &
                        ,'QMEAN_CAN_THETA_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_vpdef                                             &
                        ,'QMEAN_CAN_VPDEF_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_temp                                              &
                        ,'QMEAN_CAN_TEMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_shv                                               &
                        ,'QMEAN_CAN_SHV_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_co2                                               &
                        ,'QMEAN_CAN_CO2_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_rhos                                              &
                        ,'QMEAN_CAN_RHOS_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_prss                                              &
                        ,'QMEAN_CAN_PRSS_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_gnd_temp                                              &
                        ,'QMEAN_GND_TEMP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_gnd_shv                                               &
                        ,'QMEAN_GND_SHV_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_can_ggnd                                              &
                        ,'QMEAN_CAN_GGND_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sfcw_depth                                            &
                        ,'QMEAN_SFCW_DEPTH_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sfcw_energy                                           &
                        ,'QMEAN_SFCW_ENERGY_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sfcw_mass                                             &
                        ,'QMEAN_SFCW_MASS_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sfcw_temp                                             &
                        ,'QMEAN_SFCW_TEMP_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sfcw_fliq                                             &
                        ,'QMEAN_SFCW_FLIQ_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rshort_gnd                                            &
                        ,'QMEAN_RSHORT_GND_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_par_gnd                                               &
                        ,'QMEAN_PAR_GND_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rlong_gnd                                             &
                        ,'QMEAN_RLONG_GND_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rlongup                                               &
                        ,'QMEAN_RLONGUP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_parup                                                 &
                        ,'QMEAN_PARUP_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_nirup                                                 &
                        ,'QMEAN_NIRUP_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rshortup                                              &
                        ,'QMEAN_RSHORTUP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rnet                                                  &
                        ,'QMEAN_RNET_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_albedo                                                &
                        ,'QMEAN_ALBEDO_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_albedo_par                                            &
                        ,'QMEAN_ALBEDO_PAR_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_albedo_nir                                            &
                        ,'QMEAN_ALBEDO_NIR_PA      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_rlong_albedo                                          &
                        ,'QMEAN_RLONG_ALBEDO_PA    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_ustar                                                 &
                        ,'QMEAN_USTAR_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_tstar                                                 &
                        ,'QMEAN_TSTAR_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_qstar                                                 &
                        ,'QMEAN_QSTAR_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_cstar                                                 &
                        ,'QMEAN_CSTAR_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_carbon_ac                                             &
                        ,'QMEAN_CARBON_AC_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_carbon_st                                             &
                        ,'QMEAN_CARBON_ST_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_vapor_gc                                              &
                        ,'QMEAN_VAPOR_GC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_vapor_ac                                              &
                        ,'QMEAN_VAPOR_AC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_throughfall                                           &
                        ,'QMEAN_THROUGHFALL_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_runoff                                                &
                        ,'QMEAN_RUNOFF_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_drainage                                              &
                        ,'QMEAN_DRAINAGE_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sensible_gc                                           &
                        ,'QMEAN_SENSIBLE_GC_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sensible_ac                                           &
                        ,'QMEAN_SENSIBLE_AC_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_qthroughfall                                          &
                        ,'QMEAN_QTHROUGHFALL_PA    ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_qrunoff                                               &
                        ,'QMEAN_QRUNOFF_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_qdrainage                                             &
                        ,'QMEAN_QDRAINAGE_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_rh                                                    &
                        ,'QMSQU_RH_PA              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_cwd_rh                                                &
                        ,'QMSQU_CWD_RH_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_nep                                                   &
                        ,'QMSQU_NEP_PA             ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_rlongup                                               &
                        ,'QMSQU_RLONGUP_PA         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_parup                                                 &
                        ,'QMSQU_PARUP_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_nirup                                                 &
                        ,'QMSQU_NIRUP_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_rshortup                                              &
                        ,'QMSQU_RSHORTUP_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_rnet                                                  &
                        ,'QMSQU_RNET_PA            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_albedo                                                &
                        ,'QMSQU_ALBEDO_PA          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_ustar                                                 &
                        ,'QMSQU_USTAR_PA           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_carbon_ac                                             &
                        ,'QMSQU_CARBON_AC_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_carbon_st                                             &
                        ,'QMSQU_CARBON_ST_PA       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_vapor_gc                                              &
                        ,'QMSQU_VAPOR_GC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_vapor_ac                                              &
                        ,'QMSQU_VAPOR_AC_PA        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_sensible_gc                                           &
                        ,'QMSQU_SENSIBLE_GC_PA     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmsqu_sensible_ac                                           &
                        ,'QMSQU_SENSIBLE_AC_PA     ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (nzs,npatches).                                       !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(nzs,8)
   chnkdims(1) = int(nzs,8)
   memdims (1) = int(nzs,8)
   memsize (1) = int(nzs,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches ,8)
   chnkoffs(2) = int(sipa_index - 1 ,8)
   memdims (2) = int(csite%npatches ,8)
   memsize (2) = int(csite%npatches ,8)
   memoffs (2) = 0_8
   call hdf_getslab_r(csite%sfcwater_mass                                                  &
                     ,'SFCWATER_MASS        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%sfcwater_energy                                                &
                     ,'SFCWATER_ENERGY      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%sfcwater_depth                                                 &
                     ,'SFCWATER_DEPTH       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%sfcwater_tempk                                                 &
                     ,'SFCWATER_TEMPK       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%sfcwater_fracliq                                               &
                     ,'SFCWATER_FRACLIQ     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rshort_s                                                       &
                     ,'RSHORT_S             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rshort_s_beam                                                  &
                     ,'RSHORT_S_BEAM        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rshort_s_diffuse                                               &
                     ,'RSHORT_S_DIFFUSE     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_s                                                          &
                     ,'PAR_S                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_s_beam                                                     &
                     ,'PAR_S_BEAM           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%par_s_diffuse                                                  &
                     ,'PAR_S_DIFFUSE        ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (nzg,npatches).                                       !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims (1) = int(nzg,8)
   memsize (1) = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches ,8)
   chnkoffs(2) = int(sipa_index - 1 ,8)
   memdims (2) = int(csite%npatches ,8)
   memsize (2) = int(csite%npatches ,8)
   memoffs (2) = 0_8
   call hdf_getslab_r(csite%soil_energy                                                    &
                     ,'SOIL_ENERGY_PA       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%soil_mstpot                                                    &
                     ,'SOIL_MSTPOT_PA       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%soil_water                                                     &
                     ,'SOIL_WATER_PA        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%soil_tempk                                                     &
                     ,'SOIL_TEMPK_PA        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%soil_fracliq                                                   &
                     ,'SOIL_FRACLIQ_PA      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%rootdense                                                      &
                     ,'PATCH_ROOT_DENSITY   ',dsetrank,iparallel,.true. ,foundvar)
   if (writing_long) then
      call hdf_getslab_r(csite%dmean_soil_energy                                           &
                        ,'DMEAN_SOIL_ENERGY_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_soil_mstpot                                           &
                        ,'DMEAN_SOIL_MSTPOT_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_soil_water                                            &
                        ,'DMEAN_SOIL_WATER_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_soil_temp                                             &
                        ,'DMEAN_SOIL_TEMP_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_soil_fliq                                             &
                        ,'DMEAN_SOIL_FLIQ_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_smoist_gg                                             &
                        ,'DMEAN_SMOIST_GG_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_transloss                                             &
                        ,'DMEAN_TRANSLOSS_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%dmean_sensible_gg                                           &
                        ,'DMEAN_SENSIBLE_GG_PA ',dsetrank,iparallel,.false.,foundvar)
   end if
   if (writing_eorq) then
      call hdf_getslab_r(csite%mmean_soil_energy                                           &
                        ,'MMEAN_SOIL_ENERGY_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_soil_mstpot                                           &
                        ,'MMEAN_SOIL_MSTPOT_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_soil_water                                            &
                        ,'MMEAN_SOIL_WATER_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_soil_temp                                             &
                        ,'MMEAN_SOIL_TEMP_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_soil_fliq                                             &
                        ,'MMEAN_SOIL_FLIQ_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_smoist_gg                                             &
                        ,'MMEAN_SMOIST_GG_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_transloss                                             &
                        ,'MMEAN_TRANSLOSS_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%mmean_sensible_gg                                           &
                        ,'MMEAN_SENSIBLE_GG_PA ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (n_pft,npatches).                                     !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims (1) = int(n_pft,8)
   memsize (1) = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches,8)
   chnkoffs(2) = int(sipa_index - 1,8)
   memdims (2) = int(csite%npatches,8)
   memsize (2) = int(csite%npatches,8)
   memoffs (2) = 0_8

   call hdf_getslab_r(csite%repro                                                          &
                     ,'REPRO_PA             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%A_o_max                                                        &
                     ,'A_O_MAX              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(csite%A_c_max                                                        &
                     ,'A_C_MAX              ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (n_pft,ff_nhgt,npatches).                             !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 3
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims (1) = int(n_pft,8)
   memsize (1) = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(ff_nhgt,8)
   chnkdims(2) = int(ff_nhgt,8)
   memdims (2) = int(ff_nhgt,8)
   memsize (2) = int(ff_nhgt,8)
   chnkoffs(2) = 0_8
   memoffs (2) = 0_8
   globdims(3) = int(npatches_global,8)
   chnkdims(3) = int(csite%npatches ,8)
   chnkoffs(3) = int(sipa_index - 1 ,8)
   memdims (3) = int(csite%npatches ,8)
   memsize (3) = int(csite%npatches ,8)
   memoffs (3) = 0_8
   call hdf_getslab_r(csite%cumlai_profile                                                 &
                     ,'CUMLAI_PROFILE       ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (nzg,ndcycle,npatches).                               !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 3
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims (1) = int(nzg,8)
   memsize (1) = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs (1) = 0_8
   globdims(2) = int(ndcycle,8)
   chnkdims(2) = int(ndcycle,8)
   memdims (2) = int(ndcycle,8)
   memsize (2) = int(ndcycle,8)
   chnkoffs(2) = 0_8
   memoffs (2) = 0_8
   globdims(3) = int(npatches_global,8)
   chnkdims(3) = int(csite%npatches ,8)
   chnkoffs(3) = int(sipa_index - 1 ,8)
   memdims (3) = int(csite%npatches ,8)
   memsize (3) = int(csite%npatches ,8)
   memoffs (3) = 0_8
   if (writing_dcyc) then
      call hdf_getslab_r(csite%qmean_soil_energy                                           &
                        ,'QMEAN_SOIL_ENERGY_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_soil_mstpot                                           &
                        ,'QMEAN_SOIL_MSTPOT_PA ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_soil_water                                            &
                        ,'QMEAN_SOIL_WATER_PA  ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_soil_temp                                             &
                        ,'QMEAN_SOIL_TEMP_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_soil_fliq                                             &
                        ,'QMEAN_SOIL_FLIQ_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_smoist_gg                                             &
                        ,'QMEAN_SMOIST_GG_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_transloss                                             &
                        ,'QMEAN_TRANSLOSS_PA   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(csite%qmean_sensible_gg                                           &
                        ,'QMEAN_SENSIBLE_GG_PA ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Once upon a time, soil water was double precision, but after the Runge-Kutta      !
   ! integrator was entirely converted to double precision, we reverted back to single and !
   ! converted everything just for the Runge-Kutta.  Of course, nothing in life is certain !
   ! but death and taxes, so we kept the recipe on how to read a double precision if we    !
   ! ever change our minds again.                                                          !
   !                                                                                       !
   !   call h5dopen_f(file_id,'SOIL_WATER_PA ', dset_id, hdferr)                           !
   !   if (hdferr /= 0 ) then                                                              !
   !      call fatal_error('Dataset did not have soil water?' &                            !
   !           ,'fill_history_site','ed_init_full_history.F90')                            !
   !   endif                                                                               !
   !                                                                                       !
   !------ These lines are useful for determining data size of any given object in a set.--!
   !  call h5dget_type_f(dset_id,datatype_id,hdferr)                                       !
   !  call h5tget_size_f(datatype_id,setsize,hdferr)                                       !
   !  call h5dclose_f(dset_id  , hdferr)                                                   !
   !                                                                                       !
   !------ Check precision. ---------------------------------------------------------------!
   !   if (setsize==4_8) then  ! Single precision                                          !
   !      call hdf_getslab_r(csite%soil_water                                            & !
   !                        ,'SOIL_WATER_PA ',dsetrank,iparallel,.true.,foundvar)          !
   !   else if (setsize==8_8) then ! Double precision                                      !
   !      allocate(buff(nzg,csite%npatches))                                               !
   !      write (unit=*,fmt='(a)') '----------------------------------------------------'  !
   !      write (unit=*,fmt='(a)') '  Load 8-byte precision, and convert to 4-byte'        !
   !      write (unit=*,fmt='(a)') '----------------------------------------------------'  !
   !      call hdf_getslab_d(buff,'SOIL_WATER_PA ',dsetrank,iparallel,.true.,foundvar)     !
   !      csite%soil_water(1:nzg,1:csite%npatches) = sngl(buff(1:nzg,1:csite%npatches))    !
   !      deallocate(buff)                                                                 !
   !  else                                                                                 !
   !     call fatal_error('Soil water dataset is not real nor double?'                   & !
   !                     ,'fill_history_site','ed_init_full_history.F90')                  !
   !  end if                                                                               !
   !---------------------------------------------------------------------------------------!

   return
end subroutine fill_history_site
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine loads all state variables and partial integrations at the cohort    !
! level.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine fill_history_patch(cpatch,paco_index,ncohorts_global)
   use ed_state_vars      , only : patchtype     ! ! structure
   use grid_coms          , only : nzg           & ! intent(in)
                                 , nzs           ! ! intent(in)
   use ed_max_dims        , only : n_pft         & ! intent(in)
                                 , n_mort        & ! intent(in)
                                 , max_site      & ! intent(in)
                                 , n_dist_types  & ! intent(in)
                                 , n_radprof     ! ! intent(in)
   use hdf5
   use hdf5_coms          , only : file_id       & ! intent(inout)
                                 , dset_id       & ! intent(inout)
                                 , dspace_id     & ! intent(inout)
                                 , plist_id      & ! intent(inout)
                                 , globdims      & ! intent(inout)
                                 , chnkdims      & ! intent(inout)
                                 , chnkoffs      & ! intent(inout)
                                 , cnt           & ! intent(inout)
                                 , stride        & ! intent(inout)
                                 , memdims       & ! intent(inout)
                                 , memoffs       & ! intent(inout)
                                 , memsize       & ! intent(inout)
                                 , datatype_id   & ! intent(inout)
                                 , setsize       ! ! intent(inout)
   use ed_misc_coms       , only : ndcycle       & ! intent(in)
                                 , writing_long  & ! intent(in)
                                 , writing_eorq  & ! intent(in)
                                 , writing_dcyc  ! ! intent(in)
   use fusion_fission_coms, only : ff_nhgt       ! ! intent(in)
   use allometry          , only : dbh2ca        ! ! function
   implicit none
   !----- Interfaces. ---------------------------------------------------------------------!
#if USE_INTERF
   interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_r
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_d
      !------------------------------------------------------------------------------------!
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
         use hdf5_coms, only : memsize ! ! intent(in)
         !----- Arguments. ----------------------------------------------------------------!
         integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))          &
                                                          , intent(inout) :: buff
         character(len=*)                                 , intent(in)    :: varn
         integer                                          , intent(in)    :: dsetrank
         integer                                          , intent(in)    :: iparallel
         logical                                          , intent(in)    :: required
         logical                                          , intent(out)   :: foundvar
         !---------------------------------------------------------------------------------!
      end subroutine hdf_getslab_i
      !------------------------------------------------------------------------------------!
   end interface
#endif
   !---------------------------------------------------------------------------------------!

   !----- Arguments. ----------------------------------------------------------------------!
   type(patchtype)                  , target       :: cpatch
   integer                          , intent(in)   :: paco_index
   integer                          , intent(in)   :: ncohorts_global
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: iparallel
   integer                                         :: dsetrank
   logical                                         :: foundvar
   integer                                         :: ico
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     If this patch is empty, we don't need to do anything (actually the less we do the !
   ! better, we don't want to add segmentation violations in our beloved model!).          !
   !---------------------------------------------------------------------------------------!
   if (cpatch%ncohorts == 0) return
   !---------------------------------------------------------------------------------------!

   !----- Turn off parallel for this sub-routine. -----------------------------------------!
   iparallel = 0
   !---------------------------------------------------------------------------------------!
 

   !---------------------------------------------------------------------------------------!
   !     Reset the dimension arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   globdims(:) = 0_8
   chnkdims(:) = 0_8
   chnkoffs(:) = 0_8
   memoffs (:) = 0_8
   memdims (:) = 0_8
   memsize (:) = 1_8
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For the remainder of this sub-routine, we must load the polygon variables         !
   ! according to the dimensionality.  If you are adding a new variable, make sure that it !
   ! is placed together with variables of the same dimensions.                             !
   !                                                                                       !
   ! DSETRANK -- the rank of the variable.                                                 !
   ! GLOBDIMS -- the size of each dimension                                                !
   ! CHNKDIMS -- the size of the chunk to be read this time                                !
   ! CHNKOFFS -- the offset (number of indices to be skipped)                              !
   ! MEMDIMS  -- the dimensions of the buffer to be filled.                                !
   ! MEMOFFS  -- the offset of the memory                                                  !
   ! MEMSIZE  -- the size of variable that will house the information read.                !
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables - Integers.                                           !
   !---------------------------------------------------------------------------------------!
   dsetrank = 1
   globdims(1) = int(ncohorts_global,8)
   chnkdims(1) = int(cpatch%ncohorts,8)
   chnkoffs(1) = int(paco_index - 1 ,8)
   memdims (1) = int(cpatch%ncohorts,8)
   memsize (1) = int(cpatch%ncohorts,8)
   memoffs (1) = 0_8
   call hdf_getslab_i(cpatch%pft                                                           &
                     ,'PFT                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpatch%phenology_status                                              &
                     ,'PHENOLOGY_STATUS          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpatch%recruit_dbh                                                   &
                     ,'RECRUIT_DBH               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpatch%census_status                                                 &
                     ,'CENSUS_STATUS             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpatch%krdepth                                                       &
                     ,'KRDEPTH                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpatch%first_census                                                  &
                     ,'FIRST_CENSUS              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_i(cpatch%new_recruit_flag                                              &
                     ,'NEW_RECRUIT_FLAG          ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!


   !----- Swap phenology_status from 2 to -2. ---------------------------------------------!
   do ico = 1,cpatch%ncohorts
      if (cpatch%phenology_status(ico) == 2) cpatch%phenology_status(ico) = -2
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Single dimension variables - Real.                                               !
   !---------------------------------------------------------------------------------------!
   dsetrank = 1
   globdims(1) = int(ncohorts_global,8)
   chnkdims(1) = int(cpatch%ncohorts,8)
   chnkoffs(1) = int(paco_index - 1 ,8)
   memdims (1) = int(cpatch%ncohorts,8)
   memsize (1) = int(cpatch%ncohorts,8)
   memoffs (1) = 0_8
   call hdf_getslab_r(cpatch%nplant                                                        &
                     ,'NPLANT                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%hite                                                          &
                     ,'HITE                      ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%agb                                                           &
                     ,'AGB_CO                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%basarea                                                       &
                     ,'BA_CO                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%dagb_dt                                                       &
                     ,'DAGB_DT                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%dlnagb_dt                                                     &
                     ,'DLNAGB_DT                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%dba_dt                                                        &
                     ,'DBA_DT                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%dlnba_dt                                                      &
                     ,'DLNBA_DT                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%ddbh_dt                                                       &
                     ,'DDBH_DT                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%dlndbh_dt                                                     &
                     ,'DLNDBH_DT                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%dbh                                                           &
                     ,'DBH                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%bdead                                                         &
                     ,'BDEAD                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%bleaf                                                         &
                     ,'BLEAF                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%balive                                                        &
                     ,'BALIVE                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%broot                                                         &
                     ,'BROOT                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%bsapwooda                                                     &
                     ,'BSAPWOODA                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%bsapwoodb                                                     &
                     ,'BSAPWOODB                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%bstorage                                                      &
                     ,'BSTORAGE                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%bseeds                                                        &
                     ,'BSEEDS_CO                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lai                                                           &
                     ,'LAI_CO                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wai                                                           &
                     ,'WAI_CO                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%crown_area                                                    &
                     ,'CROWN_AREA_CO             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%cbr_bar                                                       &
                     ,'CBR_BAR                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_energy                                                   &
                     ,'LEAF_ENERGY               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_temp                                                     &
                     ,'LEAF_TEMP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_vpdef                                                    &
                     ,'LEAF_VPDEF                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_temp_pv                                                  &
                     ,'LEAF_TEMP_PV              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_hcap                                                     &
                     ,'LEAF_HCAP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_fliq                                                     &
                     ,'LEAF_FLIQ                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_water                                                    &
                     ,'LEAF_WATER                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_energy                                                   &
                     ,'WOOD_ENERGY               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_temp                                                     &
                     ,'WOOD_TEMP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_temp_pv                                                  &
                     ,'WOOD_TEMP_PV              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_hcap                                                     &
                     ,'WOOD_HCAP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_fliq                                                     &
                     ,'WOOD_FLIQ                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_water                                                    &
                     ,'WOOD_WATER                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%veg_wind                                                      &
                     ,'VEG_WIND                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lsfc_shv_open                                                 &
                     ,'LSFC_SHV_OPEN             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lsfc_shv_closed                                               &
                     ,'LSFC_SHV_CLOSED           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lsfc_co2_open                                                 &
                     ,'LSFC_CO2_OPEN             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lsfc_co2_closed                                               &
                     ,'LSFC_CO2_CLOSED           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lint_shv                                                      &
                     ,'LINT_SHV                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lint_co2_open                                                 &
                     ,'LINT_CO2_OPEN             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%lint_co2_closed                                               &
                     ,'LINT_CO2_CLOSED           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_leaf_resp                                               &
                     ,'TODAY_LEAF_RESP           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_root_resp                                               &
                     ,'TODAY_ROOT_RESP           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_gpp                                                     &
                     ,'TODAY_GPP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_gpp_pot                                                 &
                     ,'TODAY_GPP_POT             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_gpp_lightmax                                            &
                     ,'TODAY_GPP_LIGHTMAX        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_gpp_moistmax                                            &
                     ,'TODAY_GPP_MOISTMAX        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_gpp_mlmax                                               &
                     ,'TODAY_GPP_MLMAX           ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppleaf                                                 &
                     ,'TODAY_NPPLEAF             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppfroot                                                &
                     ,'TODAY_NPPFROOT            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppsapwood                                              &
                     ,'TODAY_NPPSAPWOOD          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppcroot                                                &
                     ,'TODAY_NPPCROOT            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppseeds                                                &
                     ,'TODAY_NPPSEEDS            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppwood                                                 &
                     ,'TODAY_NPPWOOD             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%today_nppdaily                                                &
                     ,'TODAY_NPPDAILY            ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_growth_resp                                              &
                     ,'LEAF_GROWTH_RESP   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%root_growth_resp                                              &
                     ,'ROOT_GROWTH_RESP   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%sapa_growth_resp                                              &
                     ,'SAPA_GROWTH_RESP   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%sapb_growth_resp                                              &
                     ,'SAPB_GROWTH_RESP   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_storage_resp                                             &
                     ,'LEAF_STORAGE_RESP       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%root_storage_resp                                             &
                     ,'ROOT_STORAGE_RESP       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%sapa_storage_resp                                             &
                     ,'SAPA_STORAGE_RESP       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%sapb_storage_resp                                             &
                     ,'SAPB_STORAGE_RESP       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%monthly_dndt                                                  &
                     ,'MONTHLY_DNDT              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%monthly_dlnndt                                                &
                     ,'MONTHLY_DLNNDT            ',dsetrank,iparallel,.true. ,foundvar)

   call hdf_getslab_r(cpatch%par_level_beam                                                &
                     ,'PAR_LEVEL_BEAM              ',dsetrank,iparallel,.false. ,foundvar)
   call hdf_getslab_r(cpatch%par_level_diffd                                               &
                     ,'PAR_LEVEL_DIFFD             ',dsetrank,iparallel,.false. ,foundvar)
   call hdf_getslab_r(cpatch%par_level_diffu                                               &
                     ,'PAR_LEVEL_DIFFU             ',dsetrank,iparallel,.false. ,foundvar)

   call hdf_getslab_r(cpatch%light_level                                                   &
                     ,'LIGHT_LEVEL               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%light_level_beam                                              &
                     ,'LIGHT_LEVEL_BEAM          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%light_level_diff                                              &
                     ,'LIGHT_LEVEL_DIFF          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%par_l                                                         &
                     ,'PAR_L                     ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%par_l_beam                                                    &
                     ,'PAR_L_BEAM                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%par_l_diffuse                                                 &
                     ,'PAR_L_DIFFUSE             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rshort_l                                                      &
                     ,'RSHORT_L                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rshort_l_beam                                                 &
                     ,'RSHORT_L_BEAM             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rshort_l_diffuse                                              &
                     ,'RSHORT_L_DIFFUSE          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rlong_l                                                       &
                     ,'RLONG_L                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rshort_w                                                      &
                     ,'RSHORT_W                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rshort_w_beam                                                 &
                     ,'RSHORT_W_BEAM             ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rshort_w_diffuse                                              &
                     ,'RSHORT_W_DIFFUSE          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%rlong_w                                                       &
                     ,'RLONG_W                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_gbh                                                      &
                     ,'LEAF_GBH                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_gbw                                                      &
                     ,'LEAF_GBW                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_gbh                                                      &
                     ,'WOOD_GBH                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%wood_gbw                                                      &
                     ,'WOOD_GBW                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%A_open                                                        &
                     ,'A_OPEN                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%A_closed                                                      &
                     ,'A_CLOSED                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%A_light                                                       &
                     ,'A_LIGHT                   ',dsetrank,iparallel,.false. ,foundvar)
   call hdf_getslab_r(cpatch%A_rubp                                                        &
                     ,'A_RUBP                    ',dsetrank,iparallel,.false. ,foundvar)
   call hdf_getslab_r(cpatch%A_co2                                                         &
                     ,'A_CO2                     ',dsetrank,iparallel,.false. ,foundvar)
   call hdf_getslab_r(cpatch%psi_open                                                      &
                     ,'PSI_OPEN                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%psi_closed                                                    &
                     ,'PSI_CLOSED                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%gsw_open                                                      &
                     ,'GSW_OPEN                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%gsw_closed                                                    &
                     ,'GSW_CLOSED                ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_gsw                                                      &
                     ,'LEAF_GSW                  ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%fsw                                                           &
                     ,'FSW                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%fsn                                                           &
                     ,'FSN                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%fs_open                                                       &
                     ,'FS_OPEN                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%water_supply                                                  &
                     ,'WATER_SUPPLY              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_maintenance                                              &
                     ,'LEAF_MAINTENANCE          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%root_maintenance                                              &
                     ,'ROOT_MAINTENANCE          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_drop                                                     &
                     ,'LEAF_DROP                 ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%leaf_respiration                                              &
                     ,'LEAF_RESPIRATION          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%root_respiration                                              &
                     ,'ROOT_RESPIRATION          ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%gpp                                                           &
                     ,'GPP                       ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%paw_avg                                                       &
                     ,'PAW_AVG                   ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%elongf                                                        &
                     ,'ELONGF                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%turnover_amp                                                  &
                     ,'TURNOVER_AMP              ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%llspan                                                        &
                     ,'LLSPAN                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%vm_bar                                                        &
                     ,'VM_BAR                    ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%sla                                                           &
                     ,'SLA                       ',dsetrank,iparallel,.true. ,foundvar)
   !----- Daily means. --------------------------------------------------------------------!
   if (writing_long) then
      call hdf_getslab_r(cpatch%dmean_nppleaf                                              &
                        ,'DMEAN_NPPLEAF_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_nppfroot                                             &
                        ,'DMEAN_NPPFROOT_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_nppsapwood                                           &
                        ,'DMEAN_NPPSAPWOOD_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_nppcroot                                             &
                        ,'DMEAN_NPPCROOT_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_nppseeds                                             &
                        ,'DMEAN_NPPSEEDS_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_nppwood                                              &
                        ,'DMEAN_NPPWOOD_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_nppdaily                                             &
                        ,'DMEAN_NPPDAILY_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_gpp                                                  &
                        ,'DMEAN_GPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_npp                                                  &
                        ,'DMEAN_NPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_resp                                            &
                        ,'DMEAN_LEAF_RESP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_root_resp                                            &
                        ,'DMEAN_ROOT_RESP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_growth_resp                                     &
                        ,'DMEAN_LEAF_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_root_growth_resp                                     &
                        ,'DMEAN_ROOT_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_sapa_growth_resp                                     &
                        ,'DMEAN_SAPA_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_sapb_growth_resp                                     &
                        ,'DMEAN_SAPB_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_storage_resp                                    &
                        ,'DMEAN_LEAF_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_root_storage_resp                                    &
                        ,'DMEAN_ROOT_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_sapa_storage_resp                                    &
                        ,'DMEAN_SAPA_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_sapb_storage_resp                                    &
                        ,'DMEAN_SAPB_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_plresp                                               &
                        ,'DMEAN_PLRESP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_energy                                          &
                        ,'DMEAN_LEAF_ENERGY_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_water                                           &
                        ,'DMEAN_LEAF_WATER_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_hcap                                            &
                        ,'DMEAN_LEAF_HCAP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_vpdef                                           &
                        ,'DMEAN_LEAF_VPDEF_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_temp                                            &
                        ,'DMEAN_LEAF_TEMP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_fliq                                            &
                        ,'DMEAN_LEAF_FLIQ_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_gsw                                             &
                        ,'DMEAN_LEAF_GSW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_leaf_gbw                                             &
                        ,'DMEAN_LEAF_GBW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wood_energy                                          &
                        ,'DMEAN_WOOD_ENERGY_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wood_water                                           &
                        ,'DMEAN_WOOD_WATER_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wood_hcap                                            &
                        ,'DMEAN_WOOD_HCAP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wood_temp                                            &
                        ,'DMEAN_WOOD_TEMP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wood_fliq                                            &
                        ,'DMEAN_WOOD_FLIQ_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wood_gbw                                             &
                        ,'DMEAN_WOOD_GBW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_fs_open                                              &
                        ,'DMEAN_FS_OPEN_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_fsw                                                  &
                        ,'DMEAN_FSW_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_fsn                                                  &
                        ,'DMEAN_FSN_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_a_open                                               &
                        ,'DMEAN_A_OPEN_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_a_closed                                             &
                        ,'DMEAN_A_CLOSED_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_a_net                                                &
                        ,'DMEAN_A_NET_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_a_light                                              &
                        ,'DMEAN_A_LIGHT_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_a_rubp                                               &
                        ,'DMEAN_A_RUBP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_a_co2                                                &
                        ,'DMEAN_A_CO2_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_psi_open                                             &
                        ,'DMEAN_PSI_OPEN_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_psi_closed                                           &
                        ,'DMEAN_PSI_CLOSED_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_water_supply                                         &
                        ,'DMEAN_WATER_SUPPLY_CO     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_light_level                                          &
                        ,'DMEAN_LIGHT_LEVEL_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_light_level_beam                                     &
                        ,'DMEAN_LIGHT_LEVEL_BEAM_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_light_level_diff                                     &
                        ,'DMEAN_LIGHT_LEVEL_DIFF_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_par_l                                                &
                        ,'DMEAN_PAR_L_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_par_l_beam                                           &
                        ,'DMEAN_PAR_L_BEAM_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_par_l_diff                                           &
                        ,'DMEAN_PAR_L_DIFF_CO       ',dsetrank,iparallel,.false.,foundvar)

      call hdf_getslab_r(cpatch%dmean_par_level_beam                                       &
                        ,'DMEAN_PAR_LEVEL_BEAM_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_par_level_diffu                                       &
                        ,'DMEAN_PAR_LEVEL_DIFFU_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_par_level_diffd                                       &
                        ,'DMEAN_PAR_LEVEL_DIFFD_CO   ',dsetrank,iparallel,.false.,foundvar)


      call hdf_getslab_r(cpatch%dmean_rshort_l                                             &
                        ,'DMEAN_RSHORT_L_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_rlong_l                                              &
                        ,'DMEAN_RLONG_L_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_sensible_lc                                          &
                        ,'DMEAN_SENSIBLE_LC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_vapor_lc                                             &
                        ,'DMEAN_VAPOR_LC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_transp                                               &
                        ,'DMEAN_TRANSP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_intercepted_al                                       &
                        ,'DMEAN_INTERCEPTED_AL_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wshed_lg                                             &
                        ,'DMEAN_WSHED_LG_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_rshort_w                                             &
                        ,'DMEAN_RSHORT_W_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_rlong_w                                              &
                        ,'DMEAN_RLONG_W_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_sensible_wc                                          &
                        ,'DMEAN_SENSIBLE_WC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_vapor_wc                                             &
                        ,'DMEAN_VAPOR_WC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_intercepted_aw                                       &
                        ,'DMEAN_INTERCEPTED_AW_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%dmean_wshed_wg                                             &
                        ,'DMEAN_WSHED_WG_CO         ',dsetrank,iparallel,.false.,foundvar)
   end if
   !----- Daily means. --------------------------------------------------------------------!
   if (writing_eorq) then
      call hdf_getslab_r(cpatch%mmean_lai                                                  &
                        ,'MMEAN_LAI_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_bleaf                                                &
                        ,'MMEAN_BLEAF_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_broot                                                &
                        ,'MMEAN_BROOT_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_bstorage                                             &
                        ,'MMEAN_BSTORAGE_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_maintenance                                     &
                        ,'MMEAN_LEAF_MAINTENANCE_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_root_maintenance                                     &
                        ,'MMEAN_ROOT_MAINTENANCE_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_drop                                            &
                        ,'MMEAN_LEAF_DROP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_cb                                                   &
                        ,'MMEAN_CB_CO               ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_gpp                                                  &
                        ,'MMEAN_GPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_npp                                                  &
                        ,'MMEAN_NPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_resp                                            &
                        ,'MMEAN_LEAF_RESP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_root_resp                                            &
                        ,'MMEAN_ROOT_RESP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_growth_resp                                     &
                        ,'MMEAN_LEAF_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_root_growth_resp                                     &
                        ,'MMEAN_ROOT_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_sapa_growth_resp                                     &
                        ,'MMEAN_SAPA_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_sapb_growth_resp                                     &
                        ,'MMEAN_SAPB_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_storage_resp                                    &
                        ,'MMEAN_LEAF_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_root_storage_resp                                    &
                        ,'MMEAN_ROOT_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_sapa_storage_resp                                    &
                        ,'MMEAN_SAPA_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_sapb_storage_resp                                    &
                        ,'MMEAN_SAPB_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_plresp                                               &
                        ,'MMEAN_PLRESP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_energy                                          &
                        ,'MMEAN_LEAF_ENERGY_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_water                                           &
                        ,'MMEAN_LEAF_WATER_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_hcap                                            &
                        ,'MMEAN_LEAF_HCAP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_vpdef                                           &
                        ,'MMEAN_LEAF_VPDEF_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_temp                                            &
                        ,'MMEAN_LEAF_TEMP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_fliq                                            &
                        ,'MMEAN_LEAF_FLIQ_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_gsw                                             &
                        ,'MMEAN_LEAF_GSW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_leaf_gbw                                             &
                        ,'MMEAN_LEAF_GBW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wood_energy                                          &
                        ,'MMEAN_WOOD_ENERGY_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wood_water                                           &
                        ,'MMEAN_WOOD_WATER_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wood_hcap                                            &
                        ,'MMEAN_WOOD_HCAP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wood_temp                                            &
                        ,'MMEAN_WOOD_TEMP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wood_fliq                                            &
                        ,'MMEAN_WOOD_FLIQ_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wood_gbw                                             &
                        ,'MMEAN_WOOD_GBW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_fs_open                                              &
                        ,'MMEAN_FS_OPEN_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_fsw                                                  &
                        ,'MMEAN_FSW_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_fsn                                                  &
                        ,'MMEAN_FSN_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_a_open                                               &
                        ,'MMEAN_A_OPEN_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_a_closed                                             &
                        ,'MMEAN_A_CLOSED_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_a_net                                                &
                        ,'MMEAN_A_NET_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_a_light                                              &
                        ,'MMEAN_A_LIGHT_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_a_rubp                                               &
                        ,'MMEAN_A_RUBP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_a_co2                                                &
                        ,'MMEAN_A_CO2_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_psi_open                                             &
                        ,'MMEAN_PSI_OPEN_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_psi_closed                                           &
                        ,'MMEAN_PSI_CLOSED_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_water_supply                                         &
                        ,'MMEAN_WATER_SUPPLY_CO     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_light_level                                          &
                        ,'MMEAN_LIGHT_LEVEL_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_light_level_beam                                     &
                        ,'MMEAN_LIGHT_LEVEL_BEAM_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_light_level_diff                                     &
                        ,'MMEAN_LIGHT_LEVEL_DIFF_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_par_l                                                &
                        ,'MMEAN_PAR_L_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_par_l_beam                                           &
                        ,'MMEAN_PAR_L_BEAM_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_par_l_diff                                           &
                        ,'MMEAN_PAR_L_DIFF_CO       ',dsetrank,iparallel,.false.,foundvar)

      call hdf_getslab_r(cpatch%mmean_par_level_beam                                       &
                        ,'MMEAN_PAR_LEVEL_BEAM_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_par_level_diffu                                       &
                        ,'MMEAN_PAR_LEVEL_DIFFU_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_par_level_diffd                                       &
                        ,'MMEAN_PAR_LEVEL_DIFFD_CO   ',dsetrank,iparallel,.false.,foundvar)

      call hdf_getslab_r(cpatch%mmean_rshort_l                                             &
                        ,'MMEAN_RSHORT_L_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_rlong_l                                              &
                        ,'MMEAN_RLONG_L_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_sensible_lc                                          &
                        ,'MMEAN_SENSIBLE_LC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_vapor_lc                                             &
                        ,'MMEAN_VAPOR_LC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_transp                                               &
                        ,'MMEAN_TRANSP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_intercepted_al                                       &
                        ,'MMEAN_INTERCEPTED_AL_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wshed_lg                                             &
                        ,'MMEAN_WSHED_LG_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_rshort_w                                             &
                        ,'MMEAN_RSHORT_W_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_rlong_w                                              &
                        ,'MMEAN_RLONG_W_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_sensible_wc                                          &
                        ,'MMEAN_SENSIBLE_WC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_vapor_wc                                             &
                        ,'MMEAN_VAPOR_WC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_intercepted_aw                                       &
                        ,'MMEAN_INTERCEPTED_AW_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_wshed_wg                                             &
                        ,'MMEAN_WSHED_WG_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppleaf                                              &
                        ,'MMEAN_NPPLEAF_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppfroot                                             &
                        ,'MMEAN_NPPFROOT_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppsapwood                                           &
                        ,'MMEAN_NPPSAPWOOD_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppcroot                                             &
                        ,'MMEAN_NPPCROOT_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppseeds                                             &
                        ,'MMEAN_NPPSEEDS_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppwood                                              &
                        ,'MMEAN_NPPWOOD_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmean_nppdaily                                             &
                        ,'MMEAN_NPPDAILY_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_gpp                                                  &
                        ,'MMSQU_GPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_npp                                                  &
                        ,'MMSQU_NPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_plresp                                               &
                        ,'MMSQU_PLRESP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_sensible_lc                                          &
                        ,'MMSQU_SENSIBLE_LC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_vapor_lc                                             &
                        ,'MMSQU_VAPOR_LC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_transp                                               &
                        ,'MMSQU_TRANSP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_sensible_wc                                          &
                        ,'MMSQU_SENSIBLE_WC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%mmsqu_vapor_wc                                             &
                        ,'MMSQU_VAPOR_WC_CO         ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
  




   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (13 ,ncohorts).                                       !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = 13_8
   chnkdims(1) = 13_8
   chnkoffs(1) = 0_8
   memdims (1) = 13_8
   memsize (1) = 13_8
   memoffs (1) = 0_8
   globdims(2) = int(ncohorts_global,8)
   chnkdims(2) = int(cpatch%ncohorts,8)
   chnkoffs(2) = int(paco_index - 1 ,8)
   memdims (2) = int(cpatch%ncohorts,8)
   memsize (2) = int(cpatch%ncohorts,8)
   memoffs (2) = 0_8
   call hdf_getslab_r(cpatch%cb                                                            &
                     ,'CB                        ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%cb_lightmax                                                   &
                     ,'CB_LIGHTMAX               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%cb_moistmax                                                   &
                     ,'CB_MOISTMAX               ',dsetrank,iparallel,.true. ,foundvar)
   call hdf_getslab_r(cpatch%cb_mlmax                                                      &
                     ,'CB_MLMAX                  ',dsetrank,iparallel,.true. ,foundvar)
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (n_mort,ncohorts).                                    !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(n_mort,8)
   chnkdims(1) = int(n_mort,8)
   chnkoffs(1) = 0_8
   memdims (1) = int(n_mort,8)
   memsize (1) = int(n_mort,8)
   memoffs (1) = 0_8
   
   globdims(2) = int(ncohorts_global,8)
   chnkdims(2) = int(cpatch%ncohorts,8)
   chnkoffs(2) = int(paco_index - 1,8)
   memdims (2) = int(cpatch%ncohorts,8)
   memsize (2) = int(cpatch%ncohorts,8)
   memoffs (2) = 0_8
   call hdf_getslab_r(cpatch%mort_rate                                                     &
                     ,'MORT_RATE_CO              ',dsetrank,iparallel,.true. ,foundvar)
   if (writing_eorq) then
      call hdf_getslab_r(cpatch%mmean_mort_rate                                            &
                        ,'MMEAN_MORT_RATE_CO        ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (n_radprof,ncohorts).                                 !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(n_radprof,8)
   chnkdims(1) = int(n_radprof,8)
   chnkoffs(1) = 0_8
   memdims (1) = int(n_radprof,8)
   memsize (1) = int(n_radprof,8)
   memoffs (1) = 0_8
   
   globdims(2) = int(ncohorts_global,8)
   chnkdims(2) = int(cpatch%ncohorts,8)
   chnkoffs(2) = int(paco_index - 1,8)
   memdims (2) = int(cpatch%ncohorts,8)
   memsize (2) = int(cpatch%ncohorts,8)
   memoffs (2) = 0_8
   call hdf_getslab_r(cpatch%rad_profile                                                   &
                     ,'RAD_PROFILE_CO              ',dsetrank,iparallel,.true. ,foundvar)
   if (writing_long) then
      call hdf_getslab_r(cpatch%dmean_rad_profile                                          &
                        ,'DMEAN_RAD_PROFILE_CO     ',dsetrank,iparallel,.false.,foundvar)
   end if
   if (writing_eorq) then
      call hdf_getslab_r(cpatch%mmean_rad_profile                                          &
                        ,'MMEAN_RAD_PROFILE_CO     ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      2-D variables, dimensions: (ndcycle,ncohorts).                                   !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 2
   globdims(1) = int(ndcycle,8)
   chnkdims(1) = int(ndcycle,8)
   chnkoffs(1) = 0_8
   memdims (1) = int(ndcycle,8)
   memsize (1) = int(ndcycle,8)
   memoffs (1) = 0_8
   globdims(2) = int(ncohorts_global,8)
   chnkdims(2) = int(cpatch%ncohorts,8)
   chnkoffs(2) = int(paco_index - 1 ,8)
   memdims (2) = int(cpatch%ncohorts,8)
   memsize (2) = int(cpatch%ncohorts,8)
   memoffs (2) = 0_8
   !----- Mean diel. ----------------------------------------------------------------------!
   if (writing_dcyc) then
      call hdf_getslab_r(cpatch%qmean_gpp                                                  &
                        ,'QMEAN_GPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_npp                                                  &
                        ,'QMEAN_NPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_resp                                            &
                        ,'QMEAN_LEAF_RESP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_root_resp                                            &
                        ,'QMEAN_ROOT_RESP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_growth_resp                                     &
                        ,'QMEAN_LEAF_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_root_growth_resp                                     &
                        ,'QMEAN_ROOT_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_sapa_growth_resp                                     &
                        ,'QMEAN_SAPA_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_sapb_growth_resp                                     &
                        ,'QMEAN_SAPB_GROWTH_RESP_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_storage_resp                                    &
                        ,'QMEAN_LEAF_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_root_storage_resp                                    &
                        ,'QMEAN_ROOT_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_sapa_storage_resp                                    &
                        ,'QMEAN_SAPA_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_sapb_storage_resp                                    &
                        ,'QMEAN_SAPB_STORAGE_RESP_CO',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_plresp                                               &
                        ,'QMEAN_PLRESP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_energy                                          &
                        ,'QMEAN_LEAF_ENERGY_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_water                                           &
                        ,'QMEAN_LEAF_WATER_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_hcap                                            &
                        ,'QMEAN_LEAF_HCAP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_vpdef                                           &
                        ,'QMEAN_LEAF_VPDEF_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_temp                                            &
                        ,'QMEAN_LEAF_TEMP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_fliq                                            &
                        ,'QMEAN_LEAF_FLIQ_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_gsw                                             &
                        ,'QMEAN_LEAF_GSW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_leaf_gbw                                             &
                        ,'QMEAN_LEAF_GBW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wood_energy                                          &
                        ,'QMEAN_WOOD_ENERGY_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wood_water                                           &
                        ,'QMEAN_WOOD_WATER_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wood_hcap                                            &
                        ,'QMEAN_WOOD_HCAP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wood_temp                                            &
                        ,'QMEAN_WOOD_TEMP_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wood_fliq                                            &
                        ,'QMEAN_WOOD_FLIQ_CO        ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wood_gbw                                             &
                        ,'QMEAN_WOOD_GBW_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_fs_open                                              &
                        ,'QMEAN_FS_OPEN_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_fsw                                                  &
                        ,'QMEAN_FSW_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_fsn                                                  &
                        ,'QMEAN_FSN_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_a_open                                               &
                        ,'QMEAN_A_OPEN_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_a_closed                                             &
                        ,'QMEAN_A_CLOSED_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_a_net                                                &
                        ,'QMEAN_A_NET_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_a_light                                              &
                        ,'QMEAN_A_LIGHT_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_a_rubp                                               &
                        ,'QMEAN_A_RUBP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_a_co2                                                &
                        ,'QMEAN_A_CO2_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_psi_open                                             &
                        ,'QMEAN_PSI_OPEN_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_psi_closed                                           &
                        ,'QMEAN_PSI_CLOSED_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_water_supply                                         &
                        ,'QMEAN_WATER_SUPPLY_CO     ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_light_level                                          &
                        ,'QMEAN_LIGHT_LEVEL_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_light_level_beam                                     &
                        ,'QMEAN_LIGHT_LEVEL_BEAM_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_light_level_diff                                     &
                        ,'QMEAN_LIGHT_LEVEL_DIFF_CO ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_par_l                                                &
                        ,'QMEAN_PAR_L_CO            ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_par_l_beam                                           &
                        ,'QMEAN_PAR_L_BEAM_CO       ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_par_l_diff                                           &
                        ,'QMEAN_PAR_L_DIFF_CO       ',dsetrank,iparallel,.false.,foundvar)

      call hdf_getslab_r(cpatch%qmean_par_level_beam                                       &
                        ,'QMEAN_PAR_LEVEL_BEAM_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_par_level_diffu                                       &
                        ,'QMEAN_PAR_LEVEL_DIFFU_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_par_level_diffd                                       &
                        ,'QMEAN_PAR_LEVEL_DIFFD_CO   ',dsetrank,iparallel,.false.,foundvar)


      call hdf_getslab_r(cpatch%qmean_rshort_l                                             &
                        ,'QMEAN_RSHORT_L_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_rlong_l                                              &
                        ,'QMEAN_RLONG_L_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_sensible_lc                                          &
                        ,'QMEAN_SENSIBLE_LC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_vapor_lc                                             &
                        ,'QMEAN_VAPOR_LC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_transp                                               &
                        ,'QMEAN_TRANSP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_intercepted_al                                       &
                        ,'QMEAN_INTERCEPTED_AL_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wshed_lg                                             &
                        ,'QMEAN_WSHED_LG_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_rshort_w                                             &
                        ,'QMEAN_RSHORT_W_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_rlong_w                                              &
                        ,'QMEAN_RLONG_W_CO          ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_sensible_wc                                          &
                        ,'QMEAN_SENSIBLE_WC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_vapor_wc                                             &
                        ,'QMEAN_VAPOR_WC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_intercepted_aw                                       &
                        ,'QMEAN_INTERCEPTED_AW_CO   ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmean_wshed_wg                                             &
                        ,'QMEAN_WSHED_WG_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_gpp                                                  &
                        ,'QMSQU_GPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_npp                                                  &
                        ,'QMSQU_NPP_CO              ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_plresp                                               &
                        ,'QMSQU_PLRESP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_sensible_lc                                          &
                        ,'QMSQU_SENSIBLE_LC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_vapor_lc                                             &
                        ,'QMSQU_VAPOR_LC_CO         ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_transp                                               &
                        ,'QMSQU_TRANSP_CO           ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_sensible_wc                                          &
                        ,'QMSQU_SENSIBLE_WC_CO      ',dsetrank,iparallel,.false.,foundvar)
      call hdf_getslab_r(cpatch%qmsqu_vapor_wc                                             &
                        ,'QMSQU_VAPOR_WC_CO         ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      3-D variables, dimensions: (n_radprof,ndcycle,ncohorts).                                 !
   !---------------------------------------------------------------------------------------!
   dsetrank    = 3
   globdims(1) = int(n_radprof,8)
   chnkdims(1) = int(n_radprof,8)
   chnkoffs(1) = 0_8
   memdims (1) = int(n_radprof,8)
   memsize (1) = int(n_radprof,8)
   memoffs (1) = 0_8
   globdims(2) = int(ndcycle,8)
   chnkdims(2) = int(ndcycle,8)
   chnkoffs(2) = 0_8
   memdims (2) = int(ndcycle,8)
   memsize (2) = int(ndcycle,8)
   memoffs (2) = 0_8
   globdims(3) = int(ncohorts_global,8)
   chnkdims(3) = int(cpatch%ncohorts,8)
   chnkoffs(3) = int(paco_index - 1,8)
   memdims (3) = int(cpatch%ncohorts,8)
   memsize (3) = int(cpatch%ncohorts,8)
   memoffs (3) = 0_8
   if (writing_dcyc) then
      call hdf_getslab_r(cpatch%qmean_rad_profile                                          &
                        ,'QMEAN_RAD_PROFILE_CO     ',dsetrank,iparallel,.false.,foundvar)
   end if
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!


  return
end subroutine fill_history_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine reads in real variables from HDF5 files.                               !
! BUFF      -- Where the dataset should be loaded to                                       !
! VARN      -- Variable name as in the HDF5 file                                           !
! DSETRANK  -- Number of dimensions.                                                       !
! IPARALLEL -- Should fortran read data in parallel?  (0 = no, 1 = yes)                    !
! REQUIRED  -- Is this variable required? (.true. = yes, .false. = no)                     !
!              If a variable is missing and required is set to true, then the model will   !
!              crash.  If it is missing but required is false, then we initialise the      !
!              variable with zeroes and print a big banner to the standard output to       !
!              warn the user.                                                              !
! FOUNDVAR  -- Output flag that tells whether the variable was found or not.               !
!------------------------------------------------------------------------------------------!
subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required,foundvar)
   
   use hdf5
   use hdf5_coms, only : file_id      & ! intent(inout)
                       , dset_id      & ! intent(inout)
                       , dspace_id    & ! intent(inout)
                       , plist_id     & ! intent(inout)
                       , filespace    & ! intent(inout)
                       , memspace     & ! intent(inout)
                       , globdims     & ! intent(inout)
                       , chnkdims     & ! intent(inout)
                       , chnkoffs     & ! intent(inout)
                       , cnt          & ! intent(inout)
                       , stride       & ! intent(inout)
                       , memdims      & ! intent(inout)
                       , memoffs      & ! intent(inout)
                       , memsize        ! intent(inout)

   use ed_misc_coms,only: suppress_h5_warnings

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))                &
                                                              , intent(inout) :: buff
   character(len=*)                                           , intent(in)    :: varn
   integer                                                    , intent(in)    :: dsetrank
   integer                                                    , intent(in)    :: iparallel
   logical                                                    , intent(in)    :: required
   logical                                                    , intent(out)   :: foundvar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                                    :: hdferr
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Try to open dataset, and save the success/failure flag.                           !
   !---------------------------------------------------------------------------------------!
   call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
   foundvar = hdferr == 0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Check whether the dataset was opened before continuing.                            !
   !---------------------------------------------------------------------------------------!
   if ((.not. foundvar) .and. required ) then
      !------------------------------------------------------------------------------------!
      !    Variable is required but wasn't found; stop the run.                            !
      !------------------------------------------------------------------------------------!
      write(unit=*,fmt=*) 'File_ID = ',file_id
      write(unit=*,fmt=*) 'Dset_ID = ',dset_id
      call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
           ,'hdf_getslab_r','ed_init_full_history.F90')
      !------------------------------------------------------------------------------------!

   else if ((.not. foundvar) .and. (.not.required) ) then
      !------------------------------------------------------------------------------------!
      !    Variable wasn't found but it wasn't required either; initialise buffer with     !
      ! zeroes, and warn the user that we are doing this.                                  !
      !------------------------------------------------------------------------------------!

      if(.not.suppress_h5_warnings)then
      write(unit=*,fmt=*) 'File_ID = ',file_id
      write(unit=*,fmt=*) 'Dset_ID = ',dset_id
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!  '
      write (unit=*,fmt='(a)') '                                                          '
      write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
      write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
      write (unit=*,fmt='(a)') ' + his may cause some of your diagnostic output related'
      write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
      write (unit=*,fmt='(a)') ''
      write (unit=*,fmt='(a)') '   This variable has been specified as:'
      write (unit=*,fmt='(a)') '   NOT ABSOUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') ''
      end if
      
      buff(:,:,:,:) = 0.
      return
      !------------------------------------------------------------------------------------!

   else
      !------------------------------------------------------------------------------------!
      !     Found the variable, read in the data, checking that every step was success-    !
      ! fully done before moving to the next.                                              !
      !------------------------------------------------------------------------------------!
      call h5dget_space_f(dset_id,filespace,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!'   &
              ,'hdf_getslab_r','ed_init_full_history.F90')
      end if
      
      call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
           chnkdims,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Couldn''t assign the hyperslab filespace for '//trim(varn)//'!' &
              ,'hdf_getslab_r','ed_init_full_history.F90')
      end if
      
      call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
         write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
         call fatal_error('Couldn''t create the hyperslab memspace for '//trim(varn)//'!'  &
                         ,'hdf_getslab_r','ed_init_full_history.F90')
      end if
      
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
           memdims,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Couldn''t assign the hyperslab filespace for '//trim(varn)//'!' &
              ,'hdf_getslab_r','ed_init_full_history.F90')
      end if
      
      if (iparallel == 1) then
         
         call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
              mem_space_id = memspace, file_space_id = filespace, &
              xfer_prp = plist_id)
         
         if (hdferr /= 0) then
            call fatal_error('Couldn''t read in hyperslab dataset for '//trim(varn)//'!'   &
                 ,'hdf_getslab_r','ed_init_full_history.F90')
         end if

      else

         call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
              mem_space_id = memspace, file_space_id = filespace )

         if (hdferr /= 0) then
            call fatal_error('Couldn''t read in hyperslab dataset for '//trim(varn)//'!'   &
                 ,'hdf_getslab_r','ed_init_full_history.F90')
         end if

      end if
      
      !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
      
      call h5sclose_f(filespace, hdferr)
      call h5sclose_f(memspace , hdferr)
      call h5dclose_f(dset_id  , hdferr)
      
   end if
   return
end subroutine hdf_getslab_r
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine reads in double precision variables from HDF5 files.                   !
! BUFF      -- Where the dataset should be loaded to                                       !
! VARN      -- Variable name as in the HDF5 file                                           !
! DSETRANK  -- Number of dimensions.                                                       !
! IPARALLEL -- Should fortran read data in parallel?  (0 = no, 1 = yes)                    !
! REQUIRED  -- Is this variable required? (.true. = yes, .false. = no)                     !
!              If a variable is missing and required is set to true, then the model will   !
!              crash.  If it is missing but required is false, then we initialise the      !
!              variable with zeroes and print a big banner to the standard output to       !
!              warn the user.                                                              !
! FOUNDVAR  -- Output flag that tells whether the variable was found or not.               !
!------------------------------------------------------------------------------------------!
subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required,foundvar)
   
   use hdf5
   use hdf5_coms, only : file_id      & ! intent(inout)
                       , dset_id      & ! intent(inout)
                       , dspace_id    & ! intent(inout)
                       , plist_id     & ! intent(inout)
                       , filespace    & ! intent(inout)
                       , memspace     & ! intent(inout)
                       , globdims     & ! intent(inout)
                       , chnkdims     & ! intent(inout)
                       , chnkoffs     & ! intent(inout)
                       , cnt          & ! intent(inout)
                       , stride       & ! intent(inout)
                       , memdims      & ! intent(inout)
                       , memoffs      & ! intent(inout)
                       , memsize      ! ! intent(inout)

   use ed_misc_coms,only: suppress_h5_warnings
   
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=8)    , dimension(memsize(1),memsize(2),memsize(3),memsize(4))                &
                                                              , intent(inout) :: buff
   character(len=*)                                           , intent(in)    :: varn
   integer                                                    , intent(in)    :: dsetrank
   integer                                                    , intent(in)    :: iparallel
   logical                                                    , intent(in)    :: required
   logical                                                    , intent(out)   :: foundvar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                                    :: hdferr
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Try to open dataset, and save the success/failure flag.                           !
   !---------------------------------------------------------------------------------------!
   call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
   foundvar = hdferr == 0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Check whether the dataset was opened before continuing.                            !
   !---------------------------------------------------------------------------------------!
   if ((.not. foundvar) .and. required ) then
      !------------------------------------------------------------------------------------!
      !    Variable is required but wasn't found; stop the run.                            !
      !------------------------------------------------------------------------------------!
      write(unit=*,fmt=*) 'File_ID = ',file_id
      write(unit=*,fmt=*) 'Dset_ID = ',dset_id
      call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
           ,'hdf_getslab_d','ed_init_full_history.F90')
      !------------------------------------------------------------------------------------!

   else if ((.not. foundvar) .and. (.not.required) ) then
      !------------------------------------------------------------------------------------!
      !    Variable wasn't found but it wasn't required either; initialise buffer with     !
      ! zeroes, and warn the user that we are doing this.                                  !
      !------------------------------------------------------------------------------------!
      if(.not.suppress_h5_warnings)then
      write(unit=*,fmt=*) 'File_ID = ',file_id
      write(unit=*,fmt=*) 'Dset_ID = ',dset_id
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!  '
      write (unit=*,fmt='(a)') '                                                          '
      write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
      write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
      write (unit=*,fmt='(a)') ' + his may cause some of your diagnostic output related'
      write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
      write (unit=*,fmt='(a)') ''
      write (unit=*,fmt='(a)') '   This variable has been specified as:'
      write (unit=*,fmt='(a)') '   NOT ABSOUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') ''
      end if
      buff(:,:,:,:) = 0.d0
      return
      !------------------------------------------------------------------------------------!

   else
      !------------------------------------------------------------------------------------!
      !     Found the variable, read in the data, checking that every step was success-    !
      ! fully done before moving to the next.                                              !
      !------------------------------------------------------------------------------------!
      call h5dget_space_f(dset_id,filespace,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!'   &
              ,'hdf_getslab_d','ed_init_full_history.F90')
      end if
      
      call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
           chnkdims,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Couldn''t assign the hyperslab filespace for '//trim(varn)//'!' &
              ,'hdf_getslab_d','ed_init_full_history.F90')
      end if
      
      call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
         write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
         call fatal_error('Couldn''t create the hyperslab memspace for '//trim(varn)//'!'  &
                         ,'hdf_getslab_d','ed_init_full_history.F90')
      end if
      
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
           memdims,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Couldn''t assign the hyperslab filespace for '//trim(varn)//'!' &
              ,'hdf_getslab_d','ed_init_full_history.F90')
      end if
      
      if (iparallel == 1) then
         
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,buff,globdims, hdferr, &
              mem_space_id = memspace, file_space_id = filespace, &
              xfer_prp = plist_id)
         
         if (hdferr /= 0) then
            call fatal_error('Couldn''t read in hyperslab dataset for '//trim(varn)//'!'   &
                 ,'hdf_getslab_d','ed_init_full_history.F90')
         end if

      else

         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,buff,globdims, hdferr, &
              mem_space_id = memspace, file_space_id = filespace )

         if (hdferr /= 0) then
            call fatal_error('Couldn''t read in hyperslab dataset for '//trim(varn)//'!'   &
                 ,'hdf_getslab_d','ed_init_full_history.F90')
         end if

      end if
      
      !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
      
      call h5sclose_f(filespace, hdferr)
      call h5sclose_f(memspace , hdferr)
      call h5dclose_f(dset_id  , hdferr)
      
   end if
   return
end subroutine hdf_getslab_d
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine reads in integer variables from HDF5 files.                            !
! BUFF      -- Where the dataset should be loaded to                                       !
! VARN      -- Variable name as in the HDF5 file                                           !
! DSETRANK  -- Number of dimensions.                                                       !
! IPARALLEL -- Should fortran read data in parallel?  (0 = no, 1 = yes)                    !
! REQUIRED  -- Is this variable required? (.true. = yes, .false. = no)                     !
!              If a variable is missing and required is set to true, then the model will   !
!              crash.  If it is missing but required is false, then we initialise the      !
!              variable with zeroes and print a big banner to the standard output to       !
!              warn the user.                                                              !
! FOUNDVAR  -- Output flag that tells whether the variable was found or not.               !
!------------------------------------------------------------------------------------------!
subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required,foundvar)
   
   use hdf5
   use hdf5_coms, only : file_id      & ! intent(inout)
                       , dset_id      & ! intent(inout)
                       , dspace_id    & ! intent(inout)
                       , plist_id     & ! intent(inout)
                       , filespace    & ! intent(inout)
                       , memspace     & ! intent(inout)
                       , globdims     & ! intent(inout)
                       , chnkdims     & ! intent(inout)
                       , chnkoffs     & ! intent(inout)
                       , cnt          & ! intent(inout)
                       , stride       & ! intent(inout)
                       , memdims      & ! intent(inout)
                       , memoffs      & ! intent(inout)
                       , memsize      ! ! intent(inout)

   use ed_misc_coms,only: suppress_h5_warnings
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer         , dimension(memsize(1),memsize(2),memsize(3),memsize(4))                &
                                                              , intent(inout) :: buff
   character(len=*)                                           , intent(in)    :: varn
   integer                                                    , intent(in)    :: dsetrank
   integer                                                    , intent(in)    :: iparallel
   logical                                                    , intent(in)    :: required
   logical                                                    , intent(out)   :: foundvar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                                    :: hdferr
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Try to open dataset, and save the success/failure flag.                           !
   !---------------------------------------------------------------------------------------!
   call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
   foundvar = hdferr == 0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Check whether the dataset was opened before continuing.                            !
   !---------------------------------------------------------------------------------------!
   if ((.not. foundvar) .and. required ) then
      !------------------------------------------------------------------------------------!
      !    Variable is required but wasn't found; stop the run.                            !
      !------------------------------------------------------------------------------------!
      write(unit=*,fmt=*) 'File_ID = ',file_id
      write(unit=*,fmt=*) 'Dset_ID = ',dset_id
      call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
           ,'hdf_getslab_i','ed_init_full_history.F90')
      !------------------------------------------------------------------------------------!

   else if ((.not. foundvar) .and. (.not.required) ) then
      !------------------------------------------------------------------------------------!
      !    Variable wasn't found but it wasn't required either; initialise buffer with     !
      ! zeroes, and warn the user that we are doing this.                                  !
      !------------------------------------------------------------------------------------!
      if(.not.suppress_h5_warnings)then
      write(unit=*,fmt=*) 'File_ID = ',file_id
      write(unit=*,fmt=*) 'Dset_ID = ',dset_id
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!  '
      write (unit=*,fmt='(a)') '                                                          '
      write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
      write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
      write (unit=*,fmt='(a)') ' + his may cause some of your diagnostic output related'
      write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
      write (unit=*,fmt='(a)') ''
      write (unit=*,fmt='(a)') '   This variable has been specified as:'
      write (unit=*,fmt='(a)') '   NOT ABSOUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') ''
      end if
      buff(:,:,:,:) = 0
      return
      !------------------------------------------------------------------------------------!

   else
      !------------------------------------------------------------------------------------!
      !     Found the variable, read in the data, checking that every step was success-    !
      ! fully done before moving to the next.                                              !
      !------------------------------------------------------------------------------------!
      call h5dget_space_f(dset_id,filespace,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!'   &
              ,'hdf_getslab_i','ed_init_full_history.F90')
      end if
      
      call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
           chnkdims,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Couldn''t assign the hyperslab filespace for '//trim(varn)//'!' &
              ,'hdf_getslab_i','ed_init_full_history.F90')
      end if
      
      call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
         write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
         call fatal_error('Couldn''t create the hyperslab memspace for '//trim(varn)//'!'  &
                         ,'hdf_getslab_i','ed_init_full_history.F90')
      end if
      
      call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
           memdims,hdferr)
      if (hdferr /= 0) then
         call fatal_error('Couldn''t assign the hyperslab filespace for '//trim(varn)//'!' &
              ,'hdf_getslab_i','ed_init_full_history.F90')
      end if
      
      if (iparallel == 1) then
         
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,buff,globdims, hdferr, &
              mem_space_id = memspace, file_space_id = filespace, &
              xfer_prp = plist_id)
         
         if (hdferr /= 0) then
            call fatal_error('Couldn''t read in hyperslab dataset for '//trim(varn)//'!'   &
                 ,'hdf_getslab_i','ed_init_full_history.F90')
         end if

      else

         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,buff,globdims, hdferr, &
              mem_space_id = memspace, file_space_id = filespace )

         if (hdferr /= 0) then
            call fatal_error('Couldn''t read in hyperslab dataset for '//trim(varn)//'!'   &
                 ,'hdf_getslab_i','ed_init_full_history.F90')
         end if

      end if
      
      !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
      
      call h5sclose_f(filespace, hdferr)
      call h5sclose_f(memspace , hdferr)
      call h5dclose_f(dset_id  , hdferr)
      
   end if
   return
end subroutine hdf_getslab_i
!==========================================================================================!
!==========================================================================================!
