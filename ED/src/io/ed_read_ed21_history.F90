!==========================================================================================!
!==========================================================================================!
!     Read the ED-2.1 history.  This will assume that we are reading from a single source, !
! like one regional run, or one polygon-of-interest collection run.  More than one grid is !
! allowed, but only if they came from the same simulation.  For a completely unstructured  !
! restart run, you may want to use ied_init_mode=5 instead (see next subroutine).          !
!------------------------------------------------------------------------------------------!
subroutine read_ed21_history_file
#if USE_HDF5
   use hdf5
#endif
   use ed_max_dims    , only : n_pft                   & ! intent(in)
                             , huge_polygon            & ! intent(in)
                             , str_len                 ! ! intent(in)
   use pft_coms       , only : SLA                     & ! intent(in)
                             , q                       & ! intent(in)
                             , qsw                     & ! intent(in)
                             , hgt_min                 & ! intent(in)
                             , min_dbh                 & ! intent(in)
                             , min_bdead               & ! intent(in)
                             , is_grass                & ! intent(in)
                             , include_pft             & ! intent(in)
                             , include_pft_ag          & ! intent(in)
                             , pft_1st_check           & ! intent(in)
                             , agf_bs                  & ! intent(in)
                             , include_these_pft       ! ! intent(in)
   use ed_misc_coms   , only : sfilin                  & ! intent(in)
                             , current_time            & ! intent(in)
                             , imonthh                 & ! intent(in)
                             , iyearh                  & ! intent(in)
                             , idateh                  & ! intent(in)
                             , igrass                  & ! intent(in)
                             , itimeh                  ! ! intent(in)
   use ed_state_vars  , only : polygontype             & ! variable type
                             , sitetype                & ! variable type
                             , patchtype               & ! variable type
                             , edtype                  & ! variable type
                             , edgrid_g                & ! variable type
                             , allocate_polygontype    & ! subroutine
                             , allocate_sitetype       & ! subroutine
                             , allocate_patchtype      ! ! subroutine
   use grid_coms      , only : ngrids                  & ! intent(in)
                             , nzg                     ! ! intent(in)
   use consts_coms    , only : pio4                    ! ! intent(in)
   use hdf5_coms      , only : file_id                 & ! intent(in)
                             , dset_id                 & ! intent(in)
                             , dspace_id               & ! intent(in)
                             , plist_id                & ! intent(in)
                             , globdims                & ! intent(in)
                             , chnkdims                & ! intent(in)
                             , chnkoffs                & ! intent(in)
                             , memdims                 & ! intent(in)
                             , memoffs                 & ! intent(in)
                             , memsize                 ! ! intent(in)
   use allometry      , only : area_indices            & ! function
                             , ed_biomass              & ! function
                             , bd2dbh                  & ! function
                             , size2bl                 & ! function
                             , dbh2h                   & ! function
                             , dbh2bd                  ! ! function
   use fuse_fiss_utils, only : terminate_cohorts       ! ! subroutine
   use disturb_coms   , only : ianth_disturb           ! ! intent(in)

   implicit none
#if USE_HDF5

   !------ Local variables. ---------------------------------------------------------------!
   type(edtype)          , pointer     :: cgrid
   type(polygontype)     , pointer     :: cpoly
   type(sitetype)        , pointer     :: csite
   type(patchtype)       , pointer     :: cpatch
   character(len=3)                    :: cgr
   character(len=str_len)              :: hnamel
   
   integer, dimension(:) , allocatable :: pysi_n
   integer, dimension(:) , allocatable :: pysi_id
   integer, dimension(:) , allocatable :: sipa_n
   integer, dimension(:) , allocatable :: sipa_id
   integer, dimension(:) , allocatable :: paco_n
   integer, dimension(:) , allocatable :: paco_id
   integer, dimension(:) , allocatable :: islakesite
   integer, dimension(:) , allocatable :: plantation
   integer                             :: year
   integer                             :: igr
   integer                             :: ipy
   integer                             :: isi
   integer                             :: is
   integer                             :: ipa
   integer                             :: ico
   integer                             :: k
   integer                             :: dset_npolygons_global
   integer                             :: dset_nsites_global
   integer                             :: dset_npatches_global
   integer                             :: dset_ncohorts_global
   integer                             :: dset_nzg
   integer                             :: hdferr
   integer                             :: ngr
   integer                             :: ifpy
   integer                             :: ipft
   integer                             :: py_index
   integer                             :: si_index
   integer                             :: pa_index
   integer                             :: dsetrank
   integer                             :: iparallel
   integer                             :: ndry_sites
   logical                             :: exists
   real, dimension(:)    , allocatable :: file_lats
   real, dimension(:)    , allocatable :: file_lons
   real                                :: minrad
   real                                :: currad
   real                                :: elim_nplant
   real                                :: elim_lai
   real                                :: salloc
   real                                :: salloci
   logical                             :: foundvar
   !----- Local constants. ----------------------------------------------------------------!
   real                  , parameter   :: tiny_biomass = 1.e-20
   !----- External function. --------------------------------------------------------------!
   real                  , external    :: dist_gc
   !---------------------------------------------------------------------------------------!



   !----- Open the HDF environment. -------------------------------------------------------!
   call h5open_f(hdferr)

   !----- Turn off automatic error printing. ----------------------------------------------!
   call h5eset_auto_f(0,hdferr)


   !----- Initialize the dimensional control variables for the H5 slabs. ------------------!
   globdims = 0_8
   chnkdims = 0_8
   chnkoffs = 0_8
   memoffs  = 0_8
   memdims  = 0_8
   memsize  = 1_8  


   !---------------------------------------------------------------------------------------!
   !     Walk the tree and pull data from the dataset.                                     !
   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids
      
      cgrid => edgrid_g(igr)


      !----- Write the grid number as a character. ----------------------------------------!
      write(cgr,'(a1,i2.2)') 'g',igr

      !------------------------------------------------------------------------------------!
      !     The file name should be the exact file that will be read in, except the part   !
      ! that lists the grid number and extension.  Because this is a history restart, only !
      ! the first SFILIN will be used so we ensure that we will continue from the same     !
      ! simulation.                                                                        !
      !------------------------------------------------------------------------------------!
      hnamel = trim(sfilin(1))//"-"//trim(cgr)//".h5"


      !------------------------------------------------------------------------------------!
      !     Check whether the file we want to open actually exists.                        !
      !------------------------------------------------------------------------------------!
      inquire(file=trim(hnamel),exist=exists)
      if (.not.exists) then
         !----- It doesn't, stop the run. -------------------------------------------------!
         write (unit=*,fmt='(a,1x,a)')    'SFILIN(1) = ',trim(sfilin(1))
         write (unit=*,fmt='(a,1x,i4.4)') 'IYEARH    = ',iyearh
         write (unit=*,fmt='(a,1x,i2.2)') 'IMONTHH   = ',imonthh
         write (unit=*,fmt='(a,1x,i2.2)') 'IDATEH    = ',idateh
         write (unit=*,fmt='(a,1x,i4.4)') 'ITIMEH    = ',itimeh
         call fatal_error ('File '//trim(hnamel)//' not found.'                            &
                          ,'read_ed21_history_file','ed_read_ed21_history.F90')
      else
         !----- File is there, read it. ---------------------------------------------------!
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a,1x,i4)') ' - Error opening HDF5 file - error - ',hdferr
            write(unit=*,fmt='(a,1x,a)') ' - File name: ',trim(hnamel)
            call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)            &
                            ,'read_ed21_history_file','ed_read_ed21_history.F90')
         end if
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Retrieve global vector sizes and mapping tree.                                !
      !------------------------------------------------------------------------------------!
      globdims    = 0_8
      chnkdims    = 0_8
      chnkoffs    = 0_8
      globdims(1) = 1_8

      call h5dopen_f(file_id,'NZG', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nzg,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nsites_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npatches_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_ncohorts_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      globdims(1) = int(dset_npolygons_global,8)

      allocate(pysi_n(dset_npolygons_global))
      allocate(pysi_id(dset_npolygons_global))
      
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
      
      globdims(1) = int(dset_nsites_global,8)
      
      allocate(sipa_n(dset_nsites_global))
      allocate(sipa_id(dset_nsites_global))
      
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
      
      globdims(1) = int(dset_npatches_global,8)
      allocate(paco_n(dset_npatches_global))
      allocate(paco_id(dset_npatches_global))
      
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

      !----- Retrieve the polygon coordinates data. ---------------------------------------!
      globdims(1) = int(dset_npolygons_global,8)
      allocate(file_lats(dset_npolygons_global))
      allocate(file_lons(dset_npolygons_global))
      
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
      !     Loop over all the input polygons.                                              !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !----- Initialise polygon map and start closes polygon with a huge distance. -----!
         py_index = 0
         minrad = 1.e20
         
         !---------------------------------------------------------------------------------!
         !     Loop over all input polygons and find the one that is the closest to this   !
         ! current polygon.                                                                !
         !---------------------------------------------------------------------------------!
         do ifpy = 1,dset_npolygons_global
            currad = dist_gc(file_lons(ifpy),cgrid%lon(ipy),file_lats(ifpy),cgrid%lat(ipy))
            if (currad <  minrad) then
               py_index = ifpy
               minrad   = currad
            end if
         end do

         !---------------------------------------------------------------------------------!
         !      A suitably close polygon has been located in the datasets.  Use these      !
         ! values, and its children values in sites, patchs and cohorts.                   !
         !---------------------------------------------------------------------------------!
         iparallel = 0
         
         !---------------------------------------------------------------------------------!
         !      POLYGON level variables.                                                   !
         !---------------------------------------------------------------------------------!
         !----- Load 1D dataset. ----------------------------------------------------------!
         dsetrank = 1
         globdims(1) = int(dset_npolygons_global,8)
         chnkdims(1) = 1_8
         chnkoffs(1) = int(py_index - 1,8)
         memdims(1)  = 1_8
         memoffs(1)  = 0_8
         memsize(1)  = 1_8
         !---- The ipy:ipy notation is needed for ifort when checking interfaces. ---------!
         call hdf_getslab_i(cgrid%load_adjacency(ipy:ipy),'LOAD_ADJACENCY '                &
                           ,dsetrank,iparallel,.true.,foundvar)
         call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.,foundvar)
         
         !----- Load the workload (2D). ---------------------------------------------------!
         dsetrank    = 2
         globdims(1) = int(13,8)
         chnkdims(1) = int(13,8)
         memdims(1)  = int(13,8)
         memsize(1)  = int(13,8)
         chnkoffs(1) = 0_8
         memoffs(1)  = 0_8
         globdims(2) = int(dset_npolygons_global,8)
         chnkdims(2) = 1_8
         chnkoffs(2) = int(py_index - 1,8)
         memdims(2)  = 1_8
         memsize(2)  = 1_8
         memoffs(2)  = 0_8
         call hdf_getslab_r(cgrid%workload(:,ipy),'WORKLOAD ',dsetrank,iparallel,.false.   &
                           ,foundvar)
         !---------------------------------------------------------------------------------!



         !----- Load the lakesite data-----------------------------------------------------!
         allocate(islakesite(pysi_n(py_index)))
         islakesite  = 0
         dsetrank    = 1_8
         globdims(1) = int(dset_nsites_global,8)
         chnkdims(1) = int(pysi_n(py_index),8)
         chnkoffs(1) = int(pysi_id(py_index) - 1,8)
         memdims(1)  = int(pysi_n(py_index),8)
         memsize(1)  = int(pysi_n(py_index),8)
         memoffs(1)  = 0_8  
         call hdf_getslab_i(islakesite,'ISLAKESITE ',dsetrank,iparallel,.false.,foundvar)

         ndry_sites = int(pysi_n(py_index))-sum(islakesite)
         !---------------------------------------------------------------------------------!


         !----- Allocate the vector of sites in the polygon. ------------------------------!
         call allocate_polygontype(cpoly,ndry_sites)
         call soil_default_fill(cgrid,igr,ipy)
         !---------------------------------------------------------------------------------!

         is = 0
         siteloop: do isi=1,pysi_n(py_index)

            if (islakesite(isi) == 0) then
               is=is+1

               !----- Reset the HDF5 auxiliary variables before moving to the next level. -!
               globdims = 0_8
               chnkdims = 0_8
               chnkoffs = 0_8
               memoffs  = 0_8
               memdims  = 0_8
               memsize  = 1_8
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      SITE level variables.                                                !
               !---------------------------------------------------------------------------!
               !----- Load 1D dataset. ----------------------------------------------------!
               dsetrank    = 1_8
               globdims(1) = int(dset_nsites_global,8)
               chnkdims(1) = int(1,8)
               chnkoffs(1) = int(pysi_id(py_index) - 2 + isi,8)
               memdims(1)  = int(1,8)
               memsize(1)  = int(1,8)
               memoffs(1)  = 0_8

               call hdf_getslab_r( cpoly%area       (is:is), 'AREA_SI '                    &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_r( cpoly%moist_f    (is:is), 'MOIST_F '                    &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_r( cpoly%moist_W    (is:is), 'MOIST_W '                    &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_r( cpoly%elevation  (is:is), 'ELEVATION '                  &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_r( cpoly%slope      (is:is), 'SLOPE '                      &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_r( cpoly%aspect     (is:is), 'ASPECT '                     &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_r( cpoly%TCI        (is:is), 'TCI '                        &
                                 , dsetrank, iparallel, .true.,foundvar)
               call hdf_getslab_i( cpoly%patch_count(is:is), 'PATCH_COUNT '                &
                                 , dsetrank, iparallel, .true.,foundvar)

               !----- Load 2D dataset. ----------------------------------------------------!
               dsetrank     = 2_8
               globdims(1)  = int(dset_nzg,8)     ! How many layers in the dataset?
               chnkdims(1)  = int(1,8)            ! We are only extracting one layer
               memdims(1)   = int(1,8)            ! We only need memory for one layer
               memsize(1)   = int(1,8)            ! On both sides
               chnkoffs(1)  = int(dset_nzg - 1,8) ! Take the top layer, not the bottom
               memoffs(1)   = 0_8
   
               globdims(2)  = int(dset_nsites_global,8)
               chnkdims(2)  = int(1,8)
               chnkoffs(2)  = int(pysi_id(py_index) - 2 + isi,8)
               memdims(2)   = int(1,8)
               memsize(2)   = int(1,8)
               memoffs(2)   = 0_8
               call hdf_getslab_i(cpoly%ntext_soil(nzg:nzg,is:is),'NTEXT_SOIL '            &
                                 ,dsetrank,iparallel,.true.,foundvar)

               !----- Now fill the soil column based on the top layer data. ---------------!
               do k=1,nzg-1
                  cpoly%ntext_soil(k,is) = cpoly%ntext_soil(nzg,is)
               end do

               csite => cpoly%site(isi)

               if (sipa_n(si_index) > 0) then

                  !----- Fill 1D polygon (site unique) level variables. -------------------!
                  call allocate_sitetype(csite,sipa_n(si_index))

                  !----- Reset the HDF5 auxiliary variables before moving to the next level. -!
                  globdims = 0_8
                  chnkdims = 0_8
                  chnkoffs = 0_8
                  memoffs  = 0_8
                  memdims  = 0_8
                  memsize  = 1_8
                  !------------------------------------------------------------------------!

                  iparallel = 0
                  
                  dsetrank = 1
                  globdims(1) = int(dset_npatches_global,8)
                  chnkdims(1) = int(csite%npatches,8)
                  chnkoffs(1) = int(sipa_id(si_index) - 1,8)
                  memdims(1)  = int(csite%npatches,8)
                  memsize(1)  = int(csite%npatches,8)
                  memoffs(1)  = 0

                  call hdf_getslab_i(csite%dist_type ,'DIST_TYPE '                         &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%age       ,'AGE '                               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%area      ,'AREA '                              &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%sum_dgd   ,'SUM_DGD '                           &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%sum_chd   ,'SUM_CHD '                           &
                                    ,dsetrank,iparallel,.true.,foundvar)

                  call hdf_getslab_r(csite%fast_soil_C       ,'FAST_SOIL_C '               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%slow_soil_C       ,'SLOW_SOIL_C '               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%fast_soil_N       ,'FAST_SOIL_N '               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%structural_soil_C ,'STRUCTURAL_SOIL_C '         &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%structural_soil_L ,'STRUCTURAL_SOIL_L '         &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N '        &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Check whether the history file is new or old.  We determine this   !
                  ! by searching for variable plantation.  In case this variable isn't     !
                  ! present, then it must be the new history.   Otherwise we correct the   !
                  ! indices for secondary forests.                                         !
                  !------------------------------------------------------------------------!
                  allocate (plantation(csite%npatches))
                  plantation(:) = 0
                  call hdf_getslab_i(plantation ,'PLANTATION '                             &
                                    ,dsetrank,iparallel,.false.,foundvar)
                  if (foundvar) then
                     do ipa=1,csite%npatches
                        select case(csite%dist_type(ipa))
                        case (2)
                           !---------------------------------------------------------------!
                           !     Secondary forests are now divided in three categories.    !
                           !---------------------------------------------------------------!
                           if (plantation(ipa) == 0 .and. ianth_disturb == 0) then
                              csite%dist_type(ipa) = 5
                           else if (plantation(ipa) == 0 .and. ianth_disturb == 1) then
                              csite%dist_type(ipa) = 6
                           else
                              csite%dist_type(ipa) = 2
                           end if
                           !---------------------------------------------------------------!
                        end select
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end if
                  deallocate(plantation)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Loop over all sites and fill the patch-level variables.            !
                  !------------------------------------------------------------------------!
                  patchloop: do ipa = 1,csite%npatches

                     !---------------------------------------------------------------------!
                     !     Reset the HDF5 auxiliary variables before moving to the next    !
                     ! level.                                                              !
                     !---------------------------------------------------------------------!
                     globdims = 0_8
                     chnkdims = 0_8
                     chnkoffs = 0_8
                     memoffs  = 0_8
                     memdims  = 0_8
                     memsize  = 1_8
                     !---------------------------------------------------------------------!

                     cpatch => csite%patch(ipa)

                     !---------------------------------------------------------------------!
                     !    Initialise patch-level variables that depend on the cohort ones. !
                     !---------------------------------------------------------------------!
                     csite%plant_ag_biomass(ipa)  = 0.0

                     pa_index = sipa_id(si_index) + ipa - 1
                     call allocate_patchtype(cpatch,paco_n(pa_index))

                     !---------------------------------------------------------------------!
                     !     Empty patches may exist, so make sure that this part is called  !
                     ! only when there are cohorts.                                        !
                     !---------------------------------------------------------------------!
                     if (cpatch%ncohorts > 0) then
                        !----- First the 1-D variables. -----------------------------------!
                        dsetrank = 1
                        globdims(1) = int(dset_ncohorts_global,8)
                        chnkdims(1) = int(cpatch%ncohorts,8)
                        chnkoffs(1) = int(paco_id(pa_index) - 1,8)
                        memdims(1)  = int(cpatch%ncohorts,8)
                        memsize(1)  = int(cpatch%ncohorts,8)
                        memoffs(1)  = 0_8

                        call hdf_getslab_r(cpatch%dbh   ,'DBH '                            &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_r(cpatch%bdead ,'BDEAD '                          &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_i(cpatch%pft   ,'PFT '                            &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_r(cpatch%nplant,'NPLANT '                         &
                                          ,dsetrank,iparallel,.true.,foundvar)

                        !------------------------------------------------------------------!
                        !    Find derived properties from Bdead.  In the unlikely case     !
                        ! that bdead is zero, then we use DBH as the starting point.  In   !
                        ! both cases we assume that plants are in allometry.               !
                        !------------------------------------------------------------------!
                        do ico=1,cpatch%ncohorts
                           ipft = cpatch%pft(ico)

                           if (igrass == 1 .and. is_grass(ipft)                            &
                                           .and. cpatch%bdead(ico)>0.0) then
                              !-- if the initial file was running with igrass = 0, bdead   !
                              ! should be nonzero.  If the new run has igrass = 1, bdead   !
                              ! is set to zero and that biomass is discarded               !
                              cpatch%dbh(ico)   = max(cpatch%dbh(ico),min_dbh(ipft))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                              cpatch%bdead(ico) = 0.0
                              
                           else if (cpatch%bdead(ico) > 0.0 .and. igrass == 0) then
                              ! grasses have bdead in both input and current run (igrass=0)
                              cpatch%bdead(ico) = max(cpatch%bdead(ico),min_bdead(ipft))
                              cpatch%dbh(ico)   = bd2dbh(ipft,cpatch%bdead(ico))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                           else 
                              ! it is either a new grass (igrass=1) in the initial file,   !
                              ! or the value for bdead is missing from the files           !
                              cpatch%dbh(ico)   = max(cpatch%dbh(ico),min_dbh(ipft))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                              cpatch%bdead(ico) = dbh2bd(cpatch%dbh  (ico),ipft)
                           end if


                           cpatch%bleaf(ico)  = size2bl(cpatch%dbh(ico),cpatch%hite(ico)   &
                                                        ,cpatch%pft(ico))

                           !----- Find the other pools. -----------------------------------!
                           salloc  = (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))
                           salloci = 1.0 / salloc
                           cpatch%balive  (ico) = cpatch%bleaf(ico) * salloc
                           cpatch%broot   (ico) = cpatch%balive(ico) * q(ipft) * salloci
                           cpatch%bsapwooda(ico) = cpatch%balive(ico) * qsw(ipft)          &
                                                 * cpatch%hite(ico) * salloci * agf_bs(ipft)
                           cpatch%bsapwoodb(ico) = cpatch%balive(ico) * qsw(ipft)          &
                                                 * cpatch%hite(ico) * salloci              &
                                                 * (1.-agf_bs(ipft))
                           cpatch%bstorage(ico) = 0.0
                           cpatch%phenology_status(ico) = 0
                           
                           
                        end do


                        !------------------------------------------------------------------!
                        !     Carbon balance variables.                                    !
                        ! MLO.  I commented this out because a restart is likely to have   !
                        !       different settings and different environment, so it's      !
                        !       likely that the system will reach a different equilibrium. !
                        !       For simplicity, we just use the typical initialisation in  !
                        !       which we give plants one year to adjust to the new         !
                        !       conditions.                                                !
                        !------------------------------------------------------------------!
                        ! dsetrank    = 2
                        ! globdims(1) = 13_8
                        ! chnkdims(1) = 13_8
                        ! chnkoffs(1) = 0_8
                        ! memdims(1)  = 13_8
                        ! memsize(1)  = 13_8
                        ! memoffs(2)  = 0_8
                        ! globdims(2) = int(dset_ncohorts_global,8)
                        ! chnkdims(2) = int(cpatch%ncohorts,8)
                        ! chnkoffs(2) = int(paco_id(pa_index) - 1,8)
                        ! memdims(2)  = int(cpatch%ncohorts,8)
                        ! memsize(2)  = int(cpatch%ncohorts,8)
                        ! memoffs(2)  = 0_8

                        ! call hdf_getslab_r(cpatch%cb    ,'CB '                           &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        ! call hdf_getslab_r(cpatch%cb_lightmax,'CB_LIGHTMAX '             &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        ! call hdf_getslab_r(cpatch%cb_moistmax,'CB_MOISTMAX '             &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        do ico = 1, cpatch%ncohorts
                           cpatch%cb          (1:12,ico) = 1.0
                           cpatch%cb_lightmax (1:12,ico) = 1.0
                           cpatch%cb_moistmax (1:12,ico) = 1.0
                           cpatch%cb          (  13,ico) = 0.0
                           cpatch%cb_lightmax (  13,ico) = 0.0
                           cpatch%cb_moistmax (  13,ico) = 0.0
                        end do
                        !------------------------------------------------------------------!



                        cohortloop: do ico=1,cpatch%ncohorts
                           !---------------------------------------------------------------!
                           !    We will now check the PFT of each cohort, so we determine  !
                           ! if this is a valid PFT.  If not, then we must decide what we  !
                           ! should do...                                                  !
                           !---------------------------------------------------------------!
                           if (.not. include_pft(cpatch%pft(ico))) then
                              select case(pft_1st_check)
                              case (0)
                                 !----- Stop the run. -------------------------------------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,' and it is not in your include_these_pft...'
                                 call fatal_error('Invalid PFT in history file'            &
                                                 ,'read_ed21_history_file'                 &
                                                 ,'ed_read_ed21_history.F90')

                              case (1)
                                 !----- Include the unexpected PFT in the list. -----------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,'... Including this PFT in your include_these_pft...'
                                 include_pft(cpatch%pft(ico))          = .true.
                                 include_these_pft(count(include_pft)) = cpatch%pft(ico)

                                 call sort_up(include_these_pft,n_pft)

                                 if (is_grass(cpatch%pft(ico))) then
                                    include_pft_ag(cpatch%pft(ico)) = .true.
                                 end if

                              case (2)
                                 !----- Ignore the unexpect PFT. --------------------------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,'... Ignoring it...'
                                 !---------------------------------------------------------!
                                 !    The way we will ignore this cohort is by setting its !
                                 ! nplant to zero, and calling the "terminate_cohorts"     !
                                 ! subroutine right after this.                            !
                                 !---------------------------------------------------------!
                                 cpatch%nplant(ico) = 0.
                              end select
                           end if

                           !---------------------------------------------------------------!
                           !     Make sure that the biomass won't lead to FPE.  This       !
                           ! should never happen when using a stable ED-2.1 version, but   !
                           ! older versions had "zombie" cohorts.  Here we ensure that     !
                           ! the model initialises with stable numbers whilst ensuring     !
                           ! that the cohorts will be eliminated.                          !
                           !---------------------------------------------------------------!
                           if (cpatch%balive(ico) > 0.            .and.                    &
                               cpatch%balive(ico) < tiny_biomass) then
                              cpatch%balive(ico) = tiny_biomass
                           end if
                           if (cpatch%bleaf(ico) > 0.                .and.                 &
                               cpatch%bleaf(ico) < tiny_biomass)     then
                              cpatch%bleaf(ico) = tiny_biomass
                           end if
                           if (cpatch%broot(ico) > 0.                .and.                 &
                               cpatch%broot(ico) < tiny_biomass)     then
                              cpatch%broot(ico) = tiny_biomass
                           end if
                           if (cpatch%bsapwooda(ico) > 0.            .and.                 &
                               cpatch%bsapwooda(ico) < tiny_biomass) then
                              cpatch%bsapwooda(ico) = tiny_biomass
                           end if
                           if (cpatch%bsapwoodb(ico) > 0.            .and.                 &
                               cpatch%bsapwoodb(ico) < tiny_biomass) then
                              cpatch%bsapwoodb(ico) = tiny_biomass
                           end if
                           if (cpatch%bdead(ico) > 0.                .and.                 &
                               cpatch%bdead(ico) < tiny_biomass)     then
                              cpatch%bdead(ico) = tiny_biomass
                           end if
                           if (cpatch%bstorage(ico) > 0.             .and.                 &
                               cpatch%bstorage(ico) < tiny_biomass)  then
                              cpatch%bstorage(ico) = tiny_biomass
                           end if
                           !---------------------------------------------------------------!




                           !----- Compute the above-ground biomass. -----------------------!
                           cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)&
                                                    ,cpatch%bsapwooda(ico),cpatch%pft(ico))
                           cpatch%basarea(ico)  = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)

                            
                           !----- Assign LAI, WAI, and CAI --------------------------------!
                           call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)          &
                                            ,cpatch%bdead(ico),cpatch%balive(ico)          &
                                            ,cpatch%dbh(ico), cpatch%hite(ico)             &
                                            ,cpatch%pft(ico), SLA(cpatch%pft(ico))         &
                                            ,cpatch%lai(ico),cpatch%wai(ico)               &
                                            ,cpatch%crown_area(ico),cpatch%bsapwooda(ico))


                           !----- Update the derived patch-level variables. ---------------!
                           csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)       &
                                                       + cpatch%agb(ico)*cpatch%nplant(ico)

                           !----- Initialise the other cohort level variables. ------------!
                           call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                        end do cohortloop

                        !------------------------------------------------------------------!
                        !    Eliminate any "unwanted" cohort (i.e., those which nplant was !
                        ! set to zero so it would be removed).                             !
                        !------------------------------------------------------------------!
                        call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)

                     end if
                  end do patchloop
               else
                  !----- This should never happen, but, just in case... -------------------!
                  call fatal_error('A site with no patches was found...'                   &
                                  ,'read_ed21_history_file','ed_read_ed21_history.F90')
               end if

               !----- Initialise the other patch-level variables. -------------------------!
               call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
            end if
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!

         deallocate(islakesite)

         !----- Initialise some site-level variables. -------------------------------------!
         call init_ed_site_vars(cpoly,cgrid%lat(ipy))
         !---------------------------------------------------------------------------------!

      end do polyloop


      !----- Initialise the other polygon-level variables. --------------------------------!
      call init_ed_poly_vars(cgrid)


      !----- Close the dataset. -----------------------------------------------------------!
      call h5fclose_f(file_id, hdferr)
      if (hdferr /= 0) then
         write (unit=*,fmt='(a,1x,a)') 'File: ',trim(hnamel)
         write (unit=*,fmt='(a)'     ) 'Problem: Failed closing the HDF5 dataset.'
         write (unit=*,fmt='(a,1x,i4)') 'HDFerr: ',hdferr
         call fatal_error('Could not close the HDF file'                                   &
                         ,'read_ed21_history_file','ed_read_ed21_history.F90')
      end if

      !------ Deallocate the temporary vectors, so no memory leak happens. ----------------!
      deallocate(file_lats)
      deallocate(file_lons)
      deallocate(paco_n   )
      deallocate(paco_id  )
      deallocate(sipa_n   )
      deallocate(sipa_id  )
      deallocate(pysi_n   )
      deallocate(pysi_id  )
      
   end do gridloop
   
   !----- Reset auto error printing to "on". ----------------------------------------------!
   call h5eset_auto_f(1,hdferr)

   !----- Close the file. -----------------------------------------------------------------!
   call h5close_f(hdferr)

#else
   call fatal_error ('You cannot restart with ED-2.1 without using HDF5...'                &
                    ,'read_ed21_history_file','ed_read_ed21_history.F90')
#endif   

   return
end subroutine read_ed21_history_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will read the ED-2.1 restart files.  The difference between this     !
! subroutine and the previous one is that here the user may provide a collection of files  !
! as the initial condition.  This may be either one regional run, one single-polygon run,  !
! a collection of single-polygon runs, and a mix of all the previous options, we must make !
! sure that every option is accounted, and that the polygon is always initialised with the !
! closest restart polygon available.                                                       !
!------------------------------------------------------------------------------------------!
subroutine read_ed21_history_unstruct
#if USE_HDF5
   use hdf5
#endif
   use ed_max_dims    , only : n_pft                   & ! intent(in)
                             , huge_polygon            & ! intent(in)
                             , str_len                 & ! intent(in)
                             , n_dist_types            & ! intent(in)
                             , maxfiles                & ! intent(in)
                             , maxlist                 ! ! intent(in)
   use pft_coms       , only : SLA                     & ! intent(in)
                             , q                       & ! intent(in)
                             , qsw                     & ! intent(in)
                             , hgt_min                 & ! intent(in)
                             , min_dbh                 & ! intent(in)
                             , min_bdead               & ! intent(in)
                             , is_grass                & ! intent(in)
                             , include_pft             & ! intent(in)
                             , include_pft_ag          & ! intent(in)
                             , pft_1st_check           & ! intent(in)
                             , include_these_pft       & ! intent(in)
                             , agf_bs                  & ! intent(in)
                             , min_cohort_size         ! ! intent(in)
   use ed_misc_coms   , only : sfilin                  & ! intent(in)
                             , current_time            & ! intent(in)
                             , imonthh                 & ! intent(in)
                             , iyearh                  & ! intent(in)
                             , idateh                  & ! intent(in)
                             , itimeh                  & ! intent(in)
                             , ied_init_mode           & ! intent(in)
                             , igrass                  & ! intent(in)
                             , max_poi99_dist          ! ! intent(in)
   use ed_state_vars  , only : polygontype             & ! variable type
                             , sitetype                & ! variable type
                             , patchtype               & ! variable type
                             , edtype                  & ! variable type
                             , edgrid_g                & ! subroutine
                             , allocate_sitetype       & ! subroutine
                             , allocate_patchtype      ! ! subroutine
   use grid_coms      , only : ngrids                  & ! intent(in)
                             , nzg                     ! ! intent(in)
   use consts_coms    , only : pio4                    ! ! intent(in)
   use hdf5_coms      , only : file_id                 & ! intent(in)
                             , dset_id                 & ! intent(in)
                             , dspace_id               & ! intent(in)
                             , plist_id                & ! intent(in)
                             , globdims                & ! intent(in)
                             , chnkdims                & ! intent(in)
                             , chnkoffs                & ! intent(in)
                             , memdims                 & ! intent(in)
                             , memoffs                 & ! intent(in)
                             , memsize                 ! ! intent(in)
   use allometry      , only : area_indices            & ! function
                             , ed_biomass              & ! function
                             , bd2dbh                  & ! function
                             , dbh2h                   & ! function
                             , dbh2bd                  & ! function
                             , size2bl                 ! ! function
   use fuse_fiss_utils, only : terminate_cohorts       ! ! subroutine
   use disturb_coms   , only : ianth_disturb           & ! intent(in)
                             , lu_rescale_file         & ! intent(in)
                             , min_patch_area          ! ! intent(in)
   use soil_coms      , only : soil                    ! ! intent(in)

   implicit none

#if (USE_HDF5)
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)          , pointer                              :: cgrid
   type(polygontype)     , pointer                              :: cpoly
   type(sitetype)        , pointer                              :: csite
   type(patchtype)       , pointer                              :: cpatch
   character(len=str_len), dimension(maxlist)                   :: full_list
   character(len=str_len), dimension(maxfiles)                  :: histo_list
   character(len=1)                                             :: vnam
   character(len=3)                                             :: cgr
   character(len=str_len)                                       :: hnamel
   integer               , dimension(maxfiles)                  :: ngridpoly
   integer               , dimension(huge_polygon)              :: pyfile_list
   integer               , dimension(huge_polygon)              :: pyindx_list
   integer               , dimension(  :)         , allocatable :: plantation
   integer               , dimension(  :)         , allocatable :: pclosest
   integer               , dimension(  :)         , allocatable :: psrcfile
   integer               , dimension(  :)         , allocatable :: pysi_n
   integer               , dimension(  :)         , allocatable :: pysi_id
   integer               , dimension(  :)         , allocatable :: sipa_n
   integer               , dimension(  :)         , allocatable :: sipa_id
   integer               , dimension(  :)         , allocatable :: paco_n
   integer               , dimension(  :)         , allocatable :: paco_id
   integer               , dimension(:,:)         , allocatable :: this_ntext
   integer               , dimension(  :)         , allocatable :: islakesite
   real                  , dimension(  :)         , allocatable :: tpoly_area
   real                  , dimension(  :)         , allocatable :: tpoly_moist_f
   real                  , dimension(  :)         , allocatable :: tpoly_moist_w
   real                  , dimension(  :)         , allocatable :: tpoly_elevation
   real                  , dimension(  :)         , allocatable :: tpoly_slope
   real                  , dimension(  :)         , allocatable :: tpoly_aspect
   real                  , dimension(  :)         , allocatable :: tpoly_TCI
   integer               , dimension(  :)         , allocatable :: tpoly_patch_count
   integer               , dimension(  :)         , allocatable :: tpoly_lsl
   integer               , dimension(:,:)         , allocatable :: tpoly_ntext_soil
   integer                                                      :: year
   integer                                                      :: igr
   integer                                                      :: ipy
   integer                                                      :: isi
   integer                                                      :: is
   integer                                                      :: ipa
   integer                                                      :: ico
   integer                                                      :: isi_best
   integer                                                      :: is_best
   integer                                                      :: isi_try
   integer                                                      :: is_try
   integer                                                      :: nsoil
   integer                                                      :: nsoil_try
   integer                                                      :: nsites_inp
   integer                                                      :: xclosest
   integer                                                      :: nflist
   integer                                                      :: nhisto
   integer                                                      :: nrescale
   integer                                                      :: k
   integer                                                      :: nf
   integer                                                      :: dset_npolygons_global
   integer                                                      :: dset_nsites_global
   integer                                                      :: dset_npatches_global
   integer                                                      :: dset_ncohorts_global
   integer                                                      :: dset_nzg
   integer                                                      :: ngr
   integer                                                      :: ifpy
   integer                                                      :: ipft
   integer                                                      :: ipya
   integer                                                      :: ipyz
   integer                                                      :: ierr
   integer                                                      :: ilu
   integer                                                      :: py_index
   integer                                                      :: si_index
   integer                                                      :: pa_index
   integer                                                      :: dsetrank,iparallel
   integer                                                      :: hdferr
   integer                                                      :: total_grid_py
   integer                                                      :: poi_minloc
   integer                                                      :: ngp1
   integer                                                      :: ndry_sites
   logical                                                      :: exists
   logical                                                      :: rescale_glob
   logical                                                      :: rescale_loc
   real                  , dimension(  :)         , allocatable :: pdist
   real                  , dimension(huge_polygon)              :: plon_list
   real                  , dimension(huge_polygon)              :: plat_list
   real                  , dimension(huge_polygon)              :: dist_rscl
   real                  , dimension(huge_polygon)              :: wlon_rscl
   real                  , dimension(huge_polygon)              :: clon_rscl
   real                  , dimension(huge_polygon)              :: elon_rscl
   real                  , dimension(huge_polygon)              :: slat_rscl
   real                  , dimension(huge_polygon)              :: clat_rscl
   real                  , dimension(huge_polygon)              :: nlat_rscl
   real                  , dimension(n_dist_types,huge_polygon) :: newarea
   real                  , dimension(n_dist_types)              :: oldarea
   real                                                         :: textdist_try
   real                                                         :: textdist_min
   real                                                         :: dummy
   real                                                         :: elim_nplant
   real                                                         :: elim_lai
   real                                                         :: salloc
   real                                                         :: salloci
   logical                                                      :: foundvar
   !----- Local constants. ----------------------------------------------------------------!
   real                                           , parameter   :: tiny_biomass = 1.e-20
   !----- External functions. -------------------------------------------------------------!
   real                                           , external    :: dist_gc 
   !---------------------------------------------------------------------------------------!



   !----- Open the HDF environment. -------------------------------------------------------!
   call h5open_f(hdferr)

   !----- Turn off automatic error printing. ----------------------------------------------!
   call h5eset_auto_f(0,hdferr)


   !---------------------------------------------------------------------------------------!
   !     Big loop over all grids.                                                          !
   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids


      !----- Retrieve all files with the specified prefix. --------------------------------!
      call ed_filelist(full_list,sfilin(igr),nflist)
      !----- Check every file and save only those that are actually history files. --------!
      call ed21_fileinfo(nflist,full_list,nhisto,histo_list)

      !----- Initialize the dimensional control variables for the H5 slabs. ---------------!
      globdims = 0_8
      chnkdims = 0_8
      chnkoffs = 0_8
      memoffs  = 0_8
      memdims  = 0_8
      memsize  = 1_8  

      !------------------------------------------------------------------------------------!
      !     First thing, we go through every file, open, and retrieve only some polygon-   !
      ! level information.  This will be used for mapping the files later.                 !
      !------------------------------------------------------------------------------------!
      ngridpoly  (:) = 0
      pyfile_list(:) = 0
      pyindx_list(:) = 0
      ipyz           = 0
      lonlatloop: do nf=1,nhisto
         hnamel = histo_list(nf)

         !----- Open the HDF5 file. -------------------------------------------------------!
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a,1x,i4)') ' - Error opening HDF5 file - error - ',hdferr
            write(unit=*,fmt='(a,1x,a)') ' - File name: ',trim(hnamel)
            call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)            &
                            ,'read_ed21_history_file','ed_read_ed21_history.F90')
         end if

         !----- Retrieve the number of polygons in this file. -----------------------------!
         call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
         
         !----- Determine the global position of these polygons. --------------------------!
         ipya = ipyz + 1
         ipyz = ipyz + dset_npolygons_global

         !----- In case we are meshing grids ----------------------------------------------!
         ngridpoly(nf) = dset_npolygons_global

         !------ Here we save the file where we can find these polygons. ------------------!
         do ipy=ipya,ipyz
            pyfile_list(ipy) = nf
            pyindx_list(ipy) = ipy-ipya+1
         end do

         !---------------------------------------------------------------------------------!
         !      Retrieve the polygon coordinates data.                                     !
         !---------------------------------------------------------------------------------!
         globdims(1) = int(dset_npolygons_global,8)

         call h5dopen_f(file_id,'LONGITUDE', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_REAL,plon_list(ipya:ipyz),globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
         
         call h5dopen_f(file_id,'LATITUDE', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_REAL,plat_list(ipya:ipyz),globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
         
         call h5fclose_f(file_id, hdferr)
         if (hdferr /= 0) then
            write (unit=*,fmt='(a,1x,a)') 'File: ',trim(hnamel)
            write (unit=*,fmt='(a)'     ) 'Problem: Failed closing the HDF5 dataset.'
            write (unit=*,fmt='(a,1x,i4)') 'HDFerr: ',hdferr
            call fatal_error('Could not close the HDF file'                                &
                            ,'read_ed21_history_unstruct','ed_read_ed21_history.F90')
         end if
      end do lonlatloop

      total_grid_py = ipyz
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If this is a simulation with anthropogenic disturbance, check and read the re- !
      ! scale file if it exists.                                                           !
      !------------------------------------------------------------------------------------!
      wlon_rscl(:) =  190.
      elon_rscl(:) = -190.
      slat_rscl(:) =  100.
      nlat_rscl(:) = -100.
      inquire(file=trim(lu_rescale_file(igr)),exist=exists)
      rescale_glob = ianth_disturb == 1 .and. exists
      nrescale = 0
      if (rescale_glob) then
         open (unit=13,file=trim(lu_rescale_file(igr)),status='old',action='read')
         read (unit=13,fmt=*)
         readrescale: do
            nrescale = nrescale + 1
            read (unit=13,fmt=*,iostat=ierr)  wlon_rscl(nrescale),elon_rscl(nrescale)      &
                                             ,slat_rscl(nrescale),nlat_rscl(nrescale)      &
                                             ,clon_rscl(nrescale),clat_rscl(nrescale)      &
                                             ,dummy,(newarea(ilu,nrescale),ilu=1           &
                                             ,n_dist_types)

            if (ierr /= 0) then
               nrescale = nrescale - 1
               exit readrescale
            end if
         end do readrescale
         rescale_glob = nrescale > 0
         close (unit=13,status='keep')
      elseif (ianth_disturb == 1 .and. len_trim(lu_rescale_file(igr)) > 0) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') '  In subroutine read_ed21_history_unstruct:'
         write (unit=*,fmt='(a,1x,i5)') '  - Grid: ',igr
         write (unit=*,fmt='(a)') '  In subroutine read_ed21_history_unstruct:'
         write (unit=*,fmt='(a)') '  - File '//trim(lu_rescale_file(igr))//                &
                                     ' wasn''t found...'
         write (unit=*,fmt='(a)') '  - Assuming no rescaling...'
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
      end if
      !------------------------------------------------------------------------------------!





      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !     Now, we will go through all polygons, and we will determine which input        !
      ! polygon is the closest one to each polygon.                                        !
      !------------------------------------------------------------------------------------!
      allocate (pclosest(cgrid%npolygons),psrcfile(cgrid%npolygons))
      allocate(pdist(total_grid_py))
      nneighloop: do ipy=1,cgrid%npolygons
         

         !----- Reset pdist to a very large number. ---------------------------------------!
         
         pdist(1:total_grid_py)   = 1.e20

         do ifpy=1,total_grid_py
            pdist(ifpy) = dist_gc(plon_list(ifpy),cgrid%lon(ipy)                           &
                                 ,plat_list(ifpy),cgrid%lat(ipy))
         end do

         select case (ied_init_mode)
         case (5)
            pclosest(ipy) = pyindx_list(minloc(pdist,dim=1))
            psrcfile(ipy) = pyfile_list(minloc(pdist,dim=1))
         case (99)
            !------------------------------------------------------------------------------!
            !      Alternative method for mixing 1 grid and POI's.  Only use the grid if   !
            ! there is NOT a POI  within a user specified resolution.  Remember, this      !
            ! assumes there is only one gridded file, and it is the first file.            !
            !------------------------------------------------------------------------------!
            ngp1          = ngridpoly(1)
            pclosest(ipy) = pyindx_list(minloc(pdist(1:ngp1),dim=1))
            psrcfile(ipy) = pyfile_list(1)
            
            poi_minloc    = minloc(pdist(ngp1+1:total_grid_py),dim=1) + ngp1
            
            if( pdist(poi_minloc) < max_poi99_dist ) then
               pclosest(ipy)  = pyindx_list(poi_minloc)
               psrcfile(ipy)  = pyfile_list(poi_minloc)
            end if
         end select
      end do nneighloop
      
      deallocate(pdist)
      !------------------------------------------------------------------------------------!
      !     Now that we have all polygons matched with their nearest neighbours, we will   !
      ! loop over all files instead of all polygons.  This is to avoid opening and closing !
      ! the files too many times.                                                          !
      !------------------------------------------------------------------------------------!
      rstfileloop: do nf=1, nhisto

         !---------------------------------------------------------------------------------!
         !    Before anything else, we check whether this file needs to be opened, i.e.,   !
         ! whether any polygon has its nearest neighbour in this file.  If not, we move    !
         ! to the next one.                                                                !
         !---------------------------------------------------------------------------------!
         if (.not. (any(psrcfile == nf))) cycle rstfileloop

         !---------------------------------------------------------------------------------!
         !     If we are here, at least one polygon closest source is in this file, open   !
         ! it and load all fill in all polygons whose nearest neighbour is in this file.   !
         !---------------------------------------------------------------------------------!
         hnamel = histo_list(nf)

         !----- Open file. ----------------------------------------------------------------!
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a,1x,i4)') ' - Error opening HDF5 file - error - ',hdferr
            write(unit=*,fmt='(a,1x,a)') ' - File name: ',trim(hnamel)
            call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)            &
                            ,'read_ed21_history_file','ed_read_ed21_history.F90')
         end if

         !---------------------------------------------------------------------------------!
         !      Retrieve global vector sizes and mapping tree.                             !
         !---------------------------------------------------------------------------------!
         globdims = 0_8
         chnkdims = 0_8
         chnkoffs = 0_8

         globdims(1) = 1_8

         call h5dopen_f(file_id,'NZG', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nzg,globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nsites_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npatches_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_ncohorts_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         globdims(1) = int(dset_npolygons_global,8)

         allocate(pysi_n(dset_npolygons_global))
         allocate(pysi_id(dset_npolygons_global))

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

         globdims(1) = int(dset_nsites_global,8)

         allocate(sipa_n(dset_nsites_global))
         allocate(sipa_id(dset_nsites_global))

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

         globdims(1) = int(dset_npatches_global,8)
         allocate(paco_n(dset_npatches_global))
         allocate(paco_id(dset_npatches_global))

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

         polyloop: do ipy = 1,cgrid%npolygons
            cpoly => cgrid%polygon(ipy)

            !----- We skip the polygon if its source polygon is not in this file. ---------!
            if (psrcfile(ipy) /= nf) cycle polyloop

            !------------------------------------------------------------------------------!
            !    Use the index corresponding to the relative position of the input polygon !
            ! in the source file to which the polygon belongs.  Use these values, and its  !
            ! children values in sites, patchs and cohorts.                                !
            !------------------------------------------------------------------------------!
            py_index = pclosest(ipy)

            !------------------------------------------------------------------------------!
            !      Retrieve the polygon coordinates data.                                  !
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     In case we seek to rescale, we must first check whether a scale for the  !
            ! current polygon.                                                             !
            !------------------------------------------------------------------------------!
            !----- Initialise distance and co-ordinates to non-sense numbers. -------------!
            dist_rscl(:) = 1.e+20 ! Initialise to a large distance and non-sense 
            if (rescale_glob) then
               neighbour: do k=1,nrescale
                  dist_rscl(k) = dist_gc(clon_rscl(k),cgrid%lon(ipy)                       &
                                        ,clat_rscl(k),cgrid%lat(ipy))
               end do neighbour
               xclosest = minloc(dist_rscl,dim=1)
               
               rescale_loc = cgrid%lon(ipy) > wlon_rscl(xclosest) .and.                    &
                             cgrid%lon(ipy) < elon_rscl(xclosest) .and.                    &
                             cgrid%lat(ipy) > slat_rscl(xclosest) .and.                    &
                             cgrid%lat(ipy) < nlat_rscl(xclosest)
            else
               rescale_loc = .false.
            end if

            iparallel = 0
            
            !------------------------------------------------------------------------------!
            !      POLYGON level variables.                                                !
            !------------------------------------------------------------------------------!
            !----- Load 1D dataset. -------------------------------------------------------!
            dsetrank = 1
            globdims(1) = int(dset_npolygons_global,8)
            chnkdims(1) = 1_8
            chnkoffs(1) = int(py_index - 1,8)
            memdims(1)  = 1_8
            memoffs(1)  = 0_8
            memsize(1)  = 1_8


            !---- The ipy:ipy notation is needed for ifort when checking interfaces. ------!
            call hdf_getslab_i(cgrid%load_adjacency(ipy:ipy),'LOAD_ADJACENCY '             &
                              ,dsetrank,iparallel,.true.,foundvar)
            call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.       &
                              ,foundvar)
            
            !----- Load the workload (2D). ------------------------------------------------!
            dsetrank    = 2
            globdims(1) = int(13,8)
            chnkdims(1) = int(13,8)
            memdims(1)  = int(13,8)
            memsize(1)  = int(13,8)
            chnkoffs(1) = 0_8
            memoffs(1)  = 0_8
            globdims(2) = int(dset_npolygons_global,8)
            chnkdims(2) = 1_8
            chnkoffs(2) = int(py_index - 1,8)
            memdims(2)  = 1_8
            memsize(2)  = 1_8
            memoffs(2)  = 0_8
            call hdf_getslab_r(cgrid%workload(:,ipy),'WORKLOAD ',dsetrank                  &
                              ,iparallel,.false.,foundvar)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check whether the input data had lakes or not.                           !
            !------------------------------------------------------------------------------!
            allocate(islakesite(pysi_n(py_index)))
            islakesite = 0
            !----- Load the lakesite data--------------------------------------------------!
            dsetrank = 1_8
            globdims(1) = int(dset_nsites_global,8)
            chnkdims(1) = int(pysi_n(py_index),8)
            chnkoffs(1) = int(pysi_id(py_index) - 1,8)
            memdims(1)  = int(pysi_n(py_index),8)
            memsize(1)  = int(pysi_n(py_index),8)
            memoffs(1)  = 0_8
            call hdf_getslab_i(islakesite,'ISLAKESITE ',dsetrank,iparallel,.false.,foundvar)
            ndry_sites = int(pysi_n(py_index))-sum(islakesite)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Site level.  Here we allocate temporary site variables that will grab    !
            ! the soil type information.  We then copy the data with the closest soil      !
            ! texture properties to the definite site, preserving the previously assigned  !
            ! area.                                                                        !
            !------------------------------------------------------------------------------!
            nsites_inp = ndry_sites
            allocate (this_ntext        (dset_nzg,nsites_inp))
            allocate (tpoly_area        (         nsites_inp))
            allocate (tpoly_moist_f     (         nsites_inp))
            allocate (tpoly_moist_w     (         nsites_inp))
            allocate (tpoly_elevation   (         nsites_inp))
            allocate (tpoly_slope       (         nsites_inp))
            allocate (tpoly_aspect      (         nsites_inp))
            allocate (tpoly_TCI         (         nsites_inp))
            allocate (tpoly_patch_count (         nsites_inp))
            allocate (tpoly_lsl         (         nsites_inp))
            allocate (tpoly_ntext_soil  (     nzg,nsites_inp))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Loop over the sites, seeking only those that are land sites.             !
            !------------------------------------------------------------------------------!
            is = 0
            siteloop1: do isi=1,pysi_n(py_index)
               if (islakesite(isi) == 0) then
                  is = is + 1

                  !------------------------------------------------------------------------!
                  !      Reset the HDF5 auxiliary variables before moving to the next      !
                  ! level.                                                                 !
                  !------------------------------------------------------------------------!
                  globdims = 0_8
                  chnkdims = 0_8
                  chnkoffs = 0_8
                  memoffs  = 0_8
                  memdims  = 0_8
                  memsize  = 1_8
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      SITE level variables.                                             !
                  !------------------------------------------------------------------------!
                  !----- Load 1D dataset. -------------------------------------------------!
                  dsetrank    = 1_8
                  globdims(1) = int(dset_nsites_global,8)
                  chnkdims(1) = int(1,8)
                  chnkoffs(1) = int(pysi_id(py_index) - 2 + isi,8)
                  memdims (1) = int(1,8)
                  memsize (1) = int(1,8)
                  memoffs (1) = 0_8

                  call hdf_getslab_r( tpoly_area       (is:is)     , 'AREA_SI '            &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_r( tpoly_moist_f    (is:is)     , 'MOIST_F '            &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_r( tpoly_moist_W    (is:is)     , 'MOIST_W '            &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_r( tpoly_elevation  (is:is)     , 'ELEVATION '          &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_r( tpoly_slope      (is:is)     , 'SLOPE '              &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_r( tpoly_aspect     (is:is)     , 'ASPECT '             &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_r( tpoly_TCI        (is:is)     , 'TCI '                &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_i( tpoly_patch_count(is:is)     , 'PATCH_COUNT '        &
                                    , dsetrank, iparallel, .true.,foundvar)
                  call hdf_getslab_i( tpoly_lsl        (is:is)     ,'LSL '                 &
                                    , dsetrank, iparallel, .true.,foundvar)
                  !------------------------------------------------------------------------!


                  !----- Load 2D dataset, currently just the soil texture. ----------------!
                  dsetrank     = 2_8
                  globdims(1)  = int(dset_nzg,8)
                  chnkdims(1)  = int(dset_nzg,8)
                  memdims(1)   = int(dset_nzg,8)
                  memsize(1)   = int(dset_nzg,8)
                  chnkoffs(1)  = 0_8
                  memoffs(1)   = 0_8
                  globdims(2)  = int(dset_nsites_global,8)
                  chnkdims(2)  = int(1,8)
                  chnkoffs(2)  = int(pysi_id(py_index) - 2 + isi,8)
                  memdims(2)   = int(1,8)
                  memsize(2)   = int(1,8)
                  memoffs(2)   = 0_8
                  call hdf_getslab_i( this_ntext     (:,is)     , 'NTEXT_SOIL '            &
                                    , dsetrank, iparallel, .true.,foundvar)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      The input file may have different number of soil layers than this !
                  ! simulation.  This is not a problem at this point because the soil maps !
                  ! don't have soil texture profiles, but it may become an issue for sites !
                  ! with different soil types along the profile.  Feel free to improve the !
                  ! code...  For the time being, we assume here that there is only one     !
                  ! soil type, so all that we need is to save one layer for each site.     !
                  !------------------------------------------------------------------------!

                  tpoly_ntext_soil(nzg,is) = this_ntext(dset_nzg,is)
                  
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do siteloop1
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Loop over all sites and fill the patch-level variables.                  !
            !------------------------------------------------------------------------------!
            siteloop2: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)

               nsoil = cpoly%ntext_soil(nzg,isi)

               !---------------------------------------------------------------------------!
               !     Loop over the sites, pick up the closest one.                         !
               !---------------------------------------------------------------------------!
               textdist_min = huge(1.)
               is_try       = 0
               do isi_try = 1, nsites_inp
                  if (islakesite(isi_try) == 0) then
                     is_try    = is_try + 1

                     nsoil_try = tpoly_ntext_soil(nzg,is_try)

                     !---------------------------------------------------------------------!
                     !     Find the "distance" between the two sites based on the sand and !
                     ! clay contents.                                                      !
                     !---------------------------------------------------------------------!
                     textdist_try = (soil(nsoil_try)%xsand - soil(nsoil)%xsand) ** 2       &
                                  + (soil(nsoil_try)%xclay - soil(nsoil)%xclay) ** 2
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Hold this site in case the "distance" is the minimum so far.    !
                     !---------------------------------------------------------------------!
                     if (textdist_try < textdist_min) then
                        isi_best     = isi_try
                        is_best      = is_try
                        textdist_min = textdist_try
                     end if
                     !---------------------------------------------------------------------!
                  end if
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Copy the variables we just loaded for the temporary site to the      !
               ! actual site, with two exceptions: area and lowest soil layer.  This is    !
               ! because we are not running a history run, so these values are likely to   !
               ! be different and we want them to be based on the user settings rather     !
               ! than the old run setting.                                                 !
               !---------------------------------------------------------------------------!
               cpoly%moist_f     (isi) = tpoly_moist_f     (is_best)
               cpoly%moist_w     (isi) = tpoly_moist_w     (is_best)
               cpoly%elevation   (isi) = tpoly_elevation   (is_best)
               cpoly%slope       (isi) = tpoly_slope       (is_best)
               cpoly%aspect      (isi) = tpoly_aspect      (is_best)
               cpoly%TCI         (isi) = tpoly_TCI         (is_best)
               cpoly%patch_count (isi) = tpoly_patch_count (is_best)
               !---------------------------------------------------------------------------!



               !----- Reset the HDF5 auxiliary variables before moving to the next level. -!
               globdims = 0_8
               chnkdims = 0_8
               chnkoffs = 0_8
               memoffs  = 0_8
               memdims  = 0_8
               memsize  = 1_8

               !----- Calculate the index of this site data in the HDF5 file. -------------!
               si_index = pysi_id(py_index) + isi_best - 1

               if (sipa_n(si_index) > 0) then

                  !----- Fill 1D polygon (site unique) level variables. -------------------!
                  call allocate_sitetype(csite,sipa_n(si_index))

                  iparallel = 0
                  
                  dsetrank = 1
                  globdims(1) = int(dset_npatches_global,8)
                  chnkdims(1) = int(csite%npatches,8)
                  chnkoffs(1) = int(sipa_id(si_index) - 1,8)
                  memdims(1)  = int(csite%npatches,8)
                  memsize(1)  = int(csite%npatches,8)
                  memoffs(1)  = 0

                  call hdf_getslab_i(csite%dist_type         ,'DIST_TYPE '                 &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%age               ,'AGE '                       &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%area              ,'AREA '                      &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%sum_dgd           ,'SUM_DGD '                   &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%sum_chd           ,'SUM_CHD '                   &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%fast_soil_C       ,'FAST_SOIL_C '               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%slow_soil_C       ,'SLOW_SOIL_C '               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%fast_soil_N       ,'FAST_SOIL_N '               &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%structural_soil_C ,'STRUCTURAL_SOIL_C '         &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%structural_soil_L ,'STRUCTURAL_SOIL_L '         &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N '        &
                                    ,dsetrank,iparallel,.true.,foundvar)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Check whether the history file is new or old.  We determine this   !
                  ! by searching for variable plantation.  In case this variable isn't     !
                  ! present, then it must be the new history.   Otherwise we correct the   !
                  ! indices for secondary forests.                                         !
                  !------------------------------------------------------------------------!
                  allocate (plantation(csite%npatches))
                  plantation(:) = 0
                  call hdf_getslab_i(plantation ,'PLANTATION '                             &
                                    ,dsetrank,iparallel,.false.,foundvar)
                  if (foundvar) then
                     do ipa=1,csite%npatches
                        select case(csite%dist_type(ipa))
                        case (2)
                           !---------------------------------------------------------------!
                           !     Secondary forests are now divided in three categories.    !
                           !---------------------------------------------------------------!
                           if (plantation(ipa) == 0 .and. ianth_disturb == 0) then
                              csite%dist_type(ipa) = 5
                           else if (plantation(ipa) == 0 .and. ianth_disturb == 1) then
                              csite%dist_type(ipa) = 6
                           else
                              csite%dist_type(ipa) = 2
                           end if
                           !---------------------------------------------------------------!
                        end select
                        !------------------------------------------------------------------!
                     end do
                     !---------------------------------------------------------------------!
                  end if
                  deallocate(plantation)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Check whether area should be re-scaled.                            !
                  !------------------------------------------------------------------------!
                  if (rescale_loc) then
                     !---------------------------------------------------------------------!
                     !     Now we loop over all land use types.                             !
                     !---------------------------------------------------------------------!
                     do ilu=1,n_dist_types
                        !---- The original (old) area. ------------------------------------!
                        oldarea(ilu) = sum(csite%area,mask=csite%dist_type == ilu)

                        !------------------------------------------------------------------!
                        !     Make sure that no area is going to be zero for a given land  !
                        ! use type when the counter part is not.                           !
                        !------------------------------------------------------------------!
                        oldarea(ilu)          = max( 0.5 * min_patch_area,oldarea(ilu))
                        newarea(ilu,xclosest) = max( 0.5 * min_patch_area                  &
                                                   , newarea(ilu,xclosest))
                        !------------------------------------------------------------------!
                     end do

                     !---- Re-scale the total areas so they are both equal to one. --------!
                     oldarea(:)          = oldarea(:)          / sum(oldarea)
                     newarea(:,xclosest) = newarea(:,xclosest)                             &
                                         / sum(newarea(:,xclosest:xclosest))
                     
                     !----- Re-scale the areas of every patch. ----------------------------!
                     do ipa=1,csite%npatches
                        ilu = csite%dist_type(ipa)
                        csite%area(ipa) = csite%area(ipa) * newarea(ilu,xclosest)          &
                                                          / oldarea(ilu)
                     end do

                     !----- Just to make sure we preserve unity. --------------------------!
                     csite%area(:) = csite%area(:) / sum(csite%area)

                  end if

                  !------------------------------------------------------------------------!
                  !     Loop over all sites and fill the patch-level variables.            !
                  !------------------------------------------------------------------------!
                  patchloop: do ipa = 1,csite%npatches
                     cpatch => csite%patch(ipa)

                     !---------------------------------------------------------------------!
                     !     Reset the HDF5 auxiliary variables before moving to the next    !
                     ! level.                                                              !
                     !---------------------------------------------------------------------!
                     globdims = 0_8
                     chnkdims = 0_8
                     chnkoffs = 0_8
                     memoffs  = 0_8
                     memdims  = 0_8
                     memsize  = 1_8
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     Initialise patch-level variables that depend on the cohort      !
                     ! ones.                                                               !
                     !---------------------------------------------------------------------!
                     csite%plant_ag_biomass(ipa)  = 0.0

                     pa_index = sipa_id(si_index) + ipa - 1
                     call allocate_patchtype(cpatch,paco_n(pa_index))

                     !---------------------------------------------------------------------!
                     !     Empty patches may exist, so make sure that this part is called  !
                     ! only when there are cohorts.                                        !
                     !---------------------------------------------------------------------!
                     if (cpatch%ncohorts > 0) then
                        !----- First the 1-D variables. -----------------------------------!
                        dsetrank = 1
                        globdims(1) = int(dset_ncohorts_global,8)
                        chnkdims(1) = int(cpatch%ncohorts,8)
                        chnkoffs(1) = int(paco_id(pa_index) - 1,8)
                        memdims(1)  = int(cpatch%ncohorts,8)
                        memsize(1)  = int(cpatch%ncohorts,8)
                        memoffs(1)  = 0_8

                        call hdf_getslab_r(cpatch%dbh             ,'DBH '                  &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_r(cpatch%bdead           ,'BDEAD '                &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_i(cpatch%pft             ,'PFT '                  &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_r(cpatch%nplant          ,'NPLANT '               &
                                          ,dsetrank,iparallel,.true.,foundvar)

                        !------------------------------------------------------------------!
                        !    Find derived properties from Bdead.  In the unlikely case     !
                        ! that bdead is zero, then we use DBH as the starting point.  In   !
                        ! both cases we assume that plants are in allometry.               !
                        !------------------------------------------------------------------!
                        do ico=1,cpatch%ncohorts
                           ipft = cpatch%pft(ico)

                           if (igrass == 1 .and. is_grass(ipft)                            &
                                           .and. cpatch%bdead(ico)>0.0) then
                              !-- if the initial file was running with igrass = 0, bdead   !
                              ! should be nonzero.  If the new run has igrass = 1, bdead   !
                              ! is set to zero and the mass is discarded                   !
                              cpatch%dbh(ico)   = max(cpatch%dbh(ico),min_dbh(ipft))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                              cpatch%bdead(ico) = 0.0
                              
                           else if (cpatch%bdead(ico) > 0.0 .and. igrass == 0) then
                              ! grasses have bdead in both input and current run (igrass=0)
                              cpatch%bdead(ico) = max(cpatch%bdead(ico),min_bdead(ipft))
                              cpatch%dbh(ico)   = bd2dbh(ipft,cpatch%bdead(ico))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                           else 
                              ! it is either a new grass (igrass=1) in the initial file,   !
                              ! or the value for bdead is missing from the files           !
                              cpatch%dbh(ico)   = max(cpatch%dbh(ico),min_dbh(ipft))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                              cpatch%bdead(ico) = dbh2bd(cpatch%dbh  (ico),ipft)
                           end if

                           cpatch%bleaf(ico)  = size2bl( cpatch%dbh (ico)                  &
                                                       , cpatch%hite(ico)                  &
                                                       , ipft )

                           !----- Find the other pools. -----------------------------------!
                           salloc  = (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))
                           salloci = 1.0 / salloc
                           cpatch%balive  (ico)  = cpatch%bleaf(ico)  * salloc
                           cpatch%broot    (ico) = cpatch%balive(ico) * q(ipft) * salloci
                           cpatch%bsapwooda(ico) = cpatch%balive(ico) * qsw(ipft)          &
                                                 * cpatch%hite(ico) * salloci * agf_bs(ipft)
                           cpatch%bsapwoodb(ico) = cpatch%balive(ico) * qsw(ipft)          &
                                                 * cpatch%hite(ico) * salloci              &
                                                 * (1.-agf_bs(ipft))
                           cpatch%bstorage(ico)  = 0.0
                           cpatch%phenology_status(ico) = 0
                           
                        end do

                        !------------------------------------------------------------------!
                        !     Carbon balance variables.                                    !
                        ! MLO.  I commented this out because a restart is likely to have   !
                        !       different settings and different environment, so it's      !
                        !       likely that the system will reach a different equilibrium. !
                        !       For simplicity, we just use the typical initialisation in  !
                        !       which we give plants one year to adjust to the new         !
                        !       conditions.                                                !
                        !------------------------------------------------------------------!
                        ! dsetrank    = 2
                        ! globdims(1) = 13_8
                        ! chnkdims(1) = 13_8
                        ! chnkoffs(1) = 0_8
                        ! memdims(1)  = 13_8
                        ! memsize(1)  = 13_8
                        ! memoffs(2)  = 0_8
                        ! globdims(2) = int(dset_ncohorts_global,8)
                        ! chnkdims(2) = int(cpatch%ncohorts,8)
                        ! chnkoffs(2) = int(paco_id(pa_index) - 1,8)
                        ! memdims(2)  = int(cpatch%ncohorts,8)
                        ! memsize(2)  = int(cpatch%ncohorts,8)
                        ! memoffs(2)  = 0_8

                        ! call hdf_getslab_r(cpatch%cb    ,'CB '                           &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        ! call hdf_getslab_r(cpatch%cb_lightmax,'CB_LIGHTMAX '             &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        ! call hdf_getslab_r(cpatch%cb_moistmax,'CB_MOISTMAX '             &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        do ico = 1, cpatch%ncohorts
                           cpatch%cb          (1:12,ico) = 1.0
                           cpatch%cb_lightmax (1:12,ico) = 1.0
                           cpatch%cb_moistmax (1:12,ico) = 1.0
                           cpatch%cb          (  13,ico) = 0.0
                           cpatch%cb_lightmax (  13,ico) = 0.0
                           cpatch%cb_moistmax (  13,ico) = 0.0
                        end do
                        !------------------------------------------------------------------!

                        cohortloop: do ico=1,cpatch%ncohorts
                           !---------------------------------------------------------------!
                           !    We will now check the PFT of each cohort, so we determine  !
                           ! if this is a valid PFT.  If not, then we must decide what we  !
                           ! should do...                                                  !
                           !---------------------------------------------------------------!
                           if (.not. include_pft(cpatch%pft(ico))) then
                              select case(pft_1st_check)
                              case (0)
                                 !----- Stop the run. -------------------------------------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,' and it is not in your include_these_pft...'
                                 call fatal_error('Invalid PFT in history file'            &
                                                 ,'read_ed21_history_file'                 &
                                                 ,'ed_read_ed21_history.F90')

                              case (1)
                                 !----- Include the unexpected PFT in the list. -----------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                      'I found a cohort with PFT=',cpatch%pft(ico)         &
                                     ,'... Including this PFT in your include_these_pft...'
                                 include_pft(cpatch%pft(ico))          = .true.
                                 include_these_pft(count(include_pft)) = cpatch%pft(ico)

                                 call sort_up(include_these_pft,n_pft)

                                 if (is_grass(cpatch%pft(ico))) then
                                    include_pft_ag(cpatch%pft(ico)) = .true.
                                 end if

                              case (2)
                                 !----- Ignore the unexpect PFT. --------------------------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,'... Ignoring it...'
                                 !---------------------------------------------------------!
                                 !    The way we will ignore this cohort is by setting its !
                                 ! nplant to zero, and calling the "terminate_cohorts"     !
                                 ! subroutine right after this.                            !
                                 !---------------------------------------------------------!
                                 cpatch%nplant(ico) = 0.
                              end select
                           end if

                           !---------------------------------------------------------------!
                           !     Make sure that the biomass won't lead to FPE.  This       !
                           ! should never happen when using a stable ED-2.1 version, but   !
                           ! older versions had "zombie" cohorts.  Here we ensure that     !
                           ! the model initialises with stable numbers whilst ensuring     !
                           ! that the cohorts will be eliminated.                          !
                           !---------------------------------------------------------------!
                           if (cpatch%balive(ico) > 0.            .and.                    &
                               cpatch%balive(ico) < tiny_biomass) then
                              cpatch%balive(ico) = tiny_biomass
                           end if
                           if (cpatch%bleaf(ico) > 0.            .and.                     &
                               cpatch%bleaf(ico) < tiny_biomass) then
                              cpatch%bleaf(ico) = tiny_biomass
                           end if
                           if (cpatch%broot(ico) > 0.            .and.                     &
                               cpatch%broot(ico) < tiny_biomass) then
                              cpatch%broot(ico) = tiny_biomass
                           end if
                           if (cpatch%bsapwooda(ico) > 0.        .and.                     &
                               cpatch%bsapwooda(ico) < tiny_biomass) then
                              cpatch%bsapwooda(ico) = tiny_biomass
                           end if
                           if (cpatch%bsapwoodb(ico) > 0.        .and.                     &
                               cpatch%bsapwoodb(ico) < tiny_biomass) then
                              cpatch%bsapwoodb(ico) = tiny_biomass
                           end if
                           if (cpatch%bdead(ico) > 0.            .and.                     &
                               cpatch%bdead(ico) < tiny_biomass) then
                              cpatch%bdead(ico) = tiny_biomass
                           end if
                           if (cpatch%bstorage(ico) > 0.            .and.                  &
                               cpatch%bstorage(ico) < tiny_biomass) then
                              cpatch%bstorage(ico) = tiny_biomass
                           end if
                           !---------------------------------------------------------------!


                           !----- Compute the above-ground biomass. -----------------------!
                           cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)&
                                                     ,cpatch%bsapwooda(ico),cpatch%pft(ico))

                           cpatch%basarea(ico)  = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)

                           !----- Assign LAI, WAI, and CAI --------------------------------!
                           call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)          &
                                            ,cpatch%bdead(ico),cpatch%balive(ico)          &
                                            ,cpatch%dbh(ico),cpatch%hite(ico)              &
                                            ,cpatch%pft(ico),SLA(cpatch%pft(ico))          &
                                            ,cpatch%lai(ico),cpatch%wai(ico)               &
                                            ,cpatch%crown_area(ico),cpatch%bsapwooda(ico))


                           !----- Update the derived patch-level variables. ---------------!
                           csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)       &
                                                       + cpatch%agb(ico)*cpatch%nplant(ico)

                           !----- Initialise the other cohort level variables. ------------!
                           call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                        end do cohortloop

                        !------------------------------------------------------------------!
                        !    Eliminate any "unwanted" cohort (i.e., those which nplant was !
                        ! set to zero so it would be removed).                             !
                        !------------------------------------------------------------------!
                        call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)

                     end if
                  end do patchloop
               else
                  !----- This should never happen, but, just in case... -------------------!
                  call fatal_error('A site with no patches was found...'                   &
                                  ,'read_ed21_history_file','ed_read_ed21_history.F90')
               end if

               !----- Initialise the other patch-level variables. -------------------------!
               call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))

            end do siteloop2


            !----- Initialise some site-level variables. ----------------------------------!
            call init_ed_site_vars(cpoly,cgrid%lat(ipy))
            !------------------------------------------------------------------------------!

            !----- Deallocate the temporary polygon and soil structure. -------------------!
            deallocate (this_ntext        )
            deallocate (tpoly_area        )
            deallocate (tpoly_moist_f     )
            deallocate (tpoly_moist_w     )
            deallocate (tpoly_elevation   )
            deallocate (tpoly_slope       )
            deallocate (tpoly_aspect      )
            deallocate (tpoly_TCI         )
            deallocate (tpoly_patch_count )
            deallocate (tpoly_lsl         )
            deallocate (tpoly_ntext_soil  )
            deallocate (islakesite        )
            !------------------------------------------------------------------------------!

         end do polyloop

         !----- Close the dataset. --------------------------------------------------------!
         call h5fclose_f(file_id, hdferr)
         if (hdferr /= 0) then
            write (unit=*,fmt='(a,1x,a)') 'File: ',trim(hnamel)
            write (unit=*,fmt='(a)'     ) 'Problem: Failed closing the HDF5 dataset.'
            write (unit=*,fmt='(a,1x,i4)') 'HDFerr: ',hdferr
            call fatal_error('Could not close the HDF file'                                &
                            ,'read_ed21_history_file','ed_read_ed21_history.F90')
         end if

         deallocate(paco_n    ,paco_id  )
         deallocate(sipa_n    ,sipa_id  )
         deallocate(pysi_n    ,pysi_id  )
      end do rstfileloop

      !----- Initialise the other polygon-level variables. --------------------------------!
      call init_ed_poly_vars(cgrid)
     
      !----- Deallocate the closest index vector. -----------------------------------------!
      deallocate(pclosest,psrcfile)
   end do gridloop
   
   !----- Turn off automatic error printing. ----------------------------------------------!
   call h5eset_auto_f(1,hdferr)

   !----- Close the HDF5 environment. -----------------------------------------------------!
   call h5close_f(hdferr)

#else
   call fatal_error ('You cannot restart with ED-2.1 without using HDF5...'                &
                    ,'read_ed21_history_unstruct','ed_read_ed21_history.F90')
#endif   


   return
end subroutine read_ed21_history_unstruct
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will read the ED-2.1 restart files via IED_INIT_MODE=10              !
! With this specification, the user wants all site data from the closest polygon.  This    !
! method also transfers soils information, geophysical information and soil moisture       !
! information along with the vegetation structure and composition from the donor polygon.  !
! Aside from doing a history restart, this is the closest method to getting exact clones   !
! of the donor polygons.                                                                   !
!==========================================================================================!
subroutine read_ed21_polyclone

#if USE_HDF5
   use hdf5    
#endif
   use ed_max_dims    , only : n_pft                   & ! intent(in)
                             , huge_polygon            & ! intent(in)
                             , str_len                 & ! intent(in)
                             , n_dist_types            & ! intent(in)
                             , maxfiles                & ! intent(in)
                             , maxlist                 ! ! intent(in)
   use pft_coms       , only : SLA                     & ! intent(in)
                             , q                       & ! intent(in)
                             , qsw                     & ! intent(in)
                             , hgt_min                 & ! intent(in)
                             , min_dbh                 & ! intent(in)
                             , min_bdead               & ! intent(in)
                             , is_grass                & ! intent(in)
                             , include_pft             & ! intent(in)
                             , include_pft_ag          & ! intent(in)
                             , pft_1st_check           & ! intent(in)
                             , include_these_pft       & ! intent(in)
                             , agf_bs                  & ! intent(in)
                             , min_cohort_size         ! ! intent(in)
   use ed_misc_coms   , only : sfilin                  & ! intent(in)
                             , current_time            & ! intent(in)
                             , imonthh                 & ! intent(in)
                             , iyearh                  & ! intent(in)
                             , idateh                  & ! intent(in)
                             , itimeh                  & ! intent(in)
                             , ied_init_mode           & ! intent(in)
                             , max_poi99_dist          & ! intent(in)
                             , igrass
   use ed_state_vars  , only : polygontype             & ! variable type
                             , sitetype                & ! variable type
                             , patchtype               & ! variable type
                             , edtype                  & ! variable type
                             , edgrid_g                & ! subroutine
                             , allocate_polygontype    & ! subroutine
                             , allocate_sitetype       & ! subroutine
                             , allocate_patchtype      ! ! subroutine
   use grid_coms      , only : ngrids                  & ! intent(in)
                             , nzg                     ! ! intent(in)
   use consts_coms    , only : pio4                    ! ! intent(in)
   use hdf5_coms      , only : file_id                 & ! intent(in)
                             , dset_id                 & ! intent(in)
                             , dspace_id               & ! intent(in)
                             , plist_id                & ! intent(in)
                             , globdims                & ! intent(in)
                             , chnkdims                & ! intent(in)
                             , chnkoffs                & ! intent(in)
                             , memdims                 & ! intent(in)
                             , memoffs                 & ! intent(in)
                             , memsize                 ! ! intent(in)
   use allometry      , only : area_indices            & ! function
                             , ed_biomass              & ! function
                             , bd2dbh                  & ! function
                             , dbh2h                   & ! function
                             , dbh2bd                  & ! function
                             , size2bl                 ! ! function
   use fuse_fiss_utils, only : terminate_cohorts       ! ! subroutine
   use disturb_coms   , only : ianth_disturb           & ! intent(in)
                             , lu_rescale_file         & ! intent(in)
                             , min_patch_area          ! ! intent(in)
   use soil_coms      , only : soil                    & ! intent(in)
                             , slzt                    &
                             , isoilcol

   implicit none

#if (USE_HDF5)
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)          , pointer                              :: cgrid
   type(polygontype)     , pointer                              :: cpoly
   type(sitetype)        , pointer                              :: csite
   type(patchtype)       , pointer                              :: cpatch
   character(len=str_len), dimension(maxlist)                   :: full_list
   character(len=str_len), dimension(maxfiles)                  :: histo_list
   character(len=1)                                             :: vnam
   character(len=3)                                             :: cgr
   character(len=str_len)                                       :: hnamel
   integer               , dimension(maxfiles)                  :: ngridpoly
   integer               , dimension(huge_polygon)              :: pyfile_list
   integer               , dimension(huge_polygon)              :: pyindx_list
   integer               , dimension(  :)         , allocatable :: plantation
   integer               , dimension(  :)         , allocatable :: pclosest
   integer               , dimension(  :)         , allocatable :: psrcfile
   integer               , dimension(  :)         , allocatable :: pysi_n
   integer               , dimension(  :)         , allocatable :: pysi_id
   integer               , dimension(  :)         , allocatable :: sipa_n
   integer               , dimension(  :)         , allocatable :: sipa_id
   integer               , dimension(  :)         , allocatable :: paco_n
   integer               , dimension(  :)         , allocatable :: paco_id
   integer               , dimension(  :)         , allocatable :: this_ntext
   integer               , dimension(  :)         , allocatable :: islakesite
   real    :: mindist
   integer :: minind
   real                  ,  dimension(:, :)        , allocatable :: this_soil_water
   real                  ,  dimension(  :)        , allocatable :: dset_slzm
   integer                , dimension(  :)        , allocatable :: slz_match
   integer                                                      :: year
   integer                                                      :: igr
   integer                                                      :: ipy
   integer                                                      :: isi
   integer                                                      :: is
   integer                                                      :: ipa
   integer                                                      :: ico
   integer                                                      :: nsoil
   integer                                                      :: nsites_inp
   integer                                                      :: xclosest
   integer                                                      :: nflist
   integer                                                      :: nhisto
   integer                                                      :: nrescale
   integer                                                      :: k,kd,km
   integer                                                      :: nf
   integer                                                      :: dset_npolygons_global
   integer                                                      :: dset_nsites_global
   integer                                                      :: dset_npatches_global
   integer                                                      :: dset_ncohorts_global
   integer                                                      :: dset_nzg
   integer                                                      :: ngr
   integer                                                      :: ifpy
   integer                                                      :: ipft
   integer                                                      :: ipya
   integer                                                      :: ipyz
   integer                                                      :: ierr
   integer                                                      :: ilu
   integer                                                      :: py_index
   integer                                                      :: si_index
   integer                                                      :: pa_index
   integer                                                      :: dsetrank,iparallel
   integer                                                      :: hdferr
   integer                                                      :: total_grid_py
   integer                                                      :: poi_minloc
   integer                                                      :: ngp1
   integer                                                      :: ndry_sites
   logical                                                      :: exists
   logical                                                      :: rescale_glob
   logical                                                      :: rescale_loc
   logical                                                      :: foundvar
   real                  , dimension(  :)         , allocatable :: pdist
   !real                  , dimension(huge_polygon)              :: pdist
   real                  , dimension(huge_polygon)              :: plon_list
   real                  , dimension(huge_polygon)              :: plat_list
   real                  , dimension(huge_polygon)              :: dist_rscl
   real                  , dimension(huge_polygon)              :: wlon_rscl
   real                  , dimension(huge_polygon)              :: clon_rscl
   real                  , dimension(huge_polygon)              :: elon_rscl
   real                  , dimension(huge_polygon)              :: slat_rscl
   real                  , dimension(huge_polygon)              :: clat_rscl
   real                  , dimension(huge_polygon)              :: nlat_rscl
   real                  , dimension(n_dist_types,huge_polygon) :: newarea
   real                  , dimension(n_dist_types)              :: oldarea
   real                                                         :: textdist_try
   real                                                         :: textdist_min
   real                                                         :: dummy
   real                                                         :: elim_nplant
   real                                                         :: elim_lai
   real                                                         :: salloc
   real                                                         :: salloci
   real                                                         :: sum_poly_area
   real                       :: zmin
   real                       :: fa
   real                       :: fb
   real                       :: te
   real                       :: t0
   real                       :: k0
   integer                    :: sc


   !----- Local constants. ----------------------------------------------------------------!
   real                                           , parameter   :: tiny_biomass = 1.e-20
   !----- External functions. -------------------------------------------------------------!
   real                                           , external    :: dist_gc 
   !---------------------------------------------------------------------------------------!

   !----- Open the HDF environment. -------------------------------------------------------!
   call h5open_f(hdferr)

   !----- Turn off automatic error printing. ----------------------------------------------!
   call h5eset_auto_f(0,hdferr)
   !---------------------------------------------------------------------------------------!
   !     Big loop over all grids.                                                          !
   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids


      !----- Retrieve all files with the specified prefix. --------------------------------!
      call ed_filelist(full_list,sfilin(igr),nflist)
      !----- Check every file and save only those that are actually history files. --------!
      call ed21_fileinfo(nflist,full_list,nhisto,histo_list)

      !----- Initialize the dimensional control variables for the H5 slabs. ---------------!
      globdims = 0_8
      chnkdims = 0_8
      chnkoffs = 0_8
      memoffs  = 0_8
      memdims  = 0_8
      memsize  = 1_8  

      !------------------------------------------------------------------------------------!
      !     First thing, we go through every file, open, and retrieve only some polygon-   !
      ! level information.  This will be used for mapping the files later.                 !
      !------------------------------------------------------------------------------------!
      ngridpoly  (:) = 0
      pyfile_list(:) = 0
      pyindx_list(:) = 0
      ipyz           = 0
      lonlatloop: do nf=1,nhisto
         hnamel = histo_list(nf)

         !----- Open the HDF5 file. -------------------------------------------------------!
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a,1x,i4)') ' - Error opening HDF5 file - error - ',hdferr
            write(unit=*,fmt='(a,1x,a)') ' - File name: ',trim(hnamel)
            call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)            &
                            ,'read_ed21_history_file','ed_read_ed21_history.F90')
         end if

         !----- Retrieve the number of polygons in this file. -----------------------------!
         call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
         
         !----- Determine the global position of these polygons. --------------------------!
         ipya = ipyz + 1
         ipyz = ipyz + dset_npolygons_global

         !----- In case we are meshing grids ----------------------------------------------!
         ngridpoly(nf) = dset_npolygons_global

         !------ Here we save the file where we can find these polygons. ------------------!
         do ipy=ipya,ipyz
            pyfile_list(ipy) = nf
            pyindx_list(ipy) = ipy-ipya+1
         end do

         !---------------------------------------------------------------------------------!
         !      Retrieve the polygon coordinates data.                                     !
         !---------------------------------------------------------------------------------!
         globdims(1) = int(dset_npolygons_global,8)

         call h5dopen_f(file_id,'LONGITUDE', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_REAL,plon_list(ipya:ipyz),globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
         
         call h5dopen_f(file_id,'LATITUDE', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_REAL,plat_list(ipya:ipyz),globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)
         
         call h5fclose_f(file_id, hdferr)
         if (hdferr /= 0) then
            write (unit=*,fmt='(a,1x,a)') 'File: ',trim(hnamel)
            write (unit=*,fmt='(a)'     ) 'Problem: Failed closing the HDF5 dataset.'
            write (unit=*,fmt='(a,1x,i4)') 'HDFerr: ',hdferr
            call fatal_error('Could not close the HDF file'                                &
                            ,'read_ed21_history_unstruct','ed_read_ed21_history.F90')
         end if
      end do lonlatloop

      total_grid_py = ipyz
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     If this is a simulation with anthropogenic disturbance, check and read the re- !
      ! scale file if it exists.                                                           !
      !------------------------------------------------------------------------------------!

      wlon_rscl(:) =  190.
      elon_rscl(:) = -190.
      slat_rscl(:) =  100.
      nlat_rscl(:) = -100.
      inquire(file=trim(lu_rescale_file(igr)),exist=exists)
      rescale_glob = ianth_disturb == 1 .and. exists
      nrescale = 0
      if (rescale_glob) then
         open (unit=13,file=trim(lu_rescale_file(igr)),status='old',action='read')
         read (unit=13,fmt=*)
         readrescale: do
            nrescale = nrescale + 1
            read (unit=13,fmt=*,iostat=ierr)  wlon_rscl(nrescale),elon_rscl(nrescale)      &
                 ,slat_rscl(nrescale),nlat_rscl(nrescale)      &
                 ,clon_rscl(nrescale),clat_rscl(nrescale)      &
                 ,dummy,(newarea(ilu,nrescale),ilu=1           &
                 ,n_dist_types)
            
            if (ierr /= 0) then
               nrescale = nrescale - 1
               exit readrescale
            end if
         end do readrescale
         rescale_glob = nrescale > 0
         close (unit=13,status='keep')
      elseif (ianth_disturb == 1 .and. len_trim(lu_rescale_file(igr)) > 0) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') '  In subroutine read_ed21_history_unstruct:'
         write (unit=*,fmt='(a,1x,i5)') '  - Grid: ',igr
         write (unit=*,fmt='(a)') '  In subroutine read_ed21_history_unstruct:'
         write (unit=*,fmt='(a)') '  - File '//trim(lu_rescale_file(igr))//                &
              ' wasn''t found...'
         write (unit=*,fmt='(a)') '  - Assuming no rescaling...'
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
      end if
      !------------------------------------------------------------------------------------!
      
      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !     Now, we will go through all polygons, and we will determine which input        !
      ! polygon is the closest one to each polygon.                                        !
      !------------------------------------------------------------------------------------!
      allocate (pclosest(cgrid%npolygons),psrcfile(cgrid%npolygons))
      allocate (pdist(total_grid_py))
      nneighloop: do ipy=1,cgrid%npolygons
         !----- Reset pdist to a very large number. ---------------------------------------!
         pdist(1:total_grid_py)   = 1.e20

         do ifpy=1,total_grid_py
            pdist(ifpy) = dist_gc(plon_list(ifpy),cgrid%lon(ipy)                           &
                                 ,plat_list(ifpy),cgrid%lat(ipy))
         end do

         pclosest(ipy) = pyindx_list(minloc(pdist,dim=1))
         psrcfile(ipy) = pyfile_list(minloc(pdist,dim=1))
         
 
      end do nneighloop
      deallocate(pdist)

      !------------------------------------------------------------------------------------!
      !     Now that we have all polygons matched with their nearest neighbours, we will   !
      ! loop over all files instead of all polygons.  This is to avoid opening and closing !
      ! the files too many times.                                                          !
      !------------------------------------------------------------------------------------!
      rstfileloop: do nf=1, nhisto

         !---------------------------------------------------------------------------------!
         !    Before anything else, we check whether this file needs to be opened, i.e.,   !
         ! whether any polygon has its nearest neighbour in this file.  If not, we move    !
         ! to the next one.                                                                !
         !---------------------------------------------------------------------------------!
         if (.not. (any(psrcfile == nf))) cycle rstfileloop

         !---------------------------------------------------------------------------------!
         !     If we are here, at least one polygon closest source is in this file, open   !
         ! it and load all fill in all polygons whose nearest neighbour is in this file.   !
         !---------------------------------------------------------------------------------!
         hnamel = histo_list(nf)

         !----- Open file. ----------------------------------------------------------------!
         call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a,1x,i4)') ' - Error opening HDF5 file - error - ',hdferr
            write(unit=*,fmt='(a,1x,a)') ' - File name: ',trim(hnamel)
            call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)            &
                            ,'read_ed21_history_file','ed_read_ed21_history.F90')
         end if

         !---------------------------------------------------------------------------------!
         !      Retrieve global vector sizes and mapping tree.                             !
         !---------------------------------------------------------------------------------!
         globdims = 0_8
         chnkdims = 0_8
         chnkoffs = 0_8

         globdims(1) = 1_8

         call h5dopen_f(file_id,'NZG', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nzg,globdims, hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nsites_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npatches_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
         call h5dget_space_f(dset_id, dspace_id, hdferr)
         call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_ncohorts_global,globdims,hdferr)
         call h5sclose_f(dspace_id, hdferr)
         call h5dclose_f(dset_id, hdferr)

         globdims(1) = int(dset_npolygons_global,8)

         allocate(pysi_n(dset_npolygons_global))
         allocate(pysi_id(dset_npolygons_global))

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

         globdims(1) = int(dset_nsites_global,8)

         allocate(sipa_n(dset_nsites_global))
         allocate(sipa_id(dset_nsites_global))

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

         globdims(1) = int(dset_npatches_global,8)
         allocate(paco_n(dset_npatches_global))
         allocate(paco_id(dset_npatches_global))

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

         dsetrank = 1
         globdims(1) = dset_nzg
         chnkdims(1) = dset_nzg
         chnkoffs(1) = 0_8
         memdims(1)  = dset_nzg 
         memsize(1)  = dset_nzg
         memoffs(1)  = 0_8
        

         allocate(this_ntext(dset_nzg))
         allocate(dset_slzm(dset_nzg))
         
         call hdf_getslab_r(dset_slzm,'SLZ',dsetrank                  &
              ,iparallel,.true.,foundvar)

         ! Calculate the mid-points of the dataset soil-layers
         do kd=1,dset_nzg-1
            dset_slzm(kd) = 0.5*(dset_slzm(kd)+dset_slzm(kd+1))
         end do
         dset_slzm(dset_nzg) = 0.5*(dset_slzm(dset_nzg)+0.0)


	 polyloop: do ipy = 1,cgrid%npolygons
	 	       cpoly => cgrid%polygon(ipy)

            !----- We skip the polygon if its source polygon is not in this file. ---------!
            if (psrcfile(ipy) /= nf) cycle polyloop

            !------------------------------------------------------------------------------!
            !    Use the index corresponding to the relative position of the input polygon !
            ! in the source file to which the polygon belongs.  Use these values, and its  !
            ! children values in sites, patchs and cohorts.                                !
            !------------------------------------------------------------------------------!
            py_index = pclosest(ipy)


            !------------------------------------------------------------------------------!
            !     In case we seek to rescale, we must first check whether a scale for the  !
            ! current polygon.                                                             !
            !------------------------------------------------------------------------------!
            !----- Initialise distance and co-ordinates to non-sense numbers. -------------!
            dist_rscl(:) = 1.e+20 ! Initialise to a large distance and non-sense 
            if (rescale_glob) then
               neighbour: do k=1,nrescale
                  dist_rscl(k) = dist_gc(clon_rscl(k),cgrid%lon(ipy)                       &
                       ,clat_rscl(k),cgrid%lat(ipy))
               end do neighbour
               xclosest = minloc(dist_rscl,dim=1)
               
               rescale_loc = cgrid%lon(ipy) > wlon_rscl(xclosest) .and.                    &
                    cgrid%lon(ipy) < elon_rscl(xclosest) .and.                    &
                    cgrid%lat(ipy) > slat_rscl(xclosest) .and.                    &
                    cgrid%lat(ipy) < nlat_rscl(xclosest)
            else
               rescale_loc = .false.
            end if
 

            iparallel = 0
            
            !------------------------------------------------------------------------------!
            !      POLYGON level variables.                                                !
            !------------------------------------------------------------------------------!
            !----- Load 1D dataset. -------------------------------------------------------!
            dsetrank = 1
            globdims(1) = int(dset_npolygons_global,8)
            chnkdims(1) = 1_8
            chnkoffs(1) = int(py_index - 1,8)
            memdims(1)  = 1_8
            memoffs(1)  = 0_8
            memsize(1)  = 1_8


            !---- The ipy:ipy notation is needed for ifort when checking interfaces. ------!
            call hdf_getslab_i(cgrid%load_adjacency(ipy:ipy),'LOAD_ADJACENCY '             &
                              ,dsetrank,iparallel,.true.,foundvar)
            call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.,foundvar)
            
            !----- Load the workload (2D). ------------------------------------------------!
            dsetrank    = 2
            globdims(1) = int(13,8)
            chnkdims(1) = int(13,8)
            memdims(1)  = int(13,8)
            memsize(1)  = int(13,8)
            chnkoffs(1) = 0_8
            memoffs(1)  = 0_8
            globdims(2) = int(dset_npolygons_global,8)
            chnkdims(2) = 1_8
            chnkoffs(2) = int(py_index - 1,8)
            memdims(2)  = 1_8
            memsize(2)  = 1_8
            memoffs(2)  = 0_8
            call hdf_getslab_r(cgrid%workload(:,ipy),'WORKLOAD ',dsetrank                  &
                              ,iparallel,.false.,foundvar)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check whether the input data had lakes or not.                           !
            !------------------------------------------------------------------------------!
            allocate(islakesite(pysi_n(py_index)))
            islakesite = 0
            !----- Load the lakesite data--------------------------------------------------!
            dsetrank = 1_8
            globdims(1) = int(dset_nsites_global,8)
            chnkdims(1) = int(pysi_n(py_index),8)
            chnkoffs(1) = int(pysi_id(py_index) - 1,8)
            memdims(1)  = int(pysi_n(py_index),8)
            memsize(1)  = int(pysi_n(py_index),8)
            memoffs(1)  = 0_8
            call hdf_getslab_i(islakesite,'ISLAKESITE ',dsetrank,iparallel,.false.,foundvar)
            ndry_sites = int(pysi_n(py_index))-sum(islakesite)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            ! Allocate the destination polygon with site level vector data                 !
            !------------------------------------------------------------------------------!
            call allocate_polygontype(cpoly,ndry_sites)
            call soil_default_fill(cgrid,igr,ipy)


            !------------------------------------------------------------------------------!
            !     Loop over the sites, seeking only those that are land sites.             !
            !------------------------------------------------------------------------------!
            is = 0
            sum_poly_area = 0.
            siteloop: do isi=1,pysi_n(py_index)
               if (islakesite(isi) == 0) then
                  is = is + 1

                  csite => cpoly%site(is)
                  
                  si_index = pysi_id(py_index)+isi-1

                  iparallel = 0
                     
                  dsetrank = 1_8
                  globdims = 0_8
                  chnkdims = 0_8
                  chnkoffs = 0_8
                  memoffs  = 0_8
                  memdims  = 0_8
                  memsize  = 1_8
                  
                  globdims(1) = int(dset_nsites_global,8)
                  chnkdims(1) = int(1,8)
                  chnkoffs(1) = int(si_index-1,8)
                  memdims(1)  = int(1,8)
                  memsize(1)  = int(1,8)
                  memoffs(1)  = 0_8

                  call hdf_getslab_i(cpoly%patch_count(is:is),'PATCH_COUNT ',dsetrank,iparallel,.true.,foundvar)  
                  call hdf_getslab_i(cpoly%sitenum(is:is),'SITENUM ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_i(cpoly%lsl(is:is),'LSL ',dsetrank,iparallel,.true.,foundvar)   
                  call hdf_getslab_i(cpoly%ncol_soil(is:is),'NCOL_SOIL ',dsetrank,iparallel,.false.,foundvar)

                  ! If this data is not available in the dataset, we should really just use
                  ! a default value.  It is probably not a good idea to use values derived from
                  ! a soils dataset earlier in the code, if we are no longer using the associated
                  ! textures.

                  if (cpoly%ncol_soil(is) == 0) then
                     write (unit=*,fmt='(a,i3)')                       &
                          'Soil color info not in ED2.1 state file, using ISOILCOL=',isoilcol
                     cpoly%ncol_soil(is)  = isoilcol
                  end if
                  

                  call hdf_getslab_r(cpoly%area(is:is),'AREA_SI ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%patch_area(is:is),'PATCH_AREA ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%elevation(is:is),'ELEVATION ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%slope(is:is),'SLOPE ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%aspect(is:is),'ASPECT ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%TCI(is:is),'TCI ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_i(cpoly%hydro_next(is:is),'HYDRO_NEXT ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_i(cpoly%hydro_prev(is:is),'HYDRO_PREV ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%moist_W(is:is),'MOIST_W ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%moist_f(is:is),'MOIST_F ',dsetrank,iparallel,.true.,foundvar)  
                  call hdf_getslab_r(cpoly%moist_tau(is:is),'MOIST_TAU ',dsetrank,iparallel,.true.,foundvar)
                  call hdf_getslab_r(cpoly%moist_zi(is:is),'MOIST_ZI ',dsetrank,iparallel,.true.,foundvar) 
                  call hdf_getslab_r(cpoly%baseflow(is:is),'BASEFLOW_SI ',dsetrank,iparallel,.true.,foundvar)

                  sum_poly_area = sum_poly_area+cpoly%area(is)

                  !----- Load 2D data soil texture
                  dsetrank     = 2_8
                  globdims(1)  = int(dset_nzg,8)
                  chnkdims(1)  = int(dset_nzg,8)
                  memdims(1)   = int(dset_nzg,8)
                  memsize(1)   = int(dset_nzg,8)
                  chnkoffs(1)  = 0_8
                  memoffs(1)   = 0_8
                  globdims(2)  = int(dset_nsites_global,8)
                  chnkdims(2)  = int(1,8)
                  chnkoffs(2)  = int(si_index-1,8)
                  memdims(2)   = int(1,8)
                  memsize(2)   = int(1,8)
                  memoffs(2)   = 0_8

                  call hdf_getslab_i( this_ntext(:), &
                       'NTEXT_SOIL ',dsetrank, iparallel, .true.,foundvar)

                  

                  !------------------------------------------------------------------------!
                  !      The input file may have different number of soil layers than this !
                  ! simulation.  This is not a problem at this point because the soil maps !
                  ! don't have soil texture profiles, but it may become an issue for sites !
                  ! with different soil types along the profile.  Feel free to improve the !
                  ! code...  For the time being, we assume here that there is only one     !
                  ! soil type, so all that we need is to save one layer for each site.     !
                  !------------------------------------------------------------------------!

                  allocate(slz_match(nzg))
                  do km=1,nzg     !km is k-model
                     mindist = 1.e10
                     minind  = -1
                     do kd=1,dset_nzg  !kd is k-dset
                        if (abs(slzt(km)-dset_slzm(kd))<mindist) then
                           mindist = abs(slzt(km)-dset_slzm(kd))
                           minind = kd
                        end if
                     end do
                     if(minind>0)then
                        slz_match(km)=minind
                     else
                        call fatal_error('Could not find soil layer match!',&
                             'ed_read_polyclone','ed_read_ed21_history.F90')
                     end if
                     cpoly%ntext_soil(km,is) = this_ntext(slz_match(km))
                  end do

                  ! ALSO, LETS RE-ASSIGN THE LSL

                  if (cpoly%lsl(is)>dset_nzg .or. cpoly%lsl(is)<1)then
                     print*,"FUNKY LSL:",cpoly%lsl(is)
                     stop
                  else
                     cpoly%lsl(is) = slz_match(cpoly%lsl(is))
                  end if
                  

                  ! We also need to set all the default properties because they were bypassed
                  ! in ed_init

                  if (sipa_n(si_index) > 0) then
                     
                     !----- Fill 1D polygon (site unique) level variables. -------------------!
                     call allocate_sitetype(csite,sipa_n(si_index))
                     
                     iparallel = 0
                     
                     dsetrank = 1
                     globdims(1) = int(dset_npatches_global,8)
                     chnkdims(1) = int(csite%npatches,8)
                     chnkoffs(1) = int(sipa_id(si_index) - 1,8)
                     memdims(1)  = int(csite%npatches,8)
                     memsize(1)  = int(csite%npatches,8)
                     memoffs(1)  = 0
                     
                     call hdf_getslab_i(csite%dist_type         ,'DIST_TYPE '                 &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%age               ,'AGE '                       &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%area              ,'AREA '                      &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%sum_dgd           ,'SUM_DGD '                   &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%sum_chd           ,'SUM_CHD '                   &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%fast_soil_C       ,'FAST_SOIL_C '               &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%slow_soil_C       ,'SLOW_SOIL_C '               &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%fast_soil_N       ,'FAST_SOIL_N '               &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%structural_soil_C ,'STRUCTURAL_SOIL_C '         &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%structural_soil_L ,'STRUCTURAL_SOIL_L '         &
                          ,dsetrank,iparallel,.true.,foundvar)
                     call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N '        &
                          ,dsetrank,iparallel,.true.,foundvar)




                     !------------------------------------------------------------------------!
                     !     Check whether the history file is new or old.  We determine this   !
                     ! by searching for variable plantation.  In case this variable isn't     !
                     ! present, then it must be the new history.   Otherwise we correct the   !
                     ! indices for secondary forests.                                         !
                     !------------------------------------------------------------------------!
                     allocate (plantation(csite%npatches))
                     plantation(:) = 0
                     call hdf_getslab_i(plantation ,'PLANTATION '                             &
                                       ,dsetrank,iparallel,.false.,foundvar)
                     if (foundvar) then
                        do ipa=1,csite%npatches
                           select case(csite%dist_type(ipa))
                           case (2)
                              !---------------------------------------------------------------!
                              !     Secondary forests are now divided in three categories.    !
                              !---------------------------------------------------------------!
                              if (plantation(ipa) == 0 .and. ianth_disturb == 0) then
                                 csite%dist_type(ipa) = 5
                              else if (plantation(ipa) == 0 .and. ianth_disturb == 1) then
                                 csite%dist_type(ipa) = 6
                              else
                                 csite%dist_type(ipa) = 2
                              end if
                              !---------------------------------------------------------------!
                           end select
                           !------------------------------------------------------------------!
                        end do
                        !---------------------------------------------------------------------!
                     end if
                     deallocate(plantation)
                     !------------------------------------------------------------------------!

                     
                     !----- Load 2D soil water
                     dsetrank     = 2_8
                     globdims(1)  = int(dset_nzg,8)
                     chnkdims(1)  = int(dset_nzg,8)
                     memdims(1)   = int(dset_nzg,8)
                     memsize(1)   = int(dset_nzg,8)
                     chnkoffs(1)  = 0_8
                     memoffs(1)   = 0_8
                     globdims(2)  = int(dset_npatches_global,8)
                     chnkdims(2)  = int(csite%npatches,8)
                     chnkoffs(2)  = int(sipa_id(si_index)-1,8)
                     memdims(2)   = int(csite%npatches,8)
                     memsize(2)   = int(csite%npatches,8)
                     memoffs(2)   = 0_8
                     
                     allocate(this_soil_water(dset_nzg,csite%npatches))
                     call hdf_getslab_r(this_soil_water,'SOIL_WATER_PA '        &
                          ,dsetrank,iparallel,.true.,foundvar)

                     ! ----------------------------------------------------------
                     ! Go through the layer centers of the model layers
                     ! Find the layer in the dataset that matches it most closely
                     ! and copy soil-water from that dataset to the model
                     ! ----------------------------------------------------------
                     do ipa=1,csite%npatches
                        do km=1,nzg
                           
!                           if( this_soil_water(slz_match(km),ipa) .gt.     &
!                                soil(cpoly%ntext_soil(km,is))%slmsts .or.  &
!                                this_soil_water(slz_match(km),ipa).le.0.0 ) then
!
!                              call fatal_error('Soil moisture is greater than porosity' &
!                                   ,'read_ed21_polyclone','ed_read_ed21_history.F90')
!                           end if
                           csite%soil_water(km,ipa) = &
                                min(this_soil_water(slz_match(km),ipa), &
                                soil(cpoly%ntext_soil(km,is))%slmsts)
                        end do
                     end do
                     deallocate(this_soil_water)

                  !------------------------------------------------------------------------!
                  !     Check whether area should be re-scaled.                            !
                  !------------------------------------------------------------------------!
                  if (rescale_loc) then
                     !---------------------------------------------------------------------!
                     !     Now we loop over all land use types.                             !
                     !---------------------------------------------------------------------!
                     do ilu=1,n_dist_types
                        !---- The original (old) area. ------------------------------------!
                        oldarea(ilu) = sum(csite%area,mask=csite%dist_type == ilu)

                        !------------------------------------------------------------------!
                        !     Make sure that no area is going to be zero for a given land  !
                        ! use type when the counter part is not.                           !
                        !------------------------------------------------------------------!
                        oldarea(ilu)          = max( 0.5 * min_patch_area,oldarea(ilu))
                        newarea(ilu,xclosest) = max( 0.5 * min_patch_area                  &
                                                   , newarea(ilu,xclosest))
                        !------------------------------------------------------------------!
                     end do

                     !---- Re-scale the total areas so they are both equal to one. --------!
                     oldarea(:)          = oldarea(:)          / sum(oldarea)
                     newarea(:,xclosest) = newarea(:,xclosest)                             &
                                         / sum(newarea(:,xclosest:xclosest))
                     
                     !----- Re-scale the areas of every patch. ----------------------------!
                     do ipa=1,csite%npatches
                        ilu = csite%dist_type(ipa)
                        csite%area(ipa) = csite%area(ipa) * newarea(ilu,xclosest)          &
                                                          / oldarea(ilu)
                     end do

                     !----- Just to make sure we preserve unity. --------------------------!
                     csite%area(:) = csite%area(:) / sum(csite%area)

                  end if

                  !------------------------------------------------------------------------!
                  !     Loop over all sites and fill the patch-level variables.            !
                  !------------------------------------------------------------------------!
                  patchloop: do ipa = 1,csite%npatches
                     cpatch => csite%patch(ipa)

                     !---------------------------------------------------------------------!
                     !     Reset the HDF5 auxiliary variables before moving to the next    !
                     ! level.                                                              !
                     !---------------------------------------------------------------------!
                     globdims = 0_8
                     chnkdims = 0_8
                     chnkoffs = 0_8
                     memoffs  = 0_8
                     memdims  = 0_8
                     memsize  = 1_8
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     Initialise patch-level variables that depend on the cohort      !
                     ! ones.                                                               !
                     !---------------------------------------------------------------------!
                     csite%plant_ag_biomass(ipa)  = 0.0

                     pa_index = sipa_id(si_index) + ipa - 1
                     call allocate_patchtype(cpatch,paco_n(pa_index))

                     !---------------------------------------------------------------------!
                     !     Empty patches may exist, so make sure that this part is called  !
                     ! only when there are cohorts.                                        !
                     !---------------------------------------------------------------------!
                     if (cpatch%ncohorts > 0) then
                        !----- First the 1-D variables. -----------------------------------!
                        dsetrank = 1
                        globdims(1) = int(dset_ncohorts_global,8)
                        chnkdims(1) = int(cpatch%ncohorts,8)
                        chnkoffs(1) = int(paco_id(pa_index) - 1,8)
                        memdims(1)  = int(cpatch%ncohorts,8)
                        memsize(1)  = int(cpatch%ncohorts,8)
                        memoffs(1)  = 0_8

                        call hdf_getslab_r(cpatch%dbh             ,'DBH '                  &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_r(cpatch%bdead           ,'BDEAD '                &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_i(cpatch%pft             ,'PFT '                  &
                                          ,dsetrank,iparallel,.true.,foundvar)
                        call hdf_getslab_r(cpatch%nplant          ,'NPLANT '               &
                                          ,dsetrank,iparallel,.true.,foundvar)

                        !------------------------------------------------------------------!
                        !    Find derived properties from Bdead.  In the unlikely case     !
                        ! that bdead is zero, then we use DBH as the starting point.  In   !
                        ! both cases we assume that plants are in allometry.               !
                        !------------------------------------------------------------------!
                        do ico=1,cpatch%ncohorts
                           ipft = cpatch%pft(ico)

                           if (igrass == 1 .and. is_grass(ipft)                            &
                                           .and. cpatch%bdead(ico)>0.0) then
                              !-- if the initial file was running with igrass = 0, bdead   !
                              ! should be nonzero.  If the new run has igrass = 1, bdead   !
                              ! is set to zero and that biomass is discarded               !
                              cpatch%dbh(ico)   = max(cpatch%dbh(ico),min_dbh(ipft))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                              cpatch%bdead(ico) = 0.0
                              
                           else if (cpatch%bdead(ico) > 0.0 .and. igrass == 0) then
                              ! grasses have bdead in both input and current run (igrass=0)
                              cpatch%bdead(ico) = max(cpatch%bdead(ico),min_bdead(ipft))
                              cpatch%dbh(ico)   = bd2dbh(ipft,cpatch%bdead(ico))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                           else 
                              ! it is either a new grass (igrass=1) in the initial file,   !
                              ! or the value for bdead is missing from the files           !
                              cpatch%dbh(ico)   = max(cpatch%dbh(ico),min_dbh(ipft))
                              cpatch%hite(ico)  = dbh2h (ipft,cpatch%dbh  (ico))
                              cpatch%bdead(ico) = dbh2bd(cpatch%dbh  (ico),ipft)
                           end if

                           cpatch%bleaf(ico)  = size2bl( cpatch%dbh (ico)                  &
                                                       , cpatch%hite(ico)                  &
                                                       , ipft )

                           !----- Find the other pools. -----------------------------------!
                           salloc  = (1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico))
                           salloci = 1.0 / salloc
                           cpatch%balive  (ico)  = cpatch%bleaf(ico) * salloc
                           cpatch%broot   (ico)  = cpatch%balive(ico) * q(ipft) * salloci
                           cpatch%bsapwooda(ico) = cpatch%balive(ico) * qsw(ipft)          &
                                                 * cpatch%hite(ico) * salloci * agf_bs(ipft)
                           cpatch%bsapwoodb(ico) = cpatch%balive(ico) * qsw(ipft)          &
                                                 * cpatch%hite(ico) * salloci              &
                                                 * (1.-agf_bs(ipft))
                           cpatch%bstorage(ico)  = 0.0
                           cpatch%phenology_status(ico) = 0
                        end do

                        !------------------------------------------------------------------!
                        !     Carbon balance variables.                                    !
                        ! MLO.  I commented this out because a restart is likely to have   !
                        !       different settings and different environment, so it's      !
                        !       likely that the system will reach a different equilibrium. !
                        !       For simplicity, we just use the typical initialisation in  !
                        !       which we give plants one year to adjust to the new         !
                        !       conditions.                                                !
                        !------------------------------------------------------------------!
                        ! dsetrank    = 2
                        ! globdims(1) = 13_8
                        ! chnkdims(1) = 13_8
                        ! chnkoffs(1) = 0_8
                        ! memdims(1)  = 13_8
                        ! memsize(1)  = 13_8
                        ! memoffs(2)  = 0_8
                        ! globdims(2) = int(dset_ncohorts_global,8)
                        ! chnkdims(2) = int(cpatch%ncohorts,8)
                        ! chnkoffs(2) = int(paco_id(pa_index) - 1,8)
                        ! memdims(2)  = int(cpatch%ncohorts,8)
                        ! memsize(2)  = int(cpatch%ncohorts,8)
                        ! memoffs(2)  = 0_8

                        ! call hdf_getslab_r(cpatch%cb    ,'CB '                           &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        ! call hdf_getslab_r(cpatch%cb_lightmax,'CB_LIGHTMAX '             &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        ! call hdf_getslab_r(cpatch%cb_moistmax,'CB_MOISTMAX '             &
                        !                   ,dsetrank,iparallel,.true.,foundvar)
                        do ico = 1, cpatch%ncohorts
                           cpatch%cb          (1:12,ico) = 1.0
                           cpatch%cb_lightmax (1:12,ico) = 1.0
                           cpatch%cb_moistmax (1:12,ico) = 1.0
                           cpatch%cb          (  13,ico) = 0.0
                           cpatch%cb_lightmax (  13,ico) = 0.0
                           cpatch%cb_moistmax (  13,ico) = 0.0
                        end do
                        !------------------------------------------------------------------!


                        cohortloop: do ico=1,cpatch%ncohorts
                           !---------------------------------------------------------------!
                           !    We will now check the PFT of each cohort, so we determine  !
                           ! if this is a valid PFT.  If not, then we must decide what we  !
                           ! should do...                                                  !
                           !---------------------------------------------------------------!
                           if (.not. include_pft(cpatch%pft(ico))) then
                              select case(pft_1st_check)
                              case (0)
                                 !----- Stop the run. -------------------------------------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,' and it is not in your include_these_pft...'
                                 call fatal_error('Invalid PFT in history file'            &
                                                 ,'read_ed21_history_file'                 &
                                                 ,'ed_read_ed21_history.F90')

                              case (1)
                                 !----- Include the unexpected PFT in the list. -----------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                      'I found a cohort with PFT=',cpatch%pft(ico)         &
                                     ,'... Including this PFT in your include_these_pft...'
                                 include_pft(cpatch%pft(ico))          = .true.
                                 include_these_pft(count(include_pft)) = cpatch%pft(ico)

                                 call sort_up(include_these_pft,n_pft)

                                 if (is_grass(cpatch%pft(ico))) then
                                    include_pft_ag(cpatch%pft(ico)) = .true.
                                 end if

                              case (2)
                                 !----- Ignore the unexpect PFT. --------------------------!
                                 write (unit=*,fmt='(a,1x,i5,1x,a)')                       &
                                       'I found a cohort with PFT=',cpatch%pft(ico)        &
                                      ,'... Ignoring it...'
                                 !---------------------------------------------------------!
                                 !    The way we will ignore this cohort is by setting its !
                                 ! nplant to zero, and calling the "terminate_cohorts"     !
                                 ! subroutine right after this.                            !
                                 !---------------------------------------------------------!
                                 cpatch%nplant(ico) = 0.
                              end select
                           end if

                           !---------------------------------------------------------------!
                           !     Make sure that the biomass won't lead to FPE.  This       !
                           ! should never happen when using a stable ED-2.1 version, but   !
                           ! older versions had "zombie" cohorts.  Here we ensure that     !
                           ! the model initialises with stable numbers whilst ensuring     !
                           ! that the cohorts will be eliminated.                          !
                           !---------------------------------------------------------------!
                           if (cpatch%balive(ico) > 0.            .and.                    &
                               cpatch%balive(ico) < tiny_biomass) then
                              cpatch%balive(ico) = tiny_biomass
                           end if
                           if (cpatch%bleaf(ico) > 0.            .and.                     &
                               cpatch%bleaf(ico) < tiny_biomass) then
                              cpatch%bleaf(ico) = tiny_biomass
                           end if
                           if (cpatch%broot(ico) > 0.            .and.                     &
                               cpatch%broot(ico) < tiny_biomass) then
                              cpatch%broot(ico) = tiny_biomass
                           end if
                           if (cpatch%bsapwooda(ico) > 0.        .and.                     &
                               cpatch%bsapwooda(ico) < tiny_biomass) then
                              cpatch%bsapwooda(ico) = tiny_biomass
                           end if
                           if (cpatch%bsapwoodb(ico) > 0.        .and.                     &
                               cpatch%bsapwoodb(ico) < tiny_biomass) then
                              cpatch%bsapwoodb(ico) = tiny_biomass
                           end if
                           if (cpatch%bdead(ico) > 0.            .and.                     &
                               cpatch%bdead(ico) < tiny_biomass) then
                              cpatch%bdead(ico) = tiny_biomass
                           end if
                           if (cpatch%bstorage(ico) > 0.            .and.                  &
                               cpatch%bstorage(ico) < tiny_biomass) then
                              cpatch%bstorage(ico) = tiny_biomass
                           end if
                           !---------------------------------------------------------------!


                           !----- Compute the above-ground biomass. -----------------------!
                           cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)&
                                                     ,cpatch%bsapwooda(ico),cpatch%pft(ico))

                           cpatch%basarea(ico)  = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)

                           !----- Assign LAI, WAI, and CAI --------------------------------!
                           call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)          &
                                            ,cpatch%bdead(ico),cpatch%balive(ico)          &
                                            ,cpatch%dbh(ico),cpatch%hite(ico)              &
                                            ,cpatch%pft(ico),SLA(cpatch%pft(ico))          &
                                            ,cpatch%lai(ico),cpatch%wai(ico)               &
                                            ,cpatch%crown_area(ico),cpatch%bsapwooda(ico))


                           !----- Update the derived patch-level variables. ---------------!
                           csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)       &
                                                       + cpatch%agb(ico)*cpatch%nplant(ico)

                           !----- Initialise the other cohort level variables. ------------!
                           call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                        end do cohortloop

                        !------------------------------------------------------------------!
                        !    Eliminate any "unwanted" cohort (i.e., those which nplant was !
                        ! set to zero so it would be removed).                             !
                        !------------------------------------------------------------------!
                        call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)

                     end if
                  end do patchloop
               else
                  !----- This should never happen, but, just in case... -------------------!
                  call fatal_error('A site with no patches was found...'                   &
                                  ,'read_ed21_history_file','ed_read_ed21_history.F90')
               end if

               !----- Initialise the other patch-level variables. -------------------------!
               call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))

               deallocate(slz_match)

            end if

         end do siteloop
         
         
         !---------------------------------------------------------------------------------!
         !     Not sure what these things do, just copying from hydrology...               !
         !---------------------------------------------------------------------------------!
         !----- Part 1. -------------------------------------------------------------------!
         Te = 0.0 
         do isi = 1,cpoly%nsites
            sc = cpoly%ntext_soil(nzg-1,isi)
            K0 = soil(sc)%slcons0
            T0 = K0 / cpoly%moist_f(isi)
            Te = Te + T0*cpoly%area(isi)
         end do
         cgrid%Te(ipy) = Te
         !----- Part 2. -------------------------------------------------------------------!
         cgrid%wbar(ipy) = 0.0
         do isi = 1,cpoly%nsites
            sc = cpoly%ntext_soil(nzg-1,isi)
            K0 = soil(sc)%slcons0
            T0 = K0 / cpoly%moist_f(isi)
            cpoly%moist_W(isi) = cpoly%TCI(isi) + log(Te) - log(T0)
            cgrid%wbar(ipy)    = cgrid%wbar(ipy) + cpoly%moist_W(isi) * cpoly%area(isi)
         end do
         !---------------------------------------------------------------------------------!
         
         ! Normalize the site area in-case not all sites were read in
         
         if (sum_poly_area>0.) then
            cpoly%area = cpoly%area/sum_poly_area
         else
            call fatal_error('Site Areas Are Nill'                                &
                 ,'read_ed21_polyclone','ed_read_ed21_history.F90')
         end if
         
         !----- Initialise some site-level variables. ----------------------------------!
         call init_ed_site_vars(cpoly,cgrid%lat(ipy))
         !------------------------------------------------------------------------------!
         
         deallocate (islakesite        )
         !------------------------------------------------------------------------------!

      end do polyloop
      
      !----- Close the dataset. --------------------------------------------------------!
      call h5fclose_f(file_id, hdferr)
      if (hdferr /= 0) then
         write (unit=*,fmt='(a,1x,a)') 'File: ',trim(hnamel)
         write (unit=*,fmt='(a)'     ) 'Problem: Failed closing the HDF5 dataset.'
         write (unit=*,fmt='(a,1x,i4)') 'HDFerr: ',hdferr
         call fatal_error('Could not close the HDF file'                                &
              ,'read_ed21_history_file','ed_read_ed21_history.F90')
      end if
      
      deallocate(paco_n    ,paco_id  )
      deallocate(sipa_n    ,sipa_id  )
      deallocate(pysi_n    ,pysi_id  )
      deallocate(this_ntext          )
      deallocate(dset_slzm           )
      
   end do rstfileloop
   
   !----- Initialise the other polygon-level variables. --------------------------------!
   call init_ed_poly_vars(cgrid)

   
   !----- Deallocate the closest index vector. -----------------------------------------!
   deallocate(pclosest,psrcfile)
end do gridloop

!----- Turn off automatic error printing. ----------------------------------------------!
call h5eset_auto_f(1,hdferr)

!----- Close the HDF5 environment. -----------------------------------------------------!
call h5close_f(hdferr)

#else
   call fatal_error ('You cannot restart with ED-2.1 without using HDF5...'                &
                    ,'read_ed21_history_unstruct','ed_read_ed21_history.F90')
#endif   


   return
 end subroutine read_ed21_polyclone



 subroutine check_rescale()


   implicit none

   




  











 end subroutine check_rescale
