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
                             , include_pft             & ! intent(in)
                             , include_pft_ag          & ! intent(in)
                             , phenology               & ! intent(in)
                             , pft_1st_check           & ! intent(in)
                             , include_these_pft       ! ! intent(in)
   use ed_misc_coms   , only : sfilin                  & ! intent(in)
                             , current_time            & ! intent(in)
                             , imonthh                 & ! intent(in)
                             , iyearh                  & ! intent(in)
                             , idateh                  & ! intent(in)
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
                             , ed_biomass              ! ! function
   use fuse_fiss_utils, only : terminate_cohorts       ! ! subroutine

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
   integer                             :: year
   integer                             :: igr
   integer                             :: ipy
   integer                             :: isi
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
   logical                             :: exists
   real, dimension(:)    , allocatable :: file_lats
   real, dimension(:)    , allocatable :: file_lons
   real                                :: minrad
   real                                :: currad
   real                                :: elim_nplant
   real                                :: elim_lai
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
      ! that lists the grid number and extension.                                          !
      !------------------------------------------------------------------------------------!
      hnamel = trim(sfilin)//"-"//trim(cgr)//".h5"


      !------------------------------------------------------------------------------------!
      !     Check whether the file we want to open actually exists.                        !
      !------------------------------------------------------------------------------------!
      inquire(file=trim(hnamel),exist=exists)
      if (.not.exists) then
         !----- It doesn't, stop the run. -------------------------------------------------!
         write (unit=*,fmt='(a,1x,a)')    'SFILIN  = ',trim(sfilin)
         write (unit=*,fmt='(a,1x,i4.4)') 'IYEARH  = ',iyearh
         write (unit=*,fmt='(a,1x,i2.2)') 'IMONTHH = ',imonthh
         write (unit=*,fmt='(a,1x,i2.2)') 'IDATEH  = ',idateh
         write (unit=*,fmt='(a,1x,i4.4)') 'ITIMEH  = ',itimeh
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
                           ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.)
         
         !----- Load the workload (2D). ---------------------------------------------------!
         dsetrank    = 2
         globdims(1) = int(13,8)
         chnkdims(1) = int(13,8)
         memdims(1)  = int(13,8)
         memsize(1)  = int(13,8)
         chnkoffs(1) = 0_8
         memoffs(1)  = 0_8
         globdims(2) = int(pysi_n(py_index),8)
         chnkdims(2) = int(pysi_n(py_index),8)
         memdims(2)  = int(pysi_n(py_index),8)
         memsize(2)  = int(pysi_n(py_index),8)
         chnkoffs(2) = 0_8
         memoffs(2)  = 0_8
         call hdf_getslab_r(cgrid%workload(:,ipy),'WORKLOAD ',dsetrank,iparallel,.false.)
         !---------------------------------------------------------------------------------!

         
         !----- Allocate the vector of sites in the polygon. ------------------------------!
         call allocate_polygontype(cpoly,pysi_n(py_index))     


         !----- Reset the HDF5 auxiliary variables before moving to the next level. -------!
         globdims = 0_8
         chnkdims = 0_8
         chnkoffs = 0_8
         memoffs  = 0_8
         memdims  = 0_8
         memsize  = 1_8

         !---------------------------------------------------------------------------------!
         !      SITE level variables.                                                      !
         !---------------------------------------------------------------------------------!
         !----- Load 1D dataset. ----------------------------------------------------------!
         dsetrank = 1_8
         globdims(1) = int(dset_nsites_global,8)
         chnkdims(1) = int(cpoly%nsites,8)
         chnkoffs(1) = int(pysi_id(py_index) - 1,8)
         memdims(1)  = int(cpoly%nsites,8)
         memsize(1)  = int(cpoly%nsites,8)
         memoffs(1)  = 0_8
         
         call hdf_getslab_r(cpoly%area       ,'AREA_SI '    ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cpoly%moist_f    ,'MOIST_F '    ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cpoly%moist_W    ,'MOIST_W '    ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cpoly%elevation  ,'ELEVATION '  ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cpoly%slope      ,'SLOPE '      ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cpoly%aspect     ,'ASPECT '     ,dsetrank,iparallel,.true.)
         call hdf_getslab_r(cpoly%TCI        ,'TCI '        ,dsetrank,iparallel,.true.)
         call hdf_getslab_i(cpoly%patch_count,'PATCH_COUNT ',dsetrank,iparallel,.true.)

         !----- Load 2D dataset. ----------------------------------------------------------!
         dsetrank     = 2_8
         globdims(1)  = int(dset_nzg,8)     ! How many layers in the dataset?
         chnkdims(1)  = int(1,8)            ! We are only extracting one layer
         memdims(1)   = int(1,8)            ! We only need memory for one layer
         memsize(1)   = int(1,8)            ! On both sides
         chnkoffs(1)  = int(dset_nzg - 1,8) ! Take the top layer, not the bottom
         memoffs(1)   = 0_8
         globdims(2)  = int(dset_nsites_global,8)
         chnkdims(2)  = int(cpoly%nsites,8)
         chnkoffs(2)  = int(pysi_id(py_index) - 1,8)
         memdims(2)   = int(cpoly%nsites,8)
         memsize(2)   = int(cpoly%nsites,8)
         memoffs(2)   = 0_8
         call hdf_getslab_i(cpoly%ntext_soil(nzg,:),'NTEXT_SOIL_SI ',dsetrank              &
                           ,iparallel,.true.)

         !---------------------------------------------------------------------------------!
         !     Loop over all sites and fill the patch-level variables.                     !
         !---------------------------------------------------------------------------------!
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Calculate the index of this site data in the HDF5 file. ----------------!
            si_index = pysi_id(py_index) + isi - 1

            if (sipa_n(si_index) > 0) then

               !---------------------------------------------------------------------------!
               !     The soil layer in this case is use defined, so take this from the     !
               ! grid level variable, and not from the dataset.                            !
               !---------------------------------------------------------------------------!
               cpoly%lsl(isi)  = cgrid%lsl(ipy)  ! Initialize lowest soil layer

               !----- Now fill the soil column based on the top layer data. ---------------!
               do k=1,nzg
                  cpoly%ntext_soil(k,isi) = cpoly%ntext_soil(nzg,isi)
               end do

               !----- Fill 1D polygon (site unique) level variables. ----------------------!
               call allocate_sitetype(csite,sipa_n(si_index))

               iparallel = 0
               
               dsetrank = 1
               globdims(1) = int(dset_npatches_global,8)
               chnkdims(1) = int(csite%npatches,8)
               chnkoffs(1) = int(sipa_id(si_index) - 1,8)
               memdims(1)  = int(csite%npatches,8)
               memsize(1)  = int(csite%npatches,8)
               memoffs(1)  = 0

               !----- Assign patch soils based off of the site level soils data. ----------!
               do k=1,nzg
                  csite%ntext_soil(k,:) = cpoly%ntext_soil(k,isi)
               end do

               call hdf_getslab_i(csite%dist_type ,'DIST_TYPE ' ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%age       ,'AGE '       ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%area      ,'AREA '      ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%sum_dgd   ,'SUM_DGD '   ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%sum_chd   ,'SUM_CHD '   ,dsetrank,iparallel,.true.)
               call hdf_getslab_i(csite%plantation,'PLANTATION ',dsetrank,iparallel,.true.)

               call hdf_getslab_r(csite%fast_soil_C       ,'FAST_SOIL_C '                  &
                                 ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%slow_soil_C       ,'SLOW_SOIL_C '                  &
                                 ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%fast_soil_N       ,'FAST_SOIL_N '                  &
                                 ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%structural_soil_C ,'STRUCTURAL_SOIL_C '            &
                                 ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%structural_soil_L ,'STRUCTURAL_SOIL_L '            &
                                 ,dsetrank,iparallel,.true.)
               call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N '           &
                                 ,dsetrank,iparallel,.true.)


               !---------------------------------------------------------------------------!
               !     Loop over all sites and fill the patch-level variables.               !
               !---------------------------------------------------------------------------!
               patchloop: do ipa = 1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !----- Initialise patch-level variables that depend on the cohort ones. -!
                  csite%lai(ipa)               = 0.0
                  csite%wpa(ipa)               = 0.0
                  csite%wai(ipa)               = 0.0
                  csite%plant_ag_biomass(ipa)  = 0.0

                  pa_index = sipa_id(si_index) + ipa - 1
                  call allocate_patchtype(cpatch,paco_n(pa_index))

                  !------------------------------------------------------------------------!
                  !     Empty patches may exist, so make sure that this part is called     !
                  ! only when there are cohorts.                                           !
                  !------------------------------------------------------------------------!
                  if (cpatch%ncohorts > 0) then
                     !----- First the 1-D variables. --------------------------------------!
                     dsetrank = 1
                     globdims(1) = int(dset_ncohorts_global,8)
                     chnkdims(1) = int(cpatch%ncohorts,8)
                     chnkoffs(1) = int(paco_id(pa_index) - 1,8)
                     memdims(1)  = int(cpatch%ncohorts,8)
                     memsize(1)  = int(cpatch%ncohorts,8)
                     memoffs(1)  = 0_8

                     call hdf_getslab_r(cpatch%dbh   ,'DBH '   ,dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%hite  ,'HITE '  ,dsetrank,iparallel,.true.)
                     call hdf_getslab_i(cpatch%pft   ,'PFT '   ,dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%nplant,'NPLANT ',dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%bdead ,'BDEAD ' ,dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%balive,'BALIVE ',dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%bleaf ,'BLEAF ' ,dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%bstorage        ,'BSTORAGE '                &
                                       ,dsetrank,iparallel,.true.)
                     call hdf_getslab_i(cpatch%phenology_status,'PHENOLOGY_STATUS '        &
                                       ,dsetrank,iparallel,.true.)

                     !----- First the 2-D variables. --------------------------------------!
                     dsetrank    = 2
                     globdims(1) = 13_8
                     chnkdims(1) = 13_8
                     chnkoffs(1) = 0_8
                     memdims(1)  = 13_8
                     memsize(1)  = 13_8
                     memoffs(2)  = 0_8
                     globdims(2) = int(dset_ncohorts_global,8)
                     chnkdims(2) = int(cpatch%ncohorts,8)
                     chnkoffs(2) = int(paco_id(pa_index) - 1,8)
                     memdims(2)  = int(cpatch%ncohorts,8)
                     memsize(2)  = int(cpatch%ncohorts,8)
                     memoffs(2)  = 0_8

                     call hdf_getslab_r(cpatch%cb    ,'CB '    ,dsetrank,iparallel,.true.)
                     call hdf_getslab_r(cpatch%cb_max,'CB_MAX ',dsetrank,iparallel,.true.)
                     
                     !----- The following variables are initialised with default values. --!
                     cpatch%dagb_dt              = 0.
                     cpatch%dba_dt               = 0.
                     cpatch%ddbh_dt              = 0.
                     cpatch%fsw                  = 1.0
                     cpatch%gpp                  = 0.0
                     cpatch%par_v                = 0.0
                     
                     cohortloop: do ico=1,cpatch%ncohorts
                        !------------------------------------------------------------------!
                        !    We will now check the PFT of each cohort, so we determine     !
                        ! if this is a valid PFT.  If not, then we must decide what we     !
                        ! should do...                                                     !
                        !------------------------------------------------------------------!
                        if (include_pft(cpatch%pft(ico)) == 0) then
                           select case(pft_1st_check)
                           case (0)
                              !----- Stop the run. ----------------------------------------!
                              write (unit=*,fmt='(a,1x,i5,1x,a)')                          &
                                    'I found a cohort with PFT=',cpatch%pft(ico)           &
                                   ,' and it is not in your include_these_pft...'
                              call fatal_error('Invalid PFT in history file'               &
                                              ,'read_ed21_history_file'                    &
                                              ,'ed_read_ed21_history.F90')

                           case (1)
                              !----- Include the unexpected PFT in the list. --------------!
                              write (unit=*,fmt='(a,1x,i5,1x,a)')                          &
                                    'I found a cohort with PFT=',cpatch%pft(ico)           &
                                   ,'... Including this PFT in your include_these_pft...'
                              include_pft(cpatch%pft(ico))        = 1
                              include_these_pft(sum(include_pft)) = cpatch%pft(ico)

                              call sort_up(include_these_pft,n_pft)

                              if (cpatch%pft(ico) == 1 .or. cpatch%pft(ico) == 5) then
                                 include_pft_ag(cpatch%pft(ico)) = 1
                              end if

                           case (2)
                              !----- Ignore the unexpect PFT. -----------------------------!
                              write (unit=*,fmt='(a,1x,i5,1x,a)')                          &
                                    'I found a cohort with PFT=',cpatch%pft(ico)           &
                                   ,'... Ignoring it...'
                              !------------------------------------------------------------!
                              !    The way we will ignore this cohort is by setting its    !
                              ! nplant to zero, and calling the "terminate_cohorts"        !
                              ! subroutine right after this.                               !
                              !------------------------------------------------------------!
                              cpatch%nplant(ico) = 0.
                           end select
                        end if

                        !----- Compute the above-ground biomass. --------------------------!
                        cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)  &
                                                    ,cpatch%bleaf(ico),cpatch%pft(ico)     &
                                                    ,cpatch%hite(ico),cpatch%bstorage(ico))

                        cpatch%basarea(ico)  = cpatch%nplant(ico) * pio4                   &
                                             * cpatch%dbh(ico) * cpatch%dbh(ico)

                        cpatch%broot(ico)    = q(cpatch%pft(ico)) * cpatch%balive(ico)     &
                                             / ( 1.0 + q(cpatch%pft(ico))                  &
                                               + qsw(cpatch%pft(ico)) * cpatch%hite(ico))
                        
                        cpatch%bsapwood(ico) = qsw(cpatch%pft(ico)) * cpatch%balive(ico)   &
                                             * cpatch%hite(ico)                            &
                                             / ( 1.0 + q(cpatch%pft(ico))                  &
                                               + qsw(cpatch%pft(ico)) * cpatch%hite(ico))
                        
                        !----- Assign LAI, WPA, and WAI -----------------------------------!
                        call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)             &
                                         ,cpatch%bdead(ico),cpatch%balive(ico)             &
                                         ,cpatch%dbh(ico), cpatch%hite(ico)                &
                                         ,cpatch%pft(ico), SLA(cpatch%pft(ico))            &
                                         ,cpatch%lai(ico),cpatch%wpa(ico), cpatch%wai(ico))


                        !----- Update the derived patch-level variables. ------------------!
                        csite%lai(ipa)  = csite%lai(ipa) + cpatch%lai(ico)
                        csite%wpa(ipa)  = csite%wpa(ipa) + cpatch%wpa(ico)
                        csite%wai(ipa)  = csite%wai(ipa) + cpatch%wai(ico)
                        csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)          &
                                                    + cpatch%agb(ico)*cpatch%nplant(ico)

                        !----- Initialise the other cohort level variables. ---------------!
                        call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                     end do cohortloop

                     !---------------------------------------------------------------------!
                     !    Eliminate any "unwanted" cohort (i.e., those which nplant was    !
                     ! set to zero so it would be removed).                                !
                     !---------------------------------------------------------------------!
                     call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)

                  end if
               end do patchloop
            else
               !----- This should never happen, but, just in case... ----------------------!
               call fatal_error('A site with no patches was found...'                      &
                               ,'read_ed21_history_file','ed_read_ed21_history.F90')
            end if

            !----- Initialise the other patch-level variables. ----------------------------!
            call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))

         end do siteloop

         !----- Initialise the other site-level variables. --------------------------------!
         call init_ed_site_vars(cpoly,cgrid%lat(ipy))

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
      deallocate(file_lats,file_lons)
      deallocate(paco_n,paco_id)
      deallocate(sipa_n,sipa_id)
      deallocate(pysi_n,pysi_id )
      
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
                             , maxfiles                & ! intent(in)
                             , maxlist                 ! ! intent(in)
   use pft_coms       , only : SLA                     & ! intent(in)
                             , q                       & ! intent(in)
                             , qsw                     & ! intent(in)
                             , hgt_min                 & ! intent(in)
                             , include_pft             & ! intent(in)
                             , include_pft_ag          & ! intent(in)
                             , phenology               & ! intent(in)
                             , pft_1st_check           & ! intent(in)
                             , include_these_pft       ! ! intent(in)
   use ed_misc_coms   , only : sfilin                  & ! intent(in)
                             , current_time            & ! intent(in)
                             , imonthh                 & ! intent(in)
                             , iyearh                  & ! intent(in)
                             , idateh                  & ! intent(in)
                             , itimeh                  & ! intent(in)
                             , ied_init_mode              
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
                             , ed_biomass              ! ! function
   use fuse_fiss_utils, only : terminate_cohorts       ! ! subroutine

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
   integer               , dimension(huge_polygon)              :: pyfile_list
   integer               , dimension(huge_polygon)              :: pyindx_list
   integer               , dimension(:)           , allocatable :: pclosest
   integer               , dimension(:)           , allocatable :: psrcfile
   integer               , dimension(:)           , allocatable :: pysi_n
   integer               , dimension(:)           , allocatable :: pysi_id
   integer               , dimension(:)           , allocatable :: sipa_n
   integer               , dimension(:)           , allocatable :: sipa_id
   integer               , dimension(:)           , allocatable :: paco_n
   integer               , dimension(:)           , allocatable :: paco_id
   integer                                                      :: year
   integer                                                      :: igr
   integer                                                      :: ipy
   integer                                                      :: isi
   integer                                                      :: ipa
   integer                                                      :: ico
   integer                                                      :: nflist
   integer                                                      :: nhisto
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
   integer                                                      :: py_index
   integer                                                      :: si_index
   integer                                                      :: pa_index
   integer                                                      :: dsetrank,iparallel
   integer                                                      :: hdferr
   logical                                                      :: exists
   real                  , dimension(huge_polygon)              :: pdist
   real                  , dimension(huge_polygon)              :: plon_list
   real                  , dimension(huge_polygon)              :: plat_list
   real                                                         :: elim_nplant
   real                                                         :: elim_lai
   real                                                         :: poi_res
   integer                                                      :: poi_minloc
   integer               , dimension(maxfiles)                 :: ngridpoly
   integer                                                      :: ngp1
   integer                                                      :: total_grid_py
   integer                                                      :: ipf
   real                                                         :: plon,plat
   real, dimension(4) :: test_array
   integer :: ipfa, ipfz

   !----- External functions. -------------------------------------------------------------!
   real                  , external                :: dist_gc ! Great circle distance.
   !---------------------------------------------------------------------------------------!

   !----- Retrieve all files with the specified prefix. -----------------------------------!
   call ed_filelist(full_list,sfilin,nflist)

   !----- Check every file and save only those that are actually history files. -----------!
   call ed21_fileinfo(nflist,full_list,nhisto,histo_list)

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
   !     First thing, we go through every file, open, and retrieve only some polygon-level !
   ! information.  This will be used for mapping the files later.                          !
   !---------------------------------------------------------------------------------------!
   pyfile_list(:) = 0
   pyindx_list(:) = 0
   ipfz           = 0
   lonlatloop: do nf=1,nhisto
      hnamel = histo_list(nf)

      !----- Open the HDF5 file. ----------------------------------------------------------!

      call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
      if (hdferr < 0) then
         write(unit=*,fmt='(a,1x,i4)') ' - Error opening HDF5 file - error - ',hdferr
         write(unit=*,fmt='(a,1x,a)') ' - File name: ',trim(hnamel)
         call fatal_error('Error opening HDF5 file - error - '//trim(hnamel)               &
                         ,'read_ed21_history_file','ed_read_ed21_history.F90')
      end if

      !----- Retrieve the number of polygons in this file. --------------------------------!
      call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      !----- Determine the global position of these polygons. -----------------------------!
      ipfa = ipfz + 1
      ipfz = ipfz + dset_npolygons_global

      !----- In case we are meshing grids -------------------------------------------------!
      ngridpoly(nf) = dset_npolygons_global
      
      !------ Here we save the file where we can find these polygons. ---------------------!
      do ipf=ipfa,ipfz
         pyfile_list(ipf) = nf
         pyindx_list(ipf) = ipf-ipfa+1
      end do

      !------------------------------------------------------------------------------------!
      !      Retrieve the polygon coordinates data.                                        !
      !------------------------------------------------------------------------------------!
      globdims(1) = int(dset_npolygons_global,8)

      call h5dopen_f(file_id,'LONGITUDE', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_IEEE_F32LE,plon_list(ipfa:ipfz),globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)

      call h5dopen_f(file_id,'LATITUDE', dset_id, hdferr)
      call h5dget_space_f(dset_id, dspace_id, hdferr)
      call h5dread_f(dset_id, H5T_IEEE_F32LE,plat_list(ipfa:ipfz),globdims, hdferr)
      call h5sclose_f(dspace_id, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5fclose_f(file_id, hdferr)
      if (hdferr /= 0) then
         write (unit=*,fmt='(a,1x,a)') 'File: ',trim(hnamel)
         write (unit=*,fmt='(a)'     ) 'Problem: Failed closing the HDF5 dataset.'
         write (unit=*,fmt='(a,1x,i4)') 'HDFerr: ',hdferr
         call fatal_error('Could not close the HDF file'                                   &
                         ,'read_ed21_history_unstruct','ed_read_ed21_history.F90')
      end if
   end do lonlatloop
   !---------------------------------------------------------------------------------------!

   total_grid_py = ipfz

   !------------------------------------------------------------------------------------!
   ! 
   ! Method for mixing 1 grid and POI's.  Only use the grid if their is NOT an POI
   ! within a user specified resolution
   ! Remember, this assumes there is only 1 grid, and it is the first file
   ! ied_init_mode=99  (Developer use only)
   !------------------------------------------------------------------------------------!
   
   poi_res = 5.0*115.0*1000;

   !---------------------------------------------------------------------------------------!
   !     Big loop over all grids.                                                          !
   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids
      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !     Now, we will go through all polygons, and we will determine which input        !
      ! polygon is the closest one to each polygon.                                        !
      !------------------------------------------------------------------------------------!
      allocate (pclosest(cgrid%npolygons),psrcfile(cgrid%npolygons))

      nneighloop: do ipy=1,cgrid%npolygons
         !----- Reset pdist to a very large number. ---------------------------------------!
         pdist(:)   = 1.e20

         do ipf=1,total_grid_py
            pdist(ipf) = dist_gc(plon_list(ipf),cgrid%lon(ipy)           &
                 ,plat_list(ipf),cgrid%lat(ipy))
         end do

         if(ied_init_mode==5) then
            pclosest(ipy) = pyindx_list(minloc(pdist,dim=1))
            psrcfile(ipy) = pyfile_list(minloc(pdist,dim=1))
         else
            ngp1 = ngridpoly(1)
            pclosest(ipy) = pyindx_list(minloc(pdist(1:ngp1),dim=1))
            psrcfile(ipy) = pyfile_list(1)
            
            poi_minloc = minloc(pdist(ngp1+1:total_grid_py),dim=1)+ngp1
            
            if( pdist(poi_minloc) < poi_res ) then

               pclosest(ipy)  = pyindx_list(poi_minloc)
               psrcfile(ipy)  = pyfile_list(poi_minloc)
            end if

         end if


      end do nneighloop
      

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
!            py_index = pyindx_list(pclosest(ipy))

            py_index = pclosest(ipy)

            !------------------------------------------------------------------------------!
            !      Retrieve the polygon coordinates data.                                        !
            !------------------------------------------------------------------------------!

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
                              ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.)
            
            !----- Load the workload (2D). ------------------------------------------------!
            dsetrank    = 2
            globdims(1) = int(13,8)
            chnkdims(1) = int(13,8)
            memdims(1)  = int(13,8)
            memsize(1)  = int(13,8)
            chnkoffs(1) = 0_8
            memoffs(1)  = 0_8
            globdims(2) = int(pysi_n(py_index),8)
            chnkdims(2) = int(pysi_n(py_index),8)
            memdims(2)  = int(pysi_n(py_index),8)
            memsize(2)  = int(pysi_n(py_index),8)
            chnkoffs(2) = 0_8
            memoffs(2)  = 0_8
            call hdf_getslab_r(cgrid%workload(:,ipy),'WORKLOAD ',dsetrank                  &
                              ,iparallel,.false.)
            !------------------------------------------------------------------------------!

            
            !----- Allocate the vector of sites in the polygon. ---------------------------!
            call allocate_polygontype(cpoly,pysi_n(py_index))     


            !----- Reset the HDF5 auxiliary variables before moving to the next level. ----!
            globdims = 0_8
            chnkdims = 0_8
            chnkoffs = 0_8
            memoffs  = 0_8
            memdims  = 0_8
            memsize  = 1_8

            !------------------------------------------------------------------------------!
            !      SITE level variables.                                                   !
            !------------------------------------------------------------------------------!
            !----- Load 1D dataset. -------------------------------------------------------!
            dsetrank = 1_8
            globdims(1) = int(dset_nsites_global,8)
            chnkdims(1) = int(cpoly%nsites,8)
            chnkoffs(1) = int(pysi_id(py_index) - 1,8)
            memdims(1)  = int(cpoly%nsites,8)
            memsize(1)  = int(cpoly%nsites,8)
            memoffs(1)  = 0_8
            
            call hdf_getslab_r(cpoly%area       ,'AREA_SI '    ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cpoly%moist_f    ,'MOIST_F '    ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cpoly%moist_W    ,'MOIST_W '    ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cpoly%elevation  ,'ELEVATION '  ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cpoly%slope      ,'SLOPE '      ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cpoly%aspect     ,'ASPECT '     ,dsetrank,iparallel,.true.)
            call hdf_getslab_r(cpoly%TCI        ,'TCI '        ,dsetrank,iparallel,.true.)
            call hdf_getslab_i(cpoly%patch_count,'PATCH_COUNT ',dsetrank,iparallel,.true.)

            !----- Load 2D dataset. -------------------------------------------------------!
            dsetrank     = 2_8
            globdims(1)  = int(dset_nzg,8)     ! How many layers in the dataset?
            chnkdims(1)  = int(1,8)            ! We are only extracting one layer
            memdims(1)   = int(1,8)            ! We only need memory for one layer
            memsize(1)   = int(1,8)            ! On both sides
            chnkoffs(1)  = int(dset_nzg - 1,8) ! Take the top layer, not the bottom
            memoffs(1)   = 0_8
            globdims(2)  = int(dset_nsites_global,8)
            chnkdims(2)  = int(cpoly%nsites,8)
            chnkoffs(2)  = int(pysi_id(py_index) - 1,8)
            memdims(2)   = int(cpoly%nsites,8)
            memsize(2)   = int(cpoly%nsites,8)
            memoffs(2)   = 0_8
            call hdf_getslab_i(cpoly%ntext_soil(nzg,:),'NTEXT_SOIL_SI ',dsetrank           &
                              ,iparallel,.true.)

            !------------------------------------------------------------------------------!
            !     Loop over all sites and fill the patch-level variables.                  !
            !------------------------------------------------------------------------------!
            siteloop: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)

               !----- Calculate the index of this site data in the HDF5 file. -------------!
               si_index = pysi_id(py_index) + isi - 1

               if (sipa_n(si_index) > 0) then

                  !------------------------------------------------------------------------!
                  !     The soil layer in this case is use defined, so take this from the  !
                  ! grid level variable, and not from the dataset.                         !
                  !------------------------------------------------------------------------!
                  cpoly%lsl(isi)  = cgrid%lsl(ipy)  ! Initialize lowest soil layer

                  !----- Now fill the soil column based on the top layer data. ------------!
                  do k=1,nzg
                     cpoly%ntext_soil(k,isi) = cpoly%ntext_soil(nzg,isi)
                  end do

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

                  !----- Assign patch soils based off of the site level soils data. -------!
                  do k=1,nzg
                     csite%ntext_soil(k,:) = cpoly%ntext_soil(k,isi)
                  end do

                  call hdf_getslab_i(csite%dist_type         ,'DIST_TYPE '                 &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%age               ,'AGE '                       &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%area              ,'AREA '                      &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%sum_dgd           ,'SUM_DGD '                   &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%sum_chd           ,'SUM_CHD '                   &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_i(csite%plantation        ,'PLANTATION '                &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%fast_soil_C       ,'FAST_SOIL_C '               &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%slow_soil_C       ,'SLOW_SOIL_C '               &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%fast_soil_N       ,'FAST_SOIL_N '               &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%structural_soil_C ,'STRUCTURAL_SOIL_C '         &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%structural_soil_L ,'STRUCTURAL_SOIL_L '         &
                                    ,dsetrank,iparallel,.true.)
                  call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N '        &
                                    ,dsetrank,iparallel,.true.)


                  !------------------------------------------------------------------------!
                  !     Loop over all sites and fill the patch-level variables.            !
                  !------------------------------------------------------------------------!
                  patchloop: do ipa = 1,csite%npatches
                     cpatch => csite%patch(ipa)

                     !---------------------------------------------------------------------!
                     !     Initialise patch-level variables that depend on the cohort      !
                     ! ones.                                                               !
                     !---------------------------------------------------------------------!
                     csite%lai(ipa)               = 0.0
                     csite%wpa(ipa)               = 0.0
                     csite%wai(ipa)               = 0.0
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
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%hite            ,'HITE '                 &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_i(cpatch%pft             ,'PFT '                  &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%nplant          ,'NPLANT '               &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%bdead           ,'BDEAD '                &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%balive          ,'BALIVE '               &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%bleaf           ,'BLEAF '                &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%bstorage        ,'BSTORAGE '             &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_i(cpatch%phenology_status,'PHENOLOGY_STATUS '     &
                                          ,dsetrank,iparallel,.true.)

                        !----- First the 2-D variables. -----------------------------------!
                        dsetrank    = 2
                        globdims(1) = 13_8
                        chnkdims(1) = 13_8
                        chnkoffs(1) = 0_8
                        memdims(1)  = 13_8
                        memsize(1)  = 13_8
                        memoffs(2)  = 0_8
                        globdims(2) = int(dset_ncohorts_global,8)
                        chnkdims(2) = int(cpatch%ncohorts,8)
                        chnkoffs(2) = int(paco_id(pa_index) - 1,8)
                        memdims(2)  = int(cpatch%ncohorts,8)
                        memsize(2)  = int(cpatch%ncohorts,8)
                        memoffs(2)  = 0_8

                        call hdf_getslab_r(cpatch%cb    ,'CB '                             &
                                          ,dsetrank,iparallel,.true.)
                        call hdf_getslab_r(cpatch%cb_max,'CB_MAX '                         &
                                          ,dsetrank,iparallel,.true.)

                        !------------------------------------------------------------------!
                        !    The following variables are initialised with default values.  !
                        !------------------------------------------------------------------!
                        cpatch%dagb_dt              = 0.
                        cpatch%dba_dt               = 0.
                        cpatch%ddbh_dt              = 0.
                        cpatch%fsw                  = 1.0
                        cpatch%gpp                  = 0.0
                        cpatch%par_v                = 0.0

                        cohortloop: do ico=1,cpatch%ncohorts
                           !---------------------------------------------------------------!
                           !    We will now check the PFT of each cohort, so we determine  !
                           ! if this is a valid PFT.  If not, then we must decide what we  !
                           ! should do...                                                  !
                           !---------------------------------------------------------------!
                           if (include_pft(cpatch%pft(ico)) == 0) then
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
                                 include_pft(cpatch%pft(ico))        = 1
                                 include_these_pft(sum(include_pft)) = cpatch%pft(ico)

                                 call sort_up(include_these_pft,n_pft)

                                 if (cpatch%pft(ico) == 1 .or. cpatch%pft(ico) == 5) then
                                    include_pft_ag(cpatch%pft(ico)) = 1
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

                           !----- Compute the above-ground biomass. -----------------------!
                           cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico)                  &
                                                       ,cpatch%balive(ico)                 &
                                                       ,cpatch%bleaf(ico),cpatch%pft(ico)  &
                                                       ,cpatch%hite(ico)                   &
                                                       ,cpatch%bstorage(ico))

                           cpatch%basarea(ico)  = cpatch%nplant(ico) * pio4                &
                                                * cpatch%dbh(ico) * cpatch%dbh(ico)

                           cpatch%broot(ico)    = q(cpatch%pft(ico)) * cpatch%balive(ico)  &
                                                / ( 1.0 + q(cpatch%pft(ico))               &
                                                  + qsw(cpatch%pft(ico))                   &
                                                  * cpatch%hite(ico))
                           
                           cpatch%bsapwood(ico) = qsw(cpatch%pft(ico))                     &
                                                * cpatch%balive(ico) * cpatch%hite(ico)    &
                                                / ( 1.0 + q(cpatch%pft(ico))               &
                                                  + qsw(cpatch%pft(ico))                   &
                                                  * cpatch%hite(ico))

                           !----- Assign LAI, WPA, and WAI --------------------------------!
                           call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)          &
                                            ,cpatch%bdead(ico),cpatch%balive(ico)          &
                                            ,cpatch%dbh(ico),cpatch%hite(ico)              &
                                            ,cpatch%pft(ico),SLA(cpatch%pft(ico))          &
                                            ,cpatch%lai(ico),cpatch%wpa(ico)               &
                                            ,cpatch%wai(ico))


                           !----- Update the derived patch-level variables. ---------------!
                           csite%lai(ipa)  = csite%lai(ipa) + cpatch%lai(ico)
                           csite%wpa(ipa)  = csite%wpa(ipa) + cpatch%wpa(ico)
                           csite%wai(ipa)  = csite%wai(ipa) + cpatch%wai(ico)
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

            end do siteloop

            !----- Initialise the other site-level variables. -----------------------------!
            call init_ed_site_vars(cpoly,cgrid%lat(ipy))

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

         deallocate(paco_n   ,paco_id  )
         deallocate(sipa_n   ,sipa_id  )
         deallocate(pysi_n   ,pysi_id  )
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
