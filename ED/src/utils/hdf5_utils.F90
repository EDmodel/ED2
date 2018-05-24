!==========================================================================================!
!==========================================================================================!
! MODULE: HDF5_UTILS
!
!> \brief Subroutines that handle the HDF5 input and output in ED2
!> \author Developed by Ryan Knox
!> \author 23 May 2018. Deleted stop commands and replaced with fatal error, to make sure
!          MPI processes are cleanly stopped. Also included 'only' in all use hdf5_coms.
!------------------------------------------------------------------------------------------!
module hdf5_utils
#if USE_HDF5
   use hdf5
#endif
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !> \brief This subroutine retrieves the HDF5 information.     
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_info_f(dsetname,ndims,dims)
      use hdf5_coms, only : dsetid_f   & ! intent(inout)
                          , dspaceid_f & ! intent(inout)
                          , fileid_f   ! ! intent(inout)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*)             , intent(in)    :: dsetname !> Dataset name
      integer                      , intent(inout) :: ndims    !> Dataset rank (in file)
      integer, dimension(*)        , intent(inout) :: dims     !> Dataset dimensions
      !----- Local variables. -------------------------------------------------------------!
#if USE_HDF5
      integer(HSIZE_T),dimension(4)                :: dimshf   !> Shuffled dimensions
      integer(HSIZE_T),dimension(4)                :: maxdims  !> Maximum # of dimension
      integer                                      :: hdferr   !> Error flag
      !------------------------------------------------------------------------------------!


      !----- Open the dataset. ------------------------------------------------------------!
      call h5dopen_f(fileid_f,trim(dsetname)//char(0), dsetid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dataset ID:     ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' File ID:        ',fileid_f
         write(unit=*,fmt='(a,1x,i18)') ' Dimension size: ',ndims
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) 'Failed opening dataset, assume dummy dimensions.'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         ndims   = 0
         dims(1) = 0
         return
      end if
      !------------------------------------------------------------------------------------!


      !----- Get data space. --------------------------------------------------------------!
      call h5dget_space_f(dsetid_f, dspaceid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dataset ID:     ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' Data space ID:  ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) 'Failed getting data space, assume dummy dims.'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         ndims   = 0
         dims(1) = 0
         return
      end if
      !------------------------------------------------------------------------------------!


      !----- Get dataset's dimension count. -----------------------------------------------!
      call h5sget_simple_extent_ndims_f(dspaceid_f,ndims,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data space ID:  ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) 'Failed getting dimension count, assume dummy ones.'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         ndims   = 0
         dims(1) = 0
         return
      end if
      !------------------------------------------------------------------------------------!


      !----- Get dataset's extent dimensions. ---------------------------------------------!
      call h5Sget_simple_extent_dims_f (dspaceid_f,dimshf,maxdims,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data space ID:  ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) 'Failed getting dimensions, assume dummy ones.'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         ndims   = 0
         dims(1) = 0
         return
      end if
      !------------------------------------------------------------------------------------!


      !----- Save dimensions.  We use int to convert from kind=8 to kind=4. ---------------!
      dims(1:ndims) = int(dimshf(1:ndims))
      !------------------------------------------------------------------------------------!



      !----- Close dataset. ---------------------------------------------------------------!
      call h5dclose_f(dsetid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dataset ID:     ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Problem closing dataset.','shdf5_info_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!

      call h5sclose_f(dspaceid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data space ID:  ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Problem closing dataset space','shdf5_info_f','hdf5_utils.F90')
      end if
#endif
      return
   end subroutine shdf5_info_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !> \brief   This subroutine opens the HDF5 file.
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_open_f(locfn,access,idelete)

      use hdf5_coms, only : fileid_f
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*), intent(in)           :: locfn   !> File name
      character(len=*), intent(in)           :: access  !> File access ('R','W','RW')
      logical         , intent(in), optional :: idelete !> Delete/overwrite file? (T|F)
      !----- Local variables. -------------------------------------------------------------!
#if USE_HDF5
      integer                                :: hdferr    !> Error flag for HDF5
      logical                                :: exists    !> File existence
      logical                                :: overwrite !> Local delete variable
      !------------------------------------------------------------------------------------!


      ! Default flag for overwriting is false.
      if (present(idelete)) then
         overwrite = idelete
      else
         overwrite = .false.
      end if

      ! Open the HDF environment
      call h5open_f(hdferr)

      ! Collect garbage
      call h5garbage_collect_f(hdferr)


      ! Check for file existence.
      inquire(file=trim(locfn),exist=exists)

      ! Decide whether to read, write, or do both.
      select case (trim(access))
      case ('R','RW')
         !----- Check that file exists before reading. ------------------------------------!
         if (.not. exists) then
            call fatal_error(' File '//trim(locfn)//' not found.','shdf5_open'             &
                            ,'hdf5_utils.f90')
         end if
         !---------------------------------------------------------------------------------!


         !----- Open file with different privileges. --------------------------------------!
         select case (trim(access))
         case ('R')
            call h5fopen_f(trim(locfn)//char(0), H5F_ACC_RDONLY_F, fileid_f, hdferr)
         case ('RW')
            call h5fopen_f(trim(locfn)//char(0), H5F_ACC_RDWR_F, fileid_f, hdferr)
         end select
         !---------------------------------------------------------------------------------!


         !----- Check for errors opening file. --------------------------------------------!
         if (hdferr < 0) then
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) '-------------------------------------------'
            write(unit=*,fmt='(a,1x,a)'  ) ' File:      ',trim(locfn )
            write(unit=*,fmt='(a,1x,i18)') ' HDF Error: ',hdferr
            write(unit=*,fmt='(a)'       ) '-------------------------------------------'
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) ''
            call fatal_error('Error opening file (R/RW).','shdf5_open_f','hdf5_utils.f90')
         end if
         !---------------------------------------------------------------------------------!

      case ('W')
         !---------------------------------------------------------------------------------!
         !     Check whether file exists or not.  In case it exists, check also whether    !
         ! the file can be overwritten.  By default files cannot be overwritten.           !
         !---------------------------------------------------------------------------------!
         if (exists .and. overwrite) then
            !------ Delete file using dirty Fortran trick. --------------------------------!
            open  (unit=99,file=trim(locfn))
            close (unit=99,status='delete')
            !------------------------------------------------------------------------------!
         elseif (exists) then
            !------ Prevent overwrite. ----------------------------------------------------!
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) '-------------------------------------------'
            write(unit=*,fmt='(a,1x,a)'  ) ' File:             ',trim(locfn )
            write(unit=*,fmt='(a,1x,l)'  ) ' Overwrite:        ',overwrite
            write(unit=*,fmt='(a,1x,l)'  ) ' OW Flag provided: ',present(idelete)
            write(unit=*,fmt='(a)'       ) '-------------------------------------------'
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) ''
            call fatal_error('Open existing file for writing with no overwriting.'         &
                            ,'shdf5_open_f','hdf5_utils.f90')
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !----- Open file and check that it opened properly. ------------------------------!
         call h5fcreate_f(trim(locfn)//char(0), H5F_ACC_TRUNC_F, fileid_f, hdferr)
         if (hdferr < 0) then
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) '-------------------------------------------'
            write(unit=*,fmt='(a,1x,a)'  ) ' File:      ',trim(locfn )
            write(unit=*,fmt='(a,1x,i18)') ' File ID:   ',fileid_f
            write(unit=*,fmt='(a,1x,i18)') ' HDF Error: ',hdferr
            write(unit=*,fmt='(a)'       ) '-------------------------------------------'
            write(unit=*,fmt='(a)'       ) ''
            write(unit=*,fmt='(a)'       ) ''
            call fatal_error('Error opening file (W).','shdf5_open_f','hdf5_utils.f90')
         end if
         !---------------------------------------------------------------------------------!
      case default
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' File:   ',trim(locfn )
         write(unit=*,fmt='(a,1x,a)'  ) ' Access: ',trim(access)
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Invalid access: it should be R, RW, or W.'                      &
                         ,'shdf5_open_f','hdf5_utils.f90')

      end select
      return
#endif
   end subroutine shdf5_open_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !> \brief Subroutine for writing variables.  It doesn't seem to be used at all in the
   !         code except for preproc (where it is defined too).  Keeping it for the time 
   !         being but possible candidate for deletion.
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_orec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara               &
                          ,ivars,rvars,cvars,dvars,lvars)
      use hdf5_coms  , only : dsetid_f   & ! intent(inout)
                            , dspaceid_f & ! intent(inout)
                            , fileid_f   & ! intent(inout)
                            , prpid_f    ! ! intent(inout)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      character(len=*)              , intent(in) :: dsetname  !> Variable label
      integer                       , intent(in) :: ndims     !> # of dimensions or rank
      integer, dimension(*)         , intent(in) :: dims      !> Dataset dimensions.
      !----- Optional arguments. Only one should appear in each call. ---------------------!
      integer                       , optional   :: ivars     !> Integer scalar
      real                          , optional   :: rvars     !> Single-precision scalar
      character(len=*)              , optional   :: cvars     !> String scalar
      real(kind=8)                  , optional   :: dvars     !> Double-precision scalar
      logical                       , optional   :: lvars     !> Logical scalar
      integer         , dimension(*), optional   :: ivara     !> Integer array
      real            , dimension(*), optional   :: rvara     !> Single-precision array
      character(len=*), dimension(*), optional   :: cvara     !> String array
      real(kind=8)    , dimension(*), optional   :: dvara     !> Double-precision array
      logical         , dimension(*), optional   :: lvara     !> Logical array
      !----- Local variables. -------------------------------------------------------------!
      integer                                    :: idx       !> Index
#if USE_HDF5
      integer                                    :: h5_type   !> Local type designator
      integer         , dimension(4)             :: dimsh     !> Dataset dimensions.
      integer(HSIZE_T), dimension(4)             :: dimshf    !> Dimensions (shuffle)
      integer(HSIZE_T), dimension(4)             :: chunksize !> Chunk size
      character(len=2)                           :: ctype     !> Variable type.
      integer                                    :: hdferr    !> Error flag
      integer                                    :: npres     !> Presence count.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Check for presence.  Only one presence should be true.                        !
      !------------------------------------------------------------------------------------!
      npres = 0
      if (present(ivars)) then
         ctype = 'is'
         npres = npres + 1
      end if
      if (present(rvars)) then
         ctype = 'rs'
         npres = npres + 1
      end if
      if (present(cvars)) then
         ctype = 'cs'
         npres = npres + 1
      end if
      if (present(dvars)) then
         ctype = 'ds'
         npres = npres + 1
      end if
      if (present(lvars)) then
         ctype = 'ls'
         npres = npres + 1
      end if
      if (present(ivara)) then
         ctype = 'ia'
         npres = npres + 1
      end if
      if (present(rvara)) then
         ctype = 'ra'
         npres = npres + 1
      end if
      if (present(cvara)) then
         ctype = 'ca'
         npres = npres + 1
      end if
      if (present(dvara)) then
         ctype = 'da'
         npres = npres + 1
      end if
      if (present(lvara)) then
         ctype = 'la'
         npres = npres + 1
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Make sure only one argument is present.                                        !
      !------------------------------------------------------------------------------------!
      if (npres /= 1) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name: ',trim(dsetname)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''ivars'' present: ',present(ivars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''rvars'' present: ',present(rvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''cvars'' present: ',present(cvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''dvars'' present: ',present(dvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''lvars'' present: ',present(lvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''ivara'' present: ',present(ivara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''rvara'' present: ',present(rvara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''cvara'' present: ',present(cvara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''dvara'' present: ',present(dvara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''lvara'' present: ',present(lvara)
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Invalid settins, one and only one of the above variables '      &
                         ,'should be present!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Make sure dimensions are reasonable.                                          !
      !------------------------------------------------------------------------------------!
      if (ndims <=0 .or. any(dims(1:ndims) <= 0)) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name: ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Number of dimensions:   ',ndims
         do idx=1,ndims
            write(unit=*,fmt='(2(a,1x,i6,1x))') ' Dimension ',idx,'; size: ',dims(idx)
         end do
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Dataset has bad dimensions!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Set dimensions. --------------------------------------------------------------!
      dimsh(1:ndims)     = dims(1:ndims)
      dimshf             = 0_8
      chunksize          = 0_8
      dimshf(1:ndims)    = int(dimsh(1:ndims),8)
      chunksize(1:ndims) = int(dimsh(1:ndims),8)
      !------------------------------------------------------------------------------------!


      !----- Create data space. -----------------------------------------------------------!
      call h5screate_simple_f(ndims,dimshf,dspaceid_f,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dimension size: ',ndims
         do idx=1,ndims
            write(unit=*,fmt='(2(a,1x,i6,1x))') ' Dimension ',idx,'; size: ',dims(idx)
         end do
         write(unit=*,fmt='(a,1x,i18)') ' Dspaceid_f:     ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed creating datasdpace!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Create property. -------------------------------------------------------------!
      call h5pcreate_f(H5P_DATASET_CREATE_F, prpid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Prpid_f: '       ,prpid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed creating property!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Set chunk sizes. -------------------------------------------------------------!
      call h5pset_chunk_f(prpid_f, ndims, chunksize, hdferr) 
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Prpid_f: '       ,prpid_f
         write(unit=*,fmt='(a,1x,i18)') ' Dimension size: ',ndims
         do idx=1,ndims
            write(unit=*,fmt='(2(a,1x,i6,1x))') ' Chunksize ',idx,'; size: ',chunksize(idx)
         end do
         write(unit=*,fmt='(a,1x,i18)') ' Dspaceid_f:     ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed setting chunk sizes!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Set shuffle. -----------------------------------------------------------------!
      call h5pset_shuffle_f(prpid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Prpid_f: '       ,prpid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed setting shuffle!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Set deflate. -----------------------------------------------------------------!
      call h5pset_deflate_f(prpid_f, 9, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Prpid_f: '       ,prpid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed setting deflate!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Define hdf5 type. ------------------------------------------------------------!
      select case (ctype)
      case ('is','ia')
         h5_type = 1
      case ('rs','ra')
         h5_type = 2
      case ('cs','ca')
         h5_type = 3
      case ('ds','da')
         h5_type = 4
      case ('ls','la')
         h5_type = 5
      end select
      !------------------------------------------------------------------------------------!


      !---- Write the dataset. ------------------------------------------------------------!
      select case (trim(ctype))
      case ('is') ! Integer, scalar
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_INTEGER             &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_INTEGER,ivars,dimshf,hdferr)

      case ('rs') ! Single-precision, scalar
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_REAL                &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_REAL,rvars,dimshf,hdferr)

      case ('cs') ! String, scalar
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_CHARACTER           &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_CHARACTER,cvars,dimshf,hdferr)

      case ('ds') ! Double-precision, scalar
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_DOUBLE              &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_DOUBLE,dvars,dimshf,hdferr)

      case ('ia') ! Integer, array
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_INTEGER             &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_INTEGER,ivara,dimshf,hdferr)

      case ('ra') ! Real, array
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_REAL                &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_REAL,rvara,dimshf,hdferr)

      case ('ca') ! String, array
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_CHARACTER           &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_CHARACTER,cvara,dimshf,hdferr)

      case ('da') ! Double-precision, array
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_DOUBLE              &
                         ,dspaceid_f,dsetid_f, hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_DOUBLE,dvara,dimshf,hdferr)

      case ('ls','la')
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a)'       ) ' There is no HDF5 Fortran API for logical vars.'
         write(unit=*,fmt='(a)'       ) ' You must change back to C or switch type.'
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Attempt to write logical variable!','shdf5_orec_f'              &
                         ,'hdf5_utils.F90')
      end select
      !------------------------------------------------------------------------------------!



      !----- Make sure everything went fine. ----------------------------------------------!
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed writing variable!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!




      !----- Close the dataspace. ---------------------------------------------------------!
      call h5sclose_f(dspaceid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dspaceid_f:     ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed closing data space!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Close the dataset. -----------------------------------------------------------!
      call h5dclose_f(dsetid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dsetid_f:       ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed closing data space!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!

      !----- Close the property. ----------------------------------------------------------!
      call h5pclose_f(prpid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Prpid_f:        ',prpid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed closing property!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!
#endif
      return
   end subroutine shdf5_orec_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !> \brief Subroutine for reading variables.
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_irec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,ivars,rvars   &
                          ,cvars,dvars,lvars)
      use hdf5_coms  , only : dsetid_f   & ! intent(inout)
                            , dspaceid_f & ! intent(inout)
                            , fileid_f   ! ! intent(inout)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      character(len=*)              , intent(in)  :: dsetname  !> Variable label
      integer                       , intent(in)  :: ndims     !> # of dimensions or rank
      integer, dimension(*)         , intent(in)  :: dims      !> Dataset dimensions.
      !----- Optional arguments. Only one should appear in each call. ---------------------!
      integer                       , optional    :: ivars     !> Integer scalar
      real                          , optional    :: rvars     !> Single-precision scalar
      character(len=*)              , optional    :: cvars     !> String scalar
      real(kind=8)                  , optional    :: dvars     !> Double-precision scalar
      logical                       , optional    :: lvars     !> Logical scalar
      integer         , dimension(*), optional    :: ivara     !> Integer array
      real            , dimension(*), optional    :: rvara     !> Single-precision array
      character(len=*), dimension(*), optional    :: cvara     !> String array
      real(kind=8)    , dimension(*), optional    :: dvara     !> Double-precision array
      logical         , dimension(*), optional    :: lvara     !> Logical array
      !----- Local variables. -------------------------------------------------------------!
      integer                                     :: idx       !> Index
#if USE_HDF5
      integer         , dimension(4)              :: dimsh     !> Dataset dimensions.
      integer(HSIZE_T), dimension(4)              :: dimshf    !> Dimensions (shuffle)
      integer(HID_T)                              :: type_id   !> Variable type ID.
      character(len=2)                            :: ctype     !> Variable type.
      integer                                     :: hdferr    !> Error flag
      integer                                     :: npres     !> Presence count.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Check for presence.  Only one presence should be true.                        !
      !------------------------------------------------------------------------------------!
      npres = 0
      if (present(ivars)) then
         ctype = 'is'
         npres = npres + 1
      end if
      if (present(rvars)) then
         ctype = 'rs'
         npres = npres + 1
      end if
      if (present(cvars)) then
         ctype = 'cs'
         npres = npres + 1
      end if
      if (present(dvars)) then
         ctype = 'ds'
         npres = npres + 1
      end if
      if (present(lvars)) then
         ctype = 'ls'
         npres = npres + 1
      end if
      if (present(ivara)) then
         ctype = 'ia'
         npres = npres + 1
      end if
      if (present(rvara)) then
         ctype = 'ra'
         npres = npres + 1
      end if
      if (present(cvara)) then
         ctype = 'ca'
         npres = npres + 1
      end if
      if (present(dvara)) then
         ctype = 'da'
         npres = npres + 1
      end if
      if (present(lvara)) then
         ctype = 'la'
         npres = npres + 1
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Make sure only one argument is present.                                        !
      !------------------------------------------------------------------------------------!
      if (npres /= 1) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name: ',trim(dsetname)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''ivars'' present: ',present(ivars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''rvars'' present: ',present(rvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''cvars'' present: ',present(cvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''dvars'' present: ',present(dvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''lvars'' present: ',present(lvars)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''ivara'' present: ',present(ivara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''rvara'' present: ',present(rvara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''cvara'' present: ',present(cvara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''dvara'' present: ',present(dvara)
         write(unit=*,fmt='(a,1x,l)'  ) ' Variable ''lvara'' present: ',present(lvara)
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Invalid settins, one and only one of the above variables '      &
                         ,'should be present!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Make sure dimensions are reasonable.                                          !
      !------------------------------------------------------------------------------------!
      if (ndims <=0 .or. any(dims(1:ndims) <= 0)) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:           ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Number of dimensions:   ',ndims
         do idx=1,ndims
            write(unit=*,fmt='(2(a,1x,i6,1x))') ' Dimension ',idx,'; size: ',dims(idx)
         end do
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Dataset has bad dimensions!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Set dimensions. --------------------------------------------------------------!
      dimsh(1:ndims)     = dims(1:ndims)
      dimshf             = 0_8
      dimshf(1:ndims)    = int(dimsh(1:ndims),8)
      !------------------------------------------------------------------------------------!


      !----- Open dataset. ----------------------------------------------------------------!
      call h5dopen_f(fileid_f,trim(dsetname)//char(0), dsetid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data set ID:    ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed opening dataset!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Retrieve data space. ---------------------------------------------------------!
      call h5dget_space_f(dsetid_f, dspaceid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data set ID:    ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' Data space ID:  ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed retrieving data space!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Check that data field argument matches data type, convert if possible. -------!
      call h5dget_type_f(dsetid_f,type_id,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data set ID:    ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed checking data type!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Read the data.                                                                !
      !------------------------------------------------------------------------------------!
      select case (ctype)
      case ('is') ! Integer, scalar
         call h5dread_f(dsetid_f, H5T_NATIVE_INTEGER ,ivars, dimshf, hdferr)

      case ('rs') ! Single-precision, scalar
         call h5dread_f(dsetid_f, H5T_NATIVE_REAL    ,rvars, dimshf, hdferr)

      case ('cs') ! Character, scalar
         call h5dread_f(dsetid_f,H5T_NATIVE_CHARACTER,cvars, dimshf, hdferr)

      case ('ds') ! Double-precision, scalar
         call h5dread_f(dsetid_f,H5T_NATIVE_DOUBLE   ,dvars, dimshf, hdferr)

      case ('ia') ! Integer, array
         call h5dread_f(dsetid_f,H5T_NATIVE_INTEGER  ,ivara, dimshf, hdferr)

      case ('ra') ! Single-precision, array
         call h5dread_f(dsetid_f,H5T_NATIVE_REAL     ,rvara, dimshf, hdferr)

      case ('ca') ! Character, array
         call h5dread_f(dsetid_f,H5T_NATIVE_CHARACTER,cvara, dimshf, hdferr)

      case ('da') ! Double-precision, array
         call h5dread_f(dsetid_f,H5T_NATIVE_DOUBLE   ,dvara, dimshf, hdferr)
      case ('ls','la')
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:   ',trim(dsetname)
         write(unit=*,fmt='(a)'       ) ' There is no HDF5 Fortran API for logical vars.'
         write(unit=*,fmt='(a)'       ) ' You must change back to C or switch type.'
         write(unit=*,fmt='(a)'       ) '-------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Attempt to read in logical variable!','shdf5_irec'              &
                         ,'hdf5_utils.F90')
      end select
      !------------------------------------------------------------------------------------!



      !----- Make sure everything went fine. ----------------------------------------------!
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:       ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data type in file:  ',type_id
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag:     ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed reading variable!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!


      !----- Close data space. ------------------------------------------------------------!
      call h5sclose_f(dspaceid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:       ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Data space ID:      ',dspaceid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag:     ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed closing data space!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!



      !----- Close dataset. ===------------------------------------------------------------!
      call h5dclose_f(dsetid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a,1x,a)'  ) ' Dataset name:       ',trim(dsetname)
         write(unit=*,fmt='(a,1x,i18)') ' Dataset ID:         ',dsetid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF error flag:     ',hdferr
         write(unit=*,fmt='(a)'       ) '------------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Failed closing dataset!','shdf5_irec_f','hdf5_utils.F90')
      end if
      !------------------------------------------------------------------------------------!

#endif
      return
   end subroutine shdf5_irec_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !> \brief Subroutine that closes HD5 file.
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_close_f()
      use hdf5_coms, only : fileid_f ! ! intent(in)
      implicit none
      !----- Local variables. -------------------------------------------------------------!
#if USE_HDF5
      integer     :: hdferr  !> Error flags
      !------------------------------------------------------------------------------------!

      !----- Close the hdf file. ----------------------------------------------------------!
      call h5fclose_f(fileid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a,1x,i18)') ' File ID:   ',fileid_f
         write(unit=*,fmt='(a,1x,i18)') ' HDF Error: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Error closing HDF5 file.','shdf5_close_f','hdf5_utils.f90')
      end if
      !------------------------------------------------------------------------------------!

      !----- Collect garbage. -------------------------------------------------------------!
      call h5garbage_collect_f(hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a,1x,i18)') ' HDF Error: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Error collecting garbage.','shdf5_close_f','hdf5_utils.f90')
      end if
      !------------------------------------------------------------------------------------!

      !----- Close the hdf environment. ---------------------------------------------------!
      call h5close_f(hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a,1x,i18)') ' HDF Error: ',hdferr
         write(unit=*,fmt='(a)'       ) '-------------------------------------------'
         write(unit=*,fmt='(a)'       ) ''
         write(unit=*,fmt='(a)'       ) ''
         call fatal_error('Error closing the HDF5 environment.'                            &
                         ,'shdf5_close_f','hdf5_utils.f90')
      end if
      !------------------------------------------------------------------------------------!
#endif
      return
   end subroutine shdf5_close_f
   !=======================================================================================!
   !=======================================================================================!
end module hdf5_utils
!==========================================================================================!
!==========================================================================================!
