!==========================================================================================!
!==========================================================================================!
!    Module hdf5_utils.  This module contains the interface with the HDF5 library.         !
!------------------------------------------------------------------------------------------!
module hdf5_utils
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will initialise the HDF5 environment.                             !
   !---------------------------------------------------------------------------------------!
   subroutine init_hdf5_env()
#if USE_HDF5
      use hdf5
      use hdf5_coms
      implicit none
      !----- Local variables. -------------------------------------------------------------!
      integer :: hdferr
      !------------------------------------------------------------------------------------!

      call h5open_f(hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') ' => HDF5 Open error #:',hdferr
         call fatal_error('Could not initialize the hdf environment'                       &
                         ,'init_hdf5_env','hdf5_utils.F90')
      end if
#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'init_hdf5_env','hdf5_utils.F90')
#endif
      return
   end subroutine init_hdf5_env
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine retrieves the file information such as dimensions.                !
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_info_f(dsetname,ndims,dims)
#if USE_HDF5
      use hdf5
      use hdf5_coms
#endif
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*)     , intent(in)    :: dsetname ! Dataset name
      integer              , intent(out)   :: ndims    ! Dataset rank (in file)
      integer, dimension(*), intent(inout) :: dims     ! Dimensions.
#if USE_HDF5
      !----- Local variables. -------------------------------------------------------------!
      integer                              :: hdferr   ! Error flag
      integer(HSIZE_T),dimension(4)        :: dimshf   ! 
      integer(HSIZE_T),dimension(4)        :: maxdims
      !------------------------------------------------------------------------------------!


      !----- Open the dataset. ------------------------------------------------------------!
      call h5dopen_f(fileid_f,trim(dsetname)//char(0), dsetid_f, hdferr)
      call h5dget_space_f(dsetid_f, dspaceid_f, hdferr)

      if (hdferr < 0) then
         write(unit=*,fmt='(a)') '=== In shdf5_info:'
         write(unit=*,fmt='(3a)') '  Variable ', trim(dsetname)                            &
                                , ' is not in the currently opened hdf5 file'
         ndims   = 0
         dims(1) = 0
         return
      end if

      !----- Get dataset's dimensions. ----------------------------------------------------!
      call h5sget_simple_extent_ndims_f(dspaceid_f, ndims, hdferr) 
      call h5sget_simple_extent_dims_f(dspaceid_f, dimshf, maxdims, hdferr  )

      dims(1:ndims) = dimshf(1:ndims)

      call h5dclose_f(dsetid_f, hdferr)
      if (hdferr /= 0) then
         call fatal_error('Could not close the dataset!','shdf5_info_f','hdf5_utils.F90')
      end if

      call h5sclose_f(dspaceid_f, hdferr)
      if (hdferr /= 0) then
         call fatal_error('Could not close the dataset!','shdf5_info_f','hdf5_utils.F90')
      end if
#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'shdf5_info_f','hdf5_utils.F90')
#endif

      return
   end subroutine shdf5_info_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This routine calls shdf5_irec or shdf5_orec to read or write a variable          !
   ! depending on whether 'action' equals 'read' or 'write'.                               !
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_io(action,ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                                 ,ivars,rvars,cvars,dvars,lvars)
#if USE_HDF5
      use hdf5
      use hdf5_coms
#endif
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*)              , intent(in)   :: dsetname
      character(len=*)              , intent(in)   :: action
      integer                       , intent(in)   :: ndims
      integer         , dimension(*), intent(in)   :: dims
      integer                       , intent(inout), optional :: ivars
      integer         , dimension(*), intent(inout), optional :: ivara
      real                          , intent(inout), optional :: rvars
      real            , dimension(*), intent(inout), optional :: rvara
      real(kind=8)                  , intent(inout), optional :: dvars
      real(kind=8)    , dimension(*), intent(inout), optional :: dvara
      character(len=*)              , intent(inout), optional :: cvars
      character(len=*), dimension(*), intent(inout), optional :: cvara
      logical                       , intent(inout), optional :: lvars
      logical         , dimension(*), intent(inout), optional :: lvara
      !------------------------------------------------------------------------------------!

#if USE_HDF5

      select case (trim(action))
      case ('READ')
         call shdf5_irec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara               &
                          ,ivars,rvars,cvars,dvars,lvars)
         
      case ('WRITE')
         call shdf5_orec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara               &
              ,ivars,rvars,cvars,dvars,lvars)

      case default
         call fatal_error('Illegal action in shdf5_io. Action should be READ or WRITE'     &
                         ,'shdf5_io','hdf5_utils.F90')

      end select
#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'shdf5_io','hdf5_utils.F90')
#endif

      return
   end subroutine shdf5_io
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will open a HDF5 file.                                            !
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_open_f(locfn,access,idelete)
#if USE_HDF5
      use hdf5
      use hdf5_coms
#endif
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*), intent(in)           :: locfn   ! file name
      character(len=*), intent(in)           :: access  ! File access ('R','W','RW')
      integer         , intent(in), optional :: idelete ! If W, delete/overwrite file if 
                                                        !    exists? 1=yes, 0=no
                                                        ! If access is R, this is ignored.
#if USE_HDF5
      !----- Local variables. -------------------------------------------------------------!
      integer          :: hdferr  ! Error flag for HDF5
      character(len=2) :: caccess ! File access ('R ','W ','RW')
      logical          :: exists  ! File existence
      !------------------------------------------------------------------------------------!

      caccess = access

      !----- Open the HDF environment. ----------------------------------------------------!
      call h5open_f(hdferr)

      !----- Collect garbage. -------------------------------------------------------------!
      call h5garbage_collect_f(hdferr)

      inquire(file=trim(locfn),exist=exists)
      
      !----- Create a new file or open an existing HDF5 file. -----------------------------!
      select case (access(1:1))
      case ('R')
         if (.not.exists) then
            call fatal_error('File '//trim(locfn)//' not found.'                           &
                            ,'shdf5_open','hdf5_utils.f90')
         else
            select case(caccess)
            case ('R')
               call h5fopen_f(trim(locfn)//char(0), H5F_ACC_RDONLY_F, fileid_f, hdferr)
            case ('RW')
               call h5fopen_f(trim(locfn)//char(0), H5F_ACC_RDWR_F, fileid_f, hdferr)
            case default
               call fatal_error('Invalid file read type access. Should be R or RW' &
                               ,'shdf5_open','hdf5_utils.f90')
            end select
            
            if (hdferr < 0) then
               write (unit=*,fmt='(a,1x,i5)') 'Error opening hdf5 file: Error =',hdferr
               call fatal_error('Error opening file '//trim(locfn)                         &
                               ,'shdf5_open','hdf5_utils.f90')
            end if
         end if

      case ('W')
         if (.not.exists) then

            call h5fcreate_f(trim(locfn)//char(0), H5F_ACC_EXCL_F, fileid_f, hdferr)
         else
            if (.not.present(idelete)) then
               call fatal_error('idelete not specified when access=W'                      &
                               ,'shdf5_open','hdf5_utils.f90')
            end if
            
            if (idelete == 0) then
               write (unit=*,fmt='(a)')    'In shdf5_open:'
               write (unit=*,fmt='(a)')    'Attempt to open an existing file for writing, '
               write (unit=*,fmt='(a,1x,i5)') 'but overwrite is disabled. idelete=',idelete
               call fatal_error('Open existing file for writing with no overwriting.'      &
                               ,'shdf5_open','hdf5_utils.f90')
            else
               !----- Avoiding system calls... --------------------------------------------!
               open  (unit=99,file=trim(locfn))
               close (unit=99,status='delete')
               call h5fcreate_f(trim(locfn)//char(0), H5F_ACC_TRUNC_F, fileid_f, hdferr)
            end if
         end if

         if(hdferr < 0) then
            write(unit=*,fmt=*) 'File name:',trim(locfn),' ',trim(access), idelete,hdferr
            call fatal_error('HDF5 file '//trim(locfn)//' create failed:'                  &
                            ,'shdf5_open','hdf5_utils.f90')
         end if
      end select
#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'shdf5_open_f','hdf5_utils.F90')

#endif

      return
   end  subroutine shdf5_open_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will write some data into a HDF5 file.                            !
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_orec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,ivars,rvars   &
                          ,cvars,dvars,lvars)
#if USE_HDF5
      use hdf5
      use hdf5_coms
#endif
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*)              , intent(in)   :: dsetname
      integer                       , intent(in)   :: ndims
      integer         , dimension(*), intent(in)   :: dims
      integer                       , intent(inout), optional :: ivars
      integer         , dimension(*), intent(inout), optional :: ivara
      real                          , intent(inout), optional :: rvars
      real            , dimension(*), intent(inout), optional :: rvara
      real(kind=8)                  , intent(inout), optional :: dvars
      real(kind=8)    , dimension(*), intent(inout), optional :: dvara
      character(len=*)              , intent(inout), optional :: cvars
      character(len=*), dimension(*), intent(inout), optional :: cvara
      logical                       , intent(inout), optional :: lvars
      logical         , dimension(*), intent(inout), optional :: lvara
#if USE_HDF5
      !----- Local variables. -------------------------------------------------------------!
      character(len=2)               :: ctype     ! Variable type: int, real, char
      integer                        :: d
      integer                        :: h5_type   ! Local type designator
      integer         , dimension(4) :: dimsh     ! Dataset dimensions.
      integer                        :: hdferr    ! Error flag
      integer(HSIZE_T), dimension(4) :: dimshf
      integer(HSIZE_T), dimension(4) :: chunksize
      !------------------------------------------------------------------------------------!

      
      
      !----- Find which data type is input. -----------------------------------------------!
      if (present(ivars)) then
         ctype='is'
      elseif (present(rvars)) then
         ctype='rs'
      elseif (present(cvars)) then
         ctype='cs'
      elseif (present(dvars)) then
         ctype='ds'
      elseif (present(lvars)) then
         ctype='ls'
      elseif (present(ivara)) then
         ctype='ia'
      elseif (present(rvara)) then
         ctype='ra'
      elseif (present(cvara)) then
         ctype='ca'
      elseif (present(dvara)) then
         ctype='da'
      elseif (present(lvara)) then
         ctype='la'
      else
         call fatal_error('Incorrect or missing data field argument in shdf5_orec.'        &
                         ,'shdf5_orec_f','hdf5_utils.F90')
      end if
      
      !----- Check dimensions and set compression chunk size. -----------------------------!
      if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
         write(unit=*,fmt='(a)') 'Dimension error in shdf5_orec:'
         write(unit=*,fmt='(a,1x,i5)') 'NDIMS =',ndims
         do d=1,ndims
            write(unit=*,fmt='(a,i2,a,1x,i5)') 'DIM(',d,') =',dims(d)
         end do
         call fatal_error('Bad dimensions.','shdf5_orec_f','hdf5_utils.F90')
      end if
      
      dimsh(1:ndims)     = dims(1:ndims)
      dimshf             = 0_8
      chunksize          = 0_8
      dimshf(1:ndims)    = int(dimsh(1:ndims),8)
      chunksize(1:ndims) = int(dimsh(1:ndims),8)
      
      !----- Prepare memory and options for the write. ------------------------------------!
      call h5screate_simple_f(ndims,dimshf,dspaceid_f,hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'NDIMS      =',ndims
         write(unit=*,fmt='(a,1x,i5)') 'DIMSHF     =',dimshf
         write(unit=*,fmt='(a,1x,i5)') 'DSPACEID_F =',dspaceid_f
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while creating dataspace!'                               &
                         ,'shdf5_orec_f','hdf5_utils.F90')
      end if

      call h5pcreate_f(H5P_DATASET_CREATE_F, prpid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'PRPID_F    =',prpid_f
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while creating property!'                                &
                         ,'shdf5_orec_f','hdf5_utils.F90')
      end if
      
      call h5pset_chunk_f(prpid_f, ndims, chunksize, hdferr) 
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'NDIMS      =',ndims
         write(unit=*,fmt='(a,1x,i5)') 'CHUNKSIZE  =',chunksize
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while setting chunk!','shdf5_orec_f','hdf5_utils.F90')
      end if
      
      call h5pset_shuffle_f(prpid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while setting shuffle!','shdf5_orec_f','hdf5_utils.F90')
      end if
      
      call h5pset_deflate_f(prpid_f, 9, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while setting deflate!','shdf5_orec_f','hdf5_utils.F90')
      end if
      
      !----- Select HDF type. -------------------------------------------------------------!
      select case(ctype(1:1))
      case ('i') 
         h5_type=1
      case ('r')
         h5_type=2
      case ('c')
         h5_type=3
      case ('d')
         h5_type=4  ! If native precision is 8 bytes, do h5_type=2
      case ('l')
         h5_type=5
      end select
      
      !----- Write the dataset. -----------------------------------------------------------!
      select case (ctype)
      case ('is')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_INTEGER,dspaceid_f   &
                         ,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_INTEGER,ivars,dimshf,hdferr)


      case ('rs')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_REAL,dspaceid_f      &
                         ,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_REAL,rvars,dimshf,hdferr)


      case ('cs')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_CHARACTER            &
                         ,dspaceid_f,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_CHARACTER,cvars,dimshf,hdferr)


      case ('ds')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_DOUBLE,dspaceid_f    &
                         ,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_DOUBLE,dvars,dimshf,hdferr)


      case ('ia')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_INTEGER,dspaceid_f   &
                         ,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_INTEGER,ivara,dimshf,hdferr)

      case ('ra')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_REAL,dspaceid_f      &
                         ,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_REAL,rvara,dimshf,hdferr)
         
      case ('ca')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_CHARACTER            &
                         ,dspaceid_f,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_CHARACTER,cvara,dimshf,hdferr)
         
      case ('da')
         call h5dcreate_f(fileid_f,trim(dsetname)//char(0),H5T_NATIVE_DOUBLE               &
                         ,dspaceid_f,dsetid_f,hdferr)
         call h5dwrite_f(dsetid_f,H5T_NATIVE_DOUBLE,dvara,dimshf,hdferr)
         
      case ('ls','la')
         write (unit=*,fmt='(a)') 'There is no HDF5 fortran API datatype for Boolean...'
         call fatal_error('Invalid variable type!','shdf5_orec_f','hdf5_utils.F90')

      end select

      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while writing to file!','shdf5_orec_f','hdf5_utils.F90')
      end if

      !------------------------------------------------------------------------------------!
      !     Close the dataset, the dataspace for the dataset, and the dataspace            !
      ! properties.                                                                        !
      !------------------------------------------------------------------------------------!
      !----- Close out the dataspace. -----------------------------------------------------!
      call h5sclose_f(dspaceid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while closing dataspace!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !----- Close out the dataset. -------------------------------------------------------!
      call h5dclose_f(dsetid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while closing dataset!','shdf5_orec_f','hdf5_utils.F90')
      end if
      !----- Close the property. ----------------------------------------------------------!
      call h5pclose_f(prpid_f, hdferr)
      if (hdferr /= 0) then
         write(unit=*,fmt='(a,1x,i5)') 'HDFERR     =',hdferr
         call fatal_error('Failed while closing property!','shdf5_orec_f','hdf5_utils.F90')
      end if

#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'shdf5_orec_f','hdf5_utils.F90')

#endif


      return
   end subroutine shdf5_orec_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will read some data from a HDF5 file.                             !
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_irec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,ivars,rvars   &
                          ,cvars,dvars,lvars)
#if USE_HDF5
      use hdf5
      use hdf5_coms
#endif
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      character(len=*)              , intent(in)   :: dsetname
      integer                       , intent(in)   :: ndims
      integer         , dimension(*), intent(in)   :: dims
      integer                       , intent(inout), optional :: ivars
      integer         , dimension(*), intent(inout), optional :: ivara
      real                          , intent(inout), optional :: rvars
      real            , dimension(*), intent(inout), optional :: rvara
      real(kind=8)                  , intent(inout), optional :: dvars
      real(kind=8)    , dimension(*), intent(inout), optional :: dvara
      character(len=*)              , intent(inout), optional :: cvars
      character(len=*), dimension(*), intent(inout), optional :: cvara
      logical                       , intent(inout), optional :: lvars
      logical         , dimension(*), intent(inout), optional :: lvara
#if USE_HDF5
      !----- Local variables. -------------------------------------------------------------!
      character(len=2)                            :: ctype     ! Variable type
      integer                                     :: d
      integer                                     :: type_id
      integer         , dimension(4)              :: dimsh     ! Dataset dimensions.
      integer                                     :: hdferr    ! Error flag
      integer(HSIZE_T), dimension(4)              :: dimshf
      !------------------------------------------------------------------------------------!
      
      
      
      
      !----- Find which data type will be read. -------------------------------------------!
      if (present(ivars)) then
         ctype='is'
      elseif (present(rvars)) then
         ctype='rs'
      elseif (present(cvars)) then
         ctype='cs'
      elseif (present(dvars)) then
         ctype='ds'
      elseif (present(lvars)) then
         ctype='ls'
      elseif (present(ivara)) then
         ctype='ia'
      elseif (present(rvara)) then
         ctype='ra'
      elseif (present(cvara)) then
         ctype='ca'
      elseif (present(dvara)) then
         ctype='da'
      elseif (present(lvara)) then
         ctype='la'
      else
         call fatal_error('Incorrect or missing data field argument in shdf5_irec.'        &
                         ,'shdf5_irec_f','hdf5_utils.F90')
      end if
      
      !----- Check dimensions and set compression chunk size. -----------------------------!
      if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
         write(unit=*,fmt='(a)') 'Dimension error in shdf5_orec:'
         write(unit=*,fmt='(a,1x,i5)') 'NDIMS =',ndims
         do d=1,ndims
            write(unit=*,fmt='(a,i2,a,1x,i5)') 'DIM(',d,') =',dims(d)
         end do
         call fatal_error('Bad dimensions.','shdf5_irec_f','hdf5_utils.F90')
      end if
      
      dimsh(1:ndims)  = dims(1:ndims)
      dimshf          = 0_8
      dimshf(1:ndims) = int(dimsh(1:ndims),8)

      !----- Open file. -------------------------------------------------------------------!
      call h5dopen_f(fileid_f,trim(dsetname)//char(0), dsetid_f, hdferr)
      call h5dget_space_f(dsetid_f, dspaceid_f, hdferr)

      !----- Check that data field arguement matches data type, convert if possible... ----!
      call h5dget_type_f(dsetid_f,type_id,hdferr)

      !----- Read the data. ---------------------------------------------------------------!
      select case (ctype)
      case ('is')
         call h5dread_f(dsetid_f, H5T_NATIVE_INTEGER, ivars, dimshf, hdferr)

      case ('rs') 
         call h5dread_f(dsetid_f, H5T_NATIVE_REAL,rvars, dimshf, hdferr )

      case ('cs') 
         call h5dread_f(dsetid_f,H5T_NATIVE_CHARACTER,cvars, dimshf, hdferr )

      case ('ds') 
         call h5dread_f(dsetid_f,H5T_NATIVE_DOUBLE,dvars, dimshf, hdferr )

      case ('ia') 
         call h5dread_f(dsetid_f,H5T_NATIVE_INTEGER,ivara, dimshf, hdferr )

      case ('ra') 
         call h5dread_f(dsetid_f,H5T_NATIVE_REAL,rvara, dimshf, hdferr )

      case ('ca') 
         call h5dread_f(dsetid_f,H5T_NATIVE_CHARACTER,cvara, dimshf, hdferr )

      case ('da') 
         call h5dread_f(dsetid_f,H5T_NATIVE_DOUBLE,dvara, dimshf, hdferr )

      case ('ls','la') 
         write (unit=*,fmt='(a)') 'There is no HDF5 fortran API datatype for Boolean...'
         call fatal_error('Invalid variable type!','shdf5_irec_f','hdf5_utils.F90')

      case default
         call fatal_error('Invalid ctype '//trim(ctype)//'!!!'                             &
                         ,'shdf5_irec_f','hdf5_utils.f90')

      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Close the dataset and the dataspace for the dataset.                           !
      !------------------------------------------------------------------------------------!
      !----- Close out the dataspace. -----------------------------------------------------!
      call h5sclose_f(dspaceid_f, hdferr)
      !----- Close out the dataset. -------------------------------------------------------!
      call h5dclose_f(dsetid_f, hdferr)
#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'shdf5_irec_f','hdf5_utils.F90')

#endif

      return
   end subroutine shdf5_irec_f
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will close a HDF5 file.                                           !
   !---------------------------------------------------------------------------------------!
   subroutine shdf5_close_f()

#if USE_HDF5
      use hdf5
      use hdf5_coms
      implicit none
      !----- Local variables. -------------------------------------------------------------!
      integer :: hdferr  ! Error flags
      !------------------------------------------------------------------------------------!
      
      !----- Close the HDF5 file. ---------------------------------------------------------!
      call h5fclose_f(fileid_f, hdferr)
      !----- Collect garbage. -------------------------------------------------------------!
      call h5garbage_collect_f(hdferr)
      !----- Close the HDF5 environment too. ----------------------------------------------!
      call h5close_f(hdferr)
#else
      call fatal_error('Cannot use hdf5 routines without compiling RAPP with HDF5.'        &
                         ,'shdf5_close_f','hdf5_utils.F90')
#endif
      return
   end  subroutine shdf5_close_f
   !=======================================================================================!
   !=======================================================================================!
end module hdf5_utils
!==========================================================================================!
!==========================================================================================!
