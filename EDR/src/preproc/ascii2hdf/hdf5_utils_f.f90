module hdf5_utils


contains
  
!===============================================================================

subroutine shdf5_info_f(dsetname,ndims,dims)

  use hdf5_coms

implicit none

character(len=*) :: dsetname ! Dataset name
integer :: dims(*)
integer(HSIZE_T),dimension(4) :: dimshf,maxdims

integer :: ndims ! Dataset rank (in file)
integer :: hdferr ! Error flag

! Open the dataset.

call h5dopen_f(fileid_f,trim(dsetname)//char(0), dsetid_f, hdferr)

call h5dget_space_f(dsetid_f, dspaceid_f, hdferr)

if (hdferr < 0) then
   print*, 'In shdf5_info:'
   print*, 'Variable ', trim(dsetname), ' is not in the currently opened hdf5 file'
   ndims   = 0
   dims(1) = 0
   return
endif

! Get dataset's dimensions

call h5sget_simple_extent_ndims_f(dspaceid_f, ndims, hdferr) 
call h5Sget_simple_extent_dims_f(dspaceid_f, dimshf, maxdims, hdferr  )

dims(1:ndims) = dimshf(1:ndims)

call h5dclose_f(dsetid_f, hdferr)
if (hdferr.ne.0) then
   print*,"COULD NOT CLOSE THE DATASET"
   stop
endif

call h5sclose_f(dspaceid_f, hdferr)
if (hdferr.ne.0) then
   print*,"COULD NOT CLOSE DATASPACE"
   stop
endif


!call fh5s_get_ndims(ndims)
!call fh5s_get_dims(dims)

!print*,'ndims: ',ndims
!print*,'dims: ',dims(1:ndims)

return
end subroutine shdf5_info_f

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                              ,ivars,rvars,cvars,dvars,lvars)
  use hdf5_coms
  implicit none

  character(len=*)           :: dsetname, action
  integer                    :: ndims, dims(*)
  integer,          optional :: ivara(*), ivars
  real,             optional :: rvara(*), rvars
  character(len=*), optional :: cvara(*), cvars
  real(kind=8),     optional :: dvara(*), dvars
  logical,          optional :: lvara(*), lvars
 
  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (trim(action) == 'READ') then
     
     call shdf5_irec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
          ,ivars,rvars,cvars,dvars,lvars)
     
  elseif (trim(action) == 'WRITE') then
     
     call shdf5_orec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
          ,ivars,rvars,cvars,dvars,lvars)
     
  else
     
     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif
  
end subroutine shdf5_io

! ========================================

subroutine shdf5_open_f(locfn,access,idelete)

  use hdf5_coms

  implicit none

  character(len=*) :: locfn     ! file name
  character(len=*) :: access    ! File access ('R','W','RW')
  integer, optional :: idelete  ! If W, delete/overwrite file if exists? 1=yes, 0=no
  ! Only needed when access='W'

  integer :: hdferr ! Error flag
  integer :: iaccess ! int access flag
  character(len=2) :: caccess ! File access ('R ','W ','RW')
  
  logical :: exists ! File existence
  
  caccess = access
  
  ! Check for existence of RAMS file.

  ! Open the HDF environment
  
  call h5open_f(hdferr)
  
  inquire(file=trim(locfn),exist=exists)
  
  ! Create a new file or open an existing RAMS file.
  if (access(1:1) == 'R') then
     if (.not.exists) then
        print*,'shdf5_open:'
        print*,'   Attempt to open a file for reading that does not exist.'
        print*,'   Filename: ',trim(locfn)
        stop 'shdf5_open: no file'
     else
        if (caccess == 'R ') then
           call h5fopen_f(trim(locfn)//char(0), H5F_ACC_RDONLY_F, fileid_f, hdferr)
        else if (caccess == 'RW') then
           call h5fopen_f(trim(locfn)//char(0), H5F_ACC_RDWR_F, fileid_f, hdferr)
        else
           print*,"INVALID FILE READ TYPE ACCESS"
           print*,"SHOULD BE READ OR READ-WRITE"
           stop
        endif
        
        if (hdferr < 0) then
           print*,'shdf5_open:'
           print*,'   Error opening hdf5 file - error -',hdferr
           print*,'   Filename: ',trim(locfn)
           stop 'shdf5_open: open error'      
        endif
     endif
  elseif (access(1:1) == 'W') then
     
     if (.not.exists) then

        call h5fcreate_f(trim(locfn)//char(0), H5F_ACC_EXCL_F, fileid_f, hdferr)
        
     else
        if(.not.present(idelete) ) then
           print*,'shdf5_open: idelete not specified when access=W'
           stop 'shdf5_open: no idelete'
        endif
        
        if(idelete == 0) then
           print*,'In shdf5_open:'
           print*,'   Attempt to open an existing file for writing, '
           print*,'      but overwrite is disabled. idelete=',idelete
           print*,'   Filename: ',trim(locfn)
           stop 'shdf5_open'
        else
           call system('rm -f '//trim(locfn)//char(0))
           
           call h5fcreate_f(trim(locfn)//char(0), H5F_ACC_TRUNC_F, fileid_f, hdferr)
           
        endif
     endif
     if(hdferr < 0) then
        print*,'HDF5 file create failed:',hdferr
        print*,'file name:',trim(locfn),' ',trim(access), idelete
        stop 'shdf5_open: bad create'
     endif
  endif
  
  return
end  subroutine shdf5_open_f

! =======================================================================

subroutine shdf5_orec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars)
  
  use hdf5_coms
  implicit none
  
  character(len=*) :: dsetname ! Variable label
  integer :: ndims             ! Number of dimensions or rank
  integer, dimension(*) :: dims ! Dataset dimensions.
  
  ! Array and scalar arguments for different types. Only specify one in each call.
  integer,          optional :: ivara(*),ivars
  real,             optional :: rvara(*),rvars
  character(len=*), optional :: cvara(*),cvars
  real(kind=8),     optional :: dvara(*),dvars
  logical,          optional :: lvara(*),lvars
  
  integer:: h5_type   ! Local type designator
  
  integer, dimension(4) :: dimsh ! Dataset dimensions.
  
  integer(HSIZE_T),dimension(4) :: dimshf,chunksize
  
  character(len=2) :: ctype    ! Variable type: int, real, char
  integer :: hdferr ! Error flag
  
  ! Find which data type is input
  if(present(ivars)) then ; ctype='is'
  elseif(present(rvars)) then ; ctype='rs'
  elseif(present(cvars)) then ; ctype='cs'
  elseif(present(dvars)) then ; ctype='ds'
  elseif(present(lvars)) then ; ctype='ls'
  elseif(present(ivara)) then ; ctype='ia'
  elseif(present(rvara)) then ; ctype='ra'
  elseif(present(cvara)) then ; ctype='ca'
  elseif(present(dvara)) then ; ctype='da'
  elseif(present(lvara)) then ; ctype='la'
  else
     print*,'Incorrect or missing data field argument in shdf5_orec'
     stop 'shdf5_orec: bad data field'
  endif
  
  ! Check dimensions and set compression chunk size
  
  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*,'Dimension error in shdf5_orec:',ndims,dims(1:ndims)
     stop 'shdf5_orec: bad dims'
  endif
  
  dimsh(1:ndims) = dims(1:ndims)
  dimshf = 0
  chunksize = 0
  dimshf(1:ndims) = dimsh(1:ndims)
  chunksize(1:ndims) = dimsh(1:ndims)
  
  ! Prepare memory and options for the write
  
  call h5screate_simple_f(ndims,dimshf,dspaceid_f,hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT CREATE DATASPACE"
     print*,ndims,dimshf,dspaceid_f,hdferr
     stop
  endif
  
  call h5pcreate_f(H5P_DATASET_CREATE_F, prpid_f, hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT CREATE PROPERTY"
     print*,prpid_f,hdferr
     stop
  endif
  
  call h5pset_chunk_f(prpid_f, ndims, chunksize, hdferr) 
  if (hdferr.ne.0) then
     print*,"COULD NOT SET CHUNK"
     print*,ndims,chunksize,hdferr
     stop
  endif
  
  call h5pset_shuffle_f(prpid_f, hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT SET SHUFFLE"
     print*,hdferr
     stop
  endif
  
  call h5pset_deflate_f(prpid_f, 9, hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT SET DEFLATE"
     print*,hdferr
     stop
  endif
  
  if (ctype(1:1) == 'i') h5_type=1
  if (ctype(1:1) == 'r') h5_type=2
  if (ctype(1:1) == 'c') h5_type=3
  if (ctype(1:1) == 'd') h5_type=4  ! If native precision is 8 bytes, do h5_type=2
  if (ctype(1:1) == 'l') h5_type=5
  
  ! Write the dataset.
  if (ctype == 'is') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_INTEGER, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_INTEGER,ivars,dimshf,hdferr)
     
  elseif (ctype == 'rs') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_REAL, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_REAL,rvars,dimshf,hdferr)
     
  elseif (ctype == 'cs') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_CHARACTER, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_CHARACTER,cvars,dimshf,hdferr)
     
  elseif (ctype == 'ds') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_DOUBLE, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_DOUBLE,dvars,dimshf,hdferr)
     
  elseif (ctype == 'ls') then
     print*,"THERE IS NO HDF5 FORTRAN API DATATYPE FOR BOOLEAN"
     print*,"YOU MUST CHANGE BACK TO C IO FOR THIS"
     print*,"STOPPING"
     stop
     
  elseif (ctype == 'ia') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_INTEGER, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_INTEGER,ivara,dimshf,hdferr)
     
  elseif (ctype == 'ra') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_REAL, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_REAL,rvara,dimshf,hdferr)
     
  elseif (ctype == 'ca') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_CHARACTER, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_CHARACTER,cvara,dimshf,hdferr)
     
  elseif (ctype == 'da') then
     
     call h5dcreate_f(fileid_f,trim(dsetname)//char(0), H5T_NATIVE_DOUBLE, dspaceid_f,dsetid_f, hdferr)
     call h5dwrite_f(dsetid_f,H5T_NATIVE_DOUBLE,dvara,dimshf,hdferr)
     
  elseif (ctype == 'la') then
     print*,"THERE IS NO HDF5 FORTRAN API DATATYPE FOR BOOLEAN"
     print*,"YOU MUST CHANGE BACK TO C IO FOR THIS"
     print*,"STOPPING"
     stop
     
  endif
  if (hdferr /= 0) then
     print*,'In shdf5_orec: hdf5 write error =',hdferr
     stop 'shdf5_orec: hdf5 write error'
  endif
  
  ! Close the dataset, the dataspace for the dataset, and the dataspace properties.
  
  !   Close out the dataspace
  call h5sclose_f(dspaceid_f, hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT CLOSE DATASPACE"
     stop
  endif
  
  !   Close the dataset
  call h5dclose_f(dsetid_f, hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT CLOSE THE DATASET"
     stop
  endif
  
  !   Close the property
  call h5pclose_f(prpid_f, hdferr)
  if (hdferr.ne.0) then
     print*,"COULD NOT CLOSE THE PROPERTY"
     stop
  endif
    
  return
end subroutine shdf5_orec_f

! ==============================================================================

subroutine shdf5_irec_f(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars)
  use hdf5_coms
  implicit none

  character(len=*) :: dsetname ! Dataset name
  integer :: ndims             ! Number of dimensions or rank
  integer, dimension(*) :: dims ! Dataset dimensions.
  
  ! Array and scalar arguments for different types. Only specify one in each call.
  integer,          optional :: ivara(*),ivars
  real,             optional :: rvara(*),rvars
  character(len=*), optional :: cvara(*),cvars
  real(kind=8),     optional :: dvara(*),dvars
  logical,          optional :: lvara(*),lvars
  
  integer:: h5_type   ! Local type designator
  
  integer, dimension(4) :: dimsh ! Dataset dimensions.
  
  integer(HSIZE_T),dimension(4) :: dimshf
  
  integer :: hdferr ! Error flag
  
  character(len=2) :: ctype
  
  
  ! Find which data type will be read
  if(present(ivars)) then ; ctype='is'
  elseif(present(rvars)) then ; ctype='rs'
  elseif(present(cvars)) then ; ctype='cs'
  elseif(present(dvars)) then ; ctype='ds'
  elseif(present(lvars)) then ; ctype='ls'
  elseif(present(ivara)) then ; ctype='ia'
  elseif(present(rvara)) then ; ctype='ra'
  elseif(present(cvara)) then ; ctype='ca'
  elseif(present(dvara)) then ; ctype='da'
  elseif(present(lvara)) then ; ctype='la'
  else
     print*,'Incorrect or missing data field argument in shdf5_irec'
     stop 'shdf5_irec: bad data field'
  endif
  
  ! Check dimensions and set compression chunk size
  
  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*,'Dimension error in shdf5_irec:',ndims,dims(1:ndims)
     stop 'shdf5_irec: bad dims'
  endif
  
  dimsh(1:ndims) = dims(1:ndims)
  dimshf = 0
  dimshf(1:ndims) = dimsh(1:ndims)
  
  call h5dopen_f(fileid_f,trim(dsetname)//char(0), dsetid_f, hdferr)
  call h5dget_space_f(dsetid_f, dspaceid_f, hdferr)
  
  if (ctype == 'is') then
     call h5dread_f(dsetid_f, H5T_NATIVE_INTEGER, ivars, dimshf, hdferr)
  elseif (ctype == 'rs') then
     call h5dread_f(dsetid_f, H5T_NATIVE_REAL,rvars, dimshf, hdferr )
  elseif (ctype == 'cs') then
     call h5dread_f(dsetid_f,H5T_NATIVE_CHARACTER,cvars, dimshf, hdferr )
  elseif (ctype == 'ds') then
     call h5dread_f(dsetid_f,H5T_NATIVE_DOUBLE,dvars, dimshf, hdferr )
  elseif (ctype == 'ls') then
     !      call h5dread_f(dsetid_f,H5T_NATIVE_HBOOL,lvars, dimsh, hdferr )
     print*,"THERE IS NO HDF5 FORTRAN API DATATYPE FOR BOOLEAN"
     print*,"YOU MUST CHANGE BACK TO C IO FOR THIS"
     print*,"STOPPING"
     stop
  elseif (ctype == 'ia') then
     call h5dread_f(dsetid_f,H5T_NATIVE_INTEGER,ivara, dimshf, hdferr )
  elseif (ctype == 'ra') then
     call h5dread_f(dsetid_f,H5T_NATIVE_REAL,rvara, dimshf, hdferr )
  elseif (ctype == 'ca') then
     call h5dread_f(dsetid_f,H5T_NATIVE_CHARACTER,cvara, dimshf, hdferr )
  elseif (ctype == 'da') then
     call h5dread_f(dsetid_f,H5T_NATIVE_DOUBLE,dvara, dimshf, hdferr )
  elseif (ctype == 'la') then
     print*,"THERE IS NO HDF5 FORTRAN API DATATYPE FOR BOOLEAN"
     print*,"YOU MUST CHANGE BACK TO C IO FOR THIS"
     print*,"STOPPING"
     stop
  endif
  
  
  call h5sclose_f(dspaceid_f, hdferr)
  call h5dclose_f(dsetid_f, hdferr)
  
  
  return
end subroutine shdf5_irec_f

! ===========================

subroutine shdf5_close_f()
  
  use hdf5_coms
  implicit none
  
  integer :: hdferr  ! Error flags
  
  ! Close the hdf file.
  
  call h5fclose_f(fileid_f, hdferr)
  
  ! Close the hdf environment too
  
  call h5close_f(hdferr)
  
  
  return
end  subroutine shdf5_close_f

!===============================================================================


end module hdf5_utils
