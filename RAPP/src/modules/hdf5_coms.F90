module hdf5_coms

#if USE_HDF5
  use hdf5
  use mod_maxdims, only : maxrank
  integer(HID_T) :: dsetid_f,dspaceid_f,fileid_f,prpid_f
  integer(HID_T) :: plist_id,file_id,dset_id,filespace,memspace,dspace_id

  integer(HSIZE_T),dimension(maxrank) :: globdims,chnkdims,chnkoffs
  integer(HSIZE_T),dimension(maxrank) :: cnt,stride

  integer(HSIZE_T),dimension(maxrank) :: memsize,memdims,memoffs !EQUIV IN MEM

  integer(HID_T)  :: datatype_id ! Datatype identifier 
  integer(SIZE_T) :: setsize  ! Datatype size
#endif

end module hdf5_coms
