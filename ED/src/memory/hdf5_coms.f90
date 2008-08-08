module hdf5_coms

  use hdf5
  
  integer(HID_T) :: dsetid_f,dspaceid_f,fileid_f,prpid_f
  integer(HID_T) :: plist_id,file_id,dset_id,filespace,memspace,dspace_id

  integer(HSIZE_T),dimension(4) :: globdims,chnkdims,chnkoffs
  integer(HSIZE_T),dimension(4) :: cnt,stride

  integer(HSIZE_T),dimension(4) :: memsize,memdims,memoffs !EQUIV IN MEM



end module hdf5_coms
