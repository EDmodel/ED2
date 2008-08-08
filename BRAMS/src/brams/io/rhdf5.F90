! =========================================================
!
! Create H5 output files from the var table
! =========================================================
subroutine anlhdf(vtype)

  use an_header
  use var_tables
  use grid_dims,only : maxgrds
  use mem_aerad, only: nwave
  use mem_cuparm, only: nclouds

  use mem_grid
  use io_params

#if USE_HDF5
  use hdf5
#endif

  implicit none

  ! This routine writes the chosen variables on the analysis file.

  character*(*) vtype

#if USE_HDF5

  character(len=128) :: anamel
  character(len=2)  :: cgrid
  character(len=25) :: subaname
  character(len=16) :: varn
  character(len=1)  :: vnam
  character(len=10) :: c0
  logical exans
  integer, save :: ncall_head=0,nvtota=0,nvtotl=0  &
       ,nvtot
  !  HDF specific data types

  integer, save :: ihdfinit = 0
  integer       :: hdferr
  integer       :: dsetrank
  integer(HID_T):: dspace_id,file_id,dset_id
  integer(HSIZE_T),dimension(4) :: datadims,maxdatadims

  
  integer :: ngr,nv,nvcnt,lenl
  real(kind=8) :: timeold
  real(kind=8),parameter :: zero = 0.0

  integer :: outyear,outmonth,outdate,outhour

  if (ioutput /= 3) return

  if (ncall_head == 0) then

     !  Find total number of fields to be written
     do ngr=1,ngrids
        do nv = 1,num_var(ngr)
           if ( vtab_r(nv,ngr)%ianal == 1) nvtota=nvtota+1
           if ( vtab_r(nv,ngr)%ilite == 1) nvtotl=nvtotl+1
        enddo
     enddo

     nvtot=max(nvtota,nvtotl)
     
     ncall_head=1
  endif
  
  timeold=time

  ! Construct the HDF file name

  if(vtype == 'INST') vnam='A'
  if(vtype == 'LITE') vnam='L'
  if(vtype == 'MEAN') vnam='M'
  if(vtype == 'BOTH') vnam='B'
  
  nvcnt=0

  do ngr=1,ngrids

     write(cgrid,'(a1,i1)') 'g',ngr

        
     call makefnam(anamel,afilout,time,iyeara,imontha,idatea,  &
                   itimea*100,vnam,cgrid,'h5')

     lenl = len_trim(anamel)

     inquire(file=anamel,exist=exans)
     if(exans.and.iclobber.eq.0) then
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*,'!!!   trying to open file name :'
        print*,'!!!       ',anamel
        print*,'!!!   but it already exists. run is ended.'
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        stop 'anlwrt'
     endif

     
     !   LETS FIRST INITIALIZE THE HDF ENVIRONMENT
     call h5open_f(hdferr)
        
     !   Open a new HDF file using IO mode H5F_ACC_TRUNC_F
     !   In this case, if a file with the same name already exists
     !   It will overwrite the data in that file, destroying all data
     !   ------------------------------------------------------------

     call h5fcreate_f(anamel, H5F_ACC_TRUNC_F, file_id, hdferr)


     !   Now we need to create HDF datasets and then put them
     !   in the file, cycle all of our variables
     !   ----------------------------------------------------

     do nv = 1,num_var(ngr)

              
        if ((vtype == 'INST' .and. vtab_r(nv,ngr)%ianal == 1) .or. &
            (vtype == 'LITE' .and. vtab_r(nv,ngr)%ilite == 1)) then

           varn= vtab_r(nv,ngr)%name
           
           datadims = 0
           maxdatadims = -1
           
           !   There are multiple possible dimensions to this
           select case (vtab_r(nv,ngr)%idim_type)
           case(2)   !   2-D
              datadims(1) = nnxp(ngr)
              datadims(2) = nnyp(ngr)
              dsetrank    = 2

           case(3)        !  3-D
              datadims(1) = nnzp(ngr)
              datadims(2) = nnxp(ngr)
              datadims(3) = nnyp(ngr)
              dsetrank    = 3
              
           case(4)    !  2-D w/ soil,patches
              datadims(1) = nzg
              datadims(2) = nnxp(ngr)
              datadims(3) = nnyp(ngr)
              datadims(4) = npatch
              dsetrank    = 4
              
           case(5)    !  2-D w/ snow,patches
              datadims(1) = nzs
              datadims(2) = nnxp(ngr)
              datadims(3) = nnyp(ngr)
              datadims(4) = npatch
              dsetrank    = 4
              
           case (6)   !   2-D w/patches
              datadims(1) = nnxp(ngr)
              datadims(2) = nnyp(ngr)
              datadims(3) = npatch
              dsetrank    = 3
              
           case (7)   !   2-D w/patches
              datadims(1) = nnxp(ngr)
              datadims(2) = nnyp(ngr)
              datadims(3) = nwave
              dsetrank    = 3
              
           case (8)   !   3-D with cloud spectrum
              datadims(1) = nnzp(ngr)
              datadims(2) = nnxp(ngr)
              datadims(3) = nnyp(ngr)
              datadims(4) = nclouds
              dsetrank    = 4
              
           case (9)   !   2-D with cloud spectrum
              datadims(1) = nnxp(ngr)
              datadims(2) = nnyp(ngr)
              datadims(3) = nclouds
              dsetrank    = 3

           end select
           
           !   Create a dataspace

           call h5screate_simple_f(dsetrank,datadims,dspace_id,hdferr)
           if (hdferr.ne.0) then
              print*,"COULD NOT CREATE DATASPACE"
              print*,trim(varn),dsetrank,datadims,dspace_id,hdferr
              stop
           endif


           !   Create a dataset, in the file, with the space just specified
           
           call h5dcreate_f(file_id,varn, H5T_NATIVE_REAL, dspace_id, &
                dset_id, hdferr)
           if (hdferr.ne.0) then
              print*,"COULD NOT CREATE DATASET"
              print*,trim(varn),dsetrank,datadims,dset_id,hdferr
              stop
           endif
        

           !   Fill that dataset, in the file, with the real data
           
           call h5dwrite_f(dset_id,H5T_NATIVE_REAL,vtab_r(nv,ngr)%var_p,datadims,hdferr)
           if (hdferr.ne.0) then
              print*,"COULD NOT WRITE TO THE DATASET"
              print*,trim(varn),dsetrank,datadims,dset_id,hdferr
              stop
           endif
           
           !   Close out the dataspace
           call h5sclose_f(dspace_id, hdferr)
           if (hdferr.ne.0) then
              print*,"COULD NOT CLOSE DATASPACE"
              print*,trim(varn),dspace_id
              stop
           endif
           
           !   Close the dataset, in the file
           call h5dclose_f(dset_id, hdferr)
           if (hdferr.ne.0) then
              print*,"COULD NOT CLOSE THE DATASET"
              print*,trim(varn),dset_id
              stop
           endif

        endif
        
     enddo
     
     !   Close our beloved HDF file
     !   --------------------------
     
     call h5fclose_f(file_id, hdferr)
     if (hdferr.ne.0) then
        print*,"COULD NOT CLOSE THE HDF FILE"
        print*,file_id
        stop
     endif
     
     

  enddo

  call h5close_f(hdferr)
  if (hdferr /= 0) then
     print*,"COULD NOT CLOSE THE HDF ENVIRONMENT"
     stop
  endif


  if(vtype == 'LITE')then
     subaname='  Analysis lite HDF write'
  elseif(vtype == 'MEAN')then
     subaname='  Averaged analysis HDF write    '
  elseif(vtype == 'BOTH')then
     subaname='  Averaged analysis lite HDF write   '
  else
     subaname='  Analysis HDF write         '
  endif

  write(c0,"(f10.1)") time
  write(*,"(/,a)") " === "//trim(adjustl(subaname))//" at Sim time "//trim(adjustl(c0))//" ==="
  write(*,"(a,/)") " === wrote file "//trim(adjustl(anamel))//" ==="

  ! Reset the time back to the original
  if(vtype == 'MEAN'.or.vtype == 'BOTH')  time=timeold

#endif


  return
end subroutine anlhdf
