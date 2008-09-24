! =========================================================
!
! Create H5 output files from the var table

! =========================================================
subroutine h5_output(vtype)

  use an_header

  use var_tables_array,only: &
       vt_info,              &
       var_table,            &
       var_table_vector,     &
       num_var

  use mem_sites,only : maxcohort

  use misc_coms, only : ffilout, &
                        sfilout, &
	                itimea,  &
	                idatea,  &
	                imontha, &
	                iyeara,  &
                        iclobber, &
                        nrec_fast, &
                        nrec_state, &
                        irec_fast, &
                        irec_state, &
                        out_time_fast, &
                        out_time_state, &
                        outstate, &
                        outfast, frqfast

  use ed_misc_coms,only: attach_metadata

  use grid_coms,only  :  ngrids, &
                        time
            
  use hdf5_coms,only : chnkdims,chnkoffs,cnt,stride,globdims

  use hdf5
  
  use ed_node_coms,only:mynum,nnodetot,recvnum,sendnum
  use max_dims, only : n_pft,n_dist_types,n_dbh,maxgrds
  use ed_state_vars,only: edgrid_g,edtype,polygontype,sitetype,patchtype

  implicit none


  include 'mpif.h'

  ! This routine writes the chosen variables on the analysis file.

  character*(*) vtype

  type(edtype),pointer      :: egrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: ipy,isi,ipa,ico

  character(len=128) :: anamel
  character(len=2)  :: cgrid
  character(len=25) :: subaname
  character(len=64) :: varn
  character(len=1)  :: vnam
  character(len=10) :: c0
  character(len=64),dimension(3) :: metadata
  logical exans
  integer, save :: ncall_head=0,nvtota=0,nvtotl=0  &
       ,nvtot
  !  HDF specific data types

  integer, save :: ihdfinit = 0
  integer       :: hdferr
  integer       :: dsetrank
  integer(HID_T):: file_id,dset_id
  integer(HID_T):: memspace,filespace
  integer(HID_T):: plist_id
  integer(HID_T) :: attr_id

  ! Attribute types
  integer(HID_T) :: aspace_id
  integer(HID_T) :: atype_id
  integer(HSIZE_T),dimension(1) :: adims
  integer :: arank
  integer(SIZE_T) :: attrlen
  
  integer :: ngr,nv,nvcnt,lenl,iptr
  real(kind=8) :: timeold
  real(kind=8),parameter :: zero = 0.0

  integer :: outyear,outmonth,outdate,outhour,outminu,outseco
  real, external :: dbh2bl
  type(var_table),pointer        :: vtinfo
  type(var_table_vector),pointer :: vtvec
  integer :: irec, nrec

  integer :: mpierror
  integer :: comm,info
  integer :: mpi_size,mpi_rank
  integer :: ping,ierr
  integer,       dimension(MPI_STATUS_SIZE) :: status
  real(kind=8)    :: dsec
  logical :: new_file = .true.

  logical,parameter :: collective_mpi = .false.


  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  
  if (nnodetot /= 1 ) then
     call MPI_COMM_SIZE(comm,mpi_size,mpierror)
     call MPI_COMM_RANK(comm,mpi_rank,mpierror)
!     call MPI_Barrier(MPI_COMM_WORLD,ierr)
  endif

  ping = 456

  timeold=time

  ! Construct the HDF file name

  select case (trim(vtype))
  case ('INST')
     vnam='I'
  case ('LITE')
     vnam='Z'
  case ('MEAN') 
     vnam='N'
  case ('BOTH') 
     vnam='O'
  case ('DAIL') 
     vnam='D'
  case ('MONT') 
     vnam='E'
  case ('YEAR') 
     vnam='Y'
  case('HIST') 
     vnam='S'  ! S for reStart, don't want confusion with RAMS' H or R files
  end select  
  nvcnt=0

  ! DOING SEQUENTIAL IO RIGHT NOW - THE NEXT 4 LINES ARE A BLOCKING PROCEDURE

  if (.not. collective_mpi) then
     if (mynum /= 1) call MPI_RECV(ping,1,MPI_INTEGER,recvnum,734,MPI_COMM_WORLD,status,ierr)
  end if


  nrec = 1
  irec = 1  
  do ngr=1,ngrids
     
     write(cgrid,'(a1,i1)') 'g',ngr

     select case (trim(vtype))
     case('DAIL')
         
        !Return the current year,month and day of the last 24hrs
        call date_add_to (iyeara,imontha,idatea,itimea*100,  &
             time-21600,'s',outyear,outmonth,outdate,outhour)
        
        call makefnam(anamel,ffilout,zero,outyear,outmonth,outdate, &
             0,vnam,cgrid,'h5 ')

     case('MONT')
         
        call date_add_to (iyeara,imontha,idatea,itimea*100,  &
             time-21600,'s',outyear,outmonth,outdate,outhour)
        
        
        call makefnam(anamel,ffilout,zero,outyear,outmonth,0, &
             0,vnam,cgrid,'h5 ')

     case('YEAR')
         
        call date_add_to (iyeara,imontha,idatea,itimea*100,  &
             time-21600,'s',outyear,outmonth,outdate,outhour)
        
        
        call makefnam(anamel,ffilout,zero,outyear,0,0, &
             0,vnam,cgrid,'h5 ')

     case('HIST')
        call makefnam(anamel,sfilout,time,iyeara,imontha,idatea,  &
             itimea*100,vnam,cgrid,'h5 ')
        
     case default
        if(nrec_fast .eq. 1) then  !! single file per output
           call makefnam(anamel,ffilout,time,iyeara,imontha,idatea,  &
                itimea*100,vnam,cgrid,'h5 ')
        else   !! group outputs
           new_file = .false.
           !! determine whether to advance out_time_fast
           call date_add_to (iyeara,imontha,idatea,itimea*100,  &
                time,'s',outyear,outmonth,outdate,outhour)
           call date_2_seconds(out_time_fast%year,out_time_fast%month, &
                out_time_fast%date,int(out_time_fast%time), &
                iyeara,imontha,idatea,itimea*100,dsec)
           if(time >= (dsec+outfast) .or. outmonth .ne. out_time_fast%month) then
              print*,"SYNC"
              out_time_fast%year  = outyear
              out_time_fast%month = outmonth
              out_time_fast%date  = outdate
              out_time_fast%time  = outhour*1.!!(3600.*int(outhour/10000)+60.*int(mod(outhour,10000)/100)+mod(outhour,100)*1.)   !! DOUBLE CHECK
              dsec = time
              new_file = .true.
           endif
           irec_fast = ((time-dsec)/frqfast) + 1
           nrec = nrec_fast
           irec = irec_fast
print*,irec,nrec,outmonth,out_time_fast,dsec,time
           !! construct filename
           call makefnam(anamel,ffilout,0.0,out_time_fast%year, &
                out_time_fast%month,out_time_fast%date,out_time_fast%time*100, &
                vnam,cgrid,'h5 ')
        endif
        
     end select

     lenl = len_trim(anamel)

     inquire(file=anamel,exist=exans)
     if(exans .and. iclobber == 0) then
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print*,'!!!   trying to open file name :'
        print*,'!!!       ',anamel
        print*,'!!!   but it already exists. run is ended.'
        print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        call fatal_error('File '//trim(anamel)//' already exists' &
                       ,'h5_output','h5_output.F90')

     endif

     !   LETS FIRST INITIALIZE THE HDF ENVIRONMENT

     call h5open_f(hdferr)
     if (hdferr /= 0) then
        print*,'HDF5 Open error #:',hdferr
        call fatal_error('Could not initialize the hdf environment' &
                        ,'h5_output','h5_output.F90')
     endif
      
     if(collective_mpi) then

        call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdferr)
        if (hdferr /= 0) &
             call fatal_error('Could not create the p-list' &
             ,'h5_output','h5_output.F90')
                 
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
! COMMENTED BECAUSE IT IS NOT IN SOME LIBRARIES...
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!         call h5pset_fapl_mpio_f(plist_id,comm,info,hdferr)
!         if (hdferr /= 0) &
!            call fatal_error('Failed using h5pset_fapl_mpi_f' &
!                            ,'h5_output','h5_output.F90')
         
         !   Open a new HDF file using IO mode H5F_ACC_TRUNC_F
         !   In this case, if a file with the same name already exists
         !   It will overwrite the data in that file, destroying all data
         !   ------------------------------------------------------------
         
         call h5fcreate_f(trim(anamel)//char(0), H5F_ACC_TRUNC_F, file_id, &
              hdferr,access_prp = plist_id)
         if (hdferr /= 0) then
            print*,"COULD NOT OPEN THE NEW HDF FILE"
            print*,trim(anamel),file_id,hdferr
            call fatal_error('Failed opening the HDF file' &
                 ,'h5_output','h5_output.F90')
         endif
         
         call h5pclose_f(plist_id,hdferr)
         if (hdferr /= 0) &
             call fatal_error('Could not close the p-list' &
                             ,'h5_output','h5_output.F90')

     else
        
        if (mynum == 1 .and. new_file) then

           call h5fcreate_f(trim(anamel)//char(0), H5F_ACC_TRUNC_F, file_id, hdferr)
           if (hdferr /= 0) then
              print*,"COULD NOT OPEN THE new HDF FILE"
              print*,trim(anamel),file_id,hdferr
              call fatal_error('Failed opening the HDF file' &
                              ,'h5_output','h5_output.F90')
           endif

        else
           
           call h5fopen_f(trim(anamel)//char(0), H5F_ACC_RDWR_F, file_id, hdferr)
           if (hdferr /= 0) then
              print*,"COULD NOT OPEN THE existing HDF FILE"
              print*,trim(anamel),file_id,hdferr
              call fatal_error('Failed opening the HDF file' &
                              ,'h5_output','h5_output.F90')
           endif
        endif

           
     endif
 
       
     if (collective_mpi) then
        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdferr)
        if (hdferr /= 0) call fatal_error('Errror at pcreate' &
                                         ,'h5_output','h5_output.f90')
        
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
! COMMENTED BECAUSE IT IS NOT IN SOME LIBRARIES...
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------!
!        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdferr)
!        if (hdferr /= 0) call fatal_error('Errror at h5pset_dxpl_mpio_f' &
!                                         ,'h5_output','h5_output.f90')
        
     endif

     
     !   Now we need to create HDF datasets and then put them
     !   in the file, cycle all of our variables
     !   ----------------------------------------------------

     varloop: do nv = 1,num_var(ngr)

!        vtinfo => vt_info(nv,ngr)

        if ((vtype == 'INST' .and. vt_info(nv,ngr)%ianal == 1) .or. &
            (vtype == 'LITE' .and. vt_info(nv,ngr)%ilite == 1) .or. &
            (vtype == 'DAIL' .and. vt_info(nv,ngr)%idail == 1) .or. &
            (vtype == 'MONT' .and. vt_info(nv,ngr)%imont == 1) .or. &
            (vtype == 'YEAR' .and. vt_info(nv,ngr)%iyear == 1) .or. &
            (vtype == 'HIST' .and. vt_info(nv,ngr)%ihist == 1)) then

           varn= vt_info(nv,ngr)%name

           ! Initialize global dimensions of the hyperslab
           
           call geth5dims(vt_info(nv,ngr)%idim_type,0,0,vt_info(nv,ngr)%var_len_global,dsetrank,varn,nrec,irec)

           !if this is a 1d vector
!           globdims(2) = nrec
!           chnkdims(2) = 1
!           chnkoffs(2) = irec

           call h5screate_simple_f(dsetrank, globdims, filespace, hdferr)
           if (hdferr /= 0) then
              call fatal_error('Could not create the first filespace' &
                                         ,'h5_output','h5_output.f90')
           end if

           if (collective_mpi) then
                 if (vt_info(nv,ngr)%dtype == 'r') then   ! real data type
                    call h5dcreate_f(file_id,varn,H5T_NATIVE_REAL, filespace, &
                         dset_id,hdferr)
                 else if (vt_info(nv,ngr)%dtype == 'i') then   ! integer data type
                    call h5dcreate_f(file_id,varn,H5T_NATIVE_INTEGER, filespace, &
                         dset_id,hdferr)
                 else if (vt_info(nv,ngr)%dtype == 'c') then   ! character data type
                    call h5dcreate_f(file_id,varn,H5T_NATIVE_CHARACTER, filespace, &
                         dset_id,hdferr)
                 else
                    print*,"YOU ARE ATTEMPTING TO WRITE AN UNDEFINED DATATYPE"
                    print*,varn,vt_info(nv,ngr)%dtype
                    stop
                    
                 endif
              else

              if (mynum == 1 .and. new_file) then
                 
                 if (vt_info(nv,ngr)%dtype == 'r') then   ! real data type
                    call h5dcreate_f(file_id,varn,H5T_NATIVE_REAL, filespace, &
                         dset_id,hdferr)
                 else if (vt_info(nv,ngr)%dtype == 'i') then   ! integer data type
                    call h5dcreate_f(file_id,varn,H5T_NATIVE_INTEGER, filespace, &
                         dset_id,hdferr)
                 else if (vt_info(nv,ngr)%dtype == 'c') then   ! character data type
                    call h5dcreate_f(file_id,varn,H5T_NATIVE_CHARACTER, filespace, &
                         dset_id,hdferr)
                 else
                    print*,"YOU ARE ATTEMPTING TO WRITE AN UNDEFINED DATATYPE"
                    print*,varn,vt_info(nv,ngr)%dtype
                    stop
                    
                 endif

                 ! REMEMBER THESE COMMANDS
                 !              h5pset_meta_block_size
                 !              h5pget_meta_block_size
                 !              h5pset_cache

                 ! If the user has decided to attach metadata
                 ! to the datasets, assign that metadata as an 
                 ! attribute here.  That attribute is a rank 1 vector
                 ! of strings.  It is only necessary to do this
                 ! on the master.
                 ! Note that descriptors max out at 64 characters
                 
                 if (attach_metadata == 1) then
                    
                    arank = 1
                    adims = 3
                    attrlen = 64
                    metadata(1) = trim('Long Name: '//trim(vt_info(nv,ngr)%lname))
                    metadata(2) = trim('Units: '//trim(vt_info(nv,ngr)%units))
                    metadata(3) = trim('Dimensions: '//trim(vt_info(nv,ngr)%dimlab))
                    
                    call h5screate_simple_f(arank,adims,aspace_id, hdferr)
                    if (hdferr /= 0) then
                       call fatal_error('Error calling h5screate_simple_f' &
                                       ,'h5_output','h5_output.F90')
                    end if
                    
                    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
                    if (hdferr /= 0) then
                       call fatal_error('Error calling h5tcopy_f' &
                                       ,'h5_output','h5_output.F90')
                    end if
                    
                    call h5tset_size_f(atype_id,attrlen,hdferr)
                    if (hdferr /= 0) then
                       call fatal_error('Error calling h5tset_size_f' &
                                       ,'h5_output','h5_output.F90')
                    end if
                 
                    call h5acreate_f( dset_id,'Metadata',atype_id,aspace_id,attr_id,hdferr)
                    if (hdferr /= 0) then
                       call fatal_error('Error calling h5acreate_f' &
                                       ,'h5_output','h5_output.F90')
                    end if
                    
                    call h5awrite_f( attr_id,atype_id,metadata,adims,hdferr )
                    if (hdferr /= 0) then
                       call fatal_error('Error calling h5awrite_f' &
                                       ,'h5_output','h5_output.F90')
                    end if
                    
                    call h5aclose_f(attr_id,hdferr)
                    call h5sclose_f(aspace_id,hdferr)

                 end if
                    
                 
              else
                 call h5dopen_f(file_id,varn,dset_id,hdferr)
              end if
           end if

           if (hdferr /= 0) then
              write (unit=*,fmt=*) 'File name:        ',trim(anamel)
              write (unit=*,fmt=*) 'Variable name:    ',trim(varn)
              write (unit=*,fmt=*) 'File ID:          ',file_id
              write (unit=*,fmt=*) 'Dataset ID:       ',dset_id
              write (unit=*,fmt=*) 'Dataset rank:     ',dsetrank
              write (unit=*,fmt=*) 'Global dimension: ',globdims
              call fatal_error('Could not create the dataset','h5_output','h5_output.F90')
           end if

           call h5sclose_f(filespace,hdferr)
           if (hdferr /= 0) then
              call fatal_error('Could not close the first filespace' &
                              ,'h5_output','h5_output.F90')
           end if


           pointerloop: do iptr = 1,vt_info(nv,ngr)%nptrs
              
              vtvec => vt_info(nv,ngr)%vt_vector(iptr)
              
              ! Set the size of the chunk and it's offset in the
              ! global dataset
              

              if (vtvec%varlen > 0 ) then
                 

                 !  Evaluate the variable output type
                 !  Resolve the dimensioning and the meta-data tags
                 !  accordingly.  See ed_state_vars.f90 for a 
                 !  description of the various datatype.
                 !  -----------------------------------------------
                 
                 ! Initialize hyperslab indexes

                 call geth5dims(vt_info(nv,ngr)%idim_type,vtvec%varlen, &
                      vtvec%globid,vt_info(nv,ngr)%var_len_global,dsetrank,varn,nrec,irec)
                 
                 ! Create the data space for the  dataset. 
                 
                 call h5screate_simple_f(dsetrank, chnkdims, memspace, hdferr)
                 if (hdferr.ne.0) then
                    write (unit=*,fmt=*) 'Chunk dimension:  ',chnkdims
                    write (unit=*,fmt=*) 'Chunk offset:     ',chnkoffs
                    write (unit=*,fmt=*) 'Global dimension: ',globdims
                    write (unit=*,fmt=*) 'Dataset rank:     ',dsetrank
                    call fatal_error('Could not create the hyperslabs memspace' &
                                    ,'h5_output','h5_output.F90')
                 end if
                 
                 
                 ! Get the hyperslab in the file
                 
                 call h5dget_space_f(dset_id,filespace,hdferr)
                 if (hdferr /= 0) then
                    call fatal_error('Could not get the hyperslabs filespace' &
                                    ,'h5_output','h5_output.F90')
                 end if
                 
                 call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
                      cnt, hdferr, stride, chnkdims)
                 if (hdferr /= 0) then
                    call fatal_error('Could not assign the hyperslabs filespace' &
                                    ,'h5_output','h5_output.F90')
                 end if
                 
                 if (collective_mpi) then
                    
                    if (vt_info(nv,ngr)%dtype .eq. 'r') then   ! real data type
                       call h5dwrite_f(dset_id,H5T_NATIVE_REAL,vtvec%var_rp,globdims, &
                            hdferr,file_space_id = filespace, mem_space_id = memspace, &
                            xfer_prp = plist_id)
                    elseif(vt_info(nv,ngr)%dtype .eq. 'i') then ! integer data type
                       
                       call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,vtvec%var_ip,globdims, &
                            hdferr,file_space_id = filespace, mem_space_id = memspace, &
                            xfer_prp = plist_id)
                    elseif(vt_info(nv,ngr)%dtype .eq. 'c') then ! character data type
                       
                       call h5dwrite_f(dset_id,H5T_NATIVE_CHARACTER,vtvec%var_cp,globdims, &
                            hdferr,file_space_id = filespace, mem_space_id = memspace, &
                            xfer_prp = plist_id)

                    end if
                    if (hdferr /= 0) then
                       write (unit=*,fmt=*) 'Variable name:    ',varn
                       write (unit=*,fmt=*) 'Global dimension: ',globdims
                       write (unit=*,fmt=*) 'Chunk dimension:  ',chnkdims
                       write (unit=*,fmt=*) 'Chunk offset:     ',chnkoffs
                       write (unit=*,fmt=*) 'Count:            ',cnt
                       write (unit=*,fmt=*) 'Stride:           ',stride
                       call fatal_error('Could not write the real hyperslab into the dataset' &
                                       ,'h5_output','h5_output.F90')
                    end if
                    
                 else
                    
                    if (vt_info(nv,ngr)%dtype .eq. 'r') then   ! real data type
                       call h5dwrite_f(dset_id,H5T_NATIVE_REAL,vtvec%var_rp,globdims, &
                            hdferr,file_space_id = filespace, mem_space_id = memspace)
                    elseif(vt_info(nv,ngr)%dtype .eq. 'i') then ! integer data type
                       
                       call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,vtvec%var_ip,globdims, &
                            hdferr,file_space_id = filespace, mem_space_id = memspace)

                    elseif(vt_info(nv,ngr)%dtype .eq. 'c') then ! character data type
                       
                       call h5dwrite_f(dset_id,H5T_NATIVE_CHARACTER,vtvec%var_cp,globdims, &
                            hdferr,file_space_id = filespace, mem_space_id = memspace)
                    end if
                    if (hdferr /= 0) then
                       write (unit=*,fmt=*) 'Variable name:    ',varn
                       write (unit=*,fmt=*) 'Global dimension: ',globdims
                       write (unit=*,fmt=*) 'Chunk dimension:  ',chnkdims
                       write (unit=*,fmt=*) 'Chunk offset:     ',chnkoffs
                       write (unit=*,fmt=*) 'Count:            ',cnt
                       write (unit=*,fmt=*) 'Stride:           ',stride
                       call fatal_error('Could not write the real hyperslab into the dataset' &
                                       ,'h5_output','h5_output.F90')
                    end if
                 end if
                 
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ! WE SWITCHED THE ORDER OF THE FOLLOWING TWO SUBROUTINES!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 call h5sclose_f(filespace,hdferr)
                 if (hdferr /= 0) then
                    call fatal_error('Could not close the hyperslabs filespace' &
                                    ,'h5_output','h5_output.F90')
                 end if
                 
                 call h5sclose_f(memspace,hdferr)
                 if (hdferr /= 0) then
                    call fatal_error('Could not close the hyperslabs memspace' &
                                    ,'h5_output','h5_output.F90')
                 end if
                 
              end if
              
           end do pointerloop
           
           call h5dclose_f(dset_id,hdferr)
           if (hdferr /= 0) then
              call fatal_error('Could not get the dataset','h5_output','h5_output.F90')
           end if
                      
        end if
        
     end do varloop

     if (collective_mpi) then
        call h5pclose_f(plist_id,hdferr)
        if (hdferr /= 0) then
           call fatal_error('could not close the plist,post write','h5_output','h5_output.F90')
        end if
     end if
     
     call h5fclose_f(file_id,hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not close the file','h5_output','h5_output.F90')
     end if
     
     call h5close_f(hdferr)
     if (hdferr /= 0) then
        call fatal_error('Could not close the hdf environment','h5_output','h5_output.F90')
     end if
     
  enddo
    
  if (.not. collective_mpi) then
     if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,734,MPI_COMM_WORLD,ierr)
!     if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
  end if

  select case (vtype)
  case ('LITE')
     subaname='  Analysis lite HDF write'
  case ('MEAN')
     subaname='  Averaged analysis HDF write    '
  case ('BOTH')
     subaname='  Averaged analysis lite HDF write   '
  case ('DAIL')
     subaname='  Daily average analysis HDF write    '
  case ('MONT')
     subaname='  Monthly average analysis HDF write   '
  case ('YEAR')
     subaname='  Annual average analysis HDF write   '
  case ('HIST')
     subaname='  History HDF write   '
  case default
     subaname='  Analysis HDF write         '
  end select
  
  if (mynum.eq.nnodetot .and. new_file) then
     write(c0,'(F10.0)') time
     write(*,"(/,a)") " === "//trim(adjustl(subaname))//" at Sim time "//trim(adjustl(c0))//" ==="
     write(*,"(a,/)") " === wrote file "//&
          &trim(adjustl(anamel))//" ==="
  end if
     
  ! Reset the time back to the original
  if(vtype == 'MEAN'.or.vtype == 'BOTH')  time=timeold
  
  
  return
end subroutine h5_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine geth5dims(idim_type,varlen,globid,var_len_global,dsetrank,varn,nrec,irec)
  
  use grid_coms,only : nzg,nzs
  use max_dims, only : n_pft,n_dist_types,n_dbh
  use hdf5_coms,only : chnkdims,chnkoffs,cnt,stride,globdims
  use fusion_fission_coms, only: ff_ndbh

  implicit none
  character(len=*) :: varn
  integer :: idim_type,varlen,globid,var_len_global,dsetrank,nrec,irec
  
  ! Initialize the size of the memory and file-space dimensioning
  
  globdims = 0
  chnkdims = 0
  chnkoffs = 0
  cnt      = 0
  stride   = 0
  ! No assumption, I am adding it to all cases.

  ! Determine the array sizes here

  select case (idim_type) 
  case(90) ! No polygon-site-patch or cohort dimension
     
     dsetrank = 1
     chnkdims(1) = varlen
     chnkoffs(1) = globid
     globdims(1) = var_len_global
     cnt(1)      = 1
     stride(1)   = 1

  case(92) ! single vector soil info
     
     dsetrank = 2
     globdims(1) = nzg
     chnkdims(1) = nzg
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
     
  case (11) ! (npolygons) 
     
     ! Vector type,dimensions set
     
     dsetrank = 1
     chnkdims(1) = varlen
     chnkoffs(1) = globid
     globdims(1) = var_len_global
     cnt(1)      = 1
     stride(1)   = 1
     
  case (12) ! (nzg,npolygons)  
     
     ! Soil column type
     dsetrank = 2
     globdims(1) = nzg
     chnkdims(1) = nzg
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
     
  case (14) ! (n_pft,npolygons)  
     
     ! Soil column type
     dsetrank = 2
     globdims(1) = n_pft
     chnkdims(1) = n_pft
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1

  case (146) ! (n_pft,n_dbh,npolygons)
     
     dsetrank = 3
     
     globdims(1) = n_pft
     chnkdims(1) = n_pft
     chnkoffs(1) = 0
     globdims(2) = n_dbh
     chnkdims(2) = n_dbh
     chnkoffs(2) = 0
     globdims(3) = var_len_global
     chnkdims(3) = varlen
     chnkoffs(3) = globid
     cnt(1:3)      = 1
     stride(1:3)   = 1

  case (15) ! (n_dist_types,npolygons)  
     
     ! Soil column type
     dsetrank = 2
     globdims(1) = n_dist_types
     chnkdims(1) = n_dist_types
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1

  case (16) ! (n_dbh,npolygons)  
     
     ! Soil column type
     dsetrank = 2
     globdims(1) = n_dbh
     chnkdims(1) = n_dbh
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
     
  case (21) !(nsites)
     
     dsetrank = 1
     chnkdims(1) = varlen
     chnkoffs(1) = globid
     globdims(1) = var_len_global
     cnt(1)      = 1
     stride(1)   = 1

  case(22) !(nzg,nsites)
              
     ! Soil column type
     dsetrank = 2
     globdims(1) = nzg
     chnkdims(1) = nzg
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1

  case (24) !(n_pft,nsites)
     
     ! PFT type
     dsetrank = 2
     globdims(1) = n_pft
     chnkdims(1) = n_pft
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1

  case (246) !(n_pft,n_dbh,nsites)
       
     dsetrank = 3

     globdims(1) = n_pft
     chnkdims(1) = n_pft
     chnkoffs(1) = 0
     globdims(2) = n_dbh
     chnkdims(2) = n_dbh
     chnkoffs(2) = 0
     globdims(3) = var_len_global
     chnkdims(3) = varlen
     chnkoffs(3) = globid
     cnt(1:3)      = 1
     stride(1:3)   = 1
     
  case (25) !(n_dist_types,nsites)
  
     dsetrank = 2
     globdims(1) = n_dist_types
     chnkdims(1) = n_dist_types
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1

  case (255) !(n_dist_types,n_dist_types,nsites)
              
     dsetrank = 3
     
     globdims(1) = n_dist_types
     chnkdims(1) = n_dist_types
     chnkoffs(1) = 0
     globdims(2) = n_dist_types
     chnkdims(2) = n_dist_types
     chnkoffs(2) = 0
     globdims(3) = var_len_global
     chnkdims(3) = varlen
     chnkoffs(3) = globid
     cnt(1:3)      = 1
     stride(1:3)   = 1

  case (28) ! (n_months,nsites)
     dsetrank = 2
     globdims(1) = 12
     chnkdims(1) = 12
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
     
  case (31) !(npatches)
     
     dsetrank = 1            
     chnkdims(1) = varlen
     chnkoffs(1) = globid
     globdims(1) = var_len_global
     cnt(1)      = 1
     stride(1)   = 1
     
  case (32) ! (nzg,npatches)
     
     ! Soil column type
     dsetrank = 2
     globdims(1) = nzg
     chnkdims(1) = nzg
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
     
  case (33) !(nzs,npatches)
     
     ! Surface water column type
     dsetrank = 2
     globdims(1) = nzs
     chnkdims(1) = nzs
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1

  case (34) !(n_pft,npatches)
     
     ! PFT type
     dsetrank = 2
     globdims(1) = n_pft
     chnkdims(1) = n_pft
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
           
  case (347) ! (n_pft,ff_ndbh,npatched)

     dsetrank = 3
     globdims(1) = n_pft
     chnkdims(1) = n_pft
     chnkoffs(1) = 0
     
     globdims(2) = ff_ndbh
     chnkdims(2) = ff_ndbh
     chnkoffs(2) = 0

     globdims(3) = var_len_global
     chnkdims(3) = varlen
     chnkoffs(3) = globid
     cnt(1:3)      = 1
     stride(1:3)   = 1

  case (36) !(n_dbh,npatches)
     
     ! DBH type
     dsetrank = 2
     globdims(1) = n_dbh
     chnkdims(1) = n_dbh
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
              
  case (41) !(ncohorts)
     
     dsetrank = 1
     chnkdims(1) = varlen
     chnkoffs(1) = globid
     globdims(1) = var_len_global
     cnt(1)      = 1
     stride(1)   = 1
     
  case (49) !(13,ncohorts)
     
     ! 13 Months type
     dsetrank = 2
     globdims(1) = 13
     chnkdims(1) = 13
     chnkoffs(1) = 0
     globdims(2) = var_len_global
     chnkdims(2) = varlen
     chnkoffs(2) = globid
     cnt(1:2)    = 1
     stride(1:2) = 1
  case default
     write (unit=*,fmt='(a)')       '--------------------------------------------------'
     write (unit=*,fmt='(a)')       ' I can''t recognize this type of variable...'
     write (unit=*,fmt='(a,1x,a)')  ' Variable: ',varn
     write (unit=*,fmt='(a,1x,i3)') ' Dimension type: ',idim_type 
     write (unit=*,fmt='(a)')       '--------------------------------------------------'
     call fatal_error ('Wrong idim_type','geth5dims','h5_output.F90')
  end select

  !!! add TIME if writing multiple observations/file
  if(nrec .gt. 1) then
     dsetrank = dsetrank + 1
     globdims(dsetrank) = nrec
     chnkdims(dsetrank) = 1
     chnkoffs(dsetrank) = irec-1
     if(chnkoffs(dsetrank) .lt. 0) chnkoffs(dsetrank) = 0
     cnt(dsetrank)      = 1
     stride(dsetrank)   = 1
  endif

  return
end subroutine geth5dims
