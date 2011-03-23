!==========================================================================================!
!==========================================================================================!
!      This sub-routine creates the output files based on the variable table specification !
! and the type of output file to be created.  ED's default format is HDF5.                 !
!------------------------------------------------------------------------------------------!
subroutine h5_output(vtype)
#if USE_HDF5
   use hdf5
#endif
   use an_header
  
   use ed_var_tables, only : vt_info               & ! intent(in)
                           , var_table             & ! intent(in)
                           , var_table_vector      & ! intent(in)
                           , num_var               ! ! intent(in)
   use ed_misc_coms , only : ffilout               & ! intent(in)
                           , sfilout               & ! intent(in)
                           , itimea                & ! intent(in)
                           , idatea                & ! intent(in)
                           , imontha               & ! intent(in)
                           , iyeara                & ! intent(in)
                           , iclobber              & ! intent(in)
                           , nrec_fast             & ! intent(in)
                           , nrec_state            & ! intent(in)
                           , irec_fast             & ! intent(in)
                           , irec_state            & ! intent(in)
                           , out_time_fast         & ! intent(in)
                           , out_time_state        & ! intent(in)
                           , outstate              & ! intent(in)
                           , outfast               & ! intent(in)
                           , frqfast               ! ! intent(in)
   use ed_misc_coms , only : attach_metadata       ! ! intent(in)
   use grid_coms    , only : ngrids                & ! intent(in)
                           , time                  ! ! intent(in)
   use hdf5_coms    , only : chnkdims              & ! intent(in)
                           , chnkoffs              & ! intent(in)
                           , cnt                   & ! intent(in)
                           , stride                & ! intent(in)
                           , globdims              ! ! intent(in)
   use ed_node_coms , only : mynum                 & ! intent(in)
                           , nnodetot              & ! intent(in)
                           , recvnum               & ! intent(in)
                           , sendnum               ! ! intent(in)
   use ed_max_dims  , only : n_pft                 & ! intent(in)
                           , n_dist_types          & ! intent(in)
                           , n_dbh                 & ! intent(in)
                           , maxgrds               & ! intent(in)
                           , str_len               ! ! intent(in)
   use ed_state_vars, only : edgrid_g              & ! structure
                           , edtype                & ! structure
                           , polygontype           & ! structure
                           , sitetype              & ! structure
                           , patchtype             & ! structure
                           , gdpy                  ! ! intent(in)
   implicit none

   !------ Include standard common blocks. ------------------------------------------------!
   include 'mpif.h'
   !------ Arguments. ---------------------------------------------------------------------!
   character(len=*)                                 , intent(in) :: vtype
   !------ Local variables. ---------------------------------------------------------------!
   type(var_table_vector)                           , pointer    :: vtvec
   character(len=str_len)                                        :: anamel
   character(len=3)                                              :: cgrid
   character(len=40)                                             :: subaname
   character(len=64)                                             :: varn
   character(len=1)                                              :: vnam
   character(len=64)                                             :: c0
   character(len=64)     , dimension(3)                          :: metadata
   integer               ,dimension(MPI_STATUS_SIZE)             :: status
   integer                                                       :: ngr
   integer                                                       :: nv
   integer                                                       :: nvcnt
   integer                                                       :: iptr
   integer                                                       :: outyear
   integer                                                       :: outmonth
   integer                                                       :: outdate
   integer                                                       :: outhour
   integer                                                       :: irec
   integer                                                       :: nrec
   integer                                                       :: mpierror
   integer                                                       :: comm
   integer                                                       :: info
   integer                                                       :: mpi_size
   integer                                                       :: mpi_rank
   integer                                                       :: ierr
   logical                                                       :: exans
   logical                                                       :: new_file
   real(kind=8)                                                  :: dsec
   !------ HDF specific data types. -------------------------------------------------------!
   integer                                                       :: hdferr
   integer                                                       :: dsetrank
   integer(HID_T)                                                :: file_id
   integer(HID_T)                                                :: dset_id
   integer(HID_T)                                                :: memspace
   integer(HID_T)                                                :: filespace
   integer(HID_T)                                                :: attr_id
   !------ Attribute types. ---------------------------------------------------------------!
   integer(HID_T)                                                :: aspace_id
   integer(HID_T)                                                :: atype_id
   integer(HSIZE_T)      , dimension(1)                          :: adims
   integer                                                       :: arank
   integer(SIZE_T)                                               :: attrlen
   !----- Locally saved variables. --------------------------------------------------------!
   integer                                          , save       :: irec_opt       = 0
   !----- Local constants. ----------------------------------------------------------------!
   logical                                          , parameter  :: collective_mpi = .false.
   real(kind=8)                                     , parameter  :: zero           = 0.0d0
   !----- External functions. -------------------------------------------------------------!
   logical                                          , external   :: isleap
   !---------------------------------------------------------------------------------------!



   !----- Initialise some variables. ------------------------------------------------------!
   comm    = MPI_COMM_WORLD
   info    = MPI_INFO_NULL
   !---------------------------------------------------------------------------------------!



   !------ Find which letter we should use to denote this type of analysis. ---------------!

   select case (trim(vtype))
   case ('INST')
      vnam='I'   ! Instantaneous analysis.
   case ('DAIL') 
      vnam='D'   ! Daily averages.
   case ('MONT') 
      vnam='E'   ! Monthly averages.
   case ('DCYC')
      vnam='Q'   ! Mean diurnal cycle
   case ('YEAR') 
      vnam='Y'   ! Yearly averages.
   case ('OPTI') 
      vnam='T'   ! Ouput for optimisation.
   case ('HIST') 
      vnam='S'   ! S for reStart, don't want confusion with RAMS' H or R files
   case ('CONT') 
      vnam='Z'   ! The first time with history start, so we don't replace the history
   end select  
   nvcnt=0


   nrec = 1
   irec = 1  

   !---------------------------------------------------------------------------------------!
   !    Loop over the grids.                                                               !
   !---------------------------------------------------------------------------------------!
   gridloop: do ngr=1,ngrids

      !----- I guess this cleans out anything that didn't finish correctly. ---------------!
      call h5garbage_collect_f(hdferr) 


      !----- Wait until the previous node has finished writing. ---------------------------!
      new_file=.true.
      if (mynum /= 1) then
         call MPI_RECV(new_file,1,MPI_LOGICAL,recvnum,3510+ngr,MPI_COMM_WORLD,status,ierr)
      end if

      !------------------------------------------------------------------------------------!
      !    If there are no polygons on this node, then the node shouldn't do anything to   !
      ! this file.  This is something that happens only in coupled model runs.             !
      !------------------------------------------------------------------------------------!
      if (gdpy(mynum,ngr)>0) then
         !----- Make the grid flag. -------------------------------------------------------!
         write(cgrid,fmt='(a1,i2.2)') 'g',ngr


         !---------------------------------------------------------------------------------!
         !     Create the file name.                                                       !
         !---------------------------------------------------------------------------------!
         select case (trim(vtype))
         case ('DAIL')
            !----- Return the year, month, and day of the previous day. -------------------!
            call date_add_to(iyeara,imontha,idatea,itimea*100,time-21600.0d0,'s',outyear   &
                            ,outmonth,outdate,outhour)
            !------------------------------------------------------------------------------!
            !     Set outhour to be zero, it was pushed backwards just to get the previous !
            ! day.                                                                         !
            !------------------------------------------------------------------------------!
            outhour = 0
            !----- Make the file name. ----------------------------------------------------!
            call makefnam(anamel,ffilout,zero,outyear,outmonth,outdate,outhour,vnam,cgrid  &
                         ,'h5 ')

         case ('MONT','DCYC')
            !----- Return the year, month, and day of the previous day. -------------------!
            call date_add_to(iyeara,imontha,idatea,itimea*100,time-21600.0d0,'s',outyear   &
                            ,outmonth,outdate,outhour)
            !------------------------------------------------------------------------------!
            !     Set outhour and outdate to zero.                                         !
            !------------------------------------------------------------------------------!
            outhour = 0
            outdate = 0
            !----- Make the file name. ----------------------------------------------------!
            call makefnam(anamel,ffilout,zero,outyear,outmonth,outdate,outhour,vnam,cgrid   &
                         ,'h5 ')

         case ('YEAR')
            !----- Return the year, month, and day of the previous day. -------------------!
            call date_add_to(iyeara,imontha,idatea,itimea*100,time-21600.0d0,'s',outyear   &
                            ,outmonth,outdate,outhour)
            !------------------------------------------------------------------------------!
            !     Set outhour, outdate, and outmonth to zero.                              !
            !------------------------------------------------------------------------------!
            outhour  = 0
            outdate  = 0
            outmonth = 0
            !----- Make the file name. ----------------------------------------------------!
            call makefnam(anamel,ffilout,zero,outyear,outmonth,outdate,outhour,vnam,cgrid  &
                         ,'h5 ')

         case ('OPTI')
            !------ Make the optimisation file. -------------------------------------------!
            new_file = .false.
            irec_opt = irec_opt + 1
            call date_add_to(iyeara,imontha,idatea,itimea*100,time-21600.0d0,'s',outyear   &
                            ,outmonth,outdate,outhour)

            if (isleap(outyear)) then
               nrec = int(366*86400/frqfast)
            else
               nrec = int(365*86400/frqfast)
            end if

            call date_add_to(iyeara,imontha,idatea,itimea*100,time-21600.0d0,'s',outyear   &
                            ,outmonth,outdate,outhour)
           
            call makefnam(anamel,ffilout,zero,outyear,0,0,0,vnam,cgrid,'h5 ')

            if (outmonth == 1 .and. outdate == 1 .and. outhour == 0) then
               new_file = .true.
               irec_opt = 1
            end if
            irec  = irec_opt
            if (irec == 1) new_file = .true.

         case ('HIST','CONT')
         
            call makefnam(anamel,sfilout,time,iyeara,imontha,idatea,itimea*100,vnam,cgrid  &
                         ,'h5 ')
            outyear  = iyeara
            outmonth = imontha
            outdate  = idatea
            outhour  = itimea*100

         case default
            if (nrec_fast == 1) then
               !----- Single file per output. ---------------------------------------------!
               call makefnam(anamel,ffilout,time,iyeara,imontha,idatea,itimea*100,vnam     &
                            ,cgrid,'h5 ')
            else   
               !----- Multiple times per file. --------------------------------------------!
               new_file = .false.
               !----- Determine whether to advance out_time_fast. -------------------------!
               call date_add_to(iyeara,imontha,idatea,itimea*100,time,'s',outyear,outmonth &
                               ,outdate,outhour)
               call date_2_seconds(out_time_fast%year,out_time_fast%month                  &
                                  ,out_time_fast%date,int(out_time_fast%time)              &
                                  ,iyeara,imontha,idatea,itimea*100,dsec)

               if (time >= (dsec+dble(outfast)) .or. outmonth /= out_time_fast%month) then
                  out_time_fast%year  = outyear
                  out_time_fast%month = outmonth
                  out_time_fast%date  = outdate
                  out_time_fast%time  = real(outhour)
                  dsec                = time
                  new_file            = .true.
               end if

               irec_fast = int(((time-dsec)/dble(frqfast))) + 1
               nrec      = nrec_fast
               irec      = irec_fast

               call makefnam(anamel,ffilout,zero,out_time_fast%year,out_time_fast%month    &
                            ,out_time_fast%date,int(out_time_fast%time),vnam,cgrid,'h5 ')
            end if           
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check whether the file exists, and stop the run if the user doesn't want    !
         ! the files to be over-written.                                                   !
         !---------------------------------------------------------------------------------!
         inquire (file=trim(anamel),exist=exans)
         if (exans .and. iclobber == 0) then
            write(unit=*,fmt='(a)'      )  '-----------------------------------------------'
            write(unit=*,fmt='(3(a,1x))') ' - File:',trim(anamel),'already exists.'
            write(unit=*,fmt='(a)'      ) ' - Based on your ICLOBBER, you don''t want'
            write(unit=*,fmt='(a)'      ) '   files to be over-written.  Stopping the run.'
            write(unit=*,fmt='(a)'      ) '-----------------------------------------------'
            call fatal_error('Attempted to over-write a file that already exists.'         &
                            ,'h5_output','h5_output.F90')
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Initialise the HDF5 environment.                                            !
         !---------------------------------------------------------------------------------!
         call h5open_f(hdferr)
         if (hdferr /= 0) then
            write(unit=*,fmt='(a,1x,i)') ' - HDF5 Open error #:',hdferr
            call fatal_error('Could not initialize the hdf environment'                    &
                            ,'h5_output','h5_output.F90')
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Either open or create the output file.                                      !
         !---------------------------------------------------------------------------------!
         if (new_file) then
            
            call h5fcreate_f(trim(anamel)//char(0), H5F_ACC_TRUNC_F, file_id, hdferr)
            if (hdferr /= 0) then
               write(unit=*,fmt='(a)'      ) '--------------------------------------------'
               write(unit=*,fmt='(a)'      ) ' Could not create the HDF5 file.'
               write(unit=*,fmt='(a,1x,i)' ) ' - File   : ',trim(anamel)
               write(unit=*,fmt='(a,1x,i)' ) ' - file_id: ',file_id
               write(unit=*,fmt='(a,1x,i)' ) ' - hdferr : ',hdferr
               write(unit=*,fmt='(a)'      ) '--------------------------------------------'
               call fatal_error('Failed creating the HDF file','h5_output','h5_output.F90')
            end if
         else
            call h5fopen_f(trim(anamel)//char(0), H5F_ACC_RDWR_F, file_id, hdferr)
            if (hdferr /= 0) then
               write(unit=*,fmt='(a)'      ) '--------------------------------------------'
               write(unit=*,fmt='(a)'      ) ' Could not open the HDF5 file.'
               write(unit=*,fmt='(a,1x,i)' ) ' - File   : ',trim(anamel)
               write(unit=*,fmt='(a,1x,i)' ) ' - file_id: ',file_id
               write(unit=*,fmt='(a,1x,i)' ) ' - hdferr : ',hdferr
               write(unit=*,fmt='(a)'      ) '--------------------------------------------'
               call fatal_error('Failed opening the HDF file','h5_output','h5_output.F90')
            end if
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !   Create HDF5 datasets and then put them in the file; we must loop over all     !
         ! variables.                                                                      !
         !---------------------------------------------------------------------------------!
         varloop: do nv = 1,num_var(ngr)
            !----- Check whether the variable goes to the output. -------------------------!
            if ((vtype == 'INST' .and. vt_info(nv,ngr)%ianal == 1) .or.                    &
                (vtype == 'OPTI' .and. vt_info(nv,ngr)%iopti == 1) .or.                    &
                (vtype == 'LITE' .and. vt_info(nv,ngr)%ilite == 1) .or.                    &
                (vtype == 'DAIL' .and. vt_info(nv,ngr)%idail == 1) .or.                    &
                (vtype == 'MONT' .and. vt_info(nv,ngr)%imont == 1) .or.                    &
                (vtype == 'DCYC' .and. vt_info(nv,ngr)%idcyc == 1) .or.                    &
                (vtype == 'YEAR' .and. vt_info(nv,ngr)%iyear == 1) .or.                    &
                (vtype == 'HIST' .and. vt_info(nv,ngr)%ihist == 1) .or.                    &
                (vtype == 'CONT' .and. vt_info(nv,ngr)%ihist == 1) ) then

               !----- Variable name. ------------------------------------------------------!
               varn= vt_info(nv,ngr)%name
               !---------------------------------------------------------------------------!

               !----- Initialize global dimensions of the hyperslab. ----------------------!
               call geth5dims(vt_info(nv,ngr)%idim_type,0,0,vt_info(nv,ngr)%var_len_global &
                             ,dsetrank,varn,nrec,irec)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Create the data set.                                                  !
               !---------------------------------------------------------------------------!
               call h5screate_simple_f(dsetrank, globdims, filespace, hdferr)
               if (hdferr /= 0 .or. globdims(1) < 1 ) then
                  write (unit=*,fmt='(a,1x,a)') ' VTYPE:    ',trim(vtype)
                  write (unit=*,fmt='(a,1x,a)') ' VAR NAME: ',trim(varn)
                  write (unit=*,fmt='(a,1x,i)') ' IDIM_TYPE:',vt_info(nv,ngr)%idim_type
                  write (unit=*,fmt='(a,1x,i)') ' VLEN_GLOB:',vt_info(nv,ngr)%var_len_global
                  write (unit=*,fmt=*)          ' DSETRANK: ',dsetrank
                  write (unit=*,fmt=*)          ' GLOBDIMS: ',globdims
                  call fatal_error('Could not create the first filespace'                  &
                                  ,'h5_output','h5_output.f90')
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Determine whether the dataset exists.                                !
               !---------------------------------------------------------------------------!
               call h5eset_auto_f(0,hdferr)
               call h5dopen_f(file_id,varn,dset_id,hdferr)

               if (hdferr < 0) then
                  call h5eset_auto_f(1,hdferr)

                  select case (vt_info(nv,ngr)%dtype)
                  case ('R','r') !----- Real variable. ------------------------------------!
                     call h5dcreate_f(file_id,varn,H5T_NATIVE_REAL, filespace,dset_id      &
                                     ,hdferr)

                  case ('D','d') !----- Double precision variable. ------------------------!
                     call h5dcreate_f(file_id,varn,H5T_NATIVE_DOUBLE, filespace,dset_id    &
                                     ,hdferr)

                  case ('I','i') !----- Integer variable. ---------------------------------!
                     call h5dcreate_f(file_id,varn,H5T_NATIVE_INTEGER, filespace,dset_id   &
                                     ,hdferr)

                  case ('C','c') !----- Character variable. -------------------------------!
                     call h5dcreate_f(file_id,varn,H5T_NATIVE_CHARACTER, filespace,dset_id &
                                     ,hdferr)

                  case default
                     write(unit=*,fmt='(a)')      '---------------------------------------'
                     write(unit=*,fmt='(a,1x,a)') ' Variable: ',trim(varn)
                     write(unit=*,fmt='(a,1x,a)') ' Type:     ',trim(vt_info(nv,ngr)%dtype)
                     write(unit=*,fmt='(a)')      '---------------------------------------'
                     call fatal_error('Attempted to write an undefined datatype'           &
                                     ,'h5_output','h5_output.F90')

                  end select



                  !------------------------------------------------------------------------!
                  !      Attached metadata if the user wants it.                           !
                  !------------------------------------------------------------------------!
                  if (attach_metadata == 1) then
                     arank       = 1
                     adims       = 3_8
                     attrlen     = 64_8
                     metadata(1) = trim('Long Name: '//trim(vt_info(nv,ngr)%lname))
                     metadata(2) = trim('Units: '//trim(vt_info(nv,ngr)%units))
                     metadata(3) = trim('Dimensions: '//trim(vt_info(nv,ngr)%dimlab))

                     call h5screate_simple_f(arank,adims,aspace_id,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Error calling h5screate_simple_f'                &
                                        ,'h5_output','h5_output.F90')
                     end if

                     call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Error calling h5tcopy_f'                         &
                                        ,'h5_output','h5_output.F90')
                     end if

                     call h5tset_size_f(atype_id,attrlen,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Error calling h5tset_size_f'                     &
                                        ,'h5_output','h5_output.F90')
                     end if

                     call h5acreate_f(dset_id,'Metadata',atype_id,aspace_id,attr_id,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Error calling h5acreate_f'                       &
                                        ,'h5_output','h5_output.F90')
                     end if

                     call h5awrite_f(attr_id,atype_id,metadata,adims,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Error calling h5awrite_f'                        &
                             ,'h5_output','h5_output.F90')
                     end if

                     call h5aclose_f(attr_id,hdferr)
                     call h5sclose_f(aspace_id,hdferr)
                  end if

                  call h5dopen_f(file_id,varn,dset_id,hdferr)

                  if (hdferr /= 0) then
                     write (unit=*,fmt=*) 'File name:           ',trim(anamel)
                     write (unit=*,fmt=*) 'Variable name:       ',trim(varn)
                     write (unit=*,fmt=*) 'File ID:             ',file_id
                     write (unit=*,fmt=*) 'Dataset ID:          ',dset_id
                     write (unit=*,fmt=*) 'Dataset rank:        ',dsetrank
                     write (unit=*,fmt=*) 'Global dimension:    ',globdims
                     write (unit=*,fmt=*) 'Chunk size           ',chnkdims
                     write (unit=*,fmt=*) 'Chunk offset         ',chnkoffs
                     write (unit=*,fmt=*) 'Vars written so far: ',nv-1
                     call fatal_error('Could not create the dataset'                       &
                                     ,'h5_output','h5_output.F90')
                  end if
               else
                  call h5eset_auto_f(1,hdferr)
               end if

               call h5sclose_f(filespace,hdferr)
               if (hdferr /= 0) then
                  call fatal_error('Could not close the first filespace'                   &
                                  ,'h5_output','h5_output.F90')
               end if
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Loop over all the pointers.                                          !
               !---------------------------------------------------------------------------!
               pointerloop: do iptr = 1,vt_info(nv,ngr)%nptrs
                  
                  vtvec => vt_info(nv,ngr)%vt_vector(iptr)
                  !------------------------------------------------------------------------!
                  !      Set the size of the chunk and it's offset in the global dataset.  !
                  !------------------------------------------------------------------------!
                  if (vtvec%varlen > 0 ) then
                     !---------------------------------------------------------------------!
                     !     Evaluate the variable output type.  Resolve the dimensioning    !
                     ! and the meta-data tags accordingly.  See ed_state_vars.f90 for a    !
                     ! description of the various datatype.                                !
                     !---------------------------------------------------------------------!
                     !----- Initialize hyperslab indices. ---------------------------------!
                     call geth5dims(vt_info(nv,ngr)%idim_type,vtvec%varlen,vtvec%globid    &
                                   ,vt_info(nv,ngr)%var_len_global,dsetrank,varn,nrec,irec)
                     
                     !----- Create the data space for the dataset. ------------------------!
                     call h5screate_simple_f(dsetrank, chnkdims, memspace, hdferr)
                     if (hdferr /= 0) then
                        write (unit=*,fmt=*) 'Chunk dimension:  ',chnkdims
                        write (unit=*,fmt=*) 'Chunk offset:     ',chnkoffs
                        write (unit=*,fmt=*) 'Global dimension: ',globdims
                        write (unit=*,fmt=*) 'Dataset rank:     ',dsetrank
                        call fatal_error('Could not create the hyperslabs memspace'        &
                                        ,'h5_output','h5_output.F90')
                     end if

                     !----- Get the hyperslab in the file. --------------------------------!
                     call h5dget_space_f(dset_id,filespace,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Could not get the hyperslabs filespace'          &
                                        ,'h5_output','h5_output.F90')
                     end if

                     call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs,cnt    &
                                               ,hdferr,stride,chnkdims)
                     if (hdferr /= 0) then
                        call fatal_error('Could not assign the hyperslabs filespace'       &
                                        ,'h5_output','h5_output.F90')
                     end if


                     !---------------------------------------------------------------------!
                     !     Choose the right pointer when writing the variable.             !
                     !---------------------------------------------------------------------!
                     select case (vt_info(nv,ngr)%dtype)
                     case ('R') !----- Real variable (vector). ----------------------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_REAL,vtvec%var_rp,globdims      &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     case ('D') !----- Double precision variable (vector). ----------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vtvec%var_dp,globdims    &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     case ('I') !----- Integer variable (vector). -------------------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,vtvec%var_ip,globdims   &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     case ('C') !----- Character variable (vector). -----------------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_CHARACTER,vtvec%var_cp,globdims &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id = memspace)

                     case ('r') !----- Real variable (scalar). ----------------------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_REAL,vtvec%sca_rp,globdims      &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     case ('d') !----- Double precision variable (scalar). ----------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,vtvec%sca_dp,globdims    &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     case ('i') !----- Integer variable (scalar). -------------------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,vtvec%sca_ip,globdims   &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     case ('c') !----- Character variable (scalar). -----------------------!
                        call h5dwrite_f(dset_id,H5T_NATIVE_CHARACTER,vtvec%sca_cp,globdims &
                                       ,hdferr,file_space_id=filespace                     &
                                       ,mem_space_id=memspace)

                     end select
                     if (hdferr /= 0) then
                        write (unit=*,fmt=*) 'Variable name:    ',varn
                        write (unit=*,fmt=*) 'Global dimension: ',globdims
                        write (unit=*,fmt=*) 'Chunk dimension:  ',chnkdims
                        write (unit=*,fmt=*) 'Chunk offset:     ',chnkoffs
                        write (unit=*,fmt=*) 'Count:            ',cnt
                        write (unit=*,fmt=*) 'Stride:           ',stride
                        call fatal_error('Could not write hyperslab into the dataset'      &
                                        ,'h5_output','h5_output.F90')
                     end if


                     call h5sclose_f(filespace,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Could not close the hyperslabs filespace'        &
                                        ,'h5_output','h5_output.F90')
                     end if

                     call h5sclose_f(memspace,hdferr)
                     if (hdferr /= 0) then
                        call fatal_error('Could not close the hyperslabs memspace'         &
                                        ,'h5_output','h5_output.F90')
                     end if
                  end if
               end do pointerloop
               !---------------------------------------------------------------------------!

               call h5dclose_f(dset_id,hdferr)
               if (hdferr /= 0) then
                  call fatal_error('Could not get the dataset','h5_output','h5_output.F90')
               end if
              
            end if
         end do varloop
         !---------------------------------------------------------------------------------!

         call h5fclose_f(file_id,hdferr)
         if (hdferr /= 0) then
            call fatal_error('Could not close the file','h5_output','h5_output.F90')
         end if
        
         call h5close_f(hdferr)
         if (hdferr /= 0) then
            call fatal_error('Could not close the hdf environment'                         &
                            ,'h5_output','h5_output.F90')
         end if

         new_file = .false.
      end if
      !------------------------------------------------------------------------------------!


      !----- Send message for the next node that it is time for it to write the output. ---!
      if (mynum < nnodetot) then
         call MPI_Send(new_file,1,MPI_LOGICAL,sendnum,3510+ngr,MPI_COMM_WORLD,ierr)
      end if
   end do gridloop

   select case (trim(vtype))
   case ('DAIL')
      subaname = '  Daily average analysis HDF write    '

   case ('MONT')
      subaname = '  Monthly average analysis HDF write   '

   case ('DCYC')
      subaname = '  Mean diurnal cycle analysis HDF write   '

   case ('YEAR')
      subaname = '  Annual average analysis HDF write   '

   case ('OPTI')
      subaname = '  Opt update   '

   case ('HIST','CONT')
      !------------------------------------------------------------------------------------!
      !      The following is helpful to those running very massive simulations.  File     !
      ! write times may be significantly long in these cases, and with automation of       !
      ! simulation queuing, the last thing you want to it try restarting partially written !
      ! files.  The generation of the CPM file ensures that its sibling HDF5 file was      !
      ! completely written.  For simplicity, this will remain commented out by default.    !
      !------------------------------------------------------------------------------------!
      if (nnodetot /= 1 .and. mynum == nnodetot) then
         !----- Write a dummy file that signals we are done. ------------------------------!
         call makefnam(anamel,sfilout,time,iyeara,imontha,idatea,itimea*100,vnam,'g00'     &
                      ,'cmp')
         open (unit=79,file=trim(anamel),form='formatted',status='unknown')
         write(unit=79,fmt='(a)') 'history write completed'
         write(unit=* ,fmt=*    ) 'Completed History Write: ',trim(anamel)
         close(unit=79,status='keep')
      end if

      select case(trim(vtype))
      case ('HIST')
         subaname = '  History HDF write   '
      case ('CONT')
         subaname = '  Resume file HDF write   '
      end select
   case default
      subaname = '  Analysis HDF write         '
   end select
     
   if (mynum == nnodetot .and. new_file) then
      write(c0,fmt='(f15.0)') sngl(time)
      write(unit=*,fmt='(a)') ' === '//trim(adjustl(subaname))//' at Sim time '            &
                                     //trim(adjustl(c0))//' ==='
   end if

   return
end subroutine h5_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine geth5dims(idim_type,varlen,globid,var_len_global,dsetrank,varn,nrec,irec)
   
   use grid_coms          , only : nzg            & ! intent(in)
                                 , nzs            ! ! intent(in)
   use ed_max_dims        , only : n_pft          & ! intent(in)
                                 , n_dist_types   & ! intent(in)
                                 , n_dbh          & ! intent(in)
                                 , n_age          & ! intent(in)
                                 , n_mort         ! ! intent(in)
   use hdf5_coms          , only : chnkdims       & ! intent(in)
                                 , chnkoffs       & ! intent(in)
                                 , cnt            & ! intent(in)
                                 , stride         & ! intent(in)
                                 , globdims       ! ! intent(in)
   use fusion_fission_coms, only : ff_ndbh        ! ! intent(in)
   use c34constants       , only : n_stoma_atts   ! ! intent(in)
   use ed_misc_coms       , only : ndcycle        ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*), intent(in)  :: varn
   integer         , intent(in)  :: idim_type
   integer         , intent(in)  :: varlen
   integer         , intent(in)  :: globid
   integer         , intent(in)  :: var_len_global
   integer         , intent(out) :: dsetrank
   integer         , intent(in)  :: nrec
   integer         , intent(in)  :: irec
   !---------------------------------------------------------------------------------------!


   !----- Initialize the size of the memory and file-space dimensioning. ------------------!
   globdims = 0_8
   chnkdims = 0_8
   chnkoffs = 0_8
   cnt      = 0_8
   stride   = 0_8

   !---------------------------------------------------------------------------------------!
   !     Determine the array sizes based on idim_type.                                     !
   !---------------------------------------------------------------------------------------!

   select case (idim_type) 
   case(90,92) ! No polygon-site-patch or cohort dimension, or a single-dimension vector
      
      dsetrank = 1
      chnkdims(1) = int(varlen,8)
      chnkoffs(1) = int(globid,8)
      globdims(1) = int(var_len_global,8)
      cnt(1)      = 1_8
      stride(1)   = 1_8
      
   case (10,11) ! (npolygons) 
      
      ! Vector type,dimensions set
      
      dsetrank = 1
      chnkdims(1) = int(varlen,8)
      chnkoffs(1) = int(globid,8)
      globdims(1) = int(var_len_global,8)
      cnt(1)      = 1_8
      stride(1)   = 1_8
      
   case (-11) ! (ndcycle,npolygons)  
      
      ! Soil column type
      dsetrank = 2
      globdims(1) = int(ndcycle,8)
      chnkdims(1) = int(ndcycle,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
      
   case (12,120) ! (nzg,npolygons)  
      
      ! Soil column type
      dsetrank = 2
      globdims(1) = int(nzg,8)
      chnkdims(1) = int(nzg,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (-12) ! (nzg,ndcycle,npolygons)
      
      dsetrank = 3
      
      globdims(1) = int(nzg,8)
      chnkdims(1) = int(nzg,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(ndcycle,8)
      chnkdims(2) = int(ndcycle,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)      = 1_8
      stride(1:3)   = 1_8
      
   case (14) ! (n_pft,npolygons)  
      
      dsetrank = 2
      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (146) ! (n_pft,n_dbh,npolygons)
      
      dsetrank = 3
      
      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(n_dbh,8)
      chnkdims(2) = int(n_dbh,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)      = 1_8
      stride(1:3)   = 1_8

   case (14567) ! (n_pft,n_dist_types,n_dbh,n_age,npolygons)
      
      dsetrank = 5
      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(n_dist_types,8)
      chnkdims(2) = int(n_dist_types,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(n_dbh,8)
      chnkdims(3) = int(n_dbh,8)
      chnkoffs(3) = 0_8
      globdims(4) = int(n_age,8)
      chnkdims(4) = int(n_age,8)
      chnkoffs(4) = 0_8
      globdims(5) = int(var_len_global,8)
      chnkdims(5) = int(varlen,8)
      chnkoffs(5) = int(globid,8)
      cnt(1:5)      = 1_8
      stride(1:5)   = 1_8

   case (15) ! (n_dist_types,npolygons)  
      
      dsetrank = 2
      globdims(1) = int(n_dist_types,8)
      chnkdims(1) = int(n_dist_types,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (155) ! (n_dist_types,n_dist_types,npolygons)  
      
      dsetrank = 3
      globdims(1) = int(n_dist_types,8)
      chnkdims(1) = int(n_dist_types,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(n_dist_types,8)
      chnkdims(2) = int(n_dist_types,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)    = 1_8
      stride(1:3) = 1_8

   case (157) ! (n_dist_types,n_age,n_age)  

      dsetrank = 3
      globdims(1) = int(n_dist_types,8)
      chnkdims(1) = int(n_dist_types,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(n_age,8)
      chnkdims(2) = int(n_age,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)    = 1_8
      stride(1:3) = 1_8

   case (16) ! (n_dbh,npolygons)  
      
      dsetrank = 2
      globdims(1) = int(n_dbh,8)
      chnkdims(1) = int(n_dbh,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (17) ! (n_age,npolygons)  
      
      dsetrank = 2
      globdims(1) = int(n_age,8)
      chnkdims(1) = int(n_age,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (18) ! (n_mort,npolygons)  
      
      dsetrank = 2
      globdims(1) = int(n_mort,8)
      chnkdims(1) = int(n_mort,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (19) ! (13 months,npolygons)  
      
      dsetrank = 2
      globdims(1) = int(13,8)
      chnkdims(1) = int(13,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (199) ! (nlai,nvars,npolygons)  
      
      dsetrank = 3
      globdims(1) = 3_8
      chnkdims(1) = 3_8
      chnkoffs(1) = 0_8
      globdims(2) = 4_8
      chnkdims(2) = 4_8
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)    = 1_8
      stride(1:3) = 1_8

   case (20,21) !(nsites)
      
      dsetrank = 1
      chnkdims(1) = int(varlen,8)
      chnkoffs(1) = int(globid,8)
      globdims(1) = int(var_len_global,8)
      cnt(1)      = 1_8
      stride(1)   = 1_8

   case(22,220) !(nzg,nsites)
               
      ! Soil column type
      dsetrank = 2
      globdims(1) = int(nzg,8)
      chnkdims(1) = int(nzg,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (24) !(n_pft,nsites)
      
      ! PFT type
      dsetrank = 2
      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (246) !(n_pft,n_dbh,nsites)
        
      dsetrank = 3

      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(n_dbh,8)
      chnkdims(2) = int(n_dbh,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)      = 1_8
      stride(1:3)   = 1_8
      
   case (25) !(n_dist_types,nsites)
   
      dsetrank = 2
      globdims(1) = int(n_dist_types,8)
      chnkdims(1) = int(n_dist_types,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (255) !(n_dist_types,n_dist_types,nsites)
               
      dsetrank = 3
      
      globdims(1) = int(n_dist_types,8)
      chnkdims(1) = int(n_dist_types,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(n_dist_types,8)
      chnkdims(2) = int(n_dist_types,8)
      chnkoffs(2) = 0_8
      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)      = 1_8
      stride(1:3)   = 1_8
      
   case (26) !(n_dbh,nsites)
   
      dsetrank = 2
      globdims(1) = int(n_dbh,8)
      chnkdims(1) = int(n_dbh,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
       
   case (27) !(n_age,nsites)
   
      dsetrank = 2
      globdims(1) = int(n_age,8)
      chnkdims(1) = int(n_age,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
       
   case (28) !(n_mort,nsites)
   
      dsetrank = 2
      globdims(1) = int(n_mort,8)
      chnkdims(1) = int(n_mort,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (29) ! (n_months,nsites)
      dsetrank = 2
      globdims(1) = 12_8
      chnkdims(1) = 12_8
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
     
   case (30,31) !(npatches)
      
      dsetrank = 1            
      chnkdims(1) = int(varlen,8)
      chnkoffs(1) = int(globid,8)
      globdims(1) = int(var_len_global,8)
      cnt(1)      = 1_8
      stride(1)   = 1_8
      
   case (-31) ! (ndcycle,npatches)
      
      ! Soil column type
      dsetrank = 2
      globdims(1) = int(ndcycle,8)
      chnkdims(1) = int(ndcycle,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
       
   case (32,320) ! (nzg,npatches)
      
      ! Soil column type
      dsetrank = 2
      globdims(1) = int(nzg,8)
      chnkdims(1) = int(nzg,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
     
   case (33) !(nzs,npatches)
      
      ! Surface water column type
      dsetrank = 2
      globdims(1) = int(nzs,8)
      chnkdims(1) = int(nzs,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (34) !(n_pft,npatches)
      
      ! PFT type
      dsetrank = 2
      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
            
   case (346) ! (n_pft,ff_ndbh,npatched)

      dsetrank = 3
      globdims(1) = int(n_pft,8)
      chnkdims(1) = int(n_pft,8)
      chnkoffs(1) = 0_8
      
      globdims(2) = int(ff_ndbh,8)
      chnkdims(2) = int(ff_ndbh,8)
      chnkoffs(2) = 0_8

      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)      = 1_8
      stride(1:3)   = 1_8

  case (316) ! (n_stoma_atts,n_pft,npatches)

      dsetrank = 3
      globdims(1) = int(n_stoma_atts,8)
      chnkdims(1) = int(n_stoma_atts,8)
      chnkoffs(1) = 0_8
      
      globdims(2) = int(n_pft,8)
      chnkdims(2) = int(n_pft,8)
      chnkoffs(2) = 0_8

      globdims(3) = int(var_len_global,8)
      chnkdims(3) = int(varlen,8)
      chnkoffs(3) = int(globid,8)
      cnt(1:3)      = 1_8
      stride(1:3)   = 1_8
      

   case (36) !(n_dbh,npatches)
      
      ! DBH type
      dsetrank = 2
      globdims(1) = int(n_dbh,8)
      chnkdims(1) = int(n_dbh,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
  
   case (37) !(n_age,npatches)
      
      ! DBH type
      dsetrank = 2
      globdims(1) = int(n_age,8)
      chnkdims(1) = int(n_age,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
  
   case (38) !(n_mort,npatches)
      
      ! Mortality type
      dsetrank = 2
      globdims(1) = int(n_mort,8)
      chnkdims(1) = int(n_mort,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
              
   case (40,41) !(ncohorts)
      
      dsetrank = 1
      chnkdims(1) = int(varlen,8)
      chnkoffs(1) = int(globid,8)
      globdims(1) = int(var_len_global,8)
      cnt(1)      = 1_8
      stride(1)   = 1_8

   case (-41) !(ndcycle,ncohorts)
      
      dsetrank = 2
      globdims(1) = int(ndcycle,8)
      chnkdims(1) = int(ndcycle,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case (416) !(16 - ncohorts (stoma data))
      
      dsetrank = 2
      globdims(1) = 16_8
      chnkdims(1) = 16_8
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
     
   case (46) !(n_dbh,ncohorts)
      
      ! DBH class type
      dsetrank = 2
      globdims(1) = int(n_dbh,8)
      chnkdims(1) = int(n_dbh,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
     
   case (47) !(n_age,ncohorts)
      
      ! Age class type
      dsetrank = 2
      globdims(1) = int(n_age,8)
      chnkdims(1) = int(n_age,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
     
   case (48) !(n_mort,ncohorts)
      
      ! Mortality type
      dsetrank    = 2
      globdims(1) = int(n_mort,8)
      chnkdims(1) = int(n_mort,8)
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8
     
   case (49) !(13,ncohorts)
      
      ! 13 Months type
      dsetrank = 2
      globdims(1) = 13_8
      chnkdims(1) = 13_8
      chnkoffs(1) = 0_8
      globdims(2) = int(var_len_global,8)
      chnkdims(2) = int(varlen,8)
      chnkoffs(2) = int(globid,8)
      cnt(1:2)    = 1_8
      stride(1:2) = 1_8

   case default
      write (unit=*,fmt='(a)')       '--------------------------------------------------'
      write (unit=*,fmt='(a)')       ' I can''t recognize this type of variable...'
      write (unit=*,fmt='(a,1x,a)')  ' Variable: ',varn
      write (unit=*,fmt='(a,1x,i3)') ' Dimension type: ',idim_type 
      write (unit=*,fmt='(a)')       '--------------------------------------------------'
      call fatal_error ('Wrong idim_type','geth5dims','h5_output.F90')
   end select

   !!! add TIME if writing multiple observations/file
   if(nrec > 1) then
      dsetrank = dsetrank + 1
      globdims(dsetrank) = int(nrec,8)
      chnkdims(dsetrank) = 1_8
      chnkoffs(dsetrank) = int(irec-1,8)
      if(chnkoffs(dsetrank) .lt. 0_8) chnkoffs(dsetrank) = 0_8
      cnt(dsetrank)      = 1_8
      stride(dsetrank)   = 1_8
   endif

   return
end subroutine geth5dims
