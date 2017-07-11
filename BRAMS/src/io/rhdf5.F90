!==========================================================================================!
!==========================================================================================!
!     This subroutine creates HDF-5 output files from the var table.                       !
!------------------------------------------------------------------------------------------!
subroutine anlhdf(vtype)

   use an_header
   use var_tables
   use grid_dims , only : maxgrds & ! intent(in)
                        , str_len ! ! intent(in)
   use mem_aerad , only : nwave   ! ! intent(in)
   use mem_cuparm, only : nclouds ! ! intent(in)

   use mem_grid
   use io_params

#if USE_HDF5
  use hdf5
#endif

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   character(len=*), intent(in) :: vtype
   !---------------------------------------------------------------------------------------!

#if USE_HDF5
   !----- Local variables -----------------------------------------------------------------!
   character(len=str_len)   :: anamel
   character(len=2)         :: cgrid
   character(len=40)        :: subaname
   character(len=16)        :: varn
   character(len=1)         :: vnam
   character(len=10)        :: c0
   integer                  :: ngr,nv,nvcnt,lenl
   integer                  :: outyear,outmonth,outdate,outhour
   integer           , save :: ncall_head=0,nvtota=0,nvtotl=0,nvtot
   logical                  :: exans
   real(kind=8)             :: timeold
   !----- Local constants -----------------------------------------------------------------!
   real(kind=8),parameter :: zero = 0.d0
   !-----  HDF specific data types --------------------------------------------------------!
   integer, save                 :: ihdfinit = 0
   integer                       :: hdferr
   integer                       :: dsetrank
   integer(HID_T)                :: dspace_id,file_id,dset_id
   integer(HSIZE_T),dimension(4) :: datadims,maxdatadims
   !---------------------------------------------------------------------------------------!

  
   !----- If this run doesn't output HDF5, there is nothing for me to do here -------------!
   if (ioutput /= 3) return

   if (ncall_head == 0) then

      !----- Find total number of fields to be written ------------------------------------!
      do ngr=1,ngrids
         do nv = 1,num_var(ngr)
            if ( vtab_r(nv,ngr)%ianal == 1) nvtota=nvtota+1
            if ( vtab_r(nv,ngr)%ilite == 1) nvtotl=nvtotl+1
         enddo
      enddo

      nvtot=max(nvtota,nvtotl)
     
      ncall_head=1
   end if
  
   timeold=time

   !----- Construct the HDF file name -----------------------------------------------------!
   if(vtype == 'INST') vnam='A'
   if(vtype == 'LITE') vnam='L'
   if(vtype == 'MEAN') vnam='M'
   if(vtype == 'BOTH') vnam='B'
  
   nvcnt=0

   do ngr=1,ngrids

      !----- Create the file name ---------------------------------------------------------!
      write(cgrid,'(a1,i1)') 'g',ngr
      call makefnam(anamel,afilout,time,iyeara,imontha,idatea,  &
                    itimea*100,vnam,cgrid,'h5')

      lenl = len_trim(anamel)
      !------------------------------------------------------------------------------------!
      !     Check whether the file already exists. In case there is such file already,     !
      ! decide between overwriting it or stopping the run.                                 !
      !------------------------------------------------------------------------------------!
      inquire(file=anamel,exist=exans)
      if(exans .and. iclobber == 0) then
         print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         print*,'!!!   trying to open file name :'
         print*,'!!!       ',anamel
         print*,'!!!   but it already exists. run is ended.'
         print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         call abort_run('Analysis file already existed and you said I can''t overwrite!'   &
                       ,'anlhdf','rhdf5.F90')
      end if

      
      !----- Let's first initialize the hdf environment -----------------------------------!
      call h5open_f(hdferr)
         
      !------------------------------------------------------------------------------------!
      !   Open a new HDF file using IO mode H5F_ACC_TRUNC_F. In this case, if a file with  !
      ! the same name already exists, it will overwrite the data in that file, destroying  !
      ! all data.                                                                          !
      !------------------------------------------------------------------------------------!
      call h5fcreate_f(anamel, H5F_ACC_TRUNC_F, file_id, hdferr)


      !------------------------------------------------------------------------------------!
      !   Now we need to create HDF datasets and then put them in the file, cycle all of   !
      ! our variables.                                                                     !
      !------------------------------------------------------------------------------------!
      do nv = 1,num_var(ngr)

               
         if ((vtype == 'INST' .and. vtab_r(nv,ngr)%ianal == 1) .or. &
             (vtype == 'LITE' .and. vtab_r(nv,ngr)%ilite == 1)) then

            varn= vtab_r(nv,ngr)%name
            
            datadims = 0
            maxdatadims = -1
            
            !----- There are multiple possible dimensions to this -------------------------!
            select case (vtab_r(nv,ngr)%idim_type)
            case (2) !----- 2-D (nxp,nyp) -------------------------------------------------!
               datadims(1) = nnxp(ngr)
               datadims(2) = nnyp(ngr)
               dsetrank    = 2

            case (3) !----- 3-D (nzp,nxp,nyp) ---------------------------------------------!
               datadims(1) = nnzp(ngr)
               datadims(2) = nnxp(ngr)
               datadims(3) = nnyp(ngr)
               dsetrank    = 3
               
            case (4) !----- 4-D (nzg,nxp,nyp,npatch) --------------------------------------!
               datadims(1) = nzg
               datadims(2) = nnxp(ngr)
               datadims(3) = nnyp(ngr)
               datadims(4) = npatch
               dsetrank    = 4
               
            case (5) !----- 4-D (nzs,nxp,nyp,npatch) --------------------------------------!
               datadims(1) = nzs
               datadims(2) = nnxp(ngr)
               datadims(3) = nnyp(ngr)
               datadims(4) = npatch
               dsetrank    = 4
               
            case (6) !----- 3-D (nxp,nyp,npatch) ------------------------------------------!
               datadims(1) = nnxp(ngr)
               datadims(2) = nnyp(ngr)
               datadims(3) = npatch
               dsetrank    = 3
              
            case (7) !----- 3-D (nxp,nyp,nwave) -------------------------------------------!
               datadims(1) = nnxp(ngr)
               datadims(2) = nnyp(ngr)
               datadims(3) = nwave
               dsetrank    = 3
               
            case (8) !----- 4-D (nzp,nxp,nyp,nclouds) -------------------------------------!
               datadims(1) = nnzp(ngr)
               datadims(2) = nnxp(ngr)
               datadims(3) = nnyp(ngr)
               datadims(4) = nclouds
               dsetrank    = 4
               
            case (9) !----- 3-D (nxp,nyp,nclouds) -----------------------------------------!
               datadims(1) = nnxp(ngr)
               datadims(2) = nnyp(ngr)
               datadims(3) = nclouds
               dsetrank    = 3
            end select
            
            !----- Create a dataspace -----------------------------------------------------!
            call h5screate_simple_f(dsetrank,datadims,dspace_id,hdferr)
            if (hdferr /= 0) then
               print*,"COULD NOT CREATE DATASPACE"
               print*,trim(varn),dsetrank,datadims,dspace_id,hdferr
               call abort_run ('Failed creating dataspace','anlhdf','rhdf5.F90')
            end if


            !----- Create a dataset, in the file, with the space just specified -----------!
            call h5dcreate_f(file_id,varn, H5T_NATIVE_REAL, dspace_id,dset_id, hdferr)
            if (hdferr /= 0) then
               print*,"COULD NOT CREATE DATASET"
               print*,trim(varn),dsetrank,datadims,dset_id,hdferr
               call abort_run ('Failed creating dataset','anlhdf','rhdf5.F90')
            end if
         

            !----- Fill that dataset, in the file, with the real data. --------------------!
            
            call h5dwrite_f(dset_id,H5T_NATIVE_REAL,vtab_r(nv,ngr)%var_p,datadims,hdferr)
            if (hdferr /= 0) then
               print*,"COULD NOT WRITE TO THE DATASET"
               print*,trim(varn),dsetrank,datadims,dset_id,hdferr
               call abort_run ('Failed writing dataset','anlhdf','rhdf5.F90')
            endif
            
            !----- Close out the dataspace ------------------------------------------------!
            call h5sclose_f(dspace_id, hdferr)
            if (hdferr /= 0) then
               print*,"COULD NOT CLOSE DATASPACE"
               print*,trim(varn),dspace_id
               call abort_run ('Failed closing dataspace','anlhdf','rhdf5.F90')
            endif
            
            !----- Close the dataset, in the file -----------------------------------------!
            call h5dclose_f(dset_id, hdferr)
            if (hdferr /= 0) then
               print*,"COULD NOT CLOSE THE DATASET"
               print*,trim(varn),dset_id
               call abort_run ('Failed closing dataset','anlhdf','rhdf5.F90')
            end if

         end if
         
      end do
      
      !----- Close our beloved HDF file ---------------------------------------------------!
      call h5fclose_f(file_id, hdferr)
      if (hdferr /= 0) then
         print*,'I COULDN''T CLOSE THE HDF FILE, WHOSE ID IS :',file_id
         call abort_run ('Failed closing HDF file','anlhdf','rhdf5.F90')
      end if
   end do

   call h5close_f(hdferr)
   if (hdferr /= 0) then
      call abort_run ('Failed closing HDF environment','anlhdf','rhdf5.F90')
   endif


   !----- Printing the banner -------------------------------------------------------------!
   select case (trim(vtype))
   case ('LITE')
      subaname='  Analysis lite HDF write'
   case ('MEAN')
      subaname='  Averaged analysis HDF write    '
   case ('BOTH')
      subaname='  Averaged analysis lite HDF write   '
   case default
      subaname='  Analysis HDF write         '
   end select

   write(c0,"(f10.1)") time
   write(*,"(/,a)") " === "//trim(adjustl(subaname))//" at Sim time "//trim(adjustl(c0))//" ==="
   write(*,"(a,/)") " === wrote file "//trim(adjustl(anamel))//" ==="

   !----- Reset the time back to the original ---------------------------------------------!
   if(vtype == 'MEAN'.or.vtype == 'BOTH')  time=timeold
#else
   if (ioutput /= 3) return
   write (unit=*,fmt='(a)') '-------------------------------------------------------------'
   write (unit=*,fmt='(a)') '    This run will stop and this time it''s your fault!'
   write (unit=*,fmt='(a)') '    You have set the analysis to output in HDF5, however you'
   write (unit=*,fmt='(a)') ' MUST compile BRAMS with the HDF5 libraries for that!!!'
   write (unit=*,fmt='(a)') '-------------------------------------------------------------'
   call abort_run('No HDF5 libraries were linked to this executable!','anlhdf','rhdf5')
#endif


   return
end subroutine anlhdf
!==========================================================================================!
!==========================================================================================!
