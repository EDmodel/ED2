!===============================================================================
! OLAM version 2.11  

! Copyright (C) 2002-2006; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================
Module hdf5_utils

Contains

subroutine shdf5_open(locfn,access,idelete)
        
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

inquire(file=trim(locfn),exist=exists)

! Create a new file or open an existing RAMS file.
if (access(1:1) == 'R') then
   if (.not.exists) then
      print*,'shdf5_open:'
      print*,'   Attempt to open a file for reading that does not exist.'
      print*,'   Filename: ',trim(locfn)
      stop 'shdf5_open: no file'
   else
      if (caccess == 'R ') iaccess = 1
      if (caccess == 'RW') iaccess = 2
      call fh5f_open(trim(locfn)//char(0), iaccess, hdferr)
      
      if (hdferr < 0) then
         print*,'shdf5_open:'
         print*,'   Error opening hdf5 file - error -',hdferr
         print*,'   Filename: ',trim(locfn)
         stop 'shdf5_open: open error'      
      endif
   endif
elseif (access(1:1) == 'W') then
   if (.not.exists) then
      iaccess=2
      call fh5f_create(trim(locfn)//char(0), iaccess, hdferr)
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
         iaccess=1
         call fh5f_create(trim(locfn)//char(0), iaccess, hdferr)
      endif
   endif
   if(hdferr < 0) then
      print*,'HDF5 file create failed:',hdferr
      print*,'file name:',trim(locfn),' ',trim(access), idelete
      stop 'shdf5_open: bad create'
   endif
endif

return
end  subroutine

!===============================================================================

subroutine shdf5_info(dsetname,ndims,dims)

implicit none

character(len=*) :: dsetname ! Dataset name
integer :: dims(*)

integer :: ndims ! Dataset rank (in file)

integer :: hdferr ! Error flag


! Open the dataset.

call fh5d_open( trim(dsetname)//char(0), hdferr)
if (hdferr < 0) then
   print*,'In h5info: call fh5d_open:',trim(dsetname),'hdf5 error =',hdferr
   stop   'shdf5_info'
endif

! Get dataset's dimensions

call fh5s_get_ndims(ndims)

!print*,'ndims: ',ndims

call fh5s_get_dims(dims)
!print*,'dims: ',dims(1:ndims)

return
end subroutine


!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars)
    
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
     
! Prepare memory and options for the write
call fh5_prepare_write(ndims, dimsh, hdferr)
if (hdferr /= 0) then
   print*,'shdf5_orec: can''t prepare requested field:',trim(dsetname)
   return
endif

if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4  ! If native precision is 8 bytes, do h5_type=2
if (ctype(1:1) == 'l') h5_type=5

! Write the dataset.
if (ctype == 'is') then
   call fh5_write(h5_type, ivars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'rs') then
   call fh5_write(h5_type, rvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'cs') then
   call fh5_write(h5_type, cvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ds') then
   call fh5_write(h5_type, dvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ls') then
   call fh5_write(h5_type, lvars, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ia') then
   call fh5_write(h5_type, ivara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ra') then
   call fh5_write(h5_type, rvara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'ca') then
   call fh5_write(h5_type, cvara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'da') then
   call fh5_write(h5_type, dvara, trim(dsetname)//char(0), hdferr)
elseif (ctype == 'la') then
   call fh5_write(h5_type, lvara, trim(dsetname)//char(0), hdferr)
endif
if (hdferr /= 0) then
   print*,'In shdf5_orec: hdf5 write error =',hdferr
   stop 'shdf5_orec: hdf5 write error'
endif

! Close the dataset, the dataspace for the dataset, and the dataspace properties.

call fh5_close_write(hdferr)

return
end subroutine

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars)
        
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
!print*,'-----:',trim(dsetname)
    
! Prepare file and memory space for the read
call fh5_prepare_read(trim(dsetname)//char(0), ndims, dimsh, hdferr)
if (hdferr < 0) then
   print*,'shdf5_irec: can''t prepare requested field:',trim(dsetname)
   return
endif

! Read data from hyperslab in the file into the hyperslab in memory.
if (ctype(1:1) == 'i') h5_type=1
if (ctype(1:1) == 'r') h5_type=2
if (ctype(1:1) == 'c') h5_type=3
if (ctype(1:1) == 'd') h5_type=4  ! If native precision is 8 bytes, do h5_type=2
if (ctype(1:1) == 'l') h5_type=5

if (ctype == 'is') then
   call fh5d_read(h5_type,ivars,hdferr)
elseif (ctype == 'rs') then
   call fh5d_read(h5_type,rvars,hdferr)
elseif (ctype == 'cs') then
   call fh5d_read(h5_type,cvars,hdferr)
elseif (ctype == 'ds') then
   call fh5d_read(h5_type,dvars,hdferr)
elseif (ctype == 'ls') then
   call fh5d_read(h5_type,lvars,hdferr)
elseif (ctype == 'ia') then
   call fh5d_read(h5_type,ivara,hdferr)
elseif (ctype == 'ra') then
   call fh5d_read(h5_type,rvara,hdferr)
elseif (ctype == 'ca') then
   call fh5d_read(h5_type,cvara,hdferr)
elseif (ctype == 'da') then
   call fh5d_read(h5_type,dvara,hdferr)
elseif (ctype == 'la') then
   call fh5d_read(h5_type,lvara,hdferr)
endif

if (hdferr /= 0) then
   print*,'shdf5_irec: call fh5d_read: hdf5 error =',hdferr
   stop
endif

! Close the dataset, the dataspace for the dataset, and the memory space.

call fh5_close_read(hdferr)

return
end subroutine

!===============================================================================

subroutine shdf5_close()
        
implicit none

integer :: hdferr  ! Error flags

! Close RAMS hdf file.

call fh5f_close(hdferr)

return
end  subroutine

end module
