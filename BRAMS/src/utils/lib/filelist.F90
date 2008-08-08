!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

! Newer version that just uses ls and C to form unique filenames

subroutine RAMS_filelist (fnames,file_prefix,nfile)
implicit none
integer :: nfile
character(len=*) :: fnames(*),file_prefix

character(len=128) :: file,command,cdir
character(len=30) :: tmpname

integer :: iflag,iprelen,nc,nf,iun

integer,parameter :: maxlist=6666,maxlen=200
character(len=maxlist*maxlen) :: filelist
character(len=256)::dir,filetemp
character(len=256)::prefix
integer,dimension(maxlist) :: indices
integer :: n,i,j

      
! this version uses nfile as flag for whether to stop if no files exist
! if nfile.ge.0, then stop

iflag=nfile

nfile = 0

iprelen=len_trim(file_prefix)
if(iprelen == 0) iprelen=len(file_prefix)
      
#if defined (PC_NT1)

   ! First change all "/" to "\" so same namelist can be used 
   !   for Unix/Linux/Windows
   
   do nc=1,iprelen
      if(file_prefix(nc:nc) == '/') file_prefix(nc:nc)='\'
   enddo

   command=  &
     'dir /b '//file_prefix(1:len_trim(file_prefix))//' >c:\temp\RAMS_filelist'
   call system(command)
 
   ! open the directory list
 
   iun=98
   open(unit=iun,file='c:\temp\RAMS_filelist',status='old',err=15)
 
   ! read through the files
   ! windows doesn't put directory names on "dir", so...

   do nc=len_trim(file_prefix),1,-1
      if(file_prefix(nc:nc).eq.'\') then
         lndir=nc
         cdir=file_prefix(1:lndir)
         goto 25
      endif
   enddo
   lndir=2
   cdir='.\'
   25 continue
 
   do nf=1,1000000
      read(iun,'(a128)',end=30,err=30) file
      fnames(nf) = cdir(1:lndir)//file
   enddo
 
   30 continue

   close(iun)

   command= 'del c:\temp\RAMS_filelist'
   call system(command)
      
#else


   ! The following subroutine uses implicit c functions to return
   ! an vector of file names that match the search prefix strings
   ! The strings contained in the string vector are indexed by the
   ! integer array "
   prefix = file_prefix(1:iprelen) // char(0)
   call filelist_c(n,indices,prefix,filelist)

   do nf=1,n

        filetemp=trim(filelist(indices(nf):indices(nf+1)-1))
        fnames(nf)= filetemp
        print*,filetemp,indices(nf),indices(nf+1),len_trim(filetemp),len_trim(fnames(nf))

    enddo

![Deprecated way]! !  ! Let C determine a unique filename
![Deprecated way]!
![Deprecated way]!    tmpname='/tmp/XXXXXX'//char(0)
![Deprecated way]!    call form_tmpname(tmpname)
![Deprecated way]! ! 
![Deprecated way]!    command = '/bin/ls -1 '//file_prefix(1:iprelen)//' > '//tmpname
![Deprecated way]!    call system(command)
![Deprecated way]!    command = 'chmod 777 '//tmpname
![Deprecated way]!    call system(command)
![Deprecated way]!
![Deprecated way]! !  ! open the directory list and read through the files
![Deprecated way]!
![Deprecated way]!    iun=98
![Deprecated way]!    open(unit=iun,file=tmpname,status='unknown',err=15)
![Deprecated way]!    rewind iun
![Deprecated way]!
![Deprecated way]!    do nf=1,1000000
![Deprecated way]!       read(iun,'(a128)',end=30,err=30) file
![Deprecated way]!       fnames(nf) = file
![Deprecated way]!    enddo
![Deprecated way]!      
![Deprecated way]!    30 continue!
![Deprecated way]! !   close(iun)
![Deprecated way]!
![Deprecated way]!    command= '/bin/rm -f '//tmpname
![Deprecated way]!    call system(command)

#endif


nfile=nf-1

if (nfile == 0) then
   print *, 'No RAMS files for prefix:',file_prefix
   if(iflag >= 0) stop 'RAMS_filelist-no_files'
endif

return 
      
15 print *, 'RAMS_filelist: Error opening temporary RAMS_filelist'
stop 'RAMS_filelist-/tmp file error : run again'
return
      
100 continue

return
end

!geodata/sst-W-0000-01-16-120000-g1.vfm
!12345678901234567890123456789012345678
