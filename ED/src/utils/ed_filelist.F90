!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################
!     This function is based on the RAMS filelist, which loads all files with the same     !
! prefix file_prefix into the fnames filelist.                                             !
!------------------------------------------------------------------------------------------!
subroutine ed_filelist(fnames,file_prefix,nfile)
 
   use max_dims, only : str_len, maxfiles

   implicit none
   integer                        , intent(out)    :: nfile
   character(len=*), dimension(*) , intent(out)    :: fnames
   character(len=str_len)         , intent(inout)  :: file_prefix

   character(len=str_len)                          :: file,command,cdir
   integer                                         :: iflag,iprelen,nc,nf,iun,lndir

   character(len=maxfiles*str_len)                 :: filelist
   character(len=str_len)                          :: dir
   character(len=str_len)                          :: prefix
   character(len=str_len)                          :: myform
   integer, dimension(maxfiles)                    :: indices
   integer                                         :: n,i,j,ierr
   logical                                         :: foundbslash

   nfile = 0
   write(myform,fmt='(a,i3.3,a)') '(a',str_len,')'

   !----- Determine the files based on the prefix -----------------------------------------!
   write (unit=*,fmt='(a,1x,a)') ' [-] filelist_f: Checking prefix: ',trim(file_prefix)

   iprelen = len_trim(file_prefix)
   if(iprelen == 0) iprelen=len(file_prefix)
         
!----- Use preprocessor tool. Windows require different way to list... --------------------!
#if defined (PC_NT1)
   !---------------------------------------------------------------------------------------!
   ! First change all "/" to "\" so same namelist can be used for Unix/Linux/Windows       !
   !---------------------------------------------------------------------------------------!
   do nc=1,iprelen
      if(file_prefix(nc:nc) == '/') file_prefix(nc:nc)='\'
   enddo
   command=  &
     'dir /b '//file_prefix(1:len_trim(file_prefix))//' >c:\temp\MEVI_filelist'
   call system(trim(command))
   
   !----- Open the directory list ---------------------------------------------------------!
   iun=66
   open(unit=iun,file='c:\temp\MEVI_filelist',status='old',iostat=ierr)
   if (ierr /= 0) call fatal_error('filelist_f: Error opening temporary MEVI_filelist'  &
                                  ,'filelist_f','filelist.F90')
   
   !----- Read through the files. Windows doesn't put directory names on "dir", so... -----!
   foundbslash=.false.
   nameloop: do nc=len_trim(file_prefix),1,-1
      if(file_prefix(nc:nc).eq.'\') then
         lndir=nc
         cdir=file_prefix(1:lndir)
         foundbslash = .true.
         exit nameloop
      endif
   end do nameloop
   if (.not. foundbslash) then
     lndir=2
     cdir='.\'
   end if
   
   readloop: do
      read(iun,'(a128)',iostat=ierr) file
      if (ierr /= 0) exit readloop
      fnames(nf) = cdir(1:lndir)//file
   end do readloop
   close(iun)

   command= 'del c:\temp\RAMS_filelist'
   call system(command)
      
#else

   !----- Appending char(0), so C can handle the characters properly. ---------------------!
   prefix = file_prefix(1:iprelen) // char(0)

   !---------------------------------------------------------------------------------------!
   !    The following subroutine uses intrinsic C functions to return an vector of file    !
   ! names that match the search prefix strings. The strings contained in the string       !
   ! vector are indexed by the integer array.                                              !
   !---------------------------------------------------------------------------------------!
   call filelist_c(n,indices,prefix,filelist)

   do nf=1,n

      fnames(nf) = trim(filelist(indices(nf):indices(nf+1)-1))
      !write (unit=*,fmt='(a,1x,i5,1x,a)') '   [-] File #: ',nf,trim(fnames(nf))
   end do

#endif
   
   nfile=nf-1

   if (nfile == 0) then
      print *, 'No INPUT files for prefix:',file_prefix
      call fatal_error('filelist_f: No files found...','filelist_f','filelist.F90')
   end if

   return
end subroutine ed_filelist
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed1_fileinfo(text,nfiles,full_list,ntype,type_list,tlon_list,tlat_list)
   use max_dims, only: str_len,maxfiles
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)                           , intent(in)  :: text   ! Type file extension
   integer                                    , intent(in)  :: nfiles ! # of input files
   character(len=str_len), dimension(nfiles)  , intent(in)  :: full_list 
   integer                                    , intent(out) :: ntype
   character(len=str_len), dimension(maxfiles), intent(out) :: type_list
   real                  , dimension(maxfiles), intent(out) :: tlon_list, tlat_list
   !----- Local variables -----------------------------------------------------------------!
   character(len=str_len)                                   :: this_file
   integer                                                  :: n,posend,posdot,posextradot
   integer                                                  :: latbeg,lonbeg,poslondot
   integer                                                  :: okdot
   !---------------------------------------------------------------------------------------!


   !----- Safety check, just to avoid assigning the wrong position... ---------------------!
   select case(text)
   case ('.site')
      okdot = 4
   case ('.pss','.css')
      okdot = 3
   end select


   !---------------------------------------------------------------------------------------!
   !     Counting how many files with the extension are there.  If there is no such file,  !
   ! return without anything else assigned...                                              !
   !---------------------------------------------------------------------------------------!
   ntype = 0
   listloop: do n=1,nfiles
      !---- Checking whether this file is the type we are looking for ---------------------!
      this_file = full_list(n)
      posdot    = index(this_file,text,.true.)
      posend    = len_trim(this_file)
            
      if (posdot == 0 .or. posend-posdot /= okdot) cycle listloop
      
      !----- Seeking for some possible dots after lon but before the file extension -------!
      posextradot = index(this_file(1:(posdot-1)),'.',.true.)

      !---- It is! Adding one more to our count and assigning the file to this count. -----!
      ntype=ntype+1
      type_list(ntype) = this_file
      
      !---- Figuring out the latitude and longitude of this file. -------------------------!
      latbeg = index(this_file,'lat',.true.) + 3
      lonbeg = index(this_file,'lon',.true.) + 3
      poslondot = index(this_file(lonbeg:(posdot-1)),'.',.false.)+lonbeg-1

      read(this_file(latbeg:lonbeg-4),fmt=*) tlat_list(ntype)


      !---- If poslondot is the same as posextradot, no extra dot was found ---------------!
      if (poslondot == posextradot) then
         read(this_file(lonbeg:posdot-1),fmt=*) tlon_list(ntype)
      else
         read(this_file(lonbeg:posextradot-1),fmt=*) tlon_list(ntype)
      end if
      
   end do listloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed1_fileinfo
!==========================================================================================!
!==========================================================================================!
