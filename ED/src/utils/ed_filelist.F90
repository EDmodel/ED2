!======================================= Change Log =======================================!
! 2.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function is based on the RAMS filelist, which loads all files with the same     !
! prefix file_prefix into the fnames filelist.                                             !
!------------------------------------------------------------------------------------------!
subroutine ed_filelist(fnames,file_prefix,nfile)
 
   use ed_max_dims, only : str_len  & ! intent(in)
                         , maxfiles ! ! intent(in)

   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                        , intent(out)    :: nfile
   character(len=*), dimension(*) , intent(inout)  :: fnames
   character(len=str_len)         , intent(inout)  :: file_prefix
   !------ Local variables. ---------------------------------------------------------------!
   character(len=maxfiles*str_len)                 :: filelist
   character(len=str_len)                          :: prefix
   character(len=str_len)                          :: myform
   integer, dimension(maxfiles)                    :: indices
   integer                                         :: iprelen
   integer                                         :: nf
   integer                                         :: n
   !---------------------------------------------------------------------------------------!



   !----- Initialise the file counter. ----------------------------------------------------!
   nfile = 0
   !---------------------------------------------------------------------------------------!



   !----- String format.  This works if str_len is between 100 and 999. -------------------!
   write(myform,fmt='(a,i3.3,a)') '(a',str_len,')'
   !---------------------------------------------------------------------------------------!





   !----- Find the length of the prefix. --------------------------------------------------!
   iprelen = len_trim(file_prefix)
   if(iprelen == 0) iprelen=len(file_prefix)
   !---------------------------------------------------------------------------------------!

   !----- Determine the files based on the prefix -----------------------------------------!
   write (unit=*,fmt='(a,1x,a )') ' [-] filelist_f: Checking prefix: ',trim(file_prefix)
   write (unit=*,fmt='(a,1x,i5)') '                 Prefix length:   ',iprelen
   !---------------------------------------------------------------------------------------!



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
   open(unit=iun,file='c:\temp\ED22_filelist',status='old',iostat=ierr)
   if (ierr /= 0) call fatal_error('filelist_f: Error opening temporary MEVI_filelist'  &
                                  ,'filelist_f','filelist.F90')
   !---------------------------------------------------------------------------------------!



   !----- Read through the files. Windows doesn't put directory names on "dir", so... -----!
   foundbslash = .false.
   nameloop: do nc =len_trim(file_prefix),1,-1
      if (file_prefix(nc:nc) == '\') then
         lndir=nc
         cdir=file_prefix(1:lndir)
         foundbslash = .true.
         exit nameloop
      end if
   end do nameloop
   if (.not. foundbslash) then
     lndir=2
     cdir='.\'
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Get the file names from the scratch file.                                          !
   !---------------------------------------------------------------------------------------!
   readloop: do
      read(unit=iun,fmt=myform,iostat=ierr) file
      if (ierr /= 0) exit readloop
      fnames(nf) = cdir(1:lndir)//file
   end do readloop
   close(unit=iun,status='delete')
   !---------------------------------------------------------------------------------------!

#else

   !----- Appending char(0), so C can handle the characters properly. ---------------------!
   prefix = file_prefix(1:iprelen) // char(0)

   !---------------------------------------------------------------------------------------!
   !    The following subroutine uses intrinsic C functions to return an vector of file    !
   ! names that match the search prefix strings. The strings contained in the string       !
   ! vector are indexed by the integer array.                                              !
   !---------------------------------------------------------------------------------------!
   call filelist_c(n,indices,prefix,filelist)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Print the first 10 files on screen, so the user can check whether it makes sense !
   ! or not.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (n > 10) then
      write (unit=*,fmt='(a,i5,a)') ' +  Showing first 10 files (out of ',n,'):'
   else
      write (unit=*,fmt='(a,i5,a)') ' +  Showing all ',n,' files:'
   end if
   do nf=1,n
      fnames(nf) = trim(filelist(indices(nf):indices(nf+1)-1))

      if (nf <= 10) then
         write (unit=*,fmt='(a,1x,i5,1x,a)') '   [-] File #: ',nf,trim(fnames(nf))
      end if
   end do
   !---------------------------------------------------------------------------------------!



#endif
   
   !------ Update the total number of files. ----------------------------------------------!
   nfile=nf-1
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Crash in case there isn't a file with the sought prefix.                         !
   !---------------------------------------------------------------------------------------!
   if (nfile == 0) then
      print *, 'No INPUT files for prefix:',file_prefix
      call fatal_error('filelist_f: No files found...','filelist_f','filelist.F90')
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_filelist
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed1_fileinfo(text,nfiles,full_list,ntype,type_list,tlon_list,tlat_list)
   use ed_max_dims, only: str_len,maxfiles,maxlist
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)                           , intent(in)  :: text   ! Type file extension
   integer                                    , intent(in)  :: nfiles ! # of input files
   character(len=str_len), dimension(maxlist) , intent(in)  :: full_list 
   integer                                    , intent(inout) :: ntype
   character(len=str_len), dimension(maxfiles), intent(inout) :: type_list
   real                  , dimension(maxfiles), intent(inout) :: tlon_list, tlat_list
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
   case ('.sss','.pss','.css','.txt')
      okdot = 3
   case ('.lu')
      okdot = 2
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






!==========================================================================================!
!==========================================================================================!
!     This subroutine will go through a list of file and retain only those that are actual !
! ED-2.1 history files.                                                                    !
!------------------------------------------------------------------------------------------!
subroutine ed21_fileinfo(nfiles,full_list,nhisto,histo_list)
   use ed_max_dims, only: str_len,maxfiles,maxlist
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                                    , intent(in)    :: nfiles ! # of input files
   character(len=str_len), dimension(maxlist) , intent(in)    :: full_list 
   integer                                    , intent(inout) :: nhisto
   character(len=str_len), dimension(maxfiles), intent(inout) :: histo_list
   !----- Local variables -----------------------------------------------------------------!
   character(len=str_len)                                     :: this_file
   integer                                                    :: n
   integer                                                    :: posend
   integer                                                    :: posdot
   integer                                                    :: posdsd
   !----- Local constants. ----------------------------------------------------------------!
   integer                                    , parameter     :: okdot = 2
   integer                                    , parameter     :: okdsd = 26
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Counting how many files with the extension are there.  If there is no such file,  !
   ! return without anything else assigned...                                              !
   !---------------------------------------------------------------------------------------!
   nhisto        = 0
   histo_list(:) = ''
   listloop: do n=1,nfiles
      !------------------------------------------------------------------------------------!
      !     Checking whether this file is the type we are looking for.  To qualify as a    !
      ! restart/history ED-2.1 file, the file must end with .h5, and must have a -S- right !
      ! before the date information.                                                       !
      !------------------------------------------------------------------------------------!
      this_file = full_list(n)
      posdot    = index(this_file,'.h5',.true.)
      posdsd    = index(this_file,'-S-',.true.)
      posend    = len_trim(this_file)

      !------------------------------------------------------------------------------------!
      !     If anything is not in place, then this file will not be used for the history   !
      ! initialisation, since it couldn't be considered an ED-2.1 history file.            !
      !------------------------------------------------------------------------------------!
      if (posdot == 0 .or. posend-posdot /= okdot) cycle listloop
      if (posdsd == 0 .or. posend-posdsd /= okdsd) cycle listloop


      !------------------------------------------------------------------------------------!
      !     If we reach this point, it means that the filename at least looks like a       !
      ! restart.  Adding one more to our count and assigning the file to this count.       !
      !------------------------------------------------------------------------------------!
      nhisto             = nhisto + 1
      histo_list(nhisto) = this_file

   end do listloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed21_fileinfo
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine copies the values defined for grid 1 when multiple grids are used   !
! but the user has defined only the first grid.                                            !
!------------------------------------------------------------------------------------------!
subroutine copy_path_from_grid_1(ngrids,varname,pathval)
   use ed_max_dims, only : undef_path & ! intent(in)
                         , str_len    & ! intent(in)
                         , maxgrds    ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                   , intent(in)    :: ngrids
   character(len=*)      , intent(in)                        :: varname
   character(len=str_len), dimension(maxgrds), intent(inout) :: pathval
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                   :: ifm
   !---------------------------------------------------------------------------------------!

   !----- No need to bother if this has a single grid or if values were defined. ----------!
   do ifm=2,ngrids
      if (trim(pathval(ifm)) == trim(undef_path)) then
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a,1x,a)')  ' - Variable:',trim(varname)
         write (unit=*,fmt='(a,1x,i5)') ' - Undefined value for grid:',ifm
         write (unit=*,fmt='(a,1x,i5)') ' - Copying value from grid:',1
         write (unit=*,fmt='(a,1x,a)')  ' - Value assigned:',trim(pathval(1))
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write (unit=*,fmt='(a)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

         pathval(ifm) = pathval(1)
      end if
   end do

   return
end subroutine copy_path_from_grid_1
!==========================================================================================!
!==========================================================================================!
