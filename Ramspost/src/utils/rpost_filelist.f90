!======================================= Change Log =======================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine loads all files with the same prefix file_prefix into the fnames   !
! filelist.  This version uses some that just uses dir if this is a Windows run, or some   !
! C code otherwise.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine rams_filelist(fnames,file_prefix,nfile)
   use rpost_dims, only : str_len  & ! intent(in)
                        , maxfiles ! ! intent(in)
   implicit none

   !------ Arguments. ---------------------------------------------------------------------!
   integer                              , intent(out)    :: nfile
   character(len=str_len), dimension(*) , intent(inout)  :: fnames
   character(len=str_len)               , intent(inout)  :: file_prefix
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



   !----- Determine the files based on the prefix -----------------------------------------!
   write (unit=*,fmt='(a,1x,a)') ' [-] filelist_f: Checking prefix: ',trim(file_prefix)
   !---------------------------------------------------------------------------------------!



   !----- Find the length of the prefix. --------------------------------------------------!
   iprelen = len_trim(file_prefix)
   if(iprelen == 0) iprelen=len(file_prefix)
   !---------------------------------------------------------------------------------------!

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
   write (unit=*,fmt='(a)') ' +  Showing first 10 files:'
   do nf=1,n
      fnames(nf) = trim(filelist(indices(nf):indices(nf+1)-1))

      if (nf <= 10) then
         write (unit=*,fmt='(a,1x,i5,1x,a)') '   [-] File #: ',nf,trim(fnames(nf))
      end if
   end do
   !---------------------------------------------------------------------------------------!


   !------ Update the total number of files. ----------------------------------------------!
   nfile=nf-1
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Crash in case there isn't a file with the sought prefix.                         !
   !---------------------------------------------------------------------------------------!
   if (nfile == 0) then
      print *, 'No INPUT files for prefix:',file_prefix
      ! call fatal_error('filelist_f: No files found...','filelist_f','filelist.F90')
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine rams_filelist
!==========================================================================================!
!==========================================================================================!
