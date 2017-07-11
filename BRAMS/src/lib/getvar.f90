!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

integer function RAMS_getvar (string,ngrd,a,b,flnm)
   use grid_dims, only : str_len
   use an_header

   implicit none

   include 'interface.h'

   real :: a(*),b(*)
   integer :: itype,ngrd,rams_c_pos
   character(len=*) :: flnm,string
   character(len=1) :: cgrid
   character(len=str_len) :: flng
   logical :: there
   integer :: ni,npts,iword

   print*,'getvar:',string

   do ni=1,nvbtab

      if(string == anal_table(ni)%string .and. ngrd == anal_table(ni)%ngrid) then
      
         write(cgrid,'(i1)') ngrd
         flng=trim(flnm)//'-g'//cgrid//'.vfm'

         inquire(file=flng,exist=there)
         if(.not.there) then
            write (unit=*,fmt='(a,1x,a)') 'File not found - ',flng
            return
         endif

         npts=anal_table(ni)%nvalues
         itype=anal_table(ni)%idim_type
         iword=anal_table(ni)%npointer

       !  print*,'gv:opening'
         call RAMS_c_open(trim(flng)//char(0),'r'//char(0))
       !  print*,'gv:opened'
         call vfirecr(10,a,npts,'LIN',b,iword)
       !  print*,'gv:vfirecr'
         call RAMS_c_close()
       !  print*,'gv:closed'

         RAMS_getvar=0
         write(unit=*,fmt='(a,1x,a)') 'getvar good:',string
         return

      endif
   enddo

   write(unit=*,fmt='(a,1x,a)') 'Variable not available in this run -',string
   RAMS_getvar=1

   return
end function RAMS_getvar
