!############################# Change Log ##################################
! 1.0.0.0
!
! 001107 MJB rams_read_header ##
!            Added routine.
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!  Mission Research Corporation / *ASTeR Division
!###########################################################################

subroutine rams_read_header (flnm)

use an_header

implicit none

character*(*) flnm

integer lenf,nv,lastchar

if(allocated(anal_table)) deallocate(anal_table)

! open analysis file and read in commons

open(10,file=flnm(1:lastchar(flnm))//'-head.txt')
read(10,*) nvbtab
allocate (anal_table(nvbtab))
do nv=1,nvbtab
   read(10,*)  anal_table(nv)%string   &
              ,anal_table(nv)%npointer  &
              ,anal_table(nv)%idim_type  &
              ,anal_table(nv)%ngrid  &
              ,anal_table(nv)%nvalues
enddo

call commio('ANAL','READ',10)
close(10)

return
end
