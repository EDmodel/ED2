!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################

subroutine master_getcflcpu()

use mem_grid
use rpara
#if defined(RAMS_MPI)
   use mpi
#endif

implicit none
include 'interface.h'
integer :: ngr,nm,k
real, save, allocatable :: buff1(:),buff2(:)
integer, save :: ncall=0
#if defined(RAMS_MPI)
integer, dimension(MPI_STATUS_SIZE) :: status
#endif
integer :: ierr

if (ncall==0) then
   allocate (buff1(maxgrds))
   allocate (buff2(maxgrds))
   ncall=1
endif

do ngr = 1,ngrids
   cflxy(ngr) = 0.
   cflz(ngr) = 0.
enddo

#if defined(RAMS_MPI)
do nm=1,nmachs
   call MPI_Recv(ptimes(machnum(nm),1),1,MPI_REAL,machnum(nm),41,MPI_COMM_WORLD,status,ierr)
   call MPI_Recv(ptimes(machnum(nm),2),1,MPI_REAL,machnum(nm),42,MPI_COMM_WORLD,status,ierr)
   call MPI_Recv(buff1,maxgrds,MPI_REAL,machnum(nm),43,MPI_COMM_WORLD,status,ierr)
   call MPI_Recv(buff2,maxgrds,MPI_REAL,machnum(nm),44,MPI_COMM_WORLD,status,ierr)
   do ngr = 1,ngrids
      cflxy(ngr) = max(cflxy(ngr),buff1(ngr))
      cflz(ngr) = max(cflz(ngr),buff2(ngr))
   enddo
enddo
#endif
return
end subroutine master_getcflcpu

!     *****************************************************************

subroutine node_putcflcpu(totcpu,totwall)

use mem_grid
use node_mod
#if defined(RAMS_MPI)
   use mpi
#endif

implicit none

#if defined(RAMS_MPI)
include 'interface.h'
#endif
real :: totcpu,totwall
integer :: ierr

#if defined(RAMS_MPI)
call MPI_Send(totcpu,1,MPI_REAL,master_num,41,MPI_COMM_WORLD,ierr)
call MPI_Send(totwall,1,MPI_REAL,master_num,42,MPI_COMM_WORLD,ierr)
call MPI_Send(cflxy,maxgrds,MPI_REAL,master_num,43,MPI_COMM_WORLD,ierr)
call MPI_Send(cflz,maxgrds,MPI_REAL,master_num,44,MPI_COMM_WORLD,ierr)
#endif
return
end subroutine node_putcflcpu

!     *****************************************************************

subroutine master_putdtsched(isendflg,isendlite,isendmean  &
                            ,isendboth,ntsend)
use mem_grid
use rpara
#if defined(RAMS_MPI)
   use mpi
#endif

implicit none

integer :: ierr
integer :: isendflg,isendlite,isendmean,isendboth,ntsend

#if defined(RAMS_MPI)
include 'interface.h'
#endif
integer :: nm

#if defined(RAMS_MPI)
call mpi_bcast(isendflg,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
call MPI_Bcast(isendlite,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
call mpi_bcast(isendmean,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
call mpi_bcast(isendboth,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

call mpi_bcast(ntsend,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

if(ntsend == 1) then
   call mpi_bcast(nnacoust,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call mpi_bcast(nndtrat,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call mpi_bcast(isched,maxsched*maxschent,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call mpi_bcast(nsubs,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call mpi_bcast(dtlongn,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call mpi_bcast(sspct,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
endif
#endif

return
end subroutine master_putdtsched

!     *****************************************************************

subroutine node_getdtsched(isendflg,isendlite,isendmean,isendboth)

use mem_grid
use node_mod
#if defined(RAMS_MPI)
   use mpi
#endif

implicit none
integer :: isendflg,isendlite,isendmean,isendboth,ntsend

#if defined(RAMS_MPI)
include 'interface.h'
#endif
integer :: ierr

#if defined(RAMS_MPI)
call MPI_Bcast(isendflg,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
call MPI_Bcast(isendlite,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
call MPI_Bcast(isendmean,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
call MPI_Bcast(isendboth,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
call MPI_Bcast(ntsend,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

if(ntsend == 1) then
   call MPI_Bcast(nnacoust,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nndtrat,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(isched,maxsched*maxschent,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nsubs,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dtlongn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(sspct,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
endif
#endif

return
end subroutine node_getdtsched




