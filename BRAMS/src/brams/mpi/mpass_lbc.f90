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

subroutine node_sendlbc()

  use mem_grid
  use node_mod

  use var_tables
  use mem_scratch
  use mem_cuparm, only: nclouds
  use grid_dims, only: maxgrds

  implicit none

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,ipos
  integer :: itype,nm,i1,i2,j1,j2,nv,mtp,mpiid

  itype=1
  !______________________
  !
  !   First, before we send anything, let's post the receives. Also, make sure
  !     any pending sends are complete.

  do nm=1,nmachs
     if (iget_paths(itype,ngrid,nm).ne.0) then

        mpiid=300000 +maxgrds*(machs(nm)-1)+ngrid
        call MPI_Irecv(node_buffs(nm)%lbc_recv_buff(1),  &
             node_buffs(nm)%nrecv*f_ndmd_size,  &
             MPI_PACKED,  &
             machs(nm),  &
             mpiid,  &
             MPI_COMM_WORLD,  &
             irecv_req(nm),ierr )
        
     endif
  enddo

  !______________________
  !
  !   Now we can actually go on to sending the stuff

  do nm=1,nmachs

     if(ipaths(1,itype,ngrid,nm).ne.0) then

        i1=ipaths(1,itype,ngrid,nm)
        i2=ipaths(2,itype,ngrid,nm)
        j1=ipaths(3,itype,ngrid,nm)
        j2=ipaths(4,itype,ngrid,nm)

  !      call par_init_put(node_buffs(nm)%lbc_send_buff(1)  &
  !                       ,node_buffs(nm)%nsend )
        ipos = 1
        call MPI_Pack(i1,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff, &
             node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(i2,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff, &
             node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(j1,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff, &
             node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(j2,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff, &
             node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(mynum,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff, &
             node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

        do nv = 1,num_var(ngrid)
           if ( vtab_r(nv,ngrid)%impt1 == 1) then
              select case (vtab_r(nv,ngrid)%idim_type)
              case (2) !---- (X,Y) --------------------------------------------------------!
                 call mkstbuff(1,mxp,myp,1,vtab_r(nv,ngrid)%var_p,scratch%vt3dp            &
                              ,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
                 call MPI_Pack(scratch%vt3dp,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff     &
                              ,node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
              case (3) !---- (Z,X,Y) ------------------------------------------------------!
                 call mkstbuff(mzp,mxp,myp,1,vtab_r(nv,ngrid)%var_p,scratch%vt3dp          &
                              ,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
                 call MPI_Pack(scratch%vt3dp,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff     &
                              ,node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
              case (9) !---- (X,Y,C) ------------------------------------------------------!
                 call mkstbuff(1,mxp,myp,nclouds,vtab_r(nv,ngrid)%var_p,scratch%vt3dp      &
                              ,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
                 call MPI_Pack(scratch%vt3dp,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff     &
                              ,node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
              case default
                 write(unit=*,fmt='(a)')       ' Invalid data type set to mpt1!!!'
                 write(unit=*,fmt='(a,1x,a)')  ' Variable name: ',vtab_r(nv,ngrid)%name
                 write(unit=*,fmt='(a,1x,i5)') ' Variable type: ',vtab_r(nv,ngrid)%idim_type
                 call abort_run ('This data type is not set up to work with LBC'           &
                                ,'node_sendlbc','mpass_lbc.f90')
              end select
           endif
        enddo
        
        mpiid=300000+maxgrds*(mchnum-1)+ngrid
        call MPI_Isend(node_buffs(nm)%lbc_send_buff,ipos-1,  &
             MPI_Packed,ipaths(5,itype,ngrid,nm),mpiid,MPI_COMM_WORLD,  &
             isend_req(nm),ierr)

     endif

  enddo


  return
end subroutine node_sendlbc

!     ****************************************************************

subroutine node_getlbc()

  use mem_grid
  use node_mod

  use var_tables
  use mem_scratch
  use mem_cuparm, only: nclouds

  implicit none

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,ipos
  integer, dimension(MPI_STATUS_SIZE) :: status1
  integer :: itype,nm,ibytes,msgid,ihostnum,i1,i2,j1,j2  &
            ,nmp,nv,node_src,mtc,mtp,nptsxy


  itype=1

  !_____________________________________________________________________
  !
  !  First, let's make sure our sends are all finished and de-allocated


  do nm=1,nmachs
     if(ipaths(1,itype,ngrid,nm).ne.0) then
        call MPI_Wait(isend_req(nm),status1,ierr)
     endif
  enddo

  !_____________________________________________________________________
  !
  !  Now, let's wait on our receives

  do nm=1,nmachs
     if (iget_paths(itype,ngrid,nm).ne.0) then
        call MPI_Wait(irecv_req(nm),status1,ierr)
     endif
  enddo

  !_____________________________________________________________________
  !
  !  We got all our stuff. Now unpack it into appropriate space.

  do nm=1,nmachs

     if (iget_paths(itype,ngrid,nm).ne.0) then



        ipos = 1
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos, &
                        i1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos, &
                        i2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos, &
                        j1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos, &
                        j2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos, &
                        node_src,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        nptsxy=(i2-i1+1)*(j2-j1+1)

        do nv = 1,num_var(ngrid)
           if ( vtab_r(nv,ngrid)%impt1 == 1) then
              select case(vtab_r(nv,ngrid)%idim_type)
              case (2) !---- (X,Y) --------------------------------------------------------!
                 mtp= nptsxy
                 call MPI_Unpack(node_buffs(nm)%lbc_recv_buff                              &
                                ,node_buffs(nm)%nrecv*f_ndmd_size,ipos,scratch%vt3dp,mtp   &
                                ,MPI_REAL,MPI_COMM_WORLD,ierr)
                 call exstbuff(1,mxp,myp,1,vtab_r(nv,ngrid)%var_p,scratch%vt3dp            &
                              ,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
              case (3) !---- (Z,X,Y) ------------------------------------------------------!
                 mtp = nnzp(ngrid) * nptsxy
                 call MPI_Unpack(node_buffs(nm)%lbc_recv_buff                              &
                                ,node_buffs(nm)%nrecv*f_ndmd_size,ipos,scratch%vt3dp,mtp   &
                                ,MPI_REAL,MPI_COMM_WORLD,ierr)
                 call exstbuff(mzp,mxp,myp,1,vtab_r(nv,ngrid)%var_p,scratch%vt3dp          &
                              ,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
              case (9) !---- (X,Y,C) ------------------------------------------------------!
                 mtp= nptsxy * nclouds
                 call MPI_Unpack(node_buffs(nm)%lbc_recv_buff                              &
                                ,node_buffs(nm)%nrecv*f_ndmd_size,ipos,scratch%vt3dp,mtp   &
                                ,MPI_REAL,MPI_COMM_WORLD,ierr)
                 call exstbuff(1,mxp,myp,nclouds,vtab_r(nv,ngrid)%var_p,scratch%vt3dp      &
                              ,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
              case default
                 write(unit=*,fmt='(a)')       ' Invalid data type set to mpt1!!!'
                 write(unit=*,fmt='(a,1x,a)')  ' Variable name: ',vtab_r(nv,ngrid)%name
                 write(unit=*,fmt='(a,1x,i5)') ' Variable type: ',vtab_r(nv,ngrid)%idim_type
                 call abort_run ('This data type is not set up to work with LBC'           &
                                ,'node_getlbc','mpass_lbc.f90')
              end select
           end if
        end do        
     end if
  end do

  return
end subroutine node_getlbc
