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

subroutine node_sendst(isflag)

  use mem_grid
  use node_mod

  use mem_scratch
  use mem_basic

  implicit none

  integer :: isflag

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,ipos
  integer :: itype,nm,i1,i2,j1,j2,mtp

  if (isflag >= 2 .and. isflag <= 4) itype=isflag
  if (isflag >= 5 .and. isflag <= 6) itype=1

  !---------------------------------------------------
  !
  !   First, before we send anything, let's post the receives


  do nm=1,nmachs
     if (iget_paths(itype,ngrid,nm).ne.0) then
        call MPI_Irecv(node_buffs(nm)%lbc_recv_buff,               &
             node_buffs(nm)%nrecv*f_ndmd_size,MPI_PACKED,             &
             machs(nm),10000+isflag,MPI_COMM_WORLD,irecv_req(nm),ierr )
     endif
  enddo

  !---------------------------------------------------
  !
  !   Now we can actually go on to sending the stuff

  do nm=1,nmachs

     if(ipaths(1,itype,ngrid,nm).ne.0) then

        i1=ipaths(1,itype,ngrid,nm)
        i2=ipaths(2,itype,ngrid,nm)
        j1=ipaths(3,itype,ngrid,nm)
        j2=ipaths(4,itype,ngrid,nm)

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

        if(isflag == 2) then
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%up  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        elseif (isflag == 3) then
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%vp  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                      node_buffs(nm)%nsend,ipos,MPI_COMM_WORLD,ierr)
        elseif (isflag == 4) then
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%pp  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                      node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        elseif(isflag == 5) then
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%up  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                      node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%vp  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        elseif (isflag == 6) then
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%wp  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
           call mkstbuff(mzp,mxp,myp,basic_g(ngrid)%pp  &
                ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtp)
           call MPI_Pack(scratch%scr1,mtp,MPI_REAL,node_buffs(nm)%lbc_send_buff,&
                node_buffs(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

        endif

        call MPI_Isend(node_buffs(nm)%lbc_send_buff,  &
             ipos-1,  &
             MPI_PACKED,ipaths(5,itype,ngrid,nm),10000+isflag,      &
             MPI_COMM_WORLD,isend_req(nm),ierr)
     endif

  enddo


  return
end subroutine node_sendst
!
!     ****************************************************************
!
subroutine node_getst(isflag)

  use mem_grid
  use node_mod

  use mem_scratch
  use mem_basic

  implicit none

  integer :: isflag

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,ipos
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer :: itype,nm,i1,i2,j1,j2,mtp,mtc,node_src,nptsxy
  !integer :: ibytes,msgid,ihostnum

  if (isflag.ge.2.and.isflag.le.4) itype=isflag
  if (isflag.ge.5.and.isflag.le.6) itype=1

  !_____________________________________________________________________
  !
  !  First, let's make sure our sends are all finished and de-allocated

  do nm=1,nmachs
     if(ipaths(1,itype,ngrid,nm).ne.0) then
        call MPI_Wait(isend_req(nm),status,ierr)
     endif
  enddo
  !_____________________________________________________________________
  !
  !  Now, let's wait on our receives

  do nm=1,nmachs
     if (iget_paths(itype,ngrid,nm).ne.0) then
        call MPI_Wait(irecv_req(nm),status,ierr)
     endif
  enddo
  !_____________________________________________________________________
  !
  !  We got all our stuff. Now unpack it into appropriate space.


  do nm=1,nmachs

     if (iget_paths(itype,ngrid,nm).ne.0) then
        ipos = 1
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
                        i1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
                        i2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
                        j1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
                        j2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
                        node_src,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

        nptsxy=(i2-i1+1)*(j2-j1+1)

        mtp=nzp * nptsxy

        if(isflag == 2) then
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%up  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        elseif(isflag == 3) then
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%vp  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        elseif(isflag == 4) then
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%pp  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        elseif(isflag == 5) then
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%up  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%vp  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        elseif(isflag == 6) then
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%wp  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
           call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
             scratch%scr1,mtp,MPI_REAL,MPI_COMM_WORLD,ierr)
           call exstbuff(mzp,mxp,myp,basic_g(ngrid)%pp  &
               ,scratch%scr1,i1-i0,i2-i0,j1-j0,j2-j0,mtc)
        endif

     endif

  enddo

return
end subroutine node_getst

!*********************************************************************

subroutine mkstbuff(n1,n2,n3,a,b,il,ir,jb,jt,ind)
  implicit none
  integer :: n1,n2,n3,il,ir,jb,jt,ind
  real :: a(n1,n2,n3),b(*)

  integer :: i,j,k

  ind=0
  do j=jb,jt
     do i=il,ir
        do k=1,n1
           ind=ind+1
           b(ind)=a(k,i,j)
        enddo
     enddo
  enddo

  return
end subroutine mkstbuff

!*********************************************************************

subroutine mkstbuffs(n1,n2,n3,a,b,il,ir,jb,jt,ind)
  implicit none
  integer :: n1,n2,n3,il,ir,jb,jt,ind
  real :: a(n1,n2,n3),b(*)

  integer :: i,j,k

  ind=0
  do j=jb,jt
     do i=il,ir
        do k=1,n3
           ind=ind+1
           b(ind)=a(i,j,k)
        enddo
     enddo
  enddo

  return
end subroutine mkstbuffs

!*********************************************************************

subroutine exstbuff(n1,n2,n3,a,b,il,ir,jb,jt,ind)
  implicit none
  integer :: n1,n2,n3,il,ir,jb,jt,ind
  real :: a(n1,n2,n3),b(*)

  integer :: i,j,k

  ind=0
  do j=jb,jt
     do i=il,ir
        do k=1,n1
           ind=ind+1
           a(k,i,j)=b(ind)
        enddo
     enddo
  enddo

  return
end subroutine exstbuff

!*********************************************************************

subroutine exstbuffs(n1,n2,n3,a,b,il,ir,jb,jt,ind)
  implicit none
  integer :: n1,n2,n3,il,ir,jb,jt,ind
  real :: a(n1,n2,n3),b(*)

  integer :: i,j,k

  ind=0
  do j=jb,jt
     do i=il,ir
        do k=1,n3
           ind=ind+1
           a(i,j,k)=b(ind)
        enddo
     enddo
  enddo

  return
end subroutine exstbuffs


