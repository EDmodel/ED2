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

subroutine node_sendnbc(ifm,icm)

  use mem_grid
  use node_mod

  use var_tables
  use mem_basic
  use grid_dims, only: maxgrds
  implicit none

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,ipos
  integer :: nm,i1,i2,j1,j2,k1,k2,ng,itype,mtp,iptr,nv
  integer :: ifm,icm
  integer :: la,lz,mpiid

  real, allocatable, save :: buffnest(:)
  integer, save :: ncall=0, membuff,membuff_extra,nvar

  itype=5

  !______________________
  !
  !   First, before we send anything, let's post the receives.

  do nm=1,nmachs
     irecv_req(nm)=0
     if (iget_paths(itype,ifm,nm).ne.0) then
        mpiid=400000+maxgrds*(machs(nm)-1)+icm
        call MPI_Irecv(node_buffs_nest(nm)%pack_recv_buff,                  &
             node_buffs_nest(nm)%nrecv*f_ndmd_size,MPI_PACKED,             &
             machs(nm),mpiid,MPI_COMM_WORLD,irecv_req(nm),ierr )
     endif
  enddo

  ! Send coarse grid points necessary for fine grid boundary interpolation
  !   to fine grid nodes. Note that even though coarse grid points are sent,
  !   ipaths is referenced by the fine grid, since all nests only have one
  !   parent, not vice versa.


  ! Compute size of buffer needed and allocate if necessary
  if(ncall == 0) then

     ncall=1
     membuff_extra=nvar*2+100
     membuff=0      

     do ng=1,ngrids
        do nm=1,nmachs
           if(ipaths(1,itype,ng,nm).ne.0) then
              i1=ipaths(1,itype,ng,nm)
              i2=ipaths(2,itype,ng,nm)
              j1=ipaths(3,itype,ng,nm)
              j2=ipaths(4,itype,ng,nm)
              k1=1
              k2=nnzp(ng)

              nvar=4 + num_scalar(ng)
              mtp=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
              membuff=max(membuff,mtp*nvar)
           endif
        enddo
     enddo

     membuff=membuff+membuff_extra
     allocate (buffnest(membuff))
     print*,'sending nesting condition: Allocate buffer for:',mynum &
              ,membuff,nvar
  endif


  do nm=1,nmachs

     isend_req(nm)=0

     if(ipaths(1,itype,ifm,nm).ne.0) then

        i1=ipaths(1,itype,ifm,nm)
        i2=ipaths(2,itype,ifm,nm)
        j1=ipaths(3,itype,ifm,nm)
        j2=ipaths(4,itype,ifm,nm)
        k1=1
        k2=nnzp(icm)

        mtp=(i2-i1+1)*(j2-j1+1)*(k2-k1+1)

  ! Put variables into buffer. All need coarse grid density weighting first.

        la = 1
        lz = mtp
        call mknest_buff(1,basic_g(icm)%uc,buffnest(la:lz)    &
            ,basic_g(icm)%dn0,mmzp(icm),mmxp(icm),mmyp(icm)   &
            ,mi0(icm),mj0(icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        la   = lz + 1
        lz   = lz + mtp
        call mknest_buff(2,basic_g(icm)%vc,buffnest(la:lz)    &
            ,basic_g(icm)%dn0,mmzp(icm),mmxp(icm),mmyp(icm)   &
            ,mi0(icm),mj0(icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        la   = lz + 1
        lz   = lz + mtp
        call mknest_buff(3,basic_g(icm)%wc,buffnest(la:lz)    &
            ,basic_g(icm)%dn0,mmzp(icm),mmxp(icm),mmyp(icm)   &
            ,mi0(icm),mj0(icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        la   = lz + 1
        lz   = lz + mtp
        call mknest_buff(4,basic_g(icm)%pc,buffnest(la:lz)    &
            ,basic_g(icm)%dn0,mmzp(icm),mmxp(icm),mmyp(icm)   &
            ,mi0(icm),mj0(icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)

        do nv=1,num_scalar(ifm)
           la   = lz + 1
           lz   = lz + mtp
           call mknest_buff(5,scalar_tab(nv,icm)%var_p,buffnest(la:lz)  &
               ,basic_g(icm)%dn0,mmzp(icm),mmxp(icm),mmyp(icm)          &
               ,mi0(icm),mj0(icm),i1,i2,j1,j2,k1,k2,mynum,nm,nv)
        end do

        iptr = lz
        ipos = 1
        call MPI_Pack(i1,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(i2,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(j1,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(j2,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(k1,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(k2,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(mynum,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(nvar,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
        call MPI_Pack(iptr,1,MPI_INTEGER,node_buffs_nest(nm)%pack_send_buff, &
             node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

        call MPI_Pack(buffnest,iptr,MPI_REAL,node_buffs_nest(nm)%pack_send_buff, &
                      node_buffs_nest(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

        mpiid=400000+maxgrds*(mchnum-1)+icm
        call MPI_Isend(node_buffs_nest(nm)%pack_send_buff,  &
             ipos - 1,  &
             MPI_PACKED,ipaths(5,itype,ifm,nm),mpiid,MPI_COMM_WORLD,isend_req(nm),ierr)
     endif

  enddo

  return
end subroutine node_sendnbc

!
!     ****************************************************************
!
subroutine mknest_buff(ivarn,ac,acs,den,m1,m2,m3,i0,j0  &
        ,i1,i2,j1,j2,k1,k2,mynum,nm,nv)
  implicit none
  integer :: ivarn,m1,m2,m3,i0,j0,i1,i2,j1,j2,k1,k2,mynum,nm,nv
  real :: ac(m1,m2,m3),acs(0:k2-k1,0:i2-i1,0:j2-j1),den(m1,m2,m3)

  integer :: i,j,k

  !     ivarn = variable types 1- u
  !                            2- v
  !                            3- w
  !                            4- p
  !                            5- scalar

  if(ivarn.eq.5) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)*den(k,i-i0,j-j0)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.1) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)  &
                   *((den(k,i-i0,j-j0)+den(k,i+1-i0,j-j0))*.5)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.2) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)  &
                   *((den(k,i-i0,j-j0)+den(k,i-i0,j+1-j0))*.5)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.3) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2-1
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)  &
                   *((den(k,i-i0,j-j0)+den(k+1,i-i0,j-j0))*.5)
           enddo
        enddo
     enddo
  elseif(ivarn.eq.4) then
     do j=j1,j2
        do i=i1,i2
           do k=k1,k2
              acs(k-k1,i-i1,j-j1)=ac(k,i-i0,j-j0)
           enddo
        enddo
     enddo
  endif

  return
end subroutine mknest_buff
!
!     ****************************************************************
!
subroutine node_getnbc(ifm,icm)

  use mem_grid
  use node_mod

  use var_tables
  use mem_scratch
  use mem_basic
  use mem_nestb

  implicit none

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,ipos
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer :: ifm,icm
  integer, dimension(maxmach) :: i1c,i2c,j1c,j2c,k1c,k2c,iptv,iptc

  integer :: itype,nm,iptr,machf,nvar,nwords,nv,nxc,nyc,nzc,mtp,la,lz
  !integer :: ibytes,msgid,ihostnum
  real, allocatable, save :: buffnest(:)
  integer, save :: ncall=0, nbuff_save=0

  itype=5

  !_____________________________________________________________________
  !
  !  First, let's make sure our sends are all finished and de-allocated
  do nm=1,nmachs
     if(ipaths(1,itype,ifm,nm).ne.0)then
        call MPI_Wait(isend_req(nm),status,ierr)
     endif
  enddo
  !_____________________________________________________________________
  !
  !  Now, let's wait on our receives

  do nm=1,nmachs
     if(iget_paths(itype,ifm,nm).ne.0)then
        call MPI_Wait(irecv_req(nm),status,ierr)
     endif
  enddo
  !_____________________________________________________________________
  !


  ! Compute size of buffer needed and allocate if necessary

  if(nbuff_nest > nbuff_save) then
     nbuff_save = nbuff_nest
     allocate (buffnest(nbuff_nest))
  endif

  !     From the fine grid nodes, get the coarse grid buffers,
  !      interpolate the boundaries, and put them in the "b" array.

  iptr=0
  do nm=1,nmachs

     if(iget_paths(itype,ifm,nm).ne.0) then
        ipos = 1
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,i1c(nm),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,i2c(nm),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,j1c(nm),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,j2c(nm),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,k1c(nm),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,k2c(nm),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,machf,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,nvar,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,nwords,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        la=1+iptr
        lz=iptr+nwords
        call MPI_Unpack(node_buffs_nest(nm)%pack_recv_buff,node_buffs_nest(nm)%nrecv*f_ndmd_size,     &
                        ipos,buffnest(la:lz),nwords,MPI_REAL,MPI_COMM_WORLD,ierr)
        iptc(nm)=1+iptr
        iptv(nm)=0

        iptr=iptr+nwords
     endif
  enddo

  !   We have all the coarse grid info. Start looping through each variable.

  do nv=1,nvar

  !            First, construct coarse grid variable in scr1.

     call azero(maxnzp*maxnxp*maxnyp,scratch%scr1)
     call azero(maxnzp*maxnxp*maxnyp,scratch%scr2)
     do nm=1,nmachs
        if(iget_paths(itype,ifm,nm).ne.0) then
           nzc=k2c(nm)-k1c(nm)+1
           nxc=i2c(nm)-i1c(nm)+1
           nyc=j2c(nm)-j1c(nm)+1
           mtp=nzc*nxc*nyc
           la=iptc(nm)+iptv(nm)
           lz=la+mtp-1
           call unmkbuff(scratch%scr1,buffnest(la:lz)  &
                ,maxnzp,maxnxp,maxnyp,nzc,nxc,nyc  &
                ,i1c(nm),i2c(nm),j1c(nm),j2c(nm),k1c(nm),k2c(nm),mynum)

           iptv(nm)=iptv(nm)+mtp
        endif
     enddo

  !            Do the actual interpolation and put stuff into the "b" array

     if(nv.eq.1) then
        call par_bintp(scratch%scr1,scratch%scr2  &
             ,basic_g(ifm)%dn0  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,1,mi0(ifm),mj0(ifm),mibcon(ifm)  &
             ,nbounds(ifm)%bux,nbounds(ifm)%buy  &
             ,nbounds(ifm)%buz,mynum)
     elseif(nv.eq.2) then
        call par_bintp(scratch%scr1,scratch%scr2  &
             ,basic_g(ifm)%dn0  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,2,mi0(ifm),mj0(ifm),mibcon(ifm)  &
             ,nbounds(ifm)%bvx,nbounds(ifm)%bvy  &
             ,nbounds(ifm)%bvz,mynum)
     elseif(nv.eq.3) then
        call par_bintp(scratch%scr1,scratch%scr2  &
             ,basic_g(ifm)%dn0  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,3,mi0(ifm),mj0(ifm),mibcon(ifm)  &
             ,nbounds(ifm)%bwx,nbounds(ifm)%bwy  &
             ,nbounds(ifm)%bwz,mynum)
     elseif(nv.eq.4) then
        call par_bintp(scratch%scr1,scratch%scr2  &
             ,basic_g(ifm)%dn0  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,4,mi0(ifm),mj0(ifm),mibcon(ifm)  &
             ,nbounds(ifm)%bpx,nbounds(ifm)%bpy  &
             ,nbounds(ifm)%bpz,mynum)
     else
        call par_bintp(scratch%scr1,scratch%scr2  &
             ,basic_g(ifm)%dn0  &
             ,maxnzp,maxnxp,maxnyp,nnzp(ifm)  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,5,mi0(ifm),mj0(ifm),mibcon(ifm)  &
             ,nbounds(ifm)%bsx(:,:,:,nv-4),nbounds(ifm)%bsy(:,:,:,nv-4)  &
             ,nbounds(ifm)%bsz(:,:,:,nv-4),mynum)
     endif


  enddo
  return
end subroutine node_getnbc
!
!     ****************************************************************
!
subroutine unmkbuff(ac,buff,max1,max2,max3,m1,m2,m3  &
     ,i1,i2,j1,j2,k1,k2,mynum)
implicit none
integer :: max1,max2,max3,m1,m2,m3,i1,i2,j1,j2,k1,k2,mynum
real :: ac(max1,max2,max3),buff(0:m1-1,0:m2-1,0:m3-1)

integer :: i,j,k

do j=j1,j2
   do i=i1,i2
      do k=k1,k2
         ac(k,i,j)=buff(k-k1,i-i1,j-j1)
      enddo
   enddo
enddo

return
end subroutine unmkbuff
!
!     ****************************************************************
!
subroutine prtlev(ac,m1,m2,m3,lev,my)
implicit none
integer :: m1,m2,m3,lev,my
real :: ac(m1,m2,m3)
integer :: j,i

print*,'======',my,'====================',m1,m2,m3
do j=m3,10,-1
   print '(62f6.0)',(ac(lev,i,j),i=10,m2)
enddo

return
end subroutine prtlev
