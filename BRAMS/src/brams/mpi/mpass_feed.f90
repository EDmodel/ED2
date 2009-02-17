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

subroutine node_sendfeed(ngr)

  use mem_grid, only: &
       nxtnest,       & ! intent(in)
       kpm,           & ! intent(in)
       nnzp,          & ! intent(in)
       nstratx,       & ! intent(in)
       nstraty          ! intent(in)

  use node_mod, only: &
       nmachs,        & ! intent(in)
       irecv_req,     & ! intent(out)
       iget_paths,    & ! intent(in)
       node_buffs,    & ! intent(inout)
       machs,         & ! intent(in)
       nbuff_feed,    & ! intent(in)
       mynum,         & ! intent(in)
       isend_req,     & ! intent(inout)
       ipaths,        & ! intent(in)
       mmzp,          & ! intent(in)
       mmxp,          & ! intent(in)
       mmyp,          & ! intent(in)
       i0,            & ! intent(in)
       j0,            & ! intent(in)
       mibcon,        & ! intent(in)
       f_ndmd_size

  use mem_basic, only: &
       basic_g          ! intent(in)

  use var_tables, only: &
       num_scalar,      & ! intent(in)
       scalar_tab         ! intent(in)

implicit none

integer :: ngr

include 'interface.h'
include 'mpif.h'

integer :: ierr,ipos
integer :: i1s,i2s,j1s,j2s,k1s,k2s,mtp,i1f,i2f,j1f,j2f,k1f,k2f,i1,i2
integer :: nm,icm,ifm,itype,itypef,nv,iptr,nvar,ibytes,msgid,ihostnum

real, save, allocatable::pbuff(:)
integer, save :: nbuff_save=0

icm=nxtnest(ngr)
ifm=ngr

itype = 6
!______________________
!
!   First, before we send anything, let's post the receives.

do nm=1,nmachs
   irecv_req(nm)=0
   if (iget_paths(itype,ifm,nm).ne.0) then
      call MPI_Irecv(node_buffs(nm)%lbc_recv_buff,  &
           node_buffs(nm)%nrecv*f_ndmd_size,MPI_PACKED,machs(nm),  &
           5500+icm,MPI_COMM_WORLD,irecv_req(nm),ierr )
   endif
enddo
!______________________

!     Allocate new temporary buffer if bigger than the old one.

if(nbuff_feed > nbuff_save) then
   print*,'Allocating feed send buffer:',mynum,nbuff_feed,nbuff_save
  allocate (pbuff(nbuff_feed))
  nbuff_save=nbuff_feed
endif


!     Feed back this fine grid's portion of the each coarse grid node

k1s=kpm(2,ifm)
k2s=kpm(nnzp(ifm)-1,ifm)

do nm=1,nmachs
   isend_req(nm)=0
   if(ipaths(1,itype,ifm,nm).ne.0) then

!            mtp=total number of coarse grid points created
!              and sent from this fine grid node

!      i1s=ipaths(1,itype,ifm,nm)
!      i2s=ipaths(2,itype,ifm,nm)
!      j1s=ipaths(3,itype,ifm,nm)
!      j2s=ipaths(4,itype,ifm,nm)
!      mtp=(i2s-i1s+1)*(j2s-j1s+1)*(k2s-k1s+1)

      itypef=7
      i1f=ipaths(1,itypef,ifm,nm)
      i2f=ipaths(2,itypef,ifm,nm)
      j1f=ipaths(3,itypef,ifm,nm)
      j2f=ipaths(4,itypef,ifm,nm)

      iptr=0
      call fdbackp(1,basic_g(ifm)%uc,pbuff(1+iptr:),mtp  &
          ,basic_g(ifm)%dn0,basic_g(ifm)%dn0u  &
          ,basic_g(ifm)%dn0v  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,icm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,mibcon(ifm) ,nstratx(ifm),nstraty(ifm),mynum,i1s,i2s)
      iptr=iptr+mtp
      call fdbackp(2,basic_g(ifm)%vc,pbuff(1+iptr:),mtp  &
          ,basic_g(ifm)%dn0,basic_g(ifm)%dn0u  &
          ,basic_g(ifm)%dn0v  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,icm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,mibcon(ifm) ,nstratx(ifm),nstraty(ifm),mynum,j1s,j2s)
      iptr=iptr+mtp
      call fdbackp(3,basic_g(ifm)%wc,pbuff(1+iptr:),mtp  &
          ,basic_g(ifm)%dn0,basic_g(ifm)%dn0u  &
          ,basic_g(ifm)%dn0v  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,icm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,mibcon(ifm) ,nstratx(ifm),nstraty(ifm),mynum,i1,i2)
      iptr=iptr+mtp
      call fdbackp(4,basic_g(ifm)%pc,pbuff(1+iptr:),mtp  &
          ,basic_g(ifm)%dn0,basic_g(ifm)%dn0u  &
          ,basic_g(ifm)%dn0v  &
          ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
          ,ifm,icm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
          ,i0,j0,mibcon(ifm) ,nstratx(ifm),nstraty(ifm),mynum,i1,i2)
      iptr=iptr+mtp

      do nv=1,num_scalar(ifm)
         call fdbackp(5,scalar_tab(nv,ifm)%var_p,pbuff(1+iptr:),mtp  &
             ,basic_g(ifm)%dn0,basic_g(ifm)%dn0u  &
             ,basic_g(ifm)%dn0v  &
             ,mmzp(ifm),mmxp(ifm),mmyp(ifm)  &
             ,ifm,icm,i1f-i0,i2f-i0,j1f-j0,j2f-j0  &
             ,i0,j0,mibcon(ifm) ,nstratx(ifm),nstraty(ifm),mynum,i1,i2)
         iptr=iptr+mtp
      enddo

!     We will send master coarse grid indices to nodes.
      ipos = 1
      call MPI_Pack(i1f,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(i2f,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(j1f,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(j2f,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(i1s,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(i2s,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(j1s,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(j2s,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(k1s,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(k2s,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(mynum,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(nvar,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Pack(iptr,1,MPI_INTEGER,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)

      call MPI_Pack(pbuff(1),iptr,MPI_REAL,node_buffs(nm)%lbc_send_buff,node_buffs(nm)%nsend*f_ndmd_size,ipos, &
           MPI_COMM_WORLD,ierr)
      call MPI_Isend(node_buffs(nm)%lbc_send_buff,  &
           ipos-1,  &
           MPI_PACKED,ipaths(5,itype,ifm,nm),5500+icm,MPI_COMM_WORLD,isend_req(nm),ierr)

   endif
enddo

return
end subroutine node_sendfeed
!
!     ****************************************************************
!
subroutine node_getfeed(icm,ifm)

  use mem_grid, only : &
       time,           & ! intent(in) ! debug
       ipm,            & ! intent(in)
       jpm,            & ! intent(in)
       kpm,            & ! intent(in)
       nnxp,           & ! intent(in)
       nnyp,           & ! intent(in)
       nnzp,           & ! intent(in)
       nstratx,        & ! intent(in)
       nstraty,        & ! intent(in)
       nrz,            & ! intent(in)
       nnsttop,        & ! intent(in)
       nnstbot           ! intent(in)
  use node_mod, only : &
       mynum,          & ! intent(in)
       nodebounds,     & ! intent(in)
       mmxp,           & ! intent(in)
       mmyp,           & ! intent(in)
       mmzp,           & ! intent(in)
       nmachs,         & ! intent(in) 
       ipaths,         & ! intent(in) 
       isend_req,      & ! intent(inout)
       irecv_req,      & ! intent(inout)
       iget_paths,     & ! intent(in)
       node_buffs,     & ! intent(inout)
       mi0,            & ! intent(in)
       mj0,            & ! intent(in)
       mibcon,         & ! intent(in)
       f_ndmd_size

  use var_tables, only : &
       num_scalar,     & ! intent(in)
       scalar_tab        ! intent(inout)
  use mem_basic, only : &
       basic_g           ! intent(inout)
  use mem_scratch1, only : &
       scratch1          ! intent(inout)

  ! Included by Alvaro L.Fazenda
  use grid_dims, only: &
       nxpmax,         &  ! intent(in)
       nypmax,         &  ! intent(in)
       nzpmax,         &  ! intent(in)
       maxgrds            ! intent(in)

implicit none

include 'interface.h'
include 'mpif.h'

integer :: ierr,ipos
integer, dimension(MPI_STATUS_SIZE) :: status
integer :: icm,ifm

integer :: i1s,i2s,j1s,j2s,k1s,k2s,mtp,i1z,i2z,j1z,j2z,k1z,k2z
integer :: i2zu,j2zv,k2zw,i2u,j2v,k2w,nfx,nfy,nfz,i1f,i2f,i3f,i4f
integer :: j1f,j2f,k1f,k2f
integer :: nm,itype,nv,iptr,nvar,machf,nwds
logical :: icall
!integer :: ibytes,msgid,ihostnum

real, save, allocatable::pbuff(:)
integer, save :: nbuff_save=0
integer :: ialerr

itype=6
!_____________________________________________________________________
!
!  First, let's make sure our sends are all finished and de-allocated

do nm=1,nmachs
   if (ipaths(1,itype,ifm,nm).ne.0 ) then
      call MPI_Wait(isend_req(nm),status,ierr)
   endif
enddo
!      print*,mynum,'done FEED send wait'
!_____________________________________________________________________
!
!  Now, let's wait on our receives

do nm=1,nmachs
   if (iget_paths(itype,ifm,nm).ne.0) then
      call MPI_Wait(irecv_req(nm),status,ierr)
   endif
enddo
!      print*,mynum,'done FEED recv wait'
!_____________________________________________________________________

!
!     Can we use existing memory for the buffers? If not, allocate new
!       temporary buffer.

!     Allocate new temporary buffer if bigger than the old one.

icall=.true.
do nm=1,nmachs
   if(iget_paths(itype,ifm,nm).ne.0) then

!      call par_assoc_buff(node_buffs(nm)%lbc_recv_buff(1)  &
!                         ,node_buffs(nm)%nrecv)

      ipos = 1
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,i1f,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,i2f,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,j1f,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,j2f,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,i1s,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,i2s,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,j1s,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,j2s,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,k1s,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,k2s,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,machf,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,nvar,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size  &
                      ,ipos,nwds,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      ! Make sure buffer for floating point info big enough
      
      if (nwds >  nbuff_save .or. .not.(allocated(pbuff))) then
         if (allocated(pbuff)) deallocate (pbuff)
         allocate (pbuff(nwds),stat=ialerr)
         if (ialerr /= 0) then
           write(unit=*,fmt='(a)') '---- FATAL Error at node_getfeed (mpass_feed.f90)'
           write(unit=*,fmt='(a)') '---- Unable to allocate pbuff'
           write(unit=*,fmt=('( 4(a,1x,i,1x) )') ) '---- mynum=',mynum,'machf=',machf,'nwds=',nwds,'nbuff_save=',nbuff_save
           stop '---- node_getfeed (mpass_feed.f90)'
         end if
         nbuff_save=nwds
      endif      
      
      call MPI_Unpack(node_buffs(nm)%lbc_recv_buff,node_buffs(nm)%nrecv*f_ndmd_size,ipos,  &
                    pbuff,nwds,MPI_REAL,MPI_COMM_WORLD,ierr)
      if(icall) then
         icall=.false.

! Set the portion of this coarse grid subdomain that will be filled
! with fine grid info to zero.

         i1z = max(ipm(2,ifm),1+mi0(icm))
         j1z = max(jpm(2,ifm),1+mj0(icm))
         k1z = kpm(2,ifm)

         i2z = min(ipm(nnxp(ifm)-1,ifm),mmxp(icm)+mi0(icm))
         j2z = min(jpm(nnyp(ifm)-1,ifm),mmyp(icm)+mj0(icm))
         k2z = kpm(nnzp(ifm)-1,ifm)

         i2zu = min(ipm(nnxp(ifm)-1-nstratx(ifm),ifm)  &
                   ,mmxp(icm)+mi0(icm))
         call zeroout(basic_g(icm)%uc,mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2zu,j1z,j2z,k1z,k2z)
         j2zv = min(jpm(nnyp(ifm)-1-nstraty(ifm),ifm)  &
                   ,mmyp(icm)+mj0(icm))
         call zeroout(basic_g(icm)%vc,mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2zv,k1z,k2z)
         k2zw = kpm(nnzp(ifm)-1-nrz(kpm(nnzp(ifm)-1,ifm)  &
                   ,ifm),ifm)
         call zeroout(basic_g(icm)%wc,mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2z,k1z,k2zw)

         call zeroout(basic_g(icm)%pc,mmzp(icm),mmxp(icm),mmyp(icm) &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2z,k1z,k2z)

         do nv=1,num_scalar(ifm)
            call zeroout(scalar_tab(nv,icm)%var_p  &
                 ,mmzp(icm),mmxp(icm),mmyp(icm)  &
                 ,mi0(icm),mj0(icm),i1z,i2z,j1z,j2z,k1z,k2z)
         enddo

      endif

      call unfpack1(pbuff,scratch1%vtu,scratch1%vtv &
           ,scratch1%vtw,scratch1%vtp &
           ,scratch1%vtscalar,num_scalar(ifm) &
           ,nodebounds(ifm,1),nodebounds(ifm,2),nodebounds(ifm,3) &
           ,nodebounds(ifm,4),2,mmzp(ifm)-1,nodebounds(icm,5) &
           ,nodebounds(icm,6),nodebounds(icm,7),nodebounds(icm,8),k1s,k2s  &
           ,i1f,i2f,j1f,j2f,2,mmzp(ifm)-1,i1s,i2s,j1s,j2s,k1s,k2s)      

   endif

enddo
! Deallocating pbuff
if (allocated(pbuff)) deallocate (pbuff)
if (.not. icall) then
    i2u = i2zu
    j2v = j2zv
    k2w = k2zw
    i1f = nodebounds(ifm,1)
    i2f = nodebounds(ifm,2)
    j1f = nodebounds(ifm,3)
    j2f = nodebounds(ifm,4)
    i1s = nodebounds(icm,5)
    i2s = nodebounds(icm,6)
    j1s = nodebounds(icm,7)
    j2s = nodebounds(icm,8)
    k1f = 2
    k2f = mmzp(ifm)-1

    ! New Sub-Routine - bugs fixed
    call unfdbackp_fixed(1,basic_g(icm)%uc,scratch1%vtu &
         ,i1s,i2s,j1f,j2f,k1f,k2f &
         ,basic_g(icm)%dn0,basic_g(icm)%dn0u &
         ,basic_g(icm)%dn0v  &
         ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
         ,i2u,j2v,k2w,i1s,i2s,j1s,j2s  &
         ,k1s,k2s  &
         ,ipm(1,ifm), nxpmax  &
         ,jpm(1,ifm), nypmax  &
         ,kpm(1,ifm), nzpmax  &
         ,mynum)

    ! New Sub-Routine by Alvaro L.Fazenda
    call unfdbackp_fixed(2,basic_g(icm)%vc,scratch1%vtv &
         ,i1f,i2f,j1s,j2s,k1f,k2f &
         ,basic_g(icm)%dn0,basic_g(icm)%dn0u &
         ,basic_g(icm)%dn0v  &
         ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
         ,i2u,j2v,k2w,i1s,i2s,j1s,j2s  &
         ,k1s,k2s  &
         ,ipm(1,ifm), nxpmax  &
         ,jpm(1,ifm), nypmax  &
         ,kpm(1,ifm), nzpmax  &
         ,mynum)

    ! New Sub-Routine by Alvaro L.Fazenda
    call unfdbackp_fixed(3,basic_g(icm)%wc,scratch1%vtw &
         ,i1f,i2f,j1f,j2f,k1s,k2s &
         ,basic_g(icm)%dn0,basic_g(icm)%dn0u &
         ,basic_g(icm)%dn0v  &
         ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
         ,i2u,j2v,k2w,i1s,i2s,j1s,j2s  &
         ,k1s,k2s  &
         ,ipm(1,ifm), nxpmax  &
         ,jpm(1,ifm), nypmax  &
         ,kpm(1,ifm), nzpmax  &
         ,mynum)

    ! New Sub-Routine by Alvaro L.Fazenda
    call unfdbackp_fixed(4,basic_g(icm)%pc,scratch1%vtp &
         ,i1f,i2f,j1f,j2f,k1f,k2f &
         ,basic_g(icm)%dn0,basic_g(icm)%dn0u &
         ,basic_g(icm)%dn0v  &
         ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
         ,i2u,j2v,k2w,i1s,i2s,j1s,j2s  &
         ,k1s,k2s  &
         ,ipm(1,ifm), nxpmax  &
         ,jpm(1,ifm), nypmax  &
         ,kpm(1,ifm), nzpmax  &
         ,mynum)

    iptr=0
    mtp = (i2f-i1f+1)*(j2f-j1f+1)*(k2f-k1f+1)

    do nv=1,num_scalar(ifm)
       ! New Sub-Routine by Alvaro L.Fazenda
       call unfdbackp_fixed(5,scalar_tab(nv,icm)%var_p  &
            ,scratch1%vtscalar(1+iptr:) &
            ,i1f,i2f,j1f,j2f,k1f,k2f &
            ,basic_g(icm)%dn0,basic_g(icm)%dn0u &
            ,basic_g(icm)%dn0v  &
            ,mmzp(icm),mmxp(icm),mmyp(icm),mi0(icm),mj0(icm)  &
            ,i2u,j2v,k2w,i1s,i2s,j1s,j2s  &
            ,k1s,k2s  &
            ,ipm(1,ifm), nxpmax  &
            ,jpm(1,ifm), nypmax  &
            ,kpm(1,ifm), nzpmax  &
            ,mynum)

       iptr=iptr+mtp
    enddo
end if

if (nnstbot(icm) == 1) then
   call botset(mmzp(icm),mmxp(icm),mmyp(icm)  &
       ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
       ,basic_g(icm)%uc,'U')
   call botset(mmzp(icm),mmxp(icm),mmyp(icm)  &
       ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
       ,basic_g(icm)%vc,'V')
   call botset(mmzp(icm),mmxp(icm),mmyp(icm)  &
       ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
       ,basic_g(icm)%pc,'P')
endif
if (nnsttop(icm) == 1) then
   call topset(mmzp(icm),mmxp(icm),mmyp(icm)  &
     ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
     ,basic_g(icm)%uc,basic_g(icm)%uc,'U')
   call topset(mmzp(icm),mmxp(icm),mmyp(icm)  &
     ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
     ,basic_g(icm)%vc,basic_g(icm)%vc,'V')
   call topset(mmzp(icm),mmxp(icm),mmyp(icm)  &
     ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
     ,basic_g(icm)%pc,basic_g(icm)%pc,'P')
endif

do nv=1,num_scalar(ifm)
   if (nnstbot(icm) == 1) then   
      call botset(mmzp(icm),mmxp(icm),mmyp(icm)  &
          ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
          ,scalar_tab(nv,icm)%var_p,'T')
   endif

   if (nnsttop(icm) == 1) then
      call topset(mmzp(icm),mmxp(icm),mmyp(icm)  &
          ,1,mmxp(icm),1,mmyp(icm),mibcon(icm)  &
          ,scalar_tab(nv,icm)%var_p  &
          ,scalar_tab(nv,icm)%var_p,'T')
   endif
enddo

return
end subroutine node_getfeed

!     ****************************************************************

subroutine zeroout(ac,m1,m2,m3,i0,j0,i1z,i2z,j1z,j2z,k1z,k2z)
implicit none
integer :: m1,m2,m3,i0,j0,i1z,i2z,j1z,j2z,k1z,k2z
real :: ac(m1,m2,m3)

integer :: i,j,k

do j=j1z,j2z
   do i=i1z,i2z
      do k=k1z,k2z
         ac(k,i-i0,j-j0)=0.
      enddo
   enddo
enddo

return
end subroutine zeroout

!     ****************************************************************

subroutine unfpack(ac,acf,i1,i2,j1,j2,k1,k2,i1f,i2f,j1f,j2f,k1f,k2f)
  implicit none

  ! Arguments:
  integer, intent(in) :: i1,i2,j1,j2,k1,k2,i1f,i2f,j1f,j2f,k1f,k2f
  real, intent(inout) :: ac(k1:k2,i1:i2,j1:j2)
  real, intent(in)    :: acf(i1f:i2f,j1f:j2f,k1f:k2f)

  ! Local Variables:
  integer :: i,j,k

  do j=j1f,j2f
     do i=i1f,i2f
        do k=k1f,k2f
           ac(k,i,j)=acf(i,j,k) 
        enddo
     enddo
  enddo

  return
end subroutine unfpack

!     ****************************************************************

subroutine unfpack1(buf,acu,acv,acw,acp,acscalar,nscalar &
     ,i1f,i2f,j1f,j2f,k1f,k2f,i1c,i2c,j1c,j2c,k1c,k2c  &    
     ,i1fp,i2fp,j1fp,j2fp,k1fp,k2fp,i1cp,i2cp,j1cp,j2cp,k1cp,k2cp)   
  implicit none

  ! Arguments:
  integer, intent(in) :: nscalar, i1f, i2f, j1f, j2f, k1f, k2f, &
       i1c, i2c, j1c, j2c, k1c, k2c, i1fp, i2fp, j1fp, j2fp, k1fp, k2fp, &
       i1cp, i2cp, j1cp, j2cp, k1cp, k2cp
  real, intent(in)    :: buf(*)
  real, intent(inout) :: acu(i1c:i2c,j1f:j2f,k1f:k2f)
  real, intent(inout) :: acv(i1f:i2f,j1c:j2c,k1f:k2f)
  real, intent(inout) :: acw(i1f:i2f,j1f:j2f,k1c:k2c)
  real, intent(inout) :: acp(i1f:i2f,j1f:j2f,k1f:k2f)
  real, intent(inout) :: acscalar(i1f:i2f,j1f:j2f,k1f:k2f,nscalar)

  ! Local Variables:
  integer :: iptr
  integer :: i

  !     ivarn = variable types 1- u
  !                            2- v
  !                            3- w
  !                            4- p
  !                            5- scalar

  iptr = 1
  call unfpack(acu,buf(iptr),i1c,i2c,j1f,j2f,k1f,k2f,i1cp,i2cp,j1fp,j2fp,k1fp,k2fp)
  iptr = iptr + (i2cp-i1cp+1)*(j2fp-j1fp+1)*(k2fp-k1fp+1)
  call unfpack(acv,buf(iptr),i1f,i2f,j1c,j2c,k1f,k2f,i1fp,i2fp,j1cp,j2cp,k1fp,k2fp)
  iptr = iptr + (i2fp-i1fp+1)*(j2cp-j1cp+1)*(k2fp-k1fp+1)
  call unfpack(acw,buf(iptr),i1f,i2f,j1f,j2f,k1c,k2c,i1fp,i2fp,j1fp,j2fp,k1cp,k2cp-1)
  iptr = iptr + (i2fp-i1fp+1)*(j2fp-j1fp+1)*(k2cp-k1cp)
  call unfpack(acp,buf(iptr),i1f,i2f,j1f,j2f,k1f,k2f,i1fp,i2fp,j1fp,j2fp,k1fp,k2fp)
  iptr = iptr + (i2fp-i1fp+1)*(j2fp-j1fp+1)*(k2fp-k1fp+1)

  do i=1,nscalar
     call unfpack(acscalar(i1f,j1f,k1f,i),buf(iptr), & 
          i1f,i2f,j1f,j2f,k1f,k2f,i1fp,i2fp,j1fp,j2fp,k1fp,k2fp)
     iptr = iptr + (i2fp-i1fp+1)*(j2fp-j1fp+1)*(k2fp-k1fp+1)
  enddo

  return
end subroutine unfpack1
!
!     ****************************************************************
!
subroutine unfdbackp(ivarn,ac,acf,den,denu,denv,m1,m2,m3,i0,j0  &
   ,ibcon,i1s,i2s,j1s,j2s,k1s,k2s,i2u,j2v,k2w,nf1,nf2,nf3,mynum)
implicit none
integer :: ivarn,m1,m2,m3,i0,j0  &
          ,ibcon,i1s,i2s,j1s,j2s,k1s,k2s,i2u,j2v,k2w,nf1,nf2,nf3,mynum
real :: ac(m1,m2,m3),acf(nf1,nf2,nf3),den(m1,m2,m3)  &
       ,denu(m1,m2,m3),denv(m1,m2,m3)

integer :: i,j,k

!     ivarn = variable types 1- u
!                            2- v
!                            3- w
!                            4- p
!                            5- scalar

if(ivarn.ge.5) then
   do j=j1s,j2s
      do i=i1s,i2s
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)/den(k,i-i0,j-j0)
         enddo
      enddo
   enddo
elseif(ivarn==1) then
   do j=j1s,j2s
      do i=i1s,i2u
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)  &
                 /denu(k,i-i0,j-j0)
         enddo
      enddo
   enddo
elseif(ivarn==2) then
   do j=j1s,j2v
      do i=i1s,i2s
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)  &
                 /denv(k,i-i0,j-j0)
         enddo
      enddo
   enddo
elseif(ivarn==3) then
   do j=j1s,j2s
      do i=i1s,i2s
         do k=k1s,k2w
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)  &
                 /(den(k,i-i0,j-j0)+den(k+1,i-i0,j-j0))
         enddo
      enddo
   enddo
else
   do j=j1s,j2s
      do i=i1s,i2s
         do k=k1s,k2s
            ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                 +acf(k-k1s+1,i-i1s+1,j-j1s+1)
         enddo
      enddo
   enddo
endif

return
end subroutine unfdbackp


!     ****************************************************************
!
subroutine unfdbackp_fixed(ivarn,ac,acf,i1,i2,j1,j2,k1,k2 &
     ,den,denu,denv,m1,m2,m3,i0,j0,i2u,j2v,k2w  &
     ,i1c,i2c,j1c,j2c,k1c,k2c  &
     ,ipm, n_ipm  &
     ,jpm, n_jpm  &
     ,kpm, n_kpm  &
     ,mynum)

  implicit none

  ! Arguments:
  integer, intent(in) :: ivarn, i1, i2, j1, j2, k1, k2, m1, m2, m3, i0, j0,  &
       i2u, j2v, k2w, i1c, i2c, j1c, j2c, k1c, k2c, mynum
  real, intent(in)    :: acf(k1:k2,i1:i2,j1:j2)
  real, intent(inout) :: ac(m1,m2,m3)
  real, intent(in)    :: den(m1,m2,m3), denu(m1,m2,m3), denv(m1,m2,m3)
  integer, intent(in) :: n_ipm, n_jpm, n_kpm
  integer, intent(in) :: ipm(n_ipm),jpm(n_jpm),kpm(n_kpm)

  ! Local Variables:
  integer :: i, j, k, jc, ic, kc

  !     ivarn = variable types 1- u
  !                            2- v
  !                            3- w
  !                            4- p
  !                            5- scalar

  if(ivarn==1) then

     do j=j1,j2
        jc = jpm(j)
        do ic=i1,i2u
           do k=k1,k2
              kc = kpm(k)
              ac(kc,ic-i0,jc-j0)=ac(kc,ic-i0,jc-j0)  &
                   + acf(k,ic,j)

           enddo
        enddo
     enddo
     do j=j1c,j2c
        do i=i1c,i2u
           do k=k1c,k2c
              ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                   /denu(k,i-i0,j-j0)

           enddo
        enddo
     enddo
  elseif(ivarn==2) then

     do jc=j1,j2v
        do i=i1,i2
           ic = ipm(i)
           do k=k1,k2
              kc = kpm(k)
              ac(kc,ic-i0,jc-j0)=ac(kc,ic-i0,jc-j0)  &
                   + acf(k,i,jc)
           enddo
        enddo
     enddo
     do j=j1c,j2v
        do i=i1c,i2c
           do k=k1c,k2c
              ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                   /denv(k,i-i0,j-j0)
           enddo
        enddo
     enddo
  elseif(ivarn==3) then
     do j=j1,j2
        jc = jpm(j)
        do i=i1,i2
           ic = ipm(i)
           do kc=k1,k2w
              ac(kc,ic-i0,jc-j0)=ac(kc,ic-i0,jc-j0)  &
                   +acf(kc,i,j)
           enddo
        enddo
     enddo
     do j=j1c,j2c
        do i=i1c,i2c
           do k=k1c,k2w
              ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                   /(den(k,i-i0,j-j0)+den(k+1,i-i0,j-j0))
           enddo
        enddo
     enddo
  else
     do j=j1,j2
        jc = jpm(j)
        do i=i1,i2
           ic = ipm(i)
           do k=k1,k2
              kc = kpm(k)
              ac(kc,ic-i0,jc-j0)=ac(kc,ic-i0,jc-j0)  &
                   +acf(k,i,j)
           enddo
        enddo
     enddo
     if (ivarn.ge.5) then
        do j=j1c,j2c
           do i=i1c,i2c
              do k=k1c,k2c
                 ac(k,i-i0,j-j0)=ac(k,i-i0,j-j0)  &
                      / den(k,i-i0,j-j0)
              enddo
           enddo
        enddo
     endif
  endif

  return
end subroutine unfdbackp_fixed


!     ****************************************************************
!
