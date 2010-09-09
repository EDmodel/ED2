!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine movenest(ifm,icm)

use mem_tend
use var_tables
use mem_grid
use mem_leaf
use mem_basic
use mem_turb
use mem_scratch

implicit none

integer :: ifm,icm

integer :: i,j,k,nf,nc,i1,i2,j1,j2,ibeg,iend,ishft,jbeg,jend,jshft,ip &
          ,idi,idj,ngr,indfm,nxy4pcm,nxy4p
real :: dtlcm          

if (gridu(ifm) .eq. 0. .and. gridv(ifm) .eq. 0.) return

dtlcm = dtlongn(icm)
idi = 0
idj = 0
dimove(ifm) = dimove(ifm) + gridu(ifm) * dtlcm / deltaxn(icm)
djmove(ifm) = djmove(ifm) + gridv(ifm) * dtlcm / deltayn(icm)

! [Considering only having the moving grid for sequential runs
! because cannot check for time-to-move on fm timesteps.  This limits
! the cm to small timesteps often.  Another difference is that movenest in
! sequential is called for a fm right after the fm is fed back.  Thus, the
! floating point move distances, dimove(ifm) and djmove(ifm), are updated
! every dtlcm, where cm is the immediate parent grid.  However, in parallel,
! dtlcm must be substituted by dtlongn(1) for all nested grids.]

!p1
write(6,557)ifm,idi,idj,ninest(ifm),njnest(ifm)  &
   ,dimove(ifm),djmove(ifm)  &
   ,dtlcm,gridu(ifm),gridv(ifm),deltaxn(nxtnest(ifm))  &
   ,deltayn(nxtnest(ifm))
557  format('1,ifm,idi,idj,ninest,njnest,dimove,djmove'  &
  ,',dtlcm,gridu,gridv,delxn,delyn',/,5i5,7f9.2)
!p2
if (abs(dimove(ifm)) .ge. 1.) then
   idi = int(dimove(ifm))
   ninest(ifm) = ninest(ifm) + idi
   dimove(ifm) = dimove(ifm) - float(idi)

   do ngr = 3,ngrids
      if (nxtnest(ngr) .eq. ifm) then
         ninest(ngr) = ninest(ngr) - idi * nstratx(ifm)
         if (ninest(ngr) .lt. 3 .or. ninest(ngr)  &
            + nnx1(ngr) / nstratx(ngr) .ge. nnx1(ifm)) then
            print*, 'Parent grid x boundary has moved too close'
            print*, 'to grid ',ngr
            print*, 'ninest = ',ninest(ngr)
            iflag = 1
         endif
      endif
   enddo
endif
!
if (abs(djmove(ifm)) .ge. 1.) then
   idj = int(djmove(ifm))
   njnest(ifm) = njnest(ifm) + idj
   djmove(ifm) = djmove(ifm) - float(idj)
!p1
write(6,559)ifm,idi,idj,ninest(ifm),njnest(ifm)  &
   ,dimove(ifm),djmove(ifm)  &
   ,dtlcm,gridu(ifm),gridv(ifm),deltaxn(nxtnest(ifm))  &
   ,deltayn(nxtnest(ifm))
559  format('3,ifm,idi,idj,ninest,njnest,dimove,djmove'  &
  ,',dtlcm,gridu,gridv,delxn,delyn',/,5i5,7f9.2)
!p2
   do ngr = 3,ngrids
      if (nxtnest(ngr) .eq. ifm) then
         njnest(ngr) = njnest(ngr) - idj * nstraty(ifm)
         if (njnest(ngr) .lt. 3 .or. njnest(ngr)  &
            + nny1(ngr) / nstraty(ngr) .ge. nny1(ifm)) then
            print*, 'Parent grid y boundary has moved too close'
            print*, 'to grid ',ngr
            print*, 'njnest = ',njnest(ngr)
            iflag=1
         endif
      endif
   enddo
endif
!
if (ninest(ifm) .lt. 3 .or. njnest(ifm) .lt. 3 .or.  &
   ninest(ifm) + nx1 / nstratx(ifm) .ge. nnx1(icm) .or.  &
   njnest(ifm) + ny1 / nstraty(ifm) .ge. nny1(icm)) then
   print*, 'Grid ',ifm,' has moved too close to cm boundary'
   print*, 'ninest,njnest = ',ninest(ifm),njnest(ifm)
   iflag=1
endif
if (abs(idi) .ge. 2 .or. abs(idj) .ge. 2) then
   print*, 'Grid ',ifm,' moving too fast ',idi,idj
   stop 'griduv'
endif

!  Shift/interpolate non-prognostic sfc variables and reference state

if (idi .ne. 0 .or. idj .ne. 0) then

   call gridset(ifm)
   CALL POLARST(NXP,NYP  &
      ,grid_g(ifm)%GLAT   ,grid_g(ifm)%GLON    & 
      ,grid_g(ifm)%FMAPU  ,grid_g(ifm)%FMAPV   &
      ,grid_g(ifm)%FMAPT  ,grid_g(ifm)%FMAPM   & 
      ,grid_g(ifm)%FMAPUI ,grid_g(ifm)%FMAPVI  &
      ,grid_g(ifm)%FMAPTI ,grid_g(ifm)%FMAPMI  )

   ibeg = 1
   iend = nnxp(ifm)
   ishft = 0
   jbeg = 1
   jend = nnyp(ifm)
   jshft = 0
   if (idi .gt. 0) then
      iend = nnxp(ifm) - nstratx(ifm)
      ishft = nstratx(ifm)
   elseif (idi .lt. 0) then
      ibeg = 1 + nstratx(ifm)
      ishft = -nstratx(ifm)
   endif
   if (idj .gt. 0) then
      jend = nnyp(ifm) - nstraty(ifm)
      jshft = nstraty(ifm)
   elseif (idj .lt. 0) then
      jbeg = 1 + nstraty(ifm)
      jshft = -nstraty(ifm)
   endif
   print*, 'jbeg,jend,ibeg,iend,ishft,jshft'
   print*,  jbeg,jend,ibeg,iend,ishft,jshft


! TOPT
   call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
      ,scratch%scr1,grid_g(icm)%topt)
   call eintp(scratch%scr1,scratch%scr2  &
      ,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
   call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
               ,scratch%scr2,grid_g(ifm)%topt)

   call transfm(nnxp(ifm),nnyp(ifm)                   &
               ,grid_g(ifm)%topt  ,grid_g(ifm)%topu   &
               ,grid_g(ifm)%topv  ,grid_g(ifm)%topm   &
               ,grid_g(ifm)%rtgt  ,grid_g(ifm)%rtgu   &
               ,grid_g(ifm)%rtgv  ,grid_g(ifm)%rtgm   &
               ,grid_g(ifm)%f13u  ,grid_g(ifm)%f13v   &
               ,grid_g(ifm)%f13t  ,grid_g(ifm)%f13m   &
               ,grid_g(ifm)%f23u  ,grid_g(ifm)%f23v   &
               ,grid_g(ifm)%f23t  ,grid_g(ifm)%f23m   &
               ,grid_g(ifm)%dxu   ,grid_g(ifm)%dxv    &
               ,grid_g(ifm)%dxt   ,grid_g(ifm)%dxm    &
               ,grid_g(ifm)%dyu   ,grid_g(ifm)%dyv    &
               ,grid_g(ifm)%dyt   ,grid_g(ifm)%dym    )

   call fmrefs3d(ifm,0)

! LEAF2 section

   nxy4pcm = nnxyp(icm) * (nzg + nzs + 3)
   nxy4p = nnxyp(ifm) * (nzg + nzs + 3)
! TGP
!!   call ae1(nnxysp(ifm),scratch%scr1(1),leaf_g(ifm)%tgp(1,1))

   do ip = 1,npatch
      do j = 1,nnyp(ifm)
         do i = 1,nnxp(ifm)
            indfm = (j-1) * nnxp(ifm) + i-1

            do k = 1,nzg + nzs + 3

!!               a(itgpn(ifm) + (ip-1)*nxy4p  &
!!                  + (k-1)*nnxyp(ifm) +indfm)  &
!!                  = valugp(1,nnxp(icm),nnyp(icm)  &
!!                  ,1,ipm(i,ifm),jpm(j,ifm)  &
!!                  ,a(itgpn(icm) + (ip-1)*nxy4pcm  &
!!                  + (k-1)*nnxyp(icm)))

            enddo
         enddo
      enddo
   enddo

!!   call fmpmoves(nnxp(ifm),nnyp(ifm),nzg+nzs+3,npatch  &
!!      ,a(itgpn(ifm)),a(iscr1)  &
!!      ,ibeg,iend,ishft,jbeg,jend,jshft,1,nzg+nzs+3)
! WGP
!!   call ae1(nnxysp(ifm),a(iscr1),a(iwgpn(ifm)))

   do ip = 1,npatch
      do j = 1,nnyp(ifm)
         do i = 1,nnxp(ifm)
            indfm = (j-1) * nnxp(ifm) + i-1

            do k = 1,nzg + nzs + 3

!!               a(iwgpn(ifm) + (ip-1)*nxy4p  &
!!                  + (k-1)*nnxyp(ifm) +indfm)  &
!!                  = valugp(1,nnxp(icm),nnyp(icm)  &
!!                  ,1,ipm(i,ifm),jpm(j,ifm)  &
!!                  ,a(iwgpn(icm) + (ip-1)*nxy4pcm  &
!!                  + (k-1)*nnxyp(icm)))

            enddo
         enddo
      enddo
   enddo

!!   call fmpmoves(nnxp(ifm),nnyp(ifm),nzg+nzs+3,npatch  &
!!      ,a(iwgpn(ifm)),a(iscr1)  &
!!      ,ibeg,iend,ishft,jbeg,jend,jshft,1,nzg+nzs+3)
! GSF
!!   call ae1(nnxysp(ifm),a(iscr1),a(igsfn(ifm)))

   do ip = 1,npatch
      do j = 1,nnyp(ifm)
         do i = 1,nnxp(ifm)
            indfm = (j-1) * nnxp(ifm) + i-1

            do k = 1,nzg + nzs + 3

!!               a(igsfn(ifm) + (ip-1)*nxy4p  &
!!                  + (k-1)*nnxyp(ifm) +indfm)  &
!!                  = valugp(1,nnxp(icm),nnyp(icm),1  &
!!                  ,ipm(i,ifm),jpm(j,ifm)  &
!!                  ,a(igsfn(icm) + (ip-1)*nxy4pcm  &
!!                  + (k-1)*nnxyp(icm)))

            enddo
         enddo
      enddo
   enddo

!!   call fmpmoves(nnxp(ifm),nnyp(ifm),nzg+nzs+3,npatch  &
!!      ,a(igsfn(ifm)),a(iscr1)  &
!!      ,ibeg,iend,ishft,jbeg,jend,jshft,1,nzg+nzs+3)
! SCHAR
!!   call ae1(nnxysp(ifm),a(iscr1),a(ischarn(ifm)))

   do ip = 1,npatch
      do j = 1,nnyp(ifm)
         do i = 1,nnxp(ifm)
            indfm = (j-1) * nnxp(ifm) + i-1

            do k = 1,nzg + nzs + 3

!!               a(ischarn(ifm) + (ip-1)*nxy4p  &
!!                  + (k-1)*nnxyp(ifm) +indfm)  &
!!                  = valugp(1,nnxp(icm),nnyp(icm),1  &
!!                  ,ipm(i,ifm),jpm(j,ifm)  &
!!                  ,a(ischarn(icm) + (ip-1)*nxy4pcm  &
!!                  + (k-1)*nnxyp(icm)))

            enddo
         enddo
      enddo
   enddo

!!   call fmpmoves(nnxp(ifm),nnyp(ifm),nzg+nzs+3,npatch  &
!!      ,a(ischarn(ifm)),a(iscr1)  &
!!      ,ibeg,iend,ishft,jbeg,jend,jshft,1,nzg+nzs+3)
! SST
!!   call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
!!      ,a(iscr1),a(itgpn(icm)+nnxyp(icm)*(nzg)))
!!   call eintp(a(iscr1),a(iscr2),1,maxnxp,maxnyp  &
!!      ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
!!   call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
!!      ,a(iscr2),a(itgpn(ifm)+nnxyp(ifm)*(nzg)))

!!   call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
!!      ,a(iscr1),a(itgpn(icm)+nnxyp(icm)*(nzg+1)))
!!   call eintp(a(iscr1),a(iscr2),1,maxnxp,maxnyp  &
!!      ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
!!   call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
!!      ,a(iscr2),a(itgpn(ifm)+nnxyp(ifm)*(nzg+1)))

!3  Prognostic 3-D atmospheric variables

   call ae1(nnxyzp(ifm),tend%ut,basic_g(ifm)%up)
   call ae1(nnxyzp(ifm),tend%vt,basic_g(ifm)%vp)
   call ae1(nnxyzp(ifm),tend%wt,basic_g(ifm)%wp)
   call ae1(nnxyzp(ifm),tend%pt,basic_g(ifm)%pp)

do nf = 1,num_scalar(ifm)
   do nc = 1,num_scalar(icm)
      if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
            call ae1(nnxyzp(ifm),scalar_tab(nf,ifm)%var_t  &
               ,scalar_tab(nf,ifm)%var_p)
         endif
      enddo
   enddo

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm)  &
      ,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
      ,maxnzp,maxnxp,maxnyp,ifm,icm  &
      ,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'u'  &
      ,basic_g(icm)%up,basic_g(ifm)%up  &
      ,basic_g(icm)%dn0u,basic_g(ifm)%dn0u  &
      ,scratch%scr1,scratch%scr2  &
      ,grid_g(ifm)%topt,scratch%vt2da  &
      ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm)  &
      ,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
      ,maxnzp,maxnxp,maxnyp,ifm,icm  &
      ,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'v'  &
      ,basic_g(icm)%vp,basic_g(ifm)%vp  &
      ,basic_g(icm)%dn0v,basic_g(ifm)%dn0v  &
      ,scratch%scr1,scratch%scr2  &
      ,grid_g(ifm)%topt,scratch%vt2da  &
      ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm)  &
      ,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
      ,maxnzp,maxnxp,maxnyp,ifm,icm  &
      ,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'w'  &
      ,basic_g(icm)%wp,basic_g(ifm)%wp  &
      ,basic_g(icm)%dn0,basic_g(ifm)%dn0  &
      ,scratch%scr1,scratch%scr2  &
      ,grid_g(ifm)%topt,scratch%vt2da  &
      ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm)  &
      ,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
      ,maxnzp,maxnxp,maxnyp,ifm,icm  &
      ,nnstbot(ifm),nnsttop(ifm),jdim,1,0,0,'t'  &
      ,basic_g(icm)%pp,basic_g(ifm)%pp  &
      ,basic_g(icm)%dn0,basic_g(ifm)%dn0  &
      ,scratch%scr1,scratch%scr2  &
      ,grid_g(ifm)%topt,scratch%vt2da  &
      ,scratch%scr1,scratch%scr1,scratch%scr1,0)

do nf = 1,num_scalar(ifm)
   do nc = 1,num_scalar(icm)
      if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then

            call fmint3(nnzp(icm),nnxp(icm),nnyp(icm)  &
               ,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
               ,maxnzp,maxnxp,maxnyp,ifm,icm  &
               ,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'t'  &
               ,scalar_tab(nc,icm)%var_p &
               ,scalar_tab(nf,ifm)%var_p  &
               ,basic_g(icm)%dn0,basic_g(ifm)%dn0  &
               ,scratch%scr1,scratch%scr2  &
               ,grid_g(ifm)%topt,scratch%vt2da  &
               ,scratch%scr1,scratch%scr1,scratch%scr1,0)

         endif
      enddo
   enddo

   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%up,tend%ut  &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%vp,tend%vt  &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%wp,tend%wt  &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%pp,tend%pt  &
               ,ibeg,iend,ishft,jbeg,jend,jshft)

do nf = 1,num_scalar(ifm)
   do nc = 1,num_scalar(icm)
      if (scalar_tab(nf,ifm)%name == scalar_tab(nc,icm)%name) then
            call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
               ,scalar_tab(nf,ifm)%var_p  &
               ,scalar_tab(nf,ifm)%var_t  &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
         endif
      enddo
   enddo

   call ae1(nnxyzp(ifm), tend%ut  , basic_g(ifm)%uc   )
   call ae1(nnxyzp(ifm), tend%vt  , basic_g(ifm)%vc   )
   call ae1(nnxyzp(ifm), tend%wt  , basic_g(ifm)%wc   )
   call ae1(nnxyzp(ifm), tend%pt  , basic_g(ifm)%pc   )
   call ae1(nnxyzp(ifm), tend%tht , basic_g(ifm)%theta)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)              &
              ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'u'    &
              ,basic_g(icm)%uc,basic_g(ifm)%uc,basic_g(icm)%dn0u,basic_g(ifm)%dn0u      &
              ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                 &
              ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)              &
              ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'v'    &
              ,basic_g(icm)%vc,basic_g(ifm)%vc,basic_g(icm)%dn0v,basic_g(ifm)%dn0v      &
              ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                 &
              ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)              &
              ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'w'    &
              ,basic_g(icm)%wc,basic_g(ifm)%wc,basic_g(icm)%dn0,basic_g(ifm)%dn0        &
              ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                 &
              ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)              &
              ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,0,0,'t'    &
              ,basic_g(icm)%pc,basic_g(ifm)%pc,basic_g(icm)%dn0,basic_g(ifm)%dn0        &
              ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                 &
              ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmint3(nnzp(icm),nnxp(icm),nnyp(icm),nnzp(ifm),nnxp(ifm),nnyp(ifm)              &
              ,maxnzp,maxnxp,maxnyp,ifm,icm,nnstbot(ifm),nnsttop(ifm),jdim,1,1,0,'t'    &
              ,basic_g(icm)%theta,basic_g(ifm)%theta,basic_g(icm)%dn0,basic_g(ifm)%dn0  &
              ,scratch%scr1,scratch%scr2,grid_g(ifm)%topt,scratch%vt2da                 &
              ,scratch%scr1,scratch%scr1,scratch%scr1,0)

   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%uc,tend%ut                   &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%vc,tend%vt                   &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%wc,tend%wt                   &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%pc,tend%pt                   &
               ,ibeg,iend,ishft,jbeg,jend,jshft)
   call fmpmove(nnzp(ifm),nnxp(ifm),nnyp(ifm),basic_g(ifm)%theta,tend%tht               &
               ,ibeg,iend,ishft,jbeg,jend,jshft)

   i1 = 2
   i2 = nnx(ifm)
   j1 = 2
   j2 = nny(ifm)

   if (ishft .lt. 0) i2 = 1 - ishft
   if (ishft .gt. 0) i1 = nnxp(ifm) - ishft
   if (jshft .lt. 0) j2 = 1 - jshft
   if (jshft .gt. 0) j2 = nnyp(ifm) - jshft

   call thermo (nnzp(ifm),nnxp(ifm),nnyp(ifm),i1,i2,j1,j2)

   call tkeinit(nzp,nxp,nyp,i1,i2,j1,j2)

endif
return
end
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fmpmove(n1,n2,n3,anew,aold,ibeg,iend,ishft,jbeg,jend,jshft)
   implicit none
   integer :: n1,n2,n3 ,ibeg,iend,ishft,jbeg,jend,jshft
   real :: anew(n1,n2,n3),aold(n1,n2,n3)

   integer :: i,j,k

   do j = jbeg,jend
      do i = ibeg,iend
         do k = 1,n1
            anew(k,i,j) = aold(k,i+ishft,j+jshft)
         end do
      end do
   end do
   return
end subroutine fmpmove
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fmpmoves(n2,n3,n4,n5,anew,aold,ibeg,iend,ishft,jbeg,jend,jshft,kbeg,kend)
   implicit none
   integer :: n2,n3,n4,n5 ,ibeg,iend,ishft,jbeg,jend,jshft,kbeg,kend
   real :: anew(n2,n3,n4,n5),aold(n2,n3,n4,n5)

   integer :: i,j,k,ip

   do ip = 1,n5
      do k = kbeg,kend
         do j = jbeg,jend
            do i = ibeg,iend
               anew(i,j,k,ip) = aold(i+ishft,j+jshft,k,ip)
            end do
         end do
      end do
   end do
   return
end subroutine fmpmoves
