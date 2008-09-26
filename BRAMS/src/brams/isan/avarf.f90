!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine makevarf(ng)

use isan_coms
use mem_grid

implicit none

integer :: ng

!---------------------------------------------------------------+
!    Interpolate model grid from isentropic/sigmaz/surface data
!---------------------------------------------------------------+

!            Vertically interpolate isentropic data to sigma-z levels

call isnsig(nnzp(ng),nnxp(ng),nnyp(ng) ,is_grids(ng)%rr_u(1,1,1)  &
        ,is_grids(ng)%rr_v(1,1,1)      ,is_grids(ng)%rr_t(1,1,1)  &
        ,is_grids(ng)%rr_r(1,1,1)      ,is_grids(ng)%rr_p(1,1,1)  &
        ,grid_g(ng)%topt(1,1),ztn(1,ng),ztop)

!            Compute Exner function on model sigma-z surfaces
!              and change relative humidity to mixing ratio.

call vshyd(nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_p(1,1,1)  &
     ,is_grids(ng)%rr_t(1,1,1)       ,is_grids(ng)%rr_r(1,1,1)  &
     ,grid_g(ng)%topt(1,1),grid_g(ng)%rtgt(1,1),ztn(1,ng))

!          Combine surface analysis with the upper air data.

call visurf(nnzp(ng),nnxp(ng),nnyp(ng) ,is_grids(ng)%rr_u(1,1,1)  &
      ,is_grids(ng)%rr_v(1,1,1)        ,is_grids(ng)%rr_t(1,1,1)  &
      ,is_grids(ng)%rr_r(1,1,1)        ,is_grids(ng)%rr_p(1,1,1)  &
      ,grid_g(ng)%topt(1,1),grid_g(ng)%rtgt(1,1),ztn(1,ng))

is_grids(ng)%rr_slp (1:nnxp(ng),1:nnyp(ng)) =rs_slp (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_sfp (1:nnxp(ng),1:nnyp(ng)) =rs_sfp (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_sft (1:nnxp(ng),1:nnyp(ng)) =rs_sft (1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_snow(1:nnxp(ng),1:nnyp(ng)) =rs_snow(1:nnxp(ng),1:nnyp(ng))
is_grids(ng)%rr_sst (1:nnxp(ng),1:nnyp(ng)) =rs_sst (1:nnxp(ng),1:nnyp(ng))

!          average the velocities to the correct points in the stagger
!             and rotate for polar stereographic transformation.

call varuv(nnzp(ng),nnxp(ng),nnyp(ng),is_grids(ng)%rr_u(1,1,1)  &
                      ,is_grids(ng)%rr_v(1,1,1))

return
end
!
!     **************************************************************
!
subroutine isnsig(n1,n2,n3,uu,vv,tt,rr,pp  &
     ,topt,zt,ztop)
   
use isan_coms
use rconstants

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: uu,vv,tt,rr,pp
real, dimension(n2,n3) :: topt
real, dimension(n1) :: zt
real :: ztop

real, allocatable :: v1(:),v2(:),v3(:),v4(:),v5(:),v6(:)
integer :: i,j,k,nki,ki,kbeg,kend
real :: rtg,wtsz

allocate (v1(n1),v2(nisn),v3(nisn),v4(nisn),v5(nisn),v6(nisn))

do j=1,n3
   do i=1,n2

      rtg=(1.-topt(i,j)/ztop)
      do k=1,n1
         v1(k)=topt(i,j)+zt(k)*rtg
      enddo

      if (guess1st.eq.'PRESS'.or.guess1st.eq.'NC') then
         nki=0
         do ki=1,nisn
!             print*,'v0:',i,j,ki,pi_s(i,j,ki),levth(ki),pi_p(i,j,ki)  &
!               ,pi_u(i,j,ki),pi_v(i,j,ki)
            if(pi_p(i,j,ki).lt.1e20.and.pi_s(i,j,ki).lt.1e20  &
                 .and.pi_u(i,j,ki).lt.1e20.and.pi_v(i,j,ki).lt.1e20)  &
                 then
               nki=nki+1
               v2(nki)=(pi_s(i,j,ki)-cp*levth(ki)  &
                    *(pi_p(i,j,ki)*p00i)**rocp)/g
!                    print*,'v2:',i,j,nki,pi_s(i,j,ki),levth(ki),pi_p(i,j,ki)
                    
               v3(nki)=pi_u(i,j,ki)
               v4(nki)=pi_v(i,j,ki)
               v5(nki)=levth(ki)
               v6(nki)=pi_r(i,j,ki)
            endif
         enddo

         call htint(nki,v3,v2,n1,uu(1,i,j),v1)
         call htint(nki,v4,v2,n1,vv(1,i,j),v1)
         call htint(nki,v5,v2,n1,tt(1,i,j),v1)
         call htint(nki,v6,v2,n1,rr(1,i,j),v1)

      endif


      kbeg=n1+1
      kend=n1+1
      do k=1,nsigz
         if(zt(k).gt.hybbot) then
            kbeg=k-1
            exit
         endif
      enddo

      do k=1,nsigz
         if(zt(k).gt.hybtop) then
            kend=k
            exit
         endif
      enddo

      do k=1,nsigz

         if(k.lt.kbeg) then
            wtsz=sigzwt
         elseif(k.gt.kend) then
            wtsz=0.
         elseif(k.ge.kbeg.and.k.le.kend) then
            wtsz= (zt(kend)-zt(k))  &
                 /(zt(kend)-zt(kbeg))  &
                 *sigzwt
         endif

         uu(k,i,j)=(1.-wtsz)*uu(k,i,j)+wtsz*ps_u(i,j,k)
         vv(k,i,j)=(1.-wtsz)*vv(k,i,j)+wtsz*ps_v(i,j,k)
         tt(k,i,j)=(1.-wtsz)*tt(k,i,j)+wtsz*ps_t(i,j,k)
         rr(k,i,j)=(1.-wtsz)*rr(k,i,j)+wtsz*ps_r(i,j,k)

         if(guess1st.eq.'RAMS') then
            pp(k,i,j)=ps_p(i,j,k)
         endif

      !   if(i==10.and.j==10) &
      !      print*,'avrf, k:',k,tt(k,i,j),rr(k,i,j), pp(k,i,j)

      enddo

   enddo
enddo

deallocate (v1,v2,v3,v4,v5,v6)

return
end

!     ****************************************************************

subroutine visurf(n1,n2,n3,up,vp,thp,rtp,pp,topt,rtgt,zt)

use isan_coms
use rconstants
use therm_lib, only: rehuil,ptrh2rvapil,virtt

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: up,vp,pp,thp,rtp
real, dimension(n2,n3) :: rtgt,topt
real, dimension(n1) :: zt

integer :: i,j,k 

real, allocatable :: v3(:)
real :: zgp,dzgobs,wt,ppp,ttt,rhm


allocate(v3(n1))

do j=1,n3
   do i=1,n2
      if(rs_qual(i,j).gt..5) then

!          Quality of this surface analysis point is okay.

!          Extend the vertical influence of the surface analysis
!            up to a distance of SFCINF above the height of the
!            objectively analyzed surface.

         do k=1,n1
            zgp=topt(i,j)+zt(k)*rtgt(i,j)
            dzgobs=abs(zgp-rs_top(i,j))
            wt=dzgobs/sfcinf
            if(wt.le.1.) then
               if(rs_u(i,j).lt.1.e10)  &
                    up(k,i,j)=up(k,i,j)*wt+rs_u(i,j)*(1.-wt)
               if(rs_v(i,j).lt.1.e10)  &
                    vp(k,i,j)=vp(k,i,j)*wt+rs_v(i,j)*(1.-wt)
               ppp=(pp(k,i,j)/cp)**cpor*p00
               ttt=thp(k,i,j)*pp(k,i,j)/cp
               rhm=rehuil(ppp,ttt,rtp(k,i,j))
               if(rs_r(i,j).lt.1.e10)  &
                    rhm=rhm*wt+rs_r(i,j)*(1.-wt)
               if(rs_t(i,j).lt.1.e10)  &
                    thp(k,i,j)=thp(k,i,j)*wt+rs_t(i,j)*(1.-wt)
               ttt=thp(k,i,j)*pp(k,i,j)/cp
               rtp(k,i,j)=ptrh2rvapil(rhm,ppp,ttt)
            endif
         enddo

         do k=1,n1
            v3(k)=topt(i,j)+zt(k)*rtgt(i,j)
         enddo
         do k=n1-1,1,-1
            pp(k,i,j)=pp(k+1,i,j)+g*(v3(k+1)-v3(k))      &
                 /((virtt(thp(k,i,j),rtp(k,i,j))         &
                  +virtt(thp(k+1,i,j),rtp(k+1,i,j)))*.5)
         enddo

      endif

   enddo
enddo

deallocate(v3)

return
end

!     ****************************************************************

subroutine vshyd(n1,n2,n3,pp,tt,rr,topt,rtg,zt)
     
use isan_coms
use rconstants
use therm_lib, only: ptrh2rvapil,virtt
implicit none
     
integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: tt,rr,pp
real, dimension(n2,n3) :: topt,rtg
real, dimension(n1) :: zt

integer :: i,j,k,lbc,kabc
real :: thmin,ppp,ttt,tmpbc,rvibc,piibc,ziibc,gd2

real, allocatable :: v3(:)

allocate(v3(n1))

if(guess1st.eq.'PRESS' .or. guess1st.eq.'NC') then

!         Find internal boundary condition for the hydrostatic
!           equation

   thmin=tt(n1-1,1,1)
   do j=1,n3
      do i=1,n2
         thmin=min(thmin,tt(n1-1,i,j))
      enddo
   enddo

   do k=nisn,1,-1
      print*,'thmins-',k,levth(k),thmin
      if(float(levth(k)).le.thmin) then
         lbc=k
         go to 121
      endif
   enddo
   stop 'vi-nobc'
121     continue

   print 1212,lbc,levth(lbc)
1212    format(//' Hydrostatic boundary condition set at level'  &
        ,2I4,' K')

endif



do j=1,n3
   do i=1,n2

!         Compute actual model heights

      do k=1,n1
         v3(k)=topt(i,j)+zt(k)*rtg(i,j)
      enddo

      if(guess1st.eq.'PRESS' .or. guess1st.eq.'NC') then

!         Find model level above boundary condition

         do k=1,n1
            if(tt(k,i,j).gt.levth(lbc))go to 141
         enddo
         stop 'vi-nbc2'
141      continue
         kabc=k

!         Integrate P to surface

         piibc=cp*(pi_p(i,j,lbc)*p00i)**rocp
         ziibc=(pi_s(i,j,lbc)-levth(lbc)*piibc)/g
         gd2=2.*g
         pp(kabc-1,i,j)=piibc+gd2*(ziibc-v3(kabc-1))  &
              /(levth(lbc)+tt(kabc-1,i,j))
         do k=kabc-2,1,-1
            pp(k,i,j)=pp(k+1,i,j)+gd2*(v3(k+1)-v3(k))  &
                 /(tt(k+1,i,j)+tt(k,i,j))
         enddo

!         Integrate P to model top

         do k=kabc,n1
            pp(k,i,j)=pp(k-1,i,j)-gd2*(v3(k)-v3(k-1))  &
                 /(tt(k,i,j)+tt(k-1,i,j))
         enddo
!
!         Compute mixing ratio from relative humidity and do final
!           hydrostatic integration.
!
         do k=1,n1
            ppp=(pp(k,i,j)/cp)**cpor*p00
            ttt=tt(k,i,j)*pp(k,i,j)/cp
            rr(k,i,j)=ptrh2rvapil(rr(k,i,j),ppp,ttt)
         enddo

         tmpbc=levth(lbc)*piibc/cp
         rvibc=ptrh2rvapil(pi_r(i,j,lbc),pi_p(i,j,lbc),tmpbc)

         pp(kabc-1,i,j)=piibc+gd2*(ziibc-v3(kabc-1))  &
              /(virtt(real(levth(lbc)),rvibc)  &
               +virtt(tt(kabc-1,i,j),rr(kabc-1,i,j)))
         do k=kabc-2,1,-1
            pp(k,i,j)=pp(k+1,i,j)+gd2*(v3(k+1)-v3(k))  &
                 /(virtt(tt(k+1,i,j),rr(k+1,i,j))  &
                  +virtt(tt(k,i,j),rr(k,i,j)))
         enddo

         do k=kabc,n1
            pp(k,i,j)=pp(k-1,i,j)-gd2*(v3(k)-v3(k-1))  &
                 /((virtt(tt(k,i,j),rr(k,i,j))         &
                   +virtt(tt(k-1,i,j),rr(k-1,i,j))))
         enddo

      elseif (guess1st.eq.'RAMS') then

         do k=1,n1
            ppp=pp(k,i,j)
            ttt=tt(k,i,j)*(ppp/p00)**rocp
            rr(k,i,j)=ptrh2rvapil(rr(k,i,j),ppp,ttt)
         enddo

         pp(1,i,j)= cp*(pp(1,i,j)/p00)**rocp
         do k=2,n1
            pp(k,i,j)=pp(k-1,i,j)-g*(v3(k)-v3(k-1))  &
                 /((virtt(tt(k,i,j),rr(k,i,j))+virtt(tt(k-1,i,j),rr(k-1,i,j)) ) *.5)
       !  if(i==10.and.j==10) &
       !     print*,'avrf222, k:',k,tt(k,i,j),rr(k,i,j), pp(k,i,j)
         enddo

      endif

   enddo
enddo

deallocate(v3)

return
end

!     ******************************************************************

subroutine varuv(n1,n2,n3,up,vp)

implicit none

integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: up,vp

integer :: i,j,k

!             Average to u,v points

do j=1,n3
   do i=1,n2-1
      do k=1,n1
         up(k,i,j)=(up(k,i,j)+up(k,i+1,j))*.5
      enddo
   enddo
enddo

do j=1,n3-1
   do i=1,n2
      do k=1,n1
         vp(k,i,j)=(vp(k,i,j)+vp(k,i,j+1))*.5
      enddo
   enddo
enddo

return
end


!     ****************************************************************

subroutine varfile_nstfeed(ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c &
                          ,nbot,ntop)

use isan_coms

implicit none

integer :: ifm,icm,n1f,n2f,n3f,n1c,n2c,n3c,nbot,ntop


!     Feed back the finer mesh to the coarser mesh.

call fdback(is_grids(icm)%rr_u   (1,1,1),is_grids(ifm)%rr_u   (1,1,1) &
           ,is_grids(icm)%rr_dn0u(1,1,1),is_grids(ifm)%rr_dn0u(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'u',rr_scr1(1))

call fdback(is_grids(icm)%rr_v   (1,1,1),is_grids(ifm)%rr_v   (1,1,1) &
           ,is_grids(icm)%rr_dn0v(1,1,1),is_grids(ifm)%rr_dn0v(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'v',rr_scr1(1))

call fdback(is_grids(icm)%rr_p  (1,1,1),is_grids(ifm)%rr_p  (1,1,1) &
           ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'p',rr_scr1(1))

call fdback(is_grids(icm)%rr_t  (1,1,1),is_grids(ifm)%rr_t  (1,1,1) &
           ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))

call fdback(is_grids(icm)%rr_r  (1,1,1),is_grids(ifm)%rr_r  (1,1,1) &
           ,is_grids(icm)%rr_dn0(1,1,1),is_grids(ifm)%rr_dn0(1,1,1) &
           ,n1c,n2c,n3c,n1f,n2f,n3f,ifm,'t',rr_scr1(1))


if(nbot == 1) then
   call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15 ,is_grids(icm)%rr_u(1,1,1),'U')
   call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15 ,is_grids(icm)%rr_v(1,1,1),'V')
   call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15 ,is_grids(icm)%rr_p(1,1,1),'P')
   call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15 ,is_grids(icm)%rr_t(1,1,1),'T')
   call botset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15 ,is_grids(icm)%rr_r(1,1,1),'T')
endif

if(ntop == 1) then
   call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15,is_grids(icm)%rr_u(1,1,1),is_grids(icm)%rr_u(1,1,1),'U')
   call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15,is_grids(icm)%rr_v(1,1,1),is_grids(icm)%rr_v(1,1,1),'V')
   call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15,is_grids(icm)%rr_p(1,1,1),is_grids(icm)%rr_p(1,1,1),'P')
   call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15,is_grids(icm)%rr_t(1,1,1),is_grids(icm)%rr_t(1,1,1),'T')
   call topset(n1c,n2c,n3c,1,n2c,1,n3c  &
        ,15,is_grids(icm)%rr_r(1,1,1),is_grids(icm)%rr_r(1,1,1),'T')
endif


return
end
