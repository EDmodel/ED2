!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine obj_anal (ctype,ng,nxp,nyp,g_lat,g_lon  &
                    ,polelat,polelon,swx,swy,delx,dely)

use isan_coms

implicit none

character(len=*) :: ctype
integer :: ng,nxp,nyp
real :: polelat,polelon,swx,swy,delx,dely
real, dimension(*) :: g_lat,g_lon
real, allocatable :: ps_scr(:,:)


integer, save :: nvar=1
real, allocatable :: pscr(:,:)

if (ctype == 'isen') then

   call obanl (nxp,nyp,nisn,nvar, pi_u, g_lat, g_lon  &
              ,nsta,nsta,        upi_u,up_lat,up_lon  &
              ,nxp,nyp,           pi_u, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nisn,nvar, pi_v, g_lat, g_lon  &
              ,nsta,nsta,        upi_v,up_lat,up_lon  &
              ,nxp,nyp,           pi_v, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nisn,nvar, pi_p, g_lat, g_lon  &
              ,nsta,nsta,        upi_p,up_lat,up_lon  &
              ,nxp,nyp,           pi_p, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nisn,nvar, pi_s, g_lat, g_lon  &
              ,nsta,nsta,        upi_s,up_lat,up_lon  &
              ,nxp,nyp,           pi_s, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nisn,nvar, pi_r, g_lat, g_lon  &
              ,nsta,nsta,        upi_r,up_lat,up_lon  &
              ,nxp,nyp,           pi_r, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)

elseif (ctype == 'sigz') then

!allocate(pscr(nxp,nyp))
!pscr(1:nxp,1:nyp)=ps_u(1:nxp,1:nyp,2)
!call ezcntr(pscr,nxp,nyp)
   call obanl (nxp,nyp,nsigz,nvar, ps_u, g_lat, g_lon  &
              ,nsta,nsta,         ups_u,up_lat,up_lon  &
              ,nxp,nyp,            ps_u, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
!pscr(1:nxp,1:nyp)=ps_u(1:nxp,1:nyp,20)
!call ezcntr(pscr,nxp,nyp)
!call clsgks
!stop
   call obanl (nxp,nyp,nsigz,nvar, ps_v, g_lat, g_lon  &
              ,nsta,nsta,         ups_v,up_lat,up_lon  &
              ,nxp,nyp,            ps_v, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nsigz,nvar, ps_t, g_lat, g_lon  &
              ,nsta,nsta,         ups_t,up_lat,up_lon  &
              ,nxp,nyp,            ps_t, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nsigz,nvar, ps_p, g_lat, g_lon  &
              ,nsta,nsta,         ups_p,up_lat,up_lon  &
              ,nxp,nyp,            ps_p, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,nsigz,nvar, ps_r, g_lat, g_lon  &
              ,nsta,nsta,         ups_r,up_lat,up_lon  &
              ,nxp,nyp,            ps_r, g_lat, g_lon  &
              ,igridfl,wvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)

elseif (ctype == 'surf') then

              !  Some surface fields do not use a first-guess, so
              !  pass in ps_scr for these

   allocate(ps_scr(nxp,nyp))
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,       sf_ur,sf_lat,sf_lon  &
              ,nxp,nyp,             rs_u, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,       sf_vr,sf_lat,sf_lon  &
              ,nxp,nyp,             rs_v, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,        sf_p,sf_lat,sf_lon  &
              ,nxp,nyp,             rs_p, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,        sf_s,sf_lat,sf_lon  &
              ,nxp,nyp,             rs_s, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,        sf_r,sf_lat,sf_lon  &
              ,nxp,nyp,             rs_r, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
              print*,'===',nssfc,maxsfc
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,        sf_t,sf_lat,sf_lon  &
              ,nxp,nyp,             rs_t, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    ps_scr, g_lat, g_lon  &
              ,maxsfc,nssfc,      sf_top,sf_lat,sf_lon  &
              ,nxp,nyp,           rs_top, g_lat, g_lon  &
              ,0,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)

! These surface variables may have gridded first guess fields, 
!    but they do not have obs now, so pass in ps_scr to obanl as obs data, 
!    set nssfc argument=0 and igridfl=4. This doesn't have to be done, 
!    but, for future compatibility...

   call obanl (nxp,nyp,1,nvar,    rs_slp, g_lat, g_lon  &
              ,maxsfc,0,          ps_scr,sf_lat,sf_lon  &
              ,nxp,nyp,           rs_slp, g_lat, g_lon  &
              ,4,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    rs_sfp, g_lat, g_lon  &
              ,maxsfc,0,          ps_scr,sf_lat,sf_lon  &
              ,nxp,nyp,           rs_sfp, g_lat, g_lon  &
              ,4,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    rs_sft, g_lat, g_lon  &
              ,maxsfc,0,          ps_scr,sf_lat,sf_lon  &
              ,nxp,nyp,           rs_sft, g_lat, g_lon  &
              ,4,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    rs_snow, g_lat, g_lon  &
              ,maxsfc,0,          ps_scr,sf_lat,sf_lon  &
              ,nxp,nyp,           rs_snow, g_lat, g_lon  &
              ,4,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
   call obanl (nxp,nyp,1,nvar,    rs_sst, g_lat, g_lon  &
              ,maxsfc,0,          ps_scr,sf_lat,sf_lon  &
              ,nxp,nyp,           rs_sst, g_lat, g_lon  &
              ,4,swvlnth(ng),respon(ng),gobsep,gobrad,gridwt(ng)  &
              ,polelat,polelon,swx,swy,delx,dely)
              
   deallocate(ps_scr)

else
   print*,'undefined type to obj_anal:', ctype
   stop 'obj_anal'
endif

return
end


!***************************************************************************

subroutine obanl (npx,npy,npz,nvar,dat_p_i    ,p_lat   ,p_lon  &
                 ,max_obs,n_obs   ,dat_obs ,dlat_obs,dlon_obs  &
                 ,nxp,nyp         ,dat_i_i    ,g_lat   ,g_lon  &
                 ,igrdflg,wvlnth,respon,gobsep,gobrad,gridwt  &
                 ,polat,polon,swx,swy,delx,dely)
                 
implicit none                 

integer ::  npx,npy,npz,nvar,max_obs,n_obs,nxp,nyp,igrdflg
real :: wvlnth,respon,gobsep,gobrad,gridwt,polat,polon,swx,swy,delx,dely
real :: dat_p_i(*),p_lat(*),p_lon(*)  &
         ,dat_obs(*),dlat_obs(*),dlon_obs(*)  &
         ,dat_i_i(*),g_lat(*),g_lon(*)
real,allocatable:: dat_qual(:,:)
real,allocatable:: obs_dat(:,:,:,:)  &
                  ,obs_xxx(:),obs_yyy(:),obs_swt(:),obs_scr(:)
                  
integer :: igns,ngs1,iqflag
real :: tm1,tm2,gamma,x4k

! allocate memory for combined pressure grid/station data array

if(igrdflg==1.or.igrdflg==2.or.igrdflg==3) then
   igns=npx*npy+max_obs
elseif(igrdflg==0.or.igrdflg==4) then
   igns=max_obs
endif


print*,'allocating in OBANL-',igns,npz,nvar,igrdflg
allocate(obs_dat(igns,npz,nvar,2))
allocate(obs_xxx(igns))
allocate(obs_yyy(igns))
allocate(obs_swt(igns))
allocate(obs_scr(igns))
allocate(dat_qual(nxp,nyp))

call timing(1,tm1)

call prebarn(npx,npy,npz,nvar  &
     ,dat_p_i,p_lat,p_lon  &
     ,max_obs,dat_obs,dlat_obs,dlon_obs,n_obs  &
     ,igns,obs_dat,obs_xxx,obs_yyy,obs_swt,ngs1  &
     ,igrdflg,gobsep,gobrad  &
     ,gridwt,obs_scr,polat,polon)
     
if(ngs1==0.and.igrdflg/=4) then
   dat_i_i(1:nxp*nyp*npz*nvar)=1.e30
   goto 10
endif


gamma=.3
call bn_parm(wvlnth,respon,gamma,x4k)

! Will we use "quality" weight?
iqflag=0
if(igrdflg==4) iqflag=1

! Get station proximity "quality" weight field.

if(iqflag==1) then
   call bn_qual(dat_qual,nxp,nyp,g_lat,g_lon  &
        ,obs_xxx,obs_yyy ,ngs1,x4k,polat,polon)
endif

if(igrdflg /= 4) then
   dat_i_i(1:nxp*nyp*npz*nvar)=0.
   call bn_pass(dat_i_i,nxp,nyp,npz,nvar,g_lat,g_lon  &
        ,igns,obs_dat,obs_xxx,obs_yyy,obs_swt,ngs1,obs_scr  &
        ,x4k,polat,polon,1,0,dat_qual)
endif

if(ngs1>=1) then
   call stainterp(dat_i_i,nxp,nyp,npz,nvar  &
        ,igns,obs_dat,obs_xxx,obs_yyy  &
        ,ngs1,polat,polon,swx,swy,delx,dely,-1.e-10)

   call bn_pass(dat_i_i,nxp,nyp,npz,nvar,g_lat,g_lon  &
        ,igns,obs_dat,obs_xxx,obs_yyy  &
        ,obs_swt,ngs1,obs_scr  &
        ,x4k*gamma,polat,polon,2,iqflag,dat_qual)
endif

10 continue

call timing(2,tm2)

print 90,ngs1,tm2-tm1,wvlnth,respon
90 format(/' Objective analysis---','  No. obs -',I5  &
         ,'  CPU time - ',F9.4,' secs',/,T31  &
         ,'Barnes parameters -WVLNTH,RESPON-'  &
         ,F8.0,F4.1)

deallocate(obs_dat,obs_xxx,obs_yyy,obs_swt,obs_scr)
deallocate(dat_qual)

RETURN
END


!***************************************************************************

subroutine prebarn (nprx,npry,nisn,nvar  &
                   ,grid,glat,glon  &
                   ,ist,ssdat,sslt,ssln,nnsta  &
                   ,igns,gsdat,gsx,gsy,gswt  &
                   ,ngs,igridfl,gobsep,gobrad  &
                   ,gridwt,scra,polat,polon)
use rconstants, only: spconkm
implicit none
integer ::  nprx,npry,nisn,nvar,nnsta,ist,igns,ngs,igridfl                  
real :: grid(nprx,npry,nisn,*),glat(nprx,npry),glon(nprx,npry)  &
         ,ssdat(ist,nisn,*),sslt(*),ssln(*),scra(*)  &
         ,gsdat(igns,nisn,*),gsx(*),gsy(*),gswt(*),gridwt
integer :: iqf(4)

integer :: n,k,nv,i,j,ns,iq
real :: gobkm,gobrd,gobsep,gobrad,polat,polon,gln,glt

ngs=0
do n=1,nnsta
   ngs=ngs+1
   call ll_xy(sslt(n),ssln(n),polat,polon,gsx(ngs),gsy(ngs))
   gswt(ngs)=1.
   do k=1,nisn
      do nv=1,nvar
         gsdat(ngs,k,nv)=ssdat(n,k,nv)
      enddo
   enddo
enddo

if(igridfl.eq.0.or.igridfl.eq.4)return

if(igridfl.eq.1) then
   do j=1,npry
      do i=1,nprx
         ngs=ngs+1
         call ll_xy(glat(i,j),glon(i,j),polat,polon,gsx(ngs),gsy(ngs))
         gswt(ngs)=gridwt
         do k=1,nisn
            do nv=1,nvar
               gsdat(ngs,k,nv)=grid(i,j,k,nv)
            enddo
         enddo
      enddo
   enddo

   return

elseif(igridfl.eq.2) then
   gobkm=gobsep*spconkm
   gobrd=gobrad*spconkm
   do i=1,nprx
      do j=1,npry
         do iq=1,4
            iqf(iq)=0
         enddo
         do ns=1,nnsta
            scra(ns)=sqrt(((glat(i,j)-sslt(ns))*spconkm)**2  &
                 +((glon(i,j)-ssln(ns))*spconkm  &
                 *cos((glat(i,j)+sslt(ns))*.5*.01745))**2)
         enddo
         do ns=1,nnsta
            if(scra(ns).lt.gobkm) goto 40
         enddo
         do ns=1,nnsta
            if(scra(ns).le.gobrd) then
               gln=glon(i,j)
               glt=glat(i,j)
               if(ssln(ns).le.gln.and.sslt(ns).ge.glt) iqf(1)=1
               if(ssln(ns).ge.gln.and.sslt(ns).ge.glt) iqf(2)=1
               if(ssln(ns).ge.gln.and.sslt(ns).le.glt) iqf(3)=1
               if(ssln(ns).le.gln.and.sslt(ns).le.glt) iqf(4)=1
            endif
         enddo
         if(iqf(1)+iqf(2)+iqf(3)+iqf(4).le.2) then
            ngs=ngs+1
            call ll_xy(glat(i,j),glon(i,j),polat,polon,gsx(ngs),gsy(ngs))
            gswt(ngs)=gridwt
            do k=1,nisn
               do nv=1,nvar
                  gsdat(ngs,k,nv)=grid(i,j,k,nv)
               enddo
            enddo
         endif

         40 continue
      enddo
   enddo
   return

endif

RETURN
END

!***************************************************************************

subroutine bn_parm(wvlnth,respon,gamma,x4k)
use rconstants, only : pi1
implicit none

real :: wvlnth,respon,gamma,x4k

real :: a,o4k,agam,agam1,e1,e2,e3
integer :: n

!       Computes the Barnes(1974) "4K" value given a
!          wavelength, response, and gamma

a=-(pi1/wvlnth)**2
o4k=(wvlnth/pi1)**2
agam=a*gamma
agam1=a*(gamma+1)
n=0
81   continue
n=n+1
e1=exp(a*o4k)
e2=exp(agam*o4k)
e3=exp(agam1*o4k)
x4k=o4k-(e1+e2-e3-respon)/(a*e1+agam*e2-agam1*e3)
if(abs(o4k-x4k).lt.(1./gamma))go to 80
o4k=x4k
if(n.lt.50)go to 81
83   continue
print 85,n,x4k
85   format(' iteration problem in barnes 4k value ',i5,e20.10)
stop 'bn_parm'

80   continue
if(x4k.le.0.)go to 83

return
end
!
!     ******************************************************************
!
subroutine bn_pass(grid,nxx,nyy,nlev,nvar,glat,glon  &
     ,ist,sdat,sxxx,syyy,swt,nsta,scra,x4k  &
     ,polat,polon,npass,iqflag,qual)
     
implicit none

integer :: nxx,nyy,nlev,nvar,ist,nsta,npass,iqflag
real :: grid(nxx,nyy,nlev,nvar),glat(nxx,nyy),glon(nxx,nyy)  &
     ,sdat(ist,nlev,nvar,*),scra(ist,*),sxxx(*),syyy(*),swt(*)  &
     ,qual(nxx,nyy)
real :: x4k,polat,polon


integer :: i,j,k,nv,n
real :: sum,sw,dfactx,dfacty,gxxx,gyyy

dfactx=-1./x4k
dfacty=-1./x4k

do j=1,nyy
   do i=1,nxx

      call ll_xy(glat(i,j),glon(i,j),polat,polon  &
           ,gxxx,gyyy)
      gxxx=gxxx*.001
      gyyy=gyyy*.001

      do n=1,nsta
         scra(n,1)=exp( max(-20.,  &
                  (gyyy-syyy(n)*.001)**2*dfacty  &
                    + (gxxx-sxxx(n)*.001)**2*dfactx) ) *swt(n)
      enddo

      do nv=1,nvar
         do k=1,nlev

            sum=0.
            sw=0.
            do n=1,nsta
               if(sdat(n,k,nv,npass).lt.1e20) then
                  sum=sum+sdat(n,k,nv,npass)*scra(n,1)
                  sw=sw+scra(n,1)
               endif
            enddo

            if(sw.gt.1e-20) then
               if(iqflag.eq.0) then
                   grid(i,j,k,nv)=grid(i,j,k,nv)+sum/sw
               else
                  if(qual(i,j).lt.1.e20) then
                     grid(i,j,k,nv)=grid(i,j,k,nv)  &
                          +sum/sw*qual(i,j)
                  else
!                            Just leave first guess field
!cc                          grid(i,j,k,nv)=grid(i,j,k,nv)
                  endif
               endif

            else
               if(npass.eq.1) then
                  grid(i,j,k,nv)=1.e30
               elseif(npass.eq.2) then
!                         Just leave first guess field
!cc                     grid(i,j,k,nv)=grid(i,j,k,nv)
               endif
            endif

         enddo
      enddo

   enddo
enddo

return
end

!     ******************************************************************

subroutine stainterp(grid,nxx,nyy,nlev,nvar  &
     ,ist,sdat,sxxx,syyy,nsta  &
     ,polat,polon,swx,swy,delx,dely,dval)
implicit none
integer, intent(in) :: nxx,nyy,nlev,nvar,ist,nsta
real,  dimension(nxx,nyy,nlev,nvar) :: grid
real,  dimension(ist,nlev,nvar,2)   :: sdat
real,  dimension(ist)               :: sxxx,syyy
real                                :: polat,polon,swx,swy,delx,dely,dval

integer :: n,nv,k,ni,nj
real :: fi,fj,dx,dy,yval1,sval

do n=1,nsta
   fj=(syyy(n)-swy)/dely
   fi=(sxxx(n)-swx)/delx
   ni=fi+1.000001
   nj=fj+1.000001
   dy=mod(fj,1.)
   dx=mod(fi,1.)

   if(ni.ge.1.and.ni.lt.nxx.and.nj.ge.1.and.nj.lt.nyy) then
      do nv=1,nvar
         do k=1,nlev
            if(sdat(n,k,nv,1).lt.1.e20) then
               yval1=grid(ni,nj,k,nv)  &
                    +(grid(ni+1,nj,k,nv)-grid(ni,nj,k,nv))*dx
               sval=yval1+(grid(ni,nj+1,k,nv)  &
                    +(grid(ni+1,nj+1,k,nv)-grid(ni,nj+1,k,nv))  &
                    *dx-yval1)*dy
               sdat(n,k,nv,2)=sdat(n,k,nv,1)-sval
            else
               sdat(n,k,nv,2)=1.e30
            endif
         enddo
      enddo
   else
      do nv=1,nvar
         do k=1,nlev
            sdat(n,k,nv,2)=dval
         enddo
      enddo
   endif
enddo

return
end


!     ******************************************************************

subroutine bn_qual(qual,nxx,nyy,glat,glon,sxxx,syyy,nsta,x4k,polat,polon)
implicit none
integer :: nxx, nyy,nsta
real :: glat(nxx,nyy),glon(nxx,nyy),sxxx(*),syyy(*),qual(nxx,nyy) &
       ,x4k,polat,polon
       
integer :: n,i,j
real :: dfactx,dfacty,gxxx,gyyy

dfactx=-1./x4k
dfacty=-1./x4k

call azero(nxx*nyy,qual)
do n=1,nsta
   do j=1,nyy
      do i=1,nxx

         call ll_xy(glat(i,j),glon(i,j),polat,polon  &
              ,gxxx,gyyy)

         qual(i,j)=max(qual(i,j),  &
            exp( max(-20.,((gyyy-syyy(n))*.001)**2*dfacty  &
              + ((gxxx-sxxx(n))*.001)**2*dfactx) ) )
      enddo
   enddo
enddo

return
end

