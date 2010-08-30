!############################# Change Log ##################################
! 2.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine latlon_ps(n1,n2,np1,np2,np3,nvar,ipos,aps,psglob,all  &
     ,glat,glon,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
implicit none
integer :: n1,n2,np1,np2,np3,nvar,ipos,idatelin,iglobew,iglobs,iglobn
real :: xswlat,xswlon,gdatdx,gdatdy
real :: aps(n1,n2,np3,*),psglob(np1+3,np2+2),all(np1,np2,np3,*)  &
     ,glon(n1,n2),glat(n1,n2)

real :: fioff,fjoff,grx,gry,gglon
integer :: np1h,ioff,joff,nv,i,j,k

np1h = np1 / 2

! Determine if psglob scratch array needs to be offset to west and/or south.

ioff = 0
joff = 0

if (iglobew .eq. 1) then
   ioff = 1
   if (iglobs .eq. 1) then
      joff = 1
   endif
endif
fioff = float(ioff)
fjoff = float(joff)

! Loop through all variable types and vertical levels

do nv = 1,nvar
   do k = 1,np3

! Fill main part of scratch array psglob from data

      do j = 1,np2
         do i = 1,np1
            psglob(i+ioff,j+joff) = all(i,j,k,nv)
         enddo
      enddo

! If all longitudes exist in input data, add 1 row on west boundary and 2
! rows on east boundary of input data in scratch array.

      if (iglobew .eq. 1) then
         do j = 1,np2
            psglob(    1,j+joff) = psglob(np1+1,j+joff)
            psglob(np1+2,j+joff) = psglob(    2,j+joff)
            psglob(np1+3,j+joff) = psglob(    3,j+joff)
          enddo
      endif

! If input data includes the south pole, add 1 row to south in scratch array.

      if (iglobs .eq. 1) then
         do i = 1,np1h
            psglob(i,1) = psglob(i+np1h,3)
         enddo
         do i = np1h+1,np1+3
            psglob(i,1) = psglob(i-np1h,3)
         enddo
      endif

! If input data includes the north pole, add 1 row to north in scratch array.

      if (iglobn .eq. 1) then
         do i = 1,np1h
            psglob(i,np2+2) = psglob(i+np1h,np2)
         enddo
         do i = np1h+1,np1+3
            psglob(i,np2+2) = psglob(i-np1h,np2)
         enddo
      endif

! Loop through all i,j points on model grid

      do j = 1,n2
         do i = 1,n1
            gry = (glat(i,j) - xswlat) / gdatdy + 1. + fjoff

            gglon = glon(i,j)
            if (idatelin .eq. 1 .and. gglon .lt. xswlon)  &
               gglon = glon(i,j) + 360.
            grx = (gglon - xswlon) / gdatdx + 1. + fioff

!     horizontally interpolate variable

            call gdtost(psglob,np1+3,np2+2,grx,gry,aps(i,j,k,nv))
            if(ipos==1) aps(i,j,k,nv)=max(aps(i,j,k,nv),0.)
         enddo
      enddo

   enddo
enddo

return
end


!c-------------lambert-conformal to polar stereographic

subroutine lambcon_ps(n1,n2,np1,np2,n3,nvar,ipos,aps,all,glat,glon  &
     ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
implicit none
integer :: n1,n2,np1,np2,n3,nvar,ipos
real :: xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon
real :: aps(n1,n2,n3,*),all(np1,np2,n3,*),glon(n1,n2),glat(n1,n2)

real :: spacing,ri,rj
integer :: i,j,k,nv

!c     real spacing
!c     spacing=  ! (meters, for AWIPS grid 212 this should be 40635.25)
!c
  spacing = gdatdx 

do i = 1, n1
   do j = 1, n2


      call ll2lc(cntlat, cntlon, glat(i,j), glon(i,j), ri, rj  &
         ,xswlat,xswlon, spacing)


      if(ri.lt.1.or.ri.gt.np1.or.rj.lt.1.or.rj.gt.np2)  &
           then
         print*,' gridpoint outside of lambert-con data coverage'
         print*,'   polar i,j,lat,lon-',i,j,glat(i,j),glon(i,j)
         print*,'   computed datafile-relative point:',ri,rj
         stop 'vi-ob'
      endif
!
!     horizontally interpolate variables
!
      do nv=1,nvar
         do k=1,n3
            call gdtost(all(1,1,k,nv),np1,np2  &
                 ,ri,rj,aps(i,j,k,nv))
            if(ipos==1) aps(i,j,k,nv)=max(aps(i,j,k,nv),0.)
         enddo
      enddo

   enddo
enddo

return
end

!***************************************************************************

! *****Subroutine to convert qlat,qlon to x,y (relative to centlat/centlon)
! on Lambert conformal Conic grid***********

subroutine ll_lc(qlat,qlon,centlat,centlon,x,y)

!     ***** Passed variables ********

!     qlat    : latitude of input point (deg)
!     qlon    : longitude of input point (deg)
!     centlat : tangent latitude (deg)
!     centlon : longitude of y axis on lcc (deg)
!     x       : x-coordinate of point on Lambert Conformal grid
!     y       : y-coordinate of point on Lambert Conformal grid

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     f
!     others are scratch variables (see code)

use rconstants, only : erad,pio180
implicit none

real :: qlat,qlon,centlat,centlon,x,y
real :: rcentlat,rcentlon,rqlat,rqlon,cone,Ffactor,rho,rho8999,theta

rcentlat = centlat * pio180
rcentlon = centlon * pio180
rqlat = qlat * pio180
rqlon = qlon * pio180

cone = sin(rcentlat)
Ffactor = cos(rcentlat) * tan(45. * pio180 + rcentlat/2.) ** cone / cone

if (qlat .le. 89.99) then
   rho = Ffactor/tan(45. * pio180 + rqlat / 2.) ** cone
else
   rho8999 = Ffactor/tan(45. * pio180 + 89.99 * pio180 / 2.) ** cone
   rho = 100. * (90. - qlat) * rho8999
endif
   
theta = cone * (rqlon - rcentlon)
x = rho * sin(theta)
y = - rho * cos(theta)

x = x * erad
y = y * erad

return
end

!***************************************************************************

! *****Subroutine to convert x,y (relative to centlon,centlat) to qlat,qlon
! on Lambert conformal Conic grid***********

subroutine lc_ll(qlat,qlon,centlat,centlon,x,y)

!     ***** Passed variables ********

!     qlat    : latitude of input point (deg)
!     qlon    : longitude of input point (deg)
!     centlat : tangent latitude (deg)
!     centlon : longitude of y axis on lcc (deg)
!     x       : x-coordinate of point on Lambert Conformal grid
!     y       : y-coordinate of point on Lambert Conformal grid

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     f
!     others are scratch variables (see code)
use rconstants, only: pio180,erad
implicit none

real, parameter :: radearth = 6370997           ! (USGS Paper 1395)

real :: qlat,qlon,centlat,centlon,x,y
real :: rcentlat,rcentlon,rqlat,rqlon,cone,Ffactor,rho,rho8999,theta

rcentlat = centlat * pio180
rcentlon = centlon * pio180
rqlat = qlat * pio180
rqlon = qlon * pio180

cone = sin(rcentlat)
Ffactor = cos(rcentlat) * tan(45. * pio180 + rcentlat/2.) ** cone / cone

if (qlat .le. 89.99) then
   rho = Ffactor/tan(45. * pio180 + rqlat / 2.) ** cone
else
   rho8999 = Ffactor/tan(45. * pio180 + 89.99 * pio180 / 2.) ** cone
   rho = 100. * (90. - qlat) * rho8999
endif
   
theta = cone * (rqlon - rcentlon)
x = rho * sin(theta)
y = - rho * cos(theta)

x = x * erad
y = y * erad

return
end
!
! *****Subroutine to convert lat/lon pair to real i,j pair on Lambert conformal Conic grid***********
!

subroutine ll2lc(dcentlat, dcentlon, dinlat, dinlon, ri, rj, doriglat, doriglon, spacing)

use rconstants, only : erad, pio180,pio4
implicit none

!     ***** Passed variables ********

!     dcentlat:tangent latitude (deg)
!     dcentlon:longitude of y axis on lcc (deg)
!     dinlat:latitude of input point (deg)
!     dinlon:longitude of input point (deg)
!     ri:place for output i-index (altered)
!     rj:place for output j-index (altered)
!     doriglat: latitude of (1,1) point on lcc grid (deg)
!     doriglon: longitude of (1,1) point on lcc grid (deg)
!     spacing: grid spacing on lcc grid (meters)

real dcentlat, dcentlon, dinlat, dinlon, ri, rj, doriglat,  &
     doriglon, spacing

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     r*: radian equivalents of d* above
!     others are scratch variables (see code)

real cone, rcentlat, rcentlon, rinlat, rinlon, roriglat,  &
     roriglon, xorig, yorig, rho, theta, Ffactor

!     ***** Parameters **************

!     radearth: mean radius of earth (m)


!     ***** Code ********************

!     Convert degrees to radians

rcentlat=dcentlat*pio180
rcentlon=dcentlon*pio180
rinlat=dinlat*pio180
rinlon=dinlon*pio180
roriglat=doriglat*pio180
roriglon=doriglon*pio180

!     Compute projection constants

cone=sin(rcentlat)
Ffactor=cos(rcentlat)*tan(pio4 + rcentlat/2.0)**cone/cone


!     Compute location of origin of lcc grid (1,1)

rho= Ffactor/tan(pio4 + roriglat/2.0)**cone
theta= cone*(roriglon-rcentlon)
xorig= rho*sin(theta)
yorig= -1.0*rho*cos(theta)


!     Compute location of input lat and lon

rho= Ffactor/tan(pio4 + rinlat/2.0)**cone
theta= cone*(rinlon-rcentlon)
ri= rho*sin(theta)
rj= -1.0*rho*cos(theta)

!     Compute location relative to origin

ri= ri - xorig
rj= rj - yorig

!     Convert to index

ri= ri*erad/spacing + 1.0
rj= rj*erad/spacing + 1.0

!     Done

return
end subroutine ll2lc


! *****Subroutine to convert real i,j pair on Lambert
!onformal Conic grid to lat/lon ***********
!

subroutine lc2ll(dcentlat, dcentlon, qlat, qlon, ri, rj,  &
                 doriglat, doriglon, spacing)

use rconstants, only : erad,eradi,pio180,pio4,halfpi,onerad
implicit none

!     ***** Passed variables ********

!     dcentlat:tangent latitude (deg)
!     dcentlon:longitude of y axis on lcc (deg)
!     qlat:output latitude (deg)
!     qlon:output longitude (deg)
!     ri:input i-index
!     rj:input j-index
!     doriglat: latitude of (1,1) point on lcc grid (deg)
!     doriglon: longitude of (1,1) point on lcc grid (deg)
!     spacing: grid spacing on lcc grid (meters)

real dcentlat, dcentlon, qlat, qlon, ri, rj, doriglat,  &
     doriglon, spacing

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     r*: radian equivalents of d* above
!     others are scratch variables (see code)

real cone, rcentlat, rcentlon, rinlat, rinlon, roriglat,  &
     roriglon, xorig, yorig, rho, theta, Ffactor, rho0


!     ***** Code ********************

!     Convert degrees to radians

rcentlat=dcentlat*pio180
rcentlon=dcentlon*pio180
roriglat=doriglat*pio180
roriglon=doriglon*pio180

!     Compute projection constants

cone=sin(rcentlat)
Ffactor=cos(rcentlat)*(tan(pio4 + rcentlat/2.0))**cone/cone

!print*, 'cone,ffactor',cone,ffactor


!     Compute location of origin of lcc grid (1,1)

rho0= Ffactor/(tan(pio4 + rcentlat/2.0))**cone
rho = Ffactor/(tan(pio4 + roriglat/2.0))**cone
theta= cone*(roriglon-rcentlon)
xorig= rho*sin(theta)
yorig= rho0 - rho*cos(theta)

!print*, 'rho,theta,xorig,yorig',rho,theta,xorig,yorig


! convert ri,rj to radians

ri = (ri - 1.) * spacing * eradi
rj = (rj - 1.) * spacing * eradi

! get ri,rj relative to rcentlat,rcentlon

ri = ri + xorig
rj = rj + yorig

!     Compute location of input lat and lon

rho = sqrt(ri ** 2 + (rho0 - rj) ** 2)
theta = atan2(ri,(rho0 - rj))
qlon = theta / cone + rcentlon
qlat = 2. * atan((Ffactor/rho) ** (1./cone)) - halfpi

! convert to degrees
 
qlat = qlat * onerad
qlon = qlon * onerad

return
end

!=======================================================================

subroutine rotate_winds(rot_type,n1,n2,n3,u1,v1,u2,v2  &
                       ,qlat,qlon,reflat1,reflon1,reflat2,reflon2)

!    Rotate 3-D wind field depending on value of rot_type
!       u1,v1           :  input winds
!       reflat1,reflon1 :  reference lat-lon of input wind projection
!       u2,v2           :  output winds
!       reflat2,reflon2 :  reference lat-lon of output wind projection
!       qlat, qlon      :  lat-lon arrays at grid points
                       
implicit none

character (len=*) :: rot_type
integer :: n1,n2,n3
real, dimension(n1,n2,n3) :: u1,v1,u2,v2
real, dimension(n1,n2) :: qlat,qlon
real :: reflat1,reflon1,reflat2,reflon2

real :: uu,vv
integer :: i,j,k

print*,'in rot:',rot_type
print*,'in rot:',n1,n2,n3,reflat1,reflon1,reflat2,reflon2
if(n1 > 1 .and. n2 > 1)print*,'in rot:',qlat(n1/2,n2/2),qlon(n1/2,n2/2)

if (rot_type == 'll_rps') then
   !   lat-lon (earth-relative) winds  to RAMS polar-stereo
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            call uevetouv(u2(i,j,k),v2(i,j,k),u1(i,j,k),v1(i,j,k)  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
         enddo
      enddo
   enddo
elseif (rot_type == 'rps_ll') then
   !   RAMS polar-stereo  to lat-lon (earth-relative) winds
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            call uvtoueve(u1(i,j,k),v1(i,j,k),u2(i,j,k),v2(i,j,k)  &
                         ,qlat(i,j),qlon(i,j),reflat1,reflon1)
         enddo
      enddo
   enddo
elseif (rot_type == 'lc_rps') then
   !   lambert-conformal-relative winds  to RAMS polar-stereo
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            call uvlc_uvll(u1(i,j,k),v1(i,j,k),uu,vv  &
                          ,qlat(i,j),qlon(i,j),reflat1,reflon1)
            call uevetouv(u2(i,j,k),v2(i,j,k),uu,vv  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
            if(i==n1/2.and.j==n2/2.and.k==8) then
               print*,'input:',u1(i,j,k),v1(i,j,k),qlat(i,j),qlon(i,j)
               print*,'earth:',uu,vv,reflat1,reflon1
               print*,'ps:',u2(i,j,k),v2(i,j,k),reflat2,reflon2
            endif
         enddo
      enddo
   enddo
elseif (rot_type == 'ps_rps') then
   !   polar-stereo-relative winds  to RAMS polar-stereo
   do j=1,n2
      do i=1,n1
         do k=1,n3
            if (u1(i,j,k) > 1.e20 .or. v1(i,j,k) > 1.e20) then
               u2(i,j,k)=1.e30
               v2(i,j,k)=1.e30
               cycle
            endif
            call uvtoueve(u1(i,j,k),v1(i,j,k),uu,vv  &
                         ,qlat(i,j),qlon(i,j),reflat1,reflon1)
            call uevetouv(u2(i,j,k),v2(i,j,k),uu,vv  &
                         ,qlat(i,j),qlon(i,j),reflat2,reflon2)
         enddo
      enddo
   enddo
endif

return
end



!=======================================================================

subroutine rotate_winds_obs(rot_type,n1,n2,u1,v1,u2,v2  &
                       ,qlat,qlon,reflat1,reflon1,reflat2,reflon2)

!    Rotate rawin or sfc wind obs depending on value of rot_type
!       u1,v1           :  input winds
!       reflat1,reflon1 :  reference lat-lon of input wind projection
!       qlat, qlon      :  lat-lon of obs
!       u2,v2           :  output winds
!       reflat2,reflon2 :  reference lat-lon of output wind projection
                       
implicit none

character (len=*) :: rot_type
integer :: n1,n2
real, dimension(n1,n2) :: u1,v1,u2,v2
real, dimension(n1) :: qlat,qlon
real :: reflat1,reflon1,reflat2,reflon2

real :: uu,vv
integer :: i,k

if (rot_type == 'll_rps') then
   !   lat-lon (earth-relative) winds  to RAMS polar-stereo
   do k=1,n2
      do i=1,n1
            call uevetouv(u2(i,k),v2(i,k),u1(i,k),v1(i,k)  &
                         ,qlat(i),qlon(i),reflat2,reflon2)
      enddo
   enddo
elseif (rot_type == 'rps_ll') then
   !   RAMS polar-stereo  to lat-lon (earth-relative) winds
   do k=1,n2
      do i=1,n1
            call uvtoueve(u1(i,k),v1(i,k),u2(i,k),v2(i,k)  &
                         ,qlat(i),qlon(i),reflat1,reflon1)
      enddo
   enddo
elseif (rot_type == 'lc_rps') then
   !   lambert-conformal-relative winds  to RAMS polar-stereo
   do k=1,n2
      do i=1,n1
            call uvlc_uvll(u1(i,k),v1(i,k),uu,vv  &
                          ,qlat(i),qlon(i),reflat1,reflon1)
            call uevetouv(u2(i,k),v2(i,k),uu,vv  &
                         ,qlat(i),qlon(i),reflat2,reflon2)
      enddo
   enddo
elseif (rot_type == 'ps_rps') then
   !   polar-stereo-relative winds  to RAMS polar-stereo
   do k=1,n2
      do i=1,n1
            call uvtoueve(u1(i,k),v1(i,k),uu,vv  &
                         ,qlat(i),qlon(i),reflat1,reflon1)
            call uevetouv(u2(i,k),v2(i,k),uu,vv  &
                         ,qlat(i),qlon(i),reflat2,reflon2)
      enddo
   enddo
endif

return
end

!***************************************************************************

subroutine uvlc_uvll(ulc,vlc,ull,vll,qlat,qlon,cntlat,cntlon)

!    Rotate lambert-conformal-relative winds to earth-relative components

implicit none
real :: ulc,vlc,ull,vll,qlat,qlon,cntlat,cntlon,x0,y0,x1,y1,angle

call ll_rotate_lc(qlat,qlon   ,cntlat,cntlon,x0,y0)
call ll_rotate_lc(qlat,qlon+.1,cntlat,cntlon,x1,y1)

angle = -atan2(y1-y0,x1-x0)
ull = ulc * cos(angle) - vlc * sin(angle)
vll = ulc * sin(angle) + vlc * cos(angle)

return
end

!***************************************************************************

subroutine uvll_uvlc(ulc,vlc,ull,vll,qlat,qlon,cntlat,cntlon)

!    Rotate earth-relative winds to lambert-conformal-relative components

implicit none
real :: ulc,vlc,ull,vll,qlat,qlon,cntlat,cntlon,x0,y0,x1,y1,angle

call ll_rotate_lc(qlat,qlon   ,cntlat,cntlon,x0,y0)
call ll_rotate_lc(qlat,qlon+.1,cntlat,cntlon,x1,y1)

angle = atan2(y1-y0,x1-x0)
ulc = ull * cos(angle) - vll * sin(angle)
vlc = ull * sin(angle) + vll * cos(angle)

return
end

!***************************************************************************

! *****Subroutine to convert qlat,qlon to x,y (relative to centlat/centlon)
! on Lambert conformal Conic grid***********

subroutine ll_rotate_lc(qlat,qlon,centlat,centlon,x,y)

!     ***** Passed variables ********

!     qlat    : latitude of input point (deg)
!     qlon    : longitude of input point (deg)
!     centlat : tangent latitude (deg)
!     centlon : longitude of y axis on lcc (deg)
!     x       : x-coordinate of point on Lambert Conformal grid
!     y       : y-coordinate of point on Lambert Conformal grid

!     ***** Local variables *********

!     cone: cone constant
!     Ffactor: another constant
!     f
!     others are scratch variables (see code)
use rconstants, only : erad, pio180
implicit none
real :: qlat,qlon,centlat,centlon,x,y

real :: rcentlat,rcentlon,rqlat,rqlon,cone,Ffactor,rho,rho8999,theta

rcentlat = centlat * pio180
rcentlon = centlon * pio180
rqlat = qlat * pio180
rqlon = qlon * pio180

cone = sin(rcentlat)
Ffactor = cos(rcentlat) * tan(45. * pio180 + rcentlat/2.) ** cone / cone

if (qlat .le. 89.99) then
   rho = Ffactor/tan(45. * pio180 + rqlat / 2.) ** cone
else
   rho8999 = Ffactor/tan(45. * pio180 + 89.99 * pio180 / 2.) ** cone
   rho = 100. * (90. - qlat) * rho8999
endif
   
theta = cone * (rqlon - rcentlon)
x = rho * sin(theta)
y = - rho * cos(theta)

x = x * erad
y = y * erad

return
end


!***************************************************************************

subroutine trueps60_ps (n1,n2,np1,np2,n3,nvar,ipos,aps,all,glat,glon  &
                       ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
implicit none
integer :: n1,n2,np1,np2,n3,nvar,ipos
real :: xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon
real :: aps(n1,n2,n3,*),all(np1,np2,n3,*),glon(n1,n2),glat(n1,n2)

real :: spacing,ri,rj
integer :: i,j,k,nv
     
! True polar stereographic to RAMS rotated polar stereographic


!-----------------------------------------------------------------------
! This code could be used if we made new ll_xy and xy_ll routines
!   that used NCEP constants, mainly
!   DATA  RERTH /6.3712E+6/,PI/3.1416/
!   DATA  SS60 /1.86603/

!  call ll_xy(xswlat,xswlon,cntlat,cntlon,xw,ys)
!  Convert distance true at 60N to distance true at polepoint
!  fact=0.5*(1.+sqrt(3.)/2.)
!  xe=xw+float(np1-1)*splamx/fact
!  yn=ys+float(np2-1)*splamy/fact
!  print*,'xw,xe,ys,yn=',xw,xe,ys,yn
!  call xy_ll(xnelat,xnelon,cntlat,cntlon,xe,yn)
!  print*,'SW,NE corner=',xswlat,xswlon,xnelat,xnelon
!
!  call xy_ll(xnwlat,xnwlon,cntlat,cntlon,xw,yn)
!  call xy_ll(xselat,xselon,cntlat,cntlon,xe,ys)
!  print*,'NW,SE corner=',xnwlat,xnwlon,xselat,xselon
!-----------------------------------------------------------------------

!real spacing
!spacing=  ! (meters, for AWIPS grid 212 this should be 40635.25)

  spacing = gdatdx

do i = 1, n1
   do j = 1, n2

      ! RAMS lat/lon coordinates on the other polar stereo coord
      call w3fb06(glat(i,j),glon(i,j),xswlat,xswlon,spacing,cntlon,ri,rj)

      if(ri.lt.1.or.ri.gt.np1.or.rj.lt.1.or.rj.gt.np2) then
         print*,' gridpoint outside true Polr Ster data coverage'
         print*,'   polar i,j,lat,lon-',i,j,glat(i,j),glon(i,j)
         print*,'   computed datafile-relative point:',ri,rj
         stop 'vi-ob in sub trueps60_ps in asti.f90'
      endif

     ! horizontally interpolate variables
      do nv=1,nvar
         do k=1,n3
            call gdtost(all(1,1,k,nv),np1,np2,ri,rj,aps(i,j,k,nv))
            if(ipos==1) aps(i,j,k,nv)=max(aps(i,j,k,nv),0.)
         enddo
      enddo

   enddo
enddo

return
end

!***************************************************************************

subroutine w3fb07 (xi,xj,alat1,alon1,dx,alonv,alat,alon)
use rconstants, only : pio180,onerad,erad,ss60
implicit none
 
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM:  W3FB07        GRID COORDS TO LAT/LON FOR GRIB
!   PRGMMR: STACKPOLE        ORG: NMC42       DATE:88-04-05
!
! ABSTRACT: CONVERTS THE COORDINATES OF A LOCATION ON EARTH GIVEN IN A
!   GRID COORDINATE SYSTEM OVERLAID ON A POLAR STEREOGRAPHIC MAP PRO-
!   JECTION TRUE AT 60 DEGREES N OR S LATITUDE TO THE
!   NATURAL COORDINATE SYSTEM OF LATITUDE/LONGITUDE
!   W3FB07 IS THE REVERSE OF W3FB06.
!   USES GRIB SPECIFICATION OF THE LOCATION OF THE GRID
!
! PROGRAM HISTORY LOG:
!   88-01-01  ORIGINAL AUTHOR:  STACKPOLE, W/NMC42
!
! USAGE:  CALL W3FB07(XI,XJ,ALAT1,ALON1,DX,ALONV,ALAT,ALON)
!   INPUT ARGUMENT LIST:
!     XI       - I COORDINATE OF THE POINT  REAL*4
!     XJ       - J COORDINATE OF THE POINT  REAL*4
!     ALAT1    - LATITUDE  OF LOWER LEFT POINT OF GRID (POINT 1,1)
!                LATITUDE <0 FOR SOUTHERN HEMISPHERE; REAL*4
!     ALON1    - LONGITUDE OF LOWER LEFT POINT OF GRID (POINT 1,1)
!                  EAST LONGITUDE USED THROUGHOUT; REAL*4
!     DX       - MESH LENGTH OF GRID IN METERS AT 60 DEG LAT
!                 MUST BE SET NEGATIVE IF USING
!                 SOUTHERN HEMISPHERE PROJECTION; REAL*4
!                   190500.0 LFM GRID,
!                   381000.0 NH PE GRID, -381000.0 SH PE GRID, ETC.
!     ALONV    - THE ORIENTATION OF THE GRID.  I.E.,
!                THE EAST LONGITUDE VALUE OF THE VERTICAL MERIDIAN
!                WHICH IS PARALLEL TO THE Y-AXIS (OR COLUMNS OF
!                THE GRID) ALONG WHICH LATITUDE INCREASES AS
!                THE Y-COORDINATE INCREASES.  REAL*4
!                   FOR EXAMPLE:
!                   255.0 FOR LFM GRID,
!                   280.0 NH PE GRID, 100.0 SH PE GRID, ETC.
!
!   OUTPUT ARGUMENT LIST:
!     ALAT     - LATITUDE IN DEGREES (NEGATIVE IN SOUTHERN HEMI.)
!     ALON     - EAST LONGITUDE IN DEGREES, REAL*4
!
!   Remarks: FORMULAE AND NOTATION LOOSELY BASED ON HOKE, HAYES,
!     AND RENNINGER'S "MAP PROJECTIONS AND GRID SYSTEMS...", MARCH 1981
!     AFGWC/TN-79/003
!
! Attributes:
!   Language: IBM VS FORTRAN
!   Machine:  NAS
!
!$$$

real ::xi,xj,alat1,alon1,dx,alonv,alat,alon,h,dxl,reflon &
      ,rebydx,ala1,rmll,alo1,polei,polej,xx,yy,r2,gi2,arccos

! Preliminary variables and redifinitions

! h = 1 for northern hemisphere; = -1 for southern

! reflon is longitude upon which the positive x-coordinate
! drawn through the pole and to the right lies
! rotated around from orientation (y-coordinate) longitude
! differently in each hemisphere

if(dx.lt.0) then
   h = -1.
   dxl = -dx
   reflon = alonv - 90.
else
   h = 1.
   dxl = dx
   reflon = alonv - 270.
endif

rebydx = erad/dxl

! Radius to lower left hand (ll) corner

ala1 =  alat1 * pio180
rmll = rebydx * cos(ala1) * ss60/(1. + h * sin(ala1))

! Use ll point info to locate pole point

alo1 = (alon1 - reflon) * pio180
polei = 1. - rmll * cos(alo1)
polej = 1. - h * rmll * sin(alo1)

! Radius to the i,j point (in grid units)

xx =  xi - polei
yy = (xj - polej) * h
r2 =  xx**2 + yy**2

! Now the magic formulae

if(r2.eq.0) then
   alat = h * 90.
   alon = reflon
else
   gi2 = (rebydx * ss60)**2
   alat = onerad * h * asin((gi2 - r2)/(gi2 + r2))
   arccos = acos(xx/sqrt(r2))
   if(yy.gt.0) then
      alon = reflon + onerad * arccos
   else
      alon = reflon - onerad * arccos
   endif
endif
!if(alon.lt.0) alon = alon + 360.

return
end

!***************************************************************************

SUBROUTINE W3FB06 (ALAT,ALON,ALAT1,ALON1,DX,ALONV,XI,XJ)
use rconstants, only : erad, pio180,ss60
implicit none
real :: ALAT,ALON,ALAT1,ALON1,DX,ALONV,XI,XJ

real :: h,dx1,reflon,rebydx,ala1,alo1,polei,polej,ala,rm,alo,dxl,rmll

! SUBPROGRAM:  W3FB06        LAT/LON TO POLA (I,J) FOR GRIB
!   PRGMMR: STACKPOLE        ORG: NMC42       DATE:88-04-05

IF(DX.LT.0) THEN
   H = -1.
   DXL = -DX
   REFLON = ALONV - 90.
ELSE
   H = 1.
   DXL = DX
   REFLON = ALONV - 270.
ENDIF

REBYDX = erad/DXL

ALA1 =  ALAT1 * pio180
RMLL = REBYDX * COS(ALA1) * SS60/(1. + H * SIN(ALA1))

ALO1 = (ALON1 - REFLON) * pio180
POLEI = 1. - RMLL * COS(ALO1)
POLEJ = 1. - H * RMLL * SIN(ALO1)

ALA =  ALAT * pio180
RM = REBYDX * COS(ALA) * SS60/(1. + H * SIN(ALA))

ALO = (ALON - REFLON) * pio180
XI = POLEI + RM * COS(ALO)
XJ = POLEJ + H * RM * SIN(ALO)

RETURN
END

!***************************************************************************

subroutine trueps60_ps_rot (n1,n2,n3,nvar,aps,glat,glon  &
                           ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon &
                           ,polelat,polelon)
implicit none
integer :: n1,n2,n3,nvar
real :: xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon,polelat,polelon    
real :: aps(n1,n2,n3,*),glon(n1,n2),glat(n1,n2)

integer :: i,j,k
real :: spacing,ue,ve

! Rotate winds from true polar stereographic to RAMS rotated polar 
!    stereographic

spacing = gdatdx

do i=1,n1
   do j=1,n2
      do k=1,n3
         call trueps60_uvtoueve(aps(i,j,k,1),aps(i,j,k,2),ue,ve  &
                               ,glat(i,j),glon(i,j),xswlat,xswlon  &
                               ,spacing,cntlon)
         call uevetouv(aps(i,j,k,1),aps(i,j,k,2),ue,ve  &
                      ,glat(i,j),glon(i,j),polelat,polelon)
      enddo
   enddo
enddo

return
end

!***************************************************************************

subroutine trueps60_uevetouv (u,v,ue,ve,qlat,qlon,xswlat,xswlon &
                             ,spacing,cntlon)
implicit none
real :: u,v,ue,ve,qlat,qlon,xswlat,xswlon,spacing,cntlon
real :: x0,y0,x1,y1,angle

call w3fb06(qlat,qlon,xswlat,xswlon,spacing,cntlon,x0,y0)
call w3fb06(qlat,qlon+0.1,xswlat,xswlon,spacing,cntlon,x1,y1)
angle = atan2(y1-y0,x1-x0)
u = ue * cos(angle) - ve * sin(angle)
v = ue * sin(angle) + ve * cos(angle)

return
end

!***************************************************************************

subroutine trueps60_uvtoueve (u,v,ue,ve,qlat,qlon,xswlat,xswlon &
                             ,spacing,cntlon)
implicit none
real :: u,v,ue,ve,qlat,qlon,xswlat,xswlon ,spacing,cntlon
real :: x0,y0,x1,y1,angle

call w3fb06(qlat,qlon,xswlat,xswlon,spacing,cntlon,x0,y0)
call w3fb06(qlat,qlon+0.1,xswlat,xswlon,spacing,cntlon,x1,y1)
angle = -atan2(y1-y0,x1-x0)
ue = u * cos(angle) - v * sin(angle)
ve = u * sin(angle) + v * cos(angle)

return
end

!***************************************************************************
!***************************************************************************
subroutine ps_ps (n1,n2,np1,np2,n3,nvar,aps,all,glat,glon  &
                 ,polelat,polelon,ppolelat,ppolelon,gglat1,gglon1,ddx,ddy)
implicit none
integer :: n1,n2,np1,np2,n3,nvar
real :: polelat,polelon,ppolelat,ppolelon,gglat1,gglon1,ddx,ddy
real :: aps(n1,n2,n3,*),all(np1,np2,n3,*),glon(n1,n2),glat(n1,n2)

integer :: i,j,k,nv
real :: xxi,xxj,xi,xj,ri,rj

! Rotate winds from rotated polar stereographic grid to RAMS rotated
!    polar stereographic grid


! The incoming polar stereo coord x,y position
call ll_xy(gglat1,gglon1,ppolelat,ppolelon,xxi,xxj)

do i=1,n1
   do j=1,n2

      ! RAMS lat/lon coordinates on the incoming polar stereo coord
      call ll_xy(glat(i,j),glon(i,j),ppolelat,ppolelon,xi,xj)
      ri = 1. + (xi - xxi) / ddx
      rj = 1. + (xj - xxj) / ddy

      if(ri.lt.1.or.ri.gt.np1.or.rj.lt.1.or.rj.gt.np2) then
         print*,' gridpoint outside Polr Ster data coverage'
         print*,'   polar i,j,lat,lon-',i,j,glat(i,j),glon(i,j)
         print*,'   computed datafile-relative point:',ri,rj
         stop 'vi-ob in sub ps_ps in asti.f90'
      endif

      ! horizontally interpolate variables
      do nv=1,nvar
         do k=1,n3
            call gdtost(all(1,1,k,nv),np1,np2,ri,rj,aps(i,j,k,nv))
         enddo
      enddo

   enddo
enddo

return
end

!***************************************************************************
subroutine ps_ps_rot (n1,n2,np1,np2,n3,nvar,aps,all,glat,glon  &
                     ,polelat,polelon,ppolelat,ppolelon)
implicit none
integer :: n1,n2,np1,np2,n3,nvar
real :: polelat,polelon,ppolelat,ppolelon
real :: aps(n1,n2,n3,*),all(np1,np2,n3,*),glon(n1,n2),glat(n1,n2)

integer :: i,j,k
real :: ue,ve

! Rotate winds from rotated polar stereographic grid to RAMS rotated
!    polar stereographic grid


do i=1,n1
   do j=1,n2
      do k=1,n3
         call uvtoueve(aps(i,j,k,1),aps(i,j,k,2),ue,ve  &
                      ,glat(i,j),glon(i,j),ppolelat,ppolelon)
         call uevetouv(aps(i,j,k,1),aps(i,j,k,2),ue,ve  &
                      ,glat(i,j),glon(i,j),polelat,polelon)
      enddo
   enddo
enddo

return
end

