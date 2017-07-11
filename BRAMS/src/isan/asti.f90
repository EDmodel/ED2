!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

SUBROUTINE ISNSTAGE ()

use isan_coms
use mem_grid

implicit none

character csuff*3
character(len=8) :: rot_type

real gxval(maxlev),gyval(maxlev)
real :: reflat1,reflon1,olat1,olat2,olon1,olon2  &
       ,olatmin,olatmax,olonmin,olonmax
integer :: ng

! ---------------------------------------------start  ngrid =1 section
if(ngrid.eq.1) then

!-----------------------------------------------------------------
!      Read rawinsondes and surface observations
!-----------------------------------------------------------------
!      Default to RAMS grid 1 for obs access and plotting
!-----------------------------------------------------------------

   call RAMS_mm (grid_g(ngrid)%glat,nnxp(1)*nnyp(1),olatmin,olatmax)
   call RAMS_mm (grid_g(ngrid)%glon,nnxp(1)*nnyp(1),olonmin,olonmax)

   olat1=olatmin-5.
   olat2=olatmax+5.
   olon1=olonmin-5.
   olon2=olonmax+5.

   nsta=0

   if(iproc_flag(natime,3).eq.1) then

      if(.not.allocated(up_p))   allocate(up_p(maxsta,maxlev))
      if(.not.allocated(up_t))   allocate(up_t(maxsta,maxlev))
      if(.not.allocated(up_z))   allocate(up_z(maxsta,maxlev))
      if(.not.allocated(up_r))   allocate(up_r(maxsta,maxlev))
      if(.not.allocated(up_uz))  allocate(up_uz(maxsta,maxlev))
      if(.not.allocated(up_vz))  allocate(up_vz(maxsta,maxlev))
      if(.not.allocated(up_ur))  allocate(up_ur(maxsta,maxlev))
      if(.not.allocated(up_vr))  allocate(up_vr(maxsta,maxlev))
      if(.not.allocated(up_zz))  allocate(up_zz(maxsta,maxlev))

      if(.not.allocated(up_lat)) allocate(up_lat(maxsta))
      if(.not.allocated(up_lon)) allocate(up_lon(maxsta))
      if(.not.allocated(up_top)) allocate(up_top(maxsta))
      if(.not.allocated(up_lp))  allocate(up_lp(maxsta))
      if(.not.allocated(up_lz))  allocate(up_lz(maxsta))
      if(.not.allocated(up_topg))allocate(up_topg(maxsta,nigrids))

      if(.not.allocated(up_chstid)) allocate (up_chstid(maxsta))

      inrawi=iproc_names(natime,2)

      call input_rawi (olat1,olat2,olon1,olon2)
           
      do ng=1,nigrids

         if(ng.gt.1) up_topg(1:maxsta,ng)=up_topg(1:maxsta,nxtnest(ng)) 

         call soundtopo (ng,nsta,maxsta,up_topg(:,ng)  &
                        ,up_lat,up_lon,up_top  &
                        ,grid_g(ng)%topt,grid_g(ng)%glat  &
                        ,grid_g(ng)%glon &
                        ,nnxp(ng),nnyp(ng),platn(ng),plonn(ng)  &
                        ,xtn(1,ng),ytn(1,ng),deltaxn(ng),deltayn(ng))
      enddo


   endif

!-----------------------------------------------------------------
!     Read surface observations
!-----------------------------------------------------------------

   nssfc=0

   if(iproc_flag(natime,4).eq.1) then

      
      if(.not.allocated(sf_u)) allocate(sf_u(maxsfc))
      if(.not.allocated(sf_v)) allocate(sf_v(maxsfc))
      if(.not.allocated(sf_ur)) allocate(sf_ur(maxsfc))
      if(.not.allocated(sf_vr)) allocate(sf_vr(maxsfc))
      if(.not.allocated(sf_p)) allocate(sf_p(maxsfc))
      if(.not.allocated(sf_t)) allocate(sf_t(maxsfc))
      if(.not.allocated(sf_s)) allocate(sf_s(maxsfc))
      if(.not.allocated(sf_r)) allocate(sf_r(maxsfc))
      if(.not.allocated(sf_lat)) allocate(sf_lat(maxsfc))
      if(.not.allocated(sf_lon)) allocate(sf_lon(maxsfc))
      if(.not.allocated(sf_top)) allocate(sf_top(maxsfc))
      if(.not.allocated(sf_scra)) allocate(sf_scra(maxsfc))
      
      if(.not.allocated(sf_chstid)) allocate (sf_chstid(maxsfc))
      if(.not.allocated(sf_date)) allocate (sf_date(maxsfc))

      insrfce=iproc_names(natime,3)

      call input_sfc (olat1,olat2,olon1,olon2)
  
   endif

endif


! end of ngrid =1 section

!-----------------------------------------------------------------
!    Interpolate pressure data to RAMS grid
!-----------------------------------------------------------------

if(igridfl.gt.0) then

   if ( guess1st.eq.'PRESS' .or. guess1st.eq.'NC') then

!         First, interpolate first guess pressure data horizontally to the 
!                RAMS polar-stereo grid
!         --------------------------------------------------------

      print*,'Allocating polar/pressure grid-'  &
           ,nnxp(ngrid),nnyp(ngrid),nprz,ngrid

      if(.not.allocated(pp_u)) allocate(pp_u(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_v)) allocate(pp_v(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_t)) allocate(pp_t(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_z)) allocate(pp_z(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_r)) allocate(pp_r(nnxp(ngrid),nnyp(ngrid),nprz))
      if(.not.allocated(pp_sglob)) allocate(pp_sglob(nprx+3,npry+2))
      
!         Rotate winds to the RAMS polar-stereo grid
!         --------------------------------------------------------

      if(inproj.eq.1) then
         rot_type='ll_rps'
         reflat1=0.
         reflon1=0.
      elseif(inproj.eq.2) then
         rot_type='lc_rps'
         reflat1=cntlat
         reflon1=cntlon
      elseif(inproj.eq.3) then
         rot_type='ps_rps'
         reflat1=cntlat
         reflon1=cntlon
      endif
         
      call rotate_winds(rot_type,nprx,npry,nprz  &
                       ,p_u,p_v,p_ur,p_vr  &
                       ,p_lat,p_lon  &
                       ,reflat1,reflon1,platn(ngrid),plonn(ngrid))

      if(inproj.eq.1) then
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_u,pp_sglob,p_ur,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_v,pp_sglob,p_vr,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_t,pp_sglob,p_t,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_z,pp_sglob,p_z,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,1  &
             ,pp_r,pp_sglob,p_r,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
              ! Surface fields 
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_slp,pp_sglob,p_slp,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sfp,pp_sglob,p_sfp,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sft,pp_sglob,p_sft,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snow,pp_sglob,p_snow,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         call latlon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sst,pp_sglob,p_sst,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,idatelin,iglobew,iglobs,iglobn)
         rot_type='ll_rps'
         reflat1=0.
         reflon1=0.
      elseif(inproj.eq.2) then
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_u,p_ur,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_v,p_vr,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_t,p_t,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_z,p_z,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,1  &
             ,pp_r,p_r,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
              ! Surface fields 
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_slp,p_slp,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sfp,p_sfp,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sft,p_sft,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snow,p_snow,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call lambcon_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sst,p_sst,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         rot_type='lc_rps'
         reflat1=cntlat
         reflon1=cntlon
      elseif(inproj.eq.3) then
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_u,p_ur,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_v,p_vr,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_t,p_t,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,0  &
             ,pp_z,p_z,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,nprz,1,1  &
             ,pp_r,p_r,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
              ! Surface fields 
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_slp,p_slp,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sfp,p_sfp,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sft,p_sft,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_snow,p_snow,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         call trueps60_ps (nnxp(ngrid),nnyp(ngrid),nprx,npry,1,1,1  &
             ,rs_sst,p_sst,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
             ,xswlat,xswlon,gdatdx,gdatdy,cntlat,cntlon)
         rot_type='tps_rps'
         reflat1=cntlat
         reflon1=cntlon
      endif
      
!         Next, interpolate vertically to the specified isentropes
!         --------------------------------------------------------

      call vterpp_i (nnxp(ngrid),nnyp(ngrid),nprz,nisn  &
                    ,pp_u,pp_v,pp_t,pp_z,pp_r  &
                    ,pi_u,pi_v,pi_p,pi_s,pi_r)

!         Last, interpolate vertically to the sigma-z levels
!         --------------------------------------------------------

      call vterpp_s (nnxp(ngrid),nnyp(ngrid),nprz,nsigz  &
                    ,pp_u,pp_v,pp_t,pp_z,pp_r  &
                    ,ps_u,ps_v,ps_p,ps_t,ps_r  &
                    ,grid_g(ngrid)%topt  &
                    ,grid_g(ngrid)%rtgt)

      if(allocated(pp_u)) deallocate(pp_u)
      if(allocated(pp_v)) deallocate(pp_v)
      if(allocated(pp_t)) deallocate(pp_t)
      if(allocated(pp_z)) deallocate(pp_z)
      if(allocated(pp_r)) deallocate(pp_r)
      if(allocated(pp_sglob)) deallocate(pp_sglob)

   elseif(guess1st.eq.'RAMS') then

      call first_RAMS(nnxp(ngrid),nnyp(ngrid),nnzp(ngrid) &
                     ,ps_u,ps_v,ps_p,ps_t,ps_r)
 
      !  print*,'RAMS first guess not working for A array removal'
      !  stop 'no RAMS first guess'
   endif

endif

!-----------------------------------------------------------------
!     Do objective analysis or interpolation to RAMS grids
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!     If only using gridded p data, fill rs_qual, then we're done...
!-----------------------------------------------------------------

if(igridfl.eq.3) then

   call strmfun (nnxp(ngrid),nnyp(ngrid),grid_g(ngrid)%topt  &
                    ,grid_g(ngrid)%rtgt)

   rs_qual(1:nnxp(ngrid),1:nnyp(ngrid)) = 0.

   goto 1200
endif

!-----------------------------------------------------------------
!           Perform upper air objective analysis/interpolation
!-----------------------------------------------------------------

if(guess1st.ne.'RAMS') then

!    Rotate obs winds to RAMS polar-stereo grid  
!        (reflat1,reflon1) doesn't matter since obs winds 
!         are earth relative
   
   if(iproc_flag(natime,3).eq.1) then
      call rotate_winds('ll_rps',maxsta,1,maxlev  &
                       ,up_uz,up_vz,up_ur,up_vr  &
                       ,up_lat,up_lon  &
                       ,reflat1,reflon1,platn(ngrid),plonn(ngrid) )
   endif

   if(iproc_flag(natime,4).eq.1) then
      call rotate_winds('ll_rps',nssfc,1,1 &
                       ,sf_u,sf_v,sf_ur,sf_vr  &
                       ,sf_lat,sf_lon  &
                       ,reflat1,reflon1,platn(ngrid),plonn(ngrid) )
   endif

!                 Isentropic coordinates
!-----------------------------------------------------------------

!         Interpolate rawindsondes vertically to isentropes

   print*,'Allocating obs/isen array-',nsta,nisn
   
   allocate(upi_u(nsta,nisn))
   allocate(upi_v(nsta,nisn))
   allocate(upi_p(nsta,nisn))
   allocate(upi_s(nsta,nisn))
   allocate(upi_r(nsta,nisn))

   if (nsta > 0) call obs_isen (maxsta,maxlev  &
                 ,up_t,up_p,up_z,up_r,up_ur,up_vr,up_zz,up_lat,up_lon  &
                 ,up_top,up_lp,up_lz  &
                 ,upi_u,upi_v,upi_p,upi_s,upi_r,nsta,up_chstid)

   call obj_anal ('isen',ngrid,nnxp(ngrid),nnyp(ngrid)  &
                 ,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
                 ,platn(ngrid),plonn(ngrid),xtn(1,ngrid),ytn(1,ngrid)  &
                 ,deltaxn(ngrid),deltayn(ngrid))
              
   deallocate (upi_u,upi_v,upi_p,upi_s,upi_r)
endif

!                 Sigma-z coordinates
!-----------------------------------------------------------------

!         Interpolate rawindsondes vertically to sigma-z

print*,'Allocating obs/sigma-z array-',nsta,nsigz

allocate(ups_u(nsta,nsigz))
allocate(ups_v(nsta,nsigz))
allocate(ups_p(nsta,nsigz))
allocate(ups_t(nsta,nsigz))
allocate(ups_r(nsta,nsigz))
   
if (nsta > 0)  call obs_sigz (maxsta,maxlev  &
              ,up_t,up_p,up_z,up_r,up_ur,up_vr,up_zz,up_lat,up_lon  &
              ,up_top,up_lp,up_lz,up_topg(1,ngrid)  &
              ,ups_u,ups_v,ups_p,ups_t,ups_r,nsta,zmn(nnzp(1)-1,1))

call obj_anal ('sigz',ngrid,nnxp(ngrid),nnyp(ngrid)  &
              ,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
              ,platn(ngrid),plonn(ngrid),xtn(1,ngrid),ytn(1,ngrid)  &
              ,deltaxn(ngrid),deltayn(ngrid))

deallocate (ups_u,ups_v,ups_p,ups_t,ups_r)

!-----------------------------------------------------------------
!          Find Montgomery streamfunction on isentropic surfaces
!-----------------------------------------------------------------

call strmfun(nnxp(ngrid),nnyp(ngrid),grid_g(ngrid)%topt  &
            ,grid_g(ngrid)%rtgt)

!-----------------------------------------------------------------
!          Perform ground surface objective analysis.
!-----------------------------------------------------------------

call obj_anal ('surf',ngrid,nnxp(ngrid),nnyp(ngrid) &
           ,grid_g(ngrid)%glat,grid_g(ngrid)%glon  &
           ,platn(ngrid),plonn(ngrid),xtn(1,ngrid),ytn(1,ngrid)  &
           ,deltaxn(ngrid),deltayn(ngrid))


! Check if "quality" of surface grid point is good

call sfcqual (nnxp(ngrid),nnyp(ngrid),grid_g(ngrid)%glat  &
             ,grid_g(ngrid)%glon  &
             ,rs_qual,sf_lat,sf_lon,sf_p,nssfc,sf_scra)

1200 continue

   
return
end

!***************************************************************************

subroutine sfcqual (nxp,nyp,glat,glon,quals3  &
                   ,sflt,sfln,pss,nssfc,scra)

use rconstants

implicit none

integer :: nxp,nyp,nssfc
real ::  glat(nxp,nyp),glon(nxp,nyp),sflt(*),sfln(*),scra(*)  &
        ,quals3(nxp,nyp),pss(*)

integer :: iqf(4),i,j,ns
real :: gobkm,gobrd

gobkm=2.*spconkm
gobrd=4.*spconkm

do i=1,nxp
   do j=1,nyp
      quals3(i,j)=0.
      iqf(1)=0
      iqf(2)=0
      iqf(3)=0
      iqf(4)=0
      do ns=1,nssfc
         scra(ns)=sqrt(((glat(i,j)-sflt(ns))*spconkm)**2  &
                 +((glon(i,j)-sfln(ns))*spconkm  &
                 *cos((glat(i,j)+sflt(ns))*.5*.01745))**2)
      enddo
      do ns=1,nssfc
         if(pss(ns).lt.1e19.and.scra(ns).lt.gobkm) then
            quals3(i,j)=1.
         endif
      enddo
      do ns=1,nssfc
         if(scra(ns).le.gobrd) then
            if(sfln(ns).le.glon(i,j).and.sflt(ns).ge.glat(i,j)) iqf(1)=1
            if(sfln(ns).ge.glon(i,j).and.sflt(ns).ge.glat(i,j)) iqf(2)=1
            if(sfln(ns).ge.glon(i,j).and.sflt(ns).le.glat(i,j)) iqf(3)=1
            if(sfln(ns).le.glon(i,j).and.sflt(ns).le.glat(i,j)) iqf(4)=1
         endif
      enddo
      if(iqf(1)+iqf(2)+iqf(3)+iqf(4).ge.3) quals3(i,j)=quals3(i,j)+2.
   enddo
enddo

return
end

!***************************************************************************

subroutine strmfun (nxp,nyp,topt,rtgt)

use isan_coms
use rconstants
use therm_lib, only : ptrh2rvapil,virtt
implicit none

integer :: nxp,nyp
real, dimension(nxp,nyp) :: topt,rtgt

real, dimension(maxsigz) :: sigzr,temp,thv

integer :: k,i,j,lbchyd,lbc,lbcp
real :: syo,po,tho,thvo,sigo,bcpr,raux

if(guess1st.ne.'RAMS') then
 
   ! Find a boundary condition as first level at or below 360K
 
   lbchyd=360
   do k=nisn,1,-1
      if(levth(k).le.lbchyd) then
         do i=1,nxp
            do j=1,nyp
               if(pi_p(i,j,k).gt.1e19) goto 50
            enddo
         enddo
         lbc=k
         print 54,lbc,levth(lbc)
         54 format(///,' isentropic hydrostatic boundary set at level',2i6,' k')
         goto 52
         50 continue
      endif
   enddo
   print 53,i,j,k
   53 format(' could not find good isentropic boundary level ',3i6)
   stop 'st3-under'
   52 continue
 
   do i=1,nxp
      do j=1,nyp
 
        syo=pi_s(i,j,lbc)
        po=pi_p(i,j,lbc)
        tho=levth(lbc)
        do k=lbc+1,nisn
           pi_s(i,j,k)=1e30
           if(pi_p(i,j,k).lt.1e19) then
               pi_s(i,j,k)=syo+cpdry*(po**rocp+pi_p(i,j,k)**rocp)  &
                    *.5/p00**rocp *(levth(k)-tho)
               syo=pi_s(i,j,k)
               po=pi_p(i,j,k)
               tho=levth(k)
            endif
!            if(i==1.and.j==1) then
!               print*,'1loop:',i,j,k,pi_s(i,j,k),pi_p(i,j,k),levth(k)
!            endif
         enddo
 
         syo=pi_s(i,j,lbc)
         po=pi_p(i,j,lbc)
         tho=levth(lbc)
         do k=lbc-1,1,-1
            pi_s(i,j,k)=1e30
            if(pi_p(i,j,k).lt.1e19) then
               pi_s(i,j,k)=syo+cpdry*(po**rocp+pi_p(i,j,k)**rocp)  &
                    *.5/p00**rocp*(levth(k)-tho)
               syo=pi_s(i,j,k)
               po=pi_p(i,j,k)
               tho=levth(k)
            endif
!            if(i==1.and.j==1) then
!               print*,'2loop:',i,j,k,pi_s(i,j,k),pi_p(i,j,k),levth(k)
!            endif
         enddo
 
      enddo
   enddo

endif

!stop'me'
bcpr=10000.
do k=nsigz,1,-1
   if(sigz(k).le.bcpr) then
      do i=1,nxp
         do j=1,nyp
            if(ps_t(i,j,k).gt.1e19) goto 60
         enddo
      enddo
      lbcp=k
      print 64,lbcp,sigz(lbcp)
      64 format(///,' Sigma-z hydrostatic boundary set at level',i6,f10.1,' m')
      goto 62
      60 continue
   endif
enddo
print 63
63 format(' Could not find good sigma-z boundary level ')
stop 'st3-unders'
62 continue

do i=1,nxp
   do j=1,nyp

      do k=1,nsigz
         temp(k)=ps_t(i,j,k)*(ps_p(i,j,k)/p00)**rocp
         sigzr(k)=topt(i,j)+sigz(k)*rtgt(i,j)
         raux = ptrh2rvapil(ps_r(i,j,k),ps_p(i,j,k),temp(k),.false.)
         thv(k)=virtt(ps_t(i,j,k),raux)
      enddo

      ps_p(i,j,lbcp)=cpdry*(ps_p(i,j,lbcp)/p00)**rocp

      thvo=thv(lbcp)
      po=ps_p(i,j,lbcp)
      sigo=sigzr(lbcp)
      do k=lbcp+1,nsigz
         ps_p(i,j,k)=1e30
         if(ps_t(i,j,k).lt.1e19.and.ps_r(i,j,k).lt.1.e19) then
            ps_p(i,j,k)=po-grav*(sigzr(k)-sigo)/((thvo+thv(k))*.5)
            thvo=thv(k)
            po=ps_p(i,j,k)
            sigo=sigzr(k)
         endif
      enddo

      thvo=thv(lbcp)
      po=ps_p(i,j,lbcp)
      sigo=sigzr(lbcp)
      do k=lbcp-1,1,-1
         ps_p(i,j,k)=1e30
         if(ps_t(i,j,k).lt.1e19.and.ps_r(i,j,k).lt.1.e19) then
            ps_p(i,j,k)=po+grav*(sigo-sigzr(k))/((thvo+thv(k))*.5)
            thvo=thv(k)
            po=ps_p(i,j,k)
            sigo=sigzr(k)
         endif
      enddo

      do k=1,nsigz
         ps_p(i,j,k)=(ps_p(i,j,k)/cpdry)**cpor*p00
      enddo

   enddo
enddo

return
end

!***************************************************************************

subroutine soundtopo (ngr,nst,m1,sndtopg,stlt,stln,topsnd,topt  &
                     ,glat,glon,n1,n2,polat,polon,swx,swy,delx,dely)
                     
implicit none                     

integer ::     ngr,nst,m1 ,n1,n2
real :: stlt(m1),stln(m1),sndtopg(m1),topsnd(m1)  &
     ,topt(n1,n2),glat(n1,n2),glon(n1,n2),sdat(2)
real :: polat,polon,swx,swy,delx,dely

integer :: ns
real :: topo,stx(1),sty(1)

do ns=1,nst
   sdat(1)=0.
   call ll_xy(stlt(ns),stln(ns),polat,polon,stx(1),sty(1))
   call stainterp(topt,n1,n2,1,1,1,sdat,stx,sty  &
                 ,1,polat,polon,swx,swy,delx,dely,-1.e30)
   topo=-sdat(2)
   if(ngr.eq.1) then
      if(topo.lt.1e20) then
         sndtopg(ns)=topo
      else
         sndtopg(ns)=topsnd(ns)
      endif
   else
      if(topo.lt.1e20) sndtopg(ns)=topo
   endif
enddo

return
end


