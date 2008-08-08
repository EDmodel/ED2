!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine eng_params

  ! eng_params:
  !  Set some constants that were formerly defined in the namelist.
  !  This set should be changed only for special tests or with code modification.

  use mem_grid, only: &
       sspct,         &
       impl,          &
       iadvl,         &
       iadvf

  use io_params, only: &
       ntopsmth,      &
       izflat
  
  implicit none
  
  sspct    = 0.  ! Sound speed fraction
  impl     = 1   ! Implicit flag for acoustic model  -  0=off, 1=on
  ntopsmth = 0   ! Number of passes for topography smoother
  izflat   = 0   ! Width of flat margin around domain (in grid points)
  iadvl    = 2   ! Order of advection - Leapfrog - 2 or 4 - only 2 now
  iadvf    = 2   ! Order of advection - Forward - 2 or 6 - only 2 now

end subroutine eng_params

!***************************************************************************

subroutine toptinit_user(n2,n3,ifm,topt,topzo)

implicit none

integer :: n2,n3,ifm,nx,ny
real, dimension(n2,n3) :: topt,topzo

integer :: i,j

!  This subroutine is the intended location for a user to customize TOPT,
!  the surface topography array.  It is called after all other types of
!  initialization of this field, so this subroutine has the last word.
!  By default the subroutine makes no change to the field.  The commented
!  lines below serve as a template for user-designed changes; the example
!  shown is the Witch of Agnesi mountain, a common test case.   Note that
!  this routine is called for each grid separately, so attention to the
!  current value of ngrid in this routine may be required by the user.


       return !<<<<<<<<<<<<<<<

!      if (itoptflg(ifm) .lt. 2) return
!
      if(ifm.eq.1)then
      print*,'tttttttttttttttttttttttttttttttttttttttt'
      print*,'Opening topography file: data/topo_AMSUL_wvln11-15-02_dz150m_g1.dat'
      print*,'tttttttttttttttttttttttttttttttttttttttt'
      print*,' '

       open(10,file='data/topo_AMSUL_wvln11-15-02_dz150m_g1.dat',&
            status='old',form='formatted')

        read(10,1998) nx,ny
        if(nx.ne.n2.or.ny.ne.n3) then
	print*,'Grades nao correspondentes nx n2 ny n3:',nx,n2,ny,n3
	print*,'STOP at ruser_traceA'
	stop
	endif
!		
        do j=ny,1,-1
            do i=1,nx
             read(10,2000) topt(i,j)
            enddo
        enddo

      endif

      if(ifm.eq.2)then
      print*,'tttttttttttttttttttttttttttttttttttttttt'
      print*,'Opening topography file: data/topo_AMSUL_wvln11-15-02_dz150m_q2003_g2.dat'
      print*,'tttttttttttttttttttttttttttttttttttttttt'
      print*,' '

       open(10,file='data/topo_AMSUL_wvln11-15-02_dz150m_g2.dat',&
            status='old',form='formatted')

        read(10,1998) nx,ny
        if(nx.ne.n2.or.ny.ne.n3) then
	print*,'Grades nao correspondentes nx n2 ny n3:',nx,n2,ny,n3
	print*,'STOP at ruser_traceA'
	stop
	endif
!		
        do j=ny,1,-1
            do i=1,nx
             read(10,2000) topt(i,j)
            enddo
        enddo

      endif

 1998   format(1x,2i6)
 2000   format(1x,f12.3)	


!
!      if(ifm.eq.1)then
!         hfwid=10000.
!         hgt=10.
!         hfwid2=hfwid**2
!         do j=1,n3
!            do i=1,n2
!               topt(i,j)=hgt*hfwid2/(hfwid2+xtn(i,1)**2)
!               topt(i,j) = 0. + float(i)
!               topzo(i,j) = 0.001 + topt(i,j)
!            enddo
!         enddo
!c         topt(5,1) = 10.
!      elseif(ifm.eq.2)then
!      endif

!if(ifm == 3) then
!   call toposmooth(n2,n3,topt,6,100,nnxp(ifm)-1,2,100)
!endif

return
end

subroutine toposmooth(n2,n3,topt,ntopsmth,i1,i2,j1,j2)

implicit none
integer :: n2,n3,ntopsmth,i1,i2,j1,j2

real, dimension(n2,n3) :: topt
real, allocatable :: vt2da(:,:)

integer :: i,j,iter

allocate (vt2da(n2,n3))
!          smooth the topography if desired.

do iter=1,ntopsmth
   do j=j1,j2
      do i=i1,i2
         vt2da(i,j)=abs(topt(i+1,j)-topt(i,j))
      enddo
   enddo
   print*,'pass:',iter, 'max dtopox:',maxval(vt2da(i1:i2,j1:j2))
   do j=j1,j2
      do i=i1,i2
         vt2da(i,j)=abs(topt(i,j+1)-topt(i,j))
      enddo
   enddo
   print*,'pass:',iter, 'max dtopoy:',maxval(vt2da(i1:i2,j1:j2))
   do j=j1,j2
      do i=i1,i2
         vt2da(i,j)=( (topt(i+1,j)+topt(i-1,j)+topt(i,j+1)+topt(i,j-1))&
                    + 4.*topt(i,j)) / 8.
      enddo
   enddo
   do j=j1,j2
      do i=i1,i2
         topt(i,j)=vt2da(i,j)
      enddo
   enddo
enddo

deallocate (vt2da)

return
end


! ****************************************************************************

subroutine sstinit_user(n2,n3,ifm,seatf)
implicit none
integer :: n2,n3,ifm,i,j
real, dimension(n2,n3) :: seatf

!  This subroutine is the intended location for a user to customize the
!  SEATP and SEATF arrays.  It is called after all other types of
!  initialization of these fields, so this subroutine has the last word.
!  By default the subroutine makes no change to the fields.  The commented
!  lines below serve as a template for user-designed changes.   Note that
!  this routine is called for each grid separately, so attention to the
!  current value of ngrid in this routine may be required by the user.

! if (ifm .eq. 1) then
!    do j = 1,n3
!       do i = 1,n2
!          seatf(i,j) =
!       enddo
!    enddo

! elseif (ifm .eq. 2) then
! endif

return
end

! ****************************************************************************

subroutine ndviinit_user(n2,n3,npat,ifm,veg_ndvif)
implicit none
integer :: n2,n3,npat,ifm,i,j
real, dimension(n2,n3,npat) :: veg_ndvif

!  This subroutine is the intended location for a user to customize the
!  NDVIP and NDVIF arrays.  It is called after all other types of
!  initialization of these fields, so this subroutine has the last word.
!  By default the subroutine makes no change to the fields.  The commented
!  lines below serve as a template for user-designed changes.   Note that
!  this routine is called for each grid separately, so attention to the
!  current value of ngrid in this routine may be required by the user.

! if (ifm .eq. 1) then

!    do j = 1,n3
!       do i = 1,n2
!          veg_ndvif(i,j,1) =
!          veg_ndvif(i,j,2) =

!          do ipat = 3,npat
!             veg(ndvif(i,j,ipat) =
!          enddo

!       enddo
!    enddo

! elseif (ifm .eq. 2) then
! endif

return
end

!*****************************************************************************

subroutine sfcinit_file_user(n2,n3,mzg,npat,ifm  &
   ,patch_area,leaf_class,soil_text)

use rconstants

use catt_start, only: CATT !INTENT(IN)

implicit none

integer :: n2,n3,mzg,npat,ifm,i,j,k,ipat

real, dimension(mzg,n2,n3,npat) :: soil_text
real, dimension(n2,n3,npat) :: patch_area,leaf_class

!  This subroutine is the intended location for a user to customize the
!  PATCH_AREA, leaf_class, and SOIL_TEXT arrays.  It is called after all
!  other types of initialization of these fields, so this subroutine has
!  the last word.  By default the subroutine makes no change to the
!  fields.  The commented lines below serve as a template for user-designed
!  changes.   Note that this routine is called for each grid separately, so
!  attention to the current value of ngrid in this routine may be required
!  by the user.

if (CATT==1) then

   !srf-06-08-2003: altera vegetacao no norte da Africa
   !srf-22-08-2003: altera vegetacao em todo dominio
   do j = 1,n3
      do i = 1,n2
         do ipat = 2,npat
            !print*,'veg=',leaf_class(i,j,ipat)
            if (nint(leaf_class(i,j,ipat)) == 10 ) then  !   semi-deserto
               
               !      if(nveg == 10 ) then                              !   semi-deserto
               ! em 22 fe agosto => rlong divergiu em regioes de semi deserto sobre a Am. do Sul
               !
               !       if(glat(i,j) .ge. 5. .and. glat(i,j) .lt. 35.) then !   5N < lat < 35N
               !        if(glon(i,j) .ge. -15. .and. glon(i,j) .le. 50.) then !  15W < lon < 50E
               
               !          schar(i,j,3,ipatch) = 13.   ! semi-deserto = cerrado
               !print*,'veg=',i,j,leaf_class(i,j,ipat)
               leaf_class(i,j,ipat)= 13.   ! semi-deserto = cerrado
               !	endif
               !       endif
            endif

         enddo
      enddo
   enddo

endif

! if (ifm .eq. 1) then

! do j = 1,n3
!    do i = 1,n2

!       patch_area(i,j,1) =         ! patch 1
!       leaf_class(i,j,1) =         ! patch 1

!       patch_area(i,j,2) =         ! patch 2
!       leaf_class(i,j,2) =         ! patch 2

!       do k = 1,nzg
!          soil_text(k,i,j,1) =     ! patch 1
!          soil_text(k,i,j,2) =     ! patch 2
!       enddo

!    enddo
! enddo

! do ipat = 3,npat
!    do j = 1,n3
!       do i = 1,n2

!          patch_area(i,j,ipat) =
!          leaf_class(i,j,ipat) =

!          do k = 1,nzg
!             soil_text(k,i,j,ipat) =
!          enddo

!       enddo
!    enddo
! enddo

! elseif (ifm .eq. 2) then
! endif

return
end

!*****************************************************************

subroutine sfcinit_nofile_user(n1,n2,n3,mzg,mzs,npat,ifm  &
   ,theta,pi0,pp,rv  &

   ,soil_water     ,soil_energy      ,soil_text       &
   ,sfcwater_mass  ,sfcwater_energy  ,sfcwater_depth  &
   ,ustar          ,tstar            ,rstar           &
   ,veg_fracarea   ,veg_lai          ,veg_tai         &
   ,veg_rough      ,veg_height       ,veg_albedo      &
   ,patch_area                                        &
   ,patch_rough    ,patch_wetind     ,leaf_class      &
   ,soil_rough     ,sfcwater_nlev    ,stom_resist     &
   ,ground_rsat    ,ground_rvap      ,veg_water       &
   ,veg_temp       ,can_rvap         ,can_temp        &
   ,veg_ndvip      ,veg_ndvic        ,veg_ndvif       &
   ,snow_mass      ,snow_depth  &

   ,rvs,prss,pis,vt2da,vt2db,glat,glon,zot,flpw)

use rconstants

implicit none

integer :: n1,n2,n3,mzg,mzs,npat,ifm,i,j,k,ipat,nveg,nsoil

real :: c1,airtemp

real, dimension(n1,n2,n3) :: theta,pi0,pp,rv
real, dimension(n2,n3)    :: rvs,prss,pis,vt2da,vt2db,glat,glon,zot  &
                            ,snow_mass, snow_depth
real, dimension(n2,n3) :: flpw

real, dimension(mzg,n2,n3,npat) :: soil_water,soil_energy,soil_text
real, dimension(mzs,n2,n3,npat) :: sfcwater_mass,sfcwater_energy  &
                                  ,sfcwater_depth

real, dimension(n2,n3,npat) :: ustar        ,tstar         ,rstar        &
                              ,veg_fracarea ,veg_lai       ,veg_tai      &
                              ,veg_rough    ,veg_height    ,veg_albedo   &
                              ,patch_area                                &
                              ,patch_rough  ,patch_wetind  ,leaf_class   &
                              ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                              ,ground_rsat  ,ground_rvap   ,veg_water    &
                              ,veg_temp     ,can_rvap      ,can_temp     &
                              ,veg_ndvip    ,veg_ndvic     ,veg_ndvif

!  This subroutine is the intended location for a user to customize the
!  primary LEAF2 arrays for which standard RAMS data files do not
!  exist.  It is called after all other types of
!  initialization of these fields, so this subroutine has the last word.
!  By default the subroutine makes no change to the
!  fields.  The commented lines below serve as a template for user-designed
!  changes.   Note that this routine is called for each grid separately, so
!  attention to the current value of ngrid in this routine may be required
!  by the user.

!  if (ifm .eq. 1) then

!     call sst_update(n2,n3,mzg,iupdsst    &
!        ,time_seatp,time_seatf,seatp,seatf,soil_energy(1,1,1,1))

!     do j = 1,n3
!        do i = 1,n2

!           patch_rough(i,j,1) =
!           can_temp(i,j,1) =
!           can_rvap(i,j,1) =

!           do ipat = 2,npat

!              nveg = nint(leaf_class(i,j,ipat))

!              soil_rough(i,j,ipat) =
!              patch_rough(i,j,ipat) =
!              veg_rough(i,j,ipat) =

!              veg_height(i,j,ipat) =
!              veg_albedo(i,j,ipat) =
!              sfcwater_nlev(i,j,ipat) =
!              stom_resist(i,j,ipat) =

!              veg_temp(i,j,ipat) =
!              can_temp(i,j,ipat) =

!              veg_water(i,j,ipat) =
!              can_rvap(i,j,ipat) =

!              do k = 1,nzg

!                 nsoil = nint(soil_text(k,i,j,ipat))
!                 soil_water(k,i,j,ipat) = max(soilcp(nsoil),slmstr(k)  &
!                    * slmsts(nsoil))

! For persistent wetlands (bogs, marshes, fens, swamps), initialize with
! saturated soil.  Currently, this corresponds to datq classes 31 and 32.

!                 if (nint(leaf_class(i,j,ipat)) .eq. 31 .or.  &
!                     nint(leaf_class(i,j,ipat)) .eq. 32) then
!                    soil_water(k,i,j,ipat) = slmsts(nsoil)
!                 endif

! By default, initialize soil thermal energy at a temperature equal to
! airtemp + stgoff(k) and with all water assumed to be liquid.  If the
! temperature is initially below 0C, this will immediately adjust to soil
! at 0C with part ice.  In order to begin with partially or totally frozen
! soil, reduce or remove the latent-heat-of-fusion term (the one with the
! factor of 3.34) from soil_energy below.  If the soil is totally frozen and the
! temperature is below zero C, the factor of 4.186 should be changed to 2.093
! to reflect the reduced heat capacity of ice compared to liquid.  These
! changes may be alternatively be done in subroutine sfcinit_user in ruser.f

!                 soil_energy(k,i,j,ipat) = (airtemp - 273.15 + stgoff(k))  &
!                    * (slcpd(nsoil) + soil_water(k,i,j,ipat) * 4.186e6)  &
!                    + soil_water(k,i,j,ipat) * 3.34e8

!              enddo

!              do k = 1,nzs

!                 sfcwater_energy(k,i,j,ipat) =
!                 sfcwater_mass(k,i,j,ipat) =
!                 sfcwater_depth(k,i,j,ipat) =

! For persistent wetlands (bogs, marshes, fens, swamps), initialize with
! 10 cm water depth.  Currently, this corresponds to datq class 31.

!                 if (nint(leaf_class(i,j,ipat)) .eq. 31) then
!                    if (k .eq. 1) then
!                       sfcwater_energy(k,i,j,ipat) =
!                       sfcwater_mass(k,i,j,ipat) =
!                       sfcwater_depth(k,i,j,ipat) =
!                    endif
!                 endif

!                 if (sfcwater_mass(k,i,j,ipat) .gt. 0.)  &
!                    sfcwater_nlev(i,j,ipat) = float(k)

!              enddo
!           enddo
!        enddo
!     enddo

!     do ipat = 2,npat

!   if (ipat >= 2) call vegndvi(n2,n3,ifm   &
!      ,patch_area  (1,1,ipat) ,leaf_class(1,1,ipat)   &
!      ,veg_fracarea(1,1,ipat) ,veg_lai   (1,1,ipat)   &
!      ,veg_tai     (1,1,ipat) ,veg_rough (1,1,ipat)   &
!      ,veg_height  (1,1,ipat) ,veg_albedo(1,1,ipat)   &
!      ,veg_ndvip   (1,1,ipat) ,veg_ndvic (1,1,ipat)   &
!      ,veg_ndvif   (1,1,ipat)                         )

!        call grndvap(n2,n3,nzg,nzs,1,n2,1,n3  &
!           ,soil_energy,soil_water,soil_text,sfcwater_energy,patch_area  &
!           ,sfcwater_nlev,ground_rsat,ground_rvap,can_temp,can_rvap,prss)

!     enddo

!  elseif(ifm .eq. 2) then
!  endif

return
end

!*****************************************************************************

subroutine bubble(m1,m2,m3,thp,rtp)

implicit none
integer :: m1,m2,m3,i,j,k
real, dimension(m1,m2,m3) :: thp,rtp
!
do j = 1,1
   do i = 17,26
     do k = 2,7
!e            thp(k,i,j) = thp(k,i,j) + 5.
!e            rtp(k,i,j) = rtp(k,i,j) * 1.2
     enddo
   enddo
enddo
return
end

!     ******************************************************************

subroutine trtend(m1,m2,m3,k1,k2,i1,i2,j1,j2,rtt,tht,wt,thp,wp,dtlt)
implicit none
integer :: m1,m2,m3,k1,k2,i1,i2,j1,j2,i,j,k

real :: dtlt
real, dimension(m1,m2,m3) :: rtt,tht,wt,thp,wp

do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      rtt(k,i,j)=rtt(k,i,j)+10./3600.
    enddo
  enddo
enddo
do j=j1,j2
  do i=i1,i2
    do k=k1,k2
      tht(k,i,j)=(310.-thp(k,i,j))/dtlt
    enddo
  enddo
enddo
do j=j1,j2
  do i=i1,i2
    do k=k1-1,k2
      wt(k,i,j)=(2.-wp(k,i,j))/dtlt
    enddo
  enddo
enddo
return
end
