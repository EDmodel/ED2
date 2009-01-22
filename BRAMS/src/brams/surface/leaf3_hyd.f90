!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

!
!     ******************************************************************
!
!  npatch = maximum number of patches in any grid cell
!  nzgp = nzg + 1
!  maxgrp = maximum allowable number of groups in domain
!  ngrps = number of groups in the domain as read from the data file
!  soil_water is the 4d soil moisture array with values in m3/m3
!  patch_area is a 3d array that contains patch fractional area of a grid cell,
!  soil_text is a 4d array that contains soil textural class
!  slz defines the depth in meters below the surface of the bottom of each
!     soil level
!  ig is the i index (x-direction) of each group
!  jg is the j index (y-direction) of each group
!  ipg is an array of (model) patch numbers in each group
!  lpg is the number of patches within each group
!  slmsts is an array of soil moisture capacities, based on soil type.
!     should this be redefined?
!  zi is a scratch array for computing local (patch) water table heights
!  fa as read in here is the fractional area of a patch relative to a grid cell
!  wi as read in here is the wetness index of a patch
!  finv is the depth (m) over which saturated hydraulic conductivity decreases
!     by a factor of e.

subroutine hydro(m2,m3,mzg,mzs,np                           &
   ,soil_water,soil_energy,soil_text                        &
   ,sfcwater_mass,sfcwater_energy,patch_area,patch_wetind)

use mem_grid
use mem_leaf
use leaf_coms

implicit none
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real, dimension(nzgmax,nstyp) :: slcons1
real, dimension(nstyp) :: slcons0,fhydraul
common/efold/slcons1,slcons0,fhydraul
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer, parameter :: maxpatch=15,maxgrp=200
integer :: m2,m3,mzg,mzs,np,ngrp,ng,l,i,j,ip
integer, save  :: icall4=0,ngrps
integer, save, dimension(maxgrp) :: ig,jg
integer, save, dimension(maxpatch) :: lpg
integer, save, dimension(maxpatch,maxgrp) :: ipg

real :: rhow,rhowi
real, dimension(maxpatch) :: zi,wateradd
real, save, dimension(maxpatch) :: fa,wi,finv

real, dimension(mzg,m2,m3,np) :: soil_energy,soil_water,soil_text
real, dimension(mzs,m2,m3,np) :: sfcwater_energy,sfcwater_mass
real, dimension(m2,m3,np)     :: patch_area,patch_wetind

real, dimension(7,10) :: profile

data rhow,rhowi/1.e3,1.e-3/


if (icall4 .ne. 4) then
   icall4 = 4
!  read in patch table values of ig, jg, lpg, and ipg
   open(91,file='patch_table',status='old',form='formatted')
   read(91,10) ngrps
   do ngrp = 1,ngrps

      read(91,10) ng,ig(ng),jg(ng),lpg(ng),(ipg(l,ng),l=1,lpg(ng))
      read(91,11) (fa(l),l=1,lpg(ng))
      read(91,12) finv(ng),(wi(l),l=1,lpg(ng))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! replace finv with value from subroutine sfcdat (use type 6 = sandy clay loam)
      finv(ng) = 1. / fhydraul(6)
      print*, 'ng,finv',ng,finv(ng)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

10         format(24i5)
11         format(20x,20f5.2)
12         format(15x,f5.3,20f5.1)

!  define patch fractional area (of a grid cell) and wetness index

      i = ig(ng)
      j = jg(ng)
      do l = 1,lpg(ng)
         ip = ipg(l,ng)
         patch_area(i,j,ip) = fa(l)
         patch_wetind(i,j,ip) = wi(l)
      enddo

   enddo
   close(91)
endif

if (ngrps .lt. 1) return

call hydrol(m2,m3,mzg,mzs,np,maxpatch,ngrps                &
   ,soil_energy,soil_water,soil_text                       &
   ,sfcwater_energy,sfcwater_mass,patch_area,patch_wetind  &
   ,slz,ig,jg,ipg,lpg,slmsts,zi,fa,wi,finv                 &
   ,rhow,rhowi,dtlt,slcons0,slcpd,wateradd,time)

return
end

!*****************************************************************************

subroutine hydrol(m2,m3,mzg,mzs,np,maxpatch,ngrps          &
   ,soil_energy,soil_water,soil_text                       &
   ,sfcwater_energy,sfcwater_mass,patch_area,patch_wetind  &
   ,slz,ig,jg,ipg,lpg,slmsts,zi,fa,wi,finv                 &
   ,rhow,rhowi,dtlt,slcons0,slcpd,wateradd,time)

use rconstants, only: cliq,cliqvlme, alli, wdns,wdnsi, tsupercool
use therm_lib, only: qwtk,qtk
implicit none
integer :: m2,m3,mzg,mzs,np,maxpatch,ngrps,ngd,i,j,k,lp,l,ip,nsoil,ibotpatch
integer, dimension(*) :: ig,jg,lpg
integer, dimension(maxpatch,*) :: ipg
integer, dimension(20) :: ksat

real, dimension(20) :: zitopm

real :: zibar,wibar,okbar,etafbar,tgfa,energysum,watersum,slcpdsum  &
       ,fracw,vol,tempktopm,fracliq,olflow,qolflow,runoff,qrunoff   &
       ,dtltopm,baseflow,ziadd,wsum,qwsum,delta_water,delta_energy  &
       ,tempk,q,capacity,add,qw,wfreeb,qwfree,rhow,rhowi,dtlt
real(kind=8) :: time

real, dimension(mzg,m2,m3,np) :: soil_energy,soil_water,soil_text
real, dimension(mzs,m2,m3,np) :: sfcwater_energy,sfcwater_mass
real, dimension(m2,m3,np)     :: patch_area,patch_wetind

real, dimension(*) :: slz,slmsts,zi,fa,wi,finv,slcons0,slcpd,wateradd

!  tgfa is the total group fractional area of the grid cell, not including
!     the group's bottomland patch
!  fa as computed here is the patch's fractional area of a group
!  wi is the patch wetness index, copied from the patch_wetind array
!  l is a consecutive counter over the patches in a group up to lpg(mzg),
!     the number in the group

do ngd = 1,ngrps

   i = ig(ngd)
   j = jg(ngd)
   lp = lpg(ngd)
   zibar = 0.
   wibar = 0.
   okbar = 0.
   etafbar = 0.
   tgfa = 0.
   energysum = 0.
   watersum = 0.
   slcpdsum = 1.e-3
   do l = 2,lp
      ip = ipg(l,ngd)
      tgfa = tgfa + patch_area(i,j,ip)
   enddo

   do l = 2,lp
      ip = ipg(l,ngd)
      zi(l) = 0.
!
!  Compute zi values starting at the bottom of the lowest soil layer and
!  summing over all saturated levels, defined here as levels that have more
!  than 95% of the full moisture capacity), but accounting for any deficit
!  in those levels.  Sum the total soil energy (soil_energy * dz * dA), water
!  content (soil_water * dz * dA), and soil heat capacity (slcpd * dz * dA)
!  in order to compute a mean temperature of transported water.

      ksat(l) = 1
      zi(l) = slz(1)
      do k = 1,mzg
         nsoil = nint(soil_text(k,i,j,ip))
         fracw = soil_water(k,i,j,ip) / slmsts(nsoil)

         if (fracw .lt. .95) go to 57

         zi(l) = zi(l) + fracw * (slz(k+1) - slz(k))
         ksat(l) = k

         vol = (slz(k+1) - slz(k)) * patch_area(i,j,ip)
         energysum = energysum + vol * soil_energy(k,i,j,ip)
         watersum = watersum + vol * soil_water(k,i,j,ip)
         slcpdsum = slcpdsum + vol * slcpd(nsoil)
      enddo
57         continue

      fa(l) = patch_area(i,j,ip) / tgfa
      wi(l) = patch_wetind(i,j,ip)

      zibar = zibar + zi(l) * fa(l)
      wibar = wibar + wi(l) * fa(l)
!
!  compute mean hydraulic conductivity K_0 and mean soil porosity etaf
!  for the group
!
      okbar = okbar + slcons0(nsoil) * fa(l)
      etafbar = etafbar + slmsts(nsoil) * fa(l)

   enddo

   call qwtk(energysum,watersum*wdns,slcpdsum,tempktopm,fracliq)

!  If there is no saturated water or if saturated soil is more than half
!  frozen, skip over soil hydrology.

   olflow = 0.
   qolflow = 0.
   runoff = 0.
   qrunoff = 0.
   ibotpatch = ipg(1,ngd)

   if (watersum .lt. 1.e-6 .or. fracliq .lt. .5) go to 60

   dtltopm = etafbar / (okbar*exp(-wibar) * exp(zibar/finv(ngd)))
   baseflow = finv(ngd) / dtltopm

   do l = 2,lp
      ip = ipg(l,ngd)

      zitopm(l) = zibar + finv(ngd) * (wi(l) - wibar)
      ziadd = dtlt * ((zitopm(l) - zi(l)) / dtltopm )
      wateradd(l) = .45 * ziadd - baseflow * dtlt

! Wateradd is the actual (positive or negative) column height (m) of pure
! water to add to each column.  It is assumed to be a constant .45 times
! ziadd in order to ensure total water conservation.

   enddo

! First, subtract water from columns (patches) that need to dry.  Compute
! internal energy (qwsum) of water that gets subtracted.

   wsum = 0.
   qwsum = 0.
   do l = 2,lp
      if (wateradd(l) .lt. 0.) then
         ip = ipg(l,ngd)
         nsoil = nint(soil_text(ksat(l),i,j,ip))
         call qwtk(soil_energy(ksat(l),i,j,ip)  &
                  ,soil_water (ksat(l),i,j,ip)*wdns  &
                  ,slcpd(nsoil),tempk,fracliq)
         delta_water = wateradd(l) / (slz(ksat(l)+1) - slz(ksat(l)))
         soil_water(ksat(l),i,j,ip) = soil_water(ksat(l),i,j,ip) + delta_water
         delta_energy = delta_water * cliqvlme * (tempk - tsupercool)
         soil_energy(ksat(l),i,j,ip) = soil_energy(ksat(l),i,j,ip) + delta_energy
         wsum = wsum + wateradd(l) * fa(l)
         qwsum = qwsum + wateradd(l) * cliqvlme * (tempk - tsupercool) * fa(l)
      endif
   enddo

! Now, add excess water to columns (patches) that need to moisten

   q = qwsum / min(-1.e-15,wsum)

   do l = 2,lp
      if (wateradd(l) .gt. 0.) then
         ip = ipg(l,ngd)
         nsoil = nint(soil_text(ksat(l),i,j,ip))

         do k = 1,mzg
            nsoil = nint(soil_text(k,i,j,ip))
            capacity = (slmsts(nsoil) - soil_water(k,i,j,ip))  &
               * (slz(k+1) - slz(k))
            add = min(wateradd(l),capacity)
            delta_water = add / (slz(k+1) - slz(k))
            soil_water(k,i,j,ip) = soil_water(k,i,j,ip) + delta_water
            soil_energy(k,i,j,ip) = soil_energy(k,i,j,ip) + delta_water * q
            wateradd(l) = wateradd(l) - add

         enddo
         if (wateradd(l) .gt. 0) then
            olflow = olflow + wateradd(l) * patch_area(i,j,ip)
         endif

      endif

   enddo

!  Add baseflow to bottomland patch

   wateradd(1) = baseflow * dtlt * tgfa / patch_area(i,j,ibotpatch)

   do k = 1,mzg
      capacity = (slmsts(nsoil) - soil_water(k,i,j,ibotpatch))  &
         * (slz(k+1) - slz(k))
      add = min(wateradd(1),capacity)
      delta_water = add / (slz(k+1) - slz(k))
      soil_water(k,i,j,ibotpatch) = soil_water(k,i,j,ibotpatch) + delta_water
      soil_energy(k,i,j,ibotpatch) = soil_energy(k,i,j,ibotpatch)  &
         + delta_water * q
      wateradd(1) = wateradd(1) - add
   enddo

   if (wateradd(1) .gt. 0) then
      qw = sfcwater_energy(1,i,j,ibotpatch)  &
         * sfcwater_mass(1,i,j,ibotpatch) + wateradd(1) * q
      sfcwater_mass(1,i,j,ibotpatch) = sfcwater_mass(1,i,j,ibotpatch)  &
         + wateradd(1) * 1.e3
      sfcwater_energy(1,i,j,ibotpatch) = qw  &
         / max(1.e-4,sfcwater_mass(1,i,j,ibotpatch))
   endif

60      continue

!  runoff from each patch - accumulate in 'runoff' and 'qrunoff' arrays

   do l = 2,lp
      ip = ipg(l,ngd)

      if (sfcwater_mass(1,i,j,ip) .gt. 0.) then
         call qtk(sfcwater_energy(1,i,j,ip),tempk,fracliq)

         if (fracliq .gt. .1) then
            qw = sfcwater_energy(1,i,j,ip) * sfcwater_mass(1,i,j,ip)
            wfreeb = sfcwater_mass(1,i,j,ip) * (fracliq - .1) / 0.9
            qwfree = wfreeb * cliq * (tempk - tsupercool)
            sfcwater_mass(1,i,j,ip) = sfcwater_mass(1,i,j,ip) - wfreeb
            sfcwater_energy(1,i,j,ip) = (qw - qwfree)  &
               / (max(1.e-4,sfcwater_mass(1,i,j,ip)))
            sfcwater_energy(1,i,j,ip) = max (0., min (6.4e5,  &
               sfcwater_energy(1,i,j,ip)))

            runoff = runoff + wfreeb * patch_area(i,j,ip)
            qrunoff = qrunoff + qwfree * patch_area(i,j,ip)

         endif
      endif

   enddo

! Add overland flow and runoff to bottomland patch.  First convert each
! to kg/m2 for bottomland patch.

   olflow = olflow * wdns / patch_area(i,j,ibotpatch)
   runoff = runoff / patch_area(i,j,ibotpatch)
   qolflow = olflow * q * wdnsi
   qrunoff = qrunoff / patch_area(i,j,ibotpatch)

   qw = sfcwater_energy(1,i,j,ibotpatch) * sfcwater_mass(1,i,j,ibotpatch)  &
      + qolflow + qrunoff
   sfcwater_mass(1,i,j,ibotpatch) = sfcwater_mass(1,i,j,ibotpatch)  &
      + olflow + runoff
   sfcwater_energy(1,i,j,ibotpatch) = qw  &
      / max(1.e-4,sfcwater_mass(1,i,j,ibotpatch))

enddo

return
end


