!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine sfcdata

use mem_grid
use mem_leaf
use leaf_coms

implicit none

integer :: k,nnn

type soil_class
   real :: slpots     ! Soil moisture potential at saturation                 [       m]
   real :: slmsts     ! Soil moisture at saturation                           [   m3/m3]
   real :: slbs       ! B exponent                                            [     n/d]
   real :: slcpd      ! Specific heat of dry soil                             [  J/m3/K]
   real :: soilcp     ! Dry soil capacity (at -3.1MPa)                        [   m3/m3]
   real :: soilwp     ! Wilting point capacity (at -1.5MPa)                   [   m3/m3]
   real :: slcons     ! hydraulic conductivity at saturation                  [     m/s]
   real :: slcons0    ! ?Dry soil density (also total soil porosity)          [   kg/m3]
   real :: soilcond0
   real :: soilcond1
   real :: soilcond2
   real :: sfldcap
   real :: xsand
   real :: xclay
   real :: xorgan
   real :: xrobulk
   real :: slden
end type soil_class

real :: romin,roorg,slfcap,refdepth,tmin,ratio,xmin
real, dimension     (nstyp) :: xsand,xclay,xorgan,xrobulk
type(soil_class) , dimension(nstyp)            :: soilparms

real, dimension (12,0:nvtyp+nvtyp_teb) :: bioparms

!  Soil Characteristics (see Clapp & Hornberger, 1978; McCumber & Pielke,
!                        1981; Pielke, 1984; Tremback & Kessler, 1985)
!
!  slpots  - saturation moisture potential (m)
!  slmsts  - saturation volumetric moisture content (m3/m3)
!  slbs    - b exponent (dimensionless)
!  slcons  - saturation soil hydraulic conductivity (m/s)
!  slcons0 - surface value for slcons (m/s)
!  slcpd   - dry soil volumetric heat capacity (J/m3/K)
!  slden   - dry soil density (kg/m3) (also total soil porosity)





!         LEAF-3 BIOPHYSICAL PARAMETERS BY LANDUSE CLASS NUMBER

data bioparms/  &
!-----------------------------------------------------------------------------
!albv_green     sr_max         veg_clump       rootdep             LEAF-3 CLASS #
!     albv_brown     tai_max        veg_frac        dead_frac      AND DESCRIPTION
!          emisv          sai            veg_ht         rcmin
!-----------------------------------------------------------------------------
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  0  Ocean
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  1  Lakes, rivers, streams
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  2  Ice cap/glacier
 .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  3  Desert, bare soil
 .14, .24, .97, 5.4, 8.0, 1.0, 1.0, .80, 20.0, 1.5, .0, 500., & !  4  Evergreen needleleaf tree
 .14, .24, .95, 5.4, 8.0, 1.0, 1.0, .80, 22.0, 1.5, .0, 500., & !  5  Deciduous needleleaf tree
 .20, .24, .95, 6.2, 7.0, 1.0,  .0, .80, 22.0, 1.5, .0, 500., & !  6  Deciduous broadleaf tree
 .12, .18, .95, 4.1, 6.5, 1.0,  .0, .90, 32.0, 2.5, .0, 285., & !  7  Evergreen broadleaf tree
 .13, .30, .96, 5.1, 4.0, 1.0,  .0, .75,   .3,  .7, .7, 100., & !  8  Short grass
 .24, .43, .96, 5.1, 5.0, 1.0,  .0, .80,  1.2, 1.0, .7, 100., & !  9  Tall grass
 .24, .24, .96, 5.1, 1.0,  .2, 1.0, .20,   .7, 1.0, .0, 500., & ! 10  Semi-desert
 .20, .24, .95, 5.1, 4.5,  .5, 1.0, .60,   .2, 1.0, .0,  50., & ! 11  Tundra
 .14, .24, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 12  Evergreen shrub
 .20, .28, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 13  Deciduous shrub
 .16, .24, .96, 6.2, 7.0, 1.0,  .5, .80, 22.0, 1.5, .0, 500., & ! 14  Mixed woodland
 .22, .40, .95, 5.1, 5.0,  .5,  .0, .85,  1.0, 1.0, .0, 100., & ! 15  Crop/mixed farming, C3 grassland
 .18, .40, .95, 5.1, 5.0,  .5,  .0, .80,  1.1, 1.0, .0, 500., & ! 16  Irrigated crop
 .12, .43, .98, 5.1, 7.0, 1.0,  .0, .80,  1.6, 1.0, .0, 500., & ! 17  Bog or marsh
 !srf  ---
 .13, .30, .96, 5.1, 6.0, 1.0,  .0, .80,  7.0, 1.5, .0, 100., & ! 18  Wooded grassland
 .20, .36, .90, 5.1, 3.6, 1.0,  .0, .74,  6.0,  .8, .0, 500., & ! 19  Urban and built up
 .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 1.5, .0, 500., & ! 20  Wetland evergreen broadleaf tree
 .16, .24, .96, 5.1, 2.0, 1.5, 1.0, .10, 20.0, 1.5, .0, 500./   ! 21  Very urban



! Soil constants
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
!                     slpots     slmsts       slbs     slcpd         soilcp        soilwp        slcons       slcons0 soilcond0 soilcond1 soilcond2       sfldcap    xsand    xclay   xorgan xrobulk    slden !
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
soilparms= (/ &
     soil_class(-0.049831046,  0.373250,  3.295000,  1465000.,  0.026183447,  0.032636854,  2.446420e-5,  0.000500000,   0.3000,   4.8000,  -2.7000,  0.132130936,  0.9200,  0.0300,  0.0500,   1200.,  1600.)  & ! 1. Sand
  ,  soil_class(-0.068643592,  0.386340,  3.796000,  1407000.,  0.041873866,  0.050698650,  1.751180e-5,  0.000600000,   0.3000,   4.6600,  -2.6000,  0.155720881,  0.8200,  0.0600,  0.1200,   1250.,  1600.)  & ! 2. Loamy Sand
  ,  soil_class(-0.155095787,  0.418940,  4.496000,  1344000.,  0.076932647,  0.090413463,  8.228700e-6,  0.000769000,   0.2900,   4.2700,  -2.3100,  0.199963447,  0.5800,  0.1000,  0.3200,   1300.,  1600.)  & ! 3. Sandy loam
  ,  soil_class(-0.659933234,  0.476050,  5.090000,  1273000.,  0.141599910,  0.163306007,  2.396240e-6,  0.000010600,   0.2700,   3.4700,  -1.7400,  0.266720155,  0.1700,  0.1300,  0.7000,   1400.,  1600.)  & ! 4. Silt loam
  ,  soil_class(-0.238341682,  0.437280,  5.797000,  1214000.,  0.126501190,  0.143377074,  4.732940e-6,  0.002200000,   0.2800,   3.6300,  -1.8500,  0.247334654,  0.4300,  0.1800,  0.3900,   1350.,  1600.)  & ! 5. Loam
  ,  soil_class(-0.121199269,  0.412650,  7.165000,  1177000.,  0.137648733,  0.152325872,  6.405180e-6,  0.001500000,   0.2800,   3.7800,  -1.9600,  0.250954745,  0.5800,  0.2700,  0.1500,   1350.,  1600.)  & ! 6. Sandy clay loam
  ,  soil_class(-0.627769194,  0.478220,  8.408000,  1319000.,  0.228171947,  0.248747504,  1.435260e-6,  0.000107000,   0.2600,   2.7300,  -1.2000,  0.333825332,  0.1000,  0.3400,  0.5600,   1500.,  1600.)  & ! 7. Silty clay loam
  ,  soil_class(-0.281968114,  0.446980,  8.342000,  1227000.,  0.192624431,  0.210137962,  2.717260e-6,  0.002200000,   0.2700,   3.2300,  -1.5600,  0.301335491,  0.3200,  0.3400,  0.3400,   1450.,  1600.)  & ! 8. Clay loam
  ,  soil_class(-0.121283019,  0.415620,  9.538000,  1177000.,  0.182198910,  0.196607427,  4.314510e-6,  0.000002167,   0.2700,   3.3200,  -1.6300,  0.286363001,  0.5200,  0.4200,  0.0600,   1450.,  1600.)  & ! 9. Sandy clay
  ,  soil_class(-0.601312179,  0.479090,  10.46100,  1151000.,  0.263228486,  0.282143846,  1.055190e-6,  0.000001033,   0.2500,   2.5800,  -1.0900,  0.360319788,  0.0600,  0.4700,  0.4700,   1650.,  1600.)  & !10. Silty clay
  ,  soil_class(-0.286417797,  0.452300,  12.14000,  1088000.,  0.253968998,  0.269618857,  1.427350e-6,  0.000001283,   0.2500,   2.4000,  -0.9600,  0.348432364,  0.2200,  0.5800,  0.2000,   1700.,  1600.)  & !11. Clay
  ,  soil_class(-0.534564359,  0.469200,  6.180000,  8740000.,  0.167047523,  0.187868805,  2.357930e-6,  0.000008000,   0.0600,   0.4600,   0.0000,  0.285709966,  0.2000,  0.2000,  0.6000,    500.,   300.)  & !12. Peat
            /)
! Thermal conductivity in J/msK
cka = 0.418684 * 0.0615
ckw = 0.418684 * 1.45
romin = 2655.0
roorg = 1300.0

slz(nzg+1) = 0.

slfcap = -10. / 3.
refdepth = -2.0

do nnn = 1,nstyp
   slcons0(nnn) = soilparms(nnn)%slcons0
   fhydraul(nnn) = log (soilparms(nnn)%slcons / soilparms(nnn)%slcons0) / refdepth

   do k = 1,nzg
      slcons1(k,nnn) = soilparms(nnn)%slcons      ! ORIGINAL form - const with depth
!     slcons1(k,nnn) = soilparms(nnn)%slcons0  &  ! TOPMODEL form - large at surface
!        * exp(slz(k) * fhydraul(nnn))            !    and exp decrease with depth
   end do

   slpots   (nnn) = soilparms(nnn)%slpots
   slmsts   (nnn) = soilparms(nnn)%slmsts
   slbs     (nnn) = soilparms(nnn)%slbs
   slcons   (nnn) = soilparms(nnn)%slcons
   slcons00 (nnn) = soilparms(nnn)%slcons0
   slcpd    (nnn) = soilparms(nnn)%slcpd
   slden    (nnn) = soilparms(nnn)%slden
   sfldcap  (nnn) = soilparms(nnn)%sfldcap
   soilcp   (nnn) = soilparms(nnn)%soilcp
   soilwp   (nnn) = soilparms(nnn)%soilwp
   soilcond0(nnn) = soilparms(nnn)%soilcond0
   soilcond1(nnn) = soilparms(nnn)%soilcond1
   soilcond2(nnn) = soilparms(nnn)%soilcond2
   xsand    (nnn) = soilparms(nnn)%xsand
   xclay    (nnn) = soilparms(nnn)%xclay
   xorgan   (nnn) = soilparms(nnn)%xorgan
   xrobulk  (nnn) = soilparms(nnn)%xrobulk


   emisg(nnn) = .98
   slfc(nnn) = slmsts(nnn) * (slfcap / slpots(nnn)) ** (-1. / slbs(nnn))

!    tmin = xsand(nnn) + xclay(nnn)
!    romean = xorgan(nnn) * roorg + tmin * romin
!    sporo(nnn) = 1.0 - xrobulk(nnn) / romean
!    ratio = (xorgan(nnn) / roorg) / (tmin / romin)
!    sorgan(nnn) = ratio / (1.0 + ratio) * (1.0 - sporo(nnn))
!    xmin = 1.0 - sporo(nnn) - sorgan(nnn)
!    ssand(nnn) = xsand(nnn) * xmin / tmin
!    sclay(nnn) = xclay(nnn) * xmin / tmin

enddo

do nnn = 1,(nvtyp+nvtyp_teb)
   albv_green(nnn) = bioparms(1,nnn)
   albv_brown(nnn) = bioparms(2,nnn)
   emisv(nnn)      = bioparms(3,nnn)
   sr_max(nnn)     = bioparms(4,nnn)
   tai_max(nnn)    = bioparms(5,nnn)
   sai(nnn)        = bioparms(6,nnn)
   veg_clump(nnn)  = bioparms(7,nnn)
   veg_frac(nnn)   = bioparms(8,nnn)
   veg_ht(nnn)     = bioparms(9,nnn)
   dead_frac(nnn)  = bioparms(11,nnn)
   rcmin(nnn)      = bioparms(12,nnn)
   glai_max(nnn)   = tai_max(nnn) - sai(nnn)

   root(1,nnn)  = 0.              ! not used
   kroot(nnn)   = nzg
   do k = nzg-1,1,-1
      if (slz(k+1) .gt. -bioparms(10,nnn)) kroot(nnn) = k
   enddo
enddo

return
end

! ****************************************************************************

subroutine snowinit(n2,n3,snow_mass,snow_depth)
implicit none
integer :: n2,n3,i,j
real, dimension(n2,n3) :: snow_mass,snow_depth

! Fill the snow_mass and snow_depth arrays with a default value of 0.  This
! default is used when snowcover information is not read from varfiles.

do j = 1,n3
   do i = 1,n2
      snow_mass(i,j) = 0.
      snow_depth(i,j) = 0.
   enddo
enddo
return
end


!*****************************************************************************

subroutine sfcinit_file(n2,n3,mzg,npat,ifm,patch_area,leaf_class,soil_text)

use mem_leaf
use rconstants

implicit none

integer :: n2,n3,mzg,npat,ifm,i,j,k,ipat

real, dimension(mzg,n2,n3,npat) :: soil_text
real, dimension(n2,n3,npat) :: patch_area,leaf_class

! This subroutine fills the arrays PATCH_AREA, leaf_class, and SOIL_TEXT
! horizontally homogeneously with default values that are defined in the
! RAMSIN namelist file.  These fields comprise the land/sea surface data
! types that are normally available on standard RAMS datasets.  The default
! values assigned here may be overridden by (1) interpolation from coarser
! grids, (2) specifying new values in subroutine sfcinit_user in the file
! ruser.f, or (3) reading data from the standard RAMS datasets.

do j = 1,n3
   do i = 1,n2

      patch_area(i,j,1) = 1. - pctlcon         ! patch 1
      leaf_class(i,j,1) = 0.                   ! patch 1

      patch_area(i,j,2) = pctlcon              ! patch 2
      leaf_class(i,j,2) = float(nvgcon)        ! patch 2

      do k = 1,mzg
         soil_text(k,i,j,1) = 0.               ! patch 1
         soil_text(k,i,j,2) = float(nslcon)    ! patch 2
      enddo

   enddo
enddo

do ipat = 3,npat
   do j = 1,n3
      do i = 1,n2

         patch_area(i,j,ipat) = 0.
         leaf_class(i,j,ipat) = leaf_class(i,j,2)

         do k = 1,mzg
            soil_text(k,i,j,ipat) = soil_text(k,i,j,2)
         enddo

      enddo
   enddo
enddo

return
end subroutine sfcinit_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine fills the primary LEAF3 arrays for which standard RAMS data files   !
! do not exist with initial default values.  Many of the initial values are horizontally   !
! homogeneous, although some depend on atmospheric conditions.  The default values         !
! assigned here may be overridden by:                                                      !
! 1. Specification from coarser grids;                                                     !
! 2. Specifying new values in subroutine sfcinit_nofile_user in the file ruser.f90.        !
!------------------------------------------------------------------------------------------!
subroutine sfcinit_nofile(n1,n2,n3,mzg,mzs,npat,ifm,theta,pi0,pp,rv,co2p,seatp,seatf       &
                         ,soil_water,soil_energy,soil_text,sfcwater_mass,sfcwater_energy   &
                         ,sfcwater_depth,ustar,tstar,rstar,cstar,veg_fracarea,veg_agb      &
                         ,veg_lai,veg_tai,veg_rough,veg_height,veg_albedo,patch_area       &
                         ,patch_rough,patch_wetind,leaf_class,soil_rough,sfcwater_nlev     &
                         ,stom_resist,ground_rsat,ground_rvap,ground_temp,ground_fliq      &
                         ,veg_water,veg_hcap,veg_energy,can_prss,can_theta,can_rvap        &
                         ,can_co2,sensible,evap,transp,gpp,plresp,resphet,veg_ndvip        &
                         ,veg_ndvic,veg_ndvif,snow_mass,snow_depth,rvv,prsv,piv,vt2da      &
                         ,vt2db,glat,glon,zot,flpw,rtgt)
   use mem_grid
   use mem_leaf
   use leaf_coms
   use io_params
   use rconstants
   use therm_lib , only : reducedpress

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)    :: n1,n2,n3,mzg,mzs,npat,ifm
   real, dimension(n1,n2,n3)      , intent(in)    :: theta,pi0,pp,rv,co2p
   real, dimension(   n2,n3)      , intent(in)    :: seatp,seatf,snow_mass,snow_depth
   real, dimension(   n2,n3)      , intent(in)    :: flpw,rtgt
   real, dimension(mzg,n2,n3,npat), intent(inout) :: soil_water,soil_energy,soil_text
   real, dimension(mzs,n2,n3,npat), intent(inout) :: sfcwater_mass,sfcwater_energy
   real, dimension(mzs,n2,n3,npat), intent(inout) :: sfcwater_depth
   real, dimension(    n2,n3,npat), intent(inout) :: ustar,tstar,rstar,cstar
   real, dimension(    n2,n3,npat), intent(inout) :: veg_fracarea,veg_agb,veg_lai,veg_tai
   real, dimension(    n2,n3,npat), intent(inout) :: veg_rough,veg_height,veg_albedo
   real, dimension(    n2,n3,npat), intent(inout) :: patch_area,patch_rough,patch_wetind
   real, dimension(    n2,n3,npat), intent(inout) :: leaf_class
   real, dimension(    n2,n3,npat), intent(inout) :: soil_rough,sfcwater_nlev,stom_resist
   real, dimension(    n2,n3,npat), intent(inout) :: ground_rsat,ground_rvap
   real, dimension(    n2,n3,npat), intent(inout) :: ground_temp,ground_fliq
   real, dimension(    n2,n3,npat), intent(inout) :: veg_water,veg_energy,veg_hcap
   real, dimension(    n2,n3,npat), intent(inout) :: can_prss,can_theta,can_rvap,can_co2
   real, dimension(    n2,n3,npat), intent(inout) :: sensible,evap,transp
   real, dimension(    n2,n3,npat), intent(inout) :: gpp,plresp,resphet
   real, dimension(    n2,n3,npat), intent(inout) :: veg_ndvip,veg_ndvic,veg_ndvif
   real, dimension(   n2,n3)      , intent(inout) :: rvv,prsv,piv,vt2da,vt2db,glat,glon,zot
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: k2
   integer                                        :: i,j,k,ipat,nveg,nsoil
   real                                           :: tsoil
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Time interpolation factor for updating SST.  Use calculation of continuous model !
   ! time.                                                                                 !
   !---------------------------------------------------------------------------------------!
   if (iupdsst == 0) then
      timefac_sst = 0.
   else
      timefac_sst = (time - ssttime1(ifm)) / (ssttime2(ifm) - ssttime1(ifm))
   end if

   jloop: do j = 1,n3
      iloop: do i = 1,n2
         k2=nint(flpw(i,j))
         piv(i,j)  = 0.5 * cpi * (pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
         prsv(i,j) = piv(i,j) ** cpor * p00
         geoht     = (zt(k2)-zm(k2-1)) * rtgt(i,j)

         atm_shv   = rv(k2,i,j) / (rv(k2,i,j) + 1.)

         patch_rough(i,j,1) = waterrough

         !---------------------------------------------------------------------------------!
         !     Canopy properties.  Copy conserved variables from lowest atmospheric grid,  !
         ! and compute pressure and temperature.                                           !
         !---------------------------------------------------------------------------------!
         can_prss(i,j,1) = reducedpress(prsv(i,j),theta(k2,i,j),atm_shv,geoht              &
                                       ,theta(k2,i,j),atm_shv,can_depth)
         can_theta(i,j,1)   = theta(k2,i,j)
         can_rvap(i,j,1)    = rv(k2,i,j)
         can_co2(i,j,1)     = co2p(k2,i,j)
         can_temp           = theta(k2,i,j) * (p00i * can_prss(i,j,1)) ** rocp

         !----- Water patch, so we set vegetation properties to zero. ---------------------!
         veg_energy(i,j,1)  = 0.0
         veg_water (i,j,1)  = 0.0
         veg_hcap  (i,j,1)  = 0.0
         
         !----- Unless we are solving water lilies or phytoplankton, these should be 0. ---!
         veg_lai   (i,j,1)  = 0.0
         veg_tai   (i,j,1)  = 0.0
         veg_agb   (i,j,1)  = 0.0

         !----- Soil properties. Except for top layer energy, everything is set to zero. --!
         soil_energy(:,i,j,1) = 0.
         soil_water (:,i,j,1) = 1.
         soil_energy(mzg,i,j,1) =  cliq * (seatp(i,j) + (seatf(i,j) - seatp(i,j))          &
                                                      * timefac_sst  - tsupercool)

         !----- Fluxes.  Initially they should be all zero. -------------------------------!
         sensible (i,j,1) = 0.0
         evap     (i,j,1) = 0.0
         transp   (i,j,1) = 0.0
         gpp      (i,j,1) = 0.0
         plresp   (i,j,1) = 0.0
         resphet  (i,j,1) = 0.0

         !----- Above-ground biomass.  This should be always 0 for water patches. ---------!
         veg_agb  (i,j,1) = 0.0

         !---------------------------------------------------------------------------------!
         !     We now loop over the land patches.  Some properties such as canopy          !
         ! properties are simply copied from the water patch.                              !
         !---------------------------------------------------------------------------------!
         patchloop: do ipat = 2,npat

            nveg = nint(leaf_class(i,j,ipat))

            soil_rough(i,j,ipat)  = zrough
            patch_rough(i,j,ipat) = max(zrough,zot(i,j))
            veg_rough(i,j,ipat)   = .13 * veg_ht(nveg)

            veg_height(i,j,ipat)  = veg_ht(nveg)
            veg_albedo(i,j,ipat)  = albv_green(nveg)
            stom_resist(i,j,ipat) = 1.e6

            !------------------------------------------------------------------------------!
            !     We cannot allow the vegetation height to be above the first level.  If   !
            ! that happens, we stop the run and ask the user to coarsen delta-z.           !
            !------------------------------------------------------------------------------!
            can_depth = max(veg_height(i,j,ipat),can_depth_min)
            if (can_depth + 6. > geoht) then
               write (unit=*,fmt='(a)') '================================================='
               write (unit=*,fmt='(a)') '   DELTA-Z is too fine, the first level is'
               write (unit=*,fmt='(a)') ' beneath or too close to the top of the canopy...'
               write (unit=*,fmt='(a)') ' Try coarsening it and run again.'
               write (unit=*,fmt='(a)') ' '
               write (unit=*,fmt='(a,1x,i5)'    ) 'I          =',i
               write (unit=*,fmt='(a,1x,i5)'    ) 'J          =',j
               write (unit=*,fmt='(a,1x,i5)'    ) 'P          =',ipat
               write (unit=*,fmt='(a,1x,i5)'    ) 'NVEG       =',nveg
               write (unit=*,fmt='(a,1x,es12.5)') 'GLON       =',glon(i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'GLAT       =',glat(i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'GEOHT      =',geoht
               write (unit=*,fmt='(a,1x,es12.5)') 'VEG_HEIGHT =',veg_height(i,j,ipat)
               write (unit=*,fmt='(a,1x,es12.5)') 'CAN_DEPTH  =',can_depth
               write (unit=*,fmt='(a)') '================================================='
               call abort_run('Delta-z is too thin!','sfcinit_nofile','leaf3_init.f90')
            end if


            veg_hcap  (i,j,ipat) = hcapveg_ref * max(veg_ht(nveg),hcapveg_hmin)
            veg_water (i,j,ipat) = 0.
            veg_energy(i,j,ipat) = veg_hcap(i,j,ipat) * can_temp

            !----- Above-ground biomass.  This is non-0 only when we run ED-2. ------------!
            veg_agb   (i,j,ipat) = 0.

            !----- Canopy air properties.  Initially we just assume same properties. ------!
            can_prss  (i,j,ipat) = can_prss (i,j,1)
            can_theta (i,j,ipat) = can_theta(i,j,1)
            can_rvap  (i,j,ipat) = can_rvap (i,j,1)
            can_co2   (i,j,ipat) = can_co2  (i,j,1)

            !----- Fluxes.  Initially they should be all zero. ----------------------------!
            sensible (i,j,ipat) = 0.0
            evap     (i,j,ipat) = 0.0
            transp   (i,j,ipat) = 0.0
            gpp      (i,j,ipat) = 0.0
            plresp   (i,j,ipat) = 0.0
            resphet  (i,j,ipat) = 0.0

            do k = 1,mzg

               nsoil = nint(soil_text(k,i,j,ipat))

               !---------------------------------------------------------------------------!
               !     For persistent wetlands (bogs, marshes, fens, swamps) and irrigated   !
               ! crops, initialize with saturated soil.  Currently, this corresponds to    !
               ! leaf classes 16, 17, and 20.  Otherwise, use the user-defined input       !
               ! profile.                                                                  !
               !---------------------------------------------------------------------------!
               select case (nint(leaf_class(i,j,ipat)))
               case (16,17,20)
                  soil_water(k,i,j,ipat) = slmsts(nsoil)
               case default
                  soil_water(k,i,j,ipat) = max(soilcp(nsoil),slmstr(k) * slmsts(nsoil))
               end select

               !---------------------------------------------------------------------------!
               !     By default, initialize soil internal energy at a temperature equal to !
               ! can_temp + stgoff(k).  If the temperature is initially below triple       !
               ! point, we assume all soil water to be frozen, otherwise we assume all     !
               ! water to be liquid.                                                       !
               !---------------------------------------------------------------------------!
               tsoil = can_temp + stgoff(k)
               if (can_temp >= t3ple) then
                  soil_energy(k,i,j,ipat) = slcpd(nsoil) * tsoil + soil_water(k,i,j,ipat)  &
                                          * cliqvlme * (tsoil - tsupercool)
               else
                  soil_energy(k,i,j,ipat) = slcpd(nsoil) * tsoil                           &
                                          + soil_water(k,i,j,ipat) * cicevlme * tsoil
               end if
            end do

            !------ Surface water, if any, will initially occupy just the first level. ----!
            do k = 1,mzs
               sfcwater_mass(k,i,j,ipat) = 0.
               sfcwater_energy(k,i,j,ipat) = 0.
               sfcwater_depth(k,i,j,ipat) = 0.
            end do

            !------------------------------------------------------------------------------!
            !    For persistent wetlands (bogs, marshes, fens, swamps), initialise with 10 !
            ! cm water depth and temperature not colder than the triple point.  Currently, !
            ! this corresponds to leaf classes 17 and  20.  Glaciers (leaf class 2) will   !
            ! be initialised with a thick layer of compact ice (currently 6 m), not warmer !
            ! than the triple point, so we make sure there will enough ice there.  The ice !
            ! will be eventually split into several layers.                                !
            !------------------------------------------------------------------------------!
            select case (nint(leaf_class(i,j,ipat)))
            case (2)
               sfcwater_depth (1,i,j,ipat) = 6.
               sfcwater_mass  (1,i,j,ipat) = idns * sfcwater_depth(1,i,j,ipat)
               sfcwater_energy(1,i,j,ipat) = cice * min(t3ple,can_temp)
            case (17,20)
               sfcwater_depth (1,i,j,ipat) = .1
               sfcwater_mass  (1,i,j,ipat) = wdns * sfcwater_depth(1,i,j,ipat) 
               sfcwater_energy(1,i,j,ipat) = cliq * (max(t3ple,can_temp) -tsupercool)
            end select

            !------------------------------------------------------------------------------!
            !    If there is initial snow information, add the information into the first  !
            ! layer.                                                                       !
            !------------------------------------------------------------------------------!
            if (snow_mass(i,j) > 0.) then
               sfcwater_energy(1,i,j,ipat) = ( sfcwater_energy(1,i,j,ipat)                 &
                                             * sfcwater_mass  (1,i,j,ipat)                 &
                                             + min(t3ple,can_temp)*cice*snow_mass(i,j) )   &
                                           / (sfcwater_mass(1,i,j,ipat) + snow_mass(i,j))
               sfcwater_mass  (1,i,j,ipat) = sfcwater_mass (1,i,j,ipat) + snow_mass(i,j)
               !---------------------------------------------------------------------------!
               !    Add a depth of 5x the equivalent of liquid depth so the extra snow is  !
               ! fresh and fluffy.                                                         !
               !---------------------------------------------------------------------------!
               sfcwater_depth(1,i,j,ipat)  = sfcwater_depth(1,i,j,ipat)                    &
                                           + snow_mass(i,j) * 5. * wdnsi
            end if

            !----- We must have either nothing or a single layer at this point. -----------!
            if (sfcwater_mass(1,i,j,ipat) > 0.) then
               sfcwater_nlev(i,j,ipat) = 1.
            else
               sfcwater_nlev(i,j,ipat) = 0.
            end if


            call vegndvi(ifm,patch_area  (i,j,ipat) ,leaf_class(i,j,ipat)                  &
                            ,veg_fracarea(i,j,ipat) ,veg_lai   (i,j,ipat)                  &
                            ,veg_tai     (i,j,ipat) ,veg_rough (i,j,ipat)                  &
                            ,veg_height  (i,j,ipat) ,veg_albedo(i,j,ipat)                  &
                            ,veg_ndvip   (i,j,ipat) ,veg_ndvic (i,j,ipat)                  &
                            ,veg_ndvif   (i,j,ipat)                                        )

            call leaf_grndvap( soil_energy(mzg,i,j,ipat) , soil_water     (mzg,i,j,ipat)   &
                             , soil_text  (mzg,i,j,ipat) , sfcwater_energy(mzs,i,j,ipat)   &
                             , sfcwater_nlev  (i,j,ipat) , can_rvap           (i,j,ipat)   &
                             , can_prss       (i,j,ipat) , ground_rsat        (i,j,ipat)   &
                             , ground_rvap    (i,j,ipat) , ground_temp        (i,j,ipat)   &
                             , ground_fliq    (i,j,ipat) )
         end do patchloop
      end do iloop
   end do jloop

   return
end subroutine sfcinit_nofile

!**********************************************************************

subroutine datp_datq(datp, datq)

  ! This subroutine maps the input datp classes to a smaller set datq
  ! which represents the full set of LEAF-2 or LEAF-3 classes for which
  ! LSP values are
  ! defined.

  ! For TEB_SPM
  use teb_spm_start, only : TEB_SPM ! INTENT(in)

  implicit none
  ! Arguments:
  real, intent(in)     :: datp
  integer, intent(out) :: datq
  ! Local variables:
!!$  integer              :: catb(0:94) 
  integer              :: catb(0:95) 
  !integer              :: catb_leaf3(0:94) ! not used


  !  Olson Global Ecosystems dataset (94 classes) mapped to LEAF-3 classes
  !  (see leaf3_document).

!!$  !-------------------------------------------!
!!$  data catb/ 0,                             & !
!!$       19, 8, 4, 5, 6, 7, 9, 3,11,16,  & !  0
!!$!!!    10, 2,17, 1, 0,12,13,14,18, 4,  & ! 10 (use for future Olson data)
!!$       !--srf--
!!$       10, 2,20, 0, 0,12,13,14,18, 4,  & ! 10 (use for current Olson data)
!!$       04, 4,14,14, 6, 6, 4, 7, 7,15,  & ! 20
!!$       15, 6, 7, 7,15,16,16,16,16, 8,  & ! 30
!!$       08, 8,18,17,17,12,12, 7,10, 3,  & ! 40
!!$       10,10,11,14,18,18,18,18,13, 6,  & ! 50
!!$       05, 4,11,12, 0, 0, 0, 0, 3, 2,  & ! 60
!!$       03,20, 0,17,17,17, 4,14, 7, 3,  & ! 70
!!$       03, 3, 3, 3, 3, 3, 8,12, 7, 6,  & ! 80
!!$       18,15,15,15                  /    ! 90
!!$  !-------------------------------------------!
!!$  !     1  2  3  4  5  6  7  8  9 10

  !-------------------------------------------!
  data catb/ 0,                             & !
       19, 8, 4, 5, 6, 7, 9, 3,11,16,  & !  0
!!!    10, 2,17, 1, 0,12,13,14,18, 4,  & ! 10 (use for future Olson data)
       !--srf--
       10, 2,20, 0, 0,12,13,14,18, 4,  & ! 10 (use for current Olson data)
       04, 4,14,14, 6, 6, 4, 7, 7,15,  & ! 20
       15, 6, 7, 7,15,16,16,16,16, 8,  & ! 30
       08, 8,18,17,17,12,12, 7,10, 3,  & ! 40
       10,10,11,14,18,18,18,18,13, 6,  & ! 50
       05, 4,11,12, 0, 0, 0, 0, 3, 2,  & ! 60
       03,20, 0,17,17,17, 4,14, 7, 3,  & ! 70
       03, 3, 3, 3, 3, 3, 8,12, 7, 6,  & ! 80
       18,15,15,15,19               /    ! 90
  !-------------------------------------------!
  !     1  2  3  4  5  6  7  8  9 10
  

  ! IF TEB_SPM:
  ! Using the Global Ecosystem legend 95 maped to RAMS type 21 (Very Urban Type)
  ! Otherwise the OGE legend 95 is maped to RAMS type 19 (Urban)
  if (TEB_SPM==1) catb(95) = 21

  ! Temporario
  if (TEB_SPM==1) catb(50) = 21

  !SRF 2008: OGE-INPE data set has values greater than 95
  if (nint(datp) >= 0. .and. nint(datp) <= 95. ) then
     datq = catb(nint(datp))
  else
     datq = catb(0)
  end if

  return
end subroutine datp_datq

!**********************************************************************

subroutine datp_datsoil(datp,datsoil)

! This subroutine maps the input datp soil classes to a smaller set datsoil
! which represents the full set of LEAF-2 classes for which soil parameter
! values are defined.

implicit none
integer datsoil,catb(0:133)
real datp

! (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our
! own missing value designated 133, to the USDA soil textural classes.  FAO
! classes [0] (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not
! used in original FAO dataset), [132] (water), and [133] (our missing value)
! are all mapped to a default class of sandy clay loam in case they happen to
! correspond to a land surface area in the landuse dataset that RAMS uses to
! define land area.  We wrote missing value class 133 to the RAMS FAO files
! whenever a negative value (which is out of range of defined values) was
! found in the original FAO dataset, which occurred in about 2.6% of the
! pixels.  For the remaining FAO classes, a cross reference table to Zobler
! soil texture classes that was provided, plus our own cross referencing table
! from Zobler to USDA classes listed below, provides the mapping from FAO to
! USDA.  In this mapping, we use only organic USDA classes and omit nonorganic
! classes since most soils do contain organic matter, and organic content
! information is not provided in the Zobler classes.

!  Zobler Class              USDA Class

!  1  Coarse                 2  Loamy sand
!  2  Medium                 4  Silt loam
!  3  Fine                   8  Clay loam
!  4  Coarse-medium          3  Sandy loam
!  5  Coarse-fine            6  Sandy clay loam
!  6  Medium-fine            7  Silty clay loam
!  7  Coarse-medium-fine     6  Sandy clay loam
!  8  Organic matter         5  Loam

!                            1  Sand (not used)
!                            9  Sandy clay (not used)
!                           10  Silty clay (not used)
!                           11  Clay (not used)
!                           12  Peat (not used)

data catb/ 6  &
         , 6, 4, 4, 7, 7, 8, 6, 4, 4, 4  &
         , 7, 4, 4, 4, 8, 4, 8, 4, 4, 8  &
         , 4, 2, 4, 4, 4, 4, 6, 8, 8, 8  &
         , 4, 8, 8, 2, 6, 4, 7, 4, 4, 3  &
         , 4, 6, 7, 4, 4, 4, 4, 4, 4, 4  &
         , 4, 4, 4, 4, 4, 4, 2, 4, 4, 2  &
         , 4, 3, 4, 2, 7, 6, 4, 4, 6, 8  &
         , 8, 7, 2, 5, 4, 5, 6, 6, 4, 2  &
         , 2, 2, 4, 6, 2, 2, 2, 2, 2, 4  &
         , 2, 2, 2, 4, 2, 4, 3, 6, 2, 7  &
         , 4, 4, 4, 8, 8, 8, 3, 7, 4, 4  &
         , 4, 3, 6, 4, 2, 4, 4, 4, 2, 2  &
         , 2, 4, 6, 4, 4, 7, 7, 6, 3, 2  &
         , 2, 6, 6 /

datsoil = catb(nint(datp))

return
end
