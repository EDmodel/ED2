!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine sfcdata_sib_driver

  use mem_grid
  use mem_leaf
  use leaf_coms

  implicit none
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  real, dimension(nzgmax,nstyp) :: slcons1
  real, dimension(nstyp) :: slcons0,fhydraul
  common/efold/slcons1,slcons0,fhydraul
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  integer :: k,nnn

  real :: romin,roorg,slfcap,refdepth,tmin,ratio,xmin
  real, dimension     (nstyp) :: xsand,xclay,xorgan,xrobulk
  real, dimension   (8,nstyp) :: soilparms
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

  data soilparms/  &
       !----------------------------------------------------------------------
       !slpots        slbs          slcons0         slden       USDA SOIL CLASS
       !      slmsts         slcons          slcpd              # AND NAME
       !----------------------------------------------------------------------
       -.121, .395,  4.05, .18e-3, .50e-3, 1465.e3, 1600.,.135  & !  1 sand
       ,-.090, .410,  4.38, .16e-3, .60e-3, 1407.e3, 1600.,.150  & !  2 loamy sand
       ,-.218, .435,  4.9 , .34e-4, .77e-3, 1344.e3, 1600.,.195  & !  3 sandy loam
       ,-.786, .485,  5.3 , .72e-5, .11e-4, 1273.e3, 1600.,.255  & !  4 silt loam
       ,-.478, .451,  5.39, .69e-5, .22e-2, 1214.e3, 1600.,.240  & !  5 loam
       ,-.299, .420,  7.12, .63e-5, .15e-2, 1177.e3, 1600.,.255  & !  6 sandy clay loam
       ,-.356, .477,  7.75, .17e-5, .11e-3, 1319.e3, 1600.,.322  & !  7 silty clay loam
       ,-.630, .476,  8.52, .24e-5, .22e-2, 1227.e3, 1600.,.325  & !  8 clay loam
       ,-.153, .426, 10.4 , .22e-5, .22e-5, 1177.e3, 1600.,.310  & !  9 sandy clay
       ,-.490, .492, 10.4 , .10e-5, .10e-5, 1151.e3, 1600.,.370  & ! 10 silty clay
       ,-.405, .482, 11.4 , .13e-5, .13e-5, 1088.e3, 1600.,.367  & ! 11 clay 
       ,-.356, .863,  7.75, .80e-5, .80e-5,  874.e3,  300.,.535/   ! 12 peat

  data  xsand  /.97,.92,.80,.57,.60,.65,.35,.48,.50,.30,.25,.20/
  data  xclay  /.03,.07,.18,.40,.35,.31,.59,.45,.42,.61,.65,.20/
  data  xorgan /.00,.01,.02,.03,.05,.04,.06,.07,.08,.09,.10,.60/

  data  xrobulk/1200.,1250.,1300.,1400.,1350.,1350.  &
       ,1500.,1450.,1450.,1650.,1700., 500./

  !         LEAF-3 BIOPHYSICAL PARAMETERS BY LANDUSE CLASS NUMBER

  data bioparms/  &
       !----------------------------------------------------------------------
       !albv_green     sr_max         veg_clump       rootdep             LEAF-3 CLASS #
       !     albv_brown     tai_max        veg_frac        dead_frac      AND DESCRIPTION
       !          emisv          sai            veg_ht         rcmin
       !----------------------------------------------------------------------
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  0  Ocean
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  1  Lakes, rivers, streams
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  2  Ice cap/glacier
       .00, .00, .00,  .0, 0.0,  .0,  .0, .00,   .0,  .0, .0,   0., & !  3  Desert, bare soil
       .14, .24, .97, 5.4, 8.0, 1.0, 1.0, .80, 20.0, 1.5, .0, 500., & !  4  Evergreen needleleaf tree
       .14, .24, .95, 5.4, 8.0, 1.0, 1.0, .80, 22.0, 1.5, .0, 500., & !  5  Deciduous needleleaf tree
       .20, .24, .95, 6.2, 7.0, 1.0,  .0, .80, 22.0, 1.5, .0, 500., & !  6  Deciduous broadleaf tree
       .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 1.5, .0, 500., & !  7  Evergreen broadleaf tree
       .21, .43, .96, 5.1, 4.0, 1.0,  .0, .75,   .3,  .7, .7, 100., & !  8  Short grass
       .24, .43, .96, 5.1, 5.0, 1.0,  .0, .80,  1.2, 1.0, .7, 100., & !  9  Tall grass
       .24, .24, .96, 5.1, 1.0,  .2, 1.0, .20,   .7, 1.0, .0, 500., & ! 10  Semi-desert
       .20, .24, .95, 5.1, 4.5,  .5, 1.0, .60,   .2, 1.0, .0,  50., & ! 11  Tundra
       .14, .24, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 12  Evergreen shrub
       .20, .28, .97, 5.1, 5.5, 1.0, 1.0, .70,  1.0, 1.0, .0, 500., & ! 13  Deciduous shrub
       .16, .24, .96, 6.2, 7.0, 1.0,  .5, .80, 22.0, 1.5, .0, 500., & ! 14  Mixed woodland
       .22, .40, .95, 5.1, 5.0,  .5,  .0, .85,  1.0, 1.0, .0, 100., & ! 15  Crop/mixed farming, C3 grassland
       .18, .40, .95, 5.1, 5.0,  .5,  .0, .80,  1.1, 1.0, .0, 500., & ! 16  Irrigated crop
       .12, .43, .98, 5.1, 7.0, 1.0,  .0, .80,  1.6, 1.0, .0, 500., & ! 17  Bog or marsh
       .20, .36, .96, 5.1, 6.0, 1.0,  .0, .80,  7.0, 1.0, .0, 100., & ! 18  Wooded grassland 
       .20, .36, .90, 5.1, 3.6, 1.0,  .0, .74,  6.0,  .8, .0, 500., & ! 19  Urban and built up
       .17, .24, .95, 4.1, 7.0, 1.0,  .0, .90, 32.0, 1.5, .0, 500., & ! 20  Wetland evergreen broadleaf tree
       .16, .24, .96, 5.1, 2.0, 1.5, 1.0, .10, 20.0, 1.5, .0, 500./   ! 21  Very urban

  ! Soil constants

  ! Thermal conductivity in J/msK
  cka = 0.418684 * 0.0615
  ckw = 0.418684 * 1.45
  romin = 2655.0
  roorg = 1300.0

  slz(nzg+1) = 0.

  slfcap = -10. / 3.
  refdepth = -2.0

  do nnn = 1,nstyp
     slcons0(nnn) = soilparms(5,nnn)
     fhydraul(nnn) = log (soilparms(4,nnn) / soilparms(5,nnn)) / refdepth

     do k = 1,nzg
        slcons1(k,nnn) = soilparms(4,nnn)   ! ORIGINAL form - const with depth
        !   slcons1(k,nnn) = soilparms(5,nnn) & !TOPMODEL form-large at surface
        !        * exp(slz(k) * fhydraul(nnn))  !  and exp decrease with depth
     enddo

     slpots(nnn) = soilparms(1,nnn)
     slmsts(nnn) = soilparms(2,nnn)
     slbs(nnn)   = soilparms(3,nnn)
     slcons(nnn) = soilparms(4,nnn)
     slcpd(nnn)  = soilparms(6,nnn)
     slden(nnn)  = soilparms(7,nnn)
     sfldcap(nnn)  = soilparms(8,nnn)

     emisg(nnn) = .98
     slfc(nnn) = slmsts(nnn) * (slfcap / slpots(nnn)) ** (-1. / slbs(nnn))
     soilcp(nnn) = 0.1 - 0.07 * xsand(nnn)

     !    tmin = xsand(nnn) + xclay(nnn)
     !    romean = xorgan(nnn) * roorg + tmin * romin
     !    sporo(nnn) = 1.0 - xrobulk(nnn) / romean
     !    ratio = (xorgan(nnn) / roorg) / (tmin / romin)
     !    sorgan(nnn) = ratio / (1.0 + ratio) * (1.0 - sporo(nnn))
     !    xmin = 1.0 - sporo(nnn) - sorgan(nnn)
     !    ssand(nnn) = xsand(nnn) * xmin / tmin
     !    sclay(nnn) = xclay(nnn) * xmin / tmin

  enddo

  do nnn = 1,nvtyp+nvtyp_teb
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

  !srf - call sfcdata_sib here (latter call at rdint.f90 and rnode.f90)
  call sfcdata_sib
  call sib_init_co2
  call sib_init_pco2ap

  return
end subroutine sfcdata_sib_driver

!**********************************************************************

!itb_usp...SiB2.5 tabular vegetation/soil data entered here.
!itb_usp...User needs to be aware that 'Morph Table' data will
!itb_usp...change with new pre-run mapper values. email
!itb_usp...Ian Baker (baker@atmos.colostate.edu) for details.

subroutine sfcdata_sib

  use mem_sib

  implicit none
  !INCLUDE 'rcommons.h'

  !itb_usp...dimension table variables here..
  integer :: nnn,result
  real,dimension(12) :: sbee,sphsat,ssatco,sporos,sslope,swopt,sskew,srspsat

  real,dimension(13) :: sz1,sz2,sfvcover,schil,ssodep,srootd,sphihalf     &
       ,stran11,stran12,stran21,stran22,sref11,sref12     &
       ,sref21,sref22,svmax0,seffcon,sgslope,sgsmin       &
       ,satheta,sbtheta,strda,strdm,strop,srespcp,sslti   & 
       ,sshti,shlti,shhti,ssoref1,ssoref2

  !...................SOIL CLASSES....................
  !
  !   The soil texture catagories are (based on the 12 USDA classes):
  !   
  !   Soil  Name            % clay  % sand
  !    1    sand               3       92
  !    2    loamy sand         5       82
  !    3    sandy loam         10      65
  !    4    silt loam          13      22
  !    5    silt               7       7
  !    6    loam               18      42
  !    7    sandy clay loam    28      58
  !    8    sandy clay         40      52
  !    9    clay loam          39      32
  !   10    silty clay loam    39      10
  !   11    silty clay         41      7
  !   12    clay               65      19
  !   
  !   Soil properties based on the approximate centroid of the soil texture 
  !   catagory within the USDA texture triangle.
  !   
  !   Modifications:
  !     Kevin Schaefer added respiration variables (Wopt, skew, RespSat) 
  !       using curve fits based on data from Raich et al., 1991 (6/19/00)
  !     Kevin Schaefer updated table using centroid % clay/sand (3/30/01)
  !   
  !   Variables:   
  !   
  !    Bee      : Clapp & Hornberge 'B' exponent
  !    phisat   : Soil tension at saturation (units)
  !    satco    : soil tension at 1/2 assimilation value (true?) units?
  !    poros    : porosity
  !    slope    : slope
  !    wopt     : optimum saturation percentage for respiration
  !    skew     : 
  !    respsat  : 

  data sbee/3.387,3.705,4.500,4.977,4.023,5.772,7.362,9.270,9.111,9.111 &
       ,9.429,13.245/

  data sphsat/-0.047,-0.064,-0.107,-0.391,-0.614,-0.214,-0.132,-0.158   &
       ,-0.289,-0.561,-0.614,-0.428/

  data ssatco/0.236E-4,0.166E-4,0.910E-5,0.200E-5,0.118E-5,0.405E-5     &
       ,0.711E-5,0.576E-5,0.285E-5,0.131E-5,0.118E-5,0.180E-5/

  data sporos/0.373,0.386,0.407,0.461,0.480,0.436,0.416,0.423,0.449     &
       ,0.476,0.480,0.465/

  data sslope/0.176,0.176,0.176,0.176,0.176,0.176,0.176,0.176,0.176     &
       ,0.176,0.176,0.176/

  data swopt/59.653,60.080,61.120,61.725,60.501,62.701,64.533,66.520    &
       ,66.363,66.363,66.675,69.920/

  data sskew/0.354,0.357,0.362,0.363,0.360,0.359,0.328,0.232,0.243      &
       ,0.243,0.221,-0.255/

  data srspsat/0.508,0.513,0.525,0.533,0.518,0.545,0.570,0.600,0.598    &
       ,0.598,0.603,0.663/


  do nnn=1,nstyp_sib
     bee_sib(nnn)    = sbee(nnn)
     phsat_sib(nnn)  = sphsat(nnn)
     satco_sib(nnn)  = ssatco(nnn)
     poros_sib(nnn)  = sporos(nnn)
     slope_sib(nnn)  = sslope(nnn)
     wopt_sib(nnn)   = swopt(nnn)
     skew_sib(nnn)   = sskew(nnn)
     respsat_sib(nnn)= srspsat(nnn)
  enddo

  !...END SOIL TABLES........

  !..................VEGETATION/BIOME TYPES............................
  !
  !     The vegetation types are:
  !     # type  Name
  !     1  C3  Tall Broadleaf-Evergreen Trees
  !                  Ref: Stanford Group, Sellers et al. (1989)
  !     2  C3  Tall Broadleaf-Deciduous Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972)
  !     3  C3  Tall Broadleaf and Needleleaf Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972)
  !     4  C3  Tall Needleleaf Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Strebel et al. (1982)
  !     5  C3  Tall Needleleaf-DECIDUOUS Trees
  !                  Ref: Klink and Willmott (1985), 
  !                       Strebel et al. (1982)
  !     6  C4  Short Vegetation, Same as Types 6, 7, 
  !                       8, and 11 (Stanford-Carnegie)
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972)
  !     7  C4  Short Vegetation: Ground Cover (Maize Optical Properties)
  !                  Ref: Klink and Willmott (1985), 
  !                       Miller (1972), Sellers (PC*)
  !     8  C4  Short Vegetation: Ground Cover (Maize Optical Properties)
  !                  Ref: Klink and Willmott (1985), 
  !                       Miller (1972), Sellers (PC*)
  !     9  C3  Short Broadleaf Shrubs with Bare Soil
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972),
  !                       Sellers (PC)
  !     10 C3  Short Ground Cover (Tundra)
  !                  Ref: Klink and Willmott (1985), 
  !                       Turner (1974), Miller (1972),
  !                       Sellers (PC)
  !     11 C4  No Vegetation (Low Latitude Desert)
  !                  Ref: Sellers (PC) and DORMAN
  !     12 C3  Agriculture (Wheat) and C3 Grasslands
  !                  Ref: Sellers and Dorman (1987), 
  !                       Turner (1974), and Dorman
  !     13 C4  ice
  !                  Ref: Sellers (PC) and Dorman
  !     
  !     * personal communication
  !     
  !     
  !............VARIABLES........................
  !
  !     z2 --- canopy top height (meters)
  !     z1 --- canopy bottom height (meters)
  !     fvcover --- fractional vegetation cover (-)
  !     chil    ---
  !     sodep   --- total soil depth (meters)
  !     rootd   --- rooting depth (meters)
  !     phi_half ---
  !     trans(2,2) --- leaf transmittance
  !                    (1,1) - SW green
  !                    (1,2) - LW green
  !                    (2,1) - SW brown
  !                    (2,2) - LW brown
  !     ref(2,2) --- leaf reflectance
  !                    (1,1) - SW green
  !                    (1,2) - LW green
  !                    (2,1) - SW brown
  !                    (2,2) - LW brown
  !     vmax0 ---
  !     effcon ---
  !     gsslope ---
  !     gsmin ---
  !     atheta ---
  !     btheta ---
  !     trda ---
  !     trdm ---
  !     trop ---
  !     respcp ---
  !     slti ---
  !     shti ---
  !     hlti ---
  !     hhti ---
  !     soref(2) --- soil reflectance
  !                     (1) - visible
  !                     (2) - nir     


  data sz2/35.00,20.00,20.00,17.00,17.00,1.00,1.00,1.00,0.500       &
       ,0.600,1.00,1.00,1.00/

  data sz1/1.00,11.500,10.00,8.500,8.500,0.100,0.100,0.100          &
       ,0.100,0.100,0.100,0.100,0.100/

  data sfvcover/0.874,0.597,0.727,0.558,0.670,0.776,0.343,0.343     &
       ,0.136,0.402,0.055,0.553,0.055/

  data schil/0.100,0.250,0.125,0.010,0.010,-0.300,-0.300,-0.300     &
       ,0.010,0.200,-0.300,-0.300,-0.300/

  data ssodep/3.500,2.000,2.000,2.000,2.000,1.500,1.500,1.500       &
       ,1.500,1.500,1.500,1.500,1.500/

  data srootd/1.500,1.500,1.500,1.500,1.000,1.000,1.000,1.000       &
       ,1.000,1.000,1.000,1.000,1.000/

  data sphihalf/-200.0,-200.0,-200.0,-200.0,-200.0,-200.0           &
       ,-200.0,-200.0,-200.0,-200.0,-200.0,-200.0,-200.0/

  data stran11/0.050,0.050,0.050,0.050,0.050,0.070,0.070,0.070      &
       ,0.050,0.070,0.070,0.070,0.070/

  data stran12/0.250,0.250,0.150,0.100,0.100,0.248,0.248,0.248      &
       ,0.250,0.248,0.248,0.248,0.248/

  data stran21/0.001,0.001,0.001,0.001,0.001,0.220,0.220,0.220      &
       ,0.001,0.220,0.220,0.220,0.220/

  data stran22/0.001,0.001,0.001,0.001,0.001,0.375,0.375,0.375      &
       ,0.001,0.375,0.375,0.375,0.375/

  data sref11/0.060,0.070,0.070,0.080,0.080,0.060,0.100,0.120       &
       ,0.120,0.120,0.120,0.080,0.140/

  data sref12/0.390,0.390,0.380,0.370,0.370,0.400,0.400,0.400       &
       ,0.400,0.400,0.400,0.400,0.420/

  data sref21/0.160,0.160,0.160,0.160,0.160,0.160,0.220,0.220       &
       ,0.220,0.220,0.220,0.160,0.220/

  data sref22/0.430,0.430,0.420,0.410,0.410,0.440,0.480,0.480       &
       ,0.480,0.480,0.480,0.480,0.320/

  data svmax0/0.100E-3,0.100E-3,0.750E-4,0.600E-4,0.100E-3          &
       ,0.300E-4,0.300E-4,0.300E-4,0.600E-4,0.600E-4          &
       ,0.300E-4,0.100E-3,0.300E-4/

  data seffcon/0.080,0.080,0.080,0.080,0.080,0.050,0.050,0.050      &
       ,0.080,0.080,0.050,0.080,0.050/

  data sgslope/9.000,9.000,9.000,9.000,9.000,4.000,4.000,4.000      &
       ,9.000,9.000,4.000,9.000,4.000/

  data sgsmin/0.010,0.010,0.010,0.010,0.010,0.040,0.040,0.040       &
       ,0.010,0.010,0.040,0.010,0.040/

  data satheta/0.980,0.980,0.980,0.980,0.980,0.800,0.800,0.800      &
       ,0.980,0.980,0.800,0.980,0.800/

  data sbtheta/0.950,0.950,0.950,0.950,0.950,0.950,0.950,0.950      &
       ,0.950,0.950,0.950,0.950,0.950/

  data strda/1.300,1.300,1.300,1.300,1.300,1.300,1.300,1.300        &
       ,1.300,1.300,1.300,1.300,1.300/

  data strdm/328.16,328.16,328.16,328.16,328.16,328.16,328.16       &
       ,328.16,328.16,328.16,328.16,328.16,328.16/

  data strop/298.16,298.16,298.16,298.16,298.16,298.16,298.16       &
       ,298.16,298.16,298.16,298.16,298.16,298.16/

  data srespcp/0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015      &
       ,0.015,0.015,0.015,0.015,0.015/

  data sslti/0.200,0.200,0.200,0.200,0.200,0.300,0.300,0.300        &
       ,0.200,0.200,0.300,0.200,0.300/

  data sshti/288.16,283.16,280.16,278.16,278.16,288.16,288.16       &
       ,288.16,283.16,278.16,288.16,281.16,288.16/

  data shlti/0.300,0.300,0.300,0.300,0.300,0.300,0.300,0.300        &
       ,0.300,0.300,0.300,0.300,0.300/

  data shhti/313.16,311.16,307.16,303.16,303.16,313.16,313.16       &
       ,313.16,313.16,303.16,313.16,308.16,313.16/

  data ssoref1/0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.150      &
       ,0.300,0.110,0.300,0.100,0.300/

  data ssoref2/0.225,0.225,0.225,0.225,0.225,0.225,0.225,0.250      &
       ,0.350,0.225,0.350,0.150,0.350/


  do nnn=1,nvtyp_sib
     z2_sib(nnn)       = sz2(nnn)
     z1_sib(nnn)       = sz1(nnn)
     fvcover_sib(nnn)  = sfvcover(nnn)
     chil_sib(nnn)     = schil(nnn)
     sodep_sib(nnn)    = ssodep(nnn)
     rootd_sib(nnn)    = srootd(nnn)
     phc_sib(nnn)      = sphihalf(nnn)
     tran_sib(nnn,1,1) = stran11(nnn)
     tran_sib(nnn,1,2) = stran12(nnn)
     tran_sib(nnn,2,1) = stran21(nnn)
     tran_sib(nnn,2,2) = stran22(nnn)
     ref_sib(nnn,1,1)  = sref11(nnn)
     ref_sib(nnn,1,2)  = sref12(nnn)
     ref_sib(nnn,2,1)  = sref21(nnn)
     ref_sib(nnn,2,2)  = sref22(nnn)
     vmax0_sib(nnn)    = svmax0(nnn)
     effcon_sib(nnn)   = seffcon(nnn)
     gslope_sib(nnn)   = sgslope(nnn)
     gsmin_sib(nnn)    = sgsmin(nnn)
     atheta_sib(nnn)   = satheta(nnn)
     btheta_sib(nnn)   = sbtheta(nnn)
     trda_sib(nnn)     = strda(nnn)
     trdm_sib(nnn)     = strdm(nnn)
     trop_sib(nnn)     = strop(nnn)
     respcp_sib(nnn)   = srespcp(nnn)
     slti_sib(nnn)     = sslti(nnn)
     shti_sib(nnn)     = sshti(nnn)
     hlti_sib(nnn)     = shlti(nnn)
     hhti_sib(nnn)     = shhti(nnn)
     soref_sib(nnn,1)  = ssoref1(nnn)
     soref_sib(nnn,2)  = ssoref2(nnn)


  enddo


  !citb_usp...now need to read in data from the morph_tab.binary file.
  !citb_usp...this data will change with different ndvi/soil/biome 
  !citb_usp...data, so it can't be hardwired in. Contact Ian Baker
  !citb_usp...baker@atmos.colostate.edu with any questions.

  open(unit=34,file='morph_tab.binary',form='unformatted',iostat=result)
  if(result > 0) then    !ERROR
     print*,'ERROR OPENING SiB MORPHOLOGICAL TABLE'
     stop
  endif

  read(34,iostat=result)laig_sib
  read(34,iostat=result)fvcg_sib
  read(34,iostat=result)a_zo_sib
  read(34,iostat=result)a_zp_sib
  read(34,iostat=result)a_rbc_sib
  read(34,iostat=result)a_rdc_sib
  read(34,iostat=result)zc_sib
  read(34,iostat=result)zlw_sib
  read(34,iostat=result)zlen_sib
  read(34,iostat=result)ltmax_sib
  read(34,iostat=result)stem_sib
  read(34,iostat=result)nd98_sib
  read(34,iostat=result)nd02_sib
  read(34,iostat=result)srmax_sib
  read(34,iostat=result)srmin_sib

  close(34)


  goto 1000

  print*,'laig:',laig_sib
  print*,'fvcg:',fvcg_sib
  do nnn=1,50
     print*,'RBC:',nnn,a_rbc_sib(6,25,nnn)
  enddo

  print*,'a_zo:',a_zo_sib
  print*,'a_zp:',a_zp_sib
  print*,'a_rbc:',a_rbc_sib
  print*,'a_rdc:',a_rdc_sib
  print*,'zc:',zc_sib
  print*,'zlw:',zlw_sib
  print*,'zlen:',zlen_sib
  print*,'ltmax:',ltmax_sib
  print*,'stem:',stem_sib
  print*,'nd98:',nd98_sib
  print*,'nd02:',nd02_sib
  print*,'srmax:',srmax_sib
  print*,'srmin:',srmin_sib



  !srf ---- debug
  print*,'SIB2 '
  do nnn=1,nvtyp_sib
     print*,z2_sib(nnn)       
     print*,   z1_sib(nnn)       
     print*,   fvcover_sib(nnn)  
     print*,   chil_sib(nnn)     
     print*,   sodep_sib(nnn)    
     print*,  rootd_sib(nnn)    
     print*,  phc_sib(nnn)      
     print*,  tran_sib(nnn,1,1) 
     print*,  tran_sib(nnn,1,2) 
     print*,  tran_sib(nnn,2,1) 
     print*,tran_sib(nnn,2,2) 
     print*,ref_sib(nnn,1,1)  
     print*,ref_sib(nnn,1,2)  
     print*,ref_sib(nnn,2,1)  
     print*,ref_sib(nnn,2,2)  
     print*,vmax0_sib(nnn)    
     print*,effcon_sib(nnn)   
     print*,gslope_sib(nnn)   
     print*,gsmin_sib(nnn)    
     print*,atheta_sib(nnn)   
     print*,btheta_sib(nnn)   
     print*,trda_sib(nnn)     
     print*,trdm_sib(nnn)     
     print*,trop_sib(nnn)     
     print*,respcp_sib(nnn)   
     print*,slti_sib(nnn)     
     print*,shti_sib(nnn)     
     print*,hlti_sib(nnn)     
     print*,hhti_sib(nnn)     
     print*,soref_sib(nnn,1)  
     print*,   soref_sib(nnn,2)  
  enddo
  !srf-debug
1000 continue




  return
end subroutine sfcdata_sib

!***************************************************************************

!itb_usp...
!itb_usp...end SiB2.5 tabular data...
!itb_usp---------------------------------------------------------------
subroutine sib_init_co2

  use sib_vars
  use mem_grid, only : nnxp, nnyp, nnzp, ngrids
  use mem_scalar, only : scalar_g
  use ref_sounding, only : maxsndg

  implicit none

  real f
  data f /1.51724e-6/  ! 1.51724e-6= 44/29/1.e+6

  integer :: i,j,k,isc,ng

  co2_init(1) = max(0., co2_init(1))

  do k=2,maxsndg
     if (co2_init(k) <= 0.) co2_init(k) = co2_init(k-1)
  enddo

  do ng=1,ngrids
     do isc = 1,N_CO2
        do j = 1,nnyp(ng)
           do i = 1,nnxp(ng)
              do k=2,nnzp(ng)
!!$                 scalar_g(isc, ng)%sclp(k, i, j) = co2_init(k-1)
                 scalar_g(isc, ng)%sclp(k, i, j) = f*co2_init(k-1)
              enddo
           enddo
        enddo
     enddo

     do isc = 1,N_CO2
        do j = 1,nnyp(ng)
           do i = 1,nnxp(ng)
              scalar_g(isc, ng)%sclp(1, i, j) = scalar_g(isc, ng)%sclp(2, i, j)
           enddo
        enddo
     enddo
  enddo

  return
end subroutine sib_init_co2

!****************************************************************************

!itb...------------------------------------------------------
!itb...initialize pco2ap (CAS CO2) to 35.0 Pa

subroutine sib_init_pco2ap

  use mem_grid, only : nnxp, nnyp, nnzp, ngrids
  use mem_sib, only : sib_brams_g

  integer :: i,j,k,isc,ng
  do ng=1,ngrids
        do j = 1,nnyp(ng)
           do i = 1,nnxp(ng)
             sib_brams_g(ng)%pco2ap(i,j) = 35.0
             sib_brams_g(ng)%rst(i,j)    = 100.0
           enddo
        enddo
   enddo

  return
end subroutine sib_init_pco2ap
