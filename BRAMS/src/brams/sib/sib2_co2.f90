
!############################# Change Log ##################################
! ----------- SiB VERSION 2.5, MODIFIED FOR USE WITH RAMS

![MLO - Added the subroutine back, but I left it commented out except for the variable declaration. This is 
!       just to avoid crashes in the other subroutines compilation. 
![MLO - No way, I'm commenting out everything!!!
! SUBROUTINE sib(mzg,mzs, dt              & ! timestep
!      ,thm,        sh_x,        ps_x     & ! atmospheric forcing
!      ,ros,        ts_x,        cupr     & ! atmospheric forcing
!      ,lspr,     dswbot,        spdm     & ! atmospheric forcing
!      ,dlwbot,    zwind,        cosz     & ! atmospheric forcing
!      ,biome_f,     cndvi,       pndvi   & ! SiB boundary conditions
! !    ,respfactor_in,                    & ! SiB boundary conditions
!      ,soiltype_f                        & ! SiB boundary conditions
!      ,ta,          sha,          tc     & ! SiB prognostic variables
! !    ,tg,           td,         www     & ! SiB prognostic variables
!      ,tempk,soil_water                  &
!      ,snow1,     snow2,      capac1     & ! SiB prognostic variables
!      ,capac2,                   rst     & ! SiB prognostic variables
!      ,      co2flx     & ! output
!      ,ustar,     rstar,       tstar     & ! output
!      ,pco2ap,    pco2c                  & ! diagnostics
!     !         ,    doy          ! forcing for mapper
!      ,igrp, jgrp                        &
!      , sclp, pco2m ,latitude,doy         &! New variables introduced
!      , fss_out,fws_out,assimn_out,respg_out,rstfac4   &  ! new diagnostics
!      , rstfac1,rstfac2,rstfac3 & ! diagnostics
!      , ect_out,eci_out,egi_out,egs_out,hc_out,hg_out   &
!        ,    w1_out       & ! diag: soil water layer 1 (fraction saturated)
!        ,    w2_out       & ! diag: soil water layer 2 (fraction saturated)
!        ,    w3_out       & ! diag: soil water layer 3 (fraction saturated)
!        ,    ww1_out      & ! diag: soil water layer 1 (kg/m^2)
!        ,    ww2_out      & ! diag: soil water layer 2 (kg/m^2)
!        ,    ww3_out      & ! diag: soil water layer 3 (kg/m^2)
!        ,    exo_out      & ! diag: total soil water in excess of saturation (m)
!        ,    ta_out       & ! diag: canopy airspace temperature (K)
!        ,    tc_out       & ! diag: canopy (vegetation) temperature (K)
!        ,    tg_out       & ! diag: ground surface temperature (K)
!        ,    td1_out      & ! diag: soil layer 1 temperature (K)
!        ,    td2_out      & ! diag: soil layer 2 temperature (K)
!        ,    td3_out      & ! diag: soil layer 3 temperature (K)
!        ,    td4_out      & ! diag: soil layer 4 temperature (K)
!        ,    td5_out      & ! diag: soil layer 5 temperature (K)
!        ,    td6_out      & ! diag: soil layer 6 temperature (K)
!        ,    ra_out       & ! diag: CAS-RAMS resistance (sec/m)
!        ,    rb_out       & ! diag: leaf sfc-CAS resistance (sec/m)
!        ,    rc_out       & ! diag: total canopy resistance (sec/m)
!        ,    rd_out       & ! diag: ground-CAS resistance (sec/m)
!        ,    roff_out     & ! diag: runoff (kg/m^2/sec)
!        ,    zlt_out      & ! diag: LAI
!        ,    green_out    & ! diag: greenness fraction (-)
!        ,    apar_out     & ! diag: absorbed fraction of PAR
!        ,    nee_out      & ! diag: net ecosystem exchange (umol/m^2/sec)
!        ,    cu_out       & ! diag: momentum transfer coefficient (-)
!        ,    ct_out       & ! diag: thermal transfer coefficient (-)
!        ,    ventmf_out   & ! diag: ventilation mass flux (kg/m^2/sec)
!        ,    pco2c_out    & ! diag: leaf chloroplast CO2 concentration (Pa)
!        ,    pco2i_out    & ! diag: leaf internal CO2 concentration (Pa)
!        ,    pco2s_out    & ! diag: leaf surface CO2 concentration (Pa)
!        ,    ea_out       & ! diag: CAS water vapor pressure (hPa)
!        ,    sha_out      & ! diag: CAS water vapor mixing ratio (kg/kg)
!        ,    em_out       & ! diag: ref level vapor pressure (hPa)
!        ,    rha_out      & ! diag: CAS relative humidity
!        ,    radvbc_out   & ! diag: radiation: visible beam (W/m^2)
!        ,    radvdc_out   & ! diag: radiation: visible diffuse (W/m^2)
!        ,    radnbc_out   & ! diag: radiation: nir beam (W/m^2)
!        ,    radndc_out   & ! diag: radiation: nir diffuse (W/m^2)
!        ,    dlwbot_out   & ! diag: radiation: longwave (W/m^2)
!        ,    cp_out       & ! diag: specific heat of dry air (J kg^-1 deg^-1)
!        ,    rho_out      & ! diag: density (kg m^-3)
!        ,    psy_out      & ! diag: psychrometric constant (hPa deg^-1)
!        ,    cupr_out     & ! diag: conv precip rate (mm/hr)
!        ,    lspr_out     ) ! diag: large-scale precip rate (mm/hr)


!   ! Now "co2flx_rams" receive data from "sib_g(ng, 1)%SRC_CO2(i, j)"


!   !itb...END OF SIB ROUTINE IS APPROX LINE 1275

!   USE mem_sib

!   !USE ref_sounding, ONLY : PS  ! Ver com Saulo

!   !itb_usp...intentionally keeping the diagnostics to a minimum right
!   !itb_usp...now. We can add more later.


!   !itb_usp...variables with the '_x' subscript were modified to avoid
!   !itb_usp...conflicts with RAMS common block variables...

!   IMPLICIT NONE

!   REAL :: sclp, pco2m  ! New variables introduced

!   !INCLUDE 'rcommons.h'


!   !-------------------------------------------------------------------
!   !     REFERENCES: Sato, N., P. J. Sellers, D. A. Randall, E. K. Schneider,
!   !          J. Shukla, J. L Kinter III, Y-T, Hou, and Albertazzi (1989)
!   !          "Effects of implementing the simple biosphere model in a general
!   !          circulation model. J. Atmos. Sci., 46, 2767-2782.

!   !                 Sellers, P. J., D. A. Randall, C. J. Collatz, J. A. Berry,
!   !          C. B. Field, D. A. Dazlich, C. Zhang, G. Collelo (1996) A revise
!   !          land-surface parameterization (SiB2) for atmospheric GCMs. Part 1:
!   !          Model formulation. (accepted by JCL)


!   !     MODIFICATIONS:
!   !        - changed VQSAT call to VNQSAT.  kwitt 10/23
!   !          - added in the prognostic stomatal conductance in addinc. changan
!   !          - moved sib diagnostics accumulation from dcontrol to bldif
!   !                dd 950202

!   !     SUBROUTINES CALLED:  VNQSAT, SNOW1, balan, VNTLAT
!   !          DELHF, DELEF, NETRAD, SIBSLV, endtem, updat2, addinc
!   !          inter2, balan, soilprop, soiltherm, begtem, rnload

!   !     FUNCTIONS CALLED:
!   !          none
!   !srf
!   !     rams input/output
!   INTEGER :: mzs,mzg,igrp,jgrp
!   REAL, DIMENSION(mzg+mzs) :: tempk
!   REAL, DIMENSION(mzg) :: soil_water




!   !itb_usp...FIXING LEN AT 1-WILL RUN AS A SINGLE POINT...
!   INTEGER, PARAMETER :: len=1
!   INTEGER, PARAMETER :: nsib=1

!   REAL, PARAMETER :: num_pi = 3.1415926
!   REAL, PARAMETER :: grav = 9.81
!   REAL, PARAMETER :: cp = 1004.
!   REAL, PARAMETER :: cv = 1952.
!   REAL, PARAMETER :: rgas = 287.
!   REAL, PARAMETER :: hltm = 2.52E6
!   REAL, PARAMETER :: delta = 0.608
!   REAL, PARAMETER :: asnow = 16.7
!   REAL, PARAMETER :: kapa = 0.2861328125
!   REAL, PARAMETER :: snomel = 3.705185e8
!   REAL, PARAMETER :: clai = 4.186*1000.0*0.2
!   REAL, PARAMETER :: cww = 4.186*1000.0*1000.
!   REAL, PARAMETER :: pr0 = 0.74
!   REAL, PARAMETER :: ribc = 3.05
!   REAL, PARAMETER :: vkrmn = 0.35
! !!$  REAL, PARAMETER :: pco2m = 35.0 ! Now it is passed as argument
!   REAL, PARAMETER :: po2m = 20900.
!   REAL, PARAMETER :: stefan = 5.67e-8
!   REAL, PARAMETER :: grav2 = grav *0.01
!   REAL, PARAMETER :: day = 24.0 * 3600.0
!   REAL, PARAMETER :: gamfac = HLTM*5417.9827/CP
!   REAL, PARAMETER :: tice = 273.16
!   REAL, PARAMETER :: snofac = hltm/ ( hltm + SNOMEL * 1.E-3 )




!   !     Argument list variables
!   !       Intent: in
!   INTEGER  nsoil    ! array and loop bounds
!   INTEGER ioffset   ! subdomain offset

!   !itb_usp...these are giving me weird compiler errors right now...
!   LOGICAL :: sibdrv         &    ! run SiB with prescribe meteorology
!        ,       forcerestore &    ! do force-restore soil thermodynamics
!        ,       dotkef       &    ! use changan instead of deardorff flux
!        !  parameterization
!        ,       louis


!   REAL  zwind, ztemp

!   !itb_usp...this one is giving me a compiler error also...
!   REAL   &  !  atmospheric forcing
!        ps_x(len)          &! surface pressure (hPa)
!        ,    bps(len)      &! (ps/1000)**kapa
!        ,    ros (len)     &! surface air density (kg/m^3)
!        ,    ts_x(len)     &! surface mixed layer air temperature
!        ,    psb(len)      &! boundary layer mass depth (hPa)
!        ,    cupr(len)     &! convective precipitation rate (mm/s)
!        ,    lspr(len)     &! stratiform precipitation rate (mm/s)
!        ,    radvbc(len)   &! surface incident visible direct beam (W/m^2)
!        ,    radnbc(len)   &! surface incident near IR direct beam (W/m^2)
!        ,    radvdc(len)   &! surface incident visible diffuse beam (W/m^2)
!        ,    radndc(len)   &! surface incident near IR diffuse beam (W/m^2)
!        ,    zb(len)       &! boundary layer thickness (m)
!        ,    spdm(len)     &! boundary layer wind speed (m/s)
!        ,    dlwbot(len)   &! surface incident longwave radiation (W/m^2)
!        ,    dswbot(len)   &! surface incident shortwave radiation (W/m^2)
!        ,    cosz(len)      ! cosine of solar zenith angle

!   REAL biome_f,soiltype_f
!   INTEGER :: biome   &! biome type, sent in as a single value
!        ,          soiltype ! soil type, sent in as a single value

!   REAL                    &!  surface parameters
!        z0d(len)           &! surface roughness length (m)
!        ,    zlt(len)      &! leaf area index
!        ,    z1(len)       &
!        ,    z2(len)       &
!        ,    cc1(len)      &
!        ,    cc2(len)      &
!        ,    dd(len)       &
!        ,    poros(len)    &     ! soil porosity
!        ,    zdepth(len,3) & ! porosity * soil hydrology model layer depths (m)
!        ,    phsat(len)    &
!        ,    bee(len)      &
!        ,    respcp(len)   &
!        ,    vmax0(len)    &
!        ,    green(len)    &
!        ,    tran(len,2,2) &
!        ,    ref(len,2,2)  &
!        ,    gmudmu(len)   &
!        ,    trop(len)     &
!        ,    phc(len)      &
!        ,    trda(len)     &
!        ,    trdm(len)     &
!        ,    slti(len)     &
!        ,    shti(len)     &
!        ,    hltii(len)    &
!        ,    hhti(len)     &
!        ,    effcon(len)   &
!        ,    binter(len)   &
!        ,    gradm(len)    &
!        ,    atheta(len)   &
!        ,    btheta(len)   &
!        ,    aparc(len)    &
!        ,    wopt(len)     &
!        ,    zm_x(len)     &
!        ,    wsat(len)     &
!        ,    vcover(len)   &  ! vegetation cover fraction
!        ,    sodep(len)    &
!        ,    rootd(len)    &
!        ,    soref(len,2)  &
!        ,    thermk(len)

!   REAL :: satco(len)      &
!        ,    slope(len)    &
!        ,    chil(len)     &
!        ,    ztdep(len,nzg_sib) & ! soil thermal model layer depths (m)
!        ,    sandfrac(len) &    ! soil texture sand fraction
!        ,    radfac(len,2,2,2)

!   REAL &    !  constants
!        dt        ! time step (s)
!   !    Intent: in/out
!   REAL &
!        tc(len)                &! canopy (vegetation) temperature (K)
!        ,    ta(len)           &! CAS temperature (K)
!        ,    sha(len)          &! CAS water vapor mixing ratio (kg/kg)
!        ,    cas_cap_heat(len) &! CAS heat capacity (J/K m^2)
!        ,    cas_cap_vap(len)  &! CAS
!        ,    cas_cap_co2(len)  &! CAS CO2 capacity (m / m^2)
!        ! will be moles air / m^2 in phosib
!        ,    tg(len)           &! surface boundary temperature (K)
!        ,    td(len,nzg_sib)   &! deep soil temperature (K)
!        ,    www(len,3)        &! soil wetness
!        ,    snow1(len)        &! vegetation snow
!        ,    snow2(len)        &! ground surface snow
!        ,    snow(len,2)       &! snow cover ( kg /m^2)
!        ,    capac1(len)       &! vegetation liquid store
!        ,    capac2(len)       &! ground surface liquid interception store
!        ,    capac(len,2)      &! liquid interception store (kg/m^2)
!        ,    rst(len)          &! stomatal resistance
!        ,    tke(len)          &! turbulent kinetic energy
!        ,    thm(len)          &! mixed layer potential temperature (K)
!        ,    sh_x(len)          ! mixed layer water vapor mixing ratio (kg/kg)
!   ! sh_x valid only when len=1

!   REAL    respfactor(len,nzg_sib+1)
!   !     real         respfactor_in(5)

!   !     Intent:out
!   REAL                   &
!        fss(len)          &! surface sensible heat flux (W/m^2)
!        ,    fws(len)     &! surface evaporation (kg/m^2/s)
!        ,    drag(len,2)  &! surface drag coefficient
!        ,    co2flx(len)  &! surface CO2 flux
!        ,    cflux(len)   &! new formulation of CO2 flux (phosib)
!        ,    cu(len)      &
!        ,    ct(len)      &
!        ,    ventmf(len)  &! ventilation mass flux
!        ,    thvgm(len)   &
!        ,    xgpp(len)     ! gross primary productivity (micromoles/m**2/s)

!   !     Local variables:
!   REAL          &
!        dtinv    &
!        ,    AUXeadem

!   REAL &
!        psy(len)         &! psychrometric 'constant'
!        ,    tha(len)    &! canopy airspace potential temperature (K)
!        ,    ea(len)     &! canopy airspace water vapor pressure (hPa)
!        ,    em(len)     &! mixed layer water vapor pressure (hPa)
!        ,    d(len)      &! dd corrected for snow covered canopy
!        ,    rbc(len)    &! cc1 corrected for snow covered canopy
!        ,    rdc(len)    &! cc2 corrected for snow covered canopy
!        ,    etmass(len) &! evapotranspiration
!        ,    totwb(len)  &! total surface and soil water at beginning of timest.
!        ,    chf(len)    &! canopy heat flux (W/m^2)
!        ,    ahf(len)    &! CAS heat flux (W/m^2)
!        ,    shf(len)    &! soil heat flux (W/m^2)
!        ,    ect(len)    &! transpiration (J)
!        ,    eci(len)    &! canopy interception evaporation (J)
!        ,    egs(len)    &! soil evaporation (J)
!        ,    egi(len)    &! ground interception evaporation (J)
!        ,    hc(len)     &
!        ,    hg(len)     &
!        ,    hs_x(len)   &
!        ,    ecdif(len)  &
!        ,    egdif(len)  &
!        ,    esdif(len)  &
!        ,    heaten(len) &
!        ,    hflux(len)  &
!        ,    tgs(len)    &! bare ground and snow surface mean temperature
!        ,    tsnow(len)  &
!        ,    czc(len)    &! canopy heat capacity
!        ,    btc(len)    &
!        ,    btg(len)    &
!        ,    btci(len)   &
!        ,    btct(len)   &
!        ,    btgs(len)   &
!        ,    bts(len)    &
!        ,    rstfac(len,4)

!   REAL      rstfac1      & ! diag: stress factor 1: leaf sfc humidity
!        ,    rstfac2      & ! diag: stress factor 2: soil moisture
!        ,    rstfac3      & ! diag: stress factor 3: temperature
!        ,    rstfac4      & ! diag: stress factor 4: combination of 1-3
!        ,    ect_out      & ! diag: transpiration flux (W/m^2)
!        ,    eci_out      & ! diag: canopy interception flux (W/m^2)
!        ,    egi_out      & ! diag: ground interception flux (W/m^2)
!        ,    egs_out      & ! diag: ground surface layer evap (W/m^2)
!        ,    hc_out       & ! diag: canopy (vegetation) sensible heat flux (W/m^2)
!        ,    hg_out       & ! diag: ground surface sensible heat flux (W/m^2)
!        ,    assimn_out   & ! diag: net carbon assimilation (umol/m^2/sec)
!        ,    respg_out    & ! diag: ground respiration flux (umol/m^2/sec)
!        ,    fss_out      & ! diag: sensible heat flux (W/m^2)
!        ,    fws_out      & ! diag: latent heat flux (W/m^2)
!        ,    w1_out       & ! diag: soil water layer 1 (fraction saturated)
!        ,    w2_out       & ! diag: soil water layer 2 (fraction saturated)
!        ,    w3_out       & ! diag: soil water layer 3 (fraction saturated)
!        ,    ww1_out      & ! diag: soil water layer 1 (kg/m^2)
!        ,    ww2_out      & ! diag: soil water layer 2 (kg/m^2)
!        ,    ww3_out      & ! diag: soil water layer 3 (kg/m^2)
!        ,    exo_out      & ! diag: total soil water in excess of saturation (kg/m^2)
!        ,    ta_out       & ! diag: canopy airspace temperature (K)
!        ,    tc_out       & ! diag: canopy (vegetation) temperature (K)
!        ,    tg_out       & ! diag: ground surface temperature (K)
!        ,    td1_out      & ! diag: soil layer 1 temperature (K)
!        ,    td2_out      & ! diag: soil layer 2 temperature (K)
!        ,    td3_out      & ! diag: soil layer 3 temperature (K)
!        ,    td4_out      & ! diag: soil layer 4 temperature (K)
!        ,    td5_out      & ! diag: soil layer 5 temperature (K)
!        ,    td6_out      & ! diag: soil layer 6 temperature (K)
!        ,    ra_out       & ! diag: CAS-RAMS resistance (sec/m)
!        ,    rb_out       & ! diag: leaf sfc-CAS resistance (sec/m)
!        ,    rc_out       & ! diag: total canopy resistance (sec/m)
!        ,    rd_out       & ! diag: ground-CAS resistance (sec/m)
!        ,    roff_out     & ! diag: runoff (kg/m^2/sec)
!        ,    zlt_out      & ! diag: LAI
!        ,    green_out    & ! diag: greenness fraction (-)
!        ,    apar_out     & ! diag: absorbed fraction of PAR
!        ,    nee_out        ! diag: net ecosystem exchange (umol/m^2/sec)

!   REAL      cu_out       & ! diag: momentum transfer coefficient (-)
!        ,    ct_out       & ! diag: thermal transfer coefficient (-)
!        ,    ventmf_out   & ! diag: ventilation mass flux (kg/m^2/sec)
!        ,    pco2c_out    & ! diag: leaf chloroplast CO2 concentration (Pa)
!        ,    pco2i_out    & ! diag: leaf internal CO2 concentration (Pa)
!        ,    pco2s_out    & ! diag: leaf surface CO2 concentration (Pa)
!        ,    ea_out       & ! diag: CAS water vapor pressure (hPa)
!        ,    sha_out      & ! diag: CAS water vapor mixing ratio (kg/kg)
!        ,    em_out       & ! diag: ref level vapor pressure (hPa)
!        ,    rha_out      & ! diag: CAS relative humidity
!        ,    radvbc_out   & ! diag: radiation: visible beam (W/m^2)
!        ,    radvdc_out   & ! diag: radiation: visible diffuse (W/m^2)
!        ,    radnbc_out   & ! diag: radiation: nir beam (W/m^2)
!        ,    radndc_out   & ! diag: radiation: nir diffuse (W/m^2)
!        ,    dlwbot_out   & ! diag: radiation: longwave (W/m^2)
!        ,    cp_out       & ! diag: specific heat of dry air (J kg^-1 deg^-1)
!        ,    rho_out      & ! diag: density (kg m^-3)
!        ,    psy_out      & ! diag: psychrometric constant (hPa deg^-1)
!        ,    cupr_out     & ! diag: conv precip rate (mm/hr)
!        ,    lspr_out       ! diag: large-scale precip rate (mm/hr)

!   REAL      rsoil(len)    &
!        ,    hr(len)       &
!        ,    wc(len)       &
!        ,    wg(len)       &
!        ,    satcap(len,2) &
!        ,    areas(len)    &
!        ,    csoil(len)    &
!        ,    gect(len)     &
!        ,    geci(len)     &
!        ,    gegs(len)     &
!        ,    gegi(len)     &
!        ,    rb(len)       &
!        ,    zmelt(len)    &
!        ,    rd(len)       &
!        ,    rds(len)      &
!        ,    ra(len)       &
!        ,    rib(len)      &! soil resistance
!        ,    rc_x(len)     &
!        ,    ecmass(len)   &
!        ,    hrr(len)      &
!        ,    bintc(len)    &
!        ,    aparkk(len)   &
!        ,    wsfws(len)    &
!        ,    wsfht(len)    &
!        ,    wsflt(len)    &
!        ,    wci(len)      &
!        ,    whs(len)      &
!        ,    omepot(len)   &
!        ,    assimpot(len)

!   REAL     &
!        assimci(len)                 &
!        ,    antemp(len)             &
!        ,    assimnp(len)            &
!        ,    wags(len)               &
!        ,    wegs(len)               &
!        ,    pfd(len)                &
!        ,    assim(len)              &
!        ,    zmstscale(len,2)        &
!        ,    zltrscale(len)          &
!        ,    zmlscale(len)           &
!        ,    drst(len)               &    ! stomatal resistance increment
!        ,    soilq10(len,nzg_sib+1)  &
!        ,    ansqr(len)              &
!        ,    soilscaleold(len)       &
!        ,    anwtc(len)              &
!        ,    hgdtg(len)              &
!        ,    hgdth(len)              &
!        ,    hgdta(len)              &
!        ,    hsdts(len)              &
!        ,    hsdtc(len)              &
!        ,    hsdth(len)              &
!        ,    hsdta(len)              &
!        ,    hcdtg(len)              &
!        ,    hcdtc(len)              &
!        ,    hcdth(len)              &
!        ,    hcdta(len)              &
!        ,    hadtg(len)              &
!        ,    hadtc(len)              &
!        ,    hadth(len)              &
!        ,    hadta(len)              &
!        ,    aag(len)                &
!        ,    aac(len)                &
!        ,    aam(len)                &
!        ,    fc(len)                 &
!        ,    fg(len)

!   !Bio these are the derivatives of the LW fluxes needed for
!   !Bio the prognostic canopy

!   REAL                   &
!        lcdtc(len)        &
!        ,    lcdtg(len)   &
!        ,    lcdts(len)   &
!        ,    lgdtg(len)   &
!        ,    lgdtc(len)   &
!        ,    lsdts(len)   &
!        ,    lsdtc(len)

!   REAL     &
!        eg(len)     &
!        ,    ec(len)     &
!        ,    es(len)     &
!        ,    egdtg(len)     &
!        ,    egdqm(len)     &
!        ,    ecdtg(len)     &
!        ,    ecdtc(len)     &
!        ,    ecdqm(len)     &
!        ,    egdtc(len)     &
!        ,    deadtg(len)     &
!        ,    deadtc(len)     &
!        ,    deadqm(len)     &
!        ,    ecdea(len)     &
!        ,    egdea(len)     &
!        ,    esdts(len)     &
!        ,    esdea(len)     &
!        ,    eadea(len)     &
!        ,    eadem(len)     &
!        ,    bbg(len)     &
!        ,    bbc(len)     &
!        ,    bbm(len)     &
!        ,    radt(len,3)     &
!        ,    dtg(len,2)   &! surface ground and snow temperature increments (K)
!        ,    dtc_x(len)    &! canopy temperature increment (K)
!        ,    dta(len)    &! CAS temperature increment (K)
!        ,    dea(len)    &! CAS moisture increment (Pa)
!        ,    dtd(len,nzg_sib)  &! deep soil temperature increments (K)
!        ,    q3l(len)   &
!        ,    q3o(len)   &
!        ,    exo(len)   &
!        ,    qqq(len,3)   &
!        ,    zmelt1(len)   &
!        ,    evt(len)   &
!        ,    eastar(len)   &
!        ,    roffo(len)   &
!        ,    zmelt2(len)   &
!        ,    rha(len)        &! canopy airspace relative humidity
!        ,    ggl(len)   &
!        ,    egmass(len)

!   REAL ::   etc(len)    &
!        ,    etg(len)    &
!        ,    etci(len)   &
!        ,    etgi(len)   &
!        ,    etct(len)   &
!        ,    etgs(len)   &
!        ,    ets(len)   &
!        ,    slamda(len,nzg_sib)   &  ! soil thermal conductivities
!        ,    shcap(len,nzg_sib)   &   ! soil heat capacities
!        ,    fac1(len)   &
!        ,    dth(len)   &   ! mixed layer potential temperature increment (K)
!        ,    dqm(len) & ! mixed layer water vapor mixing ratio increment (kg/kg)
!        ,    roff(len)&   ! total runoff (surface and subsurface)
!        ,    assimn(len)&
!        ,    soilscale(len,nzg_sib+1)&
!        ,    czh(len)   &! surface layer heat capacity
!        ,    radc3(len,2)&
!        ,    radn(len,2,2)&
!        ,    ustar(len) & ! friction velocity (m/s)
!        ,    ustaro(len) & ! friction velocity (m/s) ( for oceanic z0)
!        ,    rstar(len)   &!
!        ,    tstar(len)   &!
!        ,    cuo(len)   &! ( for oceanic z0)
!        ,    z0(len)   &! surface roughness length corrected for canopy snow (m)
!        ,    wwwtem(len,3)  &! soil wetness copy
!        ,    thgeff(len), tgeff(len)&
!        ,    shgeff(len), canex(len)&
!        ,    cuprt(len), lsprt(len) &! copies of cupr and lspr
!        ,    thmtem(len), shtem(len) & ! copies of thm and sh
!        ,    zzwind(len), zztemp(len)&
!        ,    respg(len)&
!        ,    discrim(len),discrim2(len),discrim3(len)&
!        ,    closs(len)     &! vegetation IR loss
!        ,    gloss(len)     &! ground IR loss
!        ,    sloss(len)     &! snow IR loss
!        ,    dtc4(len)      &! 1st derivative of vegetation T^4
!        ,    dtg4(len)      &! 1st derivative of ground T^4
!        ,    dts4(len)      &! 1st derivative of snow T^4
!        !itb    ! discrim factors-neil suits' programs
!        ,    pco2i(len)     &! leaf internal pCO2 (Pa)
!        ,    pco2ap(len)     ! canopy air space pCO2 (Pa)

!   REAL ::   pco2ap_old(len)& !previous time step pCO2 for cfrax eqns
!        ,    pco2c(len)     &! chloroplast pCO2 (Pa)
!        ,    pco2s(len)     &! leaf surface pCO2 (Pa)
!        ,    co2cap(len)   ! moles of air in the canopy (moles/canopy air space)

!   TYPE biome_morph_var
!      REAL zc        ! Canopy inflection height (m)
!      REAL lwidth    ! Leaf width
!      REAL llength   ! Leaf length
!      REAL laimax    ! Maximum LAI
!      REAL stems     ! Stem area index
!      REAL ndvimax   ! Maximum NDVI
!      REAL ndvimin   ! Minimum NDVI
!      REAL srmax     ! Maximum simple ratio
!      REAL srmin     ! Minimum simple ratio
!   END TYPE biome_morph_var
!   TYPE(biome_morph_var) morphtab

!   TYPE aero_var
!      REAL zo       ! Canopy roughness coeff
!      REAL zp_disp  ! Zero plane displacement
!      REAL rbc      ! RB Coefficient
!      REAL rdc      ! RC Coefficient
!   END TYPE aero_var
!   TYPE(aero_var),DIMENSION(50,50) :: aerovar ! aerodynamic interpolation tables

!   TYPE time_dep_var
!      REAL fpar    ! Canopy absorbed fraction of PAR
!      REAL lai     ! Leaf-area index
!      REAL green   ! Canopy greeness fraction of LAI
!      REAL zo      ! Canopy roughness coeff
!      REAL zp_disp ! Zero plane displacement
!      REAL rbc     ! RB Coefficient (c1)
!      REAL rdc     ! RC Coefficient (c2)
!      REAL gmudmu  ! Time-mean leaf projection
!   END TYPE time_dep_var
!   TYPE(time_dep_var) timevar


!   !itb_usp...some extra variables that will be passed in by RAMS
!   REAL  pndvi   ! past    value of ndvi for the gridcell
!   REAL  cndvi   ! current value of ndvi for the gridcell


!   !itb_usp...some extra stuff, stick it here...
!   REAL salb(len,2,2)
!   REAL tgeff4(len)
!   REAL latitude,ftime
!   INTEGER doy
!   INTEGER i, j, k, n, ksoil, l
!   INTEGER, DIMENSION(21) :: leaf_biome_map

!   integer iyear1,imonth1,idate1,ihour1

!   integer, external :: julday

!   DATA leaf_biome_map /1,13,11,4,5,2,1,6,6,9,10,9,9,3,12,12,12,7,11,1,11/

!   !-------------------------------------------------------------------
!   !DATA ubmin /.25/   ! should use ubmin=1.0 for convec case
!   REAL, PARAMETER :: ubmin = .25
!    spdm=MAX(spdm,ubmin)

!    !itb...some conversions for SiB...
!    !...pressure comes out of RAMS in Pa...
!    !!SRF      ps_x(:) = ps_x(:)/100.0
!    !...looks like ts_x (actual surface temp) need to be multiplied by 1000....
!    !      ts_x(:) = ts_x(:)*1000.0
!    !      ts_x(:) = 20.


!    DO i=1,len
!       www(i,1) = soil_water(7)
!       www(i,2) = soil_water(5)
!       www(i,3) = soil_water(1)
!  !     print*,WWW(i,:)
!       tg(i) = tempk(mzg)
!       !         do k=1,mzg-1            !srf   MZG = nzg_sib = 7
!       DO k=1,nzg_sib-1
!          td(i,k) = tempk(k)
!       ENDDO
!    ENDDO


!    !itb...some conversions...
!    !srf      cndvi = cndvi * 10.0
!    !srf      pndvi = pndvi * 10.0
!  !print*,'cndvi=',cndvi,' pndvi=',pndvi



!  !itb...assign value from sclp (lowest atm level CO2, in ppm)
!  !itb...to pco2m
!    pco2m = sclp/10.0


!  !  print*,pco2m,pco2ap,assimn,latitude,doy

!    dtinv = 1.0 / dt   ! inverse time step

!  !ogl...conversions for sib from leaf_class to biome type
!  !---patch, should be fixed so that biome_f comes in correctly---!
!  !---4,5 are left out because they map to 4,5 respectively---!
!    biome = leaf_biome_map(INT(biome_f))
!  ! biome = int(biome_f)

!    soiltype = INT(soiltype_f)

!    DO i=1,len
!       bps(i) = (ps_x(i)/1000.)**kapa
!       psb(i) = 50.0

!       snow(i,1) = snow1(i)
!       snow(i,2) = snow2(i)

!       capac(i,1) = capac1(i)
!       capac(i,2) = capac2(i)


!    ENDDO



!  !itb...print out some stuff...
!  !  PRINT*,'------------------------------------SRF'
!  ! print*,'SIB CUPR=',igrp, jgrp ,cupr
!  !  PRINT'(6g12.4)',thm,sh_x,ps_x,ros,ts_x,nzg_sib
!  !  PRINT'(5g12.4)',ta,sha,tc,tg,www(1,1)
!  !  PRINT'(5g12.4)',cupr,lspr,dswbot,spdm,dlwbot
!  !  PRINT'(5g12.4)',zwind,cosz,biome_f,cndvi,pndvi
!  !!PRINT'(5g12.4)',(respfactor(i),i=1,5) ! VER Com Saulo
!  !  PRINT'(5g12.4)',soiltype_f,www(1,2),www(1,3),snow1,snow2
!  !  PRINT'(5g12.4)',capac1,capac2,rst,latitude,doy
!  !  PRINT*,'------------------------------------SRF'
!    !srf
!    !       stop
!    !srf


!    !itb_usp...set some things here...

!    sibdrv = .TRUE.
!    forcerestore = .FALSE.
!    dotkef = .TRUE.
!    louis = .TRUE.
!    ioffset = 0

!    nsoil = nzg_sib - 1   !deep soil

!    ztemp = zwind

!    !      k = 5
!    !      do i=1, nzg_sib+1
!    !       if(i .le. 2) then
!    !         respfactor(1,i) = 0.0
!    !       else
!    !         respfactor(1,i) = respfactor_in(k)
!    !         k = k-1
!    !       endif
!    !      enddo


!    !itb_usp...first step: obtain TI BCs using biome type, etc
!    !itb_usp...this won't be exactly like was done in Owen's
!    !itb_usp...version, because I will already have a map
!    !itb_usp...of biome,soil,ndvi, etc read in-single values
!    !itb_usp...will be passed to here.


!    !itb_usp...   THESE BC'S ARE ASSIGNED IN leaf2_init.f90

!    DO i=1,len
!       bee(i)       = bee_sib(soiltype)
!       phsat(i)     = phsat_sib(soiltype)
!       satco(i)     = satco_sib(soiltype)
!       poros(i)     = poros_sib(soiltype)
!       slope(i)     = slope_sib(soiltype)
!       wopt(i)      = wopt_sib(soiltype)
!       zm_x(i)      = skew_sib(soiltype)
!       wsat(i)      = respsat_sib(soiltype)

!       !itb_usp...   THESE BC'S ARE ASSIGNED IN leaf2_init_f90

!       z2(i)        = z2_sib(biome)
!       z1(i)        = z1_sib(biome)
!       !itb_usp...will take fvcover from a table for now.
!       vcover(i)    = fvcover_sib(biome)
!       chil(i)      = chil_sib(biome)
!       sodep(i)     = sodep_sib(biome)
!       rootd(i)     = rootd_sib(biome)
!       phc(i)       = phc_sib(biome)
!       !
!       tran(i,1,1)  = tran_sib(biome,1,1)
!       tran(i,2,1)  = tran_sib(biome,2,1)
!       tran(i,1,2)  = tran_sib(biome,1,2)
!       tran(i,2,2)  = tran_sib(biome,2,2)
!       ref(i,1,1)   = ref_sib(biome,1,1)
!       ref(i,2,1)   = ref_sib(biome,2,1)
!       ref(i,1,2)   = ref_sib(biome,1,2)
!       ref(i,2,2)   = ref_sib(biome,2,2)
!       !
!       vmax0(i)     = vmax0_sib(biome)
!       effcon(i)    = effcon_sib(biome)
!       gradm(i)     = gslope_sib(biome)
!       binter(i)    = gsmin_sib(biome)
!       atheta(i)    = atheta_sib(biome)
!       btheta(i)    = btheta_sib(biome)
!       trda(i)      = trda_sib(biome)
!       trdm(i)      = trdm_sib(biome)
!       trop(i)      = trop_sib(biome)
!       respcp(i)    = respcp_sib(biome)
!       slti(i)      = slti_sib(biome)
!       hltii(i)     = hlti_sib(biome)
!       shti(i)      = shti_sib(biome)
!       hhti(i)      = hhti_sib(biome)

!       !itb_usp...taking soil reflectance from a table for now...
!       soref(i,1)   = soref_sib(biome,1)
!       soref(i,2)   = soref_sib(biome,2)

!    ENDDO

!    morphtab%zc = zc_sib(biome)
!    morphtab%lwidth = zlw_sib(biome)
!    morphtab%llength = zlen_sib(biome)
!    morphtab%laimax = ltmax_sib(biome)
!    morphtab%stems = stem_sib(biome)
!    morphtab%ndvimax = nd98_sib(biome)
!    morphtab%ndvimin = nd02_sib(biome)
!    morphtab%srmax = srmax_sib(biome)
!    morphtab%srmin = srmin_sib(biome)
!  !  PRINT*,'--',zc_sib(biome),zlw_sib(biome),srmax_sib(biome)&
!  !       , morphtab%laimax
!    !      stop

!    DO j=1,50
!       DO i=1,50
!          aerovar(i,j)%zo      = a_zo_sib(biome,i,j)
!          aerovar(i,j)%zp_disp = a_zp_sib(biome,i,j)
!          aerovar(i,j)%rbc     = a_rbc_sib(biome,i,j)
!          aerovar(i,j)%rdc     = a_rdc_sib(biome,i,j)
!       ENDDO
!    ENDDO

!    !       print*,'RBC:',a_rbc_sib
!    !       stop


!    !itb_usp......................................................



!    !itb_usp...now have to call mapper...
!    !itb_usp...will have to modify mapper for this code...
!    i=1
!  !  PRINT*,'call mapper:',vcover(i),i
!  !  PRINT*,'****************'
!  !  PRINT*,latitude
!  !  PRINT*,doy
!  !  PRINT*,pndvi
!  !  PRINT*,cndvi
!  !  PRINT*,vcover(i)
!  !  PRINT*,chil(i)
!  !  PRINT*,tran(i,1,1)
!  !  PRINT*,ref(i,1,1)
!  !  PRINT*,morphtab%zc
!  !  PRINT*,morphtab%lwidth
!  !  PRINT*,morphtab%llength
!  !  PRINT*,morphtab%laimax
!  !  PRINT*,morphtab%stems
!  !  PRINT*,morphtab%ndvimax
!  !  PRINT*,morphtab%ndvimin
!  !  PRINT*,morphtab%srmax
!  !  PRINT*,morphtab%srmin
!  !  PRINT*,'*****************'

!    !      print*,'call mapper'
!    !      do j=1,50
!    !        print*,j,laig_sib(j),fvcg_sib(j)
!    !      enddo

!  !XXXXXXXXXXXXXXXXXXXXXXXXXXXX
!  !print*,pndvi,cndvi
!    pndvi =0.6
!    cndvi =0.6
!  !XXXXXXXXXXXXXXXXXXXXXXXXXXXX
!    CALL mapper(latitude,doy,pndvi,cndvi,vcover(i)       &
!         ,chil(i),tran(i,1,1),ref(i,1,1),morphtab  &
!         ,aerovar,laig_sib,fvcg_sib,timevar)

!    DO i=1,len
!       aparc(i) = timevar%fpar
!       zlt(i) = timevar%lai
!       green(i) = timevar%green
!       z0d(i) = timevar%zo
!       dd(i) = timevar%zp_disp
!       cc1(i) = timevar%rbc
!       cc2(i) = timevar%rdc
!  !     print*,zlt(i),green(i),pndvi,cndvi



!       !itb_usp...for the time being, setting respfactor to 3.0E-6
!       !respfactor(i,:) = 3.0E-6
!       !itb_cptec...now hardwiring for RJ
!         respfactor(i,1) = 0.0
!         respfactor(i,2) = 0.0
!         respfactor(i,3) = 3.65e-7/(8.-zlt(i))
!         respfactor(i,4) = 7.62e-7/(8.-zlt(i))
!         respfactor(i,5) = 8.03e-7/(8.-zlt(i))
!         respfactor(i,6) = 4.33e-6/(8.-zlt(i))
!         respfactor(i,7) = 3.40e-6/(8.-zlt(i))


!    ENDDO

!  !  PRINT*,'after call to mapper'
!  !  PRINT*,'biome=',biome
!  !  PRINT*,'lai=',zlt
!  !  PRINT*,'green=',green
!  !  PRINT*,'z0d=',z0d
!  !  PRINT*,'dd=',dd
!  !  PRINT*,'cc1=',cc1
!  !  PRINT*,'cc2=',cc2
!    !      stop

!    !itb_usp...citb_usp...citb_usp...citb_usp...citb_usp...

!    !itb_usp...get radiation components from single SW value
!    CALL raddrv(nsib,dswbot     &
!         ,       cosz,radvbc,radvdc,radnbc,radndc)


!    !itb_usp...set up soil information
!    DO i=1,len
!       rootd(i) =     &
!            MIN(rootd(i),(sodep(i)*0.75))

!       !...soil moisture layers
!       zdepth(i,1) = 0.02 * poros(i)
!       zdepth(i,2) =                                &
!            (rootd(i) - 0.02)*poros(i)
!       zdepth(i,3) = poros(i)*sodep(i)              &
!            - ( zdepth(i,1)+zdepth(i,2) )

!       !itb...quick patch to cover some underflow problems...
!       IF(vcover(i) < zlt(i)/10.0)THEN
!          vcover(i) = vcover(i) * 10.0
!       ENDIF

!       ztdep(i,1) = 6.0 - sodep(i)
!       ztdep(i,2) = sodep(i) - rootd(i)
!       ztdep(i,3) = 8. * (rootd(i)-0.02) / 15.
!       ztdep(i,4) = 4. * (rootd(i)-0.02) / 15.
!       ztdep(i,5) = 2. * (rootd(i)-0.02) / 15.
!       ztdep(i,6) = 1. * (rootd(i)-0.02) / 15.
!    ENDDO

!    !itb_usp...now need to call rada2 to obtain albedos...
!    CALL rada2(snow,zlt,z1,z2                             &
!         ,          asnow,tg,cosz,tice,ref,tran,chil       &
!         ,          green,vcover,soref,radfac,salb,thermk      &
!         ,          tgeff4,tc,len)



!    !     some initialization, copy of soil wetness
!    DO i = 1,len
!       tsnow(i) = MIN(tg(i),tice)
!       cuprt(i) = cupr(i) * 0.001
!       lsprt(i) = lspr(i) * 0.001
!       zmelt(i) = 0.0
!       roff(i) = 0.0
!       !itb_usp...RAMS sends in water like soil, deepest layer indexed = 1
!       wwwtem(i,1) = www(i,1)/poros(i)
!       wwwtem(i,2) = www(i,2)/poros(i)
!       wwwtem(i,3) = www(i,3)/poros(i)

!       pco2ap_old(i) = pco2ap(i)

!    ENDDO

!    !     first guesses for ta and ea (see temrec 120)

!    DO I=1,len
!       THA(I) = TA(I) / BPS(I)
!    ENDDO

!    DO I=1,len
!       EA(I) = SHA(I) * PS_X(I) / (0.622 + SHA(I))
!       EM(I) = SH_X(I) * PS_X(I) / (0.622 + SH_X(I))
!    ENDDO

!    !      print'(5g12.4)',ta,sha,ps(i),bps(i)


!    !    distribute incident radiation between canopy and surface
!    CALL rnload(len, nsib,       &
!         radvbc,radvdc,radnbc,radndc,dlwbot,VCOVER,       &
!         thermk,radfac,radn,radc3)

!    DO i = 1,len
!       CANEX(i)  = 1.-( SNOW(i,1)*5.-Z1(i))/       &
!            (Z2(i)-Z1(i))
!       !         print*,canex(i),snow(i,1),z1(i),z2(i)
!       CANEX(i)  = MAX( 0.1 , CANEX(i) )
!       CANEX(i)  = MIN( 1.0 , CANEX(i) )
!       D(i)  = Z2(i) - ( Z2(i)-DD(i) ) * CANEX(i)
!  !     PRINT*,'DD', D(i)  ,Z2(i),z1(i),z0d(i),DD(i), CANEX(i)
!       Z0(i) = Z0D(i)/( Z2(i)-DD(i) ) * ( Z2(i)-D(i) )
!       RBC(i)    = CC1(i)/CANEX(i)
!       RDC(i)    = CC2(i)*CANEX(i)
!       AREAS(i)    = MIN(1. , ASNOW*SNOW(i,2))


!       SATCAP(i,1) = ZLT(i)*0.0001 * CANEX(i)

!       !
!       !c Collatz-Bounoua change satcap(2) to 0.0002
!       !      SATCAP(i,2) = 0.002
!       satcap(i,2) = 0.0002           ! lahouari
!    ENDDO

!    !    initialize energy and water budgets
!    CALL balan(1,1.0, zdepth, wwwtem, capac, cupr,     &
!         lspr, roff, etmass, totwb, radt, chf, shf,     &
!         dt, ect, eci, egs, egi, hc, hg, heaten,     &
!         hflux, snow, thm, tc, tg, tgs, td,     &
!         ps_x, kapa, nsib, len, ioffset, nsoil )

!    !         print*,'cccccccccccccccccccccccccccccccc'
!    !      print*,'3',tg(i),tc(i)
!    !         print*,'cccccccccccccccccccccccccccccccc'

!    CALL begtem(tc, tg, cp, hltm, ps_x, snomel                  &
!         ,           zlt, clai, cww, wwwtem, poros, num_pi, psy        &
!         ,           phsat, bee, czc, czh, phc                     &
!         ,           tgs, etc, etg, btc, btg, rstfac               &
!         ,           rsoil, hr,  wc, wg                            &
!         ,           snow, capac, areas, satcap, csoil, tice, grav    &
!         ,           snofac, len, nsib, forcerestore )

!    !      print*,'post begtem'
!    !      print*,'www0',www,poros(i),wwwtem
!    !         print*,'cccccccccccccccccccccccccccccccc'
!    !      print*,'4',tg(i),tc(i)
!    !         print*,'cccccccccccccccccccccccccccccccc'

!    !pl now that we have the new psy, calculate the new CAS capacities
!    !itb...PL made this max(4.,z2(:)), but we might have to boost the
!    !itb...min value upwards...
!    cas_cap_heat(:) = ros(:) * cp * MAX(4.,z2(:))
!    !itb...I think cas_cap_vap should use cv instead of cp...
!    !           cas_cap_vap(:)  = ros(:) * cp * max(4.,z2(:)) / psy(:)
!    cas_cap_vap(:)  = ros(:) * cv * MAX(4.,z2(:)) / psy(:)
!    cas_cap_co2(:)  =               MAX(4.,z2(:))  ! this goes
!                                                            ! out to phosib

!    !Bio approximate snow sfc vapor pressure and d(esnow)/dt with
!    !Bio ground surface values

!    ets(:) = etg(:)
!    bts(:) = btg(:)


!    !     CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT
!    !     LOSSES FROM CANOPY AND GROUND

!    CALL NETRAD(radc3, radt, stefan, fac1, vcover, thermk, tc,   &
!         tg, tice, dtc4, dtg4, dts4, closs, gloss, sloss,       &
!         tgeff, areas, len )
!    !         print*,'cccccccccccccccccccccccccccccccc'
!    !      print*,'5',tg(i),tc(i),rb(1)
!    !        print*, snow            !xxxxxxxxxxx
!    !         print*,'cccccccccccccccccccccccccccccccc'

!    DO I=1,len
!       THgeff(I) = Tgeff(I) / BPS(I)
!    ENDDO

!    CALL VNQSAT(1,tgeff,PS_X,SHgeff,len)

!    !     GET RESISTANCES FOR SIB
!    CALL VNTLAT(grav, tice,           &
!         pr0, ribc, vkrmn, delta, dt, tc, tg, ts_x, ps_x, zlt,           &
!         wwwtem, tgs, etc, etg,snow,                             &
!         rstfac, rsoil, hr, wc, wg, snofac,                      &
!         sh_x, z0, spdm, sha, zb, ros,cas_cap_co2,               &
!         cu, ra, thvgm, rib, ustar, rstar,tstar,                 &
!         ventmf, thm, tha, z2, d,                                &
!         fc, fg, rbc, rdc,gect,geci,gegs,gegi,                   &
!         respcp, rb, rd, rds,bps, rst, rc_x, ecmass,             &
!         ea, hrr, assimn, bintc, ta, pco2m, po2m, vmax0,         &
!         green, tran, ref,TimeVar%gmudmu, trop, trda, trdm, slti,&
!         shti, hltii, hhti, radn, effcon, binter, gradm,         &
!         atheta, btheta, aparkk, wsfws, wsfht, wsflt, wci,       &
!         whs, omepot, assimpot, assimci, antemp, assimnp,        &
!         wags, wegs, aparc, pfd, assim, td, wopt, zm_x, wsat,    &
!         soilscale, zmstscale, drst,                             &
!         soilq10, ansqr,                                         &
!         nsib, len, nsoil, forcerestore, dotkef,                 &
!         thgeff, shgeff, tke, ct, louis, zwind, ztemp,           &
!         respg, respfactor, pco2ap, pco2i, pco2c, pco2s,         &
!         co2cap,cflux)

!    !              print*,'post vntlat'
!    !         print*,'cccccccccccccccccccccccccccccccc'
!    !      print*,'6',tg(1),tc(1),rb(1)
!    !        print*, snow            !xxxxxxxxxxx
!    !         print*,'cccccccccccccccccccccccccccccccc'



!    !   this call for ustar, cu for oceanic value of z0
!    IF(louis) THEN
!       DO i = 1,len
!          zzwind(i) = z2(i)-d(i)+ zwind
!          zztemp(i) = z2(i)-d(i)+ ztemp
!       ENDDO
!       CALL VMFCALZO(PR0,RIBC,VKRMN,DELTA,GRAV               &
!            ,                PS_X,tha                                  &
!            ,                SPDM,ROS,CUo,THVGM,RIB                  &
!            ,                USTARo,zzwind,zztemp                    &
!            ,                len )
!    ELSE
!       CALL VMFCALo(PR0,RIBC,VKRMN,GRAV,SPDM                 &
!            ,              ZB,ROS,CUo,THVGM,USTARo                   &
!            ,              tha, len, thgeff, tke, dotkef )
!    ENDIF
!    DO I=1,len
!       DRAG(I,1) = ROS(I) * CU(I) * USTAR(I)
!       DRAG(I,2) = ROS(I) * CUo(I) * USTARo(I)
!       ANWTC(i) = ANTEMP(i)* TC(i)
!    ENDDO

!    !itb...calculate partial derivatives of the various heat fluxes
!    !itb...with respect to ground/canopy/snow temp, as well as
!    !itb...some other derivatives.


!    CALL DELLWF(DT,dtc4,dtg4,dts4,fac1,areas                &
!         ,       lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc, len )


!   ! PRINT*,'rb=',rb(1)

!    CALL DELHF( DT,CP,bps,ts_x,tgs,tsnow,tc,ta,ros,ra,rb,rd    &
!         ,                HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA          &
!         ,                HADTA,HADTH                                  &
!         ,                hc, hg, hs_x, fss, len)


!    CALL DELEF( DT,CP,ps_x,em,ea,ros,HRr,fc,fg                  &
!         ,                ra,rb,rd,rc_x,rsoil,snow,capac,wwwtem            &
!         ,                ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                 &
!         ,                ESDEA,EADEA,EADEM                             &
!         ,                ec,eg,es,fws,hltm,cas_cap_vap                 &
!         ,                etc,etg                                       &
!         ,                btc,btg,bts                                   &
!         ,                areas, gect,geci,gegs,gegi, psy, snofac, hr   &
!         ,                len )

!    !      print*,'after partial flux calc (dellwf,delhf,delef)'

!    !     get soil thermal properties
!    IF(.NOT.forcerestore)                                         &
!         CALL soilprop( td, tgs, slamda, shcap, wwwtem, poros, ztdep,  &
!         asnow, snow(1,2), areas, tice, snomel,         &
!         sandfrac, nsib, len, nsoil )

!    !        print*,'after soilprop'
!    !      print*,'6.4',tg(i),tc(i),dtc_x(i)


!    !     check against new code, rn derivatives may be zero

!    !print*,DT,GRAV2,CP,HLTM,tgs,tsnow
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*,td(1,nsoil),slamda(1,nsoil),num_pi
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, areas,fac1                   
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*,VENTMF,PSB,BPS,ros,psy                        
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, czh,czc,cas_cap_heat,cas_cap_vap     
!    !print*,'1ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc
!    !print*,'2ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA  
!    !print*,'3ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, HADTA,HADTH                          
!    !print*,'4ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, hc, hg, hs_x, fss                    
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, ESDEA,EADEA,EADEM                    
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, ec,eg,es,fws                         
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, etc,etg,ets                          
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, btc,btg,bts                          
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, RADT                                 
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, dtc_x, dtg, dth, dqm, dta, dea       
!    !print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!    !print*, len, sibdrv, forcerestore


!    CALL SIBSLV( DT,GRAV2,CP,HLTM,tgs,tsnow                     &
!         ,                 td(1,nsoil),slamda(1,nsoil),num_pi           &
!         ,                 areas,fac1                                   &
!         ,                 VENTMF,PSB,BPS,ros,psy                       &
!         ,                 czh,czc,cas_cap_heat,cas_cap_vap             &
!         ,                 lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc    &
!         ,                 HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA          &
!         ,                 HADTA,HADTH                                  &
!         ,                 hc, hg, hs_x, fss                            &
!         ,                 ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                &
!         ,                 ESDEA,EADEA,EADEM                            &
!         ,                 ec,eg,es,fws                                 &
!         ,                 etc,etg,ets                                  &
!         ,                 btc,btg,bts                                  &
!         ,                 RADT                                         &
!         ,                 dtc_x, dtg, dth, dqm, dta, dea               &
!         ,                 len, sibdrv, forcerestore )


!    !       print*,'after sibslv'

!    DO i = 1,len
!       radt(i,2) = (1.-areas(i))*radt(i,2)+areas(i)*radt(i,3)
!    ENDDO


!    CALL updat2(snow ,capac , snofac, ect, eci, egi,             &
!         egs, hltm, wwwtem, num_pi, czh, dtd, dtg, dtc_x, ta,dta,  &
!         dea, dt,                                                  &
!         roff, tc, td, tg, bee, poros, satco,                      &
!         slope, phsat, zdepth, ecmass, egmass,                     &
!         shf, tice, snomel, asnow, czc, csoil, chf,                &
!         hc, hg, areas, q3l, q3o,                                  &
!         qqq, zmelt1, cww, len, nsib, nsoil,                       &
!         forcerestore,etc,ea,btc                                   &
!         ,             geci,ros,cp,psy,gect,etg,btg                      &
!         ,             gegs,hr,fg,gegi,rd,rb,hcdtc,hcdta                 &
!         ,             hgdta,slamda(1,nsoil))


!    DO i = 1, len
!       EVT(i) = (55.56 *dtinv) * (ECMASS(i) + EGMASS(i))
!       wegs(i) = wegs(i) * evt(i)
!    ENDDO


!    !     get soil temperature increments
!    IF(.NOT.forcerestore)                                          &
!         CALL soiltherm( td, dtd, tgs, dtg, slamda, shcap, ztdep, dt, &
!         nsib, len, nsoil )


!    !    update prognostic variables, get total latent and sensible fluxes
!    DO i = 1,len
!       thmtem(i) = thm(i)
!       shtem(i) = sh_x(i)
!    ENDDO

!    !         print*,'cccccccccccccccccccccccccccccccc'
!    !      print*,'7.8',tg(i),tc(i)
!    !        print*, snow(:,:)            !xxxxxxxxxxx
!    !         print*,'cccccccccccccccccccccccccccccccc'

!    CALL addinc( grav2, cp, dt, hc, hg,                     &
!         ps_x, bps, ecmass,psy, ros, hltm, cas_cap_heat,      &
!         cas_cap_vap,                                         &
!         egmass, fss, fws, hflux, etmass,                     &
!         psb, td, thmtem, ts_x, shtem,                        &
!         tc, tg, ta, ea, ra, em, sha,                         &
!         dtd, dtc_x, dtg, dta, dea,                           &
!         drst, rst, bintc, len, nsib, nsoil,             &
!       fws_out,ea_out,em_out,ra_out,cp_out,rho_out,psy_out )

!    !      print*,'---------------------------------'
!    !      print*,'www0',www
!    !      print*,'wwwtem',wwwtem
!    !      print*, cuprt
!    !      print*,'---------------------------------'
!    !       print*, lsprt
!    !      print*,'---------------------------------'
!    !        print*, snow            !xxxxxxxxxxx
!    !      print*,'---------------------------------'
!    !      print*,capac
!    !      print*,'---------------------------------'
!    !       print*, wwwtem
!    !      print*,'---------------------------------'
!    !        print*,num_pi
!    !      print*,'---------------------------------'
!    !             print*,satcap   !????????????????????
!    !      print*,'---------------------------------'
!    !           print*,cww
!    !      print*,'---------------------------------'
!    !            print*,tc !xxxxxxx
!    !      print*,'---------------------------------'
!    !              print*,tg  !xxxxxxxxxxxxx
!    !      print*,'---------------------------------'
!    !                print*,clai
!    !      print*,'---------------------------------'
!    !                print*, zlt
!    !      print*,'---------------------------------'
!    !                print*,chil
!    !      print*,'---------------------------------'
!    !                print*, roff
!    !      print*,'---------------------------------'
!    !             print*,snomel
!    !      print*,'---------------------------------'
!    !           print*, zdepth
!    !      print*,'---------------------------------'
!    !            print*, ts_x
!    !      print*,'---------------------------------'
!    !             print*, tice
!    !      print*,'---------------------------------'
!    !             print*, asnow
!    !      print*,'---------------------------------'
!    !             print*, csoil
!    !      print*,'---------------------------------'
!    !             print*,satco
!    !      print*,'---------------------------------'
!    !          print*,dt
!    !      print*,'---------------------------------'
!    !          print*,vcover
!    !      print*,'---------------------------------'
!    !           print*, roffo
!    !      print*,'---------------------------------'
!    !           print*, zmelt2
!    !      print*,'---------------------------------'
!    !           print*,len, nsib, exo

!    CALL vnqsat(2, ta, ps_x, eastar, len)

!    !     inter2 replaces interc

!    CALL inter2( cuprt, lsprt, snow, capac, wwwtem , num_pi, &
!         satcap, cww, tc, tg, clai, zlt, chil, roff,           &
!         snomel, zdepth, ts_x, tice, asnow, csoil,             &
!         satco, dt, vcover, roffo, zmelt2 ,len, nsib, exo )

!    !          print*,'after inter2-done!'

!    DO i = 1,len
!       xgpp(i) = assim(i) * 1.e6
!    ENDDO
!    !     Calculate the surface flux of CO2 due to SiB (tracer T18)

!    !     The net flux to the atmosphere is given by the release of CO2
!    !     by soil respiration minus the uptake of CO2 by photosynthesis

!    !     ASSIMN is the net assimilation of CO2 by the plants (from SiB2)
!    !     respFactor*soilScale is the rate of release of CO2 by the soil

!    !     soilScale is a diagnostic of the instantaneous rate of
!    !        soil respiration (derived by Jim Collatz, similar to TEM)
!    !     respFactor is the annual total accumulation of carbon in the
!    !        previous year at each grid cell (annual total ASSIMN)
!    !        divided by the annual total of soilScale at the same grid pt.

!    !     Surface flux of CO2 used to be merely Assimn-Respg. With the
!    !     prognostic CAS, the calculation becomes
!    !
!    !     co2flux =  (CO2A - CO2M)/ra
!    !
!    !     with a temperature correction thrown in. This calculation is
!    !     performed in phosib.


!    DO i = 1,len
!       co2flx(i) = cflux(i)
!  !!$          PRINT*,i,co2flx(i)
!    ENDDO

!    !    some quantities for diagnostic output
!    DO i = 1,len
!       rha(i) = ea(i)/eastar(i)
!       zmelt(i) = zmelt1(i)+zmelt2(i)

!       !           Calculate an overall leaf conductance, which is QP2(162)
!       !           and is used as a weighting function in QP2(162 and 163)
!       ggl(i) = 1. / (rst(i) * (rb(i) + rc_x(i)))

!  !   if(igrp == 35 .and. jgrp == 35) then
!  !     print*,'inside SiB',igrp,jgrp
!  !     print'(a,4g16.6)','humidity:',ea(i),sha(i),rha(i),capac(i,1)
!  !   endif

!    ENDDO


!    !      print*,'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
!    !      print*,'www0',www,poros(i),wwwtem


!    !itb...copy soil moisture back, remember to invert for RAMS
!    DO i=1,len

!       www(i,3) = wwwtem(i,3)*poros(i)
!       www(i,2) = wwwtem(i,2)*poros(i)
!       www(i,1) = wwwtem(i,1)*poros(i)
!  !     PRINT*,'www1',www,poros(i)

!       snow1(i) = snow(i,1)
!       snow2(i) = snow(i,2)

!       capac1(i) = capac(i,1)
!       capac2(i) = capac(i,2)

!    ENDDO

!    !devolva para  o mesmo endereco de memoria
!    DO i=1,len
!       tempk(mzg) = tg(i)
!       !         do k=1,mzg-1
!       DO k=1,nzg_sib-1
!          tempk(k) = td(i,k)
!       ENDDO
!       soil_water(1)     = www(i,3)
!       soil_water(5) = www(i,2)
!       soil_water(7) = www(i,1)
!       !     print*,'www',www
!       !copia os valores para as camadas inferiores
!  !     DO k=nzg_sib-3,1,-1
!   !       soil_water(k) = www(i,1)
!    !   ENDDO
!    ENDDO


!  !SRF- Get r and T stars from heat and latent fluxes, for consistency with RAMS
!        DO i=1,len
!           tstar(i) = -    hflux(i)/(cp*ros(i)*ustar(i))  !hflux has unit W/m2
!           rstar(i) = -      fws(i)/(   ros(i)*ustar(i))  !fws has unit kg/m2s

!  !itb...paste in some diagnostics
!           rstfac1 = rstfac(i,1)
!           rstfac2 = rstfac(i,2)
!           rstfac3 = rstfac(i,3)
!           rstfac4 = rstfac(i,4)
!  !itb..also change some units...
!           ect_out = ect(i) / dt
!           eci_out = eci(i) / dt
!           egi_out = egi(i) / dt
!           egs_out = egs(i) / dt
!           hc_out  = hc(i)  / dt
!           hg_out  = hg(i)  / dt
!           assimn_out = assimn(i)*1.0e6
!           respg_out  = respg(i)*1.0e6
!           fss_out    = fss(i)
!  !         fws_out    = fws(i) * hltm
!           w1_out     = www(i,1)/poros(i)
!           w2_out     = www(i,2)/poros(i)
!           w3_out     = www(i,3)/poros(i)
!           ww1_out    = www(i,1)
!           ww2_out    = www(i,2)
!           ww3_out    = www(i,3)
!           exo_out    = exo(i)
!           ta_out     = ta(i)
!           tc_out     = tc(i)
!           tg_out     = tg(i)
!           td1_out    = td(i,1)
!           td2_out    = td(i,2)
!           td3_out    = td(i,3)
!           td4_out    = td(i,4)
!           td5_out    = td(i,5)
!           td6_out    = td(i,6)
!     !      ra_out     = ra(i)
!           rb_out     = rb(i)
!           rc_out     = rc_x(i)
!           rd_out     = rd(i)
!           roff_out   = roff(i)
!           zlt_out    = zlt(i)
!           green_out  = green(i)
!           apar_out   = aparc(i)
!           nee_out    = co2flx(i)
!           cu_out     = cu(i)
!           ct_out     = ct(i)
!           ventmf_out = ventmf(i)
!           pco2c_out  = pco2c(i)
!           pco2i_out  = pco2i(i)
!           pco2s_out  = pco2s(i)
!    !       ea_out     = ea(i)
!           sha_out    = sha(i)
!    !       em_out     = em(i)
!           rha_out    = rha(i)
!           radvbc_out = radvbc(i)
!           radvdc_out = radvdc(i)
!           radnbc_out = radnbc(i)
!           radndc_out = radndc(i)
!           dlwbot_out = dlwbot(i)
!           cupr_out = cupr(i) * 3600.
!           lspr_out = lspr(i) * 3600.
!        ENDDO


!  !                    PRINT*,'end of SiB:'
!  !                      print*,fss,   fws,assimn_out,respg_out

!  RETURN
!  END SUBROUTINE sib


!  !itb_usp..............................................................




!  SUBROUTINE raddrv(nsib,swdown                           &
!       ,       sunang,radvbc,radvdc,radnbc,radndc)


!    IMPLICIT NONE
!    !---------------------------------------------------------------------
!    !               radiation radive code to use the downward sw at bottom
!    !               and the formulation to estimate radvbc,radvdc, radndc,
!    !               radndc
!    !---------------------------------------------------------------------

!    INTEGER nsib
!    REAL swdown(nsib)
!    REAL sunang(nsib), stemp
!    REAL radvbc(nsib),radvdc(nsib)     &
!         ,      radnbc(nsib),radndc(nsib),c1,c2,c3,c4,c5,cloud,difrat   &
!         ,      vnrat

!    INTEGER i

!    C1 = 580.
!    C2 = 464.
!    C3 = 499.
!    C4 = 963.
!    C5 = 1160.

!    DO i=1,nsib
!       sunang(i) = MAX( 0.001 , sunang(i) )
!       stemp = swdown(i)
!       stemp = MAX(stemp,0.01 )
!       cloud = (c5 * sunang(i) - stemp) / (c4 * sunang(i))
!       cloud = MAX(cloud,0.)
!       cloud = MIN(cloud,1.)
!       !         cloud = max(0.58,cloud)

!       !z  use the real clouds here!
!       !         cloud = cldtot(i)
!       !         CLOUD = AMAX1(CLOUD,0.)
!       !         CLOUD = AMIN1(CLOUD,1.)

!       DIFRAT = 0.0604 / ( SUNANG(i)-0.0223 ) + 0.0683
!       IF ( DIFRAT .LT. 0. ) DIFRAT = 0.
!       IF ( DIFRAT .GT. 1. ) DIFRAT = 1.

!       DIFRAT = DIFRAT + ( 1. - DIFRAT ) * CLOUD
!       VNRAT = ( C1 - CLOUD*C2 ) / ( ( C1 - CLOUD*C3 )    &
!            + ( C1 - CLOUD*C2 ) )

!       radvbc(i) = (1.-DIFRAT)*VNRAT*stemp
!       radvdc(i) = DIFRAT*VNRAT*stemp
!       radnbc(i) = (1.-DIFRAT)*(1.-VNRAT)*stemp
!       radndc(i) = DIFRAT*(1.-VNRAT)*stemp
!       !xx         RADN(3,2) = ZLWD
!    ENDDO

!    RETURN
!  END SUBROUTINE raddrv

!  !-----------------------------------------------------------------
!  !=======================================================================
!  SUBROUTINE mapper(   &
!       lat,                &
!       DOY,                &
!       prevNDVI,           &
!       curNDVI,            &
!       fVCover,            &
!       ChiL,               &
!       LTran,              &
!       LRef,               &
!       MorphTab,           &
!       AeroVar,            &
!       LAIgrid,            &
!       fVCovergrid,        &
!       TimeVar             &
!       )
!    !=======================================================================
!    ! calculates time dependant boundary condition variables for SiB.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL lat         ! center latitude of grid cell
!    REAL curNDVI     ! FASIR NDVI values for a grid cell
!    REAL prevNDVI    ! previous month's NDVI value
!    REAL fVCover     !
!    REAL ChiL        !
!    REAL LTran(2,2)  !
!    REAL LRef(2,2)   !
!    INTEGER DOY         ! Day of Year (DOY) of ndvi input map
!    !
!    ! begin input biome dependant, physical morphology variables
!    TYPE biome_morph_var
!       REAL zc        ! Canopy inflection height (m)
!       REAL LWidth    ! Leaf width
!       REAL LLength   ! Leaf length       
!       REAL LAImax    ! Maximum LAI
!       REAL stems     ! Stem area index
!       REAL NDVImax   ! Maximum NDVI
!       REAL NDVImin   ! Minimum NDVI
!       REAL SRmax     ! Maximum simple ratio
!       REAL SRmin     ! Minimum simple ratio
!    END TYPE biome_morph_var
!    TYPE(biome_morph_var) :: MorphTab
!    !
!    ! begin input aerodynamic parameters
!    TYPE aero_var
!       REAL zo       ! Canopy roughness coeff
!       REAL zp_disp  ! Zero plane displacement
!       REAL RbC      ! RB Coefficient
!       REAL RdC      ! RC Coefficient
!    END TYPE aero_var

!    TYPE(aero_var),DIMENSION(50,50) :: AeroVar ! aerodynamic
!    !  interpolation tables

!    REAL LAIgrid(50)   ! grid of LAI values for lookup table
!    REAL fVCovergrid(50)! grid of fVCover values for interpolation table
!    !
!    ! begin time dependant, output variables
!    TYPE time_dep_var
!       REAL fPAR    ! Canopy absorbed fraction of PAR
!       REAL LAI     ! Leaf-area index
!       REAL Green   ! Canopy greeness fraction of LAI
!       REAL zo      ! Canopy roughness coeff
!       REAL zp_disp ! Zero plane displacement
!       REAL RbC     ! RB Coefficient (c1)
!       REAL RdC     ! RC Coefficient (c2)
!       REAL gmudmu  ! Time-mean leaf projection
!    END TYPE time_dep_var
!    TYPE(time_dep_var) TimeVar
!    !
!    ! begin internal variables
!    REAL prevfPAR    ! previous month's fPAR value
!    REAL, PARAMETER :: fPARmax=0.95
!    !                   ! Maximum possible FPAR corresponding to 98th percentile
!    REAL, PARAMETER :: fPARmin=0.01
!    !                   ! Minimum possible FPAR corresponding to 2nd percentile
!    !     For more information on fPARmin and fPARmax, see
!    !     Sellers et al. (1994a, pg. 3532); Los (1998, pg. 29, 37-39)
!    !

!    !-----------------------------------------------------------------------
!    ! Calculate time dependent variables
!    !-----------------------------------------------------------------------

!    !
!    ! Calculate first guess fPAR
!    ! use average of Simple Ratio (SR) and NDVI methods.
!    !
!  !  print*,'call avgapar:',prevndvi,morphtab%ndvimin,morphtab%ndvimax

!    CALL AverageAPAR (prevNDVI,                       &
!         MorphTab%NDVImin, MorphTab%NDVImax,            &
!         MorphTab%SRmin, MorphTab%SRmax,                &
!         fPARmax, fParmin, prevfPAR)

!    CALL AverageAPAR (curNDVI,                        &
!         MorphTab%NDVImin, MorphTab%NDVImax,            &
!         MorphTab%SRmin, MorphTab%SRmax,                &
!         fPARmax, fParmin, TimeVar%fPAR)
!    !
!    !
!    ! Calculate leaf area index (LAI) and greeness fraction (Green)
!    !   See S. Los et al 1998 section 4.2.
!    !
!    !   Select previous month

!    !
!    CALL laigrn (TimeVar%fPAR,                       &
!         prevfPAR,                            &
!         fPARmax,                             &
!         fVCover,                             &
!         MorphTab%stems,                      &
!         MorphTab%LAImax,                     &
!         TimeVar%Green,                       &
!         TimeVar%LAI)
!    !
!    ! Interpolate to calculate aerodynamic, time varying variables

!  !  PRINT*,'call aeroint:',timevar%lai,fvcover

!    CALL AeroInterpolate (                           &
!         TimeVar%LAI,                         &
!         fVCover,                             &
!         LAIgrid,                             &
!         fVCovergrid,                         &
!         AeroVar,                             &
!         TimeVar%zo,                          &
!         TimeVar%zp_disp,                     &
!         TimeVar%RbC,                         &
!         TimeVar%RdC)
!    !

!    ! Calculate mean leaf orientation to par flux (gmudmu)
!    CALL gmuder (lat,                                &
!         DOY,                                &
!         ChiL,                               &
!         TimeVar%gmudmu)

!    !

!    !
!    ! recalculate fPAR adjusting for Sun angle, vegetation cover fraction,
!    ! and greeness fraction, and LAI

!    CALL aparnew (TimeVar%LAI,                       &
!         TimeVar%Green,                       &
!         LTran,                               &
!         LRef,                                &
!         TimeVar%gmudmu,                      &
!         fVCover,                             &
!         TimeVar%fPAR,                        &
!         fPARmax,                             &
!         fPARmin)

!    !
!    RETURN
!  END SUBROUTINE mapper
!  !=======================================================================
!  SUBROUTINE laigrn (fPAR,fPARm,fPARmax,fVCover,stems,   &
!       LAImax,Green,LAI)
!    !=======================================================================
!    ! calculate leaf area index (LAI) and greenness fraction (Green) from fPAR.
!    ! LAI is linear with vegetation fraction and exponential with fPAR.
!    ! See Sellers et al (1994), Equations 7 through 13.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL fPAR     ! fraction of PAR absorbed by plants at current time
!    REAL fPARm    ! fraction of PAR absorbed by plants at previous time
!    REAL fPARmax  ! maximum possible FPAR corresponding to 98th percentile
!    REAL fVCover  ! vegetation cover fraction
!    REAL stems    ! stem area index for the specific biome type
!    REAL LAImax   ! maximum total leaf area index for specific biome type
!    !
!    ! begin output variables
!    REAL Green    ! greeness fraction of the total leaf area index
!    REAL LAI      ! area average total leaf area index
!    !
!    ! begin internal variables
!    REAL LAIg     ! green leaf area index at current time
!    REAL LAIgm    ! green leaf area index at previous time
!    REAL LAId     ! dead leaf area index at current time
!    !
!    ! Calculate current and previous green leaf area index (LAIg and LAIgm):
!    ! LAIg is log-linear with fPAR.  Since measured fPAR is an area average,
!    ! divide by fVCover to get portion due to vegetation.  Since fVCover can
!    ! be specified, check to assure that calculated fPAR does not exceed fPARMax.
!    !
!  !  print*,'lai_1',fpar,fvcover,fparmax
!    IF(fPAR/fVCover.GE.fPARmax) THEN
!       LAIg=LAImax
!    ELSE
!       LAIg=alog(1.-fPAR/fVCover)*LAImax/alog(1-fPARmax)
!    ENDIF
!    !
!    IF(fPARm/fVCover.GE.fPARmax) THEN
!       LAIgm=LAImax
!    ELSE
!       LAIgm=alog(1.-fPARm/fVCover)*LAImax/alog(1-fPARmax)
!    ENDIF
!    !
!    ! Calculate dead leaf area index (LAId):
!    ! If LAIg is increasing or unchanged, the vegetation is in growth mode.
!    ! LAId is then very small (very little dead matter).
!    ! If LAIg is decreasing, the peak in vegetation growth has passed and
!    ! leaves have begun to die off.  LAId is then half the change in LAIg,
!    ! assuming half the dead leaves fall off.
!    !
!    !     Growth mode dead leaf area index:
!    IF (LAIg.GE.LAIgm) LAId=0.0001
!    !
!    !     die-off (post peak growth) dead leaf area index:
!    IF (LAIg.LT.LAIgm) LAId=0.5*(LAIgm-LAIg)
!    !
!    ! Calculate area average, total leaf area index (LAI):
!    LAI=(LAIg+LAId+stems)*fVCover
!  !  print*,'laigrn1',laig,laid,stems,fvcover
!    !
!    ! Calculate greeness fraction (Green):
!    ! Greeness fraction=(green leaf area index)/(total leaf area index)
!    Green=LAIg/(LAIg+LAId+stems)
!  !  PRINT*,'end laigrn',LAI,Green,laimax

!    RETURN
!  END SUBROUTINE laigrn
!  !
!  !=======================================================================
!  SUBROUTINE Oldgmuder (Lat, DOY, ChiL, gmudmu)
!    !=======================================================================
!    ! calculates time mean leaf projection relative to the Sun.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL Lat      ! latitude in degrees
!    REAL DOY      ! day-of-year (typically middle day of the month)
!    REAL ChiL     ! leaf angle distribution factor
!    !
!    ! begin output variables
!    REAL gmudmu   ! time mean projected leaf area normal to Sun
!    !
!    ! begin internal variables
!    INTEGER itime ! time counter
!    REAL time     ! time from 0:00 Greenwhich Mean Time (GMT)
!    REAL coshr    ! cosine of the Greenwhich Meridian (GM) Hour angle
!    REAL mu       ! cosine of the Solar zenith angle
!    REAL chiv     ! dummy variable for leaf angle distribution factor
!    REAL dec      ! declination of the Sun (Solar Declination)
!    REAL sind     ! sine of the solar declination
!    REAL cosd     ! cosine of the solar declination
!    REAL pi       ! the constant pi
!    REAL pi180    ! conversion factor from degrees to radians
!    REAL aa       ! minimum possible LAI projection vs. cosine Sun angle
!    REAL bb       ! slope leaf area projection vs. cosine Sun angle
!    REAL cloud    ! (?) Cloud cover fraction
!    REAL fb       ! (?) mean cosine of Sun Angle
!    REAL swdown   ! (?) downward shortwave flux
!    REAL pardif   ! (?) PAR flux difracted into canopy by clouds
!    REAL pardir   ! (?) PAR flux directly onto canopy
!    REAL difrat   ! (?) fraction of shortwave flux defracted by clouds
!    REAL vnrat    ! (?) shortwave flux ratio of some sort
!    REAL tor      ! (?) TBD
!    REAL topint   ! (?) Total flux onto canopy adjusted for sun angle
!    REAL botint   ! (?) total PAR flux onto canopy during 24 hour period
!    !
!    ! Assign values to constants
!    DATA pi /3.141592/
!    pi180=pi/180.
!    cloud=0.5
!    !
!    ! Calculate solar declination in radians
!    dec=pi180*23.5*SIN(1.72e-2*(DOY-80))
!    !
!    ! Calculate sine and cosine of solar declination
!    sind=SIN(dec)
!    cosd=COS(dec)
!    !
!    ! Begin time loop to average leaf projection over 24 hour period
!    topint=0.
!    botint=0.
!    !
!    DO itime=1, 48, 1
!       !
!       ! Calculate time from zero Greenwhich Mean Time (GMT)
!       time=0.5*REAL(itime)
!       !
!       ! Calculate cosine of hour angle of Grenwhich Meridion (GM)
!       coshr=COS(-pi+time/24.*2.*pi)
!       !
!       ! Calculate cosine of the Sun angle (mu)
!       !     longitude=GM=0 degrees, latitude=Lat
!       mu=SIN(Lat*pi180)*sind+COS(Lat*pi180)*cosd*coshr
!       !
!       ! Ensure the cosine of Sun angle is positive, but not zero
!       !     e.g., daylight, sun angle<=89.4 degrees (about start disc set/rise)
!       mu=amax1(0.01, mu)
!       !
!       ! It looks like this code calculates the direct and difracted PAR based
!       ! on the solar constant and a cloud fraction of 0.5 at the top and
!       ! bottom of the atmosphere.  During night, mu=0.01, a constant.  These
!       ! calculations do not match the definition of G(mu)/mu described in
!       ! Bonan (1996) and Sellers (1985).
!       tor=0.7**(1./mu)
!       swdown=1375.*mu*(tor+0.271-0.294*tor)
!       difrat=0.0604/(mu-0.0223)+0.0683
!       difrat=MAX(difrat,0.)
!       difrat=MIN(difrat,1.)
!       difrat=difrat+(1.-difrat)*cloud
!       vnrat=(580.-cloud*464.)/((580.-cloud*499.)   &
!            +(580.-cloud*464.))
!       pardir=(1.-difrat)*vnrat*swdown
!       pardif=difrat*vnrat*swdown
!       topint=topint+pardir*mu+pardif*0.5
!       botint=botint+pardir+pardif
!       !
!    ENDDO
!    !
!    ! Calculate what looks like a mean value of something
!    fb=topint/botint
!    !
!    ! Calculate min and slope of LAI projection
!    chiv=ChiL
!    IF (ABS(chiv) .LE. 0.01) chiv=0.01
!    !   calculate minimum value of projected leaf area
!    aa=0.5-0.633*chiv-0.33*chiv*chiv
!    !   calculate slope of projected leaf area wrt cosine sun angle
!    bb=0.877*(1.-2.*aa)
!    !
!    ! Calculate mean projected LAI in Sun direction assuming fb approximates
!    ! the mean cosine of the sun angle
!    gmudmu=(aa+bb*fb)/fb
!    !
!    RETURN
!  END SUBROUTINE Oldgmuder
!  !
!  !=======================================================================
!  SUBROUTINE AeroInterpolate (LAI, fVCover,          &
!       LAIgrid,fVCovergrid,AeroVar,zo,zp_disp,RbC,RdC)
!    !=======================================================================
!    ! This subroutine calculates the aerodynamic parameters by bi-linear
!    ! interpolation from a lookup table of previously calculated values.
!    ! The interpolation table is a numpts x numpts LAI/fVCover grid with
!    ! LAI ranging from 0.02 to 10 and fVCover ranging from 0.01 to 1.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL LAI            ! actual area averaged LAI for interpolation
!    REAL fVCover        ! vegetation cover fraction for interpolation
!    REAL LAIgrid(50)    ! grid of LAI values for lookup table
!    REAL fVCovergrid(50)! grid of fVCover values for interpolation table
!    TYPE aero_var
!       REAL zo      ! Canopy roughness coeff
!       REAL zp_disp ! Zero plane displacement
!       REAL RbC     ! RB Coefficient
!       REAL RdC     ! RC Coefficient
!    END TYPE aero_var
!    TYPE(aero_var), dimension(50,50) :: AeroVar ! interpolation tables
!    !
!    ! begin output variables
!    REAL RbC            ! interpolated Rb coefficient
!    REAL RdC            ! interpolated Rd coefficient
!    REAL zo             ! interpolated roughness length
!    REAL zp_disp        ! interpolated zero plane displacement
!    !
!    ! begin internal variables
!    INTEGER i           ! index for LAI grid location
!    INTEGER j           ! index for fVCover grid location
!    REAL LocLAI         ! local LAI var. to prevent changing main LAI value
!    REAL LocfVCover     ! local fVCover var. to prevent changing fVCover value
!    REAL DLAI           ! grid spacing between LAI values in tables
!    REAL DfVCover       ! grid spacing between fVCover values in tables

!  !  !print*,'aerointerp:lai,fvc=',lai,fvcover
!    !
!    ! calculate grid spacing (assumed fixed)
!    DLAI=LAIgrid(2)-LAIgrid(1)
!    DfVCover=fVCovergrid(2)-fVCovergrid(1)
!    !
!    ! Assign input LAI and fVCover to local variables and make sure
!    ! they lie within the limits of the interpolation tables, assuring
!    ! the LAI and fVCover values returned from the subroutine are not modified.
!    LocLAI=MAX(LAI,0.02)
!    LocfVCover=MAX(fVCover,0.01)
!    !
!    ! determine the nearest array location for the desired LAI and fVCover
!    i=INT(LocLAI/DLAI+1)
!    j=INT(LocfVCover/DfVCover+1)
!    j=MIN(j,49)
!    !
!    ! interpolate RbC variable
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroVar(i,j)%RbC,                                 &
!         AeroVar(i+1,j)%RbC,                               &
!         AeroVar(i,j+1)%RbC,                               &
!         AeroVar(i+1,j+1)%RbC,                             &
!         RbC)


!    !
!    ! interpolate RdC variable
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroVar(i,j)%RdC,                                 &
!         AeroVar(i+1,j)%RdC,                               &
!         AeroVar(i,j+1)%RdC,                               &
!         AeroVar(i+1,j+1)%RdC,                             &
!         RdC)
!    !
!    ! interpolate roughness length'
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroVar(i,j)%zo,                                  &
!         AeroVar(i+1,j)%zo,                                &
!         AeroVar(i,j+1)%zo,                                &
!         AeroVar(i+1,j+1)%zo,                              &
!         zo)
!    !
!    ! interpolate zero plane displacement
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroVar(i,j)%zp_disp,                             &
!         AeroVar(i+1,j)%zp_disp,                           &
!         AeroVar(i,j+1)%zp_disp,                           &
!         AeroVar(i+1,j+1)%zp_disp,                         &
!         zp_disp)
!    !
!    RETURN
!  END SUBROUTINE AeroInterpolate
!  !
!  !=======================================================================
!  SUBROUTINE NewAeroInterpolate (LAI, fVCover,        &
!       LAIgrid,fVCovergrid,AeroTab,AeroVar)
!    !=======================================================================
!    ! calculates the aerodynamic parameters by bi-linear
!    ! interpolation from a lookup table of previously calculated values.
!    ! The interpolation table is a numpts x numpts LAI/fVCover grid with
!    ! LAI ranging from 0.02 to 10 and fVCover ranging from 0.01 to 1.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL LAI            ! actual area averaged LAI for interpolation
!    REAL fVCover        ! vegetation cover fraction for interpolation
!    REAL LAIgrid(50)    ! grid of LAI values for lookup table
!    REAL fVCovergrid(50)! grid of fVCover values for interpolation table
!    TYPE aero_var
!       REAL :: zo      ! Canopy roughness coeff
!       REAL :: zp_disp ! Zero plane displacement
!       REAL :: RbC     ! RB Coefficient (c1 or cc1)
!       REAL :: RdC     ! RC Coefficient (c2 or cc2)
!       REAL :: G2      ! Ratio Ra (actual) to Ra (log-linear) for momentum
!       REAL :: G3      ! Ratio Ra (actual) to Ra (log-linear) for heat transfer
!       REAL :: CORB1   ! Non-nuetral correction for Ra between Ha and z2
!       REAL :: CORB2   ! neutral value of RBB*U2^2 (RdC^2 for upper canopy)
!       REAL :: HA      ! Canopy source height
!    END TYPE aero_var
!    TYPE(aero_var) AeroTab(50,50)
!    !
!    ! begin output variables
!    TYPE(aero_var) AeroVar
!    !
!    ! begin internal variables
!    INTEGER i           ! index for LAI grid location
!    INTEGER j           ! index for fVCover grid location
!    REAL LocLAI         ! local LAI var. to prevent changing main LAI value
!    REAL LocfVCover     ! local fVCover var. to prevent changing fVCover value
!    REAL DLAI           ! grid spacing between LAI values in tables
!    REAL DfVCover       ! grid spacing between fVCover values in tables
!    !
!    ! calculate grid spacing (assumed fixed)
!    DLAI=LAIgrid(2)-LAIgrid(1)
!    DfVCover=fVCovergrid(2)-fVCovergrid(1)
!    !
!    ! Assign input LAI and fVCover to local variables and make sure
!    ! they lie within the limits of the interpolation tables, assuring
!    ! the LAI and fVCover values returned from the subroutine are not modified.
!    LocLAI=MAX(LAI,0.02)
!    LocfVCover=MAX(fVCover,0.01)
!    !
!    ! determine the nearest array location for the desired LAI and fVCover
!    i=INT(LocLAI/DLAI+1)
!    j=INT(LocfVCover/DfVCover+1)
!    j=MIN(j,49)
!    !
!    ! interpolate RbC variable
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%RbC,                                 &
!         AeroTab(i+1,j)%RbC,                               &
!         AeroTab(i,j+1)%RbC,                               &
!         AeroTab(i+1,j+1)%RbC,                             &
!         AeroVar%RbC)



!    !
!    ! interpolate RdC variable
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%RdC,                                 &
!         AeroTab(i+1,j)%RdC,                               &
!         AeroTab(i,j+1)%RdC,                               &
!         AeroTab(i+1,j+1)%RdC,                             &
!         AeroVar%RdC)
!    !
!    ! interpolate roughness length (z0)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%zo,                                  &
!         AeroTab(i+1,j)%zo,                                &
!         AeroTab(i,j+1)%zo,                                &
!         AeroTab(i+1,j+1)%zo,                              &
!         AeroVar%zo)
!    !
!    ! interpolate zero plane displacement (zp_disp)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%zp_disp,                             &
!         AeroTab(i+1,j)%zp_disp,                           &
!         AeroTab(i,j+1)%zp_disp,                           &
!         AeroTab(i+1,j+1)%zp_disp,                         &
!         AeroVar%zp_disp)
!    !
!    ! interpolate Ra ratio for momentum (G2)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%G2,                                  &
!         AeroTab(i+1,j)%G2,                                &
!         AeroTab(i,j+1)%G2,                                &
!         AeroTab(i+1,j+1)%G2,                              &
!         AeroVar%G2)
!    !
!    ! interpolate Ra ratio for Heat transfer (G3)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%G3,                                  &
!         AeroTab(i+1,j)%G3,                                &
!         AeroTab(i,j+1)%G3,                                &
!         AeroTab(i+1,j+1)%G3,                              &
!         AeroVar%G3)
!    !
!    ! interpolate Non-nuetral correction for Ra between Ha and z2 (Corb1)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%Corb1,                               &
!         AeroTab(i+1,j)%Corb1,                             &
!         AeroTab(i,j+1)%Corb1,                             &
!         AeroTab(i+1,j+1)%Corb1,                           &
!         AeroVar%Corb1)
!    !
!    ! interpolate neutral value of RBB*U2^2 (Corb2)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%Corb2,                               &
!         AeroTab(i+1,j)%Corb2,                             &
!         AeroTab(i,j+1)%Corb2,                             &
!         AeroTab(i+1,j+1)%Corb2,                           &
!         AeroVar%Corb2)
!    !
!    ! interpolate Canopy source height (HA)
!    CALL interpolate(                               &
!         LAIgrid(i),                                       &
!         LocLAI,                                           &
!         DLAI,                                             &
!         fVCovergrid(j),                                   &
!         LocfVCover,                                       &
!         DfVCover,                                         &
!         AeroTab(i,j)%HA,                                  &
!         AeroTab(i+1,j)%HA,                                &
!         AeroTab(i,j+1)%HA,                                &
!         AeroTab(i+1,j+1)%HA,                              &
!         AeroVar%HA)
!    !
!    RETURN
!  END SUBROUTINE NewAeroInterpolate
!  !
!  !=======================================================================
!  SUBROUTINE interpolate(x1, x, Dx,                  &
!       y1, y, Dy,                                        &
!       z11, z21, z12, z22, z)
!    !=======================================================================

!    IMPLICIT NONE

!    ! calculates the value of z=f(x,y) by linearly interpolating
!    ! between the 4 closest data points on a uniform grid.  The subroutine
!    ! requires a grid point (x1, y1), the grid spacing (Dx and Dy), and the
!    ! 4 closest data points (z11, z21, z12, and z22).
!    !
!    ! begin input variables
!    REAL x1  ! the x grid location of z11
!    REAL x   ! x-value at which you will interpolate z=f(x,y)
!    REAL Dx  ! grid spacing in the x direction
!    REAL y1  ! the y grid location of z11
!    REAL y   ! y-value at which you will interpolate z=f(x,y)
!    REAL Dy  ! grid spacing in the y direction
!    REAL z11 ! f(x1, y1)
!    REAL z21 ! f(x1+Dx, y1)
!    REAL z12 ! f(x1, y1+Dy)
!    REAL z22 ! f(x1+Dx, y1+Dy)
!    !
!    ! begin output variables
!    REAL z   ! f(x,y), the desired interpolated value
!    !
!    ! begin internal variables
!    REAL zp  ! z'=first interpolated value at (x, y1)
!    REAL zpp ! z''=second interpolated value at (x, Y1+Dy)
!    !
!    ! interpolate between z11 and z21 to calculate z' (zp) at (x, y1)
!    zp=z11+(x-x1)*(z21-z11)/Dx
!    !
!    ! interpolate between z12 and z22 to calculate z'' (zpp) at (x, Y1+Dy)
!    zpp=z12+(x-x1)*(z22-z12)/Dx
!    !
!    ! interpolate between zp and zpp to calculate z at (x,y)
!    z=zp+(y-y1)*(zpp-zp)/Dy
!    !
!    RETURN
!  END SUBROUTINE interpolate
!  !
!  !=======================================================================
!  SUBROUTINE srapar (ndvi, SRmin, SRmax, fPAR, fPARmax, fParmin)
!    !=======================================================================
!    ! calculates Canopy absorbed fraction of Photosynthetically
!    ! Active Radiation (fPAR) using the Simple Ratio (sr) method
!    ! (Los et al. (1998), eqn 6). This empirical method assumes a linear
!    ! relationship between fPAR and sr.
!    !----------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL ndvi    ! normalized difference vegetation index
!    REAL SRmin   ! minimum simple ratio for vegetation type
!    REAL SRmax   ! maximum simple ratio for vegetation type
!    REAL fPARmax ! Maximum possible FPAR corresponding to 98th percentile
!    REAL fPARmin ! Minimum possible FPAR corresponding to 2nd percentile
!    !
!    ! begin output variables
!    REAL fPAR    ! Canopy absorbed fraction of PAR
!    !
!    ! begin internal variables
!    REAL sr      ! simple ratio of near IR and visible radiances
!    !
!    ! Calculate simple ratio (SR)
!    sr=(1.+ndvi)/(1.-ndvi)
!    !
!    ! Insure calculated SR value falls within physical limits for veg. type
!    sr=MAX(sr,SRmin)
!    sr=MIN(sr,SRmax)
!    !
!    ! Calculate fPAR using SR method (Los et al. (1998), eqn 6)
!    fPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
!    !
!    RETURN
!  END SUBROUTINE srapar
!  !
!  !=======================================================================
!  SUBROUTINE NDVIapar (ndvi, NDVImin, NDVImax,       &
!       fPAR, fPARmax, fParmin)
!    !=======================================================================
!    ! calculates Canopy absorbed fraction of Photosynthetically
!    ! Active Radiation (fPAR) using the NDVI method
!    ! (Los et al. (1998), eqn 7). This empirical method assumes a linear
!    ! relationship between fPAR and NDVI.
!    !----------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL ndvi    ! normalized difference vegetation index
!    REAL NDVImin ! minimum NDVI for vegetation type
!    REAL NDVImax ! maximum NDVI for vegetation type
!    REAL fPARmax ! Maximum possible FPAR corresponding to 98th percentile
!    REAL fPARmin ! Minimum possible FPAR corresponding to 2nd percentile
!    !
!    ! begin output variables
!    REAL fPAR    ! Canopy absorbed fraction of PAR
!    !
!    ! begin internal variables
!    REAL LocNDVI ! local value of NDVI to prevent changes in input value
!    !
!    ! switch to local value of ndvi to prevent any changes going back to main
!    LocNDVI=NDVI
!    !
!    ! Insure calculated NDVI value falls within physical limits for veg. type
!    LocNDVI=MAX(LocNDVI,NDVImin)
!    LocNDVI=MIN(LocNDVI,NDVImax)
!    !
!    ! Calculate fPAR using NDVI method (Los et al. (1998), eqn 6)
!    fPAR=(LocNDVI-NDVImin)*(fPARmax-fPARmin)/          &
!         (NDVImax-NDVImin)+fPARmin
!    !
!    RETURN
!  END SUBROUTINE NDVIapar
!  !
!  !=======================================================================
!  SUBROUTINE AverageAPAR (ndvi, NDVImin, NDVImax, SRmin, SRmax,  &
!       fPARmax, fParmin, fPAR)
!    !=======================================================================
!    ! calculates Canopy absorbed fraction of Photosynthetically
!    ! Active Radiation (fPAR) using an average of the Simple Ratio (sr)
!    ! and NDVI methods (Los et al. (1999), eqn 5-6).  The empirical
!    ! SR method assumes a linear relationship between fPAR and SR.
!    ! The NDVI method assumes a linear relationship between fPAR and NDVI.
!    !----------------------------------------------------------------------
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL ndvi     ! normalized difference vegetation index
!    REAL NDVImin  ! minimum NDVI for vegetation type
!    REAL NDVImax  ! maximum NDVI for vegetation type
!    REAL SRmin    ! minimum NDVI for vegetation type
!    REAL SRmax    ! maximum NDVI for vegetation type
!    REAL fPARmax  ! Maximum possible FPAR corresponding to 98th percentile
!    REAL fPARmin  ! Minimum possible FPAR corresponding to 2nd percentile
!    !
!    ! begin output variables
!    REAL fPAR     ! Canopy absorbed fraction of PAR
!    !
!    ! begin internal variables
!    REAL LocNDVI  ! local value of NDVI to prevent changes in input value
!    REAL sr       ! simple ratio of near IR and visible radiances
!    REAL NDVIfPAR ! fPAR from NDVI method
!    REAL SRfPAR   ! fPAR from SR method
!    !
!    ! switch to local value of ndvi to prevent any changes going back to main
!    LocNDVI=NDVI
!    !
!    ! Insure calculated NDVI value falls within physical limits for veg. type
!    LocNDVI=MAX(LocNDVI,NDVImin)
!    LocNDVI=MIN(LocNDVI,NDVImax)
!  !  print*,'fpar1:',locndvi,ndvimin,ndvimax
!    !
!    ! Calculate simple ratio (SR)
!    sr=(1.+LocNDVI)/(1.-LocNDVI)

!    !
!    ! Calculate fPAR using SR method (Los et al. (1999), eqn 5)
!    SRfPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
!  !  print*,'fpar2:',sr,srmin,fparmax,fparmin,srmin,srmax,fparmin
!    !
!    ! Calculate fPAR using NDVI method (Los et al. (1999), eqn 6)
!    NDVIfPAR=(LocNDVI-NDVImin)*(fPARmax-fPARmin)/      &
!         (NDVImax-NDVImin)+fPARmin

!    !
!    ! take average of two methods
!    fPAR=0.5*(SRfPAR+NDVIfPAR)
!  !  print*,'fpar3:',fpar,srfpar,ndvifpar
!    !
!    RETURN
!  END SUBROUTINE AverageAPAR
!  !
!  !=======================================================================
!  SUBROUTINE aparnew (LAI,Green,LTran,LRef,gmudmu,fVCover,    &
!       fPAR, fPARmax, fPARmin)
!    !=======================================================================
!    ! recomputes the Canopy absorbed fraction of Photosynthetically
!    ! Active Radiation (fPAR), adjusting for solar zenith angle and the
!    ! vegetation cover fraction (fVCover) using a modified form of Beer's law..
!    ! See Sellers et al. Part II (1996), eqns. 9-13.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL LAI       ! Leaf Area Index
!    REAL Green     ! Greeness fraction of Leaf Area Index
!    REAL LTran(2,2)! Leaf transmittance for green/brown plants
!    REAL LRef(2,2) ! Leaf reflectance for green/brown plants
!    !                      For LTran and LRef:
!    !                        (1,1)=shortwave, green plants
!    !                        (2,1)=longwave, green plants
!    !                        (1,2)=shortwave, brown plants
!    !                        (2,2)=longwave, brown plants
!    REAL gmudmu    ! daily Time-mean canopy optical depth
!    REAL fVCover   ! Canopy cover fraction
!    REAL fPARmax   ! Maximum possible FPAR corresponding to 98th percentile
!    REAL fPARmin   ! Minimum possible FPAR corresponding to 2nd percentile
!    !
!    ! begin output variables
!    REAL fPAR      ! area average Canopy absorbed fraction of PAR
!    !
!    ! begin internal variables
!    REAL scatp     ! Canopy transmittance + reflectance coefficient wrt PAR
!    REAL PARk      ! mean canopy absorption optical depth wrt PAR
!    !
!    ! Calculate canopy transmittance + reflectance coefficient wrt PAR
!    ! transmittance + reflectance coefficient=green plants + brown plants
!    scatp=Green*(LTran(1,1)+LRef(1,1))+        &
!         (1.-Green)*(LTran(1,2)+LRef(1,2))
!    !
!    ! Calculate PAR absorption optical depth in canopy adjusting for
!    ! variance in projected leaf area wrt solar zenith angle
!    ! (Sellers et al. Part II (1996), eqn. 13b)
!    ! PAR absorption coefficient=(1-scatp)
!    PARk=SQRT(1.-scatp)*gmudmu
!    !
!    ! Calculate the new fPAR (Sellers et al. Part II (1996), eqn. 9)
!    fPAR=fVCover*(1.-EXP(-PARk*LAI/fVCover))
!    !
!    ! Ensure calculated fPAR falls within physical limits
!    fPAR=amax1(fPARmin,fPAR)
!    fPAR=amin1(fPARmax,fPAR)
!    !
!    RETURN
!  END SUBROUTINE aparnew
!  !
!  !******************************************************************
!  SUBROUTINE textclass(text)
!    !******************************************************************
!    ! Assigns soil texture classes based on the USDA texture triangle
!    ! using subroutines developed by aris gerakis
!    !******************************************************************
!    !* +-----------------------------------------------------------------------
!    !* |                         T R I A N G L E
!    !* | Main program that calls WHAT_TEXTURE, a function that classifies soil
!    !* | in the USDA textural triangle using sand and clay %
!    !* +-----------------------------------------------------------------------
!    !* | Created by: aris gerakis, apr. 98 with help from brian baer
!    !* | Modified by: aris gerakis, july 99: now all borderline cases are valid
!    !* | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
!    !* |              main program
!    !* +-----------------------------------------------------------------------
!    !* | COMMENTS
!    !* | Supply a data file with two columns, in free format: 1st column sand,
!    !* |   2nd column clay %, no header.  The output is a file with the classes.
!    !* +-----------------------------------------------------------------------
!    !* | You may use, distribute and modify this code provided you maintain
!    !* ! this header and give appropriate credit.
!    !* +-----------------------------------------------------------------------
!    !
!    ! Modifications:
!    !   Lara Prihodko customized triangle program for mapper (1/31/01)
!    !
!    IMPLICIT NONE

!    TYPE text_type
!       REAL :: clay  !  Percent clay content
!       REAL :: silt  !  Percent silt content
!       REAL :: sand  !  Percent sand content
!       INTEGER :: class  !  Soil texture class
!    END TYPE text_type

!    TYPE(text_type) :: text

!    INTEGER    :: texture, what_texture
!    REAL       :: sand, clay
!    REAL       :: silty_loam(1:7,1:2), sandy(1:7,1:2),              &
!         silty_clay_loam(1:7,1:2),                                      &
!         loam(1:7,1:2), clay_loam(1:7,1:2), sandy_loam(1:7,1:2),        &
!         silty_clay(1:7,1:2), sandy_clay_loam(1:7,1:2),                 &
!         loamy_sand(1:7,1:2), clayey(1:7,1:2), silt(1:7,1:2),           &
!         sandy_clay(1:7,1:2)
!    LOGICAL :: inpoly

!    !Initalize polygon coordinates:

!    DATA silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
!    DATA sandy/85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
!    DATA silty_clay_loam/0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 0,  &
!         0/
!    DATA loam/43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
!    DATA clay_loam/20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
!    DATA sandy_loam/50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 15,   &
!         0/
!    DATA silty_clay/0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
!    DATA sandy_clay_loam/52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35,     &
!         20, 0, 0/
!    DATA loamy_sand/70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
!    DATA clayey/20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/
!    DATA silt/0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
!    DATA sandy_clay/45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/


!    sand = 0.0
!    clay = 0.0

!    !Read input:

!    sand = text%sand
!    clay = text%clay

!    !Call function that estimates texture and put into structure:
!    text%class = what_texture (sand, clay, silty_loam, sandy,   &
!         silty_clay_loam,                                           &
!         loam, clay_loam, sandy_loam, silty_clay,                   &
!         sandy_clay_loam, loamy_sand, clayey, silt,                 &
!         sandy_clay)

!    RETURN
!  END SUBROUTINE textclass
!  !
!  !******************************************************************
!  !* +-----------------------------------------------------------------------
!  !* | WHAT TEXTURE?
!  !* | Function to classify a soil in the triangle based on sand and clay %
!  !* +-----------------------------------------------------------------------
!  !* | Created by: aris gerakis, apr. 98
!  !* | Modified by: aris gerakis, june 99.  Now check all polygons instead of
!  !* | stopping when a right solution is found.  This to cover all borderline
!  !* | cases.
!  !* +-----------------------------------------------------------------------

!  FUNCTION what_texture (sand, clay, silty_loam, sandy,       &
!       silty_clay_loam, loam, clay_loam,                          &
!       sandy_loam, silty_clay, sandy_clay_loam,                   &
!       loamy_sand, clayey, silt, sandy_clay)

!    IMPLICIT NONE

!    !Declare arguments:

!    REAL, INTENT(in) :: clay, sand, silty_loam(1:7,1:2),        &
!         sandy(1:7,1:2),                                            &
!         silty_clay_loam(1:7,1:2), loam(1:7,1:2),                   &
!         clay_loam(1:7,1:2), sandy_loam(1:7,1:2),                   &
!         silty_clay(1:7,1:2), sandy_clay_loam(1:7,1:2),             &
!         loamy_sand(1:7,1:2), clayey(1:7,1:2), silt(1:7,1:2),       &
!         sandy_clay(1:7,1:2)

!    !Declare local variables:

!    LOGICAL :: inpoly
!    INTEGER :: texture, what_texture

!    !Find polygon(s) where the point is.

!    texture = 13

!    IF (sand .GT. 0.0 .AND. clay .GT. 0.0) THEN

!       IF (inpoly(silty_loam, 6, sand, clay)) THEN
!          texture = 4
!       ENDIF
!       IF (inpoly(sandy, 3, sand, clay)) THEN
!          texture = 1
!       ENDIF
!       IF (inpoly(silty_clay_loam, 4, sand, clay)) THEN
!          texture = 10
!       ENDIF
!       IF (inpoly(loam, 5, sand, clay)) THEN
!          texture = 6
!       ENDIF
!       IF (inpoly(clay_loam, 4, sand, clay)) THEN
!          texture = 9
!       ENDIF
!       IF (inpoly(sandy_loam, 7, sand, clay)) THEN
!          texture = 3
!       ENDIF
!       IF (inpoly(silty_clay, 3, sand, clay)) THEN
!          texture = 11
!       ENDIF
!       IF (inpoly(sandy_clay_loam, 5, sand, clay)) THEN
!          texture = 7
!       ENDIF
!       IF (inpoly(loamy_sand, 4, sand, clay)) THEN
!          texture = 2
!       ENDIF
!       IF (inpoly(clayey, 5, sand, clay)) THEN
!          texture = 12
!       ENDIF
!       IF (inpoly(silt, 4, sand, clay)) THEN
!          texture = 5
!       ENDIF
!       IF (inpoly(sandy_clay, 3, sand, clay)) THEN
!          texture = 8
!       ENDIF

!    ENDIF

!    IF (sand == 100) THEN
!       texture = 1
!    ENDIF

!    IF (clay == 100) THEN
!       texture = 12
!    ENDIF

!    IF (texture == 13 ) THEN
!       texture = 13
!    ENDIF

!    what_texture = texture


!  END FUNCTION what_texture
!  !
!  !******************************************************************
!  !--------------------------------------------------------------------------
!  !                            INPOLY
!  !   Function to tell if a point is inside a polygon or not.
!  !--------------------------------------------------------------------------
!  !   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
!  !
!  !   Please feel free to use this source code for any purpose, commercial
!  !   or otherwise, as long as you don't restrict anyone else's use of
!  !   this source code.  Please give credit where credit is due.
!  !
!  !   Point-in-polygon algorithm, created especially for World-Wide Web
!  !   servers to process image maps with mouse-clickable regions.
!  !
!  !   Home for this file:  http://www.gcomm.com/develop/inpoly.c
!  !
!  !                                       6/19/95 - Bob Stein & Craig Yap
!  !                                       stein@gcomm.com
!  !                                       craig@cse.fau.edu
!  !--------------------------------------------------------------------------
!  !   Modified by:
!  !   Aris Gerakis, apr. 1998: 1.  translated to Fortran
!  !                            2.  made it work with real coordinates
!  !                            3.  now resolves the case where point falls
!  !                                on polygon border.
!  !   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
!  !   Aris Gerakis, july 1999: Now all borderline cases are valid
!  !--------------------------------------------------------------------------
!  !   Glossary:
!  !   function inpoly: true=inside, false=outside (is target point inside
!  !                    a 2D polygon?)
!  !   poly(*,2):  polygon points, [0]=x, [1]=y
!  !   npoints: number of points in polygon
!  !   xt: x (horizontal) of target point
!  !   yt: y (vertical) of target point
!  !--------------------------------------------------------------------------

!  FUNCTION inpoly (poly, npoints, xt, yt)

!    IMPLICIT NONE

!    !Declare arguments:

!    INTEGER :: npoints
!    REAL, INTENT(in)    :: poly(7, 2), xt, yt

!    !Declare local variables:

!    REAL    :: xnew, ynew, xold, yold, x1, y1, x2, y2
!    INTEGER :: i
!    LOGICAL :: inside, on_border, inpoly

!    inside = .FALSE.
!    on_border = .FALSE.

!    IF (npoints < 3)  THEN
!       inpoly = .FALSE.
!       RETURN
!    END IF

!    xold = poly(npoints,1)
!    yold = poly(npoints,2)

!    DO i = 1 , npoints
!       xnew = poly(i,1)
!       ynew = poly(i,2)

!       IF (xnew > xold)  THEN
!          x1 = xold
!          x2 = xnew
!          y1 = yold
!          y2 = ynew
!       ELSE
!          x1 = xnew
!          x2 = xold
!          y1 = ynew
!          y2 = yold
!       END IF

!       !The outer IF is the 'straddle' test and the 'vertical border' test.
!       !The inner IF is the 'non-vertical border' test and the 'north' test.

!       !The first statement checks whether a north pointing vector crosses
!       !(stradles) the straight segment.  There are two possibilities, depe-
!       !nding on whether xnew < xold or xnew > xold.  The '<' is because edge
!       !must be "open" at left, which is necessary to keep correct count when
!       !vector 'licks' a vertix of a polygon.

!       IF ((xnew < xt .AND. xt <= xold) .OR. (.NOT. xnew < xt .AND.    &
!            .NOT. xt <= xold)) THEN
!          !The test point lies on a non-vertical border:
!          IF ((yt-y1)*(x2-x1) == (y2-y1)*(xt-x1)) THEN
!             on_border = .TRUE.
!             !Check if segment is north of test point.  If yes, reverse the
!             !value of INSIDE.  The +0.001 was necessary to avoid errors due
!             !arithmetic (e.g., when clay = 98.87 and sand = 1.13):
!          ELSEIF ((yt-y1)*(x2-x1) < (y2-y1)*(xt-x1) + 0.001) THEN
!             inside = .NOT.inside ! cross a segment
!          ENDIF
!          !This is the rare case when test point falls on vertical border or
!          !left edge of non-vertical border. The left x-coordinate must be
!          !common.  The slope requirement must be met, but also point must be
!          !between the lower and upper y-coordinate of border segment.  There
!          !are two possibilities,  depending on whether ynew < yold or ynew >
!          !yold:
!       ELSEIF ((xnew == xt .OR. xold == xt) .AND. (yt-y1)*(x2-x1) ==     &
!            (y2-y1)*(xt-x1) .AND. ((ynew <= yt .AND. yt <= yold) .OR.        &
!            (.NOT. ynew < yt .AND. .NOT. yt < yold))) THEN
!          on_border = .TRUE.
!       ENDIF

!       xold = xnew
!       yold = ynew

!    ENDDO

!    !If test point is not on a border, the function result is the last state
!    !of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
!    !inside the polygon if it falls on any of its borders:

!    IF (.NOT. on_border) THEN
!       inpoly = inside
!    ELSE
!       inpoly = .TRUE.
!    ENDIF
!    !
!  END FUNCTION inpoly
!  !
!  !=======================================================================
!  SUBROUTINE FractionVegCover(NDVI,fPARmax,fPARmin,        &
!       SRmax,SRmin,fVCover)
!    !=======================================================================
!    ! calculates the vegetation cover fraction for a single pixel.
!    ! The maximum fPAR for pixel during entire year determines vegetation
!    ! cover fraction.  Maximum yearly NDVI corresponds to maximum fPAR.
!    ! Calculate fPAR from Simple Ratio using an empirical linear relationship..
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL NDVI      ! maximum FASIR NDVI value for a grid cell
!    REAL fPARmax   ! Maximum possible FPAR corresponding to 98th percentile
!    REAL fPARmin   ! Minimum possible FPAR corresponding to 2nd percentile
!    REAL SRmax     ! Maximum simple ratio for biome type
!    REAL SRmin     ! Minimum simple ratio for biome type
!    !
!    ! begin output variables
!    REAL fVCover   ! fractional vegetation cover
!    !
!    ! begin internal variables
!    REAL fPAR      ! maximum fPAR associated with maximum ndvi
!    !
!    ! The maximum fPAR for pixel during entire year determines vegetation
!    ! cover fraction.  Maximum yearly NDVI corresponds to maximum fPAR.
!    ! Calculate fPAR from Simple Ratio using an empirical linear relationship..
!    !
!    CALL srapar (NDVI,                                &
!         SRmin,                               &
!         SRmax,                               &
!         fPAR,                                &
!         fPARmax,                             &
!         fPARmin)
!    !
!    ! calculate fractional vegetation cover
!    fVCover=fPAR/fPARmax
!    !
!    RETURN
!  END SUBROUTINE FractionVegCover
!  !
!  !=======================================================================
!  SUBROUTINE SoilProperties(text,SoilVar)
!    !=======================================================================
!    ! calculates soil physical properties given sand and clay content
!    !
!    ! Modifications
!    !  Kevin Schaefer created subroutine for soil hydraulic properties (4/22/00)
!    !  Kevin Schaefer resp. variable curve fits from Raich et al., 1991 (6/19/00)
!    !  Kevin Schaefer combine code for hydraulic & respiration variables(3/30/01)
!    !
!    IMPLICIT NONE
!    !
!    ! begin Input variables
!    TYPE text_type
!       REAL :: clay  !  Percent clay content
!       REAL :: silt  !  Percent silt content
!       REAL :: sand  !  Percent sand content
!       INTEGER :: class  !  Soil texture class
!    END TYPE text_type

!    TYPE(text_type) :: text
!    !
!    ! begin Output soil property variables
!    INTEGER SoilNum   ! soil type number
!    TYPE soil_Physical
!       INTEGER SoilNum ! soil type number
!       REAL BEE     ! Wetness exponent for soil conductivity (-)
!       REAL PhiSat  ! Soil matrix potential (water tension) at Saturation (m)
!       REAL SatCo   ! Soil Hydraulic Conductivity at Saturation (m/s)
!       REAL poros   ! Soil Porosity or saturation water content (-)
!       REAL Slope   ! Cosine of mean slope
!       REAL Wopt    ! Optimal soil moisture for respiration (-)
!       REAL Skew    ! skewness exponent of soil respiration vs. wetness curve
!       REAL RespSat ! Parameter determining soil respiration at saturation (-)
!    END TYPE soil_Physical
!    !
!    TYPE(soil_Physical) SoilVar  ! time ind., soil dependant variables
!    !
!    ! begin local variables
!    REAL fclay   ! fraction of clay in soil
!    REAL fsand   ! fraction of sand in soil
!    !
!    ! calculate Soil hydraulic and thermal variables based on Klapp and
!    ! Hornberger
!    SoilVar%PhiSat=-10.*10**(1.88-0.0131*Text%Sand)/1000.
!    SoilVar%poros=0.489-0.00126*Text%Sand
!    SoilVar%SatCo=0.0070556*10**(-0.884+0.0153*Text%sand)/1000.
!    SoilVar%bee=2.91+0.159*Text%Clay
!    !
!    ! Calculate clay and sand fractions from percentages
!    fclay=Text%clay/100.
!    fsand=Text%sand/100.
!    !
!    ! Calculate soil respiration variables based on curve fits to
!    ! data shown in Raich et al. (1991)
!    SoilVar%Wopt=(-0.08*fclay**2+0.22*fclay+0.59)*100.
!    SoilVar%Skew=-2*fclay**3-0.4491*fclay**2+0.2101*fclay+0.3478
!    SoilVar%RespSat=0.25*fclay+0.5
!    !
!    ! assign value for mean slope of terrain
!    SoilVar%Slope=0.176
!    !
!    RETURN
!  END SUBROUTINE SoilProperties
!  !
!  !=======================================================================
!  SUBROUTINE gmuder (Lat, DOY, ChiL, gmudmu)
!    !=======================================================================
!    ! calculates daily time mean optical depth of canopy relative to the Sun.
!    !
!    IMPLICIT NONE
!    !
!    ! begin input variables
!    REAL Lat      ! latitude in degrees
!    REAL DOY      ! day-of-year (typically middle day of the month)
!    REAL ChiL     ! leaf angle distribution factor
!    !
!    ! begin output variables
!    REAL gmudmu   ! daily time mean canopy optical depth relative to Sun
!    REAL test     ! test variable
!    !
!    ! begin internal variables
!    REAL mumax    ! max cosine of the Solar zenith angle (noon)
!    REAL mumin    ! min cosine of the Solar zenith angle (rise/set)
!    REAL dec      ! declination of the Sun (Solar Declination)
!    REAL pi180    ! conversion factor from degrees to radians
!    REAL aa       ! minimum possible LAI projection vs. cosine Sun angle
!    REAL bb       ! slope leaf area projection vs. cosine Sun angle
!    !
!    ! Calculate conversion factor from degrees to radians
!    pi180=3.14159/180.
!    !
!    ! Calculate solar declination in degrees
!    dec=23.5*SIN(1.72e-2*(DOY-80))
!    !
!    ! Calculate maximum cosine of zenith angle corresponding to noon
!    mumax=COS((dec-lat)*pi180)
!    mumax=MAX(0.02, mumax)
!    !
!    ! Assign min cosine zenith angle corresponding to start disc set (cos(89.4))
!    mumin=0.01
!    !
!    ! The projected leaf area relative to the Sun is G(mu)=aa+bb*mu
!    ! Calculate minimum projected leaf area
!    aa=0.5-0.633*ChiL-0.33*ChiL*ChiL
!    !
!    ! Calculate slope of projected leaf area wrt cosine sun angle
!    bb=0.877*(1.-2.*aa)
!    !
!    ! Calculate mean optical depth of canopy by integrating G(mu)/mu over
!    ! all values of mu.  Since G(mu) has an analytical form, this comes to
!    gmudmu=aa*alog(mumax/mumin)/(mumax-mumin)+bb
!    !
!    RETURN
!  END SUBROUTINE gmuder
!  !
!  !-----------------------------------------------------------------

!  !
!  !=================SUBROUTINE RADA2======================================
!  !
!  SUBROUTINE RADA2(snoww,zlt,z1,z2                            &
!       ,          asnow,tg,sunang,tf,ref,tran,chil                 &
!       ,          green,vcover,soref,radfac,salb,thermk,tgeff4     &
!       ,          tc,len)

!    IMPLICIT NONE

!    !=======================================================================
!    !
!    !     CALCULATION OF ALBEDOS VIA TWO STREAM APPROXIMATION( DIRECT
!    !     AND DIFFUSE ) AND PARTITION OF RADIANT ENERGY
!    !
!    !-----------------------------------------------------------------------


!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       SALB(2,2)      SURFACE ALBEDOS
!    !       TGEFF4         EFFECTIVE SURFACE RADIATIVE TEMPERATURE (K)
!    !       RADFAC(2,2,2)  RADIATION ABSORPTION FACTORS
!    !       THERMK         CANOPY GAP FRACTION FOR TIR RADIATION
!    !
!    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!    !
!    !       ALBEDO(2,2,2)  COMPONENT REFLECTANCES
!    !
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    !
!    !     arg list declarations
!    !
!    INTEGER len
!    REAL snoww(len,2),                                   &
!         ref(len,2,2),tran(len,2,2),soref(len,2),        &
!         salb(len,2,2),radfac(len,2,2,2), asnow,         &
!         zlt(len),z1(len),z2(len),tg(len),sunang(len),   &
!         chil(len),green(len),vcover(len),               &
!         tgeff4(len),tc(len), thermk(len)
!    REAL tf

!    !     local variables
!    REAL TRANC1(2), TRANC2(2), TRANC3(2), satcap(len,2),             &
!         f(len),fmelt(len),zmew(len),albedo(len,2,2,2),              &
!         canex(len), areas(len), facs, scov, scat, chiv, aa, bb,     &
!         fac2, fac1, zkat, tg4, tc4, tgs, hh6,                       &
!         hh5, hh3, hh2, zmk, hh10, hh9, hh8, hh7, den, zp, f1, ge,   &
!         ek, epsi, power2, power1, zat, psi, hh4, hh1, fe, de, bot,  &
!         ce, be, betao, upscat, acss, extkb, proj, tran2, tran1,     &
!         reff1, reff2
!    INTEGER iwave, i, irad
!    !
!    !----------------------------------------------------------------------
!    !
!    !
!    !     MODIFICATION FOR EFFECT OF SNOW ON UPPER STOREY ALBEDO
!    !         SNOW REFLECTANCE   = 0.80, 0.40 . MULTIPLY BY 0.6 IF MELTING
!    !         SNOW TRANSMITTANCE = 0.20, 0.54
!    !
!    !
!    !-----------------------------------------------------------------------
!    !
!    !
!    DO i = 1,len                      ! loop over gridpoint
!       !        this portion is snow1 inlined
!       CANEX(i)  = 1.-( SNOWw(i,1)*5.-Z1(i))/         &
!            (Z2(i)-Z1(i))
!       CANEX(i)  = MAX( 0.1, CANEX(i) )
!       CANEX(i)  = MIN( 1.0, CANEX(i) )
!       AREAS(i)    = MIN(1.0 , ASNOW*SNOWw(i,2))
!       SATCAP(i,1) = ZLT(i)*0.0001 * CANEX(i)
!       !
!       !c Collatz-Bounoua change satcap(2) to 0.0002
!       !
!       !bl          SATCAP(i,2) = 0.002
!       SATCAP(i,2) = 0.0002             ! lahouari
!       !    end old snow1
!       F(i) = MAX(0.01746,SUNANG(i))
!       FACS  = ( TG(i)-TF ) * 0.04
!       FACS  = MAX( 0.0 , FACS)
!       FACS  = MIN( 0.4, FACS)
!       FMELT(i) = 1. - FACS
!    ENDDO
!    !
!    !-----------------------------------------------------------------------
!    !
!    DO IWAVE = 1, 2  !DO 1000 IWAVE = 1, 2
!       !
!       !----------------------------------------------------------------------
!       DO i = 1,len                      ! loop over gridpoint
!          SCOV =  MIN( 0.5, SNOWw(i,1)/SATCAP(i,1) )
!          REFF1 = ( 1. - SCOV ) * REF(i,IWAVE,1) + SCOV * ( 1.2 -           &
!               IWAVE * 0.4 ) * FMELT(i)
!          REFF2 = ( 1. - SCOV ) * REF(i,IWAVE,2) + SCOV * ( 1.2 -           &
!               IWAVE * 0.4 ) * FMELT(i)
!          TRAN1 = TRAN(i,IWAVE,1) * ( 1. - SCOV )                           &
!               + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT(i) )         &
!               * TRAN(i,IWAVE,1)
!          TRAN2 = TRAN(i,IWAVE,2) * ( 1. - SCOV )                           &
!               + SCOV * ( 1.- ( 1.2 - IWAVE * 0.4 ) * FMELT(i) ) * 0.9   &
!               * TRAN(i,IWAVE,2)
!          !
!          !----------------------------------------------------------------------
!          !
!          !     CALCULATE AVERAGE SCATTERING COEFFICIENT, LEAF PROJECTION AND
!          !     OTHER COEFFICIENTS FOR TWO-STREAM MODEL.
!          !
!          !      SCAT  (OMEGA)        : EQUATION (1,2) , SE-85
!          !      PROJ  (G(MU))        : EQUATION (13)  , SE-85
!          !      EXTKB (K, G(MU)/MU)  : EQUATION (1,2) , SE-85
!          !      ZMEW  (INT(MU/G(MU)) : EQUATION (1,2) , SE-85
!          !      ACSS  (A-S(MU))      : EQUATION (5)   , SE-85
!          !      EXTK  (K, VARIOUS)   : EQUATION (13)  , SE-85
!          !      UPSCAT(OMEGA*BETA)   : EQUATION (3)   , SE-85
!          !      BETAO (BETA-0)       : EQUATION (4)   , SE-85
!          !
!          !----------------------------------------------------------------------
!          !
!          SCAT = GREEN(i)*( TRAN1 + REFF1 ) +( 1.-GREEN(i) ) *     &
!               ( TRAN2 + REFF2)
!          CHIV = CHIL(i)
!          !
!          IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01
!          AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV
!          BB = 0.877 * ( 1. - 2. * AA )
!          !
!          PROJ = AA + BB * F(i)
!          EXTKB = ( AA + BB * F(i) ) / F(i)
!          ZMEW(i) = 1. / BB * ( 1. - AA / BB                        &
!               * LOG ( ( AA + BB ) / AA ) )
!          ACSS = SCAT / 2. * PROJ / ( PROJ + F(i) * BB )
!          ACSS = ACSS * ( 1. - F(i) * AA                               &
!               / ( PROJ + F(i) * BB ) * LOG ( ( PROJ                &
!               +   F(i) * BB + F(i) * AA ) / ( F(i) * AA ) ) )
!          !
!          UPSCAT = GREEN(i) * TRAN1 + ( 1.-GREEN(i) ) * TRAN2
!          UPSCAT = 0.5 * ( SCAT + ( SCAT - 2. * UPSCAT ) *             &
!               (( 1. - CHIV ) / 2. ) ** 2 )
!          BETAO = ( 1. + ZMEW(i) * EXTKB )                             &
!               / ( SCAT * ZMEW(i) * EXTKB ) * ACSS
!          !
!          !----------------------------------------------------------------------
!          !
!          !     Intermediate variables identified in appendix of SE-85.
!          !
!          !      BE          (B)     : APPENDIX      , SE-85
!          !      CE          (C)     : APPENDIX      , SE-85
!          !      BOT         (SIGMA) : APPENDIX      , SE-85
!          !      HH1         (H1)    : APPENDIX      , SE-85
!          !      HH2         (H2)    : APPENDIX      , SE-85
!          !      HH3         (H3)    : APPENDIX      , SE-85
!          !      HH4         (H4)    : APPENDIX      , SE-85
!          !      HH5         (H5)    : APPENDIX      , SE-85
!          !      HH6         (H6)    : APPENDIX      , SE-85
!          !      HH7         (H7)    : APPENDIX      , SE-85
!          !      HH8         (H8)    : APPENDIX      , SE-85
!          !      HH9         (H9)    : APPENDIX      , SE-85
!          !      HH10        (H10)   : APPENDIX      , SE-85
!          !      PSI         (H)     : APPENDIX      , SE-85
!          !      ZAT         (L-T)   : APPENDIX      , SE-85
!          !      EPSI        (S1)    : APPENDIX      , SE-85
!          !      EK          (S2)    : APPENDIX      , SE-85
!          !----------------------------------------------------------------------
!          !
!          BE = 1. - SCAT + UPSCAT
!          CE = UPSCAT
!          BOT = ( ZMEW(i) * EXTKB ) ** 2 + ( CE**2 - BE**2 )
!          IF ( ABS(BOT) .LE. 1.E-10) THEN
!             SCAT = SCAT* 0.98
!             BE = 1. - SCAT + UPSCAT
!             BOT = ( ZMEW(i) * EXTKB ) ** 2 + ( CE**2 - BE**2 )
!          END IF
!          DE = SCAT * ZMEW(i) * EXTKB * BETAO
!          FE = SCAT * ZMEW(i) * EXTKB * ( 1. - BETAO )
!          HH1 = -DE * BE + ZMEW(i) * DE * EXTKB - CE * FE
!          HH4 = -BE * FE - ZMEW(i) * FE * EXTKB - CE * DE
!          !
!          PSI = SQRT(BE**2 - CE**2)/ZMEW(i)
!          !
!          ZAT = ZLT(i)/VCOVER(i)*CANEX(i)
!          !
!          POWER1 = MIN( PSI*ZAT, 50.E0 )
!          POWER2 = MIN( EXTKB*ZAT, 50.E0 )
!          EPSI = EXP( - POWER1 )
!          EK = EXP ( - POWER2 )
!          !
!          ALBEDO(i,2,IWAVE,1) = SOREF(i,IWAVE)*(1.-AREAS(i))            &
!               + ( 1.2-IWAVE*0.4 )*FMELT(i) * AREAS(i)
!          ALBEDO(i,2,IWAVE,2) = SOREF(i,IWAVE)*(1.-AREAS(i))            &
!               + ( 1.2-IWAVE*0.4 )*FMELT(i) * AREAS(i)
!          GE = ALBEDO(i,2,IWAVE,1)/ALBEDO(i,2,IWAVE,2)
!          !
!          !----------------------------------------------------------------------
!          !     CALCULATION OF DIFFUSE ALBEDOS
!          !
!          !      ALBEDO(1,IWAVE,2) ( I-UP ) : APPENDIX , SE-85
!          !----------------------------------------------------------------------
!          !
!          F1 = BE - CE / ALBEDO(i,2,IWAVE,2)
!          ZP = ZMEW(i) * PSI
!          !
!          DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -                    &
!               ( BE - ZP ) * ( F1 + ZP ) * EPSI
!          HH7 = CE * ( F1 - ZP ) / EPSI / DEN
!          HH8 = -CE * ( F1 + ZP ) * EPSI / DEN
!          F1 = BE - CE * ALBEDO(i,2,IWAVE,2)
!          DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI
!          !
!          HH9 = ( F1 + ZP ) / EPSI / DEN
!          HH10 = - ( F1 - ZP ) * EPSI / DEN
!          TRANC2(IWAVE) = HH9 * EPSI + HH10 / EPSI
!          !
!          ALBEDO(i,1,IWAVE,2) =  HH7 + HH8
!          !
!          !----------------------------------------------------------------------
!          !     CALCULATION OF DIRECT ALBEDOS AND CANOPY TRANSMITTANCES.
!          !
!          !      ALBEDO(1,IWAVE,1) ( I-UP )   : EQUATION(11)   , SE-85
!          !      TRANC(IWAVE)      ( I-DOWN ) : EQUATION(10)   , SE-85
!          !
!          !----------------------------------------------------------------------
!          !
!          F1 = BE - CE / ALBEDO(i,2,IWAVE,2)
!          ZMK = ZMEW(i) * EXTKB
!          !
!          DEN = ( BE + ZP ) * ( F1 - ZP ) / EPSI -           &
!               ( BE - ZP ) * ( F1 + ZP ) * EPSI
!          HH2 = ( DE - HH1 / BOT * ( BE + ZMK ) )              &
!               * ( F1 - ZP ) / EPSI -                        &
!               ( BE - ZP ) * ( DE - CE*GE - HH1 / BOT      &
!               * ( F1 + ZMK ) ) * EK
!          HH2 = HH2 / DEN
!          HH3 = ( BE + ZP ) * (DE - CE*GE -                    &
!               HH1 / BOT * ( F1 + ZMK ))* EK -               &
!               ( DE - HH1 / BOT * ( BE + ZMK ) ) *           &
!               ( F1 + ZP ) * EPSI
!          HH3 = HH3 / DEN
!          F1 = BE - CE * ALBEDO(i,2,IWAVE,2)
!          DEN = ( F1 + ZP ) / EPSI - ( F1 - ZP ) * EPSI
!          HH5 = - HH4 / BOT * ( F1 + ZP ) / EPSI -            &
!               ( FE + CE*GE*ALBEDO(i,2,IWAVE,2) +            &
!               HH4 / BOT*( ZMK-F1 ) ) * EK
!          HH5 = HH5 / DEN
!          HH6 =   HH4 / BOT * ( F1 - ZP ) * EPSI +            &
!               ( FE + CE*GE*ALBEDO(i,2,IWAVE,2) +            &
!               HH4 / BOT*( ZMK-F1 ) ) * EK
!          HH6 = HH6 / DEN
!          TRANC1(IWAVE) = EK
!          TRANC3(IWAVE) = HH4 / BOT * EK + HH5 * EPSI + HH6 / EPSI
!          !
!          ALBEDO(i,1,IWAVE,1) = HH1 / BOT + HH2 + HH3
!          !
!          !----------------------------------------------------------------------
!          !
!          !
!          !----------------------------------------------------------------------
!          !     CALCULATION OF TERMS WHICH MULTIPLY INCOMING SHORT WAVE FLUXES
!          !     TO GIVE ABSORPTION OF RADIATION BY CANOPY AND GROUND
!          !
!          !      RADFAC      (F(IL,IMU,IV)) : EQUATION (19,20) , SE-86
!          !----------------------------------------------------------------------
!          !
!          RADFAC(i,2,IWAVE,1) = ( 1.-VCOVER(i) )                   &
!               * ( 1.-ALBEDO(i,2,IWAVE,1) ) + VCOVER(i)            &
!               * ( TRANC1(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,1) )      &
!               + TRANC3(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) ) )
!          !
!          RADFAC(i,2,IWAVE,2) = ( 1.-VCOVER(i) )                    &
!               * ( 1.-ALBEDO(i,2,IWAVE,2) ) + VCOVER(i)             &
!               *  TRANC2(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) )
!          !
!          RADFAC(i,1,IWAVE,1) = VCOVER(i)                           &
!               * ( ( 1.-ALBEDO(i,1,IWAVE,1) )                       &
!               - TRANC1(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,1) )         &
!               - TRANC3(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) ) )
!          !
!          RADFAC(i,1,IWAVE,2) = VCOVER(i)                           &
!               * ( ( 1.-ALBEDO(i,1,IWAVE,2) )                     &
!               - TRANC2(IWAVE) * ( 1.-ALBEDO(i,2,IWAVE,2) ) )
!       ENDDO
!       !
!       !----------------------------------------------------------------------
!       !     CALCULATION OF TOTAL SURFACE ALBEDOS ( SALB ) WITH WEIGHTING
!       !     FOR COVER FRACTIONS.
!       !----------------------------------------------------------------------
!       !
!       DO IRAD = 1, 2 !DO 3000 IRAD = 1, 2
!          DO i = 1,len        !  loop over gridpoint
!             SALB(i,IWAVE,IRAD) = ( 1.-VCOVER(i) )               &
!                  * ALBEDO(i,2,IWAVE,IRAD) +             &
!                  VCOVER(i) * ALBEDO(i,1,IWAVE,IRAD)
!          ENDDO
!       ENDDO  !3000 CONTINUE
!       !
!       !----------------------------------------------------------------------
!       !
!    ENDDO !1000 CONTINUE
!    !
!    !----------------------------------------------------------------------
!    !
!    !     CALCULATION OF LONG-WAVE FLUX TERMS FROM CANOPY AND GROUND
!    !
!    !----------------------------------------------------------------------
!    !
!    DO i = 1,len                  !  loop over gridpoint
!       TGS = MIN(TF,TG(i))*AREAS(i)            &
!            + TG(i)*(1.-AREAS(i))
!       TC4 = TC(i)**4
!       TG4 = TGS**4
!       !
!       ZKAT = 1./ZMEW(i) * ZLT(i) / VCOVER(i)
!       ZKAT = MIN( 50.E0 , ZKAT )
!       ZKAT = MAX( 1.E-5, ZKAT )
!       THERMK(i) = EXP(-ZKAT)
!       !
!       FAC1 =  VCOVER(i) * ( 1.-THERMK(i) )
!       FAC2 =  1.
!       TGEFF4(i) =  FAC1 * TC4                         &
!            + (1. - FAC1 ) * FAC2 * TG4
!    ENDDO

!    RETURN
!  END SUBROUTINE RADA2

!  !---------------------------------------------------------------

!  !
!  !==================SUBROUTINE RNLOAD====================================
!  !
!  SUBROUTINE rnload(len, nsib,                            &
!       radvbc,radvdc,radnbc,radndc,dlwbot,VCOVER,   &
!       thermk,radfac,radn,radc3)

!    IMPLICIT NONE
!    !
!    !=======================================================================
!    !
!    !    calculation of absorption of radiation by surface.  Note that
!    !       output from this calculation (radc3) only accounts for the
!    !       absorption of incident longwave and shortwave fluxes.  The
!    !       total net radiation calculation is performed in subroutine
!    !       netrad.
!    !
!    !=======================================================================
!    !

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       RADN(2,3)      INCIDENT RADIATION FLUXES (W M-2)
!    !       RADC3(2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len, nsib
!    REAL  radc3(len,2), radn(len,2,2), radfac(nsib,2,2,2),       &
!         radvbc(len), radvdc(len), radnbc(len), radndc(len),    &
!         dlwbot(len), VCOVER(len), thermk(len)

!    INTEGER i, iveg, iwave, irad

!    !-----------------------------------------------------------------------
!    !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
!    !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
!    !     SOURCE FOR TRANSPIRATION.
!    !
!    !      RADN        (F(IW,IMU,O)) : EQUATION (19-22) , SE-86
!    !      RADC3       (FC,FGS)      : EQUATION (21,22) , SE-86
!    !-----------------------------------------------------------------------

!    DO i = 1,len
!       RADc3(i,1) = 0.
!       RADc3(i,2) = 0.
!       radn(i,1,1) = radvbc(i)
!       radn(i,1,2) = radvdc(i)
!       radn(i,2,1) = radnbc(i)
!       radn(i,2,2) = radndc(i)
!    ENDDO

!    DO iveg=1,2
!       DO iwave=1,2
!          DO irad=1,2
!             DO i = 1,len
!                radc3(i,iveg) = radc3(i,iveg) +         &
!                     radfac(i,iveg,iwave,irad) *          &
!                     radn(i,iwave,irad)
!             ENDDO
!          ENDDO
!       ENDDO
!    ENDDO

!    !     absorbed downwelling radiation only

!    DO i = 1,len
!       RADc3(i,1) = RADc3(i,1) + dlwbot(i) *            &
!            VCOVER(i) * (1.- THERMK(i))
!       RADc3(i,2) = RADc3(i,2) + dlwbot(i) *            &
!            (1.-VCOVER(i) * (1.-THERMK(i)))
!    ENDDO

!    RETURN
!  END SUBROUTINE rnload

!  ! Function E - Defined by Alvaro
!  REAL FUNCTION E(X)
!    REAL, INTENT(IN) :: X

!    E = EXP( 21.18123 - 5418. / X ) / .622
!  END FUNCTION E

!  ! Function GE - Defined by Alvaro
!  REAL FUNCTION GE(X)
!    REAL, INTENT(IN) :: X

!    GE = EXP( 21.18123 - 5418. / X ) * 5418. / (X*X) / .622
!  END FUNCTION GE


!  !---------------------------------------------------------
!  !
!  !================SUBROUTINE BEGTEM=======================================
!  !
!  SUBROUTINE begtem(tc,tg,cpair,hlat,psur,snomel     &
!       ,               zlt,clai,cw,www,poros,pie, psy     &
!       ,               phsat,bee,ccx,cg, phc              &
!       ,               tgs,etc,etgs,getc,getgs,rstfac     &
!       ,               rsoil,hr,wc,wg,snoww,capac         &
!       ,               areas,satcap,csoil,tf,g            &
!       ,               snofac,len,nlen, forcerestore )

!    IMPLICIT NONE

!    !========================================================================
!    !
!    !     Calculation of flux potentials and constants prior to heat
!    !         flux calculations.  Corresponds to first half of TEMREC
!    !         in 1D model.
!    !
!    !========================================================================


!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       RSTFAC(2)      SOIL MOISTURE STRESS FACTOR
!    !       RSOIL          SOIL SURFACE RESISTANCE (S M-1)
!    !       HR             SOIL SURFACE RELATIVE HUMIDITY
!    !       WC             CANOPY WETNESS FRACTION
!    !       WG             GROUND WETNESS FRACTION
!    !       CCX            CANOPY HEAT CAPACITY (J M-2 K-1)
!    !       CG             GROUND HEAT CAPACITY (J M-2 K-1)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len, nlen
!    LOGICAL forcerestore
!    REAL                                                    &
!         www(len,3),rstfac(len,4),satcap(len,2),            &
!         snoww(nlen,2),capac(nlen,2),                       &
!         tc(nlen),tg(nlen),tm(nlen),hlat,psur(nlen),        &
!         zlt(nlen),poros(nlen),phsat(nlen),                 &
!         bee(nlen),ccx(len),cg(nlen),tgs(len),etc(len),     &
!         etgs(len),getgs(len),rsoil(len),                   &
!         hr(len),wc(len),wg(len),areas(len),csoil(len),     &
!         getc(len), phc(nlen), psy(len)
!    REAL                                                    &
!         cpair, snofac, snomel, cw, tf, g, clai, pie
!    !     local variables
!    INTEGER i
!    REAL                                                          &
!         e, ge, slamda, shcap, x, tsnow, rsnow, fac, psit, argg   &
!         ,    phroot(len), shcapf, shcapu, one

!    !----------------------------------------------------------------------
!    !     E(X) IS VAPOUR PRESSURE IN MBARS AS A FUNCTION OF TEMPERATURE
!    !     GE(X) IS D E(X) / D ( TEMP )
!    !----------------------------------------------------------------------

!    !E(X) = EXP( 21.18123 - 5418. / X ) / .622
!    !GE(X) = EXP( 21.18123 - 5418. / X ) * 5418.              &
!    !     / (X*X) / .622
!    ! Functions defined above.

!    one = 1.0

!    !----------------------------------------------------------------------

!    SNOFAC   = HLAT/ ( HLAT + SNOMEL * 1.E-3 )
!    !----------------------------------------------------------------------
!    !     CALCULATION OF CANOPY AND GROUND HEAT CAPACITIES.
!    !     N.B. THIS SPECIFICATION DOES NOT NECESSARILY CONSERVE ENERGY WHEN
!    !     DEALING WITH VERY LARGE SNOWPACKS.
!    !     HEAT CAPACITY OF THE SOIL, AS USED IN FORCE-RESTORE HEAT FLUX
!    !     DESCRIPTION. DEPENDENCE OF CSOIL ON POROSITY AND WETNESS.
!    !
!    !      CG          (cg)    : EQUATION (?) , CS-81
!    !----------------------------------------------------------------------
!    IF(forcerestore) THEN
!       DO i = 1,len
!          SLAMDA = ( 1.5*(1.-POROS(i)) +                         &
!               1.3*WWW(i,1)*POROS(i) ) /                    &
!               ( 0.75 + 0.65*POROS(i) -                      &
!               0.4*WWW(i,1)*POROS(i) )                      &
!               * 0.4186
!          SHCAP  = ( 0.5*(1.-POROS(i)) + WWW(i,1)*               &
!               POROS(i)) * 4.186E6
!          CSOIL(i)  = SQRT( SLAMDA * SHCAP * 86400./PIE ) * 0.5

!          !950511 adjust for different heat capacity of snow
!          !       CCX(i) = ZLT(i) * CLAI +                             &
!          !                  (0.5*SNOWw(i,1)+CAPAC(i,1))*CW
!          CCX(i) = ZLT(i)*CLAI+                               &
!               (SNOWw(i,1)+CAPAC(i,1))*CW
!          !950511 adjust for different heat capacity of snow
!          !       CG(i)  = CSOIL(i) +                                  &
!          !            MIN ( 0.025*one, (0.5 *SNOWw(i,2)+CAPAC(i,2))) * CW
!          CG(i)  = CSOIL(i) +                                  &
!               MIN ( 0.05*one, (SNOWw(i,2)+CAPAC(i,2))) * CW
!       ENDDO
!    ELSE
!       DO i = 1,len
!          !
!          !--------------------
!          ! new calculation for ground heat capacity cg - no longer force-restore
!          ! now for the top 1cm of snow and soil, with phase change in the soil
!          ! incorporated into the heat capacity from +0.5C to -0.5C

!          SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,1)*              &
!               POROS(i)) * 4.186E6
!          shcapf =  0.5 * (1. + poros(i) * (www(i,1)-1.0)) * 4.186E6
!          IF(tg(i).GE.tf+0.5) THEN
!             csoil(i) = shcapu * 0.02
!          ELSE IF(tg(i).LE.tf-0.5) THEN
!             csoil(i) = shcapf * 0.02
!          ELSE
!             csoil(i) = (0.5*(shcapu+shcapf) +                     &
!                  snomel*poros(i)*www(i,1) ) * 0.02
!          ENDIF

!          CCX(i) = ZLT(i)*CLAI+                               &
!               (0.5*SNOWw(i,1)+CAPAC(i,1))*CW
!          CG(i)  = (1.-areas(i))*CSOIL(i) +                    &
!               cw * (capac(i,2) + 0.01 * areas(i))
!          !
!       ENDDO
!    ENDIF
!    DO i = 1,len
!       ! HLAT(i)  = ( 3150.19 - 2.378 * TM(i) ) * 1000. !use constant passed in
!       PSY(i)      = CPAIR / HLAT * PSUR(i) / .622
!       !
!       !----------------------------------------------------------------------
!       !      Calculation of ground surface temperature and wetness fractions
!       !
!       !----------------------------------------------------------------------
!       !
!       TSNOW = MIN ( TF-0.01, TG(i) )
!       RSNOW = SNOWw(i,2) /                             &
!            (SNOWw(i,2)+CAPAC(i,2)+1.E-10)
!       !
!       TGS(i) = TSNOW*AREAS(i) + TG(i)*(1.-AREAS(i))
!       IF(tgs(i).LT.0.0)tgs(i) = SQRT(tgs(i))
!       !
!       ETC(i)   = E(TC(i))
!       ETGS(i)  = E(TGS(i))
!       GETC(i)  = GE(TC(i))
!       GETGS(i) = GE(TGS(i))
!       !
!       WC(i) = MIN( one,( CAPAC(i,1) +                          &
!            SNOWw(i,1))/SATCAP(i,1) )
!       WG(i) = MAX( 0.*one,  CAPAC(i,2)/SATCAP(i,2) )*0.25

!       !-----------------------------------------------------------------------
!       !     CALCULATION OF SOIL MOISTURE STRESS FACTOR.
!       !     AVERAGE SOIL MOISTURE POTENTIAL IN ROOT ZONE (LAYER-2) USED AS
!       !     SOURCE FOR TRANSPIRATION.
!       !
!       !      PHROOT      (PSI-R) : EQUATION (48) , SE-86
!       !      RSTFAC(2)  F(PSI-L) : MODIFICATION OF EQUATION (12), SE-89
!       !-----------------------------------------------------------------------
!       !
!       PHROOT(i) = PHSAT(i) *                                &
!            MAX( 0.02*one, WWW(i,2) ) ** ( - BEE(i) )
!       PHROOT(i) = MAX ( PHROOT(i), -2.E3*one )
!       RSTFAC(i,2) = 1./( 1. + EXP( 0.02*( PHC(i)-PHROOT(i)) ))
!       RSTFAC(i,2) = MAX( 0.0001*one, RSTFAC(i,2) )
!       RSTFAC(i,2) = MIN( one,     RSTFAC(i,2) )
!       !
!       !----------------------------------------------------------------------
!       !
!       !      RSOIL FUNCTION FROM FIT TO FIFE-87 DATA.  Soil surface layer
!       !         relative humidity.
!       !
!       !      RSOIL      (RSOIL) : HEISER 1992 (PERSONAL COMMUNICATION)
!       !      HR         (Fh)    : EQUATION (66) , SE-86
!       !----------------------------------------------------------------------
!       !
!       FAC = MIN( WWW(i,1), one )
!       FAC = MAX( FAC, 0.02*one  )
!       !
!       ! Collatz-Bounoua change rsoil to FIFE rsoil formulation from eq(19) SE-92
!       !
!       !cbl     RSOIL(i) =  MAX (0.1*one, 694. - FAC*1500.) + 23.6
!       rsoil(i) =  EXP(8.206 - 4.255 * fac)      ! lahouari

!       !
!       PSIT = PHSAT(i) * FAC ** (- BEE(i) )
!       ARGG = MAX(-10.*one,(PSIT*G/ (461.5*TGS(i)) ))
!       HR(i) = EXP(ARGG)
!    ENDDO
!    RETURN
!  END SUBROUTINE begtem

!  !-----------------------------------------------------------------
!  !
!  !==================SUBROUTINE NETRAD====================================
!  !
!  SUBROUTINE NETRAD(radc3, radt, stefan, fac1,             &
!       vcover, thermk, tc, tg, tf,                        &
!       dtc4, dtg4, dts4, closs, gloss, sloss,             &
!       tgeff, areas, len )
!    IMPLICIT NONE
!    !
!    !=======================================================================
!    !
!    !
!    !        CALCULATE RADT USING RADIATION FROM PHYSICS AND CURRENT
!    !        LOSSES FROM CANOPY AND GROUND
!    !
!    !
!    !=======================================================================
!    !
!    !pl bands in sib: 0.2 to 0.7 micromets are VIS, then
!    !pl               0.7 to 4.0 is NIR, above 4.0 it is thermal
!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       RADt (2)       SUM OF ABSORBED RADIATIVE FLUXES (W M-2)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len
!    REAL  radc3(len,2),radt(len,3),stefan, fac1(len),         &
!         vcover(len),thermk(len),tc(len),tg(len), tf,        &
!         areas(len), tgeff(len),                             &
!         dtc4(len), dtg4(len), dts4(len),                    &
!         closs(len), gloss(len), sloss(len)

!    !     local variables
!    INTEGER i
!    REAL tc4, tg4, ts4, tsnow, radtbar, zlwup              &
!         ,      feedfac  &!  feedback factor
!         ,      stefanx4  !  stefan-boltzmann const X 4.0
!    DATA feedfac/1.0/

!    stefanx4 = 4.0 * stefan

!    DO I=1,len
!       tsnow = MIN(tf, tg(i))
!       !itb...Joerg Kaduk assures me this is quicker than an exponent...
!       TC4 = TC(i) * TC(i) * TC(i) * TC(i)
!       TG4 = TG(i) * TG(i) * TG(i) * TG(i)
!       ts4 = tsnow * tsnow * tsnow * tsnow

!       !pl effective ground cover for thermal radiation
!       FAC1(i) =  VCOVER(i) * ( 1.-THERMK(i) )

!       !pl derivatives
!       DTC4(i) = stefanX4 * TC(i)* TC(i) * TC(i)
!       DTG4(i) = stefanX4 * TG(i)* TG(i) * TG(i)
!       Dts4(i) = stefanX4 * tsnow * tsnow * tsnow

!       CLOSS(i) =  2. * FAC1(i) * STEFAN * TC4
!       !pl canopy leaves thermal radiation loss
!       CLOSS(i) =  CLOSS(i) - FAC1(i) * STEFAN *             &
!            ( (1.-areas(i))*TG4+areas(i)*ts4)
!       !pl ground loss
!       GLOSS(i) =  STEFAN * TG4 - FAC1(i) * STEFAN * TC4
!       !pl snow loss
!       SLOSS(i) =  STEFAN * Ts4 - FAC1(i) * STEFAN * TC4
!       !pl canopy leaves net radiation
!       RADT(i,1) = RADC3(i,1) - closs(i)
!       !pl ground net radiation
!       RADT(i,2) = RADC3(i,2) - gloss(i)
!       !pl snow net radiation
!       RADT(i,3) = RADC3(i,2) - sloss(i)
!       !pl bulk, weighted net radiation from combined ground and snow
!       radtbar = areas(i)*radt(i,3) + (1.-areas(i))*radt(i,2)
!       !pl this is the exchange meant to help out exchanges between
!       !pl ground and snow

!       radt(i,2) = radtbar + (1.+feedfac)*(radt(i,2)-radtbar)
!       radt(i,3) = radtbar + (1.+feedfac)*(radt(i,3)-radtbar)
!       !pl total thermal radiation up from surface
!       zlwup = fac1(i) * tc4 +                                    &
!            (1.-fac1(i)) * (areas(i)*ts4+(1.-areas(i))*tg4)
!       !pl effective (combined) skin temperature from surface thermal radiation
!       TGEFF(i) =  ZLWUP ** 0.25
!    ENDDO

!    RETURN
!  END SUBROUTINE NETRAD

!  !-------------------------------------------------------------------------

!  !--------------------

!  SUBROUTINE VNQSAT (iflag, TQS, PQS, QSS, IM)

!    IMPLICIT NONE

!    !     Computes saturation mixing ratio or saturation vapour pressure
!    !     as a function of temperature and pressure for water vapour.

!    !     INPUT VARIABLES     TQS,PQS,iflag,im
!    !     OUTPUT VARIABLES    QSS
!    !     SUBROUTINES CALLED  (AMAX1,AMIN1)
!    !     iflag = 1 for saturation mixing ratio, otherwise for saturation
!    !             vapor pressure

!    !     Modifications:
!    !             - removed routine VHQSAT and made the call to it in c3vint.F
!    !               a call to VNQSAT adding the iflag=1 argument.  changed
!    !               subroutine name from VQSAT to VNQSAT and all calls to
!    !               VQSAT (c1subs.F, c3subs.F, comp3.F, hstatc.F, ocean.F) to
!    !               VNQSAT with the iflag=1 argument.  kwitt  10/23/91

!    !     converted to fortran 90 syntax - dd 6/17/93

!    !--------------------

!    !     argument declarations

!    INTEGER im,ic(IM), iflag
!    REAL TQS(IM), PQS(IM), QSS(IM)

!    !     local declarations

!    REAL EST(139), tq1(IM), es1(IM), es2(IM), epsinv

!    DATA EST/ 0.0031195, 0.0036135, 0.0041800, 0.0048227, 0.0055571,     &
!         0.0063934, 0.0073433, 0.0084286, 0.0096407, 0.011014,             &
!         0.012582, 0.014353, 0.016341, 0.018574, 0.021095, 0.023926,       &
!         0.027096, 0.030652, 0.034629, 0.039073, 0.044028, 0.049546,       &
!         0.055691, 0.062508, 0.070077, 0.078700, 0.088128, 0.098477,       &
!         0.10983, 0.12233, 0.13608, 0.15121, 0.16784, 0.18615, 0.20627,    &
!         0.22837, 0.25263, 0.27923, 0.30838, 0.34030, 0.37520, 0.41334,    &
!         0.45497, 0.50037, 0.54984, 0.60369, 0.66225, 0.72589, 0.79497,    &
!         0.86991, 0.95113, 1.0391, 1.1343, 1.2372, 1.3484, 1.4684,         &
!         1.5979, 1.7375, 1.8879, 2.0499, 2.2241, 2.4113, 2.6126, 2.8286,   &
!         3.0604, 3.3091, 3.5755, 3.8608, 4.1663, 4.4930, 4.8423,           &
!         5.2155, 5.6140, 6.0394, 6.4930, 6.9767, 7.4919, 8.0406, 8.6246,   &
!         9.2457, 9.9061, 10.608, 11.353, 12.144, 12.983, 13.873, 14.816,   &
!         15.815, 16.872, 17.992, 19.176, 20.428, 21.750, 23.148, 24.623,   &
!         26.180, 27.822, 29.553, 31.378, 33.300, 35.324, 37.454, 39.696,   &
!         42.053, 44.531, 47.134, 49.869, 52.741, 55.754, 58.916, 62.232,   &
!         65.708, 69.351, 73.168, 77.164, 81.348, 85.725, 90.305, 95.094,   &
!         100.10, 105.33, 110.80, 116.50, 122.46, 128.68, 135.17, 141.93,   &
!         148.99, 156.34, 164.00, 171.99, 180.30, 188.95, 197.96, 207.33,   &
!         217.08, 227.22, 237.76, 248.71/

!    EPSINV = 1./1.622

!    tq1(:) = MAX(1.00001E0,(tqs(:)-198.99999))
!    tq1(:) = MIN(138.900001E0,tq1(:))
!    ic(:) = tq1(:)

!    es1(:) = est(ic(:))
!    es2(:) = est(ic(:)+1)

!    qss(:) = ic(:)
!    qss(:) = es1(:) + (es2(:)-es1(:)) * (tq1(:)-qss(:))
!    tq1(:) = pqs(:) * epsinv
!    qss(:) = MIN(tq1(:),qss(:))

!    IF (IFLAG .EQ. 1) qss(:) = 0.622 * qss(:) / (pqs(:)-qss(:))

!    RETURN
!  END SUBROUTINE VNQSAT

!  !----------------------------------------------------------------

!  SUBROUTINE VNTLAT                                                &
!       (grav, tice,                                          &
!       pr0, ribc, vkrmn, delta, dtt, tc, gt, ts, ps, zlt,   &
!       www, tgs, etc, etg,snoww,                            &
!       rstfac, rsoil, hr, wc, wg,snofac,                    &
!       sh, z0, spdm, sha, zb, ros,cas_cap_co2,              &
!       cu, ra, thvgm, rib, ustar, rstar,tstar,              &
!       ventmf, thm, tha, z2, d,                             &
!       fc, fg, rbc, rdc,gect,geci,gegs,gegi,                &
!       respcp, rb, rd, rds, bps, rst, rc, ecmass,           &
!       ea, hrr, assimn, bintc, ta, pco2m, po2m, vmax0,      &
!       green, tran, ref, gmudmu, trop, trda, trdm, slti,    &
!       shti, hltii, hhti, radn, effcon, binter, gradm,      &
!       atheta, btheta, aparkk, wsfws, wsfht, wsflt, wci,    &
!       whs, omepot, assimpot, assimci, antemp, assimnp,     &
!       wags, wegs, aparc, pfd, assim, td, wopt, zm, wsat,   &
!       soilscale, zmstscale, drst,                          &
!       soilq10, ansqr,                                      &
!       nsib, len, nsoil, forcerestore, dotkef,              &
!       thgeff, shgeff, tke, ct, louis, zwind, ztemp,        &
!       respg, respfactor, pco2ap, pco2i, pco2c, pco2s       &
!       ,           co2cap,cflux)

!    !
!    !
!    !     - optimized subroutine phosib
!    !       dd 92.06.10
!    !
!    !
!    IMPLICIT NONE

!    !    argument list declarations
!    INTEGER nsib, len, nsoil
!    LOGICAL forcerestore, dotkef, louis
!    REAL grav, tice, pr0                                             &
!         ,    ribc, vkrmn, delta, dtt, zwind, ztemp                       &
!         ,    snofac
!    REAL tc(len), gt(len), ts(len), ps(len), zlt(len)                &
!         ,    www(len,3)                                                  &
!         ,    tgs(len), etc(len), etg(len)                                &
!         ,    rstfac(len,4), rsoil(len), hr(len)                          &
!         ,    wc(len), wg(len)                                            &
!         ,    sh(len),cflux(len)                                          &
!         ,    z0(len), spdm(len), sha(len), zb(len), ros(len)             &
!         ,    cu(len), ra(len), thvgm(len), rib(len), ustar(len)          &
!         ,    rstar,tstar                                                 &
!         ,    ventmf(len), thm(len), tha(len), z2(len), d(len)            &
!         ,    fc(len), fg(len), rbc(len), rdc(len),rds(len)               &
!         ,    respcp(len), cogr(len), cogs(len), rb(len), rd(len)         &
!         ,    bps(len), rst(len), rc(len), ecmass(len), ea(len)           &
!         ,    hrr(len), assimn(len), bintc(len), ta(len)                  &
!         ,     pco2m(len), po2m, vmax0(len), green(len), tran(len,2,2)    &
!         ,     ref(len,2,2), gmudmu(len), trop(len), trda(len)            &
!         ,     trdm(len), slti(len), shti(len), hltii(len)                &
!         ,     hhti(len), radn(len,2,2), effcon(len), binter(len)         &
!         ,     gradm(len), atheta(len), btheta(len)                       &
!         ,     aparkk(len), wsfws(len), wsfht(len), wsflt(len)            &
!         ,     wci(len), whs(len), omepot(len), assimpot(len)             &
!         ,     assimci(len), antemp(len), assimnp(len), wags(len)         &
!         ,     wegs(len), aparc(len), pfd(len), assim(len), td(nsib,5)    &
!         ,     wopt(len), zm(len), wsat(len), soilscale(len,nsoil+1)      &
!         ,     zmstscale(len,2)                                           &
!         ,     drst(len), soilq10(len,nsoil+1), ansqr(len)                &
!         ,     respg(len)                                                 &
!         ,     thgeff(len), shgeff(len), tke(len), ct(len)                &
!         ,     pco2s(len),pco2i(len)    & !added for neil suits' programs
!         ,     pco2ap(len),pco2c(len)   &! more added nsuits vars
!         ,     co2cap(len)   & ! moles of air in canopy
!         ,     snoww(len,2) & ! snow cover (veg and ground) (m)
!         ,     cas_cap_co2(len)         &
!         ,     gect(len)                &
!         ,     geci(len)                &
!         ,     gegs(len)                &
!         ,     gegi(len)

!    REAL respfactor(nsib, nsoil+1)

!    !    local variables

!    INTEGER i
!    !      LOGICAL  FRSTVM
!    REAL  zln2, ghalf, dmin, pdamp, qdamp, dttin,eps
!    REAL u2(len), cog1(len), cog2(len), vsib(len)                 &
!         ,    thsib(len), coc(len), tprcor(len)                        &
!         ,    cni(len), cuni(len), ctni(len), ctni3(len), cti(len)     &
!         ,    cui(len) , temv(len), cun(len), ctn(len)                 &
!         ,    z1z0Urt(len), z1z0Trt(len), zzwind(len), zztemp(len)     &
!         ,    epsc(len),epsg(len)


!    !
!    EPS    = 1. / SNOFAC
!    !
!    !czzggrst   calculate damping factors
!    zln2 = 6.9314718e-1
!    ghalf = 1.0257068e1
!    dttin = 3.6e3
!    dmin = 6.0e1
!    pdamp = EXP(-1.0 * zln2*(dtt*dmin)/(dttin*ghalf))
!    qdamp = 1.0 - pdamp


!    !czzggrst
!    !     FOR OCEAN POINTS, THIS IS THE ONLY CALL TO VMFCAL.  FOR VEGETATED
!    !     POINTS, ITS JUST THE FIRST CALL BEFORE AN ITERATION.
!    !
!    !      FRSTVM = .TRUE.
!    !xx
!    IF(louis) THEN
!       DO i = 1,len
!          zzwind(i) = z2(i)-d(i)+ zwind
!          zztemp(i) = z2(i)-d(i)+ ztemp
!          !            !print*,'ppppppppp',zzwind(i) ,z2(i),d(i), zwind
!       ENDDO
!       CALL VMFCALZ(VKRMN,DELTA,GRAV                       &
!            ,                SH,z0                                  &
!            ,                SPDM,SHA,ROS,CU,ct,THVGM,RIB           &
!            ,                USTAR,RSTAR,TSTAR,VENTMF,THM           &
!            ,                tha, zzwind, zztemp                    &
!            ,                cuni, cun, ctn, z1z0Urt, z1z0Trt, len )
!    ELSE
!       CALL VMFCAL(PR0,RIBC,VKRMN,DELTA,GRAV                 &
!            ,              SH,Z0,SPDM                                &
!            ,              SHA,ZB,ROS,CU,THVGM,RIB,USTAR             &
!            ,              VENTMF,THM, tha, cni, cuni, ctni          &
!            ,              ctni3, cti, cui, len                      &
!            ,              thgeff, shgeff, tke, ct, dotkef )
!    ENDIF
!    !
!    !     AERODYNAMIC RESISTANCE
!    !

!    DO I=1,len
!       RA(I)    = ROS(i) / VENTMF(i)
!       TEMV(I) = (Z2(i) - D(i)) / z0(i)
!       U2(I)     = SPDM(i) / (CUNI(i) * VKRMN)
!       !print*,'temv',u2(i),temv(I),Z2(i),D(i),z0(i)
!       temv(i) = LOG(temv(i))
!       U2(I) = U2(I) * TEMV(I)

!    ENDDO

!    !      FRSTVM = .FALSE.
!    !
!    !      DO 100 I=1,len
!    !         FC(I)=1.0
!    !         FG(I)=1.0
!    !  100 CONTINUE

!    FC(:) = 1.0
!    FG(:) = 1.0
!    !
!    CALL RBRD(tc, rbc, zlt, z2, u2              &
!         ,         rd, rb, ta,  grav, rdc, tgs       &
!         ,         len  )
!    !


!    DO I=1,len  !DO 50 I=1,len

!       !itb...here is inserted some PL prog CAS stuff...
!       epsc(i) = 1.
!       epsg(i) = 1.
!       !itb...pl says " this only makes sense for canopy leaves, since
!       !itb...there can only be water OR snow, not both. switching epsc
!       !itb...epsc to eps makes the hltm adapt to freezing/fusion.

!       IF(snoww(i,1) .GT. 0.0) epsc(i) = eps
!       IF(snoww(i,2) .GT. 0.0) epsg(i) = eps

!       RC(i) = RST(i) + RB(i) + RB(i)

!       RDS(i) = RSOIL(i) * FG(i) + RD(i)

!       GECT(i) =  (1. - WC(i)) / RC(i)
!       GECI(i) = epsc(i) * WC(i) / (RB(i) + RB(i))

!       GEGS(i) =  (1. - WG(i)) / RDS(i)
!       GEGI(i) = epsg(i) * WG(i) / RD(i)

!       COC(i) = GECT(i) + GECI(i)


!       VSIB(I)  = 1.0/RB(I) + 1.0/RD(I)
!       THSIB(I) = (GT(i)/RD(I) + TC(i)/RB(I)) /    &
!            (BPS(i)*VSIB(I))
!       ENDDO  !50 CONTINUE
!       !

!       !czzggrst calculate ecmass -- canopy evapotranspiration
!    !
!    DO i=1,len
!       ecmass(i) = (etc(i) - ea(i)) * coc(i) *              &
!            ros (i) * 0.622e0 /ps(i) * dtt
!    ENDDO


!    !pl include here a call to respsib
!    !pl pass it

!    CALL respsib(len, nsib, nsoil, wopt, zm, www, wsat,                &
!         tgs, td, forcerestore, respfactor, respg, soilscale,  &
!         zmstscale, soilq10)


!    !czzggrst
!    !
!    !     calculation of canopy conductance and photosynthesis
!    !
!    CALL phosib(pco2m,pco2ap,po2m,vmax0,tice,ps,green                 &
!         ,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda        &
!         ,          trdm,slti,shti,hltii,hhti                              &
!         ,          radn,etc                                               &
!         ,          ea,rb,ra,ts                                            &
!         ,          effcon,rstfac,binter,gradm,assimn                      &
!         ,          rst,atheta                                             &
!         ,          btheta,tgs,respcp,aparkk, len, nsib,                   &
!         omepot,assimpot,assimci,antemp,assimnp,                      &
!         wsfws,wsfht,wsflt,wci,whs,                                   &
!         wags,wegs,aparc,pfd,assim, td,                               &
!         soilscale,                                                   &
!         drst, pdamp,qdamp,ecmass,dtt,bintc,tprcor,ansqr,             &
!         nsoil, respg, pco2c, pco2i,                                  &
!         pco2s,co2cap,cflux)


!    !czzggrst block moved up (cxx-140-2000 loop)
!    !

!    DO i = 1,len
!       bintc(i) = bintc(i) * tc(i) / ( 44.6 * tprcor(i))
!       IF(ea(i).GT.etc(i)) fc(i) = 0.0
!       IF(ea(i).GT.etg(i)) fg(i) = 0.0
!       HRR(I) = HR(I)
!       IF (FG(I) .LT. 0.5) HRR(I) = 1.0
!    ENDDO
!    !


!    RETURN
!  END SUBROUTINE VNTLAT

!  !--------------------------------------------------------------

!  SUBROUTINE VMFCALZ(VKRMN,DELTA,GRAV                         &
!       ,                SH,z0                                      &
!       ,                SPDM,SHA,ROS,CU,ct,THVGM,RIB               &
!       ,                USTAR,RSTAR,TSTAR,VENTMF,THM, tha          &
!       ,                zzwind, zztemp                             &
!       ,                cuni, cun, ctn, z1z0Urt, z1z0Trt, len )

!    IMPLICIT NONE

!    !****************************************************************************
!    !    VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum
!    !    and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
!    !    and Boville(1993), and by Beljaars and Holtslag (1991).
!    !
!    !     Rerences:
!    !       Beljars and Holtslag (1991): Flux parameterization over land surfaces
!    !              for atmospheric models. J. Appl. Meteo., 30, 327-341.
!    !       Holtslag and Boville (1993): Local versus nonlocal boundary-layer
!    !              diffusion in a global climate model. J. of Climate, 6, 1825-
!    !              1842.
!    !       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
!    !              atmosphere. Boundary-Layer Meteo., 17, 187-202.
!    !       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
!    !              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
!    !              Layer parameterization, ECMWF, 59-79.
!    !
!    !     General formulation:
!    !        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.)
!    !     Transfer coefficients for mommentum and heat fluxes are:
!    !        CU = CUN*Fm, and
!    !        CT = CTN*Fh
!    !        where  CUN and CTN are nutral values of momentum and heat transfers,
!    !           and Fm and Fh are stability functions derived from surface
!    !           similarity relationships.
!    !****************************************************************************

!    INTEGER len
!    REAL  vkrmn, delta, grav                             &
!         ,    z0(len),CU(len)                                 &
!         ,    ROS(len),SH(len)                                &
!         ,    THM(len),THVGM(len),RIB(len),THGM(len)          &
!         ,    SHA(len),SPDM(len),USTAR(len)                   &
!         ,    RSTAR(len),TSTAR(len)                           &
!         ,    CUI(len),CTI(len),CT(len), tha(len)             &
!         ,    VENTMF(len), cuni(len)

!    !     local variables
!    REAL TEMV(len), wgm(len)                                           &
!         ,     bunstablM, bunstablT, cunstablM, cunstablT, bstabl, cstabl   &
!         ,    zzwind(len), zztemp(len), zrib(len), cun(len), ctn(len)       &
!         ,      z1z0U(len), z1z0Urt(len), z1z0T(len), z1z0Trt(len)          &
!         ,      fmomn(len),fheat(len),ribtemp,dm,dh

!    INTEGER i
!    !
!    !  condtants for surface flux functions, according to Holtslag and
!    !      Boville (1993, J. Climate)

!    bunstablM = 10.       ! constants for unstable function
!    bunstablT = 15
!    cunstablM = 75.
!    cunstablT = 75.
!    bstabl = 8.           ! constants for stable function
!    cstabl = 10.
!    !
!    DO I=1,len  !DO 10 I=1,len
!       zrib(i) = zzwind(i) **2 / zztemp(i)
!       !print*,'WGM:',sha(i),sh(i)
!       WGM(i)  = SHA(i) - SH(i)
!       !
!       !        SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE
!       !        WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.
!       !
!       THGM(i)  = THA(i)  - THM(i)
!       THVGM(i) = THGM(i) + THA(i) * DELTA * WGM(i)
!    ENDDO !10   CONTINUE

!    !   Ratio of reference height (zwind/ztemp) and roughness length:
!    DO i = 1, len                !for all grid points
!       z1z0U(i) = zzwind(i)/ z0(i)
!       z1z0Urt(i) = SQRT( z1z0U(i) )
!       z1z0U(i) = LOG( z1z0U(i) )
!       z1z0T(i) = zzwind(i)/ z0(i)
!       z1z0Trt(i) = SQRT( z1z0T(i) )
!       z1z0T(i) = LOG( z1z0T(i) )
!    ENDDO

!    !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

!    DO i = 1, len
!       cun(i) = VKRMN*VKRMN / (z1z0U(i)*z1z0U(i) )   !neutral Cm & Ct
!       ctn(i) = VKRMN*VKRMN / (z1z0T(i)*z1z0T(i) )
!       cuni(i) = z1z0u(i) / vkrmn
!    ENDDO
!    !
!    !   SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
!    !   RIB IS THE BULK RICHARDSON NUMBER, between reference height and surface.
!    !
!    DO I=1,len
!       TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
!       temv(i) = MAX(0.000001E0,temv(i))
!       RIB(I) = -THVGM(I) * GRAV * zrib(i) / TEMV(i)
!    ENDDO

!    !   The stability functions for momentum and heat/moisture fluxes as
!    !   derived from the surface-similarity theory by Luis (1079, 1982), and
!    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag
!    !   (1991).

!    DO I=1,len  !DO 60 I=1,len
!       IF(rib(i).GE.0.0) THEN
!          !
!          !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
!          !
!          rib(i) = MIN( rib(i), 0.5E0)
!          fmomn(i) = (1. + cstabl * rib(i) * (1.+ bstabl * rib(i)))
!          fmomn(i) = 1. / fmomn(i)
!          fmomn(i) = MAX(0.0001E0,fmomn(i))
!          fheat(i) = fmomn(i)

!       ELSE
!          !
!          !        THE UNSTABLE CASE.
!          !
!          ribtemp = ABS(rib(i))
!          ribtemp = SQRT( ribtemp )
!          dm = 1. + cunstablM * cun(i) * z1z0Urt(i) * ribtemp
!          dh = 1. + cunstablT * ctn(i) * z1z0Trt(i) * ribtemp
!          fmomn(i) = 1. - (bunstablM * rib(i) ) / dm
!          fheat(i) = 1. - (bunstablT * rib(i) ) / dh

!       END IF
!    ENDDO !60   CONTINUE

!    !   surface-air transfer coefficients for momentum CU, for heat and
!    !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

!    DO i = 1, len
!       CU(i) = CUN(i) * fmomn(i)
!       CT(i) = CTN(i) * fheat(i)
!       CUI(i) = 1. / CU(i)
!       CTI(i) = 1. / CT(i)
!    ENDDO

!    !   Ustar and ventlation mass flux: note that the ustar and ventlation
!    !   are calculated differently from the Deardoff's methods due to their
!    !   differences in define the CU and CT.

!    !print*,'VNTLAT: STAR'
!    DO I=1,len  !DO 80 I=1,len
!       USTAR(i) = SPDM(i)*SPDM(i)*CU(i)
!       !print*,'ustar:',SQRT(ustar(i)),spdm(i),cu(i)
!       USTAR(i) = SQRT( USTAR(i) )
!       RSTAR(i) = -WGM(i)*CU(i)*fmomn(i)/USTAR(i)
!       !print*,'rstar:',rstar(i),-wgm(i),cu(i),fmomn(i),ustar(i)
!       TSTAR(i) = -THGM(i)*CT(i)*fheat(i)/USTAR(i)
!       !print*,'tstar:',tstar(i),-thgm(i),ct(i),fheat(i),ustar(i)
!       VENTMF(i)= ROS(i)*CT(i)* SPDM(i)
!    ENDDO !80   CONTINUE
!    !
!    !   Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION
!    !   VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
!    !   free convection conditions.
!    !

!    RETURN
!  END SUBROUTINE VMFCALZ

!  !-----------------------------------------------------------------

!  SUBROUTINE VMFCAL(PR0,RIBC,VKRMN,DELTA,GRAV             &
!       ,                SH,z0                                  &
!       ,                SPDM,SHA,ZB,ROS,CU,THVGM,RIB           &
!       ,                USTAR,VENTMF,THM, tha                  &
!       ,                cni, cuni, ctni, ctni3, cti, cui       &
!       ,                len, thgeff, shgeff, tke, ct, dotkef )
!    IMPLICIT NONE
!    !
!    !***  VENTILATION MASS FLUX, BASED ON DEARDORFF, MWR, 1972
!    !
!    INTEGER len
!    LOGICAL dotkef
!    REAL pr0, ribc, vkrmn, delta, grav                         &
!         ,    z0(len),CU(len)                                       &
!         ,    ZB(len),ROS(len),SH(len)                              &
!         ,    THM(len),THVGM(len),RIB(len),THGM(len),CNI(len)       &
!         ,    CUNI(len),SHA(len),SPDM(len),USTAR(len),CTNI(len)     &
!         ,    CUI(len),CTI(len),CT(len), tha(len)                   &
!         ,    VENTMF(len),CTNI3(len)                                &
!         ,    thgeff(len),shgeff(len),tke(len)

!    !     local variables
!    REAL ribmax, vkrinv, athird, TEMV(len), ZDRDRF(len), CHOKE(len)   &
!         ,    wgm(len), sqrtke
!    INTEGER i
!    !
!    RIBMAX = 0.9*RIBC
!    VKRINV = 1.0/VKRMN
!    ATHIRD = 1.0/3.0
!    !
!    DO I=1,len !DO 10 I=1,len
!       WGM(i)  = SHA(i) - SH(i)
!       !
!       !        SFC-AIR DEFICITS OF MOISTURE AND POTENTIAL TEMPERATURE
!       !        WGM IS THE EFFECTIVE SFC-AIR TOTAL MIXING RATIO DIFFERENCE.
!       !
!       THGM(i)  = THA(i)  - THM(i)
!    ENDDO !10   CONTINUE
!    IF (dotkef) THEN
!       DO i = 1,len
!          thvgm(i) = (thgeff(i)-thm(i)) + thgeff(i) * delta *      &
!               ((shgeff(i)-sh(i)))
!       ENDDO
!    ELSE
!       DO i = 1,len
!          THVGM(i) = THGM(i) + THA(i) * DELTA *             &
!               WGM(i)
!       ENDDO
!    ENDIF
!    !
!    !        CUNI AND CTN1 ARE INVERSES OF THE NEUTRAL TRANSFER COEFFICIENTS
!    !        DEARDORFF EQS(33) AND (34).
!    !        PR0 IS THE TURBULENT PRANDTL NUMBER AT NEUTRAL STABILITY.
!    !
!    !
!    DO I=1,len !DO 20 I=1,len
!       CNI(i) = 0.025*ZB(i)/z0(i)
!       cni(i) = LOG(cni(i))
!       CNI(i)  = CNI(i)  * VKRINV
!       CUNI(i) = CNI(i)  + 8.4
!       CTNI(i) = CNI(i)  * PR0 + 7.3
!       CTNI3(i)= 0.3     * CTNI(i)
!       CNI(i)  = GRAV    * ZB(i)
!       !
!       !        SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
!       !        RIB IS THE BULK RICHARDSON NUMBER, DEARDORFF EQ (25)
!       !
!       TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
!       temv(i) = MAX(0.000001E0,temv(i))
!    ENDDO !20   CONTINUE
!    IF (dotkef) THEN
!       DO I=1,len
!          rib(i) = -thvgm(i)*grav*zb(i)/(thgeff(i)*tke(i))
!       ENDDO
!    ELSE
!       DO I=1,len
!          RIB(I) = -THVGM(I) * CNI(I) / TEMV(i)
!       ENDDO
!    ENDIF
!    !
!    !
!    IF (dotkef) THEN
!       DO I=1,len
!          IF(rib(i).GE.0.0) THEN
!             cu(i) = ((-0.0307*rib(i)**2)/(61.5+rib(i)**2))+0.044
!             ct(i) = ((-0.0152*rib(i)**2)/(190.4+rib(i)**2))+0.025
!          ELSE
!             cu(i) = ((-0.016*rib(i)**2)/(4.2e4+rib(i)**2))+0.044
!             ct(i) = ((-0.0195*rib(i)**2)/(2.1e4+rib(i)**2))+0.025
!          ENDIF
!       ENDDO
!       DO I=1,len
!          sqrtke = SQRT(tke(i))
!          USTAR(I) = SQRT(sqrtke*spdm(I)*CU(I))
!          VENTMF(I)=ROS(I)*CT(I)*sqrtke
!       ENDDO
!    ELSE
!       DO I=1,len  !DO 60 I=1,len
!          IF(rib(i).GE.0.0) THEN
!             !
!             !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
!             !
!             CHOKE(i)=1.0-MIN(RIB(I),RIBMAX)/RIBC
!             CU(I)=CHOKE(i)/CUNI(I)
!             CT(I)=CHOKE(i)/CTNI(I)
!          ELSE
!             !
!             !        FIRST, THE UNSTABLE CASE. DEARDORFF EQS(28), (29), AND (30)
!             !
!             ZDRDRF(i) = LOG10(-RIB(I)) - 3.5
!             CUI(I)    = CUNI(I) - 25.0 * EXP(0.26 * ZDRDRF(i)      &
!                  - 0.03 * ZDRDRF(i) * ZDRDRF(i))
!             CTI(I)    = CUI(I) + CTNI(I) - CUNI(I)
!             CU(I)     = 1.0 / MAX(CUI(I),0.5 * CUNI(I))
!             CT(I)     = 1.0 / MAX(CTI(I),CTNI3(I))
!          END IF
!       ENDDO  !60      CONTINUE
!       DO I=1,len  !DO 80 I=1,len
!          USTAR(i) =SPDM(i)*CU(i)
!          VENTMF(i)=ROS(i)*CT(i)*USTAR(i)
!       ENDDO  !80      CONTINUE
!    ENDIF

!    !
!    !     CHECK THAT VENTFC EXCEEDS TOWNSENDS(1962) FREE CONVECTION VALUE,
!    !     DEARDORFF EQ(40B)
!    !
!    IF (.NOT.dotkef) THEN
!       DO I=1,len  !DO 90 I=1,len
!          IF( rib(i).LT.0.0) THEN
!             IF( cti(i).LT.ctni3(i) )          &
!                  VENTMF(I)=MAX(VENTMF(I),       &
!                  ROS(I)*0.00186*(THVGM(I))**ATHIRD)
!          ENDIF
!       ENDDO  !90 CONTINUE
!    ENDIF
!    RETURN
!    !
!  END SUBROUTINE VMFCAL

!  !------------------------------------------------------------------

!  !cc
!  SUBROUTINE respsib(len, nsib, nsoil, wopt, zm, www, wsat,         &
!       tg, td, forcerestore, respfactor, respg, soilscale,   &
!       zmstscale, soilq10)

!    IMPLICIT NONE

!    INTEGER len, nsib, nsoil,i,l
!    REAL wopt(len), zm(len), www(len,3), wsat(len),     &
!         tg(len), td(nsib,nsoil),woptzm
!    REAL respfactor(nsib,nsoil+1)
!    ! output
!    REAL respg(len),soilscale(len,nsoil+1), zmstscale(len,2),    &
!         soilq10(len,nsoil+1),                                    &
!         ! local arrays
!         b(len,2)

!    LOGICAL forcerestore
!    !    Calculates the rate of CO2 efflux from soils, according to the "R-star"
!    !    approach of Denning et al (1996), adapted for use with the Bonan
!    !    6-layer soil thermodynamics module

!    !    Changed soil Q10 value for soil respiration from 2.0 to 2.4
!    !    following Raich and Schelsinger (1992, Tellus 44B, 81-89),
!    !    Scott Denning, 9/14/95

!    IF(.NOT.forcerestore) THEN

!       DO i = 1,len

!          ! Moisture effect on soil respiration, shallow and root zone
!          woptzm = wopt(i)**zm(i)
!          b(i,1) = (((100.*www(i,1))**zm(i)-woptzm)/                &
!               (woptzm - 100.**zm(i)))**2
!          b(i,2) = (((100.*www(i,2))**zm(i)-woptzm)/                &
!               (woptzm - 100.**zm(i)))**2
!          b(i,1) = MIN(b(i,1),10.E0)
!          b(i,2) = MIN(b(i,2),10.E0)
!          zmstscale(i,1) = 0.8*wsat(i)**b(i,1) + 0.2
!          zmstscale(i,2) = 0.8*wsat(i)**b(i,2) + 0.2

!          ! Temperature effect is a simple Q10, with a reference T of 25 C

!          ! Deepest soil layers do not respire (no carbon below root zone)
!          DO L = 1, 2
!             soilscale(i,L) = 0.0
!             soilq10(i,L) = 0.0
!          ENDDO

!          ! Layers 3 through nsoil (root zone) use WWW(2) and TD(3:nSoil)
!          DO l = 3, nsoil
!             soilQ10(i,L) = EXP(0.087547 * (td(i,L) - 298.15))
!             soilscale(i,L) = soilQ10(i,L) * zmstscale(i,2)
!          ENDDO

!          ! Surface soil uses TG and water layer 1
!          soilQ10(i,nsoil+1) = EXP(0.087547 * (tg(i) - 298.15))
!          soilscale(i,nsoil+1) = soilQ10(i,nsoil+1) * zmstscale(i,1)

!          ! Dimensionalize soil resp flux to balance annual budget
!          respg(i) = respfactor(i,1) * soilscale(i,1)

!          DO L = 2, nSoil+1
!             respg(i) = respg(i) + respfactor(i,L) * soilscale(i,L)
!          ENDDO

!       ENDDO



!    ELSE  ! (FORCERESTORE CASE ... only two soil T levels available)

!       DO i = 1,len

!          ! Moisture effect from TEM
!          woptzm = wopt(i)**zm(i)
!          b(i,2) = (((100.*www(i,2))**zm(i)-woptzm)/              &
!               (woptzm - 100.**zm(i)))**2
!          b(i,2) = MIN(b(i,2),10.E0)
!          zmstscale(i,1) = 0.8*wsat(i)**b(i,2) +0.2

!          ! Temperature effect is Q10 =2.4 from ref T of 25 C
!          soilQ10(i,1) = EXP(0.087547 * (td(i,nsoil) - 298.15))
!          soilscale(i,1) = soilQ10(i,1) * zmstscale(i,1)

!          ! Dimensionalize soil resp flux to balance annual budget
!          respg(i) = respfactor(i,1) * soilscale(i,1)

!       ENDDO


!    ENDIF



!    RETURN
!  END SUBROUTINE respsib

!  !----------------------------------------------------------------------

!  !==================SUBROUTINE PHOSIB===================================
!  !
!  SUBROUTINE PHOSIB(pco2m,pco2ap,po2m,vmax0,tf,psur,green        &
!       ,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda     &
!       ,          trdm,slti,shti,hltii,hhti,radn,etc                  &
!       ,          ea,rb,ra,tm                                         &
!       ,          effcon,rstfac,binter,gradm,assimn                   &
!       ,          rst,atheta,btheta,tgs,respcp                        &
!       ,          aparkk,len,nsib,                                    &
!       omepot,assimpot,assimci,antemp,assimnp,                   &
!       wsfws,wsfht,wsflt,wci,whs,                                &
!       wags,wegs,aparc,pfd,assim,td,                             &
!       soilscale,                                                &
!       drst,pdamp,qdamp,ecmass,dtt,bintc,tprcor,ansqr,           &
!       nsoil, respg, pco2c, pco2i,                               &
!       pco2s,co2cap,cflux)

!    IMPLICIT NONE

!    !
!    !
!    !=======================================================================
!    !
!    !     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED
!    !     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
!    !     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
!    !     BASE REFERENCE IS SE-92A
!    !
!    !                          UNITS
!    !                         -------
!    !
!    !      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
!    !      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
!    !      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
!    !      EFFCON                                   : MOL CO2 MOL QUANTA-1
!    !      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
!    !      EVAPKG                                   : KG M-2 S-1
!    !      Q                                        : KG KG-1
!    !
!    !                       CONVERSIONS
!    !                      -------------
!    !
!    !      1 MOL H2O           = 0.018 KG
!    !      1 MOL CO2           = 0.044 KG
!    !      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
!    !      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
!    !pl the next line applies to the Ci to Cs pathway
!    !      GS  (CO2)           = GS (H2O) * 1./1.6
!    !pl 44.6 is the number of moles of air per cubic meter
!    !      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
!    !      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!    !      MM  (MOLAIR/MOLH2O) = 1.611
!    !
!    !
!    !                         OUTPUT
!    !                      -------------
!    !
!    !      ASSIMN              = CANOPY NET ASSIMILATION RATE
!    !      EA                  = CANOPY AIR SPACE VAPOR PRESSURE
!    !      1/RST               = CANOPY CONDUCTANCE
!    !      PCO2I               = INTERNAL CO2 CONCENTRATION
!    !      RESPC               = CANOPY RESPIRATION
!    !      RESPG               = GROUND RESPIRATION
!    !
!    !----------------------------------------------------------------------
!    !
!    !         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
!    !         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
!    !         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
!    !         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
!    !
!    !-----------------------------------------------------------------------
!    !

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1)
!    !       RST            CANOPY RESISTANCE (S M-1)
!    !       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS
!    !
!    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!    !
!    !       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
!    !       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
!    !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
!    !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
!    !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    !     Modifications:
!    !       - gs (stomatal conductance reduced for freezing soils per Jim Collatz
!    !         dd 950221
!    !
!    !      Modified for multitasking - introduced gather/scatter indices
!    !          - DD 951206
!    !
!    !itb   Added in pco2c (chloroplast partial co2) for neil's fractionation
!    !itb   calculations
!    !itb       - IB Sep99

!    !
!    !     input arrays:
!    INTEGER len, nsib, nsoil

!    REAL vmax0(len),psur(len),green(len),gmudmu(len),        &
!         zlt(len),cas_cap_co2(len),tc(len),ta(len),           &
!         trop(len),trda(len),                                 &
!         slti(len),shti(len),hltii(len),hhti(len),            &
!         ra(len),rb(len),                                     &
!         cog1(len),cog2(len),tm(len),effcon(len),             &
!         binter(len),gradm(len),atheta(len),btheta(len),      &
!         tgs(len),respcp(len),tran(len,2,2),ref(len,2,2),     &
!         radn(len,2,2),ecmass(len),trdm(len),etc(len),        &
!         aparc(len),rst(len), cflux(len)
!    REAL pdamp, qdamp, dtt, pco2m(len), pco2ap(len), tf, po2m

!    !     output arrays:


!    REAL assimn(len),ea(len),rstfac(len,4),            &
!         pco2i(len),respc(len),respg(len),drst(len)
!    !zz new diagostics 10/14/92
!    !
!    ! output arrays
!    REAL omepot(len),assimpot(len),assimci(len),        &
!         assimnp(len),whs(len),antemp(len),             &
!         wsfws(len),wsfht(len),wsflt(len),wci(len),     &
!         wags(len),wegs(len),pfd(len),                  &
!         td(nsib,nsoil),assim(len),                     &
!         soilscale(len,nsoil+1),                        &
!         tprcor(len),bintc(len),                        &
!         ansqr(len)                                     &
!         ,pco2c(len) &!chloroplast pco2
!         ,xgah2o(len)        &
!         ,xgco2m(len)

!    !     work arrays:

!    REAL PCO2Y(len,6), EYY(len,6),assimny(len,6),    &
!         assimy(len,6)
!    REAL c3(len),c4(len),RANGE(len),gammas(len),       &
!         aparkk(len),gah2o(len),                       &
!         gbh2o(len),                                   &
!         par(len),rrkk(len),                           &
!         omss(len),vm(len),gsh2o(len),pco2s(len),      &
!         templ(len),temph(len),                        &
!         qt(len),co2s(len),scatp(len),scatg(len),      &
!         park(len),respn(len),zkc(len),                &
!         zko(len),spfy(len), co2a(len), co2m(len),co2cap(len)

!    INTEGER icconv(len),igath(len)

!    REAL soilfrz(len)
!    REAL  cwsflt, cwsfht, cwsfws,                                  &
!         ccoms, ccomc, ascitemp, dompdomc, omsci, ompci, omcci,    &
!         omcpot, omppot, sqrtin, omspot, pco2ipot, ohtp2, sttp2,   &
!         gsh2oinf, h2osrh, h2os, ecmole, h2oa, h2oi, dtti,         &
!         pco2in, pco2a, soilfrztd, soilfrztg

!    INTEGER i, ic1, ic, l



!    !pl introduce a co2 capacity
!    !pl this will basically be the mass of air under the top of the canopy (in
!    !pl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.

!    !pl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
!    !pl carbon are in m2 s / mol air

!    !pl one mole of gas occupies 22.4 cubic dm
!    !pl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
!    !pl the units for the carbon capacity are mol air /m2.
!    !pl (e.g. here 893 moles if thickness of the layer is 20m)
!    !pl this means that the units for pc02ap should be mol co2 / mol air, but
!    !pl it is also possible to keep just co2 pressure and convert


!    DO i = 1,len                !  LOOP OVER GRIDPOINT

!       TPRCOR(i) = TF*PSUR(i)*100./1.013E5
!       co2cap(i) = cas_cap_co2(i) * 44.6 * tprcor(i)/ta(i)     ! moles air / m2

!       !pl this needs to be modified as in sibslv3 to automatically use the
!       !pl thickness of the canopy air space.

!       !
!       !----------------------------------------------------------------------
!       !
!       !pl        RESPG(i) = 0. E -6 ! fixed respiration at 5 micromoles
!       !pl   no longe rused since we now have respsib
!       !
!       !----------------------------------------------------------------------
!       !
!       IF( EFFCON(i) .GT. 0.07 ) THEN
!          C3(i) = 1.
!       ELSE
!          C3(i) = 0.
!       ENDIF
!       C4(i)     = 1. - C3(i)

!       !
!       !-----------------------------------------------------------------------
!       !
!       !     CALCULATION OF CANOPY PAR USE PARAMETER.
!       !
!       !      APARKK      (PI)     : EQUATION (31) , SE-92A
!       !-----------------------------------------------------------------------
!       !
!       SCATP(I) =     GREEN(i)   *              &
!            ( TRAN(i,1,1) + REF(i,1,1) )    &
!            +( 1.-GREEN(i) ) *              &
!            ( TRAN(i,1,2) + REF(i,1,2) )
!       SCATG(i) = TRAN(i,1,1) + REF(i,1,1)
!       PARK(i) = SQRT(1.-SCATP(i)) * GMUDMU(i)
!       !
!       ! Collatz-Bounoua commented the calculation of  aparc
!       ! replaced it with theone calculated in new_mapper.
!       !
!       !b        APARC(i) = 1. - EXP ( -PARK(i)*ZLT(i) )   ! lahouari
!       !
!       APARKK(i)   = APARC(i) / PARK(i) * GREEN(i)
!       !-----------------------------------------------------------------------
!       !
!       !     Q-10 AND STRESS TEMPERATURE EFFECTS
!       !
!       !      QT          (QT)    : TABLE (2)     , SE-92A
!       !-----------------------------------------------------------------------
!       !
!       qt(i) = 0.1*( TC(i) - TROP(i) )
!       RESPN(i) = RESPCP(i) * VMAX0(i) * RSTFAC(i,2)

!       !itb...patch to prevent underflow if temp is too cool...
!       IF(TC(i) >= TRDM(i))THEN
!          RESPC(i) = RESPN(i) * 2.0**qt(i)            &
!               /( 1. + EXP( TRDA(i)*(TC(i)-TRDM(i))))
!       ELSE
!          RESPC(i) = RESPN(i) * 2.0**qt(i)
!       ENDIF

!       VM(i) = VMAX0(i) * 2.1**qt(i)
!       TEMPL(i) = 1. + EXP(SLTI(i)*(HLTIi(i)-TC(i)))
!       TEMPH(i) = 1. + EXP(SHTI(i)*(TC(i)-HHTI(i)))
!       RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))
!       VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i)    &
!            + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)
!       !
!       !-----------------------------------------------------------------------
!       !
!       !     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
!       !     COMPENSATION POINT
!       !
!       !      ZKC          (KC)     : TABLE (2)     , SE-92A
!       !      ZKO          (KO)     : TABLE (2)     , SE-92A
!       !      SPFY         (S)      : TABLE (2)     , SE-92A
!       !      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
!       !      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
!       !      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
!       !-----------------------------------------------------------------------
!       !
!       ZKC(i) = 30. * 2.1**qt(i)
!       ZKO(i) = 30000. * 1.2**qt(i)
!       SPFY(i) = 2600. * 0.57**qt(i)
!       GAMMAS(i) = 0.5 * PO2M/SPFY(i) * C3(i)
!       PFD(i)    = 4.6E-6 * GMUDMU(i)*               &
!            ( RADN(i,1,1)+RADN(i,1,2) )
!       !
!       !pl these here all go from being m/s to being mol/ (m2 sec)
!       GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
!       GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
!       GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)

!       xgah2o(i) = MAX(0.466E0, gah2o(i) )
!       xgco2m(i) = 4000.0 * vmax0(i)
!       !
!       RRKK(i)   = ZKC(i)*( 1. + PO2M/ZKO(i) ) * C3(i)        &
!            + VMAX0(i)/5.* ( 1.8**qt(i)) * C4(i)
!       PAR(i)    = pfd(i)*EFFCON(i)*( 1.-SCATG(i) )
!       soilfrztg = 1.+EXP(-1.5 * (MAX(270.0E0,tgs(i))-273.16) )
!       soilfrztd = 1.+EXP(-1.5 * (MAX(270.0E0,td (i,nsoil))-273.16) )
!       soilfrz(i) = MAX(1./soilfrztg, 1./soilfrztd)
!       soilfrz(i) = MAX( soilfrz(i), 0.05E0)
!       BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)*                &
!            RSTFAC(i,2) * soilfrz(i)
!  !     print'(a,5g16.6)','bintc:',binter(i),zlt(i),green(i),   &
!  !                rstfac(i,2),soilfrz(i)
!       OMSS(i)   = ( VMAX0(i)/2.0 ) * ( 1.8**qt(i) )         &
!            /TEMPL(i) * RSTFAC(i,2) * C3(i)       &
!            + RRKK(i) * RSTFAC(i,2) * C4(i)
!       !
!       !-----------------------------------------------------------------------
!       !
!       !     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
!       !     ASSIMILATION RATE.
!       !
!       !-----------------------------------------------------------------------


!       RANGE(i)    = PCO2M(i) * ( 1. - 1.6/GRADM(i) ) - GAMMAS(i)
!       icconv(i) = 1

!    ENDDO

!    !
!    DO IC = 1, 6 !DO 1000 IC = 1, 6
!       DO i = 1,len        ! LOOP OVER GRIDPOINT
!          PCO2Y(i,IC) = 0.
!          EYY(i,IC) = 0.
!       ENDDO
!    ENDDO        !1000 CONTINUE
!    !

!    !pl beginning of PL's setup

!    DO i=1,len
!       gah2o(i) =  1. / MAX(0.446E0,GAH2O(i))
!    ENDDO

!    DO IC = 1, 6  !DO 2000 IC = 1, 6

!       !
!       CALL       SORTIN( EYY, PCO2Y, RANGE, GAMMAS, ic,len )

!       CALL       CYCALC( APARKK, VM, ATHETA, BTHETA,par,           &
!            GAMMAS, RESPC, RRKK, OMSS, C3, C4,        &
!            PCO2Y(1,ic), assimny(1,ic), assimy(1,ic), &
!            len  )
!       !

!       DO i = 1,len

!          !pl now prognose the new CAS CO2 according to flux divergence
!          !pl we are going to do this in mol C / mol air (same as PaC/PaAir)

!          CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
!          co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)

!          CO2A(i)   = (  CO2A(i) + (dtt/co2cap(i)) *         &
!               (respg(i) - assimny(i,ic)              &
!               +co2m(i)*gah2o(i)        ) )           &
!               / (1+dtt*gah2o(i)/ co2cap(i) )

!          pco2a = co2a(i) * psur(i) * 100.

!          PCO2S(i) = PCO2A - (1.4/GBH2O(i) * ASSIMNy(i,ic)      &
!               * PSUR(i)*100.)
!          PCO2IN   = PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)   &
!               * PSUR(i)*100.)
!          EYY(i,IC) = PCO2Y(i,IC) - PCO2IN
!       ENDDO
!       !
!       IF(ic.GE.2) THEN
!          ic1 = ic-1
!          DO i = 1,len        ! LOOP OVER GRIDPOINT
!             IF(ABS(eyy(i,ic1)).GE.0.1)THEN
!                icconv(i) = ic
!             ELSE
!                eyy(i,ic) = eyy(i,ic1)
!                pco2y(i,ic) = pco2y(i,ic1)
!             ENDIF
!          ENDDO
!       ENDIF
!       !
!    ENDDO !2000 CONTINUE
!    !
!    !

!    !pl end of PL's setup

!    DO i = 1,len        ! LOOP OVER GRIDPOINT
!       icconv(i) = MIN(icconv(i),6)
!       igath(i) = i+(icconv(i)-1)*len
!    ENDDO


!    DO i = 1,len         ! LOOP OVER GRIDPOINT

!       pco2i(i) = pco2y(igath(i),1)
!       assimn(i) = assimny(igath(i),1)
!       assim(i)  = assimy(igath(i),1)



!       !        pco2i(i)  = pco2y  (i,icconv(i))
!       !        assimn(i) = assimny(i,icconv(i))
!       !  assim(i)  = assimy (i,icconv(i))

!       pco2c(i) = pco2i(i) - assimn(i)/xgco2m(i)*psur(i)*100.0
!       !print*,'pco2c:',pco2c(i),pco2i(i),assimn(i),xgco2m(i),psur(i)

!       !pl now do the real C_A forecast with the iterated fluxes.

!       CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
!       co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)

!       CO2A(i) = (CO2A(i) + (dtt/co2cap(i)) *          &
!            (respg(i) - assimn(i)               &
!            +co2m(i)*gah2o(i) ) )               &
!            / (1+dtt*gah2o(i)/co2cap(i))
!       !pl go back from molC / mol air to Pascals

!       pco2ap(i) = co2a(i) * psur(i) * 100.
!       !print*,'pco2ap',pco2ap(i),co2a(i),psur(i)

!       !itb...carbon flux between CAS and reference level
!       cflux(i) = gah2o(i)*(co2a(i)-co2m(i))*0.012
!       !print*,'cflux:',cflux(i),gah2o(i),co2a(i),co2m(i)

!    ENDDO

!    !
!    dtti = 1./dtt
!    DO i = 1,len        ! LOOP OVER GRIDPOINT
!       !czzggrst5 - new code
!       H2OI   = ETC(i)/PSUR(i)
!       H2OA   =  EA(i)/PSUR(i)
!       ECMOLE = 55.56 * ECMASS(i) * dtti  ! ecmass must be computed and passed in
!       H2OS = H2OA + ECMOLE / GBH2O(i)
!       H2OS  = MIN( H2OS, H2OI )
!       H2OS  = MAX( H2OS, 1.E-7)
!       H2OSRH = H2OS/H2OI
!       !  need qdamp and pdamp calculated and passed to here!
!       !pl        CO2S(i) = MAX(PCO2S(I),PCO2M*0.5) / (PSUR(i)*100.)

!       !pl I have relaxed this condition to 1/10 of previous. The old way made
!       !pl the CO2 on top of the leaves always at least 1/2 of the value at the
!       !pl reference level.

!       CO2S(i) = MAX(PCO2S(I),PCO2M(i)*0.05) / (PSUR(i)*100.)

!       !pl Ball-Berry equation right here !

!  !     print'(a,6g12.4)','ball-berry:',gradm(i),assimn(i),h2osrh,soilfrz(i),  &
!  !                          co2s(i),bintc(i)

!       GSH2OINF = (GRADM(i) * MAX(1.E-12,ASSIMN(i))            &
!            * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)

!       !pl this is the change in stomatal resistance

!       DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/          &
!            (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))

!  !     print'(a,6g16.6)','drst',rst(i),drst(i),qdamp,gsh2o(i),gsh2oinf,pdamp


!       !pl this is the 'would be change' if we did not use the damping factor..

!       !        rstnew = (1./gsh2oinf) * 44.6 * tprcor(i)/tc(i)
!       !        DRST(i) = rstnew - RST(i)

!       !
!       RSTFAC(i,1) = H2OS/H2OI
!       RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)
!    ENDDO
!    !
!    !Z CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
!    !-----------------------------------------------------------------------
!    !  INPUTS: PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
!    !    C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
!    !    TEMPL,RSTFAC(i,1),VM(i),ASSIM,GSH20(i),EFFCON(i),QT,GAMMAS(i),
!    !    PFD(i)
!    !

!    sttp2 = 73.**0.2
!    ohtp2 = 100.**0.2
!    DO i = 1,len
!       !-----------------------------------------------------------------------
!       ! CALCULATION OF POTENTIAL ASSIMILATION
!       !-----------------------------------------------------------------------
!       !
!       ! Make assimn a top leaf, not the canopy.
!       ASSIMNp(i) = ASSIMN(i) / APARKK(i)
!       !
!       ! Bottom stopped assim.
!       ANTEMP(i) = MAX(0.E0,ASSIMNp(i))
!       !
!       ! Potential intercellular co2.
!       PCO2IPOT = PSUR(i)*100.*(co2s(i)-(1.6*ASSIMNp(i)/         &
!            ((GRADM(i)*ANTEMP(i)/co2s(i))+BINTC(i))))
!       !
!       ! Potential rubisco limitation.
!       OMCPOT = VMAX0(i)*2.1**qt(i)*((PCO2IPOT-GAMMAS(i))/        &
!            (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
!       !
!       ! Potential light limitation.
!       OMEPOT(i) = PAR(i)*((PCO2IPOT-GAMMAS(i))/                 &
!            (PCO2IPOT+2.*GAMMAS(i))*C3(i) + C4(i))
!       !
!       ! Quad 1.
!       SQRTIN = MAX(0.E0,((OMEPOT(i)+OMCPOT)**2-                &
!            4.*ATHETA(i)*OMEPOT(i)*OMCPOT))
!       !
!       ! Quad 1. Intermediate  top leaf photosynthesis.
!       OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.*ATHETA(i))
!       !
!       ! Potential sink or pep limitation.
!       OMSPOT = (VMAX0(i)/2.0)*(1.8**qt(i))*C3(i)                  &
!            + RRKK(i)*PCO2IPOT*C4(i)
!       !
!       ! Quad 2.
!       SQRTIN=MAX(0.E0,((OMPPOT+OMSPOT)**2-4.*BTHETA(i)*          &
!            OMPPOT*OMSPOT))
!       !
!       ! Quad 2. Final Potential top leaf photosynthesis.
!       ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.*BTHETA(i))
!       !
!       !-----------------------------------------------------------------------
!       ! CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
!       !-----------------------------------------------------------------------
!       !
!       ! Stressed rubisco limitation.
!       OMCCI = VM(i)*((PCO2IPOT-GAMMAS(i))/(PCO2IPOT+RRKK(i))*C3(i)    &
!            + C4(i))
!       !
!       ! Quad 1.
!       SQRTIN = MAX(0.E0,(OMEPOT(i)+OMCCI)**2-          &
!            4.*ATHETA(i)*OMEPOT(i)*OMCCI)
!       !
!       ! Quad 1. Intermediate stress limited top leaf photosynthesis.
!       OMPCI = ((OMEPOT(i)+OMCCI)-SQRT(SQRTIN))/(2.*ATHETA(i))
!       !
!       ! Stressed sink or pep limitation.
!       OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
!       !
!       ! Quad 2.
!       SQRTIN = MAX(0.E0,(OMPCI+OMSCI)**2-4.*BTHETA(i)*OMPCI*OMSCI)
!       !
!       ! Quad 2. Final stress limited top leaf photosynthesis.
!       ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.*BTHETA(i))
!       !
!       !-----------------------------------------------------------------------
!       ! CALCULATION OF CONTROL COEFFICIENTS
!       !-----------------------------------------------------------------------
!       !
!       ! Intermediate.
!       DOMPDOMC = (OMPCI-OMEPOT(i))/                       &
!            (2.*ATHETA(i)*OMPCI-OMCCI-OMEPOT(i))
!       !
!       ! Bottom stopped final stress limited top leaf photosynthesis.
!       ASCITEMP = MAX(ASSIMCI(i),1.E-12)
!       !
!       ! Rubisco control coefficient.
!       CCOMC = (DOMPDOMC*(ASSIMCI(i)-OMSCI)/                        &
!            (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMCCI/ASCITEMP
!       !
!       ! Sink or pep control coefficient.
!       CCOMS = ((ASSIMCI(i)-OMPCI)/            &
!            (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMSCI/ASCITEMP
!       !
!       !-----------------------------------------------------------------------
!       !  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
!       !-----------------------------------------------------------------------
!       ! Canopy values (overwrites top leaf).
!       !
!       OMEPOT(i) = OMEPOT(i)*APARKK(i)
!       ASSIMPOT(i) = ASSIMPOT(i)*APARKK(i)
!       ASSIMCI(i) = ASSIMCI(i)*APARKK(i)
!       ASSIM(i) = ASSIM(i)*APARKK(i)
!       ANTEMP(i) = ANTEMP(i)*APARKK(i)
!       ANSQR(i) = ANTEMP(i)*ANTEMP(i)
!       ASSIMNp(i) = ASSIMNp(i)*APARKK(i)
!       !
!       !-----------------------------------------------------------------------
!       ! OUTPUT:WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
!       !-----------------------------------------------------------------------
!       !
!       ! Water stress.
!       WSFWS(i) = ASSIMPOT(i)*(1.-RSTFAC(i,2))*(CCOMC+CCOMS)
!       !
!       ! High temperature stress.
!       WSFHT(i) = ASSIMPOT(i)*(1.-1./TEMPH(i))*CCOMC
!       !
!       ! Low temperature stress.
!       WSFLT(i) = ASSIMPOT(i)*(1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
!       !
!       !  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
!       cwsfws = (1.-RSTFAC(i,2))*(CCOMC+CCOMS)
!       IF(cwsfws.GT.1. .OR. cwsfws.LT.0.) wsfws(i)=0.
!       cwsfht = (1.-1./TEMPH(i))*CCOMC
!       IF(cwsfht.GT.1. .OR. cwsfht.LT.0.) wsfht(i)=0.
!       cwsflt = (1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
!       IF(cwsflt.GT.1. .OR. cwsflt.LT.0.) wsflt(i)=0.

!       !
!       ! Intermediate assimilation weighted Ci.
!       WCI(i) = ANTEMP(i)*PCO2I(i)
!       !
!       ! Intermediate assimilation weighted relative humidty stress factor.
!       WHS(i) = ANTEMP(i)*RSTFAC(i,1)
!       !
!       ! Intermediate assimilation weighted stomatal conductance.
!       WAGS(i) = GSH2O(i)*ANTEMP(i)
!       !
!       ! Intermediate evaporation weighted stomatal conductance.(Step 1.
!       !   Step 2 after subroutine updat2)
!       WEGS(i) = GSH2O(i)
!       !
!    ENDDO
!    !

!    RETURN
!  END SUBROUTINE PHOSIB

!  !----------------------------------------------------------------

!  !
!  !===================SUBROUTINE SORTIN===================================
!  !
!  SUBROUTINE SORTIN( EYY, PCO2Y, RANGE, GAMMAS, IC,len )
!    !
!    !=======================================================================
!    !
!    !     ARRANGES SUCCESSIVE PCO2/ERROR PAIRS IN ORDER OF INCREASING PCO2.
!    !       ESTIMATES NEXT GUESS FOR PCO2 USING COMBINATION OF LINEAR AND
!    !       QUADRATIC FITS.
!    !
!    !=======================================================================
!    !
!    IMPLICIT NONE

!    INTEGER len, ic
!    REAL EYY(len,6), PCO2Y(len,6)   &
!         ,    RANGE(len),gammas(len)

!    !     work arrays

!    REAL eyyi1(len),eyyi2(len),eyyi3(len),eyyis(len),    &
!         eyyisp(len),pco2yis(len),pco2yisp(len),         &
!         pco2b(len),                                     &
!         pco2yi1(len),pco2yi2(len),pco2yi3(len)
!    REAL aterm, bterm, cterm, pco2yq, cc1, cc2, bc1, bc2, ac1, ac2,  &
!         pco2yl, a, b, pmin, emin, one
!    INTEGER is(len)
!    LOGICAL bitx(len)
!    INTEGER i, ix, i1, i2, i3, isp, n, l, j

!    !
!    one = 1.0
!    IF( IC .LT. 4 ) THEN
!       DO i = 1,len
!          PCO2Y(i,1) = GAMMAS(i) + 0.5*RANGE(i)
!          PCO2Y(i,2) = GAMMAS(i)                                    &
!               + RANGE(i)*( 0.5 - 0.3*SIGN(one,EYY(i,1)) )
!          PCO2Y(i,3) = PCO2Y(i,1)- (PCO2Y(i,1)-PCO2Y(i,2))          &
!               /(EYY(i,1)-EYY(i,2)+1.E-10)*EYY(i,1)
!          !
!          PMIN = MIN( PCO2Y(i,1), PCO2Y(i,2) )
!          EMIN = MIN(   EYY(i,1),   EYY(i,2) )
!          IF ( EMIN .GT. 0. .AND. PCO2Y(i,3) .GT. PMIN )            &
!               PCO2Y(i,3) = GAMMAS(i)
!       ENDDO
!    ELSE
!       !
!       N = IC - 1
!       DO l = 1,len
!          bitx(l) = ABS(eyy(l,n)).GT.0.1
!          IF(.NOT.bitx(l)) pco2y(l,ic) = pco2y(l,n)
!       ENDDO
!       DO l = 1,len
!          IF(bitx(l)) THEN
!             DO J = 2, N !DO 1000 J = 2, N
!                A = EYY(l,J)
!                B = PCO2Y(l,J)
!                DO I = J-1,1,-1 !DO 2000 I = J-1,1,-1
!                   IF(EYY(l,I) .LE. A ) GO TO 100
!                   EYY(l,I+1) = EYY(l,I)
!                   PCO2Y(l,I+1) = PCO2Y(l,I)
!                ENDDO !2000             CONTINUE
!                i = 0
!  100           CONTINUE
!                EYY(l,I+1) = A
!                PCO2Y(l,I+1) = B
!             ENDDO !1000          CONTINUE
!          ENDIF
!       ENDDO
!       !
!       !-----------------------------------------------------------------------
!       !
!       DO l = 1,len
!          IF(bitx(l)) THEN
!             PCO2B(l) = 0.
!             IS(l)    = 1
!          ENDIF
!       ENDDO

!       DO IX = 1, N !DO 3000 IX = 1, N
!          DO l = 1,len
!             IF(bitx(l)) THEN
!                IF( EYY(l,IX) .LT. 0. )  THEN
!                   PCO2B(l) = PCO2Y(l,IX)
!                   IS(l) = IX
!                ENDIF
!             ENDIF
!          ENDDO
!       ENDDO !3000    CONTINUE
!       DO l = 1,len
!          IF(bitx(l)) THEN
!             I1 = IS(l)-1
!             I1 = MAX(1, I1)
!             I1 = MIN(N-2, I1)
!             I2 = I1 + 1
!             I3 = I1 + 2
!             ISP   = IS(l) + 1
!             ISP = MIN0( ISP, N )
!             IS(l) = ISP - 1
!             eyyisp(l) = eyy(l,isp)
!             eyyis(l) = eyy(l,is(l))
!             eyyi1(l) = eyy(l,i1)
!             eyyi2(l) = eyy(l,i2)
!             eyyi3(l) = eyy(l,i3)
!             pco2yisp(l) = pco2y(l,isp)
!             pco2yis(l) = pco2y(l,is(l))
!             pco2yi1(l) = pco2y(l,i1)
!             pco2yi2(l) = pco2y(l,i2)
!             pco2yi3(l) = pco2y(l,i3)
!          ENDIF
!       ENDDO
!       !
!       DO l = 1,len
!          IF(bitx(l)) THEN

!             !itb...Neil Suits' patch to check for zero in the denominator...
!             IF(EYYis(l) /= EYYisp(l))THEN
!                PCO2YL=PCO2Yis(l)                &
!                     - (PCO2Yis(l)-PCO2Yisp(l))    &
!                     /(EYYis(l)-EYYisp(l))*EYYis(l)
!             ELSE
!                PCO2YL = PCO2Yis(l) * 1.01
!             ENDIF
!             !
!             !   METHOD USING A QUADRATIC FIT
!             !
!             AC1 = EYYi1(l)*EYYi1(l) - EYYi2(l)*EYYi2(l)
!             AC2 = EYYi2(l)*EYYi2(l) - EYYi3(l)*EYYi3(l)
!             BC1 = EYYi1(l) - EYYi2(l)
!             BC2 = EYYi2(l) - EYYi3(l)
!             CC1 = PCO2Yi1(l) - PCO2Yi2(l)
!             CC2 = PCO2Yi2(l) - PCO2Yi3(l)

!             !itb...Neil Suits' patch to prevent zero in denominator...
!             IF(BC1*AC2-AC1*BC2 /= 0.0 .AND. AC1 /= 0.0)THEN
!                BTERM = (CC1*AC2-CC2*AC1)/(BC1*AC2-AC1*BC2)
!                ATERM = (CC1-BC1*BTERM)/AC1
!                CTERM = PCO2Yi2(l)                              &
!                     -ATERM*EYYi2(l)*EYYi2(l)-BTERM*EYYi2(l)
!                PCO2YQ= CTERM
!                PCO2YQ= MAX( PCO2YQ, PCO2B(l) )
!                PCO2Y(l,IC) = ( PCO2YL+PCO2YQ)/2.
!             ELSE
!                PCO2Y(l,IC) = PCO2Y(l,IC) * 1.01
!             ENDIF

!          ENDIF
!       ENDDO
!       !
!    ENDIF
!    DO i = 1,len
!       pco2y(i,ic) = MAX(pco2y(i,ic),0.01E0)
!    ENDDO
!    !
!    RETURN
!  END SUBROUTINE SORTIN

!  !------------------------------------------------------
!  !
!  !=====================SUBROUTINE CYCALC=================================
!  !

!  SUBROUTINE CYCALC( APARKK, VM, ATHETA, BTHETA, par,     &
!       GAMMAS, RESPC, RRKK, OMSS, C3, C4,   &
!       PCO2I, ASSIMN, assim, len )
!    !
!    !=======================================================================
!    !
!    !     CALCULATION EQUIVALENT TO STEPS IN FIGURE 4 OF SE-92A
!    !     C4 CALCULATION BASED ON CO-92.
!    !
!    !=======================================================================
!    !
!    IMPLICIT NONE

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
!    !       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
!    !       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
!    !
!    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!    !
!    !       OMC            RUBISCO LIMITED ASSIMILATION (MOL M-2 S-1)
!    !       OME            LIGHT LIMITED ASSIMILATION (MOL M-2 S-1)
!    !       OMS            SINK LIMITED ASSIMILATION (MOL M-2 S-1)
!    !       CO2S           CANOPY SURFACE CO2 CONCENTRATION (MOL MOL-1)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len
!    REAL aparkk(len),vm(len),atheta(len),    &
!         btheta(len),gammas(len),par(len),   &
!         respc(len),rrkk(len),omss(len),     &
!         c3(len),c4(len),pco2i(len),         &
!         assimn(len),assim(len)

!    !    local variables
!    REAL ome, omc, omp, oms, sqrtin
!    INTEGER i

!    !-----------------------------------------------------------------------
!    !     CALCULATE ASSIMILATION RATE
!    !
!    !      OMC         (OMEGA-C): EQUATION (11) , SE-92A
!    !      OME         (OMEGA-E): EQUATION (12) , SE-92A
!    !      OMS         (OMEGA-S): EQUATION (13) , SE-92A
!    !      ASSIMN      (A-N)    : EQUATION (14,15), SE-92A
!    !-----------------------------------------------------------------------

!    DO i = 1,len
!       OMC = VM(i) *(PCO2I(i)-GAMMAS(i))/(PCO2I(i) + RRKK(i))*C3(i)   &
!            + VM(i) * C4(i)
!       OME = PAR(i)*(PCO2I(i)-GAMMAS(i))/(PCO2I(i)+2.*GAMMAS(i))*C3(i) &
!            + PAR(i) * C4(i)
!       SQRTIN= MAX( 0.E0, ( (OME+OMC)**2 - 4.*ATHETA(i)*OME*OMC ) )
!       OMP  = ( ( OME+OMC ) - SQRT( SQRTIN ) ) / ( 2.*ATHETA(i) )
!       OMS  = OMSS(i) * C3(i) + OMSS(i)*PCO2I(i) * C4(i)
!       SQRTIN= MAX( 0.E0, ( (OMP+OMS)**2 - 4.*BTHETA(i)*OMP*OMS ) )
!       ASSIM(i) = ( ( OMS+OMP ) - SQRT( SQRTIN ) ) /    &
!            ( 2.*BTHETA(i) )
!       ASSIMN(i)= ( ASSIM(i) - RESPC(i)) * APARKK(i)

!    ENDDO
!    !
!    RETURN
!  END SUBROUTINE CYCALC

!  !---------------------------------------------------------------------

!  SUBROUTINE VMFCALZO(PR0,RIBC,VKRMN,DELTA,GRAV       &
!       ,                PS,tha                             &
!       ,                SPDM,ROS,CU,THVGM,RIB              &
!       ,                USTAR,zzwind,zztemp                &
!       ,                len )

!    IMPLICIT NONE

!    !****************************************************************************
!    !    VENTILATION MASS FLUX,Ustar, and transfer coefficients for momentum
!    !    and heat fluxes, based on by Louis (1979, 1982), and revised by Holtslag
!    !    and Boville(1993), and by Beljaars and Holtslag (1991).
!    !
!    !     Rerences:
!    !       Beljars and Holtslag (1991): Flux parameterization over land surfaces
!    !              for atmospheric models. J. Appl. Meteo., 30, 327-341.
!    !       Holtslag and Boville (1993): Local versus nonlocal boundary-layer
!    !              diffusion in a global climate model. J. of Climate, 6, 1825-
!    !              1842.
!    !       Louis, J. F., (1979):A parametric model of vertical eddy fluxes in
!    !              atmosphere. Boundary-Layer Meteo., 17, 187-202.
!    !       Louis, Tiedke, and Geleyn, (1982): A short history of the PBL
!    !              parameterization at ECMWF. Proc. ECMWF Workshop on Boundary-
!    !              Layer parameterization, ECMWF, 59-79.
!    !
!    !     General formulation:
!    !        surface_flux = transfer_coef.*U1*(mean_in_regerence - mean_at_sfc.)
!    !     Transfer coefficients for mommentum and heat fluxes are:
!    !        CU = CUN*Fm, and
!    !        CT = CTN*Fh
!    !        where  CUN and CTN are nutral values of momentum and heat transfers,
!    !           and Fm and Fh are stability functions derived from surface
!    !           similarity relationships.
!    !****************************************************************************

!    INTEGER len
!    REAL pr0, ribc, vkrmn, delta, grav               &
!         ,    eve(len),CU(len)                            &
!         ,    ROS(len),SH(len),PS(len)                    &
!         ,    THa(len),THVGM(len),RIB(len)                &
!         ,    SPDM(len),USTAR(len)                        &
!         ,    CUI(len)

!    !     local variables
!    REAL TEMV(len), wgm(len)                              &
!         ,      bunstablM, cunstablM, bstabl, cstabl           &
!         ,      zzwind(len), zztemp(len), zrib(len), cun(len)  &
!         ,      z1z0U(len), z1z0Urt(len)                       &
!         ,      fmomn(len),ribtemp,dm

!    INTEGER i
!    !
!    !  condtants for surface flux functions, according to Holtslag and
!    !      Boville (1993, J. Climate)

!    bunstablM = 10.       ! constants for unstable function
!    cunstablM = 75.
!    bstabl = 8.           ! constants for stable function
!    cstabl = 10.
!    !
!    DO I=1,len

!       zrib(i) = zzwind(i) **2 / zztemp(i)

!       !   Ratio of reference height (zwind/ztemp) and roughness length:
!       z1z0U(i) = zzwind(i)/ 0.0002   ! oceanic roughness length
!       z1z0Urt(i) = SQRT( z1z0U(i) )
!       z1z0U(i) = LOG( z1z0U(i) )

!       !   Neutral surface transfers for momentum CUN and for heat/moisture CTN:

!       cun(i) = VKRMN*VKRMN / (z1z0U(i)*z1z0U(i) )   !neutral Cm & Ct
!       !
!       ! SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
!       ! RIB IS THE BULK RICHARDSON NUMBER, between reference height and surface.
!       !
!       TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
!       temv(i) = MAX(0.000001E0,temv(i))
!       RIB(I) = -THVGM(I) * GRAV * zrib(i) / TEMV(i)
!    ENDDO

!    !   The stability functions for momentum and heat/moisture fluxes as
!    !   derived from the surface-similarity theory by Luis (1079, 1982), and
!    !   revised by Holtslag and Boville(1993), and by Beljaars and Holtslag
!    !   (1991).

!    DO I=1,len  !DO 60 I=1,len
!       IF(rib(i).GE.0.0) THEN
!          !
!          !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
!          !
!          rib(i) = MIN( rib(i), 0.5E0)
!          fmomn(i) = (1. + cstabl * rib(i) * (1.+ bstabl * rib(i)))
!          fmomn(i) = 1. / fmomn(i)
!          fmomn(i) = MAX(0.0001E0,fmomn(i))

!       ELSE
!          !
!          !        THE UNSTABLE CASE.
!          !
!          ribtemp = ABS(rib(i))
!          ribtemp = SQRT( ribtemp )
!          dm = 1. + cunstablM * cun(i) * z1z0Urt(i) * ribtemp
!          fmomn(i) = 1. - (bunstablM * rib(i) ) / dm

!       END IF
!    ENDDO !60   CONTINUE

!    !   surface-air transfer coefficients for momentum CU, for heat and
!    !   moisture CT. The CUI and CTI are inversion of CU and CT respectively.

!    DO i = 1, len
!       CU(i) = CUN(i) * fmomn(i)
!       CUI(i) = 1. / CU(i)

!       !   Ustar and ventlation mass flux: note that the ustar and ventlation
!       !   are calculated differently from the Deardoff's methods due to their
!       !   differences in define the CU and CT.

!       USTAR(i) = SPDM(i)*SPDM(i)*CU(i)
!       USTAR(i) = SQRT( USTAR(i) )
!    ENDDO
!    !
!    !  Note there is no CHECK FOR VENTMF EXCEEDS TOWNSENDS(1962) FREE CONVECTION
!    !  VALUE, like DEARDORFF EQ(40B), because the above CU and CT included
!    !  free convection conditions.
!    !

!    RETURN
!  END SUBROUTINE VMFCALZO

!  !--------------------------------------------------------------------------

!  SUBROUTINE VMFCALo(PR0,RIBC,VKRMN,GRAV        &
!       ,                SPDM,ZB,ROS,CU,THVGM         &
!       ,                USTAR, tha                   &
!       ,                len, thgeff, tke, dotkef )
!    IMPLICIT NONE
!    ! subroutine is stripped from vmfcal.F and is used only to calculate
!    !   ustart and cu for oceanic values of the surface roughness length
!    !
!    !***  VENTILATION MASS FLUX, BASED ON DEARDORFF, MWR, 1972
!    !
!    INTEGER len
!    LOGICAL dotkef
!    REAL pr0, ribc, vkrmn, grav               &
!         ,    CU(len), ZB(len),ROS(len)            &
!         ,    THVGM(len),RIB(len)                  &
!         ,    SPDM(len),USTAR(len)                 &
!         ,    tha(len), thgeff(len),tke(len)

!    !     local variables
!    REAL ribmax, vkrinv, TEMV(len), ZDRDRF(len), CHOKE(len)  &
!         ,    sqrtke, cui(len), cuni(len), cni(len)
!    INTEGER i
!    !
!    RIBMAX = 0.9*RIBC
!    VKRINV = 1.0/VKRMN
!    !
!    !        CUNI AND CTN1 ARE INVERSES OF THE NEUTRAL TRANSFER COEFFICIENTS
!    !        DEARDORFF EQS(33) AND (34).
!    !        PR0 IS THE TURBULENT PRANDTL NUMBER AT NEUTRAL STABILITY.
!    !
!    DO I=1,len  !DO 20 I=1,len
!       CNI(i) = 0.025*ZB(i)/0.0002    ! oceanic surface roughness length
!       cni(i) = LOG(cni(i))
!       CNI(i)  = CNI(i)  * VKRINV
!       CUNI(i) = CNI(i)  + 8.4
!       CNI(i)  = GRAV    * ZB(i)
!       !
!       !        SURFACE TO AIR DIFFERENCE OF POTENTIAL TEMPERATURE.
!       !        RIB IS THE BULK RICHARDSON NUMBER, DEARDORFF EQ (25)
!       !
!    ENDDO !20   CONTINUE
!    IF (dotkef) THEN
!       DO I=1,len
!          rib(i) = -thvgm(i)*grav*zb(i)/(thgeff(i)*tke(i))
!       ENDDO
!       DO I=1,len
!          IF(rib(i).GE.0.0) THEN
!             cu(i) = ((-0.0307*rib(i)**2)/(61.5+rib(i)**2))+0.044
!          ELSE
!             cu(i) = ((-0.016*rib(i)**2)/(4.2e4+rib(i)**2))+0.044
!          ENDIF
!       ENDDO
!       DO I=1,len
!          sqrtke = SQRT(tke(i))
!          USTAR(I) = SQRT(sqrtke*spdm(I)*CU(I))
!       ENDDO
!    ELSE
!       DO I=1,len
!          TEMV(i) = THA(i) * SPDM(i) * SPDM(i)
!          temv(i) = MAX(0.000001E0,temv(i))
!          RIB(I) = -THVGM(I) * CNI(I) / TEMV(i)
!       ENDDO
!       DO I=1,len  !DO 60 I=1,len
!          IF(rib(i).GE.0.0) THEN
!             !
!             !        THE STABLE CASE. RIB IS USED WITH AN UPPER LIMIT
!             !
!             CHOKE(i)=1.0-MIN(RIB(I),RIBMAX)/RIBC
!             CU(I)=CHOKE(i)/CUNI(I)
!          ELSE
!             !
!             !        FIRST, THE UNSTABLE CASE. DEARDORFF EQS(28), (29), AND (30)
!             !
!             ZDRDRF(i) = LOG10(-RIB(I)) - 3.5
!             CUI(I)    = CUNI(I) - 25.0 * EXP(0.26 * ZDRDRF(i)     &
!                  - 0.03 * ZDRDRF(i) * ZDRDRF(i))
!             CU(I)     = 1.0 / MAX(CUI(I),0.5 * CUNI(I))
!          END IF
!       ENDDO !60      CONTINUE
!       DO I=1,len !DO 80 I=1,len
!          USTAR(i) =SPDM(i)*CU(i)
!       ENDDO !80      CONTINUE
!    ENDIF
!    !
!    RETURN
!  END SUBROUTINE VMFCALo

!  !--------------------------------------------------------------------
!  !
!  !===================SUBROUTINE DELLWF=====================================
!  !
!  SUBROUTINE DELLWF(DTA,dtc4,dtg4,dts4,fac1,areas    &
!       ,                   lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc   &
!       ,                   len )

!    !========================================================================
!    !
!    !     Calculation of partial derivatives of canopy and ground radiative
!    !        heat fluxes with respect to Tc, Tgs
!    !pl   Here we are doing only the long wave radiative loss, which is the
!    !pl   only radiative quantity we are trying to bring to the next time step.
!    !
!    !========================================================================

!    !------------------------------INPUT is coming from Netrad-------------
!    !
!    !       dtc4, dtg4, dts4, which are the derivatives of the LW loss
!    !
!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       LCDTC          dLC/dTC
!    !       LCDTG          dLC/dTG
!    !       LCDTS          dLC/dTS
!    !       LGDTG          dLG/dTG
!    !       LGDTC          dLG/dTC
!    !       LSDTS          dLS/dTS
!    !       LSDTC          dLS/dTC
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    IMPLICIT NONE
!    !
!    INTEGER len
!    REAL                                                          &
!         lcdtc(len),lcdtg(len),lcdts(len),lgdtg(len),lgdtc(len)   &
!         ,    lsdts(len),lsdtc(len),areas(len),fac1(len)               &
!         ,    dtc4(len),dtg4(len),dts4(len),dta
!    !     local variables
!    INTEGER i

!    DO I=1,len

!       !pl canopy leaves:
!       LCDTC(I) =   2 * DTC4(i) * fac1(i)
!       LCDTG(I) =     - DTG4(i) * fac1(i) * (1.-areas(i))
!       LCDTS(I) =     - DTS4(i) * fac1(i) * (   areas(i))

!       !pl ground:
!       LGDTG(I) =   DTG4(i)
!       LGDTC(I) = - DTC4(i) * fac1(i)

!       !pl snow:
!       LSDTS(I) =   DTS4(i)
!       LSDTC(I) = - DTC4(i) * fac1(i)
!       !
!    ENDDO

!    RETURN
!  END SUBROUTINE DELLWF

!  !-----------------------------------------------------------------
!  !
!  !==================SUBROUTINE DELEF25======================================
!  !
!  SUBROUTINE DELEF(DTA,CP,ps,em,ea,ros,HRr,fc,fg         &
!       ,                ra,rb,rd,rc,rsoil,snow,capac,www      &
!       ,                ECDTC,ECDEA,EGDTG,EGDEA,ESDTS         &
!       ,                ESDEA,EADEA,EADEM                     &
!       ,                ec,eg,es,fws,hltm,cas_cap_vap         &
!       ,                etc,etg                               &
!       ,                btc,btg,bts                           &
!       ,                areas, gect,geci,gegs,gegi, psy, snofac, hr   &
!       ,                len)
!    !========================================================================
!    !
!    !     Calculation of partial derivatives of canopy and ground latent
!    !        heat fluxes with respect to Tc, Tgs, Theta-m, and Qm.
!    !     Calculation of initial latent heat fluxes.

!    !pl the ETC, ETG and so on are the vapor pressure at temps TC, TG and so on
!    !pl the BTC, BTG are the derivatives of ETC, ETG with relation to TC, TG etc.
!    !
!    !========================================================================

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       EC             ECT + ECI
!    !       EG             EGS + EGI
!    !       ECDTC          dEC/dTC
!    !       ECDTG          dEC/dTGS
!    !       ECDQM          dEC/dQM
!    !       EGDTC          dEG/dTC
!    !       EGDTG          dEG/dTGS
!    !       EGDQM          dEG/dQM
!    !       BBC            dE/dTC
!    !       BBG            dE/dTGS
!    !       BBM            dE/dQM
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    IMPLICIT NONE
!    !
!    INTEGER len
!    REAL                                                       &
!         ECDTC(len),ECDEA(len),EGDTG(len),EGDEA(len)           &
!         ,    ESDTS(len),ESDEA(len),EADEA(len),EADEM(len)           &
!         ,    fc(len),fg(len),snow(len,2),capac(len,2),www(len,3)   &
!         ,    hrr(len),rst(len),rb(len),rd(len),rc(len),rsoil(len)  &
!         ,    ra(len),rds(len), cas_cap_vap(len)                    &
!         ,    em(len),etc(len),etg(len),ea(len)                     &
!         ,    ps(len),ec(len),eg(len),es(len)               &
!         ,    btc(len),btg(len),bts(len)                            &
!         ,    sh(len),COG1(len),COC(len),D2(len)                    &
!         ,    tc_0(len), tgs_0(len), intg(len), intc(len)           &
!         ,    COG2(len), hr(len) ,fws(len),hltm, rcp(len), rcpg(len)&
!         ,    ros(len), dta, cp, psy(len), snofac                   &
!         ,    limci(len),limgi(len),limct(len),limgs(len)
!    REAL cogr(len),cogs(len),areas(len), cpdpsy(len)
!    REAL gect(len),geci(len),gegs(len),gegi(len),resrat

!    INTEGER i
!    !
!    !     MODIFICATION FOR SOIL DRYNESS : HR=REL. HUMIDITY IN TOP LAYER
!    !
!    resrat = 0.5

!    DO I=1,len
!       hrr(I) = HR(I)
!       IF(fg(i).LT.0.5) hrr(i) = 1.0
!    ENDDO

!    DO I=1,len
!       rcp(i)  = ros(i) * cp
!       rcpg(i) = rcp(i)/psy(i)              ! this is rho * cp / gamma
!       cpdpsy(i) = cp / psy(i)
!       rds(i) = rsoil(i) + rd(i)

!       !-----------------------------------------------------------------------
!       !
!       !     CALCULATION OF SURFACE RESISTANCE COMPONENTS, SEE EQUATIONS (64,66)
!       !       OF SE-86
!       !pl the ge?? coefficients come all the way from VNTLAT and are common to
!       !pl all subroutines:
!       !pl     gect(i)  =      (1. -wc(i)) /  rc(i)
!       !pl        geci(i)  = epsc(i) * wc(i)  / (RB(I) + RB(I))
!       !pl        gegs(i)  =        (1-wg(i)) / (rds(i))
!       !pl        gegi(i)  = epsg(i) * wg(i)  /  rd(i)
!       !
!       !-----------------------------------------------------------------------

!       COC(I) =  gect(i) + geci(i)
!       COG1(i) = (gegi(i) + gegs(i)*HRR(i))
!       COG2(i) = (gegi(i) + gegs(i)       )

!       !            D2(I)   = 1.0 / RA(I) + COC(I) + COG2(I)
!       !-----------------------------------------------------------------------
!       !
!       !     FLUXES EXPRESSED IN JOULES M-2   CPL WHY ?????
!       !
!       !      ec         (EC)    : EQUATION (64) , SE-86
!       !      eg         (EG)    : EQUATION (66) , SE-86
!       !      es         (ES)    : EQUATION (66) , SE-86
!       !      ea         (EA)    : EQUATION ????
!       !-----------------------------------------------------------------------

!       !pl these are the current time step fluxes in J/m2  WHY ?????

!       !pl notice that the fluxes are already limited by the altered e*(T) values

!       ec(I)  = (etc(I) - ea(i)) * COC(I) *       &
!            ros(i)  * dta * cpdpsy(i)

!       eg(I)  = ( etg(I) * COG1(I)                &
!            - ea(i) * COG2(I)) *             &
!            ros(i) * dta * cpdpsy(i)

!       es(I)  = ((etg(I) - ea(i))/rd(i) )*                     &
!            ros(i) * dta * cpdpsy(i)/snofac
!       fws(I) = ((ea(I)  - em(i) ) / ra(i))                    &
!            * ros(i) * dta * cpdpsy(i)

!       !pl now we do the partial derivatives  these assume W/m2

!       !pl for the canopy leaves vapor pressure: W/ (m2* K)
!       ECDTC(I) =    btc(I) * COC(I)                       &
!            * ros(i) * CPDPSY(i)
!       ECDEA(I) = - COC(I) * ros(i) * CPDPSY(i)

!       !pl for ground latent heat fluxes: W/ (m2* K)
!       EGDTG(I) =   btg(I) * COG1(I)                    &
!            * ros(i) * CPDPSY(i)
!       EGDEA(I) = - (COG2(I)) * ros(i) * CPDPSY(i)

!       !pl for snow latent heat fluxes: W/ (m2* K)
!       ESDTS(I) =   btg(I) * ros(i) * CPDPSY(i)/RD(i)

!       !pl for snow latent heat fluxes: W/ (m2 * Pa)
!       ESDEA(I) = - ros(i) * CPDPSY(i)/RD(i)

!       !pl for CAS latent heat fluxes: W/ (m2* Pa)
!       EADEA(I) = ros(i) * CPDPSY(i) / ra(i)
!       EADEM(I) = - EADEA(I)

!       !PL ATTENTION !!!! DANGER !!! do not use without sibdrv = true
!       !pl these all need to be re-done for the GCM (no sibdrv)
!       !-----------------------------------------------------------------------
!       !      BBC       (dE/dTC)  : EQUATION (13) , SA-89B
!       !      BBG       (dE/dTGS) : EQUATION (13) , SA-89B
!       !      BBM       (dE/dQM)  : EQUATION (13) , SA-89B
!       !-----------------------------------------------------------------------
!       !        BBG(I) = (COG1(I) / D2(i))
!       !     *          * btg(I) * 0.622 * ps(i)
!       !     *       / ((ps(i) - etg(I)) * (ps(i) - etg(I)))
!       !        BBC(I) = (COC(I)  / D2(i))
!       !     *            * btc(I) * 0.622 * ps(i)
!       !     *       / ((ps(i) - etc(I)) * (ps(i) - etc(I)))
!       !        BBM(I) = 1.0   / (ra(I)  * D2(i))
!       !
!    ENDDO
!    !

!    !      !print*,'delef: em,ea,etg,etc=',em,ea,etg,etc
!    !      !print*,'delef: coc,cog1,cog2,dta=',coc,cog1,cog2,dta
!    !      !print*,'delef: ros, cpdpsy=',ros,cpdpsy
!    !      !print*,'delef: ec,eg,es,fws=',ec,eg,es,fws
!    RETURN
!  END SUBROUTINE DELEF

!  !------------------------------------------------------------------------
!  !
!  !===================SUBROUTINE DELHF25=====================================
!  !
!  !pl ATTENTION receiving tgs instead of tg !!!
!  SUBROUTINE DELHF(DTA,CP,bps,tm,tg,ts,tc,ta,ros,ra,rb,rd      &
!       ,                HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA         &
!       ,                HADTA,HADTH                                 &
!       ,                hc, hg, hs, fss, len )

!    !========================================================================
!    !
!    !     Calculation of partial derivatives of canopy and ground sensible
!    !        heat fluxes with respect to Tc, Tgs, and Theta-m.
!    !     Calculation of initial sensible heat fluxes.
!    !
!    !========================================================================


!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       HC             CANOPY SENSIBLE HEAT FLUX (J M-2)
!    !       HG             GROUND SENSIBLE HEAT FLUX (J M-2)
!    !       HS             SNOW   SENSIBLE HEAT FLUX (J M-2)
!    !       HA             CAS    SENSIBLE HEAT FLUX (J M-2)
!    !       HCDTC          dHC/dTC
!    !       HCDTA          dHC/dTA
!    !       HGDTG          dHG/dTG
!    !       HGDTA          dHG/dTA
!    !       HSDTS          dHS/dTS
!    !       HSDTA          dHS/dTA
!    !       HADTA          dHA/dTA
!    !       HADTH          dHA/dTH
!    !       AAC            dH/dTC
!    !       AAG            dH/dTGS
!    !       AAM            dH/dTH
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    IMPLICIT NONE
!    !
!    INTEGER len
!    REAL                                                         &
!         bps(len),tm(len),tg(len),ts(len),HCDTC(len),HCDTA(len)  &
!         ,    HGDTG(len),HGDTA(len),HSDTS(len),HSDTA(len)             &
!         ,    HADTA(len),HADTH(len),ra(len)                           &
!         ,    rb(len),rd(len),tc(len), ta(len), hc(len), hg(len)      &
!         ,    ros(len), dta, cp, fss(len), hs(len)
!    !     local variables
!    REAL D1, d1i, rai, rbi,  rdi
!    INTEGER i

!    !-----------------------------------------------------------------------
!    !
!    !     FLUXES EXPRESSED IN JOULES M-2, although in SIBSLV WE THEN WANT W/m2
!    !pl WHY ????
!    !
!    !pl if we were to keep things simple, there is no need to separate
!    !pl HG and HS, but it helps the derivatives keep clean.
!    !
!    !      HC          (HC)    : EQUATION (63) , SE-86
!    !      HG          (HG)    : EQUATION (65) , SE-86
!    !      HS          (HS)    : EQUATION (65) , SE-86
!    !      HA          (HA)    : EQUATION ???
!    !-----------------------------------------------------------------------

!    DO I=1,len
!       rai = 1.0 / ra(I)
!       rbi = 1.0 / rb(I)
!       rdi = 1.0 / rd(I)
!       !        D1 = rai + rbi + rdi
!       !        d1i = 1. / d1

!       !pl these are the current time step fluxes in J/m2
!       !pl can we change this to W/m2 ???

!       HC(I)   = CP * ros(i) * (tc(i) - ta(i)) * rbi * DTA
!       HG(I)   = CP * ros(i) * (tg(I) - ta(i)) * rdi * DTA
!       HS(I)   = CP * ros(i) * (tg(I) - ta(i)) * rdi * DTA
!       fss(I)  = CP * ros(i) * (ta(I) - tm(i)) * rai * DTA

!       !pl now we do the partial derivatives

!       !pl these are done assuming the fluxes in W/m2

!       !pl for canopy leaves sensible heat flux: W/(m2 * K)
!       HCDTC(I) =   CP * ros(i) * rbi
!       HCDTA(I) = - HCDTC(I)
!       !pl for ground and snow sensible heat fluxes: W/(m2 * K)
!       HGDTG(I) =   CP * ros(i) * rdi
!       HSDTS(I) =   HGDTG(I)
!       HGDTA(I) = - HGDTG(I)
!       HSDTA(I) = - HGDTG(I)
!       !pl for the canopy air space (CAS) sensible heat flux: W/(m2 * K)
!       HADTA(I) =   CP * ros(i) * rai
!       HADTH(I) = - HADTA(I)/bps(i)

!       !pl ATTENTION !!!! DANGER !!!!! THIS WILL NOT WORK WITHOUT sibdrv = true
!       !pl for mixed layer (ref temp if not sibdrv): YET TO BE DONE
!       !        AAG(I) = rdi * d1i
!       !        AAC(I) = rbi * d1i
!       !        AAM(I) = rai * d1i * bps(i)
!       !
!    ENDDO

!    !itb...
!    !      !print*,'delhf: rbi,rdi,rai,dta=',rbi,rdi,rai,dta
!    !      !print*,'delhf: cp, ros=',cp,ros
!    !      !print*,'delhf: tm,tc,tg,ta=',tm,tc,tg,ta
!    !      !print*,'delhf: hc,hg,hs,fss=',hc,hg,hs,fss

!    RETURN
!  END SUBROUTINE DELHF

!  !-----------------------------------------------------------------

!  SUBROUTINE soilprop( td, tg, slamda, shcap, www, poros, ztdep,  &
!       asnow, snoww, areas, tf, snomel, sandfrac,    &
!       nsib, len, nsoil )

!    IMPLICIT NONE

!    !     this subroutine calculates the soil thermal properties heat capacity
!    !         (shcap) and conductivity (slamda). Phase changes are incorporated
!    !         into the heat capacity over a range of -0.5C to 0.5C.
!    !     slamda(n) is the conductivity between layers n and n+1 / delta z
!    !     layer 1 is the deepest layer, layer nsoil is immediately below the
!    !         surface
!    !     treatment based on Bonan

!    !     argument list variables
!    INTEGER len, nsib, nsoil
!    REAL                                                       &
!         td(nsib,nsoil), tg(len), slamda(len,nsoil),         &
!         sandfrac(len),                                      &
!         shcap(len,nsoil), www(len,3), poros(len), snomel,   &
!         ztdep(nsib,nsoil), areas(len), snoww(len), asnow, tf

!    !     local variables
!    REAL                                                           &
!         shcapu, shcapf, shsnow, klamda(len,nsoil),zsnow(len),   &
!         kf, ku, ksnow, ksoil(len), kiw(len,3),kii(len,3),       &
!         shcap1, klamda1, klamda2
!    INTEGER i,n

!    !  snow density assumed 0.25 that of water to convert from Bonan's constants
!    !  for ksnow and max effective snow depth
!    ksnow = 0.085
!    DO i = 1,len
!       !         ksoil(i) = (8.8 * sandfrac(i) +
!       !     *                    2.92 * (1.-sandfrac(i)) )
!       !     *                           **(1.-poros(i))
!       ksoil(i) = 6.0**(1.-poros(i))
!    ENDDO
!    DO n = 1,3
!       DO i = 1,len
!          kiw(i,n) = 0.6**(www(i,n)*poros(i))
!          kii(i,n) = 2.2**(www(i,n)*poros(i))
!       ENDDO
!    ENDDO
!    !     heat capacity calculation and conductivity calculation

!    !     soil layers 1 and 2 ( deepest layers, use www(3) )
!    DO n = 1,2
!       DO i = 1,len
!          SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,3)*   &
!               POROS(i)) * 4.186E6
!          shcapf =  0.5 * (1. + poros(i) *            &
!               (www(i,3)-1.0)) * 4.186E6
!          ku = (ksoil(i)*kiw(i,3)-0.15)*www(i,3) + 0.15
!          kf = (ksoil(i)*kii(i,3)-0.15)*www(i,3) + 0.15
!          IF(td(i,n).GE.tf+0.5) THEN
!             shcap(i,n) = shcapu * ztdep(i,n)
!             klamda(i,n) = ku
!          ELSE IF(td(i,n).LE.tf-0.5) THEN
!             shcap(i,n) = shcapf * ztdep(i,n)
!             klamda(i,n) = kf
!          ELSE
!             shcap(i,n) = (0.5*(shcapu+shcapf) +            &
!                  poros(i)*www(i,3)*snomel) * ztdep(i,n)
!             klamda(i,n) = kf + (ku-kf)*(td(i,n)+0.5-tf)
!          ENDIF
!       ENDDO
!    ENDDO

!    !     soil layers 3,4,5 and nsoil ( intermediate layers, use www(2) )
!    DO n = 3,nsoil
!       DO i = 1,len
!          SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,2)*       &
!               POROS(i)) * 4.186E6
!          shcapf =  0.5 * (1. + poros(i) *                &
!               (www(i,2)-1.0)) * 4.186E6
!          ku = (ksoil(i)*kiw(i,2)-0.15)*www(i,2) + 0.15
!          kf = (ksoil(i)*kii(i,2)-0.15)*www(i,2) + 0.15
!          IF(td(i,n).GE.tf+0.5) THEN
!             shcap(i,n) = shcapu * ztdep(i,n)
!             klamda(i,n) = ku
!          ELSE IF(td(i,n).LE.tf-0.5) THEN
!             shcap(i,n) = shcapf * ztdep(i,n)
!             klamda(i,n) = kf
!          ELSE
!             shcap(i,n) = (0.5*(shcapu+shcapf) +           &
!                  poros(i)*www(i,2)*snomel) * ztdep(i,n)
!             klamda(i,n) = kf + (ku-kf)*(td(i,n)+0.5-tf)
!          ENDIF
!       ENDDO
!    ENDDO

!    !     soil layer nsoil ( top layer, use www(1) )
!    !     if snow covered add additional heat capacity due to snow
!    DO i = 1,len
!       SHCAPu  = ( 0.5*(1.-POROS(i)) + WWW(i,1)*       &
!            POROS(i)) * 4.186E6
!       shcapf =  0.5 * (1. + poros(i) *                &
!            (www(i,1)-1.0)) * 4.186E6
!       ku = (ksoil(i)*kiw(i,1)-0.15)*www(i,1) + 0.15
!       kf = (ksoil(i)*kii(i,1)-0.15)*www(i,1) + 0.15
!       IF(td(i,nsoil).GE.tf+0.5) THEN
!          shcap1 = shcapu
!          klamda1 = ku
!       ELSE IF(td(i,nsoil).LE.tf-0.5) THEN
!          shcap1 = shcapf
!          klamda1 = kf
!       ELSE
!          shcap1 = 0.5*(shcapu+shcapf) +                &
!               poros(i)*www(i,1)*snomel
!          klamda1 = kf + (ku-kf)*(td(i,nsoil)+0.5-tf)
!       ENDIF
!       IF(tg(i).GE.tf+0.5) THEN
!          klamda2 = ku
!       ELSE IF(tg(i).LE.tf-0.5) THEN
!          klamda2 = kf
!       ELSE
!          klamda2 = kf + (ku-kf)*(tg(i)+0.5-tf)
!       ENDIF
!       zsnow(i) = MIN( MAX(1./asnow,snoww(i)) - 0.02, 0.25E0 )
!       shsnow = (0.5 * 4.186e6 * zsnow(i) + shcap1 * 0.02) /         &
!            (zsnow(i)+0.02)
!       shcap(i,nsoil) =  shcap(i,nsoil) * (1.-areas(i)) + areas(i) *   &
!            shcap(i,nsoil)*shsnow*(ztdep(i,nsoil)+zsnow(i)+0.02) /       &
!            (shsnow*ztdep(i,nsoil)+shcap(i,nsoil) *                      &
!            (zsnow(i)+0.02)) * (ztdep(i,nsoil)+zsnow(i)+0.02)
!       klamda1 = ksnow*klamda1 * (0.02+zsnow(i)) /                 &
!            (ksnow*0.02 + klamda1*zsnow(i))

!       !    soil conductivities / delta z

!       slamda(i,nsoil) = (1.-areas(i)) * 2.*klamda(i,nsoil)*klamda2 /   &
!            (klamda2*ztdep(i,nsoil)+0.02*klamda(i,nsoil)) +              &
!            areas(i) * klamda(i,nsoil)*klamda1 /                &
!            (klamda1*ztdep(i,nsoil)*0.05 +                               &
!            klamda(i,nsoil)*(zsnow(i)+0.02))

!    ENDDO

!    DO n = 1,nsoil-1
!       DO i = 1,len
!          slamda(i,n) = 2.0 * klamda(i,n)*klamda(i,n+1) /          &
!               (klamda(i,n)*ztdep(i,n+1)+                  &
!               klamda(i,n+1)*ztdep(i,n))
!       ENDDO
!    ENDDO

!    !DO i = 1,len
!    !ENDDO

!    RETURN
!  END SUBROUTINE soilprop

!  !------------------------------------------------------------------
!  !
!  !======================SUBROUTINE SIBSLV25=================================
!  !
!  !pl ATTENTION !!! I am calling this wiht the TGS temp instead of tg
!  SUBROUTINE SIBSLV(DT,GRAV2,CP,HLTM,tg,tsnow,td,slamda,pi         &
!       ,                 areas,fac1                                     &
!       ,                 VENTMF,PSB,BPS,ros,psy                         &
!       ,                 cg,ccx,cas_cap_heat,cas_cap_vap                &
!       ,                 lcdtc,lcdtg,lcdts,lgdtg,lgdtc,lsdts,lsdtc      &
!       ,                 HCDTC,HCDTA,HGDTG,HGDTA,HSDTS,HSDTA            &
!       ,                 HADTA,HADTH                                    &
!       ,                 hc, hg, hs, fss                                &
!       ,                 ECDTC,ECDEA,EGDTG,EGDEA,ESDTS                  &
!       ,                 ESDEA,EADEA,EADEM                              &
!       ,                 ec,eg,es,fws                                   &
!       ,                 etc,etg,ets                                    &
!       ,                 btc,btg,bts                                    &
!       ,                 RADT                                           &
!       ,                 dtc, dtg, dth, dqm, dta, dea                   &
!       ,                 len, sibdrv, forcerestore  )
!    !========================================================================
!    !
!    !     Calculation of time increments in Tc, Tgs, Theta-m and Qm using an
!    !        implicit backwards method with explicit coefficients.
!    !pl   Similar to equations (10-15), SA-92B.
!    !
!    !     Longwave feedbacks are now really included

!    !pl ATTENTION !!! this is hardwired to work only with Bonan soil.
!    !pl if force-restore is wanted, we need to include approporiate if loops like
!    !pl in sibslv.F
!    !========================================================================

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       DTC            CANOPY TEMPERATURE INCREMENT (K)
!    !       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
!    !       DTH            MIXED LAYER POTENTIAL TEMPERATURE INCREMENT (K)
!    !       DQM            MIXED LAYER MIXING RATIO INCREMENT (KG KG-1)
!    !       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
!    !       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    IMPLICIT NONE
!    !
!    ! Declare parameters
!    INTEGER, PARAMETER :: sgl = SELECTED_REAL_KIND(p=6)   ! Single
!    INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(p=13)  ! Double

!    INTEGER len
!    LOGICAL sibdrv, forcerestore,ipvt

!    REAL                                                           &
!         lcdtc(len),lcdtg(len),lcdts(len),lgdtg(len),lgdtc(len)    &
!         ,    lsdts(len),lsdtc(len),areas(len)                          &
!         ,    dtc4(len),dtg4(len),dts4(len)
!    REAL                                                      &
!         VENTMF(len),PSB(len),BPS(len)                        &
!         ,    cg(len),ccx(len),z2(len),z1(len)                     &
!         ,    HCDTC(len),HCDTA(len)                                &
!         ,    HGDTG(len),HGDTA(len),HSDTS(len),HSDTA(len)          &
!         ,    HADTA(len),HADTH(len)                                &
!         ,    ECDTC(len),ECDEA(len),EGDTG(len),EGDEA(len)          &
!         ,    ESDTS(len),ESDEA(len),EADEA(len),EADEM(len)          &
!         ,    AAG(len),AAC(len)                                    &
!         ,    AAM(len),BBG(len),BBC(len),BBM(len)                  &
!         ,    FSS(len),FWS(len),td(len),RADT(len,3)                &
!         ,    tc(len), tg(len),tsnow(len), ta(len)                 &
!         ,    ec(len),eg(len),es(len),ea(len)                      &
!         ,    etc(len),etg(len),ets(len)                           &
!         ,    btc(len),btg(len),bts(len)                           &
!         ,    hc(len), hg(len), hs(len)                            &
!         ,    dtc(len), dtg(len,2), dta(len), dea(len)             &
!         ,    dth(len), dqm(len), fac1(len)           &
!         ,    dt, hltm, cp, grav2, slamda(len), pi                 &
!         ,    cas_cap_heat(len),cas_cap_vap(len),ros(len)

!    !     local variables
!    REAL timcon, tmcn2, dti, ddthl, cpdt, TEMV(len)           &
!         , cpdpsy(len), psy(len)                                   &
!         ,     AVEC(len,7,8),BVEC(len,7)
!    REAL (KIND=8) ::  dAVEC(len,7,8), dBVEC(len,7)


!    INTEGER i,j,k, error_flag

!    !pl  this routine sets up the coupled system of partial
!    !pl  differential equations described in Sato et al.,
!    !pl  with the exception that now Ta and ea are prognostic
!    !pl  variables, and so we have two more equations, reflected
!    !pl  in two more lines and two more columns as compared
!    !pl  to the old sibslv.F (used for no prognistic CAS calculations)


!    !pl          J: /variables
!    !pl J: equation/  1     2     3     4      5     6     7     8
!    !pl              TC,   TG,   TS, TREF,  EREF,   TA,   EA,  FORCING past t..l.
!    !pl 1: TC
!    !pl 2: TG
!    !pl 3: TS
!    !pl 4: TREF
!    !pl 5: EREF
!    !pl 6: TA
!    !pl 7: EA


!    TIMCON = PI  / 86400.0
!    DTI  = 1.0 / DT
!    !
!    DO I=1,len
!       TEMV(I) = GRAV2 * VENTMF(i)
!       cpdpsy(i) = cp / psy(i)
!       !
!       !     DTC EQUATION
!       !
!       AVEC(I,1,1) = ccx(I) * DTI + HCDTC(I) + ECDTC(I) + lcdtc(i)
!       AVEC(I,1,2) = LCDTG(i)
!       AVEC(I,1,3) = LCDTS(i)
!       AVEC(I,1,4) = 0.
!       AVEC(I,1,5) = 0.
!       AVEC(I,1,6) = HCDTA(I)
!       AVEC(I,1,7) = ECDEA(I)
!       AVEC(I,1,8) = RADT(i,1) - HC(I) * DTI - ec(I) * DTI
!       !
!       !     DTG EQUATION
!       !
!       AVEC(i,2,1) = lgdtc(i)
!       AVEC(i,2,2) = cg(i)   * DTI + HGDTG(I) + EGDTG(I)         &
!            + lgdtg(i) + slamda(i)
!       AVEC(i,2,3) = 0.
!       AVEC(I,2,4) = 0.
!       AVEC(I,2,5) = 0.
!       AVEC(I,2,6) = hgdta(i)
!       AVEC(i,2,7) = egdea(i)
!       AVEC(I,2,8) = RADT(i,2) - HG(I) * DTI - eg(I) * DTI         &
!            - slamda(i) * (tg(I) - td(i))
!       !
!       !     DTS EQUATION
!       !
!       AVEC(i,3,1) = lsdtc(i)
!       AVEC(i,3,2) = 0.
!       AVEC(i,3,3) = cg(i)   * DTI + HSDTS(I) + ESDTS(I)        &
!            + lsdts(i) + slamda(i)
!       AVEC(I,3,4) = 0.
!       AVEC(I,3,5) = 0.
!       AVEC(I,3,6) = hsdta(i)
!       AVEC(i,3,7) = esdea(i)
!       AVEC(I,3,8) = RADT(i,3) - HS(I) * DTI - es(I) * DTI       &
!            - slamda(i) * (tg(I) - td(i))
!       !
!       !     DTA EQUATION
!       !
!       AVEC(i,6,1) = - HCDTC(i)
!       AVEC(i,6,2) = - HGDTG(i) * (1.-areas(i))
!       AVEC(i,6,3) = - HSDTS(i) * (   areas(i))
!       AVEC(I,6,4) =   HADTH(i)
!       AVEC(I,6,5) =   0.
!       AVEC(I,6,6) = cas_cap_heat(i)   * DTI                &
!            + HADTA(I)  - HCDTA(i)                &
!            - (1.-areas(i))*HGDTA(I) - areas(i)*HSDTA(I)
!       AVEC(i,6,7) = 0.
!       AVEC(I,6,8) = HC(I) * DTI - FSS(I) * DTI         &
!            + (1.-areas(i))*HG(I) * DTI        &
!            +     areas(i) *HS(I) * DTI
!       !
!       !     DEA EQUATION
!       !
!       AVEC(i,7,1) = - ECDTC(i)
!       AVEC(i,7,2) = - EGDTG(i) * (1.-areas(i))
!       AVEC(i,7,3) = - ESDTS(i) * (   areas(i))
!       AVEC(I,7,4) = 0.
!       AVEC(I,7,5) = EADEM(i)
!       AVEC(I,7,6) = 0.
!       AVEC(i,7,7) = cas_cap_vap(i)   * DTI            &
!            + (EADEA(I)  - ECDEA(I)          &
!            -  (1.-areas(i))*EGDEA(I)        &
!            -      areas(i) *ESDEA(I))
!       AVEC(I,7,8) = (EC(I) * DTI  -  FWS(I) * DTI      &
!            +  (1.-areas(i))*EG(I) * DTI      &
!            +      areas(i) *ES(I) * DTI)
!    ENDDO
!    !

!    IF(sibdrv) THEN
!       DO i=1,len
!          !
!          !     DTHETA EQUATION
!          !
!          AVEC(I,4,1) = 0.0
!          AVEC(I,4,2) = 0.0
!          AVEC(I,4,3) = 0.0
!          AVEC(I,4,4) = 1.0
!          AVEC(I,4,5) =  0.0
!          AVEC(I,4,6) =  0.0
!          AVEC(I,4,7) =  0.0
!          AVEC(I,4,8) =  0.0
!          !
!          !     DSH EQUATION
!          !
!          AVEC(I,5,1) = 0.0
!          AVEC(I,5,2) = 0.0
!          AVEC(I,5,3) =  0.0
!          AVEC(I,5,4) =  0.0
!          AVEC(I,5,5) = 1.0
!          AVEC(I,5,6) = 0.0
!          AVEC(I,5,7) = 0.0
!          AVEC(I,5,8) = 0.0
!       ENDDO
!    ELSE
!       DO i=1,len
!          !
!          !     DTHETA EQUATION
!          !
!          AVEC(I,4,1) = -temv(I) * AAG(I) * (1.-areas(i))
!          AVEC(I,4,2) = -temv(I) * AAG(I) * areas(i)
!          AVEC(I,4,3) = -TEMV(I) * AAC(I)
!          AVEC(I,4,4) = -TEMV(I) * (AAM(I) - 1.0) +          &
!               PSB(i) * DTI
!          AVEC(I,4,5) =  0.0
!          !
!          !     DSH EQUATION
!          !
!          AVEC(I,5,1) = -TEMV(I)*BBG(I) * (1.-areas(i))
!          AVEC(I,5,2) = -TEMV(I)*BBG(I) * areas(i)
!          AVEC(I,5,3) = -TEMV(I)*BBC(I)
!          AVEC(I,5,4) =  0.0
!          AVEC(I,5,5) = -TEMV(I) * (BBM(I) - 1.0) +             &
!               PSB(i) * DTI
!          AVEC(I,4,6) = GRAV2 * FSS(i)
!          AVEC(I,5,6) = GRAV2 * FWS(i)
!       ENDDO
!    ENDIF

!    !
!    !     SOLVE 7 X 8 MATRIX EQUATION
!    !

!    CALL GAUSS(len,AVEC,7,8,BVEC)
!    !
!    DO I=1,len
!       !
!       DTC(I)   = BVEC(I,1)
!       DTG(I,1) = BVEC(I,2)   ! this is DTG
!       DTG(I,2) = BVEC(I,3)   ! this is DTS
!       DTH(i)   = BVEC(I,4)
!       DQM(i)   = BVEC(I,5)
!       DTA(i)   = BVEC(I,6)
!       DEA(i)   = BVEC(I,7)
!    ENDDO
!    RETURN
!  END SUBROUTINE SIBSLV

!  !----------------------------------------------------------------
!  !
!  !===================SUBROUTINE GAUSS=====================================
!  !
!  SUBROUTINE GAUSS(len,WORK,N,NP1,X)
!    IMPLICIT NONE
!    !========================================================================
!    !
!    !     SOLVE A LINEAR SYSTEM BY GAUSSIAN ELIMINATION.  DEVELOPED BY
!    !     DR. CHIN-HOH MOENG.  A IS THE MATRIX OF COEFFICIENTS, WITH THE
!    !     VECTOR OF CONSTANTS APPENDED AS AN EXTRA COLUMN.  X IS THE VECTOR
!    !     CONTAINING THE RESULTS.  THE INPUT MATRIX IS NOT DESTROYED.
!    !
!    !========================================================================
!    !
!    INTEGER len, n, np1
!    REAL WORK(len,N,NP1),X(len,N)
!    !     local variables
!    REAL TEMV(len,2)
!    INTEGER ii, j, i, k, kk, l
!    !
!    DO II=2,N !DO 20 II=2,N
!       DO J=II,N !DO 20 J=II,N
!          DO I=1,len !DO 10 I=1,len
!             TEMV(I,1) = WORK(I,J,II-1) / WORK(I,II-1,II-1)
!          ENDDO !10 CONTINUE
!          DO K=1,NP1 !DO 20 K=1,NP1
!             DO I=1,len !DO 20 I=1,len
!                WORK(I,J,K) = WORK(I,J,K) - TEMV(I,1) * WORK(I,II-1,K)
!             ENDDO
!          ENDDO
!       ENDDO
!    ENDDO
!    !20 CONTINUE
!    !
!    DO K=N,2,-1 !DO 45 K=N,2,-1
!       DO I=1,len !DO 30 I=1,len
!          TEMV(I,1) = WORK(I,K,NP1) / WORK(I,K,K)
!       ENDDO !30   CONTINUE
!       KK = K-1
!       DO L=KK,1,-1 !DO 45 L=KK,1,-1
!          DO I=1,len !DO 40 I=1,len
!             TEMV(I,2) = TEMV(I,1) * WORK(I,L,K)
!          ENDDO !40      CONTINUE
!          DO I=1,len !DO 45 I=1,len
!             WORK(I,L,NP1) = WORK(I,L,NP1) - TEMV(I,2)
!          ENDDO
!       ENDDO
!    ENDDO !45 CONTINUE
!    !
!    DO II=1,N !DO 50 II=1,N
!       DO I=1,len !DO 50 I=1,len
!          X(I,II) = WORK(I,II,NP1) / WORK(I,II,II)
!       ENDDO
!    ENDDO !50 CONTINUE
!    !
!    RETURN
!  END SUBROUTINE GAUSS
!  !----------------------------------------------------------------------
!  !
!  !=====================SUBROUTINE UPDAT25=================================
!  !
!  SUBROUTINE UPDAT2(snoww,capac,snofac,ect,eci,egi            &
!       ,             egs,hlat,www,pi,cg,dtd,dtg,dtc,ta,dta,dea     &
!       ,             dtt                                           &
!       ,             roff,tc,td,tg,bee, poros                      &
!       ,             satco,slope,phsat,zdepth,ecmass               &
!       ,             egmass,shf,tf,snomel,asnow                    &
!       ,             ccx,csoil,chf,hc,hg,areas                     &
!       ,             q3l,q3o,qqq,zmelt,cw, len, nsib               &
!       ,             nsoil, forcerestore,etc,ea,btc                &
!       ,             geci,ros,cp,psy,gect,etg,btg                  &
!       ,             gegs,hr,fg,gegi,rd,rb,hcdtc,hcdta             &
!       ,             hgdta,slamda)

!    IMPLICIT NONE

!    !itb====================================================================
!    !itb
!    !itb  MOVING STORAGE TERM UPDATES (SNOW, CAPAC) HERE FROM ENDTEM, WHICH
!    !itb     NO LONGER EXISTS. FLUXES PREVIOUSLY CALCULATED IN ENDTEM ARE
!    !itb     TAKEN CARE OF IN THE PROGNOSTIC C.A.S. CALCULATIONS, SO WE
!    !itb     MERELY NEED TO TAKE CARE OF STORAGE TERMS NOW.
!    !itb
!    !itb      November 2000
!    !
!    !=======================================================================
!    !
!    !     UPDATING OF ALL HYDROLOGICAL PROGNOSTIC VARIABLES.  SNOW AND
!    !        RUNOFF CALCULATIONS (SEE ALSO INTER2).  SUBROUTINES SNOW2 AND
!    !        RUN2 OF 1D MODEL ARE INLINED IN THIS CODE.
!    !
!    !=======================================================================
!    !

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       DTC            CANOPY TEMPERATURE INCREMENT (K)
!    !       DTD            DEEP SOIL TEMPERATURE INCREMENT (K)
!    !       DTG            GROUND SURFACE TEMPERATURE INCREMENT (K)
!    !       WWW(3)         GROUND WETNESS
!    !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
!    !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
!    !       ROFF           RUNOFF (MM)
!    !
!    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!    !
!    !       ECMASS         CANOPY EVAPOTRANSPIRATION (MM)
!    !       EGMASS         GROUND EVAPOTRANSPIRATION (MM)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len, nsib, nsoil, i, l, iveg, ksoil,j

!    LOGICAL forcerestore

!    !
!    !-----------------------------------------------------------
!    !
!    !              INTENT = IN VARIABLES
!    !
!    !----------------------------------------------------------

!    REAL, INTENT(in),DIMENSION(len)                ::&
!         cg              &! surface layer heat capacity (J m^-2 deg^-1)
!         ,      dta            & ! delta canopy air space (CAS) temperature (K)
!         ,      dea            & ! delta CAS vapor pressure (Pa)
!         ,      tc             & ! canopy temperature (K)
!         ,      ta              &! CAS temperature (K)
!         ,      tg             & ! ground surface temperature (K)
!         ,      bee            & ! Clapp&Hornberger 'b' exponent
!         ,      poros          & ! soil porosity (fraction)
!         ,      satco          & ! hydraulic conductivity at saturation (UNITS?)
!         ,      slope          & ! cosine of mean slope
!         ,      phsat          & ! soil tension at saturation  (UNITS?)
!         ,      ccx            & ! canopy heat capacity (J m^-2 deg^-1)
!         ,      csoil          & ! soil heat capacity (J m^-2 deg^-1)
!         ,      etc            & ! 'E-star' of the canopy - vapor pressure (Pa)
!         ,      ea             & ! CAS vapor pressure (Pa)
!         ,      btc            & ! d(E(Tc))/d(Tc) - Clausius-Clapyron
!         ,      gegs           & ! dry fraction of ground/(fg*rsoil + Rd)
!         ,      hr             & ! soil surface relative humidity
!         ,      fg             & ! flag indicating direction of vapor pressure
!         ! deficit between CAS and ground: 0 => ea>e(Tg)
!         !                                 1 => ea<e(Tg)
!         ,      gegi           & ! wet fraction of ground/Rd
!         ,      rd             & ! ground-CAS resistance
!         ,      rb             & ! leaf-CAS resistance
!         ,      hcdtc          & ! dHc/dTc
!         ,      hcdta          & ! dHc/dTa
!         ,      hgdta          & ! dHg/dTa
!         ,      slamda         & !
!         ,      etg            & ! 'E-star' of the ground surface  (Pa)
!         ,      btg            & ! d(E(tg))/d(Tg) - Clausius-Clapyron
!         ,      geci           & ! wetted fraction of canopy/2Rb
!         ,      gect           & ! dry fraction of canopy/(Rst + 2Rb)
!         ,      ros            &! air density (kg m^-3)
!         ,      psy             !


!    REAL, INTENT(in)  ::&
!         snofac         & !  ___(lat ht of vap)___     (unitless)
!         ! (lat ht vap + lat ht ice)
!         ,      hlat           & ! latent heat of vaporization of water (J kg^-1)
!         ,      pi             & ! 3.1415...
!         ,      dtt            & ! timestep (seconds)
!         ,      tf             & ! freezing temperature (273.16 K)
!         ,      snomel         & ! latent heat of fusion for ice (J m^-3)
!         ,      asnow          & ! conversion factor for kg water to depth
!         !   of snow (16.7)
!         ,      cw             & ! water heat capacity (J m^-3 deg^-1)
!         ,      cp            ! specific heat of air at const pres (J kg-1 deg-1)


!    REAL, INTENT(in),DIMENSION(nsib,3)  ::&
!         zdepth          ! soil layer depth * porosity (meters)

!    REAL, INTENT(in),DIMENSION(nsib,nsoil)         ::&
!         td            ! deep soil temperature (K)


!    !
!    !----------------------------------------------------------------------
!    !
!    !              INTENT = OUT VARIABLES
!    !
!    !----------------------------------------------------------------------

!    REAL,INTENT(out),DIMENSION(len)                ::&
!         ect             &! transpiration flux (J m^-2 for the timestep)
!         ,      eci            & ! interception flux (veg - CAS) (J m^-2)
!         ,      egi            & ! ground interception flux (J m^-2)
!         ,      egs             &! ground evaporative flux (J m^-2)
!         ,      ecmass         & ! canopy evaporation (mm)
!         ,      egmass          &! ground evaporation (mm)
!         ,      shf             &! soil heat flux (W m^-2)
!         ,      chf             &! canopy heat flux (W m^-2)
!         ,      q3l             &! 'Liston' drainage from the bottom of
!         !    soil layer 3 (mm)
!         ,      q3o             &! gravitational drainage out of soil
!         !    layer 3 (mm)
!         ,      zmelt           ! depth of melted water (m)

!    REAL,INTENT(out),DIMENSION(len,nsoil)          ::&
!         dtd             ! deep soil temperature increment (K)

!    REAL,INTENT(out),DIMENSION(len,3)              ::&
!         QQQ             ! soil layer drainage (mm m^-2 timestep)


!    !
!    !----------------------------------------------------------------------
!    !
!    !              INTENT = IN/OUT VARIABLES
!    !
!    !----------------------------------------------------------------------

!    REAL,INTENT(inout),DIMENSION(len,3)            ::&
!         www             ! soil moisture (% of saturation)

!    REAL,INTENT(inout),DIMENSION(nsib,2)           ::&
!         snoww           &! snow-interception (1-veg, 2-ground) (meters)
!         ,      capac           ! liquid interception

!    REAL,INTENT(inout),DIMENSION(len,2)            ::&
!         dtg             ! delta ground surface temperature (K)

!    REAL,INTENT(inout),DIMENSION(len)              ::&
!         roff           & ! runoff (mm)
!         ,      hc              &! canopy sensible heat flux (J m^-2)
!         ,      hg              &! ground sensible heat flux (J m^-2)
!         ,      areas           &! fractional snow coverage (0 to 1)
!         ,      dtc             ! delta vegetation temperature (K)


!    !
!    !----------------------------------------------------------------------
!    !
!    !              LOCAL VARIABLES
!    !
!    !----------------------------------------------------------------------

!    REAL,DIMENSION(len)                            ::&
!         dtgs          & ! snow/dry soil averaged temp increment (K)
!         ,      zmelt_total    &! total melting
!         ,      cctt           &
!         ,      cct            &
!         ,      ts             &
!         ,      dts            &
!         ,      flux           &
!         ,      fluxef         &
!         ,      tsnow          &
!         ,      dpdw           &
!         ,      tgs            &
!         ,      dts2           &
!         ,      cpdpsy         &
!         ,      heaten          ! energy to heat snow to ground temp (J m^-2)



!    REAL,DIMENSION(len,3)                          ::&
!         TEMW            &
!         ,    TEMWP               &
!         ,    TEMWPP

!    REAL,DIMENSION(len,2)                          ::&
!         AAA               &
!         ,    BBB              &
!         ,    CCC

!    REAL                                           ::&
!         hlati            &
!         ,    rsnow            &
!         ,    facks            &
!         ,    realc            &
!         ,    realg            &
!         ,    dpdwdz            &
!         ,    qmax            &
!         ,    qmin             &
!         ,    rdenom            &
!         ,    denom            &
!         ,    props            &
!         ,    avkmax            &
!         ,    avkmin            &
!         ,    div            &
!         ,    rsame            &
!         ,    pmin            &
!         ,    wmin            &
!         ,    pmax            &
!         ,    wmax            &
!         ,    egsdif            &
!         ,    ectdif            &
!         ,    extrak            &
!         ,    facl            &
!         ,    dtsg2            &
!         ,    dtsg3            &
!         ,    cool            &
!         ,    zmelt2            &
!         ,    dtsg            &
!         ,    heat            &
!         ,    exmelt            &
!         ,    exheat            &
!         ,    safe            &
!         ,    avheat            &
!         ,    avmelt            &
!         ,    avk            &
!         ,    pows            &
!         ,    avex            &
!         ,    freeze            &
!         ,    ecpot

!    REAL :: egpot        &
!         ,    hrr        &
!         ,    ecidif     &! actual amount of canopy moisture put into CAS(J m^-2)
!         ,    egidif     &! actual amount of ground interception moisture
!         ! put into CAS (J m^-2))
!         ,    egit       &! temporary ground heat flux holder (J m^-2)
!         ,    t1,t2      &! snow depth measures (intermediate)
!         ,    aven       &! energy difference between actual snow and areas=1
!         ,    darea      &! adjustment to areas
!         ,    arean      &! adjustment to areas
!         ,    ectmax     &! upper bound for transpiratoin (J m^-2)
!         ,    egsmax     &! upper bound for soil evaporation (J m^-2)
!         ,    dti        &! 1/timestep
!         ,    timcon     &! pi/seconds per day
!         ,    cogs1      &! non snowcovered fraction * soil humidity
!         ,    cogs2      ! non snowcovered fraction


!                                                              !
!    !----------------------------------------------------------------------
!    !
!    !----------------------------------------------------------------------
!    !

!    dti = 1./DTT
!    timcon = 3.1415926 / 86400.0

!    DO i=1,len
!       tsnow(i) = MIN(tf - 0.01, tg(i))
!       cpdpsy(i) = cp/psy(i)
!       tgs(i) = (1.0 - areas(i))*tg(i) + areas(i)*tsnow(i)
!       dtgs(i) = (1.0-areas(i))*dtg(i,1) + areas(i)*dtg(i,2)
!       rsnow = snoww(i,2) / (snoww(i,2) + capac(i,2) + 1.0E-10)
!       !pl this is the potential gradient in Pa
!       !pl this WAS realized in sibslv

!       ECPOT =     (etc(i) + btc(i)*DTC(i)) - (ea(I) + DEA(i))

!       !pl and this is the  INTERCEPTION flux in J/m2
!       ECI(i) = ECPOT * geci(i) * ros(i) * CPDPSY(i) * DTT

!       !pl and this is the TRANSPIRATION flux in J/m2
!       ECT(i) = ECPOT * gect(i) * ros(i) * CPDPSY(i) * DTT


!       !pl this is the potential gradient in Pa
!       !pl this WAS realized in sibslv25

!       EGPOT =   (etg(i) + btg(i)*DTGS(i)) - ( ea(i) + DEA(i))

!       !pl and this is the  INTERCEPTION flux in J/m2
!       EGI(i) = EGPOT * (gegi(i) * (1.-areas(i)) + areas(i)/rd(i))     &
!            * ros(i) * cpdpsy(i) * DTT

!       HRR = HR(i)
!       IF ( FG(i) .LT. .5 ) HRR = 1.
!       COGS1    =  gegs(i) * HRR * (1.-AREAS(i))
!       COGS2    =  gegs(i)       * (1.-AREAS(i))

!       !pl and this is the EVAPORATIVE flux in J/m2
!       EGS(i) =  (etg(i) + btg(i) * dtgs(i)) * COGS1     &
!            -(ea(i) +           dea(i)) * COGS2
!       EGS(i) = EGS(i) * ros(i) * CPDPSY(i) * DTT

!       !itb...make sure you don't evap more than you have...
!       EGSMAX = WWW(i,1) * ZDEPTH(i,1) * hlat * 1.e3 * 0.5
!       EGS(i) = MIN ( EGS(i), EGSMAX )
!       !itb...make sure you don't transpire more water than is in the soil
!       ECTMAX = WWW(i,2) * ZDEPTH(i,2) * hlat * 1.e3 * 0.5
!       ECT(i) = MIN ( ECT(i), ECTMAX )


!       !itb...these fluxes were all realized in sibslv. If positive, they
!       !itb...imply transfer of water vapor INTO the CAS. If negative,
!       !itb...they imply transfer OUT OF the CAS. We need to adjust
!       !itb...the various reserviors as well as the CAS vapor capacity,
!       !itb...making sure that none go to negative values.

!       !itb...the actual movement of vapor is taken care of in the
!       !itb...equations in sibslv. All we do now is adjust the surface
!       !itb...and vegetation storage reservoirs to reflect what we've
!       !itb...already added or taken out.

!       !pl this is the limitation to the ECI flux in J/m2

!       ECIdif=MAX(0.0E0,(ECI(i)-(SNOWw(i,1)+CAPAC(i,1))     &
!            *1.E3*hlat))

!       ECI(i)   =MIN(ECI(i),                               &
!            ( (SNOWw(i,1)+CAPAC(i,1))*1.E3*hlat))


!       !pl this is the EGI flux in J/m2

!       EGIdif=                                                    &
!            MAX(0.0E0,EGI(i)-(SNOWw(i,2)+CAPAC(i,2))*1.E3*hlat)   &
!            *(1.-RSNOW)


!       EGIT  =                                                    &
!            MIN(EGI(i), (SNOWw(i,2)+CAPAC(i,2))*1.E3*hlat )        &
!            *(1.-RSNOW)


!       !itb...print this stuff out, for grins
!       !        print*,'updat2: eci,ect,egi,egs,ecidif,egidif'
!       !        print*,eci(i),ect(i),egi(i),egs(i),ecidif,egidif


!       !
!       !----------------------------------------------------------------------
!       !     CALCULATION OF INTERCEPTION LOSS FROM GROUND-SNOW. IF SNOW PATCH
!       !     SHRINKS, ENERGY IS TAKEN FROM EGI TO WARM EXPOSED SOIL TO TGS.
!       !----------------------------------------------------------------------
!       !
!       T1 = SNOWw(i,2) - 1./ASNOW
!       T2 = MAX( 0.E0, T1 )
!       AVEN = EGI(i) - T2*hlat*1.E3/SNOFAC
!       IF ( (T1-T2)*EGI(i) .GT. 0. ) AVEN = EGI(i)
!       DAREA = AVEN/( (TSNOW(i)-TG(i))*CSOIL(i)            &
!            - 1./ASNOW*hlat*1.E3/SNOFAC)
!       AREAN = AREAS(i) + DAREA
!       EGIdif = EGIdif - MIN( 0.E0, AREAN )              &
!            *hlat*1.E3/(asnow*SNOFAC)*RSNOW
!       DAREA = MAX( DAREA, -AREAS(i) )
!       DAREA = MIN( 1.-AREAS(i), DAREA )
!       HEATEN(i) = (TSNOW(i)-TG(i))*CSOIL(i)*DAREA*RSNOW
!       EGIT = EGIT + ( EGI(i) - HEATEN(i) - EGIdif )*RSNOW
!       EGI(i) = EGIT


!       !---------------------------------------------------------------------
!       !     CALCULATION OF SENSIBLE HEAT FLUXES FOR THE END OF THE TIMESTEP.
!       !        SEE FIGURE (2) OF SE-86.  NOTE THAT INTERCEPTION LOSS EXCESS
!       !        ENERGIES (ECIDIF, EGIDIF) ARE ADDED.
!       !
!       !      HC          (HC)    : EQUATION (63) , SE-86
!       !      HG          (HGS)   : EQUATION (65) , SE-86
!       !----------------------------------------------------------------------
!       !
!       !        HC(i) = HC(i) + (HCDTC(i)*DTC(i)
!       !     &                +  HCDTA(i)*dta(i))*DTT + ECIdif
!       !        HG(i) = HG(i) + (HGDTC(i)*DTC(i)
!       !     &                +  HGDTA(i)*dta(i))*DTT + EGIdif

!       !itb...i've left the leaf one-sided, for now...
!       HC(i) = ( (tc(i)+dtc(i)) - (ta(i)+dta(i)) ) /rb(i)         &
!            * ros(i) * cp * DTT + ECIdif

!       !itb...ground sensible heat flux includes soil and snow by using
!       !itb...dtgs
!       HG(i) = ( (tg(i)+dtgs(i)) - (ta(i)+dta(i)) ) /rd(i)          &
!            * ros(i) * cp * DTT + EGIdif


!       CHF(i) = CCX(i) * dti * DTC(i)

!    ENDDO
!    !----------------------------------------------------------------------
!    !     CALCULATION OF STORAGE HEAT FLUXES
!    !
!    !----------------------------------------------------------------------
!    !

!    IF(forcerestore) THEN
!       DO i = 1,len
!          SHF(i) = dti * ( (1.-areas(i))*dtg(i,1) +                   &
!               areas(i)*dtg(i,2) ) * cg(i)             &
!               + TIMCON*csoil(i)*2. *( TGS(i)+dtgs(i) - TD(i,nsoil) )
!       ENDDO
!    ELSE
!       DO i = 1,len
!          !  new soil thermodynamic model
!          SHF(i) = dti * ( (1.-areas(i))*dtg(i,1) +                   &
!               areas(i)*dtg(i,2) ) * cg(i)             &
!               + slamda(i) *( TGS(i)+dtgs(i) - TD(i,nsoil) )
!          !           print*,'soil heat flux'
!          !           print*,'dti,areas,dtg1,2,slamda,tgs,dtgs,td1'
!          !           print*,dti,areas(i),dtg(i,1),dtg(i,2),slamda(i),tgs(i)
!          !     &       ,dtgs(i),td(i,nsoil)
!          !           print*,'SOIL HEAT FLUX = ',shf(i)
!          !           print*,'soil temps'
!          !           print*,'*********'
!       ENDDO
!    ENDIF

!    IF(forcerestore) THEN
!       ksoil = nsoil
!    ELSE
!       ksoil = 3
!    ENDIF
!    !
!    !----------------------------------------------------------------------
!    !    INTERCEPTION LOSSES APPLIED TO SURFACE WATER STORES.
!    !    EVAPORATION LOSSES ARE EXPRESSED IN J M-2 : WHEN DIVIDED BY
!    !    ( HLAT*1000.) LOSS IS IN M M-2. MASS TERMS ARE IN KG M-2 DT-1
!    !    INTERCEPTION AND DRAINAGE TREATED IN INTER2.
!    !
!    !      CAPAC/SNOWW(1) (M-C)   : EQUATION (3)  , SE-86
!    !      CAPAC/SNOWW(2) (M-G)   : EQUATION (4)  , SE-86
!    !----------------------------------------------------------------------
!    !
!    hlati = 1. / hlat
!    !PL HERE WE DO A CHECK FOR CONDENSATION AND MAKE SURE THAT IT ONLY
!    !PL HAPPENS TRHOUGH ECI AND EGI

!    DO i = 1,len
!       RSNOW = SNOWW(i,1)/(SNOWW(i,1)+CAPAC(i,1)+1.E-10)
!       FACKS = 1. + RSNOW * ( SNOFAC-1. )
!       IF ( (ECT(i)+ECI(i)) .LE. 0.) THEN
!          ECI(i) = ECT(i)+ECI(i)
!          ECT(i) = 0.
!          FACKS = 1. / FACKS
!       ENDIF
!       CAPAC(i,1) = CAPAC(i,1)-( 1.-RSNOW )*ECI(i)*FACKS*hlati*0.001
!       SNOWW(i,1) = SNOWW(i,1)-     RSNOW  *ECI(i)*FACKS*hlati*0.001
!       snoww(i,1) = MAX(snoww(i,1),0.0E0)
!       capac(i,1) = MAX(capac(i,1),0.0E0)
!       ECMASS(i) = ECI(i)*FACKS *hlati
!       zmelt_total(i) = 0.0
!    ENDDO
!    !      do i = 1,len
!    !         if(snoww(i,1).lt.0.0)
!    !     *      print *,'snoww1 after updat2 100 ',i,snoww(i,1)
!    !         if(capac(i,1).lt.0.0)
!    !     *      print *,'capac1 after updat2 100 ',i,capac(i,1)
!    !      enddo
!    !
!    DO i = 1,len
!       RSNOW = SNOWW(i,2)/(SNOWW(i,2)+CAPAC(i,2)+1.e-10)
!       FACKS = 1. + RSNOW * ( SNOFAC-1. )
!       IF ( (EGS(i)+EGI(i)) .LE. 0. ) THEN
!          EGI(i) = EGS(i)+EGI(i)
!          EGS(i)= 0.
!          FACKS = 1. / FACKS
!       ENDIF
!       CAPAC(i,2) = CAPAC(i,2)-( 1.-RSNOW )*EGI(i)*FACKS*hlati*0.001
!       SNOWW(i,2) = SNOWW(i,2)-     RSNOW  *EGI(i)*FACKS*hlati*0.001
!       snoww(i,2) = MAX(snoww(i,2),0.0E0)
!       capac(i,2) = MAX(capac(i,2),0.0E0)
!       EGMASS(i) = EGI(i)*FACKS *hlati
!    ENDDO
!    !      do i = 1,len
!    !         if(snoww(i,2).lt.0.0)
!    !     *      print *,'snoww2 after updat2 200 ',i,snoww(i,2)
!    !         if(capac(i,2).lt.0.0)
!    !     *      print *,'capac2 after updat2 200 ',i,capac(i,2)
!    !      enddo
!    !
!    !----------------------------------------------------------------------
!    !    DUMPING OF SMALL CAPAC VALUES ONTO SOIL SURFACE STORE
!    !----------------------------------------------------------------------
!    !
!    DO IVEG = 1, 2 !DO 1000 IVEG = 1, 2
!       DO i = 1,len
!          IF ( (SNOWW(i,iveg)+CAPAC(i,IVEG)) .LE. 0.00001 ) THEN
!             WWW(i,1) = WWW(i,1) + (SNOWW(i,IVEG)+CAPAC(i,IVEG)) /     &
!                  ZDEPTH(i,1)
!             CAPAC(i,IVEG) = 0.
!             SNOWW(i,IVEG) = 0.
!          ENDIF
!       ENDDO
!    ENDDO !1000 CONTINUE
!    !      do i = 1,len
!    !         if(www(i,1).le.0.0)print *,'www after updat2 1000 ',i,www(i,1)
!    !        print*,rsnow , snoww(i,2) , capac(i,2),  snoww(i,1) , capac(i,1)
!    !        stop
!    !      enddo
!    !
!                                                                 !
!    !=======================================================================
!    !------------------SNOW2 INLINED-------------------------------------
!    !----------------------------------------------------------------------
!    !    SNOWMELT / REFREEZE CALCULATION
!    !----------------------------------------------------------------------
!    !
!    !     CALCULATION OF SNOWMELT AND MODIFICATION OF TEMPERATURES
!    !
!    !     MODIFICATION DEALS WITH SNOW PATCHES:
!    !          TS < TF, TSNOW = TS
!    !          TS > TF, TSNOW = TF
!    !
!    !=======================================================================
!    !
!    DO IVEG = 1, 2 !DO 1002 IVEG = 1, 2
!       !
!       REALC = (2 - IVEG)*1.
!       REALG = (IVEG - 1)*1.
!       !
!       DO i = 1,len
!          CCTT(i) = REALC*CCX (i) +  REALG*CG(i)
!          CCT(i)  = REALC*CCX(i)  +  REALG*CSOIL(i)
!          TS(i)   = REALC*TC(i)   +  REALG*TG(i)
!          DTS(i)  = REALC*DTC(i)  +  REALG*DTG(i,1)
!          DTS2(i)  = REALC*DTC(i)  +  REALG*DTG(i,2)
!          FLUX(i) = REALC*CHF(i)  +  REALG*              &
!               ( (1.-areas(i))*DTG(i,1)+                   &
!               areas(i)*dtg(i,2)  )*cg(i) /DTT
!          !  fluxef moved up here to conserve energy
!          fluxef(i) = ( shf(i) - flux(i)) * realg
!          TSNOW(i) = MIN ( TF-0.01, TS(i) )
!          ZMELT(i) = 0.
!       ENDDO
!       !
!       DO i = 1,len  ! this scalar loop needs vector optimization
!          !itb   print*,'updat2:ts,dts,ts+dts,tf',ts(i),dts(i),ts(i)+dts(i),tf
!          IF ( SNOWW(i,IVEG) .GT. 0. ) GO TO 102
!          IF ( ( TS(i)+DTS(i)) .GT. TF ) GO TO 502
!          !----------------------------------------------------------------------
!          !
!          !     NO SNOW  PRESENT, SIMPLE THERMAL BALANCE WITH POSSIBLE FREEZING.
!          !
!          !----------------------------------------------------------------------
!          FREEZE = MIN ( 0.E0, (FLUX(i)*DTT-( TF-0.01 - TS(i))      &
!               *CCTT(i)))
!          SNOWW(i,IVEG) = MIN( CAPAC(i,IVEG), - FREEZE/SNOMEL )
!          ZMELT(i) = CAPAC(i,IVEG) - SNOWW(i,IVEG)
!          CAPAC(i,IVEG) = 0.
!          DTS(i) = DTS(i) + SNOWW(i,IVEG)*SNOMEL/CCTT(i)
!          GO TO 502
!          !
!          !----------------------------------------------------------------------
!          !
!          !     SNOW PRESENT
!          !
!          !---------------------------------------------------------------------
!  102     CONTINUE
!          !
!          !itb      IF ( TS(i) .LT. TF .AND. (TS(i)+DTS(i)) .LT. TF ) GO TO 502
!          IF ( TS(i) .LT. TF .AND. (TS(i)+DTS2(i)) .LT. TF ) GO TO 502
!          IF ( ts(i) .GT. TF ) GO TO 202
!          !----------------------------------------------------------------
!          !
!          !     SNOW PRESENT : TS < TF,  TS+DTS > TF
!          !
!          !------------------------------------------------------------
!          AVEX = FLUX(i) - ( TF-0.01 - TS(i) ) * CCTT(i)/DTT
!          AVMELT = ( AVEX/SNOMEL * (AREAS(i)*REALG + REALC ) )*DTT
!          ZMELT(i) = MIN( AVMELT, SNOWW(i,IVEG) )
!          SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT(i)
!          AVHEAT = AVEX*( 1.-AREAS(i) )*REALG +                        &
!               ( AVMELT-ZMELT(i) )*SNOMEL/DTT
!          AREAS(i) = MIN( 0.999E0, ASNOW*SNOWW(i,2) )
!          SAFE = MAX( ( 1.-AREAS(i)*REALG ), 1.E-8 )
!          DTS(i) = TF-0.01 - TS(i) + AVHEAT / ( CCTT(i)*SAFE )*DTT
!          GO TO 502
!          !----------------------------------------------------------------------
!          !
!          !     SNOW PRESENT AND TS > TF : GROUND ONLY.
!          !
!          !------------------------------------------------------------------
!  202     CONTINUE
!          !
!          EXHEAT = CCT(i)*( 1.001-MAX(0.1E0,AREAS(i))) * DTS(i)
!          EXMELT = FLUX(i)*DTT - EXHEAT
!          HEAT = EXHEAT
!          DTSG = EXHEAT / ( CCT(i)*(1.001-AREAS(i) ))
!          IF ( (TS(i)+DTSG) .GT. TF ) GO TO 302
!          HEAT = ( TF-0.01 - TS(i) ) * ( CCT(i)*(1.-AREAS(i)) )
!          DTSG = TF-0.01 - TS(i)
!          !
!  302     EXMELT = EXMELT + EXHEAT - HEAT
!          !
!          IF( EXMELT .LT. 0. ) GO TO 402
!          ZMELT(i) = EXMELT/SNOMEL
!          IF( ASNOW*(SNOWW(i,IVEG)-ZMELT(i)) .LT. 1. )                     &
!               ZMELT(i) = MAX( 0.E0, SNOWW(i,IVEG) - 1./ASNOW )
!          SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT(i)
!          !print*,'XX' ,snoww(i,2) , snoww(i,2) , capac(i,2)
!          EXMELT = EXMELT - ZMELT(i)*SNOMEL
!          ZMELT2 = EXMELT/ ( CCT(i)*( TS(i)-TF )*ASNOW + SNOMEL )
!          ZMELT2 = MIN( ZMELT2, SNOWW(i,IVEG) )
!          ZMELT(i) = ZMELT(i) + ZMELT2
!          SNOWW(i,IVEG) = SNOWW(i,IVEG) - ZMELT2
!          EXMELT = EXMELT - ZMELT2*( CCT(i)*( TS(i)-TF )*ASNOW + SNOMEL )
!          DTS(i)  = DTSG + EXMELT/CCT(i)
!          GO TO 502
!          !
!  402     COOL = MIN( 0.E0, TF-0.01 -(TS(i)+DTSG)) * CCT(i)             &
!               *(1.-AREAS(i))
!          DTSG2 = MAX ( COOL, EXMELT ) / ( CCT(i)*( 1.001-AREAS(i) ) )
!          EXMELT = EXMELT - DTSG2*CCT(i)*(1.-AREAS(i))
!          DTSG3 =EXMELT/CCTT(i)
!          DTS(i) = DTSG + DTSG2 + DTSG3
!          !
!  502     CONTINUE
!       ENDDO
!       !
!       DO i = 1,len
!          !itb...patch
!          IF(ZMELT(i) < 0.0 ) ZMELT(i) = 0.0
!          !itb...patch
!          WWW(i,1) = WWW(i,1) + ZMELT(i) /                     &
!               ZDEPTH(i,1)

!          IF(www(i,1).LT.0.0)zmelt_total(i) = SQRT(www(i,1))
!          !
!          DTC(i) = DTC(i)*REALG + DTS(i)*REALC
!          DTG(i,1) = DTG(i,1)*REALC + DTS(i)*REALG
!          zmelt_total(i) = zmelt_total(i) + zmelt(i)
!       ENDDO
!       !
!    ENDDO !1002 CONTINUE
!    !

!    !itb...put zmelt_total into zmelt
!    zmelt(:) = zmelt_total(:)

!    IF(forcerestore) THEN
!       DO i = 1,len
!          !960320 fluxef calculation moved up to conserve energy
!          !950511 changed cg to csoil
!          DTD(i,nsoil) = FLUXEF(i) /                              &
!               ( csoil(i) * 2.*SQRT( PI*365. ) ) * DTT
!       ENDDO
!    ENDIF
!    DO i = 1,len
!       !
!       !------------------END SNOW2  -------------------------------------
!       !----------------------------------------------------------------------
!       !    EVAPOTRANSPIRATION LOSSES APPLIED TO SOIL MOISTURE STORE.
!       !    EXTRACTION OF TRANSPIRATION LOSS FROM ROOT ZONE, SOIL EVAPORATION..
!       !
!       !      ECT         (E-DC)  : EQUATION (5,6), SE-86
!       !      EGS         (E-S)   : EQUATION (5)  , SE-86
!       !----------------------------------------------------------------------
!       !
!       !PL STEP THREE part II
!       !pl we have done the potential ECT and EGS inside ENDTEM25
!       !pl now we limit these fluxes according to 1/2 of what is
!       !pl in the soil reservoirs, WWW(i,1) for EGS and WWW(i,2) for ECT
!       !pl we 'donate' the excess to HC and HG, if any.

!       FACL   = hlati*0.001/ZDEPTH(i,2)
!       EXTRAK = ECT(i)*FACL
!       EXTRAK = MIN( EXTRAK, WWW(i,2)*0.5 )
!       ECTDIF = ECT(i) - EXTRAK/FACL
!       ECT(i)    = EXTRAK/FACL
!       HC(i)     = HC(i) + ECTDIF
!       ECMASS(i) = ECMASS(i) + ECT(i)*hlati
!       WWW(i,2) = WWW(i,2) - ECT(i)*FACL
!       !
!       FACL   = 0.001*hlati/ZDEPTH(i,1)
!       EXTRAK = EGS(i)*FACL
!       EXTRAK = MIN( EXTRAK, WWW(i,1) *0.5 )
!       EGSDIF = EGS(i) - EXTRAK/FACL
!       EGS(i)    = EXTRAK/FACL
!       HG(i)     = HG(i) + EGSDIF
!       EGMASS(i) = EGMASS(i) + EGS(i)*hlati
!       WWW(i,1) = WWW(i,1) - EGS(i)*FACL


!    ENDDO
!    !
!    !========================================================================
!    !------------------RUN2 INLINED-------------------------------------
!    !========================================================================
!    !----------------------------------------------------------------------
!    !    CALCULATION OF INTERFLOW, INFILTRATION EXCESS AND LOSS TO
!    !    GROUNDWATER .  ALL LOSSES ARE ASSIGNED TO VARIABLE 'ROFF' .
!    !----------------------------------------------------------------------
!    !
!    !
!    DO I = 1, 3  !DO 1001 I = 1, 3
!       !
!       DO l = 1,len
!          TEMW(l,I)   = MAX( 0.03E0, WWW(l,I) )
!          TEMWP(l,I)  = TEMW(l,I) ** ( -BEE(l) )
!          TEMWPP(l,I) = MIN( 1.E0, TEMW(l,I))**( 2.*BEE(l)+ 3. )
!       ENDDO
!    ENDDO !1001 CONTINUE
!    !
!    !-----------------------------------------------------------------------
!    !
!    !    CALCULATION OF GRAVITATIONALLY DRIVEN DRAINAGE FROM W(3) : TAKEN
!    !    AS AN INTEGRAL OF TIME VARYING CONDUCTIVITY.
!    !
!    !     qqq(3) (Q3) : EQUATION (62) , SE-86
!    !
!    !    QQQ(3) IS AUGMENTED BY A LINEAR LOSS TERM RECOMMENDED BY LISTON (1992)
!    !-----------------------------------------------------------------------
!    !
!    DO i = 1,len
!       POWS = 2.*BEE(i)+2.
!       qqq(i,3) = TEMW(i,3)**(-POWS) + SATCO(i)/            &
!            (ZDEPTH(i,3) )*                         &
!            SLOPE(i)*POWS*DTT
!       qqq(i,3) = qqq(i,3) ** ( 1. / POWS )
!       qqq(i,3) = - ( 1. / qqq(i,3) - WWW(i,3) ) *           &
!            ZDEPTH(i,3) / DTT
!       qqq(i,3) = MAX( 0.E0, qqq(i,3) )
!       q3o(i) = qqq(i,3) * dtt
!       qqq(i,3) = MIN( qqq(i,3), WWW(i,3)*                    &
!            ZDEPTH(i,3)/DTT )
!       !
!       Q3l(i) = 0.002*ZDEPTH(i,3)*0.5 / 86400.*          &
!            MAX(0.E0,(www(i,3)-0.01)/0.99 )
!       qqq(i,3) = qqq(i,3) + q3l(i)
!       q3l(i) = q3l(i) * dtt
!       !
!       !----------------------------------------------------------------------
!       !
!       !    CALCULATION OF INTER-LAYER EXCHANGES OF WATER DUE TO GRAVITATION
!       !    AND HYDRAULIC GRADIENT. THE VALUES OF W(X) + DW(X) ARE USED TO
!       !    CALCULATE THE POTENTIAL GRADIENTS BETWEEN LAYERS.
!       !    MODIFIED CALCULATION OF MEAN CONDUCTIVITIES FOLLOWS MILLY AND
!       !    EAGLESON (1982 ), REDUCES RECHARGE FLUX TO TOP LAYER.
!       !
!       !      DPDW           : ESTIMATED DERIVATIVE OF SOIL MOISTURE POTENTIAL
!       !                       WITH RESPECT TO SOIL WETNESS. ASSUMPTION OF
!       !                       GRAVITATIONAL DRAINAGE USED TO ESTIMATE LIKELY
!       !                       MINIMUM WETNESS OVER THE TIME STEP.
!       !
!       !      QQQ  (Q     )  : EQUATION (61) , SE-86
!       !             I,I+1
!       !            -
!       !      AVK  (K     )  : EQUATION (4.14) , ME-82
!       !             I,I+1
!       !
!       !----------------------------------------------------------------------
!       !
!       WMAX = MAX( WWW(i,1), WWW(i,2), WWW(i,3), 0.05E0 )
!       WMAX = MIN( WMAX, 1.E0 )
!       PMAX = WMAX**(-BEE(i))
!       WMIN = (PMAX-2.*poros(i)/( PHSAT(i)*                 &
!            (ZDEPTH(i,1)+2.*ZDEPTH(i,2)+ZDEPTH(i,3))))          &
!            **(-1./BEE(i))
!       WMIN = MIN( WWW(i,1), WWW(i,2), WWW(i,3), WMIN )
!       WMIN = MAX( WMIN, 0.02E0 )
!       PMIN = WMIN**(-BEE(i))
!       DPDW(i) = PHSAT(i)*( PMAX-PMIN )/( WMAX-WMIN )
!    ENDDO
!    !
!    DO I = 1, 2 !DO 2001 I = 1, 2
!       !
!       DO l = 1,len
!          RSAME = 0.
!          AVK  = TEMWP(l,I)*TEMWPP(l,I) - TEMWP(l,I+1)*TEMWPP(l,I+1)
!          DIV  = TEMWP(l,I+1) - TEMWP(l,I)
!          IF ( ABS(DIV) .LT. 1.E-6 ) RSAME = 1.
!          AVK = SATCO(l)*AVK /                                       &
!               ( ( 1. + 3./BEE(l) ) * DIV + RSAME )
!          AVKMIN = SATCO(l) * MIN( TEMWPP(l,I), TEMWPP(l,I+1) )
!          AVKMAX = SATCO(l) * MAX( TEMWPP(l,I), TEMWPP(l,I+1) )*1.01
!          AVK = MAX( AVK, AVKMIN )
!          AVK = MIN( AVK, AVKMAX )
!          !
!          !c Collatz-Bounoua change to effective hydraulic conductivity making
!          !  it 10x harder for water to move up than down if the upper soil layer
!          !  is wetter than lower soil layer.
!          !----------------------------------------------------------------------
!          IF (www(l,i) .LT. www(l,i+1)) avk = 0.1 * avk     ! lahouari
!          !
!          !----------------------------------------------------------------------
!          !     CONDUCTIVITIES AND BASE FLOW REDUCED WHEN TEMPERATURE DROPS BELOW
!          !        FREEZING.
!          !----------------------------------------------------------------------

!          TSNOW(l) = MIN ( TF-0.01, TG(l) )
!          TGS(l) = TSNOW(l)*AREAS(l) + TG(l)*(1.-AREAS(l))
!          TS(l)    = TGS(l)*(2-I) + TD(l,ksoil)*(I-1)
!          PROPS = ( TS(l)-(TF-10.) ) / 10.
!          PROPS = MAX( 0.05E0, MIN( 1.0E0, PROPS ) )
!          AVK  = AVK * PROPS
!          qqq(l,3)  = qqq(l,3) * PROPS

!          !----------------------------------------------------------------------
!          !     BACKWARD IMPLICIT CALCULATION OF FLOWS BETWEEN SOIL LAYERS.
!          !----------------------------------------------------------------------

!          DPDWDZ = DPDW(l)*2.*poros(l)/                         &
!               ( ZDEPTH(l,I) + ZDEPTH(l,I+1) )
!          AAA(l,I) = 1. + AVK*DPDWDZ*                            &
!               ( 1./ZDEPTH(l,I)+1./ZDEPTH(l,I+1) )         &
!               *DTT
!          BBB(l,I) =-AVK *   DPDWDZ * 1./ZDEPTH(l,2)*DTT
!          CCC(l,I) = AVK * ( DPDWDZ * ( WWW(l,I)-WWW(l,I+1) ) + 1. +        &
!               (I-1)*DPDWDZ*qqq(l,3)*1./ZDEPTH(l,3)*                  &
!               DTT )
!       ENDDO
!    ENDDO !2001 CONTINUE
!    !
!    DO i = 1,len
!       DENOM  = ( AAA(i,1)*AAA(i,2) - BBB(i,1)*BBB(i,2) )
!       RDENOM = 0.
!       IF ( ABS(DENOM) .LT. 1.E-6 ) RDENOM = 1.
!       RDENOM = ( 1.-RDENOM)/( DENOM + RDENOM )
!       QQQ(i,1)  = ( AAA(i,2)*CCC(i,1) - BBB(i,1)*CCC(i,2) ) * RDENOM
!       QQQ(i,2)  = ( AAA(i,1)*CCC(i,2) - BBB(i,2)*CCC(i,1) ) * RDENOM

!       !-----------------------------------------------------------------------
!       !     UPDATE WETNESS OF EACH SOIL MOISTURE LAYER DUE TO LAYER INTERFLOW
!       !        AND BASE FLOW.
!       !-----------------------------------------------------------------------

!       WWW(i,3) = WWW(i,3) - qqq(i,3)*DTT/ZDEPTH(i,3)
!       ROFF(i) = ROFF(i) + qqq(i,3) * DTT
!    ENDDO

!    DO I = 1, 2 !DO 3001 I = 1, 2
!       !
!       DO l = 1,len
!          QMAX   =  WWW(l,I)   *                       &
!               (ZDEPTH(l,I)  /DTT) * 0.5
!          QMIN   = -WWW(l,I+1) *                       &
!               (ZDEPTH(l,I+1)/DTT) * 0.5
!          QQQ(l,I) = MIN( QQQ(l,I),QMAX)
!          QQQ(l,I) = MAX( QQQ(l,I),QMIN)
!          WWW(l,I)   =   WWW(l,I)   -                  &
!               QQQ(l,I)/ZDEPTH(l,I) *DTT
!          WWW(l,I+1) =   WWW(l,I+1) +                  &
!               QQQ(l,I)/ZDEPTH(l,I+1)*DTT
!       ENDDO
!    ENDDO !3001 CONTINUE


!    !      do i = 1,len
!    !         if(www(i,1).le.0.0)
!    !      print *,'www after updat2 ',i,www(i,1)
!    !     print*,'XX' ,snoww(i,2) , snoww(i,2) , capac(i,2)
!    !      enddo
!    RETURN
!  END SUBROUTINE UPDAT2

!  !----------------------------------------------------------------
!  SUBROUTINE soiltherm( td, dtd, tgs, dtg, slamda, shcap, ztdep,   &
!       dt, nsib, len, nsoil )

!    IMPLICIT NONE

!    !     this subroutine calculates the temperature increments dtd for
!    !         the soil, based on the soil thermodynamic model of Bonan
!    !     layer 1 is the deepest layer, layer nsoil is immediately below the
!    !         surface
!    !     the time step is crank-nicholson. a tridiagonal matrix system is solved

!    !     argument list variables
!    INTEGER len, nsib, nsoil
!    REAL                                                   &
!         dtd(len,nsoil), tgs(len), slamda(len,nsoil),    &
!         shcap(len,nsoil), dtg(len), ztdep(len,nsoil),   &
!         td(nsib,nsoil), dt

!    !     local variables    a(n)*dtd(n-1)+b(n)*dtd(n)+c(n)*dtd(n+1) = r(n)
!    REAL                                            &
!         a(len,nsoil)    &! lower sub-diagonal
!         ,     b(len,nsoil)    &! diagonal
!         ,     c(len,nsoil)    &! upper sib-diagonal
!         ,     r(len,nsoil)    &! right hand side
!         ,     lamtem, rtem , dti, fac
!    INTEGER i, n

!    dti = 1. / dt   ! inverse time step

!    !     construct matrix
!    DO n = 1,nsoil
!       DO i = 1,len
!          b(i,n) = shcap(i,n) * dti
!          r(i,n) = 0.0
!       ENDDO
!    ENDDO
!    DO n = 1,nsoil-1
!       !!   !DIR$ IVDEP
!       DO i = 1,len
!          lamtem = -0.5 * slamda(i,n)
!          rtem = slamda(i,n) * (td(i,n+1) - td(i,n))
!          a(i,n+1) = lamtem
!          c(i,n) = lamtem
!          b(i,n) = b(i,n) - c(i,n)
!          b(i,n+1) = b(i,n+1) - a(i,n+1)
!          r(i,n) = r(i,n) + rtem
!          r(i,n+1) = r(i,n+1) - rtem
!       ENDDO
!    ENDDO
!    DO i = 1,len
!       r(i,nsoil) = r(i,nsoil) + slamda(i,nsoil) * (tgs(i)+dtg(i)   &
!            - td(i,nsoil))
!    ENDDO

!    !     eliminate lower diagonal
!    DO n = 1,nsoil - 1
!       !! !DIR$ IVDEP
!       DO i = 1,len
!          fac = a(i,n+1) / b(i,n)
!          b(i,n+1) = b(i,n+1) - c(i,n) * fac
!          r(i,n+1) = r(i,n+1) - r(i,n) * fac
!       ENDDO
!    ENDDO
!    !     back-substitution
!    DO i = 1,len
!       dtd(i,nsoil) = r(i,nsoil) / b(i,nsoil)
!    ENDDO
!    DO n = nsoil-1,1,-1
!       !! !DIR$ IVDEP
!       DO i = 1,len
!          dtd(i,n) = (r(i,n) - c(i,n) * dtd(i,n+1)) / b(i,n)
!       ENDDO
!    ENDDO

!    RETURN
!  END SUBROUTINE soiltherm

!  !---------------------------------------------------------------
!  !
!  !==================SUBROUTINE ADDINC25====================================
!  !
!  SUBROUTINE addinc(grav2,cp,dtt,hc,hg,                  &
!       ps, bps,ecmass,psy,rho, hltm, cas_cap_heat,      &
!       cas_cap_vap,                                     &
!       egmass,fss,fws,hflux,etmass,psb,                 &
!       td,thm,ts,sh,tc,tg,ta,ea,ra,em,sha,              &
!       dtd,dtc,dtg,dta,dea,drst,rst,bintc, len, nsib, nsoil,  &
!       fws_out,ea_out,em_out,ra_out,cp_out,rho_out,psy_out )


!    IMPLICIT NONE
!    !
!    !=======================================================================
!    !
!    !        Add prognostic variable increments to prognostic variables
!    !           and diagnose evapotranspiration and sensible heat flux.
!    !
!    !        Modified for multitasking - introduced gather/scatter indices
!    !           - DD 951206
!    !
!    !=======================================================================
!    !

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       TC             CANOPY TEMPERATURE (K)
!    !       TG             GROUND SURFACE TEMPERATURE (K)
!    !       TD             DEEP SOIL TEMPERATURE (K)
!    !       THM            MIXED LAYER POTENTIAL TEMPERATURE (K)
!    !       QM (gsh)       MIXED LAYER MIXING RATIO (KG KG-1)
!    !       ETMASS (FWS)   EVAPOTRANSPIRATION (MM)
!    !pl now FWS mm/s, not mm !
!    !       HFLUX (FSS)    SENSIBLE HEAT FLUX (W M-2)
!    !       rst (FSS)     STOMATAL RESISTANCE (S M-1)
!    !
!    !++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!    !
!    !       DTH            MIXED LAYER POTENTIAL TEMPERATURE INCREMENT (K)
!    !       DQM            MIXED LAYER MIXING RATIO INCREMENT (KG KG-1)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len, nsib, nsoil
!    REAL                                                         &
!         grav2,cp,dtt,hc(len),hg(len),bps(len),ecmass(len),     &
!         egmass(len),fss(len),fws(len),psb(len),ps(len),        &
!         ta(len),ea(len),ra(len),em(len),hflux(len),etmass(len),&
!         td(nsib,nsoil),thm(len),sh(len),tc(len),tg(len),       &
!         dtd(len,nsoil),dtc(len),dtg(len),DRST(len),            &
!         rst(len),bintc(len),psy(len), rho(len), z1(len),z2(len), &
!         sha(len),ts(len),dta(len),dea(len)
!    !     local variables
!    INTEGER i, n
!    REAL    &
!         dqm, dth, ten, hltm, cas_cap_heat, cas_cap_vap,  &
!         fws_out,ea_out,em_out,ra_out,cp_out,rho_out,psy_out

!    ten = 10.0

!    DO I=1,len
!       !   print*,dtc(i),dtg(i),dta(i),len,tc(i)

!       tc(i)   = tc(i) + dtc(i)
!       tg(i)   = tg(i) + dtg(i)
!       ta(i)   = ta(i) + dta(i)
!       ea(i)   = ea(i) + dea(i)


!       IF(tg(i).LT.0.0 .OR.   &
!            ta(i).LT.0.0 .OR.   &
!            ea(i).LT.0.0 .OR.   &
!            tc(i).LT.0.0)THEN
!          !print*,'BAD Ta OR ea VALUE:'
!          !print*,'SiB point:',i
!          !print*,'ea:',ea(i)-dea(i),dea(i),ea(i)
!          !print*,'Ta:',ta(i)-dta(i),dta(i),ta(i)
!          !print*,'Tg:',tg(i)-dtg(i),dtg(i),tg(i)
!          !print*,'Tc:',tc(i)-dtc(i),dtc(i),tc(i)
!          !print*,' '
!       ENDIF


!       !          !print*,'addinc: new tc,tg,ta,ea=',tc(1),tg(1),ta(1),ea(1)
!       !    stop

!       !pl now do the flux from the CAS to the ref level
!       !pl vidale et al. (1999) equation ??
!       !pl here we are using W/m2

!       FSS(i) = rho(i)*cp * (ta(i)-ts(i)) / ra(i)
!       hflux(i) = fss(i)

!       SH(i)  = 0.622 / ( ps(i)/em(i) -1.)
!       SHa(i) = 0.622 / ( ps(i)/ea(i) -1.)


!       !pl this is the latent heat flux from the CAS
!       !pl instead of using W/m2 we stick to kg/m^2/s
!       !pl the conversion is then done at the output, for
!       !pl instance in RAMS module scontrol
!       !pl vidale et al. (1999) equation ??

!       !pl so, here we have want W/m2
!       !pl in the next equation we need to multiply by cpdpsy



!       fws(i) = (ea(i) - em(i)) / ra(i) * cp * rho(i) / psy(i)

!  !itb_csu...diagnostics...
!       fws_out = fws(i)
!       ea_out  = ea(i)
!       em_out  = em(i)
!       ra_out  = ra(i)
!       cp_out  = cp
!       rho_out = rho(i)
!       psy_out = psy(i)

!       !pl but now let us go back to mm/s (or kg/(m2 s)) for the water flux,
!       !pl in order to keep to the (confusing) system we had before

!       fws(i) = fws(i) / hltm
!       etmass(i) = fws(i)


!       rst(i) = rst(i) + drst(i)
!       ! bintc(i)- smallest canopy stomatal conductance needs to be passed in here.
!       ! ---- c.zhang, 2/3/93
!       rst(i)=MIN( 1./bintc(i), rst(i) )
!       rst(i)=MAX( ten, rst(i) )
!       !          rst(i)=MAX( 1., rst(i) )
!    ENDDO
!    DO n = 1,nsoil
!       DO I=1,len
!          TD(i,n)    = TD(i,n) + dtd(i,n)
!       ENDDO
!    ENDDO

!    RETURN
!  END SUBROUTINE addinc

!  !----------------------------------------------------------------------
!  !
!  !=====================SUBROUTINE INTER2=================================
!  !
!  SUBROUTINE INTER2(ppc,ppl,snoww,capac,www           &
!       ,   pie,satcap,cw,tc,tg,clai,zlt,chil,roff          &
!       ,   snomel,zdepTH,tm,tf,asnow,csoil,satco,dtt,vcover,roffo   &
!       ,   zmelt, len, nsib, exo)

!    IMPLICIT NONE

!    INTEGER len, nsib,l
!    REAL snoww(nsib,2),capac(nsib,2),satcap(len,2),WWW(len,3),     &
!         ZDEPTH(nsib,3), ppc(len), ppl(len), pie, cw, tc(len),     &
!         tg(len), clai, zlt(len), chil(len), roff(len), snomel,    &
!         tm(len), tf, asnow, csoil(len), satco(len), dtt,          &
!         vcover(len), roffo(len), zmelt(len), exo(len), excess

!    !=======================================================================
!    !
!    !     CALCULATION OF  INTERCEPTION AND DRAINAGE OF RAINFALL AND SNOW
!    !     INCORPORATING EFFECTS OF PATCHY SNOW COVER AND TEMPERATURE
!    !     ADJUSTMENTS.
!    !
!    !----------------------------------------------------------------------
!    !
!    !     (1) NON-UNIFORM PRECIPITATION
!    !         CONVECTIVE PPN. IS DESCRIBED BY AREA-INTENSITY
!    !         RELATIONSHIP :-
!    !
!    !                   F(X) = A*EXP(-B*X)+C
!    !
!    !         THROUGHFALL, INTERCEPTION AND INFILTRATION
!    !         EXCESS ARE FUNCTIONAL ON THIS RELATIONSHIP
!    !         AND PROPORTION OF LARGE-SCALE PPN.
!    !         REFERENCE: SA-89B, APPENDIX.
!    !
!    !     (2) REORGANISATION OF SNOWMELT AND RAIN-FREEZE PROCEDURES.
!    !               SUBROUTINE ADJUST
!    !
!    !     (3) ADDITIONAL CALCULATION FOR PARTIAL SNOW-COVER CASE.
!    !               SUBROUTINE PATCHS
!    !
!    !     (4) REORGANISATION OF OVERLAND FLOW.
!    !         REFERENCE: SA-89B, APPENDIX.
!    !
!    !     (5) MODIFIED CALCULATION OF SOIL HEAT CAPACITY AND
!    !         CONDUCTIVITY.
!    !
!    !=======================================================================
!    !      1D MODEL SUBROUTINE PATCHS INLINED.
!    !----------------------------------------------------------------------
!    !

!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       ROFF           RUNOFF (MM)
!    !       TC             CANOPY TEMPERATURE (K)
!    !       TG             GROUND SURFACE TEMPERATURE (K)
!    !       WWW(1)         GROUND WETNESS OF SURFACE LAYER
!    !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
!    !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    !     local variables
!    REAL                                                           &
!         PCOEFS(2,2), bp, totalp, ap(len), cp(len), thru(len),    &
!         fpi(len), realc, realg, xss, xsc, chiv, aa, bb,          &
!         exrain, zload, xs, arg, tex(len), tti(len), snowwp(len),  &
!         capacp(len), spechc(len), pinf(len), ts(len), equdep,     &
!         dcap, tsd, ex, dareas, rhs, areas, snowhc, p0(len)
!    INTEGER i, iveg

!    DATA PCOEFS(1,1)/ 20. /, PCOEFS(1,2)/ .206E-8 /,                &
!         PCOEFS(2,1)/ 0.0001 /, PCOEFS(2,2)/ 0.9999 /, BP /20. /
!    !
!    !-----------------------------------------------------------------------
!    !
!    !     PREC ( PI-X )   : EQUATION (C.3), SA-89B
!    !
!    !-----------------------------------------------------------------------
!    !
!    DO i = 1,len
!       roffo(i) = 0.0
!       zmelt(i) = 0.0
!       AP(i) = PCOEFS(2,1)
!       CP(i) = PCOEFS(2,2)
!       TOTALP = (PPC(i) + PPL(i)) * dtt
!       IF( SNOWW(i,1) .GT. 0. .OR. SNOWW(i,2) .GT. 0.      &
!            .OR. TM(i) .LT. TF ) PPC(i) = 0.
!       PPL(i) = TOTALP/dtt - PPC(i)
!       IF(TOTALP.GE.1.E-8) THEN
!          AP(i) = PPC(i)*dtt/TOTALP * PCOEFS(1,1) +           &
!               PPL(i)*dtt/TOTALP * PCOEFS(2,1)
!          CP(i) = PPC(i)*dtt/TOTALP * PCOEFS(1,2) +            &
!               PPL(i)*dtt/TOTALP * PCOEFS(2,2)
!       ENDIF
!       !
!       THRU(i) = 0.
!       FPI(i)  = 0.
!       !
!       !----------------------------------------------------------------------
!       !     PRECIP INPUT INTO INTER2 IN M/SEC; TOTALP IS IN METERS
!       !----------------------------------------------------------------------
!       !
!       P0(i) = TOTALP
!    ENDDO

!    DO IVEG = 1, 2  !DO 1000 IVEG = 1, 2
!       REALC = 2. - IVEG
!       REALG = IVEG - 1.
!       !
!       DO i = 1,len
!          !
!          XSC = MAX(0.E0, CAPAC(i,IVEG) - SATCAP(i,IVEG) )
!          CAPAC(i,IVEG) = CAPAC(i,IVEG) - XSC
!          XSS = MAX(0.E0, SNOWW(i,IVEG) - SATCAP(i,IVEG) ) * REALC
!          SNOWW(i,IVEG) = SNOWW(i,IVEG) - XSS
!          ROFF(i) = ROFF(i) + XSC + XSS

!          CAPACP(i) = CAPAC(i,IVEG)
!          SNOWWP(i) = SNOWW(i,IVEG)
!          !
!          SPECHC(i) =                                                    &
!               MIN( 0.05E0, ( CAPAC(i,IVEG) + SNOWW(i,IVEG) ) ) * CW    &
!               + REALC * ZLT(i) * CLAI + REALG * CSOIL(i)
!          TS(i) = TC(i) * REALC + TG(i) * REALG
!          !
!          !----------------------------------------------------------------------
!          !    PROPORTIONAL SATURATED AREA (XS) AND LEAF DRAINAGE(TEX)
!          !
!          !     TTI ( D-D )     : EQUATION (C.4), SA-89B
!          !     XS  ( X-S )     : EQUATION (C.7), SA-89B
!          !     TEX ( D-C )     : EQUATION (C.8), SA-89B
!          !
!          !----------------------------------------------------------------------
!          !
!          CHIV = CHIL(i)
!          IF ( ABS(CHIV) .LE. 0.01 ) CHIV = 0.01
!          AA = 0.5 - 0.633 * CHIV - 0.33 * CHIV * CHIV
!          BB = 0.877 * ( 1. - 2. * AA )
!          EXRAIN = AA + BB
!          !
!          ZLOAD = CAPAC(i,IVEG) + SNOWW(i,IVEG)
!          FPI(i)=( 1.-EXP( - EXRAIN*ZLT(i)/VCOVER(i) ) )*            &
!               VCOVER(i)*REALC + REALG
!          TTI(i) = P0(i) * ( 1.-FPI(i) )
!          XS = 1.
!          IF ( P0(i) .GE. 1.E-9 ) THEN
!             ARG =  ( SATCAP(i,IVEG)-ZLOAD )/             &
!                  ( P0(i)*FPI(i)*AP(i) ) -CP(i)/AP(i)
!             IF ( ARG .GE. 1.E-9 ) THEN
!                XS = -1./BP * LOG( ARG )
!                XS = MIN( XS, 1.E0 )
!                XS = MAX( XS, 0.E0 )
!             ENDIF
!          ENDIF
!          TEX(i) = P0(i)*FPI(i) *                                &
!               ( AP(i)/BP*( 1.- EXP( -BP*XS )) + CP(i)*XS ) -    &
!               ( SATCAP(i,IVEG) - ZLOAD ) * XS
!          TEX(i) = MAX( TEX(i), 0.E0 )
!       ENDDO
!       !
!       !-------------------------------------------------------------
!       !    TOTAL THROUGHFALL (THRU) AND STORE AUGMENTATION
!       !-----------------------------------------------------------
!       !
!       IF ( IVEG .EQ. 1 ) THEN
!          !
!          !! !DIR$ INLINE
!          DO i = 1,len
!             thru(i) = TTI(i) + TEX(i)
!             PINF(i) = P0(i) - THRU(i)
!             IF( TM(i).GT.TF ) THEN
!                CAPAC(i,IVEG) = CAPAC(i,IVEG) + PINF(i)
!             ELSE
!                SNOWW(i,IVEG) = SNOWW(i,IVEG) + PINF(i)
!             ENDIF
!             !
!             CALL ADJUST(Tc(i), SPECHC(i), CAPACP(i), SNOWWP(i),      &
!                  IVEG, capac(i,1), snoww(i,1), tm(i), tf,             &
!                  snomel, www(i,1), zdepth(i,1),                       &
!                  satcap(i,1), cw, nsib, len )

!             !
!             P0(i) = THRU(i)
!          ENDDO

!       ELSE
!          !
!          !! !DIR$ INLINE
!          DO i = 1,len
!             IF ( TG(i) .GT. TF .AND. SNOWW(i,2) .GT. 0. ) THEN
!                !
!                !=============================================================
!                !------------------PATCHS INLINED-----------------------------
!                !=============================================================
!                !
!                !CALCULATION OF EFFECT OF INTERCEPTED SNOW AND RAINFALL ON GROUND
!                !PATCHY SNOWCOVER SITUATION INVOLVES COMPLEX TREATMENT TO KEEP
!                !ENERGY CONSERVED.
!                !
!                !==============================================================

!                !
!                !MARGINAL SITUATION: SNOW EXISTS IN PATCHES AT TEMPERATURE TF
!                !WITH REMAINING AREA AT TEMPERATURE TG > TF.
!                !
!                !--------------------------------------------------------------
!                !
!                PINF(i) = P0(i)
!                THRU(i) = 0.
!                SNOWHC = MIN( 0.05E0, SNOWW(i,2) ) * CW
!                areas = MIN( 1.E0,(ASNOW*SNOWW(i,2)) )
!                IF( TM(i) .LE. TF ) THEN
!                   !
!                   !-----------------------------------------------------------
!                   !     SNOW FALLING ONTO AREA
!                   !---------------------------------------------------------
!                   !
!                   RHS = TM(i)*PINF(i)*CW + TF*(SNOWHC +          &
!                        CSOIL(i)*areas)                         &
!                        + TG(i)*CSOIL(i)*(1.-areas)
!                   DAREAS = MIN( ASNOW*PINF(i), ( 1.-areas ) )
!                   EX = RHS - TF*PINF(i)*CW -                      &
!                        TF*(SNOWHC + CSOIL(i)*(areas + DAREAS))    &
!                        - TG(i)*CSOIL(i)*(1.-areas-DAREAS)
!                   IF( (areas+DAREAS) .GE. 0.999 )                 &
!                        TG(i) = TF - 0.01
!                   IF( EX .GE. 0. ) THEN
!                      !
!                      !----------------------------------------------------------
!                      !EXCESS ENERGY IS POSITIVE, SOME SNOW MELTS AND INFILTRATES
!                      !----------------------------------------------------------
!                      !
!                      ZMELT(i) = EX/SNOMEL
!                      IF( ASNOW*(SNOWW(i,2) + PINF(i) - ZMELT(i))       &
!                           .LE. 1. ) THEN
!                         ZMELT(i) = 0.
!                         IF( ASNOW*(SNOWW(i,2) + PINF(i)) .GE. 1. )      &
!                              ZMELT(i) = ( ASNOW*(SNOWW(i,2) +          &
!                              PINF(i)) - 1. ) / ASNOW
!                         ZMELT(i) = ( EX - ZMELT(i)*SNOMEL )/            &
!                              ( SNOMEL + ASNOW*CSOIL(i)*                 &
!                              (TG(i)-TF) ) + ZMELT(i)
!                      ENDIF
!                      SNOWW(i,2) =  SNOWW(i,2) + PINF(i) - ZMELT(i)
!                      WWW(i,1) = WWW(i,1) + ZMELT(i)/ZDEPTH(i,1)
!                   ELSE
!                      !
!                      !----------------------------------------------------------
!                      !EXCESS ENERGY IS NEGATIVE,
!                      !BARE GROUND COOLS TO TF, THEN WHOLE
!                      !AREA COOLS TOGETHER TO LOWER TEMPERATURE.
!                      !----------------------------------------------------------
!                      !
!                      TSD = 0.
!                      IF( (areas+DAREAS) .LE. 0.999 )               &
!                           TSD = EX/(CSOIL(i)*( 1.-areas-DAREAS))   &
!                           + TG(i)
!                      IF( TSD .LE. TF )                             &
!                           TSD = TF + ( EX - (TF-TG(i))*              &
!                           CSOIL(i)*(1.-areas-DAREAS) )        &
!                           /(SNOWHC+PINF(i)*CW+CSOIL(i))
!                      TG(i) = TSD
!                      SNOWW(i,2) = SNOWW(i,2) + PINF(i)
!                   ENDIF
!                ELSE
!                   !
!                   !-------------------------------------------------------------
!                   !     RAIN FALLING ONTO AREA
!                   !-------------------------------------------------------------
!                   !
!                   !-------------------------------------------------------------
!                   !     RAIN FALLS ONTO SNOW-FREE SECTOR FIRST.
!                   !-------------------------------------------------------------
!                   TSD = TF - 0.01
!                   IF ( areas .LT. 0.999 ) TSD =                &
!                        ( TM(i)*PINF(i)*CW +               &
!                        TG(i)*CSOIL(i) )                 &
!                        /  ( PINF(i)*CW + CSOIL(i) )
!                   TG(i) = TSD
!                   WWW(i,1)= WWW(i,1)+PINF(i)*(1.-areas)/         &
!                        ZDEPTH(i,1)
!                   !-------------------------------------------------------------
!                   !     RAIN FALLS ONTO SNOW-COVERED SECTOR NEXT.
!                   !-------------------------------------------------------------
!                   EX = ( TM(i) - TF )*PINF(i)*CW*areas
!                   DCAP = -EX / ( SNOMEL + ( TG(i)-TF )*           &
!                        CSOIL(i)*ASNOW )
!                   IF( (SNOWW(i,2) + DCAP) .GE. 0. ) THEN
!                      WWW(i,1) = WWW(i,1)+(PINF(i)*areas-DCAP)/     &
!                           ZDEPTH(i,1)
!                      SNOWW(i,2) = SNOWW(i,2) + DCAP
!                   ELSE
!                      TG(i) = ( EX - SNOMEL*SNOWW(i,2) -             &
!                           ( TG(i)-TF )*CSOIL(i)*areas ) /          &
!                           CSOIL(i) + TG(i)
!                      WWW(i,1)=WWW(i,1)+(SNOWW(i,2)+PINF(i)*         &
!                           areas)/zdepth(i,1)
!                      CAPAC(i,2) = 0.
!                      SNOWW(i,2) = 0.
!                   ENDIF
!                   !
!                ENDIF
!                !
!                !----------------------------------------------------------------
!                !---------------------END OF PATCHS -----------------------------
!                !----------------------------------------------------------------

!             ELSE
!                !
!                THRU(i) = TTI(i) + TEX(i)
!                IF ( TG(i) .LE. TF .OR. TM(i) .LE. TF )        &
!                     THRU(i) = 0.
!                PINF(i) = P0(i) - THRU(i)
!                IF( TM(i) .GT. TF )THEN
!                   CAPAC(i,IVEG) = CAPAC(i,IVEG) + PINF(i)
!                   !
!                   !-------------------------------------------------------------
!                   !
!                   !    INSTANTANEOUS OVERLAND FLOW CONTRIBUTION ( ROFF )
!                   !
!                   !     ROFF( R-I )     : EQUATION (C.13), SA-89B
!                   !
!                   !-------------------------------------------------------------
!                   !
!                   EQUDEP = SATCO(i) * DTT
!                   !
!                   XS = 1.
!                   IF ( THRU(i) .GE. 1.E-9 ) THEN
!                      ARG = EQUDEP / ( THRU(i) * AP(i) ) -CP(i)/AP(i)
!                      IF ( ARG .GE. 1.E-9 ) THEN
!                         XS = -1./BP * LOG( ARG )
!                         XS = MIN( XS, 1.E0 )
!                         XS = MAX( XS, 0.E0 )
!                      ENDIF
!                   ENDIF
!                   ROFFO(i) = THRU(i) * ( AP(i)/BP *                 &
!                        ( 1.-EXP( -BP*XS )) + CP(i)*XS ) -EQUDEP*XS
!                   ROFFO(i) = MAX ( ROFFO(i), 0.E0 )
!                   ROFF(i) = ROFF(i) + ROFFO(i)
!                   WWW(i,1) = WWW(i,1) +                            &
!                        (THRU(i) - ROFFO(i)) / ZDEPTH(i,1)
!                ELSE
!                   SNOWW(i,IVEG) = SNOWW(i,IVEG) + PINF(i)
!                ENDIF
!                !
!                !----------------------------------------------------------------
!                !
!                CALL ADJUST ( Tg(i), SPECHC(i), CAPACP(i), SNOWWP(i),   &
!                     IVEG, capac(i,1), snoww(i,1), tm(i), tf,            &
!                     snomel, www(i,1), zdepth(i,1),                      &
!                     satcap(i,1), cw, nsib, len  )
!                !
!                !----------------------------------------------------------------
!                !
!             ENDIF
!          ENDDO
!       ENDIF   ! if(iveg.eq.1)

!       !     make either all capac or all snow

!       DO i = 1,len
!          IF(capac(i,iveg).GT.snoww(i,iveg)) THEN
!             capac(i,iveg) = capac(i,iveg) + snoww(i,iveg)
!             snoww(i,iveg) = 0.0
!          ELSE
!             snoww(i,iveg) = snoww(i,iveg) + capac(i,iveg)
!             capac(i,iveg) = 0.0
!          ENDIF
!       ENDDO
!       !
!    ENDDO  !1000 CONTINUE
!    !
!    DO i = 1,len
!       exo(i) = 0.0
!    ENDDO
!    DO I = 1, 3  !DO 4001 I = 1, 3
!       DO l = 1,len
!          EXCESS = MAX(0.E0,(WWW(l,I) - 1.))
!          WWW(l,I) = WWW(l,I) - EXCESS
!          exo(l) = exo(l) + EXCESS * ZDEPTH(l,I)

!          !
!          ! Collatz-Bounoua put excess water into runoff according to
!          ! original sib2 offline code .
!          !
!          roff(l) = roff(l) + EXCESS *               &
!               ZDEPTH(l,I)    ! lahouari
!       ENDDO

!    ENDDO !4001 CONTINUE

!    RETURN
!  END SUBROUTINE INTER2

!  !-----------------------------------------------------------------
!  !
!  !======================SUBROUTINE ADJUST================================
!  !
!  SUBROUTINE ADJUST ( TS, SPECHC, CAPACP, SNOWWP, IVEG         &
!       ,     capac,snoww,tm,tf,snomel,www,zdepth                    &
!       ,     satcap,cw,nsib, len)

!    IMPLICIT NONE

!    INTEGER len, nsib, iveg
!    REAL                                                           &
!         capac(nsib,2),snoww(nsib,2),www,zdepth,satcap(len,2),     &
!         ts, spechc, capacp, snowwp, tm, tf, snomel,               &
!         cw
!    !
!    !=======================================================================
!    !
!    !     TEMPERATURE CHANGE DUE TO ADDITION OF PRECIPITATION
!    !
!    !=======================================================================
!    !
!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       TC             CANOPY TEMPERATURE (K)
!    !       WWW(1)         GROUND WETNESS OF SURFACE LAYER
!    !       CAPAC(2)       CANOPY/GROUND LIQUID INTERCEPTION STORE (M)
!    !       SNOWW(2)       CANOPY/GROUND SNOW INTERCEPTION STORE (M)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    !     local variables
!    REAL    &
!         freeze, diff, ccp, cct, tsd, tta, ttb, cca, ccb, ccc, xs
!    FREEZE = 0.
!    DIFF = ( CAPAC(1,IVEG)+SNOWW(1,IVEG) - CAPACP-SNOWWP )*CW

!    DIFF=MAX(DIFF,0.)

!    CCP = SPECHC
!    CCT = SPECHC + DIFF
!    !
!    TSD = ( TS * CCP + TM * DIFF ) / CCT
!    !
!    IF ( ( TS .GT. TF .AND. TM .LE. TF ) .OR.           &
!         ( TS .LE. TF .AND. TM .GT. TF ) )THEN
!       !
!       TTA = TS
!       TTB = TM
!       CCA = CCP
!       CCB = DIFF
!       IF ( TSD .LE. TF ) THEN
!          !
!          !----------------------------------------------------------------------
!          !    FREEZING OF WATER ON CANOPY OR GROUND
!          !----------------------------------------------------------------------
!          !
!          CCC = CAPACP * SNOMEL
!          IF ( TS .LT. TM ) CCC = DIFF * SNOMEL / CW
!          TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT
!          !
!          FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )
!          FREEZE = (MIN ( CCC, FREEZE )) / SNOMEL
!          IF(TSD .GT. TF)TSD = TF - 0.01
!          !
!       ELSE
!          !
!          !----------------------------------------------------------------------
!          !    MELTING OF SNOW ON CANOPY OR GROUND, WATER INFILTRATES.
!          !----------------------------------------------------------------------
!          !
!          CCC = - SNOWW(1,IVEG) * SNOMEL
!          IF ( TS .GT. TM ) CCC = - DIFF * SNOMEL / CW
!          !
!          TSD = ( TTA * CCA + TTB * CCB + CCC ) / CCT
!          !
!          FREEZE = ( TF * CCT - ( TTA * CCA + TTB * CCB ) )
!          FREEZE = (MAX( CCC, FREEZE )) / SNOMEL
!          IF(TSD .LE. TF)TSD = TF - 0.01
!          !
!       ENDIF
!    ENDIF
!    SNOWW(1,IVEG) = SNOWW(1,IVEG) + FREEZE
!    CAPAC(1,IVEG) = CAPAC(1,IVEG) - FREEZE
!    snoww(1,IVEG) = MAX(snoww(1,IVEG),0.0E0)
!    capac(1,IVEG) = MAX(capac(1,IVEG),0.0E0)

!    !
!    XS = MAX( 0.E0, ( CAPAC(1,IVEG) - SATCAP(1,IVEG) ) )
!    IF( SNOWW(1,IVEG) .GE. 0.0000001 ) XS = CAPAC(1,IVEG)
!    WWW = WWW + XS / ZDEPTH
!    CAPAC(1,IVEG) = CAPAC(1,IVEG) - XS
!    TS = TSD

!    !
!    RETURN
!  END SUBROUTINE ADJUST

!  !------------------------------------------------------------
!  !
!  !===================SUBROUTINE BALAN====================================
!  !
!  SUBROUTINE BALAN ( IPLACE, tau, zdepth,www         &
!       ,            capac,ppc                             &
!       ,            ppl,roff,etmass,totwb,radt            &
!       ,            chf,shf,dtt,ect,eci,egs,egi           &
!       ,            hc,hg,heaten,hflux ,snoww,thm,tc,tg,tgs,td,ps,kapa  &
!       ,             len, ioffset, nsoil )

!    IMPLICIT NONE

!    INTEGER len,  iplace, nsoil
!    INTEGER ioffset
!    REAL                                                          &
!         www(len,3),capac(len,2),zdepth(len,3),radt(len,2),       &
!         snoww(len,2), kapa, tau, ppc(len), ppl(len)              &
!         ,    roff(len), etmass(len), totwb(len), chf(len), shf(len)   &
!         ,    dtt, ect(len), eci(len), egs(len), egi(len), hc(len)     &
!         ,    hg(len), heaten(len), hflux(len), thm(len), tc(len)      &
!         ,    tg(len), tgs(len), td(len, nsoil), ps(len)

!    !     local variables
!    INTEGER i, j, igp, jgp, nerror, indxerr(len), n
!    REAL                                                        &
!         endwb(len), errorw(len), pcmeter(len), plmeter(len)    &
!         ,    emeter(len), cbal(len), gbal(len), errore(len)         &
!         ,    zlhs(len), zrhs(len), tm
!    !
!    !=======================================================================
!    !
!    !     ENERGY AND WATER BALANCE CHECK.
!    !
!    !-----------------------------------------------------------------------
!    !
!    IF( IPLACE .EQ. 1 ) THEN
!       !
!       DO i = 1,len
!          ETMASS(i) = 0.
!          ROFF(i)   = 0.
!          !
!          TOTWB(i) = WWW(i,1) * ZDEPTH(i,1)                &
!               + WWW(i,2) * ZDEPTH(i,2)                &
!               + WWW(i,3) * ZDEPTH(i,3)                &
!               + CAPAC(i,1) + CAPAC(i,2) + snoww(i,1) + snoww(i,2)
!       ENDDO
!       !
!    ELSE
!       !
!       nerror = 0
!       DO i = 1,len
!          ENDWB(i) = WWW(i,1) * ZDEPTH(i,1)                         &
!               + WWW(i,2) * ZDEPTH(i,2)                               &
!               + WWW(i,3) * ZDEPTH(i,3)                               &
!               + CAPAC(i,1) + CAPAC(i,2) + snoww(i,1) + snoww(i,2)    &
!               - (PPL(i)+PPC(i))*0.001*dtt                            &
!               + ETMASS(i)*0.001 + ROFF(i)
!          ERRORW(i)= TOTWB(i) - ENDWB(i)
!          pcmeter(i) = ppc(i) * 0.001*dtt
!          plmeter(i) = ppl(i) * 0.001*dtt
!          EMETER(i)= ETMASS(i) * 0.001
!          !itb...trying a different error check, 1.e-6 is in the
!          !itb...noise in 32-bit arithmetic, works fine in 64-bit
!          !         if(abs(errorw(i)).gt.1.e-6) then
!          IF(ABS(errorw(i)).GT.1.e-5*totwb(i))THEN
!             nerror = nerror + 1
!             indxerr(nerror) = i
!          ENDIF
!       ENDDO
!       !
!       DO j = 1,nerror
!          i = indxerr(j)
!          igp = i+ioffset
!          WRITE(6,900) tau,IGP, TOTWB(i), ENDWB(i), ERRORW(i),    &
!               WWW(i,1), WWW(i,2), WWW(i,3),                 &
!               CAPAC(i,1), CAPAC(i,2),snoww(i,1),snoww(i,2), &
!               pcmeter(i),plmeter(i), EMETER(i), ROFF(i)
!       ENDDO
!       !
!       !-----------------------------------------------------------------------
!       !
!       nerror = 0
!       DO i = 1,len
!          CBAL(i) = RADT(i,1) - CHF(i) -                           &
!               (ECT(i)+HC(i)+ECI(i) )/DTT
!          GBAL(i) = RADT(i,2) - SHF(i) - (EGS(i)+HG(i)+EGI(i)         &
!               - HEATEN(i) )/DTT
!          ZLHS(i) = RADT(i,1)+RADT(i,2) - CHF(i) - SHF(i)
!          ZRHS(i) = HFLUX(i) + (ECT(i) + ECI(i) + EGI(i) + EGS(i)    &
!               + HEATEN(i) ) /DTT
!          !
!          ERRORE(i)= ZLHS(i) - ZRHS(i)
!          IF(ABS(errore(i)).GT.1.0) THEN
!             nerror = nerror + 1
!             indxerr(nerror) = i
!          ENDIF
!       ENDDO
!       !
!       DO j = 1,nerror
!          i = indxerr(j)
!          tm = thm(i) * (0.001*ps(i))**kapa
!          igp = i+ioffset
!          WRITE(6,910) tau,IGP, ZLHS(i), ZRHS(i),            &
!               RADT(i,1), RADT(i,2),CHF(i), SHF(i),             &
!               HFLUX(i), ECT(i), ECI(i), EGI(i), EGS(i),        &
!               tm,tc(i),tg(i),tgs(i)
!          WRITE(6,911)(td(i,n),n=nsoil,1,-1)
!          WRITE(6,912) HC(i), HG(i), HEATEN(i), CBAL(i), GBAL(i)
!          WRITE(6,901)  TOTWB(i), ENDWB(i), ERRORW(i),       &
!               WWW(i,1), WWW(i,2), WWW(i,3),            &
!               CAPAC(i,1), CAPAC(i,2),snoww(i,1),snoww(i,2),   &
!               pcmeter(i),plmeter(i), EMETER(i), ROFF(i)
!       ENDDO
!    ENDIF

!  900 FORMAT(//,10X,'** WARNING: WATER BALANCE VIOLATION **  ',//,    &
!           /,1X,'TAU ', F10.2,' AT SIB POINT (I) = ',I5,                  &
!           /,1X,'BEGIN, END, DIFF  ', 2(F10.7,1X),g13.5,                   &
!           /,1X,'WWW,1-3           ', 3(F10.8,1X),                         &
!           /,1X,'CAPAC1-2,snow 1-2 ', 4(g13.5,1X),                         &
!           /,1X,'PPc,PPl, ET, ROFF ', 4(g13.5,1X) )
!  910 FORMAT(//,10X,'** WARNING: ENERGY BALANCE VIOLATION **',//,    &
!           /,1X,'TAU ', F10.2,' AT SIB POINT (I) = ',I5,                 &
!           /,1X,'RHS, LHS              ', 2G13.5,                        &
!           /,1X,'RN1, RN2, CHF, SHF, H ', 5G13.5,                        &
!           /,1X,'ECT, ECI, EGI, EGS    ', 4G13.5,                        &
!           /,1X,'TM,  TC,  TG, TGS     ', 4G13.5,                        &
!           /,1X,'HC        HG          ',  G12.5, 12X, G12.5,            &
!           /,1X,'HEATEN, C-BAL, G-BAL  ', 3G13.5 )
!  911 FORMAT(1x,'TD ', 5G13.5)
!  912 FORMAT(1X,'HC        HG          ',  G12.5, 12X, G12.5,      &
!           /,1X,'HEATEN, C-BAL, G-BAL  ', 3G13.5 )
!  901 FORMAT(10X,'WATER BALANCE'                       &
!           /,1X,'BEGIN, END, DIFF  ', 2(F10.7,1X),g13.5,    &
!           /,1X,'WWW,1-3           ', 3(F10.8,1X),          &
!           /,1X,'CAPAC1-2,snow 1-2 ', 4(g13.5,1X),          &
!           /,1X,'PPc,PPl, ET, ROFF ', 4(g13.5,1X) )
!    !
!    RETURN
!  END SUBROUTINE BALAN

!  !---------------------------------------------------------------
!  !
!  !==================SUBROUTINE RBRD=======================================
!  !
!  SUBROUTINE RBRD(tc,rbc,zlt,z2,u2         &
!       ,             rd,rb,ta,g,rdc,tgs          &
!       ,             len )

!    IMPLICIT NONE

!    !========================================================================
!    !
!    !      CALCULATION OF RB AND RD AS FUNCTIONS OF U2 AND TEMPERATURES
!    !
!    !========================================================================


!    !++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!    !
!    !       RB (GRB)       CANOPY TO CAS AERODYNAMIC RESISTANCE (S M-1)
!    !       RD (GRD)       GROUND TO CAS AERODYNAMIC RESISTANCE (S M-1)
!    !       TA (GTA)       CAS TEMPERATURE (K)
!    !
!    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    INTEGER len
!    REAL tc(len),rbc(len),zlt(len),              &
!        z2(len),u2(len),rd(len),rb(len),        &
!        ta(len),rdc(len),tgs(len),              &
!        temdif(len), g
!    INTEGER i
!    REAL fac, fih, d1(len)
!    !
!    !-----------------------------------------------------------------------
!    !      RB       (RB)       : EQUATION (A9), SE-86
!    !-----------------------------------------------------------------------
!    !
!    DO i = 1,len                      !  loop over gridpoint
!       TEMDIF(i) = MAX( 0.1E0,  TC(i)-TA(i) )
!       FAC = ZLT(i)/890.* (TEMDIF(i)*20.0)**0.25
!       RB(i)  = 1.0/(SQRT(U2(i))/RBC(i)+FAC)
!       !

!       !-----------------------------------------------------------------------
!       !      RD       (RD)       : EQUATION (A15), SE-86
!       !-----------------------------------------------------------------------
!       !
!       TEMDIF(i) = MAX( 0.1E0, TGs(i)-TA(i) )
!       FIH = SQRT( 1.+9.*G*TEMDIF(i)*Z2(i)/(TGS(i)*U2(i)*U2(i)) )
!       RD(i)  = RDC(i) / (U2(i) * FIH)

!    ENDDO
!    !
!    RETURN
!  END SUBROUTINE RBRD
!-----------------------------------------------------------------------------
SUBROUTINE co2_biosource(m1,m2,m3,ia,iz,ja,jz,ng,tend_co2,dn0,rtgt)

  USE mem_sib_co2
  use mem_grid

  IMPLICIT NONE
  INTEGER :: m1,m2,m3,ia,iz,ja,jz,ng,i,j
  real, dimension(m1,m2,m3) :: tend_co2,dn0
  real, dimension(m2,m3) :: rtgt
  real f,dz
  data f /3.666667/  !3.666... =44/12

  DO j=ja,jz
     DO i=ia,iz
        dz = rtgt(i,j)/dzt(2) ! dzt=1/(z(k)-z(k-1))
        tend_co2(2,i,j) = tend_co2(2,i,j) +  f*sib_g(ng, 1)%SRC_CO2(i, j)/(dz*dn0(2,i,j))

     ENDDO
  ENDDO
  RETURN
END SUBROUTINE co2_biosource
