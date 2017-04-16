!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
Module rconstants

   !---------------------------------------------------------------------------------------!
   ! Trigonometric constants                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: pi1        = 3.14159265358979  ! Pi                       [      ---]
   real, parameter :: pii        = 1./pi1            ! 1/Pi                     [      ---]
   real, parameter :: halfpi     = pi1/2.            ! Pi/2                     [      ---]
   real, parameter :: twopi      = pi1* 2.           ! 2 Pi                     [      ---]
   real, parameter :: pio180     = pi1/ 180.         ! Pi/180 (deg -> rad)      [      ---]
   real, parameter :: onerad     = 180. / pi1        ! 180/pi (rad -> deg)      [      ---]
   real, parameter :: pi4        = pi1 * 4.          ! 4 Pi                     [      ---]
   real, parameter :: pi4o3      = pi1 * 4. / 3.     ! 4 Pi / 3                 [      ---]
   real, parameter :: pio4       = pi1 /4.           ! Pi/4                     [      ---]
   real, parameter :: pio6       = pi1 /6.           ! Pi/6                     [      ---]
   real, parameter :: pio6i      = 6.  /pi1          ! 6/Pi                     [      ---]
   real, parameter :: sqrtpii    = 0.564189583547756 ! 1/(pi**0.5)              [      ---]
   real, parameter :: sqrthalfpi = 1.2533141373155   ! (pi/2)**0.5              [      ---]
   real, parameter :: sqrttwopi  = 2. * sqrthalfpi   ! (2*pi)**0.5              [      ---]
   real, parameter :: euler_gam  = 0.577215664901533 ! Euler's constant         [      ---]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Algebraic shortcuts                                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: srtwo     = 1.414213562373095 ! Square root of 2.         [      ---]
   real, parameter :: srthree   = 1.732050807568877 ! Square root of 3.         [      ---]
   real, parameter :: sqrt2o2   = 0.5 * srtwo       ! ½ Square root of 2.       [      ---]
   real, parameter :: srtwoi    = 1./srtwo          ! 1./ Square root of 2.     [      ---]
   real, parameter :: srthreei  = 1./srthree        ! 1./ Square root of 3.     [      ---]
   real, parameter :: onethird  = 1./3.             ! 1/3                       [      ---]
   real, parameter :: twothirds = 2./3.             ! 2/3                       [      ---]
   real, parameter :: onesixth  = 1./6.             ! 1/6                       [      ---]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Universal constants                                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: stefan      = 5.6696e-8       ! Stefan-Boltzmann constant [ W/m²/K^4]
   real, parameter :: boltzmann   = 1.3806503e-23   ! Boltzmann constant        [m²kg/s²/K]
   real, parameter :: avogrado    = 6.02252e+23     ! Avogrado number           [    #/mol]
   real, parameter :: loschmidt   = 2.68719e+25     ! Loschmidt number          [     #/m3]
   real, parameter :: loschcgs    = loschmidt*1.e-6 ! Loschmidt number          [    #/cm3]
   real, parameter :: t00         = 273.15          ! 0°C                       [       °C]
   real, parameter :: rmol        = 8.314510        ! Molar gas constant        [  J/mol/K]
   real, parameter :: rmolcgs     = rmol*1.e7       ! Molar gas constant        [erg/mol/K]
   real, parameter :: volmol      = 0.022710980     ! Molar volume at STP       [       m³]
   real, parameter :: volmoll     = volmol*1.e3     ! Molar volume at STP       [        L]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Molar masses and derived variables                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: mmdry       = 0.02897        ! Mean dry air molar mass    [   kg/mol]
   real, parameter :: mmo2        = 0.0319988      ! Mean O2 molar mass         [   kg/mol]
   real, parameter :: mmo3        = 0.0479982      ! Mean ozone molar mass      [   kg/mol]
   real, parameter :: mmh2o       = 0.01801505     ! Mean water molar mass      [   kg/mol]
   real, parameter :: mmco2       = 0.0440095      ! Mean CO2 molar mass        [   kg/mol]
   real, parameter :: mmdrycgs    = mmdry * 1.e3   ! Mean dry air molar mass    [    g/mol]
   real, parameter :: mmo2cgs     = mmo2  * 1.e3   ! Mean O2 molar mass         [    g/mol]
   real, parameter :: mmo3cgs     = mmo3  * 1.e3   ! Mean ozone molar mass      [    g/mol]
   real, parameter :: mmh2ocgs    = mmh2o * 1.e3   ! Mean water molar mass      [    g/mol]
   real, parameter :: mmco2cgs    = mmco2 * 1.e3   ! Mean CO2 molar mass        [    g/mol]
   real, parameter :: mmdoc       = mmdry/mmco2    ! mmdry/mmco2                [     ----]
   real, parameter :: mmcod       = mmco2/mmdry    ! mmco2/mmdry                [     ----]
   real, parameter :: mmdry1000   = 1000.*mmdry    ! Mean dry air molar mass    [   kg/mol]
   real, parameter :: mmcod1em6   = mmcod * 1.e-6  ! Convert ppm to kgCO2/kgair [     ----]
   real, parameter :: mmdryi      = 1./mmdry       ! 1./mmdry                   [   mol/kg]
   real, parameter :: mmh2oi      = 1./mmh2o       ! 1./mmdry                   [   mol/kg]
   real, parameter :: mmco2i      = 1./mmco2       ! 1./mmco2                   [   mol/kg]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Time conversion units                                                                 !
   !---------------------------------------------------------------------------------------!
   real, parameter :: yr_day  = 365.2425         ! # of days in a year          [   day/yr]
   real, parameter :: day_sec = 86400.           ! # of seconds in a day        [    s/day]
   real, parameter :: day_hr  = 24.              ! # of hours in a day          [   hr/day]
   real, parameter :: hr_sec  = 3600.            ! # of seconds in an hour      [     s/hr]
   real, parameter :: hr_min  = 60.              ! # of minutes in an hour      [   min/hr]
   real, parameter :: min_sec = 60.              ! # of seconds in a minute     [    s/min]
   real, parameter :: yr_sec  = yr_day * day_sec ! # of seconds in a year       [     s/yr]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! General Earth properties                                                              !
   !---------------------------------------------------------------------------------------!
   real, parameter :: vonk      = 0.40        ! Von Kármán constant             [      ---]
   real, parameter :: grav      = 9.80665     ! Gravity acceleration            [     m/s²]
   real, parameter :: gcgs      = grav * 100. ! Gravity acceleration            [    cm/s²]
   real, parameter :: gg        = .5 * grav   ! ½ g                             [     m/s²]
   real, parameter :: erad      = 6370997.    ! Earth radius                    [        m]
   real, parameter :: spcon     = pio180*erad ! One degree of latitude          [        m]
   real, parameter :: spconkm   = spcon*0.001 ! One degree of latitude          [       km]
   real, parameter :: eradi     = 1./erad     ! Inverse of Earth radius         [      1/m]
   real, parameter :: erad2     = erad*2      ! Earth diameter                  [        m]
   real, parameter :: ss60      = 1.8663      ! Polar stereo conversion to 60°  [         ]
   real, parameter :: omega     = 7.292e-5    ! Earth's rotation speed          [    rad/s]
   real, parameter :: viscos    = .15e-4      ! Viscosity coefficient           [         ]
   real, parameter :: solar     = 1.3533e3    ! Solar constant                  [     W/m²]
   real, parameter :: p00       = 1.e5        ! Reference pressure              [       Pa]
   real, parameter :: prefsea   = 101325.     ! Reference sea level pressure    [       Pa]
   real, parameter :: p00i      = 1. / p00    ! 1/p00                           [     1/Pa]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Reference for this block:                                                             !
   ! MU08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,    !
   !        third edition, Academic Press, Amsterdam, 418pp.  (Chapters 3 and 10).         !
   !                                                                                       !
   !     Air diffusion properties. These properties are temperature-dependent in reality,  !
   ! but for simplicity we assume them constants, using the value at 20°C.                 !
   !                                                                                       !
   ! Thermal diffusivity - Computed from equation on page 32 of MU08;                      !
   ! Kinematic viscosity - Computed from equation on page 32 of MU08;                      !
   ! Thermal expansion coefficient - determined by inverting the coefficient at equation   !
   !                                 10.11 (MU08).                                         !
   ! These terms could be easily made function of temperature in the future if needed be.  !
   !---------------------------------------------------------------------------------------!
   real, parameter :: th_diff0  = 1.89e-5     ! Air thermal diffusivity         [     m²/s]
   real, parameter :: dth_diff  = 0.007       ! Temperature dependency slope    [      1/K]
   real, parameter :: kin_visc0 = 1.33e-5     ! Kinematic viscosity             [     m²/s]
   real, parameter :: dkin_visc = 0.007       ! Temperature dependency slope    [      1/K]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Dry air properties                                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rdry    = rmol/mmdry    ! Gas constant for dry air (Ra)   [   J/kg/K]
   real, parameter :: rdryi   = mmdry/rmol    ! 1./Gas const. for dry air (Ra)  [   kg K/J]
   real, parameter :: cpdry   = 3.5 * rdry    ! Spec. heat at constant press.   [   J/kg/K]
   real, parameter :: cvdry   = 2.5 * rdry    ! Spec. heat at constant volume   [   J/kg/K]
   real, parameter :: cpog    = cpdry /grav   ! cp/g                            [      m/K]
   real, parameter :: rocp    = rdry  / cpdry ! Ra/cp                           [     ----]
   real, parameter :: rocv    = rdry  / cvdry ! Ra/Cv                           [     ----]
   real, parameter :: cpocv   = cpdry / cvdry ! Cp/Cv                           [     ----]
   real, parameter :: cpor    = cpdry / rdry  ! Cp/Ra                           [     ----]
   real, parameter :: cvor    = cvdry / rdry  ! Cp/Ra                           [     ----]
   real, parameter :: gocp    = grav / cpdry  ! g/Cp, dry adiabatic lapse rate  [      K/m]
   real, parameter :: gordry  = grav / rdry   ! g/Ra                            [      K/m]
   real, parameter :: cpdryi  = 1. / cpdry    ! 1/Cp                            [   kg K/J]
   real, parameter :: cpdryi4 = 4. * cpdryi   ! 4/Cp                            [   kg K/J]
   real, parameter :: p00or   = p00 / rdry    ! p0 ** (Ra/Cp)                   [   Pa^2/7]
   real, parameter :: p00k    = 26.8269579527 ! p0 ** (Ra/Cp)                   [   Pa^2/7]
   real, parameter :: p00ki   = 1. / p00k     ! p0 ** (-Ra/Cp)                  [  Pa^-2/7]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Water vapour properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rh2o    = rmol/mmh2o  ! Gas const. for water vapour (Rv)  [   J/kg/K]
   real, parameter :: cph2o   = 1859.       ! Heat capacity at const. pres.     [   J/kg/K]
   real, parameter :: cph2oi  = 1. / cph2o  ! Inverse of heat capacity          [   kg K/J]
   real, parameter :: cvh2o   = cph2o-rh2o  ! Heat capacity at const. volume    [   J/kg/K]
   real, parameter :: gorh2o  = grav / rh2o ! g/Rv                              [      K/m]
   real, parameter :: ep      = mmh2o/mmdry ! or Ra/Rv, epsilon                 [    kg/kg]
   real, parameter :: epi     = mmdry/mmh2o ! or Rv/Ra, 1/epsilon               [    kg/kg]
   real, parameter :: epim1   = epi-1.      ! that 0.61 term of virtual temp.   [    kg/kg]
   real, parameter :: toodry  = 1.e-8       ! Minimum acceptable mixing ratio.  [    kg/kg]
   real, parameter :: toowet  = 3.e-2       ! Maximum acceptable mixing ratio.  [    kg/kg]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Liquid water properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: wdns     = 1.000e3    ! Liquid water density              [    kg/m³]
   real, parameter :: wdnsi    = 1./wdns    ! Inverse of liquid water density   [    m³/kg]
   real, parameter :: cliq     = 4.186e3    ! Liquid water specific heat (Cl)   [   J/kg/K]
   real, parameter :: cliqi    = 1./cliq    ! Inverse of water heat capacity    [   kg K/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Ice properties                                                                        !
   !---------------------------------------------------------------------------------------!
   real, parameter :: idns     = 9.167e2      ! "Hard" ice density              [    kg/m³]
   real, parameter :: idnsi    = 1./idns      ! Inverse of ice density          [    m³/kg]
   real, parameter :: fdns     = 2.000e2      ! Frost density                   [    kg/m³]
   real, parameter :: fdnsi    = 1./fdns      ! Inverse of frost density        [    m³/kg]
   real, parameter :: fsdns    = 1.000e2      ! Fresh snow density              [    kg/m³]
   real, parameter :: fsdnsi   = 1./fsdns     ! Inverse of liquid water density [    m³/kg]
   real, parameter :: cice     = 2.093e3      ! Ice specific heat (Ci)          [   J/kg/K]
   real, parameter :: cicei    = 1. / cice    ! Inverse of ice heat capacity    [   kg K/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Phase change properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: t3ple     = 273.16         ! Water triple point temp. (T3)[        K]
   real, parameter :: t3plei    = 1./t3ple       ! 1./T3                        [      1/K]
   real, parameter :: es3ple    = 611.65685464   ! Vapour pressure at T3 (es3)  [       Pa]
   real, parameter :: es3plei   = 1./es3ple      ! 1./es3                       [     1/Pa]
   real, parameter :: epes3ple  = ep * es3ple    ! epsilon × es3                [ Pa kg/kg]
   real, parameter :: rh2ot3ple = rh2o * t3ple   ! Rv × T3                      [     J/kg]
   real, parameter :: alli      = 3.34e5         ! Lat. heat - fusion       (Lf)[     J/kg]
   real, parameter :: alvl3     = 2.50e6         ! Lat. heat - vaporisation (Lv)[     J/kg]
   real, parameter :: alvi3     = alli + alvl3   ! Lat. heat - sublimation  (Ls)[     J/kg]
   real, parameter :: allii     = 1.   / alli    ! 1./Lf                        [     kg/J]
   real, parameter :: aklv      = alvl3 / cpdry  ! Lv/Cp                        [        K]
   real, parameter :: akiv      = alvi3 / cpdry  ! Ls/Cp                        [        K]
   real, parameter :: lvordry   = alvl3 / rdry   ! Lv/Ra                        [        K]
   real, parameter :: lvorvap   = alvl3 / rh2o   ! Lv/Rv                        [        K]
   real, parameter :: lsorvap   = alvi3 / rh2o   ! Ls/Rv                        [        K]
   real, parameter :: lvt3ple   = alvl3 * t3ple  ! Lv × T3                      [   K J/kg]
   real, parameter :: lst3ple   = alvi3 * t3ple  ! Ls × T3                      [   K J/kg]
   real, parameter :: uiicet3   = cice * t3ple   ! u at triple point, only ice  [     J/kg]
   real, parameter :: uiliqt3   = uiicet3 + alli ! u at triple point, only liq. [     J/kg]
   real, parameter :: dcpvl     = cph2o - cliq   ! difference of sp. heat       [   J/kg/K]
   real, parameter :: dcpvi     = cph2o - cice   ! difference of sp. heat       [   J/kg/K]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The following variables are useful when defining the derivatives of theta_il.     !
   !   They correspond to L?(T) - L?' T.                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: del_alvl3 = alvl3 - dcpvl * t3ple
   real, parameter :: del_alvi3 = alvi3 - dcpvi * t3ple
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Tsupercool are defined as temperatures of supercooled liquid water (water vapour)  !
   ! that will cause the internal energy (enthalpy) to be the same as ice at 0K.  It can   !
   ! be used as an offset for temperature when defining internal energy (enthalpy). The    !
   ! next two methods of defining the internal energy for the liquid part:                 !
   !                                                                                       !
   !   Uliq = Mliq [ Cice T3 + Cliq (T - T3) + Lf]                                         !
   !   Uliq = Mliq Cliq (T - Tsupercool_liq)                                               !
   !                                                                                       !
   !   H    = Mliq [ Cice T3 + Cliq (Ts - T3) + Lv3 + (Cpv - Cliq) (Ts-T3) + Cpv (T-T3) ]  !
   !   H    = Mliq Cpv (T - Tsupercool_vap) ]                                              !
   !                                                                                       !
   !     You may be asking yourself why would we have the ice term in the internal energy  !
   ! definition. The reason is that we can think that internal energy is the amount of     !
   ! energy a parcel received to leave the 0K state to reach the current state (or if you  !
   ! prefer the inverse way, Uliq is the amount of energy the parcel would need to lose to !
   ! become solid at 0K.)                                                                  !
   !---------------------------------------------------------------------------------------!
   real, parameter :: tsupercool_liq = t3ple - (uiicet3 + alli ) * cliqi
   real, parameter :: tsupercool_vap = t3ple - (uiicet3 + alvi3) * cph2oi
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Minimum temperature for computing the condensation effect of temperature on       !
   ! theta_il, thetae_iv, and associates. Below this temperature, assuming the latent      !
   ! heats as constants becomes a really bad assumption. See :                             !
   !                                                                                       !
   ! Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential temper- !
   !    ature as a thermodynamic variable in deep atmospheric models. Mon. Wea. Rev.,      !
   !    v. 109, 1094-1102.                                                                 !
   !---------------------------------------------------------------------------------------!
   real, parameter :: ttripoli  = 253.           ! "Tripoli-Cotton" temp. (Ttr) [        K]
   real, parameter :: htripoli  = cpdry*ttripoli ! Sensible enthalpy at T=Ttr   [     J/kg]
   real, parameter :: htripolii = 1./htripoli    ! 1./htripoli                  [     kg/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Lower bounds for turbulence-related variables                                         !
   !---------------------------------------------------------------------------------------!
   real, parameter :: tkmin       = 5.e-4 ! Minimum TKE                         [     J/kg]
   real, parameter :: sigwmin     = 1.e-4 ! Minimum sigma-w                     [      m/s]
   real, parameter :: abslmomin   = 1.e-4 ! Minimum abs value of Obukhov length [        m]
   real, parameter :: ltscalemax  = 1.e5  ! Maximum Lagrangian timescale        [        s]
   real, parameter :: abswltlmin  = 1.e-4 ! Minimum abs value of Theta*         [    K m/s]
   real, parameter :: lturbmin    = 1.e-3 ! Minimum abs value of turb. lenght   [        m]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These are the lower and upper bounds in which we compute exponentials.  This is   !
   ! to avoid overflows and/or underflows when we compute exponentials.                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: lnexp_min = -38.
   real, parameter :: lnexp_max =  38.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These are the just default huge and tiny numbers that are not the actual huge or  !
   ! tiny values from Fortran intrinsic functions, so if you do any numerical operations   !
   ! you will still be fine.                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: huge_num  = 1.e+19
   real, parameter :: tiny_num  = 1.e-19
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These are useful to test numbers with a safe margin.                              !
   !---------------------------------------------------------------------------------------!
   real, parameter :: almost_one = 1.-epsilon(1.)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Tolerance for truncation errors.  "s" are the strict tolerance (literally the     !
   ! machine epsilon), whereas the "r" are the relaxed tolerance (the square root of the   !
   ! strict tolerance).                                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4), parameter :: s_tol_trunc  = epsilon(1.)
   real(kind=4), parameter :: r_tol_trunc  = 100. * s_tol_trunc
   real(kind=8), parameter :: s_tol_trunc8 = epsilon(1.d0)
   real(kind=8), parameter :: r_tol_trunc8 = 1.d2 * s_tol_trunc8
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Carbon-related unit conversions.                                                  !
   !---------------------------------------------------------------------------------------!
   real, parameter :: mol_2_umol     = 1.e6                 ! mol         => µmol
   real, parameter :: umol_2_mol     = 1.e-6                ! µmol        => mol
   real, parameter :: umol_2_kgC     = 1.20107e-8           ! µmol(CO2)   => kg(C)
   real, parameter :: Watts_2_Ein    = 4.6e-6               ! W/m2        => mol/m²/s
   real, parameter :: Ein_2_Watts    = 1./Watts_2_Ein       ! mol/m²/s    => W/m2
   real, parameter :: kgC_2_umol     = 1. / umol_2_kgC      ! kg(C)       => µmol(CO2)
   real, parameter :: kgom2_2_tonoha = 10.                  ! kg(C)/m²    => ton(C)/ha
   real, parameter :: tonoha_2_kgom2 = 0.1                  ! ton(C)/ha   => kg(C)/m²
   real, parameter :: umols_2_kgCyr  = umol_2_kgC * yr_sec  ! µmol(CO2)/s => kg(C)/yr
   real, parameter :: kgCday_2_umols = kgC_2_umol / day_sec ! kg(C)/day   => µmol(CO2)/s
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Double precision version of all constants used in Runge-Kutta.                     !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: pi18            = dble(pi1           )
   real(kind=8), parameter :: halfpi8         = dble(halfpi        )
   real(kind=8), parameter :: twopi8          = dble(twopi         )
   real(kind=8), parameter :: sqrtpii8        = dble(sqrtpii       )
   real(kind=8), parameter :: pio1808         = dble(pio180        )
   real(kind=8), parameter :: pi48            = dble(pi4           )
   real(kind=8), parameter :: pio48           = dble(pio4          )
   real(kind=8), parameter :: srtwo8          = dble(srtwo         )
   real(kind=8), parameter :: srthree8        = dble(srthree       )
   real(kind=8), parameter :: sqrt2o28        = dble(sqrt2o2       )
   real(kind=8), parameter :: srtwoi8         = dble(srtwoi        )
   real(kind=8), parameter :: srthreei8       = dble(srthreei      )
   real(kind=8), parameter :: onethird8       = dble(onethird      )
   real(kind=8), parameter :: twothirds8      = dble(twothirds     )
   real(kind=8), parameter :: onesixth8       = dble(onesixth      )
   real(kind=8), parameter :: stefan8         = dble(stefan        )
   real(kind=8), parameter :: boltzmann8      = dble(boltzmann     )
   real(kind=8), parameter :: t008            = dble(t00           )
   real(kind=8), parameter :: rmol8           = dble(rmol          )
   real(kind=8), parameter :: volmol8         = dble(volmol        )
   real(kind=8), parameter :: volmoll8        = dble(volmoll       )
   real(kind=8), parameter :: mmdry8          = dble(mmdry         )
   real(kind=8), parameter :: mmh2o8          = dble(mmh2o         )
   real(kind=8), parameter :: mmo28           = dble(mmo2          )
   real(kind=8), parameter :: mmo38           = dble(mmo3          )
   real(kind=8), parameter :: mmco28          = dble(mmco2         )
   real(kind=8), parameter :: mmdoc8          = dble(mmdoc         )
   real(kind=8), parameter :: mmcod8          = dble(mmcod         )
   real(kind=8), parameter :: mmdry10008      = dble(mmdry1000     )
   real(kind=8), parameter :: mmcod1em68      = dble(mmcod1em6     )
   real(kind=8), parameter :: mmdryi8         = dble(mmdryi        )
   real(kind=8), parameter :: mmh2oi8         = dble(mmh2oi        )
   real(kind=8), parameter :: mmco2i8         = dble(mmco2i        )
   real(kind=8), parameter :: yr_day8         = dble(yr_day        )
   real(kind=8), parameter :: day_sec8        = dble(day_sec       )
   real(kind=8), parameter :: day_hr8         = dble(day_hr        )
   real(kind=8), parameter :: hr_sec8         = dble(hr_sec        )
   real(kind=8), parameter :: min_sec8        = dble(min_sec       )
   real(kind=8), parameter :: vonk8           = dble(vonk          )
   real(kind=8), parameter :: grav8           = dble(grav          )
   real(kind=8), parameter :: erad8           = dble(erad          )
   real(kind=8), parameter :: erad28          = dble(erad2         )
   real(kind=8), parameter :: p008            = dble(p00           )
   real(kind=8), parameter :: p00i8           = dble(p00i          )
   real(kind=8), parameter :: p00k8           = dble(p00k          )
   real(kind=8), parameter :: p00ki8          = dble(p00ki         )
   real(kind=8), parameter :: rdry8           = dble(rdry          )
   real(kind=8), parameter :: rdryi8          = dble(rdryi         )
   real(kind=8), parameter :: cpdry8          = dble(cpdry         )
   real(kind=8), parameter :: cvdry8          = dble(cvdry         )
   real(kind=8), parameter :: cpog8           = dble(cpog          )
   real(kind=8), parameter :: rocp8           = dble(rocp          )
   real(kind=8), parameter :: rocv8           = dble(rocv          )
   real(kind=8), parameter :: cpocv8          = dble(cpocv         )
   real(kind=8), parameter :: cpor8           = dble(cpor          )
   real(kind=8), parameter :: cpdryi8         = dble(cpdryi        )
   real(kind=8), parameter :: cpdryi48        = dble(cpdryi4       )
   real(kind=8), parameter :: rh2o8           = dble(rh2o          )
   real(kind=8), parameter :: cph2o8          = dble(cph2o         )
   real(kind=8), parameter :: cph2oi8         = dble(cph2oi        )
   real(kind=8), parameter :: cvh2o8          = dble(cvh2o         )
   real(kind=8), parameter :: gorh2o8         = dble(gorh2o        )
   real(kind=8), parameter :: ep8             = dble(ep            )
   real(kind=8), parameter :: epi8            = dble(epi           )
   real(kind=8), parameter :: epim18          = dble(epim1         )
   real(kind=8), parameter :: toodry8         = dble(toodry        )
   real(kind=8), parameter :: wdns8           = dble(wdns          )
   real(kind=8), parameter :: wdnsi8          = dble(wdnsi         )
   real(kind=8), parameter :: cliq8           = dble(cliq          )
   real(kind=8), parameter :: cliqi8          = dble(cliqi         )
   real(kind=8), parameter :: idns8           = dble(idns          )
   real(kind=8), parameter :: idnsi8          = dble(idnsi         )
   real(kind=8), parameter :: fdns8           = dble(fdns          )
   real(kind=8), parameter :: fdnsi8          = dble(fdnsi         )
   real(kind=8), parameter :: fsdns8          = dble(fsdns         )
   real(kind=8), parameter :: fsdnsi8         = dble(fsdnsi        )
   real(kind=8), parameter :: cice8           = dble(cice          )
   real(kind=8), parameter :: cicei8          = dble(cicei         )
   real(kind=8), parameter :: t3ple8          = dble(t3ple         )
   real(kind=8), parameter :: t3plei8         = dble(t3plei        )
   real(kind=8), parameter :: es3ple8         = dble(es3ple        )
   real(kind=8), parameter :: es3plei8        = dble(es3plei       )
   real(kind=8), parameter :: epes3ple8       = dble(epes3ple      )
   real(kind=8), parameter :: alvl38          = dble(alvl3         )
   real(kind=8), parameter :: alvi38          = dble(alvi3         )
   real(kind=8), parameter :: alli8           = dble(alli          )
   real(kind=8), parameter :: allii8          = dble(allii         )
   real(kind=8), parameter :: akiv8           = dble(akiv          )
   real(kind=8), parameter :: aklv8           = dble(aklv          )
   real(kind=8), parameter :: uiicet38        = dble(uiicet3       )
   real(kind=8), parameter :: uiliqt38        = dble(uiliqt3       )
   real(kind=8), parameter :: dcpvl8          = dble(dcpvl         )
   real(kind=8), parameter :: dcpvi8          = dble(dcpvi         )
   real(kind=8), parameter :: del_alvl38      = dble(del_alvl3     )
   real(kind=8), parameter :: del_alvi38      = dble(del_alvi3     )
   real(kind=8), parameter :: tsupercool_liq8 = dble(tsupercool_liq)
   real(kind=8), parameter :: tsupercool_vap8 = dble(tsupercool_vap)
   real(kind=8), parameter :: ttripoli8       = dble(ttripoli      )
   real(kind=8), parameter :: htripoli8       = dble(htripoli      )
   real(kind=8), parameter :: htripolii8      = dble(htripolii     )
   real(kind=8), parameter :: tkmin8          = dble(tkmin         )
   real(kind=8), parameter :: sigwmin8        = dble(sigwmin       )
   real(kind=8), parameter :: abslmomin8      = dble(abslmomin     )
   real(kind=8), parameter :: ltscalemax8     = dble(ltscalemax    )
   real(kind=8), parameter :: abswltlmin8     = dble(abswltlmin    )
   real(kind=8), parameter :: lturbmin8       = dble(lturbmin      )
   real(kind=8), parameter :: th_diff08       = dble(th_diff0      )
   real(kind=8), parameter :: dth_diff8       = dble(dth_diff      )
   real(kind=8), parameter :: kin_visc08      = dble(kin_visc0     )
   real(kind=8), parameter :: dkin_visc8      = dble(dkin_visc     )
   real(kind=8), parameter :: lnexp_min8      = dble(lnexp_min     )
   real(kind=8), parameter :: lnexp_max8      = dble(lnexp_max     )
   real(kind=8), parameter :: huge_num8       = dble(huge_num      )
   real(kind=8), parameter :: tiny_num8       = dble(tiny_num      )
   real(kind=8), parameter :: almost_one8     = dble(almost_one    )
   real(kind=8), parameter :: euler_gam8      = dble(euler_gam     )
   real(kind=8), parameter :: mol_2_umol8     = dble(mol_2_umol    )
   real(kind=8), parameter :: umol_2_mol8     = dble(umol_2_mol    )
   real(kind=8), parameter :: umol_2_kgC8     = dble(umol_2_kgC    )
   real(kind=8), parameter :: Watts_2_Ein8    = dble(Watts_2_Ein   )
   real(kind=8), parameter :: Ein_2_Watts8    = dble(Ein_2_Watts   )
   !---------------------------------------------------------------------------------------!



end module rconstants
!==========================================================================================!
!==========================================================================================!
