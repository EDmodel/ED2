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
   real, parameter :: pi1       = 3.14159265358979  ! Pi                        [      ---]
   real, parameter :: pii       = 1./pi1            ! 1/Pi                      [      ---]
   real, parameter :: halfpi    = pi1/2             ! Pi/2                      [      ---]
   real, parameter :: sqrtpii   = 0.564189583547756 ! 1/(pi**0.5)               [      ---]
   real, parameter :: twopi     = pi1* 2.           ! 2 Pi                      [      ---]
   real, parameter :: pio180    = pi1/ 180.         ! Pi/180 (deg -> rad)       [      ---]
   real, parameter :: onerad    = 180. / pi1        ! 180/pi (rad -> deg)       [      ---]
   real, parameter :: pi4       = pi1 * 4.          ! 4 Pi                      [      ---]
   real, parameter :: pio4      = pi1 /4.           ! Pi/4                      [      ---]
   real, parameter :: pio6      = pi1 /6.           ! Pi/6                      [      ---]
   real, parameter :: pio6i     = 6.  /pi1          ! 6/Pi                      [      ---]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Algebraic shortcuts                                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: srtwo     = 1.414213562373095 ! Square root of 2.         [      ---]
   real, parameter :: srthree   = 1.732050807568877 ! Square root of 3.         [      ---]
   real, parameter :: srtwoi    = 1./srtwo          ! 1./ Square root of 2.     [      ---]
   real, parameter :: srthreei  = 1./srthree        ! 1./ Square root of 3.     [      ---]
   real, parameter :: onethird  = 1./3.             ! 1/3                       [      ---]
   real, parameter :: onesixth  = 1./6.             ! 1/6                       [      ---]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Universal constants                                                                   !
   !---------------------------------------------------------------------------------------!
   real, parameter :: stefan    = 5.6696e-8         ! Stefan-Boltzmann constant [ W/m²/K^4]
   real, parameter :: boltzmann = 1.3806503e-23     ! Boltzmann constant        [m²kg/s²/K]
   real, parameter :: t00       = 273.15            ! 0°C                       [       °C]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Time conversion units                                                                 !
   !---------------------------------------------------------------------------------------!
   real, parameter :: yr_day  = 365.2425 ! # of days in a year                  [   day/yr]
   real, parameter :: day_sec = 86400.   ! # of seconds in a day                [    s/day]
   real, parameter :: day_hr  = 24.      ! # of hours in a day                  [   hr/day]
   real, parameter :: hr_sec  = 3600.    ! # of seconds in an hour              [     s/hr]
   real, parameter :: min_sec = 60.      ! # of seconds in a minute             [    s/min]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! General Earth properties                                                              !
   !---------------------------------------------------------------------------------------!
   real, parameter :: vonk      = 0.40        ! Von Kármán constant             [      ---]
   real, parameter :: g         = 9.80665     ! Gravity acceleration            [     m/s²]
   real, parameter :: gg        = .5 * g      ! ½ g                             [     m/s²]
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
   ! Dry air properties                                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rgas   = 287.04    ! Gas constant for dry air (Ra)        [   J/kg/K]
   real, parameter :: cp     = 1004.     ! Specific heat at constant pressure   [   J/kg/K]
   real, parameter :: cv     = 717.      ! Specific heat at constant volume     [   J/kg/K]
   real, parameter :: cpog   = cp /g     ! cp/g                                 [      m/K]
   real, parameter :: rocp   = rgas / cp ! Ra/cp                                [     ----]
   real, parameter :: cpor   = cp / rgas ! Cp/Ra                                [     ----]
   real, parameter :: rocv   = rgas / cv ! Ra/Cv                                [     ----]
   real, parameter :: gocp   = g / cp    ! g/Cp, dry adiabatic lapse rate       [      K/m]
   real, parameter :: gordry = g / rgas  ! g/Ra                                 [      K/m]
   real, parameter :: cpi    = 1. / cp   ! 1/Cp                                 [   kg K/J]
   real, parameter :: cpi4   = 4. * cpi  ! 4/Cp                                 [   kg K/J]
   real, parameter :: p00k   = 26.870941 ! p0 ** (Ra/Cp)                        [ Pa^0.286]
   real, parameter :: p00ki  = 1. / p00k ! p0 ** (-Ra/Cp)                       [Pa^-0.286]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Water vapour properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rm     = 461.5     ! Gas constant for water vapour (Rv)   [   J/kg/K]
   real, parameter :: gorm   = g / rm    ! g/Rv                                 [      K/m]
   real, parameter :: ep     = rgas / rm ! Ra/Rv, epsilon, used to find rv      [    kg/kg]
   real, parameter :: epi    = rm / rgas ! Rv/Ra, 1/epsilon                     [    kg/kg]
   real, parameter :: rmocp  = rm / cp   ! Rv/cp                                [     ----]
   real, parameter :: toodry = 1.e-8     ! Minimum acceptable mixing ratio.     [    kg/kg]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Liquid water properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: wdns     = 1.000e3    ! Liquid water density              [    kg/m³]
   real, parameter :: wdnsi    = 1./wdns    ! Inverse of liquid water density   [    m³/kg]
   real, parameter :: cliq     = 4.186e3    ! Liquid water specific heat (Cl)   [   J/kg/K]
   real, parameter :: cliqvlme = wdns*cliq  ! Water heat capacity × water dens. [   J/m³/K]
   real, parameter :: cliqi    = 1./cliq    ! Inverse of water heat capacity    [   kg K/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Ice properties                                                                        !
   !---------------------------------------------------------------------------------------!
   real, parameter :: cice     = 2.093e3      ! Ice specific heat (Ci)          [   J/kg/K]
   real, parameter :: cicevlme = wdns * cice  ! Heat capacity × water density   [   J/m³/K]
   real, parameter :: cicei    = 1. / cice    ! Inverse of ice heat capacity    [   kg K/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Phase change properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: t3ple    = 273.16       ! Water triple point temp. (T3)   [        K]
   real, parameter :: t3plei   = 1./t3ple     ! 1./T3                           [      1/K]
   real, parameter :: es3ple   = 611.65685464 ! Vapour pressure at T3 (es3)     [       Pa]
   real, parameter :: es3plei  = 1./es3ple    ! 1./es3                          [     1/Pa]
   real, parameter :: epes3ple = ep * es3ple  ! epsilon × es3                   [ Pa kg/kg]
   real, parameter :: rmt3ple  = rm * t3ple   ! Rv × T3                         [     J/kg]
   real, parameter :: alvl     = 2.50e6       ! Latent heat - vaporisation (Lv) [     J/kg]
   real, parameter :: alvi     = 2.834e6      ! Latent heat - sublimation  (Ls) [     J/kg]
   real, parameter :: alli     = 3.34e5       ! Latent heat - fusion       (Lf) [     J/kg]
   real, parameter :: allivlme = wdns * alli  ! Latent heat × water density     [     J/m³]
   real, parameter :: alvl2    = alvl*alvl    ! Lv²                             [   J²/kg²]
   real, parameter :: alvi2    = alvi*alvi    ! Ls²                             [   J²/kg²]
   real, parameter :: allii    = 1. / alli    ! 1./Lf                           [     kg/J]
   real, parameter :: aklv     = alvl / cp    ! Lv/Cp                           [        K]
   real, parameter :: akiv     = alvi / cp    ! Ls/Cp                           [        K]
   real, parameter :: lvordry  = alvl / rgas  ! Lv/Ra                           [        K]
   real, parameter :: lvorvap  = alvl / rm    ! Lv/Rv                           [        K]
   real, parameter :: lsorvap  = alvi / rm    ! Ls/Rv                           [        K]
   real, parameter :: lvt3ple  = alvl * t3ple ! Lv × T3                         [   K J/kg]
   real, parameter :: lst3ple  = alvi * t3ple ! Ls × T3                         [   K J/kg]
   real, parameter :: cicet3   = cice * t3ple ! C_ice × T3                      [     J/kg]
   real, parameter :: cliqt3   = cliq * t3ple ! C_liquid × T3                   [     J/kg]
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
   real, parameter :: ttripoli  = 253.        ! "Tripoli-Cotton" temp. (Ttr)    [        K]
   real, parameter :: htripoli  = cp*ttripoli ! Sensible enthalpy at T=Ttr      [     J/kg]
   real, parameter :: htripolii = 1./htripoli ! 1./htripoli                     [     kg/J]
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
end module rconstants
