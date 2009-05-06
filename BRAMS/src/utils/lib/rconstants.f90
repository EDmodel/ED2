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
   real, parameter :: idns     = 9.167e2      ! "Hard" ice density              [    kg/m³]
   real, parameter :: idnsi    = 1./idns      ! Inverse of ice density          [    m³/kg]
   real, parameter :: cice     = 2.093e3      ! Ice specific heat (Ci)          [   J/kg/K]
   real, parameter :: cicevlme = wdns * cice  ! Heat capacity × water density   [   J/m³/K]
   real, parameter :: cicei    = 1. / cice    ! Inverse of ice heat capacity    [   kg K/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Phase change properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: t3ple    = 273.16        ! Water triple point temp. (T3)  [        K]
   real, parameter :: t3plei   = 1./t3ple      ! 1./T3                          [      1/K]
   real, parameter :: es3ple   = 611.65685464  ! Vapour pressure at T3 (es3)    [       Pa]
   real, parameter :: es3plei  = 1./es3ple     ! 1./es3                         [     1/Pa]
   real, parameter :: epes3ple = ep * es3ple   ! epsilon × es3                  [ Pa kg/kg]
   real, parameter :: rmt3ple  = rm * t3ple    ! Rv × T3                        [     J/kg]
   real, parameter :: alvl     = 2.50e6        ! Lat. heat - vaporisation (Lv)  [     J/kg]
   real, parameter :: alvi     = 2.834e6       ! Lat. heat - sublimation  (Ls)  [     J/kg]
   real, parameter :: alli     = 3.34e5        ! Lat. heat - fusion       (Lf)  [     J/kg]
   real, parameter :: allivlme = wdns * alli   ! Lat. heat × water density      [     J/m³]
   real, parameter :: alvl2    = alvl * alvl   ! Lv²                            [   J²/kg²]
   real, parameter :: alvi2    = alvi * alvi   ! Ls²                            [   J²/kg²]
   real, parameter :: allii    = 1.   / alli   ! 1./Lf                          [     kg/J]
   real, parameter :: aklv     = alvl / cp     ! Lv/Cp                          [        K]
   real, parameter :: akiv     = alvi / cp     ! Ls/Cp                          [        K]
   real, parameter :: lvordry  = alvl / rgas   ! Lv/Ra                          [        K]
   real, parameter :: lvorvap  = alvl / rm     ! Lv/Rv                          [        K]
   real, parameter :: lsorvap  = alvi / rm     ! Ls/Rv                          [        K]
   real, parameter :: lvt3ple  = alvl * t3ple  ! Lv × T3                        [   K J/kg]
   real, parameter :: lst3ple  = alvi * t3ple  ! Ls × T3                        [   K J/kg]
   real, parameter :: qicet3   = cice * t3ple  ! q at triple point, only ice    [     J/kg]
   real, parameter :: qliqt3   = qicet3 + alli ! q at triple point, only liquid [     J/kg]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Tsupercool is the temperature of supercooled water that will cause the energy to   !
   ! be the same as ice at 0K. It can be used as an offset for temperature when defining   !
   ! internal energy. The next two methods of defining the internal energy for the liquid  !
   ! part:                                                                                 !
   !                                                                                       !
   !   Uliq = Mliq × [ Cice × T3 + Cliq × (T - T3) + Lf]                                   !
   !   Uliq = Mliq × Cliq × (T - Tsupercool)                                               !
   !                                                                                       !
   !     You may be asking yourself why would we have the ice term in the internal energy  !
   ! definition. The reason is that we can think that internal energy is the amount of     !
   ! energy a parcel received to leave the 0K state to reach the current state (or if you  !
   ! prefer the inverse way, Uliq is the amount of energy the parcel would need to lose to !
   ! become solid at 0K.)                                                                  !
   !---------------------------------------------------------------------------------------!
   real, parameter :: tsupercool = t3ple - (qicet3+alli) * cliqi
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



   !---------------------------------------------------------------------------------------!
   !    Double precision version of constants used in Runge-Kutta.                         !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter :: alli8        = dble(alli      )
   real(kind=8), parameter :: allii8       = dble(allii     )
   real(kind=8), parameter :: alvi8        = dble(alvi      )
   real(kind=8), parameter :: alvl8        = dble(alvl      )
   real(kind=8), parameter :: cice8        = dble(cice      )
   real(kind=8), parameter :: cicei8       = dble(cicei     )
   real(kind=8), parameter :: cliq8        = dble(cliq      )
   real(kind=8), parameter :: cliqi8       = dble(cliqi     )
   real(kind=8), parameter :: cliqvlme8    = dble(cliqvlme  )
   real(kind=8), parameter :: cp8          = dble(cp        )
   real(kind=8), parameter :: cpi8         = dble(cpi       )
   real(kind=8), parameter :: day_sec8     = dble(day_sec   )
   real(kind=8), parameter :: gorvap8      = dble(gorm      )
   real(kind=8), parameter :: grav8        = dble(g         )
   real(kind=8), parameter :: idns8        = dble(idns      )
   real(kind=8), parameter :: pi18         = dble(pi1       )
   real(kind=8), parameter :: qicet38      = dble(qicet3    )
   real(kind=8), parameter :: qliqt38      = dble(qliqt3    )
   real(kind=8), parameter :: stefan8      = dble(stefan    )
   real(kind=8), parameter :: t3ple8       = dble(t3ple     )
   real(kind=8), parameter :: tsupercool8  = dble(tsupercool)
   real(kind=8), parameter :: twothirds8   = dble(twothirds )
   real(kind=8), parameter :: vonk8        = dble(vonk      )
   real(kind=8), parameter :: wdns8        = dble(wdns      )
   real(kind=8), parameter :: wdnsi8       = dble(wdnsi     )
   !---------------------------------------------------------------------------------------!

end module rconstants
