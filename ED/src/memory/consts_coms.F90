Module consts_coms

! Making the constants to match when running a coupled run.

#if defined(COUPLED)
   use rconstants, only:                                                                   &
       b_pi1        => pi1        , b_twopi      => twopi      , b_pio180     => pio180    &
     , b_pi4        => pi4        , b_pio4       => pio4       , b_srtwo      => srtwo     &
     , b_srthree    => srthree    , b_srtwoi     => srtwoi     , b_srthreei   => srthreei  &
     , b_onethird   => onethird   , b_stefan     => stefan     , b_boltzmann  => boltzmann &
     , b_t00        => t00        , b_yr_day     => yr_day     , b_day_sec    => day_sec   &
     , b_day_hr     => day_hr     , b_hr_sec     => hr_sec     , b_min_sec    => min_sec   &
     , b_vonk       => vonk       , b_g          => g          , b_erad       => erad      &
     , b_p00        => p00        , b_p00i       => p00i       , b_rgas       => rgas      &
     , b_cp         => cp         , b_cpog       => cpog       , b_rocp       => rocp      &
     , b_cpor       => cpor       , b_cpi        => cpi        , b_rm         => rm        &
     , b_ep         => ep         , b_epi        => epi        , b_toodry     => toodry    &
     , b_cliq       => cliq       , b_cliqvlme   => cliqvlme   , b_cliqi      => cliqi     &
     , b_cice       => cice       , b_cicevlme   => cicevlme   , b_cicei      => cicei     &
     , b_t3ple      => t3ple      , b_t3plei     => t3plei     , b_es3ple     => es3ple    &
     , b_es3plei    => es3plei    , b_epes3ple   => epes3ple   , b_alvl       => alvl      &
     , b_alvi       => alvi       , b_alli       => alli       , b_allivlme   => allivlme  &
     , b_allii      => allii      , b_wdns       => wdns       , b_erad2      => erad2     &
     , b_sqrtpii    => sqrtpii    , b_onesixth   => onesixth   , b_qicet3     => qicet3    &
     , b_wdnsi      => wdnsi      , b_gorm       => gorm       , b_idns       => idns      &
     , b_idnsi      => idnsi      , b_tsupercool => tsupercool , b_twothirds  => twothirds &
     , b_qliqt3     => qliqt3     , b_sqrt2o2    => sqrt2o2

   implicit none

   real, parameter :: pi1        = b_pi1        , twopi      = b_twopi
   real, parameter :: pio180     = b_pio180     , pi4        = b_pi4
   real, parameter :: pio4       = b_pio4       , srtwo      = b_srtwo
   real, parameter :: srthree    = b_srthree    , srtwoi     = b_srtwoi
   real, parameter :: srthreei   = b_srthreei   , onethird   = b_onethird
   real, parameter :: twothirds  = b_twothirds  , stefan     = b_stefan
   real, parameter :: boltzmann  = b_boltzmann  , tsupercool = b_tsupercool
   real, parameter :: t00        = b_t00        , yr_day     = b_yr_day
   real, parameter :: day_sec    = b_day_sec    , day_hr     = b_day_hr
   real, parameter :: hr_sec     = b_hr_sec     , min_sec    = b_min_sec
   real, parameter :: vonk       = b_vonk       , grav       = b_g
   real, parameter :: erad       = b_erad       , p00        = b_p00
   real, parameter :: p00i       = b_p00i       , rdry       = b_rgas
   real, parameter :: cp         = b_cp         , cpog       = b_cpog
   real, parameter :: rocp       = b_rocp       , cpor       = b_cpor
   real, parameter :: cpi        = b_cpi        , rvap       = b_rm
   real, parameter :: ep         = b_ep         , epi        = b_epi
   real, parameter :: toodry     = b_toodry     , cliq       = b_cliq
   real, parameter :: cliqvlme   = b_cliqvlme   , cliqi      = b_cliqi
   real, parameter :: cice       = b_cice       , cicevlme   = b_cicevlme
   real, parameter :: cicei      = b_cicei      , t3ple      = b_t3ple
   real, parameter :: t3plei     = b_t3plei     , es3ple     = b_es3ple
   real, parameter :: es3plei    = b_es3plei    , epes3ple   = b_epes3ple
   real, parameter :: alvl       = b_alvl       , alvi       = b_alvi
   real, parameter :: alli       = b_alli       , allivlme   = b_allivlme
   real, parameter :: allii      = b_allii      , wdns       = b_wdns
   real, parameter :: erad2      = b_erad2      , sqrtpii    = b_sqrtpii
   real, parameter :: onesixth   = b_onesixth   , qicet3     = b_qicet3
   real, parameter :: wdnsi      = b_wdnsi      , gorvap     = b_gorm
   real, parameter :: idns       = b_idns       , idnsi      = b_idnsi
   real, parameter :: qliqt3     = b_qliqt3     , sqrt2o2    = b_sqrt2o2
#else
   implicit none

   !---------------------------------------------------------------------------------------!
   ! Trigonometric constants                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: pi1       = 3.14159265358979  ! Pi                        [      ---]
   real, parameter :: twopi     = pi1* 2.           ! 2 Pi                      [      ---]
   real, parameter :: sqrtpii   = 0.564189583547756 ! 1/(pi**0.5)               [      ---]
   real, parameter :: pio180    = pi1/ 180.         ! Pi/180 (deg -> rad)       [      ---]
   real, parameter :: pi4       = pi1 * 4.          ! 4 Pi                      [      ---]
   real, parameter :: pio4      = pi1 /4.           ! Pi/4                      [      ---]
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
   real, parameter :: grav      = 9.80665     ! Gravity acceleration            [     m/s²]
   real, parameter :: erad      = 6370997.    ! Earth radius                    [        m]
   real, parameter :: erad2     = 2.*erad     ! Earth diameter                  [        m]
   real, parameter :: p00       = 1.e5        ! Reference pressure              [       Pa]
   real, parameter :: p00i      = 1. / p00    ! 1/p00                           [     1/Pa]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Dry air properties                                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rdry   = 287.04    ! Gas constant for dry air (Ra)        [   J/kg/K]
   real, parameter :: cp     = 1004.     ! Specific heat at constant pressure   [   J/kg/K]
   real, parameter :: cpog   = cp /grav  ! cp/g                                 [      m/K]
   real, parameter :: rocp   = rdry / cp ! Ra/cp                                [     ----]
   real, parameter :: cpor   = cp / rdry ! Cp/Ra                                [     ----]
   real, parameter :: cpi    = 1. / cp   ! 1/Cp                                 [   kg K/J]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Water vapour properties                                                               !
   !---------------------------------------------------------------------------------------!
   real, parameter :: rvap   = 461.5       ! Gas constant for water vapour (Rv) [   J/kg/K]
   real, parameter :: gorvap = grav / rvap ! g/Rv                               [      K/m]
   real, parameter :: ep     = rdry / rvap ! Ra/Rv, epsilon, used to find rv    [    kg/kg]
   real, parameter :: epi    = rvap / rdry ! Rv/Ra, 1/epsilon                   [    kg/kg]
   real, parameter :: toodry = 1.e-8       ! Minimum acceptable mixing ratio.   [    kg/kg]
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
   real, parameter :: alvl     = 2.50e6        ! Lat. heat - vaporisation (Lv)  [     J/kg]
   real, parameter :: alvi     = 2.834e6       ! Lat. heat - sublimation  (Ls)  [     J/kg]
   real, parameter :: alli     = 3.34e5        ! Lat. heat - fusion       (Lf)  [     J/kg]
   real, parameter :: allivlme = wdns * alli   ! Lat. heat × water density      [     J/m³]
   real, parameter :: allii    = 1./alli       ! 1/Latent heat - fusion         [     kg/J]
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

#endif

   !---------------------------------------------------------------------------------------!
   ! Unit conversion, it must be defined locally even for coupled runs.                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: umol_2_kgC      = 1.20107e-8 ! µmol(CO2) => kg(C)
   real, parameter :: kgom2_2_tonoha = 10.         ! kg(C)/m² => ton(C)/ha
   real, parameter :: tonoha_2_kgom2 = 0.1         ! ton(C)/ha => kg(C)/m²

   !---------------------------------------------------------------------------------------!
   ! Molar masses and derived variables                                                    !
   !---------------------------------------------------------------------------------------!
   real, parameter :: mmdry       = 0.02897        ! Mean dry air molar mass    [   kg/mol]
   real, parameter :: mmvap       = 0.01801505     ! Mean water molar mass      [   kg/mol]
   real, parameter :: mmco2       = 0.0440095      ! Mean CO2 molar mass        [   kg/mol]
   real, parameter :: mmdov       = mmdry/mmvap    ! mmdry/mmvap                [     ----]
   real, parameter :: mmvod       = mmvap/mmdry    ! mmvap/mmdry                [     ----]
   real, parameter :: mmdry1000   = 1000.*mmdry    ! Mean dry air molar mass    [   kg/mol]
   real, parameter :: mmdryi      = 1./mmdry       ! 1./mmdry                   [   mol/kg]
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
   real(kind=8), parameter :: gorvap8      = dble(gorvap    )
   real(kind=8), parameter :: grav8        = dble(grav      )
   real(kind=8), parameter :: hr_sec8      = dble(hr_sec    )
   real(kind=8), parameter :: idns8        = dble(idns      )
   real(kind=8), parameter :: mmdry8       = dble(mmdry     )
   real(kind=8), parameter :: mmdryi8      = dble(mmdryi    )
   real(kind=8), parameter :: pi18         = dble(pi1       )
   real(kind=8), parameter :: pio1808      = dble(pio180    )
   real(kind=8), parameter :: qicet38      = dble(qicet3    )
   real(kind=8), parameter :: qliqt38      = dble(qliqt3    )
   real(kind=8), parameter :: stefan8      = dble(stefan    )
   real(kind=8), parameter :: sqrt2o28     = dble(sqrt2o2   )
   real(kind=8), parameter :: t3ple8       = dble(t3ple     )
   real(kind=8), parameter :: tsupercool8  = dble(tsupercool)
   real(kind=8), parameter :: twopi8       = dble(twopi     )
   real(kind=8), parameter :: twothirds8   = dble(twothirds )
   real(kind=8), parameter :: umol_2_kgC8  = dble(umol_2_kgC)
   real(kind=8), parameter :: vonk8        = dble(vonk      )
   real(kind=8), parameter :: wdns8        = dble(wdns      )
   real(kind=8), parameter :: wdnsi8       = dble(wdnsi     )
   !---------------------------------------------------------------------------------------!



end Module consts_coms
