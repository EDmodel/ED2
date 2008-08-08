Module consts_coms

! Making the constants to match when running a coupled run.

#if defined(COUPLED)
   use rconstants, only:                                                                   &
        b_pi1      => pi1     , b_twopi    => twopi   , b_pio180   => pio180               &
       ,b_pi4      => pi4     , b_pio4     => pio4    , b_vonk     => vonk                 &
       ,b_stefan   => stefan  , b_yr_day   => yr_day  , b_day_sec  => day_sec              &
       ,b_day_hr   => day_hr  , b_hr_sec   => hr_sec  , b_min_sec  => min_sec              &
       ,b_g        => g       , b_erad     => erad    , b_p00      => p00                  &
       ,b_p00i     => p00i    , b_rgas     => rgas    , b_cp       => cp                   &
       ,b_cpi      => cpi     , b_cpor     => cpor    , b_rocp     => rocp                 &
       ,b_cpog     => cpog    , b_rm       => rm      , b_ep       => ep                   &
       ,b_cliq     => cliq    , b_cliq1000 => cliq1000, b_cliqi    => cliqi                &
       ,b_cice     => cice    , b_cice1000 => cice1000, b_cicei    => cicei                &
       ,b_t00      => t00     , b_t3ple    => t3ple   , b_alvl     => alvl                 &
       ,b_alli     => alli    , b_alvi     => alvi    , b_alli1000 => alli1000             &
       ,b_allii    => allii

   implicit none

   real, parameter :: pi1     = b_pi1      , twopi   = b_twopi    , pio180   = b_pio180
   real, parameter :: pi4     = b_pi4      , pio4    = b_pio4     , vonk     = b_vonk
   real, parameter :: stefan  = b_stefan   , yr_day  = b_yr_day   , day_sec  = b_day_sec
   real, parameter :: day_hr  = b_day_hr   , hr_sec  = b_hr_sec   , min_sec  = b_min_sec
   real, parameter :: grav    = b_g        , erad    = b_erad     , p00      = b_p00
   real, parameter :: p00i    = b_p00i     , rdry    = b_rgas     , cp       = b_cp
   real, parameter :: cpi     = b_cpi      , cpor    = b_cpor     , rocp     = b_rocp
   real, parameter :: cpog    = b_cpog     , rvap    = b_rm       , ep       = b_ep
   real, parameter :: cliq    = b_cliq     , cliq1000= b_cliq1000 , cliqi    = b_cliqi
   real, parameter :: cice    = b_cice     , cice1000= b_cice1000 , cicei    = b_cicei
   real, parameter :: t00     = b_t00      , t3ple   = b_t3ple    , alvl     = b_alvl
   real, parameter :: alli    = b_alli     , alvi    = b_alvi     , alli1000 = b_alli1000
   real, parameter :: allii   = b_allii    

#else
   implicit none
!--------------------------------------------------------------------------------------------------!
! Geometric constants                                                                              !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   pi1         = 3.14159265358979     ! Pi
   real, parameter ::   twopi       = pi1* 2.              ! 2 Pi
   real, parameter ::   pio180      = pi1/ 180.            ! Pi/180
   real, parameter ::   pi4         = pi1 * 4.             ! 4 Pi
   real, parameter ::   pio4        = pi1 /4.              ! Pi/4

!--------------------------------------------------------------------------------------------------!
! General constants                                                                                !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   vonk        = 0.40                 ! Von Kármán constant                ( k)
   real, parameter ::   stefan      = 5.6696e-8            ! Stefan-Boltzmann constant       (sigma)

!--------------------------------------------------------------------------------------------------!
! Time conversion units                                                                            !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   yr_day      = 365.2425             ! # of days in a year   
   real, parameter ::   day_sec     = 86400.               ! # of seconds in a day
   real, parameter ::   day_hr      = 24.                  ! # of hours in a day
   real, parameter ::   hr_sec      = 3600.                ! # of seconds in an hour
   real, parameter ::   min_sec     = 60.                  ! # of seconds in a minute

!--------------------------------------------------------------------------------------------------!
! Earth properties                                                                                 !
!--------------------------------------------------------------------------------------------------!
   real, parameter :: grav     = 9.81           ! acceleration of gravity [m/s^2]
   real, parameter :: erad     = 6370997.       ! Earth radius [m]
   real, parameter :: p00      = 1.e5           ! 1000 hPa                           (p0)
   real, parameter :: p00i     = 1. / p00       ! 1/p0

!--------------------------------------------------------------------------------------------------!
! Dry air properties                                                                               !
!--------------------------------------------------------------------------------------------------!
   real, parameter :: rdry     = 287.04         ! dry air gas constant [J/(kg K)]
   real, parameter :: cp       = 1004.          ! dry air spec heat at const P 
   real, parameter :: cpi      = 1. / cp
   real, parameter :: cpor     = cp / rdry
   real, parameter :: rocp     = rdry / cp           
   real, parameter :: cpog = cp / grav

!--------------------------------------------------------------------------------------------------!
! Water vapour properties                                                                          !
!--------------------------------------------------------------------------------------------------!
   real, parameter :: rvap     = 461.5          ! water vapor gas constant [J/(kg K)]
   real, parameter :: ep = rdry / rvap

!--------------------------------------------------------------------------------------------------!
! Liquid water properties                                                                          !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   cliq        = 4.186e3              ! Specific heat of liquid water (1 kcal)
   real, parameter ::   cliq1000    = 1000. * cliq         ! 1000*Specific heat of liquid water
   real, parameter ::   cliqi       = 1. / cliq            ! Inverse of water heat capacity


!--------------------------------------------------------------------------------------------------!
! Ice properties                                                                                   !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   cice        = 2.106e3              ! Heat capacity*ice density (J/ kg / K)
   real, parameter ::   cice1000    = 1000. * cice         ! Heat capacity*ice density (J / m³ /K)
   real, parameter ::   cicei       = 1. / cice            ! Inverse of ice heat capacity


!--------------------------------------------------------------------------------------------------!
! Phase change properties                                                                          !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   t00         = 273.15               ! 0°C                                (T0)
   real, parameter ::   t3ple       = 273.16               ! Water triple point
   real, parameter ::   alvl        = 2.50e6               ! Latent heat of vaporization        (Lv)
   real, parameter ::   alli        = 3.34e5               ! Latent heat of melting             (Lf)
   real, parameter ::   alvi        = alvl+alli            ! Latent heat of sublimation         (Ls)
   real, parameter ::   alli1000    = 1000. * alli         ! Latent heat of melting             (Lf)
   real, parameter ::   allii       = 1. / alli   

#endif

! This constants are ED specific, need to be defined locally even in a coupled run.

   real, parameter :: umol_2_kgC      = 1.20107e-8 ! µmol(CO2) => kg(C)
   real, parameter :: kgom2_2_tonoha = 10.         ! kg(C)/m² => ton(C)/ha
   real, parameter :: tonoha_2_kgom2 = 0.1         ! ton(C)/ha => kg(C)/m²

end Module consts_coms
