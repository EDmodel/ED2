!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module rconstants

!--------------------------------------------------------------------------------------------------!
! Geometric constants                                                                              !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   pi1         = 3.14159265358979     ! Pi
   real, parameter ::   pii         = 1./pi1               ! 1/Pi
   real, parameter ::   halfpi      = pi1/2                ! Pi/2
   real, parameter ::   twopi       = pi1* 2.              ! 2 Pi
   real, parameter ::   pio180      = pi1/ 180.            ! Pi/180
   real, parameter ::   onerad      = 180. / pi1           ! 180./pi
   real, parameter ::   pi4         = pi1 * 4.             ! 4 Pi
   real, parameter ::   pio4        = pi1 /4.              ! Pi/4
   real, parameter ::   pio6        = pi1 /6.              ! Pi/6
   real, parameter ::   pio6i       = 6.  /pi1             ! 6/Pi

!--------------------------------------------------------------------------------------------------!
! General constants                                                                                !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   vonk        = 0.40                 ! Von Kármán constant                ( k)
   real, parameter ::   stefan      = 5.6696e-8            ! Stefan-Boltzmann constant       (sigma)
   real, parameter ::   srtwo       = 1.414213562373095    ! Square root of 2.
   real, parameter ::   srthree     = 1.732050807568877    ! Square root of 3.
   real, parameter ::   srtwoi      = 1./srtwo             ! 1./ Square root of 2.
   real, parameter ::   srthreei    = 1./srthree           ! 1./ Square root of 3.

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
   real, parameter ::   g           = 9.80                 ! Gravity acceleration               ( g)
   real, parameter ::   gg          = .5 * g               ! ½ g
   real, parameter ::   spcon       = 111120.              ! One degree of latitude              [m]
   real, parameter ::   erad        = 6370997.             ! Earth radius                        [m]
   real, parameter ::   eradi       = 1./erad              ! Inverse of Earth radius           [1/m]
   real, parameter ::   erad2       = erad*2               ! Earth diameter [m]
   real, parameter ::   ss60        = 1.8663               ! Polar stereographic conversion to 60°
   real, parameter ::   omega       = 7.292e-5             ! Earth's rotation speed          (OMEGA)
   real, parameter ::   viscos      = .15e-4               ! Viscosity coefficient
   real, parameter ::   solar       = 1.3533e3             ! Solar constan                      (S0)
   real, parameter ::   p00         = 1.e5                 ! 1000 hPa                           (p0)
   real, parameter ::   p00i        = 1. / p00             ! 1/p0

!--------------------------------------------------------------------------------------------------!
! Dry air properties                                                                               !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   rgas        = 287.04               ! Gas constant for dry air           (Ra)
   real, parameter ::   cp          = 1004.                ! Specific heat at constant pressure (Cp)
   real, parameter ::   cv          = 717.                 ! Specific heat at constant volume   (Cv)
   real, parameter ::   cpocv       = cp / cv              ! Cp/Cv
   real, parameter ::   cvocp       = cv / cp              ! Cp/Cv
   real, parameter ::   cpog        = cp /g                ! cp/g
   real, parameter ::   rocp        = rgas / cp            ! Ra/cp 
   real, parameter ::   cpor        = cp / rgas            ! Cp/Ra
   real, parameter ::   rocv        = rgas / cv            ! Ra/Cv
   real, parameter ::   gocp        = g / cp               ! g/Cp
   real, parameter ::   gordry      = g / rgas             ! g/Ra
   real, parameter ::   cpi         = 1. / cp              ! 1/Cp
   real, parameter ::   cpi4        = 4. * cpi             ! 4/Cp
   real, parameter ::   p00k        = 26.870941            ! p0 ** (Ra/Cp)  
   real, parameter ::   p00ki       = 1. / p00k            ! p0 ** (-Ra/Cp)

!--------------------------------------------------------------------------------------------------!
! Water vapour properties                                                                          !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   rm          = 461.5                ! Gas constant for water vapour      (Rv)
   real, parameter ::   gorm        = g / rm               ! g/Ra
   real, parameter ::   ep          = rgas / rm            ! Ra/Rv
   real, parameter ::   epi         = rm / rgas            ! Rv/Ra
   real, parameter ::   eps_virt    = (rm - rgas) / rm     ! (Rv-Ra)/Rv, that 0.608 at Tv
   real, parameter ::   rmocp       = rm / cp              ! Rv/cp
   real, parameter ::   toodry      = 1.e-8                ! Minimum acceptable mixing ratio.

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
   real, parameter ::   alvl2       = 6.25e12              ! Lv²
   real, parameter ::   alvi2       = 8.032e12             ! Ls²
   real, parameter ::   alvli       = 1. / alvl            ! 1/Lv
   real, parameter ::   allii       = 1. / alli            ! 1/Lfil
   real, parameter ::   aklv        = alvl / cp            ! Lv/Cp
   real, parameter ::   akiv        = alvi / cp            ! Ls/Cp
   real, parameter ::   lvordry     = alvl / rgas          ! Lv/Ra
   real, parameter ::   lvorvap     = alvl / rm            ! Lv/Rv

!--------------------------------------------------------------------------------------------------!
!    Minimum temperature below which assuming the latent heats as constants becomes really bad.    !
!   See :                                                                                          !
!    Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential temperature as  !
!        a thermodynamic variable in deep atmospheric models. Mon. Wea. Rev., v. 109, 1094-1102.   !                                                                 !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   ttripoli    = 253.                                                          
   real, parameter ::   cp253i     = cpi / ttripoli       ! 1/(253*Cp)

!--------------------------------------------------------------------------------------------------!
! Lower bounds for turbulence-related variables                                                     
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   tkmin       = 5.e-4                ! Minimum value of TKE
   real, parameter ::   sigwmin     = 1.e-4                ! Minimum value of sigma-w
   real, parameter ::   abslmomin   = 1.e-4                ! Minimum abs value of Obukhov length
   real, parameter ::   ltscalemax  = 1.e5                 ! Maximum value of Lagrangian timescale
   real, parameter ::   abswltlmin  = 1.e-4                ! Minimum abs value of Theta*
   real, parameter ::   lturbmin    = 1.e-3                ! Minimum abs value of turb. lenght

                                                           !    ing for bisection.
!--------------------------------------------------------------------------------------------------!
! Don't know what do these variables stand for.                                                    !
!--------------------------------------------------------------------------------------------------!
   real, parameter ::   cww         = 4218.                !
   real, parameter ::   c0          = 752.55 * 4.18684e4   !
   real, parameter ::   rowt        = 1.e3                 !


end Module
