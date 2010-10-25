!==========================================================================================!
!. File: therm_lib8.f90                                                                    !
!                                                                                          !
!  Based on BRAMS-4.0.6   This file contains most functions and subroutines that deal with !
! several thermodynamic conversions that are needed in double precision.  Most of them     !
! have the equivalent in single precision in therm_lib.  These procedures were built to    !
! avoid assumptions like hydrostatic and linearisation.  Most equations could not be       !
! solved analytically, and the standard here was to use Newton's method as the default,    !
! always having bisection or, more often, the modified Regula Falsi (Illinois) method in   !
! case Newton's fails.                                                                     !
!==========================================================================================!
!==========================================================================================!
module therm_lib8
   use therm_lib, only : toler4     => toler     & ! intent(in)
                       , maxfpo4    => maxfpo    & ! intent(in)
                       , maxit4     => maxit     & ! intent(in)
                       , maxlev4    => maxlev    & ! intent(in)
                       , newthermo4 => newthermo & ! intent(in)
                       , level4     => level     & ! intent(in)
                       , vapour_on4 => vapour_on & ! intent(in)
                       , cloud_on4  => cloud_on  & ! intent(in)
                       , bulk_on4   => bulk_on   ! ! intent(in)

   !---------------------------------------------------------------------------------------!
   !     Relative tolerance for iterative methods. The smaller the value, the more         !
   ! accurate the result, but it will slow down the run.  Notice that we are using the     !
   ! tolerance that is based on the single precision...                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter ::   toler8  = dble(toler4)


   integer, parameter ::   maxfpo = maxfpo4         ! Maximum # of iterations before crash-
                                                    !   ing for false position method.

   integer, parameter ::   maxit  = maxit4          ! Maximum # of iterations before crash-
                                                    !   ing, for other methods.

   integer, parameter ::   maxlev = maxlev4         ! Maximum # of levels for adaptive     
                                                    !   quadrature methods.    

   logical, parameter ::   newthermo = newthermo4   ! Use new thermodynamics [T|F]

   !---------------------------------------------------------------------------------------!
   !   This is the "level" variable, that used to be in micphys. Since it affects more the !
   ! thermodynamics choices than the microphysics, it was moved to here.                   !
   !---------------------------------------------------------------------------------------!
   integer, parameter ::   level = level4

   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical, parameter ::   vapour_on = vapour_on4
   logical, parameter ::   cloud_on  = cloud_on4
   logical, parameter ::   bulk_on   = bulk_on4
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These constants came from the paper in which the saturation vapour pressure is    !
   ! based on:                                                                             !
   !                                                                                       !
   !  Murphy, D. M.; Koop, T., 2005: Review of the vapour pressures of ice and supercooled !
   !     water for atmospheric applications. Q. J. Royal Meteor. Soc., vol. 31, pp. 1539-  !
   !     1565 (hereafter MK05).                                                            !
   !                                                                                       !
   !  These equations give the triple point at t3ple, with vapour pressure being es3ple.   !
   !---------------------------------------------------------------------------------------!
   !----- Coefficients based on equation (7): ---------------------------------------------!
   real(kind=8), dimension(0:3), parameter :: iii_78 = (/ 9.550426d0, -5.723265d3          &
                                                        , 3.53068d0,  -7.28332d-3 /)
   !----- Coefficients based on equation (10), first fit ----------------------------------!
   real(kind=8), dimension(0:3), parameter :: l01_108 = (/ 5.4842763d1,-6.76322d3          &
                                                         ,-4.210d0    , 3.67d-4 /)
   !----- Coefficients based on equation (10), second fit ---------------------------------!
   real(kind=8), dimension(0:3), parameter :: l02_108 = (/ 5.3878d1   ,-1.33122d3          &
                                                         ,-9.44523d0  , 1.4025d-2 /)
   !----- Coefficients based on the hyperbolic tangent ------------------------------------!
   real(kind=8), dimension(2)  , parameter :: ttt_108 = (/4.15d-2     , 2.188d2 /)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     These constants came from the paper in which the saturation vapour pressure is    !
   ! based on:                                                                             !
   !                                                                                       !
   !  Flatau, P. J.; Walko, R. L.; Cotton, W. R., 1992: Polynomial fits to saturation      !
   !     vapor pressure. J. Appl. Meteor., vol. 31, pp. 1507-1513. (hereafter FWC92).      !
   !                                                                                       !
   !  These equations give the triple point at 273.004K.                                   !
   !  N.B.: The coefficients here don't seem to match those listed on FWC92, but that's    !
   !        what was on the original code...                                               !
   !---------------------------------------------------------------------------------------!
   !----- Coefficients for esat (liquid) --------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: cll8 =   &
        (/ .6105851d+03,  .4440316d+02,  .1430341d+01  &
        , .2641412d-01,  .2995057d-03,  .2031998d-05   &
        , .6936113d-08,  .2564861d-11, -.3704404d-13 /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: cii8 =   &
        (/ .6114327d+03,  .5027041d+02,  .1875982d+01  &
        , .4158303d-01,  .5992408d-03,  .5743775d-05   &
        , .3566847d-07,  .1306802d-09,  .2152144d-12 /)
   !----- Coefficients for d(esat)/dT (liquid) --------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: dll8 =   &
        (/ .4443216d+02,  .2861503d+01,  .7943347d-01  &
        , .1209650d-02,  .1036937d-04,  .4058663d-07   &
        ,-.5805342d-10, -.1159088d-11, -.3189651d-14 /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: dii8 =   &
        (/ .5036342d+02,  .3775758d+01,  .1269736d+00  &
        , .2503052d-02,  .3163761d-04,  .2623881d-06   &
        , .1392546d-08,  .4315126d-11,  .5961476d-14 /)
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!


   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour pressure as a function of   !
   ! Kelvin temperature. This expression came from MK05, equation (10).                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function eslf8(temp,l1funout,l2funout,ttfunout)
      use consts_coms, only : t008
      implicit none
      real(kind=8), intent(in)            :: temp
      real(kind=8), intent(out), optional :: l1funout,ttfunout,l2funout
      real(kind=8)                        :: l1fun,ttfun,l2fun,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         l1fun  = l01_108(0) + l01_108(1)/temp + l01_108(2)*log(temp) + l01_108(3) * temp
         l2fun  = l02_108(0) + l02_108(1)/temp + l02_108(2)*log(temp) + l02_108(3) * temp
         ttfun  = tanh(ttt_108(1) * (temp - ttt_108(2)))
         eslf8  = exp(l1fun + ttfun*l2fun)

         if (present(l1funout)) l1funout = l1fun
         if (present(l2funout)) l2funout = l2fun
         if (present(ttfunout)) ttfunout = ttfun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-8.d1,temp-t008)
         eslf8 = cll8(0) + x * (cll8(1) + x * (cll8(2) + x * (cll8(3) + x * (cll8(4)       &
                         + x * (cll8(5) + x * (cll8(6) + x * (cll8(7) + x * cll8(8)) ))))))

         if (present(l1funout)) l1funout = eslf8
         if (present(l2funout)) l2funout = eslf8
         if (present(ttfunout)) ttfunout = eslf8
      end if

      return
   end function eslf8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour pressure as a function of      !
   ! Kelvin temperature, based on MK05 equation (7).                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function esif8(temp,iifunout)
      use consts_coms, only : t008
      implicit none
      real(kind=8), intent(in)            :: temp
      real(kind=8), intent(out), optional :: iifunout
      real(kind=8)                        :: iifun,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         iifun  = iii_78(0) + iii_78(1)/temp + iii_78(2) * log(temp) + iii_78(3) * temp
         esif8  = exp(iifun)
      
         if (present(iifunout)) iifunout=iifun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-8.d1,temp-t008)
         esif8 = cii8(0) + x * (cii8(1) + x * (cii8(2) + x * (cii8(3) + x * (cii8(4)       &
                         + x * (cii8(5) + x * (cii8(6) + x * (cii8(7) + x * cii8(8)) ))))))

         if (present(iifunout)) iifunout=esif8
      end if
      return
   end function esif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation vapour pressure as a function of  Kelvin  !
   ! temperature. It chooses which phase to look depending on whether the temperature is   !
   ! below or above the triple point.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function eslif8(temp,useice)
      use consts_coms, only: t3ple8
      implicit none
      real(kind=8), intent(in)           :: temp
      logical     , intent(in), optional :: useice
      logical                            :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple8
      else 
         brrr_cold = bulk_on .and. temp < t3ple8
      end if

      if (brrr_cold) then
         eslif8 = esif8(temp) ! Ice saturation vapour pressure 
      else
         eslif8 = eslf8(temp) ! Liquid saturation vapour pressure
      end if

      return
   end function eslif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour mixing ratio as a function  !
   ! of pressure and Kelvin temperature.                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rslf8(pres,temp)
      use consts_coms, only : ep8,toodry8
      implicit none
      real(kind=8), intent(in) :: pres,temp
      real(kind=8)             :: esl

      esl   = eslf8(temp)
      rslf8 = max(toodry8,ep8*esl/(pres-esl))

      return
   end function rslf8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour mixing ratio as a function of  !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rsif8(pres,temp)
      use consts_coms, only : ep8,toodry8
      implicit none
      real(kind=8), intent(in) :: pres,temp
      real(kind=8)             :: esi

      esi   = esif8(temp)
      rsif8 = max(toodry8,ep8*esi/(pres-esi))

      return
   end function rsif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation vapour mixing ratio, over liquid or ice   !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rslif8(pres,temp,useice)
      use consts_coms, only: t3ple8,ep8
      implicit none
      real(kind=8), intent(in)           :: pres,temp
      logical     , intent(in), optional :: useice
      real(kind=8)                       :: esz
      logical                            :: brrr_cold

      !----- Checking which saturation (liquid or ice) I should use here ------------------!
      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple8
      else 
         brrr_cold = bulk_on .and. temp < t3ple8
      end if
    
      !----- Finding the saturation vapour pressure ---------------------------------------!
      if (brrr_cold) then
         esz = esif8(temp)
      else
         esz = eslf8(temp)
      end if

      rslif8 = ep8 * esz / (pres - esz)

      return
   end function rslif8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-liquid equilibrium density for vapour, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovsl8(temp)
      use consts_coms, only : rh2o8
      implicit none
      real(kind=8), intent(in) :: temp
      real(kind=8)             :: eequ
      eequ    = eslf8(temp)
      rhovsl8 = eequ / (rh2o8 * temp)
      return
   end function rhovsl8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-ice equilibrium density for vapour, as a      !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovsi8(temp)
      use consts_coms, only : rh2o8
      implicit none
      real(kind=8), intent(in) :: temp
      real(kind=8)             :: eequ
      eequ   = esif8(temp)
      rhovsi8 = eequ / (rh2o8 * temp)
      return
   end function rhovsi8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation density for vapour, as a function of tem- !
   ! perature in Kelvin. It will decide between ice-vapour or liquid-vapour based on the   !
   ! temperature.                                                                          !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovsil8(temp,useice)
      use consts_coms, only : rh2o8
      implicit none
      real(kind=8), intent(in)           :: temp
      logical     , intent(in), optional :: useice
      real(kind=8)                       :: eequ

      if (present(useice)) then
         eequ = eslif8(temp,useice)
      else
         eequ = eslif8(temp)
      end if

      rhovsil8 = eequ / (rh2o8 * temp)

      return
   end function rhovsil8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour       !
   ! pressure with respect to temperature as a function of Kelvin temperature.             !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function eslfp8(temp)
      use consts_coms, only: t008
      implicit none
      real(kind=8), intent(in) :: temp
      real(kind=8)             :: esl,l2fun,ttfun,l1prime,l2prime,ttprime,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esl     = eslf8(temp,l2funout=l2fun,ttfunout=ttfun)
         l1prime = -l01_108(1)/(temp*temp) + l01_108(2)/temp + l01_108(3)
         l2prime = -l02_108(1)/(temp*temp) + l02_108(2)/temp + l02_108(3)
         ttprime =  ttt_108(1)*(1.-ttfun*ttfun)
         eslfp8  = esl * (l1prime + l2prime*ttfun + l2fun*ttprime)
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-8.d1,temp-t008)
         eslfp8 = dll8(0) + x * (dll8(1) + x * (dll8(2) + x * (dll8(3) + x * (dll8(4)      &
                          + x * (dll8(5) + x * (dll8(6) + x * (dll8(7) + x * dll8(8)) ))))))
      end if


      return
   end function eslfp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of ice saturation vapour pressure !
   ! with respect to temperature as a function of Kelvin temperature.                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function esifp8(temp)
      use consts_coms, only: t008
      implicit none
      real(kind=8), intent(in) :: temp
      real(kind=8)             :: esi,iiprime,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esi      = esif8(temp)
         iiprime  = -iii_78(1)/(temp*temp) + iii_78(2)/temp + iii_78(3)
         esifp8   = esi * iiprime
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-8.d1,temp-t008)
         esifp8 = dii8(0) + x * (dii8(1) + x * (dii8(2) + x * (dii8(3) + x * (dii8(4)      &
                          + x * (dii8(5) + x * (dii8(6) + x * (dii8(7) + x * dii8(8)) ))))))
      end if

      return
   end function esifp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of saturation vapour pressure as  !
   ! a function of  Kelvin temperature. It chooses which phase to look depending on        !
   ! whether the temperature is below or above the triple point.                           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function eslifp8(temp,useice)
      use consts_coms, only: t3ple8
      implicit none
      real(kind=8), intent(in)           :: temp
      logical     , intent(in), optional :: useice
      logical                            :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple8
      else 
         brrr_cold = bulk_on .and. temp < t3ple8
      end if

      if (brrr_cold) then
         eslifp8 = esifp8(temp) ! d(Ice saturation vapour pressure)/dT
      else
         eslifp8 = eslfp8(temp) ! d(Liquid saturation vapour pressure)/dT
      end if

      return
   end function eslifp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour mix-  !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rslfp8(pres,temp)
      use consts_coms, only: ep8
      implicit none
      real(kind=8), intent(in) :: pres,temp
      real(kind=8)             :: desdt,esl,pdry

      esl     = eslf8(temp)
      desdt   = eslfp8(temp)
      
      pdry    = pres-esl
      rslfp8  = ep8 * pres * desdt / (pdry*pdry)

      return
   end function rslfp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of ice saturation vapour mix-     !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rsifp8(pres,temp)
      use consts_coms, only: ep8
      implicit none
      real(kind=8), intent(in) :: pres,temp
      real(kind=8)             :: desdt,esi,pdry

      esi     = esif8(temp)
      desdt   = esifp8(temp)
      
      pdry    = pres-esi
      rsifp8  = ep8 * pres * desdt / (pdry*pdry)

      return
   end function rsifp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour mix-  !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rslifp8(pres,temp,useice)
      use consts_coms, only: t3ple8
      implicit none
      real(kind=8), intent(in)           :: pres,temp
      logical     , intent(in), optional :: useice
      logical                            :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple8
      else 
         brrr_cold = bulk_on .and. temp < t3ple8
      end if

      if (brrr_cold) then
         rslifp8=rsifp8(pres,temp)
      else
         rslifp8=rslfp8(pres,temp)
      end if

      return
   end function rslifp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of vapour-liquid equilibrium density, as  !
   ! a function of temperature in Kelvin.                                                  !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovslp8(temp)
      use consts_coms, only : rh2o8
      implicit none
      real(kind=8), intent(in) :: temp
      real(kind=8)             :: es,desdt

      es       = eslf8(temp)
      desdt    = eslfp8(temp)
      rhovslp8 = (desdt-es/temp) / (rh2o8 * temp)

      return
   end function rhovslp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of vapour-ice equilibrium density, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovsip8(temp)
      use consts_coms, only : rh2o8
      implicit none
      real(kind=8), intent(in) :: temp
      real(kind=8)             :: es,desdt

      es       = esif8(temp)
      desdt    = esifp8(temp)
      rhovsip8 = (desdt-es/temp) / (rh2o8 * temp)

      return
   end function rhovsip8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of saturation density for vapour, as a    !
   ! function of temperature in Kelvin. It will decide between ice-vapour or liquid-vapour !
   ! based on the temperature.                                                             !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovsilp8(temp,useice)
      use consts_coms, only: t3ple8
      implicit none
      real(kind=8), intent(in)           :: temp
      logical     , intent(in), optional :: useice
      logical                            :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple8
      else 
         brrr_cold = bulk_on .and. temp < t3ple8
      end if

      if (brrr_cold) then
         rhovsilp8 = rhovsip8(temp)
      else
         rhovsilp8 = rhovslp8(temp)
      end if

      return
   end function rhovsilp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the temperature from the liquid sat. vapour pressure.    !
   ! This is truly the inverse of eslf, which is done iteratively since it's not a simple  !
   ! function to invert. It uses Newton's method, which should take care of most cases. In !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function tslf8(pvap)

       implicit none
      !----- Argument ---------------------------------------------------------------------!
      real(kind=8), intent(in) :: pvap      ! Saturation vapour pressure           [    Pa]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)             :: deriv     ! Function derivative                  [    Pa]
      real(kind=8)             :: fun       ! Function for which we seek a root.   [    Pa]
      real(kind=8)             :: funa      ! Smallest  guess function             [    Pa]
      real(kind=8)             :: funz      ! Largest   guess function             [    Pa]
      real(kind=8)             :: tempa     ! Smallest guess (or previous guess)   [    Pa]
      real(kind=8)             :: tempz     ! Largest guess (or new guess )        [    Pa]
      real(kind=8)             :: delta     ! Aux. var --- 2nd guess for bisection [      ]
      integer                  :: itn,itb   ! Iteration counter                    [   ---]
      logical                  :: converged ! Convergence handle                   [   ---]
      logical                  :: zside     ! Flag to check for 1-sided approach.  [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, using Bolton (1980) equation 11, giving es in Pa and T in K -----!
      tempa = (2.965d1 * log(pvap) - 5.01678d3)/(log(pvap)-2.40854d1)
      funa  = eslf8(tempa) - pvap
      deriv = eslfp8(tempa)
      !----- Copying just in case it fails at the first iteration -------------------------!
      tempz = tempa
      fun   = funa
      
      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, go with bisection ----!
         !----- Copying the previous guess ------------------------------------------------!
         tempa = tempz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         tempz = tempa - fun/deriv
         fun   = eslf8(tempz) - pvap
         deriv = eslfp8(tempz)
         
         converged = abs(tempa-tempz) < toler8 * tempz
         if (converged) then
            tslf8 = 5.d-1 * (tempa+tempz)
            return
         elseif (fun == 0.d0) then !Converged by luck!
            tslf8 = tempz
            return
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < 1.d2 * toler8 * tempa) then
            delta = 1.d2 * toler8 * tempa
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),1.d2 * toler8 * tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
            funz  = eslf8(tempz) - pvap
            zside = funa*funz < 0.d0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                          ,'tslf8','therm_lib8.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tslf8 =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tslf8-tempa) < toler8 * tslf8
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         fun       =  eslf8(tslf8) - pvap

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0.d0 ) then
            tempz = tslf8
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 5.d-1
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = tslf8
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz=funz * 5.d-1
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         call fatal_error('Temperature didn''t converge, giving up!!!'                     &
                       ,'tslf8','therm_lib8.f90')
      end if
      
      return
   end function tslf8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the temperature from the ice saturation vapour pressure. !
   ! This is truly the inverse of esif, which is done iteratively since it's not a simple  !
   ! function to invert. It uses Newton's method, which should take care of most cases. In !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function tsif8(pvap)

      implicit none
      !----- Argument ---------------------------------------------------------------------!
      real(kind=8), intent(in) :: pvap      ! Saturation vapour pressure           [    Pa]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)             :: deriv     ! Function derivative                  [    Pa]
      real(kind=8)             :: fun       ! Function for which we seek a root.   [    Pa]
      real(kind=8)             :: funa      ! Smallest  guess function             [    Pa]
      real(kind=8)             :: funz      ! Largest   guess function             [    Pa]
      real(kind=8)             :: tempa     ! Smallest guess (or previous guess)   [    Pa]
      real(kind=8)             :: tempz     ! Largest   guess (or new guess)       [    Pa]
      real(kind=8)             :: delta     ! Aux. var --- 2nd guess for bisection [      ]
      integer                  :: itn,itb   ! Iteration counter                    [   ---]
      logical                  :: converged ! Convergence handle                   [   ---]
      logical                  :: zside     ! Flag to check for one-sided approach [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, using Murphy-Koop (2005), equation 8. ---------------------------!
      tempa = (1.814625d0 * log(pvap) +6.190134d3)/(2.9120d1 - log(pvap))
      funa  = esif8(tempa) - pvap
      deriv = esifp8(tempa)
      !----- Copying just in case it fails at the first iteration -------------------------!
      tempz = tempa
      fun   = funa
      
      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, go with bisection ----!
         !----- Copying the previous guess ------------------------------------------------!
         tempa = tempz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         tempz = tempa - fun/deriv
         fun   = esif8(tempz) - pvap
         deriv = esifp8(tempz)
         
         converged = abs(tempa-tempz) < toler8 * tempz
         if (converged) then
            tsif8 = 5.d-1*(tempa+tempz)
            return
         elseif (fun == 0.d0) then
            tsif8 = tempz
            return
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < 1.d2 * toler8 * tempa) then
            delta = 1.d2 * toler8 * delta
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),1.d2 * toler8 * tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
            funz  = esif8(tempz) - pvap
            zside = funa*funz < 0.d0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                          ,'tsif8','therm_lib8.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tsif8 =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tsif8-tempa) < toler8 * tsif8
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         fun       =  esif8(tsif8) - pvap

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0. ) then
            tempz = tsif8
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 5.d-1
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = tsif8
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz=funz * 5.d-1
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         call fatal_error('Temperature didn''t converge, giving up!!!'                     &
                       ,'tsif8','therm_lib8.f90')
      end if
      
      return
   end function tsif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the temperature from the ice or liquid mixing ratio.     !
   ! This is truly the inverse of eslf and esif.                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function tslif8(pvap,useice)
      use consts_coms, only: es3ple8

      implicit none
      real(kind=8), intent(in)           :: pvap
      logical     , intent(in), optional :: useice
      logical                            :: brrr_cold
      
      !------------------------------------------------------------------------------------!
      !    Since pvap is a function of temperature only, we can check the triple point     !
      ! from the saturation at the triple point, like what we would do for temperature.    !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         brrr_cold = useice  .and. pvap < es3ple8
      else 
         brrr_cold = bulk_on .and. pvap < es3ple8
      end if


      if (brrr_cold) then
         tslif8 = tsif8(pvap)
      else
         tslif8 = tslf8(pvap)
      end if

      return
   end function tslif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the dew point temperature given the pressure and vapour    !
   ! mixing ratio. THIS IS DEWPOINT ONLY, WHICH MEANS THAT IT WILL IGNORE ICE EFFECT. For  !
   ! a full, triple-point dependent routine use DEWFROSTPOINT                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dewpoint8(pres,rsat)
      use consts_coms, only: ep8,toodry8

      implicit none
      real(kind=8), intent(in) :: pres, rsat
      real(kind=8)             :: rsatoff, pvsat
      
      rsatoff   = max(toodry8,rsat)
      pvsat     = pres*rsatoff / (ep8 + rsatoff)
      dewpoint8 = tslf8(pvsat)

      return
   end function dewpoint8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the frost point temperature given the pressure and vapour  !
   ! mixing ratio. THIS IS FROSTPOINT ONLY, WHICH MEANS THAT IT WILL IGNORE LIQUID EFFECT. !
   ! For a full, triple-point dependent routine use DEWFROSTPOINT                          !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function frostpoint8(pres,rsat)
      use consts_coms, only: ep8,toodry8

      implicit none
      real(kind=8), intent(in) :: pres, rsat
      real(kind=8)             :: rsatoff, pvsat
      
      rsatoff     = max(toodry8,rsat)
      pvsat       = pres*rsatoff / (ep8 + rsatoff)
      frostpoint8 = tsif8(pvsat)

      return
   end function frostpoint8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the saturation point temperature given the pressure and    !
   ! vapour mixing ratio. This will check whether the vapour pressure is above or below    !
   ! the triple point vapour pressure, finding dewpoint or frostpoint accordingly.         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dewfrostpoint8(pres,rsat,useice)
      use consts_coms, only: ep8,toodry8

      implicit none
      real(kind=8), intent(in)           :: pres, rsat
      logical     , intent(in), optional :: useice
      real(kind=8)                       :: rsatoff, pvsat

      rsatoff       = max(toodry8,rsat)
      pvsat         = pres*rsatoff / (ep8 + rsatoff)
      if (present(useice)) then
         dewfrostpoint8 = tslif8(pvsat,useice)
      else
         dewfrostpoint8 = tslif8(pvsat)
      end if
      return
   end function dewfrostpoint8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based on the pressure [Pa], tem-   !
   ! perature [K] and relative humidity [fraction]. IT ALWAYS ASSUME THAT RELATIVE HUMI-   !
   ! DITY IS WITH RESPECT TO THE LIQUID PHASE.  ptrh2rvapil checks which one to use        !
   ! depending on whether temperature is more or less than the triple point.               !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function ptrh2rvapl8(relh,pres,temp)
      use consts_coms, only: ep8,toodry8

      implicit none
      real(kind=8), intent(in)           :: relh, pres, temp
      real(kind=8)                       :: rsath, relhh

      rsath = max(toodry8,rslf8(pres,temp))
      relhh = min(1.d0,max(0.d0,relh))

      if (newthermo) then
         !----- Considering that Rel.Hum. is based on vapour pressure ---------------------!
         ptrh2rvapl8 = max(toodry8,ep8 * relhh * rsath / (ep8 + (1.d0-relhh)*rsath))
      else
         !----- Original RAMS way to compute mixing ratio ---------------------------------!
         ptrh2rvapl8 = max(toodry8,relhh*rsath)
      end if

      return
   end function ptrh2rvapl8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based on the pressure [Pa], tem-   !
   ! perature [K] and relative humidity [fraction]. IT ALWAYS ASSUME THAT RELATIVE HUMI-   !
   ! DITY IS WITH RESPECT TO THE ICE PHASE. ptrh2rvapil checks which one to use depending  !
   ! on whether temperature is more or less than the triple point.                         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function ptrh2rvapi8(relh,pres,temp)
      use consts_coms, only: ep8,toodry8

      implicit none
      real(kind=8), intent(in)           :: relh, pres, temp
      real(kind=8)                       :: rsath, relhh

      rsath = max(toodry8,rsif8(pres,temp))
      relhh = min(1.d0,max(0.d0,relh))

      if (newthermo) then
         !----- Considering that Rel.Hum. is based on vapour pressure ---------------------!
         ptrh2rvapi8 = max(toodry8,ep8 * relhh * rsath / (ep8 + (1.d0-relhh)*rsath))
      else
         !----- Original RAMS way to compute mixing ratio ---------------------------------!
         ptrh2rvapi8 = max(toodry8,relhh*rsath)
      end if

      return
   end function ptrh2rvapi8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based on the pressure [Pa], tem-   !
   ! perature [K] and relative humidity [fraction]. It will check the temperature to       !
   ! decide between ice or liquid saturation and whether ice should be considered.         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function ptrh2rvapil8(relh,pres,temp,useice)
      use consts_coms, only: ep8,toodry8,t3ple8

      implicit none
      real(kind=8), intent(in)           :: relh, pres, temp
      logical     , intent(in), optional :: useice
      real(kind=8)                       :: rsath, relhh
      logical                            :: brrr_cold

      !----- Checking whether I use the user or the default check for ice saturation. -----!
      if (present(useice)) then
         brrr_cold = useice .and. temp < t3ple8
      else
         brrr_cold = bulk_on .and. temp < t3ple8
      end if

      if (brrr_cold) then
         rsath = max(toodry8,rsif8(pres,temp))
      else
         rsath = max(toodry8,rslf8(pres,temp))
      end if

      relhh = min(1.d0,max(0.d0,relh))
      
      ptrh2rvapil8 = max(toodry8,ep8 * relhh * rsath / (ep8 + (1.d0-relhh)*rsath))

      return
   end function ptrh2rvapil8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio. Two important points:                              !
   ! 1. IT ALWAYS ASSUME THAT RELATIVE HUMIDITY IS WITH RESPECT TO THE LIQUID PHASE.       !
   !    If you want to switch between ice and liquid, use rehuil instead.                  !
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rehul8(pres,temp,rvpr)
      use consts_coms, only: ep8,toodry8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres     ! Air pressure                [    Pa]
      real(kind=8), intent(in)           :: temp     ! Temperature                 [     K]
      real(kind=8), intent(in)           :: rvpr     ! Vapour mixing ratio         [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: rvprsat  ! Saturation mixing ratio     [ kg/kg]
      !------------------------------------------------------------------------------------!

      rvprsat = max(toodry8,rslf8(pres,temp))
      if (newthermo) then
         !----- This is based on relative humidity being defined with vapour pressure -----!
         rehul8 = max(0.d0,rvpr*(ep8+rvprsat)/(rvprsat*(ep8+rvpr)))
      else
         !----- Original formula used by RAMS ---------------------------------------------!
         rehul8 = max(0.d0,rvpr/rvprsat)
      end if
      return
   end function rehul8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio. Two important points:                              !
   ! 1. IT ALWAYS ASSUME THAT RELATIVE HUMIDITY IS WITH RESPECT TO THE ICE PHASE.          !
   !    If you want to switch between ice and liquid, use rehuil instead.                  !
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rehui8(pres,temp,rvpr)
      use consts_coms, only: ep8,toodry8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres     ! Air pressure                [    Pa]
      real(kind=8), intent(in)           :: temp     ! Temperature                 [     K]
      real(kind=8), intent(in)           :: rvpr     ! Vapour mixing ratio         [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: rvprsat  ! Saturation mixing ratio     [ kg/kg]
      !------------------------------------------------------------------------------------!

      rvprsat = max(toodry8,rsif8(pres,temp))
      if (newthermo) then
         !----- This is based on relative humidity being defined with vapour pressure -----!
         rehui8 = max(0.d0,rvpr*(ep8+rvprsat)/(rvprsat*(ep8+rvpr)))
      else
         !----- Original formula used by RAMS ---------------------------------------------!
         rehui8 = max(0.d0,rvpr/rvprsat)
      end if
      return
   end function rehui8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio. Two important points:                              !
   ! 1. It may consider whether the temperature is above or below the freezing point       !
   !    to choose which saturation to use. It is possible to explicitly force not to use   !
   !    ice in case level is 2 or if you have reasons not to use ice (e.g. reading data    !
   !    that did not consider ice).
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rehuil8(pres,temp,rvap,useice)
      use consts_coms, only: t3ple8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres      ! Air pressure               [    Pa]
      real(kind=8), intent(in)           :: temp      ! Temperature                [     K]
      real(kind=8), intent(in)           :: rvap      ! Vapour mixing ratio        [ kg/kg]
      logical     , intent(in), optional :: useice    ! Should I consider ice?     [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: rvapsat   ! Saturation mixing ratio    [ kg/kg]
      logical                            :: brrr_cold ! I'll use ice sat. now      [   T|F]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    Checking whether I should go with ice or liquid saturation.                     !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple8
      else 
         brrr_cold = bulk_on .and. temp < t3ple8
      end if

      if (brrr_cold) then
         rehuil8 = rehui8(pres,temp,rvap)
      else
         rehuil8 = rehul8(pres,temp,rvap)
      end if

      return
   end function rehuil8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the actual temperature based on the virtual temperature and   !
   ! mixing ratio. Two notes:                                                              !
   ! 1. It will use the condensation effect in case the total mixing ratio is provided.    !
   ! 2. This can be used for virtual potential temperature, just give potential tempera-   !
   !    ture instead of temperature.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function tv2temp8(tvir,rvpr,rtot)
      use consts_coms, only: epi8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: tvir     ! Virtual temperature          [    K]
      real(kind=8), intent(in)           :: rvpr     ! Vapour mixing ratio          [kg/kg]
      real(kind=8), intent(in), optional :: rtot     ! Total mixing ratio           [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real                               :: rtothere ! Internal rtot                [kg/kg]

      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvpr
      end if

      tv2temp8 = tvir * (1.d0 + rtothere) / (1.d0 + epi8*rvpr)

      return
   end function tv2temp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the virtual temperature based on the temperature and mixing   !
   ! ratio. Two notes:                                                                     !
   ! 1. It will use the condensation in case the total mixing ratio is provided.           !
   ! 2. This can be used for virtual potential temperature, just give potential tempera-   !
   !    ture instead of temperature.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function virtt8(temp,rvpr,rtot)
      use consts_coms, only: epi8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: temp     ! Temperature                  [    K]
      real(kind=8), intent(in)           :: rvpr     ! Vapour mixing ratio          [kg/kg]
      real(kind=8), intent(in), optional :: rtot     ! Total mixing ratio           [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real                       :: rtothere ! Internal rtot, to deal with optional [kg/kg]

      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvpr
      end if

      virtt8 = temp * (1.d0 + epi8 * rvpr) / (1.d0 + rtothere)

      return
   end function virtt8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the density based on the virtual temperature and the ideal    !
   ! gas law. The condensed phase will be taken into account if the user provided both     !
   ! the vapour and the total mixing ratios.                                               !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function idealdens8(pres,temp,rvpr,rtot)
      use consts_coms, only: rdry8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres ! Pressure                         [   Pa]
      real(kind=8), intent(in)           :: temp ! Temperature                      [    K]
      real(kind=8), intent(in)           :: rvpr ! Vapour mixing ratio              [kg/kg]
      real(kind=8), intent(in), optional :: rtot ! Total mixing ratio               [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=8)                       :: tvir ! Virtual temperature              [    K]
      !------------------------------------------------------------------------------------!
      if (present(rtot)) then
        tvir = virtt8(temp,rvpr,rtot)
      else
        tvir = virtt8(temp,rvpr)
      end if

      idealdens8 = pres / (rdry8 * tvir)

      return
   end function idealdens8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the density based on the virtual temperature and the ideal    !
   ! gas law.  The only difference between this function and the one above is that here we !
   ! provide vapour and total specific mass (specific humidity) instead of mixing ratio.   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function idealdenssh8(pres,temp,qvpr,qtot)
      use consts_coms, only : rdry8 & ! intent(in)
                            , epi8  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres ! Pressure                         [   Pa]
      real(kind=8), intent(in)           :: temp ! Temperature                      [    K]
      real(kind=8), intent(in)           :: qvpr ! Vapour specific mass             [kg/kg]
      real(kind=8), intent(in), optional :: qtot ! Total water specific mass        [kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: qall ! Either qtot or qvpr...           [kg/kg]
      !------------------------------------------------------------------------------------!

      if (present(qtot)) then
        qall = qtot
      else
        qall = qvpr
      end if

      idealdenssh8 = pres / (rdry8 * temp * (1.d0 - qall + epi8 * qvpr))

      return
   end function idealdenssh8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes reduces the pressure from the reference height to the      !
   ! canopy height by assuming hydrostatic equilibrium.                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function reducedpress8(pres,thetaref,shvref,zref,thetacan,shvcan,zcan)
      use consts_coms, only : epim18    & ! intent(in)
                            , p00k8     & ! intent(in)
                            , rocp8     & ! intent(in)
                            , cpor8     & ! intent(in)
                            , cp8       & ! intent(in)
                            , grav8     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)   :: pres     ! Pressure                        [        Pa]
      real(kind=8), intent(in)   :: thetaref ! Potential temperature           [         K]
      real(kind=8), intent(in)   :: shvref   ! Vapour specific mass            [     kg/kg]
      real(kind=8), intent(in)   :: zref     ! Height at reference level       [         m]
      real(kind=8), intent(in)   :: thetacan ! Potential temperature           [         K]
      real(kind=8), intent(in)   :: shvcan   ! Vapour specific mass            [     kg/kg]
      real(kind=8), intent(in)   :: zcan     ! Height at canopy level          [         m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=8)               :: pinc     ! Pressure increment              [ Pa^(R/cp)]
      real(kind=8)               :: thvbar   ! Average virtual pot. temper.    [         K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      First we compute the average virtual potential temperature between the canopy !
      ! top and the reference level.                                                       !
      !------------------------------------------------------------------------------------!
      thvbar = 5.d-1 * ( thetaref * (1.d0 + epim18 * shvref)                               &
                       + thetacan * (1.d0 + epim18 * shvcan))

      !----- Then, we find the pressure gradient scale. -----------------------------------!
      pinc = grav8 * p00k8 * (zref - zcan) / (cp8 * thvbar)

      !----- And we can find the reduced pressure. ----------------------------------------!
      reducedpress8 = (pres**rocp8 + pinc ) ** cpor8

      return
   end function reducedpress8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the enthalpy given the pressure, temperature, vapour       !
   ! specific humidity, and height.  Currently it doesn't compute mixed phase air, but     !
   ! adding it should be straight forward (finding the inverse is another story...).       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function ptqz2enthalpy8(pres,temp,qvpr,zref)
      use consts_coms, only : ep8       & ! intent(in)
                            , grav8     & ! intent(in)
                            , t3ple8    & ! intent(in)
                            , eta3ple8  & ! intent(in)
                            , cimcp8    & ! intent(in)
                            , clmcp8    & ! intent(in)
                            , cp8       & ! intent(in)
                            , alvi8     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres  ! Pressure                                 [    Pa]
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=8), intent(in) :: qvpr  ! Vapour specific mass                     [ kg/kg]
      real(kind=8), intent(in) :: zref  ! Reference height                         [     m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=8)             :: tequ  ! Dew-frost temperature                    [     K]
      real(kind=8)             :: pequ  ! Equlibrium vapour pressure               [    Pa]
      !------------------------------------------------------------------------------------!

      !----- First, we find the equilibrium vapour pressure and dew/frost point. ----------!
      pequ = pres * qvpr / (ep8 + (1.d0 - ep8) * qvpr)
      tequ = tslif8(pequ)

      !------------------------------------------------------------------------------------!
      !     Then, based on dew/frost point, we compute the enthalpy. This accounts whether !
      ! we would have to dew or frost formation if the temperature dropped to the          !
      ! equilibrium point.  Notice that if supersaturation exists, this will still give a  !
      ! number that makes sense, similar to the internal energy of supercooled water.      !
      !------------------------------------------------------------------------------------!
      if (tequ <= t3ple8) then
         ptqz2enthalpy8 = cp8 * temp + qvpr * (cimcp8 * tequ + alvi8   ) + grav8 * zref
      else
         ptqz2enthalpy8 = cp8 * temp + qvpr * (clmcp8 * tequ + eta3ple8) + grav8 * zref
      end if

      return
   end function ptqz2enthalpy8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature given the enthalpy, pressure, vapour       !
   ! specific humidity, and reference height.  Currently it doesn't compute mixed phase    !
   ! air, but adding it wouldn't be horribly hard, though it would require some root       !
   ! finding.                                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function hpqz2temp8(enthalpy,pres,qvpr,zref)
      use consts_coms, only : ep8       & ! intent(in)
                            , grav8     & ! intent(in)
                            , t3ple8    & ! intent(in)
                            , eta3ple8  & ! intent(in)
                            , cimcp8    & ! intent(in)
                            , clmcp8    & ! intent(in)
                            , cpi8      & ! intent(in)
                            , alvi8     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: enthalpy ! Enthalpy...                           [  J/kg]
      real(kind=8), intent(in) :: pres     ! Pressure                              [    Pa]
      real(kind=8), intent(in) :: qvpr     ! Vapour specific mass                  [ kg/kg]
      real(kind=8), intent(in) :: zref     ! Reference height                      [     m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=8)             :: tequ     ! Dew-frost temperature                 [     K]
      real(kind=8)             :: pequ     ! Equlibrium vapour pressure            [    Pa]
      !------------------------------------------------------------------------------------!

      !----- First, we find the equilibrium vapour pressure and dew/frost point. ----------!
      pequ = pres * qvpr / (ep8 + (1.d0 - ep8) * qvpr)
      tequ = tslif8(pequ)

      !------------------------------------------------------------------------------------!
      !     Then, based on dew/frost point, we compute the temperature. This accounts      !
      ! whether we would have to dew or frost formation if the temperature dropped to the  !
      ! equilibrium point.  Notice that if supersaturation exists, this will still give a  !
      ! temperature that makes sense (but less than the dew/frost point), similar to the   !
      ! internal energy of supercooled water.                                              !
      !------------------------------------------------------------------------------------!
      if (tequ <= t3ple8) then
         hpqz2temp8 = cpi8 * (enthalpy - qvpr * (cimcp8 * tequ + alvi8   ) - grav8 * zref)
      else
         hpqz2temp8 = cpi8 * (enthalpy - qvpr * (clmcp8 * tequ + eta3ple8) - grav8 * zref)
      end if

      return
   end function hpqz2temp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the temperature given the potential temperature, density, and !
   ! specific humidity.  This comes from a combination of the definition of potential      !
   ! temperature and the ideal gas law, to eliminate pressure, when pressure is also       !
   ! unknown.                                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function thrhsh2temp8(theta,dens,qvpr)
      use consts_coms, only : cpocv8  & ! intent(in)
                            , p00i8   & ! intent(in)
                            , rdry8   & ! intent(in)
                            , epim18  & ! intent(in)
                            , rocv8   ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: theta    ! Potential temperature                 [     K]
      real(kind=8), intent(in) :: dens     ! Density                               [    Pa]
      real(kind=8), intent(in) :: qvpr     ! Specific humidity                     [ kg/kg]
      !------------------------------------------------------------------------------------!

      thrhsh2temp8 = theta ** cpocv8                                                       &
                   * (p00i8 * dens * rdry8 * (1.d0 + epim18 * qvpr)) ** rocv8

      return
   end function thrhsh2temp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the ice liquid potential temperature given the Exner       !
   ! function [J/kg/K], temperature [K], and liquid and ice mixing ratios [kg/kg].         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function theta_iceliq8(exner,temp,rliq,rice)
      use consts_coms, only: alvl8, alvi8, cp8, ttripoli8, htripoli8, htripolii8

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner ! Exner function                           [J/kg/K]
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=8), intent(in) :: rliq  ! Liquid mixing ratio                      [ kg/kg]
      real(kind=8), intent(in) :: rice  ! Ice mixing ratio                         [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: hh    ! Enthalpy associated with sensible heat   [  J/kg]
      real(kind=8)             :: qq    ! Enthalpy associated with latent heat     [  J/kg]
      !------------------------------------------------------------------------------------!

      !----- Finding the enthalpies -------------------------------------------------------!
      hh = cp8 * temp
      qq  = alvl8*rliq + alvi8 * rice
      
      if (newthermo) then
      
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli8) then
            theta_iceliq8 = hh * exp(-qq/hh) / exner
         else
            theta_iceliq8 = hh * exp(-qq * htripolii8) / exner
         end if
      else
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli8) then
            theta_iceliq8 = hh * hh / (exner * ( hh + qq))
         else
            theta_iceliq8 = hh * htripoli8 / (exner * ( htripoli8 + qq))
         end if
      end if

      return
   end function theta_iceliq8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the liquid potential temperature derivative with respect   !
   ! to temperature, useful in iterative methods.                                          !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dthetail_dt8(condconst,thil,exner,pres,temp,rliq,ricein)
      use consts_coms, only: alvl8, alvi8, cp8, ttripoli8,htripoli8,htripolii8,t3ple8

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      logical     , intent(in)           :: condconst  ! Condensation is constant? [   T|F]
      real(kind=8), intent(in)           :: thil       ! Ice liquid pot. temp.     [     K]
      real(kind=8), intent(in)           :: exner      ! Exner function            [J/kg/K]
      real(kind=8), intent(in)           :: pres       ! Pressure                  [    Pa]
      real(kind=8), intent(in)           :: temp       ! Temperature               [     K]
      real(kind=8), intent(in)           :: rliq       ! Liquid mixing ratio       [ kg/kg]
      real(kind=8), intent(in), optional :: ricein     ! Ice mixing ratio          [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: rice       ! Ice mixing ratio or 0.    [ kg/kg]
      real(kind=8)                       :: ldrst      ! L × d(rs)/dT × T          [  J/kg]
      real(kind=8)                       :: hh         ! Sensible heat enthalpy    [  J/kg]
      real(kind=8)                       :: qq         ! Latent heat enthalpy      [  J/kg]
      logical                            :: thereisice ! Is ice present            [   ---]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !   Checking whether I should consider ice or not.                                   !
      !------------------------------------------------------------------------------------!
      thereisice = present(ricein)
      
      if (thereisice) then
         rice=ricein
      else
         rice=0.d0
      end if
      
      !----- No condensation, dthetail_dt is a constant -----------------------------------!
      if (rliq+rice == 0.d0) then
         dthetail_dt8 = thil/temp
         return
      else
         hh    = cp8  * temp                            !----- Sensible heat enthalpy
         qq    = alvl8* rliq + alvi8 * rice             !----- Latent heat enthalpy
         !---------------------------------------------------------------------------------!
         !    This is the term L×[d(rs)/dt]×T. L may be either the vapourisation or        !
         ! sublimation latent heat, depending on the temperature and whether we are consi- !
         ! dering ice or not. Also, if condensation mixing ratio is constant, then this    !
         ! term will be always zero.                                                       !
         !---------------------------------------------------------------------------------!
         if (condconst) then
            ldrst = 0.d0
         elseif (thereisice .and. temp < t3ple8) then
            ldrst = alvi8*rsifp8(pres,temp)*temp
         else
            ldrst = alvl8*rslfp8(pres,temp)*temp  
         end if
      end if

      if (newthermo) then
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli8) then
            dthetail_dt8 = thil * (1. + (ldrst + qq)/hh) / temp
         else
            dthetail_dt8 = thil * (1. + ldrst*htripolii8) / temp
         end if
      else
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli8) then
            dthetail_dt8 = thil * (1.d0 + (ldrst + qq)/(hh+qq)) / temp
         else
            dthetail_dt8 = thil * (1.d0 + ldrst/(htripoli8 + alvl8*rliq)) / temp
         end if
      end if

      return
   end function dthetail_dt8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes temperature from the ice-liquid water potential temperature !
   !  in Kelvin, Exner function in J/kg/K, and liquid and ice mixing ratios in kg/kg.      !
   !    For now t1stguess is used only to decide whether I should use the complete case or !
   ! the 253 K to reduce the error on neglecting the changes on latent heat due to temper- !
   ! ature.                                                                                !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function thil2temp8(thil,exner,pres,rliq,rice,t1stguess)
      use consts_coms, only : cp8        & ! intent(in)
                            , cpi8       & ! intent(in)
                            , alvl8      & ! intent(in)
                            , alvi8      & ! intent(in)
                            , t008       & ! intent(in)
                            , t3ple8     & ! intent(in)
                            , ttripoli8  & ! intent(in)
                            , htripolii8 & ! intent(in)
                            , cpi48      ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: thil      ! Ice-liquid water potential temp.     [     K]
      real(kind=8), intent(in) :: exner     ! Exner function                       [J/kg/K]
      real(kind=8), intent(in) :: pres      ! Pressure                             [    Pa]
      real(kind=8), intent(in) :: rliq      ! Liquid water mixing ratio            [ kg/kg]
      real(kind=8), intent(in) :: rice      ! Ice mixing ratio                     [ kg/kg]
      real(kind=8), intent(in) :: t1stguess ! 1st. guess for temperature           [     K]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)             :: deriv     ! Function derivative 
      real(kind=8)             :: fun       ! Function for which we seek a root.
      real(kind=8)             :: funa      ! Smallest  guess function
      real(kind=8)             :: funz      ! Largest   guess function
      real(kind=8)             :: tempa     ! Smallest  guess (or previous guess in Newton)
      real(kind=8)             :: tempz     ! Largest   guess (or new guess in Newton)
      real(kind=8)             :: delta     ! Aux. var to compute 2nd guess for bisection
      integer                  :: itn,itb   ! Iteration counter
      logical                  :: converged ! Convergence handle
      logical                  :: zside     ! Flag to check for one-sided approach...
      real(kind=8)             :: til       ! Ice liquid temperature               [     K]
      !------------------------------------------------------------------------------------!


      !----- 1st. of all, check whether there is condensation. If not, theta_il = theta ---!
      if (rliq+rice == 0.d0) then
         thil2temp8 = cpi8 * thil * exner
         return
      !----- If not, check whether we are using the old thermo or the new one -------------!
      elseif (.not. newthermo) then
         til = cpi8 * thil * exner
         if (t1stguess > ttripoli8) then
            thil2temp8 = 5.d-1 * (til + sqrt(til * ( til                                   &
                                                   + cpi48 * (alvl8*rliq + alvi8*rice))))
         else
            thil2temp8 = til * ( 1.d0 + (alvl8*rliq+alvi8*rice) * htripolii8)
         end if
         return
      end if
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(5(a,1x,f11.4,1x))')                                              &
      !   'Input values: Exner =',exner,'thil=',thil,'rliq=',1000.*rliq,'rice=',1000.*rice  &
      !           ,'t1stguess=',t1stguess
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- If not, iterate: For the Newton's 1st. guess, use t1stguess ------------------!
      tempz     = t1stguess
      fun       = theta_iceliq8(exner,tempz,rliq,rice)
      deriv     = dthetail_dt8(.true.,fun,exner,pres,tempz,rliq,rice)
      fun       = fun - thil

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,2(1x,a,1x,es12.5))')            &
      !   'itn=',0,'bisection=',.false.,'tempz=',tempz-t00                                  &
      !           ,'fun=',fun,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (fun == 0.d0) then
         thil2temp8 = tempz
         converged  = .true.
         return
      else 
         tempa = tempz
         funa  = fun
      end if
      !----- Enter loop: it will probably skip when the air is not saturated --------------!
      converged = .false.
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, go to bisection ------!
         tempa = tempz
         funa  = fun
         tempz = tempa - fun/deriv
         fun   = theta_iceliq8(exner,tempz,rliq,rice)
         deriv = dthetail_dt8(.true.,fun,exner,pres,tempz,rliq,rice)
         fun   = fun - thil
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,2(1x,a,1x,es12.5))')         &
         !   'itn=',itn,'bisection=',.false.,'tempz=',tempz-t00                             &
         !           ,'fun=',fun,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         converged = abs(tempa-tempz) < toler8 * tempz
         !----- Converged, happy with that, return the average b/w the 2 previous guesses -!
         if (fun == 0.d0) then
            thil2temp8 = tempz
            converged = .true.
            return
         elseif(converged) then
            thil2temp8 = 5.d-1 * (tempa+tempz)
            return
         end if
      end do newloop
      
      !------------------------------------------------------------------------------------!
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < toler8 * tempa) then
            delta = 1.d2 * toler8 * tempa
         else 
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)), 1.d2 * toler8 * tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
            funz  = theta_iceliq8(exner,tempz,rliq,rice) - thil
            zside = funa * funz < 0.d0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,i14,1x))')    'itn  =',itn  ,'itb =',itb
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'thil =',thil ,'t1st=',t1stguess
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'exner=',exner,'pres=',pres
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'rliq =',rliq ,'rice=',rice
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'tempz=',tempz,'funz=',funz
            write (unit=*,fmt='(1(a,1x,es14.7,1x))') 'delta=',delta
            call fatal_error('Failed finding the second guess for regula falsi'            &
                          ,'thil2temp8','therm_lib8.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         thil2temp8 =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(thil2temp8-tempa) < toler8 * thil2temp8
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         fun  = theta_iceliq8(exner,tempz,rliq,rice) - thil

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,6(1x,a,1x,f11.4))')         &
         !   'itn=',itb,'bisection=',.true.                                                 &
         !           ,'temp=',thil2temp-t00,'tempa=',tempa-t00,'tempz=',tempz-t00           &
         !           ,'fun=',fun,'funa=',funa,'funz=',funz
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0.d0 ) then
            tempz = thil2temp8
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 5.d-1
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = thil2temp8
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz=funz * 5.d-1
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(a)') ' '
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      if (.not.converged) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Temperature finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'theta_il        [     K] =',thil
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Exner           [J/kg/K] =',exner
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rliq            [  g/kg] =',1.d3*rliq
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rice            [  g/kg] =',1.d3*rice
         write (unit=*,fmt='(a,1x,f12.4)' ) 't1stguess       [    °C] =',t1stguess-t008
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tempa           [    °C] =',tempa-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tempz           [    °C] =',tempz-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [     K] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [     K] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [     K] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [  ----] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [  ----] =',toler8
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                                          ,abs(thil2temp8-tempa)/thil2temp8
         write (unit=*,fmt='(a,1x,f12.4)' ) 'thil2temp8      [     K] =',thil2temp8

         call fatal_error('Temperature didn''t converge, giving up!!!'                     &
                       ,'thil2temp8','therm_lib8.f90')
      end if

      return
   end function thil2temp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the partial derivative of temperature as a function of the !
   ! saturation mixing ratio [kg/kg],, keeping pressure constant. This is based on the     !
   ! ice-liquid potential temperature equation.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dtempdrs8(exner,thil,temp,rliq,rice,rconmin)
      use consts_coms, only: alvl8, alvi8, cp8, cpi8, ttripoli8, htripolii8

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner   ! Exner function                         [J/kg/K]
      real(kind=8), intent(in) :: thil    ! Ice-liquid potential temperature (*)   [     K]
      real(kind=8), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=8), intent(in) :: rliq    ! Liquid mixing ratio                    [ kg/kg]
      real(kind=8), intent(in) :: rice    ! Ice mixing ratio                       [ kg/kg]
      real(kind=8), intent(in) :: rconmin ! Minimum non-zero cond. mixing ratio    [ kg/kg]
      !------------------------------------------------------------------------------------!
      ! (*) Thil is not used in this formulation but it may be used should you opt by      !
      !     other ways to compute theta_il, so don't remove this argument.                 !
      !------------------------------------------------------------------------------------!
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: qhydm   ! Enthalpy associated with latent heat   [  J/kg]
      real(kind=8)             :: hh      ! Enthalpy associated with sensible heat [  J/kg]
      real(kind=8)             :: rcon    ! Condensate mixing ratio                [ kg/kg]
      real(kind=8)             :: til     ! Ice-liquid temperature                 [     K]
      !------------------------------------------------------------------------------------!
            
      !----- Finding the temperature and hydrometeor terms --------------------------------!
      qhydm = alvl8 * rliq + alvi8 * rice
      rcon  = rliq+rice
      if (rcon < rconmin) then
         dtempdrs8 = 0.d0
      elseif (newthermo) then
         hh    = cp8 * temp
         !---------------------------------------------------------------------------------!
         !    Deciding how to compute, based on temperature and whether condensates exist. !
         !---------------------------------------------------------------------------------!
         if (temp > ttripoli8) then
            dtempdrs8 = - temp * qhydm / (rcon * (hh+qhydm))
         else
            dtempdrs8 = - temp * qhydm * htripolii8 / rcon
         end if
      else
         til   = cpi8 * thil * exner
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli8) then
            dtempdrs8 = - til * qhydm /( rcon * cp8 * (2.d0*temp-til))
         else
            dtempdrs8 = - til * qhydm * htripolii8 / rcon
         end if
      end if

      return
   end function dtempdrs8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the change of ice-liquid potential temperature due to      !
   ! sedimentation. The arguments are ice-liquid potential temperature, potential temper-  !
   ! ature and temperature in Kelvin, the old and new mixing ratio [kg/kg] and the old and !
   ! new enthalpy [J/kg].                                                                  !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dthil_sedimentation8(thil,theta,temp,rold,rnew,qrold,qrnew)
      use consts_coms, only: ttripoli8,cp8,alvi8,alvl8

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: thil  ! Ice-liquid potential temperature         [     K]
      real(kind=8), intent(in) :: theta ! Potential temperature                    [     K]
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=8), intent(in) :: rold  ! Old hydrometeor mixing ratio             [ kg/kg]
      real(kind=8), intent(in) :: rnew  ! New hydrometeor mixing ratio             [ kg/kg]
      real(kind=8), intent(in) :: qrold ! Old hydrometeor latent enthalpy          [  J/kg]
      real(kind=8), intent(in) :: qrnew ! New hydrometeor latent enthalpy          [  J/kg]
      !------------------------------------------------------------------------------------!

      if (newthermo) then
         dthil_sedimentation8 = - thil * (alvi8 * (rnew-rold) - (qrnew-qrold))             &
                                        / (cp8 * max(temp,ttripoli8))
      else
         dthil_sedimentation8 = - thil*thil * (alvi8*(rnew-rold) - (qrnew-qrold))          &
                                        / (cp8 * max(temp,ttripoli8) * theta)
      end if

      return
   end function dthil_sedimentation8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the ice-vapour equivalent potential temperature from        !
   !  theta_iland the total mixing ratio. This is equivalent to the equivalent potential   !
   ! temperature considering also the effects of fusion/melting/sublimation.               !
   !    In case you want to find thetae (i.e. without ice) simply provide the logical      !
   ! useice as .false. .                                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function thetaeiv8(thil,pres,temp,rvap,rtot,useice)
      use consts_coms, only : alvl8,alvi8,cp8,ep8,p008,rocp8,ttripoli8,t3ple8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: thil   ! Ice-liquid water pot. temp.   [     K]
      real(kind=8), intent(in)           :: pres   ! Pressure                      [    Pa]
      real(kind=8), intent(in)           :: temp   ! Temperature                   [     K]
      real(kind=8), intent(in)           :: rvap   ! Water vapour mixing ratio     [ kg/kg]
      real(kind=8), intent(in)           :: rtot   ! Total mixing ratio            [ kg/kg]
      logical     , intent(in), optional :: useice ! Should I use ice?             [   T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)                       :: tlcl   ! Internal LCL temperature      [     K]
      real(kind=8)                       :: plcl   ! Lifting condensation pressure [    Pa]
      real(kind=8)                       :: dzlcl  ! Thickness of layer beneath LCL[     m]
      !------------------------------------------------------------------------------------!

      if (present(useice)) then
         call lcl_il8(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      else
         call lcl_il8(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl)
      end if

      !------------------------------------------------------------------------------------!
      !     The definition of the thetae_iv is the thetae_ivs at the LCL. The LCL, in turn !
      ! is the point in which rtot = rvap = rsat, so at the LCL rliq = rice = 0.           !
      !------------------------------------------------------------------------------------!
      thetaeiv8  = thetaeivs8(thil,tlcl,rtot,0.d0,0.d0)

      return
   end function thetaeiv8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of ice-vapour equivalent potential tempera-  !
   ! ture, based on the expression used to compute the ice-vapour equivalent potential     !
   ! temperature (function thetaeiv).                                                      !
   !                                                                                       !
   !    IMPORTANT!!! This CANNOT BE USED to compute d(Thetae_ivs)/dT, because here         !
   !                 we assume that T(LCL) and saturation mixing ratio are known and       !
   !                 constants, and that the LCL pressure (actually the saturation  vapour !
   !                 pressure at the LCL) is a function of temperature. In case you want   !
   !                 d(Thetae_ivs)/dT, use the dthetaeivs_dt function instead.             !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dthetaeiv_dtlcl8(theiv,tlcl,rtot,eslcl,useice)
      use consts_coms, only : rocp8,aklv8,ttripoli8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: theiv    ! Ice-vapour equiv. pot. temp. [    K]
      real(kind=8), intent(in)           :: tlcl     ! LCL temperature              [    K]
      real(kind=8), intent(in)           :: rtot     ! Total mixing ratio (rs @ LCL)[   Pa]
      real(kind=8), intent(in)           :: eslcl    ! LCL saturation vapour press. [   Pa]
      logical     , intent(in), optional :: useice   ! Flag for considering ice     [  T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: desdtlcl ! Sat. vapour pres. deriv.     [ Pa/K]
      !------------------------------------------------------------------------------------!



      !----- Finding the derivative of rs with temperature --------------------------------!
      if (present(useice)) then
         desdtlcl = eslifp8(tlcl,useice)
      else
         desdtlcl = eslifp8(tlcl)
      end if



      !----- Finding the derivative. Depending on the temperature, use different eqn. -----!
      if (tlcl > ttripoli8) then
         dthetaeiv_dtlcl8 = theiv * (1.d0 - rocp8*tlcl*desdtlcl/eslcl - aklv8*rtot/tlcl)   &
                          / tlcl
      else
         dthetaeiv_dtlcl8 = theiv * (1.d0 - rocp8*tlcl*desdtlcl/eslcl                  )   &
                          / tlcl
      end if

      return
   end function dthetaeiv_dtlcl8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the saturation ice-vapour equivalent potential temperature  !
   ! from theta_il and the total mixing ratio (split into saturated vapour plus liquid and !
   ! ice. This is equivalent to the equivalent potential temperature considering also the  !
   ! effects of fusion/melting/sublimation, and it is done separatedly from the regular    !
   ! thetae_iv because it doesn't require iterations.                                      !
   !                                                                                       !
   !    References:                                                                        !
   !    Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential tem- !
   !        perature as a thermodynamic variable in deep atmospheric models. Mon. Wea.     !
   !        Rev., v. 109, 1094-1102. (TC81)                                                !
   !                                                                                       !
   !    Some algebra was needed to find this equation, essentially combining (TC81-26) and !
   ! (TC81-27), and the conservation of total water (TC81-16). It assumes that the divi-   !
   ! sion between the three phases is already taken care of.                               !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function thetaeivs8(thil,temp,rsat,rliq,rice)
      use consts_coms, only : aklv8, ttripoli8
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in)   :: thil     ! Ice-liquid water potential temp.    [     K]
      real(kind=8), intent(in)   :: temp     ! Temperature                         [     K]
      real(kind=8), intent(in)   :: rsat     ! Sat. water vapour mixing ratio      [ kg/kg]
      real(kind=8), intent(in)   :: rliq     ! Liquid water mixing ratio           [ kg/kg]
      real(kind=8), intent(in)   :: rice     ! Ice mixing ratio                    [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)               :: rtots    ! Saturated mixing ratio              [     K]
      !------------------------------------------------------------------------------------!

      rtots      = rsat+rliq+rice
      
      thetaeivs8 = thil * exp ( aklv8 * rtots / max(temp,ttripoli8))

      return
   end function thetaeivs8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of saturation ice-vapour equivalent          !
   ! potential temperature, based on the expression used to compute the saturation         !
   ! ice-vapour equivalent potential temperature (function thetaeivs).                     !
   !                                                                                       !
   !    IMPORTANT!!! This CANNOT BE USED to compute d(Thetae_iv)/d(T_LCL), because here    !
   !                 we assume that temperature and pressure are known and constants, and  !
   !                 that the mixing ratio is a function of temperature. In case you want  !
   !                 d(Thetae_iv)/d(T_LCL), use the dthetaeiv_dtlcl function instead.      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dthetaeivs_dt8(theivs,temp,pres,rsat,useice)
      use consts_coms, only : aklv8,alvl8,ttripoli8,htripolii8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: theivs ! Sat. ice-vap. eq. pot. temp. [      K]
      real(kind=8), intent(in)           :: temp   ! Temperature                  [      K]
      real(kind=8), intent(in)           :: pres   ! Pressure                     [     Pa]
      real(kind=8), intent(in)           :: rsat   ! Saturation mixing ratio      [  kg/kg]
      logical     , intent(in), optional :: useice ! Flag for considering ice     [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: drsdt  ! Saturated mixing ratio deriv.[kg/kg/K]
      !------------------------------------------------------------------------------------!


      !----- Finding the derivative of rs with temperature --------------------------------!
      if (present(useice)) then
         drsdt = rslifp8(pres,temp,useice)
      else
         drsdt = rslifp8(pres,temp)
      end if


      !----- Finding the derivative. Depending on the temperature, use different eqn. -----!
      if (temp > ttripoli8) then
         dthetaeivs_dt8 = theivs * (1.d0 + aklv8 * (drsdt*temp-rsat)/temp ) / temp
      else
         dthetaeivs_dt8 = theivs * (1.d0 + alvl8 * drsdt * temp * htripolii8 ) / temp
      end if

      
      return
   end function dthetaeivs_dt8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function finds the ice-liquid potential temperature from the ice-vapour equi- !
   ! valent potential temperature.                                                         !
   ! Important remarks:                                                                    !
   ! 1. If you don't want to use ice thermodynamics, simply force useice to be .false.     !
   !    Otherwise, the model will decide based on the LEVEL given by the user from their   !
   !    RAMSIN.                                                                            !
   ! 2. If rtot < rsat, then this will convert theta_e into theta, which can be thought as !
   !    a particular case.                                                                 !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function thetaeiv2thil8(theiv,pres,rtot,useice)
      use consts_coms, only : alvl8,cp8,ep8,p008,rocp8,ttripoli8,t3ple8,t008
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: theiv     ! Ice vap. equiv. pot. temp. [     K]
      real(kind=8), intent(in)           :: pres      ! Pressure                   [    Pa]
      real(kind=8), intent(in)           :: rtot      ! Total mixing ratio         [ kg/kg]
      logical     , intent(in), optional :: useice    ! Flag for ice thermo        [   T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)                       :: pvap      ! Sat. vapour pressure
      real(kind=8)                       :: theta     ! Potential temperature
      real(kind=8)                       :: deriv     ! Function derivative 
      real(kind=8)                       :: funnow    ! Function for which we seek a root.
      real(kind=8)                       :: funa      ! Smallest  guess function
      real(kind=8)                       :: funz      ! Largest   guess function
      real(kind=8)                       :: tlcla     ! Smallest  guess (or old guess)
      real(kind=8)                       :: tlclz     ! Largest   guess (or new guess)
      real(kind=8)                       :: tlcl      ! What will be the LCL temperature
      real(kind=8)                       :: es00      ! Defined as p00*rt/(epsilon + rt)
      real(kind=8)                       :: delta     ! Aux. variable (For 2nd guess).
      integer                            :: itn,itb   ! Iteration counters
      integer                            :: ii        ! Another counter
      logical                            :: converged ! Convergence handle
      logical                            :: zside     ! Aux. flag - sides for Regula Falsi
      logical                            :: brrr_cold ! Flag - considering ice thermo.
      !------------------------------------------------------------------------------------!
    
      !----- Filling the flag for ice thermo that will be always present ------------------!
      if (present(useice)) then
         brrr_cold = useice
      else
         brrr_cold = bulk_on
      end if
    
      !----- Finding es00, which is a constant --------------------------------------------!
      es00 = p008 * rtot / (ep8 + rtot)


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=36,fmt='(a)') '----------------------------------------------------------'
      !write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x))')                                  &
      !   'INPUT : it=',-1,'theiv=',theiv,'pres=',0.01*pres,'rtot=',rtot*1000.
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, we assume we are lucky and right at the LCL.                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rtot / (ep8 + rtot) 
      tlclz     = tslif8(pvap,brrr_cold)
      theta     = tlclz * (es00/pvap)**rocp8
      funnow    = thetaeivs8(theta,tlclz,rtot,0.d0,0.d0)
      deriv     = dthetaeiv_dtlcl8(funnow,tlclz,rtot,pvap,brrr_cold)
      funnow    = funnow - theiv

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')                &
      !   'NEWTON: it=',0,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Putting something in tlcla in case we never loop through Newton's method -----!
      tlcla     = tlclz
      funa      = funnow
      converged = .false.

      !----- Looping: Newton's iterative method -------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, skip to bisection ----!
         !----- Updating guesses ----------------------------------------------------------!
         tlcla   = tlclz
         funa    = funnow
         tlclz   = tlcla - funnow/deriv

         !----- Updating the function evaluation and its derivative -----------------------!
         pvap    = eslif8(tlclz,brrr_cold)
         theta   = tlclz * (es00/pvap)**rocp8
         funnow  = thetaeivs8(theta,tlclz,rtot,0.d0,0.d0)
         deriv   = dthetaeiv_dtlcl8(funnow,tlclz,rtot,pvap,brrr_cold)
         funnow  = funnow - theiv

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
         !   'NEWTON: it=',itn,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow           &
         !          ,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         converged = abs(tlcla-tlclz) < toler8 * tlclz
         if (funnow == 0.d0) then
            tlcl = tlclz
            funz = funnow
            converged = .true.
            exit newloop
         elseif (converged) then
            tlcl = 5.d-1 *(tlcla+tlclz)
            funz = funnow
            exit newloop
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         !----- Set funz, and check whether funa and funz already have opposite sign. -----!
         funz  = funnow
         zside = .true.
         if (funa*funnow > 0.d0) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler8 * tlcla) then
               delta = 1.d2 * toler8 * tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funz-funa)),1.d2 * toler8 * tlcla)
            end if
            tlclz = tlcla + delta

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')          &
            !   '2NGGSS: tt=',0,'tlclz=',tlclz-t00,'tlcla=',tlcla-t00,'pvap=',0.01*pvap     &
            !           ,'funa=',funa,'funz=',funnow,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlcla + dble((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif8(tlclz,brrr_cold)
               theta = tlclz * (es00/pvap)**rocp8
               funz  = thetaeivs8(theta,tlclz,rtot,0.d0,0.d0) - theiv

               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')       &
               !   '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz       &
               !           ,'delta=',delta
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

               zside = funa*funz < 0.d0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(a)') ' No second guess for you...'
               write (unit=*,fmt='(2(a,1x,i14,1x))')    'itn   =',itn   ,'itb   =',itb
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'theiv =',theiv ,'rtot  =',rtot
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'pres  =',pres  ,'pvap  =',pvap
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'theta =',theta ,'delta =',delta
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'tlcla =',tlcla ,'funa  =',funa
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'tlclz =',tlclz ,'funz  =',funz
               call fatal_error('Failed finding the second guess for regula falsi'         &
                             ,'thetaeiv2thil8','therm_lib8.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo

            !----- Updating the guess -----------------------------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)

            !----- Updating function evaluation -------------------------------------------!
            pvap   = eslif8(tlcl,brrr_cold)
            theta  = tlcl * (es00/pvap)**rocp8
            funnow = thetaeivs8(theta,tlcl,rtot,0.d0,0.d0) - theiv

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),1(a,1x,es11.4,1x))')          &
            !   'REGFAL: it=',itb,'tlcl =',tlcl -t00,'pvap=',0.01*pvap,'fun=',funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla) < toler8 * tlcl
            if (converged) then
               exit fpoloop
            elseif (funnow*funa < 0.d0) then 
               tlclz = tlcl
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 5.d-1
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tlcla = tlcl
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 5.d-1
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
         end do fpoloop
      end if

      if (converged) then 
         thetaeiv2thil8  = theiv * exp (- alvl8 * rtot / (cp8 * max(tlcl,ttripoli8)) )
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
         !   'ANSWER: itb=',itn,'tlcl=',tlcl-t00,'eslcl=',0.01*pvap                         &
         !          ,'thil=',thetaeiv2thil,'funa=',funa,'funz=',funz
         !write (unit=36,fmt='(a)') '-------------------------------------------------------'
         !write (unit=36,fmt='(a)') ' '
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      else
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
         write (unit=*,fmt='(a)')           ' THEIV2THIL8 failed!'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Input: '
         write (unit=*,fmt='(a,1x,f12.5)')  '    THEIV    [     K]:',theiv
         write (unit=*,fmt='(a,1x,f12.5)')  '    PRES     [    Pa]:',pres * 1.d2
         write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1.d3*rtot
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Output: '
         write (unit=*,fmt='(a,1x,i12)')    '    ITERATIONS       :',itb
         write (unit=*,fmt='(a,1x,f12.5)')  '    PVAP     [   hPa]:',pvap
         write (unit=*,fmt='(a,1x,f12.5)')  '    THETA    [     K]:',theta
         write (unit=*,fmt='(a,1x,f12.5)')  '    TLCL     [    °C]:',tlcl-t008
         write (unit=*,fmt='(a,1x,f12.5)')  '    TLCLA    [    °C]:',tlcla-t008
         write (unit=*,fmt='(a,1x,f12.5)')  '    TLCLZ    [    °C]:',tlclz-t008
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNA     [     K]:',funa
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNZ     [     K]:',funz
         write (unit=*,fmt='(a,1x,es12.5)') '    DERIV    [   ---]:',deriv
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_A    [   ---]:',abs(tlcl-tlcla)/tlcl
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_Z    [   ---]:',abs(tlcl-tlclz)/tlcl
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)

         call fatal_error('TLCL didn''t converge, gave up!','thetaeiv2thil8'               &
                         ,'therm_lib8.f90')
      end if

      return
   end function thetaeiv2thil8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine converts saturated ice-vapour equivalent potential temperature     !
   ! into temperature, given pressure in Pa, and theta_es in Kelvin. It also returns the   !
   ! potential temperature in Kelvin, and saturation vapour mixing ratio in kg/kg. As      !
   ! usual,  we seek T using Newton's method as a starting point, and if it fails, we fall !
   ! back to the modified regula falsi (Illinois method).                                  !
   !                                                                                       !
   ! OBS: In case you want to ignore ice, send useice as false. The default is to consider !
   !      when level >= 3 and to ignore otherwise.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine thetaeivs2temp8(theivs,pres,theta,temp,rsat,useice)
      use consts_coms, only : alvl8,cp8,ep8,p008,rocp8,ttripoli8,t008
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)            :: theivs     ! Sat. thetae_iv          [      K]
      real(kind=8), intent(in)            :: pres       ! Pressure                [     Pa]
      real(kind=8), intent(out)           :: theta      ! Potential temperature   [      K]
      real(kind=8), intent(out)           :: temp       ! Temperature             [      K]
      real(kind=8), intent(out)           :: rsat       ! Saturation mixing ratio [  kg/kg]
      logical     , intent(in) , optional :: useice     ! Flag for ice thermo     [    T|F]
      !----- Local variables, with other thermodynamic properties -------------------------!
      real(kind=8)                        :: exnernormi ! 1./ (Norm. Exner fctn)  [    ---]
      logical                             :: brrr_cold  ! Flag for ice thermo     [    T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)                        :: deriv      ! Function derivative 
      real(kind=8)                        :: funnow     ! Function for which we seek a root.
      real(kind=8)                        :: funa       ! Smallest  guess function
      real(kind=8)                        :: funz       ! Largest   guess function
      real(kind=8)                        :: tempa      ! Smallest  guess (or previous)
      real(kind=8)                        :: tempz      ! Largest   guess (or new)
      real(kind=8)                        :: delta      ! Aux. var. for 2nd guess finding.
      integer                             :: itn,itb    ! Iteration counters
      logical                             :: converged  ! Convergence handle
      logical                             :: zside      ! Check sides (Regula Falsi)
      !------------------------------------------------------------------------------------!
    
      !----- Setting up the ice check, in case useice is not present. ---------------------!
      if (present(useice)) then
         brrr_cold = useice
      else 
         brrr_cold = bulk_on
      end if
    
      !----- Finding the inverse of normalised Exner, which is constant in this routine ---!
      exnernormi = (p008 /pres) ** rocp8

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, no idea, guess 0°C.                                            !
      !------------------------------------------------------------------------------------!
      tempz     = t008
      theta     = tempz * exnernormi
      rsat      = rslif8(pres,tempz,brrr_cold)
      funnow    = thetaeivs8(theta,tempz,rsat,0.d0,0.d0)
      deriv     = dthetaeivs_dt8(funnow,tempz,pres,rsat,brrr_cold)
      funnow    = funnow - theivs

      !----- Saving here just in case Newton is aborted at the 1st guess ------------------!
      tempa     = tempz
      funa      = funnow

      converged = .false.
      !----- Looping ----------------------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, skip to bisection ----!
         !----- Updating guesses ----------------------------------------------------------!
         tempa   = tempz
         funa    = funnow
         
         tempz   = tempa - funnow/deriv
         theta   = tempz * exnernormi
         rsat    = rslif8(pres,tempz,brrr_cold)
         funnow  = thetaeivs8(theta,tempz,rsat,0.d0,0.d0)
         deriv   = dthetaeivs_dt8(funnow,tempz,pres,rsat,brrr_cold)
         funnow  = funnow - theivs

         converged = abs(tempa-tempz) < toler8 * tempz
         if (funnow == 0.d0) then
            converged =.true.
            temp = tempz
            exit newloop
         elseif (converged) then
            temp = 5.d-1 * (tempa+tempz)
            exit newloop
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         !----- Set funz, and check whether funa and funz already have opposite sign. -----!
         funz  = funnow
         zside = .false.
         if (funa*funnow > 0.d0) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler8 * tempa) then
               delta = 1.d2 * toler8 * tempa
            else
               delta = max(abs(funa*(tempz-tempa)/(funz-funa)), 1.d2 * toler8 * tempa)
            end if
            tempz = tempa + delta
            zgssloop: do itb=1,maxfpo
               !----- So this will be +1 -1 +2 -2 etc. ------------------------------------!
               tempz = tempz + dble((-1)**itb * (itb+3)/2) * delta
               theta = tempz * exnernormi
               rsat  = rslif8(pres,tempz,brrr_cold)
               funz  = thetaeivs8(theta,tempz,rsat,0.d0,0.d0) - theivs
               zside = funa*funz < 0.d0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside)                                                               &
               call fatal_error('Failed finding the second guess for regula falsi'         &
                             ,'thetaes2temp','therm_lib.f90')
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            if (abs(funz-funa) < toler8 * tempa) then
               temp   = 5.d-1 * (tempa+tempz)
            else
               temp   = (funz*tempa-funa*tempz)/(funz-funa)
            end if
            theta  = temp * exnernormi
            rsat   = rslif8(pres,temp,brrr_cold)
            funnow = thetaeivs8(theta,temp,rsat,0.d0,0.d0) - theivs

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa) < toler8 * temp
            if (converged) then
               exit fpoloop
            elseif (funnow*funa < 0.d0) then 
               tempz = temp
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 5.d-1
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tempa = temp
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not. zside) funz = funz * 5.d-1
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
         end do fpoloop
      end if

      if (converged) then 
         !----- Compute theta and rsat with temp just for consistency ---------------------!
         theta = temp * exnernormi
         rsat  = rslif8(pres,temp,brrr_cold)
      else
         call fatal_error('Temperature didn''t converge, I gave up!'                       &
                       ,'thetaes2temp8','therm_lib8.f90')
      end if

      return
   end subroutine thetaeivs2temp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine finds the lifting condensation level given the ice-liquid          !
   ! potential temperature in Kelvin, temperature in Kelvin, the pressure in Pascal, and   !
   ! the mixing ratio in kg/kg. The output will give the LCL temperature and pressure, and !
   ! the thickness of the layer between the initial point and the LCL.                     !
   !                                                                                       !
   !    References:                                                                        !
   !    Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential      !
   !        temperature as a thermodynamic variable in deep atmospheric models. Mon. Wea.  !
   !        Rev., v. 109, 1094-1102. (TC81)                                                !
   !    Bolton, D., 1980: The computation of the equivalent potential temperature. Mon.    !
   !        Wea. Rev., v. 108, 1046-1053. (BO80)                                           !
   !                                                                                       !
   !    Some algebra was needed to find this equation, essentially combining (TC81-26) and !
   ! (TC81-27), and the conservation of total water (TC81-16). It assumes that the divi-   !
   ! sion between the three phases is already taken care of.                               !
   !    Iterative procedure is needed, and here we iterate looking for T(LCL). Theta_il    !
   ! can be rewritten in terms of T(LCL) only, and once we know this thetae_iv becomes     !
   ! straightforward. T(LCL) will be found using Newton's method, and in the unlikely      !
   ! event it fails,we will fall back to the modified regula falsi (Illinois method).      !
   !                                                                                       !
   ! Important remarks:                                                                    !
   ! 1. TLCL and PLCL are the actual TLCL and PLCL, so in case condensation exists, they   !
   !    will be larger than the actual temperature and pressure (because one would go down !
   !    to reach the equilibrium);                                                         !
   ! 2. DZLCL WILL BE SET TO ZERO in case the LCL is beneath the starting level. So in     !
   !    case you want to force TLCL <= TEMP and PLCL <= PRES, you can use this variable    !
   !    to run the saturation check afterwards. DON'T CHANGE PLCL and TLCL here, they will !
   !    be used for conversions between theta_il and thetae_iv as they are defined here.   !
   ! 3. In case you don't want ice, simply pass useice=.false.. Otherwise let the model    !
   !    decide by itself based on the LEVEL variable.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine lcl_il8(thil,pres,temp,rtot,rvpr,tlcl,plcl,dzlcl,useice)
      use consts_coms, only: cpog8, alvl8,alvi8,cp8,ep8,p008,rocp8,ttripoli8,t3ple8,t008
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)            :: thil   ! Ice liquid pot. temp. (*)   [      K]
      real(kind=8), intent(in)            :: pres   ! Pressure                    [     Pa]
      real(kind=8), intent(in)            :: temp   ! Temperature                 [      K]
      real(kind=8), intent(in)            :: rtot   ! Total mixing ratio          [  kg/kg]
      real(kind=8), intent(in)            :: rvpr   ! Vapour mixing ratio         [  kg/kg]
      real(kind=8), intent(out)           :: tlcl   ! LCL temperature             [      K]
      real(kind=8), intent(out)           :: plcl   ! LCL pressure                [     Pa]
      real(kind=8), intent(out)           :: dzlcl  ! Sub-LCL layer thickness     [      m]
      logical     , intent(in) , optional :: useice ! Ice thermodynamics?         [    T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8) :: pvap       ! Sat. vapour pressure
      real(kind=8) :: deriv      ! Function derivative 
      real(kind=8) :: funnow     ! Function for which we seek a root.
      real(kind=8) :: funa       ! Smallest  guess function
      real(kind=8) :: funz       ! Largest   guess function
      real(kind=8) :: tlcla      ! Smallest  guess (or previous guess in Newton)
      real(kind=8) :: tlclz      ! Largest   guess (or new guess in Newton)
      real(kind=8) :: es00       ! Defined as p00*rt/(epsilon + rt)
      real(kind=8) :: delta      ! Aux. variable in case bisection is needed.
      integer      :: itn,itb    ! Iteration counters
      logical      :: converged  ! Convergence handle
      logical      :: zside      ! Aux. flag, to check sides for Regula Falsi
      !----- Other local variables --------------------------------------------------------!
      logical      :: brrr_cold ! This requires ice thermodynamics         [    T|F]
      !------------------------------------------------------------------------------------!
      ! (*) This is the most general variable. Thil is exactly theta for no condensation   !
      !     condition, and it is the liquid potential temperature if no ice is present.    !
      !------------------------------------------------------------------------------------!

      if (present(useice)) then
         brrr_cold = useice
      else 
         brrr_cold = bulk_on
      end if

      !----- Finding es00, which is a constant --------------------------------------------!
      es00 = p008 * rtot / (ep8 + rtot)


      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvpr / (ep8 + rvpr)

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a)') '----------------------------------------------------------'
      !write (unit=21,fmt='(a,1x,i5,1x,5(a,1x,f11.4,1x))')                                  &
      !   'INPUT : it=',-1,'thil=',thil,'pres=',0.01*pres,'temp=',temp-t00                  &
      !        ,'rvpr=',rvpr*1000.,'rtot=',rtot*1000.
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      tlclz     = 5.5d1 + 2.840d3 / (3.5d0 * log(temp) - log(1.d-2*pvap) - 4.805d0)

      pvap      = eslif8(tlclz,brrr_cold)

      funnow    = tlclz * (es00/pvap)**rocp8 - thil

      deriv     = (funnow+thil)*(1.d0/tlclz - rocp8*eslifp8(tlclz,brrr_cold)/pvap) 

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')                &
      !   'NEWTON: it=',0,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      tlcla     = tlclz
      funa      = funnow
      !------------------------------------------------------------------------------------!
      !     Looping: Newton's method.                                                      !
      !------------------------------------------------------------------------------------!

      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, skip to bisection ----!
         !----- Updating guesses ----------------------------------------------------------!
         tlcla  = tlclz
         funa   = funnow
         
         tlclz  = tlcla - funnow/deriv
         
         pvap   = eslif8(tlclz,brrr_cold)
         funnow = tlclz * (es00/pvap)**rocp8 - thil
         deriv  = (funnow+thil)*(1.d0/tlclz - rocp8*eslifp8(tlclz,brrr_cold)/pvap)

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
         !   'NEWTON: it=',itn,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow           &
         !          ,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !------------------------------------------------------------------------------!
         !   Convergence may happen when we get close guesses.                          !
         !------------------------------------------------------------------------------!
         converged = abs(tlcla-tlclz) < toler8 * tlclz
         if (converged) then
            tlcl = 5.d-1*(tlcla+tlclz)
            funz = funnow
            exit newloop
         elseif (funnow == 0.d0) then
            tlcl = tlclz
            funz = funnow
            converged = .true.
            exit newloop
         end if
      end do newloop

      if (.not. converged) then
         !---------------------------------------------------------------------------------!
         !     If I reached this point then it's because Newton's method failed. Using re- !
         ! gula falsi instead. First, I need to find two guesses that give me functions    !
         ! with opposite signs. If funa and funnow have opposite signs, then we are all    !
         ! set.                                                                            !
         !---------------------------------------------------------------------------------!
         if (funa*funnow < 0.d0 ) then
            funz  = funnow
            zside = .true.
         !----- They have the same sign, seeking the other guess --------------------------!
         else

            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funnow-funa) < toler8 * tlcla) then
               delta = 1.d2 * toler8 * tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funnow-funa)),1.d2 * toler8 * tlcla)
            end if
            tlclz = tlcla + delta

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=21,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')          &
            !   '2NGGSS: tt=',0,'tlclz=',tlclz-t00,'tlcla=',tlcla-t00,'pvap=',0.01*pvap     &
            !           ,'funa=',funa,'funz=',funnow,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            zside = .false.
            zgssloop: do itb=1,maxfpo

               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlcla + dble((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif8(tlclz,brrr_cold)
               funz  = tlclz * (es00/pvap)**rocp8 - thil

               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')       &
               !   '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz       &
               !           ,'delta=',delta
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               zside = funa*funz < 0.d0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(a)') ' ====== No second guess for you... ======'
               write (unit=*,fmt='(a)') ' + INPUT variables: '
               write (unit=*,fmt='(a,1x,es14.7)') 'THIL =',thil
               write (unit=*,fmt='(a,1x,es14.7)') 'TEMP =',temp
               write (unit=*,fmt='(a,1x,es14.7)') 'PRES =',pres
               write (unit=*,fmt='(a,1x,es14.7)') 'RTOT =',rtot
               write (unit=*,fmt='(a,1x,es14.7)') 'RVPR =',rvpr
               write (unit=*,fmt='(a)') ' ============ Failed guess... ==========='
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLA =',tlcla,'FUNA =',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLZ =',tlclz,'FUNC =',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'DELTA =',delta,'FUNN =',funnow
               call fatal_error('Failed finding the second guess for regula falsi'         &
                             ,'lcl_il8','therm_lib8.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo

            tlcl = (funz*tlcla-funa*tlclz)/(funz-funa)

            pvap = eslif8(tlcl,brrr_cold)

            funnow = tlcl * (es00/pvap)**rocp8 - thil

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),1(a,1x,es11.4,1x))')          &
            !   'REGFAL: it=',itb,'tlcl =',tlcl -t00,'pvap=',0.01*pvap,'fun=',funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla) < toler8*tlcl .and.  abs(tlcl-tlclz) < toler8*tlcl
            if (funnow == 0.d0 .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.d0) then 
               tlclz = tlcl
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 5.d-1
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tlcla = tlcl
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 5.d-1
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
         end do fpoloop
      end if
      !----- Finding the other LCL thermodynamic variables --------------------------------!
      if (converged) then 
         pvap  = eslif8(tlcl,brrr_cold)
         plcl  = (ep8 + rvpr) * pvap / rvpr
         dzlcl = max(cpog8*(temp-tlcl),0.d0)
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=21,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')             &
         !   'ANSWER: itb=',itn,'tlcl=',tlcl-t00,'eslcl=',0.01*pvap                         &
         !        ,'dzlcl=',dzlcl,'plcl=',plcl*0.01,'funa=',funa,'funz=',funz
         !write (unit=21,fmt='(a)') '-------------------------------------------------------'
         !write (unit=21,fmt='(a)') ' '
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      else
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' LCL Temperature didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'theta_il        [     K] =',thil
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Pressure        [   hPa] =',1.d-2*pres
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Temperature     [    °C] =',temp-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rtot            [  g/kg] =',1.d3*rtot
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rvpr            [  g/kg] =',10.d3*rvpr
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome.'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcla           [    °C] =',tlcla-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlclz           [    °C] =',tlclz-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [     K] =',funnow
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [     K] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [     K] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [  ----] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [  ----] =',toler8
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                                            ,abs(tlclz-tlcla)/tlclz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcl            [    °C] =',tlcl
         call fatal_error('TLCL didn''t converge, gave up!','lcl_il8','therm_lib8.f90')
      end if
      return
   end subroutine lcl_il8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! internal energy .   This requires double precision arguments.                         !
   !---------------------------------------------------------------------------------------!
   subroutine qtk8(q,tempk,fracliq)
      use consts_coms, only: cliqi8,cicei8,allii8,t3ple8,qicet38,qliqt38,tsupercool8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: q        ! Internal energy                     [   J/kg]
      real(kind=8), intent(out) :: tempk    ! Temperature                         [      K]
      real(kind=8), intent(out) :: fracliq  ! Liquid Fraction (0-1)               [    ---]
      !------------------------------------------------------------------------------------!


      !----- Internal energy below qwfroz, all ice  ---------------------------------------!
      if (q <= dble(qicet38)) then
         fracliq = 0.d0
         tempk   = q * cicei8 
      !----- Internal energy, above qwmelt, all liquid ------------------------------------!
      elseif (q >= dble(qliqt38)) then
         fracliq = 1.d0
         tempk   = q * cliqi8 + tsupercool8
      !----- Changing phase, it must be at freezing point ---------------------------------!
      else
         fracliq = (q-qicet38) * allii8
         tempk   = t3ple8
      endif
      !------------------------------------------------------------------------------------!

      return
   end subroutine qtk8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from inter-  !
   ! nal energy (J/m² or J/m³), mass (kg/m² or kg/m³), and heat capacity (J/m²/K or        !
   ! J/m³/K).                                                                              !
   ! This routine requires an 8-byte double precision floating point value for density.    !
   !---------------------------------------------------------------------------------------!
   subroutine qwtk8(qw,w,dryhcap,tempk,fracliq)
      use consts_coms, only: cliqi8,cliq8,cicei8,cice8,allii8,alli8,t3ple8,tsupercool8
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: qw      ! Internal energy           [  J/m²] or [  J/m³]
      real(kind=8), intent(in)  :: w       ! Density                   [ kg/m²] or [ kg/m³]
      real(kind=8), intent(in)  :: dryhcap ! Heat capacity, nonwater   [J/m²/K] or [J/m³/K]
      real(kind=8), intent(out) :: tempk   ! Temperature                           [     K]
      real(kind=8), intent(out) :: fracliq ! Liquid fraction (0-1)                 [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=8)              :: qwfroz  ! qw of ice at triple pt.   [  J/m²] or [  J/m³] 
      real(kind=8)              :: qwmelt  ! qw of liquid at triple pt.[  J/m²] or [  J/m³]
      !------------------------------------------------------------------------------------!

      !----- Converting melting heat to J/m² or J/m³ --------------------------------------!
      qwfroz = (dryhcap + w*cice8) * t3ple8
      qwmelt = qwfroz   + w*alli8
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    This is analogous to the qtk computation, we should analyse the magnitude of    !
      ! the internal energy to choose between liquid, ice, or both by comparing with our.  !
      ! know boundaries.                                                                   !
      !------------------------------------------------------------------------------------!

      !----- Internal energy below qwfroz, all ice  ---------------------------------------!
      if (qw < qwfroz) then
         fracliq = 0.d0
         tempk   = qw  / (cice8 * w + dryhcap)
      !----- Internal energy, above qwmelt, all liquid ------------------------------------!
      elseif (qw > qwmelt) then
         fracliq = 1.d0
         tempk   = (qw + w * cliq8 * tsupercool8) / (dryhcap + w*cliq8)
      !------------------------------------------------------------------------------------!
      !    We are at the freezing point.  If water mass is so tiny that the internal       !
      ! energy of frozen and melted states are the same given the machine precision, then  !
      ! we assume that water content is negligible and we impose 50% frozen for            !
      ! simplicity.                                                                        !
      !------------------------------------------------------------------------------------!
      elseif (qwfroz == qwmelt) then
         fracliq = 5.d-1
         tempk   = t3ple8
      !------------------------------------------------------------------------------------!
      !    Changing phase, it must be at freezing point.  The max and min are here just to !
      ! avoid tiny deviations beyond 0. and 1. due to floating point arithmetics.          !
      !------------------------------------------------------------------------------------!
      else
         fracliq = min(1.d0,max(0.d0,(qw - qwfroz) * allii8 / w))
         tempk = t3ple8
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine qwtk8
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib8

