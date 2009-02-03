!==========================================================================================!
!. File: therm_lib.f90                                                                            !
!                                                                                          !
!  Based on BRAMS-4.0.6   This file contains most functions and subroutines that deal with !
! several thermodynamic conversions. These procedures were built to avoid assumptions like !
! hydrostatic and linearisation. Most equations could not be solved analytically, and the  !
! standard here was to use Newton's method as the default, always having bisection or,     !
! more often, the modified Regula Falsi (Illinois) method in case Newton's fails.          !
!==========================================================================================!
!==========================================================================================!
module therm_lib
   implicit none

   !---------------------------------------------------------------------------------------!
   ! Constants that control the convergence for iterative methods                          !
   !---------------------------------------------------------------------------------------!
   real   , parameter ::   toler  = 2.* epsilon(1.) ! Relative tolerance for iterative 
                                                    !   methods. The smaller the value, the
                                                    !   more accurate the result, but it
                                                    !   will slow down the run.
   integer, parameter ::   maxfpo = 60              ! Maximum # of iterations before crash-
                                                    !   ing for false position method.

   integer, parameter ::   maxit  = 150             ! Maximum # of iterations before crash-
                                                    !   ing, for other methods.

   integer, parameter ::   maxlev = 16              ! Maximum # of levels for adaptive     
                                                    !   quadrature methods.    

   logical, parameter ::   newthermo = .true.      ! Use new thermodynamics [T|F]

   !---------------------------------------------------------------------------------------!
   !   This is the "level" variable, that used to be in micphys. Since it affects more the !
   ! thermodynamics choices than the microphysics, it was moved to here.                   !
   !---------------------------------------------------------------------------------------!
   integer, parameter ::   level = 3

   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical, parameter ::   vapour_on = .true.
   logical, parameter ::   cloud_on  = .true.
   logical, parameter ::   bulk_on   = .true.
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
   real, dimension(0:3), parameter :: iii_7 = (/ 9.550426,-5723.265, 3.53068,-0.00728332 /)
   !----- Coefficients based on equation (10), first fit ----------------------------------!
   real, dimension(0:3), parameter :: l01_10= (/54.842763,-6763.22 ,-4.210  , 0.000367   /)
   !----- Coefficients based on equation (10), second fit ---------------------------------!
   real, dimension(0:3), parameter :: l02_10= (/53.878   ,-1331.22 ,-9.44523, 0.014025   /)
   !----- Coefficients based on the hyperbolic tangent ------------------------------------!
   real, dimension(2)  , parameter :: ttt_10= (/0.0415,218.8/)
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
   real, dimension(0:8), parameter :: cll = (/ .6105851e+03,  .4440316e+02,  .1430341e+01  &
                                             , .2641412e-01,  .2995057e-03,  .2031998e-05  &
                                             , .6936113e-08,  .2564861e-11, -.3704404e-13 /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real, dimension(0:8), parameter :: cii = (/ .6114327e+03,  .5027041e+02,  .1875982e+01  &
                                             , .4158303e-01,  .5992408e-03,  .5743775e-05  &
                                             , .3566847e-07,  .1306802e-09,  .2152144e-12 /)
   !----- Coefficients for d(esat)/dT (liquid) --------------------------------------------!
   real, dimension(0:8), parameter :: dll = (/ .4443216e+02,  .2861503e+01,  .7943347e-01  &
                                             , .1209650e-02,  .1036937e-04,  .4058663e-07  &
                                             ,-.5805342e-10, -.1159088e-11, -.3189651e-14 /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real, dimension(0:8), parameter :: dii = (/ .5036342e+02,  .3775758e+01,  .1269736e+00  &
                                             , .2503052e-02,  .3163761e-04,  .2623881e-06  &
                                             , .1392546e-08,  .4315126e-11,  .5961476e-14 /)
   !---------------------------------------------------------------------------------------!



   contains
   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour pressure as a function of   !
   ! Kelvin temperature. This expression came from MK05, equation (10).                    !
   !---------------------------------------------------------------------------------------!
   real function eslf(temp,l1funout,l2funout,ttfunout)
      use consts_coms, only : t00
      implicit none
      real, intent(in)            :: temp
      real, intent(out), optional :: l1funout,ttfunout,l2funout
      real                        :: l1fun,ttfun,l2fun,x

      if(temp <= 0.0) then
         print*,"temp < 0K in eslf",temp
         stop
      endif
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         l1fun = l01_10(0) + l01_10(1)/temp + l01_10(2)*log(temp) + l01_10(3) * temp
         l2fun = l02_10(0) + l02_10(1)/temp + l02_10(2)*log(temp) + l02_10(3) * temp
         ttfun = tanh(ttt_10(1) * (temp - ttt_10(2)))
         eslf  = exp(l1fun + ttfun*l2fun)

         if (present(l1funout)) l1funout = l1fun
         if (present(l2funout)) l2funout = l2fun
         if (present(ttfunout)) ttfunout = ttfun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x    = max(-80.,temp-t00)
         eslf = cll(0) + x * (cll(1) + x * (cll(2) + x * (cll(3) + x * (cll(4)             &
                       + x * (cll(5) + x * (cll(6) + x * (cll(7) + x * cll(8)) ) ) ) ) ) )

         if (present(l1funout)) l1funout = eslf
         if (present(l2funout)) l2funout = eslf
         if (present(ttfunout)) ttfunout = eslf
      end if

      return
   end function eslf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour pressure as a function of      !
   ! Kelvin temperature, based on MK05 equation (7).                                       !
   !---------------------------------------------------------------------------------------!
   real function esif(temp,iifunout)
      use consts_coms, only : t00
      implicit none
      real, intent(in)            :: temp
      real, intent(out), optional :: iifunout
      real                        :: iifun,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         iifun = iii_7(0) + iii_7(1)/temp + iii_7(2) * log(temp) + iii_7(3) * temp
         esif  = exp(iifun)
      
         if (present(iifunout)) iifunout=iifun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-80.,temp-t00)
         esif = cii(0) + x * (cii(1) + x * (cii(2) + x * (cii(3) + x * (cii(4)             &
                       + x * (cii(5) + x * (cii(6) + x * (cii(7) + x * cii(8)) ) ) ) ) ) )

         if (present(iifunout)) iifunout=esif
      end if
      return
   end function esif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation vapour pressure as a function of  Kelvin  !
   ! temperature. It chooses which phase to look depending on whether the temperature is   !
   ! below or above the triple point.                                                      !
   !---------------------------------------------------------------------------------------!
   real function eslif(temp,useice)
      use consts_coms, only: t3ple
      implicit none
      real   , intent(in)           :: temp
      logical, intent(in), optional :: useice
      logical                       :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple
      else 
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      if (brrr_cold) then
         eslif = esif(temp) ! Ice saturation vapour pressure 
      else
         eslif = eslf(temp) ! Liquid saturation vapour pressure
      end if

      return
   end function eslif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour mixing ratio as a function  !
   ! of pressure and Kelvin temperature.                                                   !
   !---------------------------------------------------------------------------------------!
   real function rslf(pres,temp)
      use consts_coms, only : ep,toodry
      implicit none
      real, intent(in) :: pres,temp
      real             :: esl

      esl  = eslf(temp)
      rslf = max(toodry,ep*esl/(pres-esl))

      return
   end function rslf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour mixing ratio as a function of  !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real function rsif(pres,temp)
      use consts_coms, only : ep,toodry
      implicit none
      real, intent(in) :: pres,temp
      real             :: esi

      esi  = esif(temp)
      rsif = max(toodry,ep*esi/(pres-esi))

      return
   end function rsif
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation vapour mixing ratio, over liquid or ice   !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real function rslif(pres,temp,useice)
      use consts_coms, only: t3ple,ep
      implicit none
      real   , intent(in)           :: pres,temp
      logical, intent(in), optional :: useice
      real                          :: esz
    
      if (present(useice)) then
         esz = eslif(temp,useice)
      else
         esz = eslif(temp)
      end if

      rslif = ep * esz / (pres - esz)

      return
   end function rslif
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-liquid equilibrium density for vapour, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real function rhovsl(temp)
      use consts_coms, only : rvap
      implicit none
      real, intent(in) :: temp
      real             :: eequ
      eequ = eslf(temp)
      rhovsl = eequ / (rvap * temp)
      return
   end function rhovsl
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-ice equilibrium density for vapour, as a      !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real function rhovsi(temp)
      use consts_coms, only : rvap
      implicit none
      real, intent(in) :: temp
      real             :: eequ
      eequ = esif(temp)
      rhovsi = eequ / (rvap * temp)
      return
   end function rhovsi
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation density for vapour, as a function of tem- !
   ! perature in Kelvin. It will decide between ice-vapour or liquid-vapour based on the   !
   ! temperature.                                                                          !
   !---------------------------------------------------------------------------------------!
   real function rhovsil(temp,useice)
      use consts_coms, only : rvap
      implicit none
      real, intent(in)              :: temp
      logical, intent(in), optional :: useice
      real                          :: eequ

      if (present(useice)) then
         eequ = eslif(temp,useice)
      else
         eequ = eslif(temp)
      end if

      rhovsil = eequ / (rvap * temp)

      return
   end function rhovsil
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour       !
   ! pressure with respect to temperature as a function of Kelvin temperature.             !
   !---------------------------------------------------------------------------------------!
   real function eslfp(temp)
      use consts_coms, only: t00
      implicit none
      real, intent(in) :: temp
      real             :: esl,l2fun,ttfun,l1prime,l2prime,ttprime,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esl   = eslf(temp,l2funout=l2fun,ttfunout=ttfun)
         l1prime = -l01_10(1)/(temp*temp) + l01_10(2)/temp + l01_10(3)
         l2prime = -l02_10(1)/(temp*temp) + l02_10(2)/temp + l02_10(3)
         ttprime =  ttt_10(1)*(1.-ttfun*ttfun)
         eslfp = esl * (l1prime + l2prime*ttfun + l2fun*ttprime)
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-80.,temp-t00)
         eslfp = dll(0) + x * (dll(1) + x * (dll(2) + x * (dll(3) + x * (dll(4)            &
                        + x * (dll(5) + x * (dll(6) + x * (dll(7) + x * dll(8)) ) ) ) ) ) )
      end if


      return
   end function eslfp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of ice saturation vapour pressure !
   ! with respect to temperature as a function of Kelvin temperature.                      !
   !---------------------------------------------------------------------------------------!
   real function esifp(temp)
      use consts_coms, only: t00
      implicit none
      real, intent(in) :: temp
      real             :: esi,iiprime,x

      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esi     = esif(temp)
         iiprime = -iii_7(1)/(temp*temp) + iii_7(2)/temp + iii_7(3)
         esifp   = esi * iiprime
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-80.,temp-t00)
         esifp = dii(0) + x * (dii(1) + x * (dii(2) + x * (dii(3) + x * (dii(4)            &
                        + x * (dii(5) + x * (dii(6) + x * (dii(7) + x * dii(8)) ) ) ) ) ) )
      end if

      return
   end function esifp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of saturation vapour pressure as  !
   ! a function of  Kelvin temperature. It chooses which phase to look depending on        !
   ! whether the temperature is below or above the triple point.                           !
   !---------------------------------------------------------------------------------------!
   real function eslifp(temp,useice)
      use consts_coms, only: t3ple
      implicit none
      real   , intent(in)           :: temp
      logical, intent(in), optional :: useice
      logical                       :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple
      else 
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      if (brrr_cold) then
         eslifp = esifp(temp) ! d(Ice saturation vapour pressure)/dT
      else
         eslifp = eslfp(temp) ! d(Liquid saturation vapour pressure)/dT
      end if

      return
   end function eslifp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour mix-  !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real function rslfp(pres,temp)
      use consts_coms, only: ep
      implicit none
      real, intent(in) :: pres,temp
      real             :: desdt,esl,pdry

      esl    = eslf(temp)
      desdt  = eslfp(temp)
      
      pdry   = pres-esl
      rslfp  = ep * pres * desdt / (pdry*pdry)

      return
   end function rslfp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of ice saturation vapour mix-     !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real function rsifp(pres,temp)
      use consts_coms, only: ep
      implicit none
      real, intent(in) :: pres,temp
      real             :: desdt,esi,pdry

      esi    = esif(temp)
      desdt  = esifp(temp)
      
      pdry   = pres-esi
      rsifp  = ep * pres * desdt / (pdry*pdry)

      return
   end function rsifp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour mix-  !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real function rslifp(pres,temp,useice)
      use consts_coms, only: t3ple
      implicit none
      real   , intent(in)           :: pres,temp
      logical, intent(in), optional :: useice
      logical                       :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple
      else 
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      if (brrr_cold) then
         rslifp=rsifp(pres,temp)
      else
         rslifp=rslfp(pres,temp)
      end if

      return
   end function rslifp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of vapour-liquid equilibrium density, as  !
   ! a function of temperature in Kelvin.                                                  !
   !---------------------------------------------------------------------------------------!
   real function rhovslp(temp)
      use consts_coms, only : rvap
      implicit none
      real, intent(in) :: temp
      real             :: es,desdt
      es    = eslf(temp)
      desdt = eslfp(temp)
      rhovslp = (desdt-es/temp) / (rvap * temp)
      return
   end function rhovslp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of vapour-ice equilibrium density, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real function rhovsip(temp)
      use consts_coms, only : rvap
      implicit none
      real, intent(in) :: temp
      real             :: es,desdt
      es    = esif(temp)
      desdt = esifp(temp)
      rhovsip = (desdt-es/temp) / (rvap * temp)
      return
   end function rhovsip
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of saturation density for vapour, as a    !
   ! function of temperature in Kelvin. It will decide between ice-vapour or liquid-vapour !
   ! based on the temperature.                                                             !
   !---------------------------------------------------------------------------------------!
   real function rhovsilp(temp,useice)
      use consts_coms, only: t3ple
      implicit none
      real, intent(in)              :: temp
      logical, intent(in), optional :: useice
      logical                       :: brrr_cold

      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple
      else 
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      if (brrr_cold) then
         rhovsilp=rhovsip(temp)
      else
         rhovsilp=rhovslp(temp)
      end if

      return
   end function rhovsilp
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
   real function virtt(temp,rvpr,rtot)
      use consts_coms, only: epi
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: temp     ! Temperature                          [    K]
      real, intent(in)           :: rvpr     ! Vapour mixing ratio                  [kg/kg]
      real, intent(in), optional :: rtot     ! Total mixing ratio                   [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real                       :: rtothere ! Internal rtot, to deal with optional [kg/kg]

      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvpr
      end if

      virtt = temp * (1. + epi * rvpr) / (1. + rtothere)

      return
   end function virtt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates the liquid saturation mixing ratio as a function of    !
   ! pressure and Kelvin temperature for a column.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine mrsl(n1,pres,temp,rsl)
      implicit none
      integer            ,intent(in)   :: n1
      real, dimension(n1), intent(in)  :: pres,temp
      real, dimension(n1), intent(out) :: rsl
      integer                          :: n

      do n=1,n1
         rsl(n)=rslf(pres(n),temp(n))
      end do

      return
   end subroutine mrsl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates the ice saturation mixing ratio as a function of pres- !
   ! sure and Kelvin temperature for a column.                                             !
   !---------------------------------------------------------------------------------------!
   subroutine mrsi(n1,pres,temp,rsi)
      implicit none
      integer            ,intent(in)   :: n1
      real, dimension(n1), intent(in)  :: pres,temp
      real, dimension(n1), intent(out) :: rsi
      integer                          :: n

      do n=1,n1
         rsi(n)=rsif(pres(n),temp(n))
      end do

      return
   end subroutine mrsi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates the saturation mixing ratio as a function of pressure  !
   ! and Kelvin temperature for a column.                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine mrsli(n1,pres,temp,rsi,useice)
      implicit none
      integer            ,intent(in)            :: n1
      real, dimension(n1), intent(in)           :: pres,temp
      real, dimension(n1), intent(out)          :: rsi
      logical            , intent(in), optional :: useice
      integer                                   :: n

      if (present(useice)) then
         do n=1,n1
            rsi(n)=rslif(pres(n),temp(n),useice)
         end do

      else
         do n=1,n1
            rsi(n)=rslif(pres(n),temp(n))
         end do

      end if

      return
   end subroutine mrsli
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! internal energy .                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine qtk(q,tempk,fracliq)
      use consts_coms, only: cliqi,cicei,allii,alli,tsupercool,t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: q        ! Internal energy                             [   J/kg]
      real, intent(out) :: tempk    ! Temperature                                 [      K]
      real, intent(out) :: fracliq  ! Liquid Fraction (0-1)                       [    ---]
      !------------------------------------------------------------------------------------!


      !----- Negative internal energy, frozen, all ice ------------------------------------!
      if (q <= 0.) then
         fracliq = 0.
         tempk   = q * cicei + t3ple
      !----- Positive internal energy, over latent heat of melting, all liquid ------------!
      elseif (q > alli) then
         fracliq = 1.
         tempk   = q * cliqi + tsupercool
      !----- Changing phase, it must be at triple point -----------------------------------!
      else
         fracliq = q * allii
         tempk   = t3ple
      endif
      !------------------------------------------------------------------------------------!

      return
   end subroutine qtk
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature, and the fraction of liquid water from    !
   ! the internal energy. The only difference from qtk (which is actually called here) is  !
   ! that the temperature is returned in Celsius.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine qtc(q,tempc,fracliq)
      use consts_coms, only: t00
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: q        ! Internal energy                             [   J/kg]
      real, intent(out) :: tempc    ! Temperature                                 [      K]
      real, intent(out) :: fracliq  ! Liquid Fraction (0-1)                       [    ---]
      !------------------------------------------------------------------------------------!

      call qtk(q,tempc,fracliq)
      tempc = tempc - t00

      return
   end subroutine qtc
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from inter-  !
   ! nal energy (J/m² or J/m³), mass (kg/m² or kg/m³), and heat capacity (J/m²/K or        !
   ! J/m³/K).                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine qwtk(qw,w,dryhcap,tempk,fracliq)
      use consts_coms, only: cliqi,cliq,cicei,cice,allii,alli,t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: qw      ! Internal energy                   [  J/m²] or [  J/m³]
      real, intent(in)  :: w       ! Density                           [ kg/m²] or [ kg/m³]
      real, intent(in)  :: dryhcap ! Heat capacity of nonwater part    [J/m²/K] or [J/m³/K]
      real, intent(out) :: tempk   ! Temperature                                   [     K]
      real, intent(out) :: fracliq ! Liquid fraction (0-1)                         [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real              :: qwliq0  ! qw of liquid at triple point      [  J/m²] or [  J/m³]
      !------------------------------------------------------------------------------------!

      !----- Converting melting heat to J/m² or J/m³ --------------------------------------!
      qwliq0 = w * alli
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    This is analogous to the qtk computation, we should analyse the sign and        !
      ! magnitude of the internal energy to choose between liquid, ice, or both.           !
      !------------------------------------------------------------------------------------!

      !----- Negative internal energy, frozen, all ice ------------------------------------!
      if (qw < 0.) then
         fracliq = 0.
         tempk   = qw  / (cice * w + dryhcap) + t3ple
      !----- Positive internal energy, over latent heat of melting, all liquid ------------!
      elseif (qw > qwliq0) then
         fracliq = 1.
         tempk   = (qw - qwliq0) / (cliq * w + dryhcap) + t3ple
      !----- Changing phase, it must be at triple point -----------------------------------!
      else
         if(w > 0.) then
            fracliq = qw / qwliq0
         else
            fracliq = 0.0
         endif
         tempk = t3ple
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine qwtk
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from inter-  !
   ! nal energy (J/m² or J/m³), mass (kg/m² or kg/m³), and heat capacity (J/m²/K or        !
   ! J/m³/K).                                                                              !
   ! This routine requires an 8-byte double precision floating point value for density     !
   !---------------------------------------------------------------------------------------!
   subroutine qwtk8(qw,w,dryhcap,tempk,fracliq)
      use consts_coms, only: cliqi,cliq,cicei,cice,allii,alli,t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real        , intent(in)  :: qw      ! Internal energy           [  J/m²] or [  J/m³]
      real(kind=8), intent(in)  :: w       ! Density                   [ kg/m²] or [ kg/m³]
      real        , intent(in)  :: dryhcap ! Heat capacity, nonwater   [J/m²/K] or [J/m³/K]
      real        , intent(out) :: tempk   ! Temperature                           [     K]
      real        , intent(out) :: fracliq ! Liquid fraction (0-1)                 [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real              :: qwliq0  ! qw of liquid at triple point      [  J/m²] or [  J/m³]
      !------------------------------------------------------------------------------------!
     
      !----- Converting melting heat to J/m² or J/m³ --------------------------------------!
      qwliq0 = sngl(w) * alli
      !------------------------------------------------------------------------------------!
     
      !------------------------------------------------------------------------------------!
      !    This is analogous to the qtk computation, we should analyse the sign and        !
      ! magnitude of the internal energy to choose between liquid, ice, or both.           !
      !------------------------------------------------------------------------------------!
     
      !----- Negative internal energy, frozen, all ice ------------------------------------!
      if (qw < 0.) then
         fracliq = 0.
         tempk   = qw  / (cice * sngl(w) + dryhcap) + t3ple
      !----- Positive internal energy, over latent heat of melting, all liquid ------------!
      elseif (qw > qwliq0) then
         fracliq = 1.
         tempk   = (qw - qwliq0) / (cliq * sngl(w) + dryhcap) + t3ple
      !----- Changing phase, it must be at triple point -----------------------------------!
      else
         if(w > 0.d0) then
            fracliq = qw / qwliq0
         else
            fracliq = 0.0
         endif
         tempk = t3ple
      end if
      !------------------------------------------------------------------------------------!
     
      return
   end subroutine qwtk8

 end module therm_lib

