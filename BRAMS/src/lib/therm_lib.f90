!==========================================================================================!
! BRAMS-4.0.6. File: therm_lib.f90                                                                            !
!                                                                                          !
!    This file contains most functions and subroutines that deal with several thermo-      !
! dynamic conversions. These procedures were built to avoid assumptions like hydrostatic   !
! and linearisation. Most equations could not be solved analytically, and the standard     !
! here was to use Newton's method as the default, always having bisection or, more often,  !
! the modified Regula Falsi (Illinois) method in case Newton's fails.                      !
!==========================================================================================!
!==========================================================================================!
module therm_lib
   implicit none

   !---------------------------------------------------------------------------------------!
   ! Constants that control the convergence for iterative methods                          !
   !---------------------------------------------------------------------------------------!
   real   , parameter ::   toler  = 10.* epsilon(1.) ! Relative tolerance for iterative 
                                                    !   methods. The smaller the value, the
                                                    !   more accurate the result, but it
                                                    !   will slow down the run.
   integer, parameter ::   maxfpo = 80              ! Maximum # of iterations before crash-
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
   integer            ::   level

   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical            ::   vapour_on
   logical            ::   cloud_on
   logical            ::   bulk_on
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
   real(kind=4), dimension(0:3), parameter :: iii_7   = (/  9.550426, -5723.265            &
                                                         ,  3.530680,    -0.00728332 /)
   !----- Coefficients based on equation (10), first fit ----------------------------------!
   real(kind=4), dimension(0:3), parameter :: l01_10  = (/ 54.842763, -6763.220            &
                                                         , -4.210   ,     0.000367   /)
   !----- Coefficients based on equation (10), second fit ---------------------------------!
   real(kind=4), dimension(0:3), parameter :: l02_10  = (/ 53.878   , -1331.22             &
                                                         , -9.44523 , 0.014025       /)
   !----- Coefficients based on the hyperbolic tangent ------------------------------------!
   real(kind=4), dimension(2)  , parameter :: ttt_10  = (/  0.0415  ,   218.80       /)
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
   real(kind=4), dimension(0:8), parameter :: cll = (/  .6105851e+03,  .4440316e+02        &
                                                     ,  .1430341e+01,  .2641412e-01        &
                                                     ,  .2995057e-03,  .2031998e-05        &
                                                     ,  .6936113e-08,  .2564861e-11        &
                                                     , -.3704404e-13                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: cii = (/  .6114327e+03,  .5027041e+02        &
                                                     ,  .1875982e+01,  .4158303e-01        &
                                                     ,  .5992408e-03,  .5743775e-05        &
                                                     ,  .3566847e-07,  .1306802e-09        &
                                                     ,  .2152144e-12                /)
   !----- Coefficients for d(esat)/dT (liquid) --------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: dll = (/  .4443216e+02,  .2861503e+01        &
                                                     ,  .7943347e-01,  .1209650e-02        &
                                                     ,  .1036937e-04,  .4058663e-07        &
                                                     , -.5805342e-10, -.1159088e-11        &
                                                     , -.3189651e-14                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: dii = (/  .5036342e+02,  .3775758e+01        &
                                                     ,  .1269736e+00,  .2503052e-02        &
                                                     ,  .3163761e-04,  .2623881e-06        &
                                                     ,  .1392546e-08,  .4315126e-11        &
                                                     ,  .5961476e-14                /)
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour pressure as a function of   !
   ! Kelvin temperature. This expression came from MK05, equation (10).                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function eslf(temp,l1funout,l2funout,ttfunout)
      use rconstants , only : t00 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)            :: temp     ! Temperature                [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=4), intent(out), optional :: l1funout ! Function for high temperatures
      real(kind=4), intent(out), optional :: ttfunout ! Interpolation function
      real(kind=4), intent(out), optional :: l2funout ! Function for low temperatures
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                        :: l1fun    ! 
      real(kind=4)                        :: ttfun    ! 
      real(kind=4)                        :: l2fun    ! 
      real(kind=4)                        :: x        ! 
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         l1fun = l01_10(0) + l01_10(1)/temp + l01_10(2)*log(temp) + l01_10(3) * temp
         l2fun = l02_10(0) + l02_10(1)/temp + l02_10(2)*log(temp) + l02_10(3) * temp
         ttfun = tanh(ttt_10(1) * (temp - ttt_10(2)))
         eslf  = exp(l1fun + ttfun*l2fun)
         !---------------------------------------------------------------------------------!

         if (present(l1funout)) l1funout = l1fun
         if (present(l2funout)) l2funout = l2fun
         if (present(ttfunout)) ttfunout = ttfun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x    = max(-80.,temp-t00)
         eslf = cll(0) + x * (cll(1) + x * (cll(2) + x * (cll(3) + x * (cll(4)             &
                       + x * (cll(5) + x * (cll(6) + x * (cll(7) + x * cll(8)) ) ) ) ) ) )
         !---------------------------------------------------------------------------------!

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
   real(kind=4) function esif(temp,iifunout)
      use rconstants , only : t00 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)            :: temp     ! Temperature                 [    K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=4), intent(out), optional :: iifunout
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                        :: iifun
      real(kind=4)                        :: x
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then

         !----- Updated method, using MK05 ------------------------------------------------!
         iifun = iii_7(0) + iii_7(1)/temp + iii_7(2) * log(temp) + iii_7(3) * temp
         esif  = exp(iifun)
         !---------------------------------------------------------------------------------!


         if (present(iifunout)) iifunout=iifun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-80.,temp-t00)
         esif = cii(0) + x * (cii(1) + x * (cii(2) + x * (cii(3) + x * (cii(4)             &
                       + x * (cii(5) + x * (cii(6) + x * (cii(7) + x * cii(8)) ) ) ) ) ) )
         !---------------------------------------------------------------------------------!

         if (present(iifunout)) iifunout=esif
      end if
      !------------------------------------------------------------------------------------!

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
   real(kind=4) function eslif(temp,useice)
      use rconstants , only : t3ple ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      logical                            :: frozen
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         eslif = esif(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         eslif = eslf(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function eslif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour mixing ratio as a function  !
   ! of pressure and Kelvin temperature.                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rslf(pres,temp)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: esl  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esl  = eslf(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the mixing !
      ! ratio.                                                                             !
      !------------------------------------------------------------------------------------!
      rslf = max(toodry,ep*esl/(pres-esl))
      !------------------------------------------------------------------------------------!

      return
   end function rslf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour mixing ratio as a function of  !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rsif(pres,temp)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: esi  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esi  = esif(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the mixing !
      ! ratio.                                                                             !
      !------------------------------------------------------------------------------------!
      rsif = max(toodry,ep*esi/(pres-esi))
      !------------------------------------------------------------------------------------!

      return
   end function rsif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation vapour mixing ratio, over liquid or ice   !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rslif(pres,temp,useice)
      use rconstants , only : t3ple & ! intent(in)
                            , ep    ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres
      real(kind=4), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: esz
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         esz = esif(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         esz = eslf(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the mixing !
      ! ratio.                                                                             !
      !------------------------------------------------------------------------------------!
      rslif = ep * esz / (pres - esz)
      !------------------------------------------------------------------------------------!

      return
   end function rslif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation specific humidity as a function of !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function qslf(pres,temp)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: esl  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esl  = eslf(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qslf = max(toodry,ep * esl/( pres - (1.0 - ep) * esl) )
      !------------------------------------------------------------------------------------!

      return
   end function qslf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation specific humidity as a function of    !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function qsif(pres,temp)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: esi  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esi  = esif(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qsif = max(toodry,ep * esi/( pres - (1.0 - ep) * esi) )
      !------------------------------------------------------------------------------------!

      return
   end function qsif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation specific humidity, over liquid or ice     !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function qslif(pres,temp,useice)
      use rconstants , only : t3ple  & ! intent(in)
                            , ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres
      real(kind=4), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: esz
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         esz = esif(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         esz = eslf(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qslif = max(toodry, ep * esz/( pres - (1.0 - ep) * esz) )
      !------------------------------------------------------------------------------------!

      return
   end function qslif
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-liquid equilibrium density for vapour, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rhovsl(temp)
      use rconstants , only : rh2o ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: eequ ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the equilibrium (saturation) vapour pressure.                             !
      !------------------------------------------------------------------------------------!
      eequ = eslf(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsl = eequ / (rh2o * temp)
      !------------------------------------------------------------------------------------!

      return
   end function rhovsl
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-ice equilibrium density for vapour, as a      !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rhovsi(temp)
      use rconstants , only : rh2o ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: eequ ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the equilibrium (saturation) vapour pressure.                             !
      !------------------------------------------------------------------------------------!
      eequ = esif(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsi = eequ / (rh2o * temp)
      !------------------------------------------------------------------------------------!

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
   real(kind=4) function rhovsil(temp,useice)
      use rconstants , only : rh2o ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: eequ
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Pass the "useice" argument to eslif, so it may decide whether ice thermo-      !
      ! dynamics is to be used.                                                            !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         eequ = eslif(temp,useice)
      else
         eequ = eslif(temp)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsil = eequ / (rh2o * temp)
      !------------------------------------------------------------------------------------!

      return
   end function rhovsil
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour       !
   ! pressure with respect to temperature as a function of Kelvin temperature.             !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function eslfp(temp)
      use rconstants , only : t00 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------ Local variables. ------------------------------------------------------------!
      real(kind=4)             :: esl
      real(kind=4)             :: l2fun
      real(kind=4)             :: ttfun
      real(kind=4)             :: l1prime
      real(kind=4)             :: l2prime
      real(kind=4)             :: ttprime
      real(kind=4)             :: x
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide which function to use, based on the thermodynamics.                    !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esl     = eslf(temp,l2funout=l2fun,ttfunout=ttfun)
         l1prime = -l01_10(1)/(temp*temp) + l01_10(2)/temp + l01_10(3)
         l2prime = -l02_10(1)/(temp*temp) + l02_10(2)/temp + l02_10(3)
         ttprime =  ttt_10(1)*(1.-ttfun*ttfun)
         eslfp   = esl * (l1prime + l2prime*ttfun + l2fun*ttprime)
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-80.,temp-t00)
         eslfp = dll(0) + x * (dll(1) + x * (dll(2) + x * (dll(3) + x * (dll(4)            &
                        + x * (dll(5) + x * (dll(6) + x * (dll(7) + x * dll(8)) ) ) ) ) ) )
      end if
      !------------------------------------------------------------------------------------!


      return
   end function eslfp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of ice saturation vapour pressure !
   ! with respect to temperature as a function of Kelvin temperature.                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function esifp(temp)
      use rconstants , only : t00 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------ Local variables. ------------------------------------------------------------!
      real(kind=4)             :: esi
      real(kind=4)             :: iiprime
      real(kind=4)             :: x
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide which function to use, based on the thermodynamics.                    !
      !------------------------------------------------------------------------------------!
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
      !------------------------------------------------------------------------------------!

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
   real(kind=4) function eslifp(temp,useice)
      use rconstants , only : t3ple ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in)           :: temp
      !------ Local variables. ------------------------------------------------------------!
      logical     , intent(in), optional :: useice
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- d(Saturation vapour pressure)/dT for ice. ---------------------------------!
         eslifp = esifp(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- d(Saturation vapour pressure)/dT for liquid water. ------------------------!
         eslifp = eslfp(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

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
   real(kind=4) function rslfp(pres,temp)
      use rconstants , only : ep ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres  ! Pressure                                 [    Pa]
      real(kind=4), intent(in) :: temp  ! Temperature                              [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: esl   ! Partial pressure                         [    Pa]
      real(kind=4)             :: desdt ! Derivative of partial pressure of water  [  Pa/K]
      real(kind=4)             :: pdry  ! Partial pressure of dry air              [    Pa]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      esl    = eslf(temp)
      desdt  = eslfp(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial pressure of dry air. ----------------------------------------!
      pdry   = pres-esl
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of mixing ratio. ---------------------------------!
      rslfp  = ep * pres * desdt / (pdry*pdry)
      !------------------------------------------------------------------------------------!

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
   real(kind=4) function rsifp(pres,temp)
      use rconstants , only : ep ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres  ! Pressure                                 [    Pa]
      real(kind=4), intent(in) :: temp  ! Temperature                              [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: esi   ! Partial pressure                         [    Pa]
      real(kind=4)             :: desdt ! Derivative of partial pressure of water  [  Pa/K]
      real(kind=4)             :: pdry  ! Partial pressure of dry air              [    Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      esi    = esif(temp)
      desdt  = esifp(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial pressure of dry air. ----------------------------------------!
      pdry   = pres-esi
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of mixing ratio. ---------------------------------!
      rsifp  = ep * pres * desdt / (pdry*pdry)
      !------------------------------------------------------------------------------------!
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
   real(kind=4) function rslifp(pres,temp,useice)
      use rconstants , only: t3ple ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres   ! Pressure                      [    Pa]
      real(kind=4), intent(in)           :: temp   ! Temperature                   [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! May use ice thermodynamics?   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: desdt  ! Derivative of vapour pressure [  Pa/K]
      logical                            :: frozen ! Use the ice thermodynamics    [   T|F]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide whether to use liquid water of ice for saturation, based on the temper- !
      ! ature and the settings.                                                            !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the function depending on the previous check.                             !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         rslifp = rsifp(pres,temp)
      else
         rslifp = rslfp(pres,temp)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function rslifp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of vapour-liquid equilibrium density, as  !
   ! a function of temperature in Kelvin.                                                  !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rhovslp(temp)
      use rconstants , only : rh2o ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in) :: temp   ! Temperature                             [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: es     ! Vapour pressure                         [    Pa]
      real(kind=4)             :: desdt  ! Vapour pressure derivative              [  Pa/K]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      es    = eslf(temp)
      desdt = eslfp(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of saturation density . --------------------------!
      rhovslp = (desdt-es/temp) / (rh2o * temp)
      !------------------------------------------------------------------------------------!

      return
   end function rhovslp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the derivative of vapour-ice equilibrium density, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rhovsip(temp)
      use rconstants , only : rh2o ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in) :: temp   ! Temperature                             [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: es     ! Vapour pressure                         [    Pa]
      real(kind=4)             :: desdt  ! Vapour pressure derivative              [  Pa/K]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      es    = esif(temp)
      desdt = esifp(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of saturation density . --------------------------!
      rhovsip = (desdt-es/temp) / (rh2o * temp)
      !------------------------------------------------------------------------------------!

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
   real(kind=4) function rhovsilp(temp,useice)
      use rconstants , only : t3ple ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: temp   ! Temperature                   [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! May use ice thermodynamics?   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      logical                            :: frozen ! Derivative of vapour pressure [  Pa/K]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide whether to use liquid water of ice for saturation, based on the temper- !
      ! ature and the settings.                                                            !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the function depending on the previous check.                             !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         rhovsilp = rhovsip(temp)
      else
         rhovsilp = rhovslp(temp)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function rhovsilp
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
   real(kind=4) function tslf(pvap)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pvap      ! Saturation vapour pressure           [    Pa]
      !----- Local variables for iterative method. ----------------------------------------!
      real(kind=4)             :: deriv     ! Function derivative                  [    Pa]
      real(kind=4)             :: fun       ! Function for which we seek a root.   [    Pa]
      real(kind=4)             :: funa      ! Smallest  guess function             [    Pa]
      real(kind=4)             :: funz      ! Largest   guess function             [    Pa]
      real(kind=4)             :: tempa     ! Smallest guess (or previous guess)   [    Pa]
      real(kind=4)             :: tempz     ! Largest guess (new guess in Newton)  [    Pa]
      real(kind=4)             :: delta     ! Aux. var --- 2nd guess for bisection [      ]
      integer                  :: itn       ! Iteration counter                    [   ---]
      integer                  :: itb       ! Iteration counter                    [   ---]
      logical                  :: converged ! Convergence handle                   [   ---]
      logical                  :: zside     ! Flag to check for one-sided approach [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, use Bolton (1980) equation 11, giving es in Pa and T in K -------!
      tempa = (29.65 * log(pvap) - 5016.78)/(log(pvap)-24.0854)
      funa  = eslf(tempa) - pvap
      deriv = eslfp(tempa)
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration ----------------------------!
      tempz = tempa
      fun   = funa
      !------------------------------------------------------------------------------------!


      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, go with bisection -----!

         !----- Copy the previous guess. --------------------------------------------------!
         tempa = tempz
         funa  = fun
         !---------------------------------------------------------------------------------!


         !----- New guess, its function, and derivative evaluation. -----------------------!
         tempz = tempa - fun/deriv
         fun   = eslf(tempz) - pvap
         deriv = eslfp(tempz)
         !---------------------------------------------------------------------------------!


         !----- Check convergence. --------------------------------------------------------!
         converged = abs(tempa-tempz) < toler * tempz
         if (converged) then
            tslf = 0.5 * (tempa+tempz)
            return
         elseif (fun == 0.0) then 
            !----- Converged by luck. -----------------------------------------------------!
            tslf = tempz
            return
         end if
         !---------------------------------------------------------------------------------!
      end do newloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If we have reached this point, then Newton's method has failed.  Use bisection !
      ! instead.  For bisection, we need two guesses whose function has opposite signs.    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
         !----- We already have two guesses with opposite signs. --------------------------!
         funz  = fun
         zside = .true.
         !---------------------------------------------------------------------------------!
      else
         !----- Need to find the guesses with opposite signs. -----------------------------!
         if (abs(fun-funa) < 100.*toler*tempa) then
            delta = 100.*toler*tempa
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Try guesses on both sides of the first guess, increasingly further away      !
         ! until we spot a good guess.                                                     !
         !---------------------------------------------------------------------------------!
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            funz  = eslf(tempz) - pvap
            zside = funa*funz < 0.
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)')           '------------------------------------------'
            write (unit=*,fmt='(a)')           ' Failed finding the second guess:'
            write (unit=*,fmt='(a)')           '------------------------------------------'
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           ' Input:   '
            write (unit=*,fmt='(a,1x,es14.7)') ' + pvap  =',pvap
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           ' Output:  '
            write (unit=*,fmt='(a,1x,es14.7)') ' + tempa =',tempa
            write (unit=*,fmt='(a,1x,es14.7)') ' + funa  =',funa
            write (unit=*,fmt='(a,1x,es14.7)') ' + tempz =',tempz
            write (unit=*,fmt='(a,1x,es14.7)') ' + func  =',funz
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           '------------------------------------------'
            call abort_run  ('Failed finding the second guess for regula falsi'            &
                            ,'tslf','therm_lib.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tslf =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close.  If so, !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tslf-tempa) < toler * tslf
         if (converged) exit bisloop

         !------ Find the new function evaluation. ----------------------------------------!
         fun       =  eslf(tslf) - pvap
         !---------------------------------------------------------------------------------!


         !------ Define the new interval based on the intermediate value theorem. ---------!
         if (fun*funa < 0. ) then
            tempz = tslf
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method). --------!
            if (zside) funa=funa * 0.5
            !----- We have just updated zside, so we set zside to true. -------------------!
            zside = .true.
         else
            tempa = tslf
            funa   = fun
            !----- If we are updating aside again, modify zside (Illinois method). --------!
            if (.not. zside) funz=funz * 0.5
            !----- We have just updated aside, so we set zside to false. ------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         write (unit=*,fmt='(a)')           '------------------------------------------'
         write (unit=*,fmt='(a)')           ' Failed finding the solution:'
         write (unit=*,fmt='(a)')           '------------------------------------------'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' Input:   '
         write (unit=*,fmt='(a,1x,es14.7)') ' + pvap  =',pvap
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' Output:  '
         write (unit=*,fmt='(a,1x,es14.7)') ' + tempa =',tempa
         write (unit=*,fmt='(a,1x,es14.7)') ' + funa  =',funa
         write (unit=*,fmt='(a,1x,es14.7)') ' + tempz =',tempz
         write (unit=*,fmt='(a,1x,es14.7)') ' + func  =',funz
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           '------------------------------------------'
         call abort_run  ('Temperature didn''t converge, we give up!!!'                    &
                         ,'tslf','therm_lib.f90')
      end if
      
      return
   end function tslf
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
   real(kind=4) function tsif(pvap)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pvap      ! Saturation vapour pressure           [    Pa]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=4)             :: deriv     ! Function derivative                  [    Pa]
      real(kind=4)             :: fun       ! Function for which we seek a root.   [    Pa]
      real(kind=4)             :: funa      ! Smallest  guess function             [    Pa]
      real(kind=4)             :: funz      ! Largest   guess function             [    Pa]
      real(kind=4)             :: tempa     ! Smallest guess (or previous guess)   [    Pa]
      real(kind=4)             :: tempz     ! Largest guess (new guess in Newton)  [    Pa]
      real(kind=4)             :: delta     ! Aux. var --- 2nd guess for bisection [      ]
      integer                  :: itn
      integer                  :: itb       ! Iteration counter                    [   ---]
      logical                  :: converged ! Convergence handle                   [   ---]
      logical                  :: zside     ! Flag to check for one-sided approach [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, use Murphy-Koop (2005), equation 8. -----------------------------!
      tempa = (1.814625 * log(pvap) +6190.134)/(29.120 - log(pvap))
      funa  = esif(tempa) - pvap
      deriv = esifp(tempa)
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration ----------------------------!
      tempz = tempa
      fun   = funa
      !------------------------------------------------------------------------------------!


      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, go with bisection -----!

         !----- Copy the previous guess. --------------------------------------------------!
         tempa = tempz
         funa  = fun


         !----- New guess, its function, and derivative evaluation. -----------------------!
         tempz = tempa - fun/deriv
         fun   = esif(tempz) - pvap
         deriv = esifp(tempz)
         !---------------------------------------------------------------------------------!


         !----- Check convergence. --------------------------------------------------------!
         converged = abs(tempa-tempz) < toler * tempz
         if (converged) then
            tsif = 0.5*(tempa+tempz)
            return
         elseif (fun == 0.) then
            tsif = tempz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If we have reached this point, then Newton's method has failed.  Use bisection !
      ! instead.  For bisection, we need two guesses whose function has opposite signs.    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
         !----- We already have two guesses with opposite signs. --------------------------!
         funz  = fun
         zside = .true.
         !---------------------------------------------------------------------------------!
      else
         !----- Need to find the guesses with opposite signs. -----------------------------!
         if (abs(fun-funa) < 100.*toler*tempa) then
            delta = 100.*toler*delta
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Try guesses on both sides of the first guess, increasingly further away      !
         ! until we spot a good guess.                                                     !
         !---------------------------------------------------------------------------------!
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            funz  = esif(tempz) - pvap
            zside = funa*funz < 0.0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)')           '------------------------------------------'
            write (unit=*,fmt='(a)')           ' Failed finding the second guess:'
            write (unit=*,fmt='(a)')           '------------------------------------------'
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           ' Input:   '
            write (unit=*,fmt='(a,1x,es14.7)') ' + pvap  =',pvap
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           ' Output:  '
            write (unit=*,fmt='(a,1x,es14.7)') ' + tempa =',tempa
            write (unit=*,fmt='(a,1x,es14.7)') ' + funa  =',funa
            write (unit=*,fmt='(a,1x,es14.7)') ' + tempz =',tempz
            write (unit=*,fmt='(a,1x,es14.7)') ' + func  =',funz
            write (unit=*,fmt='(a)')           ' '
            write (unit=*,fmt='(a)')           '------------------------------------------'
            call abort_run  ('Failed finding the second guess for regula falsi'            &
                            ,'tsif','therm_lib.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tsif =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close.  If so, !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tsif-tempa) < toler * tsif
         if (converged) exit bisloop
         !---------------------------------------------------------------------------------!

         !------ Find the new function evaluation. ----------------------------------------!
         fun       =  esif(tsif) - pvap
         !---------------------------------------------------------------------------------!


         !------ Define the new interval based on the intermediate value theorem. ---------!
         if (fun*funa < 0. ) then
            tempz = tsif
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method). --------!
            if (zside) funa=funa * 0.5
            !----- We have just updated zside, so we set zside to true. -------------------!
            zside = .true.
         else
            tempa = tsif
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method). --------!
            if (.not. zside) funz=funz * 0.5
            !----- We have just updated aside, so we set zside to false. ------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         write (unit=*,fmt='(a)')           '------------------------------------------'
         write (unit=*,fmt='(a)')           ' Failed finding the solution:'
         write (unit=*,fmt='(a)')           '------------------------------------------'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' Input:   '
         write (unit=*,fmt='(a,1x,es14.7)') ' + pvap  =',pvap
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' Output:  '
         write (unit=*,fmt='(a,1x,es14.7)') ' + tempa =',tempa
         write (unit=*,fmt='(a,1x,es14.7)') ' + funa  =',funa
         write (unit=*,fmt='(a,1x,es14.7)') ' + tempz =',tempz
         write (unit=*,fmt='(a,1x,es14.7)') ' + func  =',funz
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           '------------------------------------------'
         call abort_run  ('Temperature didn''t converge, we give up!!!'                    &
                         ,'tsif','therm_lib.f90')
      end if
      
      return
   end function tsif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the temperature from the ice or liquid mixing ratio.     !
   ! This is truly the inverse of eslf and esif.                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function tslif(pvap,useice)
      use rconstants , only : es3ple ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pvap   ! Vapour pressure                [   Pa]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! May use ice thermodynamics     [  T|F]
      !----- Local variables. -------------------------------------------------------------!
      logical                            :: frozen ! Will use ice thermodynamics    [  T|F]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Since pvap is a function of temperature only, we can check the triple point    !
      ! from the saturation at the triple point, like what we would do for temperature.    !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. pvap < es3ple
      else 
         frozen = bulk_on .and. pvap < es3ple
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the function depending on whether we should use ice.                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         tslif = tsif(pvap)
      else
         tslif = tslf(pvap)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function tslif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the dew point temperature given the pressure and vapour    !
   ! mixing ratio. THIS IS DEW POINT ONLY, WHICH MEANS THAT IT WILL IGNORE ICE EFFECT. For !
   ! a full, triple-point dependent routine use DEWFROSTPOINT.                             !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function dewpoint(pres,rsat)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=4), intent(in) :: rsat    ! Saturation mixing ratio                [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: rsatoff ! Non-singular saturation mixing ratio   [ kg/kg]
      real(kind=4)             :: pvsat   ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !----- Make sure mixing ratio is positive. ------------------------------------------!
      rsatoff  = max(toodry,rsat)
      !------------------------------------------------------------------------------------!


      !----- Find the saturation vapour pressure. -----------------------------------------!
      pvsat    = pres * rsatoff / (ep + rsatoff)
      !------------------------------------------------------------------------------------!

      !----- Dew point is going to be the saturation temperature. -------------------------!
      dewpoint = tslf(pvsat)
      !------------------------------------------------------------------------------------!

      return
   end function dewpoint
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the frost point temperature given the pressure and vapour  !
   ! mixing ratio. THIS IS FROST POINT ONLY, WHICH MEANS THAT IT WILL IGNORE LIQUID        !
   ! EFFECT.  For a full, triple-point dependent routine use DEWFROSTPOINT.                !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function frostpoint(pres,rsat)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=4), intent(in) :: rsat    ! Saturation mixing ratio                [ kg/kg]
      !----- Local variables for iterative method. ----------------------------------------!
      real(kind=4)             :: rsatoff ! Non-singular saturation mixing ratio   [ kg/kg]
      real(kind=4)             :: pvsat   ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !----- Make sure mixing ratio is positive. ------------------------------------------!
      rsatoff    = max(toodry,rsat)
      !------------------------------------------------------------------------------------!


      !----- Find the saturation vapour pressure. -----------------------------------------!
      pvsat      = pres*rsatoff / (ep + rsatoff)
      !------------------------------------------------------------------------------------!

      !----- Frost point is going to be the saturation temperature. -----------------------!
      frostpoint = tsif(pvsat)
      !------------------------------------------------------------------------------------!

      return
   end function frostpoint
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the saturation point temperature given the pressure and    !
   ! vapour mixing ratio. This will check whether the vapour pressure is above or below    !
   ! the triple point vapour pressure, finding dewpoint or frostpoint accordingly.         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function dewfrostpoint(pres,rsat,useice)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres    ! Pressure                     [    Pa]
      real(kind=4), intent(in)           :: rsat    ! Saturation mixing ratio      [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: rsatoff ! Non-singular sat. mix. rat.  [ kg/kg]
      real(kind=4)                       :: pvsat   ! Saturation vapour pressure   [    Pa]
      !------------------------------------------------------------------------------------!


      !----- Make sure mixing ratio is positive. ------------------------------------------!
      rsatoff  = max(toodry,rsat)
      !------------------------------------------------------------------------------------!

      !----- Find the saturation vapour pressure. -----------------------------------------!
      pvsat         = pres*rsatoff / (ep + rsatoff)
      !------------------------------------------------------------------------------------!

      !----- Dew (frost) point is going to be the saturation temperature. -----------------!
      if (present(useice)) then
         dewfrostpoint = tslif(pvsat,useice)
      else
         dewfrostpoint = tslif(pvsat)
      end if
      !------------------------------------------------------------------------------------!
      return
   end function dewfrostpoint
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio (or specific humidity) based on    !
   ! the pressure [Pa], temperature [K] and relative humidity [fraction]. IT ALWAYS        !
   ! ASSUMES THAT RELATIVE HUMIDITY IS WITH RESPECT TO THE LIQUID PHASE.  Ptrh2rvapil      !
   ! checks which one to use depending on whether temperature is more or less than the     !
   ! triple point.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function ptrh2rvapl(relh,pres,temp,out_shv)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: relh    ! Relative humidity                      [    --]
      real(kind=4), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=4), intent(in) :: temp    ! Temperature                            [     K]
      logical     , intent(in) :: out_shv ! Output is specific humidity            [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=4)             :: relhh   ! Bounded relative humidity              [    --]
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.,max(0.,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure. -----------------------------------------------------!
      pvap  = relhh * eslf(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapl = max(toodry, ep * pvap / (pres - (1.0 - ep) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapl = max(toodry, ep * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function ptrh2rvapl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio (or specific humidity) based on    !
   ! the pressure [Pa], temperature [K] and relative humidity [fraction]. IT ALWAYS        !
   ! ASSUMES THAT RELATIVE HUMIDITY IS WITH RESPECT TO THE ICE PHASE.  Ptrh2rvapil         !
   ! checks which one to use depending on whether temperature is more or less than the     !
   ! triple point.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function ptrh2rvapi(relh,pres,temp,out_shv)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: relh    ! Relative humidity                      [    --]
      real(kind=4), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=4), intent(in) :: temp    ! Temperature                            [     K]
      logical     , intent(in) :: out_shv ! Output is specific humidity            [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=4)             :: relhh   ! Bounded relative humidity              [    --]
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.,max(0.,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure. -----------------------------------------------------!
      pvap  = relhh * esif(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapi = max(toodry, ep * pvap / (pres - (1.0 - ep) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapi = max(toodry, ep * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function ptrh2rvapi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based (or specific humidity) based !
   ! on the pressure [Pa], temperature [K] and relative humidity [fraction].  It checks    !
   ! the temperature to decide between ice or liquid saturation.                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function ptrh2rvapil(relh,pres,temp,out_shv,useice)
      use rconstants , only : ep      & ! intent(in)
                            , toodry  & ! intent(in)
                            , t3ple   ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: relh    ! Relative humidity            [    --]
      real(kind=4), intent(in)           :: pres    ! Pressure                     [    Pa]
      real(kind=4), intent(in)           :: temp    ! Temperature                  [     K]
      logical     , intent(in)           :: out_shv ! Output is specific humidity  [   T|F]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=4)                       :: relhh   ! Bounded relative humidity    [    --]
      logical                            :: frozen  ! Will use ice thermodynamics  [   T|F]
      !------------------------------------------------------------------------------------!


      !----- Check whether to use the user's or the default flag for ice saturation. ------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.,max(0.,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure (ice or liquid, depending on the value of frozen). ---!
      if (frozen) then
         pvap  = relhh * esif(temp)
      else
         pvap  = relhh * eslf(temp)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapil = max(toodry, ep * pvap / (pres - (1.0 - ep) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapil = max(toodry, ep * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
      return
   end function ptrh2rvapil
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio (or specific humidity). Two important points:       !
   ! 1. IT ALWAYS ASSUME THAT RELATIVE HUMIDITY IS WITH RESPECT TO THE LIQUID PHASE.       !
   !    If you want to switch between ice and liquid, use rehuil instead.                  !
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rehul(pres,temp,humi,is_shv)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres    ! Air pressure                           [    Pa]
      real(kind=4), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=4), intent(in) :: humi    ! Humidity                               [ kg/kg]
      logical     , intent(in) :: is_shv  ! Input humidity is specific humidity    [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: shv     ! Specific humidity                      [ kg/kg]
      real(kind=4)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=4)             :: psat    ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      psat = eslf (temp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehul = max(0. ,pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehul
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio (or specific humidity). Two important points:       !
   ! 1. IT ALWAYS ASSUME THAT RELATIVE HUMIDITY IS WITH RESPECT TO THE ICE PHASE.          !
   !    If you want to switch between ice and liquid, use rehuil instead.                  !
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rehui(pres,temp,humi,is_shv)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres    ! Air pressure                           [    Pa]
      real(kind=4), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=4), intent(in) :: humi    ! Humidity                               [ kg/kg]
      logical     , intent(in) :: is_shv  ! Input humidity is specific humidity    [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: shv     ! Specific humidity                      [ kg/kg]
      real(kind=4)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=4)             :: psat    ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      psat = esif (temp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehui = max(0. ,pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehui
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio (or specific humidity). Two important points:       !
   ! 1. It may consider whether the temperature is above or below the freezing point       !
   !    to choose which saturation to use. It is possible to explicitly force not to use   !
   !    ice in case level is 2 or if you have reasons not to use ice (e.g. reading data    !
   !    that did not consider ice).
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rehuil(pres,temp,humi,is_shv,useice)
      use rconstants , only : t3ple  & ! intent(in)
                            , ep     & ! intent(in)
                            , toodry ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres    ! Air pressure                 [    Pa]
      real(kind=4), intent(in)           :: temp    ! Temperature                  [     K]
      real(kind=4), intent(in)           :: humi    ! Humidity                     [ kg/kg]
      logical     , intent(in)           :: is_shv  ! Input is specific humidity   [   T|F]
       !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: shv     ! Specific humidity            [ kg/kg]
      real(kind=4)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=4)                       :: psat    ! Saturation vapour pressure   [    Pa]
      logical                            :: frozen  ! Will use ice saturation now  [   T|F]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    Check whether we should use ice or liquid saturation.                           !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      if (frozen) then
         psat = esif (temp)
      else
         psat = esif (temp)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehuil = max(0. ,pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehuil
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour pressure deficit based on pressure, temper-     !
   ! ature, and vapour mixing ratio (or specific humidity).                                !
   !                                                                                       !
   ! IMPORTANT: IT ALWAYS ASSUMES THAT VAPOUR PRESSURE DEFICIT IS WITH RESPECT TO THE      !
   !            LIQUID PHASE.  If you would like it to switch between ice and liquid, then !
   !            use vpdefil instead.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function vpdefl(pres,temp,humi,is_shv)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres    ! Air pressure                           [    Pa]
      real(kind=4), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=4), intent(in) :: humi    ! Humidity                               [ kg/kg]
      logical     , intent(in) :: is_shv  ! Input humidity is specific humidity    [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: shv     ! Specific humidity                      [ kg/kg]
      real(kind=4)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=4)             :: psat    ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      psat = eslf(temp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      vpdefl = max(0.0 , psat - pvap)
      !------------------------------------------------------------------------------------!

      return
   end function vpdefl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour pressure deficit based on pressure, temper-     !
   ! ature, and vapour mixing ratio (or specific humidity).                                !
   !                                                                                       !
   ! IMPORTANT: IT ALWAYS ASSUMES THAT VAPOUR PRESSURE DEFICIT IS WITH RESPECT TO THE      !
   !            ICE PHASE.  If you would like it to switch between ice and liquid, then    !
   !            use vpdefil instead.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function vpdefi(pres,temp,humi,is_shv)
      use rconstants , only : ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres    ! Air pressure                           [    Pa]
      real(kind=4), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=4), intent(in) :: humi    ! Humidity                               [ kg/kg]
      logical     , intent(in) :: is_shv  ! Input humidity is specific humidity    [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: shv     ! Specific humidity                      [ kg/kg]
      real(kind=4)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=4)             :: psat    ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      psat = esif(temp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      vpdefi = max(0.0 , psat - pvap)
      !------------------------------------------------------------------------------------!

      return
   end function vpdefi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour pressure deficit based on pressure, temper-     !
   ! ature, and vapour mixing ratio (or specific humidity).                                !
   !                                                                                       !
   ! IMPORTANT: This fucntion may consider whether the temperature is above or below the   !
   !            freezing point to choose which saturation to use. It is possible to        !
   !            explicitly force not to use ice in case level is 2 or if you have reasons  !
   !            not to use ice (e.g. reading data that did not consider ice).              !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function vpdefil(pres,temp,humi,is_shv,useice)
      use rconstants , only : t3ple  & ! intent(in)
                            , ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres    ! Air pressure                 [    Pa]
      real(kind=4), intent(in)           :: temp    ! Temperature                  [     K]
      real(kind=4), intent(in)           :: humi    ! Humidity                     [ kg/kg]
      logical     , intent(in)           :: is_shv  ! Input is specific humidity   [   T|F]
       !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: shv     ! Specific humidity            [ kg/kg]
      real(kind=4)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=4)                       :: psat    ! Saturation vapour pressure   [    Pa]
      logical                            :: frozen  ! Will use ice saturation now  [   T|F]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    Check whether we should use ice or liquid saturation.                           !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else 
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      if (frozen) then
         psat = esif(temp)
      else
         psat = esif(temp)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      vpdefil = max(0.0 , psat - pvap)
      !------------------------------------------------------------------------------------!

      return
   end function vpdefil
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
   real(kind=4) function tv2temp(tvir,rvap,rtot)
      use rconstants , only : epi ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)           :: tvir     ! Virtual temperature          [    K]
      real(kind=4), intent(in)           :: rvap     ! Vapour mixing ratio          [kg/kg]
      real(kind=4), intent(in), optional :: rtot     ! Total mixing ratio           [kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: rtothere ! Total or vapour mixing ratio [kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Prefer using total mixing ratio, but if it isn't provided, then use vapour    !
      ! as total (no condensation).                                                        !
      !------------------------------------------------------------------------------------!
      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvap
      end if
      !------------------------------------------------------------------------------------!


      !----- Convert using a generalised function. ----------------------------------------!
      tv2temp = tvir * (1. + rtothere) / (1. + epi*rvap)
      !------------------------------------------------------------------------------------!

      return
   end function tv2temp
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
   real(kind=4) function virtt(temp,rvap,rtot)
      use rconstants , only: epi ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)           :: temp     ! Temperature                  [    K]
      real(kind=4), intent(in)           :: rvap     ! Vapour mixing ratio          [kg/kg]
      real(kind=4), intent(in), optional :: rtot     ! Total mixing ratio           [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=4)                       :: rtothere ! Total or vapour mixing ratio [kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Prefer using total mixing ratio, but if it isn't provided, then use vapour    !
      ! as total (no condensation).                                                        !
      !------------------------------------------------------------------------------------!
      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvap
      end if
      !------------------------------------------------------------------------------------!


      !----- Convert using a generalised function. ----------------------------------------!
      virtt = temp * (1. + epi * rvap) / (1. + rtothere)
      !------------------------------------------------------------------------------------!

      return
   end function virtt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the density based on the virtual temperature and the ideal    !
   ! gas law. The condensed phase will be taken into account if the user provided both     !
   ! the vapour and the total mixing ratios.                                               !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function idealdens(pres,temp,rvap,rtot)
      use rconstants , only : rdry ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)           :: pres ! Pressure                        [    Pa]
      real(kind=4), intent(in)           :: temp ! Temperature                     [     K]
      real(kind=4), intent(in)           :: rvap ! Vapour mixing ratio             [ kg/kg]
      real(kind=4), intent(in), optional :: rtot ! Total mixing ratio              [ kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=4)                       :: tvir ! Virtual temperature             [     K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Prefer using total mixing ratio, but if it isn't provided, then use vapour    !
      ! as total (no condensation).                                                        !
      !------------------------------------------------------------------------------------!
      if (present(rtot)) then
        tvir = virtt(temp,rvap,rtot)
      else
        tvir = virtt(temp,rvap)
      end if
      !------------------------------------------------------------------------------------!


      !----- Convert using the definition of virtual temperature. -------------------------!
      idealdens = pres / (rdry * tvir)
      !------------------------------------------------------------------------------------!

      return
   end function idealdens
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the density based on the virtual temperature and the ideal    !
   ! gas law.  The only difference between this function and the one above is that here we !
   ! provide vapour and total specific mass (specific humidity) instead of mixing ratio.   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function idealdenssh(pres,temp,qvpr,qtot)
      use rconstants , only : rdry & ! intent(in)
                            , epi  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)           :: pres ! Pressure                        [    Pa]
      real(kind=4), intent(in)           :: temp ! Temperature                     [     K]
      real(kind=4), intent(in)           :: qvpr ! Vapour specific mass            [ kg/kg]
      real(kind=4), intent(in), optional :: qtot ! Total water specific mass       [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: qall ! Either qtot or qvpr...          [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Prefer using total specific humidity, but if it isn't provided, then use      !
      ! vapour phase as the total (no condensation).                                       !
      !------------------------------------------------------------------------------------!
      if (present(qtot)) then
        qall = qtot
      else
        qall = qvpr
      end if
      !------------------------------------------------------------------------------------!


      !----- Convert using a generalised function. ----------------------------------------!
      idealdenssh = pres / (rdry * temp * (1. - qall + epi * qvpr))
      !------------------------------------------------------------------------------------!

      return
   end function idealdenssh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes reduces the pressure from the reference height to the      !
   ! canopy height by assuming hydrostatic equilibrium.  For simplicity, we assume that    !
   ! R and cp are constants (in reality they are dependent on humidity).                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function reducedpress(pres,thetaref,shvref,zref,thetacan,shvcan,zcan)
      use rconstants , only : epim1    & ! intent(in)
                            , p00k     & ! intent(in)
                            , rocp     & ! intent(in)
                            , cpor     & ! intent(in)
                            , cpdry    & ! intent(in)
                            , grav     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres     ! Pressure                            [      Pa]
      real(kind=4), intent(in) :: thetaref ! Potential temperature               [       K]
      real(kind=4), intent(in) :: shvref   ! Vapour specific mass                [   kg/kg]
      real(kind=4), intent(in) :: zref     ! Height at reference level           [       m]
      real(kind=4), intent(in) :: thetacan ! Potential temperature               [       K]
      real(kind=4), intent(in) :: shvcan   ! Vapour specific mass                [   kg/kg]
      real(kind=4), intent(in) :: zcan     ! Height at canopy level              [       m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=4)             :: pinc     ! Pressure increment                  [ Pa^R/cp]
      real(kind=4)             :: thvbar   ! Average virtual pot. temperature    [       K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      First we compute the average virtual potential temperature between the canopy !
      ! top and the reference level.                                                       !
      !------------------------------------------------------------------------------------!
      thvbar = 0.5 * (thetaref * (1. + epim1 * shvref) + thetacan * (1. + epim1 * shvcan))
      !------------------------------------------------------------------------------------!



      !----- Then, we find the pressure gradient scale. -----------------------------------!
      pinc = grav * p00k * (zref - zcan) / (cpdry * thvbar)
      !------------------------------------------------------------------------------------!



      !----- And we can find the reduced pressure. ----------------------------------------!
      reducedpress = (pres**rocp + pinc ) ** cpor
      !------------------------------------------------------------------------------------!

      return
   end function reducedpress
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the Exner function [J/kg/K], given the pressure.  It       !
   ! assumes for simplicity that R and Cp are constants and equal to the dry air values.   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function press2exner(pres)
      use rconstants , only : p00i           & ! intent(in)
                            , cpdry          & ! intent(in)
                            , rocp           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres   ! Pressure                               [     Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      press2exner = cpdry * ( pres * p00i ) ** rocp
      !------------------------------------------------------------------------------------!

      return
   end function press2exner
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the pressure [Pa], given the Exner function.  Like in the  !
   ! function above, we also assume R and Cp to be constants and equal to the dry air      !
   ! values.                                                                               !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function exner2press(exner)
      use rconstants , only : p00            & ! intent(in)
                            , cpdryi         & ! intent(in)
                            , cpor           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      exner2press = p00 * ( exner * cpdryi ) ** cpor
      !------------------------------------------------------------------------------------!

      return
   end function exner2press
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the potential temperature [K], given the Exner function    !
   ! and temperature.  For simplicity we ignore the effects of humidity in R and cp and    !
   ! use the dry air values instead.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function extemp2theta(exner,temp)
      use rconstants , only : cpdry          ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=4), intent(in) :: temp   ! Temperature                            [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extemp2theta = cpdry * temp / exner
      !------------------------------------------------------------------------------------!

      return
   end function extemp2theta
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature [K], given the Exner function and          !
   ! potential temperature.  We simplify the equations by assuming that R and Cp are       !
   ! constants.                                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function extheta2temp(exner,theta)
      use rconstants , only : p00i           & ! intent(in)
                            , cpdryi         ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=4), intent(in) :: theta  ! Potential temperature                  [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extheta2temp = cpdryi * exner * theta
      !------------------------------------------------------------------------------------!

      return
   end function extheta2temp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the specific (intensive) internal energy of water [J/kg],  !
   ! given the temperature and liquid fraction.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function tl2uint(temp,fliq)
      use rconstants , only : cice           & ! intent(in)
                            , cliq           & ! intent(in)
                            , tsupercool_liq ! ! intent(in)
                            
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=4), intent(in) :: fliq  ! Fraction liquid water                    [ kg/kg]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of ice and liquid     !
      ! phases.                                                                            !
      !------------------------------------------------------------------------------------!
      tl2uint = (1.0 - fliq) * cice * temp + fliq * cliq * (temp - tsupercool_liq)
      !------------------------------------------------------------------------------------!

      return
   end function tl2uint
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the extensive internal energy of water [J/m²] or [  J/m³], !
   ! given the temperature [K], the heat capacity of the "dry" part [J/m²/K] or [J/m³/K],  !
   ! water mass [ kg/m²] or [ kg/m³], and liquid fraction [---].                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function cmtl2uext(dryhcap,wmass,temp,fliq)
      use rconstants , only : cice           & ! intent(in)
                            , cliq           & ! intent(in)
                            , tsupercool_liq ! ! intent(in)
                            
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m²/K] or [J/m³/K]
      real(kind=4), intent(in)  :: wmass   ! Water mass                [ kg/m²] or [ kg/m³]
      real(kind=4), intent(in)  :: temp    ! Temperature                           [     K]
      real(kind=4), intent(in)  :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of dry part, plus the !
      ! contribution of ice and liquid phases.                                             !
      !------------------------------------------------------------------------------------!
      cmtl2uext = dryhcap * temp + wmass * ( (1.0 - fliq) * cice * temp                    &
                                           + fliq * cliq * (temp - tsupercool_liq) )
      !------------------------------------------------------------------------------------!

      return
   end function cmtl2uext
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the specific enthalpy [J/kg] given the temperature and     !
   ! humidity (either mixing ratio or specific humidity).  If we assume that latent heat   !
   ! of vaporisation is a linear function of temperature (equivalent to assume that        !
   ! specific heats are constants and that the thermal expansion of liquids and solids are !
   ! negligible), then the saturation disappears and the enthalpy becomes a straight-      !
   ! forward state function.  In case we are accounting for the water exchange only        !
   ! (latent heat), set the specific humidity to 1.0 and multiply the result by water mass !
   ! or water flux.                                                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function tq2enthalpy(temp,humi,is_shv)
      use rconstants , only : cpdry          & ! intent(in)
                            , cph2o          & ! intent(in)
                            , tsupercool_vap ! ! intent(in)
                            
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp   ! Temperature                             [     K]
      real(kind=4), intent(in) :: humi   ! Humidity (spec. hum. or mixing ratio)   [ kg/kg]
      logical     , intent(in) :: is_shv ! Input humidity is specific humidity     [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: shv    ! Specific humidity                       [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.0)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      tq2enthalpy = (1.0 - shv) * cpdry * temp + shv * cph2o * (temp - tsupercool_vap)  
      !------------------------------------------------------------------------------------!

      return
   end function tq2enthalpy
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature [K] given the specific enthalpy and        !
   ! humidity.  If we assume that latent heat of vaporisation is a linear function of      !
   ! temperature (equivalent to assume that specific heats are constants and that the      !
   ! thermal expansion of liquid and water are negligible), then the saturation disappears !
   ! and the enthalpy becomes a straightforward state function.  In case you are looking   !
   ! at water exchange only, set the specific humidity to 1.0 and multiply the result by   !
   ! the water mass or water flux.                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function hq2temp(enthalpy,humi,is_shv)
      use rconstants , only : cpdry          & ! intent(in)
                            , cph2o          & ! intent(in)
                            , tsupercool_vap ! ! intent(in)
                            
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: enthalpy ! Specific enthalpy                     [  J/kg]
      real(kind=4), intent(in) :: humi     ! Humidity (spec. hum. or mixing ratio) [ kg/kg]
      logical     , intent(in) :: is_shv   ! Input humidity is specific humidity   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: shv      ! Specific humidity                     [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.0)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      hq2temp = ( enthalpy + shv * cph2o * tsupercool_vap )                                &
              / ( (1.0 - shv) * cpdry + shv * cph2o )
      !------------------------------------------------------------------------------------!

      return
   end function hq2temp
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of vaporisation for a given temperature.  If  !
   ! we use the definition of latent heat (difference in enthalpy between liquid and       !
   ! vapour phases), and assume that the specific heats are constants, latent heat becomes !
   ! a linear function of temperature.                                                     !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function alvl(temp)
      use rconstants , only : alvl3  & ! intent(in)
                            , dcpvl  & ! intent(in)
                            , t3ple  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvl = alvl3 + dcpvl * (temp - t3ple)
      !------------------------------------------------------------------------------------!

      return
   end function alvl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of sublimation for a given temperature.  If   !
   ! we use the definition of latent heat (difference in enthalpy between ice and vapour   !
   ! phases), and assume that the specific heats are constants, latent heat becomes a      !
   ! linear function of temperature.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function alvi(temp)
      use rconstants , only : alvi3  & ! intent(in)
                            , dcpvi  & ! intent(in)
                            , t3ple  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvi = alvi3 + dcpvi * (temp - t3ple)
      !------------------------------------------------------------------------------------!

      return
   end function alvi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the ice liquid potential temperature given the Exner       !
   ! function [J/kg/K], temperature [K], and liquid and ice mixing ratios [kg/kg].         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function theta_iceliq(exner,temp,rliq,rice)
      use rconstants , only : alvl3      & ! intent(in)
                            , alvi3      & ! intent(in)
                            , cpdry      & ! intent(in)
                            , ttripoli   & ! intent(in)
                            , htripoli   & ! intent(in)
                            , htripolii  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner ! Exner function                          [ J/kg/K]
      real(kind=4), intent(in) :: temp  ! Temperature                             [      K]
      real(kind=4), intent(in) :: rliq  ! Liquid mixing ratio                     [  kg/kg]
      real(kind=4), intent(in) :: rice  ! Ice mixing ratio                        [  kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: hh    ! Enthalpy associated with sensible heat  [   J/kg]
      real(kind=4)             :: qq    ! Enthalpy associated with latent heat    [   J/kg]
      !------------------------------------------------------------------------------------!


      !----- Find the sensible heat enthalpy (assuming dry air). --------------------------!
      hh = cpdry * temp
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the latent heat enthalpy.  If using the old thermodynamics, we use the   !
      ! latent heat at T = T3ple, otherwise we use the temperature-dependent one.          !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         qq = alvl(temp) * rliq + alvi(temp) * rice
      else
         qq = alvl3 * rliq + alvi3 * rice
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Solve the thermodynamics.  For the new thermodynamics we don't approximate    !
      ! the exponential to a linear function, nor do we impose temperature above the thre- !
      ! shold from Tripoli and Cotton (1981).                                              !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Decide how to compute, based on temperature. ------------------------------!
         theta_iceliq = hh * exp(-qq / hh) / exner
         !---------------------------------------------------------------------------------!
      else
         !----- Decide how to compute, based on temperature. ------------------------------!
         if (temp > ttripoli) then
            theta_iceliq = hh * hh / (exner * ( hh + qq))
         else
            theta_iceliq = hh * htripoli / (exner * ( htripoli + qq))
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function theta_iceliq
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the liquid potential temperature derivative with respect   !
   ! to temperature, useful in iterative methods.                                          !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function dthetail_dt(condconst,thil,exner,pres,temp,rliq,ricein)
      use rconstants , only : alvl3      & ! intent(in)
                            , alvi3      & ! intent(in)
                            , dcpvi      & ! intent(in)
                            , dcpvl      & ! intent(in)
                            , cpdry      & ! intent(in)
                            , ttripoli   & ! intent(in)
                            , htripoli   & ! intent(in)
                            , htripolii  & ! intent(in)
                            , t3ple      ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      logical     , intent(in)           :: condconst  ! Condensation is constant? [   T|F]
      real(kind=4), intent(in)           :: thil       ! Ice liquid pot. temp.     [     K]
      real(kind=4), intent(in)           :: exner      ! Exner function            [J/kg/K]
      real(kind=4), intent(in)           :: pres       ! Pressure                  [    Pa]
      real(kind=4), intent(in)           :: temp       ! Temperature               [     K]
      real(kind=4), intent(in)           :: rliq       ! Liquid mixing ratio       [ kg/kg]
      real(kind=4), intent(in), optional :: ricein     ! Ice mixing ratio          [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: rice       ! Ice mixing ratio or 0.    [ kg/kg]
      real(kind=4)                       :: ldrst      ! L × d(rs)/dT × T          [  J/kg]
      real(kind=4)                       :: rdlt       ! r × d(L)/dT × T           [  J/kg]
      real(kind=4)                       :: hh         ! Sensible heat enthalpy    [  J/kg]
      real(kind=4)                       :: qq         ! Latent heat enthalpy      [  J/kg]
      logical                            :: thereisice ! Is ice present            [   ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Check whether we should consider ice thermodynamics or not.                      !
      !------------------------------------------------------------------------------------!
      thereisice = present(ricein)
      if (thereisice) then
         rice = ricein
      else
         rice = 0.
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether the current state has condensed water.                           !
      !------------------------------------------------------------------------------------!
      if (rliq+rice == 0.) then
         !----- No condensation, so dthetail_dt is a constant. ----------------------------!
         dthetail_dt = thil/temp
         return
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Condensation exists.  Compute some auxiliary variables.                     !
         !---------------------------------------------------------------------------------!


         !---- Sensible heat enthalpy. ----------------------------------------------------!
         hh = cpdry * temp
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Find the latent heat enthalpy.  If using the old thermodynamics, we use    !
         ! the latent heat at T = T3ple, otherwise we use the temperature-dependent one.   !
         ! The term r × d(L)/dT × T is computed only when we use the new thermodynamics.   !
         !---------------------------------------------------------------------------------!
         if (newthermo) then
            qq   = alvl(temp) * rliq + alvi(temp) * rice
            rdlt = (dcpvl * rliq + dcpvi * rice ) * temp
         else
            qq   = alvl3 * rliq + alvi3 * rice
            rdlt = 0.0
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    This is the term L×[d(rs)/dt]×T. L may be either the vapourisation or        !
         ! sublimation latent heat, depending on the temperature and whether we are consi- !
         ! dering ice or not.  We still need to check whether latent heat is a function of !
         ! temperature or not.  Also, if condensation mixing ratio is constant, then this  !
         ! term will be always zero.                                                       !
         !---------------------------------------------------------------------------------!
         if (condconst) then
            ldrst = 0.
         elseif (thereisice .and. temp < t3ple) then
            if (newthermo) then
               ldrst = alvi3 * rsifp(pres,temp) * temp
            else
               ldrst = alvi(temp) * rsifp(pres,temp) * temp
            end if
         else
            if (newthermo) then
               ldrst = alvl3 * rslfp(pres,temp) * temp
            else
               ldrst = alvl(temp) * rslfp(pres,temp) * temp
            end if
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the condensed phase consistent with the thermodynamics used.              !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         dthetail_dt = thil * ( 1. + (ldrst + qq - rdlt ) / hh ) / temp
      else
         !----- Decide how to compute, based on temperature. ------------------------------!
         if (temp > ttripoli) then
            dthetail_dt = thil * ( 1. + (ldrst + qq) / (hh+qq) ) / temp
         else
            dthetail_dt = thil * ( 1. + ldrst / (htripoli + alvl3 * rliq) ) / temp
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function dthetail_dt
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
   real(kind=4) function thil2temp(thil,exner,pres,rliq,rice,t1stguess)
      use rconstants , only : cpdry      & ! intent(in)
                            , cpdryi     & ! intent(in)
                            , cpdryi4    & ! intent(in)
                            , alvl3      & ! intent(in)
                            , alvi3      & ! intent(in)
                            , t00        & ! intent(in)
                            , t3ple      & ! intent(in)
                            , ttripoli   & ! intent(in)
                            , htripolii  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: thil      ! Ice-liquid water potential temp.     [     K]
      real(kind=4), intent(in) :: exner     ! Exner function                       [J/kg/K]
      real(kind=4), intent(in) :: pres      ! Pressure                             [    Pa]
      real(kind=4), intent(in) :: rliq      ! Liquid water mixing ratio            [ kg/kg]
      real(kind=4), intent(in) :: rice      ! Ice mixing ratio                     [ kg/kg]
      real(kind=4), intent(in) :: t1stguess ! 1st. guess for temperature           [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: til       ! Ice liquid temperature               [     K]
      real(kind=4)             :: deriv     ! Function derivative 
      real(kind=4)             :: fun       ! Function for which we seek a root.
      real(kind=4)             :: funa      ! Smallest  guess function
      real(kind=4)             :: funz      ! Largest   guess function
      real(kind=4)             :: tempa     ! Smallest  guess (or previous guess in Newton)
      real(kind=4)             :: tempz     ! Largest   guess (or new guess in Newton)
      real(kind=4)             :: delta     ! Aux. var to compute 2nd guess for bisection
      integer                  :: itn       ! Iteration counter
      integer                  :: itb       ! Iteration counter
      logical                  :: converged ! Convergence flag
      logical                  :: zside     ! Flag to check for one-sided approach...
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     First we check for conditions that don't require iterative root-finding.       !
      !------------------------------------------------------------------------------------!
      if (rliq + rice == 0.) then
         !----- No condensation.  Theta_il is the same as theta. --------------------------!
         thil2temp = cpdryi * thil * exner
         return
         !---------------------------------------------------------------------------------!
      elseif (.not. newthermo) then
         !---------------------------------------------------------------------------------!
         !    There is condensation but we are using the old thermodynamics, which can be  !
         ! solved analytically.                                                            !
         !---------------------------------------------------------------------------------!
         til = cpdryi * thil * exner
         if (t1stguess > ttripoli) then
            thil2temp = 0.5                                                                &
                      * (til + sqrt(til * (til + cpdryi4 * (alvl3 * rliq + alvi3 * rice))))
         else
            thil2temp = til * ( 1. + (alvl3 * rliq + alvi3 * rice) * htripolii)
         end if
         return
         !---------------------------------------------------------------------------------!
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

      !------------------------------------------------------------------------------------!
      !      If we have reached this point, we must use root-finding method.  For the      !
      ! Newton's 1st. guess, use t1stguess.                                                !
      !------------------------------------------------------------------------------------!
      tempz     = t1stguess
      fun       = theta_iceliq(exner,tempz,rliq,rice)
      deriv     = dthetail_dt(.true.,fun,exner,pres,tempz,rliq,rice)
      fun       = fun - thil
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,2(1x,a,1x,es12.5))')            &
      !   'itn=',0,'bisection=',.false.,'tempz=',tempz-t00                                  &
      !           ,'fun=',fun,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (fun == 0.) then
         thil2temp = tempz
         converged = .true.
         return
      else 
         tempa = tempz
         funa  = fun
      end if

      !----- Enter loop: it will probably skip when the air is not saturated --------------!
      converged = .false.
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, go to bisection -------!
         tempa = tempz
         funa  = fun
         tempz = tempa - fun/deriv
         fun   = theta_iceliq(exner,tempz,rliq,rice)
         deriv = dthetail_dt(.true.,fun,exner,pres,tempz,rliq,rice)
         fun   = fun - thil
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,2(1x,a,1x,es12.5))')         &
         !   'itn=',itn,'bisection=',.false.,'tempz=',tempz-t00                             &
         !           ,'fun=',fun,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         converged = abs(tempa-tempz) < toler*tempz
         !----- Converged, happy with that, return the average b/w the 2 previous guesses -!
         if (fun == 0.) then
            thil2temp = tempz
            converged = .true.
            return
         elseif(converged) then
            thil2temp = 0.5 * (tempa+tempz)
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If we have reached this point then Newton's method failed.  Use bisection      !
      ! instead.  For bisection, We need two guesses whose function evaluations have       !
      ! opposite sign.                                                                     !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
         !----- Guesses have opposite sign. -----------------------------------------------!
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < toler*tempa) then
            delta = 100.*toler*tempa
         else 
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            funz  = theta_iceliq(exner,tempz,rliq,rice) - thil
            zside = funa*funz < 0.0
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
            call abort_run  ('Failed finding the second guess for regula falsi'            &
                          ,'thil2temp','therm_lib.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         thil2temp =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(thil2temp-tempa)< toler*thil2temp 
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         fun  = theta_iceliq(exner,tempz,rliq,rice) - thil

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,6(1x,a,1x,f11.4))')         &
         !   'itn=',itb,'bisection=',.true.                                                 &
         !           ,'temp=',thil2temp-t00,'tempa=',tempa-t00,'tempz=',tempz-t00           &
         !           ,'fun=',fun,'funa=',funa,'funz=',funz
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0. ) then
            tempz = thil2temp
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 0.5
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = thil2temp
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz=funz * 0.5
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
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rliq            [  g/kg] =',1000.*rliq
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rice            [  g/kg] =',1000.*rice
         write (unit=*,fmt='(a,1x,f12.4)' ) 't1stguess       [    °C] =',t1stguess-t00
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tempa           [    °C] =',tempa-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tempz           [    °C] =',tempz-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [     K] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [     K] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [     K] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [  ----] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [  ----] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                                            ,abs(thil2temp-tempa)/thil2temp
         write (unit=*,fmt='(a,1x,f12.4)' ) 'thil2temp       [     K] =',thil2temp

         call abort_run  ('Temperature didn''t converge, giving up!!!'                     &
                       ,'thil2temp','therm_lib.f90')
      end if

      return
   end function thil2temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the partial derivative of temperature as a function of the !
   ! saturation mixing ratio [kg/kg], keeping pressure constant.  This is based on the     !
   ! ice-liquid potential temperature equation.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function dtempdrs(exner,thil,temp,rliq,rice,rconmin)
      use rconstants , only : alvl3      & ! intent(in)
                            , alvi3      & ! intent(in)
                            , dcpvl      & ! intent(in)
                            , dcpvi      & ! intent(in)
                            , cpdry      & ! intent(in)
                            , cpdryi     & ! intent(in)
                            , ttripoli   & ! intent(in)
                            , htripolii  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner   ! Exner function                        [ J/kg/K]
      real(kind=4), intent(in) :: thil    ! Ice-liquid potential temperature (*)  [      K]
      real(kind=4), intent(in) :: temp    ! Temperature                           [      K]
      real(kind=4), intent(in) :: rliq    ! Liquid mixing ratio                   [  kg/kg]
      real(kind=4), intent(in) :: rice    ! Ice mixing ratio                      [  kg/kg]
      real(kind=4), intent(in) :: rconmin ! Min. non-zero condensate mix. ratio   [  kg/kg]
      !------------------------------------------------------------------------------------!
      ! (*) Thil is not used in this formulation but it may be used should you opt by      !
      !     other ways to compute theta_il, so don't remove this argument.                 !
      !------------------------------------------------------------------------------------!
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: qq      ! Enthalpy -- latent heat               [   J/kg]
      real(kind=4)             :: qpt     ! d(qq)/dT * T                          [   J/kg]
      real(kind=4)             :: hh      ! Enthalpy -- sensible heat             [   J/kg]
      real(kind=4)             :: rcon    ! Condensate mixing ratio               [  kg/kg]
      real(kind=4)             :: til     ! Ice-liquid temperature                [      K]
      !------------------------------------------------------------------------------------!


      !----- Find the total hydrometeor mixing ratio. -------------------------------------!
      rcon  = rliq+rice
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      If the amount of condensate is negligible, temperature does not depend on     !
      ! saturation mixing ratio.                                                           !
      !------------------------------------------------------------------------------------!
      if (rcon < rconmin) then
         dtempdrs = 0.
      else

         !---------------------------------------------------------------------------------!
         !     Find the enthalpy associated with latent heat and its derivative            !
         ! correction.  This is dependent on the thermodynamics used.                      !
         !---------------------------------------------------------------------------------!
         if (newthermo) then
            qq  = alvl(temp) * rliq + alvi(temp) * rice
            qpt = dcpvl * rliq + dcpvi * rice
         else
            qq  = alvl3 * rliq + alvi3 * rice
            qpt = 0.0
         end if
         !---------------------------------------------------------------------------------!


         !----- Find the enthalpy associated with sensible heat. --------------------------!
         hh = cpdry * temp
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Decide how to compute, based on the thermodynamics method.                   !
         !---------------------------------------------------------------------------------!
         if (newthermo) then
            dtempdrs = - temp * qq / ( rcon * (hh + qq - qpt) )
         else
            til   = cpdryi * thil * exner
            !----- Decide how to compute, based on temperature. ---------------------------!
            if (temp > ttripoli) then
               dtempdrs = - til * qq / ( rcon * cpdry * (2.*temp-til))
            else
               dtempdrs = - til * qq * htripolii / rcon
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function dtempdrs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the ice-vapour equivalent potential temperature from       !
   ! theta_iland the total mixing ratio.  This is equivalent to the equivalent potential   !
   ! temperature considering also the effects of fusion/melting/sublimation.               !
   !     In case you want to find thetae (i.e. without ice) simply set the the logical     !
   ! useice to .false. .                                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function thetaeiv(thil,pres,temp,rvap,rtot,useice)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: thil   ! Ice-liquid potential temp.    [     K]
      real(kind=4), intent(in)           :: pres   ! Pressure                      [    Pa]
      real(kind=4), intent(in)           :: temp   ! Temperature                   [     K]
      real(kind=4), intent(in)           :: rvap   ! Water vapour mixing ratio     [ kg/kg]
      real(kind=4), intent(in)           :: rtot   ! Total mixing ratio            [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! Should I use ice?             [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: tlcl   ! Internal LCL temperature      [     K]
      real(kind=4)                       :: plcl   ! Lifting condensation pressure [    Pa]
      real(kind=4)                       :: dzlcl  ! Thickness of lyr. beneath LCL [     m]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the liquid condensation level (LCL).                                      !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         call lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      else
         call lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The definition of the thetae_iv is the thetae_ivs at the LCL. The LCL, in turn !
      ! is the point in which rtot = rvap = rsat, so at the LCL rliq = rice = 0.           !
      !------------------------------------------------------------------------------------!
      thetaeiv  = thetaeivs(thil,tlcl,rtot,0.,0.)
      !------------------------------------------------------------------------------------!

      return
   end function thetaeiv
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
   !                 pressure at the LCL) is a function of temperature.  In case you want  !
   !                 d(Thetae_ivs)/dT, use the dthetaeivs_dt function instead.             !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function dthetaeiv_dtlcl(theiv,tlcl,rtot,eslcl,useice)
      use rconstants , only : rocp       & ! intent(in)
                            , cpdry      & ! intent(in)
                            , dcpvl      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: theiv    ! Ice-vap. equiv. pot. temp. [      K]
      real(kind=4), intent(in)           :: tlcl     ! LCL temperature            [      K]
      real(kind=4), intent(in)           :: rtot     ! Total mixing ratio         [  kg/kg]
      real(kind=4), intent(in)           :: eslcl    ! LCL sat. vapour pressure   [     Pa]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice   ! Flag for considering ice   [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: desdtlcl ! Sat. vapour pres. deriv.   [   Pa/K]
      real(kind=4)                       :: esterm   ! es(TLC) term               [   ----]
      real(kind=4)                       :: hhlcl    ! Enthalpy -- sensible       [   J/kg]
      real(kind=4)                       :: qqlcl    ! Enthalpy -- latent         [   J/kg]
      real(kind=4)                       :: qptlcl   ! Latent deriv. * T_LCL      [   J/kg]
      !------------------------------------------------------------------------------------!



      !----- Find the derivative of rs with temperature. ----------------------------------!
      if (present(useice)) then
         desdtlcl = eslifp(tlcl,useice)
      else
         desdtlcl = eslifp(tlcl)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Saturation term.                                                               !
      !------------------------------------------------------------------------------------!
      esterm = rocp * tlcl * desdtlcl / eslcl
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the enthalpy terms.                                                       !
      !------------------------------------------------------------------------------------!
      hhlcl  = cpdry * tlcl
      qqlcl  = alvl(tlcl) * rtot
      qptlcl = dcpvl * rtot * tlcl
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Derivative.                                                                   !
      !------------------------------------------------------------------------------------!
      dthetaeiv_dtlcl = theiv / tlcl * (1. - esterm - (qqlcl - qptlcl) / hhlcl)
      !------------------------------------------------------------------------------------!

      return
   end function dthetaeiv_dtlcl
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
   real(kind=4) function thetaeivs(thil,temp,rsat,rliq,rice)
      use rconstants , only : cpdry    ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: thil  ! Theta_il, ice-liquid water pot. temp.    [     K]
      real(kind=4), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=4), intent(in) :: rsat  ! Saturation water vapour mixing ratio     [ kg/kg]
      real(kind=4), intent(in) :: rliq  ! Liquid water mixing ratio                [ kg/kg]
      real(kind=4), intent(in) :: rice  ! Ice mixing ratio                         [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: rtots ! Saturated mixing ratio                   [     K]
      !------------------------------------------------------------------------------------!


      !------ Find the total saturation mixing ratio. -------------------------------------!
      rtots = rsat+rliq+rice
      !------------------------------------------------------------------------------------!


      !------ Find the saturation equivalent potential temperature. -----------------------!
      thetaeivs = thil * exp ( alvl(temp) * rtots / (cpdry * temp))
      !------------------------------------------------------------------------------------!

      return
   end function thetaeivs
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
   real(kind=4) function dthetaeivs_dt(theivs,temp,pres,rsat,useice)
      use rconstants , only : cpdry     & ! intent(in)
                            , dcpvl     ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: theivs ! Sat. ice-vap. eq. pot. temp. [      K]
      real(kind=4), intent(in)           :: temp   ! Temperature                  [      K]
      real(kind=4), intent(in)           :: pres   ! Pressure                     [     Pa]
      real(kind=4), intent(in)           :: rsat   ! Saturation mixing ratio      [  kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! Flag for considering ice     [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: drsdt  ! Sat. mixing ratio derivative [kg/kg/K]
      real(kind=4)                       :: hh     ! Enthalpy -- sensible         [   J/kg]
      real(kind=4)                       :: qqaux  ! Enthalpy -- sensible         [   J/kg]
      !------------------------------------------------------------------------------------!


      !----- Find the derivative of rs with temperature and associated term. --------------!
      if (present(useice)) then
         drsdt = rslifp(pres,temp,useice)
      else
         drsdt = rslifp(pres,temp)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the enthalpy terms.                                                      !
      !------------------------------------------------------------------------------------!
      hh    = cpdry * temp
      qqaux = alvl(temp) * (drsdt * temp - rsat) + dcpvl * rsat * temp
      !------------------------------------------------------------------------------------!


      !----- Find the derivative.  Depending on the temperature, use different eqn. -------!
      dthetaeivs_dt = theivs / temp * ( 1. + qqaux / hh )
      !------------------------------------------------------------------------------------!

      return
   end function dthetaeivs_dt
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
   real(kind=4) function thetaeiv2thil(theiv,pres,rtot,useice)
      use rconstants , only : ep       & ! intent(in)
                            , cpdry    & ! intent(in)
                            , p00      & ! intent(in)
                            , rocp     & ! intent(in)
                            , t3ple    & ! intent(in)
                            , t00      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: theiv     ! Ice vap. equiv. pot. temp. [     K]
      real(kind=4), intent(in)           :: pres      ! Pressure                   [    Pa]
      real(kind=4), intent(in)           :: rtot      ! Total mixing ratio         [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice    ! May I use ice thermodyn.   [   T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=4)                       :: pvap      ! Sat. vapour pressure
      real(kind=4)                       :: theta     ! Potential temperature
      real(kind=4)                       :: deriv     ! Function derivative 
      real(kind=4)                       :: funnow    ! Function for which we seek a root.
      real(kind=4)                       :: funa      ! Smallest guess function
      real(kind=4)                       :: funz      ! Largest guess function
      real(kind=4)                       :: tlcla     ! Smallest guess (Newton: old guess)
      real(kind=4)                       :: tlclz     ! Largest guess (Newton: new guess)
      real(kind=4)                       :: tlcl      ! What will be the LCL temperature
      real(kind=4)                       :: es00      ! Defined as p00*rt/(epsilon + rt)
      real(kind=4)                       :: delta     ! Aux. variable (For 2nd guess).
      integer                            :: itn       ! Iteration counters
      integer                            :: itb       ! Iteration counters
      integer                            :: ii        ! Another counter
      logical                            :: converged ! Convergence handle
      logical                            :: zside     ! Side checker for Regula Falsi
      logical                            :: frozen    ! Will use ice thermodynamics
      !----- Local constants. -------------------------------------------------------------!
      logical     , parameter            :: debug = .false.
      !------------------------------------------------------------------------------------!



      !----- Fill the flag for ice thermodynamics so it will be present. ------------------!
      if (present(useice)) then
         frozen = useice
      else
         frozen = bulk_on
      end if
      !------------------------------------------------------------------------------------!



      !----- Find es00, which is a constant. ----------------------------------------------!
      es00 = p00 * rtot / (ep+rtot)
      !------------------------------------------------------------------------------------!


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write (unit=36,fmt='(a)') '------------------------------------------------------'
         write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x))')                               &
            'INPUT : it=',-1,'theiv=',theiv,'pres=',0.01*pres,'rtot=',rtot*1000.
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, we assume we are lucky and right at the LCL.                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rtot / (ep + rtot) 
      tlclz     = tslif(pvap,frozen)
      theta     = tlclz * (es00 / pvap) ** rocp
      funnow    = thetaeivs(theta,tlclz,rtot,0.,0.)
      deriv     = dthetaeiv_dtlcl(funnow,tlclz,rtot,pvap,frozen)
      funnow    = funnow - theiv
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
             'NEWTON: it=',0,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow            &
            ,'deriv=',deriv
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Put something in tlcla in case we never loop through Newton's method. --------!
      tlcla     = tlclz
      funa      = funnow
      converged = .false.
      !------------------------------------------------------------------------------------!

      !----- Looping: Newton's iterative method -------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Update guesses. -----------------------------------------------------------!
         tlcla   = tlclz
         funa    = funnow
         tlclz   = tlcla - funnow/deriv
         !---------------------------------------------------------------------------------!


         !----- Update the function evaluation and its derivative. ------------------------!
         pvap    = eslif(tlclz,frozen)
         theta   = tlclz * (es00/pvap)**rocp
         funnow  = thetaeivs(theta,tlclz,rtot,0.,0.)
         deriv   = dthetaeiv_dtlcl(funnow,tlclz,rtot,pvap,frozen)
         funnow  = funnow - theiv
         !---------------------------------------------------------------------------------!

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (debug) then
            write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')          &
               'NEWTON: it=',itn,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow        &
                      ,'deriv=',deriv
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         converged = abs(tlcla-tlclz) < toler * tlclz
         if (funnow == 0.) then
            tlcl = tlclz
            funz = funnow
            converged = .true.
            exit newloop
         elseif (converged) then
            tlcl = 0.5*(tlcla+tlclz)
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
         funz = funnow
         zside=.true.
         if (funa*funnow > 0.) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler*tlcla) then
               delta = 100.*toler*tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funz-funa)),100.*toler*tlcla)
            end if
            tlclz = tlcla + delta

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (debug) then
               write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')       &
                  '2NGGSS: tt=',0,'tlclz=',tlclz-t00,'tlcla=',tlcla-t00,'pvap=',0.01*pvap  &
                          ,'funa=',funa,'funz=',funnow,'delta=',delta
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlcla + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif(tlclz,frozen)
               theta = tlclz * (es00/pvap)**rocp
               funz  = thetaeivs(theta,tlclz,rtot,0.,0.) - theiv

               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               if (debug) then
                  write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')    &
                     '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz    &
                             ,'delta=',delta
               end if
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

               zside = funa*funz < 0.
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
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thetaeiv2thil','therm_lib.f90')
            end if
         end if
         !---- Continue iterative method. -------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo

            !----- Update the guess. ------------------------------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)

            !----- Updating function evaluation -------------------------------------------!
            pvap   = eslif(tlcl,frozen)
            theta  = tlcl * (es00/pvap)**rocp
            funnow = thetaeivs(theta,tlcl,rtot,0.,0.) - theiv

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (debug) then
               write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),1(a,1x,es11.4,1x))')       &
                  'REGFAL: it=',itb,'tlcl =',tlcl -t00,'pvap=',0.01*pvap,'fun=',funnow
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Check for convergence. If it did, return, we found the solution.          !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla) < toler * tlcl
            if (converged) then
               exit fpoloop
            elseif (funnow*funa < 0.) then 
               tlclz = tlcl
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 0.5
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tlcla = tlcl
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 0.5
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
         end do fpoloop
      end if

      if (converged) then 
         thetaeiv2thil  = theiv * exp (- alvl(tlcl) * rtot / (cpdry * tlcl) )
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (debug) then
            write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')          &
               'ANSWER: itb=',itn,'tlcl=',tlcl-t00,'eslcl=',0.01*pvap                      &
                      ,'thil=',thetaeiv2thil,'funa=',funa,'funz=',funz
            write (unit=36,fmt='(a)') '---------------------------------------------------'
            write (unit=36,fmt='(a)') ' '
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      else
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
         write (unit=*,fmt='(a)')           ' THEIV2THIL failed!'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Input: '
         write (unit=*,fmt='(a,1x,f12.5)')  '    THEIV    [     K]:',theiv
         write (unit=*,fmt='(a,1x,f12.5)')  '    PRES     [    Pa]:',pres * 100.
         write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1000.*rtot
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Output: '
         write (unit=*,fmt='(a,1x,i12)')    '    ITERATIONS       :',itb
         write (unit=*,fmt='(a,1x,f12.5)')  '    PVAP     [   hPa]:',pvap
         write (unit=*,fmt='(a,1x,f12.5)')  '    THETA    [     K]:',theta
         write (unit=*,fmt='(a,1x,f12.5)')  '    TLCL     [    °C]:',tlcl-t00
         write (unit=*,fmt='(a,1x,f12.5)')  '    TLCLA    [    °C]:',tlcla-t00
         write (unit=*,fmt='(a,1x,f12.5)')  '    TLCLZ    [    °C]:',tlclz-t00
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNA     [     K]:',funa
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNZ     [     K]:',funz
         write (unit=*,fmt='(a,1x,es12.5)') '    DERIV    [   ---]:',deriv
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_A    [   ---]:',abs(tlcl-tlcla)/tlcl
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_Z    [   ---]:',abs(tlcl-tlclz)/tlcl
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)

         call abort_run  ('TLCL didn''t converge, qgave up!'                               &
                         ,'thetaeiv2thil','therm_lib.f90')
      end if

      return
   end function thetaeiv2thil
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
   subroutine thetaeivs2temp(theivs,pres,theta,temp,rsat,useice)
      use rconstants , only : cpdry     & ! intent(in)
                            , ep        & ! intent(in)
                            , p00       & ! intent(in)
                            , rocp      & ! intent(in)
                            , t00       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)            :: theivs     ! Sat. thetae_iv          [      K]
      real(kind=4), intent(in)            :: pres       ! Pressure                [     Pa]
      real(kind=4), intent(out)           :: theta      ! Potential temperature   [      K]
      real(kind=4), intent(out)           :: temp       ! Temperature             [      K]
      real(kind=4), intent(out)           :: rsat       ! Saturation mixing ratio [  kg/kg]
      logical     , intent(in) , optional :: useice     ! May use ice thermodyn.  [    T|F]
      !----- Local variables, with other thermodynamic properties -------------------------!
      real(kind=4)                        :: exnernormi ! 1./ (Norm. Exner func.) [    ---]
      logical                             :: frozen     ! Will use ice thermodyn. [    T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=4)                        :: deriv      ! Function derivative 
      real(kind=4)                        :: funnow     ! Current function evaluation
      real(kind=4)                        :: funa       ! Smallest guess function
      real(kind=4)                        :: funz       ! Largest guess function
      real(kind=4)                        :: tempa      ! Smallest guess (Newton: previous)
      real(kind=4)                        :: tempz      ! Largest guess (Newton: new)
      real(kind=4)                        :: delta      ! Aux. variable for 2nd guess.
      integer                             :: itn        ! Iteration counter
      integer                             :: itb        ! Iteration counter
      logical                             :: converged  ! Convergence handle
      logical                             :: zside      ! Flag for side check.
      !------------------------------------------------------------------------------------!


      !----- Set up the ice check, in case useice is not present. -------------------------!
      if (present(useice)) then
         frozen = useice
      else 
         frozen = bulk_on
      end if
      !------------------------------------------------------------------------------------!



      !----- Finding the inverse of normalised Exner, which is constant in this routine ---!
      exnernormi = (p00 /pres) ** rocp
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The 1st. guess, no idea, guess 0°C.                                            !
      !------------------------------------------------------------------------------------!
      tempz     = t00
      theta     = tempz * exnernormi
      rsat      = rslif(pres,tempz,frozen)
      funnow    = thetaeivs(theta,tempz,rsat,0.,0.)
      deriv     = dthetaeivs_dt(funnow,tempz,pres,rsat,frozen)
      funnow    = funnow - theivs
      !------------------------------------------------------------------------------------!


      !----- Copy here just in case Newton is aborted at the 1st guess. -------------------!
      tempa     = tempz
      funa      = funnow
      !------------------------------------------------------------------------------------!

      converged = .false.
      !----- Newton's method loop. --------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Updating guesses ----------------------------------------------------------!
         tempa   = tempz
         funa    = funnow
         
         tempz   = tempa - funnow/deriv
         theta   = tempz * exnernormi
         rsat    = rslif(pres,tempz,frozen)
         funnow  = thetaeivs(theta,tempz,rsat,0.,0.)
         deriv   = dthetaeivs_dt(funnow,tempz,pres,rsat,frozen)
         funnow  = funnow - theivs

         converged = abs(tempa-tempz) < toler*tempz
         if (funnow == 0.) then
            converged =.true.
            temp = tempz
            exit newloop
         elseif (converged) then
            temp = 0.5*(tempa+tempz)
            exit newloop
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If we have reached this point then it's because Newton's method failed.  Use   !
      ! bisection instead.                                                                 !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         !----- Set funz, and check whether funa and funz already have opposite sign. -----!
         funz  = funnow
         zside = .false.
         !---------------------------------------------------------------------------------!

         if (funa*funnow > 0.) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler*tempa) then
               delta = 100.*toler*tempa
            else
               delta = max(abs(funa*(tempz-tempa)/(funz-funa)),100.*toler*tempa)
            end if
            !------------------------------------------------------------------------------!

            tempz = tempa + delta
            zgssloop: do itb=1,maxfpo
               !----- So this will be +1 -1 +2 -2 etc. ------------------------------------!
               tempz = tempz + real((-1)**itb * (itb+3)/2) * delta
               theta = tempz * exnernormi
               rsat  = rslif(pres,tempz,frozen)
               funz  = thetaeivs(theta,tempz,rsat,0.,0.) - theivs
               zside = funa*funz < 0.
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thetaes2temp','therm_lib.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            if (abs(funz-funa) < toler*tempa) then
               temp   = 0.5*(tempa+tempz)
            else
               temp   = (funz*tempa-funa*tempz)/(funz-funa)
            end if
            theta  = temp * exnernormi
            rsat   = rslif(pres,temp,frozen)
            funnow = thetaeivs(theta,temp,rsat,0.,0.) - theivs

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa) < toler*temp
            if (converged) then
               exit fpoloop
            elseif (funnow*funa < 0.) then 
               tempz = temp
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 0.5
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tempa = temp
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not. zside) funz = funz * 0.5
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
         end do fpoloop
      end if

      if (converged) then 
         !----- Compute theta and rsat with temp just for consistency ---------------------!
         theta = temp * exnernormi
         rsat  = rslif(pres,temp,frozen)
      else
         call abort_run  ('Temperature didn''t converge, I gave up!'                       &
                         ,'thetaes2temp','therm_lib.f90')
      end if

      return
   end subroutine thetaeivs2temp
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
   subroutine lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      use rconstants , only : cpog     & ! intent(in)
                            , ep       & ! intent(in)
                            , p00      & ! intent(in)
                            , rocp     & ! intent(in)
                            , t3ple    & ! intent(in)
                            , t00      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)            :: thil      ! Ice liquid pot. temp. (*)[      K]
      real(kind=4), intent(in)            :: pres      ! Pressure                 [     Pa]
      real(kind=4), intent(in)            :: temp      ! Temperature              [      K]
      real(kind=4), intent(in)            :: rtot      ! Total mixing ratio       [  kg/kg]
      real(kind=4), intent(in)            :: rvap      ! Vapour mixing ratio      [  kg/kg]
      real(kind=4), intent(out)           :: tlcl      ! LCL temperature          [      K]
      real(kind=4), intent(out)           :: plcl      ! LCL pressure             [     Pa]
      real(kind=4), intent(out)           :: dzlcl     ! Sub-LCL layer thickness  [      m]
      !------------------------------------------------------------------------------------!
      ! (*) This is the most general variable. Thil is exactly theta for no condensation   !
      !     condition, and it is the liquid potential temperature if no ice is present.    !
      !------------------------------------------------------------------------------------!
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in) , optional :: useice    ! May use ice thermodyn.?  [    T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                        :: pvap      ! Sat. vapour pressure
      real(kind=4)                        :: deriv     ! Function derivative 
      real(kind=4)                        :: funnow    ! Current function evaluation
      real(kind=4)                        :: funa      ! Smallest guess function
      real(kind=4)                        :: funz      ! Largest  guess function
      real(kind=4)                        :: tlcla     ! Smallest guess (Newton: previous)
      real(kind=4)                        :: tlclz     ! Largest guess (Newton: new)
      real(kind=4)                        :: es00      ! Defined as p00*rt/(epsilon + rt)
      real(kind=4)                        :: delta     ! Aux. variable for bisection
      integer                             :: itn       ! Iteration counter
      integer                             :: itb       ! Iteration counter
      logical                             :: converged ! Convergence flag
      logical                             :: zside     ! Flag to check sides
      logical                             :: frozen    ! Will use ice thermodyn.  [    T|F]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether ice thermodynamics is the way to go.                             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice
      else 
         frozen = bulk_on
      end if
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a)') '----------------------------------------------------------'
      !write (unit=21,fmt='(a,1x,i5,1x,5(a,1x,f11.4,1x))')                                  &
      !   'INPUT : it=',-1,'thil=',thil,'pres=',0.01*pres,'temp=',temp-t00                  &
      !        ,'rvap=',rvap*1000.,'rtot=',rtot*1000.
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !----- Find es00, which is a constant. ----------------------------------------------!
      es00 = p00 * rtot / (ep+rtot)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvap / (ep + rvap)
      tlclz     = 55. + 2840. / (3.5 * log(temp) - log(0.01*pvap) - 4.805)
      pvap      = eslif(tlclz,frozen)
      funnow    = tlclz * (es00/pvap)**rocp - thil
      deriv     = (funnow+thil)*(1./tlclz - rocp*eslifp(tlclz,frozen)/pvap) 
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')                &
      !   'NEWTON: it=',0,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      tlcla     = tlclz
      funa      = funnow


      !------------------------------------------------------------------------------------!
      !      First loop: Newton's method.                                                  !
      !------------------------------------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         !----- If derivative is too small, skip Newton's and try bisection instead. ------!
         if (abs(deriv) < toler) exit newloop
         !---------------------------------------------------------------------------------!


         !----- Otherwise, update guesses. ------------------------------------------------!
         tlcla  = tlclz
         funa   = funnow
         tlclz  = tlcla - funnow/deriv
         pvap   = eslif(tlclz,frozen)
         funnow = tlclz * (es00/pvap)**rocp - thil
         deriv  = (funnow+thil)*(1./tlclz - rocp*eslifp(tlclz,frozen)/pvap)

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
         !   'NEWTON: it=',itn,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow           &
         !          ,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !---------------------------------------------------------------------------------!
         !      Check for convergence.                                                     !
         !---------------------------------------------------------------------------------!
         converged = abs(tlcla-tlclz) < toler*tlclz
         if (converged) then
            !----- Guesses are almost identical, average them. ----------------------------!
            tlcl = 0.5*(tlcla+tlclz)
            funz = funnow
            exit newloop
            !------------------------------------------------------------------------------!
         elseif (funnow == 0.) then
            !----- We've hit the answer by luck, copy the answer. -------------------------!
            tlcl = tlclz
            funz = funnow
            converged = .true.
            exit newloop
            !------------------------------------------------------------------------------!
         end if
      end do newloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Check whether Newton's method has converged.                                  !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         !---------------------------------------------------------------------------------!
         !     Newton's method has failed.  We use regula falsi instead.  First, we must   !
         ! find two guesses whose function evaluations have opposite signs.                !
         !---------------------------------------------------------------------------------!
         if (funa*funnow < 0. ) then
            !----- We already have two good guesses. --------------------------------------!
            funz  = funnow
            zside = .true.
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     We need to find another guess with opposite sign.                        !
            !------------------------------------------------------------------------------!

            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funnow-funa) < toler*tlcla) then
               delta = 100.*toler*tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funnow-funa)),100.*toler*tlcla)
            end if
            tlclz = tlcla + delta
            !------------------------------------------------------------------------------!

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
               tlclz = tlcla + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif(tlclz,frozen)
               funz  = tlclz * (es00/pvap)**rocp - thil
               !---------------------------------------------------------------------------!


               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')       &
               !   '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz       &
               !           ,'delta=',delta
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               zside = funa*funz < 0.0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(a)') ' ====== No second guess for you... ======'
               write (unit=*,fmt='(a)') ' + INPUT variables: '
               write (unit=*,fmt='(a,1x,es14.7)') 'THIL =',thil
               write (unit=*,fmt='(a,1x,es14.7)') 'TEMP =',temp
               write (unit=*,fmt='(a,1x,es14.7)') 'PRES =',pres
               write (unit=*,fmt='(a,1x,es14.7)') 'RTOT =',rtot
               write (unit=*,fmt='(a,1x,es14.7)') 'RVAP =',rvap
               write (unit=*,fmt='(a)') ' ============ Failed guess... ==========='
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLA =',tlcla,'FUNA =',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLZ =',tlclz,'FUNC =',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'DELTA =',delta,'FUNN =',funnow
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                             ,'lcl_il','therm_lib.f90')
            end if
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      We have the guesses, solve the regula falsi method.                        !
         !---------------------------------------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            !----- Update guess and function evaluation. ----------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)
            pvap   = eslif(tlcl,frozen)
            funnow = tlcl * (es00/pvap)**rocp - thil
            !------------------------------------------------------------------------------!

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),1(a,1x,es11.4,1x))')          &
            !   'REGFAL: it=',itb,'tlcl =',tlcl -t00,'pvap=',0.01*pvap,'fun=',funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Check for convergence.  If it did, return, we found the solution.         !
            ! Otherwise, we update one of the guesses.                                     !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla) < toler*tlcl .and.  abs(tlcl-tlclz) < toler*tlcl
            if (funnow == 0. .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.) then 
               tlclz = tlcl
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 0.5
               !---------------------------------------------------------------------------!


               !----- We have just updated zside, sett zside to true. ---------------------!
               zside = .true.
               !---------------------------------------------------------------------------!
            else
               tlcla = tlcl
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 0.5
               !---------------------------------------------------------------------------!


               !----- We have just updated aside, set zside to false. ---------------------!
               zside = .false. 
               !---------------------------------------------------------------------------!
            end if
         end do fpoloop
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check whether we have succeeded or not.                                        !
      !------------------------------------------------------------------------------------!
      if (converged) then 
         !----- We have found a solution, find the remaining LCL properties. --------------!
         pvap  = eslif(tlcl,frozen)
         plcl  = (ep + rvap) * pvap / rvap
         dzlcl = max(cpog*(temp-tlcl),0.)
         !---------------------------------------------------------------------------------!


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
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Pressure        [   hPa] =',0.01*pres
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Temperature     [    °C] =',temp-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rtot            [  g/kg] =',1000.*rtot
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rvap            [  g/kg] =',1000.*rvap
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome.'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcla           [    °C] =',tlcla-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlclz           [    °C] =',tlclz-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [     K] =',funnow
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [     K] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [     K] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [  ----] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [  ----] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                                            ,abs(tlclz-tlcla)/tlclz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcl            [    °C] =',tlcl
         call abort_run  ('TLCL didn''t converge, gave up!','lcl_il','therm_lib.f90')
      end if
      return
   end subroutine lcl_il
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes a consistent set of temperature and condensated phases   !
   ! mixing ratio for a given theta_il, Exner function, and total mixing ratio. This is    !
   ! very similar to the function thil2temp, except that now we don't know rliq and rice,  !
   ! and for this reason they also become functions of temperature, since they are defined !
   ! as rtot-rsat(T,p), remembering that rtot and p are known. If the air is not           !
   ! saturated, we rather use the fact that theta_il = theta and skip the hassle.          !
   ! Otherwise, we use iterative methods.  We will always try Newton's method, since it    !
   ! converges fast. The caveat is that Newton may fail, and it actually does fail very    !
   ! close to the triple point, because the saturation vapour pressure function has a      !
   ! "kink" at the triple point (continuous, but not differentiable).  If that's the case, !
   ! then we fall back to a modified regula falsi (Illinois) method, which is a mix of     !
   ! secant and bisection and will converge.                                               !
   !---------------------------------------------------------------------------------------!
   subroutine thil2tqall(thil,exner,pres,rtot,rliq,rice,temp,rvap,rsat)
      use rconstants , only : cpdry    & ! intent(in)
                            , cpdryi   & ! intent(in)
                            , t00      & ! intent(in)
                            , toodry   & ! intent(in)
                            , t3ple    ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in)    :: thil      ! Ice-liquid water potential temp.  [     K]
      real(kind=4), intent(in)    :: exner     ! Exner function                    [J/kg/K]
      real(kind=4), intent(in)    :: pres      ! Pressure                          [    Pa]
      real(kind=4), intent(in)    :: rtot      ! Total mixing ratio                [ kg/kg]
      real(kind=4), intent(out)   :: rliq      ! Liquid water mixing ratio         [ kg/kg]
      real(kind=4), intent(out)   :: rice      ! Ice mixing ratio                  [ kg/kg]
      real(kind=4), intent(inout) :: temp      ! Temperature                       [     K]
      real(kind=4), intent(out)   :: rvap      ! Water vapour mixing ratio         [ kg/kg]
      real(kind=4), intent(out)   :: rsat      ! Sat. water vapour mixing ratio    [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                :: tempa     ! Lower bound for regula falsi iteration
      real(kind=4)                :: tempz     ! Upper bound for regula falsi iteration
      real(kind=4)                :: t1stguess ! Book keeping temperature 1st guess 
      real(kind=4)                :: fun1st    ! Book keeping 1st guess function
      real(kind=4)                :: funa      ! Function evaluation at tempa
      real(kind=4)                :: funz      ! Function evaluation at tempz
      real(kind=4)                :: funnow    ! Function at this iteration.
      real(kind=4)                :: delta     ! Aux. var in case we need regula falsi.
      real(kind=4)                :: deriv     ! Derivative of this function.
      integer                     :: itn       ! Iteration counter
      integer                     :: itb       ! Iteration counter
      integer                     :: ii        ! Iteration counter
      logical                     :: converged ! Convergence handle
      logical                     :: zside     ! Aux. Flag, for two purposes:
                                               ! 1. Found a 2nd guess for regula falsi.
                                               ! 2. I retained the "zside" (T/F)
      !----- Local constants. -------------------------------------------------------------!
      logical     , parameter     :: debug = .false.
      !------------------------------------------------------------------------------------!

      t1stguess = temp

      !------------------------------------------------------------------------------------!
      !      First check: try to find temperature assuming sub-saturation and check if     !
      ! this is the case. If it is, then there is no need to go through the iterative      !
      ! loop.                                                                              !
      !------------------------------------------------------------------------------------!
      tempz  = cpdryi * thil * exner
      rsat   = max(toodry,rslif(pres,tempz))
      if (tempz >= t3ple) then
         rliq = max(0.,rtot-rsat)
         rice = 0.
      else
         rice = max(0.,rtot-rsat)
         rliq = 0.
      end if
      rvap = rtot-rliq-rice
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    If rtot < rsat, this is not saturated, we can leave the subroutine and bypass   !
      ! the iterative part.                                                                !
      !------------------------------------------------------------------------------------!
      if (rtot < rsat) then
         temp = tempz
         return
      end if

      !------------------------------------------------------------------------------------!
      !   If not, then use the temperature the user gave as first guess and solve          !
      ! iteratively.  We use the user instead of what we just found because if the air is  !
      ! saturated, then this can be too far off which may be bad for Newton's method.      !
      !------------------------------------------------------------------------------------!
      tempz = temp
      rsat   = max(toodry,rslif(pres,tempz))
      if (tempz >= t3ple) then
         rliq = max(0.,rtot-rsat)
         rice = 0.
      else
         rice = max(0.,rtot-rsat)
         rliq = 0.
      end if
      rvap = rtot-rliq-rice


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write (unit=46,fmt='(a)') '------------------------------------------------------'
         write (unit=46,fmt='(a,1x,i5,1x,5(a,1x,f11.4,1x))')                               &
            'INPUT: it=',-1,'thil=',thil,'exner=',exner,'press=',0.01*pres                 &
           ,'rtot=',1000.*rtot,'t1st=',temp-t00
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !------------------------------------------------------------------------------------!
      !     Find the function. We are seeking a temperature which is associated with the   !
      ! theta_il we provided. Thus, the function is simply the difference between the      !
      ! theta_il associated with our guess and the actual theta_il.                        !
      !------------------------------------------------------------------------------------!
      funnow = theta_iceliq(exner,tempz,rliq,rice)
      !----- Updating the derivative. -----------------------------------------------------!
      deriv  = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq,rice)
      funnow = funnow - thil
      fun1st = funnow
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x),a,1x,es11.4,1x)')                &
            'NEWTON: it=',0,'temp=',tempz-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq        &
           ,'rice=',1000.*rice,'rvap=',1000.*rvap,'fun=',funnow,'deriv=',deriv
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !------------------------------------------------------------------------------------!
      !    Now we enter at the Newton's method iterative loop. We are always going to try  !
      ! this first, because it's fast, but if it turns out to be a dangerous choice or if  !
      ! it doesn't converge fast, we will fall back to regula falsi.                       !
      !    We start by initialising the flag and copying temp to tempz, the newest guess.  !
      !------------------------------------------------------------------------------------!
      converged=.false.
      newloop: do itn=1,maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Saving previous guess. We also save the function is in case we give up on   !
         ! Newton's and switch to regula falsi.                                            !
         !---------------------------------------------------------------------------------!
         funa  = funnow
         tempa = tempz

         !----- Go to bisection if the derivative is too flat (too dangerous...) ----------!
         if (abs(deriv) < toler) exit newloop

         tempz = tempa - funnow / deriv

         !----- Finding the mixing ratios associated with this guess ----------------------!
         rsat  = max(toodry,rslif(pres,tempz))
         if (tempz >= t3ple) then
            rliq = max(0.,rtot-rsat)
            rice = 0.
         else
            rice = max(0.,rtot-rsat)
            rliq = 0.
         end if
         rvap = rtot-rliq-rice

         !----- Updating the function -----------------------------------------------------!
         funnow = theta_iceliq(exner,tempz,rliq,rice)
         !----- Updating the derivative. --------------------------------------------------!
         deriv  = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq,rice)
         funnow = funnow - thil

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (debug) then
            write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x),a,1x,es11.4,1x)')             &
               'NEWTON: it=',itn,'temp=',tempz-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq   &
              ,'rice=',1000.*rice,'rvap=',1000.*rvap,'fun=',funnow,'deriv=',deriv
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         
         converged = abs(tempa-tempz) < toler*tempz
         !---------------------------------------------------------------------------------!
         !   Convergence. The temperature will be the mid-point between tempa and tempz.   !
         ! Fix the mixing ratios and return. But first check for converged due to luck. If !
         ! the guess gives a root, then that's it. It looks unlikely, but it actually      !
         ! happens sometimes and if not checked it becomes a singularity.                  !
         !---------------------------------------------------------------------------------!
         if (funnow == 0.) then
            temp = tempz
            converged = .true.
            exit newloop
         elseif (converged) then
            temp = 0.5 * (tempa+tempz)
            rsat  = max(toodry,rslif(pres,temp))
            if (temp >= t3ple) then
               rliq = max(0.,rtot-rsat)
               rice = 0.
            else
               rice = max(0.,rtot-rsat)
               rliq = 0.
            end if
            rvap = rtot-rliq-rice
            exit newloop
         end if

      end do newloop
      !------------------------------------------------------------------------------------!

      !----- For debugging only -----------------------------------------------------------!
      itb = itn+1

      if (.not. converged) then
         !---------------------------------------------------------------------------------!
         !    If I reach this point, then it means that Newton's method failed finding the !
         ! equilibrium, so we are going to use the regula falsi instead.  If Newton's      !
         ! method didn't converge, we use tempa as one guess and now we seek a tempz with  !
         ! opposite sign.                                                                  !
         !---------------------------------------------------------------------------------!
         !----- Check funa and funnow have opposite signs. If so, we are ready to go ------!
         if (funa*funnow < 0.0) then
            funz  = funnow
            zside = .true.
         !----- Otherwise, checking whether the 1st guess had opposite sign. --------------!
         elseif (funa*fun1st < 0.0) then
            funz  = fun1st
            zside = .true.
         !---------------------------------------------------------------------------------!
         !    Looking for a guess. Extrapolate funa linearly, trying to get the -funa.  We !
         ! don't need it to be funa, just with the opposite sign.  If that's not enough,   !
         ! we keep going further... Force the guesses to be at least 1K apart              !
         !---------------------------------------------------------------------------------!
         else
            if (abs(funnow-funa) < 100.*toler*tempa) then
               delta = 0.5
            else
               delta = max(abs(funa)*abs((tempz-tempa)/(funnow-funa)),0.5)
            end if
            tempz = tempa + delta
            funz  = funa
            !----- Just to enter at least once. The 1st time tempz=tempa-2*delta ----------!
            zside = .false. 
            zgssloop: do itb=1,maxfpo
                tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
                rsat   = max(toodry,rslif(pres,tempz))
                if (tempz >= t3ple) then
                   rliq = max(0.,rtot-rsat)
                   rice = 0.
                else
                   rice = max(0.,rtot-rsat)
                   rliq = 0.
                end if
                rvap = rtot-rliq-rice
                funz = theta_iceliq(exner,tempz,rliq,rice) - thil
                zside = funa*funz < 0.0
                if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
               write (unit=*,fmt='(a)')           ' THIL2TQALL: NO SECOND GUESS FOR YOU!'
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           ' -> Input: '
               write (unit=*,fmt='(a,1x,f12.5)')  '    THETA_IL [     K]:',thil
               write (unit=*,fmt='(a,1x,f12.5)')  '    PRESS    [   hPa]:',0.01*pres
               write (unit=*,fmt='(a,1x,f12.5)')  '    EXNER    [J/kg/K]:',exner
               write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1000.*rtot
               write (unit=*,fmt='(a,1x,f12.5)')  '    T1ST     [  degC]:',t1stguess-t00
               write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPA    [  degC]:',tempa-t00
               write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPZ    [  degC]:',tempz-t00
               write (unit=*,fmt='(a,1x,f12.5)')  '    FUNNOW   [     K]:',funnow
               write (unit=*,fmt='(a,1x,f12.5)')  '    FUNA     [     K]:',funa
               write (unit=*,fmt='(a,1x,f12.5)')  '    FUNZ     [     K]:',funz
               write (unit=*,fmt='(a,1x,f12.5)')  '    DELTA    [     K]:',delta
               write (unit=*,fmt='(60a1)')        ('-',ii=1,60)

               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thil2tqall','therm_lib.f90')
            end if
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Now we loop until convergence is achieved.  One important thing to notice   !
         ! is that Newton's method fail only when T is almost T3ple, which means that ice  !
         ! and liquid should be present, and we are trying to find the saturation point    !
         ! with all ice or all liquid. This will converge but the final answer will        !
         ! contain significant error. To reduce it we redistribute the condensates between !
         ! ice and liquid conserving the total condensed mixing ratio.                     !
         !---------------------------------------------------------------------------------!
         fpoloop: do itb=itn,maxfpo
            temp = (funz*tempa-funa*tempz)/(funz-funa)
            !----- Checking whether this guess will fall outside the range ----------------!
            if (abs(temp-tempa) > abs(tempz-tempa) .or.                                    &
                abs(temp-tempz) > abs(tempz-tempa)) then
               temp = 0.5*(tempa+tempz)
            end if
            !----- Distributing vapour into the three phases ------------------------------!
            rsat   = max(toodry,rslif(pres,temp))
            rvap   = min(rtot,rsat)
            if (temp >= t3ple) then
               rliq = max(0.,rtot-rsat)
               rice = 0.
            else
               rliq = 0.
               rice = max(0.,rtot-rsat)
            end if
            !----- Updating function ------------------------------------------------------!
            funnow = theta_iceliq(exner,temp,rliq,rice) - thil

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (debug) then
               write (unit=46,fmt='(a,1x,i5,1x,10(a,1x,f11.4,1x))')                        &
                  'REGFAL: it=',itb,'temp=',temp-t00,'tempa=',tempa-t00,'tempz=',tempz-t00 &
                 ,'rsat=',1000.*rsat,'rliq=',1000.*rliq,'rice=',1000.*rice                 &
                 ,'rvap=',1000.*rvap,'fun=',funnow,'funa=',funa,'funz=',funz
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Checking for convergence or lucky guess. If it did, return, we found the  !
            ! solution. Otherwise, constrain the guesses.                                  !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa) < toler*temp .and. abs(temp-tempz) < toler*temp 
            if (funnow == 0. .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.) then 
               tempz = temp
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 0.5
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            elseif (funnow*funz < 0.) then
               tempa = temp
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 0.5
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false.
            end if
            !------------------------------------------------------------------------------!
         end do fpoloop
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Almost done... Usually when the method goes through regula falsi, it means   !
         ! that the temperature is too close to the triple point, and often all three      !
         ! phases will coexist.  The problem with the method is that it converges for      !
         ! temperature, but whenever regula falsi is called the function evaluation is     !
         ! usually far from zero.  This can be improved by finding a better partition      !
         ! between ice and liquid given the temperature and saturation mixing ratio we     !
         ! just found. So just to round these edges, we will invert the ice-liquid         !
         ! potential temperature using the set of temperature and rsat, and fiding the     !
         ! liquid mixing ratio.                                                            !
         !---------------------------------------------------------------------------------!
         if (abs(temp-t3ple) < toler*temp) then
            !----- Find rliq. -------------------------------------------------------------!
            rliq   = ( alvi(temp) * (rtot - rsat)                                          &
                     + cpdry * temp * log ( cpdryi * exner * thil / temp ) )               &
                   / ( alvi(temp) - alvl(temp) )
            !------------------------------------------------------------------------------!

            !----- Correct rliq so it is bounded, and then find rice as the remainder. ----!
            rliq   = max( 0., min( rtot - rsat,  rliq ) )
            rice   = max( 0., rtot-rsat-rliq)
            !------------------------------------------------------------------------------!


            !----- Function evaluation. ---------------------------------------------------!
            funnow = theta_iceliq(exner,temp,rliq,rice) - thil
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
         itb=itb+1
      end if
      !------------------------------------------------------------------------------------!

      if (.not. converged) then
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
         write (unit=*,fmt='(a)')           ' THIL2TQALL failed!'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Input: '
         write (unit=*,fmt='(a,1x,f12.5)')  '    THETA_IL [     K]:',thil
         write (unit=*,fmt='(a,1x,f12.5)')  '    EXNER    [J/kg/K]:',exner
         write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1000.*rtot
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Output: '
         write (unit=*,fmt='(a,1x,i12)')    '    ITERATIONS       :',itb
         write (unit=*,fmt='(a,1x,f12.5)')  '    TEMP     [    °C]:',temp-t00
         write (unit=*,fmt='(a,1x,f12.5)')  '    RVAP     [  g/kg]:',1000.*rvap
         write (unit=*,fmt='(a,1x,f12.5)')  '    RLIQ     [  g/kg]:',1000.*rliq
         write (unit=*,fmt='(a,1x,f12.5)')  '    RICE     [  g/kg]:',1000.*rice
         write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPA    [    °C]:',tempa-t00
         write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPZ    [    °C]:',tempz-t00
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNA     [     K]:',funnow
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNZ     [     K]:',funnow
         write (unit=*,fmt='(a,1x,es12.5)') '    DERIV    [   ---]:',deriv
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_A    [   ---]:',abs(temp-tempa)/temp
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_Z    [   ---]:',abs(temp-tempz)/temp
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
         call abort_run  ('Failed finding equilibrium, I gave up!','thil2tqall'            &
                         ,'therm_lib.f90')
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x))')                               &
            'ANSWER: it=',itb,'temp=',temp-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq       &
                             ,'rice=',1000.*rice,'rvap=',1000.*rvap,'funf=',funnow
         write (unit=46,fmt='(a)') '------------------------------------------------------'
         write (unit=46,fmt='(a)') ''
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      return
   end subroutine thil2tqall
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes a consistent set of temperature and condensated phases   !
   ! mixing ratio for a given theta_il, Exner function, and total mixing ratio.  This is   !
   ! very similar to the function thil2temp, except that now we don't know rliq.  Rliq     !
   ! becomes a function of temperature, since it is defined as rtot-rsat(T,p), remembering !
   ! that rtot and p are known. If the air is not saturated, we use the fact that theta_il !
   ! is theta and skip the hassle.  Otherwise, we use iterative methods.  We will always   !
   ! try Newton's method, since it converges fast.  Not always will Newton converge, and   !
   ! if that's the case we use a modified regula falsi (Illinois) method.  This method is  !
   ! a mix of secant and bisection and will always converge.                               !
   !---------------------------------------------------------------------------------------!
   subroutine thil2tqliq(thil,exner,pres,rtot,rliq,temp,rvap,rsat)
      use rconstants , only : cpdryi   & ! intent(in)
                            , toodry   ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in)    :: thil      ! Ice-liquid water potential temp.  [     K]
      real(kind=4), intent(in)    :: exner     ! Exner function                    [J/kg/K]
      real(kind=4), intent(in)    :: pres      ! Pressure                          [    Pa]
      real(kind=4), intent(in)    :: rtot      ! Total mixing ratio                [ kg/kg]
      real(kind=4), intent(out)   :: rliq      ! Liquid water mixing ratio         [ kg/kg]
      real(kind=4), intent(inout) :: temp      ! Temperature                       [     K]
      real(kind=4), intent(out)   :: rvap      ! Water vapour mixing ratio         [ kg/kg]
      real(kind=4), intent(out)   :: rsat      ! Sat. water vapour mixing ratio    [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                :: tempa     ! Lower bound for regula falsi iteration
      real(kind=4)                :: tempz     ! Upper bound for regula falsi iteration
      real(kind=4)                :: t1stguess ! Book keeping temperature 1st guess 
      real(kind=4)                :: fun1st    ! Book keeping 1st guess function
      real(kind=4)                :: funa      ! Function evaluation at tempa
      real(kind=4)                :: funz      ! Function evaluation at tempz
      real(kind=4)                :: funnow    ! Function at this iteration.
      real(kind=4)                :: delta     ! Aux. var in case we need regula falsi.
      real(kind=4)                :: deriv     ! Derivative of this function.
      integer                     :: itn       ! Iteration counter
      integer                     :: itb       ! Iteration counter
      logical                     :: converged ! Convergence handle
      logical                     :: zside     ! Aux. Flag, for two purposes:
                                               ! 1. Found a 2nd guess for regula falsi.
                                               ! 2. I retained the "zside" (T/F)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      First check: try to find temperature assuming sub-saturation and check if     !
      ! this is the case. If it is, then there is no need to go through the iterative      !
      ! loop.                                                                              !
      !------------------------------------------------------------------------------------!
      tempz = cpdryi * thil * exner
      rsat  = max(toodry,rslf(pres,tempz))
      rliq  = max(0.,rtot-rsat)
      rvap  = rtot-rliq
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    If rtot < rsat, this is not saturated, we can leave the subroutine and bypass   !
      ! the iterative part.                                                                !
      !------------------------------------------------------------------------------------!
      if (rtot < rsat) then
         temp = tempz
         return
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   If not, then use the temperature the user gave as first guess and solve          !
      ! iteratively.  We use the user instead of what we just found because if the air is  !
      ! saturated, then this can be too far off which may be bad for Newton's method.      !
      !------------------------------------------------------------------------------------!
      tempz = temp
      rsat   = max(toodry,rslf(pres,tempz))
      rliq = max(0.,rtot-rsat)
      rvap = rtot-rliq
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Finding the function. We are seeking a temperature which is associated with    !
      ! the theta_il we provided. Thus, the function is simply the difference between the  !
      ! theta_il associated with our guess and the actual theta_il.                        !
      !------------------------------------------------------------------------------------!
      funnow = theta_iceliq(exner,tempz,rliq,0.) ! Finding thil from our guess
      deriv  = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq)
      funnow = funnow - thil ! Computing the function
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Now we enter at the Newton's method iterative loop. We are always going to try  !
      ! this first, because it's fast, but if it turns out to be a dangerous choice or if  !
      ! it doesn't converge fast, we will fall back to regula falsi.                       !
      !    We start by initialising the flag and copying temp to tempz, the newest guess.  !
      !------------------------------------------------------------------------------------!
      converged=.false.
      newloop: do itn=1,maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Saving previous guess. We also save the function is in case we give up on   !
         ! Newton's and switch to regula falsi.                                            !
         !---------------------------------------------------------------------------------!
         funa  = funnow
         tempa = tempz

         !----- Go to bisection if the derivative is too flat (too dangerous...) ----------!
         if (abs(deriv) < toler) exit newloop

         tempz = tempa - funnow / deriv

         !----- Finding the mixing ratios associated with this guess ----------------------!
         rsat  = max(toodry,rslf(pres,tempz))
         rliq = max(0.,rtot-rsat)
         rvap = rtot-rliq

         !----- Updating the function -----------------------------------------------------!
         funnow = theta_iceliq(exner,tempz,rliq,0.)
         !----- Updating the derivative. --------------------------------------------------!
         deriv = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq)
         funnow = funnow - thil

         converged = abs(tempa-tempz) < toler*tempz
         !---------------------------------------------------------------------------------!
         !   Convergence. The temperature will be the mid-point between tempa and tempz.   !
         ! Fix the mixing ratios and return. But first check for converged due to luck. If !
         ! the guess gives a root, then that's it. It looks unlikely, but it actually      !
         ! happens sometimes and if not checked it becomes a singularity.                  !
         !---------------------------------------------------------------------------------!
         if (funnow == 0.) then
            temp = tempz
            converged = .true.
            exit newloop
         elseif (converged) then
            temp = 0.5 * (tempa+tempz)
            rsat  = max(toodry,rslf(pres,temp))
            rliq = max(0.,rtot-rsat)
            rvap = rtot-rliq
            exit newloop
         end if
         !---------------------------------------------------------------------------------!
      end do newloop
      !---------------------------------------------------------------------------------------!



      !---------------------------------------------------------------------------------------!
      if (.not. converged) then
         !---------------------------------------------------------------------------------!
         !    If I reach this point, then it means that Newton's method failed finding the !
         ! equilibrium, so we are going to use the regula falsi instead.  If Newton's      !
         ! method didn't converge, we use tempa as one guess and now we seek a tempz with  !
         ! opposite sign.                                                                  !
         !---------------------------------------------------------------------------------!
         !----- Check funa and funnow have opposite signs. If so, we are ready to go ------!
         if (funa*funnow < 0.0) then
            funz = funnow
            zside = .true.
         !---------------------------------------------------------------------------------!
         !    Looking for a guess. Extrapolate funa linearly, trying to get the -funa.  We !
         ! don't need it to be funa, just with the opposite sign.  If that's not enough,   !
         ! we keep going further... Force the guesses to be at least 1K apart              !
         !---------------------------------------------------------------------------------!
         else
            if (abs(funnow-funa) < toler*tempa) then
               delta = 100.*toler*tempa
            else
               delta = max(abs(funa*(tempz-tempa)/(funnow-funa)),100.*toler*tempa)
            end if
            tempz = tempa + delta
            funz  = funa
            !----- Just to enter at least once. The 1st time tempz=tempa-2*delta ----------!
            zside = .false. 
            zgssloop: do itb=1,maxfpo
               tempz = tempz + real((-1)**itb * (itb+3)/2) * delta
               rsat  = max(toodry,rslf(pres,tempz))
               rliq  = max(0.,rtot-rsat)
               rvap  = rtot-rliq
               funz  = theta_iceliq(exner,tempz,rliq,0.) - thil
               zside = funa*funz < 0.0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside)                                                               &
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thil2tqliq','rthrm.f90')
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Now we loop until convergence is achieved.                                  !
         !---------------------------------------------------------------------------------!
         fpoloop: do itb=itn,maxfpo
            temp = (funz*tempa-funa*tempz)/(funz-funa)
            !----- Distributing vapour into the three phases ------------------------------!
            rsat   = max(toodry,rslf(pres,temp))
            rvap   = min(rtot,rsat)
            rliq   = max(0.,rtot-rsat)
            !----- Updating function ------------------------------------------------------!
            funnow = theta_iceliq(exner,tempz,rliq,0.) - thil

            !------------------------------------------------------------------------------!
            !    Checking for convergence or lucky guess. If it did, return, we found the  !
            ! solution. Otherwise, constrain the guesses.                                  !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa)< toler*temp  .and. abs(temp-tempz) < toler*temp
            if (funnow == 0. .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.) then 
               tempz = temp
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 0.5
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tempa = temp
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 0.5
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if

         end do fpoloop

      end if

      if (.not. converged) call abort_run  ('Failed finding equilibrium, I gave up!'       &
                                           ,'thil2tqliq','therm_lib.f90')
      return
   end subroutine thil2tqliq
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! intensive internal energy [J/kg].                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine uint2tl(uint,temp,fliq)
      use rconstants , only : cliqi          & ! intent(in)
                            , cicei          & ! intent(in)
                            , allii          & ! intent(in)
                            , t3ple          & ! intent(in)
                            , uiicet3        & ! intent(in)
                            , uiliqt3        & ! intent(in)
                            , tsupercool_liq ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)  :: uint     ! Internal energy                     [   J/kg]
      real(kind=4), intent(out) :: temp     ! Temperature                         [      K]
      real(kind=4), intent(out) :: fliq  ! Liquid Fraction (0-1)               [    ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compare the internal energy with the reference values to decide which phase    !
      ! the water is.                                                                      !
      !------------------------------------------------------------------------------------!
      if (uint <= uiicet3) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.
         temp = uint * cicei
         !---------------------------------------------------------------------------------!
      elseif (uint >= uiliqt3) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.
         temp = uint * cliqi + tsupercool_liq
         !---------------------------------------------------------------------------------!
      else
         !----- Changing phase, it must be at freezing point ------------------------------!
         fliq = (uint - uiicet3) * allii
         temp = t3ple
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uint2tl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from         !
   ! extensive internal energy (J/m² or J/m³), water mass (kg/m² or kg/m³), and heat       !
   ! capacity (J/m²/K or J/m³/K).                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine uextcm2tl(uext,wmass,dryhcap,temp,fliq)
      use rconstants , only : cliqi          & ! intent(in)
                            , cliq           & ! intent(in)
                            , cicei          & ! intent(in)
                            , cice           & ! intent(in)
                            , allii          & ! intent(in)
                            , alli           & ! intent(in)
                            , t3ple          & ! intent(in)
                            , tsupercool_liq ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)  :: uext    ! Extensive internal energy [  J/m²] or [  J/m³]
      real(kind=4), intent(in)  :: wmass   ! Water mass                [ kg/m²] or [ kg/m³]
      real(kind=4), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m²/K] or [J/m³/K]
      real(kind=4), intent(out) :: temp    ! Temperature                           [     K]
      real(kind=4), intent(out) :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=4)              :: uefroz  ! qw of ice at triple pt.   [  J/m²] or [  J/m³] 
      real(kind=4)              :: uemelt  ! qw of liq. at triple pt.  [  J/m²] or [  J/m³]
      !------------------------------------------------------------------------------------!



      !----- Convert melting heat to J/m² or J/m³ -----------------------------------------!
      uefroz = (dryhcap + wmass * cice) * t3ple
      uemelt = uefroz   + wmass * alli
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    This is analogous to the uint2tl computation, we should analyse the magnitude   !
      ! of the internal energy to choose between liquid, ice, or both by comparing with    !
      ! the known boundaries.                                                              !
      !------------------------------------------------------------------------------------!
      if (uext < uefroz) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.
         temp = uext  / (cice * wmass + dryhcap)
         !---------------------------------------------------------------------------------!
      elseif (uext > uemelt) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.
         temp = (uext + wmass * cliq * tsupercool_liq) / (dryhcap + wmass * cliq)
         !---------------------------------------------------------------------------------!
      elseif (uefroz == uemelt) then
         !---------------------------------------------------------------------------------!
         !    We are at the freezing point.  If water mass is so tiny that the internal    !
         ! energy of frozen and melted states are the same given the machine precision,    !
         ! then we assume that water content is negligible and we impose 50% frozen for    !
         ! simplicity.                                                                     !
         !---------------------------------------------------------------------------------!
         fliq = 0.5
         temp = t3ple
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    Changing phase, it must be at freezing point.  The max and min are here just !
         ! to avoid tiny deviations beyond 0. and 1. due to floating point arithmetics.    !
         !---------------------------------------------------------------------------------!
         fliq = min(1.,max(0.,(uext - uefroz) * allii / wmass))
         temp = t3ple
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uextcm2tl
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib
!==========================================================================================!
!==========================================================================================!
