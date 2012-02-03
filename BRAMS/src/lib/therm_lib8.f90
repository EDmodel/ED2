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
                       , newthermo4 => newthermo ! ! intent(in)
  
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
   integer ::   level

   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical ::   vapour_on
   logical ::   cloud_on
   logical ::   bulk_on
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
   real(kind=8), dimension(0:8), parameter :: cll8 = (/  .6105851d+03,  .4440316d+02       &
                                                      ,  .1430341d+01,  .2641412d-01       &
                                                      ,  .2995057d-03,  .2031998d-05       &
                                                      ,  .6936113d-08,  .2564861d-11       &
                                                      , -.3704404d-13                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: cii8 = (/  .6114327d+03,  .5027041d+02       &
                                                      ,  .1875982d+01,  .4158303d-01       &
                                                      ,  .5992408d-03,  .5743775d-05       &
                                                      ,  .3566847d-07,  .1306802d-09       &
                                                      ,  .2152144d-12                /)
   !----- Coefficients for d(esat)/dT (liquid) --------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: dll8 = (/  .4443216d+02,  .2861503d+01       &
                                                      ,  .7943347d-01,  .1209650d-02       &
                                                      ,  .1036937d-04,  .4058663d-07       &
                                                      , -.5805342d-10, -.1159088d-11       &
                                                      , -.3189651d-14                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: dii8 = (/  .5036342d+02,  .3775758d+01       &
                                                      ,  .1269736d+00,  .2503052d-02       &
                                                      ,  .3163761d-04,  .2623881d-06       &
                                                      ,  .1392546d-08,  .4315126d-11       &
                                                      ,  .5961476d-14                /)
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
      use rconstants , only : t008 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)            :: temp     ! Temperature                [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=8), intent(out), optional :: l1funout ! Function for high temperatures
      real(kind=8), intent(out), optional :: ttfunout ! Interpolation function
      real(kind=8), intent(out), optional :: l2funout ! Function for low temperatures
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                        :: l1fun    ! 
      real(kind=8)                        :: ttfun    ! 
      real(kind=8)                        :: l2fun    ! 
      real(kind=8)                        :: x        ! 
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         l1fun = l01_108(0) + l01_108(1)/temp + l01_108(2)*log(temp) + l01_108(3) * temp
         l2fun = l02_108(0) + l02_108(1)/temp + l02_108(2)*log(temp) + l02_108(3) * temp
         ttfun = tanh(ttt_108(1) * (temp - ttt_108(2)))
         eslf8 = exp(l1fun + ttfun*l2fun)
         !---------------------------------------------------------------------------------!

         if (present(l1funout)) l1funout = l1fun
         if (present(l2funout)) l2funout = l2fun
         if (present(ttfunout)) ttfunout = ttfun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-8.0d1,temp-t008)
         eslf8 = cll8(0) + x * (cll8(1) + x * (cll8(2) + x * (cll8(3) + x * (cll8(4)       &
                         + x * (cll8(5) + x * (cll8(6) + x * (cll8(7) + x * cll8(8)) ))))))
         !---------------------------------------------------------------------------------!

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
      use rconstants , only : t008 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)            :: temp     ! Temperature                 [    K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=8), intent(out), optional :: iifunout
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                        :: iifun
      real(kind=8)                        :: x
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then

         !----- Updated method, using MK05 ------------------------------------------------!
         iifun = iii_78(0) + iii_78(1)/temp + iii_78(2) * log(temp) + iii_78(3) * temp
         esif8 = exp(iifun)
         !---------------------------------------------------------------------------------!


         if (present(iifunout)) iifunout=iifun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-8.d1,temp-t008)
         esif8 = cii8(0) + x * (cii8(1) + x * (cii8(2) + x * (cii8(3) + x * (cii8(4)       &
                        + x * (cii8(5) + x * (cii8(6) + x * (cii8(7) + x * cii8(8))))))))
         !---------------------------------------------------------------------------------!

         if (present(iifunout)) iifunout=esif8
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : t3ple8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      logical                            :: frozen
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         eslif8 = esif8(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         eslif8 = eslf8(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=8), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: esl  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esl  = eslf8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the mixing !
      ! ratio.                                                                             !
      !------------------------------------------------------------------------------------!
      rslf8 = max(toodry8,ep8*esl/(pres-esl))
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=8), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: esi  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esi  = esif8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the mixing !
      ! ratio.                                                                             !
      !------------------------------------------------------------------------------------!
      rsif8 = max(toodry8,ep8*esi/(pres-esi))
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : t3ple8 & ! intent(in)
                            , ep8    ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres
      real(kind=8), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: esz
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         esz = esif8(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         esz = eslf8(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the mixing !
      ! ratio.                                                                             !
      !------------------------------------------------------------------------------------!
      rslif8 = ep8 * esz / (pres - esz)
      !------------------------------------------------------------------------------------!

      return
   end function rslif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation specific humidity as a function of !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function qslf8(pres,temp)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=8), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: esl  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esl  = eslf8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qslf8 = max(toodry8,ep8 * esl/( pres - (1.d0 - ep8) * esl) )
      !------------------------------------------------------------------------------------!

      return
   end function qslf8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation specific humidity as a function of    !
   ! pressure and Kelvin temperature.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function qsif8(pres,temp)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres ! Pressure                                   [   Pa]
      real(kind=8), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: esi  ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !----- First we find the saturation vapour pressure. --------------------------------!
      esi  = esif8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qsif8 = max(toodry8,ep8 * esi/( pres - (1.d0 - ep8) * esi) )
      !------------------------------------------------------------------------------------!

      return
   end function qsif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation specific humidity, over liquid or ice     !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function qslif8(pres,temp,useice)
      use rconstants , only : t3ple8  & ! intent(in)
                            , ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres
      real(kind=8), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: esz
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         esz = esif8(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         esz = eslf8(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qslif8 = max(toodry8, ep8 * esz/( pres - (1.d0 - ep8) * esz) )
      !------------------------------------------------------------------------------------!

      return
   end function qslif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-liquid equilibrium density for vapour, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function rhovsl8(temp)
      use rconstants , only : rh2o8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: eequ ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the equilibrium (saturation) vapour pressure.                             !
      !------------------------------------------------------------------------------------!
      eequ = eslf8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsl8 = eequ / (rh2o8 * temp)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : rh2o8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: eequ ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the equilibrium (saturation) vapour pressure.                             !
      !------------------------------------------------------------------------------------!
      eequ = esif8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsi8 = eequ / (rh2o8 * temp)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : rh2o8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: eequ
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Pass the "useice" argument to eslif, so it may decide whether ice thermo-      !
      ! dynamics is to be used.                                                            !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         eequ = eslif8(temp,useice)
      else
         eequ = eslif8(temp)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsil8 = eequ / (rh2o8 * temp)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : t008 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp
      !------ Local variables. ------------------------------------------------------------!
      real(kind=8)             :: esl
      real(kind=8)             :: l2fun
      real(kind=8)             :: ttfun
      real(kind=8)             :: l1prime
      real(kind=8)             :: l2prime
      real(kind=8)             :: ttprime
      real(kind=8)             :: x
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide which function to use, based on the thermodynamics.                    !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esl     = eslf8(temp,l2funout=l2fun,ttfunout=ttfun)
         l1prime = -l01_108(1)/(temp*temp) + l01_108(2)/temp + l01_108(3)
         l2prime = -l02_108(1)/(temp*temp) + l02_108(2)/temp + l02_108(3)
         ttprime =  ttt_108(1)*(1.d0 - ttfun*ttfun)
         eslfp8  = esl * (l1prime + l2prime*ttfun + l2fun*ttprime)
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x      = max(-8.d1,temp-t008)
         eslfp8 = dll8(0) + x * (dll8(1) + x * (dll8(2) + x * (dll8(3) + x * (dll8(4)      &
                          + x * (dll8(5) + x * (dll8(6) + x * (dll8(7) + x * dll8(8))))))))
      end if
      !------------------------------------------------------------------------------------!


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
      use rconstants , only : t008 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp
      !------ Local variables. ------------------------------------------------------------!
      real(kind=8)             :: esi
      real(kind=8)             :: iiprime
      real(kind=8)             :: x
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide which function to use, based on the thermodynamics.                    !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esi     = esif8(temp)
         iiprime = -iii_78(1)/(temp*temp) + iii_78(2)/temp + iii_78(3)
         esifp8  = esi * iiprime
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x      = max(-8.d1,temp-t008)
         esifp8 = dii8(0) + x * (dii8(1) + x * (dii8(2) + x * (dii8(3) + x * (dii8(4)      &
                          + x * (dii8(5) + x * (dii8(6) + x * (dii8(7) + x * dii8(8))))))))
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : t3ple8 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=8), intent(in)           :: temp
      !------ Local variables. ------------------------------------------------------------!
      logical     , intent(in), optional :: useice
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- d(Saturation vapour pressure)/dT for ice. ---------------------------------!
         eslifp8 = esifp8(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- d(Saturation vapour pressure)/dT for liquid water. ------------------------!
         eslifp8 = eslfp8(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : ep8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres  ! Pressure                                 [    Pa]
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: esl   ! Partial pressure                         [    Pa]
      real(kind=8)             :: desdt ! Derivative of partial pressure of water  [  Pa/K]
      real(kind=8)             :: pdry  ! Partial pressure of dry air              [    Pa]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      esl    = eslf8(temp)
      desdt  = eslfp8(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial pressure of dry air. ----------------------------------------!
      pdry   = pres-esl
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of mixing ratio. ---------------------------------!
      rslfp8  = ep8 * pres * desdt / (pdry*pdry)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : ep8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres  ! Pressure                                 [    Pa]
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: esi   ! Partial pressure                         [    Pa]
      real(kind=8)             :: desdt ! Derivative of partial pressure of water  [  Pa/K]
      real(kind=8)             :: pdry  ! Partial pressure of dry air              [    Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      esi    = esif8(temp)
      desdt  = esifp8(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial pressure of dry air. ----------------------------------------!
      pdry   = pres-esi
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of mixing ratio. ---------------------------------!
      rsifp8  = ep8 * pres * desdt / (pdry*pdry)
      !------------------------------------------------------------------------------------!
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
      use rconstants , only: t3ple8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres   ! Pressure                      [    Pa]
      real(kind=8), intent(in)           :: temp   ! Temperature                   [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! May use ice thermodynamics?   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: desdt  ! Derivative of vapour pressure [  Pa/K]
      logical                            :: frozen ! Use the ice thermodynamics    [   T|F]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide whether to use liquid water of ice for saturation, based on the temper- !
      ! ature and the settings.                                                            !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the function depending on the previous check.                             !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         rslifp8 = rsifp8(pres,temp)
      else
         rslifp8 = rslfp8(pres,temp)
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : rh2o8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in) :: temp   ! Temperature                             [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: es     ! Vapour pressure                         [    Pa]
      real(kind=8)             :: desdt  ! Vapour pressure derivative              [  Pa/K]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      es    = eslf8(temp)
      desdt = eslfp8(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of saturation density . --------------------------!
      rhovslp8 = (desdt-es/temp) / (rh2o8 * temp)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : rh2o8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in) :: temp   ! Temperature                             [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: es     ! Vapour pressure                         [    Pa]
      real(kind=8)             :: desdt  ! Vapour pressure derivative              [  Pa/K]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the partial pressure of water vapour and its derivative, using temper-    !
      ! ature.                                                                             !
      !------------------------------------------------------------------------------------!
      es    = esif8(temp)
      desdt = esifp8(temp)
      !------------------------------------------------------------------------------------!


      !----- Find the partial derivative of saturation density . --------------------------!
      rhovsip8 = (desdt - es/temp) / (rh2o8 * temp)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : t3ple8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: temp   ! Temperature                   [     K]
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
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the function depending on the previous check.                             !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         rhovsilp8 = rhovsip8(temp)
      else
         rhovsilp8 = rhovslp8(temp)
      end if
      !------------------------------------------------------------------------------------!

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
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pvap      ! Saturation vapour pressure           [    Pa]
      !----- Local variables for iterative method. ----------------------------------------!
      real(kind=8)             :: deriv     ! Function derivative                  [    Pa]
      real(kind=8)             :: fun       ! Function for which we seek a root.   [    Pa]
      real(kind=8)             :: funa      ! Smallest  guess function             [    Pa]
      real(kind=8)             :: funz      ! Largest   guess function             [    Pa]
      real(kind=8)             :: tempa     ! Smallest guess (or previous guess)   [    Pa]
      real(kind=8)             :: tempz     ! Largest guess (new guess in Newton)  [    Pa]
      real(kind=8)             :: delta     ! Aux. var --- 2nd guess for bisection [      ]
      integer                  :: itn       ! Iteration counter                    [   ---]
      integer                  :: itb       ! Iteration counter                    [   ---]
      logical                  :: converged ! Convergence handle                   [   ---]
      logical                  :: zside     ! Flag to check for one-sided approach [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, use Bolton (1980) equation 11, giving es in Pa and T in K -------!
      tempa = (2.965d1 * log(pvap) - 5.01678d3)/(log(pvap)-2.40854d1)
      funa  = eslf8(tempa) - pvap
      deriv = eslfp8(tempa)
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration ----------------------------!
      tempz = tempa
      fun   = funa
      !------------------------------------------------------------------------------------!


      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, go with bisection ----!

         !----- Copy the previous guess. --------------------------------------------------!
         tempa = tempz
         funa  = fun
         !---------------------------------------------------------------------------------!


         !----- New guess, its function, and derivative evaluation. -----------------------!
         tempz = tempa - fun/deriv
         fun   = eslf8(tempz) - pvap
         deriv = eslfp8(tempz)
         !---------------------------------------------------------------------------------!


         !----- Check convergence. --------------------------------------------------------!
         converged = abs(tempa-tempz) < toler8 * tempz
         if (converged) then
            tslf8 = 5.d-1 * (tempa+tempz)
            return
         elseif (fun == 0.0d0) then 
            !----- Converged by luck. -----------------------------------------------------!
            tslf8 = tempz
            return
         end if
         !---------------------------------------------------------------------------------!
      end do newloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If we have reached this point, then Newton's method has failed.  Use bisection !
      ! instead.  For bisection, we need two guesses whose function has opposite signs.    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         !----- We already have two guesses with opposite signs. --------------------------!
         funz  = fun
         zside = .true.
         !---------------------------------------------------------------------------------!
      else
         !----- Need to find the guesses with opposite signs. -----------------------------!
         if (abs(fun-funa) < 1.d2*toler8*tempa) then
            delta = 1.d2*toler8*tempa
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),1.d2*toler8*tempa)
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Try guesses on both sides of the first guess, increasingly further away      !
         ! until we spot a good guess.                                                     !
         !---------------------------------------------------------------------------------!
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
            funz  = eslf8(tempz) - pvap
            zside = funa*funz < 0.d0
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
                            ,'tslf8','therm_lib8.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tslf8 =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close.  If so, !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tslf8-tempa) < toler8 * tslf8
         if (converged) exit bisloop

         !------ Find the new function evaluation. ----------------------------------------!
         fun       =  eslf8(tslf8) - pvap
         !---------------------------------------------------------------------------------!


         !------ Define the new interval based on the intermediate value theorem. ---------!
         if (fun*funa < 0.d0 ) then
            tempz = tslf8
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method). --------!
            if (zside) funa = funa * 5.d-1
            !----- We have just updated zside, so we set zside to true. -------------------!
            zside = .true.
         else
            tempa = tslf8
            funa   = fun
            !----- If we are updating aside again, modify zside (Illinois method). --------!
            if (.not. zside) funz = funz * 5.d-1
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
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pvap      ! Saturation vapour pressure           [    Pa]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)             :: deriv     ! Function derivative                  [    Pa]
      real(kind=8)             :: fun       ! Function for which we seek a root.   [    Pa]
      real(kind=8)             :: funa      ! Smallest  guess function             [    Pa]
      real(kind=8)             :: funz      ! Largest   guess function             [    Pa]
      real(kind=8)             :: tempa     ! Smallest guess (or previous guess)   [    Pa]
      real(kind=8)             :: tempz     ! Largest guess (new guess in Newton)  [    Pa]
      real(kind=8)             :: delta     ! Aux. var --- 2nd guess for bisection [      ]
      integer                  :: itn
      integer                  :: itb       ! Iteration counter                    [   ---]
      logical                  :: converged ! Convergence handle                   [   ---]
      logical                  :: zside     ! Flag to check for one-sided approach [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, use Murphy-Koop (2005), equation 8. -----------------------------!
      tempa = (1.814625d0 * log(pvap) +6.190134d3)/(2.9120d1 - log(pvap))
      funa  = esif8(tempa) - pvap
      deriv = esifp8(tempa)
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration ----------------------------!
      tempz = tempa
      fun   = funa
      !------------------------------------------------------------------------------------!


      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, go with bisection ----!

         !----- Copy the previous guess. --------------------------------------------------!
         tempa = tempz
         funa  = fun


         !----- New guess, its function, and derivative evaluation. -----------------------!
         tempz = tempa - fun/deriv
         fun   = esif8(tempz) - pvap
         deriv = esifp8(tempz)
         !---------------------------------------------------------------------------------!


         !----- Check convergence. --------------------------------------------------------!
         converged = abs(tempa-tempz) < toler8 * tempz
         if (converged) then
            tsif8 = 5.d-1 * (tempa+tempz)
            return
         elseif (fun == 0.d0) then
            tsif8 = tempz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If we have reached this point, then Newton's method has failed.  Use bisection !
      ! instead.  For bisection, we need two guesses whose function has opposite signs.    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         !----- We already have two guesses with opposite signs. --------------------------!
         funz  = fun
         zside = .true.
         !---------------------------------------------------------------------------------!
      else
         !----- Need to find the guesses with opposite signs. -----------------------------!
         if (abs(fun-funa) < 1.d2*toler8*tempa) then
            delta = 1.d2*toler8*delta
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),1.d2*toler8*tempa)
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Try guesses on both sides of the first guess, increasingly further away      !
         ! until we spot a good guess.                                                     !
         !---------------------------------------------------------------------------------!
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
            funz  = esif8(tempz) - pvap
            zside = funa*funz < 0.d0
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
                            ,'tsif8','therm_lib8.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tsif8 =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close.  If so, !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tsif8-tempa) < toler8 * tsif8
         if (converged) exit bisloop
         !---------------------------------------------------------------------------------!

         !------ Find the new function evaluation. ----------------------------------------!
         fun       =  esif8(tsif8) - pvap
         !---------------------------------------------------------------------------------!


         !------ Define the new interval based on the intermediate value theorem. ---------!
         if (fun*funa < 0.d0 ) then
            tempz = tsif8
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method). --------!
            if (zside) funa = funa * 5.d-1
            !----- We have just updated zside, so we set zside to true. -------------------!
            zside = .true.
         else
            tempa = tsif8
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method). --------!
            if (.not. zside) funz = funz * 5.d-1
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
      use rconstants , only : es3ple8 ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pvap   ! Vapour pressure                [   Pa]
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
         frozen = useice  .and. pvap < es3ple8
      else 
         frozen = bulk_on .and. pvap < es3ple8
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the function depending on whether we should use ice.                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         tslif8 = tsif8(pvap)
      else
         tslif8 = tslf8(pvap)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function tslif8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the dew point temperature given the pressure and vapour    !
   ! mixing ratio. THIS IS DEW POINT ONLY, WHICH MEANS THAT IT WILL IGNORE ICE EFFECT. For !
   ! a full, triple-point dependent routine use DEWFROSTPOINT.                             !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dewpoint8(pres,rsat)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=8), intent(in) :: rsat    ! Saturation mixing ratio                [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: rsatoff ! Non-singular saturation mixing ratio   [ kg/kg]
      real(kind=8)             :: pvsat   ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !----- Make sure mixing ratio is positive. ------------------------------------------!
      rsatoff   = max(toodry8,rsat)
      !------------------------------------------------------------------------------------!


      !----- Find the saturation vapour pressure. -----------------------------------------!
      pvsat     = pres * rsatoff / (ep8 + rsatoff)
      !------------------------------------------------------------------------------------!

      !----- Dew point is going to be the saturation temperature. -------------------------!
      dewpoint8 = tslf8(pvsat)
      !------------------------------------------------------------------------------------!

      return
   end function dewpoint8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the frost point temperature given the pressure and vapour  !
   ! mixing ratio. THIS IS FROST POINT ONLY, WHICH MEANS THAT IT WILL IGNORE LIQUID        !
   ! EFFECT.  For a full, triple-point dependent routine use DEWFROSTPOINT.                !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function frostpoint8(pres,rsat)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=8), intent(in) :: rsat    ! Saturation mixing ratio                [ kg/kg]
      !----- Local variables for iterative method. ----------------------------------------!
      real(kind=8)             :: rsatoff ! Non-singular saturation mixing ratio   [ kg/kg]
      real(kind=8)             :: pvsat   ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !----- Make sure mixing ratio is positive. ------------------------------------------!
      rsatoff     = max(toodry8,rsat)
      !------------------------------------------------------------------------------------!


      !----- Find the saturation vapour pressure. -----------------------------------------!
      pvsat       = pres*rsatoff / (ep8 + rsatoff)
      !------------------------------------------------------------------------------------!

      !----- Frost point is going to be the saturation temperature. -----------------------!
      frostpoint8 = tsif8(pvsat)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres    ! Pressure                     [    Pa]
      real(kind=8), intent(in)           :: rsat    ! Saturation mixing ratio      [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: rsatoff ! Non-singular sat. mix. rat.  [ kg/kg]
      real(kind=8)                       :: pvsat   ! Saturation vapour pressure   [    Pa]
      !------------------------------------------------------------------------------------!


      !----- Make sure mixing ratio is positive. ------------------------------------------!
      rsatoff  = max(toodry8,rsat)
      !------------------------------------------------------------------------------------!

      !----- Find the saturation vapour pressure. -----------------------------------------!
      pvsat         = pres*rsatoff / (ep8 + rsatoff)
      !------------------------------------------------------------------------------------!

      !----- Dew (frost) point is going to be the saturation temperature. -----------------!
      if (present(useice)) then
         dewfrostpoint8 = tslif8(pvsat,useice)
      else
         dewfrostpoint8 = tslif8(pvsat)
      end if
      !------------------------------------------------------------------------------------!
      return
   end function dewfrostpoint8
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
   real(kind=8) function ptrh2rvapl8(relh,pres,temp,out_shv)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: relh    ! Relative humidity                      [    --]
      real(kind=8), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=8), intent(in) :: temp    ! Temperature                            [     K]
      logical     , intent(in) :: out_shv ! Output is specific humidity            [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=8)             :: relhh   ! Bounded relative humidity              [    --]
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.d0,max(0.d0,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure. -----------------------------------------------------!
      pvap  = relhh * eslf8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapl8 = max(toodry8, ep8 * pvap / (pres - (1.d0 - ep8) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapl8 = max(toodry8, ep8 * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function ptrh2rvapl8
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
   real(kind=8) function ptrh2rvapi8(relh,pres,temp,out_shv)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: relh    ! Relative humidity                      [    --]
      real(kind=8), intent(in) :: pres    ! Pressure                               [    Pa]
      real(kind=8), intent(in) :: temp    ! Temperature                            [     K]
      logical     , intent(in) :: out_shv ! Output is specific humidity            [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=8)             :: relhh   ! Bounded relative humidity              [    --]
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.d0,max(0.d0,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure. -----------------------------------------------------!
      pvap  = relhh * esif8(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapi8 = max(toodry8, ep8 * pvap / (pres - (1.d0 - ep8) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapi8 = max(toodry8, ep8 * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function ptrh2rvapi8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based (or specific humidity) based !
   ! on the pressure [Pa], temperature [K] and relative humidity [fraction].  It checks    !
   ! the temperature to decide between ice or liquid saturation.                           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function ptrh2rvapil8(relh,pres,temp,out_shv,useice)
      use rconstants , only : ep8      & ! intent(in)
                            , toodry8  & ! intent(in)
                            , t3ple8   ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: relh    ! Relative humidity            [    --]
      real(kind=8), intent(in)           :: pres    ! Pressure                     [    Pa]
      real(kind=8), intent(in)           :: temp    ! Temperature                  [     K]
      logical     , intent(in)           :: out_shv ! Output is specific humidity  [   T|F]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=8)                       :: relhh   ! Bounded relative humidity    [    --]
      logical                            :: frozen  ! Will use ice thermodynamics  [   T|F]
      !------------------------------------------------------------------------------------!


      !----- Check whether to use the user's or the default flag for ice saturation. ------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.d0,max(0.d0,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure (ice or liquid, depending on the value of frozen). ---!
      if (frozen) then
         pvap  = relhh * esif8(temp)
      else
         pvap  = relhh * eslf8(temp)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapil8 = max(toodry8, ep8 * pvap / (pres - (1.d0 - ep8) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapil8 = max(toodry8, ep8 * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
      return
   end function ptrh2rvapil8
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
   real(kind=8) function rehul8(pres,temp,humi,is_shv)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres    ! Air pressure                           [    Pa]
      real(kind=8), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=8), intent(in) :: humi    ! Humidity                               [ kg/kg]
      logical     , intent(in) :: is_shv  ! Input humidity is specific humidity    [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: shv     ! Specific humidity                      [ kg/kg]
      real(kind=8)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=8)             :: psat    ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry8,humi)
      else
         shv = max(toodry8,humi) / ( 1.d0 + max(toodry8,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep8 + (1.d0 - ep8) * shv )
      psat = eslf8(temp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehul8 = max(0.d0 , pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehul8
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
   real(kind=8) function rehui8(pres,temp,humi,is_shv)
      use rconstants , only : ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres    ! Air pressure                           [    Pa]
      real(kind=8), intent(in) :: temp    ! Temperature                            [     K]
      real(kind=8), intent(in) :: humi    ! Humidity                               [ kg/kg]
      logical     , intent(in) :: is_shv  ! Input humidity is specific humidity    [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: shv     ! Specific humidity                      [ kg/kg]
      real(kind=8)             :: pvap    ! Vapour pressure                        [    Pa]
      real(kind=8)             :: psat    ! Saturation vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry8,humi)
      else
         shv = max(toodry8,humi) / ( 1.d0 + max(toodry8,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep8 + (1.d0 - ep8) * shv )
      psat = esif8(temp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehui8 = max(0.d0 , pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehui8
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
   real(kind=8) function rehuil8(pres,temp,humi,is_shv,useice)
      use rconstants , only : t3ple8  & ! intent(in)
                            , ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres    ! Air pressure                 [    Pa]
      real(kind=8), intent(in)           :: temp    ! Temperature                  [     K]
      real(kind=8), intent(in)           :: humi    ! Humidity                     [ kg/kg]
      logical     , intent(in)           :: is_shv  ! Input is specific humidity   [   T|F]
       !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: shv     ! Specific humidity            [ kg/kg]
      real(kind=8)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=8)                       :: psat    ! Saturation vapour pressure   [    Pa]
      logical                            :: frozen  ! Will use ice saturation now  [   T|F]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    Check whether we should use ice or liquid saturation.                           !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else 
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry8,humi)
      else
         shv = max(toodry8,humi) / ( 1.d0 + max(toodry8,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep8 + (1.d0 - ep8) * shv )
      if (frozen) then
         psat = esif8(temp)
      else
         psat = esif8(temp)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehuil8 = max(0.d0 ,pvap / psat)
      !------------------------------------------------------------------------------------!

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
   real(kind=8) function tv2temp8(tvir,rvap,rtot)
      use rconstants , only : epi8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: tvir     ! Virtual temperature          [    K]
      real(kind=8), intent(in)           :: rvap     ! Vapour mixing ratio          [kg/kg]
      real(kind=8), intent(in), optional :: rtot     ! Total mixing ratio           [kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: rtothere ! Total or vapour mixing ratio [kg/kg]
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
      tv2temp8 = tvir * (1.d0 + rtothere) / (1.d0 + epi8 * rvap)
      !------------------------------------------------------------------------------------!

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
   real(kind=8) function virtt8(temp,rvap,rtot)
      use rconstants , only: epi8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: temp     ! Temperature                  [    K]
      real(kind=8), intent(in)           :: rvap     ! Vapour mixing ratio          [kg/kg]
      real(kind=8), intent(in), optional :: rtot     ! Total mixing ratio           [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=8)                       :: rtothere ! Total or vapour mixing ratio [kg/kg]
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
      virtt8 = temp * (1.d0 + epi8 * rvap) / (1.d0 + rtothere)
      !------------------------------------------------------------------------------------!

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
   real(kind=8) function idealdens8(pres,temp,rvap,rtot)
      use rconstants , only : rdry8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres ! Pressure                        [    Pa]
      real(kind=8), intent(in)           :: temp ! Temperature                     [     K]
      real(kind=8), intent(in)           :: rvap ! Vapour mixing ratio             [ kg/kg]
      real(kind=8), intent(in), optional :: rtot ! Total mixing ratio              [ kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=8)                       :: tvir ! Virtual temperature             [     K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Prefer using total mixing ratio, but if it isn't provided, then use vapour    !
      ! as total (no condensation).                                                        !
      !------------------------------------------------------------------------------------!
      if (present(rtot)) then
        tvir = virtt8(temp,rvap,rtot)
      else
        tvir = virtt8(temp,rvap)
      end if
      !------------------------------------------------------------------------------------!


      !----- Convert using the definition of virtual temperature. -------------------------!
      idealdens8 = pres / (rdry8 * tvir)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : rdry8 & ! intent(in)
                            , epi8  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres ! Pressure                        [    Pa]
      real(kind=8), intent(in)           :: temp ! Temperature                     [     K]
      real(kind=8), intent(in)           :: qvpr ! Vapour specific mass            [ kg/kg]
      real(kind=8), intent(in), optional :: qtot ! Total water specific mass       [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: qall ! Either qtot or qvpr...          [ kg/kg]
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
      idealdenssh8 = pres / (rdry8 * temp * (1.d0 - qall + epi8 * qvpr))
      !------------------------------------------------------------------------------------!

      return
   end function idealdenssh8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes reduces the pressure from the reference height to the      !
   ! canopy height by assuming hydrostatic equilibrium.  For simplicity, we assume that    !
   ! R and cp are constants (in reality they are dependent on humidity).                   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function reducedpress8(pres,thetaref,shvref,zref,thetacan,shvcan,zcan)
      use rconstants , only : epim18    & ! intent(in)
                            , p00k8     & ! intent(in)
                            , rocp8     & ! intent(in)
                            , cpor8     & ! intent(in)
                            , cpdry8    & ! intent(in)
                            , grav8     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres     ! Pressure                            [      Pa]
      real(kind=8), intent(in) :: thetaref ! Potential temperature               [       K]
      real(kind=8), intent(in) :: shvref   ! Vapour specific mass                [   kg/kg]
      real(kind=8), intent(in) :: zref     ! Height at reference level           [       m]
      real(kind=8), intent(in) :: thetacan ! Potential temperature               [       K]
      real(kind=8), intent(in) :: shvcan   ! Vapour specific mass                [   kg/kg]
      real(kind=8), intent(in) :: zcan     ! Height at canopy level              [       m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=8)             :: pinc     ! Pressure increment                  [ Pa^R/cp]
      real(kind=8)             :: thvbar   ! Average virtual pot. temperature    [       K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      First we compute the average virtual potential temperature between the canopy !
      ! top and the reference level.                                                       !
      !------------------------------------------------------------------------------------!
      thvbar = 5.d-1 * ( thetaref * (1.d0 + epim18 * shvref)                               &
                       + thetacan * (1.d0 + epim18 * shvcan) )
      !------------------------------------------------------------------------------------!



      !----- Then, we find the pressure gradient scale. -----------------------------------!
      pinc   = grav8 * p00k8 * (zref - zcan) / (cpdry8 * thvbar)
      !------------------------------------------------------------------------------------!


      !----- And we can find the reduced pressure. ----------------------------------------!
      reducedpress8 = (pres**rocp8 + pinc ) ** cpor8
      !------------------------------------------------------------------------------------!

      return
   end function reducedpress8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the Exner function [J/kg/K], given the pressure.  It       !
   ! assumes for simplicity that R and Cp are constants and equal to the dry air values.   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function press2exner8(pres)
      use rconstants , only : p00i8           & ! intent(in)
                            , cpdry8          & ! intent(in)
                            , rocp8           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres   ! Pressure                               [     Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      press2exner8 = cpdry8 * ( pres * p00i8 ) ** rocp8
      !------------------------------------------------------------------------------------!

      return
   end function press2exner8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the pressure [Pa], given the Exner function.  Like in the  !
   ! function above, we also assume R and Cp to be constants and equal to the dry air      !
   ! values.                                                                               !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function exner2press8(exner)
      use rconstants , only : p008            & ! intent(in)
                            , cpdryi8         & ! intent(in)
                            , cpor8           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      exner2press8 = p008 * ( exner * cpdryi8 ) ** cpor8
      !------------------------------------------------------------------------------------!

      return
   end function exner2press8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the potential temperature [K], given the Exner function    !
   ! and temperature.  For simplicity we ignore the effects of humidity in R and cp and    !
   ! use the dry air values instead.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function extemp2theta8(exner,temp)
      use rconstants , only : cpdry8          ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=8), intent(in) :: temp   ! Temperature                            [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extemp2theta8 = cpdry8 * temp / exner
      !------------------------------------------------------------------------------------!

      return
   end function extemp2theta8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature [K], given the Exner function and          !
   ! potential temperature.  We simplify the equations by assuming that R and Cp are       !
   ! constants.                                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function extheta2temp8(exner,theta)
      use rconstants , only : p00i8           & ! intent(in)
                            , cpdryi8         ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=8), intent(in) :: theta  ! Potential temperature                  [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extheta2temp8 = cpdryi8 * exner * theta
      !------------------------------------------------------------------------------------!

      return
   end function extheta2temp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the specific internal energy of water [J/kg], given the    !
   ! temperature and liquid fraction.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function tl2uint8(temp,fliq)
      use rconstants , only : cice8           & ! intent(in)
                            , cliq8           & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)
                            
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=8), intent(in) :: fliq  ! Fraction liquid water                    [ kg/kg]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of ice and liquid     !
      ! phases.                                                                            !
      !------------------------------------------------------------------------------------!
      tl2uint8 = (1.d0 - fliq) * cice8 * temp + fliq * cliq8 * (temp - tsupercool_liq8)
      !------------------------------------------------------------------------------------!

      return
   end function tl2uint8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the internal energy of water [J/m²] or [  J/m³], given the !
   ! temperature [K], the heat capacity of the "dry" part [J/m²/K] or [J/m³/K], water mass !
   ! [ kg/m²] or [ kg/m³], and liquid fraction [---].                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function cmtl2uext8(dryhcap,wmass,temp,fliq)
      use rconstants , only : cice8           & ! intent(in)
                            , cliq8           & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)
                            
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m²/K] or [J/m³/K]
      real(kind=8), intent(in)  :: wmass   ! Mass                      [ kg/m²] or [ kg/m³]
      real(kind=8), intent(in)  :: temp    ! Temperature                           [     K]
      real(kind=8), intent(in)  :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of dry part, plus the !
      ! contribution of ice and liquid phases.                                             !
      !------------------------------------------------------------------------------------!
      cmtl2uext8 = dryhcap * temp + wmass * ( (1.d0 - fliq) * cice8 * temp                 &
                                            + fliq * cliq8 * (temp - tsupercool_liq8) )
      !------------------------------------------------------------------------------------!

      return
   end function cmtl2uext8
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
   real(kind=8) function tq2enthalpy8(temp,humi,is_shv)
      use rconstants , only : cpdry8          & ! intent(in)
                            , cph2o8          & ! intent(in)
                            , tsupercool_vap8 ! ! intent(in)
                            
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp   ! Temperature                             [     K]
      real(kind=8), intent(in) :: humi   ! Humidity (spec. hum. or mixing ratio)   [ kg/kg]
      logical     , intent(in) :: is_shv ! Input humidity is specific humidity     [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: shv    ! Specific humidity                       [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.d0)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      tq2enthalpy8 = (1.d0 - shv) * cpdry8 * temp + shv * cph2o8 * (temp - tsupercool_vap8)
      !------------------------------------------------------------------------------------!

      return
   end function tq2enthalpy8
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
   real(kind=8) function hq2temp8(enthalpy,humi,is_shv)
      use rconstants , only : cpdry8          & ! intent(in)
                            , cph2o8          & ! intent(in)
                            , tsupercool_vap8 ! ! intent(in)
                            
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: enthalpy ! Specific enthalpy                     [  J/kg]
      real(kind=8), intent(in) :: humi     ! Humidity (spec. hum. or mixing ratio) [ kg/kg]
      logical     , intent(in) :: is_shv   ! Input humidity is specific humidity   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: shv      ! Specific humidity                     [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.d0)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      hq2temp8 = ( enthalpy + shv * cph2o8 * tsupercool_vap8 )                             &
               / ( (1.d0 - shv) * cpdry8 + shv * cph2o8 )
      !------------------------------------------------------------------------------------!

      return
   end function hq2temp8
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of vaporisation for a given temperature.  If  !
   ! we use the definition of latent heat (difference in enthalpy between liquid and       !
   ! vapour phases), and assume that the specific heats are constants, latent heat becomes !
   ! a linear function of temperature.                                                     !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function alvl8(temp)
      use rconstants , only : alvl38  & ! intent(in)
                            , dcpvl8  & ! intent(in)
                            , t3ple8  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvl8 = alvl38 + dcpvl8 * (temp - t3ple8)
      !------------------------------------------------------------------------------------!

      return
   end function alvl8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of sublimation for a given temperature.  If   !
   ! we use the definition of latent heat (difference in enthalpy between ice and vapour   !
   ! phases), and assume that the specific heats are constants, latent heat becomes a      !
   ! linear function of temperature.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function alvi8(temp)
      use rconstants , only : alvi38  & ! intent(in)
                            , dcpvi8  & ! intent(in)
                            , t3ple8  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvi8 = alvi38 + dcpvi8 * (temp - t3ple8)
      !------------------------------------------------------------------------------------!

      return
   end function alvi8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the ice liquid potential temperature given the Exner       !
   ! function [J/kg/K], temperature [K], and liquid and ice mixing ratios [kg/kg].         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function theta_iceliq8(exner,temp,rliq,rice)
      use rconstants , only : alvl38      & ! intent(in)
                            , alvi38      & ! intent(in)
                            , cpdry8      & ! intent(in)
                            , ttripoli8   & ! intent(in)
                            , htripoli8   & ! intent(in)
                            , htripolii8  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner ! Exner function                          [ J/kg/K]
      real(kind=8), intent(in) :: temp  ! Temperature                             [      K]
      real(kind=8), intent(in) :: rliq  ! Liquid mixing ratio                     [  kg/kg]
      real(kind=8), intent(in) :: rice  ! Ice mixing ratio                        [  kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: hh    ! Enthalpy associated with sensible heat  [   J/kg]
      real(kind=8)             :: qq    ! Enthalpy associated with latent heat    [   J/kg]
      !------------------------------------------------------------------------------------!


      !----- Find the sensible heat enthalpy (assuming dry air). --------------------------!
      hh = cpdry8 * temp
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the latent heat enthalpy.  If using the old thermodynamics, we use the   !
      ! latent heat at T = T3ple, otherwise we use the temperature-dependent one.          !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         qq = alvl8(temp) * rliq + alvi8(temp) * rice
      else
         qq = alvl38 * rliq + alvi38 * rice
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Solve the thermodynamics.  For the new thermodynamics we don't approximate    !
      ! the exponential to a linear function, nor do we impose temperature above the thre- !
      ! shold from Tripoli and Cotton (1981).                                              !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Decide how to compute, based on temperature. ------------------------------!
         theta_iceliq8 = hh * exp(-qq / hh) / exner
         !---------------------------------------------------------------------------------!
      else
         !----- Decide how to compute, based on temperature. ------------------------------!
         if (temp > ttripoli8) then
            theta_iceliq8 = hh * hh / (exner * ( hh + qq))
         else
            theta_iceliq8 = hh * htripoli8 / (exner * ( htripoli8 + qq))
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : alvl38      & ! intent(in)
                            , alvi38      & ! intent(in)
                            , dcpvi8      & ! intent(in)
                            , dcpvl8      & ! intent(in)
                            , cpdry8      & ! intent(in)
                            , ttripoli8   & ! intent(in)
                            , htripoli8   & ! intent(in)
                            , htripolii8  & ! intent(in)
                            , t3ple8      ! ! intent(in)

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
      real(kind=8)                       :: rdlt       ! r × d(L)/dT × T           [  J/kg]
      real(kind=8)                       :: hh         ! Sensible heat enthalpy    [  J/kg]
      real(kind=8)                       :: qq         ! Latent heat enthalpy      [  J/kg]
      logical                            :: thereisice ! Is ice present            [   ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Check whether we should consider ice thermodynamics or not.                      !
      !------------------------------------------------------------------------------------!
      thereisice = present(ricein)
      if (thereisice) then
         rice = ricein
      else
         rice = 0.d0
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether the current state has condensed water.                           !
      !------------------------------------------------------------------------------------!
      if (rliq+rice == 0.d0) then
         !----- No condensation, so dthetail_dt is a constant. ----------------------------!
         dthetail_dt8 = thil/temp
         return
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Condensation exists.  Compute some auxiliary variables.                     !
         !---------------------------------------------------------------------------------!


         !---- Sensible heat enthalpy. ----------------------------------------------------!
         hh = cpdry8 * temp
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Find the latent heat enthalpy.  If using the old thermodynamics, we use    !
         ! the latent heat at T = T3ple, otherwise we use the temperature-dependent one.   !
         ! The term r × d(L)/dT × T is computed only when we use the new thermodynamics.   !
         !---------------------------------------------------------------------------------!
         if (newthermo) then
            qq   = alvl8(temp) * rliq + alvi8(temp) * rice
            rdlt = (dcpvl8 * rliq + dcpvi8 * rice ) * temp
         else
            qq   = alvl38 * rliq + alvi38 * rice
            rdlt = 0.d0
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
            ldrst = 0.d0
         elseif (thereisice .and. temp < t3ple8) then
            if (newthermo) then
               ldrst = alvi38 * rsifp8(pres,temp) * temp
            else
               ldrst = alvi8(temp) * rsifp8(pres,temp) * temp
            end if
         else
            if (newthermo) then
               ldrst = alvl38 * rslfp8(pres,temp) * temp
            else
               ldrst = alvl8(temp) * rslfp8(pres,temp) * temp
            end if
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the condensed phase consistent with the thermodynamics used.              !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         dthetail_dt8 = thil * ( 1.d0 + (ldrst + qq - rdlt ) / hh ) / temp
      else
         !----- Decide how to compute, based on temperature. ------------------------------!
         if (temp > ttripoli8) then
            dthetail_dt8 = thil * ( 1.d0 + (ldrst + qq) / (hh+qq) ) / temp
         else
            dthetail_dt8 = thil * ( 1.d0 + ldrst / (htripoli8 + alvl38 * rliq) ) / temp
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : cpdry8      & ! intent(in)
                            , cpdryi8     & ! intent(in)
                            , cpdryi48    & ! intent(in)
                            , alvl38      & ! intent(in)
                            , alvi38      & ! intent(in)
                            , t008        & ! intent(in)
                            , t3ple8      & ! intent(in)
                            , ttripoli8   & ! intent(in)
                            , htripolii8  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: thil      ! Ice-liquid water potential temp.     [     K]
      real(kind=8), intent(in) :: exner     ! Exner function                       [J/kg/K]
      real(kind=8), intent(in) :: pres      ! Pressure                             [    Pa]
      real(kind=8), intent(in) :: rliq      ! Liquid water mixing ratio            [ kg/kg]
      real(kind=8), intent(in) :: rice      ! Ice mixing ratio                     [ kg/kg]
      real(kind=8), intent(in) :: t1stguess ! 1st. guess for temperature           [     K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: til       ! Ice liquid temperature               [     K]
      real(kind=8)             :: deriv     ! Function derivative 
      real(kind=8)             :: fun       ! Function for which we seek a root.
      real(kind=8)             :: funa      ! Smallest  guess function
      real(kind=8)             :: funz      ! Largest   guess function
      real(kind=8)             :: tempa     ! Smallest  guess (or previous guess in Newton)
      real(kind=8)             :: tempz     ! Largest   guess (or new guess in Newton)
      real(kind=8)             :: delta     ! Aux. var to compute 2nd guess for bisection
      integer                  :: itn       ! Iteration counter
      integer                  :: itb       ! Iteration counter
      logical                  :: converged ! Convergence flag
      logical                  :: zside     ! Flag to check for one-sided approach...
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     First we check for conditions that don't require iterative root-finding.       !
      !------------------------------------------------------------------------------------!
      if (rliq + rice == 0.d0) then
         !----- No condensation.  Theta_il is the same as theta. --------------------------!
         thil2temp8 = cpdryi8 * thil * exner
         return
         !---------------------------------------------------------------------------------!
      elseif (.not. newthermo) then
         !---------------------------------------------------------------------------------!
         !    There is condensation but we are using the old thermodynamics, which can be  !
         ! solved analytically.                                                            !
         !---------------------------------------------------------------------------------!
         til = cpdryi8 * thil * exner
         if (t1stguess > ttripoli8) then
            thil2temp8 = 5.d-1                                                             &
                       * (til + sqrt( til                                                  &
                                    * (til + cpdryi48 * (alvl38 * rliq + alvi38 * rice))))
         else
            thil2temp8 = til * ( 1.d0 + (alvl38 * rliq + alvi38 * rice) * htripolii8)
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
      fun       = theta_iceliq8(exner,tempz,rliq,rice)
      deriv     = dthetail_dt8(.true.,fun,exner,pres,tempz,rliq,rice)
      fun       = fun - thil
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,2(1x,a,1x,es12.5))')            &
      !   'itn=',0,'bisection=',.false.,'tempz=',tempz-t00                                  &
      !           ,'fun=',fun,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (fun == 0.d0) then
         thil2temp8 = tempz
         converged = .true.
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

         converged = abs(tempa-tempz) < toler8*tempz
         !----- Converged, happy with that, return the average b/w the 2 previous guesses -!
         if (fun == 0.d0) then
            thil2temp8 = tempz
            converged  = .true.
            return
         elseif(converged) then
            thil2temp8 = 5.d-1 * (tempa+tempz)
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If we have reached this point then Newton's method failed.  Use bisection      !
      ! instead.  For bisection, We need two guesses whose function evaluations have       !
      ! opposite sign.                                                                     !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         !----- Guesses have opposite sign. -----------------------------------------------!
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < toler8 * tempa) then
            delta = 1.d2 * toler8 * tempa
         else 
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),1.d2 * toler8 * tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
            funz  = theta_iceliq8(exner,tempz,rliq,rice) - thil
            zside = funa*funz < 0.d0
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
                            ,'thil2temp8','therm_lib8.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         thil2temp8 =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(thil2temp8 - tempa) < toler8 * thil2temp8
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
            if (zside) funa = funa * 5.d-1
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa  = thil2temp8
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 5.d-1
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
         write (unit=*,fmt='(a,1x,es12.4)') 'toler8          [  ----] =',toler8
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                           , abs(thil2temp8-tempa)/thil2temp8
         write (unit=*,fmt='(a,1x,f12.4)' ) 'thil2temp8      [     K] =',thil2temp8

         call abort_run  ('Temperature didn''t converge, giving up!!!'                     &
                         ,'thil2temp8','therm_lib8.f90')
      end if

      return
   end function thil2temp8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the partial derivative of temperature as a function of the !
   ! saturation mixing ratio [kg/kg], keeping pressure constant.  This is based on the     !
   ! ice-liquid potential temperature equation.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dtempdrs8(exner,thil,temp,rliq,rice,rconmin)
      use rconstants , only : alvl38      & ! intent(in)
                            , alvi38      & ! intent(in)
                            , dcpvl8      & ! intent(in)
                            , dcpvi8      & ! intent(in)
                            , cpdry8      & ! intent(in)
                            , cpdryi8     & ! intent(in)
                            , ttripoli8   & ! intent(in)
                            , htripolii8  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner   ! Exner function                        [ J/kg/K]
      real(kind=8), intent(in) :: thil    ! Ice-liquid potential temperature (*)  [      K]
      real(kind=8), intent(in) :: temp    ! Temperature                           [      K]
      real(kind=8), intent(in) :: rliq    ! Liquid mixing ratio                   [  kg/kg]
      real(kind=8), intent(in) :: rice    ! Ice mixing ratio                      [  kg/kg]
      real(kind=8), intent(in) :: rconmin ! Min. non-zero condensate mix. ratio   [  kg/kg]
      !------------------------------------------------------------------------------------!
      ! (*) Thil is not used in this formulation but it may be used should you opt by      !
      !     other ways to compute theta_il, so don't remove this argument.                 !
      !------------------------------------------------------------------------------------!
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: qq      ! Enthalpy -- latent heat               [   J/kg]
      real(kind=8)             :: qpt     ! d(qq)/dT * T                          [   J/kg]
      real(kind=8)             :: hh      ! Enthalpy -- sensible heat             [   J/kg]
      real(kind=8)             :: rcon    ! Condensate mixing ratio               [  kg/kg]
      real(kind=8)             :: til     ! Ice-liquid temperature                [      K]
      !------------------------------------------------------------------------------------!


      !----- Find the total hydrometeor mixing ratio. -------------------------------------!
      rcon  = rliq+rice
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      If the amount of condensate is negligible, temperature does not depend on     !
      ! saturation mixing ratio.                                                           !
      !------------------------------------------------------------------------------------!
      if (rcon < rconmin) then
         dtempdrs8 = 0.d0
      else

         !---------------------------------------------------------------------------------!
         !     Find the enthalpy associated with latent heat and its derivative            !
         ! correction.  This is dependent on the thermodynamics used.                      !
         !---------------------------------------------------------------------------------!
         if (newthermo) then
            qq  = alvl8(temp) * rliq + alvi8(temp) * rice
            qpt = dcpvl8 * rliq + dcpvi8 * rice
         else
            qq  = alvl38 * rliq + alvi38 * rice
            qpt = 0.d0
         end if
         !---------------------------------------------------------------------------------!


         !----- Find the enthalpy associated with sensible heat. --------------------------!
         hh = cpdry8 * temp
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Decide how to compute, based on the thermodynamics method.                   !
         !---------------------------------------------------------------------------------!
         if (newthermo) then
            dtempdrs8 = - temp * qq / ( rcon * (hh + qq - qpt) )
         else
            til   = cpdryi8 * thil * exner
            !----- Decide how to compute, based on temperature. ---------------------------!
            if (temp > ttripoli8) then
               dtempdrs8 = - til * qq / ( rcon * cpdry8 * (2.*temp-til))
            else
               dtempdrs8 = - til * qq * htripolii8 / rcon
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function dtempdrs8
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
   real(kind=8) function thetaeiv8(thil,pres,temp,rvap,rtot,useice)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: thil   ! Ice-liquid potential temp.    [     K]
      real(kind=8), intent(in)           :: pres   ! Pressure                      [    Pa]
      real(kind=8), intent(in)           :: temp   ! Temperature                   [     K]
      real(kind=8), intent(in)           :: rvap   ! Water vapour mixing ratio     [ kg/kg]
      real(kind=8), intent(in)           :: rtot   ! Total mixing ratio            [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! Should I use ice?             [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: tlcl   ! Internal LCL temperature      [     K]
      real(kind=8)                       :: plcl   ! Lifting condensation pressure [    Pa]
      real(kind=8)                       :: dzlcl  ! Thickness of lyr. beneath LCL [     m]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the liquid condensation level (LCL).                                      !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         call lcl_il8(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      else
         call lcl_il8(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The definition of the thetae_iv is the thetae_ivs at the LCL. The LCL, in turn !
      ! is the point in which rtot = rvap = rsat, so at the LCL rliq = rice = 0.           !
      !------------------------------------------------------------------------------------!
      thetaeiv8  = thetaeivs8(thil,tlcl,rtot,0.d0,0.d0)
      !------------------------------------------------------------------------------------!

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
   !                 pressure at the LCL) is a function of temperature.  In case you want  !
   !                 d(Thetae_ivs)/dT, use the dthetaeivs_dt function instead.             !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dthetaeiv_dtlcl8(theiv,tlcl,rtot,eslcl,useice)
      use rconstants , only : rocp8       & ! intent(in)
                            , cpdry8      & ! intent(in)
                            , dcpvl8      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: theiv    ! Ice-vap. equiv. pot. temp. [      K]
      real(kind=8), intent(in)           :: tlcl     ! LCL temperature            [      K]
      real(kind=8), intent(in)           :: rtot     ! Total mixing ratio         [  kg/kg]
      real(kind=8), intent(in)           :: eslcl    ! LCL sat. vapour pressure   [     Pa]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice   ! Flag for considering ice   [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: desdtlcl ! Sat. vapour pres. deriv.   [   Pa/K]
      real(kind=8)                       :: esterm   ! es(TLC) term               [   ----]
      real(kind=8)                       :: hhlcl    ! Enthalpy -- sensible       [   J/kg]
      real(kind=8)                       :: qqlcl    ! Enthalpy -- latent         [   J/kg]
      real(kind=8)                       :: qptlcl   ! Latent deriv. * T_LCL      [   J/kg]
      !------------------------------------------------------------------------------------!



      !----- Find the derivative of rs with temperature. ----------------------------------!
      if (present(useice)) then
         desdtlcl = eslifp8(tlcl,useice)
      else
         desdtlcl = eslifp8(tlcl)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Saturation term.                                                               !
      !------------------------------------------------------------------------------------!
      esterm = rocp8 * tlcl * desdtlcl / eslcl
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the enthalpy terms.                                                       !
      !------------------------------------------------------------------------------------!
      hhlcl  = cpdry8 * tlcl
      qqlcl  = alvl8(tlcl) * rtot
      qptlcl = dcpvl8 * rtot * tlcl
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Derivative.                                                                   !
      !------------------------------------------------------------------------------------!
      dthetaeiv_dtlcl8 = theiv / tlcl * (1.d0 - esterm - (qqlcl - qptlcl) / hhlcl)
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : cpdry8    ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: thil  ! Theta_il, ice-liquid water pot. temp.    [     K]
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=8), intent(in) :: rsat  ! Saturation water vapour mixing ratio     [ kg/kg]
      real(kind=8), intent(in) :: rliq  ! Liquid water mixing ratio                [ kg/kg]
      real(kind=8), intent(in) :: rice  ! Ice mixing ratio                         [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)             :: rtots ! Saturated mixing ratio                   [     K]
      !------------------------------------------------------------------------------------!


      !------ Find the total saturation mixing ratio. -------------------------------------!
      rtots = rsat + rliq + rice
      !------------------------------------------------------------------------------------!


      !------ Find the saturation equivalent potential temperature. -----------------------!
      thetaeivs8 = thil * exp ( alvl8(temp) * rtots / (cpdry8 * temp))
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : cpdry8     & ! intent(in)
                            , dcpvl8     ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: theivs ! Sat. ice-vap. eq. pot. temp. [      K]
      real(kind=8), intent(in)           :: temp   ! Temperature                  [      K]
      real(kind=8), intent(in)           :: pres   ! Pressure                     [     Pa]
      real(kind=8), intent(in)           :: rsat   ! Saturation mixing ratio      [  kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! Flag for considering ice     [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: drsdt  ! Sat. mixing ratio derivative [kg/kg/K]
      real(kind=8)                       :: hh     ! Enthalpy -- sensible         [   J/kg]
      real(kind=8)                       :: qqaux  ! Enthalpy -- sensible         [   J/kg]
      !------------------------------------------------------------------------------------!


      !----- Find the derivative of rs with temperature. ----------------------------------!
      if (present(useice)) then
         drsdt = rslifp8(pres,temp,useice)
      else
         drsdt = rslifp8(pres,temp)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the enthalpy terms.                                                      !
      !------------------------------------------------------------------------------------!
      hh    = cpdry8 * temp
      qqaux = alvl8(temp) * (drsdt * temp - rsat) + dcpvl8 * rsat * temp
      !------------------------------------------------------------------------------------!


      !----- Find the derivative.  Depending on the temperature, use different eqn. -------!
      dthetaeivs_dt8 = theivs / temp * ( 1.d0 + qqaux / hh )
      !------------------------------------------------------------------------------------!

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
      use rconstants , only : ep8       & ! intent(in)
                            , cpdry8    & ! intent(in)
                            , p008      & ! intent(in)
                            , rocp8     & ! intent(in)
                            , t3ple8    & ! intent(in)
                            , t008      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: theiv     ! Ice vap. equiv. pot. temp. [     K]
      real(kind=8), intent(in)           :: pres      ! Pressure                   [    Pa]
      real(kind=8), intent(in)           :: rtot      ! Total mixing ratio         [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice    ! May I use ice thermodyn.   [   T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)                       :: pvap      ! Sat. vapour pressure
      real(kind=8)                       :: theta     ! Potential temperature
      real(kind=8)                       :: deriv     ! Function derivative 
      real(kind=8)                       :: funnow    ! Function for which we seek a root.
      real(kind=8)                       :: funa      ! Smallest guess function
      real(kind=8)                       :: funz      ! Largest guess function
      real(kind=8)                       :: tlcla     ! Smallest guess (Newton: old guess)
      real(kind=8)                       :: tlclz     ! Largest guess (Newton: new guess)
      real(kind=8)                       :: tlcl      ! What will be the LCL temperature
      real(kind=8)                       :: es00      ! Defined as p00*rt/(epsilon + rt)
      real(kind=8)                       :: delta     ! Aux. variable (For 2nd guess).
      integer                            :: itn       ! Iteration counters
      integer                            :: itb       ! Iteration counters
      integer                            :: ii        ! Another counter
      logical                            :: converged ! Convergence handle
      logical                            :: zside     ! Side checker for Regula Falsi
      logical                            :: frozen    ! Will use ice thermodynamics
      !------------------------------------------------------------------------------------!



      !----- Fill the flag for ice thermodynamics so it will be present. ------------------!
      if (present(useice)) then
         frozen = useice
      else
         frozen = bulk_on
      end if
      !------------------------------------------------------------------------------------!



      !----- Find es00, which is a constant. ----------------------------------------------!
      es00 = p008 * rtot / (ep8 + rtot)
      !------------------------------------------------------------------------------------!


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
      tlclz     = tslif8(pvap,frozen)
      theta     = tlclz * (es00 / pvap) ** rocp8
      funnow    = thetaeivs8(theta,tlclz,rtot,0.d0,0.d0)
      deriv     = dthetaeiv_dtlcl8(funnow,tlclz,rtot,pvap,frozen)
      funnow    = funnow - theiv
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')                &
      !   'NEWTON: it=',0,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Put something in tlcla in case we never loop through Newton's method. --------!
      tlcla     = tlclz
      funa      = funnow
      converged = .false.
      !------------------------------------------------------------------------------------!

      !----- Looping: Newton's iterative method -------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, skip to bisection ----!
         !----- Update guesses. -----------------------------------------------------------!
         tlcla   = tlclz
         funa    = funnow
         tlclz   = tlcla - funnow/deriv
         !---------------------------------------------------------------------------------!


         !----- Update the function evaluation and its derivative. ------------------------!
         pvap    = eslif8(tlclz,frozen)
         theta   = tlclz * (es00/pvap)**rocp8
         funnow  = thetaeivs8(theta,tlclz,rtot,0.d0,0.d0)
         deriv   = dthetaeiv_dtlcl8(funnow,tlclz,rtot,pvap,frozen)
         funnow  = funnow - theiv
         !---------------------------------------------------------------------------------!

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
            tlcl = 5.d-1*(tlcla+tlclz)
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
         if (funa*funnow > 0.d0) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler8*tlcla) then
               delta = 1.d2*toler8*tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funz-funa)),1.d2*toler8*tlcla)
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
               pvap  = eslif8(tlclz,frozen)
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
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thetaeiv2thil8','therm_lib8.f90')
            end if
         end if
         !---- Continue iterative method. -------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo

            !----- Update the guess. ------------------------------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)

            !----- Updating function evaluation -------------------------------------------!
            pvap   = eslif8(tlcl,frozen)
            theta  = tlcl * (es00/pvap)**rocp8
            funnow = thetaeivs8(theta,tlcl,rtot,0.d0,0.d0) - theiv

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),1(a,1x,es11.4,1x))')          &
            !   'REGFAL: it=',itb,'tlcl =',tlcl -t00,'pvap=',0.01*pvap,'fun=',funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Check for convergence. If it did, return, we found the solution.          !
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
         thetaeiv2thil8  = theiv * exp (- alvl8(tlcl) * rtot / (cpdry8 * tlcl) )
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

         call abort_run  ('TLCL didn''t converge, qgave up!'                               &
                         ,'thetaeiv2thil8','therm_lib8.f90')
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
      use rconstants , only : cpdry8     & ! intent(in)
                            , ep8        & ! intent(in)
                            , p008       & ! intent(in)
                            , rocp8      & ! intent(in)
                            , t008       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)            :: theivs     ! Sat. thetae_iv          [      K]
      real(kind=8), intent(in)            :: pres       ! Pressure                [     Pa]
      real(kind=8), intent(out)           :: theta      ! Potential temperature   [      K]
      real(kind=8), intent(out)           :: temp       ! Temperature             [      K]
      real(kind=8), intent(out)           :: rsat       ! Saturation mixing ratio [  kg/kg]
      logical     , intent(in) , optional :: useice     ! May use ice thermodyn.  [    T|F]
      !----- Local variables, with other thermodynamic properties -------------------------!
      real(kind=8)                        :: exnernormi ! 1./ (Norm. Exner func.) [    ---]
      logical                             :: frozen     ! Will use ice thermodyn. [    T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real(kind=8)                        :: deriv      ! Function derivative 
      real(kind=8)                        :: funnow     ! Current function evaluation
      real(kind=8)                        :: funa       ! Smallest guess function
      real(kind=8)                        :: funz       ! Largest guess function
      real(kind=8)                        :: tempa      ! Smallest guess (Newton: previous)
      real(kind=8)                        :: tempz      ! Largest guess (Newton: new)
      real(kind=8)                        :: delta      ! Aux. variable for 2nd guess.
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
      exnernormi = (p008 /pres) ** rocp8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The 1st. guess, no idea, guess 0°C.                                            !
      !------------------------------------------------------------------------------------!
      tempz     = t008
      theta     = tempz * exnernormi
      rsat      = rslif8(pres,tempz,frozen)
      funnow    = thetaeivs8(theta,tempz,rsat,0.d0,0.d0)
      deriv     = dthetaeivs_dt8(funnow,tempz,pres,rsat,frozen)
      funnow    = funnow - theivs
      !------------------------------------------------------------------------------------!


      !----- Copy here just in case Newton is aborted at the 1st guess. -------------------!
      tempa     = tempz
      funa      = funnow
      !------------------------------------------------------------------------------------!

      converged = .false.
      !----- Newton's method loop. --------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler8) exit newloop !----- Too dangerous, skip to bisection ----!
         !----- Updating guesses ----------------------------------------------------------!
         tempa   = tempz
         funa    = funnow
         
         tempz   = tempa - funnow/deriv
         theta   = tempz * exnernormi
         rsat    = rslif8(pres,tempz,frozen)
         funnow  = thetaeivs8(theta,tempz,rsat,0.d0,0.d0)
         deriv   = dthetaeivs_dt8(funnow,tempz,pres,rsat,frozen)
         funnow  = funnow - theivs

         converged = abs(tempa-tempz) < toler8*tempz
         if (funnow == 0.d0) then
            converged =.true.
            temp = tempz
            exit newloop
         elseif (converged) then
            temp = 5.d-1*(tempa+tempz)
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

         if (funa*funnow > 0.d0) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler8*tempa) then
               delta = 1.d2*toler8*tempa
            else
               delta = max(abs(funa*(tempz-tempa)/(funz-funa)),1.d2*toler8*tempa)
            end if
            !------------------------------------------------------------------------------!

            tempz = tempa + delta
            zgssloop: do itb=1,maxfpo
               !----- So this will be +1 -1 +2 -2 etc. ------------------------------------!
               tempz = tempz + dble((-1)**itb * (itb+3)/2) * delta
               theta = tempz * exnernormi
               rsat  = rslif8(pres,tempz,frozen)
               funz  = thetaeivs8(theta,tempz,rsat,0.d0,0.d0) - theivs
               zside = funa*funz < 0.d0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thetaes2temp8','therm_lib8.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            if (abs(funz-funa) < toler8*tempa) then
               temp   = 5.d-1*(tempa+tempz)
            else
               temp   = (funz*tempa-funa*tempz)/(funz-funa)
            end if
            theta  = temp * exnernormi
            rsat   = rslif8(pres,temp,frozen)
            funnow = thetaeivs8(theta,temp,rsat,0.d0,0.d0) - theivs

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa) < toler8*temp
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
         rsat  = rslif8(pres,temp,frozen)
      else
         call abort_run  ('Temperature didn''t converge, I gave up!'                       &
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
   subroutine lcl_il8(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      use rconstants , only : cpog8     & ! intent(in)
                            , ep8       & ! intent(in)
                            , p008      & ! intent(in)
                            , rocp8     & ! intent(in)
                            , t3ple8    & ! intent(in)
                            , t008      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)            :: thil      ! Ice liquid pot. temp. (*)[      K]
      real(kind=8), intent(in)            :: pres      ! Pressure                 [     Pa]
      real(kind=8), intent(in)            :: temp      ! Temperature              [      K]
      real(kind=8), intent(in)            :: rtot      ! Total mixing ratio       [  kg/kg]
      real(kind=8), intent(in)            :: rvap      ! Vapour mixing ratio      [  kg/kg]
      real(kind=8), intent(out)           :: tlcl      ! LCL temperature          [      K]
      real(kind=8), intent(out)           :: plcl      ! LCL pressure             [     Pa]
      real(kind=8), intent(out)           :: dzlcl     ! Sub-LCL layer thickness  [      m]
      !------------------------------------------------------------------------------------!
      ! (*) This is the most general variable. Thil is exactly theta for no condensation   !
      !     condition, and it is the liquid potential temperature if no ice is present.    !
      !------------------------------------------------------------------------------------!
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in) , optional :: useice    ! May use ice thermodyn.?  [    T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                        :: pvap      ! Sat. vapour pressure
      real(kind=8)                        :: deriv     ! Function derivative 
      real(kind=8)                        :: funnow    ! Current function evaluation
      real(kind=8)                        :: funa      ! Smallest guess function
      real(kind=8)                        :: funz      ! Largest  guess function
      real(kind=8)                        :: tlcla     ! Smallest guess (Newton: previous)
      real(kind=8)                        :: tlclz     ! Largest guess (Newton: new)
      real(kind=8)                        :: es00      ! Defined as p00*rt/(epsilon + rt)
      real(kind=8)                        :: delta     ! Aux. variable for bisection
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
      es00 = p008 * rtot / (ep8 + rtot)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvap / (ep8 + rvap)
      tlclz     = 5.5d1 + 2.840d3 / (3.5d0 * log(temp) - log(1.d-2*pvap) - 4.805d0)
      pvap      = eslif8(tlclz,frozen)
      funnow    = tlclz * (es00/pvap)**rocp8 - thil
      deriv     = (funnow+thil)*(1.d0/tlclz - rocp8*eslifp8(tlclz,frozen)/pvap)
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
         if (abs(deriv) < toler8) exit newloop
         !---------------------------------------------------------------------------------!


         !----- Otherwise, update guesses. ------------------------------------------------!
         tlcla  = tlclz
         funa   = funnow
         tlclz  = tlcla - funnow/deriv
         pvap   = eslif8(tlclz,frozen)
         funnow = tlclz * (es00/pvap)**rocp8 - thil
         deriv  = (funnow+thil)*(1.d0/tlclz - rocp8*eslifp8(tlclz,frozen)/pvap)

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
         converged = abs(tlcla-tlclz) < toler8*tlclz
         if (converged) then
            !----- Guesses are almost identical, average them. ----------------------------!
            tlcl = 5.d-1*(tlcla+tlclz)
            funz = funnow
            exit newloop
            !------------------------------------------------------------------------------!
         elseif (funnow == 0.d0) then
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
         if (funa*funnow < 0.d0 ) then
            !----- We already have two good guesses. --------------------------------------!
            funz  = funnow
            zside = .true.
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     We need to find another guess with opposite sign.                        !
            !------------------------------------------------------------------------------!

            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funnow-funa) < toler8*tlcla) then
               delta = 1.d2*toler8*tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funnow-funa)),1.d2*toler8*tlcla)
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
               tlclz = tlcla + dble((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif8(tlclz,frozen)
               funz  = tlclz * (es00/pvap)**rocp8 - thil
               !---------------------------------------------------------------------------!


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
               write (unit=*,fmt='(a,1x,es14.7)') 'RVAP =',rvap
               write (unit=*,fmt='(a)') ' ============ Failed guess... ==========='
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLA =',tlcla,'FUNA =',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLZ =',tlclz,'FUNC =',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'DELTA =',delta,'FUNN =',funnow
               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'lcl_il8','therm_lib8.f90')
            end if
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      We have the guesses, solve the regula falsi method.                        !
         !---------------------------------------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            !----- Update guess and function evaluation. ----------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)
            pvap   = eslif8(tlcl,frozen)
            funnow = tlcl * (es00/pvap)**rocp8 - thil
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
            converged = abs(tlcl-tlcla) < toler8*tlcl .and.  abs(tlcl-tlclz) < toler8*tlcl
            if (funnow == 0.d0 .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.d0) then 
               tlclz = tlcl
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 5.d-1
               !---------------------------------------------------------------------------!


               !----- We have just updated zside, sett zside to true. ---------------------!
               zside = .true.
               !---------------------------------------------------------------------------!
            else
               tlcla = tlcl
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 5.d-1
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
         pvap  = eslif8(tlcl,frozen)
         plcl  = (ep8 + rvap) * pvap / rvap
         dzlcl = max(cpog8*(temp-tlcl),0.d0)
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
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Pressure        [   hPa] =',1.d-2*pres
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Temperature     [    °C] =',temp-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rtot            [  g/kg] =',1.d3*rtot
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rvap            [  g/kg] =',1.d3*rvap
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome.'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcla           [    °C] =',tlcla-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlclz           [    °C] =',tlclz-t008
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [     K] =',funnow
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [     K] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [     K] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [  ----] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler8          [  ----] =',toler8
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                                            ,abs(tlclz-tlcla)/tlclz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcl            [    °C] =',tlcl
         call abort_run  ('TLCL didn''t converge, gave up!','lcl_il8','therm_lib8.f90')
      end if
      return
   end subroutine lcl_il8
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
   subroutine thil2tqall8(thil,exner,pres,rtot,rliq,rice,temp,rvap,rsat)
      use rconstants , only : cpdry8    & ! intent(in)
                            , cpdryi8   & ! intent(in)
                            , t008      & ! intent(in)
                            , toodry8   & ! intent(in)
                            , t3ple8    & ! intent(in)
                            , ttripoli8 ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in)    :: thil      ! Ice-liquid water potential temp.  [     K]
      real(kind=8), intent(in)    :: exner     ! Exner function                    [J/kg/K]
      real(kind=8), intent(in)    :: pres      ! Pressure                          [    Pa]
      real(kind=8), intent(in)    :: rtot      ! Total mixing ratio                [ kg/kg]
      real(kind=8), intent(out)   :: rliq      ! Liquid water mixing ratio         [ kg/kg]
      real(kind=8), intent(out)   :: rice      ! Ice mixing ratio                  [ kg/kg]
      real(kind=8), intent(inout) :: temp      ! Temperature                       [     K]
      real(kind=8), intent(out)   :: rvap      ! Water vapour mixing ratio         [ kg/kg]
      real(kind=8), intent(out)   :: rsat      ! Sat. water vapour mixing ratio    [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                :: tempa     ! Lower bound for regula falsi iteration
      real(kind=8)                :: tempz     ! Upper bound for regula falsi iteration
      real(kind=8)                :: t1stguess ! Book keeping temperature 1st guess 
      real(kind=8)                :: fun1st    ! Book keeping 1st guess function
      real(kind=8)                :: funa      ! Function evaluation at tempa
      real(kind=8)                :: funz      ! Function evaluation at tempz
      real(kind=8)                :: funnow    ! Function at this iteration.
      real(kind=8)                :: delta     ! Aux. var in case we need regula falsi.
      real(kind=8)                :: deriv     ! Derivative of this function.
      integer                     :: itn       ! Iteration counter
      integer                     :: itb       ! Iteration counter
      integer                     :: ii        ! Iteration counter
      logical                     :: converged ! Convergence handle
      logical                     :: zside     ! Aux. Flag, for two purposes:
                                               ! 1. Found a 2nd guess for regula falsi.
                                               ! 2. I retained the "zside" (T/F)
      !------------------------------------------------------------------------------------!

      t1stguess = temp

      !------------------------------------------------------------------------------------!
      !      First check: try to find temperature assuming sub-saturation and check if     !
      ! this is the case. If it is, then there is no need to go through the iterative      !
      ! loop.                                                                              !
      !------------------------------------------------------------------------------------!
      tempz  = cpdryi8 * thil * exner
      rsat   = max(toodry8,rslif8(pres,tempz))
      if (tempz >= t3ple8) then
         rliq = max(0.d0,rtot-rsat)
         rice = 0.d0
      else
         rice = max(0.d0,rtot-rsat)
         rliq = 0.d0
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
      rsat   = max(toodry8,rslif8(pres,tempz))
      if (tempz >= t3ple8) then
         rliq = max(0.d0,rtot-rsat)
         rice = 0.d0
      else
         rice = max(0.d0,rtot-rsat)
         rliq = 0.d0
      end if
      rvap = rtot-rliq-rice


      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=46,fmt='(a)') '--------------------------------------------------------'
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !------------------------------------------------------------------------------------!
      !     Find the function. We are seeking a temperature which is associated with the   !
      ! theta_il we provided. Thus, the function is simply the difference between the      !
      ! theta_il associated with our guess and the actual theta_il.                        !
      !------------------------------------------------------------------------------------!
      funnow = theta_iceliq8(exner,tempz,rliq,rice)
      !----- Updating the derivative. -----------------------------------------------------!
      deriv  = dthetail_dt8(.false.,funnow,exner,pres,tempz,rliq,rice)
      funnow = funnow - thil
      fun1st = funnow
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x),a,1x,es11.4,1x)')                  &
      !   'NEWTON: it=',0,'temp=',tempz-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq          &
      !  ,'rice=',1000.*rice,'rvap=',1000.*rvap,'fun=',funnow,'deriv=',deriv
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
         if (abs(deriv) < toler8) exit newloop

         tempz = tempa - funnow / deriv

         !----- Finding the mixing ratios associated with this guess ----------------------!
         rsat  = max(toodry8,rslif8(pres,tempz))
         if (tempz >= t3ple8) then
            rliq = max(0.d0,rtot-rsat)
            rice = 0.d0
         else
            rice = max(0.d0,rtot-rsat)
            rliq = 0.d0
         end if
         rvap = rtot-rliq-rice

         !----- Updating the function -----------------------------------------------------!
         funnow = theta_iceliq8(exner,tempz,rliq,rice)
         !----- Updating the derivative. --------------------------------------------------!
         deriv  = dthetail_dt8(.false.,funnow,exner,pres,tempz,rliq,rice)
         funnow = funnow - thil

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x),a,1x,es11.4,1x)')               &
         !   'NEWTON: it=',itn,'temp=',tempz-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq     &
         !  ,'rice=',1000.*rice,'rvap=',1000.*rvap,'fun=',funnow,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         
         converged = abs(tempa-tempz) < toler8*tempz
         !---------------------------------------------------------------------------------!
         !   Convergence. The temperature will be the mid-point between tempa and tempz.   !
         ! Fix the mixing ratios and return. But first check for converged due to luck. If !
         ! the guess gives a root, then that's it. It looks unlikely, but it actually      !
         ! happens sometimes and if not checked it becomes a singularity.                  !
         !---------------------------------------------------------------------------------!
         if (funnow == 0.d0) then
            temp = tempz
            converged = .true.
            exit newloop
         elseif (converged) then
            temp = 5.d-1 * (tempa+tempz)
            rsat  = max(toodry8,rslif8(pres,temp))
            if (temp >= t3ple8) then
               rliq = max(0.d0,rtot-rsat)
               rice = 0.d0
            else
               rice = max(0.d0,rtot-rsat)
               rliq = 0.d0
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
         if (funa*funnow < 0.d0) then
            funz  = funnow
            zside = .true.
         !----- Otherwise, checking whether the 1st guess had opposite sign. --------------!
         elseif (funa*fun1st < 0.d0 ) then
            funz  = fun1st
            zside = .true.
         !---------------------------------------------------------------------------------!
         !    Looking for a guess. Extrapolate funa linearly, trying to get the -funa.  We !
         ! don't need it to be funa, just with the opposite sign.  If that's not enough,   !
         ! we keep going further... Force the guesses to be at least 1K apart              !
         !---------------------------------------------------------------------------------!
         else
            if (abs(funnow-funa) < 1.d2 * toler8 * tempa) then
               delta = 1.d2 * toler8 * tempa
            else
               delta = max(abs(funa)*abs((tempz-tempa)/(funnow-funa)),1.d2*toler8*tempa)
            end if
            tempz = tempa + delta
            funz  = funa
            !----- Just to enter at least once. The 1st time tempz=tempa-2*delta ----------!
            zside = .false. 
            zgssloop: do itb=1,maxfpo
                tempz = tempa + dble((-1)**itb * (itb+3)/2) * delta
                rsat   = max(toodry8,rslif8(pres,tempz))
                if (tempz >= t3ple8) then
                   rliq = max(0.d0,rtot-rsat)
                   rice = 0.d0
                else
                   rice = max(0.d0,rtot-rsat)
                   rliq = 0.d0
                end if
                rvap = rtot-rliq-rice
                funz = theta_iceliq8(exner,tempz,rliq,rice) - thil
                zside = funa*funz < 0.d0
                if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
               write (unit=*,fmt='(a)')           ' THIL2TQALL: NO SECOND GUESS FOR YOU!'
               write (unit=*,fmt='(a)')           ' '
               write (unit=*,fmt='(a)')           ' -> Input: '
               write (unit=*,fmt='(a,1x,f12.5)')  '    THETA_IL [     K]:',thil
               write (unit=*,fmt='(a,1x,f12.5)')  '    PRESS    [   hPa]:',1.d-2*pres
               write (unit=*,fmt='(a,1x,f12.5)')  '    EXNER    [J/kg/K]:',exner
               write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1.d3*rtot
               write (unit=*,fmt='(a,1x,f12.5)')  '    T1ST     [  degC]:',t1stguess-t008
               write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPA    [  degC]:',tempa-t008
               write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPZ    [  degC]:',tempz-t008
               write (unit=*,fmt='(a,1x,f12.5)')  '    FUNNOW   [     K]:',funnow
               write (unit=*,fmt='(a,1x,f12.5)')  '    FUNA     [     K]:',funa
               write (unit=*,fmt='(a,1x,f12.5)')  '    FUNZ     [     K]:',funz
               write (unit=*,fmt='(a,1x,f12.5)')  '    DELTA    [     K]:',delta
               write (unit=*,fmt='(60a1)')        ('-',ii=1,60)

               call abort_run  ('Failed finding the second guess for regula falsi'         &
                               ,'thil2tqall8','therm_lib8.f90')
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
               temp = 5.d-1*(tempa+tempz)
            end if
            !----- Distributing vapour into the three phases ------------------------------!
            rsat   = max(toodry8,rslif8(pres,temp))
            rvap   = min(rtot,rsat)
            if (temp >= t3ple8) then
               rliq = max(0.d0,rtot-rsat)
               rice = 0.d0
            else
               rliq = 0.d0
               rice = max(0.d0,rtot-rsat)
            end if
            !----- Updating function ------------------------------------------------------!
            funnow = theta_iceliq8(exner,temp,rliq,rice) - thil

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=46,fmt='(a,1x,i5,1x,10(a,1x,f11.4,1x))')                          &
            !   'REGFAL: it=',itb,'temp=',temp-t00,'tempa=',tempa-t00,'tempz=',tempz-t00   &
            !  ,'rsat=',1000.*rsat,'rliq=',1000.*rliq,'rice=',1000.*rice                   &
            !  ,'rvap=',1000.*rvap,'fun=',funnow,'funa=',funa,'funz=',funz
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Checking for convergence or lucky guess. If it did, return, we found the  !
            ! solution. Otherwise, constrain the guesses.                                  !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa) < toler8*temp .and. abs(temp-tempz) < toler8*temp 
            if (funnow == 0.d0 .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.d0) then 
               tempz = temp
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa = funa * 5.d-1
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            elseif (funnow*funz < 0.d0) then
               tempa = temp
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 5.d-1
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
         if (abs(temp-t3ple8) < toler8*temp) then
            !----- Find rliq. -------------------------------------------------------------!
            rliq   = ( alvi8(temp) * (rtot - rsat)                                         &
                     + cpdry8 * temp * log ( cpdryi8 * exner * thil / temp ) )             &
                   / ( alvi8(temp) - alvl8(temp) )
            !------------------------------------------------------------------------------!

            !----- Correct rliq so it is bounded, and then find rice as the remainder. ----!
            rliq   = max( 0.d0, min( rtot - rsat,  rliq ) )
            rice   = max( 0.d0, rtot-rsat-rliq)
            !------------------------------------------------------------------------------!


            !----- Function evaluation. ---------------------------------------------------!
            funnow = theta_iceliq8(exner,temp,rliq,rice) - thil
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
         itb=itb+1
      end if
      !------------------------------------------------------------------------------------!

      if (.not. converged) then
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
         write (unit=*,fmt='(a)')           ' THIL2TQALL8 failed!'
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Input: '
         write (unit=*,fmt='(a,1x,f12.5)')  '    THETA_IL [     K]:',thil
         write (unit=*,fmt='(a,1x,f12.5)')  '    EXNER    [J/kg/K]:',exner
         write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1.d3*rtot
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(a)')           ' -> Output: '
         write (unit=*,fmt='(a,1x,i12)')    '    ITERATIONS       :',itb
         write (unit=*,fmt='(a,1x,f12.5)')  '    TEMP     [    °C]:',temp-t008
         write (unit=*,fmt='(a,1x,f12.5)')  '    RVAP     [  g/kg]:',1.d3*rvap
         write (unit=*,fmt='(a,1x,f12.5)')  '    RLIQ     [  g/kg]:',1.d3*rliq
         write (unit=*,fmt='(a,1x,f12.5)')  '    RICE     [  g/kg]:',1.d3*rice
         write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPA    [    °C]:',tempa-t008
         write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPZ    [    °C]:',tempz-t008
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNA     [     K]:',funnow
         write (unit=*,fmt='(a,1x,es12.5)') '    FUNZ     [     K]:',funnow
         write (unit=*,fmt='(a,1x,es12.5)') '    DERIV    [   ---]:',deriv
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_A    [   ---]:',abs(temp-tempa)/temp
         write (unit=*,fmt='(a,1x,es12.5)') '    ERR_Z    [   ---]:',abs(temp-tempz)/temp
         write (unit=*,fmt='(a)')           ' '
         write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
         call abort_run  ('Failed finding equilibrium, I gave up!','thil2tqall8'           &
                         ,'therm_lib8.f90')
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x))')                                 &
      !   'ANSWER: it=',itb,'funf=',funnow,'temp=',temp-t00                                &
      !  ,'rsat=',1000.*rsat,'rliq=',1000.*rliq,'rice=',1000.*rice,'rvap=',1000.*rvap
      !write (unit=46,fmt='(a)') '----------------------------------------------------------'
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      return
   end subroutine thil2tqall8
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
   subroutine thil2tqliq8(thil,exner,pres,rtot,rliq,temp,rvap,rsat)
      use rconstants , only : cpdryi8   & ! intent(in)
                            , toodry8   ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in)    :: thil      ! Ice-liquid water potential temp.  [     K]
      real(kind=8), intent(in)    :: exner     ! Exner function                    [J/kg/K]
      real(kind=8), intent(in)    :: pres      ! Pressure                          [    Pa]
      real(kind=8), intent(in)    :: rtot      ! Total mixing ratio                [ kg/kg]
      real(kind=8), intent(out)   :: rliq      ! Liquid water mixing ratio         [ kg/kg]
      real(kind=8), intent(inout) :: temp      ! Temperature                       [     K]
      real(kind=8), intent(out)   :: rvap      ! Water vapour mixing ratio         [ kg/kg]
      real(kind=8), intent(out)   :: rsat      ! Sat. water vapour mixing ratio    [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                :: tempa     ! Lower bound for regula falsi iteration
      real(kind=8)                :: tempz     ! Upper bound for regula falsi iteration
      real(kind=8)                :: t1stguess ! Book keeping temperature 1st guess 
      real(kind=8)                :: fun1st    ! Book keeping 1st guess function
      real(kind=8)                :: funa      ! Function evaluation at tempa
      real(kind=8)                :: funz      ! Function evaluation at tempz
      real(kind=8)                :: funnow    ! Function at this iteration.
      real(kind=8)                :: delta     ! Aux. var in case we need regula falsi.
      real(kind=8)                :: deriv     ! Derivative of this function.
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
      tempz = cpdryi8 * thil * exner
      rsat  = max(toodry8,rslf8(pres,tempz))
      rliq  = max(0.d0,rtot-rsat)
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
      rsat  = max(toodry8,rslf8(pres,tempz))
      rliq  = max(0.d0,rtot-rsat)
      rvap  = rtot-rliq
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Finding the function. We are seeking a temperature which is associated with    !
      ! the theta_il we provided. Thus, the function is simply the difference between the  !
      ! theta_il associated with our guess and the actual theta_il.                        !
      !------------------------------------------------------------------------------------!
      funnow = theta_iceliq8(exner,tempz,rliq,0.d0) ! Finding thil from our guess
      deriv  = dthetail_dt8(.false.,funnow,exner,pres,tempz,rliq)
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
         if (abs(deriv) < toler8) exit newloop

         tempz = tempa - funnow / deriv

         !----- Finding the mixing ratios associated with this guess ----------------------!
         rsat  = max(toodry8,rslf8(pres,tempz))
         rliq = max(0.d0,rtot-rsat)
         rvap = rtot-rliq

         !----- Updating the function -----------------------------------------------------!
         funnow = theta_iceliq8(exner,tempz,rliq,0.d0)
         !----- Updating the derivative. --------------------------------------------------!
         deriv = dthetail_dt8(.false.,funnow,exner,pres,tempz,rliq)
         funnow = funnow - thil

         converged = abs(tempa-tempz) < toler8*tempz
         !---------------------------------------------------------------------------------!
         !   Convergence. The temperature will be the mid-point between tempa and tempz.   !
         ! Fix the mixing ratios and return. But first check for converged due to luck. If !
         ! the guess gives a root, then that's it. It looks unlikely, but it actually      !
         ! happens sometimes and if not checked it becomes a singularity.                  !
         !---------------------------------------------------------------------------------!
         if (funnow == 0.d0) then
            temp = tempz
            converged = .true.
            exit newloop
         elseif (converged) then
            temp = 5.d-1 * (tempa+tempz)
            rsat  = max(toodry8,rslf8(pres,temp))
            rliq = max(0.d0,rtot-rsat)
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
         if (funa*funnow < 0.d0) then
            funz = funnow
            zside = .true.
         !---------------------------------------------------------------------------------!
         !    Looking for a guess. Extrapolate funa linearly, trying to get the -funa.  We !
         ! don't need it to be funa, just with the opposite sign.  If that's not enough,   !
         ! we keep going further... Force the guesses to be at least 1K apart              !
         !---------------------------------------------------------------------------------!
         else
            if (abs(funnow-funa) < toler8*tempa) then
               delta = 1.d2*toler8*tempa
            else
               delta = max(abs(funa*(tempz-tempa)/(funnow-funa)),1.d2*toler8*tempa)
            end if
            tempz = tempa + delta
            funz  = funa
            !----- Just to enter at least once. The 1st time tempz=tempa-2*delta ----------!
            zside = .false. 
            zgssloop: do itb=1,maxfpo
               tempz = tempz + dble((-1)**itb * (itb+3)/2) * delta
               rsat  = max(toodry8,rslf8(pres,tempz))
               rliq  = max(0.d0,rtot-rsat)
               rvap  = rtot-rliq
               funz  = theta_iceliq8(exner,tempz,rliq,0.d0) - thil
               zside = funa*funz < 0.d0
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
            rsat   = max(toodry8,rslf8(pres,temp))
            rvap   = min(rtot,rsat)
            rliq   = max(0.d0,rtot-rsat)
            !----- Updating function ------------------------------------------------------!
            funnow = theta_iceliq8(exner,tempz,rliq,0.d0) - thil

            !------------------------------------------------------------------------------!
            !    Checking for convergence or lucky guess. If it did, return, we found the  !
            ! solution. Otherwise, constrain the guesses.                                  !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa)< toler8*temp  .and. abs(temp-tempz) < toler8*temp
            if (funnow == 0.d0 .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.d0) then 
               tempz = temp
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa = funa * 5.d-1
               !----- We just updated zside, setting zside to true. -----------------------!
               zside = .true.
            else
               tempa = temp
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 5.d-1
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
            !------------------------------------------------------------------------------!
         end do fpoloop

      end if

      if (.not. converged) call abort_run  ('Failed finding equilibrium, I gave up!'       &
                                           ,'thil2tqliq8','therm_lib8.f90')
      return
   end subroutine thil2tqliq8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! intensive internal energy [J/kg].                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine uint2tl8(uint,temp,fliq)
      use rconstants , only : cliqi8          & ! intent(in)
                            , cicei8          & ! intent(in)
                            , allii8          & ! intent(in)
                            , t3ple8          & ! intent(in)
                            , uiicet38        & ! intent(in)
                            , uiliqt38        & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: uint     ! Internal energy                     [   J/kg]
      real(kind=8), intent(out) :: temp     ! Temperature                         [      K]
      real(kind=8), intent(out) :: fliq     ! Liquid Fraction (0-1)               [    ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compare the internal energy with the reference values to decide which phase    !
      ! the water is.                                                                      !
      !------------------------------------------------------------------------------------!
      if (uint <= uiicet38) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.d0
         temp    = uint * cicei8
         !---------------------------------------------------------------------------------!
      elseif (uint >= uiliqt38) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.d0
         temp    = uint * cliqi8 + tsupercool_liq8
         !---------------------------------------------------------------------------------!
      else
         !----- Changing phase, it must be at freezing point ------------------------------!
         fliq = (uint - uiicet38) * allii8
         temp    = t3ple8
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uint2tl8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from         !
   ! extensive internal energy (J/m² or J/m³), water mass (kg/m² or kg/m³), and heat       !
   ! capacity (J/m²/K or J/m³/K).                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine uextcm2tl8(uext,wmass,dryhcap,temp,fliq)
      use rconstants , only : cliqi8          & ! intent(in)
                            , cliq8           & ! intent(in)
                            , cicei8          & ! intent(in)
                            , cice8           & ! intent(in)
                            , allii8          & ! intent(in)
                            , alli8           & ! intent(in)
                            , t3ple8          & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: uext    ! Extensive internal energy [  J/m²] or [  J/m³]
      real(kind=8), intent(in)  :: wmass   ! Water mass                [ kg/m²] or [ kg/m³]
      real(kind=8), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m²/K] or [J/m³/K]
      real(kind=8), intent(out) :: temp    ! Temperature                           [     K]
      real(kind=8), intent(out) :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=8)              :: uefroz  ! qw of ice at triple pt.   [  J/m²] or [  J/m³] 
      real(kind=8)              :: uemelt  ! qw of liq. at triple pt.  [  J/m²] or [  J/m³]
      !------------------------------------------------------------------------------------!



      !----- Convert melting heat to J/m² or J/m³ -----------------------------------------!
      uefroz = (dryhcap + wmass * cice8) * t3ple8
      uemelt = uefroz   + wmass * alli8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    This is analogous to the uint2tl8 computation, we should analyse the magnitude  !
      ! of the internal energy to choose between liquid, ice, or both by comparing with    !
      ! the known boundaries.                                                              !
      !------------------------------------------------------------------------------------!
      if (uext < uefroz) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.d0
         temp = uext  / (cice8 * wmass + dryhcap)
         !---------------------------------------------------------------------------------!
      elseif (uext > uemelt) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.d0
         temp = (uext + wmass * cliq8 * tsupercool_liq8) / (dryhcap + wmass * cliq8)
         !---------------------------------------------------------------------------------!
      elseif (uefroz == uemelt) then
         !---------------------------------------------------------------------------------!
         !    We are at the freezing point.  If water mass is so tiny that the internal    !
         ! energy of frozen and melted states are the same given the machine precision,    !
         ! then we assume that water content is negligible and we impose 50% frozen for    !
         ! simplicity.                                                                     !
         !---------------------------------------------------------------------------------!
         fliq = 5.d-1
         temp = t3ple8
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    Changing phase, it must be at freezing point.  The max and min are here just !
         ! to avoid tiny deviations beyond 0. and 1. due to floating point arithmetics.    !
         !---------------------------------------------------------------------------------!
         fliq = min(1.d0,max(0.d0,(uext - uefroz) * allii8 / wmass))
         temp = t3ple8
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uextcm2tl8
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib8
!==========================================================================================!
!==========================================================================================!

