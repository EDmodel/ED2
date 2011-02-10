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
      logical                       :: brrr_cold

      !----- Checking which saturation (liquid or ice) I should use here ------------------!
      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple
      else 
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      !----- Finding the saturation vapour pressure ---------------------------------------!
      if (brrr_cold) then
         esz = esif(temp)
      else
         esz = eslf(temp)
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
      use consts_coms, only : rh2o
      implicit none
      real, intent(in) :: temp
      real             :: eequ
      eequ = eslf(temp)
      rhovsl = eequ / (rh2o * temp)
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
      use consts_coms, only : rh2o
      implicit none
      real, intent(in) :: temp
      real             :: eequ
      eequ = esif(temp)
      rhovsi = eequ / (rh2o * temp)
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
      use consts_coms, only : rh2o
      implicit none
      real, intent(in)              :: temp
      logical, intent(in), optional :: useice
      real                          :: eequ

      if (present(useice)) then
         eequ = eslif(temp,useice)
      else
         eequ = eslif(temp)
      end if

      rhovsil = eequ / (rh2o * temp)

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
      real                          :: desdt
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
      use consts_coms, only : rh2o
      implicit none
      real, intent(in) :: temp
      real             :: es,desdt

      es    = eslf(temp)
      desdt = eslfp(temp)
      rhovslp = (desdt-es/temp) / (rh2o * temp)

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
      use consts_coms, only : rh2o
      implicit none
      real, intent(in) :: temp
      real             :: es,desdt

      es    = esif(temp)
      desdt = esifp(temp)
      rhovsip = (desdt-es/temp) / (rh2o * temp)

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
   !     This function calculates the temperature from the liquid sat. vapour pressure.    !
   ! This is truly the inverse of eslf, which is done iteratively since it's not a simple  !
   ! function to invert. It uses Newton's method, which should take care of most cases. In !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real function tslf(pvap)

       implicit none
      !----- Argument ---------------------------------------------------------------------!
      real, intent(in)   :: pvap       ! Saturation vapour pressure                [    Pa]
      !----- Local variables for iterative method -----------------------------------------!
      real               :: deriv      ! Function derivative                       [    Pa]
      real               :: fun        ! Function for which we seek a root.        [    Pa]
      real               :: funa       ! Smallest  guess function                  [    Pa]
      real               :: funz       ! Largest   guess function                  [    Pa]
      real               :: tempa      ! Smallest guess (or previous guess)        [    Pa]
      real               :: tempz      ! Largest   guess (or new guess in Newton)  [    Pa]
      real               :: delta      ! Aux. var --- 2nd guess for bisection      [      ]
      integer            :: itn,itb    ! Iteration counter                         [   ---]
      logical            :: converged  ! Convergence handle                        [   ---]
      logical            :: zside      ! Flag to check for one-sided approach...   [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, using Bolton (1980) equation 11, giving es in Pa and T in K -----!
      tempa = (29.65 * log(pvap) - 5016.78)/(log(pvap)-24.0854)
      funa  = eslf(tempa) - pvap
      deriv = eslfp(tempa)
      !----- Copying just in case it fails at the first iteration -------------------------!
      tempz = tempa
      fun   = funa
      
      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, go with bisection -----!
         !----- Copying the previous guess ------------------------------------------------!
         tempa = tempz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         tempz = tempa - fun/deriv
         fun   = eslf(tempz) - pvap
         deriv = eslfp(tempz)
         
         converged = abs(tempa-tempz) < toler * tempz
         if (converged) then
            tslf = 0.5*(tempa+tempz)
            return
         elseif (fun ==0) then !Converged by luck!
            tslf = tempz
            return
         end if
      end do newloop

      !------------------------------------------------------------------------------------!
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < 100.*toler*tempa) then
            delta = 100.*toler*tempa
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            funz  = eslf(tempz) - pvap
            zside = funa*funz < 0.
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                          ,'tslf','therm_lib.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tslf =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tslf-tempa) < toler * tslf
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         fun       =  eslf(tslf) - pvap

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0. ) then
            tempz = tslf
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 0.5
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = tslf
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz=funz * 0.5
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         call fatal_error('Temperature didn''t converge, giving up!!!'                     &
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
   real function tsif(pvap)

       implicit none
      !----- Argument ---------------------------------------------------------------------!
      real, intent(in)   :: pvap       ! Saturation vapour pressure                [    Pa]
      !----- Local variables for iterative method -----------------------------------------!
      real               :: deriv      ! Function derivative                       [    Pa]
      real               :: fun        ! Function for which we seek a root.        [    Pa]
      real               :: funa       ! Smallest  guess function                  [    Pa]
      real               :: funz       ! Largest   guess function                  [    Pa]
      real               :: tempa      ! Smallest guess (or previous guess)        [    Pa]
      real               :: tempz      ! Largest   guess (or new guess in Newton)  [    Pa]
      real               :: delta      ! Aux. var --- 2nd guess for bisection      [      ]
      integer            :: itn,itb    ! Iteration counter                         [   ---]
      logical            :: converged  ! Convergence handle                        [   ---]
      logical            :: zside      ! Flag to check for one-sided approach...   [   ---]
      !------------------------------------------------------------------------------------!

      !----- First Guess, using Murphy-Koop (2005), equation 8. ---------------------------!
      tempa = (1.814625 * log(pvap) +6190.134)/(29.120 - log(pvap))
      funa  = esif(tempa) - pvap
      deriv = esifp(tempa)
      !----- Copying just in case it fails at the first iteration -------------------------!
      tempz = tempa
      fun   = funa
      
      !----- Enter Newton's method loop: --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, go with bisection -----!
         !----- Copying the previous guess ------------------------------------------------!
         tempa = tempz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         tempz = tempa - fun/deriv
         fun   = esif(tempz) - pvap
         deriv = esifp(tempz)
         
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
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < 100.*toler*tempa) then
            delta = 100.*toler*delta
         else
            delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
         end if
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            funz  = esif(tempz) - pvap
            zside = funa*funz < 0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(1(a,1x,es14.7))') 'pvap =',pvap
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                          ,'tsif','therm_lib.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         tsif =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(tsif-tempa) < toler * tsif
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         fun       =  esif(tsif) - pvap

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0. ) then
            tempz = tsif
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 0.5
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = tsif
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz=funz * 0.5
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         call fatal_error('Temperature didn''t converge, giving up!!!'                     &
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
   real function tslif(pvap,useice)
      use consts_coms, only: es3ple,alvl,alvi

      implicit none
      real   , intent(in)           :: pvap
      logical, intent(in), optional :: useice
      logical                       :: brrr_cold
      
      !------------------------------------------------------------------------------------!
      !    Since pvap is a function of temperature only, we can check the triple point     !
      ! from the saturation at the triple point, like what we would do for temperature.    !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         brrr_cold = useice  .and. pvap < es3ple
      else 
         brrr_cold = bulk_on .and. pvap < es3ple
      end if


      if (brrr_cold) then
         tslif = tsif(pvap)
      else
         tslif = tslf(pvap)
      end if

      return
   end function tslif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the dew point temperature given the pressure and vapour    !
   ! mixing ratio. THIS IS DEWPOINT ONLY, WHICH MEANS THAT IT WILL IGNORE ICE EFFECT. For  !
   ! a full, triple-point dependent routine use DEWFROSTPOINT                              !
   !---------------------------------------------------------------------------------------!
   real function dewpoint(pres,rsat)
      use consts_coms, only: ep,toodry

      implicit none
      real, intent(in) :: pres, rsat
      real             :: rsatoff, pvsat
      
      rsatoff  = max(toodry,rsat)
      pvsat    = pres*rsatoff / (ep + rsatoff)
      dewpoint = tslf(pvsat)

      return
   end function dewpoint
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the frost point temperature given the pressure and vapour  !
   ! mixing ratio. THIS IS FROSTPOINT ONLY, WHICH MEANS THAT IT WILL IGNORE LIQUID EFFECT. !
   ! For a full, triple-point dependent routine use DEWFROSTPOINT                          !
   !---------------------------------------------------------------------------------------!
   real function frostpoint(pres,rsat)
      use consts_coms, only: ep,toodry

      implicit none
      real, intent(in) :: pres, rsat
      real             :: rsatoff, pvsat
      
      rsatoff    = max(toodry,rsat)
      pvsat      = pres*rsatoff / (ep + rsatoff)
      frostpoint = tsif(pvsat)

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
   real function dewfrostpoint(pres,rsat,useice)
      use consts_coms, only: ep,toodry

      implicit none
      real   , intent(in)           :: pres, rsat
      logical, intent(in), optional :: useice
      real                          :: rsatoff, pvsat

      rsatoff       = max(toodry,rsat)
      pvsat         = pres*rsatoff / (ep + rsatoff)
      if (present(useice)) then
         dewfrostpoint = tslif(pvsat,useice)
      else
         dewfrostpoint = tslif(pvsat)
      end if
      return
   end function dewfrostpoint
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based on the pressure [Pa], tem-   !
   ! perature [K] and relative humidity [fraction]. IT ALWAYS ASSUME THAT RELATIVE HUMI-   !
   ! DITY IS WITH RESPECT TO THE LIQUID PHASE.  ptrh2rvapil checks which one to use        !
   ! depending on whether temperature is more or less than the triple point.               !
   !---------------------------------------------------------------------------------------!
   real function ptrh2rvapl(relh,pres,temp)
      use consts_coms, only: ep,toodry

      implicit none
      real   , intent(in)           :: relh, pres, temp
      real                          :: rsath, relhh

      rsath = max(toodry,rslf(pres,temp))
      relhh = min(1.,max(0.,relh))

      if (newthermo) then
         !----- Considering that Rel.Hum. is based on vapour pressure ---------------------!
         ptrh2rvapl = max(toodry,ep * relhh * rsath / (ep + (1.-relhh)*rsath))
      else
         !----- Original RAMS way to compute mixing ratio ---------------------------------!
         ptrh2rvapl = max(toodry,relhh*rsath)
      end if

      return
   end function ptrh2rvapl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based on the pressure [Pa], tem-   !
   ! perature [K] and relative humidity [fraction]. IT ALWAYS ASSUME THAT RELATIVE HUMI-   !
   ! DITY IS WITH RESPECT TO THE ICE PHASE. ptrh2rvapil checks which one to use depending  !
   ! on whether temperature is more or less than the triple point.                         !
   !---------------------------------------------------------------------------------------!
   real function ptrh2rvapi(relh,pres,temp)
      use consts_coms, only: ep,toodry

      implicit none
      real   , intent(in)           :: relh, pres, temp
      real                          :: rsath, relhh

      rsath = max(toodry,rsif(pres,temp))
      relhh = min(1.,max(0.,relh))

      if (newthermo) then
         !----- Considering that Rel.Hum. is based on vapour pressure ---------------------!
         ptrh2rvapi = max(toodry,ep * relhh * rsath / (ep + (1.-relhh)*rsath))
      else
         !----- Original RAMS way to compute mixing ratio ---------------------------------!
         ptrh2rvapi = max(toodry,relhh*rsath)
      end if

      return
   end function ptrh2rvapi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based on the pressure [Pa], tem-   !
   ! perature [K] and relative humidity [fraction]. It will check the temperature to       !
   ! decide between ice or liquid saturation and whether ice should be considered.         !
   !---------------------------------------------------------------------------------------!
   real function ptrh2rvapil(relh,pres,temp,useice)
      use consts_coms, only: ep,toodry,t3ple

      implicit none
      real   , intent(in)           :: relh, pres, temp
      logical, intent(in), optional :: useice
      real                          :: rsath, relhh
      logical                       :: brrr_cold

      !----- Checking whether I use the user or the default check for ice saturation. -----!
      if (present(useice)) then
         brrr_cold = useice .and. temp < t3ple
      else
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      if (brrr_cold) then
         rsath = max(toodry,rsif(pres,temp))
      else
         rsath = max(toodry,rslf(pres,temp))
      end if

      relhh = min(1.,max(0.,relh))
      
      ptrh2rvapil = max(toodry,ep * relhh * rsath / (ep + (1.-relhh)*rsath))

      return
   end function ptrh2rvapil
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
   real function rehul(pres,temp,rvap)
      use consts_coms, only: ep,toodry
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: pres     ! Air pressure                     [    Pa]
      real   , intent(in)           :: temp     ! Temperature                      [     K]
      real   , intent(in)           :: rvap     ! Vapour mixing ratio              [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real                          :: rvapsat  ! Saturation mixing ratio          [ kg/kg]
      !------------------------------------------------------------------------------------!

      rvapsat = max(toodry,rslf(pres,temp))
      if (newthermo) then
         !----- This is based on relative humidity being defined with vapour pressure -----!
         rehul = max(0.,rvap*(ep+rvapsat)/(rvapsat*(ep+rvap)))
      else
         !----- Original formula used by RAMS ---------------------------------------------!
         rehul = max(0.,rvap/rvapsat)
      end if
      return
   end function rehul
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
   real function rehui(pres,temp,rvap)
      use consts_coms, only: ep,toodry
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: pres     ! Air pressure                     [    Pa]
      real   , intent(in)           :: temp     ! Temperature                      [     K]
      real   , intent(in)           :: rvap     ! Vapour mixing ratio              [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real                          :: rvapsat  ! Saturation mixing ratio          [ kg/kg]
      !------------------------------------------------------------------------------------!

      rvapsat = max(toodry,rsif(pres,temp))
      if (newthermo) then
         !----- This is based on relative humidity being defined with vapour pressure -----!
         rehui = max(0.,rvap*(ep+rvapsat)/(rvapsat*(ep+rvap)))
      else
         !----- Original formula used by RAMS ---------------------------------------------!
         rehui = max(0.,rvap/rvapsat)
      end if
      return
   end function rehui
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
   real function rehuil(pres,temp,rvap,useice)
      use consts_coms, only: t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: pres      ! Air pressure                    [    Pa]
      real   , intent(in)           :: temp      ! Temperature                     [     K]
      real   , intent(in)           :: rvap      ! Vapour mixing ratio             [ kg/kg]
      logical, intent(in), optional :: useice    ! Should I consider ice?          [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real                          :: rvapsat   ! Saturation mixing ratio         [ kg/kg]
      logical                       :: brrr_cold ! I will use ice saturation now   [   T|F]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    Checking whether I should go with ice or liquid saturation.                     !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         brrr_cold = useice  .and. temp < t3ple
      else 
         brrr_cold = bulk_on .and. temp < t3ple
      end if

      if (brrr_cold) then
         rehuil = rehui(pres,temp,rvap)
      else
         rehuil = rehul(pres,temp,rvap)
      end if

      return
   end function rehuil
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
   real function tv2temp(tvir,rvap,rtot)
      use consts_coms, only: epi
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: tvir     ! Virtual temperature                  [    K]
      real, intent(in)           :: rvap     ! Vapour mixing ratio                  [kg/kg]
      real, intent(in), optional :: rtot     ! Total mixing ratio                   [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real                       :: rtothere ! Internal rtot, to deal with optional [kg/kg]

      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvap
      end if

      tv2temp = tvir * (1. + rtothere) / (1. + epi*rvap)

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
   real function virtt(temp,rvap,rtot)
      use consts_coms, only: epi
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: temp     ! Temperature                          [    K]
      real, intent(in)           :: rvap     ! Vapour mixing ratio                  [kg/kg]
      real, intent(in), optional :: rtot     ! Total mixing ratio                   [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real                       :: rtothere ! Internal rtot, to deal with optional [kg/kg]

      if (present(rtot)) then
        rtothere = rtot
      else
        rtothere = rvap
      end if

      virtt = temp * (1. + epi * rvap) / (1. + rtothere)

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
   real function idealdens(pres,temp,rvap,rtot)
      use consts_coms, only: rdry
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: pres     ! Pressure                             [   Pa]
      real, intent(in)           :: temp     ! Temperature                          [    K]
      real, intent(in)           :: rvap     ! Vapour mixing ratio                  [kg/kg]
      real, intent(in), optional :: rtot     ! Total mixing ratio                   [kg/kg]
      !----- Local variable ---------------------------------------------------------------!
      real                       :: tvir     ! Virtual temperature                  [    K]
      !------------------------------------------------------------------------------------!
      if (present(rtot)) then
        tvir = virtt(temp,rvap,rtot)
      else
        tvir = virtt(temp,rvap)
      end if

      idealdens = pres / (rdry * tvir)

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
   real function idealdenssh(pres,temp,qvpr,qtot)
      use consts_coms, only : rdry & ! intent(in)
                            , epi  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: pres ! Pressure                                 [   Pa]
      real, intent(in)           :: temp ! Temperature                              [    K]
      real, intent(in)           :: qvpr ! Vapour specific mass                     [kg/kg]
      real, intent(in), optional :: qtot ! Total water specific mass                [kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real                       :: qall ! Either qtot or qvpr...                   [kg/kg]
      !------------------------------------------------------------------------------------!

      if (present(qtot)) then
        qall = qtot
      else
        qall = qvpr
      end if

      idealdenssh = pres / (rdry * temp * (1. - qall + epi * qvpr))

      return
   end function idealdenssh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes reduces the pressure from the reference height to the      !
   ! canopy height by assuming hydrostatic equilibrium.                                    !
   !---------------------------------------------------------------------------------------!
   real function reducedpress(pres,thetaref,shvref,zref,thetacan,shvcan,zcan)
      use consts_coms, only : epim1    & ! intent(in)
                            , p00k     & ! intent(in)
                            , rocp     & ! intent(in)
                            , cpor     & ! intent(in)
                            , cp       & ! intent(in)
                            , grav     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: pres     ! Pressure                          [      Pa]
      real, intent(in)           :: thetaref ! Potential temperature             [       K]
      real, intent(in)           :: shvref   ! Vapour specific mass              [   kg/kg]
      real, intent(in)           :: zref     ! Height at reference level         [       m]
      real, intent(in)           :: thetacan ! Potential temperature             [       K]
      real, intent(in)           :: shvcan   ! Vapour specific mass              [   kg/kg]
      real, intent(in)           :: zcan     ! Height at canopy level            [       m]
      !------Local variables. -------------------------------------------------------------!
      real                       :: pinc     ! Pressure increment                [ Pa^R/cp]
      real                       :: thvbar   ! Average virtual pot. temperature  [       K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      First we compute the average virtual potential temperature between the canopy !
      ! top and the reference level.                                                       !
      !------------------------------------------------------------------------------------!
      thvbar = 0.5 * (thetaref * (1. + epim1 * shvref) + thetacan * (1. + epim1 * shvcan))

      !----- Then, we find the pressure gradient scale. -----------------------------------!
      pinc = grav * p00k * (zref - zcan) / (cp * thvbar)

      !----- And we can find the reduced pressure. ----------------------------------------!
      reducedpress = (pres**rocp + pinc ) ** cpor

      return
   end function reducedpress
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the enthalpy given the pressure, temperature, vapour       !
   ! specific humidity, and height.  Currently it doesn't compute mixed phase air, but     !
   ! adding it should be straight forward (finding the inverse is another story...).       !
   !---------------------------------------------------------------------------------------!
   real function ptqz2enthalpy(pres,temp,qvpr,zref)
      use consts_coms, only : ep       & ! intent(in)
                            , grav     & ! intent(in)
                            , t3ple    & ! intent(in)
                            , eta3ple  & ! intent(in)
                            , cimcp    & ! intent(in)
                            , clmcp    & ! intent(in)
                            , cp       & ! intent(in)
                            , alvi     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: pres  ! Pressure                               [    Pa]
      real, intent(in)           :: temp  ! Temperature                            [     K]
      real, intent(in)           :: qvpr  ! Vapour specific mass                   [ kg/kg]
      real, intent(in)           :: zref  ! Reference height                       [     m]
      !------Local variables. -------------------------------------------------------------!
      real                       :: tequ  ! Dew-frost temperature                  [     K]
      real                       :: pequ  ! Equlibrium vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!

      !----- First, we find the equilibrium vapour pressure and dew/frost point. ----------!
      pequ = pres * qvpr / (ep + (1. - ep) * qvpr)
      tequ = tslif(pequ)

      !------------------------------------------------------------------------------------!
      !     Then, based on dew/frost point, we compute the enthalpy. This accounts whether !
      ! we would have to dew or frost formation if the temperature dropped to the          !
      ! equilibrium point.  Notice that if supersaturation exists, this will still give a  !
      ! number that makes sense, similar to the internal energy of supercooled water.      !
      !------------------------------------------------------------------------------------!
      if (tequ <= t3ple) then
         ptqz2enthalpy = cp * temp + qvpr * (cimcp * tequ + alvi   ) + grav * zref
      else
         ptqz2enthalpy = cp * temp + qvpr * (clmcp * tequ + eta3ple) + grav * zref
      end if

      return
   end function ptqz2enthalpy
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature given the enthalpy, pressure, vapour       !
   ! specific humidity, and reference height.  Currently it doesn't compute mixed phase    !
   ! air, but adding it wouldn't be horribly hard, though it would require some root       !
   ! finding.                                                                              !
   !---------------------------------------------------------------------------------------!
   real function hpqz2temp(enthalpy,pres,qvpr,zref)
      use consts_coms, only : ep       & ! intent(in)
                            , grav     & ! intent(in)
                            , t3ple    & ! intent(in)
                            , eta3ple  & ! intent(in)
                            , cimcp    & ! intent(in)
                            , clmcp    & ! intent(in)
                            , cpi      & ! intent(in)
                            , alvi     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: enthalpy ! Enthalpy...                         [  J/kg]
      real, intent(in)           :: pres     ! Pressure                            [    Pa]
      real, intent(in)           :: qvpr     ! Vapour specific mass                [ kg/kg]
      real, intent(in)           :: zref     ! Reference height                    [     m]
      !------Local variables. -------------------------------------------------------------!
      real                       :: tequ  ! Dew-frost temperature                  [     K]
      real                       :: pequ  ! Equlibrium vapour pressure             [    Pa]
      !------------------------------------------------------------------------------------!

      !----- First, we find the equilibrium vapour pressure and dew/frost point. ----------!
      pequ = pres * qvpr / (ep + (1. - ep) * qvpr)
      tequ = tslif(pequ)

      !------------------------------------------------------------------------------------!
      !     Then, based on dew/frost point, we compute the temperature. This accounts      !
      ! whether we would have to dew or frost formation if the temperature dropped to the  !
      ! equilibrium point.  Notice that if supersaturation exists, this will still give a  !
      ! temperature that makes sense (but less than the dew/frost point), similar to the   !
      ! internal energy of supercooled water.                                              !
      !------------------------------------------------------------------------------------!
      if (tequ <= t3ple) then
         hpqz2temp = cpi * (enthalpy - qvpr * (cimcp * tequ + alvi   ) - grav * zref)
      else
         hpqz2temp = cpi * (enthalpy - qvpr * (clmcp * tequ + eta3ple) - grav * zref)
      end if

      return
   end function hpqz2temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the temperature given the potential temperature, density, and !
   ! specific humidity.  This comes from a combination of the definition of potential      !
   ! temperature and the ideal gas law, to eliminate pressure, when pressure is also       !
   ! unknown.                                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function thrhsh2temp(theta,dens,qvpr)
      use consts_coms, only : cpocv  & ! intent(in)
                            , p00i   & ! intent(in)
                            , rdry   & ! intent(in)
                            , epim1  & ! intent(in)
                            , rocv   ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: theta    ! Potential temperature                 [     K]
      real(kind=4), intent(in) :: dens     ! Density                               [    Pa]
      real(kind=4), intent(in) :: qvpr     ! Specific humidity                     [ kg/kg]
      !------------------------------------------------------------------------------------!

      thrhsh2temp = theta ** cpocv                                                         &
                  * (p00i * dens * rdry * (1. + epim1 * qvpr)) ** rocv

      return
   end function thrhsh2temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the ice liquid potential temperature given the Exner       !
   ! function [J/kg/K], temperature [K], and liquid and ice mixing ratios [kg/kg].         !
   !---------------------------------------------------------------------------------------!
   real function theta_iceliq(exner,temp,rliq,rice)
      use consts_coms, only: alvl, alvi, cp, ttripoli, htripoli, htripolii

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: exner   ! Exner function                       [J/kg/K]
      real, intent(in)           :: temp    ! Temperature                          [     K]
      real, intent(in)           :: rliq    ! Liquid mixing ratio                  [ kg/kg]
      real, intent(in)           :: rice    ! Ice mixing ratio                     [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real             :: hh    ! Enthalpy associated with sensible heat           [  J/kg]
      real             :: qq    ! Enthalpy associated with latent heat             [  J/kg]
      !------------------------------------------------------------------------------------!

      !----- Finding the enthalpies -------------------------------------------------------!
      hh = cp*temp
      qq  = alvl*rliq+alvi*rice
      
      if (newthermo) then
      
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli) then
            theta_iceliq = hh * exp(-qq/hh) / exner
         else
            theta_iceliq = hh * exp(-qq * htripolii) / exner
         end if
      else
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli) then
            theta_iceliq = hh * hh / (exner * ( hh + qq))
         else
            theta_iceliq = hh * htripoli / (exner * ( htripoli + qq))
         end if
      end if

      return
   end function theta_iceliq
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the liquid potential temperature derivative with respect   !
   ! to temperature, useful in iterative methods.                                          !
   !---------------------------------------------------------------------------------------!
   real function dthetail_dt(condconst,thil,exner,pres,temp,rliq,ricein)
      use consts_coms, only: alvl, alvi, cp, ttripoli,htripoli,htripolii,t3ple

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      logical, intent(in)           :: condconst  ! Condensation is constant?      [   T|F]
      real   , intent(in)           :: thil       ! Ice liquid pot. temperature    [     K]
      real   , intent(in)           :: exner      ! Exner function                 [J/kg/K]
      real   , intent(in)           :: pres       ! Pressure                       [    Pa]
      real   , intent(in)           :: temp       ! Temperature                    [     K]
      real   , intent(in)           :: rliq       ! Liquid mixing ratio            [ kg/kg]
      real   , intent(in), optional :: ricein     ! Ice mixing ratio               [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real                          :: rice       ! Ice mixing ratio or 0.         [ kg/kg]
      real                          :: ldrst      ! L × d(rs)/dT × T               [  J/kg]
      real                          :: hh         ! Sensible heat enthalpy         [  J/kg]
      real                          :: qq         ! Latent heat enthalpy           [  J/kg]
      logical                       :: thereisice ! Is ice present                 [   ---]
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !   Checking whether I should consider ice or not.                                   !
      !------------------------------------------------------------------------------------!
      thereisice = present(ricein)
      
      if (thereisice) then
         rice=ricein
      else
         rice=0.
      end if
      
      !----- No condensation, dthetail_dt is a constant -----------------------------------!
      if (rliq+rice == 0.) then
         dthetail_dt = thil/temp
         return
      else
         hh    = cp*temp                            !----- Sensible heat enthalpy
         qq    = alvl*rliq+alvi*rice                !----- Latent heat enthalpy
         !---------------------------------------------------------------------------------!
         !    This is the term L×[d(rs)/dt]×T. L may be either the vapourisation or        !
         ! sublimation latent heat, depending on the temperature and whether we are consi- !
         ! dering ice or not. Also, if condensation mixing ratio is constant, then this    !
         ! term will be always zero.                                                       !
         !---------------------------------------------------------------------------------!
         if (condconst) then
            ldrst = 0.
         elseif (thereisice .and. temp < t3ple) then
            ldrst = alvi*rsifp(pres,temp)*temp
         else
            ldrst = alvl*rslfp(pres,temp)*temp  
         end if
      end if

      if (newthermo) then
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli) then
            dthetail_dt = thil * (1. + (ldrst + qq)/hh) / temp
         else
            dthetail_dt = thil * (1. + ldrst*htripolii) / temp
         end if
      else
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli) then
            dthetail_dt = thil * (1. + (ldrst + qq)/(hh+qq)) / temp
         else
            dthetail_dt = thil * (1. + ldrst/(htripoli + alvl*rliq)) / temp
         end if
      end if

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
   real function thil2temp(thil,exner,pres,rliq,rice,t1stguess)
      use consts_coms, only: cp, cpi, alvl, alvi, t00, t3ple, ttripoli,htripolii,cpi4
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)   :: thil       ! Ice-liquid water potential temperature    [     K]
      real, intent(in)   :: exner      ! Exner function                            [J/kg/K]
      real, intent(in)   :: pres       ! Pressure                                  [    Pa]
      real, intent(in)   :: rliq       ! Liquid water mixing ratio                 [ kg/kg]
      real, intent(in)   :: rice       ! Ice mixing ratio                          [ kg/kg]
      real, intent(in)   :: t1stguess  ! 1st. guess for temperature                [     K]
      !----- Local variables for iterative method -----------------------------------------!
      real               :: deriv      ! Function derivative 
      real               :: fun        ! Function for which we seek a root.
      real               :: funa       ! Smallest  guess function
      real               :: funz       ! Largest   guess function
      real               :: tempa      ! Smallest  guess (or previous guess in Newton)
      real               :: tempz      ! Largest   guess (or new guess in Newton)
      real               :: delta      ! Aux. var to compute 2nd guess for bisection
      integer            :: itn,itb    ! Iteration counter
      logical            :: converged  ! Convergence handle
      logical            :: zside      ! Flag to check for one-sided approach...
      real               :: til        ! Ice liquid temperature                    [     K]
      !------------------------------------------------------------------------------------!


      !----- 1st. of all, check whether there is condensation. If not, theta_il = theta ---!
      if (rliq+rice == 0.) then
         thil2temp = cpi * thil * exner
         return
      !----- If not, check whether we are using the old thermo or the new one -------------!
      elseif (.not. newthermo) then
         til = cpi * thil * exner
         if (t1stguess > ttripoli) then
            thil2temp = 0.5 * (til + sqrt(til * (til + cpi4 * (alvl*rliq + alvi*rice))))
         else
            thil2temp = til * ( 1. + (alvl*rliq+alvi*rice) * htripolii)
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
      fun       = theta_iceliq(exner,tempz,rliq,rice)
      deriv     = dthetail_dt(.true.,fun,exner,pres,tempz,rliq,rice)
      fun       = fun - thil

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
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
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
            zside = funa*funz < 0
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

         call fatal_error('Temperature didn''t converge, giving up!!!'                     &
                       ,'thil2temp','therm_lib.f90')
      end if

      return
   end function thil2temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the partial derivative of temperature as a function of the !
   ! saturation mixing ratio [kg/kg],, keeping pressure constant. This is based on the     !
   ! ice-liquid potential temperature equation.                                            !
   !---------------------------------------------------------------------------------------!
   real function dtempdrs(exner,thil,temp,rliq,rice,rconmin)
      use consts_coms, only: alvl, alvi, cp, cpi, ttripoli, htripolii

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in) :: exner     ! Exner function                               [J/kg/K]
      real, intent(in) :: thil      ! Ice-liquid potential temperature (*)         [     K]
      real, intent(in) :: temp      ! Temperature                                  [     K]
      real, intent(in) :: rliq      ! Liquid mixing ratio                          [ kg/kg]
      real, intent(in) :: rice      ! Ice mixing ratio                             [ kg/kg]
      real, intent(in) :: rconmin   ! Minimum non-zero condensate mixing ratio     [ kg/kg]
      !------------------------------------------------------------------------------------!
      ! (*) Thil is not used in this formulation but it may be used should you opt by      !
      !     other ways to compute theta_il, so don't remove this argument.                 !
      !------------------------------------------------------------------------------------!
      !----- Local variables --------------------------------------------------------------!
      real             :: qhydm     ! Enthalpy associated with latent heat         [  J/kg]
      real             :: hh        ! Enthalpy associated with sensible heat       [  J/kg]
      real             :: rcon      ! Condensate mixing ratio                      [ kg/kg]
      real             :: til       ! Ice-liquid temperature                       [     K]
      !------------------------------------------------------------------------------------!
            
      !----- Finding the temperature and hydrometeor terms --------------------------------!
      qhydm = alvl*rliq+alvi*rice
      rcon  = rliq+rice
      if (rcon < rconmin) then
         dtempdrs = 0.
      elseif (newthermo) then
         hh    = cp*temp
         !---------------------------------------------------------------------------------!
         !    Deciding how to compute, based on temperature and whether condensates exist. !
         !---------------------------------------------------------------------------------!
         if (temp > ttripoli) then
            dtempdrs = - temp * qhydm / (rcon * (hh+qhydm))
         else
            dtempdrs = - temp * qhydm * htripolii / rcon
         end if
      else
         til   = cpi * thil * exner
         !----- Deciding how to compute, based on temperature -----------------------------!
         if (temp > ttripoli) then
            dtempdrs = - til * qhydm /( rcon * cp * (2.*temp-til))
         else
            dtempdrs = - til * qhydm * htripolii / rcon
         end if
      end if

      return
   end function dtempdrs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This fucntion computes the change of ice-liquid potential temperature due to      !
   ! sedimentation. The arguments are ice-liquid potential temperature, potential temper-  !
   ! ature and temperature in Kelvin, the old and new mixing ratio [kg/kg] and the old and !
   ! new enthalpy [J/kg].                                                                  !
   !---------------------------------------------------------------------------------------!
   real function dthil_sedimentation(thil,theta,temp,rold,rnew,qrold,qrnew)
      use consts_coms, only: ttripoli,cp,alvi,alvl

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in) :: thil  ! Ice-liquid potential temperature                 [     K]
      real, intent(in) :: theta ! Potential temperature                            [     K]
      real, intent(in) :: temp  ! Temperature                                      [     K]
      real, intent(in) :: rold  ! Old hydrometeor mixing ratio                     [ kg/kg]
      real, intent(in) :: rnew  ! New hydrometeor mixing ratio                     [ kg/kg]
      real, intent(in) :: qrold ! Old hydrometeor latent enthalpy                  [  J/kg]
      real, intent(in) :: qrnew ! New hydrometeor latent enthalpy                  [  J/kg]
      !------------------------------------------------------------------------------------!

      if (newthermo) then
         dthil_sedimentation = - thil * (alvi*(rnew-rold) - (qrnew-qrold))          &
                                        / (cp * max(temp,ttripoli))
      else
         dthil_sedimentation = - thil*thil * (alvi*(rnew-rold) - (qrnew-qrold))     &
                                        / (cp * max(temp,ttripoli) * theta)
      end if

      return
   end function dthil_sedimentation
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
   real function thetaeiv(thil,pres,temp,rvap,rtot,iflg,useice)
      use consts_coms, only : alvl,alvi,cp,ep,p00,rocp,ttripoli,t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: thil    ! Ice-liquid water potential temp.  [     K]
      real   , intent(in)           :: pres    ! Pressure                          [    Pa]
      real   , intent(in)           :: temp    ! Temperature                       [     K]
      real   , intent(in)           :: rvap    ! Water vapour mixing ratio         [ kg/kg]
      real   , intent(in)           :: rtot    ! Total mixing ratio                [ kg/kg]
      integer, intent(in)           :: iflg    ! Just to tell where this has been called.
      logical, intent(in), optional :: useice  ! Should I use ice?                 [   T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real                          :: tlcl    ! Internal LCL temperature          [     K]
      real                          :: plcl    ! Lifting condensation pressure     [    Pa]
      real                          :: dzlcl   ! Thickness of layer beneath LCL    [     m]
      !------------------------------------------------------------------------------------!

      if (present(useice)) then
         call lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,iflg,useice)
      else
         call lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,iflg)
      end if

      !------------------------------------------------------------------------------------!
      !     The definition of the thetae_iv is the thetae_ivs at the LCL. The LCL, in turn !
      ! is the point in which rtot = rvap = rsat, so at the LCL rliq = rice = 0.           !
      !------------------------------------------------------------------------------------!
      thetaeiv  = thetaeivs(thil,tlcl,rtot,0.,0.)

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
   !                 pressure at the LCL) is a function of temperature. In case you want   !
   !                 d(Thetae_ivs)/dT, use the dthetaeivs_dt function instead.             !
   !---------------------------------------------------------------------------------------!
   real function dthetaeiv_dtlcl(theiv,tlcl,rtot,eslcl,useice)
      use consts_coms, only : rocp,aklv,ttripoli
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: theiv     ! Ice-vapour equiv. pot. temp.   [      K]
      real   , intent(in)           :: tlcl      ! LCL temperature                [      K]
      real   , intent(in)           :: rtot      ! Total mixing ratio (rs @ LCL)  [     Pa]
      real   , intent(in)           :: eslcl     ! LCL saturation vapour pressure [     Pa]
      logical, intent(in), optional :: useice    ! Flag for considering ice       [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real                          :: desdtlcl  ! Saturated vapour pres. deriv.  [   Pa/K]
      !------------------------------------------------------------------------------------!



      !----- Finding the derivative of rs with temperature --------------------------------!
      if (present(useice)) then
         desdtlcl = eslifp(tlcl,useice)
      else
         desdtlcl = eslifp(tlcl)
      end if



      !----- Finding the derivative. Depending on the temperature, use different eqn. -----!
      if (tlcl > ttripoli) then
         dthetaeiv_dtlcl = theiv * (1. - rocp*tlcl*desdtlcl/eslcl - aklv*rtot/tlcl) / tlcl
      else
         dthetaeiv_dtlcl = theiv * (1. - rocp*tlcl*desdtlcl/eslcl                 ) / tlcl
      end if

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
   real function thetaeivs(thil,temp,rsat,rliq,rice)
      use consts_coms, only : aklv, ttripoli
      implicit none
      real, intent(in)   :: thil     ! Theta_il, ice-liquid water potential temp.  [     K]
      real, intent(in)   :: temp     ! Temperature                                 [     K]
      real, intent(in)   :: rsat     ! Saturation water vapour mixing ratio        [ kg/kg]
      real, intent(in)   :: rliq     ! Liquid water mixing ratio                   [ kg/kg]
      real, intent(in)   :: rice     ! Ice mixing ratio                            [ kg/kg]

      real               :: rtots    ! Saturated mixing ratio                      [     K]

      rtots = rsat+rliq+rice
      
      thetaeivs = thil * exp ( aklv * rtots / max(temp,ttripoli))

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
   real function dthetaeivs_dt(theivs,temp,pres,rsat,useice)
      use consts_coms, only : aklv,alvl,ttripoli,htripolii
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: theivs    ! Sat. ice-vap. eq. pot. temp.   [      K]
      real   , intent(in)           :: temp      ! Temperature                    [      K]
      real   , intent(in)           :: pres      ! Pressure                       [     Pa]
      real   , intent(in)           :: rsat      ! Saturation mixing ratio        [  kg/kg]
      logical, intent(in), optional :: useice    ! Flag for considering ice       [    T|F]
      !----- Local variables --------------------------------------------------------------!
      real                          :: drsdt     ! Saturated mixing ratio deriv.  [kg/kg/K]
      !------------------------------------------------------------------------------------!


      !----- Finding the derivative of rs with temperature --------------------------------!
      if (present(useice)) then
         drsdt = rslifp(pres,temp,useice)
      else
         drsdt = rslifp(pres,temp)
      end if


      !----- Finding the derivative. Depending on the temperature, use different eqn. -----!
      if (temp > ttripoli) then
         dthetaeivs_dt = theivs * (1. + aklv * (drsdt*temp-rsat)/temp ) / temp
      else
         dthetaeivs_dt = theivs * (1. + alvl * drsdt * temp * htripolii ) / temp
      end if

      
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
   real function thetaeiv2thil(theiv,pres,rtot,useice)
      use consts_coms, only : alvl,cp,ep,p00,rocp,ttripoli,t3ple,t00
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)           :: theiv     ! Ice vapour equiv. pot. temp.    [     K]
      real   , intent(in)           :: pres      ! Pressure                        [    Pa]
      real   , intent(in)           :: rtot      ! Total mixing ratio              [ kg/kg]
      logical, intent(in), optional :: useice    ! Flag for considering ice thermo [   T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real                          :: pvap      ! Sat. vapour pressure
      real                          :: theta     ! Potential temperature
      real                          :: deriv     ! Function derivative 
      real                          :: funnow    ! Function for which we seek a root.
      real                          :: funa      ! Smallest  guess function
      real                          :: funz      ! Largest   guess function
      real                          :: tlcla     ! Smallest  guess (or old guess in Newton)
      real                          :: tlclz     ! Largest   guess (or new guess in Newton)
      real                          :: tlcl      ! What will be the LCL temperature
      real                          :: es00      ! Defined as p00*rt/(epsilon + rt)
      real                          :: delta     ! Aux. variable (For 2nd guess).
      integer                       :: itn,itb   ! Iteration counters
      integer                       :: ii        ! Another counter
      logical                       :: converged ! Convergence handle
      logical                       :: zside     ! Aux. flag - check sides for Regula Falsi
      logical                       :: brrr_cold ! Flag - considering ice thermo.
      !------------------------------------------------------------------------------------!
    
      !----- Filling the flag for ice thermo that will be always present ------------------!
      if (present(useice)) then
         brrr_cold = useice
      else
         brrr_cold = bulk_on
      end if
    
      !----- Finding es00, which is a constant --------------------------------------------!
      es00 = p00 * rtot / (ep+rtot)


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
      pvap      = pres * rtot / (ep + rtot) 
      tlclz     = tslif(pvap,brrr_cold)
      theta     = tlclz * (es00/pvap)**rocp
      funnow    = thetaeivs(theta,tlclz,rtot,0.,0.)
      deriv     = dthetaeiv_dtlcl(funnow,tlclz,rtot,pvap,brrr_cold)
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
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Updating guesses ----------------------------------------------------------!
         tlcla   = tlclz
         funa    = funnow
         tlclz   = tlcla - funnow/deriv

         !----- Updating the function evaluation and its derivative -----------------------!
         pvap    = eslif(tlclz,brrr_cold)
         theta   = tlclz * (es00/pvap)**rocp
         funnow  = thetaeivs(theta,tlclz,rtot,0.,0.)
         deriv   = dthetaeiv_dtlcl(funnow,tlclz,rtot,pvap,brrr_cold)
         funnow  = funnow - theiv

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
         !   'NEWTON: it=',itn,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow           &
         !          ,'deriv=',deriv
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
            !write (unit=36,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')          &
            !   '2NGGSS: tt=',0,'tlclz=',tlclz-t00,'tlcla=',tlcla-t00,'pvap=',0.01*pvap     &
            !           ,'funa=',funa,'funz=',funnow,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlcla + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif(tlclz,brrr_cold)
               theta = tlclz * (es00/pvap)**rocp
               funz  = thetaeivs(theta,tlclz,rtot,0.,0.) - theiv

               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               !write (unit=36,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')       &
               !   '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz       &
               !           ,'delta=',delta
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

               zside = funa*funz < 0
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
                             ,'thetaeiv2thil','therm_lib.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo

            !----- Updating the guess -----------------------------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)

            !----- Updating function evaluation -------------------------------------------!
            pvap   = eslif(tlcl,brrr_cold)
            theta  = tlcl * (es00/pvap)**rocp
            funnow = thetaeivs(theta,tlcl,rtot,0.,0.) - theiv

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
         thetaeiv2thil  = theiv * exp (- alvl * rtot / (cp * max(tlcl,ttripoli)) )
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

         call fatal_error('TLCL didn''t converge, gave up!'                                &
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
      use consts_coms, only : alvl,cp,ep,p00,rocp,ttripoli,t00
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)            :: theivs     ! Sat. thetae_iv               [      K]
      real   , intent(in)            :: pres       ! Pressure                     [     Pa]
      real   , intent(out)           :: theta      ! Potential temperature        [      K]
      real   , intent(out)           :: temp       ! Temperature                  [      K]
      real   , intent(out)           :: rsat       ! Saturation mixing ratio      [  kg/kg]
      logical, intent(in) , optional :: useice     ! Flag for considering ice     [    T|F]
      !----- Local variables, with other thermodynamic properties -------------------------!
      real                           :: exnernormi ! 1./ (Norm. Exner function)   [    ---]
      logical                        :: brrr_cold  ! Flag for ice thermo          [    T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real                           :: deriv      ! Function derivative 
      real                           :: funnow     ! Function for which we seek a root.
      real                           :: funa       ! Smallest  guess function
      real                           :: funz       ! Largest   guess function
      real                           :: tempa      ! Smallest  guess (or previous in Newton)
      real                           :: tempz      ! Largest   guess (or new  in Newton)
      real                           :: delta      ! Aux. variable for 2nd guess finding.
      integer                        :: itn,itb    ! Iteration counters
      logical                        :: converged  ! Convergence handle
      logical                        :: zside      ! Aux. flag, check sides (Regula Falsi)
      !------------------------------------------------------------------------------------!
    
      !----- Setting up the ice check, in case useice is not present. ---------------------!
      if (present(useice)) then
         brrr_cold = useice
      else 
         brrr_cold = bulk_on
      end if
    
      !----- Finding the inverse of normalised Exner, which is constant in this routine ---!
      exnernormi = (p00 /pres) ** rocp

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, no idea, guess 0°C.                                            !
      !------------------------------------------------------------------------------------!
      tempz     = t00
      theta     = tempz * exnernormi
      rsat      = rslif(pres,tempz,brrr_cold)
      funnow    = thetaeivs(theta,tempz,rsat,0.,0.)
      deriv     = dthetaeivs_dt(funnow,tempz,pres,rsat,brrr_cold)
      funnow    = funnow - theivs

      !----- Saving here just in case Newton is aborted at the 1st guess ------------------!
      tempa     = tempz
      funa      = funnow

      converged = .false.
      !----- Looping ----------------------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Updating guesses ----------------------------------------------------------!
         tempa   = tempz
         funa    = funnow
         
         tempz   = tempa - funnow/deriv
         theta   = tempz * exnernormi
         rsat    = rslif(pres,tempz,brrr_cold)
         funnow  = thetaeivs(theta,tempz,rsat,0.,0.)
         deriv   = dthetaeivs_dt(funnow,tempz,pres,rsat,brrr_cold)
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
      !     If I reached this point then it's because Newton's method failed. Using bisec- !
      ! tion instead.                                                                      !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         !----- Set funz, and check whether funa and funz already have opposite sign. -----!
         funz  = funnow
         zside = .false.
         if (funa*funnow > 0.) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funz-funa) < toler*tempa) then
               delta = 100.*toler*tempa
            else
               delta = max(abs(funa*(tempz-tempa)/(funz-funa)),100.*toler*tempa)
            end if
            tempz = tempa + delta
            zgssloop: do itb=1,maxfpo
               !----- So this will be +1 -1 +2 -2 etc. ------------------------------------!
               tempz = tempz + real((-1)**itb * (itb+3)/2) * delta
               theta = tempz * exnernormi
               rsat  = rslif(pres,tempz,brrr_cold)
               funz  = thetaeivs(theta,tempz,rsat,0.,0.) - theivs
               zside = funa*funz < 0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside)                                                               &
               call fatal_error('Failed finding the second guess for regula falsi'         &
                             ,'thetaes2temp','therm_lib.f90')
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            if (abs(funz-funa) < toler*tempa) then
               temp   = 0.5*(tempa+tempz)
            else
               temp   = (funz*tempa-funa*tempz)/(funz-funa)
            end if
            theta  = temp * exnernormi
            rsat   = rslif(pres,temp,brrr_cold)
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
         rsat  = rslif(pres,temp,brrr_cold)
      else
         call fatal_error('Temperature didn''t converge, I gave up!'                       &
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
   !---------------------------------------------------------------------------------------!
   subroutine lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,iflg,useice)
      use consts_coms, only: cpog,alvl,alvi,cp,ep,p00,rocp,ttripoli,t3ple,t00,rdry
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)            :: thil   ! Ice liquid potential temp. (*)   [      K]
      real   , intent(in)            :: pres   ! Pressure                         [     Pa]
      real   , intent(in)            :: temp   ! Temperature                      [      K]
      real   , intent(in)            :: rtot   ! Total mixing ratio               [  kg/kg]
      real   , intent(in)            :: rvap   ! Vapour mixing ratio              [  kg/kg]
      real   , intent(out)           :: tlcl   ! LCL temperature                  [      K]
      real   , intent(out)           :: plcl   ! LCL pressure                     [     Pa]
      real   , intent(out)           :: dzlcl  ! Sub-LCL layer thickness          [      m]
      integer, intent(in)            :: iflg   ! Flag just to tell from where it was called
      logical, intent(in) , optional :: useice ! Should I use ice thermodynamics? [    T|F]
      !----- Local variables for iterative method -----------------------------------------!
      real                :: pvap       ! Sat. vapour pressure
      real                :: deriv      ! Function derivative 
      real                :: funnow     ! Function for which we seek a root.
      real                :: funa       ! Smallest  guess function
      real                :: funz       ! Largest   guess function
      real                :: tlcla      ! Smallest  guess (or previous guess in Newton)
      real                :: tlclz      ! Largest   guess (or new guess in Newton)
      real                :: es00       ! Defined as p00*rt/(epsilon + rt)
      real                :: delta      ! Aux. variable in case bisection is needed.
      integer             :: itn,itb    ! Iteration counters
      logical             :: converged  ! Convergence handle
      logical             :: zside      ! Aux. flag, to check sides for Regula Falsi
      !----- Other local variables --------------------------------------------------------!
      logical             :: brrr_cold ! This requires ice thermodynamics         [    T|F]

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
      es00 = p00 * rtot / (ep+rtot)


      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvap / (ep + rvap)

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a)') '----------------------------------------------------------'
      !write (unit=21,fmt='(a,1x,i5,1x,5(a,1x,f11.4,1x))')                                  &
      !   'INPUT : it=',-1,'thil=',thil,'pres=',0.01*pres,'temp=',temp-t00                  &
      !        ,'rvap=',rvap*1000.,'rtot=',rtot*1000.
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      tlclz     = 55. + 2840. / (3.5 * log(temp) - log(0.01*pvap) - 4.805) ! pvap in hPa.

      pvap      = eslif(tlclz,brrr_cold)

      funnow    = tlclz * (es00/pvap)**rocp - thil

      deriv     = (funnow+thil)*(1./tlclz - rocp*eslifp(tlclz,brrr_cold)/pvap) 

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
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Updating guesses ----------------------------------------------------------!
         tlcla  = tlclz
         funa   = funnow
         
         tlclz  = tlcla - funnow/deriv
         
         pvap   = eslif(tlclz,brrr_cold)
         funnow = tlclz * (es00/pvap)**rocp - thil
         deriv  = (funnow+thil)*(1./tlclz - rocp*eslifp(tlclz,brrr_cold)/pvap)

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
         converged = abs(tlcla-tlclz) < toler*tlclz
         if (converged) then
            tlcl = 0.5*(tlcla+tlclz)
            funz = funnow
            exit newloop
         elseif (funnow == 0.) then
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
         if (funa*funnow < 0. ) then
            funz  = funnow
            zside = .true.
         !----- They have the same sign, seeking the other guess --------------------------!
         else

            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funnow-funa) < toler*tlcla) then
               delta = 100.*toler*tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funnow-funa)),100.*toler*tlcla)
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
               tlclz = tlcla + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif(tlclz,brrr_cold)
               funz  = tlclz * (es00/pvap)**rocp - thil

               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')       &
               !   '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz       &
               !           ,'delta=',delta
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               zside = funa*funz < 0
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
               write (unit=*,fmt='(a,1x,i5)')     'CALL =',iflg
               write (unit=*,fmt='(a)') ' ============ Failed guess... ==========='
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLA =',tlcla,'FUNA =',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLZ =',tlclz,'FUNC =',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'DELTA =',delta,'FUNN =',funnow
               call fatal_error('Failed finding the second guess for regula falsi'         &
                             ,'lcl_il','therm_lib.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo

            tlcl = (funz*tlcla-funa*tlclz)/(funz-funa)

            pvap = eslif(tlcl,brrr_cold)

            funnow = tlcl * (es00/pvap)**rocp - thil

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
            converged = abs(tlcl-tlcla) < toler*tlcl .and.  abs(tlcl-tlclz) < toler*tlcl
            if (funnow == 0. .or. converged) then
               converged = .true.
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
      !----- Finding the other LCL thermodynamic variables --------------------------------!
      if (converged) then 
         pvap  = eslif(tlcl,brrr_cold)
         plcl  = (ep + rvap) * pvap / rvap
         dzlcl = max(cpog*(temp-tlcl),0.)
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
         write (unit=*,fmt='(a,1x,i5)'    ) 'call            [   ---] =',iflg
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
         call fatal_error('TLCL didn''t converge, gave up!','lcl_il','therm_lib.f90')
      end if
      return
   end subroutine lcl_il
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! internal energy .                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine qtk(q,tempk,fracliq)
      use consts_coms, only: cliqi,cicei,allii,t3ple,qicet3,qliqt3,tsupercool
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: q        ! Internal energy                             [   J/kg]
      real, intent(out) :: tempk    ! Temperature                                 [      K]
      real, intent(out) :: fracliq  ! Liquid Fraction (0-1)                       [    ---]
      !------------------------------------------------------------------------------------!


      !----- Internal energy below qwfroz, all ice  ---------------------------------------!
      if (q <= qicet3) then
         fracliq = 0.
         tempk   = q * cicei 
      !----- Internal energy, above qwmelt, all liquid ------------------------------------!
      elseif (q >= qliqt3) then
         fracliq = 1.
         tempk   = q * cliqi + tsupercool
      !----- Changing phase, it must be at freezing point ---------------------------------!
      else
         fracliq = (q-qicet3) * allii
         tempk   = t3ple
      endif
      !------------------------------------------------------------------------------------!

      return
   end subroutine qtk
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from inter-  !
   ! nal energy (J/m² or J/m³), mass (kg/m² or kg/m³), and heat capacity (J/m²/K or        !
   ! J/m³/K).                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine qwtk(qw,w,dryhcap,tempk,fracliq)
      use consts_coms, only: cliqi,cliq,cicei,cice,allii,alli,t3ple,tsupercool
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: qw      ! Internal energy                   [  J/m²] or [  J/m³]
      real, intent(in)  :: w       ! Density                           [ kg/m²] or [ kg/m³]
      real, intent(in)  :: dryhcap ! Heat capacity of nonwater part    [J/m²/K] or [J/m³/K]
      real, intent(out) :: tempk   ! Temperature                                   [     K]
      real, intent(out) :: fracliq ! Liquid fraction (0-1)                         [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real              :: qwfroz  ! qw of ice at triple point         [  J/m²] or [  J/m³] 
      real              :: qwmelt  ! qw of liquid at triple point      [  J/m²] or [  J/m³]
      !------------------------------------------------------------------------------------!

      !----- Converting melting heat to J/m² or J/m³ --------------------------------------!
      qwfroz = (dryhcap + w*cice) * t3ple
      qwmelt = qwfroz   + w*alli
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      !    This is analogous to the qtk computation, we should analyse the magnitude of    !
      ! the internal energy to choose between liquid, ice, or both by comparing with our.  !
      ! know boundaries.                                                                   !
      !------------------------------------------------------------------------------------!

      !----- Internal energy below qwfroz, all ice  ---------------------------------------!
      if (qw < qwfroz) then
         fracliq = 0.
         tempk   = qw  / (cice * w + dryhcap)
      !----- Internal energy, above qwmelt, all liquid ------------------------------------!
      elseif (qw > qwmelt) then
         fracliq = 1.
         tempk   = (qw + w * cliq * tsupercool) / (dryhcap + w*cliq)
      !------------------------------------------------------------------------------------!
      !    We are at the freezing point.  If water mass is so tiny that the internal       !
      ! energy of frozen and melted states are the same given the machine precision, then  !
      ! we assume that water content is negligible and we impose 50% frozen for            !
      ! simplicity.                                                                        !
      !------------------------------------------------------------------------------------!
      elseif (qwfroz == qwmelt) then
         fracliq = 0.5
         tempk   = t3ple
      !------------------------------------------------------------------------------------!
      !    Changing phase, it must be at freezing point.  The max and min are here just to !
      ! avoid tiny deviations beyond 0. and 1. due to floating point arithmetics.          !
      !------------------------------------------------------------------------------------!
      else
         fracliq = min(1.,max(0.,(qw - qwfroz) * allii / w))
         tempk   = t3ple
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine qwtk
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib
!==========================================================================================!
!==========================================================================================!
