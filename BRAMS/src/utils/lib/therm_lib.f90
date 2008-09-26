!==========================================================================================!
! therm_lib.f90                                                                            !
!                                                                                          !
!    This file contains various subroutines to handle common thermodynamic properties.     !
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
   !  These equations give the triple point at 273.16, with vapour pressure being 
   !---------------------------------------------------------------------------------------!
   !----- Coefficients based on equation (7): ---------------------------------------------!
   real, dimension(0:3), parameter :: iii_7 = (/ 9.550426,-5723.265, 3.53068,-0.00728332 /)
   !----- Coefficients based on equation (10), first fit ----------------------------------!
   real, dimension(0:3), parameter :: l01_10= (/54.842763,-6763.22 ,-4.210  , 0.000367   /)
   !----- Coefficients based on equation (10), second fit ---------------------------------!
   real, dimension(0:3), parameter :: l02_10= (/53.878   ,-1331.22 ,-9.44523, 0.014025   /)
   !----- Coefficients based on the hyperbolic tangent ------------------------------------!
   real, dimension(2)  , parameter :: ttt_10= (/0.0415,218.8/)

   contains
   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour pressure as a function of   !
   ! Kelvin temperature. This expression came from MK05, equation (10).                    !
   !---------------------------------------------------------------------------------------!
   real function eslf(temp,l1funout,l2funout,ttfunout)
      use rconstants, only : t00
      implicit none
      real, parameter             :: c0= .6105851e+03, c1= .4440316e+02, c2= .1430341e+01
      real, parameter             :: c3= .2641412e-01, c4= .2995057e-03, c5= .2031998e-05
      real, parameter             :: c6= .6936113e-08, c7= .2564861e-11, c8=-.3704404e-13
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
         eslf = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

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
      use rconstants, only : t00
      implicit none
      real, parameter             :: c0= .6114327e+03, c1= .5027041e+02, c2= .1875982e+01
      real, parameter             :: c3= .4158303e-01, c4= .5992408e-03, c5= .5743775e-05
      real, parameter             :: c6= .3566847e-07, c7= .1306802e-09, c8= .2152144e-12
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
         esif=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

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
      use rconstants, only: t3ple
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
      use rconstants, only : ep,toodry
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
      use rconstants, only : ep,toodry
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
      use rconstants, only: t3ple,ep
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
      use rconstants, only : rm
      implicit none
      real, intent(in) :: temp
      real             :: eequ
      eequ = eslf(temp)
      rhovsl = eequ / (rm * temp)
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
      use rconstants, only : rm
      implicit none
      real, intent(in) :: temp
      real             :: eequ
      eequ = esif(temp)
      rhovsi = eequ / (rm * temp)
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
      use rconstants, only : rm
      implicit none
      real, intent(in)              :: temp
      logical, intent(in), optional :: useice
      real                          :: eequ

      if (present(useice)) then
         eequ = eslif(temp,useice)
      else
         eequ = eslif(temp)
      end if

      rhovsil = eequ / (rm * temp)

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
      use rconstants, only: t00
      implicit none
      real, parameter  :: d0= .4443216e+02, d1= .2861503e+01, d2= .7943347e-01
      real, parameter  :: d3= .1209650e-02, d4= .1036937e-04, d5= .4058663e-07
      real, parameter  :: d6=-.5805342e-10, d7=-.1159088e-11, d8=-.3189651e-14
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
         eslfp=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
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
      use rconstants, only: lsorvap, t00
      implicit none
      real, parameter  :: d0= .5036342e+02, d1= .3775758e+01, d2= .1269736e+00
      real, parameter  :: d3= .2503052e-02, d4= .3163761e-04, d5= .2623881e-06
      real, parameter  :: d6= .1392546e-08, d7= .4315126e-11, d8= .5961476e-14
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
         esifp=d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
      end if

      return
   end function esifp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour mix-  !
   ! ing ratio with respect to temperature as a function of pressure and Kelvin tempera-   !
   ! ture.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real function rslfp(pres,temp)
      use rconstants, only: ep
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
      use rconstants, only: ep
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
      use rconstants, only: t3ple
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
         delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            funz  = eslf(tempz) - pvap
            zside = funa*funz < 0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            call abort_run('Failed finding the second guess for regula falsi'           &
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
         call abort_run('Temperature didn''t converge, giving up!!!'      &
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
         delta = max(abs(funa * (tempz-tempa)/(fun-funa)),100.*toler*tempa)
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
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            call abort_run('Failed finding the second guess for regula falsi'           &
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
         call abort_run('Temperature didn''t converge, giving up!!!'      &
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
      use rconstants, only: es3ple,alvl,alvi,lvt3ple,lst3ple,rmt3ple,es3plei

      implicit none
      real, intent(in) :: pvap
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
      use rconstants, only: ep,toodry

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
      use rconstants, only: ep,toodry

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
      use rconstants, only: ep,toodry

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
      use rconstants, only: ep,toodry

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
      use rconstants, only: ep,toodry

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
      use rconstants, only: ep,toodry,t3ple

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
      use rconstants, only: ep,toodry
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
      use rconstants, only: ep,toodry
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
      use rconstants, only: t3ple
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
      use rconstants, only: epi
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
      use rconstants, only: epi
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
   !     This fucntion computes the ice liquid potential temperature given the Exner       !
   ! function [J/kg/K], temperature [K], and liquid and ice mixing ratios [kg/kg].         !
   !---------------------------------------------------------------------------------------!
   real function theta_iceliq(exner,temp,rliq,rice)
      use rconstants, only: alvl, alvi, cp, ttripoli, htripoli, htripolii

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in) :: exner ! Exner function                                   [J/kg/K]
      real, intent(in) :: temp  ! Temperature                                      [     K]
      real, intent(in) :: rliq  ! Liquid mixing ratio                              [ kg/kg]
      real, intent(in) :: rice  ! Ice mixing ratio                                 [ kg/kg]
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
            theta_iceliq = hh * exp(qq/hh) / exner
         else
            theta_iceliq = hh * exp(qq * htripolii) / exner
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
   real function dthetail_dt(exner,pres,temp,rsat,rliq,ricein)
      use rconstants, only: alvl, alvi, cp, ttripoli,htripoli,htripolii,t3ple

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)           :: exner       ! Exner function                   [J/kg/K]
      real, intent(in)           :: pres        ! Pressure                         [    Pa]
      real, intent(in)           :: temp        ! Temperature                      [     K]
      real, intent(in)           :: rsat        ! Saturation mixing ratio          [ kg/kg]
      real, intent(in)           :: rliq        ! Liquid mixing ratio              [ kg/kg]
      real, intent(in), optional :: ricein      ! Ice mixing ratio                 [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real                       :: rice        ! Ice mixing ratio or 0.           [ kg/kg]
      real                       :: thil        ! Ice-liquid potential temperature [     K]
      real                       :: ldrst       ! L × d(rs)/dT × T                 [  J/kg]
      real                       :: hh          ! Sensible heat enthalpy           [  J/kg]
      real                       :: qq          ! Latent heat enthalpy             [  J/kg]
      logical                    :: thereisice  ! Is ice present                   [   ---]
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
         dthetail_dt = cp/exner
         return
      else
         thil  = theta_iceliq(exner,temp,rliq,rice) !----- Ice liquid potential temp.
         hh    = cp*temp                            !----- Sensible heat enthalpy
         qq    = alvl*rliq+alvi*rice                !----- Latent heat enthalpy
         !---------------------------------------------------------------------------------!
         !    This is the term L×[d(rs)/dt]×T. L may be either the vapourisation or        !
         ! sublimation latent heat, depending on the temperature and whether we are consi- !
         ! dering ice or not.                                                              !
         !---------------------------------------------------------------------------------!
         if (thereisice .and. temp < t3ple) then
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
   real function thil2temp(thil,exner,rliq,rice,t1stguess)
      use rconstants, only: cp, cpi, alvl, alvi, t00, t3ple, ttripoli,htripolii,cpi4
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)   :: thil       ! Ice-liquid water potential temperature    [     K]
      real, intent(in)   :: exner      ! Exner function                            [J/kg/K]
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
      real               :: alpha      ! Part inside the exponential in theta_il formula
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

      !----- If not, iterate: For the Newton's 1st. guess, ignore rliq and rice -----------!
      tempz     = t1stguess
      alpha     = - (alvl*rliq + alvi*rice)/ (cp*max(tempz,ttripoli))
      fun       = cp * tempz/ exner * exp(alpha) - thil
      if (tempz >= ttripoli) then
         deriv     = (fun+thil)*(1.-alpha)/tempz
      else
         deriv     = (fun+thil)/tempz
      end if

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,3(1x,a,1x,es12.5))')            &
      !   'itn=',0,'bisection=',.false.,'tempz=',tempz-t00,'alpha=',alpha                   &
      !           ,'fun=',fun,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      tempa = tempz
      funa  = fun
      !----- Enter loop: it will probably skip when the air is not saturated --------------!
      converged = .false.
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, go to bisection -------!
         tempa = tempz
         funa  = fun
         tempz = tempa - fun/deriv
         alpha     = - (alvl*rliq + alvi*rice)/ (cp*max(tempz,ttripoli))
         fun       = cp * tempz/ exner * exp(alpha) - thil
         if (tempz >= ttripoli) then
            deriv     = (fun+thil)*(1.-alpha)/tempz
         else
            deriv     = (fun+thil)/tempz
         end if

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,f11.4,3(1x,a,1x,es12.5))')         &
         !   'itn=',itn,'bisection=',.false.,'tempz=',tempz-t00,'alpha=',alpha              &
         !           ,'fun=',fun,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         converged = abs(tempa-tempz)/tempz < toler .or. fun == 0.
         !----- Converged, happy with that, return the average b/w the 2 previous guesses -!
         if (converged) then
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
         delta = -funa * (tempz-tempa)/(fun-funa)
         tempz = tempa + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
            alpha     = - (alvl*rliq + alvi*rice)/ (cp*max(tempz,ttripoli))
            funz      = cp * tempz/ exner * exp(alpha) - thil
            zside     = funa*funz < 0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            write (unit=*,fmt='(2(a,1x,es14.7))') 'delta=',delta,'alpha=',alpha
            call abort_run('Failed finding the second guess for regula falsi'           &
                          ,'thil2temp','therm_lib.f90')
         end if
      end if


      bisloop: do itb=itn,maxfpo
         thil2temp =  (funz*tempa-funa*tempz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(thil2temp-tempa)/thil2temp < toler
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         alpha     = - (alvl*rliq + alvi*rice)/ (cp*max(thil2temp,ttripoli))
         fun       =  cp*thil2temp/exner*exp(alpha)-thil

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,a,1x,es12.5,6(1x,a,1x,f11.4))')         &
         !   'itn=',itb,'bisection=',.true.,'alpha=',alpha                                  &
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

         call abort_run('Temperature didn''t converge, giving up!!!'      &
                       ,'thil2temp','therm_lib.f90')
      end if

      return
   end function thil2temp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the derivative of temperature as a function of the satur-  !
   ! ation mixing ratio [kg/kg], based on the ice-liquid potential temperature equation.   !
   !---------------------------------------------------------------------------------------!
   real function dtempdrs(exner,thil,temp,rliq,rice,rconmin)
      use rconstants, only: alvl, alvi, cp, cpi, ttripoli, htripolii

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
      use rconstants, only: ttripoli,cpi

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in) :: thil  ! Ice-liquid potential temperature                 [     K]
      real, intent(in) :: theta ! Potential temperature (*)                        [     K]
      real, intent(in) :: temp  ! Temperature                                      [     K]
      real, intent(in) :: rold  ! Old hydrometeor mixing ratio                     [ kg/kg]
      real, intent(in) :: rnew  ! New hydrometeor mixing ratio                     [ kg/kg]
      real, intent(in) :: qrold ! Old hydrometeor latent enthalpy                  [  J/kg]
      real, intent(in) :: qrnew ! New hydrometeor latent enthalpy                  [  J/kg]
      !------------------------------------------------------------------------------------!


      if (newthermo) then
         dthil_sedimentation = - thil * (2820.*(rnew-rold) - cpi * (qrnew-qrold))          &
                                        / (max(temp,ttripoli))
      else
         dthil_sedimentation = - thil*thil * (2820.*(rnew-rold) - cpi * (qrnew-qrold))     &
                                        / (max(temp,ttripoli) * theta)
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
   !---------------------------------------------------------------------------------------!
   real function thetaeiv(thil,pres,temp,rvap,rtot,tlclout)
      use rconstants, only : alvl,alvi,cp,ep,p00,rocp,ttripoli,t3ple
      implicit none
      real, intent(in)            :: thil    ! Ice-liquid water potential temp.    [     K]
      real, intent(in)            :: pres    ! Pressure                            [    Pa]
      real, intent(in)            :: temp    ! Temperature                         [     K]
      real, intent(in)            :: rvap    ! Water vapour mixing ratio           [ kg/kg]
      real, intent(in)            :: rtot    ! Total mixing ratio                  [ kg/kg]
      real, intent(out), optional :: tlclout ! Lifting condensation temperature    [     K]
      !----- Arguments --------------------------------------------------------------------!
      !----- Local variables for iterative method -----------------------------------------!
      real               :: tlcl       ! TLCL current guess
      real               :: pvap       ! Sat. vapour pressure
      real               :: deriv      ! Function derivative 
      real               :: funnow     ! Function for which we seek a root.
      real               :: funa       ! Smallest  guess function
      real               :: funz       ! Largest   guess function
      real               :: tlcla      ! Smallest  guess (or previous guess in Newton)
      real               :: tlclz      ! Largest   guess (or new guess in Newton)
      real               :: es00       ! Defined as p00*rt/(epsilon + rt)
      real               :: delta      ! Aux. variable in case bisection is needed.
      integer            :: itn,itb    ! Iteration counters
      logical            :: converged  ! Convergence handle
      logical            :: zside      ! Aux. flag, to check sides for Regula Falsi
      !------------------------------------------------------------------------------------!
    
      !----- Finding es00, which is a constant --------------------------------------------!
      es00 = p00 * rtot / (ep+rtot)

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvap / (ep + rvap) 
      tlclz     = 55. + 2840. / (3.5 * log(temp) - log(0.01*pvap) - 4.805) ! pvap in hPa.
      pvap      = eslif(tlclz)
      funnow    = tlclz * (es00/pvap)**rocp - thil
      if (tlclz >= t3ple) then
         deriv  = (funnow+thil)*(cp*tlclz+ep*alvl)/(cp*tlclz*tlclz)
      else
         deriv  = (funnow+thil)*(cp*tlclz+ep*alvi)/(cp*tlclz*tlclz)
      end if
      converged = abs(funnow) > tlclz*epsilon(1.)
      !----- Looping ----------------------------------------------------------------------!
      if (converged) then
         tlcl      = tlclz
      else
         tlcla     = tlclz
         funa      = funnow
         newloop: do itn=1,maxfpo/6
            if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection --!
            !----- Updating guesses -------------------------------------------------------!
            tlcla   = tlclz
            funa    = funnow
            
            tlclz   = tlcla - funnow/deriv

            pvap = eslif(tlclz)
            funnow = tlclz * (es00/pvap)**rocp - thil
            if (tlclz >= t3ple) then
               deriv  = (funnow+thil)*(cp*tlclz+ep*alvl)/(cp*tlclz*tlclz)
            else
               deriv  = (funnow+thil)*(cp*tlclz+ep*alvi)/(cp*tlclz*tlclz)
            end if
            !------------------------------------------------------------------------------!
            !   Convergence may happen when we are lucky or get close guesses.             !
            !------------------------------------------------------------------------------!
            !----- 1. We are lucky. -------------------------------------------------------!
            converged = abs(funnow) < epsilon(1.)*tlclz
            if (converged) then
               tlcl = tlclz
               exit newloop
            end if
            !----- 2. The two guesses are really close ------------------------------------!
            converged = abs(tlcla-tlclz)/tlclz < toler
            if (converged) then
               tlcl = 0.5*(tlcla+tlclz)
               exit newloop
            end if
         end do newloop
      end if

      if (.not. converged) then
         !---------------------------------------------------------------------------------!
         !     If I reached this point then it's because Newton's method failed. Using bi- !
         ! section instead. First, I need to find two guesses that give me functions with  !
         ! opposite signs. If funa and funnow have opposite signs, then we are all set.    !                                                                      !
         !---------------------------------------------------------------------------------!
         if (funa*funnow < 0. ) then
            funz  = funnow
            zside = .true.
         !----- They have the same sign, seeking the other guess --------------------------!
         else
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            delta = -funa*(tlclz-tlcla)/(funnow-funa)
            tlclz = tlcla + delta
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlcla + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif(tlclz)
               funz  = tlclz * (es00/pvap)**rocp - thil
               zside = funa*funz < 0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(a)') ' No second guess for you...'
               write (unit=*,fmt='(2(a,1x,es14.7))') 'tlcla=',tlcla,'funa=',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'tlclz=',tlclz,'func=',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'delta=',delta,'funn=',funnow
               call abort_run('Failed finding the second guess for regula falsi'           &
                             ,'thetaeiv','therm_lib.f90')
            end if
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            tlcl = (funz*tlcla-funa*tlclz)/(funz-funa)
            pvap = eslif(tlcl)
            funnow = tlcl * (es00/pvap)**rocp - thil

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla)/tlcl < toler
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
         thetaeiv  = thil * exp ( alvl * rtot / (cp * max(tlcl,ttripoli)) )
         if (present(tlclout)) tlclout = tlcl
      else
         call abort_run('Thetae_iv didn''t converge, gave up!','thetaeiv','therm_lib.f90')
      end if

      return
   end function thetaeiv
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
      use rconstants, only : aklv, ttripoli
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
   !    This function finds the ice-liquid potential temperature from the ice-vapour equi- !
   ! valent potential temperature. It is also an iterative method since we need the LCL    !
   ! temperature and we don't have it. We find T(LCL) using (you guessed it right!) an     !
   ! iterative method, using, as usual Newton's method as a starting point, and if it      !
   ! fails, fall back to the modified regula falsi (Illinois method).                      !
   !---------------------------------------------------------------------------------------!
   real function thetaeiv2thil(thetaeiv,pres,rtot,tlclout)
      use rconstants, only : alvl,alvi,cp,ep,p00,rocp,ttripoli,t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)            :: thetaeiv ! Ice vapour equivalent pot. temp.   [     K]
      real, intent(in)            :: pres     ! Pressure                           [    Pa]
      real, intent(in)            :: rtot     ! Total mixing ratio                 [ kg/kg]
      real, intent(out), optional :: tlclout  ! Temperature @ lifting cond. level  [     K]
      !----- Local variables for iterative method -----------------------------------------!
      real               :: pvap       ! Sat. vapour pressure
      real               :: deriv      ! Function derivative 
      real               :: funnow     ! Function for which we seek a root.
      real               :: funa       ! Smallest  guess function
      real               :: funz       ! Largest   guess function
      real               :: tlcla      ! Smallest  guess (or previous guess in Newton)
      real               :: tlclz      ! Largest   guess (or new guess in Newton)
      real               :: tlcl       ! What will be the LCL temperature
      real               :: es00       ! Defined as p00*rt/(epsilon + rt)
      real               :: delta      ! Aux. variable in case bisection is needed.
      integer            :: itn,itb    ! Iteration counters
      logical            :: converged  ! Convergence handle
      logical            :: zside      ! Aux. flag, to check sides for Regula Falsi
      !------------------------------------------------------------------------------------!
    
      !----- Finding es00, which is a constant --------------------------------------------!
      es00 = p00 * rtot / (ep+rtot)

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, we assume we are lucky and right at the LCL.                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rtot / (ep + rtot) 
      tlclz     = tslif(pvap) 
      funnow    = tlclz * (es00/pvap)**rocp * exp(alvl*rtot/(cp*max(tlclz,ttripoli)))      &
                - thetaeiv
      if (tlclz >= t3ple) then
         deriv  = (funnow+thetaeiv)*(cp*tlclz-alvl*ep-alvl*rtot)/(cp*tlclz*tlclz)
      elseif (tlclz >= ttripoli) then
         deriv  = (funnow+thetaeiv)*(cp*tlclz-alvi*ep-alvl*rtot)/(cp*tlclz*tlclz)
      else
         deriv  = (funnow+thetaeiv)*(cp*tlclz-alvi*ep)/(cp*tlclz*tlclz)
      end if
      tlcla     = tlclz
      funa      = funnow
      converged = .false.
      !----- Looping ----------------------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Updating guesses ----------------------------------------------------------!
         tlcla   = tlclz
         funa    = funnow
         
         tlclz   = tlcla - funnow/deriv

         pvap = eslif(tlclz)
         funnow = tlclz * (es00/pvap)**rocp * exp(alvl*rtot/(cp*max(tlclz,ttripoli)))      &
                - thetaeiv
         if (tlclz >= t3ple) then
            deriv  = (funnow+thetaeiv)*(cp*tlclz-alvl*ep-alvl*rtot)/(cp*tlclz*tlclz)
         elseif (tlclz >= ttripoli) then
            deriv  = (funnow+thetaeiv)*(cp*tlclz-alvi*ep-alvl*rtot)/(cp*tlclz*tlclz)
         else
            deriv  = (funnow+thetaeiv)*(cp*tlclz-alvi*ep)/(cp*tlclz*tlclz)
         end if

         converged = abs(tlcla-tlclz)/tlclz < toler
         if (converged) then
            tlcl = 0.5*(tlcla+tlclz)
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
            delta = -funa*(tlclz-tlcla)/(funz-funa)
            tlclz = tlcla + delta
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlclz + delta
               pvap  = eslif(tlclz)
               funz  = tlclz*(es00/pvap)**rocp * exp(alvl*rtot/(cp*max(tlclz,ttripoli)))   &
                     - thetaeiv
               zside = funa*funz < 0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside)                                                               &
               call abort_run('Failed finding the second guess for regula falsi'           &
                             ,'thetaeiv2thil','rthrm.f90')
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            tlcl = (funz*tlcla-funa*tlclz)/(funz-funa)
            pvap = eslif(tlcl)
            funnow = tlcl*(es00/pvap)**rocp * exp(alvl*rtot/(cp*max(tlcl,ttripoli)))       &
                   - thetaeiv

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla)/tlcl < toler
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
         thetaeiv2thil  = thetaeiv * exp (- alvl * rtot / (cp * max(tlcl,ttripoli)) )
         if (present(tlclout)) tlclout = tlcl
      else
         call abort_run('Theta_il didn''t converge, gave up!'                              &
                       ,'thetaeiv2thil','therm_lib.f90')
      end if

      return
   end function thetaeiv2thil
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the equivalent potential temperature, given pressure in  !
   ! Pa, temperature in Kelvin, vapour mixing ratio in kg/kg.                              !
   !---------------------------------------------------------------------------------------!
   subroutine thetae(pres,temp,rvap,the)
      use rconstants, only : alvl,cp,ep,p00,rocp,ttripoli
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)            :: temp    ! Temperature                         [     K]
      real, intent(in)            :: pres    ! Pressure                            [    Pa]
      real, intent(in)            :: rvap    ! Water vapour mixing ratio           [ kg/kg]
      real, intent(out)           :: the     ! Equivalent potential temperature    [     K]
      !----- Local variables for iterative method -----------------------------------------!
      real               :: theta      ! Potential temperature 
      real               :: tlcl       ! TLCL current guess
      real               :: pvap       ! Sat. vapour pressure
      real               :: deriv      ! Function derivative 
      real               :: funnow     ! Function for which we seek a root.
      real               :: funa       ! Smallest  guess function
      real               :: funz       ! Largest   guess function
      real               :: tlcla      ! Smallest  guess (or previous guess in Newton)
      real               :: tlclz      ! Largest   guess (or new guess in Newton)
      real               :: es00       ! Defined as p00*rt/(epsilon + rt)
      real               :: delta      ! Aux. variable in case bisection is needed.
      integer            :: itn,itb    ! Iteration counters
      logical            :: converged  ! Convergence handle
      logical            :: zside      ! Aux. flag, to check sides for Regula Falsi
      !------------------------------------------------------------------------------------!
    
      !----- Finding es00 and theta, which are constants ----------------------------------!
      es00  = p00 * rvap / (ep+rvap)
      theta = temp * (p00 / pres) ** rocp

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvap / (ep + rvap) 
      tlclz     = 55. + 2840. / (3.5 * log(temp) - log(0.01*pvap) - 4.805) ! pvap in hPa.
      pvap      = eslf(tlclz)
      funnow    = tlclz * (es00/pvap)**rocp - theta
      deriv     = (funnow+theta)*(cp*tlclz-ep*alvl)/(cp*tlclz*tlclz)
      tlcla     = tlclz
      funa      = funnow
      converged = .false.
      !----- Looping ----------------------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         if (abs(deriv) < toler) exit newloop !----- Too dangerous, skip to bisection -----!
         !----- Updating guesses ----------------------------------------------------------!
         tlcla   = tlclz
         funa    = funnow
         
         tlclz   = tlcla - funnow/deriv

         pvap = eslf(tlclz)
         funnow = tlclz * (es00/pvap)**rocp - theta
         deriv  = (funnow+theta)*(cp*tlclz-ep*alvl)/(cp*tlclz*tlclz)

         converged = abs(tlcla-tlclz)/tlclz < toler
         if (converged) then
            tlcl = 0.5*(tlcla+tlclz)
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
         if (funa*funnow > 0.) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            delta = max(abs(funa*(tlclz-tlcla)/(funz-funa)),100.*toler*tlcla)
            tlclz = tlcla + delta
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So this will be +1 -1 +2 -2 etc. ------------------------------------!
               tlclz = tlclz + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslf(tlclz)
               funz  = tlclz * (es00/pvap)**rocp - theta
               zside = funa*funz < 0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside)                                                               &
               call abort_run('Failed finding the second guess for regula falsi'           &
                             ,'thetae','therm_lib.f90')
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            tlcl = (funz*tlcla-funa*tlclz)/(funz-funa)
            pvap = eslf(tlcl)
            funnow = tlcl * (es00/pvap)**rocp - theta

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla)/tlcl < toler
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
         the  = theta * exp ( alvl * rvap / (cp * max(tlcl,ttripoli)) )
      else
         call abort_run('Thetae didn''t converge, gave up!','thetae','therm_lib.f90')
      end if

      return
   end subroutine thetae
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine converts saturated equivalent potential temperature into tempera-  !
   ! ture, given pressure in Pa, and theta_es in Kelvin. It also returns potential temper- !
   ! ature in Kelvin, and saturation vapour mixing ratio in kg/kg. Since this is equiva-   !
   ! lent potential temperature, it will never consider ice. As usual,  we seek T using    !
   ! Newton's method as a starting point, and if it fails, we fall back to the modified    !
   ! regula falsi (Illinois method).                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine the2t(thes,pres,theta,temp,rsat)
      use rconstants, only : alvl,cp,ep,p00,rocp,ttripoli,t00
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)   :: thes       ! Sat. equivalent potential temperature    [      K]
      real, intent(in)   :: pres       ! Pressure                                 [     Pa]
      real, intent(out)  :: theta      ! Potential temperature                    [      K]
      real, intent(out)  :: temp       ! Temperature                              [      K]
      real, intent(out)  :: rsat       ! Saturation mixing ratio                  [  kg/kg]
      !----- Local variables, with other thermodynamic properties -------------------------!
      real               :: exnernormi ! 1./ (Normalised Exner function)          [    ---]
      real               :: drsdt      ! d(rsat)/dT, keeping pressure constant    [kg/kg/K]
      !----- Local variables for iterative method -----------------------------------------!
      real               :: deriv      ! Function derivative 
      real               :: funnow     ! Function for which we seek a root.
      real               :: funa       ! Smallest  guess function
      real               :: funz       ! Largest   guess function
      real               :: tempa      ! Smallest  guess (or previous guess in Newton)
      real               :: tempz      ! Largest   guess (or new guess in Newton)
      real               :: delta      ! Aux. variable in case bisection is needed.
      integer            :: itn,itb    ! Iteration counters
      logical            :: converged  ! Convergence handle
      logical            :: zside      ! Aux. flag, to check sides for Regula Falsi
      !------------------------------------------------------------------------------------!
    
      !----- Finding the inverse of normalised Exner, which is constant in this routine ---!
      exnernormi = (p00 /pres) ** rocp

      !------------------------------------------------------------------------------------!
      !     The 1st. guess, no idea, guess 0°C.                                            !
      !------------------------------------------------------------------------------------!
      tempz     = t00
      theta     = tempz * exnernormi
      rsat      = rslf(pres,tempz)
      drsdt     = rslfp(pres,tempz)
      funnow    = theta * exp(alvl*rsat / (cp*max(tempz,ttripoli))) - thes
      if (tempz >= ttripoli) then
         deriv  = (funnow+thes)*(cp*tempz-alvl*(rsat - drsdt*tempz))/(cp*tempz*tempz)
      else
         deriv  = (funnow+thes)*(cp*ttripoli+alvl*drsdt*tempz)/(cp*tempz*ttripoli)
      end if
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
         rsat    = rslf(pres,tempz)
         drsdt   = rslfp(pres,tempz)
         funnow  = theta * exp(alvl*rsat / (cp*max(tempz,ttripoli))) - thes
         if (tempz >= ttripoli) then
            deriv  = (funnow+thes)*(cp*tempz-alvl*(rsat - drsdt*tempz))/(cp*tempz*tempz)
         else
            deriv  = (funnow+thes)*(cp*ttripoli+alvl*drsdt*tempz)/(cp*tempz*ttripoli)
         end if

         converged = abs(tempa-tempz)/tempz < toler
         if (converged) then
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
         funz = funnow
         if (funa*funnow > 0.) then
            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            delta = max(abs(funa*(tempz-tempa)/(funz-funa)),100.*toler*tempa)
            tempz = tempa + delta
            zside = .false.
            zgssloop: do itb=1,maxfpo
               !----- So this will be +1 -1 +2 -2 etc. ------------------------------------!
               tempz = tempz + real((-1)**itb * (itb+3)/2) * delta
               theta = tempz * exnernormi
               rsat  = rslf(pres,tempz)
               funz  = theta * exp(alvl*rsat / (cp*max(tempz,ttripoli))) - thes
               zside = funa*funz < 0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside)                                                               &
               call abort_run('Failed finding the second guess for regula falsi'           &
                             ,'the2t','therm_lib.f90')
         end if
         !---- Continue iterative method --------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            temp   = (funz*tempa-funa*tempz)/(funz-funa)
            theta  = temp * exnernormi
            rsat   = rslf(pres,temp)
            funnow = theta * exp(alvl*rsat / (cp*max(temp,ttripoli))) - thes

            !------------------------------------------------------------------------------!
            !    Checking for convergence. If it did, return, we found the solution.       !
            ! Otherwise, constrain the guesses.                                            !
            !------------------------------------------------------------------------------!
            converged = abs(temp-tempa)/temp < toler
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
               if (.not.zside) funz = funz * 0.5
               !----- We just updated aside, setting zside to false -----------------------!
               zside = .false. 
            end if
         end do fpoloop
      end if

      if (converged) then 
         !----- Compute theta and rsat with temp just for consistency ---------------------!
         theta = temp * exnernormi
         rsat  = rslf(pres,temp)
      else
         call abort_run('Temperature didn''t converge, I gave up!'                         &
                       ,'the2t','therm_lib.f90')
      end if

      return
   end subroutine the2t
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
      use rconstants, only: cliqi,cicei,allii,alli,tsupercool,t3ple
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
      use rconstants, only: t00
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
      use rconstants, only: cliqi,cliq,cicei,cice,allii,alli,t3ple
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: qw      ! Internal energy                   [  J/m²] or [  J/m³]
      real, intent(in)  :: w       ! Density                           [ kg/m²] or [ kg/m³]
      real, intent(in)  :: dryhcap ! Heat capacity of nonwater part    [J/m²/K] or [J/m³/K]
      real, intent(out) :: tempk   ! Temperature                                   [     K]
      real, intent(out) :: fracliq ! Liquid fraction (0-1)                         [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real              :: qwliq0  ! qw of liquid at triple point      [  J/m²] or [  J/m³]
      real              :: ch2ow   ! heat capacity of water            [  J/m²] or [  J/m³]
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
         fracliq = qw / qwliq0
         tempk = t3ple
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine qwtk
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib
!==========================================================================================!
!==========================================================================================!
