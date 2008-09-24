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
   integer, parameter ::   level=3

   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical, parameter ::   vapour_on=.true.
   logical, parameter ::   cloud_on=.true.
   logical, parameter ::   bulk_on=.true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These constants came from the paper in which the saturation vapour pressure is    !
   ! based on:                                                                             !
   !                                                                                       !
   !  Murphy, D. M.; Koop, T., 2005: Review of the vapour pressures of ice and supercooled !
   !     water for atmospheric applications. Q. J. Royal Meteor. Soc., vol. 31, pp. 1539-  !
   !     1565 (hereafter MK05).                                                            !
   !                                                                                       !
   !  These equations give the triple point at 273.16.                                     !
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
      use consts_coms, only : t00
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
      use consts_coms, only : t00
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
      use consts_coms, only: t00
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

