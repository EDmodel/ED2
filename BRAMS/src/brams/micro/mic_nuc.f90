!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine cldnuc(m1)

   use micphys
   use micro_coms, only : cldnuctab

   implicit none

   !----- Argument ------------------------------------------------------------------------!
   integer, intent(in) :: m1
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,jtemp,jw,jconcen,iw,iconc
   real    :: rnuc,excessrv,rcnew,tab,concen_tab,cxadd,tairc_nuc,w_nuc,rjw,wtw2,wtw1
   real    :: concen_nuc,rjconcen,wtconcen2,wtconcen1
   !---------------------------------------------------------------------------------------!
   k1cnuc = lpw
   k2cnuc = 1

   select case (jnmb(1))
   !----- Cloud number specified in parm(1) -----------------------------------------------!
   case (1,4)
      rnuc = parm(1) * emb0(1)
      do k = lpw,m1-1
         excessrv = rvap(k) - 1.0001 * rvlsair(k)
         rcnew = 0.
         if (excessrv > 0.) then
            rcnew = min(rnuc,.5*excessrv)
            rx(k,1) = rx(k,1) + rcnew
            rvap(k) = rvap(k) - rcnew
            k2cnuc = k
            cx(k,1) = min(parm(1),rx(k,1) / emb0(1))
         elseif (k2cnuc == 1) then
            k1cnuc = k + 1
         end if
      end do

   !----- Cloud number predicted from ccn field -------------------------------------------!
   case (5:7)


      do k = lpw,m1-1
         excessrv = rvap(k) - 1.0001 * rvlsair(k)

         if (excessrv > 0.) then

            tairc_nuc = min(30.,max(-30.,tairc(k)))
            jtemp     = nint(.1 * (tairc_nuc + 30.)) + 1

            w_nuc = min(99.99,max(.010001,vertvelo(k)))
            rjw   = 2. * log10(100. * w_nuc) + 1.
            jw    = int(rjw)
            wtw2  = rjw - real(jw)
            wtw1  = 1. - wtw2

            select case (jnmb(1))
            case (5)
               concen_nuc = cparm    !----- ccn concen const value specified in cparm -----!
            case (6)
              concen_nuc = cparm     !----- Temporarily 6 = 5 -----------------------------!
              ! concen_nuc = prof(k) !----- ccn concen specified vertical profile ---------!
            case(7)
               concen_nuc = cccnx(k) !----- ccn concen predicted in cccnp -----------------!
            case default
               call abort_run('ICLOUD set to value greater than 7!!!'                      &
                             ,'cldnuc','mic_nuc.f90')
            end select

            concen_nuc = min(9999.e6,max(10.001e6,concen_nuc))
            rjconcen   = 2. * log10(1.e-7 * concen_nuc) + 1.
            jconcen    = int(rjconcen)
            wtconcen2  = rjconcen - float(jconcen)
            wtconcen1  = 1. - wtconcen2

            tab = wtconcen1 * (wtw1 * cldnuctab(jw  ,jconcen  ,jtemp)                      &
                            +  wtw2 * cldnuctab(jw+1,jconcen  ,jtemp))                     &
                + wtconcen2 * (wtw1 * cldnuctab(jw  ,jconcen+1,jtemp)                      &
                            +  wtw2 * cldnuctab(jw+1,jconcen+1,jtemp))

            concen_tab = concen_nuc * tab
            !------------------------------------------------------------------------------!
            !    Nucleate cloud droplets only if concen_tab is greater than existing       !
            ! cloud concentration.                                                         !
            !------------------------------------------------------------------------------!
            if (concen_tab > cx(k,1)) then
               cxadd = concen_tab - cx(k,1)
               if (cxadd > excessrv / emb0(1)) cxadd = excessrv / emb0(1)
               cx(k,1) = cx(k,1) + cxadd
               rx(k,1) = rx(k,1) + excessrv
               k2cnuc = k
            end if

         elseif (k2cnuc.eq.1) then
            k1cnuc = k + 1
         end if
      end do

   case default
      call abort_run ('ICLOUD not allowed to be 2 or 3','cldnuc','mic_nuc.f90')
   end select

   return
end subroutine cldnuc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine icenuc(m1,ngr,dtlt)

   use rconstants
   use micphys
   use micro_coms, only: ssi0
   use therm_lib , only: rehuil,rehui

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: m1,ngr
   real   , intent(in) :: dtlt
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,idnc,itc,irhhz,ithz
   real                :: dn1,fraccld,ridnc,wdnc2,tc,ritc,wtc2
   real                :: pbvi,ptvi,pdvi,ptotvi,fracifn,cldnuc,cldnucr,rhhz,haznuc
   real                :: rirhhz,wrhhz2,thz,rithz,wthz2,frachaz,ssi,diagni
   real                :: vapnuc,vapnucr,availvap
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Implement Paul's immersion freezing of rain here.  This would replace Mike's homo- !
   ! geneous freezing of rain which was in h03.                                            !
   !---------------------------------------------------------------------------------------!
   do k = k1(1),k1(2)
   !----- Homogeneous ice nucleation of cloud droplets ------------------------------------!
      dn1 = dnfac(1) * emb(k,1) ** pwmasi(1)

      fraccld = 0.

      if (rx(k,1) > rxmin(1) .and. tairc(k) <= -30.01) then
         ridnc = max(1.,min(real(ndnc-1),dn1 / ddnc))
         idnc  = int(ridnc)
         wdnc2 = ridnc - real(idnc)

         tc = max(-49.99,tairc(k))
         ritc = (tc + 50.00) / dtc + 1.0
         itc = int(ritc)
         wtc2 = ritc - float(itc)
         fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,ngr)  &
                 +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,ngr)  &
                 + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,ngr)  &
                 +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,ngr)

      end if

      !------------------------------------------------------------------------------------!
      !     Heterogeneous contact ice nucleation of cloud droplets by diffusiophoresis,    !
      ! thermophoresis, and Brownian motion (transport of IN).                             !
      !------------------------------------------------------------------------------------!
      call contnuc (rx(k,1),cx(k,1),tx(k,1),vap(k,1),press(k),dynvisc(k) ,thrmcon(k)       &
                   ,tair(k),tairc(k),pbvi,ptvi,pdvi,ptotvi,dn1,dtlt,rxmin(1))

      !------------------------------------------------------------------------------------!
      ! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction              !
      !::later   ptotvi = ptotvi * fracifn                                                 !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but probably do not want   !
      ! sink of ccinp here.                                                                !
      !------------------------------------------------------------------------------------!
      cldnuc = ptotvi + max(0.,fraccld * cx(k,1) - cx(k,3))
      !cldnucr = cldnuc * emb(k,1)
      cldnucr = min(rx(k,1),ptotvi * emb(k,1) + fraccld * rx(k,1))

      rx(k,3) = rx(k,3) + cldnucr
      rx(k,1) = rx(k,1) - cldnucr
      cx(k,3) = cx(k,3) + cldnuc
      cx(k,1) = cx(k,1) - cldnuc

   end do

   !---------------------------------------------------------------------------------------!
   ! DEMOTT'S NEW SCHEME: In 4.3 and beyond, assume that it gives #/KG.                    !
   !---------------------------------------------------------------------------------------!
   !----- Homogeneous nucleation of haze --------------------------------------------------!
   k1pnuc = lpw
   k2pnuc = 1

   do k = lpw,m1-1
      !------------------------------------------------------------------------------------!
      ! THERMODYNAMIC DILEMMA: Shouldn't it be rehuil instead of rehul here, since we      !
      !     seek haze nucleation at cold temperatures?                                     !
      !------------------------------------------------------------------------------------!
      rhhz = rehuil(press(k),tair(k),rvap(k))
      haznuc = 0.
      if (rhhz > 0.82 .and. tairc(k) <= -35.01) then
         rirhhz  = min(0.1799,rhhz-0.82) / drhhz + 1.0
         irhhz   = int(rirhhz)
         wrhhz2  = rirhhz - float(irhhz)
         thz     = max(-59.99,tairc(k))
         rithz   = (thz + 60.00) / dthz + 1.0
         ithz    = int(rithz)
         wthz2   = rithz - float(ithz)
         frachaz = (1.-wrhhz2) * (1.-wthz2) * frachz(irhhz  ,ithz  )  &
                 +     wrhhz2  * (1.-wthz2) * frachz(irhhz+1,ithz  )  &
                 + (1.-wrhhz2) *     wthz2  * frachz(irhhz  ,ithz+1)  &
                 +     wrhhz2  *     wthz2  * frachz(irhhz+1,ithz+1)
         frachaz = 1. - exp(-frachaz * dtlt)
         !----- OPTION 1 ------------------------------------------------------------------!
         haznuc = frachaz * 300.e6
         ! OPTION 2
         ! haznuc = frachaz * caero(k)
      end if

      !----- Meyers -  no cloud aerosol source or sink here. ------------------------------!
      !------------------------------------------------------------------------------------!
      !     Heterogeneous nucleation by deposition condensation freezing with deposition   !
      ! nuclei.  In 4.3 and beyond, assume that it gives #/kg.                             !
      !------------------------------------------------------------------------------------!
      ssi = min(ssi0,rehui(press(k),tair(k),rvap(k)) - 1.)
      if (ssi > 0. .and. tairc(k) <= -5.) then
         fracifn = exp(12.96 * (ssi - ssi0))
      else
         fracifn = 0.
      end if

      !----- Diagnose maximum number of IFN to activate based on ipris --------------------!
      select case (ipris)
      case (5)
         diagni = fracifn * 1.e5
      case (6)
         diagni = fracifn * rhoa(k) ** 5.4 * 1.e5
      case (7)
         diagni = fracifn * cifnx(k)
      end select
      !----- orig Meyers formula:     +      diagni = exp(6.269 + 12.96 * ssi) ------------!

      !------------------------------------------------------------------------------------!
      !     Combine nucleation types, and limit amounts vapnuc is #/kg_air and vapnucr is  !
      ! kg/kg_air.                                                                         !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Begin Mike's section for limiting number of crystals nucleated by number of    !
      ! ice crystals present already.                                                      !
      !------------------------------------------------------------------------------------!
      vapnuc  = max(0.,haznuc + diagni - cx(k,3))
      vapnucr = vapnuc * emb0(3)
      if (vapnucr > 0.) then
         availvap = .5 * (rvap(k) - rvisair(k))
         if (vapnucr > availvap) vapnucr = min(vapnucr, max(0.,availvap))
      end if
      vapnuc = vapnucr / emb0(3)

      rx(k,3) = rx(k,3) + vapnucr
      cx(k,3) = cx(k,3) + vapnuc

      if (rx(k,3) > rxmin(3)) k2pnuc = k
      if (k2pnuc == 1 .and. rx(k,3) < rxmin(3)) k1pnuc = k + 1

   end do

   !---------------------------------------------------------------------------------------!
   !     Here Mike has the habit diagnosis. option 1 is to use habit at cloud top, option  !
   ! 2 is to use new habit at each level. Need to consider other options.  How about       !
   ! method of formation? My question about how much of habit is due to existing ice       !
   ! structure, and how much is due to current growth environment (temp and supsat).       !
   ! Relevant supsat is wrt liquid?                                                        !
   !---------------------------------------------------------------------------------------!
   return
end subroutine icenuc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Heterogeneous contact ice nucleation of cloud droplets by diffusiophoresis, thermo- !
! phoresis, and Brownian motion (transport of IN):                                         !
!                                                                                          !
!  ana   = # IN per kg available for contact freezing (from Meyers et al. 1992             !
!          where ana was interpreted as # per m^3)                                         !
!  akn   = Knudsen number (Walko et al. 1995, Eq. 58)                                      !
!  dfar  = aerosol diffusivity (Pruppacher and Klett Eq. 12-15)                            !
!  f1    = "function 1" (Walko et al. 1995 Eq. 55) multiplied by delta t                   !
!           but now cld concen in #/kg_air so (pvbi, ptvi, pdvi) all per kg_air            !
!  f2    = "function 2" (Walko et al. 1995 Eq. 56)                                         !
!  ft    = "function ft" (Walko et al. 1995 Eq. 57)                                        !
!  pbvi  = Brownian motion nucleation amount this timestep [#/kg_air]                      !
!  ptvi  = Thermophoretic nucleation amount this timestep [#/kg_air]                       !
!  pdvi  = Diffusiophoretic nucleation amount this timestep [#/kg_air],                    !
!          reformulated to use vapor diffusion directly.  Factor of 1.2                    !
!          is (1+sigma_va x_a) from Pruppacher and Klett Eq. 12-102                        !
!          divided by .622, the molecular weight ratio between water and air.              !
!------------------------------------------------------------------------------------------!
subroutine contnuc (rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi    &
                   ,dn1,dtlt,rxmin)

   use micro_coms, only : &
           aka        & ! aerosol thermal conductivity
          ,raros      & ! aerosol radius = 3.e-7 m from Cotton et al. (1986)
          ,boltzo6pi  & ! Boltzmann constant / (6*pi)
          ,w95_58     & ! all constants of Walko et al. 1995 equation 58.
          ,ticenucmin ! ! Minimum temperature for ice to nucleate (C)
   
   use rconstants, only : &
           t00        & ! zero Celsius
          ,twopi      ! ! 2*pi
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)  :: rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc,dn1,dtlt,rxmin
   real, intent(out) :: pbvi,ptvi,pdvi,ptotvi
   !----- Local variables -----------------------------------------------------------------!
   real              :: ana,akn,dfar,f1,f2,ft
   !---------------------------------------------------------------------------------------!
   ptotvi = 0.

   if (tx <= ticenucmin .and. rx > rxmin) then
      ana = exp(4.11 - 0.262 * (tx-t00))
      akn = w95_58 * tair / (press * raros)
      dfar = boltzo6pi * tair * (1.+ akn) / (raros * dynvisc)

      f1 = twopi * dn1 * cx * ana * dtlt
      f2 = thrmcon * (tair - tx) / press

      ft = 0.4 * (1. + 1.45 * akn + 0.4 * akn * exp(-1. / akn))                            &
         * (thrmcon + 2.5 * akn * aka) / ((1. + 3. * akn)                                  &
         * (2. * thrmcon + 5. * aka * akn + aka))

      pbvi   = f1 * dfar
      ptvi   = f1 * f2 * ft
      pdvi   = 1.2 * ana * vap
      ptotvi = max(0.,pbvi + ptvi + pdvi)

   end if
   return
end subroutine contnuc
!==========================================================================================!
!==========================================================================================!
