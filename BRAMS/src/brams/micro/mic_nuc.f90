!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine cldnuc(m1,k1cnuc,k2cnuc,flpw,rv,wp,i,j)

use micphys

implicit none

integer :: m1,k1cnuc,k2cnuc,i,j
real :: flpw
real, dimension(m1) :: rv,wp

integer :: k,jtemp,jw,jconcen,iw,iconc,lpw
real :: rnuc,excessrv,rcnew,tab,concen_tab,cxadd  &
   ,tairc_nuc,w_nuc,rjw,wtw2,wtw1,concen_nuc,rjconcen,wtconcen2,wtconcen1

real, dimension(9,7,7) :: cldnuctab
data ((cldnuctab(iw,iconc,1),iw=1,9),iconc=1,7)/                 &
  .307,  .520,  .753,  .919,  .990,  .990,  .990,  .990,  .990,  &
  .230,  .426,  .643,  .860,  .969,  .990,  .990,  .990,  .990,  &
  .164,  .336,  .552,  .777,  .940,  .990,  .990,  .990,  .990,  &
  .098,  .254,  .457,  .701,  .892,  .979,  .990,  .990,  .990,  &
  .045,  .145,  .336,  .614,  .822,  .957,  .990,  .990,  .990,  &
  .018,  .073,  .206,  .426,  .672,  .877,  .969,  .990,  .990,  &
  .008,  .027,  .085,  .206,  .280,  .336,  .336,  .336,  .906/
  
data ((cldnuctab(iw,iconc,2),iw=1,9),iconc=1,7)/                 &
  .230,  .426,  .643,  .860,  .969,  .990,  .990,  .990,  .990,  &
  .164,  .336,  .552,  .777,  .930,  .990,  .990,  .990,  .990,  &
  .112,  .254,  .457,  .701,  .877,  .974,  .990,  .990,  .990,  &
  .073,  .184,  .365,  .583,  .822,  .949,  .990,  .990,  .990,  &
  .038,  .112,  .254,  .489,  .727,  .906,  .982,  .990,  .990,  &
  .015,  .054,  .145,  .365,  .614,  .841,  .957,  .990,  .990,  &
  .005,  .018,  .073,  .184,  .395,  .614,  .800,  .940,  .990/
  
data ((cldnuctab(iw,iconc,3),iw=1,9),iconc=1,7)/                 &
  .164,  .336,  .552,  .800,  .949,  .990,  .990,  .990,  .990,  &
  .128,  .254,  .457,  .701,  .892,  .979,  .990,  .990,  .990,  &
  .085,  .184,  .365,  .583,  .822,  .949,  .990,  .990,  .990,  &
  .054,  .128,  .280,  .489,  .727,  .906,  .982,  .990,  .990,  &
  .027,  .085,  .206,  .395,  .643,  .841,  .963,  .990,  .990,  &
  .012,  .038,  .112,  .280,  .520,  .777,  .930,  .990,  .990,  &
  .004,  .015,  .054,  .145,  .365,  .614,  .822,  .949,  .990/
  
data ((cldnuctab(iw,iconc,4),iw=1,9),iconc=1,7)/                 &
  .145,  .280,  .489,  .727,  .919,  .990,  .990,  .990,  .990,  &
  .098,  .206,  .395,  .614,  .841,  .963,  .990,  .990,  .990,  &
  .063,  .145,  .307,  .520,  .753,  .919,  .990,  .990,  .990,  &
  .038,  .098,  .230,  .426,  .643,  .860,  .963,  .990,  .990,  &
  .022,  .063,  .164,  .336,  .552,  .777,  .930,  .990,  .990,  &
  .010,  .027,  .085,  .230,  .426,  .701,  .877,  .974,  .990,  &
  .003,  .012,  .038,  .112,  .280,  .552,  .777,  .940,  .990/
  
data ((cldnuctab(iw,iconc,5),iw=1,9),iconc=1,7)/                 &
  .112,  .230,  .457,  .701,  .892,  .982,  .990,  .990,  .990,  &
  .073,  .164,  .336,  .552,  .800,  .940,  .990,  .990,  .990,  &
  .054,  .112,  .254,  .457,  .672,  .877,  .979,  .990,  .990,  &
  .032,  .085,  .184,  .365,  .583,  .800,  .940,  .990,  .990,  &
  .018,  .045,  .128,  .254,  .457,  .701,  .892,  .979,  .990,  &
  .008,  .022,  .073,  .184,  .365,  .614,  .822,  .949,  .990,  &
  .003,  .010,  .032,  .098,  .230,  .489,  .727,  .906,  .979/
  
data ((cldnuctab(iw,iconc,6),iw=1,9),iconc=1,7)/                 &
  .098,  .206,  .395,  .643,  .860,  .974,  .990,  .990,  .990,  &
  .063,  .145,  .307,  .520,  .753,  .930,  .990,  .990,  .990,  &
  .045,  .098,  .206,  .395,  .643,  .841,  .963,  .990,  .990,  &
  .027,  .063,  .145,  .307,  .520,  .753,  .919,  .990,  .990,  &
  .015,  .038,  .098,  .230,  .426,  .643,  .841,  .963,  .990,  &
  .007,  .018,  .054,  .145,  .307,  .552,  .777,  .919,  .990,  &
  .003,  .008,  .027,  .073,  .206,  .395,  .672,  .860,  .969/
  
data ((cldnuctab(iw,iconc,7),iw=1,9),iconc=1,7)/                 &
  .098,  .206,  .365,  .614,  .841,  .969,  .990,  .990,  .990,  &
  .054,  .128,  .280,  .489,  .727,  .906,  .990,  .990,  .990,  &
  .038,  .085,  .184,  .365,  .583,  .822,  .957,  .990,  .990,  &
  .022,  .063,  .128,  .280,  .457,  .701,  .892,  .982,  .990,  &
  .012,  .038,  .085,  .184,  .365,  .583,  .822,  .949,  .990,  &
  .005,  .018,  .045,  .128,  .280,  .489,  .727,  .892,  .979,  &
  .002,  .007,  .022,  .063,  .164,  .365,  .614,  .822,  .949/

k1cnuc = nint(flpw)
k2cnuc = 1

if (jnmb(1) == 1 .or. jnmb(1) == 4) then

! cloud number specified in parm(1)   

   rnuc = parm(1) * emb0(1)

   lpw = nint(flpw)

   do k = lpw,m1-1
      excessrv = rv(k) - 1.0001 * rvlsair(k)
      rcnew = 0.
      if (excessrv > 0.) then
         rcnew = min(rnuc,.5*excessrv)
         rx(k,1) = rx(k,1) + rcnew
         rv(k) = rv(k) - rcnew
         k2cnuc = k
         cx(k,1) = min(parm(1),rx(k,1) / emb0(1))
      elseif (k2cnuc == 1) then
         k1cnuc = k + 1
      endif
   enddo

elseif (jnmb(1) >= 5) then

! cloud number predicted from ccn field

   lpw = nint(flpw)

   do k = lpw,m1-1
      excessrv = rv(k) - 1.0001 * rvlsair(k)

      if (excessrv > 0.) then

          tairc_nuc = tairc(k)
          if (tairc_nuc < -30.) then
             tairc_nuc = -30.
          elseif (tairc_nuc > 30.) then
             tairc_nuc = 30.
          endif
          jtemp = nint(.1 * (tairc_nuc + 30.)) + 1

          w_nuc = wp(k)
          if (w_nuc < .010001) then 
             w_nuc = .010001
          elseif (w_nuc > 99.99) then
             w_nuc = 99.99
          endif
          rjw = 2. * log10(100. * w_nuc) + 1.
          jw = int(rjw)
          wtw2 = rjw - float(jw)
          wtw1 = 1. - wtw2

          if (jnmb(1) == 5) then
             concen_nuc = cparm     ! ccn concen const value specified in cparm
          elseif (jnmb(1) == 6) then
            ! concen_nuc = prof(k)   ! ccn concen specified vertical profile
          elseif (jnmb(1) == 7) then
             concen_nuc = cccnx(k)  ! ccn concen predicted in cccnp
          else
             print*, 'icloud set to value greater than 7: ',icloud
             print*, 'stopping model '
             stop 'icloud'
          endif  

          if (concen_nuc < 10.001e6) then
             concen_nuc = 10.001e6
          elseif (concen_nuc > 9999.e6) then
             concen_nuc = 9999.e6
          endif
          rjconcen = 2. * log10(1.e-7 * concen_nuc) + 1.
          jconcen = int(rjconcen)
          wtconcen2 = rjconcen - float(jconcen)
          wtconcen1 = 1. - wtconcen2

          tab = wtconcen1 * (wtw1 * cldnuctab(jw  ,jconcen  ,jtemp)   &
                          +  wtw2 * cldnuctab(jw+1,jconcen  ,jtemp))  &
              + wtconcen2 * (wtw1 * cldnuctab(jw  ,jconcen+1,jtemp)   &
                          +  wtw2 * cldnuctab(jw+1,jconcen+1,jtemp))

          concen_tab = concen_nuc * tab
          
! Nucleate cloud droplets only if concen_tab > existing cloud concentration

          if (concen_tab > cx(k,1)) then
             cxadd = concen_tab - cx(k,1)
             if (cxadd > excessrv / emb0(1)) cxadd = excessrv / emb0(1)
             cx(k,1) = cx(k,1) + cxadd
             rx(k,1) = rx(k,1) + excessrv
             k2cnuc = k
          endif

      elseif (k2cnuc.eq.1) then
         k1cnuc = k + 1
      endif

   enddo

else
   print*, 'icloud not allowed to be 2 or 3'
   print*, 'stopping model '
   stop 'icloud'
endif

return
end subroutine cldnuc

!******************************************************************************

subroutine icenuc(m1,kc1,kc2,k1pnuc,k2pnuc,flpw,ngr,rv,dn0,dtlt,i,j)

use rconstants
use micphys

implicit none

integer :: m1,kc1,kc2,k1pnuc,k2pnuc,ngr,i,j,k,idnc,itc,irhhz,ithz
real :: dn1,fraccld,ridnc,dtlt,ssi0,wdnc2,tc,ritc,wtc2  &
       ,pbvi,ptvi,pdvi,ptotvi,fracifn,cldnuc,cldnucr,rhhz,haznuc  &
       ,rirhhz,wrhhz2,thz,rithz,wthz2,frachaz,ssi,diagni  &
       ,vapnuc,vapnucr,availvap
real :: flpw
integer :: lpw
real, dimension(m1) :: rv,dn0

! Define ssi0 to be maximum supersaturation with respect to ice for
! determining total number of IFN that can nucleate in Meyers' formula
data ssi0/0.40/
save

!
! implement paul's immersion freezing of rain here.  This would
! replace mike's homogeneous freezing of rain which was in h03.
!

do k = kc1,kc2

!  Homogeneous ice nucleation of cloud droplets

! define dn locally from emb

      dn1 = dnfac(1) * emb(k,1) ** pwmasi(1)

   fraccld = 0.

   if (rx(k,1) > rxmin(1) .and. tairc(k) <= -30.01) then

      ridnc = max(1.,min(float(ndnc-1),dn1 / ddnc))
      idnc = int(ridnc)
      wdnc2 = ridnc - float(idnc)

      tc = max(-49.99,tairc(k))
      ritc = (tc + 50.00) / dtc + 1.0
      itc = int(ritc)
      wtc2 = ritc - float(itc)
      fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,ngr)  &
              +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,ngr)  &
              + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,ngr)  &
              +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,ngr)

   endif

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)

   call contnuc (rx(k,1),cx(k,1),tx(k,1),vap(k,1),press(k)  &
      ,dynvisc(k),thrmcon(k),tair(k),tairc(k)  &
      ,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt,i,k)

! progIFN: Scale ptotvi returned from contnuc by prognosed IFN fraction

!::later   ptotvi = ptotvi * fracifn

! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
! probably do not want sink of ccinp here.

   cldnuc = ptotvi + max(0.,fraccld * cx(k,1) - cx(k,3))
!         cldnucr = cldnuc * emb(k,1)
   cldnucr = min(rx(k,1),ptotvi * emb(k,1) + fraccld * rx(k,1))

   rx(k,3) = rx(k,3) + cldnucr
   rx(k,1) = rx(k,1) - cldnucr
   cx(k,3) = cx(k,3) + cldnuc
   cx(k,1) = cx(k,1) - cldnuc

enddo

! DEMOTT'S NEW SCHEME: In 4.3 and beyond, assume that it gives #/KG

!  Homogeneous nucleation of haze

k1pnuc = 2
k2pnuc = 1

lpw = nint(flpw)

do k = lpw,m1-1
   rhhz = rv(k) / rvlsair(k)
   haznuc = 0.
   if (rhhz .gt. 0.82 .and. tairc(k) .le. -35.01) then
      rirhhz = min(0.1799,rhhz-0.82) / drhhz + 1.0
      irhhz = int(rirhhz)
      wrhhz2 = rirhhz - float(irhhz)
      thz = max(-59.99,tairc(k))
      rithz = (thz + 60.00) / dthz + 1.0
      ithz = int(rithz)
      wthz2 = rithz - float(ithz)
      frachaz = (1.-wrhhz2) * (1.-wthz2) * frachz(irhhz  ,ithz  )  &
              +     wrhhz2  * (1.-wthz2) * frachz(irhhz+1,ithz  )  &
              + (1.-wrhhz2) *     wthz2  * frachz(irhhz  ,ithz+1)  &
              +     wrhhz2  *     wthz2  * frachz(irhhz+1,ithz+1)
      frachaz = 1. - exp(-frachaz * dtlt)
! OPTION 1
      haznuc = frachaz * 300.e6
! OPTION 2
!           haznuc = frachaz * caero(k)
   endif

! meyers -  no cld aerosol source or sink here

!  Heterogeneous nucleation by deposition condensation freezing
!  with deposition nuclei.  In 4.3 and beyond, assume that it gives #/kg.

   ssi = min(ssi0,rv(k) / rvisair(k) - 1.)

   if (ssi .gt. 0. .and. tairc(k) .le. -5.) then
      fracifn = exp(12.96 * (ssi - ssi0))
   else
      fracifn = 0.
   endif

! Diagnose maximum number of IFN to activate based on ipris

   if (ipris .eq. 5) then
      diagni = fracifn * 1.e5
   elseif (ipris .eq. 6) then
      diagni = fracifn * dn0(k) ** 5.4 * 1.e5
   elseif (ipris .eq. 7) then
      diagni = fracifn * cifnx(k)
   endif

! orig Meyers formula:     +      diagni = exp(6.269 + 12.96 * ssi)

!  Combine nucleation types, and limit amounts
! vapnuc is #/kg_air and vapnucr is kg/kg_air

! BEGIN MIKE'S SECTION FOR LIMITING NUMBER OF CRYSTALS NUCLEATED
! BY NUMBER OF ICE CRYSTALS PRESENT ALREADY

   vapnuc = max(0.,haznuc + diagni - cx(k,3))
   vapnucr = vapnuc * emb0(3)
   if (vapnucr .gt. 0.) then
      availvap = .5 * (rv(k) - rvisair(k))
      if (vapnucr .gt. availvap) then
         vapnucr = min(vapnucr, max(0.,availvap))
      endif
   endif
   vapnuc = vapnucr / emb0(3)

   rx(k,3) = rx(k,3) + vapnucr
   cx(k,3) = cx(k,3) + vapnuc

   if (rx(k,3) > rxmin(3)) k2pnuc = k
   if (k2pnuc .eq. 1 .and. rx(k,3) < rxmin(3)) k1pnuc = k + 1

enddo

! here mike has the habit diagnosis. option 1 is to use habit
! at cloud top, option 2 is to use new habit at each level.
! need to consider other options.  how about method of formation?
! my question about how much of habit is due to existing ice
! structure, and how much is due to current growth environment
! (temp and supsat). relevant supsat is wrt liquid?

return
end

!******************************************************************************

subroutine contnuc (rx,cx,tx,vap,press  &
   ,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt,i,k)

implicit none

integer :: i,k
real :: rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi  &
       ,dn1,dtlt,aka,raros,ana,akn,dfar,f1,f2,ft
data aka,raros/5.39e-3,3.e-7/

!  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
!  phoresis, thermophoresis, and Brownian motion (transport of IN)
!
!  ana   = # IN per kg available for contact freezing (from Meyers et al. 1992
!          where ana was interpreted as # per m^3)
!  akn   = Knudsen number (Walko et al. 1995, Eq. 58)
!          [2.28e-5 = mfp * p00 / 293.15]
!  raros = aerosol radius = 3.e-7 m from Cotton et al. (1986)
!  dfar  = aerosol diffusivity (Pruppacher and Klett Eq. 12-15)
!          [7.32e-25 = Boltzmann constant / (6 pi)]
!  f1    = "function 1" (Walko et al. 1995 Eq. 55) multiplied by delta t
!           but now cld concen in #/kg_air so (pvbi, ptvi, pdvi) all per kg_air
!  f2    = "function 2" (Walko et al. 1995 Eq. 56)
!  ft    = "function ft" (Walko et al. 1995 Eq. 57)
!  pbvi  = Brownian motion nucleation amount this timestep [#/kg_air]
!  ptvi  = Thermophoretic nucleation amount this timestep [#/kg_air]
!  pdvi  = Diffusiophoretic nucleation amount this timestep [#/kg_air],
!          reformulated to use vapor diffusion directly.  Factor of 1.2
!          is (1+sigma_va x_a) from Pruppacher and Klett Eq. 12-102
!          divided by .622, the molecular weight ratio between water and air.

   ptotvi = 0.

   if (tx .le. -2. .and. rx .gt. 1.e-10) then

      ana = exp(4.11 - 0.262 * tx)
      akn = 2.28e-5 * tair / (press * raros)
      dfar = 7.32e-25 * tair * (1.+ akn) / (raros * dynvisc)
      f1 = 6.28318 * dn1 * cx * ana * dtlt
      f2 = thrmcon * (tairc - tx) / press
      ft = 0.4 * (1. + 1.45 * akn + 0.4 * akn * exp(-1. / akn))  &
         * (thrmcon + 2.5 * akn * aka)  &
         / ((1. + 3. * akn)  &
         * (2. * thrmcon + 5. * aka * akn + aka))
      pbvi = f1 * dfar
      ptvi = f1 * f2 * ft
      pdvi = 1.2 * ana * vap
      ptotvi = max(0.,pbvi + ptvi + pdvi)

   endif
   return
   end
