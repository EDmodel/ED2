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
!    This subroutine build the haze nucleation table.                                      !
!------------------------------------------------------------------------------------------!
subroutine haznuc()

   use micphys, only : &
        nthz,          & !intent(in)
        dthz,          & !intent(in)
        nrhhz,         & !intent(in)
        drhhz,         & !intent(in)
        frachz           !intent(out)
   use rconstants, only : pio6,onethird
   implicit none

   !----- Local Variables: ----------------------------------------------------------------!
   integer :: ithz,irhhz,k
   real :: denccn,gnuccn,dnccn,ddccn,rhhz,c1hz,c2hz,c3hz,bhz,dm,sum,dccn,y,dum,thz
   real :: gammln

   denccn = 1.769
   gnuccn = 1.
   dnccn =   .075E-4
   ddccn = .005e-4
   do ithz = 1,nthz
      thz = -60. + dthz * float(ithz - 1)
      do irhhz = 1,nrhhz
         rhhz = 0.82 + drhhz * float(irhhz - 1)
         c1hz = (pio6 * denccn) ** (-onethird)
         c2hz = -14.65 - 1.045 * thz
         c3hz = -492.35 - 8.34 * thz - 0.0608 * thz ** 2
         bhz = min(38., max(-38., c2hz + c3hz * (1. - rhhz)))
         dm = c1hz * 10 ** (-bhz/6.)

         sum = 0.
         dccn = 0.
         do k=1,200
            dccn = dccn + ddccn
            y=dccn / dnccn
            !------------------------------------------------------------------------------!
            !  This IF is needed to avoid underflow...                                     !
            !------------------------------------------------------------------------------!
            if (abs(dccn / dm) < 1.e-4) then
               dum=0.
            else
               dum=min(50., (dccn / dm) ** 6)
            end if
            sum = sum + y ** (gnuccn-1.) * exp(-y) * (1. - exp(-dum))
         end do
         frachz(irhhz,ithz) = sum*ddccn/(exp(gammln(gnuccn))*dnccn)
      end do
   end do

   return
end subroutine haznuc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine produces the table for homogeneous freezing of cloud droplets. Most  !
! internal computation was switched to double precision to avoid overflow/underflow.       !
!------------------------------------------------------------------------------------------!
subroutine homfrzcl(dtlt,ngr)

   use micphys, only: &
        ntc,          & !INTENT(IN)
        ndnc,         & !INTENT(IN)
        dtc8,         & !INTENT(IN)
        ddnc8,        & !INTENT(IN)
        fracc           !INTENT(OUT)

   implicit none

   !----- Arguments: ----------------------------------------------------------------------!
   integer, intent(in) :: ngr 
   real   , intent(in) :: dtlt
   !----- Local Variables: ----------------------------------------------------------------!
   integer                 :: itc,k,idnc
   real(kind=8)            :: ajlso,dnc,osum,dc,v1,tc,y,dfracc,expvar,dtlt8
   real                    :: gammln
   real(kind=8), parameter :: gnuc=1., ddc=0.5d-6
   real(kind=8), parameter :: tinyexp=-38. ! A number small enough to assume exp to be zero
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Make table for homogeneous freezing of cloud droplets. Need gnuc = gnu(1) ???     !
   !---------------------------------------------------------------------------------------!
   dtlt8 = dble(dtlt)
   do itc = 1,ntc
      tc = -50. + dtc8*dble(itc-1)
      y = -(606.3952+tc*(52.6611+tc*(1.7439+tc*(.0265+tc*1.536e-4))))
      ajlso = 1.e6 * 10. ** y
      do idnc = 1,ndnc
         dnc  = ddnc8 * dble(idnc)
         osum = 0.
         dc   = 0.
         do k = 1,2000 !MLO - What is this 2000? 
            dc = dc + ddc
            v1 = 0.523599 * dc ** 3
            expvar=-ajlso * v1 * dtlt8
            if (expvar > tinyexp) then
              osum = osum + (dc / dnc) ** (gnuc - 1.) * exp(-dc / dnc) * (1. - exp(expvar))
            else
              osum = osum + (dc / dnc) ** (gnuc - 1.) * exp(-dc / dnc)
            end if
         end do
         dfracc = osum * ddc / (exp(dble(gammln(sngl(gnuc)))) * dnc)
         fracc(idnc,itc,ngr) = sngl(dfracc)
      end do
   end do
   return
end subroutine homfrzcl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     As before, sedimentation is not yet designed to transfer hydrometeor mass between    !
! grids in the case where a nested grid does not reach the top and/or bottom of the model  !
! domain.  Thus, vertical nested grid boundaries should be avoided where sedimentation     !
! occurs.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine mksedim_tab(m1,zm,dzt,pcpfillc,pcpfillr,sfcpcp)

   use micphys, only: &
            nembfall  & ! intent(in)
           ,maxkfall  & ! intent(in)
           ,nhcat     & ! intent(in)
           ,sedtime0  & ! intent(out)
           ,sedtime1  & ! intent(out)
           ,dispemb0  & ! intent(out)
           ,dispemb0i & ! intent(out)
           ,cfvt      & ! intent(in)
           ,emb0      & ! intent(in)
           ,cfmas     & ! intent(in)
           ,cfmasi    & ! intent(in)
           ,pwvt      & ! intent(in)
           ,pwmasi    & ! intent(in)
           ,dispemb1  & ! intent(out)
           ,ch2       & ! intent(out)
           ,emb1      & ! intent(in)
           ,gnu       & ! intent(in)
           ,pwmas       ! intent(in)

   use micro_coms, only : lcat_lhcat
   implicit none

   !----- Arguments: ----------------------------------------------------------------------!
   integer                                    , intent(in)  :: m1
   real, dimension(m1                        ), intent(in)  :: zm,dzt
   real, dimension(m1,maxkfall,nembfall,nhcat), intent(out) :: pcpfillc,pcpfillr
   real, dimension(   maxkfall,nembfall,nhcat), intent(out) :: sfcpcp
   !----- Local Constant: -----------------------------------------------------------------!
   integer, parameter    :: nbin=50
   !----- Local Variables: ----------------------------------------------------------------!
   integer               :: iembs,lcat,lhcat,k,kkf,ibin,kk,jbin
   real                  :: dmbodn,diam0,diam1,fac1,fac3,sumc,sumr,diam,fac2,fac4
   real                  :: disp,ztopnew,zbotnew,fallin,delzsfc,dispemb,dispmax,dispmx
   real, dimension(nbin) :: cbin,rbin,reldisp
   !----- Functions -----------------------------------------------------------------------!
   real, external        :: gammln,gammp
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Because timestep may now be variable in time, define sedtime0 and sedtime1 here   !
   ! as 0.1 seconds and 3000 seconds.  The former is supposed to be less than 0.7 of the   !
   ! shortest timestep on any grid (sqrt(dn0i) never exceeds 0.7) and the latter is the    !
   ! longest timestep expected to ever be used (300 seconds) times a factor of 2 for the   !
   ! maximum inverse of rtgt times a factor of 5 for the largest value of sqrt(dn0i).      !
   !---------------------------------------------------------------------------------------!
   sedtime0 = .1
   sedtime1 = 3000.
   dispmax  = 500.

   !----- Loop over hydrometeor categories ------------------------------------------------!

   do lhcat = 1,nhcat
      lcat = lcat_lhcat(lhcat)

      dispemb0(lhcat) = sedtime0 * cfvt(lhcat)                                             &
                      * (emb0(lcat) * cfmasi(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))
      
      dispemb0i(lhcat) = 1. / dispemb0(lhcat)

      dispemb1(lhcat) = sedtime1 * cfvt(lhcat)                                             &
                      * (emb1(lcat) * cfmasi(lhcat)) ** (pwvt(lhcat) * pwmasi(lhcat))

      !----- Bob (10/24/00):  Limit dispemb1 to a maximum of dispmax ----------------------!
      if (dispemb1(lhcat) > dispmax) dispemb1(lhcat) = dispmax

      ch2(lhcat) = real(nembfall-1) / log10(dispemb1(lhcat) * dispemb0i(lhcat))

      !------------------------------------------------------------------------------------!
      !    Loop over bins, filling them with fractional number, fractional mass, and dis-  !
      ! placement quotient relative to emb.                                                !
      !------------------------------------------------------------------------------------!
      dmbodn = (exp(gammln(gnu(lcat) + pwmas(lhcat)) - gammln(gnu(lcat)))) ** pwmasi(lhcat)
      diam0 = 0.06 * dmbodn
      diam1 = 1.0 * dmbodn
      fac1 = gammp(gnu(lcat),diam0)
      fac3 = gammp(gnu(lcat) + pwmas(lhcat),diam0)
      sumc = 0.
      sumr = 0.

      do jbin = 1,nbin
         diam = diam0 * (diam1 / diam0) ** (real(jbin)/real(nbin))
         fac2 = gammp(gnu(lcat),diam)
         fac4 = gammp(gnu(lcat) + pwmas(lhcat),diam)
         cbin(jbin) = fac2 - fac1
         rbin(jbin) = fac4 - fac3
         fac1 = fac2
         fac3 = fac4
         sumc = sumc + cbin(jbin)
         sumr = sumr + rbin(jbin)
         reldisp(jbin) = diam ** pwvt(lhcat)
      end do

      do jbin = 1,nbin
         cbin(jbin) = cbin(jbin) / sumc
         rbin(jbin) = rbin(jbin) / sumr
      end do

      !----- Loop over displacement distance for size emb. --------------------------------!

      do iembs = 1,nembfall
         dispemb = dispemb0(lhcat) * (dispemb1(lhcat) * dispemb0i(lhcat))                  &
                                  ** (real(iembs-1) / real(nembfall-1))

         !---------------------------------------------------------------------------------!
         !     Zero out concentration and mass fill arrays and surface precip array before !
         ! accumulation.                                                                   !
         !---------------------------------------------------------------------------------!
         do k = 1,m1
            do kkf = 1,maxkfall
               pcpfillc(k,kkf,iembs,lhcat) = 0.
               pcpfillr(k,kkf,iembs,lhcat) = 0.
            end do
            if (k <= maxkfall) sfcpcp(k,iembs,lhcat) = 0.
         end do

         !----- Loop over vertical grid index. --------------------------------------------!
         do k = 2,m1-1
            !----- Bob (10/24/00):  Limit disp distance to (maxkfall-1) levels ------------!
            dispmx = dispmax
            if (k > maxkfall) then
               dispmx = min(dispmx,zm(k-1) - zm(k-maxkfall))
            endif

            !----- Loop over bins ---------------------------------------------------------!

            do ibin = 1,nbin
               disp = max(dispmx,dispemb * reldisp(ibin))

               ztopnew = zm(k)   - disp
               zbotnew = zm(k-1) - disp

               !----- Loop over grid cells that a parcel falls into. ----------------------! 

               parcelloop: do kkf = 1,min(k-1,maxkfall)

                  kk = k + 1 - kkf
                  if (zbotnew > zm(kk)) exit parcelloop

                  if (ztopnew <= zm(kk-1)) then
                     fallin = 0.
                  else
                     fallin = dzt(kk) * (min(zm(kk),ztopnew) - max(zm(kk-1),zbotnew))
                  end if

                  pcpfillc(k,kkf,iembs,lhcat) = pcpfillc(k,kkf,iembs,lhcat)                &
                                              + fallin * cbin(ibin)
                  pcpfillr(k,kkf,iembs,lhcat) = pcpfillr(k,kkf,iembs,lhcat)                &
                                              + fallin * rbin(ibin)
               end do parcelloop

               !----- Compute surface precipitation. --------------------------------------!
               if (zbotnew < 0.) then
                  delzsfc = min(0.,ztopnew) - zbotnew
                  if (k <= maxkfall) then
                     sfcpcp(k,iembs,lhcat) = sfcpcp(k,iembs,lhcat) + delzsfc * rbin(ibin)
                  end if
               end if
            end do
         end do
      end do
   end do

   return
end subroutine mksedim_tab
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine tabmelt()

   use micphys, only: &
           ncat       & ! intent(in)
          ,nhcat      & ! intent(in)
          ,gnu        & ! intent(in)
          ,rmlttab    & ! intent(out)
          ,ninc       & ! intent(in)
          ,enmlttab   & ! intent(out)
          ,ndns       & ! intent(in)
          ,shedtab    & ! intent(out)
          ,cfmas      & ! intent(in)
          ,pwmas      & ! intent(in)
          ,cfvt       & ! intent(in)
          ,pwvt       & ! intent(in)
          ,shapefac   ! ! intent(in)
        
   use micro_coms, only : &
           lcat_lhcat     & ! intent(in)
          ,dmean          & ! intent(in)
          ,vk             ! ! intent(in)

   implicit none

   !----- Local Constants -----------------------------------------------------------------!
   integer, parameter     :: nbins=500
   !----- Local Variables -----------------------------------------------------------------!
   integer                :: lhcat,lcat,ndns1,ibin,inc,iter,idns
   real                   :: dn,gammaa,totfmg,totmass,vtx,fre,totqm,qmgoal,qmnow,totmdqdt
   real                   :: deltat,pliqmass,picemass,critmass
   real, dimension(nbins) :: db,fmg,pmass,binmass,dqdt,q
   !----- Functions -----------------------------------------------------------------------!
   real, external         :: gammln
   !---------------------------------------------------------------------------------------!

   do lhcat = 1,nhcat
      lcat = lcat_lhcat(lhcat)

      dn     = dmean(lcat) / gnu(lcat)
      gammaa = exp(gammln(gnu(lcat)))

      rmlttab(1)           = 0.0
      rmlttab(ninc)        = 1.0
      enmlttab(1,lhcat)    = 0.0
      enmlttab(ninc,lhcat) = 1.0

      ndns1 = 1
      if (lcat == 7) ndns1 = ndns

      do idns = 1,ndns1
         shedtab(1,idns) = 0.0
         shedtab(ninc,idns) = 0.0

         if (ndns1 > 1) dn = 1.e-3 * real(idns) / gnu(lcat)

         totfmg = 0.
         totmass = 0.
         do ibin = 1,nbins
            db(ibin)      = 0.02 * dn * (real(ibin) - 0.5)
            fmg(ibin)     = (db(ibin)/dn) ** (gnu(lcat)-1.) / (dn*gammaa) *                &
                            exp(-db(ibin)/dn)
            totfmg        = totfmg + fmg(ibin)
            q(ibin)       = 0.
            pmass(ibin)   = cfmas(lhcat) * db(ibin) ** pwmas(lhcat)
            binmass(ibin) = pmass(ibin) * fmg(ibin)
            totmass       = totmass + binmass(ibin)
            vtx           = cfvt(lhcat) * db(ibin) ** pwvt(lhcat)
            fre           = (1.0 + 0.229 * sqrt(vtx*db(ibin)/vk)) * shapefac(lhcat)
            dqdt(ibin)    = db(ibin) ** (1. - pwmas(lhcat)) * fre
         end do
         totqm = totmass * 80.

         do inc = 2,ninc-1
            qmgoal = totqm * real(inc-1) / real(ninc-1)
            do iter = 1,2
               qmnow = 0.
               totmdqdt = 0.
               do ibin = 1,nbins
                  if(q(ibin) < 79.9999)then
                     totmdqdt = totmdqdt + binmass(ibin) * dqdt(ibin)
                  endif
                  qmnow = qmnow + q(ibin) * binmass(ibin)
               end do
               deltat = max(0.,(qmgoal - qmnow) / totmdqdt)
               do ibin = 1,nbins
                  q(ibin) = min(80.,q(ibin) + dqdt(ibin) * deltat)
               end do
            end do

            !  For the current inc value (representing total liquid fraction), compute
            !  melted mixing ratio (rmlttab) and number (enmlttab) from totally-melted
            !  bins and compute shedded mixing ratio (shedtab) from partially-melted bins.

            if(idns == 7)then
               rmlttab(inc) = 0.0
               do ibin = 1,nbins
                  if(q(ibin) > 79.9)then
                     rmlttab(inc) = rmlttab(inc) + binmass(ibin)
                  end if
               end do
               rmlttab(inc) = rmlttab(inc) / totmass
            end if

            if(idns == 7 .or. ndns1 == 1)then
               enmlttab(inc,lhcat) = 0.0
               do ibin = 1,nbins
                  if(q(ibin) .gt. 79.9)then
                     enmlttab(inc,lhcat) = enmlttab(inc,lhcat)  &
                          + fmg(ibin)
                  end if
               end do
               enmlttab(inc,lhcat) = enmlttab(inc,lhcat) / totfmg
            end if

            if(lcat == 7)then
               shedtab(inc,idns) = 0.0
               !                  do ibin = kbin,nbins
               do ibin = 1,nbins
                  if(q(ibin) <= 79.9)then
                     pliqmass = pmass(ibin) * q(ibin) / 80.
                     picemass = pmass(ibin) - pliqmass
                     critmass = .268e-3 + .1389 * picemass
                     shedtab(inc,idns) = shedtab(inc,idns)  &
                          + max(0.0, pliqmass - critmass) * fmg(ibin)
                  end if
               end do
               shedtab(inc,idns) = shedtab(inc,idns) / totmass
            endif

         enddo
      enddo
   enddo
   return
end subroutine tabmelt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine constructs some useful collection tables.                             !
!------------------------------------------------------------------------------------------!
subroutine mkcoltb

   use micphys, only: &
           nhcat      & ! intent(in)
          ,gnu        & ! intent(in)
          ,pwmas      & ! intent(in)
          ,emb0       & ! intent(in)
          ,cfmas      & ! intent(in)
          ,emb1       & ! intent(in)
          ,ipairc     & ! intent(in)
          ,ipairr     & ! intent(in)
          ,pwvt       & ! intent(in)
          ,nembc      & ! intent(in)
          ,cfvt       & ! intent(in)
          ,coltabc    & ! intent(out)
          ,coltabr    ! ! intent(out)
   use micro_coms, only : lcat_lhcat
   implicit none

   !----- Local Constant: -----------------------------------------------------------------!
   integer, parameter :: ndx=20
   !----- Local Variables: ----------------------------------------------------------------!
   integer              :: ihx,ix,ihy,iy,iemby,iembx,idx
   real                 :: gxm,dnminx,dnmaxx,dxlo,dxhi,gyn,gyn1,gyn2,gynp,gynp1,gynp2,gym
   real                 :: dnminy,dnmaxy,dny,vny,dnx,ans
   real, dimension(ndx) :: dx,fx,gx
   !----- Function ------------------------------------------------------------------------!
   real, external       :: gammln,xj
   !---------------------------------------------------------------------------------------!

   do ihx = 1,nhcat
      ix = lcat_lhcat(ihx)

      gxm    = exp(gammln(gnu(ix)) - gammln(gnu(ix) + pwmas(ihx)))
      dnminx = ((emb0(ix) / cfmas(ihx)) * gxm) ** (1. / pwmas(ihx))
      dnmaxx = ((emb1(ix) / cfmas(ihx)) * gxm) ** (1. / pwmas(ihx))
      dxlo   = .01 * dnminx
      dxhi   = 10. * dnmaxx

      do ihy = 1,nhcat

         iy = lcat_lhcat(ihy)

         if (ipairc(ihx,ihy) > 0 .or. ipairr(ihx,ihy) > 0) then
            gyn = exp(gammln(gnu(iy)))
            gyn1 = exp(gammln(gnu(iy) + 1.)) / gyn
            gyn2 = exp(gammln(gnu(iy) + 2.)) / gyn
            gynp = exp(gammln(gnu(iy) + pwvt(ihy))) / gyn
            gynp1 = exp(gammln(gnu(iy) + pwvt(ihy) + 1.)) / gyn
            gynp2 = exp(gammln(gnu(iy) + pwvt(ihy) + 2.)) / gyn

            gym = exp(gammln(gnu(iy)) - gammln(gnu(iy) + pwmas(ihy)))
            dnminy = ((emb0(iy) / cfmas(ihy)) * gym) ** (1. /pwmas(ihy))
            dnmaxy = ((emb1(iy) / cfmas(ihy)) * gym) ** (1. /pwmas(ihy))

            do iemby = 1,nembc
               dny = dnminy * (dnmaxy / dnminy) ** (real(iemby-1)/real(nembc-1))
               vny = cfvt(ihy) * dny ** pwvt(ihy)
               do iembx = 1,nembc

                  dnx = dnminx * (dnmaxx / dnminx) ** (real(iembx-1)/ real(nembc-1))
                  do idx = 1,ndx
                     dx(idx) = dxlo * (dxhi / dxlo)** (real(idx-1) / real(ndx-1))
                     fx(idx) = xj(dx(idx),cfvt(ihx),pwvt(ihx),cfvt(ihy),pwvt(ihy),vny,dnx  &
                                 ,dny,gnu(ix),gnu(iy),gyn1,gyn2,gynp,gynp1,gynp2)
                     gx(idx) = fx(idx) * cfmas(ihx)* dx(idx) ** pwmas(ihx)

                  end do
                  if (ipairc(ihx,ihy) > 0) then
                     call avint(dx,fx,ndx,dxlo,dxhi,ans)
                     coltabc(iembx,iemby,ipairc(ihx,ihy))= -log10(max(1.e-30,ans))
                  end if
                  if (ipairr(ihx,ihy) .gt. 0) then
                     call avint(dx,gx,ndx,dxlo,dxhi,ans)
                     coltabr(iembx,iemby,ipairr(ihx,ihy))= -log10(max(1.e-30,ans))
                  end if
               end do
            end do
         end if
      end do
   end do
   return
end subroutine mkcoltb






!==========================================================================================!
!==========================================================================================!
real function xj(dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu,gyn1,gyn2,gynp,gynp1,gynp2)

   implicit none

   !----- Arguments: ----------------------------------------------------------------------!
   real, intent(in) :: dx,cvx,pvx,cvy,pvy,vny,dnx,dny,xnu,ynu,gyn1,gyn2,gynp,gynp1,gynp2
   !----- Local Variables: ----------------------------------------------------------------!
   real(kind=8) :: dx8,cvx8,pvx8,cvy8,pvy8,vny8,dnx8,dny8,xnu8,ynu8,gyn18,gyn28,gynp8
   real(kind=8) :: gynp18,gynp28
   real(kind=8) :: dnxi8,rdx8,vx8,xj8
   real         :: dxy,ynup
   real(kind=8) :: glxnu8,glynu8,gpynu8,gqynu8,gpynua18,gqynua18,gpynua28,gqynua28
   real(kind=8) :: gpynup8,gqynup8,gpynupa18,gqynupa18,gpynupa28,gqynupa28
   !----- External functions --------------------------------------------------------------!
   real, external :: gammln,gammp,gammq
   !---------------------------------------------------------------------------------------!

   dx8     = dble(dx   )
   pvx8    = dble(pvx  )
   cvy8    = dble(cvy  )
   pvy8    = dble(pvy  )
   vny8    = dble(vny  )
   dnx8    = dble(dnx  )
   dny8    = dble(dny  )
   xnu8    = dble(xnu  ) 
   ynu8    = dble(ynu  ) 
   gyn18   = dble(gyn1 ) 
   gyn28   = dble(gyn2 ) 
   gynp8   = dble(gynp )
   gynp18  = dble(gynp1)
   gynp28  = dble(gynp2)
           
   dnxi8 = 1. / dnx8
   rdx8  = dx8 * dnxi8
   vx8   = cvx8 * dx8 ** pvx8
   dxy   = (sngl(vx8) / cvy) ** (1. / pvy) / dny
   ynup  = ynu + pvy 
   

   if (rdx8 < 38.) then
      glxnu8    = dble(gammln(xnu))
      glynu8    = dble(gammln(ynu))
      gpynu8    = dble(gammp(ynu,dxy))
      gqynu8    = dble(gammq(ynu,dxy))
      gpynua18  = dble(gammp(ynu+1,dxy))
      gqynua18  = dble(gammq(ynu+1,dxy))
      gpynua28  = dble(gammp(ynu+2,dxy))
      gqynua28  = dble(gammq(ynu+2,dxy))
      gpynup8   = dble(gammp(ynu,dxy))
      gqynup8   = dble(gammq(ynu,dxy))
      gpynupa18 = dble(gammp(ynup+1,dxy))
      gqynupa18 = dble(gammq(ynup+1,dxy))
      gpynupa28 = dble(gammp(ynup+2,dxy))
      gqynupa28 = dble(gammq(ynup+2,dxy))
     
      
      xj8 = exp(-rdx8-glxnu8-glynu8) * rdx8**(xnu8-1.) * dnxi8                             &
          * (vx8 * ( dx8*dx8*(gpynu8-gqynu8))                                              &
                   + 2.*dx8*dny8*gyn18*(gpynua18-gqynua18)                                 &
                   + dny8*dny8*gyn28*(gpynua28-gqynua28)                                   &
                   - vny8 * (dx8*dx8*gynp8*(gpynup8-gqynup8)                               &
                            + 2.*dx8*dny8*gynp18*(gpynupa18-gqynupa18)                     &
                            + dny8*dny8*gynp28*(gpynupa28-gqynupa28)))
      !----- Making sure underflow won't happen. If overflow happens, let it complain... --!
      if (abs(xj8) < 1.d-30)  xj8 = sign(1.d-30,xj8)
   else
      xj8 = 0.
   end if
   
   xj = sngl(xj8)
   return
end function xj
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine creates the (auto collection?) table. It works in CGS units...        !
!------------------------------------------------------------------------------------------!
subroutine make_autotab()

   use micphys, only: &
           nd1cc      & ! intent(in)
          ,nd1cr      & ! intent(in)
          ,nr2cr      & ! intent(in)
          ,nr2rr      & ! intent(in)
          ,gnu        & ! intent(in)
          ,nd2cr      & ! intent(in)
          ,nd2rr      & ! intent(in)
          ,d1min      & ! intent(out)
          ,d1max      & ! intent(out)
          ,d1ecc      & ! intent(out)
          ,d1ecr      & ! intent(out)
          ,r2min      & ! intent(out)
          ,r2max      & ! intent(out)
          ,r2ecr      & ! intent(out)
          ,r2err      & ! intent(out)
          ,d2min      & ! intent(out)
          ,d2max      & ! intent(out)
          ,r1tabcc    & ! intent(out)
          ,c1tabcc    & ! intent(out)
          ,c2tabcc    & ! intent(out)
          ,r1tabcr    & ! intent(out)
          ,c1tabcr    & ! intent(out)
          ,c2tabrr      ! intent(out)

   implicit none

   !----- Local constants -----------------------------------------------------------------!
   integer, parameter :: ibins=36
   integer, parameter :: ithresh=15
   !----- Local Variables -----------------------------------------------------------------!
   integer                        :: i,k,id1cc,id1cr,ir2cr,id2cr,ir2rr,id2rr
   real                           :: r2,en1,en2,en1i,en1i2,d1,r1,d2minx,d2ecr,d2
   real                           :: sum1,sum10,sun10,sun1,sun20,sum20,sun2,sum2,d2err
   real, dimension(ibins+1)       :: x,diam
   real, dimension(ibins)         :: ank0,amk0,ank,amk,ank1,amk1,ank2,amk2
   real, dimension(ibins,ibins,3) :: akbarx
   real, dimension(ibins,ibins)   :: akbar
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Read in mass grid x(k+1)=2*x(k), diameters (diam) and collection kernel kbar          !
   ! GNF: kernels for ice with cloud?                                                      !
   !---------------------------------------------------------------------------------------!
   call micro_data(x,diam,akbar,ibins)

   call azero (ibins*ibins*3,akbarx)
   do i=1,ibins
      do k=1,ibins
         if(i > ithresh .and. k > ithresh) then
            akbarx(i,k,3) = akbar(i,k)
         elseif(i <= ithresh .and. k <= ithresh) then
            akbarx(i,k,1) = akbar(i,k)
         else
            akbarx(i,k,2) = akbar(i,k)
         end if
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !     d1min and d1max are equivalent to dmb0 and dmb1, but may have different values.   !
   !---------------------------------------------------------------------------------------!
   d1min = 4.e-4
   d1max = 50.e-4
   d1ecc = log10 (d1max / d1min) / float(nd1cc-1)
   d1ecr = log10 (d1max / d1min) / float(nd1cr-1)

   r2min = .01e-6
   r2max = 20.e-6
   r2ecr = log10 (r2max / r2min) / float(nr2cr-1)
   r2err = log10 (r2max / r2min) / float(nr2rr-1)

   d2min = 1.e-2
   d2max = 1.

   !---------------------------------------------------------------------------------------!
   ! Start 1 cc loop for dm1, dn1, and dn2.                                                !
   !---------------------------------------------------------------------------------------!
   r2 = .01e-06
   en1 = 100.
   en2 = 1.e-6
   en1i = 1. / en1
   en1i2 = en1i ** 2

   do id1cc = 1,nd1cc

      d1 = d1min + (d1max - d1min) * float(id1cc-1) / float(nd1cc-1)
      r1 = en1 * .5236 * d1 ** 3

      call initg(r1,r2,en1,en2,gnu(1),gnu(2),diam,x,amk0,ank0,ank1,amk1,ank2,amk2,ibins    &
                ,ithresh)
      call sumn(ank0,amk0,1,ithresh,ibins,sun10,sum10)
      call sumn(ank0,amk0,ithresh+1,ibins,ibins,sun20,sum20)
      call sxy(x,amk0,ank0,amk,ank,akbarx(1:ibins,1:ibins,1))
      call sumn(ank,amk,1,ithresh,ibins,sun1,sum1)
      call sumn(ank,amk,ithresh+1,ibins,ibins,sun2,sum2)

      r1tabcc(id1cc) = max(0.,(sum10-sum1) * en1i2)
      c1tabcc(id1cc) = max(0.,(sun10-sun1) * en1i2)
      c2tabcc(id1cc) = max(0.,(sun2-sun20) * en1i2)

   end do

   !---------------------------------------------------------------------------------------!
   ! Start 3 cr loops for dm1 and dn1.                                                     !
   !---------------------------------------------------------------------------------------!
   do id1cr = 1,nd1cr
      d1 = d1min * 10. ** (d1ecr * float(id1cr-1))
      r1 = en1 * .5236 * d1 ** 3

      do ir2cr = 1,nr2cr
         r2 = r2min * 10. ** (r2ecr * float(ir2cr-1))
         d2minx = max(d2min,(r2 / (.1 * .5236)) ** .333333)
         d2ecr = alog10(d2max / d2minx) / float(nd2cr-1)

         do id2cr = 1,nd2cr
            d2 = d2minx * 10. ** (d2ecr * float(id2cr-1))
            en2 = r2 / (.5236 * d2 ** 3)

            call initg(r1,r2,en1,en2,gnu(1),gnu(2),diam,x,amk0,ank0,ank1,amk1,ank2,amk2    &
                      ,ibins,ithresh)
            call sumn(ank0,amk0,1,ithresh,ibins,sun10,sum10)
            call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,2))
            call sumn(ank,amk,1,ithresh,ibins,sun1,sum1)

            r1tabcr(id1cr,ir2cr,id2cr) = alog10(max(1.e-20,(sum10-sum1) * en1i))
            c1tabcr(id1cr,ir2cr,id2cr) = alog10(max(1.e-20,(sun10-sun1) * en1i))

         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   ! Start 2 rr loops for dn2.                                                             !
   !---------------------------------------------------------------------------------------!
   d1 = 4.e-4
   r1 = en1 * .5236 * d1 ** 3

   do ir2rr = 1,nr2rr
      r2 = r2min * 10. ** (r2err * float(ir2rr-1))
      d2minx = max(d2min,(r2 / (.1 * .5236)) ** .333333)
      d2err = alog10(d2max / d2minx) / float(nd2rr-1)

      do id2rr = 1,nd2rr
         d2 = d2minx * 10. ** (d2err * float(id2rr-1))
         en2 = r2 / (.5236 * d2 ** 3)

         call initg(r1,r2,en1,en2,gnu(1),gnu(2),diam,x,amk0,ank0,ank1,amk1,ank2,amk2,ibins &
                   ,ithresh)
         call sumn(ank0,amk0,ithresh+1,ibins,ibins,sun20,sum20)
         call sxy(x,amk0,ank0,amk,ank,akbarx(1,1,3))
         call sumn(ank,amk,ithresh+1,ibins,ibins,sun2,sum2)

         c2tabrr(ir2rr,id2rr) = alog10(max(1.e-25,sun20-sun2))

      end do
   end do

   return
end subroutine make_autotab
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine was adapted to double precision to avoid underflow and overflow.      !
!------------------------------------------------------------------------------------------!
subroutine sxy(x,amkd,ankd,amk,ank,akbar)
   implicit none
  
   !----- Local Constants -----------------------------------------------------------------!
   integer     , parameter :: ibins=36
   real(kind=8), parameter :: ap=1.062500000
   !-----  Arguments: ---------------------------------------------------------------------!
   real, dimension(ibins+1)    , intent(in)  :: x
   real, dimension(ibins)      , intent(in)  :: amkd, ankd
   real, dimension(ibins)      , intent(out) :: amk, ank
   real, dimension(ibins,ibins), intent(in)  :: akbar
   !-----  Local Variables ----------------------------------------------------------------!
   integer                              :: i,ik,k,l
   real(kind=8), dimension(ibins+1)     :: x8
   real(kind=8), dimension(ibins)       :: xave,am2,am3,am4,psi,f,amkd8,ankd8
   real(kind=8), dimension(ibins,ibins) :: akbar8
   real(kind=8)                         :: dm,dn,sm1,sm2,sm3,sm4,sm5
   real(kind=8)                         :: sn1,sn2,sn3,sn4,dm4,dm2
   !---------------------------------------------------------------------------------------!

   x8     = dble(x)
   amkd8  = dble(amkd)
   ankd8  = dble(ankd)
   akbar8 = dble(akbar)

   do l=1,ibins
      if(ankd8(l) > 0.)then
         xave(l)=amkd8(l)/ankd8(l)
      else
         xave(l)=0.
      endif
   enddo

   do k=1,ibins
      !------------------------------------------------------------------------------------!
      !     Calculation of the 2nd, 3rd, and 4th moments of the mass distribution based on !
      ! equation 8 in reference.                                                           !
      !------------------------------------------------------------------------------------!
      am2(k) = ap*xave(k)*amkd8(k)
      am3(k) = ap*ap*xave(k)*am2(k)
      am4(k) = ap*ap*ap*xave(k)*am3(k)

      !------------------------------------------------------------------------------------!
      !     These functions come out of the linear approximations used to integrate over   !
      ! partial bins.  they are defined:                                                   !
      !      psi(k) = nk(k+1)                                                              !
      !        f(k) = nk(k)                                                                !
      ! where nk is the distribution function.  see equation 13 in reference.              !
      !------------------------------------------------------------------------------------!
      psi(k) = 2./x8(k)*(amkd8(k)/x8(k)-ankd8(k))
      f(k)   = 2./x8(k)*(2.*ankd8(k)-amkd8(k)/x8(k)) 

      !----- Zeroing the tendencies on the moments. ---------------------------------------!
      sm1=0.
      sm2=0.
      sm3=0.
      sm4=0.
      sm5=0.
      sn1=0.
      sn2=0.
      sn3=0.
      sn4=0.

      !----- Calculation of tendencies on moments -----------------------------------------!
      do i=k,ibins
         dm = akbar8(i,k)*(  am2(k)*ankd8(i)+amkd8(k)*amkd8(i))
         dn = akbar8(i,k)*(ankd8(k)*amkd8(i)+amkd8(k)*ankd8(i))

         sm5 = sm5+dm
         sn4 = sn4+dn
      end do

      if(k > 1)then
         sm3 = akbar8(k-1,k-1)*(am2(k-1)*ankd8(k-1)+amkd8(k-1)*amkd8(k-1))
         sn2 = akbar8(k-1,k-1)*ankd8(k-1)*amkd8(k-1)
         dn  = sn2
         dm  = sm3
      end if

      do i=1,k-1
         dm4=akbar8(k,i)*(ankd8(k)*am2(i)+amkd8(k)*amkd8(i))
         sm4=sm4+dm4

         if(xave(k) >= x8(k)) then
            dm2=akbar8(k,i)*(4.*x8(k)*x8(k)*psi(k)*amkd8(i)                                &
                 +0.5*x8(k)*(4.*psi(k)+f(k))*am2(i)                                        &
                 -(psi(k)-f(k))*am3(i)                                                     &
                 -0.5/x8(k)*(psi(k)-f(k)*am4(i)))
            sm2=sm2+dm2
            dn=akbar8(k,i)*(2.*x8(k)*psi(k)*amkd8(i)                                       &
                 +0.5*(f(k)*am2(i))                                                        &
                 -0.5/x8(k)*((psi(k)-f(k))*am3(i)))
            sn3=sn3+dn
         end if
      end do

      do i=1,k-2
         ik=k-1
         if(xave(ik) >= x8(ik))then

            dm=akbar8(ik,i)*(4.*x8(ik)*x8(ik)*psi(ik)*amkd8(i)                             &
                 +0.5*x8(ik)*(4.*psi(ik)+f(ik))*am2(i)                                     &
                 -(psi(ik)-f(ik))*am3(i)                                                   &
                 -0.5/x8(ik)*(psi(ik)-f(ik))*am4(i))
            sm1=sm1+dm
            dn=akbar8(ik,i)*(2.*x8(ik)*psi(ik)*amkd8(i)                                    &
                 +0.5*(f(ik)*am2(i))                                                       &
                 -0.5/x8(ik)*(psi(ik)-f(ik))*am3(i))
            sn1=sn1+dn

         endif
      enddo

      amk(k)=sngl(amkd8(k)+sm1-sm2+sm3+sm4-sm5)
      ank(k)=sngl(ankd8(k)+sn1+sn2-sn3-sn4)

   enddo
   return
end subroutine sxy

!******************************************************************************

subroutine micro_data(x,diam,akbar,ibins)
  use rconstants, only : pio6,pio6i,onethird
  implicit none

  !Arguments:
  integer, intent(in)                       :: ibins
  real, dimension(ibins+1), intent(out)     :: x, diam
  real, dimension(ibins,ibins), intent(out) :: akbar

  !Local Variables:
  integer :: l,i,j,kount,n
  real :: p,ap
  real, dimension(36,36) :: aabar

  data (aabar( 1,n),n=1, 1) /-.47757E-01  /
  data (aabar( 2,n),n=1, 2) /-.26460E+00,-.47965E-01 /
  data (aabar( 3,n),n=1, 3) /-.82258E+00,-.26760E+00,-.20453E-01 /
  data (aabar( 4,n),n=1, 4) /-.19050E+01,-.82072E+00,-.11992E+00, .78909E-01 /
  data (aabar( 5,n),n=1, 5) /-.39171E+01,-.18915E+01,-.33270E+00, .41936E+00  &
       ,.34801E+00 /
  data (aabar( 6,n),n=1, 6) /-.76415E+01,-.38808E+01,-.73737E+00, .14121E+01  &
       ,.18851E+01,.99793E+00 /
  data (aabar( 7,n),n=1, 7) /-.14595E+02,-.75638E+01,-.14861E+01, .33598E+01  &
       ,.61219E+01, .54314E+01, .24751E+01 /
  data (aabar( 8,n),n=1, 8) /-.27720E+02,-.14442E+02,-.28741E+01, .69895E+01  &
       ,.14394E+02, .17479E+02, .13500E+02, .57110E+01 /
  data (aabar( 9,n),n=1, 9) /-.52737E+02,-.27428E+02,-.54729E+01, .13703E+02  &
       ,.29792E+02, .40971E+02, .43267E+02, .31185E+02, .12630E+02 /
  data (aabar(10,n),n=1,10) /-.10083E+03,-.52188E+02,-.10391E+02, .26218E+02  &
       ,.58283E+02, .84686E+02, .10128E+03, .99726E+02, .69014E+02, .27176E+02 /
  data (aabar(11,n),n=1,11) /-.19396E+03,-.99799E+02,-.19790E+02, .49801E+02  &
       ,.11143E+03, .16558E+03, .20922E+03, .23326E+03, .22039E+03, .14858E+03  &
       ,.57396E+02 /
  data (aabar(12,n),n=1,12) /-.37536E+03,-.19200E+03,-.37896E+02, .94692E+02  &
       ,.21165E+03, .31650E+03, .40896E+03, .48169E+03, .51524E+03, .47402E+03  &
       ,.31389E+03, .11962E+03 /
  data (aabar(13,n),n=1,13) /-.73047E+03,-.37164E+03,-.73015E+02, .18089E+03  &
       ,.40253E+03, .60115E+03, .78166E+03, .94143E+03, .10638E+04, .11078E+04  &
       ,.10008E+04, .65436E+03, .24691E+03 /
  data (aabar(14,n),n=1,14) /-.14285E+04,-.72333E+03,-.14152E+03, .34764E+03  &
       ,.76925E+03, .11434E+04, .14846E+04, .17993E+04, .20789E+04, .22870E+04  &
       ,.23385E+04, .20854E+04, .13509E+04, .50600E+03 /
  data (aabar(15,n),n=1,15) /-.41365E+04,-.20869E+04,-.40697E+03, .99310E+03  &
       ,.21878E+04, .32394E+04, .41995E+04, .51084E+04, .59888E+04, .68297E+04  &
       ,.75528E+04, .79583E+04, .76785E+04, .62489E+04, .76776E+03 /
  data (aabar(16,n),n=1,16) / .63760E+04, .64739E+04, .65970E+04, .67516E+04  &
       ,.69451E+04, .71861E+04, .74835E+04, .78448E+04, .82709E+04, .87453E+04  &
       ,.92111E+04, .95276E+04, .94079E+04, .83797E+04, .26045E+04, .89777E+03 /
  data (aabar(17,n),n=1,17) / .62974E+04, .63746E+04, .64717E+04, .65934E+04  &
       ,.67457E+04, .69355E+04, .71702E+04, .74571E+04, .78005E+04, .81957E+04  &
       ,.86163E+04, .89879E+04, .91399E+04, .87394E+04, .46530E+04, .26045E+04  &
       ,.89777E+03 /
  data (aabar(18,n),n=1,18) / .62353E+04, .62963E+04, .63729E+04, .64689E+04  &
       ,.65889E+04, .67383E+04, .69233E+04, .71502E+04, .74238E+04, .77446E+04  &
       ,.81009E+04, .84538E+04, .87067E+04, .86514E+04, .59471E+04, .46530E+04  &
       ,.26045E+04, .89777E+03 /
  data (aabar(19,n),n=1,19) / .61862E+04, .62344E+04, .62949E+04, .63707E+04  &
       ,.64653E+04, .65831E+04, .67290E+04, .69080E+04, .71250E+04, .73819E+04  &
       ,.76742E+04, .79815E+04, .82491E+04, .83524E+04, .66125E+04, .59471E+04  &
       ,.46530E+04, .26045E+04, .89777E+03 /
  data (aabar(20,n),n=1,20) / .61474E+04, .61855E+04, .62334E+04, .62932E+04  &
       ,.63679E+04, .64608E+04, .65759E+04, .67172E+04, .68887E+04, .70932E+04  &
       ,.73291E+04, .75856E+04, .78311E+04, .79911E+04, .68735E+04, .66125E+04  &
       ,.59471E+04, .46530E+04, .26045E+04, .89777E+03 /
  data (aabar(21,n),n=1,21) / .61166E+04, .61468E+04, .61847E+04, .62320E+04  &
       ,.62910E+04, .63644E+04, .64552E+04, .65668E+04, .67023E+04, .68644E+04  &
       ,.70531E+04, .72625E+04, .74738E+04, .76415E+04, .69140E+04, .68735E+04  &
       ,.66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
  data (aabar(22,n),n=1,22) / .60923E+04, .61162E+04, .61462E+04, .61836E+04  &
       ,.62303E+04, .62883E+04, .63600E+04, .64481E+04, .65553E+04, .66836E+04  &
       ,.68338E+04, .70027E+04, .71786E+04, .73330E+04, .68498E+04, .69140E+04  &
       ,.68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04, .89777E+03 /
  data (aabar(23,n),n=1,23) / .60730E+04, .60919E+04, .61157E+04, .61453E+04  &
       ,.61823E+04, .62281E+04, .62848E+04, .63545E+04, .64392E+04, .65408E+04  &
       ,.66601E+04, .67953E+04, .69391E+04, .70729E+04, .67447E+04, .68498E+04  &
       ,.69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04, .26045E+04  &
       ,.89777E+03 /
  data (aabar(24,n),n=1,24) / .60577E+04, .60727E+04, .60915E+04, .61150E+04  &
       ,.61443E+04, .61806E+04, .62254E+04, .62805E+04, .63475E+04, .64279E+04  &
       ,.65225E+04, .66304E+04, .67467E+04, .68590E+04, .66311E+04, .67447E+04  &
       ,.68498E+04, .69140E+04, .68735E+04, .66125E+04, .59471E+04, .46530E+04  &
       ,.26045E+04, .89777E+03 /
  data (aabar(25,n),n=1,25) / .77967E+04, .78122E+04, .78316E+04, .78560E+04  &
       ,.78863E+04, .79242E+04, .79713E+04, .80294E+04, .81008E+04, .81878E+04  &
       ,.82924E+04, .84158E+04, .85571E+04, .87104E+04, .86265E+04, .88325E+04  &
       ,.90719E+04, .93363E+04, .95996E+04, .98007E+04, .98157E+04, .94274E+04  &
       ,.83361E+04, .63023E+04, .57988E+03 /
  data (aabar(26,n),n=1,26) / .69349E+04, .69458E+04, .69595E+04, .69766E+04  &
       ,.69979E+04, .70244E+04, .70573E+04, .70978E+04, .71473E+04, .72072E+04  &
       ,.72788E+04, .73623E+04, .74565E+04, .75566E+04, .74715E+04, .76064E+04  &
       ,.77647E+04, .79435E+04, .81311E+04, .82983E+04, .83827E+04, .82640E+04  &
       ,.77406E+04, .65488E+04, .15807E+04, .51662E+03 /
  data (aabar(27,n),n=1,27) / .61704E+04, .61781E+04, .61877E+04, .61997E+04  &
       ,.62147E+04, .62333E+04, .62562E+04, .62843E+04, .63186E+04, .63598E+04  &
       ,.64086E+04, .64648E+04, .65271E+04, .65912E+04, .65100E+04, .65961E+04  &
       ,.66976E+04, .68135E+04, .69382E+04, .70574E+04, .71390E+04, .71194E+04  &
       ,.68816E+04, .62379E+04, .26526E+04, .14083E+04, .46025E+03 /
  data (aabar(28,n),n=1,28) / .54916E+04, .54971E+04, .55038E+04, .55123E+04  &
       ,.55228E+04, .55357E+04, .55517E+04, .55712E+04, .55949E+04, .56232E+04  &
       ,.56562E+04, .56938E+04, .57345E+04, .57747E+04, .57001E+04, .57533E+04  &
       ,.58161E+04, .58880E+04, .59661E+04, .60426E+04, .61008E+04, .61062E+04  &
       ,.59940E+04, .56500E+04, .31742E+04, .23632E+04, .12546E+04, .41004E+03 /
  data (aabar(29,n),n=1,29) / .48886E+04, .48924E+04, .48971E+04, .49031E+04  &
       ,.49104E+04, .49195E+04, .49306E+04, .49441E+04, .49604E+04, .49797E+04  &
       ,.50020E+04, .50269E+04, .50530E+04, .50774E+04, .50108E+04, .50422E+04  &
       ,.50792E+04, .51212E+04, .51667E+04, .52111E+04, .52447E+04, .52486E+04  &
       ,.51861E+04, .49913E+04, .32935E+04, .28279E+04, .21054E+04, .11177E+04  &
       ,.36530E+03 /
  data (aabar(30,n),n=1,30) / .43524E+04, .43551E+04, .43585E+04, .43626E+04  &
       ,.43678E+04, .43741E+04, .43818E+04, .43912E+04, .44024E+04, .44155E+04  &
       ,.44304E+04, .44467E+04, .44631E+04, .44771E+04, .44188E+04, .44361E+04  &
       ,.44561E+04, .44786E+04, .45022E+04, .45241E+04, .45384E+04, .45339E+04  &
       ,.44893E+04, .43663E+04, .31847E+04, .29342E+04, .25193E+04, .18757E+04  &
       ,.99579E+03, .32545E+03 /
  data (aabar(31,n),n=1,31) / .38756E+04, .38775E+04, .38799E+04, .38828E+04  &
       ,.38864E+04, .38908E+04, .38961E+04, .39026E+04, .39102E+04, .39191E+04  &
       ,.39290E+04, .39395E+04, .39494E+04, .39568E+04, .39066E+04, .39149E+04  &
       ,.39241E+04, .39340E+04, .39435E+04, .39507E+04, .39516E+04, .39392E+04  &
       ,.39006E+04, .38129E+04, .29707E+04, .28372E+04, .26141E+04, .22445E+04  &
       ,.16710E+04, .88715E+03, .28994E+03 /
  data (aabar(32,n),n=1,32) / .30106E+04, .30118E+04, .30132E+04, .30149E+04  &
       ,.30171E+04, .30197E+04, .30229E+04, .30266E+04, .30309E+04, .30357E+04  &
       ,.30408E+04, .30456E+04, .30491E+04, .30494E+04, .30032E+04, .30013E+04  &
       ,.29981E+04, .29929E+04, .29844E+04, .29706E+04, .29480E+04, .29111E+04  &
       ,.28504E+04, .27503E+04, .20956E+04, .19717E+04, .17926E+04, .15245E+04  &
       ,.11243E+04, .55892E+03, .22135E+03, .00000E+00 /
  data (aabar(33,n),n=1,33) / .23888E+04, .23895E+04, .23903E+04, .23914E+04  &
       ,.23927E+04, .23943E+04, .23962E+04, .23983E+04, .24007E+04, .24033E+04  &
       ,.24057E+04, .24075E+04, .24077E+04, .24049E+04, .23645E+04, .23582E+04  &
       ,.23497E+04, .23382E+04, .23225E+04, .23007E+04, .22699E+04, .22258E+04  &
       ,.21613E+04, .20655E+04, .15572E+04, .14495E+04, .13057E+04, .11057E+04  &
       ,.82055E+03, .41732E+03, .17636E+03, .00000E+00, .00000E+00 /
  data (aabar(34,n),n=1,34) / .18955E+04, .18959E+04, .18964E+04, .18971E+04  &
       ,.18979E+04, .18988E+04, .18999E+04, .19011E+04, .19024E+04, .19036E+04  &
       ,.19045E+04, .19047E+04, .19033E+04, .18990E+04, .18647E+04, .18567E+04  &
       ,.18462E+04, .18326E+04, .18145E+04, .17905E+04, .17581E+04, .17138E+04  &
       ,.16525E+04, .15662E+04, .11695E+04, .10771E+04, .95995E+03, .80547E+03  &
       ,.59516E+03, .30433E+03, .13287E+03, .00000E+00, .00000E+00, .00000E+00 /
  data (aabar(35,n),n=1,35) / .15041E+04, .15044E+04, .15047E+04, .15051E+04  &
       ,.15056E+04, .15061E+04, .15067E+04, .15074E+04, .15080E+04, .15084E+04  &
       ,.15085E+04, .15079E+04, .15058E+04, .15011E+04, .14725E+04, .14642E+04  &
       ,.14536E+04, .14399E+04, .14221E+04, .13988E+04, .13682E+04, .13274E+04  &
       ,.12725E+04, .11976E+04, .88682E+03, .80897E+03, .71340E+03, .59221E+03  &
       ,.43360E+03, .22074E+03, .97241E+02, .00000E+00, .00000E+00, .00000E+00  &
       ,.00000E+00 /
  data (aabar(36,n),n=1,36) / .11936E+04, .11938E+04, .11940E+04, .11942E+04  &
       ,.11945E+04, .11948E+04, .11951E+04, .11954E+04, .11957E+04, .11957E+04  &
       ,.11954E+04, .11943E+04, .11921E+04, .11876E+04, .11640E+04, .11563E+04  &
       ,.11464E+04, .11337E+04, .11174E+04, .10963E+04, .10689E+04, .10330E+04  &
       ,.98554E+03, .92214E+03, .67808E+03, .61344E+03, .53582E+03, .44015E+03  &
       ,.31885E+03, .16089E+03, .70536E+02, .00000E+00, .00000E+00, .00000E+00  &
       ,.00000E+00, .00000E+00 /



  ! calculating the mass categories (c.g.s) with the lowest diameter
  ! of 3.125 microns and mass doubling every bin
  !
  diam(1)=1.5625*2.e-04
  x(1)=pio6*diam(1)**3.
  do l=2,ibins+1
     x(l)=2.*x(l-1)
     diam(l)=(pio6i*x(l))**onethird
  enddo

  l=1
  p=2.0**(1/l)
  ap=0.5+(p+1.0)*(p+1.0)/(8.0*p)

  !
  ! long's collection kernel as calculated by myself (1-36) with (1-36)
  ! weighted for (x+y)
  !
  kount=0
  do i=1,ibins
     do j=1,i
        kount=kount+1
        akbar(i,j)=aabar(i,j)
        if(akbar(i,j).lt.0.) akbar(i,j)=0.
        akbar(36,i)=0.
     enddo
  enddo

  do j=1,ibins
     do i=1,j
        akbar(i,j)=akbar(j,i)
     enddo
  enddo

  return
end subroutine micro_data

!******************************************************************************

subroutine initg(r1,r2,n1,n2,gnu1,gnu2,diam,x,amk,ank  &
     ,ank1,amk1,ank2,amk2,ibins,ithresh)
  use rconstants, only: pi1,pio6i,onethird
  implicit none

  ! Arguments:
  integer, intent(in) :: ibins, ithresh
  real, intent(in) :: r1,r2,n1,n2,gnu1,gnu2
  real, dimension(ibins+1), intent(in) :: x,diam
  real, dimension(ibins), intent(out) :: amk,ank,amk1,ank1,amk2,ank2

  ! Local Variables:
  integer :: i
  real :: dn1,dn2,trunc,trunc1,fac1,fac2 &
       ,gamp,dmean,sum,sumn,xntot,xr3
  real :: gammln,gammp
  
  ! * *
  ! * Initial double gamma distribution: n(D) = n1(D) + n2(D)
  ! * *

  ! * *
  ! * gamma spectrum
  ! * *


  gamp = exp(gammln(gnu1))
  dmean = (pio6i*r1/n1)**onethird  !mass mean diam
  dn1 = dmean * (exp(gammln(gnu1) - gammln(gnu1+3.))) ** onethird

  do i=1,ibins
     ank1(i)=0.
     amk1(i)=0.
     ank2(i)=0.
     amk2(i)=0.
  enddo

  sum=0.
  sumn=0.
  do i=1,ithresh
     fac1=gammp(gnu1,diam(i)/dn1)
     fac2=gammp(gnu1,diam(i+1)/dn1)
     trunc1=fac2-fac1
     ank1(i)=n1*trunc1

     fac1=gammp(gnu1+3.,(diam(i)/dn1))
     fac2=gammp(gnu1+3.,(diam(i+1)/dn1))
     trunc=fac2-fac1
     amk1(i)=r1*trunc

     sum=sum+amk1(i)
     sumn=sumn+ank1(i)
  enddo

  ! * Scale to exactly r1,n1

  !      do i=1,ithresh
  !        ank1(i)=ank1(i)*n1/sumn
  !        amk1(i)=amk1(i)*r1/sum
  !      enddo

  gamp = exp(gammln(gnu2))
  dmean = (pio6i * r2 /n2) ** onethird   !mass mean diam
  dn2 = dmean * (exp(gammln(gnu2) - gammln(gnu2+3.))) ** onethird

  sum=0.
  sumn=0.

  do i=ithresh+1,ibins
     fac1=gammp(gnu2,diam(i)/dn2)
     fac2=gammp(gnu2,diam(i+1)/dn2)
     trunc=fac2-fac1
     ank2(i)=n2*trunc

     fac1=gammp(gnu2+3.,(diam(i)/dn2))
     fac2=gammp(gnu2+3.,(diam(i+1)/dn2))
     trunc=fac2-fac1
     amk2(i)=r2*trunc

     sum=sum+amk2(i)
     sumn=sumn+ank2(i)

  enddo

  ! * Scale to exactly r2,n2

  !      do i=ithresh+1,ibins
  !        ank2(i)=ank2(i)*n2/sumn
  !        amk2(i)=amk2(i)*r2/sum
  !      enddo

  do i=1,ibins
     ank(i)=ank1(i)+ank2(i)
     amk(i)=amk1(i)+amk2(i)
  enddo

  xntot  = 0.0
  xr3=0.

  do i=1,ibins
     xr3       = xr3 + amk(i)
     xntot     = xntot + ank(i)
  enddo
end subroutine initg

!******************************************************************************

subroutine sumn(ank,amk,imin,imax,ibins,sun,sum)

  implicit none

  ! Arguments:
  integer, intent(in) :: imin,imax,ibins
  real, dimension(ibins), intent(in) :: ank,amk
  real, intent(out) :: sun,sum

  !Local Variables:
  integer :: i

  sum=0.
  sun=0.

  do i=imin,imax
     sun=sun+ank(i)
     sum=sum+amk(i)
  enddo

  return
end subroutine sumn

!******************************************************************************

subroutine tabhab()

  use micphys, only: jhabtab !INTENT(OUT)
       
  implicit none

  ! Local Variables:
  integer, parameter :: nhab=0
  integer :: it,is

  !if (nhab .eq.  0) print*,'VARIABLE HABIT PREDICTION'
  !if (nhab .eq.  3) print*,'ASSUMED HABIT IS COLUMNS'
  !if (nhab .eq.  8) print*,'ASSUMED HABIT IS HEX PLATES'
  !if (nhab .eq.  9) print*,'ASSUMED HABIT IS DENDRITES'
  !if (nhab .eq. 10) print*,'ASSUMED HABIT IS NEEDLES'
  !if (nhab .eq. 11) print*,'ASSUMED HABIT IS ROSETTES'
  !c    if (nhab .eq.  x) print*,'ASSUMED HABIT IS SPHERES'

  ! nt is temp, ns = satur (liq)

  do it = 1,31
     do is = 1,100
        if (nhab .eq. 0) then
           if (it .ge. 0 .and. it .le. 2) then
              if (is .le. 95) then
                 jhabtab(it,is,1) = 3
                 jhabtab(it,is,2) = 4
              else
                 jhabtab(it,is,1) = 8
                 jhabtab(it,is,2) = 12
              endif
           else if(it .gt. 2 .and. it .le. 4) then
              if (is .lt. 90) then
                 jhabtab(it,is,1) = 3
                 jhabtab(it,is,2) = 4
              else
                 jhabtab(it,is,1) = 8
                 jhabtab(it,is,2) = 12
              endif
           else if(it .gt. 4 .and. it .le. 6) then
              if (is .lt. 85) then
                 jhabtab(it,is,1) = 3
                 jhabtab(it,is,2) = 4
              else
                 jhabtab(it,is,1) = 10
                 jhabtab(it,is,2) = 14
              endif
           else if(it .gt. 6 .and. it .le. 9) then
              if (is .lt. 90) then
                 jhabtab(it,is,1) = 3
                 jhabtab(it,is,2) = 4
              else
                 jhabtab(it,is,1) = 10
                 jhabtab(it,is,2) = 14
              endif
           else if(it .gt. 9 .and. it .le. 22) then
              if (is .lt. 90) then
                 jhabtab(it,is,1) = 8
                 jhabtab(it,is,2) = 12
              else
                 jhabtab(it,is,1) = 9
                 jhabtab(it,is,2) = 13
              endif
           elseif(it .gt. 22 .and. it .le. 30) then
              if (is .lt. 80) then
                 jhabtab(it,is,1) = 3
                 jhabtab(it,is,2) = 4
              else
                 jhabtab(it,is,1) = 10
                 jhabtab(it,is,2) = 14
              endif
           elseif(it .gt. 30) then
              if (is .lt. 90) then
                 jhabtab(it,is,1) = 3
                 jhabtab(it,is,2) = 4
              else
                 jhabtab(it,is,1) = 11
                 jhabtab(it,is,2) = 15
              endif
           endif
        else
           jhabtab(it,is,1) = nhab
           jhabtab(it,is,2) = nhab + 4
           if (nhab .eq. 3) jhabtab(it,is,2) = 4
        endif
     enddo
  enddo
  return
end subroutine tabhab
