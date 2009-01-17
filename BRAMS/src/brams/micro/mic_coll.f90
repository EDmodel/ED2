!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!  08/31/08 - MLO - Adjusted the routine to OLAM-style.                                    !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine getict(lcat)

   use micphys

   implicit none
   !----- Argument ------------------------------------------------------------------------!
   integer, intent(in) :: lcat
   !----- Local variables -----------------------------------------------------------------!
   integer :: k
   real    :: rict,rictmm
   !---------------------------------------------------------------------------------------!

   mainloop: do k = k1(lcat),k2(lcat)

      if (rx(k,lcat) < rxmin(lcat)) cycle mainloop

      rict         = dict(lcat) * (log(emb(k,lcat)) - emb0log(lcat)) + 1.
      rictmm       = max(rictmin,min(rictmax,rict))
      ict1(k,lcat) = int(rictmm)
      ict2(k,lcat) = ict1(k,lcat) + 1
      wct2(k,lcat) = rictmm - real(ict1(k,lcat))
      wct1(k,lcat) = 1.0 - wct2(k,lcat)

   end do mainloop
   return
end subroutine getict
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    I think this subroutine converts cloud droplets into rain drops due to collision and  !
! coalescence of cloud droplets (aka autoconversion), following Meyers et al. (1997),      !
! section 2.1.1.                                                                           !
!------------------------------------------------------------------------------------------!
subroutine auto_accret(dtlt)

   use micphys

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in) :: dtlt
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,id1cc,id1cr,id1crn,ir2cr,id2cr,ir2rr,id2rr
   real    :: dtlt3,dtlt6,dmb1cgs,dmb2cgs,r2cgs,en1cgs,ad1,ar2,d2minx,ad2
   real    :: bd1,br2,bd2,d2e,bd1cc,bd1cr,br2cr,bd2cr,br2rr,bd2rr,wd1cr,wr2dr
   real    :: wr2rr,wd2rr,tm1cc,tn1cc,tn2cc,tm1cr,tn1cr,tn2rr,en1cgs_2
   real    :: um1cc,un1cc,un2cc,um1cr,un1cr,un2rr,um2,un1,wr2cr
   !---------------------------------------------------------------------------------------!

   dtlt3 = 1.e3 * dtlt
   dtlt6 = 1.e6 * dtlt

   mainloop: do k = k1(1),k2(1)
      if (rx(k,1) < rxmin(1)) cycle mainloop

      !----- This subroutine works in cgs units, so convert inputs from mks ---------------!

      dmb1cgs = 100. * (emb(k,1) * cfmasi(1)) ** pwmasi(1)  ! cm
      dmb2cgs = 100. * (emb(k,2) * cfmasi(2)) ** pwmasi(2)  ! cm
      r2cgs   = 1.e-3 * rx(k,2) * rhoa(k)
      en1cgs  = 1.e-6 * cx(k,1) * rhoa(k)

      ad1     = max(d1min,min(d1max,dmb1cgs))
      ar2     = max(r2min,min(r2max,r2cgs))
      d2minx  = max(d2min,(r2cgs / (.1 * .5236)) ** pwmasi(2))
      ad2     = max(d2minx,min(d2max,dmb2cgs))

      bd1    = alog10(ad1/d1min)
      br2    = alog10(ar2/r2min)
      bd2    = alog10(ad2/d2minx)
      d2e    = alog10(d2max/d2minx)

      bd1cc = real(nd1cc-1) * (ad1 - d1min) / (d1max - d1min) + 1.
      bd1cr = bd1 / d1ecr + 1.
      br2cr = br2 / r2ecr + 1.
      bd2cr = bd2 / d2e * real(nd2cr-1) + 1.
      br2rr = br2 / r2err + 1.
      bd2rr = bd2 / d2e * real(nd2rr-1) + 1.

      id1cc  = nint(bd1cc)
      id1cr  =  int(bd1cr)
      id1crn = nint(bd1cr)
      ir2cr  =  int(br2cr)
      id2cr  = nint(bd2cr)
      ir2rr  =  int(br2rr)
      id2rr  =  int(bd2rr)

      wd1cr = bd1cr - real(id1cr)
      wr2cr = br2cr - real(ir2cr)
      wr2rr = br2rr - real(ir2rr)
      wd2rr = bd2rr - real(id2rr)

      tm1cc =                            r1tabcc(id1cc)

      tn1cc =                            c1tabcc(id1cc)

      tn2cc =                            c2tabcc(id1cc)

      tm1cr = (1.-wd1cr) * ((1.-wr2cr) * r1tabcr(id1cr  ,ir2cr  ,id2cr)                    &
            +                   wr2cr  * r1tabcr(id1cr  ,ir2cr+1,id2cr))                   &
            +     wd1cr  * ((1.-wr2cr) * r1tabcr(id1cr+1,ir2cr  ,id2cr)                    &
            +                   wr2cr  * r1tabcr(id1cr+1,ir2cr+1,id2cr))

      tn1cr =               (1.-wr2cr) * c1tabcr(id1crn,ir2cr  ,id2cr)                     &
            +                   wr2cr  * c1tabcr(id1crn,ir2cr+1,id2cr)

      tn2rr = (1.-wd2rr) * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr  )                          &
            +                   wr2rr  * c2tabrr(ir2rr+1,id2rr  ))                         &
            +     wd2rr  * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr+1)                          &
            +                   wr2rr  * c2tabrr(ir2rr+1,id2rr+1))

      en1cgs_2 = en1cgs ** 2

      um1cc = tm1cc * en1cgs_2      * dtlt3
      un1cc = tn1cc * en1cgs_2      * dtlt6
      un2cc = tn2cc * en1cgs_2      * dtlt6
      um1cr = 10. ** tm1cr * en1cgs * dtlt3
      un1cr = 10. ** tn1cr * en1cgs * dtlt6
      un2rr = 10. ** tn2rr          * dtlt6

      !------------------------------------------------------------------------------------!
      !     The above values are amounts in kg/m^3 or #/m^3 converted in the present time- !
      ! step, but must still be corrected for the effect of density on fall velocity.      !
      ! Thus, they must be multiplied by (rhoi ** .5) which fall velocity is proportional  !
      ! to.  Also, since rxfer and enxfer are in units of kg/kg and #/kg, respectively,    !
      ! the above transfer amounts must also be multiplied by rhoi. Together, these        !
      ! factors make (rhoi ** 1.5).                                                        !
      !------------------------------------------------------------------------------------!
      um2 = min(rx(k,1),(um1cc + um1cr) * rhoi(k))
      un1 = min(cx(k,1)*rhoa(k),(un1cc + un1cr))

      rxfer(k,1,2)  =  rxfer(k,1,2) + um2
      qrxfer(k,1,2) = qrxfer(k,1,2) + um2 * qx(k,1)
      enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
      enxfer(k,1,2) = enxfer(k,1,2) + un2cc

   end do mainloop
   return
end subroutine auto_accret
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine effxy(m1)

   use micphys
   use rconstants, only : cicet3,t00
   use micro_coms, only : ticegrowth

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)       :: m1
   !----- Local variables -----------------------------------------------------------------!
   integer                   :: k
   real                      :: dmr
   logical            , save :: firstcall7 = .true.
   !---------------------------------------------------------------------------------------!

   !----- 1 = rp,rs,ra,rg,rh --------------------------------------------------------------!

   if (firstcall7 .and. availcat(2) .and. availcat(3)) then
      firstcall7 = .false.
      do k = 2,m1-1
         eff(k,1) = 1.0
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !----- 2 = cs,ca -----------------------------------------------------------------------!
   if (availcat(2) .or. availcat(3)) then
      do k = k1(1),k2(1)
         !---------------------------------------------------------------------------------!
         ! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:                          !
         !  close to curve for 404 microns.  Replace with auto_accret eventually.          !
         !---------------------------------------------------------------------------------!
         if (emb(k,1) > 9.e-13) then
            eff(k,2) = min(1.,30. * (emb(k,1) - 9.e-13) ** .15)
         else
            eff(k,2) = 0.
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !----- 3 = rr --------------------------------------------------------------------------!
   if (availcat(2)) then
      rainloop: do k = k1(2),k2(2)

         if (rx(k,2) < rxmin(2)) cycle rainloop

         !----- Rain breakup (temporary; eventually combine with autoconv/accret). --------!
         if (emb(k,2) < .113e-6) then
            eff(k,3) = 1.0
         elseif (emb(k,2) > .158e-5) then
            eff(k,3) = -5.0
         else
            eff(k,3) = 2. - exp(.1326e7 * (emb(k,2) - .113e-6))
         end if
      end do rainloop
   end if
   !---------------------------------------------------------------------------------------!


   !---- 4 = pp,ps,pa ---------------------------------------------------------------------!
   if (availcat(5)) then
      do k = k1(3),k2(3)
         if (abs(tx(k,3)-ticegrowth) <= 2.) then
            eff(k,4) = 1.4
         else
            eff(k,4) = min(0.2,10. ** (0.035 * (tx(k,3)-t00) - 0.7))
         endif

      enddo
      !------------------------------------------------------------------------------------!


      !----- 5 = ss,sa --------------------------------------------------------------------!
      do k = k1(4),k2(4)
         if (abs(tx(k,4)-ticegrowth) <= 2.) then
            eff(k,5) = 1.4
         else
            eff(k,5) = min(0.2,10. ** (0.035 * (tx(k,4)-t00) - 0.7))
         end if
      end do
      !------------------------------------------------------------------------------------!


      !----- 6 = aa -----------------------------------------------------------------------!
      aggrloop: do k = k1(5),k2(5)

         if (rx(k,5) < rxmin(5)) cycle aggrloop

         if (abs(tx(k,5)-ticegrowth) <= 2.) then
            eff(k,6) = 1.4
         elseif (tx(k,5) >= -1.) then
            eff(k,6) = 1.
         else
            eff(k,6) = min(0.2,10. ** (0.035 * (tx(k,5)-t00) - 0.7))
         end if
      end do aggrloop
      !------------------------------------------------------------------------------------!

   end if

   !----- 7 = pg,sg,ag,gg,gh --------------------------------------------------------------!
   if (availcat(6)) then
      graupelloop: do k = k1(6),k2(6)
         if (rx(k,6) < rxmin(6)) cycle graupelloop
         if (qr(k,6) > rx(k,6)*cicet3) then
            eff(k,7) = 1.0
         else
            eff(k,7) = min(0.2,10. ** (0.035 * (tx(k,6)-t00) - 0.7))
         end if
      end do graupelloop
   end if
   !---------------------------------------------------------------------------------------!

   !----- 8 = ph,sh,ah,gh -----------------------------------------------------------------!
   if (availcat(7)) then
      hailloop:do k = k1(7),k2(7)
         if (rx(k,7) < rxmin(7)) cycle hailloop

         if (qr(k,7) > rx(k,7)*cicet3) then
            eff(k,8) = 1.0
         else
            eff(k,8) = min(0.2,10. ** (0.035 * (tx(k,7)-t00) - 0.7))
         end if
      end do hailloop
   end if
   !---------------------------------------------------------------------------------------!


   !----- 9 = cg,ch -----------------------------------------------------------------------!
   if (availcat(2) .or. availcat(3)) then
      do k = k1(1),k2(1)

         !---------------------------------------------------------------------------------!
         ! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:                          !
         !  close to curves for 142 and 305 microns.  Replace with auto_accret eventually. !
         !---------------------------------------------------------------------------------!
         if (emb(k,1) > 3.4e-14) then
            eff(k,9) = min(1.,1426. * (emb(k,1) - 3.4e-14) ** .28)
         else
            eff(k,9) = 0.
         end if
      end do
   end if
   !---------------------------------------------------------------------------------------!


   !---- 10 = hh (trial) ------------------------------------------------------------------!
   if (availcat(7)) then
      do k = k1(7),k2(7)
         eff(k,10) = max(0.,.1 + .005 * (tx(k,7)-t00))
      enddo
   endif
   !---------------------------------------------------------------------------------------!

   return
end subroutine effxy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine cols(mx,mc1)

   use micphys

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: mx,mc1
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ipc,k
   real                :: colnum,tabval
   !---------------------------------------------------------------------------------------!

   mainloop: do k = k1(mx),k2(mx)
      if(rx(k,mx) < rxmin(mx)) cycle mainloop


      ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

      tabval = wct1(k,mx) ** 2              * coltabc(ict1(k,mx),ict1(k,mx),ipc)           &
             + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)           &
             + wct2(k,mx) ** 2              * coltabc(ict2(k,mx),ict2(k,mx),ipc)

      colnum = colfacc(k) * eff(k,mc1) * cx(k,mx) ** 2 * 10. ** (-tabval)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(0.5 * cx(k,mx),colnum)
   end do mainloop

   return
end subroutine cols
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine col3344(mx,mz,mc1)

   use micphys

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: mx,mz,mc1
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,ip,ipc
   real :: c1,tabvalx,colamt,tabvaln,colnum
   !---------------------------------------------------------------------------------------!

   mainloop: do k = k1(mx),k2(mx)

      if(rx(k,mx) < rxmin(mx)) cycle mainloop

      ip  = ipairr(jhcat(k,mx),jhcat(k,mx))
      ipc = ipairc(jhcat(k,mx),jhcat(k,mx))
      c1  = eff(k,mc1) * cx(k,mx) ** 2

      tabvalx  = wct1(k,mx) ** 2               * coltabr(ict1(k,mx),ict1(k,mx),ip)         &
               + 2. * wct1(k,mx) * wct2(k,mx)  * coltabr(ict1(k,mx),ict2(k,mx),ip)         &
               + wct2(k,mx) ** 2               * coltabr(ict2(k,mx),ict2(k,mx),ip)

      colamt          = min(rx(k,mx),colfacr2(k) * c1 * 10. ** (-tabvalx))
      rxfer(k,mx,mz)  = rxfer(k,mx,mz)  + colamt
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colamt * qx(k,mx)

      if (.not. progncat(mz)) cycle mainloop

      tabvaln = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)         &
              + 2. * wct1(k,mx) * wct2(k,mx)  * coltabc(ict1(k,mx),ict2(k,mx),ipc)         &
              + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

      colnum          = min(0.5 * cx(k,mx),colfacc2(k) * c1 * 10. ** (-tabvaln))
      enxfer(k,mx,mz) = enxfer(k,mx,mz) + colnum
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnum

   end do mainloop
   return
end subroutine col3344
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine col3443(mx,my,mz)

   use micphys

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer :: mx,my,mz
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,jhcatx,jhcaty,ipxy,ipyx,ipc
   real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum
   !---------------------------------------------------------------------------------------!

   mainloop: do k = k1(mx),k2(mx)  ! k1,k2 now same for lcat mx & my (aka 3 & 4)

      if(rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle mainloop


      jhcatx = jhcat(k,mx)
      jhcaty = jhcat(k,my)
      ipxy   = ipairr(jhcatx,jhcaty)
      ipyx   = ipairr(jhcaty,jhcatx)
      ipc    = ipairc(jhcatx,jhcaty)
      c1     = eff(k,4) * cx(k,mx) * cx(k,my)

      tabvalx = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)             &
              + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)             &
              + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)             &
              + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)
      rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

      tabvaly = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)             &
              + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)             &
              + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)             &
              + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)
      rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

      rxfer(k,mx,mz)  = rxfer(k,mx,mz)  + rcx
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

      rxfer(k,my,mz)  = rxfer(k,my,mz)  + rcy
      qrxfer(k,my,mz) = qrxfer(k,my,mz) + rcy * qx(k,my)

      tabvaln = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)              &
              + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)              &
              + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)              &
              + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)
      colnum = c1 * colfacc(k) * 10. ** (-tabvaln)

      if (cx(k,mx) > cx(k,my)) then
         enxfer(k,my,mz) = min(cx(k,my),colnum)
         enxfer(k,mx,mx) = min(cx(k,mx),colnum)
      else
         enxfer(k,mx,mz) = min(cx(k,mx),colnum)
         enxfer(k,my,my) = min(cx(k,my),colnum)
      end if

      !---- also loss for aerosol... ------------------------------------------------------!

   end do mainloop

   return
end subroutine col3443
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine col1(mx,my,mz,mc4,j1,j2)

   use micphys

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: mx,my,mz,mc4,j1,j2
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,ipxy,ipc
   real                :: c1,tabvalx,rcx,tabvaln,colnum
   !---------------------------------------------------------------------------------------!


   mainloop: do k = j1,j2
      
      if(rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle mainloop

      
      ipxy = ipairr(jhcat(k,mx),jhcat(k,my))
      ipc  = ipairc(jhcat(k,mx),jhcat(k,my))
      c1   = eff(k,mc4) * cx(k,mx) * cx(k,my)

      tabvalx = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)             &
              + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)             &
              + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)             &
              + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)


      rcx             = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))
      rxfer(k,mx,mz)  = rxfer(k,mx,mz) + rcx
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

      if (.not. progncat(mx)) cycle mainloop

      tabvaln = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)              &
              + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)              &
              + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)              &
              + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

      colnum          = c1 * colfacc(k) * 10. ** (-tabvaln)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))
   end do mainloop

   return
end subroutine col1
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine col2(mx,my,mz,mc2,j1,j2,dtlt)

   use rconstants
   use micphys
   use therm_lib, only : qtc
   use micro_coms, only : alpha_coll2,beta_coll2

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: mx,my,mz,mc2,j1,j2
   real   , intent(in) :: dtlt
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,jhcatx,jhcaty,ipxy,ipyx,ipc,it
   real                :: c1,c2,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum0,colnum,rcoal
   real                :: qrcx,qrcy,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice,area
   real                :: cn13,cn24,sip,rsip,qrsip,rfinlz,xtoz
   !---------------------------------------------------------------------------------------!


   mainloop: do k = j1,j2
      
      if(rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle mainloop
      
      jhcatx = jhcat(k,mx)
      jhcaty = jhcat(k,my)
      ipxy   = ipairr(jhcatx,jhcaty)
      ipyx   = ipairr(jhcaty,jhcatx)
      ipc    = ipairc(jhcatx,jhcaty)
      c2     = cx(k,mx) * cx(k,my)
      c1     = eff(k,mc2) * c2

      tabvalx = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)             &
              + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)             &
              + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)             &
              + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)
      rcx     = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))


      tabvaly = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)             &
              + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)             &
              + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)             &
              + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)
      rcy     = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

      if (progncat(mx) .or. progncat(my)) then 
         tabvaln = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)           &
                 + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)           &
                 + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)           &
                 + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

         colnum0 = c2 * colfacc(k) * 10. ** (-tabvaln)
         colnum  = colnum0 * eff(k,mc2)
      else
         !---------------------------------------------------------------------------------!
         ! MLO - Not sure about that, but it may happen that we need colnum and the above  !
         !       if is the only place in which colnum is defined in this sub-routine, so   !
         !       I need to add something, here, so I put 0.                                !
         !       Code dilemma, should the entire IF be removed???                          !
         !---------------------------------------------------------------------------------!
         colnum0=0.
         colnum=0.
      end if

      rcoal  = rcx + rcy
      qrcx   = rcx * qx(k,mx)
      qrcy   = rcy * qx(k,my)
      qrcoal = qrcx + qrcy
      qcoal  = qrcoal / max(1.e-13,rcoal)

      call qtc(qcoal,tcoal,fracliq)

      coalliq = rcoal * fracliq
      coalice = rcoal - coalliq

      !------------------------------------------------------------------------------------!
      !    Secondary ice production: cn24 is the number fraction of collected cloud drop-  !
      ! lets larger than 24 microns and is obtained from an incomplete gamma function      !
      ! table.  cn13 is the fraction of collected cloud droplets smaller than 13 microns.  !
      ! Area is cross section area of collecting ice per m^3 of atmospheric volume.        !
      !------------------------------------------------------------------------------------!
      if (tcoal > -8. .and. tcoal < -3.) then

         area = cx(k,my) * rhoa(k) * sipfac(jhcaty) * emb(k,my) ** (2.*pwmasi(jhcaty))
         it   = nint(emb(k,mx) / emb1(1) * ngam)
         cn13 = colnum * gamsip13(it) / (area * dtlt)
         cn24 = min(cx(k,mx)*rhoa(k),colnum0) * gamsip24(it)
         sip  = 9.1e-10 * cn24 * cn13 ** .93
         if (tcoal < -5.) then
            sip = onethird * (tcoal + 8.) * sip
         else
            sip = -0.5 * (tcoal + 3.) * sip
         end if

         rsip  = sip * emb0(3) * rhoi(k)
         qrsip = qcoal * rsip

         rcoal  = rcoal - rsip
         qrcoal = qrcoal - qrsip

         enxfer(k,mx,3) = enxfer(k,mx,3) + sip
         rxfer (k,mx,3) = rxfer (k,mx,3) + rsip
         qrxfer(k,mx,3) = qrxfer(k,mx,3) + qrsip

      end if

      !------------------------------------------------------------------------------------!
      !     ALWAYS NEED (ALPHA + BETA) >= 1 but in the (rare) case that fracliq may be a   !
      ! little larger than fracx due to collected liquid being above the triple point,     !
      ! need (ALPHA + BETA) to be at least 1.1 or 1.2, or need ALPHA itself to be at least !
      ! 1.0.                                                                               !
      !------------------------------------------------------------------------------------!
      rfinlz = min(rcoal, alpha_coll2(jhcaty) * coalliq + beta_coll2(jhcaty) * rcx)

      xtoz = min(rcx,rfinlz)

      rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
      rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
      if (my /= mz) rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz

      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
      qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
      if (my /= mz) qrxfer(k,my,mz) = qrxfer(k,my,mz) + qx(k,my) * (rfinlz - xtoz)

      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))
      if (my /= mz) then
         enxfer(k,my,mz) = enxfer(k,my,mz) + (rfinlz-xtoz) * min(colnum,cx(k,my))          &
                         / (1.e-20 + rcy)
      end if

      !----- BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y ----------------!

      !----- also include loss of aerosol -------------------------------------------------!
   end do mainloop

   return
end subroutine col2
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine col3(mx,my,mz,j1,j2)

   use micphys
   use therm_lib, only : qtc
   use micro_coms, only : alpha_coll3,beta_coll3

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: mx,my,mz,j1,j2
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,ipxy,ipyx,ipc,jhcaty
   real                :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum
   real                :: colnumx,colnumy,coalnum,rcoal,qrcx,qrcy,qrcoal,qcoal,fracliq
   real                :: coalliq,coalice,xtoz,rfinlz,tcoal,cfinlz
   !---------------------------------------------------------------------------------------!
   mainloop: do k = j1,j2

      if(rx(k,mx) < rxmin(mx) .or. rx(k,my) < rxmin(my)) cycle mainloop

      jhcaty = jhcat(k,my)
      ipxy   = ipairr(jhcat(k,mx),jhcaty)
      ipyx   = ipairr(jhcaty,jhcat(k,mx))
      ipc    = ipairc(jhcat(k,mx),jhcaty)
      c1     = eff(k,1) * cx(k,mx) * cx(k,my)

      tabvalx = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)             &
              + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)             &
              + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)             &
              + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)
      rcx     = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

      tabvaly = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)             &
              + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)             &
              + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)             &
              + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)
      rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

      if (progncat(mx)) then
         tabvaln = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)           &
                 + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)           &
                 + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)           &
                 + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

         colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
         colnumx = min(cx(k,mx),colnum)
         colnumy = min(cx(k,my),colnum)
         coalnum = min(colnumx,colnumy)
      end if

      rcoal  = rcx + rcy
      qrcx   = rcx * qx(k,mx)
      qrcy   = rcy * qx(k,my)
      qrcoal = qrcx + qrcy
      qcoal  = qrcoal / (1.e-20 + rcoal)

      call qtc(qcoal,tcoal,fracliq)

      coalliq = rcoal * fracliq
      coalice = rcoal - coalliq

      if (fracliq >= .99) then

         rxfer(k,my,mx)  = rxfer(k,my,mx) + rcy
         qrxfer(k,my,mx) = qrxfer(k,my,mx) + qrcy
         if (progncat(mx)) enxfer(k,my,my) = enxfer(k,my,my) + colnumy
      else

         rfinlz = min(rcoal,alpha_coll3(jhcaty) * coalliq + beta_coll3(jhcaty) * rcx)

         xtoz = min(rcx,rfinlz)

         rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
         rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
         if (my /= mz) rxfer(k,my,mz) = rxfer(k,my,mz) + rfinlz - xtoz

         !----- NEED TO USE QCOAL TO TRANSFER Q? ------------------------------------------!
         qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
         qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
         if (my .ne. mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
            + qx(k,my) * (rfinlz - xtoz)

         if (progncat(mx)) then
            if (my == mz) then
               enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
            elseif (colnumy >= colnumx) then
               cfinlz          = coalnum * rfinlz / (rcoal + 1.e-20)
               enxfer(k,mx,mz) = enxfer(k,mx,mz) + cfinlz
               enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx - cfinlz
               enxfer(k,my,my) = enxfer(k,my,my) + colnumy
            else
               cfinlz          = coalnum * rfinlz / (rcoal + 1.e-20)
               enxfer(k,my,mz) = enxfer(k,my,mz) + cfinlz
               enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
               enxfer(k,my,my) = enxfer(k,my,my) + colnumy - cfinlz
            end if
         end if

      end if
   end do mainloop

   !----- also include loss of aerosol ----------------------------------------------------!
   return
end subroutine col3
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine colxfers()

   use micphys

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer                               :: k,lcat,jcat,kd1,kd2
   !---------------------------------------------------------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   !    All rxfer values are nonnegative.                                                  !
   !---------------------------------------------------------------------------------------!
   do lcat = 1,7
      if (availcat(lcat)) then
         kd1 = k1(lcat)
         kd2 = k2(lcat)

         do k = kd1,kd2
            rloss(k)  = 0.
            enloss(k) = 0.
         end do

         do jcat = 1,7
            !----- Change this to include enxfer of the same categories -------------------!
            if (availcat(jcat)) then
               if (lcat /= jcat) then
                  do k = kd1,kd2
                     rloss(k) = rloss(k) + rxfer(k,lcat,jcat)
                  end do
               endif
               do k = kd1,kd2
                  enloss(k) = enloss(k) + enxfer(k,lcat,jcat)
               end do
            end if
         end do

         do k = kd1,kd2
            rloss(k) = min(1.,rx(k,lcat) / max(1.e-20,rloss(k)))
            enloss(k) = min(1.,cx(k,lcat) / max(1.e-10,enloss(k)))
         end do

         do jcat = 1,7
            if (availcat(jcat)) then
               if (lcat /= jcat) then
                  do k = kd1,kd2
                     rxfer(k,lcat,jcat) = rxfer(k,lcat,jcat)*rloss(k)
                     qrxfer(k,lcat,jcat)=qrxfer(k,lcat,jcat)*rloss(k)
                  end do
               end if
               do k = kd1,kd2
                  enxfer(k,lcat,jcat) = enxfer(k,lcat,jcat)*enloss(k)
               end do
            end if
         end do
      end if
   end do

   do lcat = 1,7

      if (availcat(lcat)) then

         kd1 = k1(lcat)
         kd2 = k2(lcat)

         do jcat = 1,7
            if (availcat(jcat) .and. lcat /= jcat) then
               do k = kd1,kd2
                  rx(k,lcat) = rx(k,lcat) - rxfer(k,lcat,jcat)
                  rx(k,jcat) = rx(k,jcat) + rxfer(k,lcat,jcat)
                  qr(k,lcat) = qr(k,lcat) - qrxfer(k,lcat,jcat)
                  qr(k,jcat) = qr(k,jcat) + qrxfer(k,lcat,jcat)
                  cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,jcat)
                  cx(k,jcat) = cx(k,jcat) + enxfer(k,lcat,jcat)
               end do
            end if
         end do

         if (progncat(lcat)) then
            do k = kd1,kd2
               cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,lcat)
            end do
         end if

      end if
   end do
   return
end subroutine colxfers


