!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine getict(k1,k2,lcat,i,j,mynum)

   use micphys

   implicit none

   integer :: k1,k2,lcat,i,j,k,mynum
   real :: rict,rictmm

   do k = k1,k2

      if (rx(k,lcat) >= rxmin(lcat)) then

      rict = dict(lcat) * (log(emb(k,lcat)) - emb0log(lcat)) + 1.
      rictmm = max(rictmin,min(rictmax,rict))
      ict1(k,lcat) = int(rictmm)
      ict2(k,lcat) = ict1(k,lcat) + 1
      wct2(k,lcat) = rictmm - float(ict1(k,lcat))
      wct1(k,lcat) = 1.0 - wct2(k,lcat)

      endif

   enddo
   return
end subroutine getict

!******************************************************************************

subroutine auto_accret(m1,k1,k2,dn0,dtlt,i,j)

  use micphys

  implicit none

  integer :: m1,i,j,k,k1,k2,id1cc,id1cr,id1crn,ir2cr,id2cr,ir2rr,id2rr
  real :: dtlt3,dtlt6,dmb1cgs,dmb2cgs,r2cgs,en1cgs,ad1,ar2,d2minx,ad2  &
         ,bd1,br2,bd2,d2e,bd1cc,bd1cr,br2cr,bd2cr,br2rr,bd2rr,wd1cr,wr2dr  &
         ,wr2rr,wd2rr,tm1cc,tn1cc,tn2cc,tm1cr,tn1cr,tn2rr,en1cgs_2  &
         ,um1cc,un1cc,un2cc,um1cr,un1cr,un2rr,um2,un1,dtlt,cfmasi1,cfmasi2  &
         ,pwmasi1,pwmasi2,wr2cr

  real, dimension(m1) :: dn0

  dtlt3 = 1.e3 * dtlt
  dtlt6 = 1.e6 * dtlt
  cfmasi1 = 1. / cfmas(1)
  cfmasi2 = 1. / cfmas(2)
  pwmasi1 = 1. / pwmas(1)
  pwmasi2 = 1. / pwmas(2)

  do k = k1,k2
  if(rx(k,1) >= rxmin(1)) then

  ! This subroutine works in cgs units, so convert inputs from mks

     dmb1cgs = 100. * (emb(k,1) * cfmasi1) ** pwmasi1
     dmb2cgs = 100. * (emb(k,2) * cfmasi2) ** pwmasi2
     r2cgs = 1.e-3 * rx(k,2) * dn0(k)
     en1cgs = 1.e-6 * cx(k,1) * dn0(k)

     ad1 = max(d1min,min(d1max,dmb1cgs))
     ar2 = max(r2min,min(r2max,r2cgs))
      d2minx = max(d2min,(r2cgs / (.1 * .5236)) ** pwmasi2)
     ad2 = max(d2minx,min(d2max,dmb2cgs))

     bd1 = alog10(ad1/d1min)
     br2 = alog10(ar2/r2min)
     bd2 = alog10(ad2/d2minx)
        d2e =  alog10(d2max / d2minx)

     bd1cc = float(nd1cc-1) * (ad1 - d1min) / (d1max - d1min) + 1.
     bd1cr = bd1 / d1ecr + 1.
     br2cr = br2 / r2ecr + 1.
     bd2cr = bd2 / d2e * float(nd2cr-1) + 1.
     br2rr = br2 / r2err + 1.
     bd2rr = bd2 / d2e * float(nd2rr-1) + 1.

  !         id1cc  =  int(bd1cc)
     id1cc  =  nint(bd1cc)
     id1cr  =  int(bd1cr)
     id1crn = nint(bd1cr)
     ir2cr  =  int(br2cr)
     id2cr  = nint(bd2cr)
     ir2rr  =  int(br2rr)
     id2rr  =  int(bd2rr)

     wd1cr = bd1cr - float(id1cr)
     wr2cr = br2cr - float(ir2cr)
     wr2rr = br2rr - float(ir2rr)
     wd2rr = bd2rr - float(id2rr)

     tm1cc =                            r1tabcc(id1cc)

     tn1cc =                            c1tabcc(id1cc)

     tn2cc =                            c2tabcc(id1cc)

     tm1cr = (1.-wd1cr) * ((1.-wr2cr) * r1tabcr(id1cr  ,ir2cr  ,id2cr)  &
           +                   wr2cr  * r1tabcr(id1cr  ,ir2cr+1,id2cr))  &
           +     wd1cr  * ((1.-wr2cr) * r1tabcr(id1cr+1,ir2cr  ,id2cr)  &
           +                   wr2cr  * r1tabcr(id1cr+1,ir2cr+1,id2cr))

     tn1cr =               (1.-wr2cr) * c1tabcr(id1crn,ir2cr  ,id2cr)  &
           +                   wr2cr  * c1tabcr(id1crn,ir2cr+1,id2cr)

     tn2rr = (1.-wd2rr) * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr  )  &
           +                   wr2rr  * c2tabrr(ir2rr+1,id2rr  ))  &
           +     wd2rr  * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr+1)  &
           +                   wr2rr  * c2tabrr(ir2rr+1,id2rr+1))

     en1cgs_2 = en1cgs ** 2

     um1cc = tm1cc * en1cgs_2 * dtlt3
     un1cc = tn1cc * en1cgs_2 * dtlt6
     un2cc = tn2cc * en1cgs_2 * dtlt6
     um1cr = 10. ** tm1cr * en1cgs * dtlt3
     un1cr = 10. ** tn1cr * en1cgs * dtlt6
     un2rr = 10. ** tn2rr * dtlt6

  ! The above values are amounts in kg/m^3 or #/m^3 converted in the
  ! present timestep, but must still be corrected for the effect of
  ! density on fall velocity.  Thus, they must be multiplied by
  ! (dn0i ** .5) which fall velocity is proportional to.  Also, since
  ! rxfer and enxfer are in units of kg/kg and #/kg, respectively, the
  ! above transfer amounts must also be multiplied by dn0i.  Together,
  ! these factors make (dn0i ** 1.5).

     um2 = min(rx(k,1),(um1cc + um1cr) * dn0i(k))
     un1 = min(cx(k,1)*dn0(k),(un1cc + un1cr))

     rxfer(k,1,2)  =  rxfer(k,1,2) + um2
     qrxfer(k,1,2) = qrxfer(k,1,2) + um2 * qx(k,1)
     enxfer(k,1,1) = enxfer(k,1,1) + un1 - un2cc
     enxfer(k,1,2) = enxfer(k,1,2) + un2cc

  ! no collis breakup yet - do not use next line but use col(2,2) in 3d micro

  !cc         enxfer(k,2,2) = enxfer(k,2,2) + un2rr

  ! aerosol loss here?

  endif
  enddo
  return
end subroutine auto_accret

!******************************************************************************

subroutine effxy(m1,k1,k2,i,j)

   use micphys

   implicit none

   integer :: m1,i,j,k,ncall7
   integer, dimension(10) :: k1,k2
   real :: dmr
   data ncall7/0/
   save

   !     1 = rp,rs,ra,rg,rh

   if (ncall7 .eq. 0 .and. jnmb(2) .ge. 1 .and. jnmb(3) .ge. 1) then
      ncall7 = 7
      do k = 2,m1-1
         eff(k,1) = 1.0
      enddo
   endif

   !     2 = cs,ca

   if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
      do k = k1(1),k2(1)

   ! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
   !  close to curve for 404 microns.  Replace with auto_accret eventually.

         if (emb(k,1) .gt. 9.e-13) then
            eff(k,2) = min(1.,30. * (emb(k,1) - 9.e-13) ** .15)
         else
            eff(k,2) = 0.
         endif
      enddo
   endif

   !     3 = rr

   if (jnmb(2) .ge. 1) then
      do k = k1(2),k2(2)

         if (rx(k,2) >= rxmin(2)) then

   ! rain breakup (old)

   !            dmr = dn(k,2) * gnu2
   !            if (dmr .lt. .0006) then
   !               eff(k,3) = 1.0
   !            elseif (dmr .gt. .001446) then
   !               eff(k,3) = -5.0
   !            else
   !               eff(k,3) = exp(2300. * (dmr - .0006))
   !            endif

   ! rain breakup (new - temporary; eventually combine with autoconv/accret

         if (emb(k,2) .lt. .113e-6) then
            eff(k,3) = 1.0
         elseif (emb(k,2) .gt. .158e-5) then
            eff(k,3) = -5.0
         else
            eff(k,3) = 2. - exp(.1326e7 * (emb(k,2) - .113e-6))
         endif

         endif

      enddo

   endif

   !     4 = pp,ps,pa

   if (jnmb(5) .ge. 1) then
      do k = k1(3),k2(3)
         if (abs(tx(k,3)+14.) .le. 2.) then
            eff(k,4) = 1.4
         else
            eff(k,4) = min(0.2,10. ** (0.035 * tx(k,3) - 0.7))
         endif

      enddo

   !     5 = ss,sa

      do k = k1(4),k2(4)
         if (abs(tx(k,4)+14.) .le. 2.) then
            eff(k,5) = 1.4
         else
            eff(k,5) = min(0.2,10. ** (0.035 * tx(k,4) - 0.7))
         endif
      enddo

   !     6 = aa

      do k = k1(5),k2(5)

         if (rx(k,5) >= rxmin(5)) then

         if (abs(tx(k,5)+14.) .le. 2.) then
            eff(k,6) = 1.4
         elseif (tx(k,5) .ge. -1.) then
            eff(k,6) = 1.
         else
            eff(k,6) = min(0.2,10. ** (0.035 * tx(k,5) - 0.7))
         endif

         endif

      enddo
   endif

   !     7 = pg,sg,ag,gg,gh

   if (jnmb(6) .ge. 1) then
      do k = k1(6),k2(6)
         if (qr(k,6) .gt. 0.) then
            eff(k,7) = 1.0
         else
            eff(k,7) = min(0.2,10. ** (0.035 * tx(k,6) - 0.7))
         endif
      enddo
   endif

   !     8 = ph,sh,ah,gh

   if (jnmb(7) .ge. 1) then
      do k = k1(7),k2(7)

         if (rx(k,7) >= rxmin(7)) then

         if (qr(k,7) .gt. 0.) then
            eff(k,8) = 1.0
         else
            eff(k,8) = min(0.2,10. ** (0.035 * tx(k,7) - 0.7))
         endif

         endif

      enddo
   endif

   !     9 = cg,ch

   if (jnmb(2) .ge. 1 .or. jnmb(3) .ge. 1) then
      do k = k1(1),k2(1)


   ! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
   !  close to curves for 142 and 305 microns.  Replace with auto_accret eventually.

         if (emb(k,1) .gt. 3.4e-14) then
            eff(k,9) = min(1.,1426. * (emb(k,1) - 3.4e-14) ** .28)
         else
            eff(k,9) = 0.
         endif
      enddo
   endif

   !     10 = hh (trial)

   if (jnmb(7) .ge. 1) then
      do k = k1(7),k2(7)
         eff(k,10) = max(0.,.1 + .005 * tx(k,7))
      enddo
   endif

   return
end subroutine effxy

!*********************************************************************************

subroutine cols(m1,mx,mc1,k1,k2,i,j)

   use micphys

   implicit none

   integer :: ipc,m1,mx,mc1,k1,k2,i,j,k
   real :: colnum,tabval

   do k = k1,k2
   if(rx(k,mx) >= rxmin(mx)) then
      ipc = ipairc(jhcat(k,mx),jhcat(k,mx))

      tabval  &
      = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
      + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
      + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

      colnum = colfacc(k) * eff(k,mc1) * cx(k,mx) ** 2 * 10. ** (-tabval)
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(0.5 * cx(k,mx),colnum)
   endif
   enddo
   return
end subroutine cols

!******************************************************************************

subroutine col3344(m1,mx,mz,mc1,k1,k2,i,j)

   use micphys

   implicit none

   integer :: m1,mx,mz,mc1,k1,k2,i,j,k,ip,ipc
   real :: c1,tabvalx,colamt,tabvaln,colnum

   do k = k1,k2
   if(rx(k,mx) >= rxmin(mx)) then
      ip = ipairr(jhcat(k,mx),jhcat(k,mx))
      ipc = ipairc(jhcat(k,mx),jhcat(k,mx))
      c1 = eff(k,mc1) * cx(k,mx) ** 2

      tabvalx  &
       = wct1(k,mx) ** 2               * coltabr(ict1(k,mx),ict1(k,mx),ip)  &
       + 2. * wct1(k,mx) * wct2(k,mx) * coltabr(ict1(k,mx),ict2(k,mx),ip)  &
       + wct2(k,mx) ** 2               * coltabr(ict2(k,mx),ict2(k,mx),ip)

      colamt = min(rx(k,mx),colfacr2(k) * c1 * 10. ** (-tabvalx))
      rxfer(k,mx,mz) = rxfer(k,mx,mz) + colamt
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + colamt * qx(k,mx)

      if (jnmb(mz) >= 5) then

      tabvaln  &
      = wct1(k,mx) ** 2               * coltabc(ict1(k,mx),ict1(k,mx),ipc)  &
      + 2. * wct1(k,mx) * wct2(k,mx) * coltabc(ict1(k,mx),ict2(k,mx),ipc)  &
      + wct2(k,mx) ** 2               * coltabc(ict2(k,mx),ict2(k,mx),ipc)

      colnum = min(0.5 * cx(k,mx),colfacc2(k) * c1 * 10. ** (-tabvaln))
      enxfer(k,mx,mz) = enxfer(k,mx,mz) + colnum
      enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnum

      endif
   endif
   enddo
   return
end

!******************************************************************************

subroutine col3443(m1,mx,my,mz,k1,k2,i,j)

   use micphys

   implicit none

   integer :: m1,mx,my,mz,k1,k2,i,j,k,jhcatx,jhcaty,ipxy,ipyx,ipc
   real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum

   do k = k1,k2
   if(rx(k,mx) >= rxmin(mx) .and. rx(k,my) >= rxmin(my)) then
      jhcatx = jhcat(k,mx)
      jhcaty = jhcat(k,my)
      ipxy = ipairr(jhcatx,jhcaty)
      ipyx = ipairr(jhcaty,jhcatx)
      ipc  = ipairc(jhcatx,jhcaty)
      c1 = eff(k,4) * cx(k,mx) * cx(k,my)

      tabvalx  &
        = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
        + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
        + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
        + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)
      rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

      tabvaly  &
        = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
        + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
        + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
        + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)
      rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

      rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

      rxfer(k,my,mz) = rxfer(k,my,mz) + rcy
      qrxfer(k,my,mz) = qrxfer(k,my,mz) + rcy * qx(k,my)

      tabvaln  &
          = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
          + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
          + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
          + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)
      colnum = c1 * colfacc(k) * 10. ** (-tabvaln)

         if (cx(k,mx) .gt. cx(k,my)) then
            enxfer(k,my,mz) = min(cx(k,my),colnum)
            enxfer(k,mx,mx) = min(cx(k,mx),colnum)
         else
            enxfer(k,mx,mz) = min(cx(k,mx),colnum)
            enxfer(k,my,my) = min(cx(k,my),colnum)
         endif

   ! also loss for aerosol

   endif
   enddo
   return
end subroutine col3443

!******************************************************************************

subroutine col1(m1,mx,my,mz,mc4,k1,k2,i,j)

   use micphys

   implicit none

   integer :: m1,mx,my,mz,mc4,k1,k2,i,j,k,ipxy,ipc
   real :: c1,tabvalx,rcx,tabvaln,colnum

   do k = k1,k2
   if(rx(k,mx) >= rxmin(mx) .and. rx(k,my) >= rxmin(my)) then
      ipxy = ipairr(jhcat(k,mx),jhcat(k,my))
      ipc  = ipairc(jhcat(k,mx),jhcat(k,my))
      c1 = eff(k,mc4) * cx(k,mx) * cx(k,my)

      tabvalx  &
        = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
        + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
        + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
        + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

      rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))
      rxfer(k,mx,mz) = rxfer(k,mx,mz) + rcx
      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + rcx * qx(k,mx)

      if (jnmb(mx) >= 5) then
         tabvaln  &
           = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
           + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
           + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
           + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

         colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
         enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))

   ! also loss for aerosol

      endif

   endif
   enddo
   return
end subroutine col1

!******************************************************************************

subroutine col2(m1,mx,my,mz,mc2,k1,k2,dn0,dtlt,i,j)

   use rconstants
   use micphys

   implicit none

   integer :: m1,mx,my,mz,mc2,k1,k2,i,j,k,jhcatx,jhcaty,ipxy,ipyx,ipc,it
   real :: c1,c2,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum0,colnum,rcoal  &
          ,qrcx,qrcy,qrcoal,qcoal,fracliq,tcoal,coalliq,coalice,area,cn13,cn24  &
          ,sip,rsip,qrsip,rfinlz,xtoz,dtlt

   real, dimension(m1) :: dn0

   real, dimension(15) ::  alpha,beta
   !            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
   data alpha /00.,00.,00., 1., 1., 1., 1.,00.,00.,00.,00., 1., 1., 1., 1./
   data beta  /00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,1.2,1.1,1.1,1.3/

   do k = k1,k2
   if(rx(k,mx) >= rxmin(mx) .and. rx(k,my) >= rxmin(my)) then
      jhcatx = jhcat(k,mx)
      jhcaty = jhcat(k,my)
      ipxy = ipairr(jhcatx,jhcaty)
      ipyx = ipairr(jhcaty,jhcatx)
      ipc  = ipairc(jhcatx,jhcaty)
      c2 = cx(k,mx) * cx(k,my)
      c1 = eff(k,mc2) * c2

      tabvalx  &
        = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
        + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
        + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
        + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

      rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

      tabvaly  &
        = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
        + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
        + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
        + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)

      rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

      if (jnmb(mx) >= 5 .or. jnmb(my) >= 5) then

         tabvaln  &
          = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
          + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
          + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
          + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

         colnum0 = c2 * colfacc(k) * 10. ** (-tabvaln)
         colnum = colnum0 * eff(k,mc2)
      else
      ![MLO - Not sure about that, but it may happen that we need colnum and the above
      !       if is the only place in which colnum is defined in this sub-routine, so
      !       I need to add something, here, so I put 0.
         colnum0=0.
         colnum=0.
      endif

      rcoal = rcx + rcy
      qrcx = rcx * qx(k,mx)
      qrcy = rcy * qx(k,my)
      qrcoal = qrcx + qrcy
      qcoal = qrcoal / (1.e-13 + rcoal)

      call qtc(qcoal,tcoal,fracliq)

      coalliq = rcoal * fracliq
      coalice = rcoal - coalliq

   ! secondary ice production: cn24 is the number fraction of collected cloud
   ! droplets larger than 24 microns and is obtained from an incomplete gamma
   ! function table.  cn13 is the fraction of collected cloud droplets
   ! smaller than 13 microns.  area is cross section area of collecting ice
   ! per m^3 of atmospheric volume.

      if (tcoal .gt. -8. .and. tcoal .lt. -3.) then

         area = cx(k,my) * dn0(k) * sipfac(jhcaty) * emb(k,my)  &
            ** (2.*pwmasi(jhcaty))
         it = nint(emb(k,mx) / emb1(1) * 5000.)
         cn13 = colnum * gamsip13(it) / (area * dtlt)
         cn24 = min(cx(k,mx)*dn0(k),colnum0) * gamsip24(it)
         sip = 9.1e-10 * cn24 * cn13 ** .93
         if (tcoal .lt. -5.) then
            sip = 0.33333 * (tcoal + 8.) * sip
         else
            sip = -0.5 * (tcoal + 3.) * sip
         endif

         rsip = sip * emb0(3) * dn0i(k)
         qrsip = qcoal * rsip

         rcoal = rcoal - rsip
         qrcoal = qrcoal - qrsip

         enxfer(k,mx,3) = enxfer(k,mx,3) + sip
         rxfer(k,mx,3) = rxfer(k,mx,3) + rsip
         qrxfer(k,mx,3) = qrxfer(k,mx,3) + qrsip

      endif

   ! ALWAYS NEED (ALPHA + BETA) .GE. 1 but in the (rare) case that
   ! fracliq may be a little larger than fracx due to collected
   ! liquid being above 0C, need (ALPHA + BETA) to be at least 1.1
   ! or 1.2, or need ALPHA itself to be at least 1.0.

      rfinlz = min(rcoal,  &
         alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

      xtoz = min(rcx,rfinlz)

      rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
      rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
      if (my .ne. mz) rxfer(k,my,mz) = rxfer(k,my,mz)  &
         + rfinlz - xtoz

      qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
      qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
      if (my .ne. mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
         + qx(k,my) * (rfinlz - xtoz)

      enxfer(k,mx,mx) = enxfer(k,mx,mx) + min(colnum,cx(k,mx))
      if (my .ne. mz) enxfer(k,my,mz) = enxfer(k,my,mz)  &
         + (rfinlz - xtoz) * min(colnum,cx(k,my)) / (1.e-20 + rcy)

   ! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X = Y

   ! also include loss of aerosol

   endif
   enddo
   return
end subroutine col2

!******************************************************************************

subroutine col3(m1,mx,my,mz,k1,k2,i,j)

   use micphys

   implicit none

   integer :: m1,mx,my,mz,k1,k2,i,j,k,ipxy,ipyx,ipc,jhcaty
   real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum,colnumx,colnumy,coalnum  &
          ,rcoal,qrcx,qrcy,qrcoal,qcoal,fracliq,coalliq,coalice,xtoz  &
          ,rfinlz,tcoal,cfinlz
   real, dimension(15) :: alpha,beta

   !            1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
   data alpha /00.,00., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1./
   data beta  /00.,00., 2., 2., 2., 1., 0., 2., 2., 2., 2., 2., 2., 2., 2./

   do k = k1,k2
   if(rx(k,mx) >= rxmin(mx) .and. rx(k,my) >= rxmin(my)) then
      jhcaty = jhcat(k,my)
      ipxy = ipairr(jhcat(k,mx),jhcaty)
      ipyx = ipairr(jhcaty,jhcat(k,mx))
      ipc  = ipairc(jhcat(k,mx),jhcaty)
      c1 = eff(k,1) * cx(k,mx) * cx(k,my)

      tabvalx  &
        = wct1(k,mx) * wct1(k,my) * coltabr (ict1(k,mx),ict1(k,my),ipxy)  &
        + wct2(k,mx) * wct1(k,my) * coltabr (ict2(k,mx),ict1(k,my),ipxy)  &
        + wct1(k,mx) * wct2(k,my) * coltabr (ict1(k,mx),ict2(k,my),ipxy)  &
        + wct2(k,mx) * wct2(k,my) * coltabr (ict2(k,mx),ict2(k,my),ipxy)

      rcx = min(rx(k,mx),c1 * colfacr(k) * 10. ** (-tabvalx))

      tabvaly  &
        = wct1(k,my) * wct1(k,mx) * coltabr (ict1(k,my),ict1(k,mx),ipyx)  &
        + wct2(k,my) * wct1(k,mx) * coltabr (ict2(k,my),ict1(k,mx),ipyx)  &
        + wct1(k,my) * wct2(k,mx) * coltabr (ict1(k,my),ict2(k,mx),ipyx)  &
        + wct2(k,my) * wct2(k,mx) * coltabr (ict2(k,my),ict2(k,mx),ipyx)

      rcy = min(rx(k,my),c1 * colfacr(k) * 10. ** (-tabvaly))

      if (jnmb(mx) >= 5) then
         tabvaln  &
          = wct1(k,mx) * wct1(k,my) * coltabc (ict1(k,mx),ict1(k,my),ipc)  &
          + wct2(k,mx) * wct1(k,my) * coltabc (ict2(k,mx),ict1(k,my),ipc)  &
          + wct1(k,mx) * wct2(k,my) * coltabc (ict1(k,mx),ict2(k,my),ipc)  &
          + wct2(k,mx) * wct2(k,my) * coltabc (ict2(k,mx),ict2(k,my),ipc)

         colnum = c1 * colfacc(k) * 10. ** (-tabvaln)
         colnumx = min(cx(k,mx),colnum)
         colnumy = min(cx(k,my),colnum)
         coalnum = min(colnumx,colnumy)
      endif

      rcoal = rcx + rcy
      qrcx = rcx * qx(k,mx)
      qrcy = rcy * qx(k,my)
      qrcoal = qrcx + qrcy
      qcoal = qrcoal / (1.e-20 + rcoal)

      call qtc(qcoal,tcoal,fracliq)

      coalliq = rcoal * fracliq
      coalice = rcoal - coalliq

      if (fracliq .ge. .99) then

         rxfer(k,my,mx) = rxfer(k,my,mx) + rcy
         qrxfer(k,my,mx) = qrxfer(k,my,mx) + qrcy
         if (jnmb(mx) >= 5)  &
            enxfer(k,my,my) = enxfer(k,my,my) + colnumy
      else

         rfinlz = min(rcoal,  &
            alpha(jhcaty) * coalliq + beta(jhcaty) * rcx)

         xtoz = min(rcx,rfinlz)

         rxfer(k,mx,mz) = rxfer(k,mx,mz) + xtoz
         rxfer(k,mx,my) = rxfer(k,mx,my) + rcx - xtoz
         if (my .ne. mz) rxfer(k,my,mz) = rxfer(k,my,mz)  &
            + rfinlz - xtoz

   ! NEED TO USE QCOAL TO TRANSFER Q?

         qrxfer(k,mx,mz) = qrxfer(k,mx,mz) + qx(k,mx) * xtoz
         qrxfer(k,mx,my) = qrxfer(k,mx,my) + qx(k,mx) * (rcx - xtoz)
         if (my .ne. mz) qrxfer(k,my,mz) = qrxfer(k,my,mz)  &
            + qx(k,my) * (rfinlz - xtoz)

         if (jnmb(mx) >= 5) then
            if (my .eq. mz) then
               enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
            elseif (colnumy .ge. colnumx) then
               cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
               enxfer(k,mx,mz) = enxfer(k,mx,mz) + cfinlz
               enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx - cfinlz
               enxfer(k,my,my) = enxfer(k,my,my) + colnumy
            else
               cfinlz = coalnum * rfinlz / (rcoal + 1.e-20)
               enxfer(k,my,mz) = enxfer(k,my,mz) + cfinlz
               enxfer(k,mx,mx) = enxfer(k,mx,mx) + colnumx
               enxfer(k,my,my) = enxfer(k,my,my) + colnumy - cfinlz
            endif
         endif

      endif
   endif
   enddo

   ! also include loss of aerosol

   return
end subroutine col3

!******************************************************************************

subroutine colxfers(m1,k1,k2,i,j,rloss,enloss)

   use micphys

   implicit none

   integer :: m1,i,j,k,lcat,kd1,kd2,jcat
   integer, dimension(10) :: k1,k2
   real, dimension(m1) :: rloss,enloss

   !  All rxfer values are nonnegative.

   do lcat = 1,7
      if (jnmb(lcat) .ge. 1) then
         kd1 = k1(lcat)
         kd2 = k2(lcat)

         do k = kd1,kd2
            rloss(k) = 0.
            enloss(k) = 0.
         enddo

         do jcat = 1,7
   ! change this to include enxfer of the same categories
            if (jnmb(jcat) .ge. 1) then
               if (lcat .ne. jcat) then
                  do k = kd1,kd2
                     rloss(k) = rloss(k) + rxfer(k,lcat,jcat)
                  enddo
               endif
               do k = kd1,kd2
                  enloss(k) = enloss(k) + enxfer(k,lcat,jcat)
               enddo
            endif
         enddo

         do k = kd1,kd2
            rloss(k) = min(1.,rx(k,lcat) / max(1.e-20,rloss(k)))
            enloss(k) = min(1.,cx(k,lcat) / max(1.e-10,enloss(k)))
         enddo

         do jcat = 1,7
            if (jnmb(jcat) .ge. 1) then
               if (lcat .ne. jcat) then
                  do k = kd1,kd2
                     rxfer(k,lcat,jcat) = rxfer(k,lcat,jcat)*rloss(k)
                     qrxfer(k,lcat,jcat)=qrxfer(k,lcat,jcat)*rloss(k)
                  enddo
               endif
               do k = kd1,kd2
                  enxfer(k,lcat,jcat) = enxfer(k,lcat,jcat)*enloss(k)
               enddo
            endif
         enddo
      endif
   enddo

   do lcat = 1,7

      if (jnmb(lcat) .ge. 1) then

         kd1 = k1(lcat)
         kd2 = k2(lcat)

         do jcat = 1,7
            if (jnmb(jcat) .ge. 1 .and. lcat .ne. jcat) then
               do k = kd1,kd2
                  rx(k,lcat) = rx(k,lcat) - rxfer(k,lcat,jcat)
                  rx(k,jcat) = rx(k,jcat) + rxfer(k,lcat,jcat)
                  qr(k,lcat) = qr(k,lcat) - qrxfer(k,lcat,jcat)
                  qr(k,jcat) = qr(k,jcat) + qrxfer(k,lcat,jcat)
                  cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,jcat)
                  cx(k,jcat) = cx(k,jcat) + enxfer(k,lcat,jcat)
               enddo
            endif
         enddo

         if (jnmb(lcat) >= 5) then
            do k = kd1,kd2
               cx(k,lcat) = cx(k,lcat) - enxfer(k,lcat,lcat)
            enddo
         endif

      endif
   enddo
   return
end subroutine colxfers


