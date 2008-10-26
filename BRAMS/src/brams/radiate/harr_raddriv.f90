!===============================================================================
! OLAM version 2.12  

! Copyright (C) 2002-2006; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================
!MLO - Adapted to ED-BRAMS 1.4
subroutine harr_raddriv(m1,m2,m3,nclouds,ifm,if_adap,time,deltat,ia,iz,ja,jz,nadd_rad    &
                       ,iswrtyp,ilwrtyp,icumfdbk,flpw,topt,glat,rtgt,pi0,pp,rho,theta,rv &
                       ,rshort,rlong,fthrd,rlongup,cosz                                  &
                       ,albedt,rshort_top,rshortup_top,rlongup_top,fthrd_lw              &
                       ,sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h                               &
                       ,con_c,con_r,con_p,con_s,con_a,con_g,con_h                        &
                       ,cuprliq,cuprice,mynum)

  use mem_harr,        only: mg, mb
  use mem_grid,        only: zm, zt
  use rconstants,      only: cpor, p00i, stefan, cp, cpi, p00, hr_sec
  use micphys,         only: ncat
  use mem_leaf,        only: isfcl

  implicit none
  integer                     , intent(in)    :: m1,m2,m3,nclouds,mynum,ifm,if_adap
  integer                     , intent(in)    :: ia,iz,ja,jz
  integer                     , intent(in)    :: nadd_rad,iswrtyp,ilwrtyp,icumfdbk
  real(kind=8)                , intent(in)    :: time
  real                        , intent(in)    :: deltat

  real, dimension(m1,m2,m3)   , intent(in)    :: pi0,pp,rho,theta,rv
  real, dimension(m1,m2,m3)   , intent(in)    :: sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h
  real, dimension(m1,m2,m3)   , intent(in)    :: con_c,con_r,con_p,con_s,con_a,con_g,con_h

  real, dimension(m2,m3)      , intent(in)    :: topt,glat,flpw,rtgt
  
  real, dimension(m1,m2,m3)   , intent(inout) :: fthrd,fthrd_lw
  real, dimension(m2,m3)      , intent(inout) :: rshort,rlong,rlongup,cosz,albedt &
                                                ,rshort_top,rshortup_top,rlongup_top

  !MLO variables cumulus feedback
  ! This is the liquid water that came from cumulus parameterization
  real, dimension(m1,m2,m3,nclouds), intent(in) :: cuprliq,cuprice

!------ Local arrays -----------------------------------------------------------------------!
  integer :: ka,nrad,koff,icld
  integer :: jhcat(m1,ncat)  ! hydrom category table with ice habits

  real, dimension(m1)       :: tairk   ! air temperature [K]
  real, dimension(m1)       :: rhov    ! vapor density [kg_vap/m^3]
  real, dimension(m1)       :: press   ! atmospheric pressure [Pa]
  real, dimension(m1,ncat)  :: rx      ! hydrom bulk spec dens [kg_hyd/kg_air]
  real, dimension(m1,ncat)  :: cx      ! hydrom bulk number [num_hyd/kg_air]
  real, dimension(m1,ncat)  :: embharr ! hydrom mean particle mass [kg/particle]

  integer, dimension(m1)       :: ns,nt ! Classes for moisture and temperature => habit table

  !Dimension nrad
  real, allocatable, dimension(:) :: rl  ! vapor density of all radiation levels (kg/m^3)
  real, allocatable, dimension(:) :: dzl ! delta-z (m) of all radiation levels
  real, allocatable, dimension(:) :: dl  ! air density of all radiation levels (kg/m^3)
  real, allocatable, dimension(:) :: pl  ! pressure (Pa)
  real, allocatable, dimension(:) :: o3l ! stores the calculated ozone profile (g/m^3)
  real, allocatable, dimension(:) :: vp  ! vapor pressure (Pa)

  !Dimension (nrad,3)
  real, allocatable, dimension(:,:) :: u ! path-length for gases (H_2O, CO_2, O_3)  (Pa)
  
  !Dimension(nrad,mb)
  real, allocatable, dimension(:,:) :: tp   ! optical depth of hydrometeors (m^-1)
  real, allocatable, dimension(:,:) :: omgp ! Single scatter albedo of hydrometeors
  real, allocatable, dimension(:,:) :: gp   ! Asymmetry factor of hydrometeors
  
  !Dimension(nrad)
  real, allocatable, dimension(:) :: zml   ! heights of W points of all radiation levels (m)
  real, allocatable, dimension(:) :: ztl   ! heights of T points of all radiation levels (m)
  real, allocatable, dimension(:) :: tl    ! temperature (K)
  real, allocatable, dimension(:) :: flxus ! Total upwelling s/w flux (W/m^2)
  real, allocatable, dimension(:) :: flxds ! Total downwelling s/w flux (W/m^2)
  real, allocatable, dimension(:) :: flxul ! Total upwelling l/w flux (W/m^2)
  real, allocatable, dimension(:) :: flxdl ! Total downwelling l/w flux (W/m^2)

  !Dimension(nrad,6)
  real, allocatable, dimension(:,:) :: fu  ! upwelling fluxes for pseudo-bands (W/m^2)
  real, allocatable, dimension(:,:) :: fd  ! downwelling fluxes for pseudo-bands (W/m^2)
  ! Set activation flags for gases of importance:  
  ! Flag = 1: gas active;    Flag = 0: gas not active

  integer, save :: ngass(mg)=(/1, 1, 1/)  ! Flags for (H2O, CO2, O3) for shortwave
  integer, save :: ngast(mg)=(/1, 1, 1/)  ! Flags for (H2O, CO2, O3) for longwave

  integer :: i,j,k,ib,ig,kk,ik,krad,mcat

  real :: rmix
  real :: dzl9,rvk0,rvk1
  logical :: first_with_ed

!MLO - Auxiliary variables for the stupid cloud 
  real, dimension(m1) :: rcl_parm,rpl_parm
  
  first_with_ed= time < dble(deltat) .and. isfcl == 5

  ! Copy surface and vertical-column values from model to radiation memory space
  ! In this loop, (k-koff) ranges from 2 to m1 + 1 - nint(flpw(i,j))
  do j=ja,jz
     do i=ia,iz
        ka=nint(flpw(i,j))
        koff=ka-2
        nrad=m1-1-koff+nadd_rad

        ![MLO - Allocating the variables that have nrad as a dimension
        allocate(rl(nrad),dzl(nrad),dl(nrad),pl(nrad),o3l(nrad),vp(nrad))
        allocate(u(nrad,3),tp(nrad,mb),omgp(nrad,mb),gp(nrad,mb))
        allocate(zml(nrad),ztl(nrad),tl(nrad))
        allocate(flxus(nrad),flxds(nrad),flxul(nrad),flxdl(nrad))
        allocate(fu(nrad,6),fd(nrad,6))
        !MLO]

        do k = ka,m1-1
           tairk(k)    = theta(k,i,j) * (pi0(k,i,j)+pp(k,i,j)) * cpi
           rhov(k)     = max(0.,rv(k,i,j)) * rho(k,i,j)
           press(k)    = p00 * (cpi * (pi0(k,i,j)+pp(k,i,j)) ) ** cpor
           dl(k-koff)  = rho(k,i,j)
           pl(k-koff)  = press(k)
           tl(k-koff)  = tairk(k)
           rl(k-koff)  = rhov (k)
        enddo
!----- Computing the heights, considering whether we are using terrain following or adaptive coordinate
        if (if_adap == 1) then
           do k=ka,m1-1
             zml(k-koff) = zm(k)
             ztl(k-koff) = zt(k)
           end do
        else
           do k=ka,m1-1
             zml(k-koff) = topt(i,j) + zm(k) * rtgt(i,j)
             ztl(k-koff) = topt(i,j) + zt(k) * rtgt(i,j)
           end do
        end if

!----- Now I fill the cloud droplets mixing ratio based on the tendencies. The idea is very simple
!      I will assume that a cloud existed wherever the moistening term was positive, and I assumed
!      a on/off cloud, and assuming that all contribution for moistening was once a cloud.
        rcl_parm = 0.
        rpl_parm = 0.
        do icld=1,nclouds
           do k=ka,m1-1
              rcl_parm(k) = rcl_parm(k)+ cuprliq(k,i,j,icld)
              rpl_parm(k) = rpl_parm(k)+ cuprice(k,i,j,icld)
           end do
        end do


        ! Fill surface values
        if (if_adap == 1) then
          zml(1) = zm(1+koff)
          ztl(1) = zt(1+koff)
        else
          zml(1) = topt(i,j) + zm(1) * rtgt(i,j)
          ztl(1) = topt(i,j) + zt(1) * rtgt(i,j)
        end if
        pl (1) = pl(2) + (zml(1) - ztl(2)) / (ztl(2) - ztl(3)) * (pl(2) - pl(3))

        call rad_mclat(m1,nrad,koff,glat(i,j),rtgt(i,j),dl,pl,rl,tl,o3l,zml,ztl,dzl)

!MLO - ED doesn't initialize rlongup before the radiation is called,
!      so I included extrapolation for the first time only.
        if (first_with_ed) then
          tl(1) = tl(2)
        else
          tl(1) = sqrt(sqrt(rlongup(i,j) / stefan))
        end if
        dl(1) = dl(2)
        rl(1) = rl(2)

        ! zero out scratch arrays

        u (:,:) = 0.0
        fu(:,:) = 0.0
        fd(:,:) = 0.0

        ! Fill arrays rx, cx, and emb with hydrometeor properties

        call cloudprep_rad(m1,ka,mcat,jhcat,tairk,rhov,rx,cx,embharr,ns,nt  &
             ,sh_c(1:m1,i,j),sh_r(1:m1,i,j),sh_p(1:m1,i,j),sh_s(1:m1,i,j)   &
             ,sh_a(1:m1,i,j),sh_g(1:m1,i,j),sh_h(1:m1,i,j)                  &
             ,con_c(1:m1,i,j),con_r(1:m1,i,j),con_p(1:m1,i,j)               &
             ,con_s(1:m1,i,j),con_a(1:m1,i,j),con_g(1:m1,i,j),con_h(1:m1,i,j))

        ! Fill hydrometeor optical property arrays [tp, omgp, gp]
        call cloud_opt(m1,ka,nrad,koff,mcat,icumfdbk,jhcat,dzl,rx           &
                      ,cx,embharr,tp,omgp,gp,sngl(time),rho(1:m1,i,j)       &
                      ,ns,nt,rcl_parm,rpl_parm,mynum)

        ! Get the path lengths for the various gases...

        call path_lengths(nrad,u,rl,dzl,dl,o3l,vp,pl)

        do k = 1,nrad
           if (rl(k) <   0. .or.   dl(k) <   0. .or.  &
               pl(k) <   0. .or.  o3l(k) <   0. .or.   tl(k) < 160.) then   
               ! tl(k) < 160.: This is -113 C, which is much colder than the 
               ! Vostok, Antarctica world record and should also
               ! be colder than any atmospheric temperature

              print*, 'Temperature too low or negative value of'
              print*, 'density, vapor, pressure, or ozone'
              print*, 'before calling Harrington radiation'
              print*, 'at grid: ',ifm,' at (k,i,j) = ',k,i,j
              print*, 'stopping model'
              print*, 'rad: k, rl(k), dl(k), pl(k), o3l(k), tl(k)'
              do kk=1,nrad
                print'(i3,5g15.6)', kk, rl(kk), dl(kk), pl(kk), o3l(kk), tl(kk)
              enddo
              stop '+++++ STOP! harr_raddriv (harr_raddriv.f90)'
           else
!              write(unit=76,fmt='(a,1x,i5,1x,5(a,1x,es13.6,1x))') &
!                'k=',k,'rl=',rl(k),'dl=',dl(k),'pl=',pl(k),'o3l=',o3l(k),'tl=',tl(k)
           endif

        enddo

        ! Harrington shortwave scheme (valid only if cosz > .03)


        if (iswrtyp == 3) then
! First I zero out rshort, which will remain zero if it is between dusk and dawn.
          rshort(i,j)       = 0.0
          rshort_top(i,j)   = 0.0
          rshortup_top(i,j) = 0.0

          if (cosz(i,j) > 0.03) then
             flxus(:) = 0.0
             flxds(:) = 0.0
             call harr_swrad(nrad,albedt(i,j),cosz(i,j),sngl(time),   &
                              u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxus,flxds,ngass,mynum)

              rshort(i,j) = flxds(1)
              rshort_top(i,j) = flxds(nrad)
              rshortup_top(i,j) = flxus(nrad)

              do k = ka,m1-1
                 krad = k - koff
                 fthrd(k,i,j) = fthrd(k,i,j)  &
                    + (flxds(krad) - flxds(krad-1) + flxus(krad-1) - flxus(krad))  &
                    / (dl(krad) * dzl(krad) * cp)
              end do
              !? rshort_diffuse(i,j) = rshort(i,j) * diffuse_fraction
           end if
        end if

        ! Harrington longwave scheme

        if (ilwrtyp == 3) then

           flxul(:) = 0.0
           flxdl(:) = 0.0

           call harr_lwrad(nrad,u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxul,flxdl,ngast,mynum)

           rlong(i,j) = flxdl(1)
           rlongup(i,j) = flxul(1)
           rlongup_top(i,j) = flxul(nrad)

           do k = ka,m1-1
              krad = k - koff
              fthrd(k,i,j) = fthrd(k,i,j)  &
                 + (flxdl(krad) - flxdl(krad-1) + flxul(krad-1) - flxul(krad))  &
                 / (dl(krad) * dzl(krad) * cp)
              fthrd_lw(k,i,j) = fthrd_lw(k,i,j)  &
                 + (flxdl(krad) - flxdl(krad-1) + flxul(krad-1) - flxul(krad))  &
                 / (dl(krad) * dzl(krad) * cp)
           end do

        end if

![MLO - Deallocating nrad dependent matrices
        deallocate(rl,dzl,dl,pl,o3l,vp)
        deallocate(u,tp,omgp,gp)
        deallocate(zml,ztl,tl,flxus,flxds,flxul,flxdl)
        deallocate(fu,fd)
!MLO]

     end do
  end do
  return
end subroutine harr_raddriv

!******************************************************************************

subroutine cloudprep_rad(m1,ka,mcat,jhcat,tairk,rhov,rx,cx,embharr        &
                        ,ns,nt,sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h                   &
                        ,con_c,con_r,con_p,con_s,con_a,con_g,con_h)

   ! This subroutine was developed from parts of subroutine MIC_COPY in
   ! omic_driv.f90 and subroutines EACH_COLUMN and ENEMB in omic_misc.f90.

   ! Arrays rx and cx are the bulk mass and number of hydrometeors PER KG OF AIR, 
   ! NOT PER M^3 AS IN THE ORIGINAL MICROPHYSICS SUBROUTINES.

   use rconstants, only: t00
   
   use micphys,    only: ncat,jnmb,jhabtab,emb0,emb1,emb2,rxmin,parm
   use therm_lib , only: level,rhovsil
   implicit none
   integer, intent(in)  :: m1
   integer, intent(in)  :: ka

   integer, intent(out) :: mcat            ! # of active hydrom categories (0,1, or 7)
   integer, intent(out) :: jhcat(m1,ncat)  ! hydrom category table with ice habits

   real, intent(in), dimension(m1)   :: tairk       ! air temperature [K]
   real, intent(in), dimension(m1)   :: rhov       ! water vapor density [kg_vap/m³]

   real, intent(out)    :: rx (m1,ncat)     ! hydrom bulk spec dens [kg_hyd/kg_air]
   real, intent(out)    :: cx (m1,ncat)     ! hydrom bulk number [num_hyd/kg_air]
   real, intent(out)    :: embharr(m1,ncat) ! hydrom mean particle mass [kg/particle]

   integer               :: k
   integer               :: icat
   integer               :: ihcat
   integer, intent(out), dimension(m1) :: ns ! Moisture class for habit table
   integer, intent(out), dimension(m1) :: nt ! Temperature class for habit table

   real                 :: rhovsilair
   real                 :: relhum
   real                 :: tairc
   real                 :: parmi

   real, dimension(m1), intent(in) :: sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h
   real, dimension(m1), intent(in) :: con_c,con_r,con_p,con_s,con_a,con_g,con_h

   real, parameter :: offset=1.e-16

   ! If level <= 1, there is no condensate of any type in this simulation.
   if (level <= 1) then
      ! Set mcat to 0 and return
      mcat = 0
      return
   endif

   ! If level = 2, cloud water is the only form of condensate that may exist. 
   if (level == 2) then
      ! Set mcat to 1
      mcat = 1
      ! Set jnmb flag for cloud water to 1
      jnmb(1) = 1
      ! In OLAM, with level = 2, cloud number concentration is specified in cparm
      ! or parm(1).  Diagnose cloud droplet mean mass.
      parmi = 1. / parm(1)
      do k = ka,m1-1
         rx(k,1) = sh_c(k)
         embharr(k,1) = rx(k,1) * parmi
         cx(k,1) = parm(1)
         jhcat(k,1) = 1
      end do
      return
   endif

   ! If level = 3, up to 7 forms of condensate that may exist.
   if (level == 3) then
      ! Set mcat to 7.
      mcat = 7
      ! Zero out microphysics scratch arrays for the present i,j column
      rx(:,:) = 0.
      cx(:,:) = 0.
      ! Copy hydrometeor bulk mass and number concentration from main model arrays
      ! to microphysics column arrays rx and cx

      ! Cloud water
      if (jnmb(1) >= 1) then
         do k = ka,m1-1
            if (sh_c(k) >= rxmin(1)) then
               ! If cloud bulk density is sufficiently abundant, copy to rx.
               rx(k,1) = sh_c(k)
               ! If cloud water number concentration is prognosed, copy to cx.
               if (jnmb(1) >= 5) cx(k,1) = con_c(k)
            end if
         end do
      end if
  
      ! Rain
      if (jnmb(2) >= 1) then
         do k = ka,m1-1
            if (sh_r(k) >= rxmin(2)) then
               ! If rain bulk density is sufficiently abundant, copy to rx,
               rx(k,2) = sh_r(k)
               ! If rain water number concentration is prognosed, copy to cx.
               if (jnmb(2) >= 5) cx(k,2) = con_r(k)
            end if
         end do
      end if

      ! Pristine ice
      if (jnmb(3) >= 1) then
         do k = ka,m1-1
            if (sh_p(k) >= rxmin(3)) then
               ! If pristine ice bulk density is sufficiently abundant, copy to rx.
               rx(k,3) = sh_p(k)
               ! If pristine ice number concentration is prognosed, copy to cx.
               if (jnmb(3) >= 5) cx(k,3) = con_p(k)
            end if
         end do
      end if

      ! Snow
      if (jnmb(4) >= 1) then
         do k = ka,m1-1
            if (sh_s(k) >= rxmin(4)) then
               ! If snow bulk density is sufficiently abundant, copy to rx.
               rx(k,4) = sh_s(k)
               ! If snow number concentration is prognosed, copy to cx.
               if (jnmb(4) >= 5) cx(k,4) = con_s(k)
            end if
         end do
      end if

      ! Aggregates
      if (jnmb(5) >= 1) then
         do k = ka,m1-1
            if (sh_a(k) >= rxmin(5)) then
               ! If aggregates bulk density is sufficiently abundant, copy to rx.
               rx(k,5) = sh_a(k)
               ! If aggregates number concentration is prognosed, copy to cx.
               if (jnmb(5) >= 5) cx(k,5) = con_a(k)
            end if
         end do
      end if

     ! Graupel
      if (jnmb(6) >= 1) then
         do k = ka,m1-1
            if (sh_g(k) >= rxmin(6)) then
               ! If graupel bulk density is sufficiently abundant, copy to rx,
               rx(k,6) = sh_g(k)
               ! If graupel number concentration is prognosed, copy to cx.
               if (jnmb(6) >= 5) cx(k,6) = con_g(k)
            end if
         end do
      end if

      ! Hail
      if (jnmb(7) >= 1) then
         do k = ka,m1-1
            if (sh_h(k) >= rxmin(7)) then
               ! If hail bulk density is sufficiently abundant, copy to rx,
               rx(k,7) = sh_h(k)
               ! If hail number concentration is prognosed, copy to cx.
               if (jnmb(7) >= 5) cx(k,7) = con_h(k)
            end if
         end do
      end if

      ! Diagnose pristine ice and snow habits from atmospheric temperature and humidity.
      ! This section of code copied or adapted from subroutines THRMSTR in omic_vap.f90 
      ! and EACH_COLUMN in omic_misc.f90
      do k = ka,m1-1
         tairc = tairk(k) - t00
         rhovsilair = rhovsil(tairk(k))
         relhum = min(1.,rhov(k) / rhovsilair)

         ns(k) = max(1,nint(100. * relhum))
         nt(k) = max(1,min(31,-nint(tairc)))

         jhcat(k,1) = 1
         jhcat(k,2) = 2
         jhcat(k,3) = jhabtab(nt(k),ns(k),1)
         jhcat(k,4) = jhabtab(nt(k),ns(k),2)
         jhcat(k,5) = 5
         jhcat(k,6) = 6
         jhcat(k,7) = 7
      enddo

      ! Loop over all hydrometeor categories
      do icat = 1,ncat
      ! Evaluate hydrometeor mean mass emb and concentration cx   
         
      ! This section of code was developed from subroutine enemb in omic_misc.f90
      ! by removing parts that are not needed for radiation calculations.  
      ! Arrays rx and cx are the bulk mass and number of hydrometeors PER KG OF AIR, 
      ! NOT PER M^3 AS IN THE ORIGINAL SUBROUTINE MIC_COPY.
         if (jnmb(icat) == 2) then

            do k = ka,m1-1
               ihcat = jhcat(k,icat)
               embharr(k,icat) = emb2(ihcat)
               cx(k,icat) = rx(k,icat) / max(offset,embharr(k,icat))
            enddo

         elseif (jnmb(icat) == 4) then

            parmi = 1. / max(parm(icat),offset)
            do k = ka,m1-1
               embharr(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat) * parmi))
               cx(k,icat) = rx(k,icat) / max(offset,embharr(k,icat))
            enddo

         elseif (jnmb(icat) >= 5) then

            do k = ka,m1-1
               embharr(k,icat) = max(emb0(icat),min(emb1(icat),rx(k,icat)  &
                           / max(rxmin(icat),cx(k,icat))))
               cx(k,icat) = rx(k,icat) / max(offset,embharr(k,icat))
            enddo

         endif

      enddo
      

   endif

   return
end subroutine cloudprep_rad

!******************************************************************************

subroutine cloud_opt(m1,ka,nrad,koff,mcat,icumfdbk,jhcat,dzl,rx,cx,embharr &
                    ,tp,omgp,gp,time,rho,ns,nt,rcl_parm,rpl_parm,mynum)

  use mem_harr, only: mb, nb, ocoef, bcoef, gcoef, nsolb
  use micphys, only: ncat, jnmb, pwmasi, dnfac,rxmin,cfmas,pwmas,emb0,emb1,gnu,jhabtab
  use rconstants, only: p00i, rocp, hr_sec

  ! computing properties of spherical liquid water and irregular ice
  ! using fits to adt theory
  !
  ! ib .......... band number
  ! mb .......... maximum number of bands
  ! nb .......... total number of bands
  ! mza.......... number of vertical levels
  ! dzl ......... delta z in each level (m)
  ! dn .......... characteristic diameter (m)
  ! embharr...... mean hydrometeor mass (kg)
  ! rx .......... hydrometeor mixing ratio (kg/kg)
  ! cx .......... hydrometeor concentration (#/kg)
  ! ocoef ....... scattering albedo fit coefficients
  ! bcoef ....... extinction fit coefficients
  ! gcoef ....... asymmetry fit coefficients
  ! ncog ........ number of fit coefficients (omega and asym)
  ! ncb ......... number of fit coefficients (extinction)
  ! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category
  !                 numbers as a function of 15 microphysics category numbers
  ! rho ......... model air density (kg/m^3)
  ! dnfac ....... factor for computing dn from emb
  ! pwmasi ...... inverse of power used in mass power law

  ! MLO - New variables included to describe the stupid cloud which was drawn from cumulus parameterization
  ! rcl_parm .... cloud droplets mixing ratio due to deep convection in the parameterized cumulus;

  implicit none

  integer, intent(in) :: m1
  integer, intent(in) :: ka
  integer, intent(in) :: nrad
  integer, intent(in) :: koff
  integer, intent(in) :: mcat
  integer, intent(in) :: icumfdbk

  real, intent(in) :: time

  integer, intent(in) :: jhcat(m1,ncat)

  real, intent(in) :: dzl(nrad)
  real, intent(in) :: rx(m1,ncat)
  real, intent(in) :: cx(m1,ncat)
  real, intent(in) :: embharr(m1,ncat)

  real, intent(out) :: tp(nrad,mb)   ! optical depth
  real, intent(out) :: omgp(nrad,mb) ! scattering albedo
  real, intent(out) :: gp(nrad,mb)   ! asymmetry parameter (assumes all 
                                       !   particles are spherical)
  real, dimension(m1), intent(in) :: rho

  ! Adding parameterized rain properties:
  integer, dimension(m1), intent(in) :: ns,nt    ! Moisture and temperature categories
  real   , dimension(m1), intent(in) :: rcl_parm
  real   , dimension(m1), intent(in) :: rpl_parm

  integer, intent(in) :: mynum


  integer ib,iz,krc
  integer icat,k,ihcat,krad

  real :: dn  ! hydrometeor characteristic diameter [microns]
  real :: ext
  real :: om
  real :: gg

  ! MLO : variables to account for the parameterized rain
  real, parameter :: parmi_cloud=1./.3e9                       ! This came from microphysics default count for cloud droplets [#/kg].
                                                               !    I only need the reciprocal, so I will only provide this instead.
  real, parameter :: parmi_prist=1./.1e4                       ! This came from microphysics default count for pristine ice   [#/kg].
                                                               !    This is rarely used since pristine ice is prognosed,  but it's 
                                                               !    probably okay for the cumulus cloud (I hope...). This will be 
                                                               !    used only at the radiation, and the cloud is parametrised anyway...

  real            :: emb2_cloud,cx_cloud,glg_cloud,glgm_cloud,dnfac_cloud
  real            :: emb2_prist,cx_prist,glg_prist,glgm_prist,dnfac_prist

  real, parameter :: offset=1.e-16
  real, external  :: gammln

  ! Use the following dn limiters (rather than those in microphysics) for 
  ! consistency with fit coefficients

  real, parameter, dimension(7) :: dnmin = (/   1.,   10.,   1.,  125.,   10.,   10.,   10./)
  real, parameter, dimension(7) :: dnmax = (/1000.,10000., 125.,10000.,10000.,10000.,10000./)

  ! Array kradcat maps RAMS/OLAM microphysics hydrometeor categories to those
  ! represented in Harrington radiation code according to the following numbering:

  !     Harrington radiation code             Microphysics
  ! ----------------------------------------------------------------
  !  1:   cloud drops                 1.  cloud drops
  !  2:   rain                        2.  rain
  !  3:   pristine ice columns        3.  pristine ice columns
  !  4:   pristine ice rosettes       4.  snow columns
  !  5:   pristine ice plates         5.  aggregates
  !  6:   snow columns                6.  graupel
  !  7:   snow rosettes               7.  hail
  !  8:   snow plates                 8.  pristine ice hexagonal plates
  !  9:   aggregates columns          9.  pristine ice dendrites
  !  10:  aggregates rosettes        10.  pristine ice needles
  !  11:  aggregates plates          11.  pristine ice rosettes
  !  12:  graupel                    12.  snow hexagonal plates
  !  13:  hail                       13.  snow dendrites
  !                                  14.  snow needles
  !                                  15.  snow rosettes

  integer, parameter, dimension(15) :: kradcat = (/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/)

  ! Initialize arrays to zero prior to summation over any hydrometeor species

  tp    (:,:) = 0.0
  omgp  (:,:) = 0.0
  gp    (:,:) = 0.0

  ! Loop over active (mcat) hydrometeor categories

  do icat = 1,mcat
     if (jnmb(icat) > 0) then

        do k = ka,m1-1
           krad = k - koff

           if (rx(k,icat) > rxmin(icat)) then

              ihcat = jhcat(k,icat)
              krc = kradcat(ihcat)
              dn = 1.e6 * dnfac(ihcat) * embharr(k,icat) ** pwmasi(ihcat)
              dn = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
              do ib = 1,nb
                 ext = cx(k,icat) * rho(k) * dzl(krad)  &
                    * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)
!                 write (unit=60+mynum,fmt='(5(a,1x,i4,1x),6(a,1x,es12.5,1x))') 'icat=',icat,'ihcat=',ihcat,'k=',k,'krc=',krc,'ib=',ib &
!                                   ,'cx=',cx(k,icat),'rho=',rho(k),'dzl=',dzl(k),'dn=',dn,'b0=',bcoef(1,ib,krc),'b2=',bcoef(2,ib,krc)

                 om = ocoef(1,ib,krc)  &
                    + ocoef(2,ib,krc) * exp(ocoef(3,ib,krc) * dn)  &
                    + ocoef(4,ib,krc) * exp(ocoef(5,ib,krc) * dn)

                 gg = gcoef(1,ib,icat)  &
                    + gcoef(2,ib,icat) * exp(gcoef(3,ib,icat) * dn)  &
                    + gcoef(4,ib,icat) * exp(gcoef(5,ib,icat) * dn)

                 if(ib <= nsolb)then
                    if(icat /= 1)then
                       gg = gg * 1.08
                    else
                       gg = gg * 1.13
                    endif
                 else
                    if(icat /= 1)then
                       gg = gg * 2.9
                    else
                       gg = gg * 2.1
                    endif
                 endif

                 tp(krad,ib)   = tp(krad,ib)   +           ext
                 omgp(krad,ib) = omgp(krad,ib) +      om * ext
                 gp(krad,ib)   = gp(krad,ib)   + gg * om * ext

              enddo
           endif
        enddo

     endif
  enddo
! MLO - Adding effects of parameterized rain. It is a very simple parameterization.
!       It assumes that all the parameterized condensed water was cloud droplet . It is crude, but 
!       it is just to give some guess...
  if (icumfdbk == 1) then
   !MLO - Now I add the parameterized liquid water

     do k=ka,m1-1
       krad=k-koff
       icat  = 1
       ihcat = 1
       krc   = kradcat(ihcat)
       if (rcl_parm(k) > rxmin(icat)) then
          emb2_cloud  = max(emb0(icat),min(emb1(icat),rcl_parm(k) * parmi_cloud))
          cx_cloud    = rcl_parm(k) / max(offset,emb2_cloud)
          glg_cloud   = gammln(gnu(icat))
          glgm_cloud  = gammln(gnu(icat) + pwmas(ihcat))
          dnfac_cloud = (exp(glg_cloud - glgm_cloud)/cfmas(ihcat)) ** pwmasi(ihcat)
          dn = 1.e6 * dnfac_cloud * emb2_cloud ** pwmasi(ihcat)
          dn = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
   ! Do the same I would do for the resolved clouds
          do ib = 1,nb
             ext = cx_cloud * rho(k) * dzl(krad)  &
                * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)

             om = ocoef(1,ib,krc)  &
                + ocoef(2,ib,krc) * exp(ocoef(3,ib,krc) * dn)  &
                + ocoef(4,ib,krc) * exp(ocoef(5,ib,krc) * dn)

             gg = gcoef(1,ib,krc)  &
                + gcoef(2,ib,krc) * exp(gcoef(3,ib,krc) * dn)  &
                + gcoef(4,ib,krc) * exp(gcoef(5,ib,krc) * dn)

             if(ib <= nsolb)then
                gg = gg * 1.13
             else
                gg = gg * 2.1
             endif

             tp(krad,ib)   = tp(krad,ib)   +           ext
             omgp(krad,ib) = omgp(krad,ib) +      om * ext
             gp(krad,ib)   = gp(krad,ib)   + gg * om * ext

          end do
        end if

       !----- Adding the parametrised ice as pristine ice. 
       icat  = 3
       if (rpl_parm(k) > rxmin(icat)) then
          ihcat = jhabtab(nt(k),ns(k),1)
          krc   = kradcat(ihcat)
          emb2_prist = max(emb0(icat),min(emb1(icat),rpl_parm(k) * parmi_prist))
          cx_prist = rpl_parm(k) / max(offset,emb2_prist)
          glg_prist = gammln(gnu(icat))
          glgm_prist = gammln(gnu(icat) + pwmas(icat))
          dnfac_prist = (exp(glg_prist - glgm_prist)/cfmas(ihcat)) ** pwmasi(ihcat)
          dn = 1.e6 * dnfac_prist * emb2_prist ** pwmasi(ihcat)
          dn = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
   ! Do the same I would do for the resolved clouds
          do ib = 1,nb
             ext = cx_prist * rho(k) * dzl(krad)  &
                * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)

             om = ocoef(1,ib,krc)  &
                + ocoef(2,ib,krc) * exp(ocoef(3,ib,krc) * dn)  &
                + ocoef(4,ib,krc) * exp(ocoef(5,ib,krc) * dn)

             gg = gcoef(1,ib,krc)  &
                + gcoef(2,ib,krc) * exp(gcoef(3,ib,krc) * dn)  &
                + gcoef(4,ib,krc) * exp(gcoef(5,ib,krc) * dn)

             if(ib <= nsolb)then
                gg = gg * 1.13
             else
                gg = gg * 2.1
             endif

             tp(krad,ib)   = tp(krad,ib)   +           ext
             omgp(krad,ib) = omgp(krad,ib) +      om * ext
             gp(krad,ib)   = gp(krad,ib)   + gg * om * ext

          end do
        end if


     end do
  end if
  
  ! Combine the optical properties....

  do ib = 1,nb
     do k = ka,m1-1
        krad = k - koff
        if (tp(krad,ib) > 0.0) then
           gp(krad,ib) = min(0.9999,gp(krad,ib) / omgp(krad,ib))
           omgp(krad,ib) = min(0.9999,omgp(krad,ib) / tp(krad,ib))
        else
           omgp(krad,ib) = 0.0
           gp(krad,ib) = 0.0
        endif
     enddo
  enddo

  return
end subroutine cloud_opt

!******************************************************************************

subroutine path_lengths(nrad,u,rl,dzl,dl,o3l,vp,pl)

  ! Get the path lengths for the various gases...

  use rconstants, only: g,ep

  implicit none

  integer :: nrad

  real, intent(out) :: u  (nrad,3)
  real, intent(out) :: vp (nrad)

  real, intent(in)  :: rl (nrad)
  real, intent(in)  :: dzl(nrad)
  real, intent(in)  :: dl (nrad)
  real, intent(in)  :: o3l(nrad)
  real, intent(in)  :: pl (nrad)

  real, parameter :: eps_rad = 1.e-15,rvmin=1.e-6

  real :: rvk0,rvk1,dzl9,rmix
  integer :: k
  real, parameter :: co2_mixing_ratio = 360.0e-6 * 44.011 / 28.966 ! [kg/kg]

  u(1,1) = .5 * (rl(2) + rl(1)) * g * dzl(1)
  u(1,2) = .5 * (dl(2) + dl(1)) * co2_mixing_ratio * g * dzl(1)
  u(1,3) = o3l(1) * g * dzl(1)

  rvk0=rl(1)
  do k = 2,nrad
     dzl9   = g * dzl(k)
     rmix = rl(k) / dl(k)
     vp(k)  = pl(k) * rmix / (ep + rmix)
     rvk1=(rl(k)+rvmin)
![MLO Trying the original formulation
!     u(k,1) = dzl9 * (rvk1-rvk0) / (log(rvk1 / rvk0) + eps_rad)
!     u(k,2) = dzl9 * (dl(k)-dl(k-1))/(log(dl(k)/dl(k-1))+eps_rad) * co2_mixing_ratio
!     u(k,3) = dzl9 * 0.5 * (o3l(k)+o3l(k-1))
!     rvk0=rvk1
!MLO]
     u(k,1) = 0.5 * dzl9 * (rl(k) + rl(k-1))
     u(k,2) = 0.5 * dzl9 * (dl(k) + dl(k-1)) * co2_mixing_ratio
     u(k,3) = 0.5 * dzl9 * (o3l(k) + o3l(k-1))
  enddo

  return
end subroutine path_lengths
