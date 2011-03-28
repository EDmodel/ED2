!==========================================================================================!
!==========================================================================================!
! BRAMS 4.0.6 - Based on OLAM version 2.12                                                 !
!                                                                                          !
! Copyright (C) 2002-2006; All Rights Reserved;                                            !
! Duke University, Durham, North Carolina, USA                                             !
!                                                                                          !
! Portions of this software are copied or derived from the RAMS software                   !
! package.  The following copyright notice pertains to RAMS and its derivatives,           !
! including OLAM:                                                                          !
!                                                                                          !
!     !-----------------------------------------------------------------------------!      !
!     ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; !      !
!     ! Colorado State University Research Foundation ; ATMET, LLC                  !      !
!     !                                                                             !      !
!     ! This software is free software; you can redistribute it and/or modify it    !      !
!     ! under the terms of the GNU General Public License as published by the Free  !      !
!     ! Software Foundation; either version 2 of the License, or (at your option)   !      !
!     ! any later version.                                                          !      !
!     !                                                                             !      !
!     ! This software is distributed in the hope that it will be useful, but        !      !
!     ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  !      !
!     ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    !      !
!     ! for more details.                                                           !      !
!     !                                                                             !      !
!     ! You should have received a copy of the GNU General Public License along     !      !
!     ! with this program; if not, write to the Free Software Foundation, Inc.,     !      !
!     ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA                     !      !
!     ! (http://www.gnu.org/licenses/gpl.html)                                      !      !
!     !-----------------------------------------------------------------------------!      !
!                                                                                          !
! It is requested that any scientific publications based on application of OLAM            !
! include the following acknowledgment:  "OLAM was developed at the                        !
! Edmund T. Pratt Jr. School of Engineering, Duke University."                             !
!                                                                                          !
! For additional information, including published references, please contact               !
! the software authors, Robert L. Walko (robert.walko@duke.edu)                            !
! or Roni Avissar (avissar@duke.edu).                                                      !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This is the main driver for Harrington et al. (2000) radiation scheme, called during  !
! the integration.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine harr_raddriv(m1,m2,m3,nclouds,ncrad,ifm,if_adap,time,deltat,ia,iz,ja,jz         &
                       ,nadd_rad,iswrtyp,ilwrtyp,icumfdbk,flpw,topt,glat,rtgt,pi0,pp,rho   &
                       ,theta,rv,co2p,rshort,rshort_diffuse,rlong,fthrd,rlongup,cosz       &
                       ,albedt,rshort_top,rshortup_top,rlongup_top,fthrd_lw                &
                       ,sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h                                 &
                       ,con_c,con_r,con_p,con_s,con_a,con_g,con_h                          &
                       ,cuprliq,cuprice,cuparea,cupierr,mynum)

   use mem_harr,        only: mg, mb, mpb
   use mem_grid,        only: zm, zt
   use rconstants,      only: cpor, p00i, stefan, cp, cpi, p00, hr_sec, toodry
   use micphys,         only: ncat,rxmin
   use mem_leaf,        only: isfcl
   use mem_radiate,     only: rad_cosz_min
   use harr_coms,       only: rl,dzl,dl,pl,co2l,o3l,vp,u,tp,omgp,gp,zml,ztl,tl             &
                             ,flxus,flxds,flxul,flxdl,fu,fd,zero_harr_met_scratch          &
                             ,zero_harr_flx_scratch,nradmax,tairk,rhoi,rhoe,rhov,press     &
                             ,rcl_parm,rpl_parm,area_parm,flx_diff
 
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                            , intent(in)    :: m1,m2,m3,nclouds,ncrad
   integer                            , intent(in)    :: mynum,ifm,if_adap
   integer                            , intent(in)    :: ia,iz,ja,jz
   integer                            , intent(in)    :: nadd_rad,iswrtyp,ilwrtyp,icumfdbk
   real(kind=8)                       , intent(in)    :: time
   real                               , intent(in)    :: deltat
   real, dimension(m2,m3)             , intent(in)    :: topt,glat,flpw,rtgt
   real, dimension(m1,m2,m3)          , intent(in)    :: pi0,pp,rho,theta,rv,co2p
   real, dimension(m1,m2,m3)          , intent(in)    :: sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h
   real, dimension(m1,m2,m3)          , intent(in)    :: con_c,con_r,con_p,con_s,con_a
   real, dimension(m1,m2,m3)          , intent(in)    :: con_g,con_h
   real, dimension(m2,m3,nclouds)     , intent(in)    :: cuparea
   real, dimension(m2,m3,nclouds)     , intent(in)    :: cupierr
   real, dimension(m1,m2,m3,nclouds)  , intent(in)    :: cuprliq
   real, dimension(m1,m2,m3,nclouds)  , intent(in)    :: cuprice
   real, dimension(m2,m3)             , intent(inout) :: rshort,rlong,rlongup,cosz,albedt
   real, dimension(m2,m3)             , intent(inout) :: rshort_top,rshortup_top
   real, dimension(m2,m3)             , intent(inout) :: rshort_diffuse
   real, dimension(m2,m3)             , intent(inout) :: rlongup_top
   real, dimension(m1,m2,m3)          , intent(inout) :: fthrd,fthrd_lw
   !------ Local arrays -------------------------------------------------------------------!
   integer                                            :: ka
   integer                                            :: nrad
   integer                                            :: koff
   integer                                            :: icld
   integer                                            :: i,j,k,ib,ig,kk,ik,krad,mcat
   real                                               :: area_csky
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Copy surface and vertical-column values from model to radiation memory space.  In  !
   ! this loop, (k-koff) ranges from 2 to m1 + 1 - nint(flpw(i,j)).                        !
   !---------------------------------------------------------------------------------------!
   jloop: do j=ja,jz
      iloop: do i=ia,iz

         !----- Flush all scratch arrays to zero ------------------------------------------!
         call zero_harr_met_scratch(m1,ncat,nradmax,mg,mb,mpb,ncrad)

         ka=nint(flpw(i,j))
         koff=ka-2
         nrad=m1-1-koff+nadd_rad


         do k = ka,m1-1
            tairk(k)     = theta(k,i,j) * (pi0(k,i,j)+pp(k,i,j)) * cpi
            rhoe(k)      = rho(k,i,j)
            rhov(k)      = max(toodry,rv(k,i,j)) * rhoe(k)
            rhoi(k)      = 1./rho(k,i,j)
            press(k)     = p00 * (cpi * (pi0(k,i,j)+pp(k,i,j)) ) ** cpor
            dl(k-koff)   = rho(k,i,j)
            pl(k-koff)   = press(k)
            tl(k-koff)   = tairk(k)
            rl(k-koff)   = rhov (k)
            co2l(k-koff) = co2p(k,i,j) * rhoe(k)
         end do

         !---------------------------------------------------------------------------------!
         !     Find the heights, considering whether we are using terrain following or     !
         ! adaptive coordinate.                                                            !
         !---------------------------------------------------------------------------------!
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


         !----- Fill surface values -------------------------------------------------------!
         if (if_adap == 1) then
           zml(1) = zm(1+koff)
           ztl(1) = zt(1+koff)
         else
           zml(1) = topt(i,j) + zm(1+koff) * rtgt(i,j)
           ztl(1) = topt(i,j) + zt(1+koff) * rtgt(i,j)
         end if
         pl(1) = pl(2) + (zml(1) - ztl(3)) / (ztl(2) - ztl(3)) * (pl(2) - pl(3))
         tl(1) = sqrt(sqrt(rlongup(i,j) / stefan))
         dl(1)   = dl(2)
         rl(1)   = rl(2)
         co2l(1) = co2l(2)

         !------ Filling upper levels if the domain doesn't go too far up. ----------------!
         call rad_mclat(m1,nrad,koff,glat(i,j),rtgt(i,j))

         !----- Fill arrays rxharr, cxharr, and embharr with hydrometeor properties. ------!
         call cloudprep_rad(m1,ka,mcat,sh_c(:,i,j),sh_r(:,i,j),sh_p(:,i,j),sh_s(:,i,j)     &
                           ,sh_a(:,i,j),sh_g(:,i,j),sh_h(:,i,j),con_c(:,i,j)               &
                           ,con_r(:,i,j),con_p(:,i,j),con_s(:,i,j),con_a(:,i,j)            &
                           ,con_g(:,i,j),con_h(:,i,j))

         !---------------------------------------------------------------------------------!
         !     Here we will account the parametrised clouds.  The way we consider it is by !
         ! calling the radiation solver several times, for each cloud and an additional    !
         ! time for non cumulus region.  In case the user doesn't want the cumulus feed-   !
         ! back, then the local nclouds is 0.                                              !
         !---------------------------------------------------------------------------------!
         area_csky = 1.
         cldloop: do icld=ncrad,0,-1
            !----- Flush all scratch arrays to zero ---------------------------------------!
            call zero_harr_flx_scratch(m1,ncat,nradmax,mg,mb,mpb,ncrad)

            !----- Assign the area of this cloud and the related parameters. --------------!
            if (icld == 0) then
               area_parm = area_csky
               do k=ka,m1-1
                  rcl_parm(k-koff) = 0.
                  rpl_parm(k-koff) = 0.
               end do
            else
               !----- No need to compute this cloud if it didn't happen. ------------------!
               if (cupierr(i,j,icld) /= 0.) cycle cldloop

               area_parm = cuparea(i,j,icld)
               area_csky = area_csky - area_parm
               do k=ka,m1-1
                  rcl_parm(k-koff) = cuprliq(k,i,j,icld)
                  rpl_parm(k-koff) = cuprice(k,i,j,icld)
               end do
            end if

            !----- Fill hydrometeor optical property arrays [tp, omgp, gp] ----------------!
            call cloud_opt(m1,ka,nrad,koff,mcat,icld,sngl(time),mynum)

            !----- Get the path lengths for the various gases... --------------------------!
            call path_lengths(nrad)

            !------------------------------------------------------------------------------!
            !    Sanity check. tl(k) < 160.: This is -113 C, which is much colder than the !
            ! Vostok, Antarctica world record and should also be colder than any atmo-     !
            ! spheric temperature.                                                         !
            !------------------------------------------------------------------------------!
            do k = 1,nrad
               if (rl(k)   <   0. .or.   dl(k) <   0. .or.   pl(k) <   0. .or.             &
                   co2l(k) <   0. .or.  o3l(k) <   0. .or.   tl(k) < 160.      ) then   

                   write (unit=*,fmt='(a)') '============================================='
                   write (unit=*,fmt='(a)') ' ERROR - harr_raddriv!!!'
                   write (unit=*,fmt='(a)') '         The model is about to stop!'
                   write (unit=*,fmt='(2(a,1x,i5,1x))') ' - Node:',mynum,' Grid: ',ifm
                   write (unit=*,fmt='(3(a,1x,i5,1x))') ' - k = ',k,' i = ',i,' j = ',j
                   write (unit=*,fmt='(a)') ' - Either the temperature is too low, or some'
                   write (unit=*,fmt='(a)') '   negative density, mixing ratio, '
                   write (unit=*,fmt='(a)') '   or pressure was detected!'
                   write (unit=*,fmt='(a)') ' - Sanity check at Harrington:'
                   write (unit=*,fmt='(a)') '---------------------------------------------'
                   write (unit=*,fmt='(a3,1x,6(a12,1x))') &
                      'LEV','  MIX. RATIO','     DENSITY','    PRESSURE','        CO_2'    &
                           ,'       OZONE',' TEMPERATURE'
                   do kk=1,nrad
                      write (unit=*,fmt='(i3,1x,6(es12.3,1x))')                            &
                                         kk,rl(kk),dl(kk),pl(kk),co2l(kk),o3l(kk),tl(kk)
                   enddo
                   write (unit=*,fmt='(a)') '---------------------------------------------'
                   write (unit=*,fmt='(a)') ' '
                   write (unit=*,fmt='(a)') '============================================='
                   call abort_run ('Weird thermodynamic values, caught at radiation'       &
                                  ,'harr_raddriv','harr_raddriv.f90')
               end if

            end do

            !------------------------------------------------------------------------------!
            !    Harrington shortwave scheme (valid only if cosz > cosz_min)               !
            !------------------------------------------------------------------------------!
            if (iswrtyp == 3) then
               !---------------------------------------------------------------------------!
               !    First I flush rshort, rshort_top and rshortup_top to zero, and they    !
               ! should remain zero if it is between dusk and dawn.                        !
               !---------------------------------------------------------------------------!
               if (cosz(i,j) > rad_cosz_min) then
                  call harr_swrad(nrad,albedt(i,j),cosz(i,j),sngl(time),mynum)
                  rshort      (i,j)   = rshort(i,j)         + flxds(1)    * area_parm
                  rshort_top  (i,j)   = rshort_top  (i,j)   + flxds(nrad) * area_parm
                  rshortup_top(i,j)   = rshortup_top(i,j)   + flxus(nrad) * area_parm
                  rshort_diffuse(i,j) = rshort_diffuse(i,j) + flx_diff    * area_parm
                  do k = ka,m1-1
                     krad = k - koff
                     fthrd(k,i,j) = fthrd(k,i,j)                                           &
                                  + ( (flxds(krad)-flxds(krad-1)                           &
                                    + flxus(krad-1)-flxus(krad))                           &
                                    / (dl(krad) * dzl(krad) * cp)) * area_parm
                  end do
               end if
            end if
            !------------------------------------------------------------------------------!
            !    Harrington longwave scheme.                                               !
            !------------------------------------------------------------------------------!
            if (ilwrtyp == 3) then

               call harr_lwrad(nrad,mynum)

               rlong      (i,j) = rlong      (i,j) + flxdl(1)    * area_parm
               rlongup    (i,j) = rlongup    (i,j) + flxul(1)    * area_parm
               rlongup_top(i,j) = rlongup_top(i,j) + flxul(nrad) * area_parm

               do k = ka,m1-1
                  krad = k - koff
                  fthrd(k,i,j)    = fthrd(k,i,j)                                           &
                                  + ( (flxdl(krad)-flxdl(krad-1)                           &
                                    + flxul(krad-1)-flxul(krad))                           &
                                    / (dl(krad) * dzl(krad) * cp)) * area_parm
                  fthrd_lw(k,i,j) = fthrd_lw(k,i,j)                                        &
                                  + ( (flxdl(krad)-flxdl(krad-1)                           &
                                    + flxul(krad-1)-flxul(krad))                           &
                                    / (dl(krad) * dzl(krad) * cp)) * area_parm
               end do
            end if
         end do cldloop
         !---------------------------------------------------------------------------------!

         if (ilwrtyp == 3 .and. rlong(i,j) > 600.) then
            call print_longwave(mynum,nrad,m1,i,j,time,rlong(i,j),rlongup(i,j))
            call abort_run('Non-sense radiation','harr_raddriv','harr_raddriv.f90')
         end if

      end do iloop
   end do jloop

   return
end subroutine harr_raddriv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine determines the hydrometeor basic radiation properties. This subrout-  !
! ine was developed from parts of subroutine MIC_COPY in mic_driv.f90 and subroutines      !
! EACH_COLUMN and ENEMB in mic_misc.f90.                                                   !
!------------------------------------------------------------------------------------------!
subroutine cloudprep_rad(m1,ka,mcat,sh_c,sh_r,sh_p,sh_s,sh_a,sh_g,sh_h,con_c,con_r,con_p   &
                        ,con_s,con_a,con_g,con_h)
   use rconstants , only : t00       ! ! intent(in)
   use micphys    , only : ncat      & ! intent(in)
                          ,jnmb      & ! intent(in)
                          ,jhabtab   & ! intent(in)
                          ,availcat  & ! intent(in)
                          ,progncat  & ! intent(in)
                          ,emb0      & ! intent(in)
                          ,emb1      & ! intent(in)
                          ,emb2      & ! intent(in)
                          ,rxmin     & ! intent(in)
                          ,parm      & ! intent(in)
                          ,cfemb0    & ! intent(in)
                          ,cfen0     & ! intent(in)
                          ,pwemb0    & ! intent(in)
                          ,pwen0     ! ! intent(in)
   use therm_lib  , only : level     & ! intent(in)
                          ,rhovsil   ! ! function
   use harr_coms  , only : jhcatharr & ! intent(out) - hydrom cat. table with ice habits
                          ,tairk     & ! intent(in)  - Air temperature             [     K]
                          ,rhov      & ! intent(in)  - Water vapour density        [ kg/m³]
                          ,rhoe      & ! intent(in)  - Air density                 [ kg/m³]
                          ,rhoi      & ! intent(in)  - Specific volume             [ m³/kg]
                          ,rxharr    & ! intent(out) - hydrom bulk spec dens       [ kg/kg]
                          ,cxharr    & ! intent(out) - hydrom bulk number          [  #/kg]
                          ,embharr   & ! intent(out) - hydrom mean particle mass   [  kg/#]
                          ,nsharr    & ! intent(out) - Moisture class for habit table
                          ,ntharr    ! ! intent(out) - Temperature class for habit table
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                    , intent(in)  :: m1
   integer                    , intent(in)  :: ka
   real   , dimension(m1)     , intent(in)  :: sh_c    ! Cloud mixing ratio        [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: sh_r    ! Rain mixing ratio         [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: sh_p    ! Pristine ice mixing ratio [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: sh_s    ! Snow mixing ratio         [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: sh_a    ! Aggregates mixing ratio   [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: sh_g    ! Graupel mixing ratio      [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: sh_h    ! Hail mixing ratio         [ kg/kg]
   real   , dimension(m1)     , intent(in)  :: con_c   ! Cloud concentration       [  #/kg]
   real   , dimension(m1)     , intent(in)  :: con_r   ! Rain concentration        [  #/kg]
   real   , dimension(m1)     , intent(in)  :: con_p   ! Pristine ice concent.     [  #/kg]
   real   , dimension(m1)     , intent(in)  :: con_s   ! Snow concentration        [  #/kg]
   real   , dimension(m1)     , intent(in)  :: con_a   ! Aggregates concentration  [  #/kg]
   real   , dimension(m1)     , intent(in)  :: con_g   ! Graupel concentration     [  #/kg]
   real   , dimension(m1)     , intent(in)  :: con_h   ! Hail concentration        [  #/kg]
   integer                    , intent(out) :: mcat    ! # of active hydrom categories
   !----- Local variables -----------------------------------------------------------------!
   integer                                  :: k,icat,ihcat
   real                                     :: rhovsilair,relhum,tairc,parmi
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!

   select case (level)
   case (0,1)
      !------------------------------------------------------------------------------------!
      !    If level <= 1, there is no condensate of any type in this simulation.           !
      !  Set mcat to 0 and return                                                          !
      !------------------------------------------------------------------------------------!
      mcat = 0
      return
   
   case (2)
      !------------------------------------------------------------------------------------!
      !    If level = 2, cloud water is the only form of condensate that may exist.        !
      ! Set mcat to 1                                                                      !
      !------------------------------------------------------------------------------------!
      mcat = 1
      !----- Set jnmb flag for cloud water to 1. ------------------------------------------!
      jnmb(1)     = 1
      availcat(1) = .true.
      progncat(1) = .false.
      !------------------------------------------------------------------------------------!
      !     In OLAM, with level = 2, cloud number concentration is specified in cparm      !
      ! or parm(1).  Diagnose cloud droplet mean mass.                                     !
      !------------------------------------------------------------------------------------!
      parmi = 1. / parm(1)
      do k = ka,m1-1
         rxharr(k,1)    = sh_c(k)
         embharr(k,1)   = rxharr(k,1) * parmi
         cxharr(k,1)    = parm(1)
         jhcatharr(k,1) = 1
      end do
      return
   
   case (3)
      !------------------------------------------------------------------------------------!
      !    If level = 3, up to 7 forms of condensate that may exist.                       !
      ! Set mcat to 7.                                                                     !
      !------------------------------------------------------------------------------------!
      mcat = 7


      !------------------------------------------------------------------------------------!
      !    Copy hydrometeor bulk mass and number concentration from main model arrays      !
      ! to microphysics column arrays rx and cx                                            !
      !------------------------------------------------------------------------------------!

      !----- Cloud water ------------------------------------------------------------------!
      if (availcat(1)) then
         do k = ka,m1-1
            if (sh_c(k) >= rxmin(1)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,1) = sh_c(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(1)) cxharr(k,1) = con_c(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------ Rain ------------------------------------------------------------------------!
      if (availcat(2)) then
         do k = ka,m1-1
            if (sh_r(k) >= rxmin(2)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,2) = sh_r(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(2)) cxharr(k,2) = con_r(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !----- Pristine ice -----------------------------------------------------------------!
      if (availcat(3)) then
         do k = ka,m1-1
            if (sh_p(k) >= rxmin(3)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,3) = sh_p(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(3)) cxharr(k,3) = con_p(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !----- Snow -------------------------------------------------------------------------!
      if (availcat(4)) then
         do k = ka,m1-1
            if (sh_s(k) >= rxmin(4)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,4) = sh_s(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(4)) cxharr(k,4) = con_s(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !----- Aggregates -------------------------------------------------------------------!
      if (availcat(5)) then
         do k = ka,m1-1
            if (sh_a(k) >= rxmin(5)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,5) = sh_a(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(5)) cxharr(k,5) = con_a(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !----- Graupel ----------------------------------------------------------------------!
      if (availcat(6)) then
         do k = ka,m1-1
            if (sh_g(k) >= rxmin(6)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,6) = sh_g(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(6)) cxharr(k,6) = con_g(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !----- Hail -------------------------------------------------------------------------!
      if (availcat(7)) then
         do k = ka,m1-1
            if (sh_h(k) >= rxmin(7)) then
               !----- If bulk density is sufficiently abundant, copy to rx. ---------------!
               rxharr(k,7) = sh_h(k)
               !----- If number concentration is prognosed, copy to cx. -------------------!
               if (progncat(7)) cxharr(k,7) = con_h(k)
            end if
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Diagnose pristine ice and snow habits from atmospheric temperature and          !
      ! humidity.  This section of code copied or adapted from subroutines THRMSTR in      !
      ! mic_vap.f90 and EACH_COLUMN in mic_misc.f90.                                       !
      !------------------------------------------------------------------------------------!
      do k = ka,m1-1
         tairc = tairk(k) - t00
         rhovsilair = rhovsil(tairk(k))
         relhum = min(1.,rhov(k) / rhovsilair)

         nsharr(k) = max(1,nint(100. * relhum))
         ntharr(k) = max(1,min(31,-nint(tairc)))

         jhcatharr(k,1) = 1
         jhcatharr(k,2) = 2
         jhcatharr(k,3) = jhabtab(ntharr(k),nsharr(k),1)
         jhcatharr(k,4) = jhabtab(ntharr(k),nsharr(k),2)
         jhcatharr(k,5) = 5
         jhcatharr(k,6) = 6
         jhcatharr(k,7) = 7
      end do
      !------------------------------------------------------------------------------------!



      !----- Loop over all hydrometeor categories -----------------------------------------!
      do icat = 1,ncat
         !---------------------------------------------------------------------------------!
         !    Evaluate hydrometeor mean mass emb and concentration cx. This section of     !
         ! code was developed from subroutine enemb in mic_misc.f90 by removing parts that !
         ! are not needed for radiation calculations.                                      !
         !---------------------------------------------------------------------------------!
         select case (jnmb(icat))
         case (2)
            do k = ka,m1-1
               ihcat = jhcatharr(k,icat)
               if (rxharr(k,icat) >= rxmin(icat)) then
                  embharr(k,icat) = emb2(ihcat)
                  cxharr(k,icat)  = rxharr(k,icat) / embharr(k,icat)
               end if
            end do

         case (3)
            do k = ka,m1-1
               ihcat       = jhcatharr(k,icat)
               if (rxharr(k,icat) >= rxmin(icat)) then
                  embharr(k,icat) = cfemb0(ihcat)                                          &
                                  * (rhoe(k) * rxharr(k,icat)) ** pwemb0(ihcat)
                  cxharr(k,icat)  = cfen0(ihcat) * rhoi(k)                                 &
                                  * (rhoe(k) * rxharr(k,icat)) ** pwen0(ihcat)
               end if
            end do

         case (4)
            parmi = 1. / parm(icat)
            do k = ka,m1-1
               if (rxharr(k,icat) >= rxmin(icat)) then
                  embharr(k,icat) = max(emb0(icat),min(emb1(icat),rxharr(k,icat) * parmi))
                  cxharr(k,icat)  = rxharr(k,icat) / embharr(k,icat)
               end if
            end do

         case (5:7)
            do k = ka,m1-1
               if (rxharr(k,icat) >= rxmin(icat)) then
                  embharr(k,icat) = max(emb0(icat)                                         &
                                       ,min(emb1(icat)                                     &
                                       ,rxharr(k,icat)/max(1.e-9,cxharr(k,icat))))
                  cxharr(k,icat) = rxharr(k,icat) / embharr(k,icat)
               end if
            end do
         end select
      end do
   end select

   return
end subroutine cloudprep_rad
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute properties of spherical liquid water and irregular ice   !
! using fits to adt theory. It will also provide a crude first guess of properties of      !
! parametrised clouds, assuming all liquid is cloud droplets and all ice is pristine ice.  !
!------------------------------------------------------------------------------------------!
subroutine cloud_opt(m1,ka,nrad,koff,mcat,icld,time,mynum)

   use mem_harr   , only : mb           & ! intent(in)
                         , nb           & ! intent(in)
                         , ocoef        & ! intent(in)
                         , bcoef        & ! intent(in)
                         , gcoef        & ! intent(in)
                         , nsolb        ! ! intent(in)
   use micphys    , only : ncat         & ! intent(in)
                         , jnmb         & ! intent(in)
                         , pwmasi       & ! intent(in)
                         , dnfac        & ! intent(in)
                         , rxmin        & ! intent(in)
                         , cfmas        & ! intent(in)
                         , pwmas        & ! intent(in)
                         , emb0         & ! intent(in)
                         , emb1         & ! intent(in)
                         , emb2         & ! intent(in)
                         , gnu          & ! intent(in) 
                         , jhabtab      & ! intent(in)
                         , parm         & ! intent(in) 
                         , availcat     ! ! intent(in)
   use rconstants , only : p00i         & ! intent(in)
                         , rocp         & ! intent(in) 
                         , hr_sec       & ! intent(in)
                         , lnexp_min    ! ! intent(in)
   use harr_coms  , only : jhcatharr    & ! intent(in)
                          ,dzl          & ! intent(in)
                          ,dl           & ! intent(in)
                          ,rxharr       & ! intent(in)
                          ,cxharr       & ! intent(in)
                          ,embharr      & ! intent(in)
                          ,tp           & ! intent(in)
                          ,omgp         & ! intent(in)
                          ,gp           & ! intent(in)
                          ,rhoe         & ! intent(in)
                          ,nsharr       & ! intent(out)
                          ,ntharr       & ! intent(out)
                          ,rcl_parm     & ! intent(out)
                          ,rpl_parm     & ! intent(out)
                          ,parmi_cloud  & ! intent(in)
                          ,parmi_prist  & ! intent(in)
                          ,dnmin        & ! intent(in)
                          ,dnmax        & ! intent(in)
                          ,sacoef       & ! intent(in)
                          ,lacoef       & ! intent(in)
                          ,kradcat      ! ! intent(in)

   !---------------------------------------------------------------------------------------!
   ! TABLE OF CONTENTS:                                                                    !
   !                                                                                       !
   ! ib .......... band number                                                             !
   ! mb .......... maximum number of bands                                                 !
   ! nb .......... total number of bands                                                   !
   ! mza.......... number of vertical levels                                               !
   ! dzl ......... delta z in each level (m)                                               !
   ! dn .......... characteristic diameter (m)                                             !
   ! embharr...... mean hydrometeor mass (kg)                                              !
   ! rx .......... hydrometeor mixing ratio (kg/kg)                                        !
   ! cx .......... hydrometeor concentration (#/kg)                                        !
   ! ocoef ....... scattering albedo fit coefficients                                      !
   ! bcoef ....... extinction fit coefficients                                             !
   ! gcoef ....... asymmetry fit coefficients                                              !
   ! ncog ........ number of fit coefficients (omega and asym)                             !
   ! ncb ......... number of fit coefficients (extinction)                                 !
   ! kradcat ..... cross-reference table giving Jerry's 13 hydrometeor category            !
   !                 numbers as a function of 15 microphysics category numbers             !
   ! rhoe ........ model air density (kg/m^3)                                              !
   ! dnfac ....... factor for computing dn from emb                                        !
   ! pwmasi ...... inverse of power used in mass power law                                 !
   ! rcl_parm .... cloud droplets mixing ratio due to parametrised cumuli                  !
   ! rpl_parm .... pristine ice mixing ratio due to parametrised cumuli.                   !
   ! dn .......... hydrometeor characteristic diameter.                                    !
   ! nt .......... Temperature category for ice habits.                                    !
   ! ns .......... Moisture category for ice habits.                                       !
   !                                                                                       !
   ! tp .......... optical depth                                                           !
   ! omgp ........ scattering albedo                                                       !
   ! gp .......... assymetry parameter (assumes all particles are spherical)               !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)  :: m1,ka,nrad,koff,mcat,icld,mynum
   real                          , intent(in)  :: time
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: ib,iz,krc,icat,k,ihcat,krad
   real                                        :: dn,ext,om,gg
   real                                        :: emb_cloud,cx_cloud,glg_cloud
   real                                        :: glgm_cloud,dnfac_cloud
   real                                        :: emb_prist,cx_prist,glg_prist
   real                                        :: glgm_prist,dnfac_prist
   !----- Functions -----------------------------------------------------------------------!
   real                          , external    :: gammln
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Bulk microphysics first: loop over active (mcat) hydrometeor categories.          !
   !---------------------------------------------------------------------------------------!
   do icat = 1,mcat
      if (availcat(icat)) then

         do k = ka,m1-1
            krad = k - koff

            if (rxharr(k,icat) > rxmin(icat)) then
               ihcat = jhcatharr(k,icat)
               krc = kradcat(ihcat)
               dn  = 1.e6 * dnfac(ihcat) * embharr(k,icat) ** pwmasi(ihcat)
               dn  = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
               do ib = 1,nb
                  ext = cxharr(k,icat) * rhoe(k) * dzl(krad)                               &
                      * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)


                  om = ocoef(1,ib,krc)                                                     &
                     + ocoef(2,ib,krc) * exp(max(-60.0,ocoef(3,ib,krc)*dn))                &
                     + ocoef(4,ib,krc) * exp(max(-60.0,ocoef(5,ib,krc)*dn))

                  gg = gcoef(1,ib,icat)                                                    &
                     + gcoef(2,ib,icat) * exp(max(-60.0,gcoef(3,ib,icat)*dn))              &
                     + gcoef(4,ib,icat) * exp(max(-60.0,gcoef(5,ib,icat)*dn))

                  if (ib <= nsolb) then
                     gg = gg * sacoef(icat)
                  else
                     gg = gg * lacoef(icat)
                  end if

                  tp(krad,ib)   = tp(krad,ib)   +           ext
                  omgp(krad,ib) = omgp(krad,ib) +      om * ext
                  gp(krad,ib)   = gp(krad,ib)   + gg * om * ext
               end do
            end if
         end do
      end if
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Adding effects of parameterized rain. It is a very simple parameterization. It    !
   ! assumes that all the parameterized condensed water was cloud droplet or pristine ice. !
   ! It is really a 1st guess, don't expect wonderful results...                           !
   !---------------------------------------------------------------------------------------!
   if (icld /= 0) then

      !----- Liquid water -----------------------------------------------------------------!
      do k=ka,m1-1
         krad=k-koff
         icat  = 1
         ihcat = 1
         krc   = kradcat(ihcat)
         if (rcl_parm(k) > rxmin(icat)) then
            !----- Finding concentration assuming known diameter --------------------------!
            emb_cloud   = max(emb0(icat),min(emb1(icat),rcl_parm(k) * parmi_cloud))
            cx_cloud    = rcl_parm(k) / emb_cloud
            glg_cloud   = gammln(gnu(icat))
            glgm_cloud  = gammln(gnu(icat) + pwmas(ihcat))
            dnfac_cloud = (exp(glg_cloud - glgm_cloud)/cfmas(ihcat)) ** pwmasi(ihcat)
            dn          = 1.e6 * dnfac_cloud * emb_cloud ** pwmasi(ihcat)
            dn          = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
            !----- Now we find the optical properties as if it were resolved cloud. -------!
            do ib = 1,nb
               ext = cx_cloud * rhoe(k) * dzl(krad)                                        &
                  * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)

               om = ocoef(1,ib,krc)                                                        &
                  + ocoef(2,ib,krc) * exp(ocoef(3,ib,krc) * dn)                            &
                  + ocoef(4,ib,krc) * exp(ocoef(5,ib,krc) * dn)

               gg = gcoef(1,ib,icat)                                                       &
                  + gcoef(2,ib,icat) * exp(gcoef(3,ib,icat) * dn)                          &
                  + gcoef(4,ib,icat) * exp(gcoef(5,ib,icat) * dn)
               if(ib <= nsolb)then
                  gg = gg * sacoef(icat)
               else
                  gg = gg * lacoef(icat)
               endif

               tp(krad,ib)   = tp(krad,ib)   +           ext
               omgp(krad,ib) = omgp(krad,ib) +      om * ext
               gp(krad,ib)   = gp(krad,ib)   + gg * om * ext

            end do
         end if

         !----- Adding the parametrised ice as pristine ice -------------------------------!
         icat  = 3
         if (rpl_parm(k) > rxmin(icat)) then
            ihcat = jhabtab(ntharr(k),nsharr(k),1)
            krc   = kradcat(ihcat)
            emb_prist   = max(emb0(icat),min(emb1(icat),rpl_parm(k) * parmi_prist))
            cx_prist    = rpl_parm(k) / emb_prist
            glg_prist   = gammln(gnu(icat))
            glgm_prist  = gammln(gnu(icat) + pwmas(icat))
            dnfac_prist = (exp(glg_prist - glgm_prist)/cfmas(ihcat)) ** pwmasi(ihcat)
            dn = 1.e6 * dnfac_prist * emb_prist ** pwmasi(ihcat)
            dn = max(dnmin(icat),min(dnmax(icat),dn))  ! dn units are microns
            !----- Now we find the optical properties as if it were resolved cloud. -------!
            do ib = 1,nb
               ext = cx_prist * rhoe(k) * dzl(krad)                                        &
                  * bcoef(1,ib,krc) * dn ** bcoef(2,ib,krc)

               om = ocoef(1,ib,krc)                                                        &
                  + ocoef(2,ib,krc) * exp(max(-60.0,ocoef(3,ib,krc)*dn))                   &
                  + ocoef(4,ib,krc) * exp(max(-60.0,ocoef(5,ib,krc)*dn))

               gg = gcoef(1,ib,icat)                                                       &
                  + gcoef(2,ib,icat) * exp(max(-60.0,gcoef(3,ib,icat)*dn))                 &
                  + gcoef(4,ib,icat) * exp(max(-60.0,gcoef(5,ib,icat)*dn))

               if (ib <= nsolb) then
                  gg = gg * sacoef(icat)
               else
                  gg = gg * lacoef(icat)
               end if

               tp(krad,ib)   = tp(krad,ib)   +           ext
               omgp(krad,ib) = omgp(krad,ib) +      om * ext
               gp(krad,ib)   = gp(krad,ib)   + gg * om * ext
            end do
         end if
      end do
   end if
   
   !----- Combine the optical properties.... ----------------------------------------------!
   do ib = 1,nb
      do k = ka,m1-1
         krad = k - koff
         if (tp(krad,ib) > 0.0) then
            gp(krad,ib) = min(0.9999,gp(krad,ib) / omgp(krad,ib))
            omgp(krad,ib) = min(0.9999,omgp(krad,ib) / tp(krad,ib))
         else
            omgp(krad,ib) = 0.0
            gp(krad,ib) = 0.0
         end if
      end do
   end do

   return
end subroutine cloud_opt
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Get the path lengths for the various gases...                                         !
!------------------------------------------------------------------------------------------!
subroutine path_lengths(nrad)
   use rconstants , only : grav  & ! intent(in)
                         , ep    ! ! intent(in)
   use mem_harr   , only : mg    ! ! intent(in)
   use harr_coms  , only : u     & ! intent(out)
                         , rl    & ! intent(in)
                         , dzl   & ! intent(in)
                         , dl    & ! intent(in)
                         , co2l  & ! intent(in)
                         , o3l   & ! intent(in)
                         , vp    & ! intent(out)
                         , pl    ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)    :: nrad
   !----- Local variables -----------------------------------------------------------------!
   real                                        :: dzl9
   real                                        :: rh2obot,rh2otop
   real                                        :: rco2bot,rco2top
   real                                        :: ro3bot , ro3top
   integer                                     :: k
   !----- Constants  ----------------------------------------------------------------------!
   real, parameter :: eps_rad = 1.e-15
   real, parameter :: rvmin=1.e-6
   !---------------------------------------------------------------------------------------!

   rh2otop = max(rl(1) / dl(1), rvmin)
   rco2top = co2l(1) / dl(1)
   ro3top  = o3l(1)  / dl(1)

   do k = 2,nrad
      rh2obot = rh2otop
      rco2bot = rco2top
      ro3bot  = ro3top

      rh2otop = max(rl(k) / dl(k), rvmin)
      rco2top = co2l(k) / dl(k)
      ro3top  = o3l(k)  / dl(k)
      vp(k)   = pl(k) * rh2otop  / (ep + rh2otop)
      u(k,1)  = 0.5 * (rh2obot + rh2otop) * (pl(k-1) - pl(k))
      u(k,2)  = 0.5 * (rco2bot + rco2top) * (pl(k-1) - pl(k))
      u(k,3)  = 0.5 * (ro3bot  + ro3top ) * (pl(k-1) - pl(k))
   end do

   u(1,1) = u(2,1)
   u(1,2) = u(2,2)
   u(1,3) = u(2,3)

   return
end subroutine path_lengths
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine print_longwave(mynum,nrad,m1,i,j,time,rlong,rlongup)
   use micphys  , only: ncat
   use harr_coms, only: rl,dl,pl,co2l, o3l,tl,rcl_parm,rpl_parm,rxharr,cxharr,flxus,flxds  &
                       ,flxul,flxdl,embharr
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer     , intent(in) :: mynum,nrad,m1,i,j
   real        , intent(in) :: rlong,rlongup
   real(kind=8), intent(in) :: time
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k,icat
   integer             :: m
   !---------------------------------------------------------------------------------------!
   m = mynum+30

   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   write (unit=m,fmt='(a)')               ' NON-SENSE LONGWAVE RADIATION: '
   write (unit=m,fmt='(3(a,1x,i5),a)')    ' + On node: ',mynum,' (i,j)=(',i,',',j,')'
   write (unit=m,fmt='(a,1x,es12.5)')     ' + RLONG  : ',rlong
   write (unit=m,fmt='(a,1x,es12.5)')     ' + RLONGUP: ',rlongup
   write (unit=m,fmt='(a,1x,es12.5)')     ' + TIME   : ',time
   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   write (unit=m,fmt='(a3,1x,10(a12,1x))') 'LEV','  MIX. RATIO','     DENSITY'             &
                                                ,'    PRESSURE','        CO_2'             &
                                                ,'       OZONE',' TEMPERATURE'             &
                                                ,'       FLXUS','       FLXDS'             &
                                                ,'       FLXUL','       FLXDL'

   do k=1,nrad
      write (unit=m,fmt='(i3,1x,10(es12.3,1x))') k, rl(k), dl(k), pl(k)                    &
                                                  , co2l(k), o3l(k), tl(k)                 &
                                                  , flxus(k), flxds(k), flxul(k), flxdl(k)
   end do
   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   write (unit=m,fmt='(a3,1x,7(a12,1x))') 'LEV','      RCLOUD','       RRAIN'              &
                                               ,'   RPRISTINE','       RSNOW'              &
                                               ,' RAGGREGATES','    RGRAUPEL'              &
                                               ,'       RHAIL'
   do k=1,m1
      write (unit=m,fmt='(i3,1x,7(es12.3,1x))') k, (rxharr(k,icat),icat=1,ncat)
   end do
   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   write (unit=m,fmt='(a3,1x,7(a12,1x))') 'LEV','      CCLOUD','       CRAIN'              &
                                               ,'   CPRISTINE','       CSNOW'              &
                                               ,' CAGGREGATES','    CGRAUPEL'              &
                                               ,'       CHAIL'
   do k=1,m1
      write (unit=m,fmt='(i3,1x,7(es12.3,1x))') k, (cxharr(k,icat),icat=1,ncat)
   end do
   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   write (unit=m,fmt='(a3,1x,7(a12,1x))') 'LEV','      ECLOUD','       ERAIN'              &
                                               ,'   EPRISTINE','       ESNOW'              &
                                               ,' EAGGREGATES','    EGRAUPEL'              &
                                               ,'       EHAIL'
   do k=1,m1
      write (unit=m,fmt='(i3,1x,7(es12.3,1x))') k, (embharr(k,icat),icat=1,ncat)
   end do
   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   write (unit=m,fmt='(a3,1x,2(a12,1x))') 'LEV','    RCL_PARM','    RPL_PARM'
   do k=1,m1
      write (unit=m,fmt='(i3,1x,2(es12.3,1x))') k, rcl_parm(k), rpl_parm(k)
   end do
   write (unit=m,fmt='(80a)')             ('-',k=1,80)
   return 
end subroutine print_longwave
!==========================================================================================!
!==========================================================================================!
