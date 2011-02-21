!==================================== Change Log ==========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! 08/31/2008 - Changing thermodynamic calls, agreeing with the new thermodynamic library.  !
!              Also making Kuo as the deepest cloud in case the user opted by it.          !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!

subroutine kuo_cupar_driver()

   use mem_tend
   use mem_cuparm
   use mem_basic
   use mem_grid
   use node_mod

   implicit none

   ! If I reached here, then I am at the deepest convection
   integer, parameter :: icld=1

   
   select case (if_cuinv)
   case (0) !----- This is the "direct" cumulus parametrisation ---------------------------!

      !----- Call the main subroutine -----------------------------------------------------!
      call conpar( mzp,mxp,myp,ia,iz,ja,jz,ibcon                                           &
                 , basic_g(ngrid)%up                  , basic_g(ngrid)%vp                  &
                 , basic_g(ngrid)%wp                  , basic_g(ngrid)%theta               &
                 , basic_g(ngrid)%pp                  , basic_g(ngrid)%pi0                 &
                 , basic_g(ngrid)%dn0                 , basic_g(ngrid)%rv                  &
                 , cuparm_g(ngrid)%thsrc (:,:,:,icld) , cuparm_g(ngrid)%rtsrc (:,:,:,icld) &
                 , grid_g(ngrid)%rtgt                 , cuparm_g(ngrid)%conprr(  :,:,icld) &
                 , grid_g(ngrid)%flpw                 )

   case (1) !----- This is the cumulus inversion method -----------------------------------!
      !------------------------------------------------------------------------------------!
      !     Check cumulus inversion tendencies and see if they are usable. If so, put in   !
      ! thsrc,rtscr,conprr arrays.                                                         !
      !------------------------------------------------------------------------------------!
      call cu_inv_tend(mzp,mxp,myp,ia,iz,ja,jz                                             &
             ,cuparm_g(ngrid)%thsrc    (:,:,:,icld) ,cuparm_g(ngrid)%thsrcp                &
             ,cuparm_g(ngrid)%thsrcf                ,cuparm_g(ngrid)%rtsrc    (:,:,:,icld) &
             ,cuparm_g(ngrid)%rtsrcp                ,cuparm_g(ngrid)%rtsrcf                &
             ,cuparm_g(ngrid)%conprr   (  :,:,icld) ,cuparm_g(ngrid)%conprrp               &
             ,cuparm_g(ngrid)%conprrf               )
   end select

   return
end subroutine kuo_cupar_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine cu_inv_tend(m1,m2,m3,ia,iz,ja,jz,thsrc,thsrcp,thsrcf,rtsrc,rtsrcp,rtsrcf,conprr &
                      ,conprrp,conprrf)

   use mem_cuparm
   use mem_grid

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)  :: m1,m2,m3,ia,iz,ja,jz
   real, dimension(   m2,m3), intent(in)  :: conprrp,conprrf
   real, dimension(m1,m2,m3), intent(in)  :: thsrcp,thsrcf,rtsrcp,rtsrcf
   real, dimension(   m2,m3), intent(out) :: conprr
   real, dimension(m1,m2,m3), intent(out) :: thsrc,rtsrc
   !----- Local variables -----------------------------------------------------------------!
   integer                                :: k,i,j
   real                                   :: tfact,grwt
   !---------------------------------------------------------------------------------------!

   thsrc(:,:,:) = 0.
   rtsrc(:,:,:) = 0.
   conprr(:,:)  = 0.

   if (time < tcu_beg .or. time > tcu_end ) return

   grwt = wt_cu_grid(ngrid)/tnudcu

   if ( (cutime2-cutime1) <= CU_TIL ) then
      !----- If past and future files are good for interpolation... -----------------------!
      tfact= (time-cutime1)/(cutime2-cutime1) * grwt
      do j=ja,jz
         do i=ia,iz
            do k=2,m1-1
               thsrc(k,i,j) = thsrcp(k,i,j) + tfact * (thsrcf(k,i,j)-thsrcp(k,i,j))
               rtsrc(k,i,j) = rtsrcp(k,i,j) + tfact * (rtsrcf(k,i,j)-rtsrcp(k,i,j))
            end do
            conprr(i,j)     = conprrp(i,j)  + tfact * (conprrf(i,j)-conprrp(i,j))
         end do
      end do

   elseif ( abs(time-cutime1) <= CU_TEL ) then
      !----- Past file close enough... ----------------------------------------------------!
      do j=ja,jz
         do i=ia,iz
            do k=2,m1-1
               thsrc(k,i,j) = thsrcp(k,i,j) * grwt
               rtsrc(k,i,j) = rtsrcp(k,i,j) * grwt
            end do
            conprr(i,j) = conprrp(i,j) * grwt
         end do
      end do

   elseif ( abs(time-cutime2) <= CU_TEL ) then
      !----- Future file close enough... --------------------------------------------------!
      do j=ja,jz
         do i=ia,iz
            do k=2,m1-1
               thsrc(k,i,j) = thsrcf(k,i,j) * grwt 
               rtsrc(k,i,j) = rtsrcf(k,i,j) * grwt
            end do
            conprr(i,j) = conprrf(i,j) * grwt
         end do
      end do

   end if

   return
end subroutine cu_inv_tend
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine conpar(m1,m2,m3,ia,iz,ja,jz,ibcon,up,vp,wp,theta,pp,pi0,dn0,rv,thsrc,rtsrc,rtgt &
                 ,conprr,flpw)

   use conv_coms
   use mem_grid
   use mem_cuparm
   use rconstants

   implicit none

   integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,lpw
   real :: up(m1,m2,m3),vp(m1,m2,m3),wp(m1,m2,m3),theta(m1,m2,m3)  &
            ,pp(m1,m2,m3),pi0(m1,m2,m3),dn0(m1,m2,m3),rv(m1,m2,m3)  &
            ,thsrc(m1,m2,m3),rtsrc(m1,m2,m3),conprr(m2,m3)  &
            ,rtgt(m2,m3)
   real, dimension(m2,m3) :: flpw

   integer :: icpcnt=0,i1,i2,j1,j2,i,j,k,iprtfrq,iqmax,jqmax,kqmax
   real :: dthmax

   !
   !        FLAG TO CONTROL PRINTOUT
   !          ICPRTFL=0 - NO PRINTOUT
   !                  1 - BRIEF INFO ON UP/DOWN DRAFT AND ITERATIONS
   !                  2 - 1 PLUS MODEL TENDENCIES
   !                  3 - 2 PLUS FINAL CONVECTIVE STRUCTURE
   !                  4 - 3 PLUS UP/DOWN DRAFT AND ENVIRONMENT

   icprtfl=0
   iprtfrq=8
   icpltfl=0
   icpcnt=icpcnt+1
   if(mod(icpcnt-iprtfrq+1,iprtfrq).eq.0) then
     icprtfl=1
   endif

   i1 = ia
   i2 = iz
   j1 = ja
   j2 = jz

   !  If variable initialization, on a coarse grid, not in a global simulation,
   !  do not run convective parameterization in the lateral boundary region.

   !  This is commented out for now since it is annoying to compute i1,i2,j1,j2
   !    when in parallel. This means we are now computing conv tendencies in the
   !    coarse grid nudging boundary regions.
   !           We will see if this causes problems ...

   !      IF (INITIAL .EQ. 2 .AND. nxtnest(ngrid) .EQ. 0 .and.
   !     +    nhemgrd2 .le. 1) THEN

   !         if (iand(ibcon,1) .gt. 0) i1 = 1 + nupts
   !         if (iand(ibcon,2) .gt. 0) i2 = nxp - nupts
   !         if (iand(ibcon,4) .gt. 0) j1 = 1 + nupts * jdim
   !         if (iand(ibcon,8) .gt. 0) j2 = max (1,nyp - nupts)

   !      ENDIF

   dthmax=0.

   do j = j1,j2
      do i = i1,i2

         do k = 1,m1
            ucon(k)      = up(k,i,j)
            vcon(k)      = vp(k,i,j)
            wcon(k)      = wp(k,i,j)
            thtcon(k)    = theta(k,i,j)
            picon(k)     = (pp(k,i,j)+pi0(k,i,j))
            tmpcon(k)    = thtcon(k)*picon(k)/cp
            dncon(k)     = dn0(k,i,j)
            prcon(k)     = (picon(k)/cp)**cpor*p00
            rvcon(k)     = rv(k,i,j)
            zcon(k)      = zt(k) *rtgt(i,j)
            zzcon(k)     = zm(k) *rtgt(i,j)
            thsrc(k,i,j) = 0.
            rtsrc(k,i,j) = 0.
         end do
         conprr(i,j)=0.
         wconmin=wcldbs
         contim=confrq

         lpw = nint(flpw(i,j))

         call cu_environ(lpw,m1-1)

         if(igo /= 0)  call kuocp()
         !----- igo may change in kuocp, need two ifs here --------------------------------!
         if(igo /= 0) then

           call cp2mod(lpw,m1)

           do k=lpw,m1-1
              thsrc(k,i,j)=ftcon(k)
              rtsrc(k,i,j)=frcon(k)
           end do

           conprr(i,j)=cprecip

           do k= lpw,m1-1
              if(thsrc(k,i,j).gt.dthmax) then
                dthmax=thsrc(k,i,j)
                iqmax=i
                jqmax=j
                kqmax=k
              end if
           end do
         end if
      end do
   end do

   return
end subroutine conpar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine cu_environ(k1,k2)

   use conv_coms
   use rconstants
   use therm_lib, only: thetaeiv,thetaeivs2temp,virtt,lcl_il

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)  :: k1,k2
   !----- Local variables -----------------------------------------------------------------!
   real, dimension(nkp) :: hz
   real                 :: wcpmax,themax,thlll,tlll,plll,rlll,zlll,dzlll,dzdd,abe
   real                 :: thdu,tdu,rdsu,znz
   integer              :: k,nkmid,nk
   logical              :: bypass
   !---------------------------------------------------------------------------------------!

   !----- Basic constants -----------------------------------------------------------------!
   dzlow=200.
   dzhigh=500.
   zmid=3000.
   cdzmin=3000.

   !----- Compute moist static energy profile ---------------------------------------------!
   do k=k1,k2
     hz(k)=cp*tmpcon(k)+grav*zcon(k)+alvl*rvcon(k)
   enddo

   !---------------------------------------------------------------------------------------!
   !     Check for conditional instability and any upward motion greater than WCONMIN      !
   ! under ZMID.                                                                           !
   !---------------------------------------------------------------------------------------!
   igo=0
   coinsloop1: do k=k1,k2
      if(hz(k) > hz(k+1))then
        igo=1
        exit coinsloop1
      end if
   end do coinsloop1
   if(igo == 0) return

   igo=0
   wcpmax=-1.e10
   coinsloop2: do k=k1,k2
      if(zcon(k) > zmid) exit coinsloop2
      wcpmax=max(wcpmax,wcon(k))
   end do coinsloop2

   if(wcpmax > 0.0 .and. wcpmax > wconmin) igo=1
   if(igo == 0) return

   !----- Interpolate model sounding (environment) to higher resolution grid. -------------!
   nkmid = zmid/dzlow+1
   zc(1) = 0.

   do k=2,nkmid
     zc(k) = zc(k-1) + dzlow
   end do

   do k=nkmid+1,nkp
     zc(k) = zc(k-1) + dzhigh
   end do

   ze(1) = 0.
   do k=2,nkp
     ze(k)= (zc(k)+zc(k-1)) * .5
   end do

   !----- Find model top on convective grid -----------------------------------------------!
   znz=zcon(k2)
   bypass=.false.
   toploop: do k=nkp,1,-1
     if (ze(k) < znz) then
        bypass=.true.
        exit toploop
     end if 
   end do toploop
   if (.not. bypass) call abort_run('Envir stop 12','cu_environ','kuo_cupar_driver.f90')
   kmt=k

   !----- Do actual interpolation ---------------------------------------------------------!
   nk=k2-k1+1
   call htint(nk,ucon,zcon(k1),kmt,upe,ze)
   call htint(nk,vcon,zcon(k1),kmt,vpe,ze)
   call htint(nk,wcon,zzcon(k1),kmt,wpe,ze)
   call htint(nk,thtcon,zcon(k1),kmt,the,ze)
   call htint(nk,rvcon,zcon(k1),kmt,rve,ze)
   do k=1,kmt
      rve(k)=max(rve(k),1e-8)
   end do

   !----- Compute theta v, theta e, and get pressure profile ------------------------------!
   pke(1)=picon(1)
   do k=1,kmt
     thve(k)=virtt(the(k),rve(k))
   end do

   do k=2,kmt
     pke(k)=pke(k-1)-grav*2.*(ze(k)-ze(k-1))/(thve(k)+thve(k-1))
   end do

   do k=1,kmt
     te(k)=the(k)*pke(k)/cp
     pe(k)=(pke(k)/cp)**cpor*p00
     rhoe(k)=pe(k)/(rdry*virtt(te(k),rve(k)))
   end do
   do k=1,kmt
     thee(k)=thetaeiv(the(k),pe(k),te(k),rve(k),rve(k),3,.false.)
   end do


   !----- Find the main source level of the updraft. --------------------------------------!
   !----- First test - any inversion below 1.2 km -----------------------------------------!
   bypass = .false.
   uploop: do k=3,nkmid
     if(te(k) > te(k-1) .and. te(k) > te(k+1) .and. ze(k) <= 1200.) then
        kcon=k
        bypass = .true.
        exit uploop
     end if
   end do uploop

   !----- If there isn't an inversion, use the level of highest theta_E. ------------------!
   if (.not. bypass) then
     themax=0.
     do k=2,nkmid
       if(thee(k).gt.themax)then
         themax=thee(k)
         kcon=k
       endif
     enddo
   end if

   !----- Find the LCL of a layer average around the source level -------------------------!
   tlll  = (te(kcon)+te(kcon+1)+te(kcon-1))/3.
   plll  = pe(kcon)
   rlll  = (rve(kcon)+rve(kcon+1)+rve(kcon-1))/3.
   zlll  = ze(kcon)
   thlll = tlll * (p00/plll)**rocp
   call lcl_il(thlll,plll,tlll,rlll,rlll,tlcl,plcl,dzlcl,1,.false.)
   if (dzlcl == 0.) then
      tlcl = tlll
      plcl = plll
   end if


   !----- Find the closest level on the convective grid to the LCL ------------------------!
   dzlll=1e20
   do k=1,kmt
      dzdd=abs(ze(k)-(zlll+dzlcl))
      if(dzdd < dzlll)then
        dzlll=dzdd
        klcl=k
      end if
   end do


   !----- If there is not upward motion at the lcl, no convection (must be > wconmin) -----!
   if(wpe(klcl) < 0.0 .or. wpe(klcl) < wconmin) then
     igo=0
     return
   endif

   !---------------------------------------------------------------------------------------!
   !     Locate equilibrium temperature level of an unentrained parcel. compute initial    !
   ! ABE. If ABE is less than 0, no convection.                                            !
   !---------------------------------------------------------------------------------------!
   theu(klcl) = the(kcon)*exp(alvl*rve(kcon)/(cp*max(tlcl,ttripoli)))

   bypass = .false.
   eqloop: do k=klcl,kmt
     if(theu(klcl) <= thve(k)) then
        bypass = .true.
        exit eqloop
     end if 
   end do eqloop
   if (.not. bypass) then
      print*,'convection above model top:',klcl,theu(klcl),thve(k)
      ketl=kmt-2
   end if
   ketl=k
   if(ze(ketl)-ze(klcl)  < cdzmin)then
     igo=0
     return
   end if

   abe=0.
   do k=klcl,ketl
     call thetaeivs2temp(theu(klcl),pe(k),thdu,tdu,rdsu,.false.)
     abe=abe+(virtt(thdu,rdsu)-thve(k))/thve(k)*(zc(k)-zc(k-1))
   end do
   if(abe <= 0.)then
     igo=0
     return
   end if

   return
end subroutine cu_environ
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine kuocp()

   use conv_coms
   use rconstants
   use therm_lib, only : thetaeivs2temp

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer :: k,idownd,klfs,kdiv,kdet,kover,kcoolh,kheat
   real    :: supplyw,anegl,apos,anegh,dddt,dzdiv,wtlfs,wtlcl,wtdiv,wtgnd,bkuo,zdetr,dzdet
   real    :: vhint,vmint,vdint,avgmin,avtdiff,overmax,factr,heatmx,coolhi,c1
   logical :: bypass
   !---------------------------------------------------------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   ! Downdraft flag - 0 - no downdrafts                                                    !
   !                  1 - simple downdraft model                                           !
   !---------------------------------------------------------------------------------------!
   idownd=1

   do k=1,nkp
      ftcon(k)=0.
      frcon(k)=0.
   end do

   !---------------------------------------------------------------------------------------!
   !    Compute vertical moisture convergence into the cloud layer. Vertical flux out      !
   ! cloud top is assumed small.                                                           !
   !---------------------------------------------------------------------------------------!
   supplyw=rhoe(klcl)*rve(klcl)*(wpe(klcl)+wpe(klcl-1))*.5
   supply=supplyw

   if(supply <= 0.) then
     igo=0
     return
   end if

   !---------------------------------------------------------------------------------------!
   !    This is the cloud model.  Updraft is constant THETA e and saturated with respect   !
   ! to water.  There is no ice.  Cloud top is one level above ETL. THETA e of the updraft !
   !---------------------------------------------------------------------------------------!
   theu(klcl)=the(kcon)*exp(alvl*rve(kcon)/(cp*tlcl))

   !----- Equilibrium Temperature Level of the source level air. --------------------------!
   igo    = 0
   bypass = .false.
   srcloop: do k=klcl,kmt
     call thetaeivs2temp(theu(klcl),pe(k),thu(k),tu(k),rsu(k),.false.)
     if(thu(k) > the(k) .and. igo == 0) then
       igo=1
       klfc=k
     end if
     if(thu(k).le.the(k).and.igo.eq.1) then
        bypass=.true.
        exit srcloop
     end if
   end do srcloop
   if(igo == 0) return
   if (.not. bypass) then
      print *,' Convection beyond model top - THup, THenv ',THU(KMT),THE(KMT)
      k=kmt-1
   end if

   ketl = min(k,kmt)
   kct  = min(ketl+1,kmt)
   call thetaeivs2temp(theu(klcl),pe(kct),thu(kct),tu(kct),rsu(kct),.false.)
   do k=1,klcl-1
     thu(k) = the(k)
   end do

   !---------------------------------------------------------------------------------------!
   !    If the cloud is not at least CDZMIN deep or cloud top is under 500 mb, no          !
   ! convection.                                                                           !
   !---------------------------------------------------------------------------------------!
   if (ze(ketl)-ze(klfc) < cdzmin .or. pe(kct) > 50000.) then
      igo=0
      return
   end if

   !---------------------------------------------------------------------------------------!
   !    Require the positive area be 50% greater than the negative area below the LFC and  !
   ! 5% greater in total.                                                                  !
   !---------------------------------------------------------------------------------------!
   anegl=0.
   do k=klcl,klfc-1
      anegl = anegl+(thu(k)-the(k))*(zc(k)-zc(k-1))
   end do

   apos=0.
   do k=klfc,ketl-1
      apos = apos+(thu(k)-the(k))*(zc(k)-zc(k-1))
   end do

   anegh=0.
   do k=ketl,kct
      anegh = anegh+(thu(k)-the(k))*(zc(k)-zc(k-1))
   end do

   if(apos < abs(anegl)*1.5 .or. apos < abs(anegl+anegh)*1.05) then
      igo=0
      return
   end if

   if (idownd == 1) then
      !------------------------------------------------------------------------------------!
      !    The downdraft model - starts at THETA e minimum (LFS). Downdraft is 2 degrees   !
      ! colder than environment at cloud base increasing to 5 degrees colder at the        !
      ! ground. Find LFS as THETA e minimum                                                !
      !------------------------------------------------------------------------------------!
      bypass = .false.
      downloop: do k=kct,2,-1
        if(thee(k) < thee(k+1) .and. thee(k) < thee(k-1)) then
           bypass=.true.
           exit downloop
        end if
      end do downloop
      if (.not. bypass) k=2
      klfs=k
      if(klfs <= klcl) klfs=klcl+1
      thd(klfs) = the(klfs)

      !------------------------------------------------------------------------------------!
      !     Limit dd deficit at the ground to the maximum of positive temperature          !
      ! difference of updraft if less than 2.5 degrees.                                    !
      !------------------------------------------------------------------------------------!
      dddt=0.
      do k=klcl,kct
        dddt=max(dddt,thu(k)-the(k))
      end do
      if(dddt > 2.5) dddt=5.

      thd(2)    = the(2)   - dddt
      thd(klcl) = the(klcl)- dddt*.2
      do k=klcl,klfs
        thd(k) = thd(klcl) + (thd(klfs)-thd(klcl)) / (ze(klfs)-ze(klcl)) * (ze(k)-ze(klcl))
      end do
      do k=3,klcl-1
        thd(k) = thd(2) + (thd(klcl)-thd(2)) / (ze(klcl)-ze(2)) * (ze(k)-ze(2))
      enddo

      !------------------------------------------------------------------------------------!
      !    Now we need to weight the downdraft relative to the updraft. Assume that the dd !
      ! weight is zero at the LFS, 1/2 of updraft at cloud base, and equal to the updraft  !
      ! at cloud base at the ground.                                                       !
      !------------------------------------------------------------------------------------!
      dzdiv = 1.e20
      do k=1,kmt
         if(abs(ze(k)-800.) < dzdiv) then
            kdiv  = k
            dzdiv = abs(ze(k)-800.)
         end if
      end do
      kdiv=max(min(klcl,kdiv),2)
      if(kdiv == klcl) kdiv=klcl-1

      do k=1,nkp
        wtd(k)=0.
      end do
      wtlfs=0.
      wtlcl=.1
      wtdiv=.2
      wtgnd=1.
      do k=klcl+1,klfs
        wtd(k) = wtlcl + (wtlfs-wtlcl) / (ze(klfs)-ze(klcl)) * (ze(k)-ze(klcl))
      end do
      do k=kdiv,klcl
        wtd(k) = wtdiv + (wtlcl-wtdiv) / (ze(klcl)-ze(kdiv)) * (ze(k)-ze(kdiv))
      end do
      do k=2,kdiv-1
        wtd(k) = wtgnd + (wtdiv-wtgnd) / (ze(kdiv)-ze(2)) * (ze(k)-ze(2))
      end do

   else
      do k=1,nkp
        wtd(k)=0.
      end do
      do k=2,klcl-1
        thu(k)=the(k)
      end do
   end if


   !----- Compute infamous b parameter.  Use Fritsch/Chappell's precipitation efficiency. -!
   envshr = sqrt((upe(kct)-upe(klfc))**2+(vpe(kct)-vpe(klfc))**2)/(ze(kct)-ze(klfc))*1e3
   if(envshr > 1.35) then
     preff = 1.591-.639*envshr+.0953*envshr**2-.00496*envshr**3
   else
     preff = .9
   end if
   bkuo=1.-preff

   !----- Vertical profiles of convective heating and moistening. -------------------------!
   do k=2,kmt
     vheat(k)=0.
     vmois(k)=0.
     vmdry(k)=0.
   end do

   !----- Find the weighted THETA to use for the convection. ------------------------------!
   do k=2,kct
      thcon(k) = wtd(k)*thd(k)+(1.-wtd(k))*thu(k)
   end do

   !----- Heating profile is difference between convective THETAs and environment. --------!
   do k=2,kct
     vheat(k) = thcon(k)-the(k)
   end do

   !---------------------------------------------------------------------------------------!
   !    Moisture profile is difference between vapor's of updraft and environment in the   !
   ! cloud layer.  Below cloud base, air is dried by SUPPLY.  Downdrafts are assumed to    !
   ! have no effect on this.                                                               !
   !---------------------------------------------------------------------------------------!
   zdetr = .66667*ze(kct)
   dzdet = 1000000.
   do k=klcl,kct
      if(abs(ze(k)-zdetr) < dzdet) then
         dzdet = abs(ze(k)-zdetr)
         kdet  = k
      end if
   end do

   do k=kdet,kct
     vmois(k)=1.
   end do

   do k=2,klcl-1
     vmdry(k)=rve(k)
   end do

   vhint=0.
   vmint=0.
   vdint=0.
   do k=2,kmt
      vhint=vhint+vheat(k)*(zc(k)-zc(k-1))
      vmint=vmint+vmois(k)*(zc(k)-zc(k-1))
      vdint=vdint+vmdry(k)*(zc(k)-zc(k-1))
   end do

   !---------------------------------------------------------------------------------------!
   !     If VHINT is less than 0, there is more negative area than positive area.  No      !
   ! convection allowed.                                                                   !
   !---------------------------------------------------------------------------------------!
   if(vhint <= 0.) then
     igo=0
     return
   end if

   !---------------------------------------------------------------------------------------!
   !     Also require that there is a minimum average temperature difference between the   !
   ! updraft and environment from the LFC to the ETL.  This eliminates the cases where     !
   ! VHINT is very small and the heating and cooling rates get astronomically large.       !
   !---------------------------------------------------------------------------------------!
   avgmin  = .10
   avtdiff = 0.
   do k = klfc,ketl-1
      avtdiff = avtdiff+(thcon(k)-the(k))
   end do
   avtdiff=avtdiff/max(1,ketl-klfc)
   if(avtdiff < avgmin) then
      igo=0
      return
   end if

   !----- Heating and moistening rates. ---------------------------------------------------!

   rateloop: do
      do k=2,kmt
        ftcon(k) = alvl*preff*supply*vheat(k) /(pke(k)*rhoe(k)*vhint)
      end do
      do k=klcl,kct
        frcon(k)=bkuo*supply*vmois(k)/(rhoe(k)*vmint)
      end do
      do k=2,klcl-1
        frcon(k)=-supply*vmdry(k)/(rhoe(k)*vdint)
      end do

      do k=klfc,ketl-1
        qvct1(k)=the(k)+contim*ftcon(k)
      end do
      overmax=0.
      do k=klfc,ketl-1
        if(qvct1(k)-thu(k) > overmax)then
          overmax=(qvct1(k)-thu(k))/(ftcon(k)*contim)
          kover=k
        endif
      end do

      if(overmax > 0.) then
        factr=1.-overmax
        supply=factr*supply
      else
         exit rateloop
      end if
   end do rateloop
   cprecip = preff*supply

   if (icprtfl > 0) then
      coolhi=100000.
      heatmx=-10000.
      kcoolh=0
      kheat=0
      do k=ketl,kct
        if(ftcon(k) < coolhi) then
          coolhi=ftcon(k)
          kcoolh=k
        end if
      end do
      do k=klcl,kct
        if(ftcon(k) > heatmx) then
          heatmx=ftcon(k)
          kheat=k
        end if
      end do

      c1=86400.
   end if

   return
end subroutine kuocp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine cp2mod(k1,k2)

   use conv_coms
   use mem_scratch
   use rconstants

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: k1,k2
   !----- Local variables -----------------------------------------------------------------!
   integer             :: k
   real                :: tftc,tftm,tfrc,tfrm,ftres,frres
   !----- Functions -----------------------------------------------------------------------!
   real, external :: ssum
   !---------------------------------------------------------------------------------------!


   !----- Compute integrated heating and moistening tendencies. ---------------------------!
   do k=2,kmt
     qvct1(k) = rhoe(k)*ftcon(k)*pke(k)
     qvct2(k) = rhoe(k)*alvl*frcon(k)
     qvct3(k) = (zc(k)-zc(k-1))*qvct1(k)
     qvct4(k) = (zc(k)-zc(k-1))*qvct2(k)
   end do
   tftc=ssum(kmt-1,qvct3(2),1)
   tfrc=ssum(kmt-1,qvct4(2),1)

   !----- Transfer tendencies to model grid. ----------------------------------------------!
   call vertmap2(qvct1,zc,kmt,vctr5,zzcon,k2)
   call vertmap2(qvct2,zc,kmt,vctr6,zzcon,k2)
   do k=k1,k2
     vctr5(k)=vctr5(k)*(zzcon(k)-zzcon(k-1))
     vctr6(k)=vctr6(k)*(zzcon(k)-zzcon(k-1))
   end do


   !---------------------------------------------------------------------------------------!
   !    Make sure the transfer from the convective grid to the model grid happened         !
   ! correctly.                                                                            !
   !---------------------------------------------------------------------------------------!
   tftm=ssum(k2-k1+1,vctr5(k1),1)
   tfrm=ssum(k2-k1+1,vctr6(k1),1)
   ftres=tftm-tftc
   frres=tfrm-tfrc
   if(abs(ftres) > .01*abs(tftc)) then
     print*,' energy error in grid tranfser in convective param.'
     print*,' tftm,tftc ',tftm,tftc
   end if


   !----- Change energy tendencies to temperature and mixing ratio tendencies. ------------!
   do k=k1,k2
     ftcon(k) = vctr5(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*picon(k))
     frcon(k) = vctr6(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*alvl)
   end do

return
end subroutine cp2mod
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine vertmap2(datin,zin,n3in,datout,zout,n3out)

   implicit none

   integer :: n3in,n3out
   real, dimension(n3in) :: datin,zin
   real, dimension(n3out) :: datout,zout(n3out)

   real, allocatable :: qvct(:),vctr(:)
   integer :: k,l
   real :: dzlft
   
   !---------------------------------------------------------------------------------------!
   !    This routine assumes that output vertical structure will be lower than input       !
   ! vertical structure.                                                                   !
   !---------------------------------------------------------------------------------------!


   !----- Transfer quantity from input grid levels to output grid. ------------------------!
   allocate(qvct(n3in),vctr(n3out))

   do k=1,n3out
      vctr(k)=0.
   end do

   dzlft=0.
   l=2
   do k=2,n3out
      if(dzlft.ne.0.) then
         if(zin(l) .gt. zout(k)) then
            vctr(k)=vctr(k)+datin(l)*(zout(k)-zout(k-1))
            dzlft=zin(l)-zout(k)
            go to 61
         else
            vctr(k)=vctr(k)+datin(l)*dzlft
            l=l+1
            if(l > n3in) exit
            dzlft=0.
         endif
      endif
   60 continue
      if(zin(l) <= zout(k)) then
         vctr(k)=vctr(k)+datin(l)*(zin(l)-zin(l-1))
         l=l+1
         if(l > n3in) exit
         dzlft=0.
         if(zin(l-1) == zout(k)) go to 61
         go to 60
      else
         vctr(k)=vctr(k)+datin(l)*(zout(k)-zin(l-1))
         dzlft=zin(l)-zout(k)
      endif
   61 continue
   end do


   !----- Change energy tendencies to temperature and mixing ratio tendencies. ------------!
   do k=2,n3out
     datout(k)=vctr(k)/(zout(k)-zout(k-1))
   end do

   deallocate(qvct,vctr)

   return
end subroutine vertmap2
!==========================================================================================!
!==========================================================================================!

