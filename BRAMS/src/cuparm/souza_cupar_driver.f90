!==========================================================================================!
!                                                                                          !
! Programmed by  Enio Pereira de Souza                                                     !
! Adapted    by  Alvaro Luiz Fazenda    (for V.5.04)                                       !
!                                                                                          !
!==========================================================================================!
subroutine souza_cupar_driver()
   ! USE Modules for 5.0
   use mem_basic
   use mem_micro
   use mem_grid
   use mem_turb
   use mem_tend
   use node_mod, only : MXP,   &   ! INTENT(IN)
        MYP,                   &   ! INTENT(IN)
        MZP,                   &   ! INTENT(IN)
        IA,                    &   ! INTENT(IN)
        IZ,                    &   ! INTENT(IN)
        JA,                    &   ! INTENT(IN)
        JZ,                    &   ! INTENT(IN)
        I0,                    &   ! INTENT(IN)  ! Rever função
        J0                         ! INTENT(IN)  ! Rever função

   use mem_cuparm, only : cuparm_g,nclouds       ! INTENT(IN)

   implicit none

   integer :: I, J
   integer :: icld
   icld = nclouds ! Just to make it similar to other methods

   call shcupar(mzp,mxp,myp,ia,iz,ja,jz,i0,j0,                   &
        basic_g(ngrid)%wp, basic_g(ngrid)%theta,   &
        basic_g(ngrid)%pp, basic_g(ngrid)%pi0,     &
        basic_g(ngrid)%dn0, basic_g(ngrid)%rv,     &
        cuparm_g(ngrid)%thsrc(:,:,:,icld),                       &
        cuparm_g(ngrid)%rtsrc(:,:,:,icld),                       &
        cuparm_g(ngrid)%upmf(:,:,icld), grid_g(ngrid)%rtgt, &
        turb_g(ngrid)%sflux_t,                              &
        turb_g(ngrid)%sflux_r,                              &
        turb_g(ngrid)%vkh,                                &
        micro_g(ngrid)%rcp)
   return
end subroutine souza_cupar_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine shcupar(m1,m2,m3,ia,iz,ja,jz,i0,j0,wp,theta,pp,pi0,dn0,rv,thsrcsh,rtsrcsh,shmf  &
                  ,rtgt,tfz,qfz,khv,rcloud)

   use conv_coms, only :         &
        wcon,                    &   ! intent(out)
        thtcon,                  &   ! intent(out)
        picon,                   &   ! intent(out)
        dncon,                   &   ! intent(out)
        zcon,                    &   ! intent(out)
        zzcon,                   &   ! intent(out)
        icprtfl,                 &   ! intent(out)   ! maybe local variable
        icpltfl,                 &   ! intent(out)   ! maybe local variable
        igo,                     &   ! intent(out)   ! maybe local variable
        klcl                         ! intent(in)
   use shcu_vars_const, only :   &
        entf,                    &   ! intent(out)
        alhf,                    &   ! intent(out)
        qvcon,                   &   ! intent(out)
        akvd,                    &   ! intent(out)
        cl_con,                  &   ! intent(out)
        dtdt,                    &   ! intent(in)
        wc,                      &   ! intent(in)
        drdt                         ! intent(in)

   use mem_grid, only :          &
        time,                    &   ! intent(in)
        ZT,                      &   ! intent(in)
        ZM                           ! intent(in/out)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1,m2,m3,ia,iz,ja,jz,i0,j0
   real, dimension(m1,m2,m3), intent(in)    :: wp,theta,pp,pi0,dn0,rv,khv,rcloud
   real, dimension(   m2,m3), intent(in)    :: rtgt,tfz,qfz
   real, dimension(m1,m2,m3), intent(inout) :: thsrcsh,rtsrcsh
   real, dimension   (m2,m3), intent(inout) :: shmf
   !----- Local variables -----------------------------------------------------------------!
   integer                                :: icpcnt = 0
   integer                                :: iprtfrq, i, j, k
   !---------------------------------------------------------------------------------------!

   icprtfl=0
   iprtfrq=8
   icpltfl=0
   icpcnt=icpcnt+1
   if(mod(icpcnt-iprtfrq+1,iprtfrq).eq.0) then
      icprtfl=1
   endif

   jloop: do j=ja,jz
      iloop: do i=ia,iz
         !----- Def. of the enthalpy flux ENTF (K*m/s) and latent energy flux (m/s)*(g/g) -!
         entf=min(.50,tfz(i,j))
         alhf=min(.00028,qfz(i,j))

         do k=1,m1
            wcon(k)        = wp(k,i,j)
            thtcon(k)      = theta(k,i,j)
            picon(k)       = (pp(k,i,j)+pi0(k,i,j))
            dncon(k)       = dn0(k,i,j) 
            !----- At this point, we change rv by qv as required by the basic equations. --!
            qvcon(k)       = rv(k,i,j)/(1+rv(k,i,j))
            zcon(k)        = zt(k) *rtgt(i,j)
            zzcon(k)       = zm(k) *rtgt(i,j)
            akvd(k)        = khv(k,i,j)
            cl_con(k)      = rcloud(k,i,j)
            thsrcsh(k,i,j) = 0.
            rtsrcsh(k,i,j) = 0.
         end do

         if(entf <= 0.0) cycle iloop
         igo=1

         call shcu_env(m1-1) 

         if(igo /= 0) call cl_top
         if(igo /= 0) call w_shallow(i,j,time)
         if(igo /= 0) call sh_rates
         if(igo /= 0) then
            call sh2mod(m1)
            do k=2,m1-1
               thsrcsh(k,i,j)=dtdt(k)
               rtsrcsh(k,i,j)=drdt(k) 
            end do
            shmf(i,j) = dncon(klcl)*wc(klcl)
         end if
      end do iloop
   end do jloop
   return
end subroutine SHCUPAR
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine shcu_env(nz)  

   use conv_coms, only : nkp,    &   ! intent(in)  ! parameter
        zc,                      &   ! intent(out)
        ze,                      &   ! intent(out)
        pke,                     &   ! intent(out)
        the,                     &   ! intent(in)
        thve,                    &   ! intent(out)
        te,                      &   ! intent(out)
        pe,                      &   ! intent(out)
        rhoe,                    &   ! intent(out)
        dzlow,                   &   ! intent(out) ! talvez var.local
        dzhigh,                  &   ! intent(out) ! talvez var.local
        zmid,                    &   ! intent(out) ! talvez var.local
        cdzmin,                  &   ! intent(out) ! rever função.
        zcon,                    &   ! intent(in)
        kmt,                     &   ! intent(out)
        wcon,                    &   ! intent(in/out)
        zzcon,                   &   ! intent(in)
        wpe,                     &   ! intent(in/out)
        thtcon,                  &   ! intent(in)
        picon,                   &   ! intent(in)
        kcon,                    &   ! intent(out)  ! talvez var.local
        tlcl,                    &   ! intent(in/out)
        plcl,                    &   ! intent(in/out)
        dzlcl,                   &   ! intent(in/out)
        klcl,                    &   ! intent(out)
        igo                          ! intent(out)
   use shcu_vars_const, only : qve,    &   ! intent(in/out)
        dse,                     &   ! intent(out)
        uhe,                     &   ! intent(out)
        evaps,                   &   ! intent(out) ! maybe local var.?
        qvse,                    &   ! intent(out) ! maybe local var.?
        uhes,                    &   ! intent(out)
        rhe,                     &   ! intent(out) ! **
        gamma,                   &   ! intent(out)
        uhc,                     &   ! intent(out)
        delz,                    &   ! intent(out)
        dldzby2,                 &   ! intent(out) ! maybe local var.?
        dsc,                     &   ! intent(out)
        dsc0,                    &   ! intent(out)
        qvc,                     &   ! intent(out)
        wlc,                     &   ! intent(out)
        qvcon,                   &   ! intent(in/out)
        akvd,                    &   ! intent(in/out)
        cl_con,                  &   ! intent(in/out)
        cl_pe,                   &   ! intent(in/out)
        g,                       &   ! intent(in)  ! parameter
        cp,                      &   ! intent(in)  ! parameter
        cpr,                     &   ! intent(in)  ! parameter
        p00,                     &   ! intent(in)  ! parameter
        r,                       &   ! intent(in)  ! parameter
        kzi,                     &   ! intent(out) ! maybe local var.?
        alvl,                    &   ! intent(in)  ! parameter
        akvde                        ! intent(in/out)
   use therm_lib, only : lcl_il
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer             , intent(in) :: NZ
   !----- Local variables -----------------------------------------------------------------!
   real, dimension(nkp)             :: dlamb
   !----- Basic constants -----------------------------------------------------------------!
   real                             :: const1,const2,es00,epslon,ummeps,ta0,const3
   real                             :: c0,dlamb0,zref,znz,thlll,tlll,plll,zlll,rlll,dzlll
   real                             :: dzdd,dthv
   integer                          :: nkmid, k
   logical                          :: bypass
   !---------------------------------------------------------------------------------------!

   !----- All variables above have been defined in USE shcu_vars_const --------------------!
   dzlow  = 200.
   dzhigh = 500.
   zmid   = 4000.
   const1 = 17223003.15
   const2 = 29.65
   es00   = 611.2
   epslon = .622
   ummeps = .378
   ta0    = 273.15
   const3 = 17.67
   cdzmin = 3000.
   c0     = 0.0
   dlamb0 = 1.e-06
   zref   = 700.

   !----- Interpolate model sounding (environment) to higher resolution grid --------------!
   nkmid=zmid/dzlow+1
   zc(1)=0.
   do k=2,nkmid
      zc(k)=zc(k-1)+dzlow
   end do
   do k=nkmid+1,nkp
      zc(k)=zc(k-1)+dzhigh
   end do

   ze(1)=0.
   do k=2,nkp
      ze(k)=(zc(k)+zc(k-1))*.5
   enddo

   !----- Find model top on convective grid -----------------------------------------------!
   bypass=.false.
   znz=zcon(nz)
   do k=nkp,1,-1
      if(ze(k).lt.znz) bypass=.true.
   enddo
   if (.not.bypass) call abort_run ('envir stop 12','shcu_env','souza_cupar_driver.f90')
   kmt=k

   !----- Do actual interpolation ---------------------------------------------------------!
   call htint(nz,wcon,zzcon,kmt,wpe,ze)
   call htint(nz,thtcon,zcon,kmt,the,ze)
   call htint(nz,qvcon,zcon,kmt,qve,ze)
   call htint(nz,akvd,zcon,kmt,akvde,ze)
   call htint(nz,cl_con,zcon,kmt,cl_pe,ze)

   do k=1,kmt
      qve(k)=max(qve(k),1e-8)
   end do

   !----- Compute theta v, theta e, and get pressure profile ------------------------------!
   pke(1)=picon(1)
   do k=1,kmt
      thve(k)=the(k)*(1.+.61*qve(k))
   end do

   do k=2,kmt
      pke(k)=pke(k-1)-g*2.*(ze(k)-ze(k-1))/(thve(k)+thve(k-1))
   end do
   do k=1,kmt
      te(k)=the(k)*pke(k)/cp
      pe(k)=(pke(k)/cp)**cpr*p00
      rhoe(k)=pe(k)/(r*te(k)*(1.+.61*qve(k)))
   end do

   !---------------------------------------------------------------------------------------!
   ! EPS - Find the main source level of the updraft. We will assume that parcels start    ! 
   !       from the first layer defined by k=2.                                            !
   !---------------------------------------------------------------------------------------!
   kcon=2
   !----- Find the lcl of a layer average around the source level -------------------------!
   tlll  = (te(kcon)+te(kcon-1))/2.      
   plll  = pe(kcon)
   rlll  = (qve(kcon)+qve(kcon-1))/2.      
   zlll  = ze(kcon)
   thlll = tlll*(p00/plll)**(r/cp)
   call lcl_il(thlll,plll,tlll,rlll,rlll,tlcl,plcl,dzlcl,2,.false.)
   if (dzlcl == 0.) then
      tlcl = tlll
      plcl = plll
   end if

   !----- Find the closest level on the convective grid to the LCL ------------------------!
   dzlll=1e20
   do k=1,kmt
      dzdd=abs(ze(k)-(zlll+dzlcl))
      if(dzdd.lt.dzlll)then
         dzlll=dzdd
         klcl=k
      endif
   end do

   !---------------------------------------------------------------------------------------!
   !    Determination of Zi, the PBL top as the level where the virtual potential tempera- !
   ! ture THVE increases by more than DTHV as compared to the previous level.              !
   !---------------------------------------------------------------------------------------!
   dthv=0.5
   ziloop: do k=3,kmt
      if((thve(k)-(thve(k-1))).gt.dthv) then
         kzi=k
         exit ziloop
      end if
   end do ziloop
   if(kzi < klcl) then
      igo=0
      return
   endif

   !----- EPS - Determination of environment variables ------------------------------------!
   do k=1,kmt
      dse(k)   = cp*te(k)+g*ze(k)
      uhe(k)   = dse(k)+alvl*qve(k)
      evaps(k) = es00*exp(const3*(te(k)-ta0)/(te(k)-const2))
      qvse(k)  = epslon*evaps(k)/(pe(k)-ummeps*evaps(k))
      uhes(k)  = dse(k)+alvl*qvse(k)
      rhe(k)   = qve(k)/qvse(k)
      gamma(k) = const1*pe(k)*qvse(k)**2/evaps(k)
      gamma(k) = gamma(k)/((te(k)-const2)*(te(k)-const2))
   end do

   !---------------------------------------------------------------------------------------!
   !    Calculating the cloud moist static energy profile. We assume that entrainment      !
   ! occurs only above the LCL.                                                            !
   !---------------------------------------------------------------------------------------!
   uhc(1)=uhe(1)
   uhc(2)=uhe(2)
   delz(2)=ze(2)-ze(1)

   !---------------------------------------------------------------------------------------!
   !    The parcel rises without mixing till its LCL and starts being entrained from the   !
   ! LCL up to the cloud top.                                                              !
   !---------------------------------------------------------------------------------------!
   if (klcl >= 3) then
      do k=3,klcl
         uhc(k)=uhe(2)
         delz(k)=ze(k)-ze(k-1)
      end do

      do k=klcl+1,kmt
         delz(k)    = ze(k)-ze(k-1)
         dlamb(k)   = exp(log(dlamb0)+2.3*(ze(k)-ze(1))/zref)
         dldzby2(k) = min(1.,dlamb(k))
         dldzby2(k) = dldzby2(k)*delz(k)/2
         uhc(k)     = (uhc(k-1)-dldzby2(k)*(uhc(k-1)-uhe(k)-uhe(k-1)))/(1+dldzby2(k))
      enddo
   else
      igo=0
      return
   end if

   !----- Calculating the in-cloud variables ----------------------------------------------!
   do k=1,kmt
      dsc(k)  = dse(k)+(uhc(k)-uhes(k))/(1+gamma(k))
      dsc0(k) = dse(k)+(uhe(2)-uhes(k))/(1+gamma(k))
      qvc(k)  = qvse(k)+gamma(k)*(uhc(k)-uhes(k))/(alvl*(1+gamma(k)))
      wlc(k)  = 0.0
   end do

   do k=klcl+1,kmt    
      wlc(k) = wlc(k-1)-(qvc(k)-qvc(k-1))-dlamb(k)                                         &
             * (qvc(k)-qve(k))*delz(k)+(c0-dlamb(k))*wlc(k-1)*delz(k)

      wlc(k) = min(qvc(klcl),wlc(k))
      wlc(k) = max(.1e-12,wlc(k))
   end do
   return
end subroutine shcu_env
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sh2mod(m1)

   use conv_coms, only :       &
        kmt,                   &   ! intent(in)
        qvct1,                 &   ! intent(out)
        rhoe,                  &   ! intent(in)
        pke,                   &   ! intent(in)
        qvct2,                 &   ! intent(out)
        qvct3,                 &   ! intent(out)
        zc,                    &   ! intent(in)
        qvct4,                 &   ! intent(out)
        zzcon,                 &   ! intent(in)
        dncon,                 &   ! intent(in)
        picon                      ! intent(in)
   use shcu_vars_const, only : &
        DTDT,                  &   ! intent(in/out)
        ALVL,                  &   ! intent(in) ! parameter
        DRDT                       ! intent(in/out)

   use mem_scratch, only : VCTR5,  & ! INTENT(IN/OUT)
        VCTR6                        ! INTENT(IN/OUT)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: m1
   !----- Local variables -----------------------------------------------------------------!
   integer             :: K, TFTC, TFRC, TFTM, TFRM, FTRES, FRRES
   !----- External Function to Sum a array ------------------------------------------------!
   real, external :: ssum
   !---------------------------------------------------------------------------------------!

   !----- Compute integrated heating and moistening tendencies ----------------------------!
   do k=2,kmt
      qvct1(k) = rhoe(k)*dtdt(k)*pke(k)
      qvct2(k) = rhoe(k)*alvl*drdt(k)
      qvct3(k) = (zc(k)-zc(k-1))*qvct1(k)
      qvct4(k) = (zc(k)-zc(k-1))*qvct2(k)
   end do
   tftc=ssum(kmt-1,qvct3(2),1)
   tfrc=ssum(kmt-1,qvct4(2),1)

   !----- Transfer tendencies to model grid -----------------------------------------------!
   call vertmap2(qvct1,zc,kmt,vctr5,zzcon,m1-1)
   call vertmap2(qvct2,zc,kmt,vctr6,zzcon,m1-1)

   do k=2,m1-1
      vctr5(k)=vctr5(k)*(zzcon(k)-zzcon(k-1))
      vctr6(k)=vctr6(k)*(zzcon(k)-zzcon(k-1))
   end do

   !----- Make sure the transfer from the convective grid to the model happened fine ------!
   tftm=ssum(m1-2,vctr5(2),1)
   tfrm=ssum(m1-2,vctr6(2),1)
   ftres=tftm-tftc
   frres=tfrm-tfrc
   if (abs(ftres).gt..01*abs(tftc)) then
      print*,' Energy error in grid tranfser in convective param.'
      print*,' TFTM,TFTC ',tftm,tftc
      print*,'ERRO SHCU'
   end if
  
   !----- Change energy tendencies to temperature and mixing ratio tendencies. ------------!
   do k=2,m1-1
      dtdt(k)=vctr5(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*picon(k))
      drdt(k)=vctr6(k)/((zzcon(k)-zzcon(k-1))*dncon(k)*alvl)
   end do

   return
end subroutine sh2mod
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine cl_top()

   use conv_coms, only :       &
        nkp,                   &   ! intent(in) ! parameter
        kmt,                   &   ! intent(in)
        te,                    &   ! intent(in)
        ze,                    &   ! intent(in)
        klcl,                  &   ! intent(in)
        igo                        ! intent(out)
   use shcu_vars_const, only : &
        dscv,                  &   ! intent(out)
        dsc,                   &   ! intent(in)
        cp,                    &   ! intent(in) ! parameter
        qvc,                   &   ! intent(in)
        wlc,                   &   ! intent(in)
        dsev,                  &   ! intent(out)
        dse,                   &   ! intent(in)
        qve,                   &   ! intent(in)
        dsc0v,                 &   ! intent(out) ! maybe local var.?
        dsc0,                  &   ! intent(in)
        dsc0vm,                &   ! intent(out) ! maybe local var.?
        entf,                  &   ! intent(in)
        g,                     &   ! intent(in)  ! parameter
        cape,                  &   ! intent(out)
        delz,                  &   ! intent(in)
        ktop                       ! intent(out)

   implicit none

   !----- Local variables -----------------------------------------------------------------!
   real, dimension(nkp) :: dscvm, dsevm, tem, emp
   real                 :: gammad, buoy1, buoy2
   integer              :: k, l
   !---------------------------------------------------------------------------------------!

   gammad=.00976

   !---------------------------------------------------------------------------------------!
   !     Determination of the cloud top based on integrated cloud buoyancy. First,         !
   ! determination of virtual static energy profiles                                       !
   !---------------------------------------------------------------------------------------!
   do K=1,KMT
      !------------------------------------------------------------------------------------!
      !     This should use the actual virtual temperature, but that would require finding !
      ! the derivative following virtt so everything would become consistent. Left at the  !
      ! "to-do" list for now...                                                            !
      !------------------------------------------------------------------------------------!
      dscv(k)  = dsc(k)+cp*te(k)*(.608*qvc(k)-wlc(k))
      dsev(k)  = dse(k)+.608*cp*te(k)*qve(k)
      dsc0v(k) = dsc0(k)+.608*cp*te(k)*qve(k)
   end do

   do k=2,kmt
      dscvm(k)=(dscv(k)+dscv(k-1))/2
      dsevm(k)=(dsev(k)+dsev(k-1))/2
      dsc0vm(k)=(dsc0v(k)+dsc0v(k-1))/2
      tem(k)=(te(k)+te(k-1))/2
   end do
   !---------------------------------------------------------------------------------------!
   !    Determination of the integrated bouyancy between the surface and the LCL (see      !
   ! Albrecht et al. 1986, Eq. A5). ENTF is the enthalpy flux at surface in (K*m/s).       !
   !---------------------------------------------------------------------------------------!
   if(entf > 0.0) then
      buoy1 = entf*(1+.608*qve(1))
      buoy1 = g*ze(klcl)*buoy1/te(1)
      buoy1 = 1.3333*buoy1**.6667
   else
      buoy1=.0
      igo=0
      return
   end if

   !---------------------------------------------------------------------------------------!
   !    Checking if the parcel is able to sustain positive buoyancy one level above the    !
   ! LCL.                                                                                  !
   !---------------------------------------------------------------------------------------!
   cape  = .0
   buoy2 = .0
   l     = klcl+1
   buoy2 = buoy2+gammad*(dscvm(l)-dsevm(l))*delz(l)/tem(l)
                           
   if((BUOY1+BUOY2) <= 0.0) then
      IGO=0
      return
   end if

   !---------------------------------------------------------------------------------------!
   !     Integration continues till the level of zero buoyancy is found. The cloud top is  !
   ! assumed to be one level below.                                                        !
   !---------------------------------------------------------------------------------------!
   zerobuoyloop: do k=klcl+2,kmt                 
      buoy2=buoy2+gammad*(dscvm(k)-dsevm(k))*delz(k)/tem(k)
      if((buoy1+buoy2).le.0.0) then
         ktop=k-1
         exit zerobuoyloop
      endif
   end do zerobuoyloop

   do k=klcl,ktop
      emp(k)=dsc0vm(k)-dsevm(k)
      emp(k)=max(0.,emp(k))
      cape=cape+gammad*emp(k)*delz(k)/tem(k)
   end do
   return
end subroutine CL_TOP
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine W_SHALLOW(IP,JP,TIME)
   !---------------------------------------------------------------------------------------!
   !    The vertical velocity at cloud base is calculated according to the heat engine     !
   ! framework as defined by Renno' and Ingersoll, 1996 Eq(42)                             !
   !---------------------------------------------------------------------------------------!

   use conv_coms, only :       &
        te,                    &   ! intent(in)
        klcl,                  &   ! intent(in)
        igo,                   &   ! intent(out)
        rhoe,                  &   ! intent(in)
        ze                         ! intent(in)
   use shcu_vars_const, only : &
        ktop,                  &   ! intent(in)
        efic,                  &   ! intent(out) ! maybe local var.?
        cape,                  &   ! intent(in)
        cp,                    &   ! intent(in)
        entf,                  &   ! intent(in)  ! parameter
        alvl,                  &   ! intent(in)  ! parameter
        alhf,                  &   ! intent(in)
        dcape,                 &   ! intent(out) ! maybe local var.?
        tcape,                 &   ! intent(out) ! maybe local var.?
        wc                         ! intent(out)

   implicit none

   !---- Arguments ------------------------------------------------------------------------!
   integer     , intent(in) :: ip, jp
   real(kind=8), intent(in) :: time
   !----- Local variables -----------------------------------------------------------------!
   real :: tcold, thot, fin, sigwshb
   integer :: icount, k
   !---------------------------------------------------------------------------------------!

   tcold=te(klcl)
   icount=1
   do k=klcl+1,ktop
      tcold=tcold+te(k)
      icount=icount+1
   end do

   !---------------------------------------------------------------------------------------!
   !    TCOLD is the temperature of the environment, averaged between the first level and  !
   ! the cloud top.                                                                        !
   !---------------------------------------------------------------------------------------!
   tcold = tcold/float(icount)
   thot  = te(2)
   efic  = (thot-tcold)/thot

   if(efic <= 0.0 .or. cape < 20.0) then
      igo=0 
      return
   end if

   !---------------------------------------------------------------------------------------!
   !    The effective vertical velocity at cloud base is calculated  according to the heat !
   ! engine framework as deffined by Renno' and Ingersoll, 1996 Eq.(34)                    !
   !---------------------------------------------------------------------------------------!
   fin=rhoe(2)*(cp*entf+alvl*alhf)

   if(fin <= 50.0) then
      igo=0
      return
   end if

   dcape=0.
   tcape=cape+dcape     
   if(tcape < 50.0) then
      igo=0
      return
   end if

   sigwshb = efic*fin/(rhoe(klcl)*tcape)
   if(sigwshb <= 0.0) then
      igo=0
      return
   end if
   !---------------------------------------------------------------------------------------!
   !    Calculating the effetive vertical velocity inside the cloud interpolating from     !
   ! SIGWSHB in the cloud base to zero on the cloud top.                                   !
   !---------------------------------------------------------------------------------------!
   do k=klcl,ktop
      wc(k)=sigwshb*((ze(ktop)-ze(k))/(ze(ktop)-ze(klcl)))
   end do
   return
end subroutine w_shallow
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine SH_RATES

   use conv_coms, only :       &
           nkp,                & ! intent(in) ! parameter
           kmt,                & ! intent(in)
           klcl,               & ! intent(in)
           ze,                 & ! intent(in)
           pke                 ! ! intent(in)
   use shcu_vars_const, only : &
           drdt,               & ! intent(out)
           dtdt,               & ! intent(out)
           ktop,               & ! intent(in)
           wc,                 & ! intent(in/out)
           dsc,                & ! intent(in)
           alvl,               & ! intent(in)  ! parameter
           wlc,                & ! intent(in)
           dse,                & ! intent(in)
           qvc,                & ! intent(in)
           qve,                & ! intent(in)
           dqdt,               & ! intent(out) ! maybe local var.?
           alhf,               & ! intent(in)
           delz                ! ! intent(in)

   implicit none

   !----- Local variables -----------------------------------------------------------------!
   real, dimension(nkp) :: wssc, wqsc, dsdt
   real                 :: wrsup, wrbase
   integer              :: k, isubcl
   !---------------------------------------------------------------------------------------!

   do k=1,kmt
      wssc(k)=0.
      wqsc(k)=0.
      dsdt(k)=0.
      drdt(k)=0.
      dtdt(k)=0.
   end do

   !----- Calculating the transports w's' and w'r' ----------------------------------------!
   do k=klcl+1,ktop
      wssc(k) = wc(k)*(dsc(k)-alvl*wlc(k)-dse(k))
      wqsc(k) = wc(k)*(qvc(k)+wlc(k)-qve(k))
   end do

   !---------------------------------------------------------------------------------------!
   !    Calculating heating and moistening rates due to shallow non-precipitating cumulus. !
   !---------------------------------------------------------------------------------------!
   do k=klcl+1,ktop-1
      dqdt(k) = -(wqsc(k+1)-wqsc(k-1))/(ze(k+1)-ze(k-1))
      drdt(k) = dqdt(k)/(1-qve(k))**2
      dsdt(k) = -(wssc(k+1)-wssc(k-1))/(ze(k+1)-ze(k-1))
      dtdt(k) = dsdt(k)/pke(k)
   end do

   !---------------------------------------------------------------------------------------!
   !     ISUBCL is a flag corresponding to removal (=1) or not (=0) of moisture from the   !
   ! boudary layer by shallow clouds.                                                      !
   !---------------------------------------------------------------------------------------!
   isubcl=0
   if(isubcl == 0)return

   !---------------------------------------------------------------------------------------!
   !     For calculating the removal of moisture from the PBL, we consider that the latent !
   ! heat flux is linearly interpolated between the surface and the cloud base.            !
   !---------------------------------------------------------------------------------------!
   wrsup=alhf
   wrbase=wc(klcl)*(qvc(klcl)-qve(klcl))

   do k=1,klcl
      wc(k)  =(wrbase*ze(k)+(wrsup*(ze(klcl)-ze(k))))/(qve(k)*ze(klcl))
      wqsc(k)= wc(k)*(qvc(k)-qve(k))
   end do

   do k=2,klcl
      drdt(k)=-(wqsc(k)-wqsc(k-1))/delz(k)
   end do

   return
end subroutine sh_rates
!==========================================================================================!
!==========================================================================================!
