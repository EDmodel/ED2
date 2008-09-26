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
!     This part used to be in micro, but Harrington radiation needs to have the tables     !
! filled before the microphysics scheme is called by the 1st time, so this is now called   !
! right before the timeloop.                                                               !
!------------------------------------------------------------------------------------------!
subroutine micro_1st()
   use therm_lib  , only: bulk_on
   use mem_radiate, only: ilwrtyp,iswrtyp,icumfdbk
   use micro_coms , only: pcp_tab, nullify_sedimtab, alloc_sedimtab
   use node_mod   , only: mmzp
   use mem_grid   , only: ngrids
   use micphys    , only: maxkfall,nembfall,nhcat
   implicit none
   !----- Local variable ------------------------------------------------------------------!
   integer       :: ifm
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     This should be needed if one of the following situations occur                    !
   ! 1. Bulk microphysics is activated;                                                    !
   ! 2. Harrington radiation (either SW or LW) will be run with cumulus feedback;          !
   !---------------------------------------------------------------------------------------!
   if (bulk_on .or. ((iswrtyp == 3 .or. ilwrtyp == 3) .and. icumfdbk == 1)) then
      call micinit()
      call make_autotab()
      call haznuc()
      call tabmelt()
      call tabhab()

      !----- Allocating the precipitation table -------------------------------------------!
      allocate(pcp_tab(ngrids))
      do ifm=1,ngrids
         call nullify_sedimtab(pcp_tab(ifm))
         call alloc_sedimtab(pcp_tab(ifm),mmzp(ifm),maxkfall,nembfall,nhcat)
      end do

   end if

   return
end subroutine micro_1st
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine micro_master(dn01dntop)

   use therm_lib , only : &
           bulk_on        ! ! intent(in)

   use micphys   , only : &
           cparm          & !intent(in)
          ,rparm          & !intent(in)
          ,pparm          & !intent(in)
          ,sparm          & !intent(in)
          ,aparm          & !intent(in)
          ,gparm          & !intent(in)
          ,hparm          & !intent(in)
          ,icloud         & !intent(in)
          ,irain          & !intent(in)
          ,ipris          & !intent(in)
          ,isnow          & !intent(in)
          ,iaggr          & !intent(in) 
          ,igraup         & !intent(in)
          ,ihail          & !intent(in)
          ,parm           & !intent(out)
          ,nhcat          & !intent(in)
          ,cfmas          & !intent(out)
          ,pwmas          & !intent(out)
          ,cfvt           & !intent(out)
          ,pwvt           & !intent(out)
          ,ipairc         & !intent(out)
          ,ipairr         & !intent(out)
          ,ncat           & !intent(in)
          ,emb0           & !intent(out)
          ,emb1           & !intent(out)
          ,emb2           & !intent(out)
          ,rxmin          & !intent(in)
          ,cxmin          & !intent(out)
          ,mkcoltab       & !intent(in)
          ,coltabfn       & !intent(in)
          ,gnu            & !intent(inout)
          ,npairc         & !intent(in)
          ,coltabc        & !intent(inout)
          ,nembc          & !intent(inout)
          ,npairr         & !intent(in)
          ,coltabr        & !intent(inout)
          ,progncat       & !intent(in)
          ,shapefac       ! !intent(out)

   use micro_coms, only : &
          dstprms         & ! intent(in)
         ,lcat_lhcat      & ! intent(in)
         ,jpairr          & ! intent(in)
         ,jpairc          ! ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   real             , intent(in)  :: dn01dntop
   !----- Local variables -----------------------------------------------------------------!
   integer                        :: lhcat,khcat,lcat,nd1,nd2,nip,ilcat,ilhcat,idum
   character(len=80)              :: dataline,cname
   real                           :: pwen0_tmp, glg, glgm, cfen0_tmp
   !----- Functions -----------------------------------------------------------------------!
   real, external                 :: gammln
   !---------------------------------------------------------------------------------------!

   !----- Initialize parm with prescribed values ------------------------------------------!
   parm(1) = cparm
   parm(2) = rparm
   parm(3) = 1.e-4 !----- Used only for cumulus feedback on radiation ---------------------!
   parm(4) = sparm
   parm(5) = aparm
   parm(6) = gparm
   parm(7) = hparm

   if (icloud <= 1) parm(1) = dstprms(8,1)
   if (irain  == 1) parm(2) = dstprms(8,2)
   if (isnow  == 1) parm(4) = dstprms(8,4)
   if (iaggr  == 1) parm(5) = dstprms(8,5)
   if (igraup == 1) parm(6) = dstprms(8,6)
   if (ihail  == 1) parm(7) = dstprms(8,7)

   !----- Define several parameters from above data list ----------------------------------!
   do lhcat=1,nhcat
      shapefac(lhcat) = dstprms(1,lhcat)
      cfmas(lhcat)    = dstprms(2,lhcat)
      pwmas(lhcat)    = dstprms(3,lhcat)
      cfvt (lhcat)    = dstprms(4,lhcat)
      pwvt (lhcat)    = dstprms(5,lhcat)

      do khcat=1,nhcat
         ipairc(lhcat,khcat) = jpairc(lhcat,khcat)
         ipairr(lhcat,khcat) = jpairr(lhcat,khcat)
      end do
   end do

   do lcat=1,ncat
      emb0 (lcat) = cfmas(lcat) * dstprms(6,lcat) ** pwmas(lcat)
      emb1 (lcat) = cfmas(lcat) * dstprms(7,lcat) ** pwmas(lcat)
      rxmin(lcat) = dstprms(9,lcat)
   enddo

   do lhcat=1,nhcat
      lcat = lcat_lhcat(lhcat)
      if (lcat == 1 .and. progncat(1)) then
         emb2 (lhcat) = cfmas(lhcat) * dstprms(8,1) ** pwmas(lhcat)
      else
         emb2 (lhcat) = cfmas(lhcat) * parm(lcat) ** pwmas(lhcat)
      end if

      !------------------------------------------------------------------------------------!
      !    Initialising the minimum cx for each category. This is just a lower bound so we !
      ! don't need to be too picky.                                                        !
      !------------------------------------------------------------------------------------!
      !select case (jnmb(lcat))
      !case (0)
      !!----- 2: Prescribed diameter -------------------------------------------------------!
      !case (2)
      !   cxmin(lhcat) = rxmin(lcat) / emb2(lhcat)
      !!----- 3: Prescribed y intercept ----------------------------------------------------!
      !case (3)
      !   pwen0_tmp = 1. / (pwmas(lhcat) + 1.)
      !   glg       = gammln(gnu(lcat))
      !   glgm      = gammln(gnu(lcat) + pwmas(lhcat))
      !   cfen0_tmp = parm(lcat) * (exp(glg - glgm) / parm(lcat)) ** pwen0_tmp
      !   cxmin(lhcat) = (cfen0_tmp/dn01dntop) * (dn01dntop*rxmin(lcat))** pwen0_tmp
      !!----- 4: Prescribed concentration (I guess cxmin could be simply parm(lcat)... -----!
      !case (4)
      !   cxmin(lhcat) = rxmin(lcat)                                                        &
      !               / max(emb0(lcat),min(emb1(lcat),rxmin(lcat)/parm(lcat)))
      !!----- 5 and beyond are prognostic, use the prescribed only for a lower bound... ----!
      !case (5,6,7)
      !   select case (lcat)
      !   case (1)
      !      cxmin(lhcat) = rxmin(lcat)                                                     &
      !                   / max(emb0(lcat),min(emb1(lcat),rxmin(lcat)/dstprms(8,1)))
      !   case default
      !      cxmin(lhcat) = rxmin(lcat)                                                     &
      !                   / max(emb0(lcat),min(emb1(lcat),rxmin(lcat)/parm(lcat)))
      !   end select
      !end select
      cxmin(lhcat) = 1.e-9
   !---------------------------------------------------------------------------------------!
   end do

   if (.not. bulk_on) return

   cname=coltabfn(1:len_trim(coltabfn))

   select case (mkcoltab)
   case (1) !----- Make collection table and write to file --------------------------------!
      call mkcoltb()

      open(unit=91,file=cname,form='formatted',status='replace',action='write')
      rewind(unit=91)
      write (unit=91,fmt='(a5,4(1x,a14))')                                                 &
         adjustr('lcat'),adjustr('gnu'),adjustr('emb0'),adjustr('emb1'),adjustr('rxmin')
      do lcat = 1,ncat
         write (unit=91,fmt='(i5,4(1x,es14.7))') lcat,gnu(lcat),emb0(lcat),emb1(lcat)      &
                                                ,rxmin(lcat)
      end do
      write (unit=91,fmt='(a)') ' '

      write (unit=91,fmt='(a5,6(1x,a14))') adjustr('lhcat'),adjustr('cfmas')               &
                                          ,adjustr('pwmas'),adjustr('cfvt')                &
                                          ,adjustr('pwvt'),adjustr('emb2'),adjustr('cxmin')
      do lhcat = 1,nhcat
         write (unit=91,fmt='(i5,6(1x,es14.7))') lhcat,cfmas(lhcat),pwmas(lhcat)           &
                                                ,cfvt(lhcat),pwvt(lhcat),emb2(lhcat)       &
                                                ,cxmin(lhcat)
      enddo
      write (unit=91,fmt='(a)') ' '

      do nip=1,npairc
         write (unit=91,fmt='(a,1x,i4)') 'ipairc',nip
         do nd2=1,nembc
            write (unit=91,fmt='(i4,20(1x,es14.7))') nd2,(coltabc(nd1,nd2,nip),nd1=1,nembc)
         end do
      end do
      write (unit=91,fmt='(a)') ' '

      do nip=1,npairr
         write (unit=91,fmt='(a,1x,i4)') 'ipairr',nip
         do nd2=1,nembc
            write (unit=91,fmt='(i4,20(1x,es14.7))') nd2,(coltabr(nd1,nd2,nip),nd1=1,nembc)
         end do
      end do



   case (0) !----- Read collection table --------------------------------------------------!
      open (unit=91,file=trim(cname),form='formatted',status='old',action='read')

      read (unit=91,fmt=*)
      do ilcat = 1,ncat
         read (unit=91,fmt=*) lcat,gnu(lcat),emb0(lcat),emb1(lcat),rxmin(lcat)
      end do
      read (unit=91,fmt=*)

      read (unit=91,fmt=*)
      do ilhcat = 1,nhcat
         read (unit=91,fmt=*) lhcat,cfmas(lhcat),pwmas(lhcat),cfvt(lhcat),pwvt(lhcat)       &
                            ,emb2(lhcat),cxmin(lhcat)
      end do
      read (unit=91,fmt=*)

      do nip=1,npairc
         read (unit=91,fmt=*)
         read (unit=91,fmt=*) (idum,(coltabc(nd1,nd2,nip),nd1=1,nembc),nd2=1,nembc)
      end do
      read (unit=91,fmt=*)

      do nip=1,npairr
         read (unit=91,fmt=*)
         read (unit=91,fmt=*) (idum,(coltabr(nd1,nd2,nip),nd1=1,nembc),nd2=1,nembc)
      end do

   end select

   close(unit=91,status='keep')

   return
end subroutine micro_master
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine assign initial energy for hydrometeors and initial concentration of   !
! CCNs in case they are needed.                                                            !
!------------------------------------------------------------------------------------------!
subroutine initqin(n1,n2,n3,q2,q6,q7,pi0,pp,theta,dn0,cccnp,cifnp)

   use micphys, only : &
        exner,         & ! intent(out)
        tair,          & ! intent(out)
        icloud,        & ! intent(in)
        irain,         & ! intent(in)
        igraup,        & ! intent(in)
        ihail,         & ! intent(in)
        ipris            ! intent(in)
   use rconstants, only : cpi,t3ple,tsupercool !intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer,                   intent(in)  :: n1,n2,n3
   real, dimension(n1,n2,n3), intent(out) :: q2, q6, q7, cifnp, cccnp
   real, dimension(n1,n2,n3), intent(in)  :: pi0, pp, theta, dn0
   !----- Local Variables -----------------------------------------------------------------!
   integer :: i,j,k
   !---------------------------------------------------------------------------------------!


   !----- Initialize Q2, Q6, Q7, CCN, IFN. ------------------------------------------------!
   do j = 1,n3
      do i = 1,n2
         do k = 1,n1
            exner(k) = pi0(k,i,j) + pp(k,i,j)
            tair(k)  = theta(k,i,j) * exner(k) * cpi

            if (irain  >= 1) q2(k,i,j)    = tair(k) - tsupercool
            if (igraup >= 1) q6(k,i,j)    = 0.5 * min(0.,tair(k) - t3ple)
            if (ihail  >= 1) q7(k,i,j)    = 0.5 * min(0.,tair(k) - t3ple)
            !----- Making up something, but not sure about this one. ----------------------!
            if (icloud == 7) cccnp(k,i,j) = 6.66e9
            if (ipris  == 7) cifnp(k,i,j) = 1.e5 * dn0(k,i,j) ** 5.4
         end do
      end do
   end do
   return
end subroutine initqin

!******************************************************************************

subroutine jnmbinit()

   use therm_lib, only : &
           level         ! ! intent(in)
   use micphys  , only : &
           ncat          & ! intent(in)
          ,icloud        & ! intent(in)
          ,irain         & ! intent(in)
          ,ipris         & ! intent(in)
          ,isnow         & ! intent(in)
          ,iaggr         & ! intent(in)
          ,igraup        & ! intent(in)
          ,ihail         & ! intent(in)
          ,jnmb          & ! intent(out)
          ,availcat      & ! intent(out)
          ,progncat      ! ! intent(out)
   implicit none
   !----- Local Variable ------------------------------------------------------------------!
   integer :: lcat
   !---------------------------------------------------------------------------------------!

   select case (level)
   case (0,1) !----- No condensation, all hydrometeors are off ----------------------------!
      jnmb(1:7) = 0
   case (2)   !----- Cloud condensation only, diagnostic concentration on -----------------!
      jnmb(1)   = 4
      jnmb(2:7) = 0
   case (3)
      !----- User-definition --------------------------------------------------------------!
      jnmb(1) = icloud
      jnmb(2) = irain
      jnmb(3) = ipris
      jnmb(4) = isnow
      jnmb(5) = iaggr
      jnmb(6) = igraup
      jnmb(7) = ihail

      !----- i???? equals to one is the "default" combination -----------------------------!
      if (icloud == 1) jnmb(1) = 4
      if (irain  == 1) jnmb(2) = 2
      if (ipris  == 1) jnmb(3) = 5
      if (isnow  == 1) jnmb(4) = 2
      if (iaggr  == 1) jnmb(5) = 2
      if (igraup == 1) jnmb(6) = 2
      if (ihail  == 1) jnmb(7) = 2

      !------------------------------------------------------------------------------------!
      !     Redefining the prognostic species. If any of the hydrometeors other than cloud !
      ! and pristine ice is set prognostic, then turn all the hydrometeors other than      !
      ! cloud and pristine ice that are available into prognostic mode.                    !
      !------------------------------------------------------------------------------------!
      if (irain == 5 .or. isnow == 5 .or. iaggr == 5 .or. igraup == 5 .or. ihail == 5) then
         if (irain  >= 1) jnmb(2) = 5
         if (isnow  >= 1) jnmb(4) = 5
         if (iaggr  >= 1) jnmb(5) = 5
         if (igraup >= 1) jnmb(6) = 5
         if (ihail  >= 1) jnmb(7) = 5
      end if
   end select

   
   !----- Saving these logical flags to speed up checks -----------------------------------!
   do lcat = 1,ncat
      availcat(lcat) = jnmb(lcat) >= 1
      progncat(lcat) = jnmb(lcat) >= 5
   end do

   return
end subroutine jnmbinit
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine defines all parameters that were not initialised at micro_master, so  !
! this subroutine is called on each node, right before the time step loop, from mic_1st.   !
!------------------------------------------------------------------------------------------!
subroutine micinit()

   use micphys, only : &
           parm        & ! intent(out)
          ,cparm       & ! intent(in)
          ,rparm       & ! intent(in)
          ,sparm       & ! intent(in)
          ,aparm       & ! intent(in)
          ,gparm       & ! intent(in)
          ,hparm       & ! intent(in)
          ,icloud      & ! intent(in)
          ,irain       & ! intent(in)
          ,isnow       & ! intent(in)
          ,iaggr       & ! intent(in)
          ,igraup      & ! intent(in)
          ,ihail       & ! intent(in)
          ,dps         & ! intent(out)
          ,dps2        & ! intent(out)
          ,rictmin     & ! intent(out)
          ,rictmax     & ! intent(out)
          ,nembc       & ! intent(in)
          ,nhcat       & ! intent(in)
          ,cfden       & ! intent(out)
          ,cfmas       & ! intent(in)
          ,pwden       & ! intent(out)
          ,pwmas       & ! intent(in)
          ,emb0log     & ! intent(out)
          ,emb0        & ! intent(in)
          ,emb1log     & ! intent(out)
          ,emb1        & ! intent(in)
          ,pwmasi      & ! intent(out)
          ,pwen0       & ! intent(out)
          ,pwemb0      & ! intent(out)
          ,pwvt        & ! intent(in)
          ,gnu         & ! intent(in)
          ,jnmb        & ! intent(in)
          ,cfemb0      & ! intent(out)
          ,cfen0       & ! intent(out)
          ,dnfac       & ! intent(out)
          ,vtfac       & ! intent(out)
          ,cfvt        & ! intent(in)
          ,frefac1     & ! intent(out)
          ,shapefac    & ! intent(in)
          ,frefac2     & ! intent(out)
          ,sipfac      & ! intent(out)
          ,cfmasft     & ! intent(out)
          ,dict        & ! intent(out)
          ,dpsmi       & ! intent(out)
          ,gamm        & ! intent(inout)
          ,gamn1       & ! intent(out)
          ,ngam        & ! intent(in)
          ,gam         & ! intent(out)
          ,gaminc      & ! intent(out)
          ,gamsip13    & ! intent(out)
          ,gamsip24    & ! intent(out)
          ,cfmasi      & ! intent(out)
          ,pwmasi      & ! intent(out)
          ,ch3         & ! intent(out)
          ,cdp1        & ! intent(out)
          ,pwvtmasi    ! ! intent(out)
   use rconstants, only: pio6i
   use micro_coms, only : lcat_lhcat
   
   implicit none

   !----- Local Variables: ----------------------------------------------------------------!
   integer :: lhcat,lcat,ia
   real :: c1,glg,glg1,glg2,glgm,glgc,glgmv,flngi,dpsi,embsip,dnsip
   real :: gammln,gammp,gammq
   real :: aux_loop1, aux_loop2
   !---------------------------------------------------------------------------------------!



   !----- Initialize arrays based on microphysics namelist parameters ---------------------!

   dps = 125.e-6
   dps2 = dps ** 2
   rictmin = 1.0001
   rictmax = 0.9999 * real(nembc)

   do lhcat = 1,nhcat
      lcat = lcat_lhcat(lhcat)

      cfden(lhcat)  = cfmas(lhcat) * pio6i
      pwden(lhcat)  = pwmas(lhcat) - 3.
      emb0log(lcat) = log(emb0(lcat))
      emb1log(lcat) = log(emb1(lcat))

      !------------------------------------------------------------------------------------!
      !    Define coefficients [vtfac, frefac1, frefac2] used for terminal velocity and    !
      ! Reynolds number.                                                                   !
      !------------------------------------------------------------------------------------!
      cfmasi(lhcat) = 1. / cfmas(lhcat)
      pwmasi(lhcat) = 1. / pwmas(lhcat)
      pwen0(lhcat)  = 1. / (pwmas(lhcat) + 1.)
      pwemb0(lhcat) = pwmas(lhcat) / (pwmas(lhcat) + 1.)
      c1            = 1.5 + .5 * pwvt(lhcat)

      glg   = gammln(gnu(lcat))
      glg1  = gammln(gnu(lcat) + 1.)
      glg2  = gammln(gnu(lcat) + 2.)
      glgm  = gammln(gnu(lcat) + pwmas(lhcat))
      glgc  = gammln(gnu(lcat) + c1)
      glgmv = gammln(gnu(lcat) + pwmas(lhcat) + pwvt(lhcat))

      if (jnmb(lcat) == 3) then
         cfemb0(lhcat) = cfmas(lhcat) * exp(glgm - glg)                                    &
                      ** pwen0(lhcat) * (1. / parm(lcat)) ** pwemb0(lhcat)
         cfen0(lhcat) = parm(lcat) * (exp(glg - glgm) / parm(lcat)) ** pwen0(lhcat)
      end if

      dnfac(lhcat) = (cfmasi(lhcat) * exp(glg - glgm)) ** pwmasi(lhcat)

      vtfac(lhcat) = cfvt(lhcat) * exp(glgmv - glgm)                                       &
                   * (cfmasi(lhcat) * exp(glg - glgm)) ** (pwvt(lhcat) *pwmasi(lhcat))

      frefac1(lhcat) = shapefac(lhcat) * exp(glg1 - glg)                                  &
                     * (cfmasi(lhcat) * exp(glg - glgm)) ** pwmasi(lhcat)

      frefac2(lhcat) = shapefac(lhcat) * 0.229 * sqrt(cfvt(lcat))                          &
                     * (cfmasi(lhcat) * exp(glg - glgm)) ** (pwmasi(lhcat) * c1)           &
                     * exp(glgc - glg)

      sipfac(lhcat) = .785 * exp(glg2 - glg)                                               &
                    * (cfmasi(lhcat) * exp(glg - glgm)) ** (2. * pwmasi(lhcat))

      cfmasft(lhcat) = cfmas(lhcat)                                                        &
                     * exp(gammln(gnu(lcat) + pwmas(lhcat)) - gammln(gnu(lcat)))

      dict(lcat) = float(nembc-1) / (emb1log(lcat) - emb0log(lcat))

      dpsmi(lhcat) = 1. / (cfmas(lhcat) * dps ** pwmas(lhcat))
      if (lhcat <= 4) gamm(lhcat)  = exp(glg)
      if (lhcat <= 4) gamn1(lhcat) = exp(glg1)


      ch3(lhcat)      = pwvt(lhcat)   * pwmasi(lhcat)
      cdp1(lhcat)     = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
      pwvtmasi(lhcat) = pwvt(lhcat)   * pwmasi(lhcat)

      ! gam1   :  the integral of the pristine distribution from dps to infty
      ! gam2   :  the integral of the snow dist. from 0 to dps
      ! gam3   :  values of the exponential exp(-dps/dn)

   end do

   flngi = 1. / float(ngam)

   ! ALF
   aux_loop1 = dps * 1.e6
   aux_loop2 = emb1(1) * flngi

   do ia=1,ngam
      dpsi = aux_loop1 / float(ia)

      gam(ia,1) = gammq(gnu(3) + 1., dpsi)
      gam(ia,2) = gammp(gnu(4) + 1., dpsi)
      gam(ia,3) = exp(-dpsi)

      gaminc(ia,1)=gammq(gnu(3),dpsi)
      gaminc(ia,2)=gammp(gnu(4),dpsi)

      embsip       = aux_loop2 * float(ia)
      dnsip        = dnfac(1) * embsip ** pwmasi(1)
      gamsip13(ia) = gammp(gnu(1),13.e-6/dnsip)
      gamsip24(ia) = gammq(gnu(1),24.e-6/dnsip)
   end do

   return
end subroutine micinit
!==========================================================================================!
!==========================================================================================!

