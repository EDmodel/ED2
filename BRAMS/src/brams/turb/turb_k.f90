!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine diffuse()

  ! +-----------------------------------------------------------------+
  ! \	this routine is the subdriver to compute tendencies due to    \
  ! \	  subgrid-scale turbulence.				      \
  ! +-----------------------------------------------------------------+

  use mem_tend, only:    &
       tend                   ! %tket, %epst, %ut, %vt, %wt

  use mem_basic, only:   &
       basic_g                !  %up(IN), %vp(IN)

  use var_tables, only:  &
       num_scalar,       &    ! INTENT(IN)
       scalar_tab             ! %var_p, %var_t

  use mem_turb, only:    &
       idiffk,           &    !INTENT(IN)
       turb_g,           &    ! %tkep, %hkm, %vkh
       xkhkm

  use mem_grid, only:     &
       if_adap,           &  !           INTENT(IN)
       jdim,              &  !           INTENT(IN)
       ngrid,             &  !           INTENT(IN)
       grid_g,            &  ! %rtgt     INTENT(IN)
       dzm,               &  !           INTENT(IN)
       dzt,               &  !           INTENT(IN)
       npatch,            &  !           INTENT(IN)
       nstbot,            &  !           INTENT(IN)
       nscl,              &  !           INTENT(IN)
       naddsc,            &  !           INTENT(IN)
       dtlt               !  !           INTENT(IN)

  use mem_leaf, only:     &
       leaf_g                  !INTENT(IN)

  use mem_micro, only:    &
       micro_g                 !%rcp

  use mem_scratch, only:  &
       scratch,           &    ! %vt3da, %vt3db, %vt3dc, %vt3dd, %vt3de,
                               ! %vt3df, %vt3dg, %vt3dh, %vt3di, %vt3dj,
                               ! %vt3dn, %scr2
       vctr34                  !

  use node_mod, only:     &
       mxp,               &  !INTENT(IN)
       myp,               &  !INTENT(IN)
       mzp,               &  !INTENT(IN)
       ia,                &  !INTENT(IN)
       iz,                &  !INTENT(IN)
       ja,                &  !INTENT(IN)
       jz,                &  !INTENT(IN)
       ia_1,              &  !INTENT(IN)
       ja_1,              &  !INTENT(IN)
       iz1,               &  !INTENT(IN)
       jz1,               &  !INTENT(IN)
       ibcon,             &  !INTENT(IN)
       mynum,             &  !INTENT(IN)
       mi0,               &  !INTENT(IN)
       mj0,               &  !INTENT(IN)
       mmyp,              &  !INTENT(IN)
       mmxp,              &  !INTENT(IN)
       mmzp,              &  !INTENT(IN)
       izu,               &  !INTENT(IN)
       jzv,               &  !INTENT(IN)
       ia1,               &  !INTENT(IN)
       ja1,               &  !INTENT(IN)
       iz_1,              &  !INTENT(IN)
       jz_1                  !INTENT(IN)

  use ke_coms, only:      &
       alf_eps,           &
       alf_tke

  use micphys, only:      &
       level

  use mem_turb_scalar, only:   &
       turb_s

  use catt_start, only:        &
       CATT           ! intent(in)
  
  use mem_cuparm, only : &
       nnqparm

  ! Mass flux
  use mem_mass, only: mass_g, imassflx

  implicit none

  !local variables:
  integer :: mxyzp,ind,n
  real :: s1,s2,s3
  real, pointer, dimension(:) :: vkh_p,hkh_p
  integer :: i,j,k,ksf

  !###########################################################################
  ! CATT
  ! srf - fev-2003 - Large Scale Forcing for shallow and deep cumulus
  real, pointer :: lsfcupar_p(:,:,:)
    
  !#########################################################################

  ! CATT
  ! Nullifing pointer to Large Scale Forcing for GRELL CUPAR - Not used
  nullify(lsfcupar_p)
  !

  mxyzp = mxp * myp * mzp
  ![MLO - Nullifying and allocating pointers for the coefficients
  nullify(vkh_p)
  nullify(hkh_p)

  scratch%vt3dg = 0. !LFR - Avoiding overflow in strain

  if (if_adap == 0) then

     call strain(mzp,mxp,myp,ia,iz,ja,jz                       &
          ,ia_1,ja_1,iz1,jz1,jdim                                &
          ,basic_g(ngrid)%up (1,1,1) ,basic_g(ngrid)%vp (1,1,1)  &
          ,basic_g(ngrid)%wp (1,1,1) ,scratch%vt3da     (1)      &
          ,scratch%vt3db     (1)     ,scratch%vt3dc     (1)      &
          ,scratch%vt3dd     (1)     ,scratch%vt3de     (1)      &
          ,scratch%vt3df     (1)     ,scratch%vt3dg     (1)      &
          ,scratch%vt3dh     (1)     ,scratch%vt3di     (1)      &
          ,scratch%vt3dn     (1)     ,scratch%scr2      (1)      &
          ,idiffk(ngrid))

  else

     call strain_adap(mzp,mxp,myp,ia,iz,ja,jz                  &
          ,ia_1,ja_1,iz1,jz1,jdim                                &
          ,grid_g(ngrid)%flpu (1,1)   ,grid_g(ngrid)%flpv (1,1)    &
          ,grid_g(ngrid)%flpw (1,1)   ,basic_g(ngrid)%up (1,1,1)  &
          ,basic_g(ngrid)%vp (1,1,1) ,basic_g(ngrid)%wp (1,1,1)  &
          ,scratch%vt3da     (1)     ,scratch%vt3db     (1)      &
          ,scratch%vt3dc     (1)     ,scratch%vt3dd     (1)      &
          ,scratch%vt3de     (1)     ,scratch%vt3df     (1)      &
          ,scratch%vt3dg     (1)     ,scratch%vt3dh     (1)      &
          ,scratch%vt3di     (1)     ,scratch%vt3dn     (1)      &
          ,scratch%scr2      (1)     ,idiffk(ngrid)              &
          ,grid_g(ngrid)%dxm (1,1)   ,grid_g(ngrid)%dxt (1,1)    &
          ,grid_g(ngrid)%dxu (1,1)   ,grid_g(ngrid)%dxv (1,1)    &
          ,grid_g(ngrid)%dym (1,1)   ,grid_g(ngrid)%dyt (1,1)    &
          ,grid_g(ngrid)%dyu (1,1)   ,grid_g(ngrid)%dyv (1,1)    &
          ,dzm,dzt)

  endif

  if (level <= 1) call azero(mxyzp,scratch%vt3dp(:))
  if (level >= 2) call ae1  (mxyzp,scratch%vt3dp(:),micro_g(ngrid)%rcp(:,:,:))

  call bruvais(mzp,mxp,myp,ia,iz,ja,jz                          &
       ,basic_g(ngrid)%theta (1,1,1) ,basic_g(ngrid)%rtp (1,1,1)  &
       ,basic_g(ngrid)%rv    (1,1,1) ,scratch%vt3dp(1)            &
       ,basic_g(ngrid)%pp    (1,1,1) ,basic_g(ngrid)%pi0 (1,1,1)  &
       ,scratch%vt3dj        (1)     ,grid_g(ngrid)%rtgt (1,1)    &
       ,grid_g(ngrid)%flpw    (1,1)   )

  if (idiffk(ngrid) <= 3 .or. idiffk(ngrid) == 7) then
     call mxdefm(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim            &
          ,scratch%vt3dh      (1)     ,scratch%vt3di      (1)    &
          ,scratch%vt3dj      (1)     ,scratch%vt3dk      (1)    &
          ,scratch%scr1       (1)     ,scratch%scr2       (1)    &
          ,basic_g(ngrid)%dn0 (1,1,1) ,grid_g(ngrid)%rtgt (1,1)  &
          ,grid_g(ngrid)%dxt  (1,1)   ,grid_g(ngrid)%dyt  (1,1)  &
          ,grid_g(ngrid)%flpw  (1,1)   ,mynum  )

     ! CATT
     if (CATT == 1) then
        !srf------
        !coef de difusao horizontal diferente (dum4 e ivt3dp) para tracers
        call mxdefm_tracer(mzp,mxp,myp,ia,iz,ja,jz  &
             ,ibcon,jdim,scratch%vt3dh(1),scratch%scr3(1) &
             ,basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%dxt(1,1),&
             grid_g(ngrid)%dyt(1,1),grid_g(ngrid)%flpw(1,1),mynum)

        !srf------
     endif

  endif


! Original Mellor-Yamada Level 2½, following Helfand and Labraga (1988)
  if (idiffk(ngrid) == 1) then
     call tkemy(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim,mi0(ngrid),mj0(ngrid)  &
          ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
          ,scratch%vt3dh        (1)     ,scratch%vt3di        (1)      &
          ,scratch%vt3dj        (1)     ,scratch%scr1         (1)      &
          ,grid_g(ngrid)%rtgt   (1,1)   ,basic_g(ngrid)%theta (1,1,1)  &
          ,basic_g(ngrid)%dn0   (1,1,1) ,basic_g(ngrid)%up    (1,1,1)  &
          ,basic_g(ngrid)%vp    (1,1,1) ,basic_g(ngrid)%wp    (1,1,1)  &
          ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
          ,turb_g(ngrid)%sflux_w(1,1)   ,turb_g(ngrid)%sflux_t(1,1),vctr34 &
          ,grid_g(ngrid)%flpw    (1,1)   ,grid_g(ngrid)%flpu    (1,1)   &
          ,grid_g(ngrid)%flpv    (1,1))
! Deardoff (1980) LES scheme
  elseif (idiffk(ngrid) == 4) then
     call mxtked(mzp,mxp,myp,ia,iz,ja,jz  &
          ,ibcon,jdim  &
          ,turb_g(ngrid)%tkep   (1,1,1) ,tend%tket            (1)      &
          ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
          ,basic_g(ngrid)%wp    (1,1,1) ,basic_g(ngrid)%rtp   (1,1,1)  &
          ,basic_g(ngrid)%rv    (1,1,1) ,basic_g(ngrid)%theta (1,1,1)  &
          ,scratch%vt3da        (1)     ,scratch%vt3dc        (1)      &
          ,scratch%vt3dh        (1)     ,scratch%vt3dj        (1)      &
          ,scratch%scr1         (1)     ,scratch%scr2         (1)      &
          ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
          ,turb_g(ngrid)%sflux_w(1,1)   ,turb_g(ngrid)%sflux_t(1,1)    &
          ,grid_g(ngrid)%dxt    (1,1)   ,grid_g(ngrid)%rtgt   (1,1)    &
          ,grid_g(ngrid)%flpw    (1,1)   )
  ![STC..................................................
  !_STC    Note: from subroutines TKESCL, TKEEPS :
  !_STC           VT3DI=Ke
  !_STC           SCR1=Km
  !_STC           VT3DH = Kh
  !_STC           SCR2 = SCR1 = Km
  !_STC..................................................
  !_STC............................................................
  !_STC............................................................
  !_STC Call to subroutine tkescl for E-l closure
  !_STC (S. Trini Castelli)
  !_STC............................................................
  elseif (idiffk(ngrid) == 5) then
     call tkescl(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
          ,turb_g(ngrid)%tkep(1,1,1),tend%tket(1)  &
          ,turb_g(ngrid)%epsp(1,1,1)  &
          ,scratch%vt3da(1),scratch%vt3dc(1)  &
          ,scratch%vt3dh(1),scratch%vt3di(1)  &
          ,scratch%vt3dj(1),scratch%scr1(1)  &
          ,scratch%scr2(1) ,grid_g(ngrid)%rtgt(1,1)  &
          ,scratch%vt3dd(1),scratch%vt3de(1),grid_g(ngrid)%dxt(1,1)  &
          ,leaf_g(ngrid)%ustar(1,1,1),leaf_g(ngrid)%patch_area(1,1,1) &
          ,grid_g(ngrid)%flpw(1,1),basic_g(ngrid)%dn0(1,1,1)  )
  !_STC............................................................
  !_STC Call to subroutine tkeeps for E-eps closure
  !_STC (S. Trini Castelli)
  !_STC............................................................
  elseif (idiffk(ngrid) == 6) then
     call tkeeps(mzp,mxp,myp,npatch,ia,iz,ja,jz  &
          ,turb_g(ngrid)%tkep(1,1,1),tend%tket(1)  &
          ,turb_g(ngrid)%epsp(1,1,1),tend%epst(1)  &
          ,scratch%vt3da(1),scratch%vt3dc(1)  &
          ,scratch%vt3dh(1),scratch%vt3di(1)  &
          ,scratch%vt3dj(1),scratch%scr1(1)  &
          ,scratch%scr2(1) ,grid_g(ngrid)%rtgt(1,1)  &
          ,leaf_g(ngrid)%ustar(1,1,1),leaf_g(ngrid)%patch_area(1,1,1) &
          ,grid_g(ngrid)%flpw(1,1),basic_g(ngrid)%dn0(1,1,1)  )
![MLO -> Nananishi and Niino (2004), scheme based on Mellor-Yamada Level 2½
  elseif (idiffk(ngrid) == 7) then
     if(level >= 1) then
       call nakanishi(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim,level         &
          ,turb_g(ngrid)%tkep       (1,1,1)   ,tend%tket                (1)     &
          ,scratch%vt3dd            (1)       ,scratch%vt3de            (1)     &
          ,scratch%vt3dh            (1)       ,scratch%vt3di            (1)     &
          ,scratch%vt3dj            (1)       ,scratch%scr1             (1)     &
          ,grid_g(ngrid)%rtgt       (1,1)     ,basic_g(ngrid)%theta     (1,1,1) &
          ,basic_g(ngrid)%rv        (1,1,1)   ,basic_g(ngrid)%dn0       (1,1,1) &
          ,basic_g(ngrid)%up        (1,1,1)   ,basic_g(ngrid)%vp        (1,1,1) &
          ,leaf_g(ngrid)%veg_rough  (1,1,1)   ,leaf_g(ngrid)%patch_rough(1,1,1) &
          ,leaf_g(ngrid)%tstar      (1,1,1)   ,leaf_g(ngrid)%ustar      (1,1,1) &
          ,leaf_g(ngrid)%patch_area (1,1,1)   ,turb_g(ngrid)%sflux_u    (1,1)   &
          ,turb_g(ngrid)%sflux_v    (1,1)     ,turb_g(ngrid)%sflux_t    (1,1)   &
          ,grid_g(ngrid)%flpu       (1,1)     ,grid_g(ngrid)%flpv       (1,1)   &
          ,grid_g(ngrid)%flpw       (1,1)     ,turb_g(ngrid)%kpbl       (1,1)   &
          ,turb_g(ngrid)%pblhgt     (1,1)     ,turb_g(ngrid)%lmo        (1,1)   &
          ,turb_g(ngrid)%ltscale    (1,1,1)   ,turb_g(ngrid)%sigw       (1,1,1))
     else
       call azero(mxyzp,scratch%vt3dp(1))
       call nakanishi(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim,level         &
          ,turb_g(ngrid)%tkep       (1,1,1)   ,tend%tket                (1)     &
          ,scratch%vt3dd            (1)       ,scratch%vt3de            (1)     &
          ,scratch%vt3dh            (1)       ,scratch%vt3di            (1)     &
          ,scratch%vt3dj            (1)       ,scratch%scr1             (1)     &
          ,grid_g(ngrid)%rtgt       (1,1)     ,basic_g(ngrid)%theta     (1,1,1) &
          ,scratch%vt3dp            (1)       ,basic_g(ngrid)%dn0       (1,1,1) &
          ,basic_g(ngrid)%up        (1,1,1)   ,basic_g(ngrid)%vp        (1,1,1) &
          ,leaf_g(ngrid)%veg_rough  (1,1,1)   ,leaf_g(ngrid)%patch_rough(1,1,1) &
          ,leaf_g(ngrid)%tstar      (1,1,1)   ,leaf_g(ngrid)%ustar      (1,1,1) &
          ,leaf_g(ngrid)%patch_area (1,1,1)   ,turb_g(ngrid)%sflux_u    (1,1)   &
          ,turb_g(ngrid)%sflux_v    (1,1)     ,turb_g(ngrid)%sflux_t    (1,1)   &
          ,grid_g(ngrid)%flpu       (1,1)     ,grid_g(ngrid)%flpv       (1,1)   &
          ,grid_g(ngrid)%flpw       (1,1)     ,turb_g(ngrid)%kpbl       (1,1)   &
          ,turb_g(ngrid)%pblhgt     (1,1)     ,turb_g(ngrid)%lmo        (1,1)   &
          ,turb_g(ngrid)%ltscale    (1,1,1)   ,turb_g(ngrid)%sigw       (1,1,1))
     end if
  endif
!MLO]

  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%scr1 (1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%flpw(1,1))
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%scr2 (1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%flpw(1,1))
  call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%vt3dh(1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%flpw(1,1))


  ! CATT
  if (CATT == 1) then
     call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%scr3(1) &
          ,basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%flpw(1,1))

     !srf----

  endif

  !_STC ....... boundary conditions even on Ke diffusion coefficient
  if(idiffk(ngrid) ==  5 .or. idiffk(ngrid) == 6) &
       call klbnd(mzp,mxp,myp,ibcon,jdim  &
       ,scratch%vt3di(1),basic_g(ngrid)%dn0(1,1,1),grid_g(ngrid)%flpw(1,1))

  !bob  swap new hkm, vkm, and vkh with past time level:  lagged K's have
  !bob  internal lateral boundary values from neighboring nodes

  ind = 0
  do j = 1,mmyp(ngrid)
     do i = 1,mmxp(ngrid)
        do k = 1,mmzp(ngrid)
           ind = ind + 1
           s1 = scratch%scr2(ind)
           s2 = scratch%scr1(ind)
           s3 = scratch%vt3dh(ind)
           scratch%scr2(ind) = turb_g(ngrid)%hkm(k,i,j)
           scratch%scr1(ind) = turb_g(ngrid)%vkm(k,i,j)
           scratch%vt3dh(ind) = turb_g(ngrid)%vkh(k,i,j)
           !! also for vt3di = K(tke) ?????    22 March 02
           !!         scratch%vt3di(ind) = turb_g(ngrid)%vke(k,i,j)
           turb_g(ngrid)%hkm(k,i,j) = s1
           turb_g(ngrid)%vkm(k,i,j) = s2
           turb_g(ngrid)%vkh(k,i,j) = s3
        enddo
     enddo
  enddo

  if (if_adap == 0) then

     call diffvel(mzp,mxp,myp,ia,iz,ja,jz,jdim,ia_1,ja_1             &
          ,ia1,ja1,iz_1,jz_1,iz1,jz1,izu,jzv,idiffk(ngrid)             &
          ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
          ,basic_g(ngrid)%wp    (1,1,1) ,tend%ut              (1)      &
          ,tend%vt              (1)     ,tend%wt              (1)      &
          ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
          ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
          ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
          ,scratch%vt3dg        (1)     ,scratch%vt3dj        (1)      &
          ,scratch%vt3dk        (1)     ,scratch%vt3dl        (1)      &
          ,scratch%vt3dm        (1)     ,scratch%vt3dn        (1)      &
          ,scratch%vt3do        (1)     ,grid_g(ngrid)%rtgu   (1,1)    &
          ,grid_g(ngrid)%rtgv   (1,1)   ,grid_g(ngrid)%rtgt   (1,1)    &
          ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
          ,turb_g(ngrid)%sflux_w(1,1)   ,basic_g(ngrid)%dn0   (1,1,1)  &
          ,basic_g(ngrid)%dn0u  (1,1,1) ,basic_g(ngrid)%dn0v  (1,1,1)  &
          ,scratch%scr1         (1)     ,scratch%scr2         (1),ibcon,mynum)

  else

     call diffvel_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim                  &
          ,iz1,jz1,izu,jzv,idiffk(ngrid)                               &
          ,basic_g(ngrid)%up    (1,1,1) ,basic_g(ngrid)%vp    (1,1,1)  &
          ,basic_g(ngrid)%wp    (1,1,1) ,tend%ut              (1)      &
          ,tend%vt              (1)     ,tend%wt              (1)      &
          ,scratch%vt3da        (1)     ,scratch%vt3db        (1)      &
          ,scratch%vt3dc        (1)     ,scratch%vt3dd        (1)      &
          ,scratch%vt3de        (1)     ,scratch%vt3df        (1)      &
          ,scratch%vt3dg        (1)     ,scratch%vt3dj        (1)      &
          ,scratch%vt3dk        (1)     ,scratch%vt3dl        (1)      &
          ,scratch%vt3dm        (1)     ,scratch%vt3dn        (1)      &
          ,scratch%vt3do        (1)     ,grid_g(ngrid)%aru    (1,1,1)  &
          ,grid_g(ngrid)%arv    (1,1,1) ,grid_g(ngrid)%arw    (1,1,1)  &
          ,grid_g(ngrid)%volu   (1,1,1) ,grid_g(ngrid)%volv   (1,1,1)  &
          ,grid_g(ngrid)%volw   (1,1,1) ,grid_g(ngrid)%flpu    (1,1)    &
          ,grid_g(ngrid)%flpv    (1,1)   ,grid_g(ngrid)%flpw    (1,1)    &
          ,turb_g(ngrid)%sflux_u(1,1)   ,turb_g(ngrid)%sflux_v(1,1)    &
          ,turb_g(ngrid)%sflux_w(1,1)   ,basic_g(ngrid)%dn0   (1,1,1)  &
          ,basic_g(ngrid)%dn0u  (1,1,1) ,basic_g(ngrid)%dn0v  (1,1,1)  &
          ,scratch%scr1         (1)     ,scratch%scr2         (1)      &
          ,grid_g(ngrid)%topma  (1,1)   ,ibcon,mynum)

  endif


  ! Convert momentum K's to scalar K's, if necessary

  if (idiffk(ngrid) <= 3 .or. idiffk(ngrid) == 7) then
     do ind = 1,mxyzp
        scratch%scr2(ind) = scratch%scr2(ind) * xkhkm(ngrid)
     enddo
  elseif (idiffk(ngrid) == 4) then
     do ind = 1,mxyzp
        scratch%vt3di(ind) = 2. * scratch%scr1(ind)
     enddo
  endif


  !- CATT
  if (CATT == 1) then
     ind = 0
     do j = 1,mmyp(ngrid)
        do i = 1,mmxp(ngrid)
           do k = 1,mmzp(ngrid)
              ind = ind + 1

              s1 = scratch%scr3(ind)
              scratch%scr3(ind)  = turb_s(ngrid)%hksc(k,i,j)
              turb_s(ngrid)%hksc(k,i,j) = s1

              ! salva o atual coef para o proximo passo no tempo.
           enddo
        enddo
     enddo

     if (idiffk(ngrid) <= 3 .or. idiffk(ngrid) == 7) then
        ind = 0
        do j = 1,mmyp(ngrid)
           do i = 1,mmxp(ngrid)
              do k = 1,mmzp(ngrid)
                 ind = ind + 1
                 scratch%scr3(ind)  =  scratch%scr3(ind)  * xkhkm(ngrid)
              enddo
           enddo
        enddo

     endif
  endif

  do n = 1,num_scalar(ngrid)

     call azero(mxp*myp,scratch%vt2da(:))
     if (nstbot == 1) then
        if (scalar_tab(n,ngrid)%name == 'THP') then
           call atob(mxp*myp,turb_g(ngrid)%sflux_t(1,1),scratch%vt2da(1))
        elseif (scalar_tab(n,ngrid)%name == 'RTP') then
           call atob(mxp*myp,turb_g(ngrid)%sflux_r(1,1),scratch%vt2da(1))
        endif
     endif

     ! 3/10/01 - Define ksf below, the "K scalar flag", to let subroutine diffsclr
     ! know which vertical K is being passed to it.  If diffsclr sees that it's
     ! a different K from the previous one, diffsclr will re-compute the tridiff
     ! matrix coefficients.  In order to use vertical scalar K's other than
     ! vt3dh and vt3di, use ksf = 3, ksf = 4, etc. for each different K.

     !_STC..................................................
     !_STC Corrections to account for the new idiffk options
     !_STC for E-l and E-eps closure. Isotropy hypothesis.
     !_STC (S. Trini Castelli)
     !_STC..................................................

     if (scalar_tab(n,ngrid)%name == 'TKEP') then
        vkh_p => scratch%vt3di
        hkh_p => scratch%scr2
        if (idiffk(ngrid) >= 4 .and. idiffk(ngrid) /= 7) hkh_p => scratch%vt3di
        ksf = 1
     elseif (scalar_tab(n,ngrid)%name == 'EPSP') then
        vkh_p => scratch%vt3di
        hkh_p => scratch%scr2
        if (idiffk(ngrid) >= 4 .and. idiffk(ngrid) /= 7)  hkh_p => scratch%vt3di
        ksf = 3
        ! Convert Ktke to Keps; it will be converted back after use below
        call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_EPS/ALF_TKE)
        call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_EPS/ALF_TKE)
     else
        vkh_p => scratch%vt3dh
        hkh_p => scratch%scr2
        if (idiffk(ngrid) >= 4 .and. idiffk(ngrid) /= 7) hkh_p => scratch%vt3dh
        ksf = 2
     endif


     ! CATT
     if (CATT == 1) then
        !srf----------------- Hor. Diffusion Coef for tracers
        if(n > (num_scalar(ngrid) - NADDSC)) then ! if(n .gt. nscl-NADDSC)
           if (idiffk(ngrid) /= 4) then
              hkh_p => scratch%scr3  !hkh_pscr3 => scratch%vt3dp(1) !ivt3dp
           endif
        endif
        !srf--------------------
     endif

     if (if_adap == 0) then


        call diffsclr(mzp,mxp,myp,ia,iz,ja,jz,jdim                    &
             ,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,n,ksf               &
             ,scalar_tab(n,ngrid)%var_p  ,scalar_tab(n,ngrid)%var_t   &
             ,scratch%vt3da(1)                                        &
             ,scratch%vt3db      (1)     ,scratch%vt3df      (1)      &
             ,scratch%vt3dg      (1)     ,scratch%vt3dj      (1)      &
             ,scratch%vt3dk      (1)     ,scratch%vt3do      (1)      &
             ,scratch%vt3dc      (1)     ,scratch%vt3dd      (1)      &
             ,scratch%vt3dl      (1)     ,scratch%vt3dm      (1)      &
             ,scratch%vt2db      (1)     ,grid_g(ngrid)%rtgt (1,1)    &
             ,scratch%vt2da      (1)     ,basic_g(ngrid)%dn0 (1,1,1)  &
             ,vkh_p                      ,hkh_p                       ) !&
        ! large and subgrid scale for GRELL CUPAR (CATT) not used
        !srf - large and subgrid scale for GRELL CUPAR (CATT)
        !!,lsfcupar_p,n)
     else

        call diffsclr_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim,n,ksf        &
             ,grid_g(ngrid)%flpw(1,1)     ,scalar_tab(n,ngrid)%var_p  &
             ,scalar_tab(n,ngrid)%var_t  ,scratch%vt3da     (1)      &
             ,scratch%vt3dc      (1)     ,scratch%vt3df     (1)      &
             ,scratch%vt3dg      (1)     ,scratch%vt3dj     (1)      &
             ,scratch%vt3dk      (1)     ,scratch%vt3dl     (1)      &
             ,scratch%vt3dm      (1)     ,scratch%vt3do     (1)      &
             ,scratch%vt2da      (1)     ,scratch%vt2db     (1)      &
             ,basic_g(ngrid)%dn0 (1,1,1) ,vkh_p                      &
             ,hkh_p                      ,grid_g(ngrid)%aru (1,1,1)  &
             ,grid_g(ngrid)%arv  (1,1,1) ,grid_g(ngrid)%arw (1,1,1)  &
             ,grid_g(ngrid)%volt (1,1,1) ,scratch%vt3db     (1)      &
             ,grid_g(ngrid)%dxu  (1,1)   ,grid_g(ngrid)%dyv (1,1)    &
             ,grid_g(ngrid)%topma(1,1)                               )
        ! CATT not ready for Shaved-ETA
        !srf - large and subgrid scale for GRELL CUPAR (CATT)
        !!,lsfcupar_p,n)
     endif

     if (scalar_tab(n,ngrid)%name == 'EPSP') then
        call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_TKE/ALF_EPS)
        call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_TKE/ALF_EPS)
     endif

  enddo


  !----------------------------------------------------------------------------------------!
  !    Integrating average turbulence parameters for mass flux.                            !
  !----------------------------------------------------------------------------------------!
  if (imassflx == 1) then
     select case (idiffk(ngrid))
     case (1,4,5,6)
        call prepare_tke_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt              &
             ,turb_g(ngrid)%tkep    (1,1,1) ,mass_g(ngrid)%tkepb    (1,1,1))
     case (7)
        call prepare_tke_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt              &
             ,turb_g(ngrid)%tkep    (1,1,1) ,mass_g(ngrid)%tkepb    (1,1,1))
        call prepare_turb_to_mass(mzp,mxp,myp,ia,iz,ja,jz,dtlt             &
             ,turb_g(ngrid)%sigw    (1,1,1) ,turb_g(ngrid)%ltscale  (1,1,1)&
             ,mass_g(ngrid)%sigwb   (1,1,1) ,mass_g(ngrid)%ltscaleb (1,1,1))
     end select
  end if

  return
end subroutine diffuse
