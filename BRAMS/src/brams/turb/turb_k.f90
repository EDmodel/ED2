!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This routine is the subdriver to compute tendencies due to subgrid-scale turbulence. !
!------------------------------------------------------------------------------------------!
subroutine diffuse()
   use mem_tend       , only : tend         ! ! intent(inout)
   use mem_basic      , only : basic_g      ! ! intent(in)
   use var_tables     , only : num_scalar   & ! intent(in)
                             , scalar_tab   ! ! intent(inout)
   use mem_turb       , only : idiffk       & ! intent(in)
                             , ibruvais     & ! intent(in)
                             , turb_g       & ! intent(inout)
                             , xkhkm        ! ! intent(in)
   use mem_grid       , only : if_adap      & ! intent(in)
                             , jdim         & ! intent(in)
                             , ngrid        & ! intent(in)
                             , grid_g       & ! intent(in)
                             , dzm          & ! intent(in)
                             , dzt          & ! intent(in)
                             , npatch       & ! intent(in)
                             , nstbot       & ! intent(in)
                             , nscl         & ! intent(in)
                             , naddsc       & ! intent(in)
                             , dtlt         ! ! intent(in)
   use mem_leaf       , only : leaf_g       ! ! intent(in)
   use mem_micro      , only : micro_g      ! ! intent(in)
   use mem_scratch    , only : scratch      & ! intent(out)
                             , vctr34       ! ! intent(out)
   use node_mod       , only : mxp          & ! intent(in)
                             , myp          & ! intent(in)
                             , mzp          & ! intent(in)
                             , ia           & ! intent(in)
                             , iz           & ! intent(in)
                             , ja           & ! intent(in)
                             , jz           & ! intent(in)
                             , ia_1         & ! intent(in)
                             , ja_1         & ! intent(in)
                             , iz1          & ! intent(in)
                             , jz1          & ! intent(in)
                             , ibcon        & ! intent(in)
                             , mynum        & ! intent(in)
                             , mi0          & ! intent(in)
                             , mj0          & ! intent(in)
                             , mmyp         & ! intent(in)
                             , mmxp         & ! intent(in)
                             , mmzp         & ! intent(in)
                             , izu          & ! intent(in)
                             , jzv          & ! intent(in)
                             , ia1          & ! intent(in)
                             , ja1          & ! intent(in)
                             , iz_1         & ! intent(in)
                             , jz_1         ! ! intent(in)
   use ke_coms        , only : alf_eps      & ! intent(in)
                             , alf_tke      ! ! intent(in)
   use mem_turb_scalar, only : turb_s       ! ! intent(inout)
   use catt_start     , only : catt         ! ! intent(in)
   use mem_mass       , only : imassflx     ! ! intent(in)
   use therm_lib      , only : vapour_on    ! ! intent(in)

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: mxyzp,ind,n,i,j,k,ksf
   logical                             :: stc_closure,lesprog_closure
   real                                :: s1,s2,s3,s4
   real    , dimension(:) , pointer    :: vkh_p,hkh_p
   !---------------------------------------------------------------------------------------!

   !----- Defining # of grid points for this grid. ----------------------------------------!
   mxyzp = mxp * myp * mzp

   !----- Checking whether we are calling one of the Trini-Castelli -----------------------!
   stc_closure = idiffk(ngrid) == 5 .or. idiffk(ngrid) == 6

   !----- Checking whether we are using one of LES closures with prognostic TKE -----------!
   lesprog_closure = stc_closure .or. idiffk(ngrid) == 4


   !----- Nullifying and allocating pointers for the coefficients -------------------------!
   nullify(vkh_p)
   nullify(hkh_p)

   !---------------------------------------------------------------------------------------!
   !    Coping the vapour and total mixing ratio to scratch vectors. This way the Brunt-   !
   ! -Väisälä frequency and Nakanishi-Niino subroutines will also work for runs that       !
   !  don't use water.                                                                     !
   !  - vt3dp -> water vapour mixing ratio;                                                !
   !  - vt3dq -> total water substance mixing ratio.                                       !
   !---------------------------------------------------------------------------------------!
   call azero2(mxyzp,scratch%vt3dp,scratch%vt3dq)
   if (vapour_on) then
     call atob(mxyzp,basic_g(ngrid)%rv,scratch%vt3dp)
     call atob(mxyzp,basic_g(ngrid)%rtp,scratch%vt3dq)
   end if
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Finding the wind gradients and divergence.                                        !
   !---------------------------------------------------------------------------------------!
   if (if_adap == 0) then
      call strain(mzp,mxp,myp,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jdim                           &
             ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        &
             ,scratch%vt3da            ,scratch%vt3db            ,scratch%vt3dc            &
             ,scratch%vt3dd            ,scratch%vt3de            ,scratch%vt3df            &
             ,scratch%vt3dg            ,scratch%vt3dh            ,scratch%vt3di            &
             ,scratch%vt3dn            ,scratch%scr2             ,idiffk(ngrid)            )
   else
      call strain_adap(mzp,mxp,myp,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jdim                      &
             ,grid_g(ngrid)%flpu       ,grid_g(ngrid)%flpv       ,grid_g(ngrid)%flpw       &
             ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        &
             ,scratch%vt3da            ,scratch%vt3db            ,scratch%vt3dc            &
             ,scratch%vt3dd            ,scratch%vt3de            ,scratch%vt3df            &
             ,scratch%vt3dg            ,scratch%vt3dh            ,scratch%vt3di            &
             ,scratch%vt3dn            ,scratch%scr2             ,idiffk(ngrid)            &
             ,grid_g(ngrid)%dxm        ,grid_g(ngrid)%dxt        ,grid_g(ngrid)%dxu        &
             ,grid_g(ngrid)%dxv        ,grid_g(ngrid)%dym        ,grid_g(ngrid)%dyt        &
             ,grid_g(ngrid)%dyu        ,grid_g(ngrid)%dyv        ,dzm ,dzt                 )
   end if
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Finding the Brunt-Väisälä frequency.                                              !
   !---------------------------------------------------------------------------------------!
   call bruvais(ibruvais,mzp,mxp,myp,ia,iz,ja,jz                                           &
             ,basic_g(ngrid)%pi0       ,basic_g(ngrid)%pp        ,basic_g(ngrid)%theta     &
             ,scratch%vt3dq            ,scratch%vt3dp            ,grid_g(ngrid)%rtgt       &
             ,grid_g(ngrid)%flpw       ,scratch%vt3dj                                      )

   !---------------------------------------------------------------------------------------!
   !    Finding the horizontal diffusivity coefficients for heat and momentum. Also, for   !
   ! idiffk=2 or idiffk=3, this will compute the vertical ones.                            !
   !---------------------------------------------------------------------------------------!
   if (.not. lesprog_closure) then
      call mxdefm(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim                                       &
             ,scratch%vt3dh            ,scratch%vt3di            ,scratch%vt3dj            &
             ,scratch%vt3dk            ,scratch%scr1             ,scratch%scr2             &
             ,basic_g(ngrid)%dn0       ,grid_g(ngrid)%rtgt       ,grid_g(ngrid)%dxt        &
             ,grid_g(ngrid)%dyt        ,grid_g(ngrid)%flpw       ,turb_g(ngrid)%akscal     &
             ,mynum                    )

      !------------------------------------------------------------------------------------!
      !    If CATT is on, we need to find specific horizontal diffusion coefficients for   !
      ! tracers.                                                                           !
      !------------------------------------------------------------------------------------!
      if (CATT == 1) then
         call mxdefm_tracer(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim                             &
             ,scratch%vt3dh            ,scratch%scr3             ,basic_g(ngrid)%dn0       &
             ,grid_g(ngrid)%dxt        ,grid_g(ngrid)%dyt        ,grid_g(ngrid)%flpw       &
             ,turb_g(ngrid)%akscal     ,mynum                                              )
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Calling the unique subroutine to find the vertical diffusion coefficients of heat,   !
   ! momentum and energy. This will not be called for idiffk=2 or idiffk=3 because the     !
   ! vertical coefficient has already been found. In addition, these schemes will also     !
   ! determine the TKE tendency.                                                           !
   !---------------------------------------------------------------------------------------!
   select case (idiffk(ngrid))
   case (1) !----- Original Mellor-Yamada Level 2½: Helfand and Labraga (1988). -----------!
      call tkemy(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim,mi0(ngrid),mj0(ngrid)                  &
             ,turb_g(ngrid)%tkep       ,tend%tket                ,scratch%vt3dh            &
             ,scratch%vt3di            ,scratch%vt3dj            ,scratch%scr1             &
             ,grid_g(ngrid)%rtgt       ,basic_g(ngrid)%theta     ,basic_g(ngrid)%dn0       &
             ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        &
             ,turb_g(ngrid)%sflux_u    ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_w    &
             ,turb_g(ngrid)%sflux_t    ,vctr34                   ,grid_g(ngrid)%flpw       &
             ,grid_g(ngrid)%flpu       ,grid_g(ngrid)%flpv       ,turb_g(ngrid)%sigw       )

   case (4) !----- Deardoff (1980) LES scheme. --------------------------------------------!
      call mxtked(mzp,mxp,myp,ia,iz,ja,jz,ibcon,jdim                                       &
             ,turb_g(ngrid)%tkep       ,tend%tket                ,basic_g(ngrid)%up        &
             ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        ,basic_g(ngrid)%rtp       &
             ,basic_g(ngrid)%rv        ,basic_g(ngrid)%theta     ,scratch%vt3da            &
             ,scratch%vt3dc            ,scratch%vt3dh            ,scratch%vt3dj            &
             ,scratch%scr1             ,scratch%scr2             ,turb_g(ngrid)%sflux_u    &
             ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_w    ,turb_g(ngrid)%sflux_t    &
             ,grid_g(ngrid)%dxt        ,grid_g(ngrid)%rtgt       ,grid_g(ngrid)%flpw       )

   case (5) !----- Trini Castelli tkescl for E-l closure. ---------------------------------!
      call tkescl(mzp,mxp,myp,npatch,ia,iz,ja,jz                                           &
             ,turb_g(ngrid)%tkep       ,tend%tket                ,turb_g(ngrid)%epsp       &
             ,scratch%vt3da            ,scratch%vt3dc            ,scratch%vt3dh            &
             ,scratch%vt3di            ,scratch%vt3dj            ,scratch%scr1             &
             ,scratch%scr2             ,grid_g(ngrid)%rtgt       ,scratch%vt3dd            &
             ,scratch%vt3de            ,grid_g(ngrid)%dxt        ,leaf_g(ngrid)%ustar      &
             ,leaf_g(ngrid)%patch_area ,grid_g(ngrid)%flpw       ,basic_g(ngrid)%dn0       )

   case (6) !----- Trini Castelli subroutine TKE-eps for E-eps closure. -------------------!
      call tkeeps(mzp,mxp,myp,npatch,ia,iz,ja,jz                                           &
             ,turb_g(ngrid)%tkep       ,tend%tket                ,turb_g(ngrid)%epsp       &
             ,tend%epst                ,scratch%vt3da            ,scratch%vt3dc            &
             ,scratch%vt3dh            ,scratch%vt3di            ,scratch%vt3dj            &
             ,scratch%scr1             ,scratch%scr2             ,grid_g(ngrid)%rtgt       &
             ,leaf_g(ngrid)%ustar      ,leaf_g(ngrid)%patch_area ,grid_g(ngrid)%flpw       &
             ,basic_g(ngrid)%dn0                                                           )

   case (7) !----- Nananishi and Niino (2004), scheme based on Mellor-Yamada Level 2½. ----!
      call nakanishi(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim                           &
              ,turb_g(ngrid)%tkep       ,tend%tket                ,scratch%vt3dd           &
              ,scratch%vt3de            ,scratch%vt3dh            ,scratch%vt3di           &
              ,scratch%vt3dj            ,scratch%scr1             ,grid_g(ngrid)%rtgt      &
              ,basic_g(ngrid)%theta     ,scratch%vt3dp            ,scratch%vt3dq           &
              ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp       &
              ,leaf_g(ngrid)%veg_rough  ,leaf_g(ngrid)%patch_rough,leaf_g(ngrid)%tstar     &
              ,leaf_g(ngrid)%ustar      ,leaf_g(ngrid)%patch_area ,turb_g(ngrid)%sflux_u   &
              ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_t    ,grid_g(ngrid)%flpu      &
              ,grid_g(ngrid)%flpv       ,grid_g(ngrid)%flpw       ,turb_g(ngrid)%kpbl      &
              ,turb_g(ngrid)%pblhgt     ,turb_g(ngrid)%lmo        ,turb_g(ngrid)%ltscale   &
              ,turb_g(ngrid)%sigw                                                          )

   case (8) !----- Nananishi and Niino (2004), scheme based on Mellor-Yamada Level 2½. ----!
      call nakanishi(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim                           &
              ,turb_g(ngrid)%tkep       ,tend%tket                ,scratch%vt3dd           &
              ,scratch%vt3de            ,scratch%vt3dh            ,scratch%vt3di           &
              ,scratch%vt3dj            ,scratch%scr1             ,grid_g(ngrid)%rtgt      &
              ,basic_g(ngrid)%theta     ,scratch%vt3dp            ,scratch%vt3dq           &
              ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp       &
              ,leaf_g(ngrid)%veg_rough  ,leaf_g(ngrid)%patch_rough,leaf_g(ngrid)%tstar     &
              ,leaf_g(ngrid)%ustar      ,leaf_g(ngrid)%patch_area ,turb_g(ngrid)%sflux_u   &
              ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_t    ,grid_g(ngrid)%flpu      &
              ,grid_g(ngrid)%flpv       ,grid_g(ngrid)%flpw       ,turb_g(ngrid)%kpbl      &
              ,turb_g(ngrid)%pblhgt     ,turb_g(ngrid)%lmo        ,turb_g(ngrid)%ltscale   &
              ,turb_g(ngrid)%sigw                                                          )
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    We now compute the boundary condition for the various mixing coefficients.         !
   !---------------------------------------------------------------------------------------!
   call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%scr1 ,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%scr2 ,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%vt3dh,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   if (CATT == 1) then
      call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%scr3,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   end if
   if (stc_closure) then
      call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%vt3di,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Bob: Swap new hkm, vkm, and vkh with past time level:  Lagged K's have internal       !
   !      lateral boundary values from neighboring nodes.                                  !
   !---------------------------------------------------------------------------------------!
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
            turb_g(ngrid)%hkm(k,i,j) = s1
            turb_g(ngrid)%vkm(k,i,j) = s2
            turb_g(ngrid)%vkh(k,i,j) = s3
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Compute fluxes using the computed K's and previous strains.                        !
   !---------------------------------------------------------------------------------------!
   if (if_adap == 0) then
      call diffvel(mzp,mxp,myp,ia,iz,ja,jz,jdim,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,izu    &
             ,jzv                      ,idiffk(ngrid)            ,basic_g(ngrid)%up        &
             ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        ,tend%ut                  &
             ,tend%vt                  ,tend%wt                  ,scratch%vt3da            &
             ,scratch%vt3db            ,scratch%vt3dc            ,scratch%vt3dd            &
             ,scratch%vt3de            ,scratch%vt3df            ,scratch%vt3dg            &
             ,scratch%vt3dj            ,scratch%vt3dk            ,scratch%vt3dl            &
             ,scratch%vt3dm            ,scratch%vt3dn            ,scratch%vt3do            &
             ,grid_g(ngrid)%rtgu       ,grid_g(ngrid)%rtgv       ,grid_g(ngrid)%rtgt       &
             ,turb_g(ngrid)%sflux_u    ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_w    &
             ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%dn0u      ,basic_g(ngrid)%dn0v      &
             ,scratch%scr1             ,scratch%scr2             ,ibcon                    &
             ,mynum                                                                        )
   else
      call diffvel_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim,iz1,jz1,izu,jzv,idiffk(ngrid)         &
             ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        &
             ,tend%ut                  ,tend%vt                  ,tend%wt                  &
             ,scratch%vt3da            ,scratch%vt3db            ,scratch%vt3dc            &
             ,scratch%vt3dd            ,scratch%vt3de            ,scratch%vt3df            &
             ,scratch%vt3dg            ,scratch%vt3dj            ,scratch%vt3dk            &
             ,scratch%vt3dl            ,scratch%vt3dm            ,scratch%vt3dn            &
             ,scratch%vt3do            ,grid_g(ngrid)%aru        ,grid_g(ngrid)%arv        &
             ,grid_g(ngrid)%arw        ,grid_g(ngrid)%volu       ,grid_g(ngrid)%volv       &
             ,grid_g(ngrid)%volw       ,grid_g(ngrid)%flpu       ,grid_g(ngrid)%flpv       &
             ,grid_g(ngrid)%flpw       ,turb_g(ngrid)%sflux_u    ,turb_g(ngrid)%sflux_v    &
             ,turb_g(ngrid)%sflux_w    ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%dn0u      &
             ,basic_g(ngrid)%dn0v      ,scratch%scr1             ,scratch%scr2             &
             ,grid_g(ngrid)%topma      ,ibcon                    ,mynum                    )
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert momentum K's to scalar K's, if necessary.                                 !
   !---------------------------------------------------------------------------------------!
   if (.not. lesprog_closure) then
      do ind = 1,mxyzp
         scratch%scr2(ind) = scratch%scr2(ind) * xkhkm(ngrid)
      end do
   else if (.not. stc_closure) then
      do ind = 1,mxyzp
         scratch%vt3di(ind) = 2. * scratch%scr1(ind)
      end do
   end if
   !----- If this is a CATT run, also update the trace gas K ------------------------------!
   if (CATT == 1) then
      ind = 0
      do j = 1,mmyp(ngrid)
         do i = 1,mmxp(ngrid)
            do k = 1,mmzp(ngrid)
               ind = ind + 1
               s1 = scratch%scr3(ind)
               scratch%scr3(ind)  = turb_s(ngrid)%hksc(k,i,j)
               turb_s(ngrid)%hksc(k,i,j) = s1
            end do
         end do
      end do
      if (.not. lesprog_closure) then
         ind = 0
         do j = 1,mmyp(ngrid)
            do i = 1,mmxp(ngrid)
               do k = 1,mmzp(ngrid)
                  ind = ind + 1
                  scratch%scr3(ind)  =  scratch%scr3(ind)  * xkhkm(ngrid)
               end do
            end do
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   Computing the tendency due to turbulent diffusion for all scalars.                  !
   !---------------------------------------------------------------------------------------!
   do n = 1,num_scalar(ngrid)

      call azero(mxp*myp,scratch%vt2da)
      if (nstbot == 1) then
         select case (trim(scalar_tab(n,ngrid)%name))
         case ('THP')
            call atob(mxp*myp,turb_g(ngrid)%sflux_t,scratch%vt2da)
         case ('RTP')
            call atob(mxp*myp,turb_g(ngrid)%sflux_r,scratch%vt2da)
         case ('CO2P')
            call atob(mxp*myp,turb_g(ngrid)%sflux_c,scratch%vt2da)
         end select
      end if

      !------------------------------------------------------------------------------------!
      !     Define ksf below, the "K scalar flag", to let subroutine diffsclr know which   !
      ! vertical K is being passed to it.  If diffsclr sees that it's a different K from   !
      ! the previous one, diffsclr will re-compute the tridiff matrix coefficients. In     !
      ! order to use vertical scalar K's other than vt3dh and vt3di, use ksf = 3, ksf = 4, !
      ! etc. for each different K.                                                         !
      !------------------------------------------------------------------------------------!
      select case (trim(scalar_tab(n,ngrid)%name))
      case ('TKEP')
         vkh_p => scratch%vt3di
         hkh_p => scratch%scr2
         if (lesprog_closure) hkh_p => scratch%vt3di
         ksf = 1
      case ('EPSP')
         vkh_p => scratch%vt3di
         hkh_p => scratch%scr2
         if (lesprog_closure)  hkh_p => scratch%vt3di
         ksf = 3
         ! Convert Ktke to Keps; it will be converted back after use below
         call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_EPS/ALF_TKE)
         call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_EPS/ALF_TKE)
      case default
         vkh_p => scratch%vt3dh
         hkh_p => scratch%scr2
         if (lesprog_closure) hkh_p => scratch%vt3dh
         ksf = 2
      end select

      !----- CATT -------------------------------------------------------------------------!
      if (CATT == 1 .and. (n > num_scalar(ngrid) - naddsc) .and. (.not. lesprog_closure))  &
      then 
          hkh_p => scratch%scr3  !hkh_pscr3 => scratch%vt3dp(1) !ivt3dp
      end if

      if (if_adap == 0) then
         call diffsclr(mzp, mxp, myp, ia, iz, ja, jz, jdim, ia_1, ja_1, ia1, ja1           &
             ,iz_1, jz_1, iz1, jz1, n, ksf                                                 &
             ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t,scratch%vt3da            &
             ,scratch%vt3db            ,scratch%vt3df            ,scratch%vt3dg            &
             ,scratch%vt3dj            ,scratch%vt3dk            ,scratch%vt3do            &
             ,scratch%vt3dc            ,scratch%vt3dd            ,scratch%vt3dl            &
             ,scratch%vt3dm            ,scratch%vt2db            ,grid_g(ngrid)%rtgt       &
             ,scratch%vt2da            ,basic_g(ngrid)%dn0       ,vkh_p                    &
             ,hkh_p                    )
      else
         call diffsclr_adap(mzp,mxp,myp, ia,iz,ja,jz, jdim,n,ksf ,grid_g(ngrid)%flpw       &
             ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t,scratch%vt3da            &
             ,scratch%vt3dc            ,scratch%vt3df            ,scratch%vt3dg            &
             ,scratch%vt3dj            ,scratch%vt3dk            ,scratch%vt3dl            &
             ,scratch%vt3dm            ,scratch%vt3do            ,scratch%vt2da            &
             ,scratch%vt2db            ,basic_g(ngrid)%dn0       ,vkh_p                    &
             ,hkh_p                    ,grid_g(ngrid)%aru        ,grid_g(ngrid)%arv        &
             ,grid_g(ngrid)%arw        ,grid_g(ngrid)%volt       ,scratch%vt3db            &
             ,grid_g(ngrid)%dxu        ,grid_g(ngrid)%dyv        ,grid_g(ngrid)%topma      )
      end if

      if (scalar_tab(n,ngrid)%name == 'EPSP') then
         call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_TKE/ALF_EPS)
         call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_TKE/ALF_EPS)
      end if

   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Integrating average turbulence parameters for mass flux.                           !
   !---------------------------------------------------------------------------------------!
   call prepare_timeavg_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,ngrid,idiffk(ngrid))
   !---------------------------------------------------------------------------------------!

   return
end subroutine diffuse
!==========================================================================================!
!==========================================================================================!
