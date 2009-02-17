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
! This subroutine is the same as diffuse, except that it has some machine-specific         !
! optimisations.                                                                           !
!------------------------------------------------------------------------------------------!
subroutine diffuse_brams31()
   use mem_tend       , only : tend         ! ! intent(inout)
   use mem_basic      , only : basic_g      ! ! intent(in)
   use var_tables     , only : num_scalar   & ! intent(in)
                             , scalar_tab   ! ! intent(inout)
   use mem_turb       , only : idiffk       & ! intent(in)
                             , ibruvais     & ! intent(in)
                             , turb_g       & ! intent(inout)
                             , xkhkm        ! ! intent(in)
   use mem_grid       , only : jdim         & ! intent(in)
                             , ngrid        & ! intent(in)
                             , grid_g       & ! intent(in)
                             , zt           & ! intent(in)
                             , npatch       & ! intent(in)
                             , nstbot       & ! intent(in)
                             , dtlt         ! ! intent(in)
   use mem_leaf       , only : leaf_g       ! ! intent(in)
   use mem_micro      , only : micro_g      ! ! intent(in)
   use mem_scratch    , only : scratch      & ! intent(out)
                             , vctr1        & ! intent(out)
                             , vctr2        & ! intent(out)
                             , vctr3        & ! intent(out)
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
   use mem_mass       , only : imassflx     ! ! intent(in)
   use therm_lib      , only : vapour_on    ! ! intent(in)
   use mem_opt        , only : jstep        & ! intent(in)
                             , istep        & ! intent(in)
                             , opt          ! ! intent(out)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: mxyzp,ind,n,i,j,k,ksf
   integer                             :: htint_i, htint_j, iia, jja, iiz, jjz
   integer                             :: iistep, jjstep
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
   call strain(mzp,mxp,myp,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jdim                              &
             ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        ,basic_g(ngrid)%wp        &
             ,scratch%vt3da            ,scratch%vt3db            ,scratch%vt3dc            &
             ,scratch%vt3dd            ,scratch%vt3de            ,scratch%vt3df            &
             ,scratch%vt3dg            ,scratch%vt3dh            ,scratch%vt3di            &
             ,scratch%vt3dn            ,scratch%scr2             ,idiffk(ngrid)            )


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
             ,grid_g(ngrid)%dyt        ,grid_g(ngrid)%flpw       ,mynum                    )
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
     call nakanishi(mzp, mxp, myp, npatch, ia, iz, ja, jz, jdim                            &
             ,turb_g(ngrid)%tkep       ,tend%tket                ,scratch%vt3dd            &
             ,scratch%vt3de            ,scratch%vt3dh            ,scratch%vt3di            &
             ,scratch%vt3dj            ,scratch%scr1             ,grid_g(ngrid)%rtgt       &
             ,basic_g(ngrid)%theta     ,scratch%vt3dp            ,scratch%vt3dq            &
             ,basic_g(ngrid)%dn0       ,basic_g(ngrid)%up        ,basic_g(ngrid)%vp        &
             ,leaf_g(ngrid)%veg_rough  ,leaf_g(ngrid)%patch_rough,leaf_g(ngrid)%tstar      &
             ,leaf_g(ngrid)%ustar      ,leaf_g(ngrid)%patch_area ,turb_g(ngrid)%sflux_u    &
             ,turb_g(ngrid)%sflux_v    ,turb_g(ngrid)%sflux_t    ,grid_g(ngrid)%flpu       &
             ,grid_g(ngrid)%flpv       ,grid_g(ngrid)%flpw       ,turb_g(ngrid)%kpbl       &
             ,turb_g(ngrid)%pblhgt     ,turb_g(ngrid)%lmo        ,turb_g(ngrid)%ltscale    &
             ,turb_g(ngrid)%sigw                                                           )
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    We now compute the boundary condition for the various mixing coefficients.         !
   !---------------------------------------------------------------------------------------!
   call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%scr1 ,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%scr2 ,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
   call klbnd(mzp,mxp,myp,ibcon,jdim,scratch%vt3dh,basic_g(ngrid)%dn0,grid_g(ngrid)%flpw)
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
   call diffvel(mzp,mxp,myp,ia,iz,ja,jz,jdim,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,izu       &
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
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Calculating Index and Weights for interpolation (htint).  Setting indexes and     !
   ! weights for HTINT with Optmization.                                                   !
   !---------------------------------------------------------------------------------------!

   !----- Loop for a block of spatial local region. ---------------------------------------!
   jjstep = min(jstep,(jz-ja))
   iistep = min(istep,(iz-ia))
   do jja = ja, jz, jjstep
   
      jjz = MIN(jja+jjstep-1, jz) !----- Calculating last element.
      do iia = ia, iz, iistep
         iiz = MIN(iia+iistep-1, iz) !----- Calculating last element.

         !------ Compute the Indexes and Weights for x_dir - ALF --------------------------!
         do j=jja,jjz
            do i=iia,iiz
               do k=1,mzp
                  vctr1(k)=grid_g(ngrid)%topt(i-1,j) + zt(k)*grid_g(ngrid)%rtgt(i-1,j)
                  vctr2(k)=grid_g(ngrid)%topt(i  ,j) + zt(k)*grid_g(ngrid)%rtgt(i  ,j)
                  vctr3(k)=grid_g(ngrid)%topt(i+1,j) + zt(k)*grid_g(ngrid)%rtgt(i+1,j)
               end do

               htint_i = i-iia+1
               htint_j = j-jja+1

               call htint_index(mzp, vctr1, mzp, vctr2                                     &
                               , opt%ind1_x_a  (:,htint_i,htint_j)                         &
                               , opt%ind2_x_a  (:,htint_i,htint_j)                         &
                               , opt%weight_x_a(:,htint_i,htint_j)                         )

               call htint_index(mzp, vctr3, mzp, vctr2                                     &
                               , opt%ind1_x_b  (:,htint_i,htint_j)                         &
                               , opt%ind2_x_b  (:,htint_i,htint_j)                         &
                               , opt%weight_x_b(:,htint_i,htint_j)                         )

               !----- Compute the Indexes and Weights for y_dir - ALF ---------------------!
               do k=1,mzp !m1
                  vctr1(k)=grid_g(ngrid)%topt(i,j-jdim)+zt(k)*grid_g(ngrid)%rtgt(i,j-jdim)
                  vctr3(k)=grid_g(ngrid)%topt(i,j+jdim)+zt(k)*grid_g(ngrid)%rtgt(i,j+jdim)
               end do

               call htint_index(mzp, vctr1, mzp, vctr2                                     &
                               , opt%ind1_y_a  (:,htint_i,htint_j)                         &
                               , opt%ind2_y_a  (:,htint_i,htint_j)                         &
                               , opt%weight_y_a(:,htint_i,htint_j)                         )

               call htint_index(mzp, vctr3, mzp, vctr2                                     &
                               , opt%ind1_y_b  (:,htint_i,htint_j)                         &
                               , opt%ind2_y_b  (:,htint_i,htint_j)                         &
                               , opt%weight_y_b(:,htint_i,htint_j)                         )
   
            end do
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !   Computing the tendency due to turbulent diffusion for all scalars.            !
         !---------------------------------------------------------------------------------!
         do n = 1,num_scalar(ngrid)

            call azero(mxp*myp,scratch%vt2da)
            if (nstbot == 1) then
               select case (trim(scalar_tab(n,ngrid)%name))
               case ('THP')
                  call atob(mxp*myp,turb_g(ngrid)%sflux_t,scratch%vt2da)
               case ('RTP')
                  call atob(mxp*myp,turb_g(ngrid)%sflux_r,scratch%vt2da)
               end select
            end if

            !------------------------------------------------------------------------------!
            !     Define ksf below, the "K scalar flag", to let subroutine diffsclr know   !
            ! which vertical K is being passed to it.  If diffsclr sees that it's a        !
            ! different K from the previous one, diffsclr will re-compute the tridiff      !
            ! matrix coefficients. In order to use vertical scalar K's other than vt3dh    !
            ! and vt3di, use ksf = 3, ksf = 4, etc. for each different K.                  !
            !------------------------------------------------------------------------------!
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

            call diffsclr_brams31(mzp,mxp,myp,iia,iiz,jja,jjz,jdim                         &
                 ,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,n,ksf                                &
             ,scalar_tab(n,ngrid)%var_p,scalar_tab(n,ngrid)%var_t,scratch%vt3da            &
             ,scratch%vt3db            ,scratch%vt3df            ,scratch%vt3dg            &
             ,scratch%vt3dj            ,scratch%vt3dk            ,scratch%vt3do            &
             ,scratch%vt3dc            ,scratch%vt3dd            ,scratch%vt3dl            &
             ,scratch%vt3dm            ,scratch%vt2db            ,grid_g(ngrid)%rtgt       &
             ,scratch%vt2da            ,basic_g(ngrid)%dn0       ,vkh_p                    &
             ,hkh_p                                                                        )
   
            if (scalar_tab(n,ngrid)%name == 'EPSP') then
               call ae1t0 (mxyzp,vkh_p,vkh_p,ALF_TKE/ALF_EPS)
               call ae1t0 (mxyzp,hkh_p,hkh_p,ALF_TKE/ALF_EPS)
            end if
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Integrating average turbulence parameters for mass flux.                           !
   !---------------------------------------------------------------------------------------!
   call prepare_timeavg_driver(mzp,mxp,myp,ia,iz,ja,jz,dtlt,ngrid,idiffk(ngrid))
   !---------------------------------------------------------------------------------------!


   return
end subroutine diffuse_brams31
!==========================================================================================!
!==========================================================================================!
