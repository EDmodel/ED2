!==========================================================================================!
!==========================================================================================!
!    This subroutine will solve the within canopy two-stream radiation for shortwave,      !
! considering the light absorption by leaves, and branches if the user wants so.  In both  !
! cases we consider the clumpiness effect (i.e., the fact that leaves and branches aren't  !
! randomly distributed).  The crown area will be also considered when the user wants it.   !
!------------------------------------------------------------------------------------------!
subroutine sw_twostream_clump(salb,scosz,scosaoi,ncoh,pft,TAI,canopy_area                  &
                             ,PAR_beam_flip,PAR_diffuse_flip,SW_abs_beam_flip              &
                             ,SW_abs_diffuse_flip,DW_vislo_beam,DW_vislo_diffuse           &
                             ,UW_vishi_beam,UW_vishi_diffuse,DW_nirlo_beam                 &
                             ,DW_nirlo_diffuse,UW_nirhi_beam,UW_nirhi_diffuse              &
                             ,beam_level,diff_level,lambda_coh,lambda_tot)

   use ed_max_dims          , only : n_pft                   ! ! intent(in) 
   use pft_coms             , only : clumping_factor         & ! intent(in) 
                                   , phenology               ! ! intent(in) 
   use canopy_radiation_coms, only : diffuse_backscatter_nir & ! intent(in)
                                   , diffuse_backscatter_vis & ! intent(in)
                                   , leaf_scatter_nir        & ! intent(in)
                                   , leaf_scatter_vis        & ! intent(in)
                                   , visible_fraction_dir    & ! intent(in)
                                   , visible_fraction_dif    & ! intent(in)
                                   , cosz_min8               ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, dimension(ncoh)     , intent(in)    :: pft
   integer                      , intent(in)    :: ncoh
   real(kind=8), dimension(ncoh), intent(in)    :: TAI
   real(kind=8), dimension(ncoh), intent(in)    :: canopy_area
   real                         , intent(in)    :: salb
   real                         , intent(in)    :: scosz
   real                         , intent(in)    :: scosaoi
   real, dimension(ncoh)        , intent(out)   :: PAR_beam_flip
   real, dimension(ncoh)        , intent(out)   :: PAR_diffuse_flip
   real, dimension(ncoh)        , intent(out)   :: SW_abs_beam_flip
   real, dimension(ncoh)        , intent(out)   :: SW_abs_diffuse_flip
   real                         , intent(out)   :: UW_vishi_beam
   real                         , intent(out)   :: UW_vishi_diffuse
   real                         , intent(out)   :: UW_nirhi_beam
   real                         , intent(out)   :: UW_nirhi_diffuse
   real                         , intent(out)   :: DW_vislo_beam
   real                         , intent(out)   :: DW_vislo_diffuse
   real                         , intent(out)   :: DW_nirlo_beam
   real                         , intent(out)   :: DW_nirlo_diffuse
   real(kind=8), dimension(ncoh), intent(out)   :: beam_level
   real(kind=8), dimension(ncoh), intent(out)   :: diff_level
   real(kind=8), dimension(ncoh), intent(out)   :: lambda_coh
   real(kind=8)                 , intent(out)   :: lambda_tot
   !----- Local variables -----------------------------------------------------------------!
   integer     , dimension(2*ncoh)        :: indx
   integer                                :: il
   integer                                :: ipft
   integer                                :: ncoh2
   integer                                :: iband
   integer                                :: i
   integer                                :: j
   integer                                :: ind
   real(kind=8), dimension(n_pft)         :: leaf_scatter
   real(kind=8), dimension(n_pft)         :: diffuse_backscatter
   real(kind=8), dimension(ncoh)          :: expkl_top
   real(kind=8), dimension(ncoh)          :: expkl_bot
   real(kind=8), dimension(ncoh)          :: expamk_top
   real(kind=8), dimension(ncoh)          :: expamk_bot
   real(kind=8), dimension(ncoh)          :: expapk_top
   real(kind=8), dimension(ncoh)          :: expapk_bot
   real(kind=8), dimension(ncoh)          :: A_top
   real(kind=8), dimension(ncoh)          :: A_bot
   real(kind=8), dimension(ncoh)          :: B_top
   real(kind=8), dimension(ncoh)          :: B_bot
   real(kind=8), dimension(ncoh)          :: C_top
   real(kind=8), dimension(ncoh)          :: C_bot
   real(kind=8), dimension(ncoh)          :: F_top
   real(kind=8), dimension(ncoh)          :: F_bot
   real(kind=8), dimension(ncoh)          :: G_top
   real(kind=8), dimension(ncoh)          :: G_bot
   real(kind=8), dimension(ncoh)          :: H_top
   real(kind=8), dimension(ncoh)          :: H_bot
   real(kind=8), dimension(ncoh)          :: beam_bot
   real(kind=8), dimension(ncoh)          :: beam_bot_crown
   real(kind=8), dimension(ncoh+1)        :: upward_vis_beam
   real(kind=8), dimension(ncoh+1)        :: upward_vis_diffuse
   real(kind=8), dimension(ncoh+1)        :: upward_nir_beam
   real(kind=8), dimension(ncoh+1)        :: upward_nir_diffuse
   real(kind=8), dimension(ncoh+1)        :: downward_nir_beam
   real(kind=8), dimension(ncoh+1)        :: downward_nir_diffuse
   real(kind=8), dimension(ncoh+1)        :: downward_vis_beam
   real(kind=8), dimension(ncoh+1)        :: downward_vis_diffuse
   real(kind=8), dimension(2*ncoh)        :: mastervec_beam
   real(kind=8), dimension(2*ncoh)        :: masveccp_beam
   real(kind=8), dimension(2*ncoh)        :: mastervec_diffuse
   real(kind=8), dimension(2*ncoh)        :: masveccp_diffuse
   real(kind=8), dimension(2*ncoh,2)      :: matal
   real(kind=8), dimension(2*ncoh,5)      :: mastermat
   real(kind=8), dimension(2*ncoh,2*ncoh) :: masmatcp
   real(kind=8)                           :: alb
   real(kind=8)                           :: cosz
   real(kind=8)                           :: cosaoi
   real(kind=8)                           :: lambda
   real(kind=8)                           :: beam_backscatter
   real(kind=8)                           :: eta
   real(kind=8)                           :: zeta
   real(kind=8)                           :: exk
   real(kind=8)                           :: exki
   real(kind=8)                           :: zetai
   real(kind=8)                           :: d
   real(kind=8)                           :: rhoo
   real(kind=8)                           :: sigma
   real(kind=8)                           :: source_bot
   real(kind=8)                           :: source_top
   real(kind=8), dimension(ncoh)          :: eff_tai
   !---------------------------------------------------------------------------------------!

   
   !----- Convert input variable to double precision. -------------------------------------!
   alb    = dble(salb)
   cosz   = max(cosz_min8,dble(scosz))
   cosaoi = max(cosz_min8,dble(scosaoi))

   !----- Calculate factors common for NIR, PAR. ------------------------------------------!
   ncoh2      = 2*ncoh
   lambda     = 5.d-1/cosaoi
   lambda_tot = 0.0d0
   do il=1,ncoh
      ipft           = pft(il)
      lambda_tot     = lambda_tot + clumping_factor(ipft)
      lambda_coh(il) = lambda * clumping_factor(ipft) / canopy_area(il)
      eff_tai(il)    = clumping_factor(ipft)*TAI(il)
   end do
   lambda_tot = lambda_tot * lambda / dble(ncoh)
   beam_backscatter = (5.d-1 + cosz) * (1.0d0 - cosz*log(1.0d0+1.0d0/cosz))
  
   !----- Loop over bands (currently Visible and near infrared). --------------------------!
   bandloop: do iband = 1,2
      select case(iband)
      case (1) !----- Visible (or PAR). ---------------------------------------------------!
         do ipft = 1,n_pft
            leaf_scatter(ipft)        = dble(leaf_scatter_vis(ipft))
            diffuse_backscatter(ipft) = dble(diffuse_backscatter_vis(ipft))
         end do
      case (2) !----- Near infrared (or NIR). ---------------------------------------------!
         do ipft = 1,n_pft
            leaf_scatter(ipft)        = dble(leaf_scatter_nir(ipft))
            diffuse_backscatter(ipft) = dble(diffuse_backscatter_nir(ipft))
         end do
      end select
      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     Calculate more factors for this band.  We start by calculating the forcings,   !
      ! using the effective tree area index, considering the leaf clumpiness and, if the   !
      ! user asked for it, the branchwood area index (in case the user doesn't want,       !
      ! WAI=0.).  We distinguish between the light underneath a crown (beam_bot_crown) and !
      ! the total light passing through a layer (beam_bot).  Also, we consider that the    !
      ! canopy may be partially open in case the user wants this distinction (otherwise    !
      ! canopy_area is one.).                                                              !
      !------------------------------------------------------------------------------------!
      !----- Start with the tallest cohort, moving downwards. -----------------------------!
      beam_bot_crown(ncoh) = exp(-lambda*eff_tai(ncoh)/canopy_area(ncoh))
      beam_level(ncoh)     = exp(-5.d-1*lambda*eff_tai(ncoh)/canopy_area(ncoh))
      beam_bot(ncoh)       = (1.d0-canopy_area(ncoh))                                      &
                           + canopy_area(ncoh)*beam_bot_crown(ncoh)
      do il=ncoh-1,1,-1
         beam_bot_crown(il) = beam_bot(il+1) * exp(-lambda*eff_tai(il)/canopy_area(il))
         beam_bot(il)       = beam_bot(il+1)*(1.d0-canopy_area(il))                        & 
                            + canopy_area(il)*beam_bot_crown(il)
         beam_level(il)     = beam_bot(il+1)                                               &
                            * exp(-5.d-1*lambda*eff_tai(il)/canopy_area(il))               &
                            * canopy_area(il)                                              &
                            + (1.d0-canopy_area(il)) * beam_bot(il+1)
      end do

      do il=1,ncoh
         ipft  = pft(il)
         eta   = clumping_factor(ipft)                                                     &
               * (1.0d0 - (1.0d0-diffuse_backscatter(ipft)) * leaf_scatter(ipft))
         zeta  = leaf_scatter(ipft) * diffuse_backscatter(ipft) * clumping_factor(ipft)
         exk   = sqrt(eta*eta - zeta*zeta)
         exki  = 1.0d0/exk
         zetai = 1.0d0/zeta

         !----- Sources. ------------------------------------------------------------------!
         source_bot = clumping_factor(ipft)*lambda*leaf_scatter(ipft)                      &
                    * beam_backscatter * beam_bot_crown(il)
         source_top = source_bot * exp(lambda*eff_tai(il))

         !----- Forcing coefficients. -----------------------------------------------------!
         rhoo  = - (zeta + eta + clumping_factor(ipft)*lambda)                             &
               * clumping_factor(ipft) * lambda                                            &
               * leaf_scatter(ipft) * beam_backscatter                                     &
               * beam_bot_crown(il)
         sigma = clumping_factor(ipft) * lambda

         !----- Calculate exponentials only once. -----------------------------------------!
         expkl_bot(il)  = 1.0d0
         expkl_top(il)  = exp(exk*TAI(il))
         expamk_bot(il) = 1.0d0
         expamk_top(il) = exp((sigma-exk)*TAI(il))
         expapk_bot(il) = 1.0d0
         expapk_top(il) = exp((sigma+exk)*TAI(il))
         A_bot(il)      = -source_bot*zetai
         A_top(il)      = -source_top*zetai                                                &
                        + 5.d-1*zetai*(eta*exki-1.0d0)*expkl_top(il)*rhoo/(sigma-exk)      &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 5.d-1*zetai*(eta*exki+1.0d0) / expkl_top(il)*rhoo/(sigma+exk)    &
                        * (expapk_top(il)-expapk_bot(il))
         B_bot(il)      = 5.d-1*zetai*(eta*exki-1.0d0)
         B_top(il)      = 5.d-1*zetai*(eta*exki-1.0d0)*expkl_top(il)
         C_bot(il)      = -5.d-1*zetai*(eta*exki+1.0d0)
         C_top(il)      = -5.d-1*zetai*(eta*exki+1.0d0)/expkl_top(il)
         F_bot(il)      = 0.0d0
         F_top(il)      = 5.d-1*exki*expkl_top(il)*rhoo/(sigma-exk)                        &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 5.d-1*exki/expkl_top(il)*rhoo/(sigma+exk)                        &
                        * (expapk_top(il)-expapk_bot(il))
         G_bot(il)      = 5.d-1*exki
         G_top(il)      = 5.d-1*exki*expkl_top(il)
         H_bot(il)      = -5.d-1*exki
         H_top(il)      = -5.d-1*exki/expkl_top(il)
      end do

      
      !----- Organize the matrix coefficients. --------------------------------------------!
      do j=1,ncoh2
         do i=1,ncoh2
            masmatcp(i,j)=0.0d0
         end do
         mastervec_beam(j)    = 0.0d0
         masveccp_beam(j)     = 0.0d0
         mastervec_diffuse(j) = 0.0d0
         masveccp_diffuse(j)  = 0.0d0
      end do

      masmatcp(1,1)        = G_top(ncoh)
      masmatcp(1,2)        = H_top(ncoh)
      mastervec_beam(1)    = -F_top(ncoh)
      mastervec_diffuse(1) = 1.0d0
      masveccp_beam(1)     = mastervec_beam(1)
      masveccp_diffuse(1)  = mastervec_diffuse(1)

      do i=2,ncoh2-2,2
         masmatcp(i,i-1) = G_bot(nint(real(ncoh2-i+2)*0.5))
         masmatcp(i,i)   = H_bot(nint(real(ncoh2-i+2)*0.5))

         masmatcp(i,i+1)      = -G_top(nint(real(ncoh2-i)*0.5))
         masmatcp(i,i+2)      = -H_top(nint(real(ncoh2-i)*0.5))
         mastervec_beam(i)    = -F_bot(nint(real(ncoh2-i+2)*0.5))                          &
                              + F_top(nint(real(ncoh2-i)*0.5))
         mastervec_diffuse(i) = 0.0d0
         masveccp_beam(i)     = mastervec_beam(i)
         masveccp_diffuse(i)  = mastervec_diffuse(i)
      end do

      do i=3,ncoh2-1,2
         masmatcp(i,i-2)      = B_bot(nint(real(ncoh2-i+3)*0.5))
         masmatcp(i,i-1)      = C_bot(nint(real(ncoh2-i+3)*0.5))
         masmatcp(i,i)        = -B_top(nint(real(ncoh2-i+1)*0.5))
         masmatcp(i,i+1)      = -C_top(nint(real(ncoh2-i+1)*0.5))
         mastervec_beam(i)    = -A_bot(nint(real(ncoh2-i+3)*0.5))                          &
                              + A_top(nint(real(ncoh2-i+1)*0.5))
         masveccp_beam(i)     = mastervec_beam(i)
         mastervec_diffuse(i) = 0.0d0
         masveccp_diffuse(i)  = mastervec_diffuse(i)
      end do
      masmatcp(ncoh2,ncoh2-1)  = B_bot(1)-alb*G_bot(1)
      masmatcp(ncoh2,ncoh2)    = C_bot(1)-alb*H_bot(1)
      mastervec_beam(ncoh2)    = -A_bot(1)+alb*beam_bot(1)
      masveccp_beam(ncoh2)     = mastervec_beam(ncoh2)
      mastervec_diffuse(ncoh2) = 0.0d0
      masveccp_diffuse(ncoh2)  = mastervec_diffuse(ncoh2)
      
      !----- Prep for inversion. ----------------------------------------------------------!
      mastermat(1,1) = 0.d0
      mastermat(1,2) = 0.d0
      mastermat(1,3) = masmatcp(1,1)
      mastermat(1,4) = masmatcp(1,2)
      mastermat(1,5) = 0.d0

      do i=2,ncoh2-2,2
         mastermat(i,1) = 0.d0
         mastermat(i,2) = masmatcp(i,i-1)
         mastermat(i,3) = masmatcp(i,i)
         mastermat(i,4) = masmatcp(i,i+1)
         mastermat(i,5) = masmatcp(i,i+2)
      end do

      do i=3,ncoh2-1,2
         mastermat(i,1) = masmatcp(i,i-2)
         mastermat(i,2) = masmatcp(i,i-1)
         mastermat(i,3) = masmatcp(i,i)
         mastermat(i,4) = masmatcp(i,i+1)
         mastermat(i,5) = 0.d0
      end do

      mastermat(ncoh2,1) = 0.d0
      mastermat(ncoh2,2) = masmatcp(ncoh2,ncoh2-1)
      mastermat(ncoh2,3) = masmatcp(ncoh2,ncoh2)
      mastermat(ncoh2,4) = 0.d0
      mastermat(ncoh2,5) = 0.d0
      
      !----- Invert the matrix. -----------------------------------------------------------!
      call bandec(mastermat,ncoh2,ncoh2,2,2,matal,indx,d)

      !----- Backsubstitute for beam and diffuse. -----------------------------------------!
      call banbks(mastermat,ncoh2,ncoh2,2,2,matal,indx,mastervec_beam)
      call banbks(mastermat,ncoh2,ncoh2,2,2,matal,indx,mastervec_diffuse)
      
      !----- Improve the solution. --------------------------------------------------------!
      call mprove(masmatcp,mastermat,matal,ncoh2,ncoh2,5,2,indx                            &
                 ,masveccp_beam,mastervec_beam)
      call mprove(masmatcp,mastermat,matal,ncoh2,ncoh2,5,2,indx                            &
                 ,masveccp_diffuse,mastervec_diffuse)

      select case (iband)
      case (1) !---- Visible (or PAR) band. -----------------------------------------------!
         do i=3,ncoh2-1,2
            ind                     = nint(real(ncoh2-i+1)*0.5)+1
            upward_vis_beam(ind)    = A_bot(ind) + B_bot(ind) * mastervec_beam(i-2)        &
                                    + C_bot(ind) * mastervec_beam(i-1)
            upward_vis_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2)                  &
                                    + C_bot(ind) * mastervec_diffuse(i-1)
         end do

         do i=2,ncoh2-2,2
            ind                       = nint(real(ncoh2-i)*0.5)+1
            downward_vis_beam(ind)    = beam_bot(ind) + F_bot(ind)                         &
                                      + H_bot(ind) * mastervec_beam(i)                     &
                                      + G_bot(ind) * mastervec_beam(i-1)
            downward_vis_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)                  &
                                      + G_bot(ind) * mastervec_diffuse(i-1)
         end do

         upward_vis_beam(ncoh+1)      = B_top(ncoh) * mastervec_beam(1)                    &
                                      + C_top(ncoh) * mastervec_beam(2) + A_top(ncoh)
         upward_vis_diffuse(ncoh+1)   = B_top(ncoh) * mastervec_diffuse(1)                 &
                                      + C_top(ncoh) * mastervec_diffuse(2)
         downward_vis_beam(ncoh+1)    = 1.0d0
         downward_vis_diffuse(ncoh+1) = 1.0d0
         downward_vis_beam(1)         = G_bot(1) * mastervec_beam(ncoh2-1)                 &
                                      + H_bot(1) * mastervec_beam(ncoh2) + F_bot(1)        &
                                      + beam_bot(1)
         downward_vis_diffuse(1)      = G_bot(1) * mastervec_diffuse(ncoh2-1)              &
                                      + H_bot(1) * mastervec_diffuse(ncoh2)
         upward_vis_beam(1)           = alb * downward_vis_beam(1)
         upward_vis_diffuse(1)        = alb * downward_vis_diffuse(1)
      case (2)

         do i=3,ncoh2-1,2
            ind                     = nint(real(ncoh2-i+1)*0.5)+1
            upward_nir_beam(ind)    = A_bot(ind) + B_bot(ind) * mastervec_beam(i-2)        &
                                    + C_bot(ind) * mastervec_beam(i-1)
            upward_nir_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2)                  &
                                    + C_bot(ind) * mastervec_diffuse(i-1)
         end do

         do i=2,ncoh2-2,2
            ind                       = nint(real(ncoh2-i)*0.5)+1
            downward_nir_beam(ind)    = beam_bot(ind) + F_bot(ind)                         &
                                      + H_bot(ind) * mastervec_beam(i)                     &
                                      + G_bot(ind) * mastervec_beam(i-1)
            downward_nir_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)                  &
                                      + G_bot(ind) * mastervec_diffuse(i-1)
         end do

         upward_nir_beam(ncoh+1)      = B_top(ncoh) * mastervec_beam(1)                    &
                                      + C_top(ncoh) * mastervec_beam(2) + A_top(ncoh)
         upward_nir_diffuse(ncoh+1)   = B_top(ncoh) * mastervec_diffuse(1)                 &
                                      + C_top(ncoh) * mastervec_diffuse(2)
         downward_nir_beam(ncoh+1)    = 1.0d0
         downward_nir_diffuse(ncoh+1) = 1.0d0
         downward_nir_beam(1)         = G_bot(1) * mastervec_beam(ncoh2-1)                 &
                                      + H_bot(1) * mastervec_beam(ncoh2) + F_bot(1)        &
                                      + beam_bot(1)
         downward_nir_diffuse(1)      = G_bot(1) * mastervec_diffuse(ncoh2-1)              &
                                      + H_bot(1) * mastervec_diffuse(ncoh2)
         upward_nir_beam(1)           = alb*downward_nir_beam(1)
         upward_nir_diffuse(1)        = alb*downward_nir_diffuse(1)
      end select
   end do bandloop
   
   do il=1,ncoh
      diff_level(il)          = downward_vis_diffuse(il+1)
      PAR_beam_flip(il)       = visible_fraction_dir                                       &
                              * sngl( downward_vis_beam(il+1)-downward_vis_beam(il)        &
                                    + upward_vis_beam(il)-upward_vis_beam(il+1))
      PAR_diffuse_flip(il)    = visible_fraction_dif                                       &
                              * sngl( downward_vis_diffuse(il+1)-downward_vis_diffuse(il)  &
                                    + upward_vis_diffuse(il)-upward_vis_diffuse(il+1))
      SW_abs_beam_flip(il)    = PAR_beam_flip(il) + (1.0 - visible_fraction_dir)           &
                              * sngl( downward_nir_beam(il+1)-downward_nir_beam(il)        &
                                    + upward_nir_beam(il)-upward_nir_beam(il+1))
      SW_abs_diffuse_flip(il) = PAR_diffuse_flip(il) + (1.0 - visible_fraction_dif)        &
                              * sngl(downward_nir_diffuse(il+1) - downward_nir_diffuse(il) &
                                    + upward_nir_diffuse(il) - upward_nir_diffuse(il+1))

      !----- Ensure that we don't get any negative radiation... ---------------------------!
      PAR_beam_flip(il)       = max(0.0,PAR_beam_flip(il)       )
      PAR_diffuse_flip(il)    = max(0.0,PAR_diffuse_flip(il)    )
      SW_abs_beam_flip(il)    = max(0.0,SW_abs_beam_flip(il)    )
      SW_abs_diffuse_flip(il) = max(0.0,SW_abs_diffuse_flip(il) )
   end do
   !---------------------------------------------------------------------------------------!





   !----- Copy to the output variables. ---------------------------------------------------!
   DW_vislo_beam    = max(0.0,sngl(downward_vis_beam(1)))       * visible_fraction_dir
   DW_vislo_diffuse = max(0.0,sngl(downward_vis_diffuse(1)))    * visible_fraction_dif
   UW_vishi_beam    = max(0.0,sngl(upward_vis_beam(ncoh+1)))    * visible_fraction_dir
   UW_vishi_diffuse = max(0.0,sngl(upward_vis_diffuse(ncoh+1))) * visible_fraction_dif
   DW_nirlo_beam    = max(0.0,sngl(downward_nir_beam(1)))       * (1.-visible_fraction_dir)
   DW_nirlo_diffuse = max(0.0,sngl(downward_nir_diffuse(1)))    * (1.-visible_fraction_dif)
   UW_nirhi_beam    = max(0.0,sngl(upward_nir_beam(ncoh+1)))    * (1.-visible_fraction_dir)
   UW_nirhi_diffuse = max(0.0,sngl(upward_nir_diffuse(ncoh+1))) * (1.-visible_fraction_dif)
   !---------------------------------------------------------------------------------------!

   return
end subroutine sw_twostream_clump
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will solve the within canopy two-stream radiation for shortwave,      !
! considering the light absorption by leaves, and branches if the user wants so.  In both  !
! cases we consider the clumpiness effect (i.e., the fact that leaves and branches aren't  !
! randomly distributed).   In this case, we solve the canopy layer-by-layer rather than    !
! cohort-by-cohort.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine sw_twostream_layer(salb,scosz,scosaoi,ncoh,pft,TAI,canopy_area,hgttop,hgtbot    &
                             ,PAR_beam_flip,PAR_diffuse_flip,SW_abs_beam_flip              &
                             ,SW_abs_diffuse_flip,DW_vislo_beam,DW_vislo_diffuse           &
                             ,UW_vishi_beam,UW_vishi_diffuse,DW_nirlo_beam                 &
                             ,DW_nirlo_diffuse,UW_nirhi_beam,UW_nirhi_diffuse              &
                             ,beam_level,diff_level,lambda_coh,lambda_tot)

   use ed_max_dims          , only : n_pft                   ! ! intent(in)
   use pft_coms             , only : clumping_factor         & ! intent(in)
                                   , phenology               ! ! intent(in)
   use canopy_radiation_coms, only : diffuse_backscatter_nir & ! intent(in)
                                   , diffuse_backscatter_vis & ! intent(in)
                                   , leaf_scatter_nir        & ! intent(in)
                                   , leaf_scatter_vis        & ! intent(in)
                                   , visible_fraction_dir    & ! intent(in)
                                   , visible_fraction_dif    & ! intent(in)
                                   , cosz_min8               ! ! intent(in)
   use canopy_layer_coms    , only : dzcan8                  & ! intent(in)
                                   , zztop8                  & ! intent(in)
                                   , zzmid8                  & ! intent(in)
                                   , zzbot8                  & ! intent(in)
                                   , zztop0i8                & ! intent(in)
                                   , ehgti8                  & ! intent(in)
                                   , ncanlyr                 & ! intent(in)
                                   , ncanlyrp1               & ! intent(in)
                                   , ncanlyrt2               & ! intent(in)
                                   , indx                    & ! intent(inout)
                                   , populated               & ! intent(inout)
                                   , dzcanpop                & ! intent(inout)
                                   , matal                   & ! intent(inout)
                                   , mastermat               & ! intent(inout)
                                   , masmatcp                & ! intent(inout)
                                   , opencan8                & ! intent(inout)
                                   , expkl_top               & ! intent(inout)
                                   , expkl_bot               & ! intent(inout)
                                   , expamk_top              & ! intent(inout)
                                   , expamk_bot              & ! intent(inout)
                                   , expapk_top              & ! intent(inout)
                                   , expapk_bot              & ! intent(inout)
                                   , A_top                   & ! intent(inout)
                                   , A_bot                   & ! intent(inout)
                                   , B_top                   & ! intent(inout)
                                   , B_bot                   & ! intent(inout)
                                   , C_top                   & ! intent(inout)
                                   , C_bot                   & ! intent(inout)
                                   , F_top                   & ! intent(inout)
                                   , F_bot                   & ! intent(inout)
                                   , G_top                   & ! intent(inout)
                                   , G_bot                   & ! intent(inout)
                                   , H_top                   & ! intent(inout)
                                   , H_bot                   & ! intent(inout)
                                   , beam_bot                & ! intent(inout)
                                   , beam_mid                & ! intent(inout)
                                   , beam_bot_crown          & ! intent(inout)
                                   , upward_vis_beam         & ! intent(inout)
                                   , upward_vis_diffuse      & ! intent(inout)
                                   , upward_nir_beam         & ! intent(inout)
                                   , upward_nir_diffuse      & ! intent(inout)
                                   , downward_nir_beam       & ! intent(inout)
                                   , downward_nir_diffuse    & ! intent(inout)
                                   , downward_vis_beam       & ! intent(inout)
                                   , downward_vis_diffuse    & ! intent(inout)
                                   , mastervec_beam          & ! intent(inout)
                                   , masveccp_beam           & ! intent(inout)
                                   , mastervec_diffuse       & ! intent(inout)
                                   , masveccp_diffuse        & ! intent(inout)
                                   , PAR_beam_layer          & ! intent(inout)
                                   , PAR_diffuse_layer       & ! intent(inout)
                                   , SW_abs_beam_layer       & ! intent(inout)
                                   , SW_abs_diffuse_layer    & ! intent(inout)
                                   , zero_canopy_layer       ! ! subroutine
   use consts_coms          , only : tiny_num8               ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer     , dimension(ncoh) , intent(in)  :: pft
   integer                       , intent(in)  :: ncoh
   real(kind=8), dimension(ncoh) , intent(in)  :: TAI
   real(kind=8), dimension(ncoh) , intent(in)  :: canopy_area
   real(kind=8), dimension(ncoh) , intent(in)  :: hgttop
   real(kind=8), dimension(ncoh) , intent(in)  :: hgtbot
   real(kind=4)                  , intent(in)  :: salb
   real(kind=4)                  , intent(in)  :: scosz
   real(kind=4)                  , intent(in)  :: scosaoi
   real(kind=4), dimension(ncoh) , intent(out) :: PAR_beam_flip
   real(kind=4), dimension(ncoh) , intent(out) :: PAR_diffuse_flip
   real(kind=4), dimension(ncoh) , intent(out) :: SW_abs_beam_flip
   real(kind=4), dimension(ncoh) , intent(out) :: SW_abs_diffuse_flip
   real(kind=4)                  , intent(out) :: UW_vishi_beam
   real(kind=4)                  , intent(out) :: UW_vishi_diffuse
   real(kind=4)                  , intent(out) :: UW_nirhi_beam
   real(kind=4)                  , intent(out) :: UW_nirhi_diffuse
   real(kind=4)                  , intent(out) :: DW_vislo_beam
   real(kind=4)                  , intent(out) :: DW_vislo_diffuse
   real(kind=4)                  , intent(out) :: DW_nirlo_beam
   real(kind=4)                  , intent(out) :: DW_nirlo_diffuse
   real(kind=8), dimension(ncoh) , intent(out) :: beam_level
   real(kind=8), dimension(ncoh) , intent(out) :: diff_level
   real(kind=8), dimension(ncoh) , intent(out) :: lambda_coh
   real(kind=8)                  , intent(out) :: lambda_tot
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: il
   integer                                     :: ico
   integer                                     :: jco
   integer                                     :: ipft
   integer                                     :: iband
   integer                                     :: i
   integer                                     :: j
   integer                                     :: ind
   integer                                     :: nactlyr
   integer                                     :: nactlyrp1
   integer                                     :: nactlyrt2
   integer     , dimension(ncoh)               :: kapartial
   integer     , dimension(ncoh)               :: kafull
   integer     , dimension(ncoh)               :: kzpartial
   integer     , dimension(ncoh)               :: kzfull
   integer     , dimension(ncoh)               :: nlyr_coh
   real(kind=8), dimension(n_pft)              :: leaf_scatter
   real(kind=8), dimension(n_pft)              :: diffuse_backscatter
   real(kind=8)                                :: alb
   real(kind=8)                                :: cosz
   real(kind=8)                                :: cosaoi
   real(kind=8)                                :: lambda
   real(kind=8)                                :: tad
   real(kind=8)                                :: this_cai
   real(kind=8)                                :: this_tai
   real(kind=8)                                :: this_eff_tai
   real(kind=8)                                :: this_ext
   real(kind=8)                                :: layer_ext
   real(kind=8)                                :: beam_backscatter
   real(kind=8)                                :: eta
   real(kind=8)                                :: zeta
   real(kind=8)                                :: exk
   real(kind=8)                                :: exki
   real(kind=8)                                :: zetai
   real(kind=8)                                :: d
   real(kind=8)                                :: rhoo
   real(kind=8)                                :: sigma
   real(kind=8)                                :: iota
   real(kind=8)                                :: source_bot
   real(kind=8)                                :: source_top
   !----- Variables that must be allocated every time. ------------------------------------!
   real(kind=8), dimension(:,:)  , allocatable :: extinct_full
   real(kind=8), dimension(:,:)  , allocatable :: extinct_half
   real(kind=8), dimension(:,:)  , allocatable :: canfrac_lyr
   real(kind=8), dimension(:,:)  , allocatable :: tai_lyr
   real(kind=8), dimension(:,:)  , allocatable :: eff_tai_lyr
   !------ External functions. ------------------------------------------------------------!
   real(kind=4)                  , external    :: sngloff
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These must be allocated and de-allocated every time as the number of cohorts is   !
   ! unbounded and may change.                                                             !
   !---------------------------------------------------------------------------------------!
   allocate(extinct_full      (ncoh,ncanlyr))
   allocate(extinct_half      (ncoh,ncanlyr))
   allocate(canfrac_lyr       (ncoh,ncanlyr))
   allocate(tai_lyr           (ncoh,ncanlyr))
   allocate(eff_tai_lyr       (ncoh,ncanlyr))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Reset the temporary variables.                                                    !
   !---------------------------------------------------------------------------------------!
   call zero_canopy_layer('sw_twostream_layer')
   extinct_full          (:,:) = 0.d0
   extinct_half          (:,:) = 0.d0
   canfrac_lyr           (:,:) = 0.d0
   tai_lyr               (:,:) = 0.d0
   eff_tai_lyr           (:,:) = 0.d0
   !---------------------------------------------------------------------------------------!



   !----- Convert input variable to double precision. -------------------------------------!
   alb    = dble(salb)
   cosz   = max(cosz_min8,dble(scosz))
   cosaoi = max(cosz_min8,dble(scosaoi))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the crown area for each layer, and determine the layer bounds for each       !
   ! cohort.                                                                               !
   !---------------------------------------------------------------------------------------!
   do ico = 1,ncoh
      !------ Find the layer bounds. ------------------------------------------------------!
      kapartial (ico) = min(ncanlyr,floor  ((hgtbot(ico) * zztop0i8)**ehgti8) + 1)
      kafull    (ico) = min(ncanlyr,ceiling((hgtbot(ico) * zztop0i8)**ehgti8) + 1)
      kzpartial (ico) = min(ncanlyr,ceiling((hgttop(ico) * zztop0i8)**ehgti8))
      kzfull    (ico) = min(ncanlyr,floor  ((hgttop(ico) * zztop0i8)**ehgti8))


      !------------------------------------------------------------------------------------!
      !     Find the actual crown area that each cohort has in each layer.  The allometry  !
      ! value should be regarded as the potential value if no other cohort is present.     !
      ! In case the total crown area for a given layer exceeds one, then we squeeze the    !
      ! plants in that layer so the cohorts don't sit on top of each other.                !
      !------------------------------------------------------------------------------------!
      do il = kapartial(ico),kzpartial(ico)
         populated(il)        = .true.
         canfrac_lyr (ico,il) = canopy_area(ico)
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Find the open canopy fraction.  This is done in three steps:                       !
   ! 1. Find the closed canopy fraction.                                                   !
   ! 2. Check the size. If it exceeds one, scale down the canopy area for each layer and   !
   !    force it to be one.                                                                !
   ! 3. Set the open canopy area to be 1- closed canopy area.                              !
   !---------------------------------------------------------------------------------------!
   do il=1,ncanlyr
      opencan8(il) = sum(canfrac_lyr(:,il))
      if (opencan8(il) > 1.d0) then
         canfrac_lyr(:,il) = canfrac_lyr(:,il) / opencan8(il)
         opencan8(il)      = 1.d0
      end if
      opencan8(il) = 1.d0 - opencan8(il)
   end do
   !---------------------------------------------------------------------------------------!


   lambda        = 5.d-1/cosaoi
   lambda_coh(:) = 0.d0 
   do ico=1,ncoh
      !------ Alias for current PFT. ------------------------------------------------------!
      ipft = pft(ico)

      !------ Find the "tree" area density. -----------------------------------------------!
      tad              = tai(ico) / (hgttop(ico) - hgtbot(ico))

      !------------------------------------------------------------------------------------!
      !    Integrate the extinction coefficients for partial layers.  Here we scale the    !
      ! leaf area density to the layer (partial or full), so the total leaf area index is  !
      ! preserved.  The only special case is when the top and bottom of the crown are in   !
      ! the same layer, in which case we simply use the leaf area index.  If the user is   !
      ! running with branch thermodynamics, then it should be total (leaf + branch) area   !
      ! index instead of leaf area index.                                                  !
      !------------------------------------------------------------------------------------!
      if (kapartial(ico) == kzpartial(ico)) then
         il                    = kapartial(ico)
         tai_lyr      (ico,il) = tai(ico)
         eff_tai_lyr  (ico,il) = clumping_factor(ipft) * tai_lyr (ico,il)
         this_ext              = lambda * eff_tai_lyr(ico,il) / canfrac_lyr(ico,il)
         extinct_full (ico,il) = exp( - this_ext)
         extinct_half (ico,il) = exp( - 5.d-1  * this_ext )
         lambda_coh      (ico) = lambda_coh(ico) + this_ext 
      else
         !------ Fully vegetated layers. --------------------------------------------------!
         do il=kafull(ico),kzfull(ico)
            tai_lyr      (ico,il) = tad * dzcan8(il)
            eff_tai_lyr  (ico,il) = clumping_factor(ipft) * tai_lyr (ico,il)
            this_ext              = lambda * eff_tai_lyr(ico,il) / canfrac_lyr(ico,il)
            extinct_full (ico,il) = exp( - this_ext)
            extinct_half (ico,il) = exp( - 5.d-1 * this_ext)
            lambda_coh      (ico) = lambda_coh(ico) + this_ext 
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !  Partial layer at the bottom of the crown.  This shouldn't be done if the       !
         ! partial and full indices are the same, otherwise we double account the layer.   !
         !---------------------------------------------------------------------------------!
         if (kapartial(ico) /= kafull(ico)) then
            il                    = kapartial(ico)
            tai_lyr      (ico,il) = tad * (zztop8(il)-hgtbot(ico))
            eff_tai_lyr  (ico,il) = clumping_factor(ipft) * tai_lyr (ico,il)
            this_ext              = lambda * eff_tai_lyr(ico,il) / canfrac_lyr(ico,il)
            extinct_full (ico,il) = exp( - this_ext)
            extinct_half (ico,il) = exp( - 5.d-1 * this_ext)
            lambda_coh      (ico) = lambda_coh(ico) + this_ext 
         end if
         if (kzpartial(ico) /= kzfull(ico)) then
            il                    = kzpartial(ico)
            tai_lyr      (ico,il) = tad * (hgttop(ico) - zzbot8(il))
            eff_tai_lyr  (ico,il) = clumping_factor(ipft) * tai_lyr (ico,il)
            this_ext              = lambda * eff_tai_lyr(ico,il) / canfrac_lyr(ico,il)
            extinct_full (ico,il) = exp( - this_ext)
            extinct_half (ico,il) = exp( - 5.d-1 * this_ext)
            lambda_coh      (ico) = lambda_coh(ico) + this_ext 
         end if
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !------ Find the average extinction coefficients. --------------------------------------!
   lambda_tot       = sum(lambda_coh(:)) / sum(tai(:))
   lambda_coh       = lambda_coh(:) / tai(:)
   !---------------------------------------------------------------------------------------!




   !------ Find the beam backscattering coefficient. --------------------------------------!
   beam_backscatter = (5.d-1 + cosz) * (1.0d0 - cosz*log(1.0d0+1.0d0/cosz))
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Pack the information so we skip empty layers in the middle, thus avoinding         !
   ! singularities.  From this point on we must refer to nactlyr, not ncanlyr.             !
   !---------------------------------------------------------------------------------------!
   !----- Find out how many active layers we have. ----------------------------------------!
   nactlyr   = count(populated)
   nactlyrp1 = nactlyr+1
   nactlyrt2 = 2 * nactlyr
   !----- Squeeze the arrays. -------------------------------------------------------------!
   opencan8(1:nactlyr) = pack(opencan8, populated)
   dzcanpop(1:nactlyr) = pack(dzcan8  , populated)
   if (nactlyrp1 < ncanlyr) then
      opencan8(nactlyrp1:ncanlyr) = 0.d0
      dzcanpop(nactlyrp1:ncanlyr) = 0.d0
   end if
   do ico = 1,ncoh
      tai_lyr      (ico,1:nactlyr) = pack(tai_lyr      (ico,:), populated)
      eff_tai_lyr  (ico,1:nactlyr) = pack(eff_tai_lyr  (ico,:), populated)
      canfrac_lyr  (ico,1:nactlyr) = pack(canfrac_lyr  (ico,:), populated)
      extinct_full (ico,1:nactlyr) = pack(extinct_full (ico,:), populated)
      extinct_half (ico,1:nactlyr) = pack(extinct_half (ico,:), populated)
      if (nactlyrp1 < ncanlyr) then
         tai_lyr      (ico,nactlyrp1:ncanlyr) = 0.d0
         eff_tai_lyr  (ico,nactlyrp1:ncanlyr) = 0.d0
         canfrac_lyr  (ico,nactlyrp1:ncanlyr) = 0.d0
         extinct_full (ico,nactlyrp1:ncanlyr) = 0.d0
         extinct_half (ico,nactlyrp1:ncanlyr) = 0.d0
      end if
      !----- Count the number of layers that have leaves/branches for this cohort. --------!
      nlyr_coh       (ico) = count(canfrac_lyr(ico,:) > 0.d0)
   end do
   !---------------------------------------------------------------------------------------!



   !----- Loop over bands (currently Visible and near infrared). --------------------------!
   bandloop: do iband = 1,2
      select case(iband)
      case (1) !----- Visible (or PAR). ---------------------------------------------------!
         do ipft = 1,n_pft
            leaf_scatter(ipft)        = dble(leaf_scatter_vis(ipft))
            diffuse_backscatter(ipft) = dble(diffuse_backscatter_vis(ipft))
         end do
      case (2) !----- Near infrared (or NIR). ---------------------------------------------!
         do ipft = 1,n_pft
            leaf_scatter(ipft)        = dble(leaf_scatter_nir(ipft))
            diffuse_backscatter(ipft) = dble(diffuse_backscatter_nir(ipft))
         end do
      end select
      
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     Calculate more factors for this band.  We start by calculating the forcings,   !
      ! using the effective tree area index, considering the leaf clumpiness and, if the   !
      ! user asked for it, the branchwood area index (in case the user doesn't want,       !
      ! WAI=0.).  We distinguish between the light underneath a crown (beam_bot_crown) and !
      ! the total light passing through a layer (beam_bot).                                !
      !------------------------------------------------------------------------------------!
      !----- Start with the highest layer, moving downwards. ------------------------------!
      beam_bot_crown(nactlyr) = sum(canfrac_lyr(:,nactlyr) * extinct_full(:,nactlyr))      &
                              / (1.d0 - opencan8(nactlyr) )
      beam_bot      (nactlyr) = opencan8(nactlyr)                                          &
                              + sum(canfrac_lyr(:,nactlyr) * extinct_full(:,nactlyr))
      beam_mid      (nactlyr) = opencan8(nactlyr)                                          &
                              + sum(canfrac_lyr(:,nactlyr) * extinct_half(:,nactlyr))
      do il = nactlyr-1,1,-1
         beam_bot_crown  (il) = beam_bot(il+1)                                             &
                              * sum(canfrac_lyr(:,il) * extinct_full(:,il))                &
                              / (1.d0 - opencan8(il))
         beam_bot        (il) = beam_bot(il+1)                                             &
                              * ( opencan8(il)                                             &
                                + sum(canfrac_lyr(:,il) * extinct_full(:,il)) )
         beam_mid        (il) = beam_bot(il+1)                                             &
                              * ( opencan8(il)                                             &
                                + sum(canfrac_lyr(:,il) * extinct_half(:,il)) )
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      srclyrloop: do il=1,nactlyr
         !---------------------------------------------------------------------------------!
         !     Find the terms by averaging them using the effective TAI as the weighting   !
         ! factor.  This is not the best way of doing it, but it is the simplest I could   !
         ! think, because it retains the fact that eta and zeta are constants within each  !
         ! layer.  In addition, we average sigma and iota, coefficients related to         !
         ! extinction and transmittance.                                                   !
         !---------------------------------------------------------------------------------!
         eta          = 0.d0
         zeta         = 0.d0
         sigma        = 0.d0
         iota         = 0.d0
         this_tai     = sum(tai_lyr(:,il))
         this_eff_tai = sum(eff_tai_lyr(:,il))
         do ico=1,ncoh
            ipft       = pft(ico)
            eta        = eta  + eff_tai_lyr(ico,il)                                        &
                              * ( 1.0d0 - (1.0d0 - diffuse_backscatter(ipft))              &
                                * leaf_scatter(ipft))
            zeta       = zeta + eff_tai_lyr(ico,il)                                        &
                              * leaf_scatter(ipft) * diffuse_backscatter(ipft)
            sigma      = sigma + eff_tai_lyr(ico,il) * lambda
            iota       = iota  + eff_tai_lyr(ico,il) * lambda * leaf_scatter(ipft)         &
                               * beam_backscatter
         end do
         !----- Normalise eta, zeta, and sigma (but retain clumping factor in them). ------!
         eta   =   eta / this_tai
         zeta  =  zeta / this_tai
         sigma = sigma / this_tai
         iota  =  iota / this_tai
         
         !----- Find derived properties. --------------------------------------------------!
         exk   = sqrt(eta*eta - zeta*zeta)
         exki  = 1.0d0/exk
         zetai = 1.0d0/zeta

         !----- Sources. ------------------------------------------------------------------!
         source_bot = iota * beam_bot_crown(il)
         source_top = source_bot * exp(lambda * this_eff_tai)

         !----- Forcing coefficients. -----------------------------------------------------!
         rhoo  = - (zeta + eta + sigma) * source_bot


         !----- Calculate exponentials only once. -----------------------------------------!
         expkl_bot (il) = 1.d0
         expkl_top (il) = exp(       exk  * this_tai)
         expamk_bot(il) = 1.d0
         expamk_top(il) = exp((sigma-exk) * this_tai)
         expapk_bot(il) = 1.d0
         expapk_top(il) = exp((sigma+exk) * this_tai)


         !---------------------------------------------------------------------------------!
         !     Find the terms for the integration.                                         !
         !---------------------------------------------------------------------------------!
         A_bot(il)      = -source_bot*zetai
         A_top(il)      = -source_top*zetai                                                &
                        + 5.d-1*zetai*(eta*exki-1.0d0)*expkl_top(il)*rhoo/(sigma-exk)      &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 5.d-1*zetai*(eta*exki+1.0d0) / expkl_top(il)*rhoo/(sigma+exk)    &
                        * (expapk_top(il)-expapk_bot(il))
         B_bot(il)      = 5.d-1*zetai*(eta*exki-1.0d0)
         B_top(il)      = 5.d-1*zetai*(eta*exki-1.0d0)*expkl_top(il)
         C_bot(il)      = -5.d-1*zetai*(eta*exki+1.0d0)
         C_top(il)      = -5.d-1*zetai*(eta*exki+1.0d0)/expkl_top(il)
         F_bot(il)      = 0.0d0
         F_top(il)      = 5.d-1*exki*expkl_top(il)*rhoo/(sigma-exk)                        &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 5.d-1*exki/expkl_top(il)*rhoo/(sigma+exk)                        &
                        * (expapk_top(il)-expapk_bot(il))
         G_bot(il)      = 5.d-1*exki
         G_top(il)      = 5.d-1*exki*expkl_top(il)
         H_bot(il)      = -5.d-1*exki
         H_top(il)      = -5.d-1*exki/expkl_top(il)
         !---------------------------------------------------------------------------------!
      end do srclyrloop
      !------------------------------------------------------------------------------------!




      !----- Initialise the matrix coefficients. ------------------------------------------!
      do j=1,nactlyrt2
         do i=1,nactlyrt2
            masmatcp(i,j) = 0.0d0
         end do
         mastervec_beam(j)    = 0.0d0
         masveccp_beam(j)     = 0.0d0
         mastervec_diffuse(j) = 0.0d0
         masveccp_diffuse(j)  = 0.0d0
      end do
      !------------------------------------------------------------------------------------!

      masmatcp(1,1)        = G_top(nactlyr)
      masmatcp(1,2)        = H_top(nactlyr)
      mastervec_beam(1)    = -F_top(nactlyr)
      mastervec_diffuse(1) = 1.0d0
      masveccp_beam(1)     = mastervec_beam(1)
      masveccp_diffuse(1)  = mastervec_diffuse(1)

      do i=2,nactlyrt2-2,2
         masmatcp(i,i-1)      = G_bot(nint(real(nactlyrt2-i+2)*0.5))
         masmatcp(i,i)        = H_bot(nint(real(nactlyrt2-i+2)*0.5))

         masmatcp(i,i+1)      = - G_top(nint(real(nactlyrt2-i)*0.5))
         masmatcp(i,i+2)      = - H_top(nint(real(nactlyrt2-i)*0.5))
         mastervec_beam(i)    = - F_bot(nint(real(nactlyrt2-i+2)*0.5))                     &
                                + F_top(nint(real(nactlyrt2-i)*0.5))
         mastervec_diffuse(i) = 0.0d0
         masveccp_beam(i)     = mastervec_beam(i)
         masveccp_diffuse(i)  = mastervec_diffuse(i)
      end do

      do i=3,nactlyrt2-1,2
         masmatcp(i,i-2)      = B_bot(nint(real(nactlyrt2-i+3)*0.5))
         masmatcp(i,i-1)      = C_bot(nint(real(nactlyrt2-i+3)*0.5))
         masmatcp(i,i)        = - B_top(nint(real(nactlyrt2-i+1)*0.5))
         masmatcp(i,i+1)      = - C_top(nint(real(nactlyrt2-i+1)*0.5))
         mastervec_beam(i)    = - A_bot(nint(real(nactlyrt2-i+3)*0.5))                     &
                                + A_top(nint(real(nactlyrt2-i+1)*0.5))
         masveccp_beam(i)     = mastervec_beam(i)
         mastervec_diffuse(i) = 0.0d0
         masveccp_diffuse(i)  = mastervec_diffuse(i)
      end do
      masmatcp(nactlyrt2,nactlyrt2-1) = B_bot(1) - alb * G_bot(1)
      masmatcp(nactlyrt2,nactlyrt2)   = C_bot(1) - alb * H_bot(1)
      mastervec_beam(nactlyrt2)       = - A_bot(1)+alb * beam_bot(1)
      masveccp_beam(nactlyrt2)        = mastervec_beam(nactlyrt2)
      mastervec_diffuse(nactlyrt2)    = 0.0d0
      masveccp_diffuse(nactlyrt2)     = mastervec_diffuse(nactlyrt2)
      
      !----- Prep for inversion. ----------------------------------------------------------!
      mastermat(1,1) = 0.d0
      mastermat(1,2) = 0.d0
      mastermat(1,3) = masmatcp(1,1)
      mastermat(1,4) = masmatcp(1,2)
      mastermat(1,5) = 0.d0

      do i=2,nactlyrt2-2,2
         mastermat(i,1) = 0.d0
         mastermat(i,2) = masmatcp(i,i-1)
         mastermat(i,3) = masmatcp(i,i)
         mastermat(i,4) = masmatcp(i,i+1)
         mastermat(i,5) = masmatcp(i,i+2)
      end do

      do i=3,nactlyrt2-1,2
         mastermat(i,1) = masmatcp(i,i-2)
         mastermat(i,2) = masmatcp(i,i-1)
         mastermat(i,3) = masmatcp(i,i)
         mastermat(i,4) = masmatcp(i,i+1)
         mastermat(i,5) = 0.d0
      end do

      mastermat(nactlyrt2,1) = 0.d0
      mastermat(nactlyrt2,2) = masmatcp(nactlyrt2,nactlyrt2-1)
      mastermat(nactlyrt2,3) = masmatcp(nactlyrt2,nactlyrt2)
      mastermat(nactlyrt2,4) = 0.d0
      mastermat(nactlyrt2,5) = 0.d0
      
      !----- Invert the matrix. -----------------------------------------------------------!
      call bandec(mastermat,ncanlyrt2,nactlyrt2,2,2,matal,indx,d)

      !----- Backsubstitute for beam and diffuse. -----------------------------------------!
      call banbks(mastermat,ncanlyrt2,nactlyrt2,2,2,matal,indx,mastervec_beam)
      call banbks(mastermat,ncanlyrt2,nactlyrt2,2,2,matal,indx,mastervec_diffuse)
      
      !----- Improve the solution. --------------------------------------------------------!
      call mprove(masmatcp,mastermat,matal,ncanlyrt2,nactlyrt2,5,2,indx                    &
                 ,masveccp_beam,mastervec_beam)
      call mprove(masmatcp,mastermat,matal,ncanlyrt2,nactlyrt2,5,2,indx                    &
                 ,masveccp_diffuse,mastervec_diffuse)

      select case (iband)
      case (1) !---- Visible (or PAR) band. -----------------------------------------------!
         do i=3,nactlyrt2-1,2
            ind                     = nint(real(nactlyrt2-i+1)*0.5)+1
            upward_vis_beam(ind)    = A_bot(ind) + B_bot(ind) * mastervec_beam(i-2)        &
                                    + C_bot(ind) * mastervec_beam(i-1)
            upward_vis_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2)                  &
                                    + C_bot(ind) * mastervec_diffuse(i-1)
         end do

         do i=2,nactlyrt2-2,2
            ind                       = nint(real(nactlyrt2-i)*0.5)+1
            downward_vis_beam(ind)    = beam_bot(ind) + F_bot(ind)                         &
                                      + H_bot(ind) * mastervec_beam(i)                     &
                                      + G_bot(ind) * mastervec_beam(i-1)
            downward_vis_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)                  &
                                      + G_bot(ind) * mastervec_diffuse(i-1)
         end do

         upward_vis_beam(nactlyr+1)      = B_top(nactlyr) * mastervec_beam(1)              &
                                         + C_top(nactlyr) * mastervec_beam(2)              &
                                         + A_top(nactlyr)
         upward_vis_diffuse(nactlyr+1)   = B_top(nactlyr) * mastervec_diffuse(1)           &
                                         + C_top(nactlyr) * mastervec_diffuse(2)
         downward_vis_beam(nactlyr+1)    = 1.0d0
         downward_vis_diffuse(nactlyr+1) = 1.0d0
         downward_vis_beam(1)            = G_bot(1) * mastervec_beam(nactlyrt2-1)          &
                                         + H_bot(1) * mastervec_beam(nactlyrt2) + F_bot(1) &
                                         + beam_bot(1)
         downward_vis_diffuse(1)         = G_bot(1) * mastervec_diffuse(nactlyrt2-1)       &
                                         + H_bot(1) * mastervec_diffuse(nactlyrt2)
         upward_vis_beam(1)              = alb * downward_vis_beam(1)
         upward_vis_diffuse(1)           = alb * downward_vis_diffuse(1)
      case (2)

         do i=3,nactlyrt2-1,2
            ind                     = nint(real(nactlyrt2-i+1)*0.5)+1
            upward_nir_beam(ind)    = A_bot(ind) + B_bot(ind) * mastervec_beam(i-2)        &
                                    + C_bot(ind) * mastervec_beam(i-1)
            upward_nir_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2)                  &
                                    + C_bot(ind) * mastervec_diffuse(i-1)
         end do

         do i=2,nactlyrt2-2,2
            ind                       = nint(real(nactlyrt2-i)*0.5)+1
            downward_nir_beam(ind)    = beam_bot(ind) + F_bot(ind)                         &
                                      + H_bot(ind) * mastervec_beam(i)                     &
                                      + G_bot(ind) * mastervec_beam(i-1)
            downward_nir_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)                  &
                                      + G_bot(ind) * mastervec_diffuse(i-1)
         end do

         upward_nir_beam(nactlyr+1)      = B_top(nactlyr) * mastervec_beam(1)              &
                                         + C_top(nactlyr) * mastervec_beam(2)              &
                                         + A_top(nactlyr)
         upward_nir_diffuse(nactlyr+1)   = B_top(nactlyr) * mastervec_diffuse(1)           &
                                         + C_top(nactlyr) * mastervec_diffuse(2)
         downward_nir_beam(nactlyr+1)    = 1.0d0
         downward_nir_diffuse(nactlyr+1) = 1.0d0
         downward_nir_beam(1)            = G_bot(1) * mastervec_beam(nactlyrt2-1)          &
                                         + H_bot(1) * mastervec_beam(nactlyrt2) + F_bot(1) &
                                         + beam_bot(1)
         downward_nir_diffuse(1)         = G_bot(1) * mastervec_diffuse(nactlyrt2-1)       &
                                         + H_bot(1) * mastervec_diffuse(nactlyrt2)
         upward_nir_beam(1)              = alb*downward_nir_beam(1)
         upward_nir_diffuse(1)           = alb*downward_nir_diffuse(1)
      end select
   end do bandloop
   
   do il=1,nactlyr
      PAR_beam_layer(il)       = visible_fraction_dir                                      &
                               * sngl( downward_vis_beam(il+1)-downward_vis_beam(il)       &
                                     + upward_vis_beam(il)-upward_vis_beam(il+1))
      PAR_diffuse_layer(il)    = visible_fraction_dif                                      &
                               * sngl( downward_vis_diffuse(il+1)-downward_vis_diffuse(il) &
                                     + upward_vis_diffuse(il)-upward_vis_diffuse(il+1))
      SW_abs_beam_layer(il)    = PAR_beam_layer(il) + (1.0 - visible_fraction_dir)         &
                               * sngl( downward_nir_beam(il+1)-downward_nir_beam(il)       &
                                     + upward_nir_beam(il)-upward_nir_beam(il+1))
      SW_abs_diffuse_layer(il) = PAR_diffuse_layer(il) + (1.0 - visible_fraction_dif)      &
                               * sngl( downward_nir_diffuse(il+1)                          &
                                     - downward_nir_diffuse(il)                            &
                                     + upward_nir_diffuse(il)                              &
                                     - upward_nir_diffuse(il+1) )

      !----- Ensure that we don't get any negative radiation... ---------------------------!
      PAR_beam_layer(il)       = max(0.0,PAR_beam_layer(il)       )
      PAR_diffuse_layer(il)    = max(0.0,PAR_diffuse_layer(il)    )
      SW_abs_beam_layer(il)    = max(0.0,SW_abs_beam_layer(il)    )
      SW_abs_diffuse_layer(il) = max(0.0,SW_abs_diffuse_layer(il) )
   end do
   
   !----- Copying to the output variables. ------------------------------------------------!
   DW_vislo_beam    = max(0.0,sngl(downward_vis_beam(1)))                                  &
                    * visible_fraction_dir
   DW_vislo_diffuse = max(0.0,sngl(downward_vis_diffuse(1)))                               &
                    * visible_fraction_dif
   UW_vishi_beam    = max(0.0,sngl(upward_vis_beam(nactlyr+1)))                            &
                    * visible_fraction_dir
   UW_vishi_diffuse = max(0.0,sngl(upward_vis_diffuse(nactlyr+1)))                         &
                    * visible_fraction_dif
   DW_nirlo_beam    = max(0.0,sngl(downward_nir_beam(1)))                                  &
                    * (1.-visible_fraction_dir)
   DW_nirlo_diffuse = max(0.0,sngl(downward_nir_diffuse(1)))                               &
                    * (1.-visible_fraction_dif)
   UW_nirhi_beam    = max(0.0,sngl(upward_nir_beam(nactlyr+1)))                            &
                    * (1.-visible_fraction_dir)
   UW_nirhi_diffuse = max(0.0,sngl(upward_nir_diffuse(nactlyr+1)))                         &
                    * (1.-visible_fraction_dif)

   !---------------------------------------------------------------------------------------!
   !     Integrate the total amount of energy for each cohort.                             !
   !---------------------------------------------------------------------------------------!
   PAR_beam_flip       (:) = 0.0
   PAR_diffuse_flip    (:) = 0.0
   SW_abs_beam_flip    (:) = 0.0
   SW_abs_diffuse_flip (:) = 0.0
   do il = 1, nactlyr
      this_eff_tai = sum(eff_tai_lyr(:,il))
      this_tai     = sum(tai_lyr    (:,il))
      do ico = 1,ncoh
         PAR_beam_flip(ico)       = PAR_beam_flip(ico) + PAR_beam_layer(il)                &
                                  * sngloff(eff_tai_lyr(ico,il) / this_eff_tai,tiny_num8)
         PAR_diffuse_flip(ico)    = PAR_diffuse_flip(ico) + PAR_diffuse_layer(il)          &
                                  * sngloff(tai_lyr(ico,il) / this_tai,tiny_num8)
         SW_abs_beam_flip(ico)    = SW_abs_beam_flip(ico) + SW_abs_beam_layer(il)          &
                                  * sngloff(eff_tai_lyr(ico,il) / this_eff_tai,tiny_num8)
         SW_abs_diffuse_flip(ico) = SW_abs_diffuse_flip(ico) + SW_abs_diffuse_layer(il)    &
                                  * sngloff(tai_lyr(ico,il) / this_tai,tiny_num8)
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the light level for beam and diffuse radiation.                              !
   !---------------------------------------------------------------------------------------!
   beam_level(:) = 0.0
   diff_level(:) = SW_abs_diffuse_flip(:)
   do ico = 1,ncoh
      do il=kapartial(ico),kzpartial(ico)
         beam_level(ico) = beam_level(ico)                                                 &
                         + sngloff( beam_mid(il) * dzcanpop(il)                            &
                                  / ( zztop8(kzpartial(ico)) - zzbot8(kapartial(ico)) )    &
                                  , tiny_num8 )
      end do
   end do
   !---------------------------------------------------------------------------------------!





   !------ De-allocate the temporary structures. ------------------------------------------!
   deallocate(extinct_full        )
   deallocate(extinct_half        )
   deallocate(canfrac_lyr         )
   deallocate(tai_lyr             )
   deallocate(eff_tai_lyr         )
   !---------------------------------------------------------------------------------------!

   return
end subroutine sw_twostream_layer
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the matrix pseudo-inverse.                                       !
!------------------------------------------------------------------------------------------!
subroutine bandec(a,ntot,nuse,m1,m2,al,indx,d)
   use consts_coms, only : tiny_num8
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                              , intent(in)    :: ntot
   integer                              , intent(in)    :: nuse
   integer                              , intent(in)    :: m1
   integer                              , intent(in)    :: m2
   integer     , dimension(ntot)        , intent(inout) :: indx
   real(kind=8), dimension(ntot,m1+m2+1), intent(inout) :: a
   real(kind=8), dimension(ntot,m1)     , intent(inout) :: al
   real(kind=8)                         , intent(out)   :: d
   !----- Local variables -----------------------------------------------------------------!
   integer                                              :: i
   integer                                              :: j
   integer                                              :: k
   integer                                              :: l
   integer                                              :: mm
   real(kind=8)                                         :: tvar
   real(kind=8)                                         :: dum
   !---------------------------------------------------------------------------------------!


   mm=m1+m2+1
   l=m1
   do i=1,m1
      do j=m1+2-i,mm
         a(i,j-l)=a(i,j)
      end do

      l=l-1

      do j=mm-l,mm
         a(i,j)=0.d0
      end do
   end do

   d=1.0d0

   l=m1
   do k=1,nuse
      dum=a(k,1)
      i=k
      if (l < nuse) l = l + 1
      do j=k+1,l
         if (abs(a(j,1)) > abs(dum)) then
            dum = a(j,1)
            i   = j
         end if
      end do
      indx(k) = i
      if (dum == 0.0d0) a(k,1) = tiny_num8
      if (i /= k) then
         d=-d
         do j=1,mm
            tvar   = a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = tvar
         end do
      end if

      do i=k+1,l
         dum       = a(i,1) / a(k,1)
         al(k,i-k) = dum
         do j=2,mm
            a(i,j-1) = a(i,j) - dum * a(k,j)
         end do

         a(i,mm) = 0.0d0

      end do
   end do
   return
end subroutine bandec
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine will perform the back substitution of the linear system.                 !
!------------------------------------------------------------------------------------------!
subroutine banbks(a,ntot,nuse,m1,m2,al,indx,b)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                              , intent(in)    :: ntot
   integer                              , intent(in)    :: nuse
   integer                              , intent(in)    :: m1
   integer                              , intent(in)    :: m2
   integer     , dimension(ntot)        , intent(in)    :: indx
   real(kind=8), dimension(ntot,m1+m2+1), intent(in)    :: a
   real(kind=8), dimension(ntot,m1)     , intent(in)    :: al
   real(kind=8), dimension(ntot)        , intent(inout) :: b
   !----- Local variables -----------------------------------------------------------------!
   integer                                              :: i
   integer                                              :: k
   integer                                              :: l
   integer                                              :: mm
   real(kind=8)                                         :: dum
   real(kind=8)                                         :: tvar
   !---------------------------------------------------------------------------------------!

   mm = m1 + m2 + 1
   l  = m1

   do k=1,nuse
      i=indx(k)

      if (i /= k) then
         tvar = b(k)
         b(k) = b(i)
         b(i) = tvar
      end if

      if(l < nuse) l = l + 1

      do i=k+1,l
         b(i) = b(i) - al(k,i-k) * b(k)
      end do
   end do

   l=1
   do i=nuse,1,-1
      dum = b(i)
      do k=2,l
         dum = dum - a(i,k) * b(k+i-1)
      end do
      b(i) = dum / a(i,1)
      if (l < mm) l = l + 1
   end do

   return
end subroutine banbks
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine will perform add a higher order term to improve the solution.          !
!------------------------------------------------------------------------------------------!
subroutine mprove(a,alud,matal,ntot,nuse,np,npp,indx,b,x)
   implicit none
   !----- Local constants -----------------------------------------------------------------!
   integer     , parameter :: nmax=100
   !----- Arguments -----------------------------------------------------------------------!
   integer                           , intent(in)    :: ntot
   integer                           , intent(in)    :: nuse
   integer                           , intent(in)    :: np
   integer                           , intent(in)    :: npp
   integer     , dimension(ntot)     , intent(in)    :: indx
   real(kind=8), dimension(ntot,ntot), intent(in)    :: a
   real(kind=8), dimension(ntot)     , intent(in)    :: b
   real(kind=8), dimension(ntot,np)  , intent(in)    :: alud
   real(kind=8), dimension(ntot,npp) , intent(in)    :: matal
   real(kind=8), dimension(ntot)     , intent(inout) :: x
   !----- Local variables -----------------------------------------------------------------!
   integer                                           :: i
   integer                                           :: j
   real(kind=8), dimension(ntot)                     :: r
   real(kind=8)                                      :: sdp
   !---------------------------------------------------------------------------------------!


   do i=1,nuse
      sdp = -b(i)
      do j=1,nuse
         sdp = sdp + a(i,j) * x(j)
      end do
      r(i) = sdp
   end do

   call banbks(alud,ntot,nuse,npp,npp,matal,indx,r)

   do i=1,nuse
      x(i) = x(i) - r(i)
   end do

   return
end subroutine mprove
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will solve the long wave radiation within the canopy.  This solves   !
! the two stream approach considering the size distribution of cohorts, acknowledging the  !
! effect of leaves, and, if that's the user's will, the branches.                          !
!------------------------------------------------------------------------------------------!
subroutine lw_twostream(ncoh,semgs,sT_grnd, pft,TAI,canopy_area,T_veg,lw_v_surf,lw_v_incid &
                       ,downward_lw_below_surf,downward_lw_below_incid                     &
                       ,upward_lw_below_surf,upward_lw_below_incid,upward_lw_above_surf    &
                       ,upward_lw_above_incid)
   use canopy_radiation_coms , only : emis_v          & ! intent(in)
                                    , mubar           ! ! intent(in)
   use pft_coms              , only : clumping_factor ! ! intent(in)
   use consts_coms           , only : stefan8         & ! intent(in)
                                    , tiny_num8       ! ! intent(in)
   use rk4_coms              , only : tiny_offset     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                      , intent(in)  :: ncoh        ! # of cohorts
   real                         , intent(in)  :: semgs       ! Gnd. emmissivity
   real                         , intent(in)  :: st_grnd     ! Gnd. temperature
   integer     , dimension(ncoh), intent(in)  :: pft         ! Plant functional type
   real(kind=8), dimension(ncoh), intent(in)  :: TAI         ! "Tree" Area Index
   real(kind=8), dimension(ncoh), intent(in)  :: T_veg       ! Vegetation temperature
   real(kind=8), dimension(ncoh), intent(in)  :: canopy_area ! canopy area
   real        , dimension(ncoh), intent(out) :: lw_v_surf
   real        , dimension(ncoh), intent(out) :: lw_v_incid
   real                         , intent(out) :: downward_lw_below_surf
   real                         , intent(out) :: upward_lw_above_surf
   real                         , intent(out) :: upward_lw_below_surf
   real                         , intent(out) :: downward_lw_below_incid
   real                         , intent(out) :: upward_lw_above_incid
   real                         , intent(out) :: upward_lw_below_incid
   !----- Local variables. ----------------------------------------------------------------!
   integer     , dimension(2*ncoh)            :: indx
   integer                                    :: ncoh2
   integer                                    :: il
   integer                                    :: i
   integer                                    :: j
   integer                                    :: ind
   real(kind=8), dimension(2*ncoh,2*ncoh)     :: masmatcp
   real(kind=8), dimension(2*ncoh,5)          :: mastermat
   real(kind=8), dimension(2*ncoh,2)          :: matal
   real(kind=8), dimension(2*ncoh)            :: mastervec_surf
   real(kind=8), dimension(2*ncoh)            :: mastervec_incid
   real(kind=8), dimension(ncoh+1)            :: explai
   real(kind=8), dimension(ncoh+1)            :: exmlai
   real(kind=8), dimension(ncoh+1)            :: downward_lw_incid
   real(kind=8), dimension(ncoh+1)            :: downward_lw_surf
   real(kind=8), dimension(ncoh+1)            :: upward_lw_incid
   real(kind=8), dimension(ncoh+1)            :: upward_lw_surf
   real(kind=8), dimension(ncoh)              :: source_lw
   real(kind=8), dimension(ncoh)              :: forcing_lw
   real(kind=8), dimension(ncoh)              :: A_dw
   real(kind=8), dimension(ncoh)              :: B_dw
   real(kind=8), dimension(ncoh)              :: C_dw
   real(kind=8), dimension(ncoh)              :: D_dw
   real(kind=8), dimension(ncoh)              :: A_uw
   real(kind=8), dimension(ncoh)              :: B_uw
   real(kind=8), dimension(ncoh)              :: C_uw
   real(kind=8), dimension(ncoh)              :: D_uw
   real(kind=8), dimension(ncoh)              :: E_uw
   real(kind=8), dimension(ncoh)              :: F_uw
   real(kind=8)                               :: emgs
   real(kind=8)                               :: T_grnd
   real(kind=8)                               :: zeta
   real(kind=8)                               :: eta
   real(kind=8)                               :: exk
   real(kind=8)                               :: zetai
   real(kind=8)                               :: exki
   real(kind=8)                               :: d
   real(kind=8)                               :: lw_v_surf_tmp
   real(kind=8)                               :: lw_v_incid_tmp
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external                     :: sngloff
   !---------------------------------------------------------------------------------------!


   !----- Convert some variables to double precision. -------------------------------------!
   emgs   = dble(semgs)
   t_grnd = dble(st_grnd)

   ncoh2  = 2*ncoh

   do il=1,ncoh
      zeta           = 2.0d0 * (1.0d0 - emis_v(pft(il))) / (3.0d0 * mubar)
      eta            = (2.0d0 + emis_v(pft(il)))/(3.0d0 * mubar)
      exk            = sqrt(eta*eta-zeta*zeta)
      exki           = 1.0d0 / exk
      zetai          = 1.0d0 / zeta
      source_lw(il)  = emis_v(pft(il)) * stefan8 *T_veg(il)**4
      forcing_lw(il) = -(zeta + eta) * source_lw(il)
      explai(il)     = exp( exk *TAI(il))
      exmlai(il)     = exp(-exk *TAI(il))


      !------------------------------------------------------------------------------------!
      !     Coefficient of lambda1 (and minus the coefficient of lambda2) for the bottom   !
      ! of a layer, downwelling radiation.                                                 !
      !------------------------------------------------------------------------------------!
      A_dw(il) = 5.d-1 * exki

      !----- Coefficient of lambda1, top of layer, downwelling radiation. -----------------!
      B_dw(il) = 5.d-1*exki*explai(il)

      !----- Coefficient of lambda2, top of layer, downwelling radiation. -----------------!
      C_dw(il) = -5.d-1*exki*exmlai(il)

      !----- Term of downwelling radiation not multiplying a lambda. ----------------------!
      D_dw(il) = 5.d-1*(exki**2)*forcing_lw(il) * (explai(il) + exmlai(il) - 2.0d0)

      A_uw(il) =   5.d-1 * zetai * (eta * exki - 1.0d0)
      B_uw(il) = - 5.d-1 * zetai * (eta * exki + 1.0d0)
      C_uw(il) = - source_lw(il) * zetai
      D_uw(il) = A_uw(il) * explai(il)
      E_uw(il) = B_uw(il) * exmlai(il)
      F_uw(il) = -source_lw(il) * zetai                                                    &
               + 5.d-1 * zetai * (eta*exki - 1.0d0) * explai(il)                           &
               * (forcing_lw(il) * exki * (1.0d0-exmlai(il)))                                 &
               - 5.d-1 * zetai * (eta*exki + 1.0d0) * exmlai(il)                           &
               * (forcing_lw(il) * exki * (explai(il)-1.0d0))
   end do


   do j=1,ncoh2
      do i=1,ncoh2
         masmatcp(i,j) = 0.0d0
      end do
      mastervec_surf(j)  = 0.0d0
      mastervec_incid(j) = 0.0d0
   end do

   !----- Vector is of the form: (lambda_N, lambda_{N-1},...,lambda_1). -------------------!

   masmatcp(1,1)      = B_dw(ncoh)
   masmatcp(1,2)      = C_dw(ncoh)
   mastervec_surf(1)  =-D_dw(ncoh)
   mastervec_incid(1) = 1.0d0

   do i=2,ncoh2-2,2
      ind = nint(real(ncoh2-i)*0.5)
      masmatcp(i,i-1)    = -A_dw(ind+1)
      masmatcp(i,i)      =  A_dw(ind+1)
      masmatcp(i,i+1)    =  B_dw(ind)
      masmatcp(i,i+2)    =  C_dw(ind)
      mastervec_surf(i)  = -D_dw(ind)
      mastervec_incid(i) = 0.0d0
   end do

   do i=3,ncoh2-1,2
      ind = nint(real(ncoh2-i+1)*0.5)
      masmatcp(i,i-2)    = -A_uw(ind+1)
      masmatcp(i,i-1)    = -B_uw(ind+1)
      masmatcp(i,i)      =  D_uw(ind)
      masmatcp(i,i+1)    =  E_uw(ind)
      mastervec_surf(i)  =  C_uw(ind+1) - F_uw(ind)
      mastervec_incid(i) =  0.0d0
   end do
   masmatcp(ncoh2,ncoh2-1) = A_uw(1) - (1.d0 - emgs) * A_dw(1)
   masmatcp(ncoh2,ncoh2)   = B_uw(1) + (1.d0 - emgs) * A_dw(1)
   mastervec_surf(ncoh2)   = emgs * stefan8 * T_grnd**4 - C_uw(1)
   mastervec_incid(ncoh2)  = 0.0d0

   mastermat(1,1) = 0.d0
   mastermat(1,2) = 0.d0
   mastermat(1,3) = masmatcp(1,1)
   mastermat(1,4) = masmatcp(1,2)
   mastermat(1,5) = 0.d0

   do i=2,ncoh2-2,2
      mastermat(i,1) = 0.d0
      mastermat(i,2) = masmatcp(i,i-1)
      mastermat(i,3) = masmatcp(i,i)
      mastermat(i,4) = masmatcp(i,i+1)
      mastermat(i,5) = masmatcp(i,i+2)
   end do

   do i=3,ncoh2-1,2
      mastermat(i,1) = masmatcp(i,i-2)
      mastermat(i,2) = masmatcp(i,i-1)
      mastermat(i,3) = masmatcp(i,i)
      mastermat(i,4) = masmatcp(i,i+1)
      mastermat(i,5) = 0.d0
   end do

   mastermat(ncoh2,1) = 0.d0
   mastermat(ncoh2,2) = masmatcp(ncoh2,ncoh2-1)
   mastermat(ncoh2,3) = masmatcp(ncoh2,ncoh2)
   mastermat(ncoh2,4) = 0.d0
   mastermat(ncoh2,5) = 0.d0
   
   !----- Invert matrix. ------------------------------------------------------------------!
   call bandec(mastermat,ncoh2,ncoh2,2,2,matal,indx,d)

   !----- Backsubstitute for contributions of ground and vegetation. ----------------------!
   call banbks(mastermat,ncoh2,ncoh2,2,2,matal,indx,mastervec_surf)

   !----- Backsubstitute for contribution of incident longwave at canopy top. -------------!
   call banbks(mastermat,ncoh2,ncoh2,2,2,matal,indx,mastervec_incid)

   do i=3,ncoh2-1,2
      ind = nint(real(ncoh2-i+1)*0.5)
      upward_lw_surf(ind+1)  = masmatcp(i,i) * mastervec_surf(i)                           &
                             + masmatcp(i,i+1) * mastervec_surf(i+1) + F_uw(ind)
      upward_lw_incid(ind+1) = masmatcp(i,i) * mastervec_incid(i)                          &
                             + masmatcp(i,i+1) * mastervec_incid(i+1)
   end do

   do i=2,ncoh2-2,2
      ind = nint(real(ncoh2-i)*0.5)
      downward_lw_surf(ind+1)  = masmatcp(i,i+1) * mastervec_surf(i+1)                     &
                               + masmatcp(i,i+2) * mastervec_surf(i+2) + D_dw(ind)
      downward_lw_incid(ind+1) = masmatcp(i,i+1) * mastervec_incid(i+1)                    &
                               + masmatcp(i,i+2) * mastervec_incid(i+2)
   end do

   upward_lw_surf(ncoh+1)    = D_uw(ncoh) * mastervec_surf(1)                              &
                             + E_uw(ncoh) * mastervec_surf(2) + F_uw(ncoh)
   upward_lw_incid(ncoh+1)   = D_uw(ncoh) * mastervec_incid(1)                             &
                             + E_uw(ncoh) * mastervec_incid(2)
   downward_lw_surf(ncoh+1)  = 0.0d0
   downward_lw_incid(ncoh+1) = 1.0d0
   downward_lw_surf(1)       = A_dw(1) * (mastervec_surf(ncoh2-1)  - mastervec_surf(ncoh2))
   downward_lw_incid(1)      = A_dw(1) * (mastervec_incid(ncoh2-1) - mastervec_incid(ncoh2))
   upward_lw_surf(1)         = (1.0d0-emgs) * downward_lw_surf(1)                          &
                             + emgs * stefan8 * T_grnd**4
   upward_lw_incid(1)        = (1.0d0-emgs) * downward_lw_incid(1)

   do il = 1,ncoh
      lw_v_surf_tmp  = downward_lw_surf(il+1)  - downward_lw_surf(il)                      &
                     + upward_lw_surf(il)  - upward_lw_surf(il+1)
      lw_v_incid_tmp = downward_lw_incid(il+1) - downward_lw_incid(il)                     &
                     + upward_lw_incid(il) - upward_lw_incid(il+1)
      lw_v_surf(il)  = sngloff(lw_v_surf_tmp , tiny_num8)
      lw_v_incid(il) = sngloff(lw_v_incid_tmp, tiny_num8)
   end do

   downward_lw_below_surf  = sngloff(downward_lw_surf(1)    , tiny_num8)
   downward_lw_below_incid = sngloff(downward_lw_incid(1)   , tiny_num8)
   upward_lw_below_surf    = sngloff(upward_lw_surf(1)      , tiny_num8)
   upward_lw_below_incid   = sngloff(upward_lw_incid(1)     , tiny_num8)
   upward_lw_above_surf    = sngloff(upward_lw_surf(ncoh+1) , tiny_num8)
   upward_lw_above_incid   = sngloff(upward_lw_incid(ncoh+1), tiny_num8)

   return
end subroutine lw_twostream
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will solve the long wave radiation within the canopy.  This solves   !
! the two stream approach considering the size distribution of cohorts, acknowledging the  !
! effect of leaves, and, if that's the user's will, the branches.                          !
!------------------------------------------------------------------------------------------!
subroutine lw_twostream_layer(ncoh,semgs,sT_grnd, pft,TAI,canopy_area,hgttop,hgtbot,T_veg  &
                             ,lw_v_surf,lw_v_incid,downward_lw_below_surf                  &
                             ,downward_lw_below_incid,upward_lw_below_surf                 &
                             ,upward_lw_below_incid,upward_lw_above_surf                   &
                             ,upward_lw_above_incid)
   use canopy_radiation_coms , only : emis_v            & ! intent(in)
                                    , mubar             ! ! intent(in)
   use pft_coms              , only : clumping_factor   ! ! intent(in)
   use consts_coms           , only : stefan8           & ! intent(in)
                                    , tiny_num8         ! ! intent(in)
   use rk4_coms              , only : tiny_offset       ! ! intent(in)
   use canopy_layer_coms     , only : ncanlyr           & ! intent(in)
                                    , ncanlyrp1         & ! intent(in)
                                    , ncanlyrt2         & ! intent(in)
                                    , zztop0i8          & ! intent(in)
                                    , ehgti8            & ! intent(in)
                                    , dzcan8            & ! intent(in)
                                    , zztop8            & ! intent(in)
                                    , zzbot8            & ! intent(in)
                                    , indx              & ! intent(out)
                                    , populated         & ! intent(out)
                                    , mastervec_surf    & ! intent(out)
                                    , mastervec_incid   & ! intent(out)
                                    , explai            & ! intent(out)
                                    , exmlai            & ! intent(out)
                                    , downward_lw_incid & ! intent(out)
                                    , downward_lw_surf  & ! intent(out)
                                    , upward_lw_incid   & ! intent(out)
                                    , upward_lw_surf    & ! intent(out)
                                    , source_lw         & ! intent(out)
                                    , forcing_lw        & ! intent(out)
                                    , A_dw              & ! intent(out)
                                    , B_dw              & ! intent(out)
                                    , C_dw              & ! intent(out)
                                    , D_dw              & ! intent(out)
                                    , A_uw              & ! intent(out)
                                    , B_uw              & ! intent(out)
                                    , C_uw              & ! intent(out)
                                    , D_uw              & ! intent(out)
                                    , E_uw              & ! intent(out)
                                    , F_uw              & ! intent(out)
                                    , lw_v_surf_layer   & ! intent(out)
                                    , lw_v_incid_layer  & ! intent(out)
                                    , matal             & ! intent(out)
                                    , mastermat         & ! intent(out)
                                    , masmatcp          & ! intent(out)
                                    , zero_canopy_layer ! ! subroutine

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                      , intent(in)  :: ncoh        ! # of cohorts
   real                         , intent(in)  :: semgs       ! Gnd. emmissivity
   real                         , intent(in)  :: st_grnd     ! Gnd. temperature
   integer     , dimension(ncoh), intent(in)  :: pft         ! Plant functional type
   real(kind=8), dimension(ncoh), intent(in)  :: TAI         ! "Tree" Area Index
   real(kind=8), dimension(ncoh), intent(in)  :: T_veg       ! Vegetation temperature
   real(kind=8), dimension(ncoh), intent(in)  :: canopy_area ! canopy area
   real(kind=8), dimension(ncoh), intent(in)  :: hgttop      ! "Tree" Area Index
   real(kind=8), dimension(ncoh), intent(in)  :: hgtbot      ! "Tree" Area Index
   real        , dimension(ncoh), intent(out) :: lw_v_surf
   real        , dimension(ncoh), intent(out) :: lw_v_incid
   real                         , intent(out) :: downward_lw_below_surf
   real                         , intent(out) :: upward_lw_above_surf
   real                         , intent(out) :: upward_lw_below_surf
   real                         , intent(out) :: downward_lw_below_incid
   real                         , intent(out) :: upward_lw_above_incid
   real                         , intent(out) :: upward_lw_below_incid
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: il
   integer                                    :: ico
   integer                                    :: ipft
   integer                                    :: i
   integer                                    :: j
   integer                                    :: ind
   integer                                    :: nactlyr
   integer                                    :: nactlyrp1
   integer                                    :: nactlyrt2
   integer     , dimension(ncoh)              :: kapartial
   integer     , dimension(ncoh)              :: kafull
   integer     , dimension(ncoh)              :: kzpartial
   integer     , dimension(ncoh)              :: kzfull
   integer     , dimension(ncoh)              :: nlyr_coh
   real(kind=8)                               :: emgs
   real(kind=8)                               :: T_grnd
   real(kind=8)                               :: zeta
   real(kind=8)                               :: eta
   real(kind=8)                               :: exk
   real(kind=8)                               :: zetai
   real(kind=8)                               :: exki
   real(kind=8)                               :: d
   real(kind=8)                               :: lw_v_surf_tmp
   real(kind=8)                               :: lw_v_incid_tmp
   real(kind=8)                               :: tad
   real(kind=8)                               :: this_tai
   !----- Variables that must be allocated every time. ------------------------------------!
   real(kind=8), dimension(:,:) , allocatable :: tai_lyr
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external                     :: sngloff
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These must be allocated and de-allocated every time as the number of cohorts is   !
   ! unbounded and may change.                                                             !
   !---------------------------------------------------------------------------------------!
   allocate(tai_lyr           (ncoh,ncanlyr))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Reset the temporary variables.                                                    !
   !---------------------------------------------------------------------------------------!
   call zero_canopy_layer('lw_twostream_layer')
   tai_lyr               (:,:) = 0.d0
   !---------------------------------------------------------------------------------------!



   !----- Convert some variables to double precision. -------------------------------------!
   emgs   = dble(semgs)
   t_grnd = dble(st_grnd)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the crown area for each layer, and determine the layer bounds for each       !
   ! cohort.                                                                               !
   !---------------------------------------------------------------------------------------!
   do ico = 1,ncoh
      !------ Find the layer bounds. ------------------------------------------------------!
      kapartial (ico) = min(ncanlyr,floor  ((hgtbot(ico) * zztop0i8)**ehgti8) + 1)
      kafull    (ico) = min(ncanlyr,ceiling((hgtbot(ico) * zztop0i8)**ehgti8) + 1)
      kzpartial (ico) = min(ncanlyr,ceiling((hgttop(ico) * zztop0i8)**ehgti8))
      kzfull    (ico) = min(ncanlyr,floor  ((hgttop(ico) * zztop0i8)**ehgti8))

      !------ Find the "tree" area density. -----------------------------------------------!
      tad             = tai(ico) / (hgttop(ico) - hgtbot(ico))


      !------------------------------------------------------------------------------------!
      !     Here we scale the leaf area density to the layer (partial or full), so the     !
      ! total leaf area index is preserved.  The only special case is when the top and     !
      ! bottom of the crown are in the same layer, in which case we simply use the leaf    !
      ! area index.  If the user is running with branch thermodynamics, then it should be  !
      ! total (leaf + branch) area index instead of leaf area index.                       !
      !------------------------------------------------------------------------------------!
      if (kapartial(ico) == kzpartial(ico)) then
         il                 = kapartial(ico)
         populated     (il) = .true.
         tai_lyr   (ico,il) = tai(ico)
      else 
         !------ Start with the fully vegetated layers. -----------------------------------!
         do il = kafull(ico),kzfull(ico)
            populated(il)         = .true.
            tai_lyr      (ico,il) = tad * dzcan8(il)
         end do

         !---------------------------------------------------------------------------------!
         !  Partial layer at the bottom of the crown.  This shouldn't be done if the       !
         ! partial and full indices are the same, otherwise we double account the layer.   !
         !---------------------------------------------------------------------------------!
         if (kapartial(ico) /= kafull(ico)) then
            il                    = kapartial(ico)
            populated        (il) = .true.
            tai_lyr      (ico,il) = tad * (zztop8(il)-hgtbot(ico))
         end if
         if (kzpartial(ico) /= kzfull(ico)) then
            il                    = kzpartial(ico)
            populated        (il) = .true.
            tai_lyr      (ico,il) = tad * (hgttop(ico) - zzbot8(il))
         end if
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Pack the information so we skip empty layers in the middle, thus avoinding         !
   ! singularities.  From this point on we must refer to nactlyr, not ncanlyr.             !
   !---------------------------------------------------------------------------------------!
   nactlyr   = count(populated)
   nactlyrp1 = nactlyr+1
   nactlyrt2 = 2 * nactlyr
   !----- Squeeze the arrays. -------------------------------------------------------------!
   do ico = 1,ncoh
      tai_lyr      (ico,1:nactlyr) = pack(tai_lyr      (ico,:), populated)
      if (nactlyrp1 < ncanlyr) then
         tai_lyr      (ico,nactlyrp1:ncanlyr) = 0.d0
      end if
      !----- Count the number of layers that have leaves/branches for this cohort. --------!
      nlyr_coh       (ico) = count(tai_lyr(ico,:) > 0.d0)
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   srclyrloop: do il=1,nactlyr
      !------------------------------------------------------------------------------------!
      !     Find the terms by averaging them using TAI as the weighting factor.  This is   !
      ! not the best way of doing it, but it is the simplest I could think, because it     !
      ! retains the fact that eta and zeta are constants within each layer.                !
      !------------------------------------------------------------------------------------!
      this_tai      = sum(tai_lyr(:,il))
      eta           = 0.d0
      zeta          = 0.d0
      source_lw(il) = 0.d0
      do ico = 1,ncoh
         ipft = pft(ico)
         eta           = eta  + tai_lyr(ico,il) * (2.d0 + emis_v(ipft)) / (3.d0 * mubar)
         zeta          = zeta + tai_lyr(ico,il)                                            &
                              * (2.d0 * (1.d0 - emis_v(ipft))) / (3.d0 * mubar)
         source_lw(il) = source_lw(il)                                                     &
                       + tai_lyr(ico,il) * emis_v(ipft) * stefan8 *T_veg(ico)**4
      end do
      !----- Normalise eta, zeta, and the LW source for this layer. -----------------------!
      eta           = eta           / this_tai
      zeta          = zeta          / this_tai
      source_lw(il) = source_lw(il) / this_tai


      !----- Find derived properties. -----------------------------------------------------!
      exk            = sqrt(eta*eta-zeta*zeta)
      exki           = 1.0d0 / exk
      zetai          = 1.0d0 / zeta
      forcing_lw(il) = -(zeta + eta) * source_lw(il)
      explai(il)     = exp( exk * this_tai)
      exmlai(il)     = exp(-exk * this_tai)


      !------------------------------------------------------------------------------------!
      !     Coefficient of lambda1 (and minus the coefficient of lambda2) for the bottom   !
      ! of a layer, downwelling radiation.                                                 !
      !------------------------------------------------------------------------------------!
      A_dw(il) = 5.0d-1 * exki

      !----- Coefficient of lambda1, top of layer, downwelling radiation. -----------------!
      B_dw(il) = 5.0d-1*exki*explai(il)

      !----- Coefficient of lambda2, top of layer, downwelling radiation. -----------------!
      C_dw(il) = -5.0d-1*exki*exmlai(il)

      !----- Term of downwelling radiation not multiplying a lambda. ----------------------!
      D_dw(il) = 5.d-1*(exki**2)*forcing_lw(il) * (explai(il) + exmlai(il) - 2.0d0)

      A_uw(il) =   5.d-1 * zetai * (eta * exki - 1.0d0)
      B_uw(il) = - 5.d-1 * zetai * (eta * exki + 1.0d0)
      C_uw(il) = - source_lw(il) * zetai
      D_uw(il) = A_uw(il) * explai(il)
      E_uw(il) = B_uw(il) * exmlai(il)
      F_uw(il) = -source_lw(il) * zetai                                                       &
               + 5.d-1 * zetai * (eta*exki - 1.0d0) * explai(il)                           &
               * (forcing_lw(il) * exki * (1.0d0-exmlai(il)))                                 &
               - 5.d-1 * zetai * (eta*exki + 1.0d0) * exmlai(il)                           &
               * (forcing_lw(il) * exki * (explai(il)-1.0d0))
   end do srclyrloop
   !---------------------------------------------------------------------------------------!



   !----- Initialise the matrix coefficients. ---------------------------------------------!
   do j=1,nactlyrt2
      do i=1,nactlyrt2
         masmatcp(i,j) = 0.0d0
      end do
      mastervec_surf(j)  = 0.0d0
      mastervec_incid(j) = 0.0d0
   end do
   !---------------------------------------------------------------------------------------!




   !----- Vector is of the form: (lambda_N, lambda_{N-1},...,lambda_1). -------------------!
   masmatcp(1,1)      = B_dw(nactlyr)
   masmatcp(1,2)      = C_dw(nactlyr)
   mastervec_surf(1)  =-D_dw(nactlyr)
   mastervec_incid(1) = 1.0d0

   do i=2,nactlyrt2-2,2
      ind = nint(real(nactlyrt2-i)*0.5)
      masmatcp(i,i-1)    = -A_dw(ind+1)
      masmatcp(i,i)      =  A_dw(ind+1)
      masmatcp(i,i+1)    =  B_dw(ind)
      masmatcp(i,i+2)    =  C_dw(ind)
      mastervec_surf(i)  = -D_dw(ind)
      mastervec_incid(i) = 0.0d0
   end do

   do i=3,nactlyrt2-1,2
      ind = nint(real(nactlyrt2-i+1)*0.5)
      masmatcp(i,i-2)    = -A_uw(ind+1)
      masmatcp(i,i-1)    = -B_uw(ind+1)
      masmatcp(i,i)      =  D_uw(ind)
      masmatcp(i,i+1)    =  E_uw(ind)
      mastervec_surf(i)  =  C_uw(ind+1) - F_uw(ind)
      mastervec_incid(i) =  0.0d0
   end do
   masmatcp(nactlyrt2,nactlyrt2-1) = A_uw(1) - (1.d0 - emgs) * A_dw(1)
   masmatcp(nactlyrt2,nactlyrt2)   = B_uw(1) + (1.d0 - emgs) * A_dw(1)
   mastervec_surf(nactlyrt2)       = emgs * stefan8 * T_grnd**4 - C_uw(1)
   mastervec_incid(nactlyrt2)      = 0.0d0

   mastermat(1,1) = 0.d0
   mastermat(1,2) = 0.d0
   mastermat(1,3) = masmatcp(1,1)
   mastermat(1,4) = masmatcp(1,2)
   mastermat(1,5) = 0.d0

   do i=2,nactlyrt2-2,2
      mastermat(i,1) = 0.d0
      mastermat(i,2) = masmatcp(i,i-1)
      mastermat(i,3) = masmatcp(i,i)
      mastermat(i,4) = masmatcp(i,i+1)
      mastermat(i,5) = masmatcp(i,i+2)
   end do

   do i=3,nactlyrt2-1,2
      mastermat(i,1) = masmatcp(i,i-2)
      mastermat(i,2) = masmatcp(i,i-1)
      mastermat(i,3) = masmatcp(i,i)
      mastermat(i,4) = masmatcp(i,i+1)
      mastermat(i,5) = 0.d0
   end do

   mastermat(nactlyrt2,1) = 0.d0
   mastermat(nactlyrt2,2) = masmatcp(nactlyrt2,nactlyrt2-1)
   mastermat(nactlyrt2,3) = masmatcp(nactlyrt2,nactlyrt2)
   mastermat(nactlyrt2,4) = 0.d0
   mastermat(nactlyrt2,5) = 0.d0
   
   !----- Invert matrix. ------------------------------------------------------------------!
   call bandec(mastermat,ncanlyrt2,nactlyrt2,2,2,matal,indx,d)

   !----- Backsubstitute for contributions of ground and vegetation. ----------------------!
   call banbks(mastermat,ncanlyrt2,nactlyrt2,2,2,matal,indx,mastervec_surf)

   !----- Backsubstitute for contribution of incident longwave at canopy top. -------------!
   call banbks(mastermat,ncanlyrt2,nactlyrt2,2,2,matal,indx,mastervec_incid)

   do i=3,nactlyrt2-1,2
      ind = nint(real(nactlyrt2-i+1)*0.5)
      upward_lw_surf(ind+1)  = masmatcp(i,i) * mastervec_surf(i)                                  &
                      + masmatcp(i,i+1) * mastervec_surf(i+1) + F_uw(ind)
      upward_lw_incid(ind+1) = masmatcp(i,i) * mastervec_incid(i)                                 &
                      + masmatcp(i,i+1) * mastervec_incid(i+1)
   end do

   do i=2,nactlyrt2-2,2
      ind = nint(real(nactlyrt2-i)*0.5)
      downward_lw_surf(ind+1)  = masmatcp(i,i+1) * mastervec_surf(i+1)                              &
                      + masmatcp(i,i+2) * mastervec_surf(i+2) + D_dw(ind)
      downward_lw_incid(ind+1) = masmatcp(i,i+1) * mastervec_incid(i+1)                             &
                      + masmatcp(i,i+2) * mastervec_incid(i+2)
   end do

   upward_lw_surf(nactlyr+1)    = D_uw(nactlyr) * mastervec_surf(1)                        &
                                + E_uw(nactlyr) * mastervec_surf(2) + F_uw(nactlyr)
   upward_lw_incid(nactlyr+1)   = D_uw(nactlyr) * mastervec_incid(1)                       &
                                + E_uw(nactlyr) * mastervec_incid(2)
   downward_lw_surf(nactlyr+1)  = 0.0d0
   downward_lw_incid(nactlyr+1) = 1.0d0
   downward_lw_surf(1)       = A_dw(1) * (mastervec_surf(nactlyrt2-1)  - mastervec_surf(nactlyrt2))
   downward_lw_incid(1)      = A_dw(1) * (mastervec_incid(nactlyrt2-1) - mastervec_incid(nactlyrt2))
   upward_lw_surf(1)       = (1.0d0-emgs) * downward_lw_surf(1) + emgs * stefan8 * T_grnd**4
   upward_lw_incid(1)      = (1.0d0-emgs) * downward_lw_incid(1)

   do il = 1,nactlyr
      lw_v_surf_tmp  = downward_lw_surf(il+1)  - downward_lw_surf(il)  + upward_lw_surf(il)  - upward_lw_surf(il+1)
      lw_v_incid_tmp = downward_lw_incid(il+1) - downward_lw_incid(il) + upward_lw_incid(il) - upward_lw_incid(il+1)
      lw_v_surf_layer(il)  = sngloff(lw_v_surf_tmp , tiny_num8)
      lw_v_incid_layer(il) = sngloff(lw_v_incid_tmp, tiny_num8)
   end do

   downward_lw_below_surf  = sngloff(downward_lw_surf(1)       , tiny_num8)
   downward_lw_below_incid = sngloff(downward_lw_incid(1)      , tiny_num8)
   upward_lw_below_surf    = sngloff(upward_lw_surf(1)         , tiny_num8)
   upward_lw_below_incid   = sngloff(upward_lw_incid(1)        , tiny_num8)
   upward_lw_above_surf    = sngloff(upward_lw_surf(nactlyr+1) , tiny_num8)
   upward_lw_above_incid   = sngloff(upward_lw_incid(nactlyr+1), tiny_num8)

   !---------------------------------------------------------------------------------------!
   !     Integrate the total amount of energy for each cohort.                             !
   !---------------------------------------------------------------------------------------!
   lw_v_surf (:) = 0.0
   lw_v_incid(:) = 0.0
   do il = 1, nactlyr
      this_tai = sum(tai_lyr(:,il))
      do ico = 1, ncoh
         lw_v_surf (ico) = lw_v_surf(ico) + lw_v_surf_layer(il)                            &
                         * sngloff(tai_lyr(ico,il) / this_tai, tiny_num8)
         lw_v_incid(ico) = lw_v_incid(ico) + lw_v_incid_layer(il)                          &
                         * sngloff(tai_lyr(ico,il) / this_tai, tiny_num8)
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !------ De-allocate the temporary structure. -------------------------------------------!
   deallocate(tai_lyr    )
   !---------------------------------------------------------------------------------------!

   return
end subroutine lw_twostream_layer
!==========================================================================================!
!==========================================================================================!
