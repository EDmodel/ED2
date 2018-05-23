!==========================================================================================!
!==========================================================================================!
!    This subroutine will solve the within canopy two-stream radiation for shortwave,      !
! considering the light absorption by leaves, and branches if the user wants so.  In both  !
! cases we consider the clumpiness effect (i.e., the fact that leaves and branches aren't  !
! randomly distributed).  The crown area will be also considered when the user wants it.   !
!------------------------------------------------------------------------------------------!
subroutine old_sw_two_stream (salbedo_par,salbedo_nir,scosaoi,ncoh,pft,lai,wai,canopy_area &
                             ,radprof_flip,par_beam_flip,par_diff_flip,sw_abs_beam_flip    &
                             ,sw_abs_diff_flip,dw_parlo_beam,dw_parlo_diff                 &
                             ,uw_parhi_diff,dw_nirlo_beam,dw_nirlo_diff,uw_nirhi_diff      &
                             ,par_level_beam,par_level_diffd,par_level_diffu               &
                             ,light_level,light_beam_level                                 &
                             ,light_diff_level)

   use ed_max_dims          , only : n_pft                   & ! intent(in)
                                   , n_radprof               ! ! intent(in)
   use rk4_coms             , only : tiny_offset             ! ! intent(in)
   use canopy_radiation_coms, only : leaf_backscatter_nir    & ! intent(in)
                                   , leaf_backscatter_vis    & ! intent(in)
                                   , leaf_scatter_nir        & ! intent(in)
                                   , leaf_scatter_vis        & ! intent(in)
                                   , wood_backscatter_nir    & ! intent(in)
                                   , wood_backscatter_vis    & ! intent(in)
                                   , wood_scatter_nir        & ! intent(in)
                                   , wood_scatter_vis        & ! intent(in)
                                   , clumping_factor         & ! intent(in)
                                   , phi1                    & ! intent(in)
                                   , phi2                    & ! intent(in)
                                   , mu_bar                  & ! intent(in)
                                   , par_beam_norm           & ! intent(in)
                                   , par_diff_norm           & ! intent(in)
                                   , nir_beam_norm           & ! intent(in)
                                   , nir_diff_norm           & ! intent(in)
                                   , cosz_min8               ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, dimension(ncoh)       , intent(in)    :: pft
   integer                        , intent(in)    :: ncoh
   real(kind=8), dimension(ncoh)  , intent(in)    :: lai
   real(kind=8), dimension(ncoh)  , intent(in)    :: wai
   real(kind=8), dimension(ncoh)  , intent(in)    :: canopy_area
   real                           , intent(in)    :: salbedo_par
   real                           , intent(in)    :: salbedo_nir
   real                           , intent(in)    :: scosaoi
   real, dimension(n_radprof,ncoh), intent(inout) :: radprof_flip
   real, dimension(ncoh)          , intent(out)   :: par_beam_flip
   real, dimension(ncoh)          , intent(out)   :: par_diff_flip
   real, dimension(ncoh)          , intent(out)   :: sw_abs_beam_flip
   real, dimension(ncoh)          , intent(out)   :: sw_abs_diff_flip
   real                           , intent(out)   :: uw_parhi_diff
   real                           , intent(out)   :: uw_nirhi_diff
   real                           , intent(out)   :: dw_parlo_beam
   real                           , intent(out)   :: dw_parlo_diff
   real                           , intent(out)   :: dw_nirlo_beam
   real                           , intent(out)   :: dw_nirlo_diff
   real(kind=8), dimension(ncoh)  , intent(out)   :: par_level_beam
   real(kind=8), dimension(ncoh)  , intent(out)   :: par_level_diffu
   real(kind=8), dimension(ncoh)  , intent(out)   :: par_level_diffd
   real(kind=8), dimension(ncoh)  , intent(out)   :: light_level
   real(kind=8), dimension(ncoh)  , intent(out)   :: light_beam_level
   real(kind=8), dimension(ncoh)  , intent(out)   :: light_diff_level
   !----- Local variables -----------------------------------------------------------------!
   integer     , dimension(2*ncoh)                :: indx
   integer                                        :: il
   integer                                        :: ipft
   integer                                        :: ncoh2
   integer                                        :: iband
   integer                                        :: i
   integer                                        :: j
   integer                                        :: ind
   real(kind=8), dimension(ncoh)                  :: expkl_top
   real(kind=8), dimension(ncoh)                  :: expkl_bot
   real(kind=8), dimension(ncoh)                  :: expamk_top
   real(kind=8), dimension(ncoh)                  :: expamk_bot
   real(kind=8), dimension(ncoh)                  :: expapk_top
   real(kind=8), dimension(ncoh)                  :: expapk_bot
   real(kind=8), dimension(ncoh)                  :: a_top
   real(kind=8), dimension(ncoh)                  :: a_bot
   real(kind=8), dimension(ncoh)                  :: b_top
   real(kind=8), dimension(ncoh)                  :: b_bot
   real(kind=8), dimension(ncoh)                  :: c_top
   real(kind=8), dimension(ncoh)                  :: c_bot
   real(kind=8), dimension(ncoh)                  :: f_top
   real(kind=8), dimension(ncoh)                  :: f_bot
   real(kind=8), dimension(ncoh)                  :: g_top
   real(kind=8), dimension(ncoh)                  :: g_bot
   real(kind=8), dimension(ncoh)                  :: h_top
   real(kind=8), dimension(ncoh)                  :: h_bot
   real(kind=8), dimension(ncoh)                  :: beam_bot
   real(kind=8), dimension(ncoh)                  :: beam_bot_crown
   real(kind=8), dimension(ncoh)                  :: tai
   real(kind=8), dimension(ncoh)                  :: eff_tai
   real(kind=8), dimension(ncoh)                  :: lambda
   real(kind=8), dimension(ncoh)                  :: beam_backscatter
   real(kind=8), dimension(ncoh)                  :: cohort_scatter_vis
   real(kind=8), dimension(ncoh)                  :: cohort_backscatter_vis
   real(kind=8), dimension(ncoh)                  :: cohort_scatter_nir
   real(kind=8), dimension(ncoh)                  :: cohort_backscatter_nir
   real(kind=8), dimension(ncoh)                  :: cohort_scatter
   real(kind=8), dimension(ncoh)                  :: cohort_backscatter
   real(kind=8), dimension(ncoh)                  :: cohort_clumping
   real(kind=8), dimension(ncoh+1)                :: upward_par_beam
   real(kind=8), dimension(ncoh+1)                :: upward_par_diff
   real(kind=8), dimension(ncoh+1)                :: upward_nir_beam
   real(kind=8), dimension(ncoh+1)                :: upward_nir_diff
   real(kind=8), dimension(ncoh+1)                :: downward_nir_beam
   real(kind=8), dimension(ncoh+1)                :: downward_nir_diff
   real(kind=8), dimension(ncoh+1)                :: downward_par_beam
   real(kind=8), dimension(ncoh+1)                :: downward_par_diff
   real(kind=8), dimension(2*ncoh)                :: mastervec_beam
   real(kind=8), dimension(2*ncoh)                :: masveccp_beam
   real(kind=8), dimension(2*ncoh)                :: mastervec_diff
   real(kind=8), dimension(2*ncoh)                :: masveccp_diff
   real(kind=8), dimension(2*ncoh,2)              :: matal
   real(kind=8), dimension(2*ncoh,5)              :: mastermat
   real(kind=8), dimension(2*ncoh,2*ncoh)         :: masmatcp
   real(kind=8)                                   :: albedo
   real(kind=8)                                   :: cosaoi
   real(kind=8)                                   :: eta
   real(kind=8)                                   :: zeta
   real(kind=8)                                   :: iota
   real(kind=8)                                   :: exk
   real(kind=8)                                   :: exki
   real(kind=8)                                   :: zetai
   real(kind=8)                                   :: d
   real(kind=8)                                   :: rhoo
   real(kind=8)                                   :: psi
   real(kind=8)                                   :: sigma
   real(kind=8)                                   :: source_bot
   real(kind=8)                                   :: source_top
   real(kind=8)                                   :: beam_top
   real(kind=8)                                   :: diff_top
   real(kind=8)                                   :: weight_leaf
   real(kind=8)                                   :: weight_wood
   real(kind=8)                                   :: proj_area
   real(kind=8)                                   :: snglscat_alb
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)                  , external       :: sngloff
   !---------------------------------------------------------------------------------------!

   
   !----- Convert input variable to double precision. -------------------------------------!
   cosaoi = max(cosz_min8,dble(scosaoi))

   !----- Calculate factors common for NIR, PAR. ------------------------------------------!
   ncoh2      = 2*ncoh
   lambda     = 5.d-1/cosaoi
   do il=1,ncoh
      ipft = pft(il)

      !------------------------------------------------------------------------------------!
      !     We find the optical depth of the direct beam (lambda), following CLM technical !
      ! note (equation 3.3 and text after equation 3.2).                                   !
      !------------------------------------------------------------------------------------!
      proj_area      = phi1(ipft) + phi2(ipft) * cosaoi
      lambda(il)     = proj_area / cosaoi
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the backscatter for direct beam radiation, following equations (3.14)    !
      ! and (3.15) of the CLM technical note.  The forward scattering was omitted on both  !
      ! equations because they are different for PAR and NIR, but in both cases they       !
      ! cancel out.                                                                        !
      !------------------------------------------------------------------------------------!
      snglscat_alb         = 5.d-1 * proj_area / (phi2(ipft) * cosaoi + proj_area)         &
                           * ( 1.d0 - phi1(ipft) * cosaoi                                  &
                             / (phi2(ipft) * cosaoi + proj_area)                           &
                             * log (1.d0  + (phi2(ipft) * cosaoi + proj_area)              &
                                          / (phi1(ipft) * cosaoi)))
      beam_backscatter(il) = ( 1.d0 + mu_bar(ipft) * lambda(il) ) * snglscat_alb           &
                           / ( mu_bar(ipft) * lambda(il) )
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the effective TAI for this layer, and the net clumping factor.            !
      !------------------------------------------------------------------------------------!
      tai                   (il) = lai(il) + wai(il)
      eff_tai               (il) = clumping_factor(ipft) * lai(il) + wai(il)
      cohort_clumping       (il) = eff_tai(il) / tai(il)
      !------------------------------------------------------------------------------------!




      !----- Find the scattering coefficients. --------------------------------------------!
      weight_leaf                = clumping_factor(ipft) * lai(il) / eff_tai(il)
      weight_wood                = 1.d0 - weight_leaf
      cohort_scatter_vis    (il) = leaf_scatter_vis(ipft)     * weight_leaf                &
                                 + wood_scatter_vis(ipft)     * weight_wood
      cohort_scatter_nir    (il) = leaf_scatter_nir(ipft)     * weight_leaf                &
                                 + wood_scatter_nir(ipft)     * weight_wood
      cohort_backscatter_vis(il) = leaf_backscatter_vis(ipft) * weight_leaf                &
                                 + wood_backscatter_vis(ipft) * weight_wood
      cohort_backscatter_nir(il) = leaf_backscatter_nir(ipft) * weight_leaf                &
                                 + wood_backscatter_nir(ipft) * weight_wood
      
   end do

   !----- Loop over bands (currently Visible and near infrared). --------------------------!
   bandloop: do iband = 1,2
      select case(iband)
      case (1) !----- Visible (or PAR). ---------------------------------------------------!
         do il = 1,ncoh
            cohort_scatter    (il) = cohort_scatter_vis    (il)
            cohort_backscatter(il) = cohort_backscatter_vis(il)
         end do
         beam_top = par_beam_norm
         diff_top = par_diff_norm
         albedo   = dble(salbedo_par)
      case (2) !----- Near infrared (or NIR). ---------------------------------------------!
         do il = 1,ncoh
            cohort_scatter    (il) = cohort_scatter_nir    (il)
            cohort_backscatter(il) = cohort_backscatter_nir(il)
         end do
         beam_top = nir_beam_norm
         diff_top = nir_diff_norm
         albedo   = dble(salbedo_nir)
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
      beam_bot_crown(ncoh)  = beam_top * exp( - lambda(ncoh) * eff_tai(ncoh)               &
                                              / canopy_area(ncoh) )
!!      par_level_beam(ncoh)  = beam_top * exp( - 5.d-1 * lambda(ncoh) * eff_tai(ncoh)       &
!!                                              / canopy_area(ncoh) )
      beam_bot(ncoh)        = beam_top * (1.d0 - canopy_area(ncoh))                        &
                            + canopy_area(ncoh) * beam_bot_crown(ncoh)
      do il=ncoh-1,1,-1
         beam_bot_crown(il) = beam_bot(il+1) * exp( - lambda(il) * eff_tai(il)             &
                                                    / canopy_area(il) )
         beam_bot(il)       = beam_bot(il+1) * (1.d0 - canopy_area(il))                    & 
                            + canopy_area(il) * beam_bot_crown(il)
!!         par_beam_level(il) = beam_bot(il+1)                                               &
 !!                           * exp(- 5.d-1 * lambda(il) * eff_tai(il) / canopy_area(il))    &
 !!                           * canopy_area(il)                                              &
 !!                           + (1.d0-canopy_area(il)) * beam_bot(il+1)
      end do

      do il=1,ncoh
         ipft  = pft(il)
         eta   = cohort_clumping(il)                                                       &
               * (1.0d0 - ( 1.0d0 - cohort_backscatter(il)) * cohort_scatter(il))          &
               / mu_bar(ipft)
         zeta  = cohort_scatter(il) * cohort_backscatter(il) * cohort_clumping(il)         &
               / mu_bar(ipft)
         psi   = cohort_clumping(il) * lambda(il)
         iota  = psi * cohort_scatter(il) * beam_backscatter(il)


         !----- Find derived properties. --------------------------------------------------!
         exk   = sqrt(eta*eta - zeta*zeta)
         exki  = 1.0d0/exk
         zetai = 1.0d0/zeta

         sigma = psi * psi + zeta * zeta - eta * eta

         !----- Sources. ------------------------------------------------------------------!
         source_bot = iota * beam_bot_crown(il)
         source_top = source_bot * exp(lambda(il) * eff_tai(il))

         !----- Forcing coefficients. -----------------------------------------------------!
         rhoo  = - (zeta + eta + psi) * source_bot

         !----- Calculate exponentials only once. -----------------------------------------!
         expkl_bot (il) = 1.0d0
         expkl_top (il) = exp(exk * tai(il))
         expamk_bot(il) = 1.0d0
         expamk_top(il) = exp((psi-exk) * tai(il))
         expapk_bot(il) = 1.0d0
         expapk_top(il) = exp((psi+exk) * tai(il))
         a_bot(il)      = -source_bot*zetai
         a_top(il)      = -source_top*zetai                                                &
                        + 5.d-1*zetai*(eta*exki-1.0d0)*expkl_top(il)*rhoo/(psi-exk)        &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 5.d-1*zetai*(eta*exki+1.0d0) / expkl_top(il)*rhoo/(psi+exk)      &
                        * (expapk_top(il)-expapk_bot(il))
         b_bot(il)      = 5.d-1*zetai*(eta*exki-1.0d0)
         b_top(il)      = 5.d-1*zetai*(eta*exki-1.0d0)*expkl_top(il)
         c_bot(il)      = -5.d-1*zetai*(eta*exki+1.0d0)
         c_top(il)      = -5.d-1*zetai*(eta*exki+1.0d0)/expkl_top(il)
         f_bot(il)      = 0.0d0
         f_top(il)      = 5.d-1*exki*expkl_top(il)*rhoo/(psi-exk)                          &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 5.d-1*exki/expkl_top(il)*rhoo/(psi+exk)                          &
                        * (expapk_top(il)-expapk_bot(il))
         g_bot(il)      = 5.d-1*exki
         g_top(il)      = 5.d-1*exki*expkl_top(il)
         h_bot(il)      = -5.d-1*exki
         h_top(il)      = -5.d-1*exki/expkl_top(il)
      end do

      
      !----- Organize the matrix coefficients. --------------------------------------------!
      do j=1,ncoh2
         do i=1,ncoh2
            masmatcp(i,j)=0.0d0
         end do
         mastervec_beam(j)    = 0.0d0
         masveccp_beam(j)     = 0.0d0
         mastervec_diff(j) = 0.0d0
         masveccp_diff(j)  = 0.0d0
      end do

      masmatcp(1,1)        = g_top(ncoh)
      masmatcp(1,2)        = h_top(ncoh)
      mastervec_beam(1)    = -f_top(ncoh)
      mastervec_diff(1) = diff_top
      masveccp_beam(1)     = mastervec_beam(1)
      masveccp_diff(1)  = mastervec_diff(1)

      do i=2,ncoh2-2,2
         masmatcp(i,i-1) = g_bot(nint(real(ncoh2-i+2)*0.5))
         masmatcp(i,i)   = h_bot(nint(real(ncoh2-i+2)*0.5))

         masmatcp(i,i+1)      = -g_top(nint(real(ncoh2-i)*0.5))
         masmatcp(i,i+2)      = -h_top(nint(real(ncoh2-i)*0.5))
         mastervec_beam(i)    = -f_bot(nint(real(ncoh2-i+2)*0.5))                          &
                              + f_top(nint(real(ncoh2-i)*0.5))
         mastervec_diff(i) = 0.0d0
         masveccp_beam(i)     = mastervec_beam(i)
         masveccp_diff(i)  = mastervec_diff(i)
      end do

      do i=3,ncoh2-1,2
         masmatcp(i,i-2)      = b_bot(nint(real(ncoh2-i+3)*0.5))
         masmatcp(i,i-1)      = c_bot(nint(real(ncoh2-i+3)*0.5))
         masmatcp(i,i)        = -b_top(nint(real(ncoh2-i+1)*0.5))
         masmatcp(i,i+1)      = -c_top(nint(real(ncoh2-i+1)*0.5))
         mastervec_beam(i)    = -a_bot(nint(real(ncoh2-i+3)*0.5))                          &
                              + a_top(nint(real(ncoh2-i+1)*0.5))
         masveccp_beam(i)     = mastervec_beam(i)
         mastervec_diff(i) = 0.0d0
         masveccp_diff(i)  = mastervec_diff(i)
      end do
      masmatcp(ncoh2,ncoh2-1)  =  b_bot(1) - albedo * g_bot(1)
      masmatcp(ncoh2,ncoh2)    =  c_bot(1) - albedo * h_bot(1)
      mastervec_beam(ncoh2)    = -a_bot(1) + albedo * beam_bot(1)
      masveccp_beam(ncoh2)     = mastervec_beam(ncoh2)
      mastervec_diff(ncoh2) = 0.0d0
      masveccp_diff(ncoh2)  = mastervec_diff(ncoh2)
      
      !----- Prepare for inversion. -------------------------------------------------------!
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
      call banbks(mastermat,ncoh2,ncoh2,2,2,matal,indx,mastervec_diff)
      
      !----- Improve the solution. --------------------------------------------------------!
      call mprove(masmatcp,mastermat,matal,ncoh2,ncoh2,5,2,indx                            &
                 ,masveccp_beam,mastervec_beam)
      call mprove(masmatcp,mastermat,matal,ncoh2,ncoh2,5,2,indx                            &
                 ,masveccp_diff,mastervec_diff)

      select case (iband)
      case (1) 
         !---- Visible (or PAR) band. -----------------------------------------------------!
         do i=3,ncoh2-1,2
            ind                     = nint(real(ncoh2-i+1)*0.5)+1
            upward_par_beam(ind)    = a_bot(ind) + b_bot(ind) * mastervec_beam(i-2)        &
                                    + c_bot(ind) * mastervec_beam(i-1)
            upward_par_diff(ind)    = b_bot(ind) * mastervec_diff(i-2)                     &
                                    + c_bot(ind) * mastervec_diff(i-1)
         end do

         do i=2,ncoh2-2,2
            ind                       = nint(real(ncoh2-i)*0.5)+1
            downward_par_beam(ind)    = beam_bot(ind) + f_bot(ind)                         &
                                      + h_bot(ind) * mastervec_beam(i)                     &
                                      + g_bot(ind) * mastervec_beam(i-1)
            downward_par_diff(ind)    = h_bot(ind) * mastervec_diff(i)                     &
                                      + g_bot(ind) * mastervec_diff(i-1)
         end do

         upward_par_beam(ncoh+1)      = b_top(ncoh) * mastervec_beam(1)                    &
                                      + c_top(ncoh) * mastervec_beam(2) + a_top(ncoh)
         upward_par_diff(ncoh+1)      = b_top(ncoh) * mastervec_diff(1)                    &
                                      + c_top(ncoh) * mastervec_diff(2)
         downward_par_beam(ncoh+1)    = beam_top
         downward_par_diff(ncoh+1) = diff_top
         downward_par_beam(1)         = g_bot(1) * mastervec_beam(ncoh2-1)                 &
                                      + h_bot(1) * mastervec_beam(ncoh2) + f_bot(1)        &
                                      + beam_bot(1)
         downward_par_diff(1)         = g_bot(1) * mastervec_diff(ncoh2-1)                 &
                                      + h_bot(1) * mastervec_diff(ncoh2)
         upward_par_beam(1)           = albedo * downward_par_beam(1)
         upward_par_diff(1)           = albedo * downward_par_diff(1)
      case (2)

         !---- Near infra-red. ------------------------------------------------------------!
         do i=3,ncoh2-1,2
            ind                     = nint(real(ncoh2-i+1)*0.5)+1
            upward_nir_beam(ind)    = a_bot(ind) + b_bot(ind) * mastervec_beam(i-2)        &
                                    + c_bot(ind) * mastervec_beam(i-1)
            upward_nir_diff(ind)    = b_bot(ind) * mastervec_diff(i-2)                     &
                                    + c_bot(ind) * mastervec_diff(i-1)
         end do

         do i=2,ncoh2-2,2
            ind                       = nint(real(ncoh2-i)*0.5)+1
            downward_nir_beam(ind)    = beam_bot(ind) + f_bot(ind)                         &
                                      + h_bot(ind) * mastervec_beam(i)                     &
                                      + g_bot(ind) * mastervec_beam(i-1)
            downward_nir_diff(ind)    = h_bot(ind) * mastervec_diff(i)                     &
                                      + g_bot(ind) * mastervec_diff(i-1)
         end do

         upward_nir_beam(ncoh+1)      = b_top(ncoh) * mastervec_beam(1)                    &
                                      + c_top(ncoh) * mastervec_beam(2) + a_top(ncoh)
         upward_nir_diff(ncoh+1)      = b_top(ncoh) * mastervec_diff(1)                    &
                                      + c_top(ncoh) * mastervec_diff(2)
         downward_nir_beam(ncoh+1)    = beam_top
         downward_nir_diff(ncoh+1)    = diff_top
         downward_nir_beam(1)         = g_bot(1) * mastervec_beam(ncoh2-1)                 &
                                      + h_bot(1) * mastervec_beam(ncoh2) + f_bot(1)        &
                                      + beam_bot(1)
         downward_nir_diff(1)         = g_bot(1) * mastervec_diff(ncoh2-1)                 &
                                      + h_bot(1) * mastervec_diff(ncoh2)
         upward_nir_beam(1)           = albedo * downward_nir_beam(1)
         upward_nir_diff(1)           = albedo * downward_nir_diff(1)
      end select
   end do bandloop
   
   do il=1,ncoh
      !------------------------------------------------------------------------------------!
      !     Find the light level of the visible band, which is the one that the plants     !
      ! use, and the total light level, which includes the near-infrared.                  !
      !------------------------------------------------------------------------------------!
      par_level_diffd(il)          = 5.d-1 * ( downward_par_diff(il  )                     &
                                 + downward_par_diff(il+1) ) / (par_diff_norm+par_beam_norm)

      par_level_diffu(il)          = 5.d-1 * ( upward_par_diff(il  )                     &
                                 + upward_par_diff(il+1) ) / (par_diff_norm+par_beam_norm)

      par_level_beam(il)          = 5.d-1 * ( downward_par_beam   (il  )                   &
                                        + downward_par_beam   (il+1) ) / (par_beam_norm+par_diff_norm)

      light_level(il)         = 5.d-1 * ( downward_par_diff(il  )                          &
                                        + downward_par_beam   (il  )                       &
                                        + downward_nir_diff(il  )                          &
                                        + downward_nir_beam   (il  )                       &
                                        + downward_par_diff(il+1)                          &
                                        + downward_par_beam   (il+1)                       &
                                        + downward_nir_diff(il+1)                          &
                                        + downward_nir_beam   (il+1) )

      light_beam_level(il)    = 5.d-1 * ( downward_par_beam   (il  )                       &
                                        + downward_nir_beam   (il  )                       &
                                        + downward_par_beam   (il+1)                       &
                                        + downward_nir_beam   (il+1) )                     &
                              / ( par_beam_norm + nir_beam_norm )

      light_diff_level(il)    = 5.d-1 * ( downward_par_diff(il  )                          &
                                        + downward_nir_diff(il  )                          &
                                        + downward_par_diff(il+1)                          &
                                        + downward_nir_diff(il+1) )                        &
                              / ( par_diff_norm + nir_beam_norm )
      !------------------------------------------------------------------------------------!

      par_beam_flip(il)    = sngloff(  downward_par_beam       (il+1)                      &
                                     - downward_par_beam       (il  )                      &
                                     + upward_par_beam         (il  )                      &
                                     - upward_par_beam         (il+1)                      &
                                     , tiny_offset )
      par_diff_flip(il)    = sngloff(  downward_par_diff (il+1)                            &
                                     - downward_par_diff (il  )                            &
                                     + upward_par_diff   (il  )                            &
                                     - upward_par_diff   (il+1)                            &
                                     , tiny_offset )

      radprof_flip(1,il)   = sngloff(  downward_par_beam (il), tiny_offset)
      radprof_flip(2,il)   = sngloff(  upward_par_beam   (il), tiny_offset)
      radprof_flip(3,il)   = sngloff(  downward_par_diff (il), tiny_offset)
      radprof_flip(4,il)   = sngloff(  upward_par_diff   (il), tiny_offset)
      radprof_flip(5,il)   = sngloff(  downward_nir_beam (il), tiny_offset)
      radprof_flip(6,il)   = sngloff(  upward_nir_beam   (il), tiny_offset)
      radprof_flip(7,il)   = sngloff(  downward_nir_diff (il), tiny_offset)
      radprof_flip(8,il)   = sngloff(  upward_nir_diff   (il), tiny_offset)


      sw_abs_beam_flip(il) = par_beam_flip(il) + sngloff(  downward_nir_beam   (il+1)      &
                                                         - downward_nir_beam   (il  )      &
                                                         + upward_nir_beam     (il  )      &
                                                         - upward_nir_beam     (il+1)      &
                                                         , tiny_offset )

      sw_abs_diff_flip(il) = par_diff_flip(il) + sngloff(  downward_nir_diff(il+1)         &
                                                         - downward_nir_diff(il  )         &
                                                         + upward_nir_diff  (il  )         &
                                                         - upward_nir_diff  (il+1)         &
                                                         , tiny_offset )

      !----- Ensure that we don't get any negative radiation... ---------------------------!
      par_beam_flip   (il)  = max(0.0,par_beam_flip   (il) )
      par_diff_flip   (il)  = max(0.0,par_diff_flip   (il) )
      sw_abs_beam_flip(il)  = max(0.0,sw_abs_beam_flip(il) )
      sw_abs_diff_flip(il)  = max(0.0,sw_abs_diff_flip(il) )
   end do
   !---------------------------------------------------------------------------------------!





   !----- Copy to the output variables. ---------------------------------------------------!
   dw_parlo_beam    = max(0.0,sngl(downward_par_beam        (1) ) )
   dw_parlo_diff    = max(0.0,sngl(downward_par_diff        (1) ) )
   uw_parhi_diff    = max(0.0,sngl(upward_par_beam(ncoh+1) + upward_par_diff(ncoh+1) ) )
   dw_nirlo_beam    = max(0.0,sngl(downward_nir_beam        (1) ) )
   dw_nirlo_diff    = max(0.0,sngl(downward_nir_diff        (1) ) )
   uw_nirhi_diff    = max(0.0,sngl(upward_nir_beam(ncoh+1) + upward_nir_diff(ncoh+1) ) )
   !---------------------------------------------------------------------------------------!

   return
end subroutine old_sw_two_stream
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
subroutine old_lw_two_stream(semgs,sT_grnd,rlong_top4,ncoh, pft,lai,wai,leaf_temp          &
                            ,wood_temp,radprof_flip,tir_flip,dw_tirlo,uw_tirlo,uw_tirhi)
   use ed_max_dims           , only : n_radprof            ! ! intent(in)
   use canopy_radiation_coms , only : leaf_emiss_tir       & ! intent(in)
                                    , wood_emiss_tir       & ! intent(in)
                                    , leaf_backscatter_tir & ! intent(in)
                                    , wood_backscatter_tir & ! intent(in)
                                    , mu_bar               ! ! intent(in)
   use consts_coms           , only : stefan8              & ! intent(in)
                                    , tiny_num8            ! ! intent(in)
   use rk4_coms              , only : tiny_offset          ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                , intent(in)    :: ncoh        ! # of cohorts
   real                                   , intent(in)    :: semgs       ! Gnd. emmissivity
   real                                   , intent(in)    :: st_grnd     ! Gnd. temperature
   real(kind=4)                           , intent(in)    :: rlong_top4
   integer     , dimension(ncoh)          , intent(in)    :: pft         ! Plant func. type
   real(kind=8), dimension(ncoh)          , intent(in)    :: lai         ! Leaf Area Index
   real(kind=8), dimension(ncoh)          , intent(in)    :: wai         ! Wood Area Index
   real(kind=8), dimension(ncoh)          , intent(in)    :: leaf_temp   ! Leaf temperature
   real(kind=8), dimension(ncoh)          , intent(in)    :: wood_temp   ! Leaf temperature
   real(kind=4), dimension(n_radprof,ncoh), intent(inout) :: radprof_flip
   real(kind=4), dimension(ncoh)          , intent(out)   :: tir_flip
   real(kind=4)                           , intent(out)   :: dw_tirlo
   real(kind=4)                           , intent(out)   :: uw_tirlo
   real(kind=4)                           , intent(out)   :: uw_tirhi
   !----- Local variables. ----------------------------------------------------------------!
   integer     , dimension(2*ncoh)                        :: indx
   integer                                                :: ncoh2
   integer                                                :: il
   integer                                                :: i
   integer                                                :: j
   integer                                                :: ind
   integer                                                :: ipft
   real(kind=8), dimension(2*ncoh,2*ncoh)                 :: masmatcp
   real(kind=8), dimension(2*ncoh,5)                      :: mastermat
   real(kind=8), dimension(2*ncoh,2)                      :: matal
   real(kind=8), dimension(2*ncoh)                        :: mastervec_surf
   real(kind=8), dimension(2*ncoh)                        :: mastervec_incid
   real(kind=8), dimension(ncoh+1)                        :: explai
   real(kind=8), dimension(ncoh+1)                        :: exmlai
   real(kind=8), dimension(ncoh+1)                        :: downward_lw_incid
   real(kind=8), dimension(ncoh+1)                        :: downward_lw_surf
   real(kind=8), dimension(ncoh+1)                        :: upward_lw_incid
   real(kind=8), dimension(ncoh+1)                        :: upward_lw_surf
   real(kind=8), dimension(ncoh)                          :: source_lw
   real(kind=8), dimension(ncoh)                          :: forcing_lw
   real(kind=8), dimension(ncoh)                          :: A_dw
   real(kind=8), dimension(ncoh)                          :: B_dw
   real(kind=8), dimension(ncoh)                          :: C_dw
   real(kind=8), dimension(ncoh)                          :: D_dw
   real(kind=8), dimension(ncoh)                          :: A_uw
   real(kind=8), dimension(ncoh)                          :: B_uw
   real(kind=8), dimension(ncoh)                          :: C_uw
   real(kind=8), dimension(ncoh)                          :: D_uw
   real(kind=8), dimension(ncoh)                          :: E_uw
   real(kind=8), dimension(ncoh)                          :: F_uw
   real(kind=8)                                           :: emgs
   real(kind=8)                                           :: T_grnd
   real(kind=8)                                           :: zeta
   real(kind=8)                                           :: eta
   real(kind=8)                                           :: exk
   real(kind=8)                                           :: zetai
   real(kind=8)                                           :: exki
   real(kind=8)                                           :: d
   real(kind=8)                                           :: lw_v_surf_tmp
   real(kind=8)                                           :: lw_v_incid_tmp
   real(kind=8)                                           :: cohort_emis
   real(kind=8)                                           :: cohort_backscatter
   real(kind=8)                                           :: cohort_temp_four
   real(kind=8)                                           :: cohort_tai
   real        , dimension(ncoh)                          :: lw_v_surf
   real        , dimension(ncoh)                          :: lw_v_incid
   real                                                   :: downward_lw_below_surf
   real                                                   :: upward_lw_above_surf
   real                                                   :: upward_lw_below_surf
   real                                                   :: downward_lw_below_incid
   real                                                   :: upward_lw_above_incid
   real                                                   :: upward_lw_below_incid
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external                                 :: sngloff
   !---------------------------------------------------------------------------------------!


   !----- Convert some variables to double precision. -------------------------------------!
   emgs   = dble(semgs)
   t_grnd = dble(st_grnd)

   ncoh2  = 2*ncoh

   do il=1,ncoh
      ipft               = pft(il)
      cohort_tai         = lai(il) + wai(il)
      cohort_emis        = ( leaf_emiss_tir(ipft) * lai(il)                                &
                           + wood_emiss_tir(ipft) * wai(il) )                              &
                           / cohort_tai
      cohort_backscatter = ( leaf_backscatter_tir(ipft) * lai(il)                          &
                           + wood_backscatter_tir(ipft) * wai(il) )                        &
                         / cohort_tai
      cohort_temp_four   = ( leaf_temp(il)**4 * leaf_emiss_tir(ipft) * lai(il)             &
                           + wood_temp(il)**4 * wood_emiss_tir(ipft) * wai(il) )           &
                         / ( leaf_emiss_tir(ipft) * lai(il)                                &
                           + wood_emiss_tir(ipft) * wai(il) )
   
      zeta           = (1.d0 - cohort_emis) * cohort_backscatter / mu_bar(ipft)
      eta            = (1.d0 - (1.d0 - cohort_backscatter) * (1.d0 - cohort_emis))         &
                     / mu_bar(ipft)
      exk            = sqrt(eta * eta - zeta * zeta)
      exki           = 1.0d0 / exk
      zetai          = 1.0d0 / zeta
      source_lw (il) = cohort_emis * stefan8 * cohort_temp_four
      forcing_lw(il) = - (zeta + eta) * source_lw(il)
      explai    (il) = exp( exk * cohort_tai)
      exmlai    (il) = exp(-exk * cohort_tai)


      !------------------------------------------------------------------------------------!
      !     Coefficient of lambda1 (and minus the coefficient of lambda2) for the bottom   !
      ! of a layer, downwelling radiation.                                                 !
      !------------------------------------------------------------------------------------!
      a_dw(il) = 5.d-1 * exki

      !----- Coefficient of lambda1, top of layer, downwelling radiation. -----------------!
      b_dw(il) = 5.d-1*exki*explai(il)

      !----- Coefficient of lambda2, top of layer, downwelling radiation. -----------------!
      c_dw(il) = -5.d-1*exki*exmlai(il)

      !----- Term of downwelling radiation not multiplying a lambda. ----------------------!
      d_dw(il) = 5.d-1*(exki**2)*forcing_lw(il) * (explai(il) + exmlai(il) - 2.0d0)

      a_uw(il) =   5.d-1 * zetai * (eta * exki - 1.0d0)
      b_uw(il) = - 5.d-1 * zetai * (eta * exki + 1.0d0)
      c_uw(il) = - source_lw(il) * zetai
      d_uw(il) = a_uw(il) * explai(il)
      e_uw(il) = b_uw(il) * exmlai(il)
      f_uw(il) = -source_lw(il) * zetai                                                    &
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

   masmatcp(1,1)      = b_dw(ncoh)
   masmatcp(1,2)      = c_dw(ncoh)
   mastervec_surf(1)  =-d_dw(ncoh)
   mastervec_incid(1) = dble(rlong_top4)

   do i=2,ncoh2-2,2
      ind = nint(real(ncoh2-i)*0.5)
      masmatcp(i,i-1)    = -a_dw(ind+1)
      masmatcp(i,i)      =  a_dw(ind+1)
      masmatcp(i,i+1)    =  b_dw(ind)
      masmatcp(i,i+2)    =  c_dw(ind)
      mastervec_surf(i)  = -d_dw(ind)
      mastervec_incid(i) = 0.0d0
   end do

   do i=3,ncoh2-1,2
      ind = nint(real(ncoh2-i+1)*0.5)
      masmatcp(i,i-2)    = -a_uw(ind+1)
      masmatcp(i,i-1)    = -b_uw(ind+1)
      masmatcp(i,i)      =  d_uw(ind)
      masmatcp(i,i+1)    =  e_uw(ind)
      mastervec_surf(i)  =  c_uw(ind+1) - f_uw(ind)
      mastervec_incid(i) =  0.0d0
   end do
   masmatcp(ncoh2,ncoh2-1) = a_uw(1) - (1.d0 - emgs) * a_dw(1)
   masmatcp(ncoh2,ncoh2)   = b_uw(1) + (1.d0 - emgs) * a_dw(1)
   mastervec_surf(ncoh2)   = emgs * stefan8 * t_grnd**4 - c_uw(1)
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
                             + masmatcp(i,i+1) * mastervec_surf(i+1) + f_uw(ind)
      upward_lw_incid(ind+1) = masmatcp(i,i) * mastervec_incid(i)                          &
                             + masmatcp(i,i+1) * mastervec_incid(i+1)
   end do

   do i=2,ncoh2-2,2
      ind = nint(real(ncoh2-i)*0.5)
      downward_lw_surf(ind+1)  = masmatcp(i,i+1) * mastervec_surf(i+1)                     &
                               + masmatcp(i,i+2) * mastervec_surf(i+2) + d_dw(ind)
      downward_lw_incid(ind+1) = masmatcp(i,i+1) * mastervec_incid(i+1)                    &
                               + masmatcp(i,i+2) * mastervec_incid(i+2)
   end do

   upward_lw_surf(ncoh+1)    = d_uw(ncoh) * mastervec_surf(1)                              &
                             + e_uw(ncoh) * mastervec_surf(2) + f_uw(ncoh)
   upward_lw_incid(ncoh+1)   = d_uw(ncoh) * mastervec_incid(1)                             &
                             + e_uw(ncoh) * mastervec_incid(2)
   downward_lw_surf(ncoh+1)  = 0.0d0
   downward_lw_incid(ncoh+1) = dble(rlong_top4)
   downward_lw_surf(1)       = a_dw(1) * (mastervec_surf(ncoh2-1)  - mastervec_surf(ncoh2))
   downward_lw_incid(1)      = a_dw(1) * (mastervec_incid(ncoh2-1) - mastervec_incid(ncoh2))
   upward_lw_surf(1)         = (1.0d0-emgs) * downward_lw_surf(1)                          &
                             + emgs * stefan8 * t_grnd**4
   upward_lw_incid(1)        = (1.0d0-emgs) * downward_lw_incid(1)

   do il = 1,ncoh
      lw_v_surf_tmp       = downward_lw_surf(il+1)  - downward_lw_surf(il)                 &
                          + upward_lw_surf(il)  - upward_lw_surf(il+1)
      lw_v_incid_tmp      = downward_lw_incid(il+1) - downward_lw_incid(il)                &
                          + upward_lw_incid(il) - upward_lw_incid(il+1)
      lw_v_surf(il)       = sngloff(lw_v_surf_tmp , tiny_num8)
      lw_v_incid(il)      = sngloff(lw_v_incid_tmp, tiny_num8)
      tir_flip(il)        = lw_v_surf(il) + lw_v_incid(il)
      radprof_flip( 9,il) = sngloff( downward_lw_surf(il) + downward_lw_incid(il)          &
                                   , tiny_num8 )
      radprof_flip(10,il) = sngloff( upward_lw_surf  (il) + upward_lw_incid  (il)          &
                                   , tiny_num8 )
   end do

   downward_lw_below_surf  = sngloff(downward_lw_surf(1)    , tiny_num8)
   downward_lw_below_incid = sngloff(downward_lw_incid(1)   , tiny_num8)
   upward_lw_below_surf    = sngloff(upward_lw_surf(1)      , tiny_num8)
   upward_lw_below_incid   = sngloff(upward_lw_incid(1)     , tiny_num8)
   upward_lw_above_surf    = sngloff(upward_lw_surf(ncoh+1) , tiny_num8)
   upward_lw_above_incid   = sngloff(upward_lw_incid(ncoh+1), tiny_num8)

   dw_tirlo = downward_lw_below_surf + downward_lw_below_incid
   uw_tirlo = upward_lw_below_surf   + upward_lw_below_incid
   uw_tirhi = upward_lw_above_surf   + upward_lw_above_incid

   return
end subroutine old_lw_two_stream
!==========================================================================================!
!==========================================================================================!
