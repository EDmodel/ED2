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
                                   , visible_fraction_dif    ! ! intent(in)
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
   real(kind=8), dimension(ncoh), intent(inout) :: beam_level
   real(kind=8), dimension(ncoh), intent(inout) :: diff_level
   real(kind=8), dimension(ncoh), intent(inout) :: lambda_coh
   real(kind=8)                 , intent(inout) :: lambda_tot
   !----- Local variables -----------------------------------------------------------------!
   integer     , dimension(2*ncoh)        :: indx
   integer                                :: il,ipft,ncoh2,iband,i,j,ind
   real(kind=8), dimension(n_pft)         :: leaf_scatter
   real(kind=8), dimension(n_pft)         :: diffuse_backscatter
   real(kind=8), dimension(ncoh)          :: expkl_top,expkl_bot,expamk_top,expamk_bot
   real(kind=8), dimension(ncoh)          :: expapk_top,expapk_bot,A_top,A_bot,B_top
   real(kind=8), dimension(ncoh)          :: B_bot,C_top,C_bot,F_top,F_bot,G_top,G_bot
   real(kind=8), dimension(ncoh)          :: H_top,H_bot,beam_bot,beam_bot_crown
   real(kind=8), dimension(ncoh+1)        :: upward_vis_beam,upward_vis_diffuse
   real(kind=8), dimension(ncoh+1)        :: upward_nir_beam, upward_nir_diffuse
   real(kind=8), dimension(ncoh+1)        :: downward_nir_beam, downward_nir_diffuse
   real(kind=8), dimension(ncoh+1)        :: downward_vis_beam,downward_vis_diffuse
   real(kind=8), dimension(2*ncoh)        :: mastervec_beam,masveccp_beam
   real(kind=8), dimension(2*ncoh)        :: mastervec_diffuse,masveccp_diffuse
   real(kind=8), dimension(2*ncoh,2)      :: matal
   real(kind=8), dimension(2*ncoh,5)      :: mastermat
   real(kind=8), dimension(2*ncoh,2*ncoh) :: masmatcp  
   real(kind=8)                           :: alb,cosz,cosaoi,lambda
   real(kind=8)                           :: beam_backscatter,eta,zeta
   real(kind=8)                           :: exk,exki,zetai
   real(kind=8)                           :: d,rhoo,sigma,source_bot,source_top
   real(kind=8), dimension(ncoh)          :: eff_tai
   !---------------------------------------------------------------------------------------!

   
   !----- Convert input variable to double precision. -------------------------------------!
   alb    = dble(salb)
   cosz   = max(3.d-2,dble(scosz))
   cosaoi = max(3.d-2,dble(scosaoi))

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
            leaf_scatter(ipft)        = dble(leaf_scatter_nir)
            diffuse_backscatter(ipft) = dble(diffuse_backscatter_nir)
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
         beam_level(il)     = beam_level(il+1)                                             &
                            * exp(-5.d-1*lambda*eff_tai(il)/canopy_area(il))               &
                            * canopy_area(il)                                              &
                            + (1.d0-canopy_area(il)) * beam_level(il+1)
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
                        + 0.5d0*zetai*(eta*exki-1.0d0)*expkl_top(il)*rhoo/(sigma-exk)      &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 0.5d0*zetai*(eta*exki+1.0d0) / expkl_top(il)*rhoo/(sigma+exk)    &
                        * (expapk_top(il)-expapk_bot(il))
         B_bot(il)      = 0.5d0*zetai*(eta*exki-1.0d0)
         B_top(il)      = 0.5d0*zetai*(eta*exki-1.0d0)*expkl_top(il)
         C_bot(il)      = -0.5d0*zetai*(eta*exki+1.0d0)
         C_top(il)      = -0.5d0*zetai*(eta*exki+1.0d0)/expkl_top(il)
         F_bot(il)      = 0.0d0
         F_top(il)      = 0.5d0*exki*expkl_top(il)*rhoo/(sigma-exk)                        &
                        * (expamk_top(il)-expamk_bot(il))                                  &
                        - 0.5d0*exki/expkl_top(il)*rhoo/(sigma+exk)                        &
                        * (expapk_top(il)-expapk_bot(il))
         G_bot(il)      = 0.5d0*exki
         G_top(il)      = 0.5d0*exki*expkl_top(il)
         H_bot(il)      = -0.5d0*exki
         H_top(il)      = -0.5d0*exki/expkl_top(il)
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
      call bandec(mastermat,ncoh2,2,2,matal,indx,d)

      !----- Backsubstitute for beam and diffuse. -----------------------------------------!
      call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_beam)
      call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_diffuse)
      
      !----- Improve the solution. --------------------------------------------------------!
      call mprove(masmatcp,mastermat,matal,ncoh2,5,2,indx,masveccp_beam,mastervec_beam)
      call mprove(masmatcp,mastermat,matal,ncoh2,5,2,indx,masveccp_diffuse                 &
                 ,mastervec_diffuse)

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
   
   !----- Copying to the output variables. ------------------------------------------------!
   DW_vislo_beam    = max(0.0,sngl(downward_vis_beam(1)))       * visible_fraction_dir
   DW_vislo_diffuse = max(0.0,sngl(downward_vis_diffuse(1)))    * visible_fraction_dif
   UW_vishi_beam    = max(0.0,sngl(upward_vis_beam(ncoh+1)))    * visible_fraction_dir
   UW_vishi_diffuse = max(0.0,sngl(upward_vis_diffuse(ncoh+1))) * visible_fraction_dif
   DW_nirlo_beam    = max(0.0,sngl(downward_nir_beam(1)))       * (1.-visible_fraction_dir)
   DW_nirlo_diffuse = max(0.0,sngl(downward_nir_diffuse(1)))    * (1.-visible_fraction_dif)
   UW_nirhi_beam    = max(0.0,sngl(upward_nir_beam(ncoh+1)))    * (1.-visible_fraction_dir)
   UW_nirhi_diffuse = max(0.0,sngl(upward_nir_diffuse(ncoh+1))) * (1.-visible_fraction_dif)
   
   return
end subroutine sw_twostream_clump
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the matrix pseudo-inverse.                                       !
!------------------------------------------------------------------------------------------!
subroutine bandec(a,n,m1,m2,al,indx,d)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                           , intent(in)    :: n,m1,m2
   integer, dimension(n)             , intent(out)   :: indx
   real(kind=8), dimension(n,m1+m2+1), intent(inout) :: a
   real(kind=8), dimension(n,m1)     , intent(inout) :: al
   real(kind=8)                      , intent(out)   :: d
   !----- Local variables -----------------------------------------------------------------!
   integer                                           :: i,j,k,l,mm
   real(kind=8)                                      :: tvar
   real(kind=8)                                      :: dum 
   !----- Local constants -----------------------------------------------------------------!
   real(kind=8)                      , parameter     :: tiny_offset=1.0d-20
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
   do k=1,n
      dum=a(k,1)
      i=k
      if(l.lt.n)l=l+1
      do j=k+1,l
         if (abs(a(j,1)) > abs(dum)) then
            dum=a(j,1)
            i=j
         end if
      end do
      indx(k)=i
      if(dum == 0.0d0) a(k,1)=tiny_offset
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
subroutine banbks(a,n,m1,m2,al,indx,b)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                           , intent(in)  :: n,m1,m2
   integer     , dimension(n)        , intent(in)  :: indx
   real(kind=8), dimension(n,m1+m2+1), intent(in)  :: a
   real(kind=8), dimension(n,m1)     , intent(in)  :: al
   real(kind=8), dimension(n)        , intent(out) :: b
   !----- Local variables -----------------------------------------------------------------!
   integer                                         :: i,k,l,mm
   real(kind=8)                                    :: dum,tvar
   !---------------------------------------------------------------------------------------!

   mm=m1+m2+1
   l=m1
   do k=1,n
      i=indx(k)
      if (i /= k) then
         tvar = b(k)
         b(k) = b(i)
         b(i) = tvar
      end if
      if(l < n)l=l+1
      do i=k+1,l
         b(i) = b(i) - al(k,i-k) * b(k)
      end do
   end do

   l=1
   do i=n,1,-1
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
subroutine mprove(a,alud,matal,n,np,npp,indx,b,x)
   implicit none
   !----- Local constants -----------------------------------------------------------------!
   integer     , parameter :: nmax=100
   !----- Arguments -----------------------------------------------------------------------!
   integer     , dimension(n)    , intent(in)    :: indx
   integer                       , intent(in)    :: n,np,npp
   real(kind=8), dimension(n,n)  , intent(in)    :: a
   real(kind=8), dimension(n)    , intent(in)    :: b
   real(kind=8), dimension(n,np) , intent(in)    :: alud
   real(kind=8), dimension(n,npp), intent(in)    :: matal
   real(kind=8), dimension(n)    , intent(inout) :: x
   !----- Local variables -----------------------------------------------------------------!
   integer                                       :: i,j
   real(kind=8), dimension(n)                    :: r
   real(kind=8)                                  :: sdp
   !---------------------------------------------------------------------------------------!


   do i=1,n
      sdp = -b(i)
      do j=1,n
         sdp = sdp + a(i,j) * x(j)
      end do
      r(i) = sdp
   end do

   call banbks(alud,n,npp,npp,matal,indx,r)

   do i=1,n
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
   use consts_coms           , only : stefan8         ! ! intent(in)
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
   real(kind=8), dimension(ncoh+1)            :: DW_incid
   real(kind=8), dimension(ncoh+1)            :: DW_surf
   real(kind=8), dimension(ncoh+1)            :: UW_incid
   real(kind=8), dimension(ncoh+1)            :: UW_surf
   real(kind=8), dimension(ncoh)              :: source
   real(kind=8), dimension(ncoh)              :: forcing
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
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter                     :: negligible = 1.d-32
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external                     :: sngloff
   !---------------------------------------------------------------------------------------!


   !----- Converting some variables to double precision. ----------------------------------!
   emgs   = dble(semgs)
   t_grnd = dble(st_grnd)

   ncoh2  = 2*ncoh

   do il=1,ncoh
      zeta        = 2.0d0 * (1.0d0 - emis_v(pft(il))) / (3.0d0 * mubar)
      eta         = (2.0d0 + emis_v(pft(il)))/(3.0d0 * mubar)
      exk         = sqrt(eta*eta-zeta*zeta)
      exki        = 1.0d0 / exk
      zetai       = 1.0d0 / zeta
      source(il)  = emis_v(pft(il)) * stefan8 *T_veg(il)**4
      forcing(il) = -(zeta + eta) * source(il)
      !explai(il)  = exp( exk * clumping_factor(pft(il))*TAI(il))
      !exmlai(il)  = exp(-exk * clumping_factor(pft(il))*TAI(il))
      explai(il)  = exp( exk *TAI(il))
      exmlai(il)  = exp(-exk *TAI(il))


      !------------------------------------------------------------------------------------!
      !     Coefficient of lambda1 (and minus the coefficient of lambda2) for the bottom   !
      ! of a layer, downwelling radiation.                                                 !
      !------------------------------------------------------------------------------------!
      A_dw(il) = 0.5d0 * exki

      !----- Coefficient of lambda1, top of layer, downwelling radiation. -----------------!
      B_dw(il) = 0.5d0*exki*explai(il)

      !----- Coefficient of lambda2, top of layer, downwelling radiation. -----------------!
      C_dw(il) = -0.5d0*exki*exmlai(il)

      !----- Term of downwelling radiation not multiplying a lambda. ----------------------!
      D_dw(il) = 0.5d0*(exki**2)*forcing(il) * (explai(il) + exmlai(il) - 2.0d0)

      A_uw(il) =   0.5d0 * zetai * (eta * exki - 1.0d0)
      B_uw(il) = - 0.5d0 * zetai * (eta * exki + 1.0d0)
      C_uw(il) = - source(il) * zetai
      D_uw(il) = A_uw(il) * explai(il)
      E_uw(il) = B_uw(il) * exmlai(il)
      F_uw(il) = -source(il) * zetai                                                       &
               + 0.5d0 * zetai * (eta*exki - 1.0d0) * explai(il)                           &
               * (forcing(il) * exki * (1.0d0-exmlai(il)))                                 &
               - 0.5d0 * zetai * (eta*exki + 1.0d0) * exmlai(il)                           &
               * (forcing(il) * exki * (explai(il)-1.0d0))
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
   call bandec(mastermat,ncoh2,2,2,matal,indx,d)

   !----- Backsubstitute for contributions of ground and vegetation. ----------------------!
   call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_surf)

   !----- Backsubstitute for contribution of incident longwave at canopy top. -------------!
   call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_incid)

   do i=3,ncoh2-1,2
      ind = nint(real(ncoh2-i+1)*0.5)
      UW_surf(ind+1)  = masmatcp(i,i) * mastervec_surf(i)                                  &
                      + masmatcp(i,i+1) * mastervec_surf(i+1) + F_uw(ind)
      UW_incid(ind+1) = masmatcp(i,i) * mastervec_incid(i)                                 &
                      + masmatcp(i,i+1) * mastervec_incid(i+1)
   end do

   do i=2,ncoh2-2,2
      ind = nint(real(ncoh2-i)*0.5)
      DW_surf(ind+1)  = masmatcp(i,i+1) * mastervec_surf(i+1)                              &
                      + masmatcp(i,i+2) * mastervec_surf(i+2) + D_dw(ind)
      DW_incid(ind+1) = masmatcp(i,i+1) * mastervec_incid(i+1)                             &
                      + masmatcp(i,i+2) * mastervec_incid(i+2)
   end do

   UW_surf(ncoh+1)  = D_uw(ncoh) * mastervec_surf(1)                                       &
                    + E_uw(ncoh) * mastervec_surf(2) + F_uw(ncoh)
   UW_incid(ncoh+1) = D_uw(ncoh) * mastervec_incid(1)                                      &
                    + E_uw(ncoh) * mastervec_incid(2)
   DW_surf(ncoh+1)  = 0.0d0
   DW_incid(ncoh+1) = 1.0d0
   DW_surf(1)       = A_dw(1) * (mastervec_surf(ncoh2-1)  - mastervec_surf(ncoh2))
   DW_incid(1)      = A_dw(1) * (mastervec_incid(ncoh2-1) - mastervec_incid(ncoh2))
   UW_surf(1)       = (1.0d0-emgs) * DW_surf(1) + emgs * stefan8 * T_grnd**4
   UW_incid(1)      = (1.0d0-emgs) * DW_incid(1)

   do il = 1,ncoh
      lw_v_surf_tmp  = DW_surf(il+1)  - DW_surf(il)  + UW_surf(il)  - UW_surf(il+1)
      lw_v_incid_tmp = DW_incid(il+1) - DW_incid(il) + UW_incid(il) - UW_incid(il+1)
      lw_v_surf(il)  = sngloff(lw_v_surf_tmp , negligible)
      lw_v_incid(il) = sngloff(lw_v_incid_tmp, negligible)
   end do

   downward_lw_below_surf  = sngloff(DW_surf(1)      , negligible)
   downward_lw_below_incid = sngloff(DW_incid(1)     , negligible)
   upward_lw_below_surf    = sngloff(UW_surf(1)      , negligible)
   upward_lw_below_incid   = sngloff(UW_incid(1)     , negligible)
   upward_lw_above_surf    = sngloff(UW_surf(ncoh+1) , negligible)
   upward_lw_above_incid   = sngloff(UW_incid(ncoh+1), negligible)

   return
end subroutine lw_twostream
!==========================================================================================!
!==========================================================================================!
