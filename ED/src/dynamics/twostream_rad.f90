!==========================================================================================!
!==========================================================================================!
!     This sub-routine solves the long-wave radiation method using the two-stream model.   !
! We consider the finite crown area when computing the layer transmittance.                !
! Contrary to short-wave radiation, we don't normalise long-wave radiation.                !
!                                                                                          !
! References:                                                                              !
!                                                                                          !
! Liou, K. N., 2002: An introduction to atmospheric radiation. Academic Press, Intl.       !
!    Geophys. Series, vol. 84.  Chapter 6., p. 257-347. (L02)                              !
!                                                                                          !
! Oleson, K. W., and co-authors, 2004: Technical description of version 4.0 of the         !
!    community land model (CLM). NCAR Technical note NCAR/TN-478+STR. 257pp. (CLM10)       !
!                                                                                          !
! Sellers, P. J., 1985: Canopy reflectance, photosynthesis and transpiration. Intl. J.     !
!    of remote sensing, 6, 1335-1372. (S85)                                                !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine lw_two_stream(grnd_emiss4,grnd_temp4,rlong_top4,ncoh,pft,lai,wai,cai,leaf_temp  &
                        ,wood_temp,radprof_flip,tir_flip,dw_tirlo,uw_tirlo,uw_tirhi)
   use ed_max_dims          , only : n_pft                   & ! intent(in)
                                   , n_radprof               ! ! intent(in)
   use rk4_coms             , only : tiny_offset             ! ! intent(in)
   use canopy_radiation_coms, only : clumping_factor         & ! intent(in)
                                   , orient_factor           & ! intent(in)
                                   , leaf_emiss_tir          & ! intent(in)
                                   , wood_emiss_tir          & ! intent(in)
                                   , phi1                    & ! intent(in)
                                   , phi2                    & ! intent(in)
                                   , mu_bar                  & ! intent(in)
                                   , leaf_backscatter_tir    & ! intent(in)
                                   , wood_backscatter_tir    ! ! intent(in)
   use consts_coms          , only : stefan8                 ! ! intent(in)
   use ed_misc_coms         , only : current_time            ! ! intent(in)
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4)                              , intent(in)    :: grnd_emiss4
   real(kind=4)                              , intent(in)    :: grnd_temp4
   real(kind=4)                              , intent(in)    :: rlong_top4
   integer                                   , intent(in)    :: ncoh
   integer     , dimension(ncoh)             , intent(in)    :: pft
   real(kind=8), dimension(ncoh)             , intent(in)    :: lai
   real(kind=8), dimension(ncoh)             , intent(in)    :: wai
   real(kind=8), dimension(ncoh)             , intent(in)    :: cai
   real(kind=8), dimension(ncoh)             , intent(in)    :: leaf_temp
   real(kind=8), dimension(ncoh)             , intent(in)    :: wood_temp
   real(kind=4), dimension(n_radprof,ncoh)   , intent(inout) :: radprof_flip
   real(kind=4), dimension(ncoh)             , intent(out)   :: tir_flip
   real(kind=4)                              , intent(out)   :: dw_tirlo
   real(kind=4)                              , intent(out)   :: uw_tirlo
   real(kind=4)                              , intent(out)   :: uw_tirhi
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                   :: ipft
   integer                                                   :: nsiz
   integer                                                   :: i
   integer                                                   :: im1
   integer                                                   :: ip1
   integer                                                   :: i2
   integer                                                   :: i2p1
   integer                                                   :: i2m1
   integer                                                   :: i2p2
   logical                                                   :: sing
   real(kind=8), dimension(ncoh+1)                           :: black
   real(kind=8), dimension(ncoh+1)                           :: expl_plus
   real(kind=8), dimension(ncoh+1)                           :: expl_minus
   real(kind=8), dimension(ncoh+1)                           :: gamm_plus
   real(kind=8), dimension(ncoh+1)                           :: gamm_minus
   real(kind=8), dimension(ncoh+1)                           :: down
   real(kind=8), dimension(ncoh+1)                           :: up
   real(kind=8), dimension(2*ncoh+2,2*ncoh+2)                :: mmat
   real(kind=8), dimension(2*ncoh+2)                         :: yvec
   real(kind=8), dimension(2*ncoh+2)                         :: xvec
   real(kind=8), dimension(ncoh+1)                           :: elai
   real(kind=8), dimension(ncoh+1)                           :: etai
   real(kind=8), dimension(ncoh+1)                           :: mu
   real(kind=8), dimension(ncoh+1)                           :: leaf_weight
   real(kind=8), dimension(ncoh+1)                           :: wood_weight
   real(kind=8), dimension(ncoh+1)                           :: beta
   real(kind=8), dimension(ncoh+1)                           :: epsil
   real(kind=8), dimension(ncoh+1)                           :: iota
   real(kind=8), dimension(ncoh+1)                           :: lambda
   real(kind=8)                                              :: temiss_four
   real(kind=8)                                              :: iota_g
   real(kind=8)                                              :: black_g
   real(kind=8)                                              :: down_sky
   !----- External functions. -------------------------------------------------------------!
   real(kind=8)                              , external      :: eifun8
   real(kind=4)                              , external      :: sngloff
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the size of the solution matrix and array.                                   !
   !---------------------------------------------------------------------------------------!
   nsiz = 2*ncoh + 2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert some properties to double precision.                                      !
   !---------------------------------------------------------------------------------------!
   iota_g   = 1.d0 - dble(grnd_emiss4)
   black_g  = stefan8 * dble(grnd_temp4) ** 4
   down_sky = dble(rlong_top4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert some leaf variables to double precision, then find some general           !
   ! properties of leaves/branches (i.e., properties that do not depend on which band we   !
   ! are solving):  the extinction coefficient lambda, the layer transmittance             !
   ! coefficients tau_beam and tau_diff, and the weight of LAI and WAI relative to the     !
   ! total area.                                                                           !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      ipft = pft(i)
      !----- Area indices, correct LAI by clumping factor. --------------------------------!
      elai(i) = clumping_factor(ipft) * lai(i)
      etai(i) = elai(i) + wai(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    The weighting factors are defined by the relative LAI and WAI areas.            !
      !------------------------------------------------------------------------------------!
      leaf_weight(i)  =  elai(i) / etai(i)
      wood_weight(i)  = 1.d0 -  leaf_weight(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     mu is the corrected mu_bar to account for finite crown area, if it is          !
      ! provided, otherwise, it is just the same as mu_bar.                                !
      !------------------------------------------------------------------------------------!
      mu(i) = - etai(i) / log( 1.d0 - cai(i)                                               &
                             + cai(i) * exp ( - etai(i) / ( cai(i) * mu_bar(ipft) ) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Black is the emission of a black body at the same temperature of the layer.   !
      !  The temperature is the weighted average where the weights are the product between !
      !  the area and the emissivity.                                                      !
      !------------------------------------------------------------------------------------!
      temiss_four = ( leaf_weight(i) * leaf_emiss_tir(ipft) * leaf_temp(i) ** 4            &
                    + wood_weight(i) * wood_emiss_tir(ipft) * wood_temp(i) ** 4 )          &
                  / ( leaf_weight(i) * leaf_emiss_tir(ipft)                                &
                    + wood_weight(i) * wood_emiss_tir(ipft) )
      black(i) = cai(i) * stefan8 * temiss_four
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      iota is the single scattering albedo.  Because we assume absorptivity to be   !
      ! the same as emissivity, we use the latter to define it.                            !
      !------------------------------------------------------------------------------------!
      iota(i) = 1.d0 - ( leaf_weight(i) * leaf_emiss_tir(ipft)                             &
                       + wood_weight(i) * wood_emiss_tir(ipft) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Beta is the backward scattering of diffuse radiation. Epsil = 1 - 2*beta.      !
      !------------------------------------------------------------------------------------!
      beta (i) = leaf_weight(i) * leaf_backscatter_tir(ipft)                               &
               + wood_weight(i) * wood_backscatter_tir(ipft)
      epsil(i) = 1.d0 - 2.d0 * beta(i)
      !------------------------------------------------------------------------------------!



      !----- lambda is the coefficient associated with the optical depth. -----------------!
      lambda(i) = sqrt( ( 1.d0 - epsil(i) * iota(i) ) * ( 1.d0 - iota(i) ) ) / mu(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    gamm_plus and gamm_minus are the coefficients that relate upwelling and down-   !
      ! welling radiation.                                                                 !
      !------------------------------------------------------------------------------------!
      gamm_plus (i) = 5.d-1 * ( 1.d0 + sqrt( ( 1.d0 -            iota(i) )                 &
                                           / ( 1.d0 - epsil(i) * iota(i) ) ) )
      gamm_minus(i) = 5.d-1 * ( 1.d0 - sqrt( ( 1.d0 -            iota(i) )                 &
                                           / ( 1.d0 - epsil(i) * iota(i) ) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    expl_plus and expl_minus are the transmitivity of diffuse light.                !
      !------------------------------------------------------------------------------------!
      expl_plus (i) = exp(   lambda(i) * etai(i) )
      expl_minus(i) = exp( - lambda(i) * etai(i) )
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Define some boundary conditions for the vector properties above.                  !
   !---------------------------------------------------------------------------------------!
   i              = ncoh+1
   elai       (i) = 0.d0
   etai       (i) = 0.d0
   leaf_weight(i) = 5.d-1
   wood_weight(i) = 5.d-1
   mu         (i) = 1.d0
   black      (i) = 0.d0
   iota       (i) = 1.d0
   beta       (i) = 0.d0
   epsil      (i) = 1.d0 - 2.d0 * beta(i)
   lambda     (i) = 0.d0
   gamm_plus  (i) = 1.d0
   gamm_minus (i) = 0.d0
   expl_plus  (i) = 1.d0
   expl_minus (i) = 1.d0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fill in the right hand side vector (Y) and the matrix (X)                         !
   !---------------------------------------------------------------------------------------!
   !------ Initialise the vector and the matrix.  The matrix is a sparse one... -----------!
   mmat(:,:)         = 0.d0
   yvec(:)           = 0.d0
   !------ Add the bottom and top boundary conditions. ------------------------------------!
   mmat(1,1)         = (gamm_minus(1) - iota_g * gamm_plus (1)) * expl_minus(1)
   mmat(1,2)         = (gamm_plus (1) - iota_g * gamm_minus(1)) * expl_plus (1)
   mmat(nsiz,nsiz-1) = gamm_plus (ncoh+1)
   mmat(nsiz,nsiz  ) = gamm_minus(ncoh+1)
   yvec(1)           = (1.d0 - iota_g) * black_g - (1.d0 - iota_g) * black(1)
   yvec(nsiz)        = down_sky - black(ncoh+1)
   do i=1,ncoh
      !----- Find auxiliary indices. ------------------------------------------------------!
      ip1  =  i + 1
      i2   =  2 * i
      i2m1 = i2 - 1
      i2p1 = i2 + 1
      i2p2 = i2 + 2
      !------------------------------------------------------------------------------------!

      !----- Make elements of C vector. ---------------------------------------------------!
      yvec(i2  ) = black(ip1) - black(i)
      yvec(i2p1) = black(ip1) - black(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Make elements of the tri-diagonal A array.                                     !
      !------------------------------------------------------------------------------------!
      mmat(i2  ,i2m1) =   gamm_plus   (i)
      mmat(i2  ,i2  ) =   gamm_minus  (i)
      mmat(i2  ,i2p1) = - gamm_plus (ip1) * expl_minus(ip1)
      mmat(i2  ,i2p2) = - gamm_minus(ip1) * expl_plus (ip1)
      mmat(i2p1,i2m1) =   gamm_minus  (i)
      mmat(i2p1,i2  ) =   gamm_plus   (i)
      mmat(i2p1,i2p1) = - gamm_minus(ip1) * expl_minus(ip1)
      mmat(i2p1,i2p2) = - gamm_plus (ip1) * expl_plus (ip1)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !       Solve the linear system.  In the future we could use a band diagonal solver,    !
   ! which is a lot cheaper than the regular Gauss elimination, but for the time being, we !
   ! go with a tested method.                                                              !
   !---------------------------------------------------------------------------------------!
   call lisys_solver8(nsiz,mmat,yvec,xvec,sing)
   if (sing) then
      call fatal_error('LW radiation failed... The matrix is singular!'                    &
                      ,'lw_two_stream','twostream_rad.f90')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Copy the solution to the vectors, using the properties:                           !
   !  Un = Un(P)                                                                           !
   !  Dn = Dn(P)                                                                           !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh+1

      !----- Auxiliary indices. -----------------------------------------------------------!
      i2   =  2 * i
      i2m1 = i2 - 1
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Retrieve the downward and upward diffuse (hemispheric) radiation, by using    !
      ! the solutions for full TAI                                                         !
      !------------------------------------------------------------------------------------!
      down(i) = xvec(i2m1) * gamm_plus (i) * expl_minus(i)                                 &
              + xvec  (i2) * gamm_minus(i) * expl_plus (i)                                 &
              + black  (i)
      up  (i) = xvec(i2m1) * gamm_minus(i) * expl_minus(i)                                 &
              + xvec  (i2) * gamm_plus (i) * expl_plus (i)                                 &
              + black  (i)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !------ Save the fluxes reaching the surface and leaving the top. ----------------------!
   dw_tirlo = sngloff(down     (1), tiny_offset)
   uw_tirlo = sngloff(up       (1), tiny_offset)
   uw_tirhi = sngloff(up  (ncoh+1), tiny_offset)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Save the radiation fluxes to the output variable.                                 !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      ip1 = i + 1
      radprof_flip ( 9,i) = sngloff(down(i), tiny_offset)
      radprof_flip (10,i) = sngloff(up  (i), tiny_offset)
      tir_flip        (i) = sngloff(down(ip1) - down(i) + up(i) - up(ip1), tiny_offset)
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine lw_two_stream
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine solves the short-wave radiation method using the two-stream model.  !
! We consider the finite crown area when computing the transmittance of both direct and    !
! diffuse radiation.  References:                                                          !
!                                                                                          !
! Liou, K. N., 2002: An introduction to atmospheric radiation. Academic Press, Intl.       !
!    Geophys. Series, vol. 84.  Chapter 6., p. 257-347. (L02)                              !
!                                                                                          !
! Oleson, K. W., and co-authors, 2004: Technical description of the community land model   !
!   (CLM). NCAR Technical note NCAR/TN-461+STR. 186pp. (CLM04)                             !
!                                                                                          !
! Sellers, P. J., 1985: Canopy reflectance, photosynthesis and transpiration. Intl. J.     !
!    of remote sensing, 6, 1335-1372. (S85)                                                !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine sw_two_stream(grnd_alb_par4,grnd_alb_nir4,cosaoi4,ncoh,pft,lai,wai,cai          &
                        ,radprof_flip,par_beam_flip,par_diff_flip,sw_abs_beam_flip         &
                        ,sw_abs_diff_flip,dw_parlo_beam,dw_parlo_diff                      &
                        ,uw_parhi_diff,dw_nirlo_beam,dw_nirlo_diff,uw_nirhi_diff           &
                        ,par_beam_level,par_diff_level,light_level,light_beam_level        &
                        ,light_diff_level)
   use ed_max_dims          , only : n_pft                   & ! intent(in)
                                   , n_radprof               ! ! intent(in)
   use rk4_coms             , only : tiny_offset             ! ! intent(in)
   use canopy_radiation_coms, only : clumping_factor         & ! intent(in)
                                   , orient_factor           & ! intent(in)
                                   , phi1                    & ! intent(in)
                                   , phi2                    & ! intent(in)
                                   , mu_bar                  & ! intent(in)
                                   , leaf_scatter_vis        & ! intent(in)
                                   , leaf_scatter_nir        & ! intent(in)
                                   , leaf_backscatter_nir    & ! intent(in)
                                   , leaf_backscatter_vis    & ! intent(in)
                                   , wood_scatter_nir        & ! intent(in)
                                   , wood_scatter_vis        & ! intent(in)
                                   , wood_backscatter_nir    & ! intent(in)
                                   , wood_backscatter_vis    & ! intent(in)
                                   , par_beam_norm           & ! intent(in)
                                   , par_diff_norm           & ! intent(in)
                                   , nir_beam_norm           & ! intent(in)
                                   , nir_diff_norm           & ! intent(in)
                                   , cosz_min8               ! ! intent(in)
   use ed_misc_coms         , only : current_time            ! ! intent(in)

   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4)                              , intent(in)     :: grnd_alb_par4
   real(kind=4)                              , intent(in)     :: grnd_alb_nir4
   real(kind=4)                              , intent(in)     :: cosaoi4
   integer                                   , intent(in)     :: ncoh
   integer     , dimension(ncoh)             , intent(in)     :: pft
   real(kind=8), dimension(ncoh)             , intent(in)     :: lai
   real(kind=8), dimension(ncoh)             , intent(in)     :: wai
   real(kind=8), dimension(ncoh)             , intent(in)     :: cai
   real(kind=4), dimension(n_radprof,ncoh)   , intent(inout)  :: radprof_flip
   real(kind=4), dimension(ncoh)             , intent(out)    :: par_beam_flip
   real(kind=4), dimension(ncoh)             , intent(out)    :: par_diff_flip
   real(kind=4), dimension(ncoh)             , intent(out)    :: sw_abs_beam_flip
   real(kind=4), dimension(ncoh)             , intent(out)    :: sw_abs_diff_flip
   real(kind=4)                              , intent(out)    :: uw_parhi_diff
   real(kind=4)                              , intent(out)    :: uw_nirhi_diff
   real(kind=4)                              , intent(out)    :: dw_parlo_beam
   real(kind=4)                              , intent(out)    :: dw_parlo_diff
   real(kind=4)                              , intent(out)    :: dw_nirlo_beam
   real(kind=4)                              , intent(out)    :: dw_nirlo_diff
   real(kind=8), dimension(ncoh)             , intent(out)    :: par_beam_level
   real(kind=8), dimension(ncoh)             , intent(out)    :: par_diff_level
   real(kind=8), dimension(ncoh)             , intent(out)    :: light_level
   real(kind=8), dimension(ncoh)             , intent(out)    :: light_beam_level
   real(kind=8), dimension(ncoh)             , intent(out)    :: light_diff_level
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                    :: ipft
   integer                                                    :: iband
   integer                                                    :: nsiz
   integer                                                    :: i
   integer                                                    :: ip1
   integer                                                    :: i2
   integer                                                    :: i2p1
   integer                                                    :: i2m1
   integer                                                    :: i2p2
   logical                                                    :: sing
   real(kind=4), dimension(ncoh)                              :: nir_beam_flip
   real(kind=4), dimension(ncoh)                              :: nir_diff_flip
   real(kind=8), dimension(ncoh+1)                            :: expl_plus
   real(kind=8), dimension(ncoh+1)                            :: expl_minus
   real(kind=8), dimension(ncoh+1)                            :: expm0_minus
   real(kind=8), dimension(ncoh+1)                            :: gamm_plus
   real(kind=8), dimension(ncoh+1)                            :: gamm_minus
   real(kind=8), dimension(ncoh+1)                            :: delta
   real(kind=8), dimension(ncoh+1)                            :: upsilon
   real(kind=8), dimension(ncoh+1)                            :: down0
   real(kind=8), dimension(ncoh+1)                            :: down
   real(kind=8), dimension(ncoh+1)                            :: up
   real(kind=8), dimension(2*ncoh+2,2*ncoh+2)                 :: mmat
   real(kind=8), dimension(2*ncoh+2)                          :: yvec
   real(kind=8), dimension(2*ncoh+2)                          :: xvec
   real(kind=8), dimension(n_pft)                             :: leaf_scatter
   real(kind=8), dimension(n_pft)                             :: wood_scatter
   real(kind=8), dimension(n_pft)                             :: leaf_backscatter
   real(kind=8), dimension(n_pft)                             :: wood_backscatter
   real(kind=8), dimension(ncoh+1)                            :: elai
   real(kind=8), dimension(ncoh+1)                            :: etai
   real(kind=8), dimension(ncoh+1)                            :: mu
   real(kind=8), dimension(ncoh+1)                            :: mu0
   real(kind=8), dimension(ncoh+1)                            :: leaf_weight
   real(kind=8), dimension(ncoh+1)                            :: wood_weight
   real(kind=8), dimension(ncoh+1)                            :: beta0
   real(kind=8), dimension(ncoh+1)                            :: beta
   real(kind=8), dimension(ncoh+1)                            :: epsil
   real(kind=8), dimension(ncoh+1)                            :: epsil0
   real(kind=8), dimension(ncoh+1)                            :: iota
   real(kind=8), dimension(ncoh+1)                            :: lambda
   real(kind=8), dimension(ncoh+1)                            :: iota_ratio
   real(kind=8), dimension(ncoh+1)                            :: proj_area
   real(kind=8), dimension(ncoh+1)                            :: a_aux
   real(kind=8), dimension(ncoh+1)                            :: s_aux
   real(kind=8)                                               :: iota_g
   real(kind=8)                                               :: iota_g_par
   real(kind=8)                                               :: iota_g_nir
   real(kind=8)                                               :: down_sky
   real(kind=8)                                               :: down0_sky
   real(kind=8)                                               :: czen
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)                              , external       :: sngloff
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the size of the solution matrix and array.                                   !
   !---------------------------------------------------------------------------------------!
   nsiz = 2*ncoh + 2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert cosine of angle of incidence to double precision.                         !
   !---------------------------------------------------------------------------------------!
   czen       = max(cosz_min8,dble(cosaoi4))
   iota_g_par = dble(grnd_alb_par4)
   iota_g_nir = dble(grnd_alb_nir4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the light level variables.                                             !
   !---------------------------------------------------------------------------------------!
   light_level     (:) = 0.d0
   light_beam_level(:) = 0.d0
   light_diff_level(:) = 0.d0
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Here we loop over the two bands (band 1 is PAR and band 2 is NIR), and solve the  !
   ! same set of equations.                                                                !
   !---------------------------------------------------------------------------------------!
   bandloop: do iband = 1,2
      select case (iband)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Visible (PAR).                                                              !
         !---------------------------------------------------------------------------------!
         leaf_scatter    (:) = leaf_scatter_vis    (:)
         wood_scatter    (:) = wood_scatter_vis    (:)
         leaf_backscatter(:) = leaf_backscatter_vis(:)
         wood_backscatter(:) = wood_backscatter_vis(:)
         down_sky            = par_diff_norm
         down0_sky           = par_beam_norm
         iota_g              = iota_g_par
         !---------------------------------------------------------------------------------!
      case (2)
         !---------------------------------------------------------------------------------!
         !     Near-Infrared (NIR).                                                        !
         !---------------------------------------------------------------------------------!
         leaf_scatter    (:) = leaf_scatter_nir    (:)
         wood_scatter    (:) = wood_scatter_nir    (:)
         leaf_backscatter(:) = leaf_backscatter_nir(:)
         wood_backscatter(:) = wood_backscatter_nir(:)
         down_sky            = nir_diff_norm
         down0_sky           = nir_beam_norm
         iota_g              = iota_g_nir
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find some general properties of leaves/branches and also direct radiation      !
      ! properties.                                                                        !
      !------------------------------------------------------------------------------------!
      directloop: do i=1,ncoh
         ipft      = pft(i)

         !----- Area indices, correct LAI by clumping factor. -----------------------------!
         elai(i) = clumping_factor(ipft) * lai(i)
         etai(i) = elai(i) + wai(i)
         !---------------------------------------------------------------------------------!



         !----- The weighting factors are defined by the relative LAI and WAI areas. ------!
         leaf_weight(i) = elai(i) / etai(i)
         wood_weight(i) = 1.d0 - leaf_weight(i)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     We find the inverse optical depth of the direct radiation (mu0), following  !
         ! CLM10 (equation 3.3 and text after equation 3.3).                               !
         !---------------------------------------------------------------------------------!
         proj_area(i) = phi1(ipft) + phi2(ipft) * czen
         mu0      (i) = - etai(i)                                                          &
                        / log( ( 1.d0 - cai(i) )                                           &
                             + cai(i) * exp( - proj_area(i) * etai(i)                      &
                                             / ( cai(i) * czen ) ) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     The inverse optical depth of the diffuse radiation (mu), same as in long-   !
         ! wave radiation.                                                                 !
         !---------------------------------------------------------------------------------!
         mu (i) = - etai(i) / log( ( 1.d0 - cai(i) )                                       &
                                 + cai(i) * exp( - etai(i) / mu_bar(ipft) ) )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the backscatter for direct beam radiation, following CLM10, equations !
         ! (3.14) and (3.15).  The forward scattering was omitted on both equations        !
         ! because they are different for PAR and NIR, but in both cases they would cancel !
         ! out.  epsil0 = 1 - 2*beta0, no special definition...                            !
         !---------------------------------------------------------------------------------!
         iota_ratio(i) = 1.d0 / ( 2.d0 * ( 1.d0 + phi2(ipft) * mu0(i) ) )                  &
                       * ( 1.d0 - phi1(ipft) * mu0(i) / ( 1.d0 + phi2(ipft) * mu0(i) )     &
                                * log ( ( 1.d0  + ( phi1(ipft) + phi2(ipft) ) * mu0(i) )   &
                                      / ( phi1(ipft) * mu0(i) ) ) )
         beta0     (i) = iota_ratio(i) * ( 1.d0 + mu0(i) / mu(i) )
         epsil0    (i) = 1.d0 - 2.d0 * beta0(i)
         !---------------------------------------------------------------------------------!


         !----- Expm0_minus is the transmissivity of direct radiation. --------------------!
         expm0_minus(i) = exp( - etai  (i) / mu0 (i) )
         !---------------------------------------------------------------------------------!
      end do directloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Define some boundary conditions for the vector properties above.               !
      !------------------------------------------------------------------------------------!
      i              = ncoh+1
      elai       (i) = 0.d0
      etai       (i) = 0.d0
      leaf_weight(i) = 5.d-1
      wood_weight(i) = 5.d-1
      proj_area  (i) = 5.d-1
      mu0        (i) = czen / proj_area(i)
      mu         (i) = 1.d0
      iota_ratio (i) = 5.d-1 * ( 1.d0 - 5.d-1 * mu0(i)                                     &
                                      * log ( 1.d0  / ( 5.d-1 * mu0(i) ) + 1.d0 ) )
      beta0      (i) = iota_ratio(i) * ( mu0(i) + mu(i) ) / mu(i)
      epsil0     (i) = 1.d0 - 2.d0 * beta0(i)
      expm0_minus(i) = 1.d0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the direct radiation profile using the exponential attenuation curve.     !
      !------------------------------------------------------------------------------------!
      down0(ncoh+1) = down0_sky
      do i = ncoh,1,-1
         ip1      = i + 1
         down0(i) = down0(ip1) * expm0_minus(i) 
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the diffuse radiation properties.                                         !
      !------------------------------------------------------------------------------------!
      diffuseloop: do i=1,ncoh
         !----- Scattering coefficient. ---------------------------------------------------!
         iota      (i) = leaf_weight(i) * leaf_scatter(ipft)                               &
                       + wood_weight(i) * wood_scatter(ipft)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Beta is the backward scattering of diffuse radiation. Epsil = 1 - 2*beta.   !
         !---------------------------------------------------------------------------------!
         beta  (i) = leaf_weight(i) * leaf_backscatter(ipft)                               &
                   + wood_weight(i) * wood_backscatter(ipft)
         epsil (i) = 1.d0 - 2.d0 * beta(i)
         !---------------------------------------------------------------------------------!



         !----- lambda is the coefficient associated with the optical depth. --------------!
         lambda(i) = sqrt( ( 1.d0 - epsil(i) * iota(i) ) * ( 1.d0 - iota(i) ) ) / mu(i)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find some auxiliary variables to determine the right hand side.             !
         !---------------------------------------------------------------------------------!
         a_aux  (i) = - ( ( 1.d0 - epsil(i) * iota(i) ) * iota(i) / mu(i)                  &
                        + epsil0(i) * iota(i) / mu0(i) ) * down0(i+1) / mu0(i)
         s_aux  (i) = - ( ( 1.d0 - iota (i) ) * epsil0(i) * iota(i) / mu(i)                &
                        + iota(i) / mu0(i) ) * down0(i+1) / mu0(i)
         delta  (i) = ( a_aux(i) + s_aux(i) ) * mu0(i) * mu0(i)                            &
                    / ( 2.d0 * ( 1.d0 - lambda(i) * lambda(i) * mu0(i) * mu0(i) ) )
         upsilon(i) = ( a_aux(i) - s_aux(i) ) * mu0(i) * mu0(i)                            &
                    / ( 2.d0 * ( 1.d0 - lambda(i) * lambda(i) * mu0(i) * mu0(i) ) )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    gamm_plus and gamm_minus are the coefficients that relate upwelling and      !
         ! downwelling radiation.                                                          !
         !---------------------------------------------------------------------------------!
         gamm_plus (i) = 5.d-1 * ( 1.d0 + sqrt( ( 1.d0 -            iota(i) )              &
                                              / ( 1.d0 - epsil(i) * iota(i) ) ) )
         gamm_minus(i) = 5.d-1 * ( 1.d0 - sqrt( ( 1.d0 -            iota(i) )              &
                                              / ( 1.d0 - epsil(i) * iota(i) ) ) )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    expl_plus and expl_minus are the transmitivity of diffuse light.             !
         !---------------------------------------------------------------------------------!
         expl_plus  (i) = exp(   lambda(i) * etai(i) )
         expl_minus (i) = exp( - lambda(i) * etai(i) )
         !---------------------------------------------------------------------------------!
      end do diffuseloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Define some boundary conditions for the vector properties above.               !
      !------------------------------------------------------------------------------------!
      i              = ncoh+1
      iota       (i) = 1.d0
      beta       (i) = 0.d0
      epsil      (i) = 1.d0 - 2.d0 * beta(i)
      lambda     (i) = 0.d0
      a_aux      (i) = 0.d0
      s_aux      (i) = 0.d0
      delta      (i) = 0.d0
      upsilon    (i) = 0.d0
      a_aux      (i) = -  epsil0(i) * down0_sky / ( mu0(i) * mu0(i) )
      s_aux      (i) = -  iota  (i) * down0_sky / ( mu0(i) * mu0(i) )
      delta      (i) = 5.d-1 * ( a_aux(i) + s_aux(i) ) * mu0(i) * mu0(i)
      upsilon    (i) = 5.d-1 * ( a_aux(i) - s_aux(i) ) * mu0(i) * mu0(i)
      gamm_plus  (i) = 1.d0
      gamm_minus (i) = 0.d0
      expl_plus  (i) = 1.d0
      expl_minus (i) = 1.d0
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Fill in the right hand side vector (Y) and the matrix (X)                      !
      !------------------------------------------------------------------------------------!
      !------ Initialise the vector and the matrix.  The matrix is a sparse one... --------!
      mmat(:,:)         = 0.d0
      yvec(:)           = 0.d0
      !------ Add the bottom and top boundary conditions. ---------------------------------!
      mmat(1,1)         = (gamm_minus(1) - iota_g * gamm_plus (1)) * expl_minus(1)
      mmat(1,2)         = (gamm_plus (1) - iota_g * gamm_minus(1)) * expl_plus (1)
      mmat(nsiz,nsiz-1) = gamm_plus (ncoh+1)
      mmat(nsiz,nsiz  ) = gamm_minus(ncoh+1)
      yvec(1)           = iota_g * down0(1)                                                &
                        - ( upsilon(1) - iota_g * delta(1) ) * expm0_minus(1)
      yvec(nsiz)        = down_sky - delta(ncoh+1)
      do i=1,ncoh
         !----- Find auxiliary indices. ---------------------------------------------------!
         ip1  =  i + 1
         i2   =  2 * i
         i2m1 = i2 - 1
         i2p1 = i2 + 1
         i2p2 = i2 + 2
         !---------------------------------------------------------------------------------!

         !----- Make elements of C vector. ------------------------------------------------!
         yvec(i2  ) = delta  (ip1) * expm0_minus(ip1) - delta  (i)
         yvec(i2p1) = upsilon(ip1) * expm0_minus(ip1) - upsilon(i)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Make elements of the tri-diagonal A array.                                  !
         !---------------------------------------------------------------------------------!
         mmat(i2  ,i2m1) =   gamm_plus   (i)
         mmat(i2  ,i2  ) =   gamm_minus  (i)
         mmat(i2  ,i2p1) = - gamm_plus (ip1) * expl_minus(ip1)
         mmat(i2  ,i2p2) = - gamm_minus(ip1) * expl_plus (ip1)
         mmat(i2p1,i2m1) =   gamm_minus  (i)
         mmat(i2p1,i2  ) =   gamm_plus   (i)
         mmat(i2p1,i2p1) = - gamm_minus(ip1) * expl_minus(ip1)
         mmat(i2p1,i2p2) = - gamm_plus (ip1) * expl_plus (ip1)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Solve the linear system.  In the future we could use a tridiagonal solver,   !
      ! which is a lot cheaper than the regular Gauss elimination, but for the time being, !
      ! we go with a tested method.                                                        !
      !------------------------------------------------------------------------------------!
      call lisys_solver8(nsiz,mmat,yvec,xvec,sing)
      if (sing) then
         call fatal_error('SW radiation failed... The matrix is singular!'                 &
                         ,'sw_two_stream','twostream_rad.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Copy the solution to the vectors, using the properties:                        !
      !  Un = Un(P)                                                                        !
      !  Dn = Dn(P)                                                                        !
      !------------------------------------------------------------------------------------!
      do i=1,ncoh+1

         !----- Auxiliary indices. --------------------------------------------------------!
         i2   =  2 * i
         i2m1 = i2 - 1
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Retrieve the downward and upward diffuse (hemispheric) radiation, by using !
         ! the solutions for full TAI                                                      !
         !---------------------------------------------------------------------------------!
         down(i) = xvec(i2m1) * gamm_plus (i) * expl_minus (i)                             &
                 + xvec  (i2) * gamm_minus(i) * expl_plus  (i)                             &
                              + delta     (i) * expm0_minus(i)
         up  (i) = xvec(i2m1) * gamm_minus(i) * expl_minus (i)                             &
                 + xvec  (i2) * gamm_plus (i) * expl_plus  (i)                             &
                              + upsilon   (i) * expm0_minus(i)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Integrate the light levels.                                                    !
      !------------------------------------------------------------------------------------!
      do i=1,ncoh
         ip1 = i + 1
         light_level     (i) = light_level     (i) + 5.d-1 * ( down (i) + down (ip1)       &
                                                             + down0(i) + down0(ip1) )
         light_beam_level(i) = light_beam_level(i) + 5.d-1 * ( down0(i) + down0(ip1) )
         light_diff_level(i) = light_diff_level(i) + 5.d-1 * ( down (i) + down (ip1) )
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Save the fluxes that we will use outside this sub-routine.                     !
      !------------------------------------------------------------------------------------!
      select case (iband)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Visible (PAR).                                                              !
         !---------------------------------------------------------------------------------!

         !------ Integrate the visible light levels. --------------------------------------!
         do i=1,ncoh
            ip1 = i + 1
            par_diff_level(i) = 5.d-1 * ( down (i) + down (ip1) ) / par_diff_norm
            par_beam_level(i) = 5.d-1 * ( down0(i) + down0(ip1) ) / par_beam_norm
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Save the radiation fluxes to the output variable.                           !
         !---------------------------------------------------------------------------------!
         do i=1,ncoh
            ip1 = i + 1
            par_beam_flip   (i) = sngloff( down0(ip1) - down0  (i), tiny_offset)
            par_diff_flip   (i) = sngloff( down (ip1) - down   (i)                         &
                                         + up     (i) - up   (ip1), tiny_offset)
            radprof_flip ( 1,i) = sngloff(down0(i), tiny_offset)
            radprof_flip ( 2,i) = 0.0
            radprof_flip ( 3,i) = sngloff(down(i), tiny_offset)
            radprof_flip ( 4,i) = sngloff(up  (i), tiny_offset)
         end do
         !---------------------------------------------------------------------------------!



         !------ Save the fluxes reaching the surface and leaving the top. ----------------!
         dw_parlo_beam = sngloff(down0          (1), tiny_offset)
         dw_parlo_diff = sngloff(down           (1), tiny_offset)
         uw_parhi_diff = sngloff(up        (ncoh+1), tiny_offset)
         !---------------------------------------------------------------------------------!

      case (2)
         !---------------------------------------------------------------------------------!
         !     Near-Infrared (NIR).                                                        !
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Save the radiation fluxes to the output variable.                           !
         !---------------------------------------------------------------------------------!
         do i=1,ncoh
            ip1 = i + 1
            nir_beam_flip   (i) = sngloff( down0(ip1) - down0  (i), tiny_offset)
            nir_diff_flip   (i) = sngloff( down (ip1) - down   (i)                         &
                                         + up     (i) - up   (ip1), tiny_offset)
            radprof_flip ( 5,i) = sngloff( down0  (i), tiny_offset)
            radprof_flip ( 6,i) = 0.0
            radprof_flip ( 7,i) = sngloff( down   (i), tiny_offset)
            radprof_flip ( 8,i) = sngloff( up     (i), tiny_offset)
         end do
         !---------------------------------------------------------------------------------!



         !------ Save the fluxes reaching the surface and leaving the top. ----------------!
         dw_nirlo_beam = sngloff(down0          (1), tiny_offset)
         dw_nirlo_diff = sngloff(down           (1), tiny_offset)
         uw_nirhi_diff = sngloff(up        (ncoh+1), tiny_offset)
         !---------------------------------------------------------------------------------!

      end select
   end do bandloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Normalise the light levels.                                                       !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      light_beam_level(i) = light_beam_level(i) / ( par_beam_norm + nir_beam_norm )
      light_diff_level(i) = light_diff_level(i) / ( par_diff_norm + nir_diff_norm )
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Total normalised radiation.                                                       !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      sw_abs_beam_flip(i) = par_beam_flip(i) + nir_beam_flip(i)
      sw_abs_diff_flip(i) = par_diff_flip(i) + nir_diff_flip(i)
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine sw_two_stream
!==========================================================================================!
!==========================================================================================!
