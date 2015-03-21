!==========================================================================================!
!==========================================================================================!
!     This sub-routine solves the long-wave radiation method using the multiple-scatter-   !
! ing model.  We consider the finite crown area when computing the layer transmittance.    !
! Contrary to short-wave radiation, we don't normalise long-wave radiation.                !
!                                                                                          !
! References:                                                                              !
!                                                                                          !
! Zhao, W., R. J. Qualls, 2006: Modeling of long-wave and net radiation energy             !
!    distribution within a homogeneous plant canopy via multiple scattering processes.     !
!    Water Resources Res., 42, W08435, doi: 10.1029/2005WR004581. (ZQ06)                   !
!                                                                                          !
! Oleson, K. W., and co-authors, 2004: Technical description of the community land model   !
!   (CLM). NCAR Technical note NCAR/TN-461+STR. 186pp. (CLM04)                             !
!                                                                                          !
! Sellers, P. J., 1985: Canopy reflectance, photosynthesis and transpiration. Intl. J.     !
!    of remote sensing, 6, 1335-1372. (S85)                                                !
!                                                                                          !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine lw_multiple_scatter(grnd_emiss4,grnd_temp4,rlong_top4,ncoh,pft,lai,wai,cai      &
                              ,leaf_temp,wood_temp,radprof_flip,tir_flip                   &
                              ,dw_tirlo,uw_tirlo,uw_tirhi)
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
   integer                                                   :: ip1
   integer                                                   :: im1
   integer                                                   :: i2
   integer                                                   :: i2p1
   integer                                                   :: i2m1
   integer                                                   :: i2p2
   logical                                                   :: sing
   real(kind=8), dimension(ncoh)                             :: locetai
   real(kind=8), dimension(ncoh)                             :: elai
   real(kind=8), dimension(ncoh)                             :: etai
   real(kind=8), dimension(ncoh)                             :: tai
   real(kind=8), dimension(ncoh)                             :: leaf_weight
   real(kind=8), dimension(ncoh)                             :: wood_weight
   real(kind=8), dimension(ncoh)                             :: source_lw
   real(kind=8), dimension(0:ncoh+1)                         :: tau      ! tau_i
   real(kind=8), dimension(0:ncoh+1)                         :: omt      ! 1-tau_i
   real(kind=8), dimension(0:ncoh+1)                         :: omr      ! 1-r_i
   real(kind=8), dimension(0:ncoh+1)                         :: r        ! r_i
   real(kind=8), dimension(0:ncoh+1)                         :: epsil    ! epsilon_i
   real(kind=8), dimension(0:ncoh+1)                         :: ome      ! 1-epsilon_i
   real(kind=8), dimension(  ncoh+1)                         :: lwd0
   real(kind=8), dimension(0:ncoh  )                         :: lwu0
   real(kind=8), dimension(  ncoh+1)                         :: lwd
   real(kind=8), dimension(0:ncoh  )                         :: lwu
   real(kind=8), dimension(2*ncoh+2)                         :: lwvec
   real(kind=8), dimension(2*ncoh+2)                         :: cvec
   real(kind=8), dimension(2*ncoh+2,2*ncoh+2)                :: amat
   real(kind=8)                                              :: temiss_four
   real(kind=8)                                              :: grnd_emiss
   real(kind=8)                                              :: grnd_temp
   real(kind=8)                                              :: rlong_top
   real(kind=8)                                              :: ext_diff1
   real(kind=8)                                              :: ext_diff2
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
   grnd_emiss = dble(grnd_emiss4)
   grnd_temp  = dble(grnd_temp4 )
   rlong_top  = dble(rlong_top4 )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert some leaf variables to double precision, then find some general           !
   ! properties of leaves/branches (i.e., properties that do not depend on which band we   !
   ! are solving):  the extinction coefficient lambda, the layer transmittance             !
   ! coefficients tau_beam and tau_diff, and the weight of LAI and WAI relative to the     !
   ! total area.                                                                           !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      ipft      = pft(i)
      !----- Area indices, correct LAI by clumping factor. --------------------------------!
      elai(i) = clumping_factor(ipft) * lai(i)
      tai(i)  =  lai(i) + wai(i)
      etai(i) = elai(i) + wai(i)
      !------------------------------------------------------------------------------------!



      !----- LOCETAI is the local exposed LAI. --------------------------------------------!
      locetai(i) = etai(i) / cai(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    The weighting factors are defined by the relative LAI and WAI areas.            !
      !------------------------------------------------------------------------------------!
      leaf_weight(i)  =  lai(i) / tai(i)
      wood_weight(i)  = 1.d0 -  leaf_weight(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The layer transmittance coefficient for diffuse radiation is defined in a      !
      ! similar way from equation (2) of ZQ06.  We must integrate the contribution coming  !
      ! from each hemispheric direction.  Again, we must consider the finite crown area,   !
      ! and use the same extinction function we use for shortwave radiation, which comes   !
      ! from CLM04 equation 3.3.                                                           !
      !     The integral can be solved analytically for our case, with the help of the     !
      ! handy http://www.integrals.com website, of course ;-).                             !
      !------------------------------------------------------------------------------------!
      ext_diff1   = phi1(ipft) * locetai(i)
      ext_diff2   = phi2(ipft) * locetai(i)
      tau(i)      = (1.d0 - cai(i))                                                        &
                  - cai(i) * exp(- ext_diff1 - ext_diff2)                                  &
                  * ( ext_diff1 * ext_diff1 * exp(ext_diff1) * eifun8(-ext_diff1)          &
                    + (ext_diff1 - 1.d0) )
      !------------------------------------------------------------------------------------!



      !----- Backward scattering of diffuse radiation. ------------------------------------!
      r(i)     = leaf_weight(i) * leaf_backscatter_tir(ipft)                               &
               + wood_weight(i) * wood_backscatter_tir(ipft)
      !------------------------------------------------------------------------------------!



      !----- Layer emissivity. ------------------------------------------------------------!
      epsil(i) = leaf_weight(i) * leaf_emiss_tir(ipft)                                     &
               + wood_weight(i) * wood_emiss_tir(ipft)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Thermal emission.  The temperature is the weighted average where the weights  !
      ! are the product between the area and the emissivity.                               !
      !------------------------------------------------------------------------------------!
      temiss_four  = ( leaf_weight(i) * leaf_emiss_tir(ipft) * leaf_temp(i) ** 4           &
                     + wood_weight(i) * wood_emiss_tir(ipft) * wood_temp(i) ** 4 )         &
                   / ( leaf_weight(i) * leaf_emiss_tir(ipft)                               &
                     + wood_weight(i) * wood_emiss_tir(ipft) )
      source_lw(i) = epsil(i) * stefan8 * temiss_four
      !------------------------------------------------------------------------------------!


   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Define some boundary conditions for transmission, emission, and backscattering,   !
   ! and the complement of the three properties.  The boundary conditions come from ZQ06,  !
   ! equations 16-22.                                                                      !
   !---------------------------------------------------------------------------------------!
   tau       (0) = 0.d0
   tau  (ncoh+1) = 1.d0
   r         (0) = 1.d0
   r    (ncoh+1) = 0.d0
   epsil     (0) = grnd_emiss
   epsil(ncoh+1) = 0.d0
   omr(:)        = 1.d0 - r    (:)
   omt(:)        = 1.d0 - tau  (:)
   ome(:)        = 1.d0 - epsil(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Find the top diffuse radiation, using equations (19) and (22) from ZQ06.           !
   !---------------------------------------------------------------------------------------!
   lwd0     (ncoh+1) = rlong_top
   lwu0     (     0) = epsil(0) * stefan8 * grnd_temp ** 4
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fill in the right hand side vector (C) and the matrix (A) described by            !
   ! equations (25) and ("25 1/2") of ZQ06.                                                !
   !---------------------------------------------------------------------------------------!
   !------ Initialise the vector and the matrix.  The matrix is a sparse one... -----------!
   amat(:,:)       = 0.d0
   cvec(:)         = 0.d0
   !------ Add the edge values first. -----------------------------------------------------!
   amat(1,1)       = 1.d0
   amat(nsiz,nsiz) = 1.d0
   cvec(1)         = lwu0(0)
   cvec(nsiz)      = lwd0(ncoh+1)
   do i=1,ncoh
      !----- Find auxiliary indices. ------------------------------------------------------!
      ip1  =  i + 1
      im1  =  i - 1
      i2   =  2 * i
      i2p1 = i2 + 1
      i2m1 = i2 - 1
      i2p2 = i2 + 2
      !------------------------------------------------------------------------------------!

      !----- Make elements of C vector. ---------------------------------------------------!
      cvec(i2)   = ( 1.d0 - r(im1) * ome(im1) * omt(im1) * r(i)   * ome(i)   * omt(i))     &
                   * omt(i) * source_lw(i)                                                
      cvec(i2p1) = ( 1.d0 - r(i)   * ome(i)   * omt(i) * r(ip1) * ome(ip1) * omt(ip1))     &
                   * omt(i) * source_lw(i)                                                
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Make elements of the tri-diagonal A array.                                     !
      !------------------------------------------------------------------------------------!
      amat(i2  ,i2m1) = - (tau(i) + omr(i) * ome(i) * omt(i))
      amat(i2  ,i2  ) = - r(im1) * ome(im1) * omt(im1) * (tau(i) + omt(i) * ome(i) * omr(i))
      amat(i2  ,i2p1) = 1.d0 - r(im1) * ome(im1) * omt(im1) * r(i)   * ome(i)   * omt(i)
      amat(i2p1,i2  ) = 1.d0 - r(i)   * ome(i)   * omt(i) * r(ip1) * ome(ip1) * omt(ip1)
      amat(i2p1,i2p1) = - r(ip1) * ome(ip1) * omt(ip1) * (tau(i) + omt(i) * ome(i) * omr(i))
      amat(i2p1,i2p2) = - (tau(i) + omt(i) * ome(i) * omr(i))
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !       Solve the linear system.  In the future we could use a tridiagonal solver,      !
   ! which is a lot cheaper than the regular Gauss elimination, but for the time being, we !
   ! go with a tested method.                                                              !
   !---------------------------------------------------------------------------------------!
   call lisys_solver8(nsiz,amat,cvec,lwvec,sing)
   if (sing) then
      call fatal_error('LW radiation failed... The matrix is singular!'                    &
                      ,'lw_multiple_scatter','multiple_scatter.f90')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Copy the solution without considering the multiple scattering to the vectors.     !
   !---------------------------------------------------------------------------------------!
   lwu0(0)      = lwvec(1)
   lwd0(ncoh+1) = lwvec(nsiz)
   do i=1,ncoh

      !----- Auxiliary indices. -----------------------------------------------------------!
      i2   = i * 2
      i2p1 = i * 2 + 1
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Retrieve the downward and upward diffuse (hemispheric) radiation, using       !
      ! equation (24) from ZQ06.                                                           !
      !------------------------------------------------------------------------------------!
      lwd0(i) = lwvec(i2)
      lwu0(i) = lwvec(i2p1)
      !------------------------------------------------------------------------------------!

   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the diffuse radiation when multiple scattering is considered.  This is based !
   ! on equations (8) and (9) of ZQ06, with the difference that we find LWd_i instead of   !
   ! LWd_i+1.                                                                              !
   !---------------------------------------------------------------------------------------!
   !----- Downward diffuse (hemispheric) radiation. ---------------------------------------!
   do i=1,ncoh+1
      im1 = i - 1
      ip1 = i + 1
      lwd(i) = ( lwd0(i) + r(i) * ome(i) * omt(i) * lwu0(im1) )                            &
             / ( 1.d0 - r(im1) * ome(im1) * omt(im1) * r(i)   * ome(i)   * omt(i)   )
   end do
   !----- Upward diffuse (hemispheric) radiation. -----------------------------------------!
   do i=0,ncoh
      im1 = i - 1
      ip1 = i + 1
      lwu(i) = ( lwu0(i) + r(i) * ome(i) * omt(i) * lwd0(ip1) )                            &
             / ( 1.d0 - r(i)   * ome(i)   * omt(i) * r(ip1) * ome(ip1) * omt(ip1)   )
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Save the fluxes that we will use outside this sub-routine.                        !
   !---------------------------------------------------------------------------------------!



   !------ Save the fluxes reaching the surface and leaving the top. ----------------------!
   dw_tirlo = sngloff(lwd          (1), tiny_offset)
   uw_tirlo = sngloff(lwu          (0), tiny_offset)
   uw_tirhi = sngloff(lwu       (ncoh), tiny_offset)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Save the radiation fluxes to the output variable.                                 !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      im1 = i - 1
      ip1 = i + 1
      radprof_flip( 9,i) = sngloff(lwd  (i), tiny_offset)
      radprof_flip(10,i) = sngloff(lwu(im1), tiny_offset)
      tir_flip       (i) = sngloff(lwd(ip1) - lwd(i) + lwu(im1) - lwu(i), tiny_offset)
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine lw_multiple_scatter
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine solves the short-wave radiation method using the multiple-scatter-  !
! ing model.  We consider the finite crown area when computing the transmittance of both   !
! direct and diffuse radiation.  References:                                               !
!                                                                                          !
! Zhao, W., R. J. Qualls, 2005: A multiple-layer canopy scattering model to simulate       !
!    shortwave radiation distribution within a homogeneous plant canopy.  Water Resources  !
!    Res., 41, W08409, doi: 10.1029/2005WR004016. (ZQ05)                                   !
!                                                                                          !
! Oleson, K. W., and co-authors, 2004: Technical description of the community land model   !
!   (CLM). NCAR Technical note NCAR/TN-461+STR. 186pp. (CLM04)                             !
!                                                                                          !
! Sellers, P. J., 1985: Canopy reflectance, photosynthesis and transpiration. Intl. J.     !
!    of remote sensing, 6, 1335-1372. (S85)                                                !
!                                                                                          !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine sw_multiple_scatter(grnd_alb_par4,grnd_alb_nir4,cosaoi4,ncoh,pft,lai,wai,cai    &
                              ,radprof_flip,par_beam_flip,par_diff_flip,sw_abs_beam_flip   &
                              ,sw_abs_diff_flip,dw_parlo_beam,dw_parlo_diff                &
                              ,uw_parhi_diff,dw_nirlo_beam,dw_nirlo_diff,uw_nirhi_diff     &
                              ,par_level_beam,par_level_diffd,par_level_diffu              &
                              ,light_level,light_beam_level  &
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

   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4)                              , intent(in)    :: grnd_alb_par4
   real(kind=4)                              , intent(in)    :: grnd_alb_nir4
   real(kind=4)                              , intent(in)    :: cosaoi4
   integer                                   , intent(in)    :: ncoh
   integer     , dimension(ncoh)             , intent(in)    :: pft
   real(kind=8), dimension(ncoh)             , intent(in)    :: lai
   real(kind=8), dimension(ncoh)             , intent(in)    :: wai
   real(kind=8), dimension(ncoh)             , intent(in)    :: cai
   real(kind=4), dimension(n_radprof,ncoh)   , intent(inout) :: radprof_flip
   real(kind=4), dimension(ncoh)             , intent(out)   :: par_beam_flip
   real(kind=4), dimension(ncoh)             , intent(out)   :: par_diff_flip
   real(kind=4), dimension(ncoh)             , intent(out)   :: sw_abs_beam_flip
   real(kind=4), dimension(ncoh)             , intent(out)   :: sw_abs_diff_flip
   real(kind=4)                              , intent(out)   :: uw_parhi_diff
   real(kind=4)                              , intent(out)   :: uw_nirhi_diff
   real(kind=4)                              , intent(out)   :: dw_parlo_beam
   real(kind=4)                              , intent(out)   :: dw_parlo_diff
   real(kind=4)                              , intent(out)   :: dw_nirlo_beam
   real(kind=4)                              , intent(out)   :: dw_nirlo_diff
   real(kind=8), dimension(ncoh)             , intent(out)   :: par_level_beam
   real(kind=8), dimension(ncoh)             , intent(out)   :: par_level_diffu
   real(kind=8), dimension(ncoh)             , intent(out)   :: par_level_diffd
   real(kind=8), dimension(ncoh)             , intent(out)   :: light_level
   real(kind=8), dimension(ncoh)             , intent(out)   :: light_beam_level
   real(kind=8), dimension(ncoh)             , intent(out)   :: light_diff_level
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                 :: ipft
   integer                                                 :: iband
   integer                                                 :: nsiz
   integer                                                 :: i
   integer                                                 :: ip1
   integer                                                 :: im1
   integer                                                 :: i2
   integer                                                 :: i2p1
   integer                                                 :: i2m1
   integer                                                 :: i2p2
   logical                                                 :: sing
   real(kind=8)                                            :: alb_par
   real(kind=8)                                            :: alb_nir
   real(kind=8)                                            :: mu
   real(kind=4), dimension(ncoh)                           :: nir_beam_flip
   real(kind=4), dimension(ncoh)                           :: nir_diff_flip
   real(kind=8), dimension(ncoh)                           :: locetai
   real(kind=8), dimension(ncoh)                           :: etai
   real(kind=8), dimension(ncoh)                           :: elai
   real(kind=8), dimension(ncoh)                           :: leaf_weight
   real(kind=8), dimension(ncoh)                           :: wood_weight
   real(kind=8), dimension(ncoh)                           :: beam_backscatter
   real(kind=8), dimension(ncoh)                           :: lambda
   real(kind=8), dimension(0:ncoh+1)                       :: tau_beam ! tau_i(psi)
   real(kind=8), dimension(0:ncoh+1)                       :: tau_diff ! tau_i
   real(kind=8), dimension(0:ncoh+1)                       :: omt_beam ! 1-tau_i(psi)
   real(kind=8), dimension(0:ncoh+1)                       :: omt_diff ! 1-tau_i
   real(kind=8), dimension(0:ncoh+1)                       :: omr_beam ! 1-r_i(psi)
   real(kind=8), dimension(0:ncoh+1)                       :: omr_diff ! 1-r_i
   real(kind=8), dimension(0:ncoh+1)                       :: r_beam   ! r_i(psi)
   real(kind=8), dimension(0:ncoh+1)                       :: r_diff   ! r_i
   real(kind=8), dimension(0:ncoh+1)                       :: alpha    ! alpha_i
   real(kind=8), dimension(0:ncoh+1)                       :: oma      ! 1-alpha_i
   real(kind=8), dimension(  ncoh+1)                       :: beam_down
   real(kind=8), dimension(  ncoh+1)                       :: swd0
   real(kind=8), dimension(0:ncoh  )                       :: swu0
   real(kind=8), dimension(  ncoh+1)                       :: swd
   real(kind=8), dimension(0:ncoh  )                       :: swu
   real(kind=8), dimension(2*ncoh+2)                       :: swvec
   real(kind=8), dimension(2*ncoh+2)                       :: cvec
   real(kind=8), dimension(2*ncoh+2,2*ncoh+2)              :: amat
   real(kind=8)                                            :: proj_area
   real(kind=8)                                            :: ext_diff1
   real(kind=8)                                            :: ext_diff2
   real(kind=8)                                            :: snglscat_alb
   !----- External functions. -------------------------------------------------------------!
   real(kind=8)                              , external    :: eifun8
   real(kind=4)                              , external    :: sngloff
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the size of the solution matrix and array.                                   !
   !---------------------------------------------------------------------------------------!
   nsiz = 2*ncoh + 2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert cosine of angle of incidence to double precision.                         !
   !---------------------------------------------------------------------------------------!
   mu      = max(cosz_min8,dble(cosaoi4))
   alb_par = dble(grnd_alb_par4)
   alb_nir = dble(grnd_alb_nir4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert some leaf variables to double precision, then find some general           !
   ! properties of leaves/branches (i.e., properties that do not depend on which band we   !
   ! are solving):  the extinction coefficient lambda, the layer transmittance             !
   ! coefficients tau_beam and tau_diff, and the weight of LAI and WAI relative to the     !
   ! total area.                                                                           !
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      ipft      = pft(i)
      !----- Area indices, correct LAI by clumping factor. --------------------------------!
      elai(i) = clumping_factor(ipft) * lai(i)
      etai(i) = elai(i) + wai(i)
      !------------------------------------------------------------------------------------!



      !----- LOCETAI is the local exposed LAI. --------------------------------------------!
      locetai(i) = etai(i) / cai(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    The weighting factors are defined by the relative LAI and WAI areas.            !
      !------------------------------------------------------------------------------------!
      leaf_weight(i) = elai(i) / etai(i)
      wood_weight(i) = 1.d0 - leaf_weight(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     We find the optical depth of the direct beam (lambda), following CLM04         !
      ! (equation 3.3 and text after equation 3.2).                                        !
      !------------------------------------------------------------------------------------!
      proj_area     = phi1(ipft) + phi2(ipft) * mu
      lambda    (i) = proj_area / mu
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The layer transmittance coefficient for direct radiation is defined from       !
      ! equation (1) of ZQ05, with the difference that the crown area is taken into        !
      ! account.                                                                           !
      !------------------------------------------------------------------------------------!
      tau_beam(i) = (1.d0 - cai(i)) + cai(i) * exp (- lambda(i) * locetai(i))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     The layer transmittance coefficient for diffuse radiation is defined in a      !
      ! similar way from equation (2) of ZQ06.  We must integrate the contribution coming  !
      ! from each hemispheric direction.  Again, we must consider the finite crown area.   !
      !     The integral can be solved analytically for our case, with the help of the     !
      ! handy http://www.integrals.com website, of course ;-).                             !
      !------------------------------------------------------------------------------------!
      ext_diff1   = phi1(ipft) * locetai(i)
      ext_diff2   = phi2(ipft) * locetai(i)
      tau_diff(i) = (1.d0 - cai(i))                                                        &
                  - cai(i) * exp(- ext_diff1 - ext_diff2)                                  &
                  * ( ext_diff1 * ext_diff1 * exp(ext_diff1) * eifun8(-ext_diff1)          &
                    + (ext_diff1 - 1.d0) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the backscatter for direct beam radiation, following CLM04, equations    !
      ! (3.14) and (3.15).  The forward scattering was omitted on both equations because   !
      ! they are different for PAR and NIR, but in both cases they would cancel out.       !
      !------------------------------------------------------------------------------------!
      snglscat_alb        = 5.d-1 * proj_area / (phi2(ipft) * mu + proj_area)              &
                          * ( 1.d0 - phi1(ipft) * mu / (phi2(ipft) * mu + proj_area)       &
                            * log (1.d0  + (phi2(ipft) * mu + proj_area)                   &
                                         / (phi1(ipft) * mu)))
      beam_backscatter(i) = ( 1.d0 + mu_bar(ipft) * lambda(i) ) * snglscat_alb             &
                             / ( mu_bar(ipft) * lambda(i) )
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Define some boundary conditions for transmission, and the complement of           !
   ! transmission.  The boundary conditions come from ZQ05, equations 32-37.               !
   !---------------------------------------------------------------------------------------!
   tau_beam(0)      = 0.d0
   tau_beam(ncoh+1) = 1.d0
   tau_diff(0)      = 0.d0
   tau_diff(ncoh+1) = 1.d0
   omt_beam(:)      = 1.d0 - tau_beam(:)
   omt_diff(:)      = 1.d0 - tau_diff(:)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the light level variables.                                             !
   !---------------------------------------------------------------------------------------!
   light_level     (:) = 0.d0
   light_beam_level(:) = 0.d0
   light_diff_level(:) = 0.d0
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Initialise the PAR level variables.                                             !
   !---------------------------------------------------------------------------------------!
   par_level_diffu     (:) = 0.d0
   par_level_diffd     (:) = 0.d0
   par_level_beam      (:) = 0.d0
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Here we loop over the two bands (band 1 is PAR and band 2 is NIR), following the  !
   ! algorithm proposed by ZQ05.                                                           !
   !---------------------------------------------------------------------------------------!
   bandloop: do iband = 1,2
      select case (iband)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Visible (PAR).                                                              !
         !---------------------------------------------------------------------------------!
         do i=1,ncoh
            !----- Alias for PFT of this layer. -------------------------------------------!
            ipft = pft(i)

            !----- Absorptance: the fraction that is not scattered. -----------------------!
            alpha(i)    = 1.d0 - leaf_weight(i) * leaf_scatter_vis(ipft)                   &
                               - wood_weight(i) * wood_scatter_vis(ipft)

            !----- Backward scattering of diffuse radiation. ------------------------------!
            r_diff(i)   = leaf_weight(i) * leaf_backscatter_vis(ipft)                      &
                        + wood_weight(i) * wood_backscatter_vis(ipft)

            !----- Backward scattering of direct radiation. -------------------------------!
            r_beam(i)   = beam_backscatter(i)

         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the boundary conditions following ZQ05, equations (32-37).             !
         !---------------------------------------------------------------------------------!
         alpha (0)      = 1.d0 - alb_par
         alpha (ncoh+1) = 0.d0
         r_beam(0)      = 1.d0
         r_diff(0)      = 1.d0
         r_beam(ncoh+1) = 0.d0
         r_diff(ncoh+1) = 0.d0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the top beam and diffuse radiation.  For the direct radiation, we       !
         ! follow a notation that is slightly different from ZQ05 notation, in which       !
         ! the beam radiation at level i is the amount that reaches the bottom of the      !
         ! layer.                                                                          !
         !---------------------------------------------------------------------------------!
         beam_down(ncoh+1) = par_beam_norm
         swd0     (ncoh+1) = par_diff_norm
         !---------------------------------------------------------------------------------!


      case (2)
         !---------------------------------------------------------------------------------!
         !     Near-Infrared (NIR).                                                        !
         !---------------------------------------------------------------------------------!
         do i=1,ncoh
            !----- Alias for PFT of this layer. -------------------------------------------!
            ipft = pft(i)

            !----- Absorptance: the fraction that is not scattered. -----------------------!
            alpha(i)    = 1.d0 - leaf_weight(i) * leaf_scatter_nir(ipft)                   &
                               - wood_weight(i) * wood_scatter_nir(ipft)

            !----- Backward scattering of diffuse radiation. ------------------------------!
            r_diff(i)   = leaf_weight(i) * leaf_backscatter_nir(ipft)                      &
                        + wood_weight(i) * wood_backscatter_nir(ipft)

            !----- Backward scattering of direct radiation. -------------------------------!
            r_beam(i)   = beam_backscatter(i)

         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the boundary conditions following ZQ05, equations (32-37).             !
         !---------------------------------------------------------------------------------!
         alpha (0)      = 1.d0 - alb_nir
         alpha (ncoh+1) = 0.d0
         r_beam(0)      = 1.d0
         r_diff(0)      = 1.d0
         r_beam(ncoh+1) = 0.d0
         r_diff(ncoh+1) = 0.d0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the top beam and diffuse radiation.  For the direct radiation, we       !
         ! follow a notation that is slightly different from ZQ05 notation, in which       !
         ! the beam radiation at level i is the amount that reaches the bottom of the      !
         ! layer.                                                                          !
         !---------------------------------------------------------------------------------!
         beam_down(ncoh+1) = nir_beam_norm
         swd0     (ncoh+1) = nir_diff_norm
         !---------------------------------------------------------------------------------!
      end select



      !----- Find the complement of the absorptance and scattering. -----------------------!
      oma     (:) = 1.d0 - alpha (:)
      omr_diff(:) = 1.d0 - r_diff(:)
      omr_beam(:) = 1.d0 - r_beam(:)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the profile of downward direct radiation using the typical atenuation      !
      ! function.  We define the beam_down vector in a slightly different way from ZQ05,   !
      ! because here we will define S_i as the radiation that reaches the bottom of the    !
      ! layer, whereas they defined it as the radiation that reaches the top.              !
      !------------------------------------------------------------------------------------!
      do i=ncoh,1,-1
         beam_down(i) = beam_down(i+1) * tau_beam(i)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Fill in the right hand side vector (C) and the matrix (A) described by         !
      ! equations (41) and (42) of ZQ05.                                                   !
      !------------------------------------------------------------------------------------!
      !------ Initialise the vector and the matrix.  The matrix is a sparse one... --------!
      amat(:,:)       = 0.d0
      cvec(:)         = 0.d0
      !------ Add the edge values first. --------------------------------------------------!
      amat(1,1)       = 1.d0
      amat(nsiz,nsiz) = 1.d0
      cvec(1)         = oma(0) * beam_down(1)
      cvec(nsiz)      = swd0(ncoh+1)
      do i=1,ncoh
         !----- Find auxiliary indices. ---------------------------------------------------!
         ip1  =  i + 1
         im1  =  i - 1
         i2   =  2 * i
         i2p1 = i2 + 1
         i2m1 = i2 - 1
         i2p2 = i2 + 2
         !---------------------------------------------------------------------------------!

         !----- Make elements of C vector. ------------------------------------------------!
         cvec(i2)      = ( 1.d0                                                            &
                         - r_diff(im1) * oma(im1) * omt_diff(im1)                          &
                         * r_diff(i)   * oma(i)   * omt_diff(i)   )                        &
                         * oma(i) * omt_beam(i) * r_beam(i) * beam_down(ip1)
         cvec(i2p1)    = ( 1.d0                                                            &
                         - r_diff(i)   * oma(i)   * omt_diff(i)                            &
                         * r_diff(ip1) * oma(ip1) * omt_diff(ip1)   )                      &
                         * oma(i) * omt_beam(i) * omr_beam(i) * beam_down(ip1)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Make elements of the tri-diagonal A array.                                  !
         !---------------------------------------------------------------------------------!
         amat(i2  ,i2m1)  = - (tau_diff(i) + omt_diff(i) * oma(i) * omr_diff(i))
         amat(i2  ,i2  )  = - r_diff(im1) * oma(im1) * omt_diff(im1)                       &
                          * (tau_diff(i) + omt_diff(i) * oma(i) * omr_diff(i) )
         amat(i2  ,i2p1)  = 1.d0 - r_diff(im1) * oma(im1) * omt_diff(im1)                  &
                                 * r_diff(i)   * oma(i)   * omt_diff(i)
         amat(i2p1,i2  )  = 1.d0 - r_diff(i)   * oma(i)   * omt_diff(i)                    &
                                 * r_diff(ip1) * oma(ip1) * omt_diff(ip1)
         amat(i2p1,i2p1)  = - r_diff(ip1) * oma(ip1) * omt_diff(ip1)                       &
                            * (tau_diff(i) + omt_diff(i) * oma(i) * omr_diff(i))
         amat(i2p1,i2p2)  = - (tau_diff(i) + omt_diff(i) * oma(i) * omr_diff(i))
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Solve the linear system.  In the future we could use a tridiagonal solver,   !
      ! which is a lot cheaper than the regular Gauss elimination, but for the time being, !
      ! we go with a tested method.                                                        !
      !------------------------------------------------------------------------------------!
      call lisys_solver8(nsiz,amat,cvec,swvec,sing)
      if (sing) then
         call fatal_error('SW radiation failed... The matrix is singular!'                 &
                         ,'sw_multiple_scatter','multiple_scatter.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Copy the solution without considering the multiple scattering to the vectors.  !
      !------------------------------------------------------------------------------------!
      swu0(0)      = swvec(1)
      swd0(ncoh+1) = swvec(nsiz)
      do i=1,ncoh

         !----- Auxiliary indices. --------------------------------------------------------!
         i2   = i * 2
         i2p1 = i * 2 + 1
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Retrieve the downward and upward diffuse (hemispheric) radiation, using    !
         ! equation (40) from ZQ05.                                                        !
         !---------------------------------------------------------------------------------!
         swd0(i) = swvec(i2)
         swu0(i) = swvec(i2p1)
         !---------------------------------------------------------------------------------!

      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the diffuse radiation when multiple scattering is considered.  This is    !
      ! based on equations (24) and (25) of ZQ05, with the difference that we find SWd_i   !
      ! instead of SWd_i+1.                                                                !
      !------------------------------------------------------------------------------------!
      !----- Downward diffuse (hemispheric) radiation. ------------------------------------!
      do i=1,ncoh+1
         im1 = i - 1
         ip1 = i + 1
         swd(i) = ( swd0(i) + r_diff(i) * oma(i) * omt_diff(i) * swu0(im1) )               &
                / ( 1.d0 - r_diff(im1) * oma(im1) * omt_diff(im1)                          &
                         * r_diff(i)   * oma(i)   * omt_diff(i)   )
      end do
      !----- Upward diffuse (hemispheric) radiation. --------------------------------------!
      do i=0,ncoh
         im1 = i - 1
         ip1 = i + 1
         swu(i) = ( swu0(i) + r_diff(i) * oma(i) * omt_diff(i) * swd0(ip1) )               &
                / ( 1.d0 - r_diff(i)   * oma(i)   * omt_diff(i)                            &
                         * r_diff(ip1) * oma(ip1) * omt_diff(ip1)   )
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Integrate the light levels.                                                    !
      !------------------------------------------------------------------------------------!
      do i=1,ncoh
         ip1 = i + 1
         light_level     (i) = light_level     (i)                                         &
                             + 5.d-1 * (swd(i) + swd(ip1) + beam_down(i) + beam_down(ip1))
         light_beam_level(i) = light_beam_level(i) + 5.d-1 * (beam_down(i) + beam_down(ip1))
         light_diff_level(i) = light_diff_level(i) + 5.d-1 * (swd(i) + swd(ip1))
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
         ! NEEDS TO BE CHECKED (PARTICULARLY THE UPWARD)
         ! THIS SHOULD BE THE LEVEL (COHORT) CENTERED FLUX OF PAR
         do i=1,ncoh
            ip1 = i + 1
            im1 = i - 1
            par_level_diffd(i) = 5.d-1 * (swd(i) + swd(ip1)) / (par_diff_norm + par_beam_norm)
            par_level_diffu(i) = 5.d-1 * (swu(i) + swu(im1)) / (par_diff_norm + par_beam_norm)
            par_level_beam (i) = 5.d-1 * (beam_down(i) + beam_down(ip1)) / (par_diff_norm+par_beam_norm)
         end do
         !---------------------------------------------------------------------------------!



         !------ Save the fluxes reaching the surface and leaving the top. ----------------!
         dw_parlo_beam = sngloff(beam_down      (1), tiny_offset)
         dw_parlo_diff = sngloff(swd            (1), tiny_offset)
         uw_parhi_diff = sngloff(swu         (ncoh), tiny_offset)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Save the radiation fluxes to the output variable.                           !
         !---------------------------------------------------------------------------------!
         do i=1,ncoh
            im1 = i - 1
            ip1 = i + 1
            par_beam_flip  (i) = sngloff(beam_down(ip1) - beam_down(i)        ,tiny_offset)
            par_diff_flip  (i) = sngloff(swd(ip1) - swd(i) + swu(im1) - swu(i),tiny_offset)
            radprof_flip (1,i) = sngloff(beam_down  (i),tiny_offset)
            radprof_flip (2,i) = 0.0
            radprof_flip (3,i) = sngloff(swd        (i),tiny_offset)
            radprof_flip (4,i) = sngloff(swu      (im1),tiny_offset)
         end do
         !---------------------------------------------------------------------------------!

      case (2)
         !---------------------------------------------------------------------------------!
         !     Near-Infrared (NIR).                                                        !
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Save the radiation fluxes to the output variable.                           !
         !---------------------------------------------------------------------------------!
         do i=1,ncoh
            im1 = i - 1
            ip1 = i + 1
            nir_beam_flip  (i) = sngloff(beam_down(ip1) - beam_down(i)        ,tiny_offset)
            nir_diff_flip  (i) = sngloff(swd(ip1) - swd(i) + swu(im1) - swu(i),tiny_offset)
            radprof_flip (5,i) = sngloff(beam_down  (i),tiny_offset)
            radprof_flip (6,i) = 0.0
            radprof_flip (7,i) = sngloff(swd        (i),tiny_offset)
            radprof_flip (8,i) = sngloff(swu      (im1),tiny_offset)
         end do
         !---------------------------------------------------------------------------------!



         !------ Save the fluxes reaching the surface and leaving the top. ----------------!
         dw_nirlo_beam = sngloff(beam_down      (1), tiny_offset)
         dw_nirlo_diff = sngloff(swd            (1), tiny_offset)
         uw_nirhi_diff = sngloff(swu         (ncoh), tiny_offset)
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
end subroutine sw_multiple_scatter
!==========================================================================================!
!==========================================================================================!
