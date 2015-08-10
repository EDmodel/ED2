!==========================================================================================!
!==========================================================================================!
!     Some constants, happening mostly at rk4_derivs.f90.  These used to be at ed_commons, !
! moved it here so they stay together. Still need some explanation of what these variables !
! represent.                                                                               !
!------------------------------------------------------------------------------------------!
module canopy_air_coms

   use consts_coms, only : twothirds & ! intent(in)
                         , vonk      ! ! intent(in)

   !=======================================================================================!
   !=======================================================================================!
   !     Parameter that will be read in the namelist.                                      !
   !---------------------------------------------------------------------------------------!
   integer :: icanturb      ! Which canopy turbulence we will use?
                            ! -1. is the original ED-2.0 scheme
                            !  0. is the default scheme
                            !  1. will rescale reference wind-speed if it
                            !     the reference height is (inapropriately)
                            !     below the canopy
                            !  2. (recommended) uses the method of Massman 
                            !     1997 and the bulk Richardson number of 
                            !     instability.  This method will not work 
                            !     when zref<h

   integer :: isfclyrm      ! Surface layer model (used to compute ustar, tstar,...)
                            !  1 - Louis, 1979: Boundary-Layer Meteor., 17, 187-202.
                            !      This is the ED-2.0 default, also used in (B)RAMS
                            !  2 - Oncley and Dudhia, 1995: Mon. Wea. Rev., 123, 3344-3357.
                            !      This is used in MM5 and WRF.
                            !  3 - Beljaars and Holtslag, 1991: J. Appl. Meteor., 30, 
                            !      328-341.
                            !  4 - BH91, using OD95 to find zeta.


   integer :: ied_grndvap   ! Methods to find the ground -> canopy conductance:
                            ! In all cases the beta term is modified so it approaches
                            !    zero as soil moisture goes to dry air soil. 
                            !  0. Modified Lee Pielke (1992), adding field capacity, but 
                            !     using beta factor without the square, like in 
                            !     Noilhan and Planton (1989).  This is the closest to
                            !     the original ED-2.0
                            !  1. Test # 1 of Mahfouf and Noilhan (1991)
                            !  2. Test # 2 of Mahfouf and Noilhan (1991)
                            !  3. Test # 3 of Mahfouf and Noilhan (1991)
                            !  4. Test # 4 of Mahfouf and Noilhan (1991)

   real :: leaf_maxwhc      !   Maximum amount of water that can stay on top of the leaf
                            ! surface.  If this amount is reached, the leaf stops collect-
                            ! ing water, thus increasing the through fall fraction.  This 
                            ! value is given in kg/[m2 leaf], so it will be always scaled
                            ! by LAI.
   !---------------------------------------------------------------------------------------!

   !----- Minimum speed for stars [m/s]. --------------------------------------------------!
   real         :: ubmin
   !----- Minimum speed for conductances [m/s]. -------------------------------------------!
   real         :: ugbmin
   !----- Minimum Ustar [m/s]. ------------------------------------------------------------!
   real         :: ustmin
   !----- Used by OD95 and BH91. ----------------------------------------------------------!
   real   :: gamm        ! Gamma for momentum.
   real   :: gamh        ! Gamma for heat.
   real   :: tprandtl    ! Turbulent Prandtl number.
   real   :: vh2vr       ! vegetation roughness:vegetation height ratio
   real   :: vh2dh       ! displacement height:vegetation height ratio
   real   :: ribmax      ! Maximum bulk Richardson number (ignored when ISFCLYRM = 1)
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!
   !     Parameters used in Euler and Runge-Kutta.                                         !
   !---------------------------------------------------------------------------------------!
   !----- Exponential wind atenuation factor [dimensionless]. -----------------------------!
   real         :: exar
   !----- Scaling factor of Tree Area Index, for computing wtveg [dimensionless]. ---------!
   real         :: covr
   !----- Some parameters that were used in ED-2.0, added here for some tests. ------------!
   real         :: ez


   !----- Double precision version of some of these variables (for Runge-Kutta). ----------!
   real(kind=8) :: exar8
   real(kind=8) :: ustmin8
   real(kind=8) :: ugbmin8
   real(kind=8) :: ubmin8
   real(kind=8) :: ez8
   real(kind=8) :: vh2vr8
   real(kind=8) :: vh2dh8
   real(kind=8) :: rasveg_min8
   real(kind=8) :: taumin8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Parameters for Massman (1997) and Massman and Weil (1999) canopy turbulence       !
   ! closures.                                                                             !
   !                                                                                       !
   ! Massman, W. J., 1997: An analytical one-dimensional model of momentum transfer by     !
   !    vegetation of arbitrary structure.  Boundary-Layer Meteorol., 83, 407-421.         !
   !                                                                                       !
   ! Massman, W. J., and J. C. Weil, 1999: An analytical one-dimension second-order clos-  !
   !    ure model turbulence statistics and the Lagrangian time scale within and above     !
   !    plant canopies of arbitrary structure.  Boundary-Layer Meteorol., 91, 81-107.      !
   !                                                                                       !
   ! Wohlfahrt, G., and A. Cernusca, 2002: Momentum transfer by a mountain meadow canopy:  !
   !    a simulation analysis based on Massman's (1997) model.  Boundary-Layer Meteorol.,  !
   !    103, 391-407.
   !---------------------------------------------------------------------------------------!
   !----- Fluid drag coefficient for turbulent flow in leaves at the top. -----------------!
   real(kind=4)  :: cdrag0
   !----- Values from re-fit of the data used by WC02. ------------------------------------!
   real(kind=4)  :: cdrag1
   real(kind=4)  :: cdrag2
   real(kind=4)  :: cdrag3
   !----- Sheltering factor of fluid drag at the top of the canopy. -----------------------!
   real(kind=4)  :: pm0
   !----- Surface drag parameters (Massman 1997). -----------------------------------------!
   real(kind=4)  :: c1_m97
   real(kind=4)  :: c2_m97
   real(kind=4)  :: c3_m97
   !----- Eddy diffusivity due to Von Karman Wakes in gravity flows. ----------------------!
   real(kind=4)  :: kvwake
   !---------------------------------------------------------------------------------------!
   !     Alpha factors to produce the profile of sheltering factor and within canopy drag, !
   ! as suggested by Massman (1997) and Massman and Weil (1999).                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4)  :: alpha_m97
   real(kind=4)  :: alpha_mw99
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for Massman and Weil (1999).                                           !
   !  Gamma and nu are the parameters that close equation 10 in Massman and Weil (1999).   !
   !  VERY IMPORTANT: If you mess with gamma, you must recompute nu!                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4), dimension(3) :: gamma_mw99
   real(kind=4), dimension(3) :: nu_mw99
   !----- Parameter to represent the Roughness sublayer effect. ---------------------------!
   real(kind=4)               :: infunc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters for CLM, at equation 5.103 of CLM-4 techical note.                     !
   !     Oleson, K. W., et al.; Technical description of version 4.0 of the community land !
   !        model (CLM) NCAR Technical Note NCAR/TN-478+STR, Boulder, CO, April 2010.      !
   !---------------------------------------------------------------------------------------!
   real(kind=4)               :: cs_dense0
   real(kind=4)               :: gamma_clm4
   !---------------------------------------------------------------------------------------!



   !----- Double precision version of all variables above. --------------------------------!
   real(kind=8)                            :: dz_m978
   real(kind=8)                            :: cdrag08
   real(kind=8)                            :: cdrag18
   real(kind=8)                            :: cdrag28
   real(kind=8)                            :: cdrag38
   real(kind=8)                            :: pm08
   real(kind=8)                            :: c1_m978
   real(kind=8)                            :: c2_m978
   real(kind=8)                            :: c3_m978
   real(kind=8)                            :: kvwake8
   real(kind=8)                            :: alpha_m97_8
   real(kind=8)                            :: alpha_mw99_8
   real(kind=8), dimension(3)              :: gamma_mw99_8
   real(kind=8), dimension(3)              :: nu_mw99_8
   real(kind=8)                            :: infunc_8
   real(kind=8)                            :: cs_dense08
   real(kind=8)                            :: gamma_clm48
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Parameters for the aerodynamic resistance between the leaf (flat surface) and    !
   ! wood (kind of cylinder surface), and the canopy air space.  These are the A, B, n,    !
   ! and m parameters that define the Nusselt number for forced and free convection, at    !
   ! equations 10.7 and 10.9.  The numbers are defined in appendix A, table A.5.           !
   !                                                                                       !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !---------------------------------------------------------------------------------------!
   real(kind=4)            :: aflat_turb  ! A (forced convection), turbulent flow
   real(kind=4)            :: aflat_lami  ! A (forced convection), laminar   flow
   real(kind=4)            :: nflat_turb  ! n (forced convection), turbulent flow
   real(kind=4)            :: nflat_lami  ! n (forced convection), laminar   flow

   real(kind=4)            :: bflat_turb  ! B (free   convection), turbulent flow
   real(kind=4)            :: bflat_lami  ! B (free   convection), laminar   flow
   real(kind=4)            :: mflat_turb  ! m (free   convection), turbulent flow
   real(kind=4)            :: mflat_lami  ! m (free   convection), laminar   flow

   real(kind=4)            :: ocyli_turb  ! intercept (forced convection), turbulent flow
   real(kind=4)            :: ocyli_lami  ! intercept (forced convection), laminar   flow
   real(kind=4)            :: acyli_turb  ! A (forced convection), turbulent flow
   real(kind=4)            :: acyli_lami  ! A (forced convection), laminar   flow
   real(kind=4)            :: ncyli_turb  ! n (forced convection), turbulent flow
   real(kind=4)            :: ncyli_lami  ! n (forced convection), laminar   flow

   real(kind=4)            :: bcyli_turb  ! B (free   convection), turbulent flow
   real(kind=4)            :: bcyli_lami  ! B (free   convection), laminar   flow
   real(kind=4)            :: mcyli_turb  ! m (free   convection), turbulent flow
   real(kind=4)            :: mcyli_lami  ! m (free   convection), laminar   flow

   real(kind=8)            :: aflat_turb8 ! A (forced convection), turbulent flow
   real(kind=8)            :: aflat_lami8 ! A (forced convection), laminar   flow
   real(kind=8)            :: nflat_turb8 ! n (forced convection), turbulent flow
   real(kind=8)            :: nflat_lami8 ! n (forced convection), laminar   flow

   real(kind=8)            :: bflat_turb8 ! B (free   convection), turbulent flow
   real(kind=8)            :: bflat_lami8 ! B (free   convection), laminar   flow
   real(kind=8)            :: mflat_turb8 ! m (free   convection), turbulent flow
   real(kind=8)            :: mflat_lami8 ! m (free   convection), laminar   flow

   real(kind=8)            :: ocyli_turb8 ! intercept (forced convection), turbulent flow
   real(kind=8)            :: ocyli_lami8 ! intercept (forced convection), laminar   flow
   real(kind=8)            :: acyli_turb8 ! A (forced convection), turbulent flow
   real(kind=8)            :: acyli_lami8 ! A (forced convection), laminar   flow
   real(kind=8)            :: ncyli_turb8 ! n (forced convection), turbulent flow
   real(kind=8)            :: ncyli_lami8 ! n (forced convection), laminar   flow

   real(kind=8)            :: bcyli_turb8 ! B (free   convection), turbulent flow
   real(kind=8)            :: bcyli_lami8 ! B (free   convection), laminar   flow
   real(kind=8)            :: mcyli_turb8 ! m (free   convection), turbulent flow
   real(kind=8)            :: mcyli_lami8 ! m (free   convection), laminar   flow

   real(kind=8)            :: beta_lami8  ! Correction term for Nusselt #, laminar flow
   real(kind=8)            :: beta_turb8  ! Correction term for Nusselt #, turbulent flow
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Parameters for surface layer models.                                             !
   !---------------------------------------------------------------------------------------!
   !----- Louis (1979) model. -------------------------------------------------------------!
   real   :: bl79        ! b prime parameter
   real   :: csm         ! C* for momentum (eqn. 20, not co2 char. scale)
   real   :: csh         ! C* for heat (eqn.20, not co2 char. scale)
   real   :: dl79        ! ???
   !----- Oncley and Dudhia (1995) model. -------------------------------------------------!
   real   :: beta_s      ! Beta for the stable case
   !----- Beljaars and Holtslag (1991) model. ---------------------------------------------!
   real   :: abh91       ! -a from equation  (28) and (32)
   real   :: bbh91       ! -b from equation  (28) and (32)
   real   :: cbh91       !  c from equations (28) and (32)
   real   :: dbh91       !  d from equations (28) and (32)
   real   :: ebh91       ! - factor multiplying a*zeta in equation (32)
   real   :: fbh91       ! exponent in equation (32)
   real   :: cod         ! c/d
   real   :: bcod        ! b*c/d
   real   :: fm1         ! f-1
   real   :: ate         ! a * e
   real   :: atetf       ! a * e * f
   real   :: z0moz0h     ! z0(M)/z0(h)
   real   :: z0hoz0m     ! z0(M)/z0(h)
   !----- Modified CLM (2004) model. ------------------------------------------------------!
   real   :: beta_vs     ! Beta for the very stable case (CLM eq. 5.30)
   real   :: chim        ! CLM coefficient for very unstable case (momentum)
   real   :: chih        ! CLM coefficient for very unstable case (heat)
   real   :: zetac_um    ! critical zeta below which it becomes very unstable (momentum)
   real   :: zetac_uh    ! critical zeta below which it becomes very unstable (heat)
   real   :: zetac_sm    ! critical zeta above which it becomes very stable   (momentum)
   real   :: zetac_sh    ! critical zeta above which it becomes very stable   (heat)
   real   :: zetac_umi   ! 1. / zetac_umi
   real   :: zetac_uhi   ! 1. / zetac_uhi
   real   :: zetac_smi   ! 1. / zetac_smi
   real   :: zetac_shi   ! 1. / zetac_shi
   real   :: zetac_umi16 ! 1/(-zetac_umi)^(1/6)
   real   :: zetac_uhi13 ! 1/(-zetac_umi)^(1/6)
   real   :: psimc_um    ! psim evaluation at zetac_um
   real   :: psihc_uh    ! psih evaluation at zetac_uh
   !---------------------------------------------------------------------------------------!

   !----- Double precision of all these variables. ----------------------------------------!
   real(kind=8)   :: bl798
   real(kind=8)   :: csm8
   real(kind=8)   :: csh8
   real(kind=8)   :: dl798
   real(kind=8)   :: beta_s8
   real(kind=8)   :: gamm8
   real(kind=8)   :: gamh8
   real(kind=8)   :: ribmax8
   real(kind=8)   :: tprandtl8
   real(kind=8)   :: abh918
   real(kind=8)   :: bbh918
   real(kind=8)   :: cbh918
   real(kind=8)   :: dbh918
   real(kind=8)   :: ebh918
   real(kind=8)   :: fbh918
   real(kind=8)   :: cod8
   real(kind=8)   :: bcod8
   real(kind=8)   :: fm18
   real(kind=8)   :: ate8
   real(kind=8)   :: atetf8
   real(kind=8)   :: z0moz0h8
   real(kind=8)   :: z0hoz0m8
   real(kind=8)   :: beta_vs8
   real(kind=8)   :: chim8
   real(kind=8)   :: chih8
   real(kind=8)   :: zetac_um8
   real(kind=8)   :: zetac_uh8
   real(kind=8)   :: zetac_sm8
   real(kind=8)   :: zetac_sh8
   real(kind=8)   :: zetac_umi8
   real(kind=8)   :: zetac_uhi8
   real(kind=8)   :: zetac_smi8
   real(kind=8)   :: zetac_shi8
   real(kind=8)   :: zetac_umi168
   real(kind=8)   :: zetac_uhi138
   real(kind=8)   :: psimc_um8
   real(kind=8)   :: psihc_uh8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   Tuneable parameters that will be set in ed_params.f90.                              !
   !---------------------------------------------------------------------------------------!
   

   !---------------------------------------------------------------------------------------!
   !    Minimum leaf water content to be considered.  Values smaller than this will be     !
   ! flushed to zero.  This value is in kg/[m2 leaf], so it will be always scaled by LAI.  !
   !---------------------------------------------------------------------------------------!
   real :: leaf_drywhc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Parameter to set the minimum vegetation-level aerodynamic conductance.             !
   !---------------------------------------------------------------------------------------!
   real         :: gbhmos_min
   real(kind=8) :: gbhmos_min8
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      This is the minimum vegetation height.  [m]                                      !
   !---------------------------------------------------------------------------------------!
   real         :: veg_height_min
   real(kind=8) :: veg_height_min8
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      This is the minimum canopy depth that is used to calculate the heat and moisture !
   ! storage capacity in the canopy air [m].                                               !
   !---------------------------------------------------------------------------------------!
   real         :: minimum_canopy_depth
   real(kind=8) :: minimum_canopy_depth8
   !---------------------------------------------------------------------------------------!

   !=======================================================================================!
   !=======================================================================================!
   !   Soil conductance terms, from:                                                       !
   !                                                                                       !
   ! Passerat de Silans, A., 1986: Transferts de masse et de chaleur dans un sol stratifié !
   !     soumis à une excitation amtosphérique naturelle. Comparaison: Modèles-expérience. !
   !     Thesis, Institut National Polytechnique de Grenoble. (P86)                        !
   !                                                                                       !
   ! retrieved from:                                                                       !
   ! Mahfouf, J. F., J. Noilhan, 1991: Comparative study of various formulations of        !
   !     evaporation from bare soil using in situ data. J. Appl. Meteorol., 30, 1354-1365. !
   !     (MN91)                                                                            !
   !                                                                                       !
   !     Please notice that the values are inverted because we compute conductance, not    !
   ! resistance.                                                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: ggsoil0
   real(kind=4) :: kksoil
   real(kind=8) :: ggsoil08
   real(kind=8) :: kksoil8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for momentum.            !
   !---------------------------------------------------------------------------------------!
   real function psim(zeta,stable)
      use consts_coms, only : halfpi   & ! intent(in)
                            , onesixth ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: xx
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psim = - beta_s * zeta 
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            psim = abh91 * zeta                                                            &
                 + bbh91 * (zeta - cod) * exp(max(-38.,-dbh91 * zeta))                     &
                 + bcod
         case (4) !----- CLM (2004) (including neglected terms). --------------------------!
            if (zeta > zetac_sm) then
               !----- Very stable case. ---------------------------------------------------!
               psim = (1.0 - beta_vs) * log(zeta * zetac_smi)                              &
                    + (1.0 - beta_s ) * zetac_sm - zeta
            else
               !----- Normal stable case. -------------------------------------------------!
               psim = - beta_s * zeta
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            xx   = sqrt(sqrt(1.0 - gamm * zeta))
            psim = log(0.125 * (1.0+xx) * (1.0+xx) * (1.0 + xx*xx)) - 2.0*atan(xx) + halfpi
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               psim = log(zeta * zetac_umi)                                                &
                    + 6.0 * chim * ((- zeta) ** (-onesixth) - zetac_umi16)                 &
                    + psimc_um
            else
               !----- Normal unstable case. -----------------------------------------------!
               xx   = sqrt(sqrt(1.0 - gamm * zeta))
               psim = log(0.125 * (1.0+xx) * (1.0+xx) * (1.0 + xx*xx))                     &
                    - 2.0*atan(xx) + halfpi
            end if
         end select
      end if
      return
   end function psim
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for heat (and vapour,    !
   ! and carbon dioxide too.)                                                              !
   !---------------------------------------------------------------------------------------!
   real function psih(zeta,stable)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: yy
      !----- External functions. ----------------------------------------------------------!
      real   , external   :: cbrt
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psih = - beta_s * zeta 
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            psih = 1.0 - (1.0 + ate * zeta)**fbh91                                         &
                 + bbh91 * (zeta - cod) * exp(max(-38.,-dbh91 * zeta)) + bcod
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sh) then
               !----- Very stable case. ---------------------------------------------------!
               psih = (1.0 - beta_vs) * log(zeta * zetac_shi)                              &
                    + (1.0 - beta_s ) * zetac_sh - zeta
            else
               !----- Normal stable case. -------------------------------------------------!
               psih = - beta_s * zeta 
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            yy   = sqrt(1.0 - gamh * zeta)
            psih = log(0.25 * (1.0+yy) * (1.0+yy))
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               psih = log(zeta * zetac_uhi)                                                &
                    + 3.0 * chih * (1./cbrt(-zeta) - zetac_uhi13)                          &
                    + psihc_uh
            else
               !----- Normal unstable case. -----------------------------------------------!
               yy   = sqrt(1.0 - gamh * zeta)
               psih = log(0.25 * (1.0+yy) * (1.0+yy))
            end if
         end select
      end if
      return
   end function psih
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for momentum.            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function psim8(zeta,stable)
      use consts_coms, only : halfpi8   & ! intent(in)
                            , onesixth8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: zeta   ! z/L, z = height, and L = Obukhov length   [ ---]
      logical     , intent(in) :: stable ! Flag... This surface layer is stable      [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: xx
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psim8 = - beta_s8 * zeta 
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            psim8 = abh918 * zeta                                                          &
                  + bbh918 * (zeta - cod8) * exp(max(-3.8d1,-dbh918 * zeta))               &
                  + bcod8
         case (4) !----- CLM (2004) (including neglected terms). --------------------------!
            if (zeta > zetac_sm8) then
               !----- Very stable case. ---------------------------------------------------!
               psim8 = (1.d0 - beta_vs8) * log(zeta * zetac_smi8)                          &
                     + (1.d0 - beta_s8 ) * zetac_sm8 - zeta
            else
               !----- Normal stable case. -------------------------------------------------!
               psim8 = - beta_s8 * zeta
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            xx   = sqrt(sqrt(1.d0 - gamm8 * zeta))
            psim8 = log(1.25d-1 * (1.d0+xx) * (1.d0+xx) * (1.d0 + xx*xx))                  &
                  - 2.d0*atan(xx) + halfpi8
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um8) then
               !----- Very unstable case. -------------------------------------------------!
               psim8 = log(zeta * zetac_umi8)                                              &
                     + 6.d0 * chim8 * ((-zeta) ** (-onesixth8) - zetac_umi168)             &
                     + psimc_um
            else
               !----- Normal unstable case. -----------------------------------------------!
               xx   = sqrt(sqrt(1.d0 - gamm8 * zeta))
               psim8 = log(1.25d-1 * (1.d0+xx) * (1.d0+xx) * (1.d0 + xx*xx))               &
                     - 2.d0*atan(xx) + halfpi8
            end if
         end select
      end if
      return
   end function psim8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the stability  correction function for heat (and vapour,    !
   ! and carbon dioxide too.)                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function psih8(zeta,stable)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: zeta   ! z/L, z = height, and L = Obukhov length   [ ---]
      logical     , intent(in) :: stable ! Flag... This surface layer is stable      [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: yy
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external   :: cbrt8
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            psih8 = - beta_s8 * zeta 
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            psih8 = 1.d0 - (1.d0 + ate8 * zeta)**fbh918                                    &
                  + bbh918 * (zeta - cod8) * exp(max(-3.8d1,-dbh918 * zeta)) + bcod8
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sh8) then
               !----- Very stable case. ---------------------------------------------------!
               psih8 = (1.d0 - beta_vs8) * log(zeta * zetac_shi8)                          &
                     + (1.d0 - beta_s8 ) * zetac_sh8 - zeta
            else
               !----- Normal stable case. -------------------------------------------------!
               psih8 = - beta_s8 * zeta 
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            yy    = sqrt(1.d0 - gamh8 * zeta)
            psih8 = log(2.5d-1 * (1.d0+yy) * (1.d0+yy))
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               psih8 = log(zeta * zetac_uhi8)                                              &
                     + 3.d0 * chih8 * (1.d0/cbrt8(-zeta) - zetac_uhi138)                   &
                     + psihc_uh8
            else
               !----- Normal unstable case. -----------------------------------------------!
               yy    = sqrt(1.d0 - gamh8 * zeta)
               psih8 = log(2.5d-1 * (1.d0+yy) * (1.d0+yy))
            end if
         end select
      end if
      return
   end function psih8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! momentum with respect to zeta.                                                        !
   !---------------------------------------------------------------------------------------!
   real function dpsimdzeta(zeta,stable)
      use consts_coms, only : onesixth ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: xx
      !----- External functions. ----------------------------------------------------------!
      real   , external   :: cbrt
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsimdzeta = - beta_s 
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            dpsimdzeta = abh91 + bbh91 * (1.0 - dbh91 * zeta + cbh91)                      &
                               * exp(max(-38.,-dbh91 * zeta))
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sm) then
               !----- Very stable case. ---------------------------------------------------!
               dpsimdzeta = (1.0 - beta_vs) / zeta - 1.0
            else
               !----- Normal stable case. -------------------------------------------------!
               dpsimdzeta = - beta_s 
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            xx         = sqrt(sqrt(1.0 - gamm * zeta))
            dpsimdzeta = - gamm / (xx * (1.0+xx) * (1.0 + xx*xx)) 
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               dpsimdzeta = (1.0 - chim * (-zeta)**onesixth) / zeta
            else
               !----- Normal unstable case. -----------------------------------------------!
               xx         = sqrt(sqrt(1.0 - gamm * zeta))
               dpsimdzeta = - gamm / (xx * (1.0+xx) * (1.0 + xx*xx))
            end if
         end select
      end if
      return
   end function dpsimdzeta
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! heat/moisture/CO2 with respect to zeta.                                               !
   !---------------------------------------------------------------------------------------!
   real function dpsihdzeta(zeta,stable)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: zeta   ! z/L, z is the height, and L the Obukhov length [ ---]
      logical, intent(in) :: stable ! Flag... This surface layer is stable           [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: yy
      !----- External functions. ----------------------------------------------------------!
      real   , external   :: cbrt
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsihdzeta = - beta_s
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            dpsihdzeta = - atetf * (1.0 + ate * zeta)**fm1                                 &
                         + bbh91 * (1.0 - dbh91 * zeta + cbh91)                            &
                         * exp(max(-38.,-dbh91 * zeta))
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sh) then
               !----- Very stable case. ---------------------------------------------------!
               dpsihdzeta = (1.0 - beta_vs) / zeta - 1.0
            else
               !----- Normal stable case. -------------------------------------------------!
               dpsihdzeta = - beta_s
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            yy   = sqrt(1.0 - gamh * zeta)
            dpsihdzeta = -gamh / (yy * (1.0 + yy))
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um) then
               !----- Very unstable case. -------------------------------------------------!
               dpsihdzeta = (1.0 + chih / cbrt(zeta)) / zeta
            else
               !----- Normal unstable case. -----------------------------------------------!
               yy   = sqrt(1.0 - gamh * zeta)
               dpsihdzeta = -gamh / (yy * (1.0 + yy))
            end if
         end select
      end if
      return
   end function dpsihdzeta
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! momentum with respect to zeta.                                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dpsimdzeta8(zeta,stable)
      use consts_coms, only : halfpi8   & ! intent(in)
                            , onesixth8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: zeta   ! z/L, z = height, and L = Obukhov length   [ ---]
      logical     , intent(in) :: stable ! Flag... This surface layer is stable      [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: xx
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external   :: cbrt8
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsimdzeta8 = - beta_s8
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            dpsimdzeta8 = abh918                                                           &
                        + bbh918 * (1.d0 - dbh918 * zeta + cbh918)                         &
                        * exp(max(-3.8d1,-dbh918 * zeta))
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sm8) then
               !----- Very stable case. ---------------------------------------------------!
               dpsimdzeta8 = (1.d0 - beta_vs8) / zeta - 1.d0
            else
               !----- Normal stable case. -------------------------------------------------!
               dpsimdzeta8 = - beta_s8
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            xx          = sqrt(sqrt(1.d0 - gamm8 * zeta))
            dpsimdzeta8 = - gamm8 / (xx * (1.d0+xx) * (1.d0 + xx*xx)) 
         case (4)   !----- Modified CLM (2004). -------------------------------------------!
            if (zeta < zetac_um8) then
               !----- Very unstable case. -------------------------------------------------!
               dpsimdzeta8 = (1.d0 - chim8 * (-zeta)**onesixth8) / zeta
            else
               !----- Normal unstable case. -----------------------------------------------!
               xx          = sqrt(sqrt(1.d0 - gamm8 * zeta))
               dpsimdzeta8 = - gamm8 / (xx * (1.d0+xx) * (1.d0 + xx*xx)) 
            end if
         end select
      end if
      return
   end function dpsimdzeta8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the derivative of the stability correction function for     !
   ! heat/moisture/CO2 with respect to zeta.                                               !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function dpsihdzeta8(zeta,stable)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: zeta   ! z/L, z = height, and L = Obukhov length   [ ---]
      logical     , intent(in) :: stable ! Flag... This surface layer is stable      [ T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: yy
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external   :: cbrt8
      !------------------------------------------------------------------------------------!
      if (stable) then
         select case (isfclyrm)
         case (2) !----- Oncley and Dudhia (1995). ----------------------------------------!
            dpsihdzeta8 = - beta_s8
         case (0,3) !----- Beljaars and Holtslag (1991). ----------------------------------!
            dpsihdzeta8 = - atetf8 * (1.d0 + ate8 * zeta)**fm18                            &
                          + bbh918 * (1.d0 - dbh918 * zeta + cbh918)                       &
                          * exp(max(-3.8d1,-dbh918 * zeta))
         case (4) !----- CLM (2004). ------------------------------------------------------!
            if (zeta > zetac_sh8) then
               !----- Very stable case. ---------------------------------------------------!
               dpsihdzeta8 = (1.d0 - beta_vs8) / zeta - 1.d0
            else
               !----- Normal stable case. -------------------------------------------------!
               dpsihdzeta8 = - beta_s8
            end if
         end select
      else
         select case (isfclyrm)
         case (0,2,3) !----- Oncley and Dudhia (1995) and Beljaars and Holtslag (1991). ---!
            yy          = sqrt(1.d0 - gamh8 * zeta)
            dpsihdzeta8 = -gamh8 / (yy * (1.d0 + yy))
         case (4)   !----- CLM (2004) (including neglected terms). ------------------------!
            if (zeta < zetac_um8) then
               !----- Very unstable case. -------------------------------------------------!
               dpsihdzeta8 = (1.d0 + chih8 / cbrt8(zeta)) / zeta
            else
               !----- Normal unstable case. -----------------------------------------------!
               yy   = sqrt(1.d0 - gamh8 * zeta)
               dpsihdzeta8 = -gamh8 / (yy * (1.d0 + yy))
            end if
         end select
      end if
      return
   end function dpsihdzeta8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the value of zeta for a given Richardson number, reference    !
   ! height and the roughness scale.  This is solved by using the definition of Obukhov    !
   ! length scale as stated in Louis (1979) equation (10), modified to define z/L rather   !
   ! than L.  The solution is found  iteratively since it's not a simple function to       !
   ! invert.  It tries to use Newton's method, which should take care of most cases.  In   !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real function zoobukhov(rib,zstar,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
      use therm_lib, only : toler  & ! intent(in)
                          , maxfpo & ! intent(in)
                          , maxit  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real   , intent(in) :: rib       ! Bulk Richardson number                    [   ---]
      real   , intent(in) :: zstar     ! Reference height - displacement height    [     m]
      real   , intent(in) :: rough     ! Roughness length scale                    [     m]
      real   , intent(in) :: zoz0m     ! zstar/roughness(momentum)                 [   ---]
      real   , intent(in) :: lnzoz0m   ! ln[zstar/roughness(momentum)]             [   ---]
      real   , intent(in) :: zoz0h     ! zstar/roughness(heat)                     [   ---]
      real   , intent(in) :: lnzoz0h   ! ln[zstar/roughness(heat)]                 [   ---]
      logical, intent(in) :: stable    ! Flag... This surface layer is stable      [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real                :: ribuse    ! Richardson number to use                  [   ---]
      real                :: fm        ! lnzoz0 - psim(zeta) + psim(zeta0)         [   ---]
      real                :: fh        ! lnzoz0 - psih(zeta) + psih(zeta0)         [   ---]
      real                :: dfmdzeta  ! d(fm)/d(zeta)                             [   ---]
      real                :: dfhdzeta  ! d(fh)/d(zeta)                             [   ---]
      real                :: z0moz     ! Roughness(momentum) / Reference height    [   ---]
      real                :: zeta0m    ! Roughness(momentum) / Obukhov length      [   ---]
      real                :: z0hoz     ! Roughness(heat) / Reference height        [   ---]
      real                :: zeta0h    ! Roughness(heat) / Obukhov length          [   ---]
      real                :: zetaa     ! Smallest guess (or previous guess)        [   ---]
      real                :: zetaz     ! Largest guess (or new guess in Newton's)  [   ---]
      real                :: deriv     ! Function Derivative                       [   ---]
      real                :: fun       ! Function for which we seek a root.        [   ---]
      real                :: funa      ! Smallest guess function.                  [   ---]
      real                :: funz      ! Largest guess function.                   [   ---]
      real                :: delta0    ! Aux. var --- 2nd guess for bisection      [   ---]
      real                :: delta     ! Aux. var --- 2nd guess for bisection      [   ---]
      real                :: coeff     ! RiB * zstar / (Pr * (zstar - z0))         [   ---]
      real                :: zetamin   ! Minimum zeta for stable case.             [   ---]
      real                :: zetamax   ! Maximum zeta for unstable case.           [   ---]
      real                :: zetasmall ! Number sufficiently close to zero         [   ---]
      integer             :: itb       ! Iteration counters                        [   ---]
      integer             :: itn       ! Iteration counters                        [   ---]
      integer             :: itp       ! Iteration counters                        [   ---]
      logical             :: converged ! Flag... The method converged!             [   T|F]
      logical             :: zside     ! Flag... I'm on the z-side.                [   T|F]
      !------------------------------------------------------------------------------------!



      !----- Define some values that won't change during the iterative method. ------------!
      z0moz = 1. / zoz0m
      z0hoz = 1. / zoz0h
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     First thing, check whether this is a stable case and we are running methods 2  !
      ! or 4.  In these methods, there is a singularity that must be avoided.              !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (2,4)
         ribuse = min(rib, (1.0 - toler) * tprandtl / beta_s)

         !---------------------------------------------------------------------------------!
         !    Stable case, use Oncley and Dudhia, we can solve it analytically.            !
         !---------------------------------------------------------------------------------!
         if (stable .and. isfclyrm == 2) then
            zoobukhov = ribuse * zstar * min(lnzoz0m,lnzoz0h)                              &
                      / ( (zstar-rough) * (tprandtl - beta_s * ribuse) )
            return
         end if
         !---------------------------------------------------------------------------------!
      case default
         ribuse = rib
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Define the coefficient Ri * zstar / [Pr * (zstar-z0)]                          !
      !------------------------------------------------------------------------------------!
      coeff = ribuse * zstar / (tprandtl * (zstar - rough))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the bulk Richardson number is zero or almost zero, then we rather just      !
      ! assign z/L to be the one similar to Oncley and Dudhia (1995).  This saves time and !
      ! also avoids the risk of having zeta with the opposite sign.                        !
      !------------------------------------------------------------------------------------!
      zetasmall = coeff * min(lnzoz0m,lnzoz0h)
      if (ribuse <= 0. .and. zetasmall > - z0moz0h * toler) then
         zoobukhov = zetasmall
         return
      elseif (ribuse > 0. .and. zetasmall < z0moz0h * toler) then
         zoobukhov = zetasmall / (1.0 - beta_s * ribuse / tprandtl)
         return
      else
         zetamin    =  toler
         zetamax    = -toler
      end if
      !------------------------------------------------------------------------------------!



      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(5(a,1x,f11.4,1x),a,l1)')                                        &
      !   'Input values: Rib =',rib,'zstar=',star,'rough=',rough,'zoz0=',zoz0              &
      !           ,'lnzoz0=',lnzoz0,'stable=',stable
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !------------------------------------------------------------------------------------!
      !     First guess, use Oncley and Dudhia (1995) approximation for unstable case.     !
      ! We won't use the stable case to avoid FPE or zeta with opposite sign when          !
      ! Ri is too positive.                                                                !
      !------------------------------------------------------------------------------------!
      zetaa = zetasmall
      !------------------------------------------------------------------------------------!



      !----- Find the function and its derivative. ----------------------------------------!
      zeta0m   = zetaa * z0moz
      zeta0h   = zetaa * z0hoz
      fm       = lnzoz0m - psim(zetaa,stable) + psim(zeta0m,stable)
      fh       = lnzoz0h - psih(zetaa,stable) + psih(zeta0h,stable)
      dfmdzeta = z0moz * dpsimdzeta(zeta0m,stable) - dpsimdzeta(zetaa,stable)
      dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaa,stable)
      funa     = coeff * fm * fm / fh - zetaa
      deriv    = coeff * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / (fh * fh) - 1.
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration. ---------------------------!
      zetaz = zetaa
      fun   = funa
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
      !   '1STGSS: itn=',0,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm          &
      !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !----- Enter Newton's method loop. --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1, maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Newton's method converges fast when it's on the right track, but there are  !
         ! cases in which it becomes ill-behaved.  Two situations are known to cause       !
         ! trouble:                                                                        !
         ! 1.  If the derivative is tiny, the next guess can be too far from the actual    !
         !     answer;                                                                     !
         ! 2.  For this specific problem, when zeta is too close to zero.  In this case    !
         !     the derivative will tend to infinity at this point and Newton's method is   !
         !     not going to perform well and can potentially enter in a weird behaviour or !
         !     lead to the wrong answer.  In any case, so we rather go with bisection.     !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler) then
            exit newloop
         elseif(stable .and. (zetaz - fun/deriv < zetamin)) then
            exit newloop
         elseif((.not. stable) .and. (zetaz - fun/deriv > zetamax)) then
            exit newloop
         end if

         !----- Copy the previous guess ---------------------------------------------------!
         zetaa = zetaz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         zetaz    = zetaa - fun/deriv
         zeta0m   = zetaz * z0moz
         zeta0h   = zetaz * z0hoz
         fm       = lnzoz0m - psim(zetaz,stable) + psim(zeta0m,stable)
         fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
         dfmdzeta = z0moz * dpsimdzeta(zeta0m,stable) - dpsimdzeta(zetaz,stable)
         dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaz,stable)
         fun      = coeff * fm * fm / fh - zetaz
         deriv    = coeff * (2. * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / (fh * fh) - 1.
         !---------------------------------------------------------------------------------!

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'NEWTON: itn=',itn,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm        &
         !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         converged = abs(zetaz-zetaa) < toler * abs(zetaz)

         if (converged) then
            zoobukhov = 0.5 * (zetaa+zetaz)
            return
         elseif (fun == 0.0) then !---- Converged by luck. --------------------------------!
            zoobukhov = zetaz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If we reached this point then it's because Newton's method failed or it has    !
      ! become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a  !
      ! fancier bisection.  For this we need two guesses, and the guesses must have        !
      ! opposite signs.                                                                    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.0) then
         funz  = fun
         zside = .true. 
      else
         if (abs(fun-funa) < 100. * toler * abs(zetaa)) then
            if (stable) then
               delta = max(0.5 * abs(zetaa-zetamin),100. * toler * abs(zetaa))
            else
               delta = max(0.5 * abs(zetaa-zetamax),100. * toler * abs(zetaa))
            end if
         else
            if (stable) then
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,100. * toler * abs(zetaa)                                       &
                          ,0.5 * abs(zetaa-zetamin))
            else
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,100. * toler * abs(zetaa)                                       &
                          ,0.5 * abs(zetaa-zetamax))
            end if
         end if
         if (stable) then
            zetaz = max(zetamin,zetaa + delta)
         else
            zetaz = min(zetamax,zetaa + delta)
         end if
         zside = .false.
         zgssloop: do itp=1,maxfpo
            if (stable) then
               zetaz    = max(zetamin,zetaa + real((-1)**itp * (itp+3)/2) * delta)
            else
               zetaz    = min(zetamax,zetaa + real((-1)**itp * (itp+3)/2) * delta)
            end if
            zeta0m   = zetaz * z0moz
            zeta0h   = zetaz * z0hoz
            fm       = lnzoz0m - psim(zetaz,stable) + psim(zeta0m,stable)
            fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
            funz     = coeff * fm * fm / fh - zetaz
            zside    = funa * funz < 0.0
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                 &
            !   '2NDGSS: itp=',itp,'zside=',zside,'zetaa=',zetaa,'zetaz=',zetaz             &
            !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(a)') '    No second guess for you...'
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zstar  =',zstar  ,'rough  =',rough
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'lnzoz0m=',lnzoz0m,'lnzoz0h=',lnzoz0h
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'rib    =',rib    ,'ribuse =',ribuse
            write (unit=*,fmt='(1(a,1x,l1,1x))')     'stable =',stable
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'fun    =',fun    ,'delta  =',delta
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaa  =',zetaa  ,'funa   =',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaz  =',zetaz  ,'funz   =',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                            ,'zoobukhov','canopy_air_coms.f90')
         end if
      end if

      !----- Now we are ready to start the regula falsi method. ---------------------------!
      bisloop: do itb=itn,maxfpo
         zoobukhov = (funz*zetaa-funa*zetaz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(zoobukhov-zetaa) < toler * abs(zoobukhov)
         if (converged) exit bisloop

         !------ Update function evaluation. ----------------------------------------------!
         zeta0m   = zoobukhov * z0moz
         zeta0h   = zoobukhov * z0hoz
         fm       = lnzoz0m - psim(zoobukhov,stable) + psim(zeta0m,stable)
         fh       = lnzoz0h - psih(zoobukhov,stable) + psih(zeta0h,stable)
         fun      = coeff * fm * fm / fh - zoobukhov
         !---------------------------------------------------------------------------------!

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'REGULA: itn=',itb,'bisection=',.true.,'zetaa=',zetaa,'zetaz=',zetaz,'fun=',fun   &
         !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

         !------ Define new interval based on the intermediate value theorem. -------------!
         if (fun*funa < 0. ) then
            zetaz = zoobukhov
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 0.5
            !----- We just updated zside, set zside to true. ------------------------------!
            zside = .true.
         else
            zetaa = zoobukhov
            funa  = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 0.5
            !----- We just updated aside, set aside to true. ------------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged .or.                                                              &
          (stable .and. zoobukhov < 0.0) .or. (.not. stable .and. zoobukhov > 0.0)) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Zeta finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rib             [   ---] =',rib
         write (unit=*,fmt='(a,1x,f12.4)' ) 'ribuse          [   ---] =',ribuse
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zstar           [     m] =',zstar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rough           [     m] =',rough
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0m           [   ---] =',zoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0m         [   ---] =',lnzoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0h           [   ---] =',zoz0h
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0h         [   ---] =',lnzoz0h
         write (unit=*,fmt='(a,1x,l1)'    ) 'stable          [   T|F] =',stable
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaa           [   ---] =',zetaa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaz           [   ---] =',zetaz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [   ---] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fm              [   ---] =',fm
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fh              [   ---] =',fh
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [   ---] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [   ---] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [   ---] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [   ---] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [   ---] ='                   &
                                                            ,abs(zetaz-zetaa)/abs(zetaz)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoobukhov       [   ---] =',zoobukhov
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

         call fatal_error('Zeta didn''t converge, giving up!!!'                            &
                         ,'zoobukhov','canopy_air_coms.f90')
      end if

      return
   end function zoobukhov
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the value of zeta for a given Richardson number, reference    !
   ! height and the roughness scale.  This is solved by using the definition of Obukhov    !
   ! length scale as stated in Louis (1979) equation (10), modified to define z/L rather   !
   ! than L.  The solution is found  iteratively since it's not a simple function to       !
   ! invert.  It tries to use Newton's method, which should take care of most cases.  In   !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function zoobukhov8(rib,zstar,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
      use therm_lib8, only : toler8 & ! intent(in)
                           , maxfpo & ! intent(in)
                           , maxit  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: rib       ! Bulk Richardson number               [   ---]
      real(kind=8), intent(in) :: zstar     ! Reference height - displacement hgt. [     m]
      real(kind=8), intent(in) :: rough     ! Roughness length scale               [     m]
      real(kind=8), intent(in) :: zoz0m     ! zstar/roughness(momentum)            [   ---]
      real(kind=8), intent(in) :: lnzoz0m   ! ln[zstar/roughness(momentum)]        [   ---]
      real(kind=8), intent(in) :: zoz0h     ! zstar/roughness(heat)                [   ---]
      real(kind=8), intent(in) :: lnzoz0h   ! ln[zstar/roughness(heat)]            [   ---]
      logical     , intent(in) :: stable    ! Flag... This surface layer is stable [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: ribuse    ! Richardson number to use             [   ---]
      real(kind=8)             :: fm        ! lnzoz0 - psim(zeta) + psim(zeta0)    [   ---]
      real(kind=8)             :: fh        ! lnzoz0 - psih(zeta) + psih(zeta0)    [   ---]
      real(kind=8)             :: dfmdzeta  ! d(fm)/d(zeta)                        [   ---]
      real(kind=8)             :: dfhdzeta  ! d(fh)/d(zeta)                        [   ---]
      real(kind=8)             :: z0moz     ! Roughness(momentum) / Reference hgt. [   ---]
      real(kind=8)             :: zeta0m    ! Roughness(momentum) / Obukhov length [   ---]
      real(kind=8)             :: z0hoz     ! Roughness(heat) / Reference height   [   ---]
      real(kind=8)             :: zeta0h    ! Roughness(heat) / Obukhov length     [   ---]
      real(kind=8)             :: zetaa     ! Smallest guess (or previous guess)   [   ---]
      real(kind=8)             :: zetaz     ! Largest guess (new guess in Newton)  [   ---]
      real(kind=8)             :: deriv     ! Function Derivative                  [   ---]
      real(kind=8)             :: fun       ! Function for which we seek a root.   [   ---]
      real(kind=8)             :: funa      ! Smallest guess function.             [   ---]
      real(kind=8)             :: funz      ! Largest guess function.              [   ---]
      real(kind=8)             :: delta0    ! Aux. var --- 2nd guess for bisection [   ---]
      real(kind=8)             :: delta     ! Aux. var --- 2nd guess for bisection [   ---]
      real(kind=8)             :: coeff     ! RiB * zstar / (Pr * (zstar - z0))    [   ---]
      real(kind=8)             :: zetamin   ! Minimum zeta for stable case.        [   ---]
      real(kind=8)             :: zetamax   ! Maximum zeta for unstable case.      [   ---]
      real(kind=8)             :: zetasmall ! Zeta dangerously close to zero       [   ---]
      integer                  :: itn       ! Iteration counters                   [   ---]
      integer                  :: itb       ! Iteration counters                   [   ---]
      integer                  :: itp       ! Iteration counters                   [   ---]
      logical                  :: converged ! Flag... The method converged!        [   T|F]
      logical                  :: zside     ! Flag... I'm on the z-side.           [   T|F]
      !------------------------------------------------------------------------------------!



      !----- Define some values that won't change during the iterative method. ------------!
      z0moz = 1.d0 / zoz0m
      z0hoz = 1.d0 / zoz0h
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     First thing, check whether this is a stable case and we are running methods 2  !
      ! or 4.  In these methods, there is a singularity that must be avoided.              !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (2,4)
         ribuse = min(rib, (1.d0 - toler8) * tprandtl8 / beta_s8 )

         !---------------------------------------------------------------------------------!
         !    Stable case, use Oncley and Dudhia, we can solve it analytically.            !
         !---------------------------------------------------------------------------------!
         if (stable .and. isfclyrm == 2) then
            zoobukhov8 = ribuse * zstar * min(lnzoz0m,lnzoz0h)                             &
                      / ( (zstar-rough) * (tprandtl8 - beta_s8 * ribuse) )
            return
         end if
         !---------------------------------------------------------------------------------!
      case default
         ribuse = rib
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Define the coefficient Ri * zstar / [Pr * (zstar-z0)]                          !
      !------------------------------------------------------------------------------------!
      coeff = ribuse * zstar / (tprandtl8 * (zstar - rough))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the bulk Richardson number is zero or almost zero, then we rather just      !
      ! assign z/L to be the one similar to Oncley and Dudhia (1995).  This saves time and !
      ! also avoids the risk of having zeta with the opposite sign.                        !
      !------------------------------------------------------------------------------------!
      zetasmall = coeff * min(lnzoz0m,lnzoz0h)
      if (ribuse <= 0.d0 .and. zetasmall > - z0moz0h * toler8) then
         zoobukhov8 = zetasmall
         return
      elseif (ribuse > 0.d0 .and. zetasmall < z0moz0h * toler8) then
         zoobukhov8 = zetasmall / (1.d0 - beta_s8 * ribuse / tprandtl8)
         return
      else
         zetamin    =  toler8
         zetamax    = -toler8
      end if
      !------------------------------------------------------------------------------------!



      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(5(a,1x,f11.4,1x),a,l1)')                                        &
      !   'Input values: Rib =',rib,'zstar=',star,'rough=',rough,'zoz0=',zoz0              &
      !           ,'lnzoz0=',lnzoz0,'stable=',stable
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!



      !------------------------------------------------------------------------------------!
      !     First guess, use Oncley and Dudhia (1995) approximation for unstable case.     !
      ! We won't use the stable case to avoid FPE or zeta with opposite sign when          !
      ! Ri is too positive.                                                                !
      !------------------------------------------------------------------------------------!
      zetaa = zetasmall
      !------------------------------------------------------------------------------------!



      !----- Find the function and its derivative. ----------------------------------------!
      zeta0m   = zetaa * z0moz
      zeta0h   = zetaa * z0hoz
      fm       = lnzoz0m - psim8(zetaa,stable) + psim8(zeta0m,stable)
      fh       = lnzoz0h - psih8(zetaa,stable) + psih8(zeta0h,stable)
      dfmdzeta = z0moz * dpsimdzeta8(zeta0m,stable) - dpsimdzeta8(zetaa,stable)
      dfhdzeta = z0hoz * dpsihdzeta8(zeta0h,stable) - dpsihdzeta8(zetaa,stable)
      funa     = coeff * fm * fm / fh - zetaa
      deriv    = coeff * (2.d0 * fm * dfmdzeta * fh - fm * fm * dfhdzeta) / (fh * fh) - 1.d0
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration. ---------------------------!
      zetaz = zetaa
      fun   = funa
      !------------------------------------------------------------------------------------!



      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
      !   '1STGSS: itn=',0,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm          &
      !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Enter Newton's method loop. --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1, maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Newton's method converges fast when it's on the right track, but there are  !
         ! cases in which it becomes ill-behaved.  Two situations are known to cause       !
         ! trouble:                                                                        !
         ! 1.  If the derivative is tiny, the next guess can be too far from the actual    !
         !     answer;                                                                     !
         ! 2.  For this specific problem, when zeta is too close to zero.  In this case    !
         !     the derivative will tend to infinity at this point and Newton's method is   !
         !     not going to perform well and can potentially enter in a weird behaviour or !
         !     lead to the wrong answer.  In any case, so we rather go with bisection.     !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler8) then
            exit newloop
         elseif(stable .and. (zetaz - fun/deriv < zetamin)) then
            exit newloop
         elseif((.not. stable) .and. (zetaz - fun/deriv > zetamax)) then
            exit newloop
         end if

         !----- Copying the previous guess ------------------------------------------------!
         zetaa = zetaz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         zetaz = zetaa - fun/deriv

         zeta0m   = zetaz * z0moz
         zeta0h   = zetaz * z0hoz
         fm       = lnzoz0m - psim8(zetaz,stable) + psim8(zeta0m,stable)
         fh       = lnzoz0h - psih8(zetaz,stable) + psih8(zeta0h,stable)
         dfmdzeta = z0moz * dpsimdzeta8(zeta0m,stable) - dpsimdzeta8(zetaz,stable)
         dfhdzeta = z0hoz * dpsihdzeta8(zeta0h,stable) - dpsihdzeta8(zetaz,stable)
         fun      = coeff * fm * fm / fh - zetaz
         deriv    = coeff * (2.d0 * fm * dfmdzeta * fh - fm * fm * dfhdzeta)               &
                  / (fh * fh) - 1.d0

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'NEWTON: itn=',itn,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm        &
         !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         converged = abs(zetaz-zetaa) < toler8 * abs(zetaz)

         if (converged) then
            zoobukhov8 = 5.d-1 * (zetaa+zetaz)
            return
         elseif (fun == 0.d0) then !---- Converged by luck. -------------------------------!
            zoobukhov8 = zetaz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If we reached this point then it's because Newton's method failed or it has    !
      ! become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a  !
      ! fancier bisection.  For this we need two guesses, and the guesses must have        !
      ! opposite signs.                                                                    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         funz  = fun
         zside = .true. 
      else
         if (abs(fun-funa) < 1.d2 * toler8 * abs(zetaa)) then
            if (stable) then
               delta = max(5.d-1 * abs(zetaa-zetamin),1.d2 * toler8 * abs(zetaa))
            else
               delta = max(5.d-1 * abs(zetaa-zetamax),1.d2 * toler8 * abs(zetaa))
            end if
         else
            if (stable) then
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,1.d2 * toler8 * abs(zetaa)                                      &
                          ,5.d-1 * abs(zetaa-zetamin))
            else
               delta = max(abs(funa * (zetaz-zetaa)/(fun-funa))                            &
                          ,1.d2 * toler8 * abs(zetaa)                                      &
                          ,5.d-1 * abs(zetaa-zetamax))
            end if
         end if
         if (stable) then
            zetaz = max(zetamin,zetaa + delta)
         else
            zetaz = min(zetamax,zetaa + delta)
         end if
         zside = .false.
         zgssloop: do itp=1,maxfpo
            if (stable) then
               zetaz    = max(zetamin,zetaa + dble((-1)**itp * (itp+3)/2) * delta)
            else
               zetaz    = min(zetamax,zetaa + dble((-1)**itp * (itp+3)/2) * delta)
            end if
            zeta0m   = zetaz * z0moz
            zeta0h   = zetaz * z0hoz
            fm       = lnzoz0m - psim8(zetaz,stable) + psim8(zeta0m,stable)
            fh       = lnzoz0h - psih8(zetaz,stable) + psih8(zeta0h,stable)
            funz     = coeff * fm * fm / fh - zetaz
            zside    = funa * funz < 0.d0
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                 &
            !   '2NDGSS: itp=',itp,'zside=',zside,'zetaa=',zetaa,'zetaz=',zetaz             &
            !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(a)') '    No second guess for you...'
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zstar  =',zstar  ,'rough  =',rough
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'lnzoz0m=',lnzoz0m,'lnzoz0h=',lnzoz0h
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'rib    =',rib    ,'ribuse=',ribuse
            write (unit=*,fmt='(1(a,1x,l1,1x))')     'stable =',stable
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'fun    =',fun   ,'delta  =',delta
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaa  =',zetaa ,'funa   =',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaz  =',zetaz ,'funz   =',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                            ,'zoobukhov8','canopy_air_coms.f90')
         end if
      end if

      !----- Now we are ready to start the regula falsi method. ---------------------------!
      bisloop: do itb=itn,maxfpo
         zoobukhov8 = (funz*zetaa-funa*zetaz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(zoobukhov8-zetaa) < toler8 * abs(zoobukhov8)
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         zeta0m   = zoobukhov8 * z0moz
         zeta0h   = zoobukhov8 * z0hoz
         fm       = lnzoz0m - psim8(zoobukhov8,stable) + psim8(zeta0m,stable)
         fh       = lnzoz0h - psih8(zoobukhov8,stable) + psih8(zeta0h,stable)
         fun      = coeff * fm * fm / fh - zoobukhov8

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'REGULA: itn=',itb,'bisection=',.true.,'zetaa=',zetaa,'zetaz=',zetaz,'fun=',fun   &
         !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0.d0 ) then
            zetaz = zoobukhov8
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 5.d-1
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            zetaa = zoobukhov8
            funa  = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 5.d-1
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged .or.                                                              &
          (stable .and. zoobukhov8 < 0.d0) .or. (.not. stable .and. zoobukhov8 > 0.d0)) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Zeta finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rib             [   ---] =',rib
         write (unit=*,fmt='(a,1x,f12.4)' ) 'ribuse          [   ---] =',ribuse
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zstar           [     m] =',zstar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rough           [     m] =',rough
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0m           [   ---] =',zoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0m         [   ---] =',lnzoz0m
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0h           [   ---] =',zoz0h
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0h         [   ---] =',lnzoz0h
         write (unit=*,fmt='(a,1x,l1)'    ) 'stable          [   T|F] =',stable
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (zoobukhov8 values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaa           [   ---] =',zetaa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaz           [   ---] =',zetaz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [   ---] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fm              [   ---] =',fm
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fh              [   ---] =',fh
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [   ---] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [   ---] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [   ---] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [   ---] =',toler8
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [   ---] ='                   &
                                                            ,abs(zetaz-zetaa)/abs(zetaz)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoobukhov8      [   ---] =',zoobukhov8
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

         call fatal_error('Zeta didn''t converge, giving up!!!'                            &
                         ,'zoobukhov8','canopy_air_coms.f90')
      end if

      return
   end function zoobukhov8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the value of zeta for a given Richardson number, reference    !
   ! height and the roughness scale.  This is solved by using the definition of Obukhov    !
   ! length scale as stated in Louis (1979) equation (10), modified to define z/L rather   !
   ! than L.  The solution is found  iteratively since it's not a simple function to       !
   ! invert.  It tries to use Newton's method, which should take care of most cases.  In   !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real function zoobukhov_ustar(rib,zstar,rough,zoz0h,lnzoz0h,kuoustar,stable)
      use therm_lib  , only : toler  & ! intent(in)
                            , maxfpo & ! intent(in)
                            , maxit  ! ! intent(in)
      use consts_coms, only : lnexp_min ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: rib       ! Bulk Richardson number               [   ---]
      real(kind=4), intent(in) :: zstar     ! Ref. height - displacement height    [     m]
      real(kind=4), intent(in) :: rough     ! Roughness length scale               [     m]
      real(kind=4), intent(in) :: zoz0h     ! zstar/roughness(heat)                [   ---]
      real(kind=4), intent(in) :: lnzoz0h   ! ln[zstar/roughness(heat)]            [   ---]
      real(kind=4), intent(in) :: kuoustar  ! k * u / u*                           [   ---]
      logical     , intent(in) :: stable    ! Flag... This surface layer is stable [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: fh        ! lnzoz0 - psih(zeta) + psih(zeta0)    [   ---]
      real(kind=4)             :: dfhdzeta  ! d(fh)/d(zeta)                        [   ---]
      real(kind=4)             :: z0hoz     ! Roughness(heat) / Reference height   [   ---]
      real(kind=4)             :: zeta0h    ! Roughness(heat) / Obukhov length     [   ---]
      real(kind=4)             :: zetaa     ! Smallest guess (or previous guess)   [   ---]
      real(kind=4)             :: zetae     ! Current guess                        [   ---]
      real(kind=4)             :: zetaz     ! Largest guess (or Newton's guess)    [   ---]
      real(kind=4)             :: deriv     ! Function Derivative                  [   ---]
      real(kind=4)             :: fun       ! Function for which we seek a root.   [   ---]
      real(kind=4)             :: funa      ! Smallest guess function.             [   ---]
      real(kind=4)             :: funz      ! Largest guess function.              [   ---]
      real(kind=4)             :: delta0    ! Aux. var --- 2nd guess for bisection [   ---]
      real(kind=4)             :: delta     ! Aux. var --- 2nd guess for bisection [   ---]
      real(kind=4)             :: coeff     ! RiB*zstar/(Pr*(zstar-z0))*(k*u/u*)^2 [   ---]
      real(kind=4)             :: coeffi    ! 1/coeff                              [   ---]
      real(kind=4)             :: zetamin   ! Minimum zeta for stable case.        [   ---]
      real(kind=4)             :: zetamax   ! Maximum zeta for unstable case.      [   ---]
      real(kind=4)             :: zetasmall ! Number sufficiently close to zero    [   ---]
      integer                  :: itb       ! Iteration counters                   [   ---]
      integer                  :: itn       ! Iteration counters                   [   ---]
      integer                  :: itp       ! Iteration counters                   [   ---]
      logical                  :: converged ! Flag... The method converged!        [   T|F]
      logical                  :: zside     ! Flag... I'm on the z-side.           [   T|F]
      !------------------------------------------------------------------------------------!



      !----- Define some values that won't change during the iterative method. ------------!
      z0hoz = 1. / zoz0h
      coeff = rib * zstar / (tprandtl * (zstar - rough)) * kuoustar * kuoustar
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the bulk Richardson number is zero or almost zero, then we rather just      !
      ! assign z/L to be the one similar to Oncley and Dudhia (1995).  This saves time and !
      ! also avoids the risk of having zeta with the opposite sign.                        !
      !------------------------------------------------------------------------------------!
      zetasmall = coeff / lnzoz0h
      if (rib <= 0. .and. zetasmall > - toler) then
         zoobukhov_ustar = zetasmall
         return
      elseif (rib > 0. .and. zetasmall < toler) then
         zoobukhov_ustar = zetasmall / (1.0 - beta_s * rib / tprandtl)
         return
      else
         zetamin    =  toler
         zetamax    = -toler
         coeffi     = 1. / coeff
      end if
      !------------------------------------------------------------------------------------!



      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(60a1)') ('-',itn=1,60)
      !write(unit=89,fmt='(5(a,1x,f11.4,1x),a,l1)')                                        &
      !   'Input values: Rib =',rib,'zstar=',star,'rough=',rough,'zoz0=',zoz0              &
      !           ,'lnzoz0=',lnzoz0,'stable=',stable
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !------------------------------------------------------------------------------------!
      !     First guess, use Oncley and Dudhia (1995) approximation for unstable case.     !
      ! We won't use the stable case to avoid FPE or zeta with opposite sign when          !
      ! Ri is too positive.                                                                !
      !------------------------------------------------------------------------------------!
      zetaa = zetasmall
      !------------------------------------------------------------------------------------!



      !----- Find the function and its derivative. ----------------------------------------!
      zeta0h   = zetaa * z0hoz
      fh       = lnzoz0h - psih(zetaa,stable) + psih(zeta0h,stable)
      dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaa,stable)
      !------ Force derivative to be zero if the function is ill-behaved. -----------------!
      if (fh == 0.0) then
         funa  = lnexp_min
         deriv = 0.0
      else
         funa  = log(coeffi * zetaa * fh)
         deriv =  ( 1.0 / zetaa + dfhdzeta / fh )
      end if
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration. ---------------------------!
      zetaz = zetaa
      fun   = funa
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
      !   '1STGSS: itn=',0,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm          &
      !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !----- Enter Newton's method loop. --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1, maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Newton's method converges fast when it's on the right track, but there are  !
         ! cases in which it becomes ill-behaved.  Two situations are known to cause       !
         ! trouble:                                                                        !
         ! 1.  If the derivative is tiny, the next guess can be too far from the actual    !
         !     answer;                                                                     !
         ! 2.  For this specific problem, when zeta is too close to zero.  In this case    !
         !     the derivative will tend to infinity at this point and Newton's method is   !
         !     not going to perform well and can potentially enter in a weird behaviour or !
         !     lead to the wrong answer.  In any case, so we rather go with bisection.     !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler) then
            exit newloop
         elseif(stable .and. (zetaz - fun/deriv < zetamin)) then
            exit newloop
         elseif((.not. stable) .and. (zetaz - fun/deriv > zetamax)) then
            exit newloop
         end if

         !----- Copy the previous guess ---------------------------------------------------!
         zetaa = zetaz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         zetaz    = zetaa - fun/deriv
         zeta0h   = zetaz * z0hoz
         fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
         dfhdzeta = z0hoz * dpsihdzeta(zeta0h,stable) - dpsihdzeta(zetaz,stable)
         !------ Force derivative to be zero if the function is ill-behaved. --------------!
         if (zetaz == 0.0 .or. fh == 0.0) then
            fun   = lnexp_min
            deriv = 0.0
         else
            fun      = log(coeffi * zetaz * fh)
            deriv =  ( 1.0 / zetaz + dfhdzeta / fh )
         end if
         !------------------------------------------------------------------------------------!



         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'NEWTON: itn=',itn,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fm=',fm        &
         !  ,'fh=',fh,'dfmdzeta=',dfmdzeta,'dfhdzeta=',dfhdzeta,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         converged = abs(zetaz-zetaa) < toler * abs(zetaz)

         if (converged) then
            zoobukhov_ustar = 0.5 * (zetaa+zetaz)
            return
         elseif (fun == 0.0) then !---- Converged by luck. --------------------------------!
            zoobukhov_ustar = zetaz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If we reached this point then it's because Newton's method failed or it has    !
      ! become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a  !
      ! fancier bisection.  For this we need two guesses, and the guesses must have        !
      ! opposite signs.                                                                    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.0) then
         funz  = fun
         zside = .true. 
      else
         !----- zetaa becomes the minimum absolute value for searching. -------------------!
         if (stable) then
            zetaa  = zetamin
         else
            zetaa  = zetamax
         end if
         zeta0h = zetaa * z0hoz
         fh     = lnzoz0h - psih(zetaa,stable) + psih(zeta0h,stable)
         if (fh == 0.0) then
            funa = lnexp_min
         else
            funa = log ( coeffi * zetaa * fh )
         end if
         !---------------------------------------------------------------------------------!

         !------ Use a multiplicative function. -------------------------------------------!
         delta = 10.
         zetaz = zetaa
         zside = .false.
         zgssloop: do itp=1,maxfpo
            zetaz = zetaz * delta
            zeta0h   = zetaz * z0hoz
            fh       = lnzoz0h - psih(zetaz,stable) + psih(zeta0h,stable)
            if (fh == 0.0) then
               funz  = lnexp_min 
            else
               funz  = log ( coeffi * zetaz * fh )
            end if
            zside    = funa * funz < 0.0
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                 &
            !   '2NDGSS: itp=',itp,'zside=',zside,'zetaa=',zetaa,'zetaz=',zetaz             &
            !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(a)') '    No second guess for you...'
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zstar  =',zstar   ,'rough  =',rough
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'kuoust =',kuoustar,'coeff  =',coeff
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'lnzoz0h=',lnzoz0h ,'rib    =',rib
            write (unit=*,fmt='(1(a,1x,l1,1x))')     'stable =',stable
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'fun    =',fun     ,'delta  =',delta
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaa  =',zetaa   ,'funa   =',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaz  =',zetaz   ,'funz   =',funz
            call fatal_error('Failed finding the second guess for regula falsi'            &
                            ,'zoobukhov_ustar','canopy_air_coms.f90')
         end if
      end if

      !----- Now we are ready to start the regula falsi method. ---------------------------!
      bisloop: do itb=itn,maxfpo
         zetae = (funz*zetaa-funa*zetaz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(zetae-zetaa) < toler * abs(zetae)
         if (converged) exit bisloop

         !------ Update function evaluation. ----------------------------------------------!
         zeta0h   = zetae * z0hoz
         fh       = lnzoz0h - psih(zetae,stable) + psih(zeta0h,stable)
         if (fh == 0.0) then
            fun = lnexp_min
         else
            fun = log(coeffi * zetae * fh)
         end if
         !---------------------------------------------------------------------------------!

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,7(1x,a,1x,es12.5))')                       &
         !   'REGULA: itn=',itb,'bisection=',.true.,'zetaa=',zetaa,'zetaz=',zetaz,'fun=',fun   &
         !  ,'funa=',funa,'funz=',funz,'fm=',fm,'fh=',fh
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

         !------ Define new interval based on the intermediate value theorem. -------------!
         if (fun*funa < 0. ) then
            zetaz = zetae
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 0.5
            !----- We just updated zside, set zside to true. ------------------------------!
            zside = .true.
         else
            zetaa = zetae
            funa  = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 0.5
            !----- We just updated aside, set aside to true. ------------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged .or.                                                              &
          (stable .and. zetae < 0.0) .or. (.not. stable .and. zetae > 0.0)) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Zeta finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rib             [   ---] =',rib
         write (unit=*,fmt='(a,1x,f12.4)' ) 'kuoustar        [   ---] =',kuoustar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zstar           [     m] =',zstar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rough           [     m] =',rough
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0h           [   ---] =',zoz0h
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0h         [   ---] =',lnzoz0h
         write (unit=*,fmt='(a,1x,l1)'    ) 'stable          [   T|F] =',stable
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaa           [   ---] =',zetaa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaz           [   ---] =',zetaz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [   ---] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fh              [   ---] =',fh
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [   ---] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [   ---] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [   ---] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [   ---] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [   ---] ='                   &
                                                            ,abs(zetaz-zetaa)/abs(zetaz)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetae           [   ---] =',zetae
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

         call fatal_error('Zeta didn''t converge, giving up!!!'                            &
                         ,'zoobukhov_ustar','canopy_air_coms.f90')
      else
         zoobukhov_ustar = zetae
      end if

      return
   end function zoobukhov_ustar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the value of zeta for a given Richardson number, reference    !
   ! height and the roughness scale.  This is solved by using the definition of Obukhov    !
   ! length scale as stated in Louis (1979) equation (10), modified to define z/L rather   !
   ! than L.  The solution is found  iteratively since it's not a simple function to       !
   ! invert.  It tries to use Newton's method, which should take care of most cases.  In   !
   ! the unlikely case in which Newton's method fails, switch back to modified Regula      !
   ! Falsi method (Illinois).                                                              !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function zoobukhov_ustar8(rib,zstar,rough,zoz0h,lnzoz0h,kuoustar,stable)
      use therm_lib8 , only : toler8     & ! intent(in)
                            , maxfpo     & ! intent(in)
                            , maxit      ! ! intent(in)
      use consts_coms, only : lnexp_min8 ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: rib       ! Bulk Richardson number               [   ---]
      real(kind=8), intent(in) :: zstar     ! Reference height - displacement hgt. [     m]
      real(kind=8), intent(in) :: rough     ! Roughness length scale               [     m]
      real(kind=8), intent(in) :: zoz0h     ! zstar/roughness(heat)                [   ---]
      real(kind=8), intent(in) :: lnzoz0h   ! ln[zstar/roughness(heat)]            [   ---]
      real(kind=8), intent(in) :: kuoustar  ! k * u / u*                           [   ---]
      logical     , intent(in) :: stable    ! Flag... This surface layer is stable [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: fh        ! lnzoz0 - psih(zeta) + psih(zeta0)    [   ---]
      real(kind=8)             :: dfhdzeta  ! d(fh)/d(zeta)                        [   ---]
      real(kind=8)             :: z0hoz     ! Roughness(heat) / Reference height   [   ---]
      real(kind=8)             :: zeta0h    ! Roughness(heat) / Obukhov length     [   ---]
      real(kind=8)             :: zetaa     ! Smallest guess (or previous guess)   [   ---]
      real(kind=8)             :: zetae     ! Current guess                        [   ---]
      real(kind=8)             :: zetaz     ! Largest guess (new guess in Newton)  [   ---]
      real(kind=8)             :: deriv     ! Function Derivative                  [   ---]
      real(kind=8)             :: fun       ! Function for which we seek a root.   [   ---]
      real(kind=8)             :: funa      ! Smallest guess function.             [   ---]
      real(kind=8)             :: funz      ! Largest guess function.              [   ---]
      real(kind=8)             :: delta0    ! Aux. var --- 2nd guess for bisection [   ---]
      real(kind=8)             :: delta     ! Aux. var --- 2nd guess for bisection [   ---]
      real(kind=8)             :: coeff     ! RiB*zstar/(Pr*(zstar-z0))*(k*u/u*)^2 [   ---]
      real(kind=8)             :: coeffi    ! 1.d0 / coeff                         [   ---]
      real(kind=8)             :: zetamin   ! Minimum zeta for stable case.        [   ---]
      real(kind=8)             :: zetamax   ! Maximum zeta for unstable case.      [   ---]
      real(kind=8)             :: zetasmall ! Zeta dangerously close to zero       [   ---]
      integer                  :: itn       ! Iteration counters                   [   ---]
      integer                  :: itb       ! Iteration counters                   [   ---]
      integer                  :: itp       ! Iteration counters                   [   ---]
      logical                  :: converged ! Flag... The method converged!        [   T|F]
      logical                  :: zside     ! Flag... I'm on the z-side.           [   T|F]
      !----- Local constants. -------------------------------------------------------------!
      logical, parameter       :: debug=.false. ! Print debugging information?
      !------------------------------------------------------------------------------------!



      !----- Define some values that won't change during the iterative method. ------------!
      z0hoz = 1.d0 / zoz0h
      coeff = rib * zstar / (tprandtl8 * (zstar - rough)) * kuoustar * kuoustar
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If the bulk Richardson number is zero or almost zero, then we rather just      !
      ! assign z/L to be the one similar to Oncley and Dudhia (1995).  This saves time and !
      ! also avoids the risk of having zeta with the opposite sign.                        !
      !------------------------------------------------------------------------------------!
      zetasmall = coeff / lnzoz0h
      if (rib <= 0.d0 .and. zetasmall > - toler8) then
         zoobukhov_ustar8 = zetasmall
         return
      elseif (rib > 0.d0 .and. zetasmall < toler8) then
         zoobukhov_ustar8 = zetasmall / (1.d0 - beta_s8 * rib / tprandtl8)
         return
      else
         zetamin    =  toler8
         zetamax    = -toler8
         coeffi     = 1.d0 / coeff
      end if
      !------------------------------------------------------------------------------------!



      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write(unit=89,fmt='(60a1)') ('-',itn=1,60)
         write(unit=89,fmt='(6(a,1x,f11.4,1x),a,l1)')                                      &
            'Input values: Rib =',rib,'zstar=',zstar,'rough=',rough,'zoz0h=',zoz0h         &
                    ,'lnzoz0h=',lnzoz0h,'kuoustar=',kuoustar,'stable=',stable
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!



      !------------------------------------------------------------------------------------!
      !     First guess, use Oncley and Dudhia (1995) approximation for unstable case.     !
      ! We won't use the stable case to avoid FPE or zeta with opposite sign when          !
      ! Ri is too positive.                                                                !
      !------------------------------------------------------------------------------------!
      zetaa = zetasmall
      !------------------------------------------------------------------------------------!



      !----- Find the function and its derivative. ----------------------------------------!
      zeta0h   = zetaa * z0hoz
      fh       = lnzoz0h - psih8(zetaa,stable) + psih8(zeta0h,stable)
      dfhdzeta = z0hoz * dpsihdzeta8(zeta0h,stable) - dpsihdzeta8(zetaa,stable)
      if (fh == 0.d0) then
         funa  =  lnexp_min8
         deriv =  0.d0
      else
         funa  = log(coeffi * zetaa * fh)
         deriv = 1.d0 / zetaa + dfhdzeta / fh
      end if
      !------------------------------------------------------------------------------------!


      !----- Copy just in case it fails at the first iteration. ---------------------------!
      zetaz = zetaa
      fun   = funa
      !------------------------------------------------------------------------------------!



      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (debug) then
         write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,5(1x,a,1x,es12.5))')                    &
            '1STGSS: itn=',0,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fh=',fh       &
                   ,'dfhdzeta=',dfhdzeta,'deriv=',deriv
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- Enter Newton's method loop. --------------------------------------------------!
      converged = .false.
      newloop: do itn = 1, maxfpo/6
         !---------------------------------------------------------------------------------!
         !     Newton's method converges fast when it's on the right track, but there are  !
         ! cases in which it becomes ill-behaved.  Two situations are known to cause       !
         ! trouble:                                                                        !
         ! 1.  If the derivative is tiny, the next guess can be too far from the actual    !
         !     answer;                                                                     !
         ! 2.  For this specific problem, when zeta is too close to zero.  In this case    !
         !     the derivative will tend to infinity at this point and Newton's method is   !
         !     not going to perform well and can potentially enter in a weird behaviour or !
         !     lead to the wrong answer.  In any case, so we rather go with bisection.     !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler8) then
            exit newloop
         elseif(stable .and. (zetaz - fun/deriv < zetamin)) then
            exit newloop
         elseif((.not. stable) .and. (zetaz - fun/deriv > zetamax)) then
            exit newloop
         end if

         !----- Copying the previous guess ------------------------------------------------!
         zetaa = zetaz
         funa  = fun
         !----- New guess, its function and derivative evaluation -------------------------!
         zetaz    = zetaa - fun/deriv
         zeta0h   = zetaz * z0hoz
         fh       = lnzoz0h - psih8(zetaz,stable) + psih8(zeta0h,stable)
         dfhdzeta = z0hoz * dpsihdzeta8(zeta0h,stable) - dpsihdzeta8(zetaz,stable)
         if (fh == 0.d0) then
            fun   =  lnexp_min8
            deriv =  0.d0
         else
            fun   = log(coeffi * zetaz * fh)
            deriv = 1.d0 / zetaz + dfhdzeta / fh
         end if
         !---------------------------------------------------------------------------------!


         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (debug) then
            write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,5(1x,a,1x,es12.5))')                 &
               'NEWTON: itn=',itn,'bisection=',.false.,'zetaz=',zetaz,'fun=',fun,'fh=',fh  &
              ,'dfhdzeta=',dfhdzeta,'deriv=',deriv
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         converged = abs(zetaz-zetaa) < toler8 * abs(zetaz)

         if (converged) then
            zoobukhov_ustar8 = 5.d-1 * (zetaa+zetaz)
            return
         elseif (fun == 0.d0) then !---- Converged by luck. -------------------------------!
            zoobukhov_ustar8 = zetaz
            return
         end if
      end do newloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If we reached this point then it's because Newton's method failed or it has    !
      ! become too dangerous.  We use the Regula Falsi (Illinois) method, which is just a  !
      ! fancier bisection.  For this we need two guesses, and the guesses must have        !
      ! opposite signs.                                                                    !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.d0) then
         funz  = fun
         zside = .true. 
      else
         !----- zetaa becomes the minimum absolute value for searching. -------------------!
         if (stable) then
            zetaa  = zetamin
         else
            zetaa  = zetamax
         end if
         zeta0h = zetaa * z0hoz
         fh     = lnzoz0h - psih8(zetaa,stable) + psih8(zeta0h,stable)
         if (fh == 0.d0) then
            funa = lnexp_min8
         else
            funa = log ( coeffi * zetaa * fh )
         end if
         !---------------------------------------------------------------------------------!

         !------ Use a multiplicative function. -------------------------------------------!
         delta = 1.d1
         zetaz = zetaa
         zside = .false.
         zgssloop: do itp=1,maxfpo
            zetaz = zetaz * delta
            zeta0h   = zetaz * z0hoz
            fh       = lnzoz0h - psih8(zetaz,stable) + psih8(zeta0h,stable)
            if (fh == 0.d0) then
               funz  = lnexp_min8
            else
               funz  = log ( coeffi * zetaz * fh )
            end if
            zside    = funa * funz < 0.d0
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (debug) then
               write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,6(1x,a,1x,es12.5))')              &
                  '2NDGSS: itp=',itp,'zside=',zside,'zetaa=',zetaa,'zetaz=',zetaz          &
                 ,'funa=',funa,'funz=',funz,'fh=',fh,'delta=',delta
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(a)') '    No second guess for you...'
            write (unit=*,fmt='(a)') '=================================================='
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zstar  =',zstar   ,'rough  =',rough
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'kuoust =',kuoustar,'coeff  =',coeff
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'lnzoz0h=',lnzoz0h ,'rib    =',rib
            write (unit=*,fmt='(1(a,1x,l1,1x))')     'stable =',stable
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'fun    =',fun     ,'delta  =',delta
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaa  =',zetaa   ,'funa   =',funa
            write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'zetaz  =',zetaz   ,'funz   =',funz
            write (unit=*,fmt='(a)') '=================================================='
            call fatal_error('Failed finding the second guess for regula falsi'            &
                            ,'zoobukhov_ustar8','canopy_air_coms.f90')
         end if
      end if

      !----- Now we are ready to start the regula falsi method. ---------------------------!
      bisloop: do itb=itn,maxfpo
         zetae = (funz*zetaa-funa*zetaz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = abs(zetae-zetaa) < toler8 * abs(zetae)
         if (converged) exit bisloop

         !------ Finding the new function -------------------------------------------------!
         zeta0h   = zetae * z0hoz
         fh       = lnzoz0h - psih8(zetae,stable) + psih8(zeta0h,stable)
         if (fh == 0.d0) then
            fun  =  lnexp_min8
         else
            fun  = log(coeffi * zetae * fh)
         end if


         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (debug) then
            write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,6(1x,a,1x,es12.5))')                 &
               'REGULA: itn=',itb,'bisection=',.true.,'zetaa=',zetaa,'zetaz=',zetaz        &
               ,'fun=',fun,'funa=',funa,'funz=',funz,'fh=',fh
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0.d0 ) then
            zetaz = zetae
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 5.d-1
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            zetaa = zetae
            funa  = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 5.d-1
            !----- We just updated aside, setting aside to true. --------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged .or.                                                              &
          (stable .and. zetae < 0.d0) .or. (.not. stable .and. zetae > 0.d0)) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Zeta finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rib             [   ---] =',rib
         write (unit=*,fmt='(a,1x,f12.4)' ) 'kuoustar        [   ---] =',kuoustar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zstar           [     m] =',zstar
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rough           [     m] =',rough
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zoz0h           [   ---] =',zoz0h
         write (unit=*,fmt='(a,1x,f12.4)' ) 'lnzoz0h         [   ---] =',lnzoz0h
         write (unit=*,fmt='(a,1x,l1)'    ) 'stable          [   T|F] =',stable
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (zoobukhov8 values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaa           [   ---] =',zetaa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetaz           [   ---] =',zetaz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [   ---] =',fun
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fh              [   ---] =',fh
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [   ---] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [   ---] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [   ---] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [   ---] =',toler8
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [   ---] ='                   &
                                                            ,abs(zetaz-zetaa)/abs(zetaz)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'zetae           [   ---] =',zetae
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

         call fatal_error('Zeta didn''t converge, giving up!!!'                            &
                         ,'zoobukhov_ustar8','canopy_air_coms.f90')
      else

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (debug) then
            write(unit=89,fmt='(a,1x,i5,1x,a,1x,l1,1x,6(1x,a,1x,es12.5))')                 &
               'ANSWER: itn=',itb,'converged=',.true.,'zetaa=',zetaa,'zetaz=',zetaz        &
               ,'fun=',fun,'funa=',funa,'funz=',funz,'fh=',fh
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         zoobukhov_ustar8 = zetae
      end if

      return
   end function zoobukhov_ustar8
   !=======================================================================================!
   !=======================================================================================!
end module canopy_air_coms
!==========================================================================================!
!==========================================================================================!
