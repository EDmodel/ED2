!==========================================================================================!
!==========================================================================================!
!   Module ed_nbg_init
!
!> \brief   Sub-routines that will initialise ED-2 with near-bare ground conditions
!> \details A group of sub-routines that provide initial near-bare ground conditions
!>          for three cases: standard ED simulation, "big-leaf" ED, and test with 
!>          fixed layers (for code testing only).
!> \author  11 Oct 2017 Turned file into a module, Manfredo di Porcia e Brugnera
!> \author  18 Oct 2017 Minor changes to include bark as an active tissue, and to 
!           define initial storage pool as a function of balive, Marcos Longo
!------------------------------------------------------------------------------------------!
module ed_nbg_init
   contains

   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine initializes a near-bare ground polygon.                          !
   !---------------------------------------------------------------------------------------!
   subroutine near_bare_ground_init(cgrid)
      use ed_state_vars  , only : edtype              & ! structure
                                , polygontype         & ! structure
                                , sitetype            & ! structure
                                , allocate_sitetype   ! ! subroutine
      use grid_coms      , only : nzg                 ! ! intent(in)
      use ed_misc_coms   , only : ied_init_mode       ! ! intent(in)
      use physiology_coms, only : n_plant_lim         ! ! intent(in)
      use ed_type_init   , only : init_ed_patch_vars  & ! subroutine
                                , init_ed_site_vars   & ! subroutine
                                , init_ed_poly_vars   ! ! subroutine
      use decomp_coms    , only : decomp_scheme       & ! intent(in)
                                , agf_fsc             & ! intent(in)
                                , agf_stsc            & ! intent(in)
                                , c2n_fast_0          & ! intent(in)
                                , c2n_structural      & ! intent(in)
                                , c2n_slow            & ! intent(in)
                                , f0_msc              & ! intent(in)
                                , f0_ssc              & ! intent(in)
                                , f0_psc              & ! intent(in)
                                , nbg_nlim_fsc        & ! intent(in)
                                , nbg_nlim_stsc       & ! intent(in)
                                , nbg_nlim_ssc        ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)      , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype) , pointer :: cpoly
      type(sitetype)    , pointer :: csite
      integer                     :: ipy
      integer                     :: isi
      !------------------------------------------------------------------------------------!


      !----- Big loop ---------------------------------------------------------------------!
      do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            !----- We start with a single patch per site, always with primary vegetation. -!
            csite%npatches = 1
            call allocate_sitetype(csite,1)
            csite%dist_type          (1) = 3
            csite%age                (1) = 0.0
            csite%area               (1) = 1.0
            csite%fbeam              (1) = 1.0
            csite%light_type         (1) = 1


            !------------------------------------------------------------------------------!
            !     Someone that uses the nitrogen model should check whether this is neces- !
            ! sary or not.  If nitrogen limitation is off, then we start all carbon and    !
            ! nitrogen pools with zeroes, otherwise we initialise with the former default  !
            ! values.  Values have been updated so C and N pools necessarily start         !
            ! following stoichiometry.                                                     !
            !------------------------------------------------------------------------------!
            select case (n_plant_lim)
            case (0)
               csite%fast_grnd_C       (1) = 0.0
               csite%fast_soil_C       (1) = 0.0
               csite%structural_grnd_C (1) = 0.0
               csite%structural_soil_C (1) = 0.0
               csite%structural_grnd_L (1) = 0.0
               csite%structural_soil_L (1) = 0.0
               csite%microbial_soil_C  (1) = 0.0
               csite%slow_soil_C       (1) = 0.0
               csite%passive_soil_C    (1) = 0.0
               csite%fast_grnd_N       (1) = 0.0
               csite%fast_soil_N       (1) = 0.0
               csite%structural_grnd_N (1) = 0.0
               csite%structural_soil_N (1) = 0.0
               csite%mineralized_soil_N(1) = 0.0
            case (1)
               csite%fast_grnd_C        (1) =        agf_fsc   * nbg_nlim_fsc
               csite%fast_soil_C        (1) = (1.0 - agf_fsc)  * nbg_nlim_fsc
               csite%structural_grnd_C  (1) =        agf_stsc  * nbg_nlim_stsc
               csite%structural_soil_C  (1) = (1.0 - agf_stsc) * nbg_nlim_stsc
               csite%structural_grnd_L  (1) = csite%structural_grnd_C(1)
               csite%structural_soil_L  (1) = csite%structural_soil_C(1)
               select case (decomp_scheme)
               case (5)
                  csite%microbial_soil_C(1) = f0_msc * nbg_nlim_ssc
                  csite%slow_soil_C     (1) = f0_ssc * nbg_nlim_ssc
                  csite%passive_soil_C  (1) = f0_psc * nbg_nlim_ssc
               case default
                  csite%microbial_soil_C(1) = 0.0
                  csite%slow_soil_C     (1) = nbg_nlim_ssc
                  csite%passive_soil_C  (1) = 0.0
               end select
               csite%fast_grnd_N       (1) = csite%fast_grnd_C       (1)   / c2n_fast_0
               csite%fast_soil_N       (1) = csite%fast_soil_C       (1)   / c2n_fast_0
               csite%structural_grnd_N (1) = csite%structural_grnd_C (1)   / c2n_structural
               csite%structural_soil_N (1) = csite%structural_soil_C (1)   / c2n_structural
               csite%mineralized_soil_N(1) = ( csite%microbial_soil_C(1)                   &
                                             + csite%slow_soil_C     (1)                   &
                                             + csite%passive_soil_C  (1) ) / c2n_slow
            end select
            !------------------------------------------------------------------------------!

            csite%sum_dgd            (1) = 0.0
            csite%sum_chd            (1) = 0.0
            csite%plant_ag_biomass   (1) = 0.

            !----- We now populate the cohorts with near bare ground condition. -----------!
            select case (ied_init_mode)
            case (-8)
               call init_cohorts_by_layers(csite,cpoly%lsl(isi),1,csite%npatches           &
                                          ,nzg,cpoly%ntext_soil(:,isi))
            case (-1,0)
               call init_nbg_cohorts(csite,cpoly%lsl(isi),1,csite%npatches                 &
                                    ,nzg,cpoly%ntext_soil(:,isi))
            end select

            !----- Initialise the patches now that cohorts are there. ---------------------!
            call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
         end do
         !----- Initialise some site-level variables. -------------------------------------!
         call init_ed_site_vars(cpoly)
      end do

      !----- Last, but not the least, the polygons. ---------------------------------------!
      call init_ed_poly_vars(cgrid)

      return
   end subroutine near_bare_ground_init
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine assigns a near-bare ground (NBG) state for some patches.         !
   !---------------------------------------------------------------------------------------!
   subroutine init_nbg_cohorts(csite,lsl,ipa_a,ipa_z,mzg,ntext_soil)
      use ed_state_vars      , only : edtype              & ! structure
                                    , polygontype         & ! structure
                                    , sitetype            & ! structure
                                    , patchtype           & ! structure
                                    , allocate_sitetype   & ! subroutine
                                    , allocate_patchtype  ! ! subroutine
      use ed_max_dims        , only : n_pft               ! ! intent(in)
      use ed_misc_coms       , only : ied_init_mode       ! ! intent(in)
      use pft_coms           , only : q                   & ! intent(in)
                                    , qsw                 & ! intent(in)
                                    , qbark               & ! intent(in)
                                    , sla                 & ! intent(in)
                                    , hgt_min             & ! intent(in)
                                    , include_pft         & ! intent(in)
                                    , init_density        & ! intent(in)
                                    , include_pft_pt      & ! intent(in)
                                    , include_pft_fp      & ! intent(in)
                                    , include_pft_ag      & ! intent(in)
                                    , f_bstorage_init     & ! intent(in)
                                    , agf_bs              ! ! intent(in)
      use consts_coms        , only : t3ple               & ! intent(in)
                                    , pio4                & ! intent(in)
                                    , almost_zero         & ! intent(in)
                                    , kgom2_2_tonoha      & ! intent(in)
                                    , tonoha_2_kgom2      ! ! intent(in)
      use physiology_coms    , only : iddmort_scheme      ! ! intent(in)
      use allometry          , only : h2dbh               & ! function
                                    , size2bd             & ! function
                                    , size2bl             & ! function
                                    , size2bt             & ! function
                                    , size2xb             & ! function
                                    , ed_balive           & ! function
                                    , ed_biomass          & ! function
                                    , area_indices        ! ! subroutine
      use fuse_fiss_utils    , only : sort_cohorts        ! ! subroutine
      use ed_type_init       , only : init_ed_cohort_vars & ! subroutine
                                    , init_ed_patch_vars  & ! subroutine
                                    , init_ed_site_vars   & ! subroutine
                                    , init_ed_poly_vars   ! ! subroutine

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target     :: csite             ! Current site
      integer               , intent(in) :: lsl               ! Lowest soil level
      integer               , intent(in) :: ipa_a             ! 1st patch
      integer               , intent(in) :: ipa_z             ! Last patch
      integer               , intent(in) :: mzg               ! Number of soil layers
      integer,dimension(mzg), intent(in) :: ntext_soil        ! Soil texture profile
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)       , pointer    :: cpatch            ! Current patch
      real                               :: bdeadx            ! Heartwood biomass
      integer                            :: ipa               ! Patch number
      integer                            :: ico               ! Cohort counter
      integer                            :: mypfts            ! Number of included PFTs
      integer                            :: ipft              ! PFT counter
      !------------------------------------------------------------------------------------!

      !----- Patch loop. ------------------------------------------------------------------!
      patchloop: do ipa=ipa_a,ipa_z
         cpatch => csite%patch(ipa)

         !----- Check which initialisation we are going to use. ---------------------------!
         select case (ied_init_mode)
         case (-1) !------ True bare ground simulation (absolute desert). -----------------!
            mypfts = 0
            return
         case ( 0) !------ Nearly bare ground simulation (start with a few seedlings). ----!

            !------------------------------------------------------------------------------!
            !     Decide how many cohorts to allocate.                                     !
            !------------------------------------------------------------------------------!
            select case (csite%dist_type(ipa))
            case (1)     !---- Pasture. ---------------------------------------------------!
               mypfts = count(include_pft_pt)
            case (2)     !---- Forest plantation. -----------------------------------------!
               mypfts = count(include_pft_fp)
            case (8)     !---- Cropland. --------------------------------------------------!
               mypfts = count(include_pft_ag)
            case default !---- Logged or non-managed forests. -----------------------------!
               mypfts = count(include_pft)
            end select
            !------------------------------------------------------------------------------!
         end select

         !----- Perform cohort allocation. ------------------------------------------------!
         call allocate_patchtype(cpatch,mypfts)

         !---------------------------------------------------------------------------------!
         !    Here we loop over PFTs rather than assigning the cohorts, so we ensure to    !
         ! only include the PFTs that should be  included (i.e., patches that should not   !
         ! happen in agricultural patches will be skipped, and agricultural PFTs will not  !
         ! be included in the forest sites.                                                !
         !---------------------------------------------------------------------------------!
         ico = 0
         pftloop: do ipft = 1,n_pft
            select case (csite%dist_type(ipa))
            case (1)
               !----- Pasture, only the ones allowed are included. ------------------------!
               if (include_pft_pt(ipft)) then
                  ico = ico + 1
               else
                  cycle pftloop
               end if
               !---------------------------------------------------------------------------!
            case (2)
               !----- Forest plantation, only the ones allowed are included. --------------!
               if (include_pft_fp(ipft)) then
                  ico = ico + 1
               else
                  cycle pftloop
               end if
               !---------------------------------------------------------------------------!
            case (8)
               !----- Cropland, only the ones allowed are included. -----------------------!
               if (include_pft_ag(ipft)) then
                  ico = ico + 1
               else
                  cycle pftloop
               end if
               !---------------------------------------------------------------------------!
            case default
               !----- Primary or secondary vegetation, only allowed PFTs are included. ----!
               if (include_pft(ipft)) then
                  ico = ico + 1
               else
                  cycle pftloop
               end if
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!

            !----- The PFT is the plant functional type. ----------------------------------!
            cpatch%pft(ico)              = ipft
            !------------------------------------------------------------------------------!
            !     Define the near-bare ground state using the standard minimum height and  !
            ! minimum plant density.  We assume all NBG PFTs to have leaves fully flushed, !
            ! but with no storage biomass.  We then compute the other biomass quantities   !
            ! using the standard allometry for this PFT.                                   !
            !------------------------------------------------------------------------------!
            cpatch%nplant(ico)           = init_density(ipft)
            cpatch%height(ico)           = hgt_min(ipft)
            cpatch%phenology_status(ico) = 0
            cpatch%dbh(ico)              = h2dbh(cpatch%height(ico),ipft)
            bdeadx                       = size2bd(cpatch%dbh(ico),cpatch%height(ico),ipft)
            cpatch%bdeada(ico)           =       agf_bs(ipft)  * bdeadx
            cpatch%bdeadb(ico)           = (1. - agf_bs(ipft)) * bdeadx
            cpatch%sla(ico)              = sla(ipft)
            cpatch%bleaf(ico)            = size2bl(cpatch%dbh(ico),cpatch%height(ico)      &
                                                  ,cpatch%sla(ico),ipft)
            cpatch%broot(ico)            = q(ipft) * cpatch%bleaf(ico)
            cpatch%bsapwooda(ico)        = agf_bs(ipft) * cpatch%bleaf(ico)                &
                                         * qsw(ipft) * cpatch%height(ico)
            cpatch%bsapwoodb(ico)        = (1.-agf_bs(ipft)) * cpatch%bleaf(ico)           &
                                         * qsw(ipft) * cpatch%height(ico)
            cpatch%bbarka(ico)           = agf_bs(ipft) * cpatch%bleaf(ico)                &
                                         * qbark(ipft) * cpatch%height(ico)
            cpatch%bbarkb(ico)           = (1.-agf_bs(ipft)) * cpatch%bleaf(ico)           &
                                         * qbark(ipft) * cpatch%height(ico)
            cpatch%balive(ico)           = ed_balive(cpatch,ico)
            cpatch%bstorage(ico)         = max(almost_zero,f_bstorage_init(ipft))          &
                                         * cpatch%balive(ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Initialise the carbon balance.  For initial conditions, we always assume !
            ! storage biomass for the previous months so the scale is correct (carbon      !
            ! balance is given in kgC/pl).  The current month carbon balance must be ini-  !
            ! tialised consistently with the iddmort_scheme set by the user.               !
            !------------------------------------------------------------------------------!
            cpatch%cb         (1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_lightmax(1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_moistmax(1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_mlmax   (1:12,ico) = cpatch%bstorage(ico)
            select case (iddmort_scheme)
            case (0)
               !------ Storage is not accounted. ------------------------------------------!
               cpatch%cb         (13,ico) = 0.0
               cpatch%cb_lightmax(13,ico) = 0.0
               cpatch%cb_moistmax(13,ico) = 0.0
               cpatch%cb_mlmax   (13,ico) = 0.0
               !---------------------------------------------------------------------------!
            case (1)
               !------ Storage is accounted. ----------------------------------------------!
               cpatch%cb         (13,ico) = cpatch%bstorage(ico)
               cpatch%cb_lightmax(13,ico) = cpatch%bstorage(ico)
               cpatch%cb_moistmax(13,ico) = cpatch%bstorage(ico)
               cpatch%cb_mlmax   (13,ico) = cpatch%bstorage(ico)
               !---------------------------------------------------------------------------!
            end select
            cpatch%cbr_bar          (ico) = 1.0
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !    Find derived properties: leaf and wood area indices, above-ground         !
            ! biomass, basal area, timber stocks and bark thickness.                       !
            !------------------------------------------------------------------------------!
            call area_indices(cpatch, ico)
            cpatch%agb    (ico) = ed_biomass(cpatch, ico)
            cpatch%basarea(ico) = pio4 * cpatch%dbh(ico)*cpatch%dbh(ico)
            cpatch%btimber(ico) = size2bt(cpatch%dbh(ico),cpatch%height(ico)               &
                                         ,cpatch%bdeada(ico),cpatch%bsapwooda(ico)         &
                                         ,cpatch%bbarka(ico),cpatch%pft(ico))
            cpatch%thbark (ico) = size2xb(cpatch%dbh(ico),cpatch%height(ico)               &
                                         ,cpatch%bbarka(ico),cpatch%bbarkb(ico)            &
                                         ,cpatch%sla(ico),cpatch%pft(ico))
            !------------------------------------------------------------------------------!


            !----- Initialize other cohort-level variables. -------------------------------!
            call init_ed_cohort_vars(cpatch,ico,lsl,mzg,ntext_soil)
            !------------------------------------------------------------------------------!

            !----- Update total patch-level above-ground biomass --------------------------!
            csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                      &
                                        + cpatch%nplant(ico) * cpatch%agb(ico)
            !------------------------------------------------------------------------------!
         end do pftloop

         !---------------------------------------------------------------------------------!
         !     Because initial heights may not be constant, we must sort the cohorts.      !
         !---------------------------------------------------------------------------------!
         call sort_cohorts(cpatch)
         !---------------------------------------------------------------------------------!
      end do patchloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine init_nbg_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine assigns a near-bare ground (NBG) state for some patches.         !
   !---------------------------------------------------------------------------------------!
   subroutine init_cohorts_by_layers(csite,lsl,ipa_a,ipa_z,mzg,ntext_soil)
      use ed_state_vars      , only : edtype              & ! structure
                                    , polygontype         & ! structure
                                    , sitetype            & ! structure
                                    , patchtype           & ! structure
                                    , allocate_sitetype   & ! subroutine
                                    , allocate_patchtype  ! ! subroutine
      use ed_max_dims        , only : n_pft               ! ! intent(in)
      use pft_coms           , only : q                   & ! intent(in)
                                    , qsw                 & ! intent(in)
                                    , qbark               & ! intent(in)
                                    , sla                 & ! intent(in)
                                    , include_pft         & ! intent(in)
                                    , include_these_pft   & ! intent(in)
                                    , f_bstorage_init     & ! intent(in)
                                    , agf_bs              ! ! intent(in)
      use consts_coms        , only : t3ple               & ! intent(in)
                                    , pio4                & ! intent(in)
                                    , almost_zero         & ! intent(in)
                                    , kgom2_2_tonoha      & ! intent(in)
                                    , tonoha_2_kgom2      ! ! intent(in)
      use physiology_coms    , only : iddmort_scheme      ! ! intent(in)
      use allometry          , only : h2dbh               & ! function
                                    , size2bd             & ! function
                                    , size2bl             & ! function
                                    , size2bt             & ! function
                                    , size2xb             & ! function
                                    , ed_balive           & ! function
                                    , ed_biomass          & ! function
                                    , area_indices        ! ! subroutine
      use fuse_fiss_utils    , only : sort_cohorts        ! ! subroutine
      use ed_type_init       , only : init_ed_cohort_vars & ! subroutine
                                    , init_ed_patch_vars  & ! subroutine
                                    , init_ed_site_vars   & ! subroutine
                                    , init_ed_poly_vars   ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target     :: csite      ! Current site
      integer               , intent(in) :: lsl        ! Lowest soil level
      integer               , intent(in) :: ipa_a      ! 1st patch assigned with NBG stat
      integer               , intent(in) :: ipa_z      ! Last patch assigned with NBG state
      integer               , intent(in) :: mzg        ! Number of soil layers
      integer,dimension(mzg), intent(in) :: ntext_soil ! Soil texture profile
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch        ! Current patch
      integer                     :: ipa           ! Patch number
      integer                     :: ico           ! Cohort counter
      integer                     :: ipft          ! PFT counter
      real                        :: bdeadx        ! Total heartwood biomass
      real                        :: height        ! Cohort initial height
      !----- Local constants. -------------------------------------------------------------!
      integer        , parameter  :: nlayers = 8   ! # of cohort layers to be included.
      real           , parameter  :: dheight = 1.5 ! height interval.
      real           , parameter  :: lai0    = 1.  ! Initial LAI for each layer.
      real           , parameter  :: h0      = 1.5 ! Height of the lowest cohort.
      !------------------------------------------------------------------------------------!

      !----- Patch loop. ------------------------------------------------------------------!
      patchloop: do ipa=ipa_a,ipa_z
         cpatch => csite%patch(ipa)

         if (count (include_pft) /= 1) then
            call fatal_error('Multi-layer run cannot be run with more than 1 PFT...'       &
                            ,'init_cohorts_by_layers','ed_nbg_init.f90')
         end if

         !----- Assigning the PFT. --------------------------------------------------------!
         ipft = include_these_pft(1)

         !----- Perform cohort allocation. ------------------------------------------------!
         call allocate_patchtype(cpatch,nlayers)

         !---------------------------------------------------------------------------------!
         !    Here we loop over layers and assign the new cohort.                          !
         !---------------------------------------------------------------------------------!
         height = h0 - dheight
         layerloop: do ico = 1,nlayers
            height = height + dheight

            !----- The PFT is the plant functional type. ----------------------------------!
            cpatch%pft(ico)              = ipft

            !------------------------------------------------------------------------------!
            !     Define the initial state in such a way that the LAI is always the        !
            ! initial LAI.  We then compute the other biomass quantities using the stand-  !
            ! ard allometry for this PFT.                                                  !
            !------------------------------------------------------------------------------!
            cpatch%height(ico)           = height
            cpatch%phenology_status(ico) = 0
            cpatch%dbh(ico)              = h2dbh(cpatch%height(ico),ipft)
            bdeadx                       = size2bd(cpatch%dbh(ico),cpatch%height(ico),ipft)
            cpatch%bdeada(ico)           =       agf_bs(ipft)  * bdeadx
            cpatch%bdeadb(ico)           = (1. - agf_bs(ipft)) * bdeadx
            cpatch%sla(ico)              = sla(ipft)
            cpatch%bleaf(ico)            = size2bl(cpatch%dbh(ico),cpatch%height(ico)      &
                                                  ,cpatch%sla(ico),ipft)
            cpatch%broot(ico)            = q(ipft) * cpatch%bleaf(ico)
            cpatch%bsapwooda(ico)        = agf_bs(ipft) * cpatch%bleaf(ico)                &
                                         * qsw(ipft) * cpatch%height(ico)
            cpatch%bsapwoodb(ico)        = (1.-agf_bs(ipft)) * cpatch%bleaf(ico)           &
                                         * qsw(ipft) * cpatch%height(ico)
            cpatch%bbarka(ico)           = agf_bs(ipft) * cpatch%bleaf(ico)                &
                                         * qbark(ipft) * cpatch%height(ico)
            cpatch%bbarkb(ico)           = (1.-agf_bs(ipft)) * cpatch%bleaf(ico)           &
                                         * qbark(ipft) * cpatch%height(ico)
            cpatch%balive(ico)           = ed_balive(cpatch,ico)
            cpatch%bstorage(ico)         = max(almost_zero,f_bstorage_init(ipft))          &
                                         * cpatch%balive(ico)

            !------------------------------------------------------------------------------!
            !     Initialise the carbon balance.  For initial conditions, we always assume !
            ! storage biomass for the previous months so the scale is correct (carbon      !
            ! balance is given in kgC/pl).  The current month carbon balance must be ini-  !
            ! tialised consistently with the iddmort_scheme set by the user.               !
            !------------------------------------------------------------------------------!
            cpatch%cb         (1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_lightmax(1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_moistmax(1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_mlmax   (1:12,ico) = cpatch%bstorage(ico)
            select case (iddmort_scheme)
            case (0)
               !------ Storage is not accounted. ------------------------------------------!
               cpatch%cb         (13,ico) = 0.0
               cpatch%cb_lightmax(13,ico) = 0.0
               cpatch%cb_moistmax(13,ico) = 0.0
               cpatch%cb_mlmax   (13,ico) = 0.0
               !---------------------------------------------------------------------------!
            case (1)
               cpatch%cb         (13,ico) = cpatch%bstorage(ico)
               cpatch%cb_lightmax(13,ico) = cpatch%bstorage(ico)
               cpatch%cb_moistmax(13,ico) = cpatch%bstorage(ico)
               cpatch%cb_mlmax   (13,ico) = cpatch%bstorage(ico)
            end select
            cpatch%cbr_bar          (ico) = 1.0
            !------------------------------------------------------------------------------!

            !----- NPlant is defined such that the cohort LAI is equal to LAI0
            cpatch%nplant(ico)           = lai0 / (cpatch%bleaf(ico) * cpatch%sla(ico))




            !------------------------------------------------------------------------------!
            !    Find derived properties: leaf and wood area indices, above-ground         !
            ! biomass, basal area, timber stocks and bark thickness.                       !
            !------------------------------------------------------------------------------!
            call area_indices(cpatch, ico)
            cpatch%agb    (ico) = ed_biomass(cpatch, ico)
            cpatch%basarea(ico) = pio4 * cpatch%dbh(ico)*cpatch%dbh(ico)
            cpatch%btimber(ico) = size2bt(cpatch%dbh(ico),cpatch%height(ico)               &
                                         ,cpatch%bdeada(ico),cpatch%bsapwooda(ico)         &
                                         ,cpatch%bbarka(ico),cpatch%pft(ico))
            cpatch%thbark (ico) = size2xb(cpatch%dbh(ico),cpatch%height(ico)               &
                                         ,cpatch%bbarka(ico),cpatch%bbarkb(ico)            &
                                         ,cpatch%sla(ico),cpatch%pft(ico))
            !------------------------------------------------------------------------------!

            !----- Initialize other cohort-level variables. -------------------------------!
            call init_ed_cohort_vars(cpatch,ico,lsl,mzg,ntext_soil)
            !------------------------------------------------------------------------------!

            !----- Update total patch-level above-ground biomass --------------------------!
            csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                      &
                                        + cpatch%nplant(ico) * cpatch%agb(ico)
            !------------------------------------------------------------------------------!
         end do layerloop
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Because initial heights may not be constant, we must sort the cohorts.      !
         !---------------------------------------------------------------------------------!
         call sort_cohorts(cpatch)
      end do patchloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine init_cohorts_by_layers
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine initializes a near-bare ground 'big-leaf' ed run.                !
   !---------------------------------------------------------------------------------------!
   subroutine near_bare_ground_big_leaf_init(cgrid)
      use ed_state_vars      , only : edtype              & ! structure
                                    , polygontype         & ! structure
                                    , sitetype            & ! structure
                                    , patchtype           & ! structure
                                    , allocate_sitetype   & ! subroutine
                                    , allocate_patchtype  ! ! subroutine
      use ed_misc_coms       , only : ied_init_mode       ! ! intent(in)
      use physiology_coms    , only : n_plant_lim         ! ! intent(in)
      use grid_coms          , only : nzg                 ! ! intent(in)
      use pft_coms           , only : q                   & ! intent(in)
                                    , qsw                 & ! intent(in)
                                    , qbark               & ! intent(in)
                                    , sla                 & ! intent(in)
                                    , agf_bs              & ! intent(in)
                                    , dbh_bigleaf         & ! intent(in)
                                    , hgt_max             & ! intent(in)
                                    , include_pft         & ! intent(in)
                                    , include_these_pft   & ! intent(in)
                                    , f_bstorage_init     & ! intent(in)
                                    , init_density        ! ! intent(in)
      use consts_coms        , only : t3ple               & ! intent(in)
                                    , pio4                & ! intent(in)
                                    , almost_zero         & ! intent(in)
                                    , kgom2_2_tonoha      & ! intent(in)
                                    , tonoha_2_kgom2      ! ! intent(in)
      use physiology_coms    , only : iddmort_scheme      ! ! intent(in)
      use allometry          , only : size2bd             & ! function
                                    , size2bl             & ! function
                                    , size2bt             & ! function
                                    , size2xb             & ! function
                                    , ed_balive           & ! function
                                    , ed_biomass          & ! function
                                    , area_indices        ! ! subroutine
      use ed_type_init       , only : init_ed_cohort_vars & ! subroutine
                                    , init_ed_patch_vars  & ! subroutine
                                    , init_ed_site_vars   & ! subroutine
                                    , init_ed_poly_vars   ! ! subroutine
      use decomp_coms        , only : decomp_scheme       & ! intent(in)
                                    , agf_fsc             & ! intent(in)
                                    , agf_stsc            & ! intent(in)
                                    , c2n_fast_0          & ! intent(in)
                                    , c2n_structural      & ! intent(in)
                                    , c2n_slow            & ! intent(in)
                                    , f0_msc              & ! intent(in)
                                    , f0_ssc              & ! intent(in)
                                    , f0_psc              & ! intent(in)
                                    , nbg_nlim_fsc        & ! intent(in)
                                    , nbg_nlim_stsc       & ! intent(in)
                                    , nbg_nlim_ssc        ! ! intent(in)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)      , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype) , pointer :: cpoly             ! Current polygon
      type(sitetype)    , pointer :: csite             ! Current site
      type(patchtype)   , pointer :: cpatch            ! Current patch
      real                        :: bdeadx            ! Total heartwood biomass
      integer                     :: ipy               ! Patch counter
      integer                     :: ico               ! Cohort counter
      integer                     :: isi               ! Site counter
      integer                     :: ipa               ! Patch counter
      integer                     :: mypfts            ! Number of included PFTs
      integer                     :: ipft              ! PFT counter
      !------------------------------------------------------------------------------------!


      !----- Big loop ---------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !-- Figure out how many patches (1 patch per pft), always primary vegetation. -!
            select case (ied_init_mode)
            case (-1) !------ True bare ground simulation (absolute desert). --------------!
               mypfts = 1
            case ( 0) !------ Nearly bare ground simulation (start with a few seedlings). -!
               mypfts = count(include_pft)
            end select
            csite%npatches = mypfts

            call allocate_sitetype(csite,mypfts)

            !----- Patch loop  ------------------------------------------------------------!
            patchloop: do ipa=1, csite%npatches
               cpatch => csite%patch(ipa)
               ipft = include_these_pft(ipa)

               csite%dist_type          (ipa) = 3
               csite%age                (ipa) = 0.0
               csite%area               (ipa) = 1.0 / mypfts
               csite%fbeam              (ipa) = 1.0
               csite%light_type         (ipa) = 1


               !---------------------------------------------------------------------------!
               !     Someone that uses the nitrogen model should check whether this is     !
               ! necessary or not.  If nitrogen limitation is off, then we start all       !
               ! carbon and nitrogen pools with zeroes, otherwise we initialise with the   !
               ! former default values.  Values have been updated so C and N pools         !
               ! necessarily start following stoichiometry.                                !
               !---------------------------------------------------------------------------!

               select case (n_plant_lim)
               case (0)
                  csite%fast_grnd_C       (ipa) = 0.0
                  csite%fast_soil_C       (ipa) = 0.0
                  csite%structural_grnd_C (ipa) = 0.0
                  csite%structural_soil_C (ipa) = 0.0
                  csite%structural_grnd_L (ipa) = 0.0
                  csite%structural_soil_L (ipa) = 0.0
                  csite%microbial_soil_C  (ipa) = 0.0
                  csite%slow_soil_C       (ipa) = 0.0
                  csite%passive_soil_C    (ipa) = 0.0
                  csite%fast_grnd_N       (ipa) = 0.0
                  csite%fast_soil_N       (ipa) = 0.0
                  csite%structural_grnd_N (ipa) = 0.0
                  csite%structural_soil_N (ipa) = 0.0
                  csite%mineralized_soil_N(ipa) = 0.0
               case (1)
                  csite%fast_grnd_C        (ipa) =        agf_fsc   * nbg_nlim_fsc
                  csite%fast_soil_C        (ipa) = (1.0 - agf_fsc)  * nbg_nlim_fsc
                  csite%structural_grnd_C  (ipa) =        agf_stsc  * nbg_nlim_stsc
                  csite%structural_soil_C  (ipa) = (1.0 - agf_stsc) * nbg_nlim_stsc
                  csite%structural_grnd_L  (ipa) = csite%structural_grnd_C(ipa)
                  csite%structural_soil_L  (ipa) = csite%structural_soil_C(ipa)
                  select case (decomp_scheme)
                  case (5)
                     csite%microbial_soil_C(ipa) = f0_msc * nbg_nlim_ssc
                     csite%slow_soil_C     (ipa) = f0_ssc * nbg_nlim_ssc
                     csite%passive_soil_C  (ipa) = f0_psc * nbg_nlim_ssc
                  case default
                     csite%microbial_soil_C(ipa) = 0.0
                     csite%slow_soil_C     (ipa) = nbg_nlim_ssc
                     csite%passive_soil_C  (ipa) = 0.0
                  end select
                  csite%fast_grnd_N       (ipa) = csite%fast_grnd_C       (ipa)            &
                                                / c2n_fast_0
                  csite%fast_soil_N       (ipa) = csite%fast_soil_C       (ipa)            &
                                                / c2n_fast_0
                  csite%structural_grnd_N (ipa) = csite%structural_grnd_C (ipa)            &
                                                / c2n_structural
                  csite%structural_soil_N (ipa) = csite%structural_soil_C (ipa)            &
                                                / c2n_structural
                  csite%mineralized_soil_N(ipa) = ( csite%microbial_soil_C(ipa)            &
                                                + csite%slow_soil_C     (ipa)              &
                                                + csite%passive_soil_C  (ipa) )            &
                                                / c2n_slow
               end select
               !---------------------------------------------------------------------------!

               csite%sum_dgd            (ipa) = 0.0
               csite%sum_chd            (ipa) = 0.0
               csite%plant_ag_biomass   (ipa) = 0.

               !---------------------------------------------------------------------------!
               !     We now populate the cohorts with near bare ground condition.  In case !
               ! of a desert, we initialise with an empty patch, otherwise we add one co-  !
               ! hort per patch.                                                           !
               !---------------------------------------------------------------------------!
               select case (ied_init_mode)
               case (-1)
                  call allocate_patchtype(cpatch,0)
               case default
                  ico = 1

                  call allocate_patchtype(cpatch,1)
                  !----- The PFT is the plant functional type. ----------------------------!
                  cpatch%pft(ico)              = ipft

                  !------------------------------------------------------------------------!
                  !     Define the near-bare ground state using the standard minimum       !
                  ! height and minimum plant density.  We assume all NBG PFTs to have      !
                  ! leaves fully flushed, but with no storage biomass.  We then compute    !
                  ! the other biomass quantities using the standard allometry for this     !
                  ! PFT.                                                                   !
                  !------------------------------------------------------------------------!
                  cpatch%nplant(ico)           = init_density(ipft)
                  cpatch%dbh(ico)              = dbh_bigleaf(ipft)
                  cpatch%height(ico)           = hgt_max(ipft)
                  cpatch%phenology_status(ico) = 0
                  bdeadx                       = size2bd(cpatch%dbh(ico)                   &
                                                        ,cpatch%height(ico),ipft)
                  cpatch%bdeada(ico)           =       agf_bs(ipft)  * bdeadx
                  cpatch%bdeadb(ico)           = (1. - agf_bs(ipft)) * bdeadx
                  cpatch%sla(ico)              = sla(ipft)
                  cpatch%bleaf(ico)            = size2bl(cpatch%dbh(ico)                   &
                                                        ,cpatch%height(ico)                &
                                                        ,cpatch%sla(ico),ipft)
                  cpatch%broot(ico)            = q(ipft) * cpatch%bleaf(ico)
                  cpatch%bsapwooda(ico)        = agf_bs(ipft) * cpatch%bleaf(ico)          &
                                               * qsw(ipft) * cpatch%height(ico)
                  cpatch%bsapwoodb(ico)        = (1.-agf_bs(ipft)) * cpatch%bleaf(ico)     &
                                               * qsw(ipft) * cpatch%height(ico)
                  cpatch%bbarka(ico)           = agf_bs(ipft) * cpatch%bleaf(ico)          &
                                               * qbark(ipft) * cpatch%height(ico)
                  cpatch%bbarkb(ico)           = (1.-agf_bs(ipft)) * cpatch%bleaf(ico)     &
                                               * qbark(ipft) * cpatch%height(ico)
                  cpatch%balive(ico)           = ed_balive(cpatch,ico)
                  cpatch%bstorage(ico)         = max(almost_zero,f_bstorage_init(ipft))    &
                                               * cpatch%balive(ico)

                  !------------------------------------------------------------------------!
                  !     Initialise the carbon balance.  For initial conditions, we always  !
                  ! assume storage biomass for the previous months so the scale is correct !
                  ! (carbon balance is given in kgC/pl).  The current month carbon balance !
                  ! must be initialised consistently with the iddmort_scheme set by the    !
                  ! user.                                                                  !
                  !------------------------------------------------------------------------!
                  cpatch%cb         (1:12,ico) = cpatch%bstorage(ico)
                  cpatch%cb_lightmax(1:12,ico) = cpatch%bstorage(ico)
                  cpatch%cb_moistmax(1:12,ico) = cpatch%bstorage(ico)
                  cpatch%cb_mlmax   (1:12,ico) = cpatch%bstorage(ico)
                  select case (iddmort_scheme)
                  case (0)
                     !------ Storage is not accounted. ------------------------------------!
                     cpatch%cb         (13,ico) = 0.0
                     cpatch%cb_lightmax(13,ico) = 0.0
                     cpatch%cb_moistmax(13,ico) = 0.0
                     cpatch%cb_mlmax   (13,ico) = 0.0
                     !---------------------------------------------------------------------!
                  case (1)
                     cpatch%cb         (13,ico) = cpatch%bstorage(ico)
                     cpatch%cb_lightmax(13,ico) = cpatch%bstorage(ico)
                     cpatch%cb_moistmax(13,ico) = cpatch%bstorage(ico)
                     cpatch%cb_mlmax   (13,ico) = cpatch%bstorage(ico)
                  end select
                  cpatch%cbr_bar          (ico) = 1.0
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !    Find derived properties: leaf and wood area indices, above-ground   !
                  ! biomass, basal area, timber stocks and bark thickness.                 !
                  !------------------------------------------------------------------------!
                  call area_indices(cpatch, ico)
                  cpatch%agb    (ico) = ed_biomass(cpatch, ico)
                  cpatch%basarea(ico) = pio4 * cpatch%dbh(ico)*cpatch%dbh(ico)
                  cpatch%btimber(ico) = size2bt(cpatch%dbh(ico),cpatch%height(ico)         &
                                               ,cpatch%bdeada(ico),cpatch%bsapwooda(ico)   &
                                               ,cpatch%bbarka(ico),cpatch%pft(ico))
                  cpatch%thbark (ico) = size2xb(cpatch%dbh(ico),cpatch%height(ico)         &
                                               ,cpatch%bbarka(ico),cpatch%bbarkb(ico)      &
                                               ,cpatch%sla(ico),cpatch%pft(ico))
                  !------------------------------------------------------------------------!


                  !----- Initialize other cohort-level variables. -------------------------!
                  call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi),nzg                   &
                                          ,cpoly%ntext_soil(:,isi))
                  !------------------------------------------------------------------------!


                  !----- Update total patch-level above-ground biomass --------------------!
                  csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                &
                                              + cpatch%nplant(ico) * cpatch%agb(ico)
                  !------------------------------------------------------------------------!
               end select
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!

            !----- Initialise the patches now that cohorts are there. ---------------------!
            call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!



         !----- Initialise some site-level variables. -------------------------------------!
         call init_ed_site_vars(cpoly)
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      !----- Last, but not the least, the polygons. ---------------------------------------!
      call init_ed_poly_vars(cgrid)
      !------------------------------------------------------------------------------------!

      return
   end subroutine near_bare_ground_big_leaf_init
   !=======================================================================================!
   !=======================================================================================!
end module ed_nbg_init
!==========================================================================================!
!==========================================================================================!
