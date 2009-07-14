!==========================================================================================!
!==========================================================================================!
!      This subroutine initializes a near-bare ground polygon.                             !
!------------------------------------------------------------------------------------------!
subroutine bare_ground_init(cgrid)
   use ed_state_vars , only: edtype           & ! structure
                           , polygontype      & ! structure
                           , sitetype         & ! structure
                           ,allocate_sitetype ! ! subroutine

   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   integer                     :: ipy
   integer                     :: isi
   !---------------------------------------------------------------------------------------!


   !----- Big loop ------------------------------------------------------------------------!
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         !----- We start with a single patch per site, always with primary vegetation. ----!
         csite%npatches = 1
         call allocate_sitetype(csite,1)
         csite%dist_type          (1) = 3
         csite%age                (1) = 0.0
         csite%area               (1) = 1.0
         csite%fast_soil_C        (1) = 0.2
         csite%slow_soil_C        (1) = 0.01
         csite%structural_soil_C  (1) = 10.0
         csite%structural_soil_L  (1) = csite%structural_soil_C (1)
         csite%mineralized_soil_N (1) = 1.0
         csite%fast_soil_N        (1) = 1.0
         csite%sum_dgd            (1) = 0.0
         csite%sum_chd            (1) = 0.0
         csite%plantation         (1) = 0
         csite%plant_ag_biomass   (1) = 0.

         !----- We now populate the cohorts with near bare ground condition. --------------!
         call init_nbg_cohorts(csite,cpoly%lsl(isi),cpoly%met(isi)%atm_tmp,1,csite%npatches)

         !----- Initialise the patches now that cohorts are there. ------------------------!
         call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
      end do
      !----- Once all patches are set, then we can assign initial values for sites. -------!
      call init_ed_site_vars(cpoly,cgrid%lat(ipy))
   end do
   
   !----- Last, but not the least, the polygons. ------------------------------------------!
   call init_ed_poly_vars(cgrid)

   return
end subroutine bare_ground_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine assigns a near-bare ground (NBG) state for some patches.            !
!------------------------------------------------------------------------------------------!
subroutine init_nbg_cohorts(csite,lsl,atm_tmp,ipa_a,ipa_z)
   use ed_state_vars      , only : edtype             & ! structure
                                 , polygontype        & ! structure
                                 , sitetype           & ! structure
                                 , patchtype          & ! structure
                                 , allocate_sitetype  & ! subroutine
                                 , allocate_patchtype ! ! subroutine
   use max_dims           , only : n_pft              ! ! intent(in)
   use pft_coms           , only : q                  & ! intent(in)
                                 , qsw                & ! intent(in)
                                 , sla                & ! intent(in)
                                 , hgt_min            & ! intent(in)
                                 , include_pft        & ! intent(in)
                                 , include_these_pft  & ! intent(in)
                                 , include_pft_ag     & ! intent(in)
                                 , init_density       ! ! intent(in)
   use consts_coms        , only : t3ple              ! ! intent(in)
   use ed_therm_lib       , only : calc_hcapveg       ! ! function
   use allometry          , only : h2dbh              & ! function
                                 , dbh2bd             & ! function
                                 , dbh2bl             & ! function
                                 , ed_biomass         & ! function
                                 , area_indices       ! ! subroutine
   use fuse_fiss_utils , only : sort_cohorts    ! ! subroutine

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target     :: csite   ! Current site
   integer               , intent(in) :: lsl     ! Lowest soil level
   integer               , intent(in) :: ipa_a   ! 1st patch to be assigned with NBG state
   integer               , intent(in) :: ipa_z   ! Last patch to be assigned with NBG state
   real                  , intent(in) :: atm_tmp ! Atmospheric temperature
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch  ! Current patch
   integer                            :: ipa     ! Patch number
   integer                            :: ico     ! Cohort counter
   integer                            :: mypfts  ! Number of PFTs to be included.
   integer                            :: ipft    ! PFT counter
   !---------------------------------------------------------------------------------------!


   !----- Patch loop. ---------------------------------------------------------------------!
   patchloop: do ipa=ipa_a,ipa_z
      cpatch => csite%patch(ipa)

      ! Decide how many cohorts to allocate
      select case (csite%dist_type(ipa))
      case (1)   ! Agriculture
         mypfts = sum(include_pft_ag)
      case (2,3) ! Secondary or primary forest
         mypfts= sum(include_pft)
      end select

      !----- Perform cohort allocation. ---------------------------------------------------!
      call allocate_patchtype(cpatch,mypfts)
      
      !------------------------------------------------------------------------------------!
      !    Here we loop over PFTs rather than assigning the cohorts, so we ensure to only  !
      ! include the PFTs that should be  included (i.e., patches that should not happen in !
      ! agricultural patches will be skipped, and agricultural PFTs will not be included   !
      ! in the forest sites.                                                               !
      !------------------------------------------------------------------------------------!
      ico = 0
      pftloop: do ipft = 1,n_pft
         select case (csite%dist_type(ipa))
         case (1)
            !----- Agriculture, only the ones allowed are included. -----------------------!
            if (include_pft_ag(ipft) == 1) then
               ico = ico + 1
            else
               cycle pftloop
            end if
         case (2,3)
            !----- Forest, only the ones allowed are included. ----------------------------!
            if (include_pft(ipft) == 1) then
               ico = ico + 1
            else
               cycle pftloop
            end if
         end select

         !----- The PFT is the plant functional type. -------------------------------------!
         cpatch%pft(ico)              = ipft
   
         !---------------------------------------------------------------------------------!
         !     Define the near-bare ground state using the standard minimum height and     !
         ! minimum plant density.  We assume all NBG PFTs to have leaves fully flushed,    !
         ! but with no storage biomass.  We then compute the other biomass quantities      !
         ! using the standard allometry for this PFT.                                      !
         !---------------------------------------------------------------------------------!
         cpatch%nplant(ico)           = init_density(ipft)
         cpatch%hite(ico)             = hgt_min(ipft)
         cpatch%phenology_status(ico) = 0
         cpatch%bstorage(ico)         = 0.0
         cpatch%dbh(ico)              = h2dbh(cpatch%hite(ico),ipft)
         cpatch%bdead(ico)            = dbh2bd(cpatch%dbh(ico),cpatch%hite(ico),ipft)
         cpatch%bleaf(ico)            = dbh2bl(cpatch%dbh(ico),ipft)
         cpatch%sla(ico)              = sla(ipft)
         cpatch%balive(ico)           = cpatch%bleaf(ico)                                  &
                                      * ( 1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico) )

         !----- Find the initial area indices (LAI, WPA, WAI). ----------------------------!
         call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bdead(ico)          &
                          ,cpatch%balive(ico),cpatch%dbh(ico), cpatch%hite(ico)            &
                          ,cpatch%pft(ico),cpatch%sla(ico),cpatch%lai(ico)                 &
                          ,cpatch%wpa(ico),cpatch%wai(ico))

         
         !----- Initialize other cohort-level variables. ----------------------------------!
         call init_ed_cohort_vars(cpatch,ico,lsl)
         
         !---------------------------------------------------------------------------------!
         !     Set the initial vegetation thermodynamic properties.  We assume the veget-  !
         ! ation to be with no condensed/frozen water in their surfaces, and the temper-   !
         ! ature to be the same as the canopy air space.  Then we find the internal energy !
         ! and heat capacity.                                                              !
         !---------------------------------------------------------------------------------!
         cpatch%veg_water(ico)  = 0.0
         cpatch%veg_fliq(ico)   = 0.0
         cpatch%veg_temp(ico)   = atm_tmp
         cpatch%hcapveg(ico)    = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico)         &
                                              ,cpatch%balive(ico),cpatch%nplant(ico)       &
                                              ,cpatch%hite(ico),cpatch%pft(ico)            &
                                              ,cpatch%phenology_status(ico))
         cpatch%veg_energy(ico) = cpatch%hcapveg(ipa) * cpatch%veg_temp(ico)

         !----- Update total patch-level above-ground biomass -----------------------------!
         csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa) + cpatch%nplant(ico)    &
                                     * ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)     &
                                                 ,cpatch%bleaf(ico),cpatch%pft(ico)        &
                                                 ,cpatch%hite(ico),cpatch%bstorage(ico))
      end do pftloop
      
      !------------------------------------------------------------------------------------!
      !     Since initial heights may not be constant, we must sort the cohorts.           !
      !------------------------------------------------------------------------------------!
      call sort_cohorts(cpatch)
   end do patchloop

   return
end subroutine init_nbg_cohorts
!==========================================================================================!
!==========================================================================================!
