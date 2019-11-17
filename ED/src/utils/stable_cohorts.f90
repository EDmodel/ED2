!==========================================================================================!
!==========================================================================================!
!  Module stable_cohorts
!
!> \brief   Sub-routine and functions that flag whether cohorts are resolvable or not
!> \details These routines and functions define two logical flags (leaf_resolvable and 
!>          wood_resolvable) that will decide whether cohorts can be solved in the water,
!>          energy, and CO2 budget.  These flags ensure that the cohort is consistently
!>          skipped or resolved.
!------------------------------------------------------------------------------------------!
module stable_cohorts
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will go through all cohorts in a polygon, and assign them as      !
   ! either numerically safe or unsafe.  The idea of to use a unique test to decide which  !
   ! cohorts we will solve in photosynthesis, radiation, and the energy balance (RK4 or    !
   ! Euler).                                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine flag_stable_cohorts(cgrid)
      use ed_state_vars         , only : edtype          & ! structure
                                       , polygontype     & ! structure
                                       , sitetype        & ! structure
                                       , patchtype       ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target   :: cgrid         ! Current grid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer  :: cpoly         ! Current polygon
      type(sitetype)   , pointer  :: csite         ! Current site
      type(patchtype)  , pointer  :: cpatch        ! Current patch
      integer                     :: ipy           ! Polygon index
      integer                     :: isi           ! Site index
      integer                     :: ipa           ! Patch index
      integer                     :: ico           ! Cohort index
      !------------------------------------------------------------------------------------!

      polyloop: do ipy=1, cgrid%npolygons

         cpoly => cgrid%polygon(ipy)
         siteloop: do isi=1,cpoly%nsites

            csite => cpoly%site(isi)
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               cohortloop: do ico=1, cpatch%ncohorts

                  !----- Check whether we can resolve this cohort. ------------------------!
                  call is_resolvable(csite,ipa,ico)
                  !------------------------------------------------------------------------!

               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
         end do siteloop
      end do polyloop

      return
   end subroutine flag_stable_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks whether cohort leaves and wood thermodynamics can or can- !
   ! not be solved.  This was taken out of flag_stable_cohorts so we can call it in other  !
   ! places and change the cohort status as soon as the conditions change.  We will skip   !
   ! the cohort if any of these conditions happens:                                        !
   ! 1. The cohort leaf biomass is much less than what it would have when all              !
   !    leaves are fully flushed (leaves only);                                            !
   ! 2. The cohort is buried in snow or under water                                        !
   ! 3. The cohort is extremely sparse.                                                    !
   ! 4. The user doesn't want to solve wood thermodynamics (wood only).                    !
   !---------------------------------------------------------------------------------------!
   subroutine is_resolvable(csite,ipa,ico)
      use ed_state_vars  , only : sitetype               & ! structure
                                , patchtype              ! ! structure
      use pft_coms       , only : C2B                    & ! intent(in)
                                , q                      & ! intent(in)
                                , qsw                    & ! intent(in)
                                , is_grass               & ! intent(in)
                                , veg_hcap_min           & ! intent(in)
                                , leaf_rwc_min           & ! intent(in)
                                , wood_rwc_min           & ! intent(in)
                                , leaf_water_cap         & ! intent(in)
                                , wood_water_cap         & ! intent(in)
                                , hgt_min                ! ! intent(in)
      use ed_max_dims    , only : n_pft                  ! ! intent(in)
      use physiology_coms, only : istomata_scheme        & ! function
                                , plant_hydro_scheme     ! ! function
      use allometry      , only : dbh2sf                 & ! function
                                , size2bl                & ! function
                                , size2bd                ! ! function
      use plant_hydro    , only : rwc2psi                & ! subroutine
                                , rwc2tw                 & ! subroutine
                                , twi2twe                ! ! subroutine
      use ed_therm_lib   , only : update_veg_energy_cweh ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target     :: csite      ! Current site
      integer               , intent(in) :: ipa        ! Patch index
      integer               , intent(in) :: ico        ! Cohort index
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch             ! Current patch
      integer                     :: ipft               ! Cohort PFT
      logical                     :: exposed            ! Cohort is above snow     [   T|F]
      logical                     :: leaf_enough        ! Cohort has enough leaves [   T|F]
      logical                     :: wood_enough        ! Cohort has enough wood   [   T|F]
      logical                     :: hydro_req          ! Plant hydr. requirement  [   T|F]
      logical                     :: old_leaf_res       ! Leaf was resolvable      [   T|F]
      logical                     :: old_wood_res       ! Wood was resolvable      [   T|F]
      logical                     :: update_hydro       ! Update hydraulic props   [   T|F]
      real                        :: old_leaf_hcap      ! Old leaf heat capacity   [J/kg/K]
      real                        :: old_wood_hcap      ! Old wood heat capacity   [J/kg/K]
      real                        :: old_leaf_water_im2 ! Old leaf internal water  [ kg/m2]
      real                        :: old_wood_water_im2 ! Old wood internal water  [ kg/m2]
      real                        :: c_leaf             ! Leaf water capacity      [  kg/m]
      real                        :: c_stem             ! Stem water capacity      [  kg/m]
      real                        :: sap_frac           ! Sapwood fraction         [   ---]
      real                        :: bleafhydro         ! Hydro leaf biomass       [kgC/pl]
      real                        :: broothydro         ! Hydro fine-root biomass  [kgC/pl]
      real                        :: bsapwhydro         ! Hydro sapwood biomass    [kgC/pl]
      real                        :: bdeadhydro         ! Hydro heartwood biomass  [kgC/pl]
      !------------------------------------------------------------------------------------!



      !------ Handy aliases. --------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      ipft   =  cpatch%pft(ico)
      !------------------------------------------------------------------------------------!


      !------ Start by assuming that no update in internal water is needed. ---------------!
      update_hydro = .false.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Save previous condition, to avoid singularities for leaves and wood coming     !
      ! out of senescence.                                                                 !
      !------------------------------------------------------------------------------------!
      old_leaf_res       = cpatch%leaf_resolvable(ico)
      old_wood_res       = cpatch%wood_resolvable(ico)
      old_leaf_hcap      = cpatch%leaf_hcap      (ico)
      old_wood_hcap      = cpatch%wood_hcap      (ico)
      old_leaf_water_im2 = cpatch%leaf_water_im2 (ico)
      old_wood_water_im2 = cpatch%wood_water_im2 (ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 1.  Check for cohort height relative to snow/water depth.  If the cohort is buried !
      !     in snow or has drowned in the standing water, we can't solve it.               !
      !------------------------------------------------------------------------------------!
      exposed      = cpatch%hite(ico)  > csite%total_sfcw_depth(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 2.   Check whether this cohort is not extremely sparse.  Wood area heat capacity   !
      !      is always set to zero when branch thermodynamics is turned off, so this will  !
      !      always be .false. in this case.                                               !
      !------------------------------------------------------------------------------------!
      leaf_enough  = cpatch%leaf_hcap(ico) > veg_hcap_min(ipft)
      wood_enough  = cpatch%wood_hcap(ico) > veg_hcap_min(ipft)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 3.   In case dynamic plant hydraulics is active, wood must be resolved whenever    !
      !      leaf is resolved.                                                             !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         hydro_req = .false.
      case default
         hydro_req = leaf_enough
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Save the tests in the cohort variable, so the checks are done consistently.    !
      !------------------------------------------------------------------------------------!
      cpatch%leaf_resolvable(ico) = exposed .and. leaf_enough
      cpatch%wood_resolvable(ico) = exposed .and. ( wood_enough .or. hydro_req )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      If plants are emerging from non-resolvable status, we must ensure that        !
      ! internal water is bounded.                                                         !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Impose saturated conditions. ----------------------------------------------!
         if (cpatch%wood_resolvable(ico) .and. (.not. old_wood_res)) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%wood_rwc(ico) = 1.0
            update_hydro         = .true.
            !------------------------------------------------------------------------------!
         end if
         if (cpatch%leaf_resolvable(ico) .and. (.not. old_leaf_res)) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%leaf_rwc(ico) = 1.0
            update_hydro         = .true.
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      case default
         if (cpatch%wood_resolvable(ico) .and. (.not. old_wood_res)) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%wood_rwc(ico) = max(wood_rwc_min(ipft),min(1.0,cpatch%wood_rwc(ico)))
            update_hydro         = .true.
            !------------------------------------------------------------------------------!
         end if
         if (cpatch%leaf_resolvable(ico) .and. (.not. old_leaf_res)) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%leaf_rwc(ico) = max(leaf_rwc_min(ipft),min(1.0,cpatch%leaf_rwc(ico)))
            update_hydro         = .true.
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      If plants are emerging from non-resolvable status, we must ensure that        !
      ! internal water is bounded.                                                         !
      !------------------------------------------------------------------------------------!
      if (update_hydro) then
         !----- Convert water potential to relative water content. ------------------------!
         call rwc2psi(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico),cpatch%pft(ico)            &
                     ,cpatch%leaf_psi(ico),cpatch%wood_psi(ico))
         !---------------------------------------------------------------------------------!


         !----- Convert relative water content to total water content. --------------------!
         call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                             &
                    ,cpatch%bleaf(ico),cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico)         &
                    ,cpatch%bdeada(ico),cpatch%bdeadb(ico),cpatch%broot(ico)               &
                    ,cpatch%dbh(ico),cpatch%pft(ico),cpatch%leaf_water_int(ico)            &
                    ,cpatch%wood_water_int(ico))
         !---------------------------------------------------------------------------------!


         !----- Convert total water content (kgW/plant) to extensive (kgW/m2). ---------------!
         call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)                &
                     ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)                        &
                     ,cpatch%wood_water_im2(ico))
         !---------------------------------------------------------------------------------!


         !----- Update internal energy. ---------------------------------------------------!
         call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap             &
                                    ,old_leaf_water_im2,old_wood_water_im2,.true.)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Flag cohorts as small or large.  "Small" cohorts are only flagged when plant   !
      ! hydraulics is on and stomatal conductance is calculated using methods based on     !
      ! leaf water potential.                                                              !
      !------------------------------------------------------------------------------------!
      select case (istomata_scheme)
      case (0)
         !---------------------------------------------------------------------------------!
         !     Stomatal conductance is not based on potential.  Ignore small/large cohort  !
         ! flag.                                                                           !
         !---------------------------------------------------------------------------------!
         cpatch%is_small(ico) = .false.
         !---------------------------------------------------------------------------------!
      case default
         select case (plant_hydro_scheme)
         case (0)
            !------------------------------------------------------------------------------!
            !     Plant hydraulics is not dynamic.  No need to distinguish small and large !
            ! cohorts.                                                                     !
            !------------------------------------------------------------------------------!
            cpatch%is_small(ico) = .false.
            !------------------------------------------------------------------------------!
         case default
            !------------------------------------------------------------------------------!
            !      Check the relative magnitude of leaf and sapwood water pool.  In case   !
            ! the cohort is a small tree/liana or a grass, then we flag them so we treat   !
            ! their leaf and wood potential to be the same.  Cohorts are considered too    !
            ! small when c_leaf is greater than half of c_stem, or if the cohort has not   !
            ! grown since recruitment.  These are arbitrary thresholds.  Users are         !
            ! welcome to modify this term in case leaf_psi shows strong oscillations from  !
            ! each timestep to another.                                                    !
            !------------------------------------------------------------------------------!
            bleafhydro = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
            broothydro = q(ipft) * bleafhydro
            sap_frac   = dbh2sf(cpatch%dbh(ico),ipft)
            bsapwhydro = qsw(ipft) * cpatch%hite(ico) * bleafhydro
            bdeadhydro = size2bd(cpatch%dbh(ico),cpatch%hite(ico),ipft)
            bsapwhydro = ( bsapwhydro + bdeadhydro ) * sap_frac
            !----- Find leaf and stem capacities. -----------------------------------------!
            c_leaf               = leaf_water_cap(ipft) * C2B * bleafhydro
            c_stem               = wood_water_cap(ipft) * C2B * (broothydro  + bsapwhydro)
            cpatch%is_small(ico) = is_grass(ipft)                     .or.                 &
                                   c_leaf           >  (0.5 * c_stem) .or.                 &
                                   cpatch%hite(ico) == hgt_min(ipft)
            cpatch%is_small(ico) = .false.
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
      return
   end subroutine is_resolvable
   !=======================================================================================!
   !=======================================================================================!
end module stable_cohorts
!==========================================================================================!
!==========================================================================================!
