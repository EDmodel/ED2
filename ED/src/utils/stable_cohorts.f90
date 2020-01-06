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
                  call is_resolvable(csite,ipa,ico,.false.,.false.,'flag_stable_cohorts')
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
   !                                                                                       !
   !    The one exception is during cohort fusion (except during initialisation).  In this !
   ! case we force cohorts to be resolvable to ensure water and energy changes are         !
   ! accounted                                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine is_resolvable(csite,ipa,ico,is_initial,force_resolvable,called_from)
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
      use ed_misc_coms   , only : frqsumi                & ! intent(in)
                                , frqsum                 & ! intent(in)
                                , current_time           ! ! intent(in)
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
      type(sitetype)  , target     :: csite              ! Current site
      integer         , intent(in) :: ipa                ! Patch index
      integer         , intent(in) :: ico                ! Cohort index
      logical         , intent(in) :: is_initial         ! Initial assignment?
      logical         , intent(in) :: force_resolvable   ! Impose resolvable?
      character(len=*), intent(in) :: called_from        ! Parent subroutine
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch             ! Current patch
      integer                      :: ipft               ! Cohort PFT
      logical                      :: exposed            ! Cohort is above snow     [   T|F]
      logical                      :: leaf_enough        ! Cohort has enough leaves [   T|F]
      logical                      :: wood_enough        ! Cohort has enough wood   [   T|F]
      logical                      :: hydro_req          ! Plant hydr. requirement  [   T|F]
      logical                      :: old_leaf_res       ! Leaf was resolvable      [   T|F]
      logical                      :: old_wood_res       ! Wood was resolvable      [   T|F]
      logical                      :: update_hydro       ! Update hydraulic props   [   T|F]
      logical                      :: is_pheneffect      ! Phenology effect updated [   T|F]
      logical                      :: switch_leaf_on     ! Leaf became resolvable   [   T|F]
      logical                      :: switch_wood_on     ! Wood became resolvable   [   T|F]
      logical                      :: switch_leaf_off    ! Leaf became unresolvable [   T|F]
      logical                      :: switch_wood_off    ! Wood became unresolvable [   T|F]
      real                         :: mid_leaf_hcap      ! Mid leaf heat capacity   [J/kg/K]
      real                         :: mid_wood_hcap      ! Mid wood heat capacity   [J/kg/K]
      real                         :: mid_leaf_water     ! Mid leaf surface water   [ kg/m2]
      real                         :: mid_wood_water     ! Mid wood surface water   [ kg/m2]
      real                         :: mid_leaf_water_im2 ! Mid leaf internal water  [ kg/m2]
      real                         :: mid_wood_water_im2 ! Mid wood internal water  [ kg/m2]
      real                         :: c_leaf             ! Leaf water capacity      [  kg/m]
      real                         :: c_stem             ! Stem water capacity      [  kg/m]
      real                         :: sap_frac           ! Sapwood fraction         [   ---]
      real                         :: bleafhydro         ! Hydro leaf biomass       [kgC/pl]
      real                         :: broothydro         ! Hydro fine-root biomass  [kgC/pl]
      real                         :: bsapwhydro         ! Hydro sapwood biomass    [kgC/pl]
      real                         :: bdeadhydro         ! Hydro heartwood biomass  [kgC/pl]
      !----- Local constants. -------------------------------------------------------------!
      logical          , parameter :: is_debug    = .false.
      character(len=18), parameter :: stable_file = 'stable_details.txt'
      character(len= 8), parameter :: fmtc='(a,1x,a)'
      character(len=10), parameter :: fmti='(a,1x,i14)'
      character(len=10), parameter :: fmtl='(a,14x,l1)'
      character(len=12), parameter :: fmtf='(a,1x,f14.3)'
      character(len=13), parameter :: fmte='(a,1x,es14.7)'
      character(len=27), parameter :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !----- Saved variables. -------------------------------------------------------------!
      logical          , save      :: first_time = .true.
      !------------------------------------------------------------------------------------!



      !------ Handy aliases. --------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      ipft   =  cpatch%pft(ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we are debugging, create the debugging file.                           !
      !------------------------------------------------------------------------------------!
      if (is_debug .and. first_time) then
         first_time = .false.
         open (unit=53,file=stable_file,status='replace',action='write')
         write(unit=53,fmt='(a)') '-------------------------------------------------------'
         write(unit=53,fmt='(a)') ' Subroutine is_resolvable'
         write(unit=53,fmt='(a)') '-------------------------------------------------------'
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt='(a)') ' '
         close(unit=53,status='keep')
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Save original "resolvable" condition.  These are used to decide whether to     !
      ! update heat capacity, internal water, and the budget terms.                        !
      !------------------------------------------------------------------------------------!
      old_leaf_res       = cpatch%leaf_resolvable(ico)
      old_wood_res       = cpatch%wood_resolvable(ico)
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
      if (force_resolvable) then
         cpatch%leaf_resolvable(ico) = .true.
         cpatch%wood_resolvable(ico) = .true.
      else
         cpatch%leaf_resolvable(ico) = exposed .and. leaf_enough
         cpatch%wood_resolvable(ico) = exposed .and. ( wood_enough .or. hydro_req )
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Save logical tests to decide whether there was any switch in resolvable        !
      ! statuses.                                                                          !
      !------------------------------------------------------------------------------------!
      switch_leaf_on  = cpatch%leaf_resolvable(ico) .and. (.not. old_leaf_res)
      switch_wood_on  = cpatch%wood_resolvable(ico) .and. (.not. old_wood_res)
      switch_leaf_off = (.not. cpatch%leaf_resolvable(ico)) .and. old_leaf_res
      switch_wood_off = (.not. cpatch%wood_resolvable(ico)) .and. old_wood_res
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Update phenology effect in the case the leaf and/or wood status have changed.  !
      !------------------------------------------------------------------------------------!
      is_pheneffect = .false.
      if (.not. is_initial) then
         !----- Leaves became resolvable. -------------------------------------------------!
         if (switch_leaf_on) then
            csite%wbudget_pheneffect(ipa) = csite%wbudget_pheneffect(ipa)                  &
                                          + frqsumi * ( cpatch%leaf_water    (ico)         &
                                                      + cpatch%leaf_water_im2(ico) )
            csite%ebudget_pheneffect(ipa) = csite%ebudget_pheneffect(ipa)                  &
                                          + frqsumi * cpatch%leaf_energy(ico)
            is_pheneffect = .true.
         end if
         !----- Wood became resolvable. ---------------------------------------------------!
         if (switch_wood_on) then
            csite%wbudget_pheneffect(ipa) = csite%wbudget_pheneffect(ipa)                  &
                                          + frqsumi * ( cpatch%wood_water    (ico)         &
                                                      + cpatch%wood_water_im2(ico) )
            csite%ebudget_pheneffect(ipa) = csite%ebudget_pheneffect(ipa)                  &
                                          + frqsumi * cpatch%wood_energy(ico)
            is_pheneffect = .true.
         end if
         !----- Leaves are no longer resolvable. ------------------------------------------!
         if (switch_leaf_off) then
            csite%wbudget_pheneffect(ipa) = csite%wbudget_pheneffect(ipa)                  &
                                          - frqsumi * ( cpatch%leaf_water    (ico)         &
                                                      + cpatch%leaf_water_im2(ico) )
            csite%ebudget_pheneffect(ipa) = csite%ebudget_pheneffect(ipa)                  &
                                          - frqsumi * cpatch%leaf_energy(ico)
            is_pheneffect = .true.
         end if
         !----- Wood is no longer resolvable. ---------------------------------------------!
         if (switch_wood_off) then
            csite%wbudget_pheneffect(ipa) = csite%wbudget_pheneffect(ipa)                  &
                                          - frqsumi * ( cpatch%wood_water    (ico)         &
                                                      + cpatch%wood_water_im2(ico) )
            csite%ebudget_pheneffect(ipa) = csite%ebudget_pheneffect(ipa)                  &
                                          - frqsumi * cpatch%wood_energy(ico)
            is_pheneffect = .true.
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Before we move on, we save the heat capacity and water (surface and internal). !
      ! In case we must update internal water, we will need these to adjust the budget     !
      ! variables.                                                                         !
      !------------------------------------------------------------------------------------!
      mid_leaf_hcap      = cpatch%leaf_hcap      (ico)
      mid_wood_hcap      = cpatch%wood_hcap      (ico)
      mid_leaf_water     = cpatch%leaf_water     (ico)
      mid_wood_water     = cpatch%wood_water     (ico)
      mid_leaf_water_im2 = cpatch%leaf_water_im2 (ico)
      mid_wood_water_im2 = cpatch%wood_water_im2 (ico)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide whether to update hydrology.                                            !
      !------------------------------------------------------------------------------------!
      update_hydro = ( switch_leaf_on .or. switch_wood_on ) .and. (.not. force_resolvable)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      If plants are emerging from non-resolvable status, we must ensure that        !
      ! internal water is bounded.                                                         !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Impose saturated conditions. ----------------------------------------------!
         if (switch_wood_on) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%wood_rwc(ico) = 1.0
            !------------------------------------------------------------------------------!
         end if
         if (switch_leaf_on) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%leaf_rwc(ico) = 1.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      case default
         if (switch_wood_on) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%wood_rwc(ico) = max(wood_rwc_min(ipft),min(1.0,cpatch%wood_rwc(ico)))
            !------------------------------------------------------------------------------!
         end if
         if (switch_wood_on) then
            !----- Make sure relative water content is bounded. ---------------------------!
            cpatch%leaf_rwc(ico) = max(leaf_rwc_min(ipft),min(1.0,cpatch%leaf_rwc(ico)))
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
         call update_veg_energy_cweh(csite,ipa,ico,mid_leaf_hcap,mid_wood_hcap             &
                                    ,mid_leaf_water,mid_wood_water,mid_leaf_water_im2      &
                                    ,mid_wood_water_im2,.true.)
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




      !------------------------------------------------------------------------------------!
      !     Print debugging information.                                                   !
      !------------------------------------------------------------------------------------!
      if (is_pheneffect .and. is_debug) then
         open (unit=53,file=stable_file,status='old',position='append',action='write')
         write(unit=53,fmt='(a)') '--------------------------------------------------------'
         write(unit=53,fmt=fmtt )  ' TIME       = ', current_time%year                     &
                                                   , current_time%month                    &
                                                   , current_time%date                     &
                                                   , current_time%time
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt=fmtc ) ' CALLED_FROM = ', called_from
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt=fmti ) ' PATCH       = ', ipa
         write(unit=53,fmt=fmti ) ' COHORT      = ', ico
         write(unit=53,fmt=fmti ) ' PFT         = ', cpatch%pft(ico)
         write(unit=53,fmt=fmti ) ' PHEN_STATUS = ', cpatch%phenology_status(ico)
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt=fmtf ) ' DBH         = ', cpatch%dbh   (ico)
         write(unit=53,fmt=fmtf ) ' HITE        = ', cpatch%hite  (ico)
         write(unit=53,fmt=fmtf ) ' NPLANT_HA   = ', cpatch%nplant(ico) * 10000.
         write(unit=53,fmt=fmtf ) ' LAI         = ', cpatch%lai   (ico)
         write(unit=53,fmt=fmtf ) ' WAI         = ', cpatch%wai   (ico)
         write(unit=53,fmt=fmtf ) ' ELONGF      = ', cpatch%elongf(ico)
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt=fmtl ) ' FORCE_RES   = ', force_resolvable
         write(unit=53,fmt=fmtl ) ' LEAF_ON     = ', switch_leaf_on
         write(unit=53,fmt=fmtl ) ' LEAF_OFF    = ', switch_leaf_off
         write(unit=53,fmt=fmtl ) ' WOOD_ON     = ', switch_wood_on
         write(unit=53,fmt=fmtl ) ' WOOD_OFF    = ', switch_wood_off
         write(unit=53,fmt=fmtl ) ' UP_HYDRO    = ', update_hydro
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt=fmte ) ' EB_PHEN     = ', csite%ebudget_pheneffect(ipa) * frqsum
         write(unit=53,fmt=fmte ) ' EB_WCAP     = ', csite%ebudget_wcapeffect(ipa) * frqsum
         write(unit=53,fmt=fmte ) ' EB_HCAP     = ', csite%ebudget_hcapeffect(ipa) * frqsum
         write(unit=53,fmt=fmte ) ' WB_PHEN     = ', csite%wbudget_pheneffect(ipa) * frqsum
         write(unit=53,fmt=fmte ) ' WB_WCAP     = ', csite%wbudget_wcapeffect(ipa) * frqsum
         write(unit=53,fmt='(a)') '========================================================'
         write(unit=53,fmt='(a)') ' '
         write(unit=53,fmt='(a)') ' '
         close(unit=53,status='keep')
      end if
      !------------------------------------------------------------------------------------!
      return
   end subroutine is_resolvable
   !=======================================================================================!
   !=======================================================================================!
end module stable_cohorts
!==========================================================================================!
!==========================================================================================!
