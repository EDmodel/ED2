!==========================================================================================!
!==========================================================================================!
!    Module rk4_misc, it contains multiple utility functions that are needed by RK4.       !
!------------------------------------------------------------------------------------------!
module rk4_misc
  contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is called before the sanity check, and updates the diagnostic     !
   ! variables, namely the temperature and liquid fraction of leaf water, soil layers and  !
   ! temporary snow/pond layers.                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine update_diagnostic_vars(initp, csite,ipa,ibuff)
      use rk4_coms              , only : rk4eps                & ! intent(in)
                                       , rk4site               & ! intent(in)
                                       , ibranch_thermo        & ! intent(in)
                                       , rk4tiny_sfcw_mass     & ! intent(in)
                                       , rk4min_sfcw_mass      & ! intent(in)
                                       , rk4min_virt_water     & ! intent(in)
                                       , rk4min_can_shv        & ! intent(in)
                                       , rk4max_can_shv        & ! intent(in)
                                       , rk4min_veg_lwater     & ! intent(in)
                                       , rk4min_veg_temp       & ! intent(in)
                                       , rk4max_veg_temp       & ! intent(in)
                                       , rk4min_soil_temp      & ! intent(in)
                                       , rk4max_soil_temp      & ! intent(in)
                                       , rk4aux                & ! structure
                                       , rk4min_sfcw_temp      & ! intent(in)
                                       , rk4max_sfcw_temp      & ! intent(in)
                                       , tiny_offset           & ! intent(in)
                                       , rk4patchtype          ! ! structure
      use ed_state_vars         , only : sitetype              & ! structure
                                       , patchtype             ! ! structure
      use soil_coms             , only : soil8                 & ! intent(in)
                                       , dslz8                 & ! intent(in)
                                       , dslzi8                & ! intent(in)
                                       , soil_rough8           & ! intent(in)
                                       , ny07_eq04_a8          & ! intent(in)
                                       , ny07_eq04_m8          & ! intent(in)
                                       , matric_potential8     ! ! function
      use grid_coms             , only : nzg                   & ! intent(in)
                                       , nzs                   ! ! intent(in)
      use therm_lib8            , only : uextcm2tl8            & ! subroutine
                                       , uint2tl8              & ! subroutine
                                       , tl2uint8              & ! function
                                       , cmtl2uext8            & ! function
                                       , thetaeiv8             & ! function
                                       , rehuil8               & ! function
                                       , qslif8                & ! function
                                       , hq2temp8              & ! function
                                       , extemp2theta8         & ! function
                                       , thil2tqall8           ! ! function
      use consts_coms           , only : t3ple8                & ! intent(in)
                                       , cpdry8                & ! intent(in)
                                       , cph2o8                & ! intent(in)
                                       , wdns8                 & ! intent(in)
                                       , fsdns8                & ! intent(in)
                                       , fsdnsi8               & ! intent(in)
                                       , rdryi8                & ! intent(in)
                                       , rdry8                 & ! intent(in)
                                       , rmol8                 & ! intent(in)
                                       , ep8                   & ! intent(in)
                                       , epim18                & ! intent(in)
                                       , cph2o_tscvap8         & ! intent(in)
                                       , cpdiff_epim18         & ! intent(in)
                                       , cpor8                 & ! intent(in)
                                       , mmdry8                & ! intent(in)
                                       , mmdryi8               & ! intent(in)
                                       , toodry8               ! ! intent(in)
      use plant_hydro           , only : om_buff_d             & ! intent(in)
                                       , op_buff_d             ! ! intent(in)
      use canopy_struct_dynamics, only : canopy_turbulence8    ! ! subroutine
      use ed_therm_lib          , only : ed_grndvap8           ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: initp
      type(sitetype)     , target     :: csite
      integer            , intent(in) :: ipa
      integer            , intent(in) :: ibuff
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer :: cpatch
      integer                          :: ico
      integer                          :: k
      integer                          :: ksn
      integer                          :: nsoil
      logical                          :: ok_can_enthalpy
      logical                          :: ok_can_theta
      logical                          :: ok_can_shv
      logical                          :: ok_ground
      logical                          :: ok_sfcw
      logical                          :: ok_leaf
      logical                          :: ok_wood
      logical                          :: ok_slwater
      logical                          :: ok_sltemp
      real(kind=8)                     :: int_sfcw_energy
      real(kind=8)                     :: int_virt_energy
      real(kind=8)                     :: energy_tot
      real(kind=8)                     :: wmass_tot
      real(kind=8)                     :: veg_temp
      real(kind=8)                     :: veg_fliq
      real(kind=8)                     :: hcapdry_tot
      real(kind=8)                     :: rk4min_veg_water
      real(kind=8)                     :: rk4min_leaf_water
      real(kind=8)                     :: rk4min_wood_water
      real(kind=8)                     :: rk4min_leaf_water_im2
      real(kind=8)                     :: rk4max_leaf_water_im2
      real(kind=8)                     :: rk4min_wood_water_im2
      real(kind=8)                     :: rk4max_wood_water_im2
      real(kind=8)                     :: wgt_leaf
      real(kind=8)                     :: wgt_wood
      real(kind=8)                     :: bulk_sfcw_dens
      !------------------------------------------------------------------------------------!

      !----- Then we define some logicals to make the code cleaner. -----------------------!
      ok_can_shv      = initp%can_shv      >= rk4min_can_shv                      .and.    &
                        initp%can_shv      <= rk4max_can_shv
      ok_can_enthalpy = initp%can_enthalpy >= rk4aux(ibuff)%rk4min_can_enthalpy   .and.    &
                        initp%can_enthalpy <= rk4aux(ibuff)%rk4max_can_enthalpy
      !------------------------------------------------------------------------------------!



      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !     Here we convert enthalpy into temperature, potential temperature and find some !
      ! derived quantities.  We only calculate these quantities when enthalpy and specific !
      ! humidity are reasonable.                                                           !
      !------------------------------------------------------------------------------------!
      if (ok_can_shv .and. ok_can_enthalpy) then


         !----- Update the canopy air space heat capacity at constant pressure. -----------!
         initp%can_cp = (1.d0 - initp%can_shv) * cpdry8 + initp%can_shv * cph2o8
         !---------------------------------------------------------------------------------!

         !----- Find the canopy air temperature. ------------------------------------------!
         initp%can_temp  = hq2temp8(initp%can_enthalpy,initp%can_shv,.true.)
         !---------------------------------------------------------------------------------!


         !----- Find the new potential temperature. ---------------------------------------!
         initp%can_theta = extemp2theta8(initp%can_exner,initp%can_temp)
         !---------------------------------------------------------------------------------!


         !----- Check whether the potential temperature makes sense or not. ---------------!
         ok_can_theta = initp%can_theta >= rk4aux(ibuff)%rk4min_can_theta .and.            &
                        initp%can_theta <= rk4aux(ibuff)%rk4max_can_theta
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Compute the other canopy air parameters only if the potential temperature  !
         ! makes sense.  Sometimes enthalpy makes sense even though the temperature is     !
         ! bad, because can_shv is way off in the opposite direction of temperature.       !
         !---------------------------------------------------------------------------------!
         if (ok_can_theta) then
            !------------------------------------------------------------------------------!
            !     Find the derived humidity variables.                                     !
            !------------------------------------------------------------------------------!
            initp%can_rhv   = rehuil8(initp%can_prss,initp%can_temp,initp%can_shv,.true.)
            initp%can_ssh   = qslif8(initp%can_prss,initp%can_temp)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      else
         !----- Either enthalpy or specific humidity is screwed, reject theta too. --------!
         ok_can_theta = .false.
         !---------------------------------------------------------------------------------!
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




      !----- Update soil temperature and liquid water fraction. ---------------------------!
      ok_slwater = .true.
      ok_sltemp  = .true.
      do k = rk4site%lsl, nzg
         nsoil = rk4site%ntext_soil(k)
         !----- Check whether soil water is fine. -----------------------------------------!
         ok_slwater = ok_slwater                                                .and.      &
                      initp%soil_water(k) >= rk4aux(ibuff)%rk4min_soil_water(k) .and.      &
                      initp%soil_water(k) <= rk4aux(ibuff)%rk4max_soil_water(k)
         if (ok_slwater) then
            initp%soil_mstpot(k) = matric_potential8(nsoil,initp%soil_water(k))
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      We find the temperature in any case.  If screwy soil water leads to screwy !
         ! temperature, the sanity check will fail twice.                                  !
         !---------------------------------------------------------------------------------!
         call uextcm2tl8(initp%soil_energy(k),initp%soil_water(k)*wdns8,soil8(nsoil)%slcpd &
                        ,initp%soil_tempk(k),initp%soil_fracliq(k))
         ok_sltemp = ok_sltemp                                   .and.                     &
                     initp%soil_tempk(k) >= rk4min_soil_temp     .and.                     &
                     initp%soil_tempk(k) <= rk4max_soil_temp
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Update surface water temperature and liquid water fraction, remembering that    !
      ! inside the RK4 integration, surface water energy is in J/m2. The abs is            !
      ! necessary because surface mass may indeed become too negative during the           !
      ! integration process andif it happens, we want the step to be rejected.             !
      !------------------------------------------------------------------------------------!
      ok_sfcw                = .true.
      ksn                    = initp%nlev_sfcwater
      initp%total_sfcw_depth = 0.d0
      initp%total_sfcw_mass  = 0.d0
      sfcwloop: do k=1,ksn
         if (initp%sfcwater_mass(k) < rk4min_sfcw_mass) then 
            !------------------------------------------------------------------------------!
            !     Surface water mass doesn't make sense.  Skip because the step is going   !
            ! to be rejected and we may not be able to compute the temperature.            !
            !------------------------------------------------------------------------------!
            ok_sfcw = .false.
            cycle sfcwloop
            !------------------------------------------------------------------------------!
         elseif (initp%flag_sfcwater == 1) then
            !------------------------------------------------------------------------------!
            !     Water layer is too thin to be computationally stable, apply the quick    !
            ! heat exchange between the soil top layer and the temporary surface water.    !
            ! Find the total internal energy of the combined pool (top soil layer plus the !
            ! thin temporary surface water).  The units of soil properties are J/m3 for    !
            ! the internal energy, and m3/m3 for soil water, whilst the temporary surface  !
            ! water has units of J/m2 for internal energy and kg/m2 for mass.  We use the  !
            ! standard for the temporary surface water.                                    !
            !------------------------------------------------------------------------------!
            nsoil       = rk4site%ntext_soil(nzg)
            energy_tot  = initp%sfcwater_energy(k) + initp%soil_energy(nzg) * dslz8(nzg)
            wmass_tot   = initp%sfcwater_mass  (k)                                         &
                        + initp%soil_water (nzg) * dslz8(nzg) * wdns8
            hcapdry_tot = soil8(nsoil)%slcpd * dslz8(nzg)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the equilibrium temperature and liquid/ice partition.   Because we !
            ! are assuming thermal equilibrium, the temperature and liquid fraction of the !
            ! attempted layer is the same as the average temperature of the augmented      !
            ! pool.                                                                        !
            !------------------------------------------------------------------------------!
            call uextcm2tl8(energy_tot,wmass_tot,hcapdry_tot                               &
                           ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
            initp%soil_tempk  (nzg) = initp%sfcwater_tempk  (k)
            initp%soil_fracliq(nzg) = initp%sfcwater_fracliq(k) 
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Re-compute the internal energy of the temporary layer, using the temper-  !
            ! ature and fraction of liquid water distribution we have just found, keeping  !
            ! the mass constant.                                                           !
            !------------------------------------------------------------------------------!
            initp%sfcwater_energy(k) = initp%sfcwater_mass(k)                              &
                                     * tl2uint8( initp%sfcwater_tempk(k)                   &
                                               , initp%sfcwater_fracliq(k) )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Re-calculate the top soil internal energy, by removing the attempted sur- !
            ! face water energy from the total energy, and converting it back to J/m3.     !
            ! The total amount of water does not need to be re-calculated at this time.    !
            !------------------------------------------------------------------------------!
            initp%soil_energy(nzg)  = (energy_tot - initp%sfcwater_energy(k)) * dslzi8(nzg)
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !      Convert surface water energy from extensive quantity (J/m2) to inten-   !
            ! sive quantity (J/kg), then update the temperature and liquid water fraction. !
            ! We only do this when there is enough mass, because otherwise the amount of   !
            ! energy is too small that the code would give some unreasonable results.      !
            ! Also, this is the minimum amount of water needed for a layer to exist, any-  !
            ! thing less than that would be entirely absorbed by the soil layer.           !
            !------------------------------------------------------------------------------!
            int_sfcw_energy = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
            call uint2tl8(int_sfcw_energy,initp%sfcwater_tempk(k)                          &
                         ,initp%sfcwater_fracliq(k))
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !----- Aggregate total surface water (mass and depth). ---------------------------!
         initp%total_sfcw_depth = initp%total_sfcw_depth + initp%sfcwater_depth(k)
         initp%total_sfcw_mass  = initp%total_sfcw_mass  + initp%sfcwater_mass (k)
         !---------------------------------------------------------------------------------!
      end do sfcwloop
      !------------------------------------------------------------------------------------!
      !    For non-existent layers of temporary surface water, we copy the temperature and !
      ! liquid fraction from the layer beneath.                                            !
      !------------------------------------------------------------------------------------!
      aboveloop: do k=ksn+1,nzs
         if (k == 1) then
            initp%sfcwater_tempk  (k) = initp%soil_tempk  (nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk  (k) = initp%sfcwater_tempk  (k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do aboveloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the fraction of the canopy covered in snow.  I could not find any       !
      ! reference for the original method (commented out), so I implemented the method     !
      ! used in CLM-4, which is based on:                                                  !
      !                                                                                    !
      ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover  !
      !    fraction and its evaluation over large North American river basins,             !
      !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                         !
      !------------------------------------------------------------------------------------!
      ! initp%snowfac = min(9.9d-1,initp%total_sfcw_depth/initp%veg_height)
      if (initp%total_sfcw_mass > rk4tiny_sfcw_mass) then
         bulk_sfcw_dens = max( fsdns8, min( wdns8                                          &
                                          , initp%total_sfcw_mass/initp%total_sfcw_depth ) )
         initp%snowfac  = max( 0.d0, min( 9.9d-1                                           &
                             , tanh( initp%total_sfcw_depth                                &
                                   / ( ny07_eq04_a8 * soil_rough8                          &
                                     * (bulk_sfcw_dens * fsdnsi8) ** ny07_eq04_m8 ) ) ) )
      else
         initp%snowfac  = 0.d0
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update the virtual pool temperature and liquid water fraction, only when the   !
      ! total mass of the virtual pool makes sense and is not zero.  When the amount of    !
      ! water is too small, we impose the same temperature as the exposed surface (either  !
      ! the top soil layer or the top temporary surface water layer), because the amount   !
      ! of energy is too small to get an acurate temperature anyway, and all the water and !
      ! energy are going to be absorbed by the top soil layer.                             !
      !------------------------------------------------------------------------------------!
      if (initp%virtual_water < rk4min_virt_water) then
         continue
      elseif (abs(initp%virtual_water) > rk4tiny_sfcw_mass) then
         int_virt_energy = initp%virtual_energy / initp%virtual_water
         call uint2tl8(int_virt_energy,initp%virtual_tempk,initp%virtual_fracliq)
      elseif (ksn == 0) then
         initp%virtual_tempk   = initp%soil_tempk(nzg)
         initp%virtual_fracliq = initp%soil_fracliq(nzg)
      else
         initp%virtual_tempk   = initp%sfcwater_tempk(ksn)
         initp%virtual_fracliq = initp%sfcwater_fracliq(ksn)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compute the ground temperature and specific humidity.                          !
      !------------------------------------------------------------------------------------!
      ok_ground = initp%soil_tempk(nzg) >= rk4min_soil_temp                     .and.      &
                  initp%soil_tempk(nzg) <= rk4max_soil_temp                     .and.      &
                  initp%soil_water(nzg) >= rk4aux(ibuff)%rk4min_soil_water(nzg) .and.      &
                  initp%soil_water(nzg) <= rk4aux(ibuff)%rk4max_soil_water(nzg)
      if (ksn > 0) then
         ok_ground = ok_ground                                     .and.                   &
                     initp%sfcwater_tempk(ksn) >= rk4min_sfcw_temp .and.                   &
                     initp%sfcwater_tempk(ksn) <= rk4max_sfcw_temp
      end if
      if (ok_ground) then
         k = max(1,ksn)
         call ed_grndvap8(ksn,initp%soil_water(nzg),initp%soil_tempk(nzg)                  &
                         ,initp%soil_fracliq(nzg),initp%sfcwater_tempk(k)                  &
                         ,initp%snowfac,initp%can_prss                                     &
                         ,initp%can_shv,initp%ground_shv,initp%ground_ssh                  &
                         ,initp%ground_temp,initp%ground_fliq,initp%ggsoil)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Initialise the leaf and wood checks here.                                       !
      !------------------------------------------------------------------------------------!
      ok_leaf = .true.
      ok_wood = .true.
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Now we update leaf and branch properties, based on which kind of branch        !
      ! thermodynamics we're using.                                                        !
      !------------------------------------------------------------------------------------!
      select case(ibranch_thermo)
      case (1)

         !---------------------------------------------------------------------------------!
         !     The combined case.  Here we assume the leaves and wood are eternally in     !
         ! thermal equilibrium, and that the total water is evenly spread over branches    !
         ! and leaves.  We then proceed to reconstruct the leaf and wood internal energy.  !
         ! This, of course, if the step produced decent values, otherwise we bypass and    !
         ! reject the  step.                                                               !
         !---------------------------------------------------------------------------------!
         vegloop: do ico=1,cpatch%ncohorts
            !----- Check whether the cohort is safe... ------------------------------------!
            if (initp%veg_resolvable(ico)) then

               !---------------------------------------------------------------------------!
               !     Find the weighting factors for leaves and branches, so we put the     !
               ! right amount of each on top of each surface.                              !
               !---------------------------------------------------------------------------!
               if (initp%leaf_resolvable(ico) .and. initp%wood_resolvable(ico)) then
                  !------------------------------------------------------------------------!
                  !    Both leaves and branchwood are solved, split according to LAI/WAI.  !
                  !------------------------------------------------------------------------!
                  wgt_leaf = initp%lai(ico) / initp%tai(ico)
                  wgt_wood = 1.d0 - wgt_leaf
               elseif (initp%leaf_resolvable(ico)) then
                  wgt_leaf = 1.d0
                  wgt_wood = 0.d0
               else
                  wgt_leaf = 0.d0
                  wgt_wood = 1.d0
               end if
               !---------------------------------------------------------------------------!



               !----- Find the minimum vegetation water for this cohort. ------------------!
               rk4min_veg_water      = rk4min_veg_lwater * initp%tai(ico)
               rk4min_leaf_water_im2 = rk4aux(ibuff)%rk4min_leaf_water_im2(ico)
               rk4max_leaf_water_im2 = rk4aux(ibuff)%rk4max_leaf_water_im2(ico)
               rk4min_wood_water_im2 = rk4aux(ibuff)%rk4min_wood_water_im2(ico)
               rk4max_wood_water_im2 = rk4aux(ibuff)%rk4max_wood_water_im2(ico)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Update leaf temperature and liquid fraction only if leaf water makes  !
               ! sense.                                                                    !
               !---------------------------------------------------------------------------!
               if ( ( initp%veg_water     (ico) < rk4min_veg_water      ) .or. &
                    ( initp%leaf_water_im2(ico) < rk4min_leaf_water_im2 ) .or. &
                    ( initp%leaf_water_im2(ico) > rk4max_leaf_water_im2 ) .or. &
                    ( initp%wood_water_im2(ico) < rk4min_wood_water_im2 ) .or. &
                    ( initp%wood_water_im2(ico) > rk4max_wood_water_im2 ) ) then

                  !------------------------------------------------------------------------!
                  !     Water is too negative, we must reject the step.  To be safe, we    !
                  ! cheat  and put the total water to both pools.                          !
                  !------------------------------------------------------------------------!
                  initp%leaf_water    (ico) = initp%veg_water    (ico)
                  initp%wood_water    (ico) = initp%veg_water    (ico)
                  !------------------------------------------------------------------------!



                  !----- Flag the step as bad so it doesn't find turbulence parameters. ---!
                  ok_leaf = .false.
                  ok_wood = .false.
                  !------------------------------------------------------------------------!



                  !----- Move on to the next cohort. --------------------------------------!
                  cycle vegloop
                  !------------------------------------------------------------------------!
               else
                  !----- Find the temperature and liquid fraction. ------------------------!
                  call uextcm2tl8(initp%veg_energy(ico)                                    &
                                 ,initp%veg_water(ico) + initp%veg_water_im2(ico)          &
                                 ,initp%veg_hcap(ico),veg_temp,veg_fliq)
                  !------------------------------------------------------------------------!



                  if (veg_temp < rk4min_veg_temp .or. veg_temp > rk4max_veg_temp) then
                     !---------------------------------------------------------------------!
                     !     Temperature is off, we must reject the step.  To be safe, we    !
                     ! cheat and put the bad temperature in both pools.                    !
                     !---------------------------------------------------------------------!
                     initp%leaf_temp(ico) = veg_temp
                     initp%wood_temp(ico) = veg_temp
                     initp%leaf_fliq(ico) = veg_fliq
                     initp%wood_fliq(ico) = veg_fliq
                     !---------------------------------------------------------------------!

                     !---------------------------------------------------------------------!
                     !      Flag the step as bad so it doesn't find turbulence parameters. !
                     !---------------------------------------------------------------------!
                     ok_leaf = .false.
                     ok_wood = .false.
                     !---------------------------------------------------------------------!



                     !----- Move on to the next cohort. -----------------------------------!
                     cycle vegloop
                     !---------------------------------------------------------------------!
                  else
                     !----- Copy veg temperature and liquid fraction to leaves and wood. --!
                     initp%leaf_temp(ico) = veg_temp
                     initp%wood_temp(ico) = veg_temp
                     initp%leaf_fliq(ico) = veg_fliq
                     initp%wood_fliq(ico) = veg_fliq
                     !---------------------------------------------------------------------!


                     !----- Split vegetation water according to LAI and WAI. --------------!
                     initp%leaf_water(ico) = initp%veg_water(ico) * wgt_leaf
                     initp%wood_water(ico) = initp%veg_water(ico) * wgt_wood
                     !---------------------------------------------------------------------!


                     !----- Find lead and wood internal energy. ---------------------------!
                     initp%leaf_energy(ico) = cmtl2uext8( initp%leaf_hcap     (ico)        &
                                                        , initp%leaf_water    (ico)        &
                                                        + initp%leaf_water_im2(ico)        &
                                                        , initp%leaf_temp     (ico)        &
                                                        , initp%leaf_fliq     (ico) )
                     initp%wood_energy(ico) = cmtl2uext8( initp%wood_hcap     (ico)        &
                                                        , initp%wood_water    (ico)        &
                                                        + initp%wood_water_im2(ico)        &
                                                        , initp%wood_temp     (ico)        &
                                                        , initp%wood_fliq     (ico) )

                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Compute the leaf intercellular specific humidity, assumed to be !
                     ! at saturation.                                                      !
                     !---------------------------------------------------------------------!
                     initp%lint_shv(ico) = qslif8(initp%can_prss,initp%leaf_temp(ico))
                     !---------------------------------------------------------------------!
                  end if
               end if
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     For plants with minimal foliage or very sparse patches, fix the leaf  !
               ! and wood temperatures to the canopy air space and force leaf_water and    !
               ! wood_water to be zero.                                                    !
               !---------------------------------------------------------------------------!
               initp%leaf_temp(ico)   = initp%can_temp 
               initp%wood_temp(ico)   = initp%can_temp
               initp%leaf_water(ico)  = 0.d0
               initp%wood_water(ico)  = 0.d0
               initp%veg_water(ico)   = 0.d0
               if (initp%can_temp == t3ple8) then
                  initp%leaf_fliq(ico) = 5.d-1
                  initp%wood_fliq(ico) = 5.d-1
               elseif (initp%can_temp > t3ple8) then
                  initp%leaf_fliq(ico) = 1.d0
                  initp%wood_fliq(ico) = 1.d0
               else
                  initp%leaf_fliq(ico) = 0.d0
                  initp%wood_fliq(ico) = 0.d0
               end if
               initp%leaf_energy(ico) = cmtl2uext8( initp%leaf_hcap     (ico)              &
                                                  , initp%leaf_water    (ico)              &
                                                  + initp%leaf_water_im2(ico)              &
                                                  , initp%leaf_temp     (ico)              &
                                                  , initp%leaf_fliq     (ico) )

               initp%wood_energy(ico) = cmtl2uext8( initp%wood_hcap     (ico)              &
                                                  , initp%wood_water    (ico)              &
                                                  + initp%wood_water_im2(ico)              &
                                                  , initp%wood_temp     (ico)              &
                                                  , initp%wood_fliq     (ico) )
               initp%veg_energy(ico)  = initp%leaf_energy(ico) + initp%wood_energy(ico)
               !---------------------------------------------------------------------------!





               !---------------------------------------------------------------------------!
               !     Compute the leaf intercellular specific humidity, assumed to be at    !
               ! saturation.                                                               !
               !---------------------------------------------------------------------------!
               initp%lint_shv(ico) = qslif8(initp%can_prss,initp%leaf_temp(ico))
               !---------------------------------------------------------------------------!
            end if
         end do vegloop
      case default
         !---------------------------------------------------------------------------------!
         !     Loop over cohorts to update the leaf temperature and liquid fraction.  This !
         ! is done here only if leaves aren't being solved together with branches as a     !
         ! combined pool (either being the only vegetation pool or treated as a separate   !
         ! pool).                                                                          !
         !---------------------------------------------------------------------------------!
         leafloop: do ico=1,cpatch%ncohorts
            !----- Check whether the leaves of this cohort are safe... --------------------!
            if (initp%leaf_resolvable(ico)) then

               !----- Find the minimum leaf water for this cohort. ------------------------!
               rk4min_leaf_water     = rk4min_veg_lwater * initp%lai(ico)
               rk4min_leaf_water_im2 = rk4aux(ibuff)%rk4min_leaf_water_im2(ico) * om_buff_d
               rk4max_leaf_water_im2 = rk4aux(ibuff)%rk4max_leaf_water_im2(ico) * op_buff_d
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Update leaf temperature and liquid fraction only if leaf water makes  !
               ! sense.                                                                    !
               !---------------------------------------------------------------------------!
               if ( ( initp%leaf_water    (ico) < rk4min_leaf_water     ) .or.             &
                    ( initp%leaf_water_im2(ico) < rk4min_leaf_water_im2 ) .or.             &
                    ( initp%leaf_water_im2(ico) > rk4max_leaf_water_im2 ) )                &
               then
                  ok_leaf = .false.
                  cycle leafloop
               else
                  call uextcm2tl8( initp%leaf_energy   (ico)                               &
                                 , initp%leaf_water    (ico)                               &
                                 + initp%leaf_water_im2(ico)                               &
                                 , initp%leaf_hcap     (ico)                               &
                                 , initp%leaf_temp     (ico)                               &
                                 , initp%leaf_fliq     (ico) )
                  if (initp%leaf_temp(ico) < rk4min_veg_temp .or.                          &
                      initp%leaf_temp(ico) > rk4max_veg_temp) then
                     ok_leaf = .false.
                     cycle leafloop
                  else
                     !---------------------------------------------------------------------!
                     !     Compute the leaf intercellular specific humidity, assumed to be !
                     ! at saturation.                                                      !
                     !---------------------------------------------------------------------!
                     initp%lint_shv(ico) = qslif8(initp%can_prss,initp%leaf_temp(ico))
                     !---------------------------------------------------------------------!
                  end if
               end if
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     For plants with minimal foliage or very sparse patches, fix the leaf  !
               ! temperature to the canopy air space and force leaf_water to be zero.      !
               !---------------------------------------------------------------------------!
               initp%leaf_temp(ico)   = initp%can_temp
               initp%leaf_water(ico)  = 0.d0
               initp%leaf_energy(ico) = cmtl2uext8( initp%leaf_hcap     (ico)              &
                                                  , initp%leaf_water    (ico)              &
                                                  + initp%leaf_water_im2(ico)              &
                                                  , initp%leaf_temp     (ico)              &
                                                  , initp%leaf_fliq     (ico) )

               if (initp%leaf_temp(ico) == t3ple8) then
                  initp%leaf_fliq(ico) = 5.d-1
               elseif (initp%leaf_temp(ico) > t3ple8) then
                  initp%leaf_fliq(ico) = 1.d0
               else
                  initp%leaf_fliq(ico) = 0.d0
               end if
               !---------------------------------------------------------------------------!





               !---------------------------------------------------------------------------!
               !     Compute the leaf intercellular specific humidity, assumed to be at    !
               ! saturation.                                                               !
               !---------------------------------------------------------------------------!
               initp%lint_shv(ico) = qslif8(initp%can_prss,initp%leaf_temp(ico))
               !---------------------------------------------------------------------------!
            end if
         end do leafloop
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Loop over cohorts to update the wood temperature and liquid fraction.  This !
         ! is done only when branches are being solved by themselves.                      !
         !---------------------------------------------------------------------------------!
         woodloop: do ico=1,cpatch%ncohorts
            !----- Check whether the leaves of this cohort are safe... --------------------!
            if (initp%wood_resolvable(ico)) then

               !----- Find the minimum leaf water for this cohort. ------------------------!
               rk4min_wood_water     = rk4min_veg_lwater * initp%wai(ico)
               rk4min_wood_water_im2 = rk4aux(ibuff)%rk4min_wood_water_im2(ico)            &
                                     * (1.d0-rk4eps)
               rk4max_wood_water_im2 = rk4aux(ibuff)%rk4max_wood_water_im2(ico)            &
                                     * (1.d0+rk4eps)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Update wood temperature and liquid fraction only if wood water makes  !
               ! sense.                                                                    !
               !---------------------------------------------------------------------------!
               if ( ( initp%wood_water    (ico) < rk4min_wood_water     ) .or.             &
                    ( initp%wood_water_im2(ico) < rk4min_wood_water_im2 ) .or.             &
                    ( initp%wood_water_im2(ico) > rk4max_wood_water_im2 ) )                &
               then
                  ok_wood = .false.
                  cycle woodloop
               else
                  call uextcm2tl8( initp%wood_energy   (ico)                               &
                                 , initp%wood_water    (ico)                               &
                                 + initp%wood_water_im2(ico)                               &
                                 , initp%wood_hcap     (ico)                               &
                                 , initp%wood_temp     (ico)                               &
                                 , initp%wood_fliq     (ico) )
                  if (initp%wood_temp(ico) < rk4min_veg_temp .or.                          &
                      initp%wood_temp(ico) > rk4max_veg_temp) then
                     ok_wood = .false.
                     cycle woodloop
                  end if
               end if
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     In case we are not solving branches, or for cohorts that are very     !
               ! sparse, fix the wood temperature to the canopy air space and force        !
               ! wood_water to be zero.                                                    !
               !---------------------------------------------------------------------------!
               initp%wood_temp(ico)   = initp%can_temp
               initp%wood_water(ico)  = 0.d0
               initp%wood_energy(ico) = cmtl2uext8( initp%wood_hcap     (ico)              &
                                                  , initp%wood_water    (ico)              &
                                                  + initp%wood_water_im2(ico)              &
                                                  , initp%wood_temp     (ico)              &
                                                  , initp%wood_fliq     (ico) )

               if (initp%wood_temp(ico) == t3ple8) then
                  initp%wood_fliq(ico) = 5.d-1
               elseif (initp%wood_temp(ico) > t3ple8) then
                  initp%wood_fliq(ico) = 1.d0
               else
                  initp%wood_fliq(ico) = 0.d0
               end if
               !---------------------------------------------------------------------------!
            end if
         end do woodloop
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Update vegetation properties, so everything is neat and consistent.  This   !
         ! can be done for non-resolvable and unrealistic cohorts, it's just a sum and the !
         ! time step is about to be rejected anyway.                                       !
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts
            initp%veg_energy   (ico) = initp%leaf_energy   (ico) + initp%wood_energy   (ico)
            initp%veg_water    (ico) = initp%leaf_water    (ico) + initp%wood_water    (ico)
            initp%veg_water_im2(ico) = initp%leaf_water_im2(ico) + initp%wood_water_im2(ico)
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Compute canopy turbulence properties. ----------------------------------------!
      if (ok_leaf .and. ok_wood .and. ok_can_shv .and. ok_can_theta .and. ok_ground) then
         call canopy_turbulence8(csite,initp,ipa,ibuff)
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine update_diagnostic_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine updates canopy air density (mass-based and molar-based), the CAS !
   ! capacities, and the integrates the CAS density effect.                                !
   !---------------------------------------------------------------------------------------!
   subroutine update_density_vars(initp,prevp)
      use rk4_coms              , only : checkbudget           & ! intent(in)
                                       , rk4patchtype          ! ! structure
      use therm_lib8            , only : idealdenssh8          & ! function
                                       , idealdmolsh8          ! ! function
      use canopy_struct_dynamics, only : can_whccap8           ! ! subroutine
      use budget_utils          , only : ddens_dt_effect8      ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: initp
      type(rk4patchtype) , target     :: prevp
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Update density (mass- and molar-based).                                        !
      !------------------------------------------------------------------------------------!
      initp%can_rhos  = idealdenssh8(initp%can_prss,initp%can_temp,initp%can_shv)
      initp%can_dmol  = idealdmolsh8(initp%can_prss,initp%can_temp,initp%can_shv)
      !------------------------------------------------------------------------------------!



      !----- Update the canopy air space capacities. --------------------------------------!
      call can_whccap8(initp%can_rhos,initp%can_dmol,initp%can_depth                       &
                      ,initp%wcapcan ,initp%hcapcan ,initp%ccapcan                         &
                      ,initp%wcapcani,initp%hcapcani,initp%ccapcani  )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Integrate the density effect.                                                  !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         initp%co2budget_denseffect = initp%co2budget_denseffect                           &
                                    + ddens_dt_effect8( initp%can_depth, prevp%can_dmol    &
                                                      , initp%can_dmol , initp%can_co2     )
         initp%ebudget_denseffect   = initp%ebudget_denseffect                             &
                                    + ddens_dt_effect8( initp%can_depth, prevp%can_rhos    &
                                                      , initp%can_rhos , initp%can_enthalpy)
         initp%wbudget_denseffect   = initp%wbudget_denseffect                             &
                                    + ddens_dt_effect8( initp%can_depth, prevp%can_rhos    &
                                                      , initp%can_rhos , initp%can_shv     )
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_density_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine performs the following tasks:                                      !
   ! 1. Check how many layers of temporary water or snow we have, and include the virtual  !
   !    pools at the topmost if needed;                                                    !
   ! 2. Force thermal equilibrium between topmost soil layer and a single snow/water layer !
   !    if the layer is too thin;                                                          !
   ! 3. Compute the amount of mass each layer has, and redistribute them accordingly.      !
   ! 4. Percolates excessive liquid water if needed.                                       !
   !---------------------------------------------------------------------------------------!
   subroutine adjust_sfcw_properties(nzg,nzs,initp,hdid,csite,ipa)

      use ed_state_vars , only : sitetype              ! ! structure
      use rk4_coms      , only : rk4patchtype          & ! structure
                               , rk4site               & ! intent(in)
                               , rk4min_sfcw_mass      & ! intent(in)
                               , rk4min_virt_water     & ! intent(in)
                               , rk4water_stab_thresh  & ! intent(in)
                               , rk4tiny_sfcw_mass     & ! intent(in)
                               , rk4tiny_sfcw_depth    & ! intent(in)
                               , rk4min_can_shv        & ! intent(in)
                               , rk4snowmin            & ! intent(in)
                               , ipercol               & ! intent(in)
                               , dtrk4i                & ! intent(in)
                               , tiny_offset           & ! intent(in)
                               , checkbudget           ! ! intent(in)
      use ed_state_vars , only : sitetype              & ! structure
                               , patchtype             ! ! structure
      use ed_misc_coms  , only : current_time          & ! intent(in)
                               , fast_diagnostics      ! ! intent(in)
      use soil_coms     , only : soil8                 & ! intent(in)
                               , dslz8                 & ! intent(in)
                               , dslzi8                & ! intent(in)
                               , thick                 ! ! intent(in)
      use consts_coms   , only : t3ple8                & ! intent(in)
                               , wdns8                 & ! intent(in)
                               , wdnsi8                & ! intent(in)
                               , fdnsi8                & ! intent(in)
                               , uiliqt38              & ! intent(in)
                               , wdnsi8                & ! intent(in)
                               , fdnsi8                & ! intent(in)
                               , fsdnsi8               ! ! intent(in)
      use therm_lib8    , only : uint2tl8              & ! subroutine
                               , uextcm2tl8            & ! subroutine
                               , tl2uint8              & ! function
                               , tq2enthalpy8          ! ! function

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype)     , target     :: initp      !< Current RK4 patch
      type(sitetype)         , target     :: csite      !< Current site
      real(kind=8)           , intent(in) :: hdid       !< Time step
      integer                , intent(in) :: ipa        !< Current patch ID
      integer                , intent(in) :: nzg        !< Number of soil layers
      integer                , intent(in) :: nzs        !< Max. number of sfc. water layers
      !----- Local variables --------------------------------------------------------------!
      integer                             :: kold
      integer                             :: newlayers
      integer                             :: nlayers
      integer                             :: ksn
      integer                             :: ksnnew
      integer                             :: k
      !----- Control variables ------------------------------------------------------------!
      real(kind=8)                        :: hdidi
      real(kind=8)                        :: wtold
      real(kind=8)                        :: wtnew
      real(kind=8), dimension(nzs)        :: newsfcw_mass
      real(kind=8), dimension(nzs)        :: newsfcw_energy
      real(kind=8), dimension(nzs)        :: newsfcw_depth
      real(kind=8)                        :: wdiff
      real(kind=8)                        :: sum_sfcw_mass
      real(kind=8)                        :: sum_sfcw_energy
      real(kind=8)                        :: sum_sfcw_depth
      real(kind=8)                        :: energy_free
      real(kind=8)                        :: wmass_free
      real(kind=8)                        :: depth_free
      real(kind=8)                        :: tempk_free
      real(kind=8)                        :: fracliq_free
      real(kind=8)                        :: energy_latent
      real(kind=8)                        :: energy_available
      real(kind=8)                        :: wmass_available
      real(kind=8)                        :: depth_available
      real(kind=8)                        :: energy_needed
      real(kind=8)                        :: wmass_needed
      real(kind=8)                        :: depth_needed
      real(kind=8)                        :: tempk_needed
      real(kind=8)                        :: fracliq_needed
      real(kind=8)                        :: wmass_perc
      real(kind=8)                        :: energy_perc
      real(kind=8)                        :: depth_perc
      real(kind=8)                        :: i_energy_try
      real(kind=8)                        :: energy_try
      real(kind=8)                        :: wmass_try
      real(kind=8)                        :: depth_try
      real(kind=8)                        :: temp_try
      real(kind=8)                        :: fliq_try
      real(kind=8)                        :: energy_runoff
      real(kind=8)                        :: wmass_runoff
      real(kind=8)                        :: energy_tot
      real(kind=8)                        :: wmass_tot
      real(kind=8)                        :: hcapdry_tot
      real(kind=8)                        :: wmass_room
      real(kind=8)                        :: energy_room
      real(kind=8)                        :: snden
      real(kind=8)                        :: sndenmin
      real(kind=8)                        :: sndenmax
      real(kind=8)                        :: Cr               ! snow waterholding capacity
      real(kind=8)                        :: gi               ! Partial density of ice
      integer                             :: nsoil
      !----- Variables used for the water and energy budget. ------------------------------!
      character(len=1)                    :: flag_first
      character(len=1)                    :: flag_stable
      character(len=1)                    :: flag_percol
      character(len=1)                    :: flag_second
      integer                             :: flag_sfcw_beg
      integer                             :: flag_sfcw_end
      integer                             :: nlev_sfcw_beg
      integer                             :: nlev_sfcw_end
      real(kind=8)                        :: wmass_cas_beg
      real(kind=8)                        :: wmass_cas_end
      real(kind=8)                        :: enthalpy_cas_beg
      real(kind=8)                        :: enthalpy_cas_end
      real(kind=8)                        :: wmass_virtual_beg
      real(kind=8)                        :: energy_virtual_beg
      real(kind=8)                        :: wmass_sfcw_beg
      real(kind=8)                        :: energy_sfcw_beg
      real(kind=8)                        :: wmass_soil_beg
      real(kind=8)                        :: energy_soil_beg
      real(kind=8)                        :: wmass_total_beg
      real(kind=8)                        :: energy_total_beg
      real(kind=8)                        :: wmass_virtual_end
      real(kind=8)                        :: energy_virtual_end
      real(kind=8)                        :: wmass_sfcw_end
      real(kind=8)                        :: energy_sfcw_end
      real(kind=8)                        :: wmass_soil_end
      real(kind=8)                        :: energy_soil_end
      real(kind=8)                        :: wmass_total_end
      real(kind=8)                        :: energy_total_end
      real(kind=8)                        :: wmass_total_rch
      real(kind=8)                        :: energy_total_rch
      !----- Constants --------------------------------------------------------------------!
      real(kind=8)           , parameter  :: Crmin      = 3.d-2
      real(kind=8)           , parameter  :: Crmax      = 1.d-1
      real(kind=8)           , parameter  :: ge         = 2.d2
      real(kind=8)           , parameter  :: sfcw_toler = 1.d-6
      character(len= 9)      , parameter  :: fmtc       = '(a,13x,a)'
      character(len=10)      , parameter  :: fmti       = '(a,1x,i14)'
      character(len=13)      , parameter  :: fmtf       = '(a,1x,es14.7)'
      character(len=27)      , parameter  :: fmtt       = '(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !----- External function. -----------------------------------------------------------!
      real                   , external   :: sngloff
      !------------------------------------------------------------------------------------!


      !----- Find the inverse of the time step. -------------------------------------------!
      hdidi      = 1.d0 / hdid
      !------------------------------------------------------------------------------------!


      !----- Copy the original number of temporary surface water layers to ksn. -----------!
      ksn       = initp%nlev_sfcwater
      !------------------------------------------------------------------------------------!


      !----- Copy the soil type at the topmost level to nsoil. ----------------------------!
      nsoil     = rk4site%ntext_soil(nzg)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Determine the total amount of temporary surface water available as well as    !
      ! derived properties.                                                                !
      !------------------------------------------------------------------------------------!
      sum_sfcw_mass   = 0.d0
      sum_sfcw_energy = 0.d0
      sum_sfcw_depth  = 0.d0
      do k=1,ksn
         sum_sfcw_mass   = sum_sfcw_mass   + initp%sfcwater_mass  (k)
         sum_sfcw_energy = sum_sfcw_energy + initp%sfcwater_energy(k)
         sum_sfcw_depth  = sum_sfcw_depth  + initp%sfcwater_depth (k)
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Initialise the budget variables.                                              !
      !------------------------------------------------------------------------------------!
      flag_sfcw_beg      = initp%flag_sfcwater
      nlev_sfcw_beg      = initp%nlev_sfcwater
      wmass_cas_beg      = initp%can_shv      * initp%wcapcan
      enthalpy_cas_beg   = initp%can_enthalpy * initp%hcapcan
      wmass_virtual_beg  = initp%virtual_water
      energy_virtual_beg = initp%virtual_energy
      wmass_sfcw_beg     = sum_sfcw_mass
      energy_sfcw_beg    = sum_sfcw_energy
      wmass_soil_beg     = initp%soil_water(nzg)  * dslz8(nzg) * wdns8
      energy_soil_beg    = initp%soil_energy(nzg) * dslz8(nzg)
      wmass_total_beg    = wmass_virtual_beg  + wmass_sfcw_beg  + wmass_soil_beg           &
                         + wmass_cas_beg
      energy_total_beg   = energy_virtual_beg + energy_sfcw_beg + energy_soil_beg          &
                         + enthalpy_cas_beg
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check the total amount of water that has just fallen plus the amount that is   !
      ! already sitting over the top soil layer.  We must do this as the first step be-    !
      ! cause we may want to eliminate this water by binding it to the top soil layer in   !
      ! case there is too little water.                                                    !
      !------------------------------------------------------------------------------------!
      if (initp%virtual_water < rk4min_virt_water .or. sum_sfcw_mass < rk4min_sfcw_mass )  &
      then
         !---------------------------------------------------------------------------------!
         !     Either the virtual layer or the temporary surface water has too negative    !
         ! mass, so this step doesn't make sense.  We quit the sub-routine here so the     !
         ! sanity check can reject this step.                                              !
         !---------------------------------------------------------------------------------!
         return
         !---------------------------------------------------------------------------------!
      elseif ((initp%virtual_water + sum_sfcw_mass) < 0.d0) then
         !---------------------------------------------------------------------------------!
         !     The mass of the potential new temporary surface water is within bounds but  !
         ! it is negative.  This can only happen if the layer evaporated more than what it !
         ! should, so we condense some of the water back from the canopy air space.  If it !
         ! is going to deplete the canopy air space specific humidity too much, then we    !
         ! leave the remainder to be stolen from the top soil, but this is dangerous be-   !
         ! cause the soil may be too dry too.                                              !
         !---------------------------------------------------------------------------------!
         wmass_needed  = - (initp%virtual_water  + sum_sfcw_mass  )
         energy_needed = - (initp%virtual_energy + sum_sfcw_energy)
         depth_needed  = - (initp%virtual_depth  + sum_sfcw_depth )
         !---------------------------------------------------------------------------------!



         !----- Find the amount available at the canopy air space. ------------------------!
         wmass_available = initp%wcapcan * (initp%can_shv - 5.d0 * rk4min_can_shv)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Check whether we canopy air space can provide the water needed.             !
         !---------------------------------------------------------------------------------!
         if ( wmass_available >= wmass_needed) then
            !------------------------------------------------------------------------------!
            !     Canopy air space has sufficient water.                                   !
            !------------------------------------------------------------------------------!
            flag_first = 'A'



            !----- Find the temperature associated with the water needed. -----------------!
            call uint2tl8(energy_needed/wmass_needed,tempk_needed,fracliq_needed)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Water vapour is at higher energetic state than ice or liquid water.  Find !
            ! the excess energy (enthalpy really) associated with condensation.            !
            !------------------------------------------------------------------------------!
            energy_latent = wmass_needed * tq2enthalpy8(tempk_needed,1.d0,.true.)          &
                          - energy_needed
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    There is enough water vapour. The transfer will require phase change, so  !
            ! the energy transfer will be a latent heat flux.  Remove the water from the   !
            ! canopy air space.  The energy lost by the canopy air space to pad the miss-  !
            ! ing water at the virtual+temporary surface water layer must go somewhere, so !
            ! we add it to the soil because it is the closest to the pounding layer.       !
            !------------------------------------------------------------------------------!
            initp%can_shv          = initp%can_shv - wmass_needed  * initp%wcapcani
            initp%can_enthalpy     = initp%can_enthalpy                                    &
                                   - ( energy_needed + energy_latent ) * initp%hcapcani
            initp%soil_energy(nzg) = initp%soil_energy(nzg) + energy_latent * dslzi8(nzg)
            initp%avg_vapor_gc     = initp%avg_vapor_gc - wmass_needed * hdidi
            !------------------------------------------------------------------------------!


            !----- Set free water to zero. ------------------------------------------------!
            wmass_free  = 0.d0
            energy_free = 0.d0
            depth_free  = 0.d0
            !------------------------------------------------------------------------------!

         elseif (wmass_available > 0.d0) then
            !------------------------------------------------------------------------------!
            !     Canopy air space does not have sufficient water.  Take all the water     !
            ! vapour that can be taken, then figure out what to do with the remaining      !
            ! demand.                                                                      !
            !------------------------------------------------------------------------------!
            flag_first = 'B'



            !----- Find the temperature associated with the water needed. -----------------!
            call uint2tl8(energy_needed/wmass_needed,tempk_needed,fracliq_needed)
            !------------------------------------------------------------------------------!



            !----- Find the energy and depth available to reduce debt. --------------------!
            energy_available = wmass_available * energy_needed / wmass_needed
            depth_available  = wmass_available *                                           &
                 ( fracliq_needed * wdnsi8 + (1.d0-fracliq_needed * fdnsi8))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Water vapour is at higher energetic state than ice or liquid water.  Find !
            ! the excess energy (enthalpy really) associated with condensation.            !
            !------------------------------------------------------------------------------!
            energy_latent = wmass_available * tq2enthalpy8(tempk_needed,1.d0,.true.)       &
                          - energy_available
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    There is not enough water vapour. Dry down to the minimum, and correct    !
            ! energy.  Since there is so little negative mass needed, we include the       !
            ! latent heat associated with this condensation to the soil, because otherwise !
            ! we could end up with energy and water mass with opposite signs.              !
            !------------------------------------------------------------------------------!
            initp%can_shv          = initp%can_shv                                         &
                                   - wmass_available  * initp%wcapcani
            initp%can_enthalpy     = initp%can_enthalpy                                    &
                                   - ( energy_available + energy_latent ) * initp%hcapcani
            initp%soil_energy(nzg) = initp%soil_energy(nzg)                                &
                                   + energy_latent * dslzi8(nzg)
            initp%avg_vapor_gc     = initp%avg_vapor_gc  - wmass_available * hdidi
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Set free water to the remaining demand (and thus negative).               !
            !------------------------------------------------------------------------------!
            wmass_free  = wmass_available  - wmass_needed 
            energy_free = energy_available - energy_needed
            depth_free  = depth_available  - depth_needed
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !    There is not any water vapour.  The demand remains the same (but it must  !
            ! be negative).                                                                !
            !------------------------------------------------------------------------------!
            flag_first  = 'C'
            wmass_free  = - wmass_needed
            energy_free = - energy_needed
            depth_free  = - depth_needed
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!

         !----- Reset both the temporary surface water and the virtual layer. -------------!
         initp%virtual_water      = 0.d0
         initp%virtual_energy     = 0.d0
         initp%virtual_depth      = 0.d0
         initp%sfcwater_mass  (:) = 0.d0
         initp%sfcwater_energy(:) = 0.d0
         initp%sfcwater_depth (:) = 0.d0
         !---------------------------------------------------------------------------------!


         !----- Set ksnnew to zero to force all free water to go to the soil. -------------!
         ksnnew                   = 0
         !---------------------------------------------------------------------------------!

      elseif ((initp%virtual_water + sum_sfcw_mass) < rk4tiny_sfcw_mass) then
         !---------------------------------------------------------------------------------!
         !     The mass of the potential new temporary surface water is within bounds but  !
         ! it is too small to be maintained.  We add both the virtual mass and the total   !
         ! surface water and dump in the free water, but set ksnnew to zero so all the     !
         ! water is infiltrated in the top soil layer.                                     !
         !---------------------------------------------------------------------------------!
         flag_first  = 'D'
         wmass_free  = initp%virtual_water  + sum_sfcw_mass
         energy_free = initp%virtual_energy + sum_sfcw_energy
         depth_free  = initp%virtual_depth  + sum_sfcw_depth
         !---------------------------------------------------------------------------------!


         !----- Reset both the temporary surface water and the virtual layer. -------------!
         initp%virtual_water      = 0.d0
         initp%virtual_energy     = 0.d0
         initp%virtual_depth      = 0.d0
         initp%sfcwater_mass  (:) = 0.d0
         initp%sfcwater_energy(:) = 0.d0
         initp%sfcwater_depth (:) = 0.d0
         !---------------------------------------------------------------------------------!


         !----- Set ksnnew to zero to force all free water to go to the soil. -------------!
         ksnnew                   = 0
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     The mass of the potential new temporary surface water is within bounds and  !
         ! could create at least one layer.  If there is already a temporary surface water !
         ! or snow layer, the new amount is initially put there, otherwise, we attempt to  !
         ! create the first layer.                                                         !
         !---------------------------------------------------------------------------------!
         flag_first           = 'E'
         wmass_free           = initp%virtual_water
         energy_free          = initp%virtual_energy
         depth_free           = initp%virtual_depth
         initp%virtual_water  = 0.d0
         initp%virtual_energy = 0.d0
         initp%virtual_depth  = 0.d0
         ksnnew               = max(ksn,1)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Update the prognostic and diagnostic variables by adding the free standing      !
      ! water.  Then we check the size of the temporary surface water layers, and update   !
      ! the temperature in a way that ensure the layer stability.  During this process, we !
      ! update the total temporary surface water mass, energy, and depth, which will be    !
      ! used later in the sub-routine.                                                     !
      !------------------------------------------------------------------------------------!
      sum_sfcw_mass   = 0.d0
      sum_sfcw_energy = 0.d0
      sum_sfcw_depth  = 0.d0
      flag_stable     = 'Z'
      flag_percol     = 'Z'
      do k = ksnnew,1,-1
         !---------------------------------------------------------------------------------!
         !    Find the potential mass, energy, and depth of the temporary layer if all the !
         ! free water became part of this layer.                                           !
         !---------------------------------------------------------------------------------!
         energy_try = initp%sfcwater_energy(k) + energy_free
         wmass_try  = initp%sfcwater_mass  (k) + wmass_free
         depth_try  = initp%sfcwater_depth (k) + depth_free
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    In case this is a single layer, and a very thin one, we may have a hard time !
         ! achieving numerical stability.  We can treat this case like the leaf case, in   !
         ! the sense that the water sitting on the top of the surface is in thermal        !
         ! equilibrium with the surface.                                                   !
         !---------------------------------------------------------------------------------!
         if (ksnnew == 1 .and. wmass_try < rk4water_stab_thresh) then
            !---- Flag layer as not stable. -----------------------------------------------!
            flag_stable = 'A'
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the total internal energy of the combined pool (top soil layer plus !
            ! the thin temporary surface water).  The units of soil properties are J/m3    !
            ! for the internal energy, and m3/m3 for soil water, whilst the temporary sur- !
            ! face water has units of J/m2 for internal energy and kg/m2 for mass.  We use !
            ! the standard for the temporary surface water.                                !
            !------------------------------------------------------------------------------!
            energy_tot  = energy_try * dslzi8(nzg)          + initp%soil_energy(nzg)
            wmass_tot   = wmass_try  * dslzi8(nzg) * wdnsi8 + initp%soil_water (nzg)
            hcapdry_tot = soil8(nsoil)%slcpd
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the equilibrium temperature and liquid/ice partition.   Because we !
            ! are assuming thermal equilibrium, the temperature and liquid fraction of the !
            ! attempted layer is the same as the average temperature of the augmented      !
            ! pool.                                                                        !
            !------------------------------------------------------------------------------!
            call uextcm2tl8(energy_tot,wmass_tot,hcapdry_tot,temp_try,fliq_try)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Re-compute the internal energy of the temporary layer, using the temper-  !
            ! ature and fraction of liquid water distribution we have just found, keeping  !
            ! the mass constant.                                                           !
            !------------------------------------------------------------------------------!
            energy_try = wmass_try * tl2uint8(temp_try,fliq_try)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Re-calculate the top soil internal energy, by removing the attempted sur- !
            ! face water energy from the total energy, and converting it back to J/m3.     !
            ! The total amount of water does not need to be re-calculated at this time.    !
            !------------------------------------------------------------------------------!
            initp%soil_energy(nzg)  = energy_tot - energy_try * dslzi8(nzg)
            !------------------------------------------------------------------------------!
         else
            !---- Flag layer as stable. ---------------------------------------------------!
            flag_stable = 'B'
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Layer is computationally stable, find temperature and liquid fraction   !
            ! of the attempted layer.                                                      !
            !------------------------------------------------------------------------------!
            i_energy_try = energy_try / wmass_try
            call uint2tl8(i_energy_try,temp_try,fliq_try)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Determine a first guess for the amount of mass that can be lost from this   !
         ! layer through percolation (wmass_perc).                                         !
         !---------------------------------------------------------------------------------!
         select case (ipercol)
         case (0)
            !------------------------------------------------------------------------------!
            !     Original method, from LEAF-3.  Shed liquid in excess of a 1:9            !
            ! liquid-to-ice ratio through percolation.                                     !
            !------------------------------------------------------------------------------!
            wmass_perc  = max(0.d0, wmass_try * (fliq_try - 1.d-1) / 9.d-1)
            !------------------------------------------------------------------------------!
         case (1,2)
            !------------------------------------------------------------------------------!
            !    Alternative "free" water calculation.                                     !
            !    Anderson (1976), NOAA Tech Report NWS 19.                                 !
            !------------------------------------------------------------------------------!
            gi          = wmass_try/max(rk4tiny_sfcw_depth,depth_try) * (1.d0 - fliq_try)
            Cr          = max(Crmin, Crmin + (Crmax - Crmin) * (ge - gi) / ge)
            wmass_perc  = max(0.d0,wmass_try * (fliq_try - Cr / (1.d0 + Cr)))
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Determinte whether the layer beneath the current one is another temporary   !
         ! surface water/snow layer, or the top soil layer.  In case it is the latter, we  !
         ! must check whether there is enough room for the percolate water to infiltrate   !
         ! (i.e., the soil will not become super-saturated), in which case we must reduce  !
         ! the total amount of percolation.                                                !
         !---------------------------------------------------------------------------------!
         if (k == 1) then
            !------------------------------------------------------------------------------!
            !     Compute the available "room" for water at the top soil layer.  We must   !
            ! multiply by density and depth to make sure that the units match.             !
            !------------------------------------------------------------------------------!
            wmass_room = max(0.d0, soil8(nsoil)%soilbp - initp%soil_water(nzg))            &
                       * wdns8 * dslz8(nzg) 
            wmass_perc = min(wmass_perc,wmass_room)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Re-calculate the total water mass and energy of this temporary surface      !
         ! water.  Here we must check whether the soil layer would be with too little      !
         ! mass, and if that is the case, we will eliminate the layer by forcing the tiny  !
         ! left-over to go to the layer beneath.                                           !
         !---------------------------------------------------------------------------------!
         if ( (wmass_try - wmass_perc) > rk4tiny_sfcw_mass) then
            if (k == 1) flag_percol = 'A'
            !------------------------------------------------------------------------------!
            !      Enough mass to keep this layer.                                         !
            !------------------------------------------------------------------------------!
            !----- Compute internal energy and depth associated with percolated water. ----!
            energy_perc = wmass_perc * tl2uint8(temp_try,1.d0)
            depth_perc  = wmass_perc * wdnsi8
            !----- Find the new water mass and energy for this layer. ---------------------!
            initp%sfcwater_mass  (k) = wmass_try  - wmass_perc
            initp%sfcwater_energy(k) = energy_try - energy_perc
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Calculate density and depth of snow.  Start with the difference of      !
            ! depths, but then we adjust it because the loss through percolation changes   !
            ! the ratio between ice and liquid in this layer                               !
            !------------------------------------------------------------------------------!
            initp%sfcwater_depth (k) = depth_try  - depth_perc
            snden    = initp%sfcwater_mass(k)                                              &
                     / max(rk4tiny_sfcw_depth,initp%sfcwater_depth(k))
            sndenmax = wdns8
            sndenmin = max(3.d1, 2.d2 * (wmass_free + wmass_perc) / initp%sfcwater_mass(k) )
            snden    = min(sndenmax, max(sndenmin,snden))
            initp%sfcwater_depth (k) = initp%sfcwater_mass(k) / snden
            !------------------------------------------------------------------------------!
         else
            if (k == 1) flag_percol = 'B'

            !------------------------------------------------------------------------------!
            !      The layer would be too small, eliminate mass from this layer and send   !
            ! all mass to the layer beneath as percolated water.                           !
            !------------------------------------------------------------------------------!
            initp%sfcwater_mass  (k) = 0.d0
            initp%sfcwater_energy(k) = 0.d0
            initp%sfcwater_depth (k) = 0.d0
            wmass_perc               = wmass_try
            energy_perc              = energy_try
            depth_perc               = depth_try
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Integrate the total temporary surface water properties.                     !
         !---------------------------------------------------------------------------------!
         sum_sfcw_mass   = sum_sfcw_mass   + initp%sfcwater_mass  (k)
         sum_sfcw_energy = sum_sfcw_energy + initp%sfcwater_energy(k)
         sum_sfcw_depth  = sum_sfcw_depth  + initp%sfcwater_depth (k)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     The water available for the layer beneath is going to be the total per-     !
         ! colated water.                                                                  !
         !---------------------------------------------------------------------------------!
         wmass_free  = wmass_perc
         energy_free = energy_perc
         depth_free  = depth_perc
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     There may be some free standing water left.  We dump what we can into the top  !
      ! soil layer.  Then, in case there is still some water to be removed, we evaporate   !
      ! what is left.                                                                      !
      !------------------------------------------------------------------------------------!
      wmass_runoff  = 0.d0
      energy_runoff = 0.d0
      if (wmass_free > 0.d0) then
         !---------------------------------------------------------------------------------!
         !     No water is needed.                                                         !
         !---------------------------------------------------------------------------------!
         wmass_needed  = 0.d0
         energy_needed = 0.d0
         depth_needed  = 0.d0
         !---------------------------------------------------------------------------------!


         !----- Find the amount of water (and energy) the top soil layer can receive. -----!
         wmass_room  = max(0.d0, soil8(nsoil)%soilbp - initp%soil_water(nzg))              &
                     * wdns8 * dslz8(nzg)
         energy_room = energy_free * wmass_room / wmass_free
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Check whether the top soil layer can absorb all the standing water.        !
         !---------------------------------------------------------------------------------!
         if (wmass_room >= wmass_free) then

            !------------------------------------------------------------------------------!
            !     There is enough space in the top soil layer for the remaining water, put !
            ! all the free water there.                                                    !
            !------------------------------------------------------------------------------!
            flag_second = 'A'
            !------------------------------------------------------------------------------!


            !------ Update soil water and internal energy. --------------------------------!
            initp%soil_water(nzg)  = initp%soil_water(nzg)  + wmass_free  * dslzi8(nzg)    &
                                   * wdnsi8
            initp%soil_energy(nzg) = initp%soil_energy(nzg) + energy_free * dslzi8(nzg)
            !------------------------------------------------------------------------------!


            !----- All remaining water was taken care of, set free water to zero. ---------!
            wmass_free  = 0.d0
            energy_free = 0.d0
            depth_free  = 0.d0
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     There is not enough space in the top soil layer for the remaining water, !
            ! put what we can there, and eliminate the remaining as runoff.                !
            !------------------------------------------------------------------------------!
            flag_second = 'B'
            !------------------------------------------------------------------------------!


            !----- Send excess free water (and associated energy) to runoff. --------------!
            wmass_runoff  = wmass_free  - wmass_room
            energy_runoff = energy_free - energy_room
            !------------------------------------------------------------------------------!



            !----- Find the temperature associated with the free water. -------------------!
            call uint2tl8(energy_free/wmass_free,tempk_free,fracliq_free)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Dump what we can dump on the top soil layer.                             !
            !------------------------------------------------------------------------------!
            initp%soil_water(nzg)  = initp%soil_water(nzg)  + wmass_room  * dslzi8(nzg)    &
                                   * wdnsi8
            initp%soil_energy(nzg) = initp%soil_energy(nzg) + energy_room * dslzi8(nzg)
            !------------------------------------------------------------------------------!



            !----- Compute runoff for output. ---------------------------------------------!
            if (fast_diagnostics) then
               !---------------------------------------------------------------------------!
               !      There is no need to divide wfreeb and qwfree  by time step,          !
               ! which will be done in subroutine normalize_averaged_vars.                 !
               !---------------------------------------------------------------------------!
               csite%runoff       (ipa) = csite%runoff(ipa)                                &
                                        + sngloff(wmass_runoff * dtrk4i,tiny_offset)
               csite%fmean_runoff (ipa) = csite%fmean_runoff(ipa)                          &
                                        + sngloff(wmass_runoff         ,tiny_offset)
               csite%fmean_qrunoff(ipa) = csite%fmean_qrunoff(ipa)                         &
                                        + sngloff(energy_runoff        ,tiny_offset)
            end if
            if (checkbudget) then
               !---------------------------------------------------------------------------!
               !      To make sure that the previous values of wbudget_loss2runoff         !
               ! and ebudget_loss2runoff are accumulated to the next time step.            !
               !---------------------------------------------------------------------------!
               initp%wbudget_loss2runoff = initp%wbudget_loss2runoff + wmass_runoff
               initp%ebudget_loss2runoff = initp%ebudget_loss2runoff + energy_runoff
               initp%wbudget_storage     = initp%wbudget_storage     - wmass_runoff
               initp%ebudget_storage     = initp%ebudget_storage     - energy_runoff
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !----- All remaining water was taken care of, set free water to zero. ---------!
            wmass_free  = 0.d0
            energy_free = 0.d0
            depth_free  = 0.d0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      elseif (wmass_free < 0.d0) then
         !---------------------------------------------------------------------------------!
         !     Free water is negative.  We must take water from other sources to make up   !
         ! the difference.                                                                 !
         !---------------------------------------------------------------------------------!
         wmass_needed  = - wmass_free
         energy_needed = - energy_free
         depth_needed  = - depth_free
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Set free water to zero (as there is no excess, but demand for water).       !
         !---------------------------------------------------------------------------------!
         wmass_free  = 0.d0
         energy_free = 0.d0
         depth_free  = 0.d0
         !---------------------------------------------------------------------------------!


         !------ Find the amount of water that the soil can provide. ----------------------!
         wmass_available  = max(0.d0,initp%soil_water(nzg) - soil8(nsoil)%soilcp)          &
                          * wdns8 * dslz8(nzg)
         energy_available = energy_free * wmass_available / wmass_free
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Ideally, we take the water needed from the top soil layer.  Check whether   !
         ! the amount of available water is sufficient to clear the debt.                  !
         !---------------------------------------------------------------------------------!
         if (wmass_available >= wmass_needed) then
            !------------------------------------------------------------------------------!
            !     There is enough space in the top soil layer to correct remaining         !
            ! negative water, get all the water needed there.                              !
            !------------------------------------------------------------------------------!
            flag_second = 'C'
            !------------------------------------------------------------------------------!


            !----- Update soil water and internal energy. ---------------------------------!
            initp%soil_water(nzg)  = initp%soil_water(nzg)  - wmass_needed  * dslzi8(nzg)  &
                                   * wdnsi8
            initp%soil_energy(nzg) = initp%soil_energy(nzg) - energy_needed * dslzi8(nzg)
            !------------------------------------------------------------------------------!



            !----- All water demand was taken care of, set needed water to zero. ----------!
            wmass_needed           = 0.d0
            energy_needed          = 0.d0
            depth_needed           = 0.d0
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     There is not enough space in the top soil layer to correct remaining     !
            ! negative water, get all the water we can from the top soil and condense the  !
            ! remaining.                                                                   !
            !------------------------------------------------------------------------------!
            flag_second = 'D'
            !------------------------------------------------------------------------------!


            !----- Add the water that can come from the soil. -----------------------------!
            wmass_needed  = wmass_needed  - wmass_available
            energy_needed = energy_needed - energy_available
            !------------------------------------------------------------------------------!



            !----- Find the temperature associated with the water needed. -----------------!
            call uint2tl8(energy_needed/wmass_needed,tempk_needed,fracliq_needed)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Water vapour is at higher energetic state than ice or liquid water.  Find !
            ! the excess energy (enthalpy really) associated with condensation.            !
            !------------------------------------------------------------------------------!
            energy_latent = wmass_needed * tq2enthalpy8(tempk_needed,1.d0,.true.)          &
                          - energy_needed
            !------------------------------------------------------------------------------!


            !----- Dump what we can dump on the top soil layer. ---------------------------!
            initp%soil_water(nzg)  = initp%soil_water(nzg)                                 &
                                   - wmass_available  * dslzi8(nzg) * wdnsi8
            initp%soil_energy(nzg) = initp%soil_energy(nzg)                                &
                                   - ( energy_available - energy_latent ) * dslzi8(nzg)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Condense the remaining, and hope for the best.                           !
            !------------------------------------------------------------------------------!
            !----- Update the canopy air space properties. --------------------------------!
            initp%can_shv      = initp%can_shv - wmass_needed   * initp%wcapcani
            initp%can_enthalpy = initp%can_enthalpy                                        &
                               - ( energy_needed + energy_latent ) * initp%hcapcani
            !----- Update the fluxes. -----------------------------------------------------!
            initp%avg_vapor_gc = initp%avg_vapor_gc  - wmass_needed   * hdidi
            !------------------------------------------------------------------------------!




            !----- All water demand was taken care of, set needed water to zero. ----------!
            wmass_needed  = 0.d0
            energy_needed = 0.d0
            depth_needed  = 0.d0
            !------------------------------------------------------------------------------!

         end if
         !---------------------------------------------------------------------------------!
      else
         !----- No water transfer needed (wmass_free is exactly zero). --------------------!
         flag_second = 'E'
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check the total amount of mass in the temporary surface water/snow, and adjust !
      ! the number of layer accordingly.                                                   !
      !------------------------------------------------------------------------------------!
      if (sum_sfcw_mass <= rk4tiny_sfcw_mass) then
         !----- Not enough water in the temporary surface water, eliminate all layers. ----!
         initp%nlev_sfcwater = 0
         !---------------------------------------------------------------------------------!


         !----- Update the flag for temporary surface water. ------------------------------!
         initp%flag_sfcwater = 0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      The total mass should be either zero or greater than rk4tiny_sfcw_mass,    !
         ! but, just in case, we add any remaining energy to the top soil layer.           !
         !---------------------------------------------------------------------------------!
         initp%soil_water(nzg)  = initp%soil_water(nzg)  + sum_sfcw_mass   * dslzi8(nzg)   &
                                                         * wdnsi8
         initp%soil_energy(nzg) = initp%soil_energy(nzg) + sum_sfcw_energy * dslzi8(nzg)
         !---------------------------------------------------------------------------------!

         !----- Loop all layers and re-set all extensive variables to zero. ---------------!
         do k = 1, nzs
            initp%sfcwater_mass  (k)  = 0.d0
            initp%sfcwater_energy(k)  = 0.d0
            initp%sfcwater_depth (k)  = 0.d0
         end do
         !---------------------------------------------------------------------------------!

      elseif (sum_sfcw_mass < rk4water_stab_thresh) then


         !----- Not much water in the temporary surface water, impose a single layer. -----!
         initp%nlev_sfcwater = 1
         !---------------------------------------------------------------------------------!


         !----- Update the flag for temporary surface water. ------------------------------!
         initp%flag_sfcwater = 1
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    If the total amount of temporary surface water is not enough to make it      !
         ! stable, we impose it to have a single layer with all the ponding/snow in there. !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_mass  (1) = sum_sfcw_mass
         initp%sfcwater_energy(1) = sum_sfcw_energy
         initp%sfcwater_depth (1) = sum_sfcw_depth
         do k=2,nzs
            initp%sfcwater_mass   (k) = 0.d0
            initp%sfcwater_energy (k) = 0.d0
            initp%sfcwater_depth  (k) = 0.d0
         end do
         !---------------------------------------------------------------------------------!


      else
         !----- Update the flag for temporary surface water. ------------------------------!
         initp%flag_sfcwater = 2
         !---------------------------------------------------------------------------------!


         !---- Check whether there is enough snow for a new layer. ------------------------!
         nlayers   = ksnnew
         newlayers = 1
         do k = 1,ksnnew
            !------------------------------------------------------------------------------!
            !     Check whether the layer as is meet the minimum requirements to stand as  !
            ! a new layer by itself.                                                       !
            !------------------------------------------------------------------------------!
            if ( initp%sfcwater_mass(k)   >=  rk4snowmin                       .and.       &
                 rk4snowmin * fsdnsi8     <= initp%sfcwater_depth(k)           .and.       &
                 initp%sfcwater_energy(k) <  initp%sfcwater_mass(k) * uiliqt38       ) then
               newlayers = newlayers + 1
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!

         !----- Newlayers is the new number of temporary surface water layers. ------------!
         newlayers = min(newlayers, nzs, nlayers + 1)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Decide how many temporary surface water or snowpack layers to keep.         !
         !---------------------------------------------------------------------------------!
         if (newlayers == 1) then
            !----- Keep all temporary surface water in one layer. -------------------------!
            newsfcw_mass  (1) = sum_sfcw_mass
            newsfcw_energy(1) = sum_sfcw_energy
            newsfcw_depth (1) = sum_sfcw_depth
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     Snowpack.  Ensure that the layer thickness of every layer is reasonable, !
            ! and exchange energy and mass.                                                !
            !------------------------------------------------------------------------------!
            kold  = 1
            wtnew = 1.d0
            wtold = 1.d0
            do k = 1,newlayers
               !---- Find the new mass of this layer. -------------------------------------!
               newsfcw_mass  (k) = sum_sfcw_mass * thick(k,newlayers)
               newsfcw_energy(k) = 0.d0
               newsfcw_depth (k) = 0.d0
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the properties of this new layer.                                !
               !---------------------------------------------------------------------------!
               find_layer: do

                  !----- Difference between old and new snow ------------------------------!
                  wdiff = wtnew * newsfcw_mass(k) - wtold * initp%sfcwater_mass(kold)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Check how to redistribute mass, energy and depth.                   !
                  !------------------------------------------------------------------------!
                  if (wdiff > 0.d0) then
                     !---------------------------------------------------------------------!
                     !     Take wtold fraction of energy and depth from former "kold"      !
                     ! layer.                                                              !
                     !---------------------------------------------------------------------!
                     newsfcw_energy(k) = newsfcw_energy(k)                                 &
                                       + wtold * initp%sfcwater_energy(kold)
                     newsfcw_depth(k)  = newsfcw_depth(k)                                  &
                                       + wtold * initp%sfcwater_depth(kold)
                     !---------------------------------------------------------------------!


                     !----- Update weights. -----------------------------------------------!
                     wtnew  = wtnew - wtold * initp%sfcwater_mass(kold) / newsfcw_mass(k)
                     wtold  = 1.d0
                     !---------------------------------------------------------------------!

                     !----- Shift "kold" layer upwards. -----------------------------------!
                     kold   = kold + 1
                     if (kold > nlayers) exit find_layer
                     !---------------------------------------------------------------------!
                  else
                     !---------------------------------------------------------------------!
                     !      Energy and depth are going to be a mix of new and old layers.  !
                     !---------------------------------------------------------------------!
                     newsfcw_energy(k) = newsfcw_energy(k) + wtnew * newsfcw_mass(k)       &
                                       * initp%sfcwater_energy(kold)                       &
                                       / max(rk4tiny_sfcw_mass,initp%sfcwater_mass(kold))
                     newsfcw_depth(k)  = newsfcw_depth(k)  + wtnew * newsfcw_mass(k)       &
                                       * initp%sfcwater_depth(kold)                        &
                                       / max(rk4tiny_sfcw_mass,initp%sfcwater_mass(kold))
                     !---------------------------------------------------------------------!



                     !----- Update weights. -----------------------------------------------!
                     wtold = wtold - wtnew * newsfcw_mass(k)                               &
                                   / max(rk4tiny_sfcw_mass,initp%sfcwater_mass(kold))
                     wtnew = 1.d0
                     !---------------------------------------------------------------------!


                     !----- Move to the next new layer. -----------------------------------!
                     exit find_layer
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end do find_layer
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!

         !----- Update the water/snow layer prognostic properties. ------------------------!
         initp%nlev_sfcwater = newlayers
         do k = 1,newlayers
            initp%sfcwater_mass  (k) = newsfcw_mass  (k)
            initp%sfcwater_energy(k) = newsfcw_energy(k)
            initp%sfcwater_depth (k) = newsfcw_depth (k)
         end do
         do k = newlayers + 1, nzs
            initp%sfcwater_mass   (k) = 0.d0
            initp%sfcwater_energy (k) = 0.d0
            initp%sfcwater_depth  (k) = 0.d0
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Determine the total amount of temporary surface water available as well as    !
      ! derived properties.                                                                !
      !------------------------------------------------------------------------------------!
      sum_sfcw_mass   = 0.d0
      sum_sfcw_energy = 0.d0
      sum_sfcw_depth  = 0.d0
      do k=1,initp%nlev_sfcwater
         sum_sfcw_mass   = sum_sfcw_mass   + initp%sfcwater_mass  (k)
         sum_sfcw_energy = sum_sfcw_energy + initp%sfcwater_energy(k)
         sum_sfcw_depth  = sum_sfcw_depth  + initp%sfcwater_depth (k)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Compute the budget variables after the adjustments.                           !
      !------------------------------------------------------------------------------------!
      flag_sfcw_end      = initp%flag_sfcwater
      nlev_sfcw_end      = initp%nlev_sfcwater
      wmass_cas_end      = initp%can_shv      * initp%wcapcan
      enthalpy_cas_end   = initp%can_enthalpy * initp%hcapcan
      wmass_virtual_end  = initp%virtual_water
      energy_virtual_end = initp%virtual_energy
      wmass_sfcw_end     = sum_sfcw_mass
      energy_sfcw_end    = sum_sfcw_energy
      wmass_soil_end     = initp%soil_water(nzg)  * dslz8(nzg) * wdns8
      energy_soil_end    = initp%soil_energy(nzg) * dslz8(nzg)
      wmass_total_end    = wmass_virtual_end  + wmass_sfcw_end  + wmass_soil_end           &
                         + wmass_cas_end + wmass_runoff
      energy_total_end   = energy_virtual_end + energy_sfcw_end + energy_soil_end          &
                         + enthalpy_cas_end + energy_runoff
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Check whether energy and mass are conserved.                                  !
      !------------------------------------------------------------------------------------!
      wmass_total_rch  = 2.d0 * abs(wmass_total_end  - wmass_total_beg)                    &
                       / (abs(wmass_total_end ) + abs(wmass_total_beg ))
      energy_total_rch = 2.d0 * abs(energy_total_end - energy_total_beg)                   &
                       / (abs(energy_total_end) + abs(energy_total_beg))
      if (wmass_total_rch > sfcw_toler .or. energy_total_rch > sfcw_toler) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') '    Water or energy conservation was violated!!!'
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt=fmtt ) ' Time                    = ', current_time%year         &
                                                               , current_time%month        &
                                                               , current_time%date         &
                                                               , current_time%time
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' - Top soil layer: '
         write (unit=*,fmt=fmtf ) '   + Depth               = ',dslz8 (nzg)
         write (unit=*,fmt=fmtf ) '   + 1/Depth             = ',dslzi8(nzg)
         write (unit=*,fmt=fmtf ) '   + Density             = ',wdns8
         write (unit=*,fmt=fmtf ) '   + 1/Density           = ',wdnsi8
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' - Initial conditions: '
         write (unit=*,fmt=fmtf ) '   + Total water mass    = ',wmass_total_beg
         write (unit=*,fmt=fmtf ) '   + CAS mass            = ',wmass_cas_beg
         write (unit=*,fmt=fmtf ) '   + Virtual mass        = ',wmass_virtual_beg
         write (unit=*,fmt=fmtf ) '   + Ponding/snow mass   = ',wmass_sfcw_beg
         write (unit=*,fmt=fmtf ) '   + Soil mass           = ',wmass_soil_beg
         write (unit=*,fmt=fmtf ) '   + Total energy        = ',energy_total_beg
         write (unit=*,fmt=fmtf ) '   + CAS enthalpy        = ',enthalpy_cas_beg
         write (unit=*,fmt=fmtf ) '   + Virtual energy      = ',energy_virtual_beg
         write (unit=*,fmt=fmtf ) '   + Ponding/snow energy = ',energy_sfcw_beg
         write (unit=*,fmt=fmtf ) '   + Soil energy         = ',energy_soil_beg
         write (unit=*,fmt=fmti ) '   + # Sfc. water layers = ',nlev_sfcw_beg
         write (unit=*,fmt=fmti ) '   + Sfc. water flag     = ',flag_sfcw_beg
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' - Final conditions: '
         write (unit=*,fmt=fmtf ) '   + Total water mass    = ',wmass_total_end
         write (unit=*,fmt=fmtf ) '   + CAS mass            = ',wmass_cas_end
         write (unit=*,fmt=fmtf ) '   + Virtual mass        = ',wmass_virtual_end
         write (unit=*,fmt=fmtf ) '   + Ponding/snow mass   = ',wmass_sfcw_end
         write (unit=*,fmt=fmtf ) '   + Soil mass           = ',wmass_soil_end
         write (unit=*,fmt=fmtf ) '   + Runoff mass         = ',wmass_runoff
         write (unit=*,fmt=fmtf ) '   + Total energy        = ',energy_total_end
         write (unit=*,fmt=fmtf ) '   + CAS enthalpy        = ',enthalpy_cas_end
         write (unit=*,fmt=fmtf ) '   + Virtual energy      = ',energy_virtual_end
         write (unit=*,fmt=fmtf ) '   + Ponding/snow energy = ',energy_sfcw_end
         write (unit=*,fmt=fmtf ) '   + Soil energy         = ',energy_soil_end
         write (unit=*,fmt=fmtf ) '   + Runoff energy       = ',energy_runoff
         write (unit=*,fmt=fmti ) '   + # Sfc. water layers = ',nlev_sfcw_end
         write (unit=*,fmt=fmti ) '   + Sfc. water flag     = ',flag_sfcw_end
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' - Exchange cases: '
         write (unit=*,fmt=fmtc ) '   + First exchange      = ',flag_first
         write (unit=*,fmt=fmtc ) '   + Second exchange     = ',flag_second
         write (unit=*,fmt=fmtc ) '   + Stable flag         = ',flag_stable
         write (unit=*,fmt=fmtc ) '   + Percolation flag    = ',flag_percol
         write (unit=*,fmt='(a)') ' '
          write (unit=*,fmt='(a)') ' - Unresolved water: '
         write (unit=*,fmt=fmtf ) '   + Free mass           = ',wmass_free
         write (unit=*,fmt=fmtf ) '   + Free energy         = ',energy_free
         write (unit=*,fmt=fmtf ) '   + Needed mass         = ',wmass_needed
         write (unit=*,fmt=fmtf ) '   + Needed energy       = ',energy_needed
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' - Relative error: '
         write (unit=*,fmt=fmtf ) '   + Tolerance           = ',sfcw_toler
         write (unit=*,fmt=fmtf ) '   + Total water mass    = ',wmass_total_rch
         write (unit=*,fmt=fmtf ) '   + Total energy        = ',energy_total_rch
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         call fatal_error('Energy or water is not being conserved!!!'                      &
                         ,'adjust_sfcw_properties','rk4_misc.f90')
      end if

      return
   end subroutine adjust_sfcw_properties
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will ensure that top soil properties are within the allowed       !
   ! range.  We currently test this only for the top soil layer because this is the most   !
   ! likely to cause problems, but if it does happen in other layers, we could easily      !
   ! extend this for all layers.  Depending on its derivative, soil moisture can go under  !
   ! the minimum soil moisture possible (soilcp) or above the bubbling point (soilbp).     !
   ! Both are bad things, and if the value is way off-bounds, then we leave it like that   !
   ! so the step can be rejected.   However, if the value is just slightly off(*) these    !
   ! limits, we make a small exchange of moisture with the neighbouring layer.  This will  !
   ! prevent the soil to go outside the range in those double precision => single          !
   ! precision => double precision conversion.                                             !
   !                                                                                       !
   ! (*) slightly off is defined as outside the range but within the desired accuracy      !
   !     (rk4eps).                                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine adjust_topsoil_properties(initp,hdid)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4site              & ! intent(in)
                                      , rk4eps               & ! intent(in)
                                      , rk4tiny_sfcw_mass    & ! intent(in)
                                      , rk4min_can_shv       ! ! intent(in)
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            ! ! structure
      use consts_coms          , only : t3ple8               & ! intent(in)
                                      , wdns8                & ! intent(in)
                                      , fdnsi8               & ! intent(in)
                                      , toodry8              & ! intent(in)
                                      , wdnsi8               ! ! intent(in)
      use therm_lib8           , only : uextcm2tl8           & ! subroutine
                                      , uint2tl8             & ! subroutine
                                      , tl2uint8             & ! function
                                      , tq2enthalpy8         ! ! function
      use grid_coms            , only : nzg                  ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , dslzi8               & ! intent(in)
                                      , dslz8                ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype)     , target     :: initp  ! Integration buffer
      real(kind=8)           , intent(in) :: hdid   ! Time step 
      !----- Local variables --------------------------------------------------------------!
      integer                             :: kt
      integer                             :: kb
      integer                             :: kw
      integer                             :: nstop
      integer                             :: nsbeneath
      real(kind=8)                        :: hdidi
      real(kind=8)                        :: shctop
      real(kind=8)                        :: shcbeneath
      real(kind=8)                        :: water_needed
      real(kind=8)                        :: energy_needed
      real(kind=8)                        :: depth_needed
      real(kind=8)                        :: energy_excess
      real(kind=8)                        :: water_excess
      real(kind=8)                        :: depth_excess
      real(kind=8)                        :: water_available
      real(kind=8)                        :: energy_available
      real(kind=8)                        :: water_room
      real(kind=8)                        :: energy_room
      real(kind=8)                        :: virtual_tempk
      real(kind=8)                        :: virtual_fracliq
      real(kind=8)                        :: virtual_latent
      logical                             :: slightlymoist
      logical                             :: slightlydry
      !------------------------------------------------------------------------------------!

      !----- Inverse of time increment ----------------------------------------------------!
      hdidi = 1.d0 / hdid
      
      !----- Defining some aliases that will be often used during the integration. --------!
      kt            = nzg
      nstop         = rk4site%ntext_soil(kt)

      !----- Check whether we are just slightly off. --------------------------------------!
      slightlymoist = initp%soil_water(kt) > soil8(nstop)%soilbp 
      slightlydry   = initp%soil_water(kt) < soil8(nstop)%soilcp

      !------------------------------------------------------------------------------------!
      !     If we reached this point, then we may be slightly off-track.  It is very       !
      ! likely that we will need top soil layer temperature and liquid fraction, so we     !
      ! find them now...                                                                   !
      !------------------------------------------------------------------------------------!
      shctop     = soil8(nstop)%slcpd
      call uextcm2tl8(initp%soil_energy(kt),initp%soil_water(kt)*wdns8,shctop              &
                     ,initp%soil_tempk(kt),initp%soil_fracliq(kt))


      !------------------------------------------------------------------------------------!
      !      If we access this IF statement, then the soil is slightly dry.  Here we will  !
      ! look for water in adjacent environments and "steal" the amount of water that the   !
      ! top soil layer needs to be exactly at soilcp.                                      !
      !------------------------------------------------------------------------------------!
      if (slightlydry) then

         !---------------------------------------------------------------------------------!
         !     Now we find how much water we need.  Since we will need to exchange with    !
         ! other environments, find it in kg/m2, the standard units.                       !
         !---------------------------------------------------------------------------------!
         water_needed  = (soil8(nstop)%soilcp - initp%soil_water(kt)) * dslz8(kt) * wdns8

         !---------------------------------------------------------------------------------!
         !    First, we check whether there is a temporary surface water.  If so, then we  !
         ! "steal" water from this layer, even if it is not enough.                        !
         !---------------------------------------------------------------------------------!
         if (initp%nlev_sfcwater > 0) then
            sfcwsrc: do
               water_available = initp%sfcwater_mass(1)
               if (water_available > water_needed) then
                  !------------------------------------------------------------------------!
                  !    The surface layer has enough water to make up the difference, find  !
                  ! the energy and depth that will be removed from this layer.             !
                  !------------------------------------------------------------------------!
                  energy_needed = water_needed * initp%sfcwater_energy(1) / water_available
                  depth_needed  = water_needed * initp%sfcwater_depth(1)  / water_available

                  !------------------------------------------------------------------------!
                  !     Add the water and energy into the top layer, remove it from the    !
                  ! surface, and quit.                                                     !
                  !------------------------------------------------------------------------!
                  initp%soil_water(kt)     = initp%soil_water(kt)                          &
                                           + water_needed * dslzi8(kt) * wdnsi8
                  initp%soil_energy(kt)    = initp%soil_energy(kt)                         &
                                           + energy_needed * dslzi8(kt)

                  initp%sfcwater_depth(1)  = initp%sfcwater_depth(1)  - depth_needed
                  initp%sfcwater_mass(1)   = initp%sfcwater_mass(1)   - water_needed
                  initp%sfcwater_energy(1) = initp%sfcwater_energy(1) - energy_needed
                  return

               elseif (water_available > 0.d0) then
                  !------------------------------------------------------------------------!
                  !     This water will not be enough, but we use all this water, elimi-   !
                  ! nate the layer,  and seek for other sources to fill the remainder.     !
                  !------------------------------------------------------------------------!
                  water_needed = water_needed - water_available
                  initp%soil_water(kt)     = initp%soil_water(kt)                          &
                                           + initp%sfcwater_mass(1) * dslzi8(kt) * wdnsi8
                  initp%soil_energy(kt)    = initp%soil_energy(kt)                         &
                                           + initp%sfcwater_energy(1) * dslzi8(kt)
                  !----- If there were more layers, then we move them down. ---------------!
                  do kw=1,initp%nlev_sfcwater-1
                     initp%sfcwater_mass  (kw) = initp%sfcwater_mass  (kw+1)
                     initp%sfcwater_energy(kw) = initp%sfcwater_energy(kw+1)
                     initp%sfcwater_depth (kw) = initp%sfcwater_depth (kw+1)
                  end do
                  
                  !----- Remove the top layer. --------------------------------------------!
                  kw=initp%nlev_sfcwater
                  initp%sfcwater_mass  (kw) = 0.d0
                  initp%sfcwater_energy(kw) = 0.d0
                  initp%sfcwater_depth(kw)  = 0.d0
                  initp%nlev_sfcwater       = initp%nlev_sfcwater - 1
                  !----- If no surface water layer is left, look for another source... ----!
                  if (initp%nlev_sfcwater == 0) exit sfcwsrc
               end if
            end do sfcwsrc
         end if

         !---------------------------------------------------------------------------------!
         !     If we hit this point, we aren't done yet and we no longer have temporary    !
         ! surface water.  The next candidate is the virtual layer.                        !
         !---------------------------------------------------------------------------------!
         water_available = initp%virtual_water
         if (water_available > water_needed) then
            !------------------------------------------------------------------------------!
            !    Virtual layer has enough water to solve the problem.  Find the energy     !
            ! associated with the water that will be transferred to the top layer, and     !
            ! move it to there too.                                                        !
            !------------------------------------------------------------------------------!
            energy_needed = water_needed * initp%virtual_energy / water_available
            depth_needed  = water_needed * initp%virtual_depth  / water_available

            !------------------------------------------------------------------------------!
            !     Add the water and energy into the top layer, remove it from the surface, !
            ! and quit.                                                                    !
            !------------------------------------------------------------------------------!
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  + water_needed  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)

            initp%virtual_depth   = initp%virtual_depth  - depth_needed
            initp%virtual_water   = initp%virtual_water  - water_needed
            initp%virtual_energy  = initp%virtual_energy - energy_needed
            return
         elseif (water_available > 0.d0) then
            !------------------------------------------------------------------------------!
            !     This water will not be enough, but we use it and seek for other sources  !
            ! to fill the remainder.                                                       !
            !------------------------------------------------------------------------------!
            water_needed = water_needed - water_available
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  + initp%virtual_water * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt)                                  &
                                  + initp%virtual_energy * dslzi8(kt)

            !----- Remove the virtual layer. ----------------------------------------------!
            initp%virtual_water   = 0.d0
            initp%virtual_energy  = 0.d0
            initp%virtual_depth   = 0.d0
         end if

         !---------------------------------------------------------------------------------!
         !    If we hit this point the temporary layers were not enough.  We now look for  !
         ! water in the layer immediately beneath the top soil layer.  We now find the     !
         ! amount of water this layer can donate.  We will also need energy associated     !
         ! with the water that is donate, so we compute the temperature and liquid         !
         ! fraction.                                                                       !
         !---------------------------------------------------------------------------------!
         kb              = nzg -1 
         nsbeneath       = rk4site%ntext_soil(kb)
         water_available = (initp%soil_water(kb)-soil8(nsbeneath)%soilcp)                  &
                         * wdns8 * dslz8(kb)
         shcbeneath      = soil8(nsbeneath)%slcpd
         call uextcm2tl8(initp%soil_energy(kb),initp%soil_water(kb)*wdns8,shcbeneath       &
                        ,initp%soil_tempk(kb),initp%soil_fracliq(kb))
         !---------------------------------------------------------------------------------!



         if (water_available > water_needed) then
            !------------------------------------------------------------------------------!
            !    The layer beneath has enough water.  Find the energy associated with this !
            ! water.                                                                       !
            !------------------------------------------------------------------------------!
            energy_needed = water_needed                                                   &
                          * tl2uint8(initp%soil_tempk(kb),initp%soil_fracliq(kb))
            !------------------------------------------------------------------------------!


            !----- Update water and energy in both layers. --------------------------------!
            initp%soil_water (kt) = initp%soil_water (kt)                                  &
                                  + water_needed  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)
            initp%soil_water (kb) = initp%soil_water (kb)                                  &
                                  - water_needed  * dslzi8(kb) * wdnsi8
            initp%soil_energy(kb) = initp%soil_energy(kb) - energy_needed * dslzi8(kb)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    We must also update the soil fluxes.                                      !
            !------------------------------------------------------------------------------!
            initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) + water_needed * hdidi
            initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) - water_needed * hdidi

            !------ Leave the subroutine. -------------------------------------------------!
            return

         elseif (water_available > 0.d0) then
            !------------------------------------------------------------------------------!
            !    The water in that layer will not be enough, but we extract all that we    !
            ! can to reduce the amount we still need.                                      !
            !------------------------------------------------------------------------------!
            water_needed = water_needed - water_available
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find the energy associated with the water that will be transferred.       !
            !------------------------------------------------------------------------------!
            energy_available      = water_available                                        &
                                  * tl2uint8(initp%soil_tempk(kb),initp%soil_fracliq(kb))
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  + water_available  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt) + energy_available * dslzi8(kt)

            initp%soil_water(kb)  = initp%soil_water(kb)                                   &
                                  - water_available  * dslzi8(kb) * wdnsi8
            initp%soil_energy(kb) = initp%soil_energy(kb) - energy_available * dslzi8(kb)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    We must also update the soil fluxes.                                      !
            !------------------------------------------------------------------------------!
            initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) + water_available * hdidi
            initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) - water_available * hdidi
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     If we hit this point, our final hope is to trap some water from the canopy  !
         ! air space.  The water available is not everything above the minimum mixing      !
         ! ratio because we don't want to risk having the canopy air space crashing.       !
         !---------------------------------------------------------------------------------!
         water_available = initp%wcapcan * (initp%can_shv - 5.d0 * rk4min_can_shv)
         if (water_available > water_needed) then
            !------------------------------------------------------------------------------!
            !    There is enough water vapour. The transfer will require phase change, so  !
            ! the energy transfer will contain a latent heat flux of vaporisation.         !
            !------------------------------------------------------------------------------!
            energy_needed = water_needed * tq2enthalpy8(initp%soil_tempk(kt),1.d0,.true.)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Add the water and energy into the top layer, remove it from the canopy   !
            ! air space and quit.                                                          !
            !------------------------------------------------------------------------------!
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  + water_needed  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt) + energy_needed * dslzi8(kt)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Remove mass and energy from the canopy air space.                        !
            !------------------------------------------------------------------------------!
            initp%can_shv      = initp%can_shv                                             &
                               - water_needed  * initp%wcapcani
            initp%can_enthalpy = initp%can_enthalpy                                        &
                               - energy_needed * initp%hcapcani
            initp%avg_vapor_gc = initp%avg_vapor_gc                                        &
                               - water_needed  * hdidi
            !------------------------------------------------------------------------------!

            return

         elseif (water_available > 0.d0) then
            !------------------------------------------------------------------------------!
            !     This is a critically dry situation and the only reason we won't start    !
            ! crying is because we don't have enough water to waste in tears... As the     !
            ! final desperate act, we will trap any water in excess of the minimum         !
            ! moisture from the canopy air space.  Even if this is a tiny amount, it may   !
            ! be enough to just avoid the model to crash at the sanity check.              !
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the enthalpy associated with the partial condensation.              !
            !------------------------------------------------------------------------------!
            energy_available = water_available                                             &
                             * tq2enthalpy8(initp%soil_tempk(kt),1.d0,.true.)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Add the water and energy into the top layer, remove it from the canopy   !
            ! air space and quit.                                                          !
            !------------------------------------------------------------------------------!
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  + water_available  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt) + energy_available * dslzi8(kt)

            initp%can_shv         = initp%can_shv                                          &
                                  - water_available  * initp%wcapcani
            initp%can_enthalpy    = initp%can_enthalpy                                     &
                                  - energy_available * initp%hcapcani
            initp%avg_vapor_gc    = initp%avg_vapor_gc                                     &
                                  - water_available  * hdidi
            !------------------------------------------------------------------------------!

            return
         end if

      !------------------------------------------------------------------------------------!
      !      If we access this ELSEIF part, then the soil is slightly above saturation.    !
      !------------------------------------------------------------------------------------!
      elseif (slightlymoist) then
         !---------------------------------------------------------------------------------!
         !     Now we find how much water the top soil layer needs to withdraw.  Since we  !
         ! will need to exchange with other environments, find it in kg/m2, the standard   !
         ! units.  Also, find the total energy that must go away with the water.           !
         !---------------------------------------------------------------------------------!
         water_excess  = (initp%soil_water(kt) - soil8(nstop)%soilbp) * dslz8(kt) * wdns8
         energy_excess = water_excess                                                      &
                       * tl2uint8(initp%soil_tempk(kt),initp%soil_fracliq(kt))
         depth_excess  = water_excess * ( initp%soil_fracliq(kt) * wdnsi8                  &
                                        + (1.d0 - initp%soil_fracliq(kt)) * fdnsi8)
         !---------------------------------------------------------------------------------!

         if (initp%nlev_sfcwater > 0) then
            !------------------------------------------------------------------------------!
            !    If there is already a temporary surface water layer, we simply dump the   !
            ! excess of water in the first layer.                                          !
            !------------------------------------------------------------------------------!
            initp%soil_water(kt)     = initp%soil_water(kt)                                &
                                     - water_excess * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt)    = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
            
            initp%sfcwater_mass(1)   = initp%sfcwater_mass(1)   + water_excess
            initp%sfcwater_energy(1) = initp%sfcwater_energy(1) + energy_excess
            initp%sfcwater_depth(1)  = initp%sfcwater_depth(1)  + depth_excess
            !------------------------------------------------------------------------------!

            return
         elseif (initp%virtual_water + water_excess > rk4tiny_sfcw_mass) then
            !------------------------------------------------------------------------------!
            !     If the virtual layer will have some significant mass after adding the    !
            ! water excess, we transfer the water to there.  It will likely become the new !
            ! temporary surface water layer.                                               !
            !------------------------------------------------------------------------------!
            initp%soil_water(kt)     = initp%soil_water(kt)                                &
                                     - water_excess * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt)    = initp%soil_energy(kt) - energy_excess * dslzi8(kt)
            
            initp%virtual_water      = initp%virtual_water     + water_excess
            initp%virtual_energy     = initp%virtual_energy    + energy_excess
            initp%virtual_depth      = initp%virtual_depth     + depth_excess
            !------------------------------------------------------------------------------!

            return
         end if

         !---------------------------------------------------------------------------------!
         !    If we hit this point, then the amount of water is small but it can't go to   !
         ! the preferred destinations.  Adding on virtual pool wouldn't help because the   !
         ! water would be sent back to the top soil layer in adjust_sfcw_properties call.  !
         ! Thus the excess now becomes the excess of water in the layer, plus any water    !
         ! left in the virtual layer...                                                    !
         !    Anyway, we first eliminate any water left in the virtual layer by boiling it !
         ! to the atmosphere.  This is a tiny amount and even if supersaturation occurs,   !
         ! it shouldn't be enough to cause trouble.                                        !
         !---------------------------------------------------------------------------------!
         if (initp%virtual_water > rk4eps * rk4eps * rk4tiny_sfcw_mass) then 

            !------------------------------------------------------------------------------!
            !      Find the associated temperature of the virtual water.  The remaining    !
            ! water will be boiled.  The boiling will eliminate the virtual layer, so the  !
            ! latent heat must be taken from somewhere else (top soil layer in this case). !
            !------------------------------------------------------------------------------!
            call uint2tl8(initp%virtual_energy/initp%virtual_water                         &
                         ,virtual_tempk,virtual_fracliq)
            virtual_latent = initp%virtual_water * tq2enthalpy8(virtual_tempk,1.d0,.true.) &
                           - initp%virtual_energy
            !------------------------------------------------------------------------------!


            !----- Correct the canopy air space and fluxes. -------------------------------!
            initp%can_shv          = initp%can_shv                                         &
                                   + initp%virtual_water  * initp%wcapcani
            initp%can_enthalpy     = initp%can_enthalpy                                    &   
                                   + ( initp%virtual_energy + virtual_latent )             &
                                   * initp%hcapcani
            initp%soil_energy(nzg) = initp%soil_energy(nzg)                                &
                                   - virtual_latent * dslz8(nzg)
            initp%avg_vapor_gc     = initp%avg_vapor_gc + initp%virtual_water  * hdidi
            !------------------------------------------------------------------------------!

            
            !----- Say goodbye to the virtual layer... ------------------------------------!
            initp%virtual_energy = 0.d0
            initp%virtual_water  = 0.d0
            initp%virtual_depth  = 0.d0
            !------------------------------------------------------------------------------!
         elseif (initp%virtual_water > 0.d0) then
            !------------------------------------------------------------------------------!
            !    The amount of water is so small that round-off errors are bound to become !
            ! more important than the error of not conserving energy and water.  Simply    !
            ! extinguish the virtual layer.                                                !
            !------------------------------------------------------------------------------!
            initp%virtual_energy = 0.d0
            initp%virtual_water  = 0.d0
            initp%virtual_depth  = 0.d0
         endif 

         !---------------------------------------------------------------------------------!
         !     Back to the top soil layer, we still need to decide where we should send    !
         ! the excess...  First we try to send to the layer beneath.  This transfer of     !
         ! water implies that some energy is also transferred. Since soil layers are not   !
         ! pure water, we must find the actual amount of energy associated with the water  !
         ! transfer.  So first we compute the temperature and liquid fraction of the layer !
         ! beneath the top.                                                                !
         !---------------------------------------------------------------------------------!
         kb         = nzg-1
         nsbeneath  = rk4site%ntext_soil(kb)
         shcbeneath = soil8(nsbeneath)%slcpd
         call uextcm2tl8(initp%soil_energy(kb),initp%soil_water(kb)*wdns8,shcbeneath       &
                        ,initp%soil_tempk(kb),initp%soil_fracliq(kb))
         water_room = (soil8(nsbeneath)%soilbp-initp%soil_water(kb)) * wdns8 * dslz8(kb)

         if (water_room > water_excess) then
            !------------------------------------------------------------------------------!
            !    The layer beneath still has some room for this water excess, send the     !
            ! water down one level.                                                        !
            !------------------------------------------------------------------------------!
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  - water_excess  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt)                                  &
                                  - energy_excess * dslzi8(kt)
            initp%soil_water(kb)  = initp%soil_water(kb)                                   &
                                  + water_excess  * dslzi8(kb) * wdnsi8
            initp%soil_energy(kb) = initp%soil_energy(kb)                                  &
                                  + energy_excess * dslzi8(kb)
            !----- Update the fluxes too... -----------------------------------------------!
            initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) - water_excess * hdidi
            initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) + water_excess * hdidi
            !------------------------------------------------------------------------------!

            return
         elseif (water_room > 0.d0) then
            !------------------------------------------------------------------------------!
            !   The layer beneath the top layer can't take all water, but we send all that !
            ! we can to that layer so we reduce the amount that will be "boiled".  The     !
            ! remaining will go to the canopy air space.  Even if some supersaturation     !
            ! happens, the excess won't be too much to cause the run to crash.             !
            !------------------------------------------------------------------------------!
            water_excess = water_excess - water_room
            energy_room  = water_room                                                      &
                          * tl2uint8(initp%soil_tempk(kt),initp%soil_fracliq(kt))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    The layer beneath still has some room for this water excess, send the     !
            ! water down one level.                                                        !
            !------------------------------------------------------------------------------!
            initp%soil_water (kt) = initp%soil_water (kt)                                  &
                                  - water_room  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt)                                  &
                                  - energy_room * dslzi8(kt)
            initp%soil_water (kb) = initp%soil_water (kb)                                  &
                                  + water_room  * dslzi8(kb) * wdnsi8
            initp%soil_energy(kb) = initp%soil_energy(kb)                                  &
                                  + energy_room * dslzi8(kb)
            !----- Update the fluxes too... -----------------------------------------------!
            initp%avg_smoist_gg(kt) = initp%avg_smoist_gg(kt) - water_room * hdidi
            initp%avg_smoist_gg(kb) = initp%avg_smoist_gg(kb) + water_room * hdidi
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     The water that is about to leave will do it through "boiling".  Find the !
            ! latent heat associated with this phase change.                               !
            !------------------------------------------------------------------------------!
            energy_excess = water_excess * tq2enthalpy8(initp%soil_tempk(kt),1.d0,.true.)
            !------------------------------------------------------------------------------!



            !----- Send the water and energy to the canopy. -------------------------------!
            initp%soil_water(kt)  = initp%soil_water(kt)                                   &
                                  - water_excess  * dslzi8(kt) * wdnsi8
            initp%soil_energy(kt) = initp%soil_energy(kt)                                  &
                                  - energy_excess * dslzi8(kt)
            initp%can_shv         = initp%can_shv                                          &
                                  + water_excess  * initp%wcapcani
            initp%can_enthalpy    = initp%can_enthalpy                                     &
                                  + energy_excess * initp%hcapcani
            initp%avg_vapor_gc    = initp%avg_vapor_gc                                     &
                                  + water_excess  * hdidi
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine adjust_topsoil_properties
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will ensure that leaf water is positively defined.  Depending on   !
   ! its derivative, it can go under zero, in which case we must correct the derivatives   !
   ! rather than forcing it to be zero.  This guarantees mass conservation.  Likewise, if  !
   ! in the end of the step the leaf water is over the maximum, we remove the excess       !
   ! through shedding.  After this is checked, we then update the remaining leaf proper-   !
   ! ties, namely the temperature and liquid water fraction.                               !
   !---------------------------------------------------------------------------------------!
   subroutine adjust_veg_properties(initp,hdid,csite,ipa,ibuff)
      use rk4_coms             , only : rk4patchtype       & ! structure
                                      , rk4site            & ! structure
                                      , rk4aux             & ! structure
                                      , rk4min_veg_lwater  & ! intent(in)
                                      , rk4leaf_drywhc     & ! intent(in)
                                      , rk4leaf_maxwhc     & ! intent(in)
                                      , print_detailed     ! ! intent(in)
      use grid_coms            , only : nzg                ! ! intent(in)
      use ed_state_vars        , only : sitetype           & ! structure
                                      , patchtype          ! ! structure
      use ed_misc_coms         , only : fast_diagnostics   ! ! intent(in)
      use consts_coms          , only : t3ple8             & ! intent(in)
                                      , wdns8              & ! intent(in)
                                      , wdnsi8             & ! intent(in)
                                      , fdnsi8             ! ! intent(in)
      use therm_lib8           , only : uextcm2tl8         & ! subroutine
                                      , tl2uint8           & ! function
                                      , tq2enthalpy8       ! ! function
      use soil_coms            , only : soil8              & ! intent(in)
                                      , dslzi8             & ! intent(in)
                                      , dslz8              ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype)     , target     :: initp  ! Integration buffer
      type(sitetype)         , target     :: csite  ! Current site
      integer                , intent(in) :: ipa    ! Current patch ID
      real(kind=8)           , intent(in) :: hdid   ! Time step 
      integer                , intent(in) :: ibuff  ! Multithread ID
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer    :: cpatch
      integer                             :: ico
      integer                             :: ksn
      integer                             :: kt
      integer                             :: nstop
      real(kind=8)                        :: shctop
      real(kind=8)                        :: rk4min_leaf_water
      real(kind=8)                        :: rk4min_wood_water
      real(kind=8)                        :: rk4min_leaf_water_im2
      real(kind=8)                        :: rk4min_wood_water_im2
      real(kind=8)                        :: rk4max_leaf_water_im2
      real(kind=8)                        :: rk4max_wood_water_im2
      real(kind=8)                        :: min_leaf_water
      real(kind=8)                        :: max_leaf_water
      real(kind=8)                        :: min_wood_water
      real(kind=8)                        :: max_wood_water
      real(kind=8)                        :: leaf_water_hint
      real(kind=8)                        :: leaf_water_uint
      real(kind=8)                        :: leaf_water_zint
      real(kind=8)                        :: leaf_excess
      real(kind=8)                        :: leaf_demand
      real(kind=8)                        :: leaf_qexcess
      real(kind=8)                        :: leaf_qdemand
      real(kind=8)                        :: leaf_wshed
      real(kind=8)                        :: leaf_qwshed
      real(kind=8)                        :: leaf_dwshed
      real(kind=8)                        :: leaf_dew
      real(kind=8)                        :: leaf_qdew
      real(kind=8)                        :: leaf_boil
      real(kind=8)                        :: leaf_qboil
      real(kind=8)                        :: leaf_wshed_tot
      real(kind=8)                        :: leaf_qwshed_tot
      real(kind=8)                        :: leaf_dwshed_tot
      real(kind=8)                        :: leaf_dew_tot
      real(kind=8)                        :: leaf_qdew_tot
      real(kind=8)                        :: leaf_boil_tot
      real(kind=8)                        :: leaf_qboil_tot
      real(kind=8)                        :: wood_water_hint
      real(kind=8)                        :: wood_water_uint
      real(kind=8)                        :: wood_water_zint
      real(kind=8)                        :: wood_excess
      real(kind=8)                        :: wood_qexcess
      real(kind=8)                        :: wood_demand
      real(kind=8)                        :: wood_qdemand
      real(kind=8)                        :: wood_wshed
      real(kind=8)                        :: wood_qwshed
      real(kind=8)                        :: wood_dwshed
      real(kind=8)                        :: wood_dew
      real(kind=8)                        :: wood_qdew
      real(kind=8)                        :: wood_boil
      real(kind=8)                        :: wood_qboil
      real(kind=8)                        :: wood_wshed_tot
      real(kind=8)                        :: wood_qwshed_tot
      real(kind=8)                        :: wood_dwshed_tot
      real(kind=8)                        :: wood_dew_tot
      real(kind=8)                        :: wood_qdew_tot
      real(kind=8)                        :: wood_boil_tot
      real(kind=8)                        :: wood_qboil_tot
      real(kind=8)                        :: soil_demand
      real(kind=8)                        :: soil_qdemand
      real(kind=8)                        :: soil_excess
      real(kind=8)                        :: soil_qexcess
      real(kind=8)                        :: soil_water_uint
      real(kind=8)                        :: old_leaf_energy
      real(kind=8)                        :: old_leaf_water
      real(kind=8)                        :: old_leaf_water_im2
      real(kind=8)                        :: old_leaf_temp
      real(kind=8)                        :: old_leaf_fliq
      real(kind=8)                        :: old_wood_energy
      real(kind=8)                        :: old_wood_water
      real(kind=8)                        :: old_wood_water_im2
      real(kind=8)                        :: old_wood_temp
      real(kind=8)                        :: old_wood_fliq
      real(kind=8)                        :: hdidi
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      
      !----- Inverse of time increment ----------------------------------------------------!
      hdidi = 1.d0 / hdid

      !------------------------------------------------------------------------------------!
      !     If we reached this point, then we may be slightly off-track.  It is very       !
      ! likely that we will need top soil layer temperature and liquid fraction, so we     !
      ! find them now.                                                                     !
      !------------------------------------------------------------------------------------!
      kt              = nzg
      nstop           = rk4site%ntext_soil(kt)
      shctop          = soil8(nstop)%slcpd
      call uextcm2tl8(initp%soil_energy(kt),initp%soil_water(kt)*wdns8,shctop              &
                     ,initp%soil_tempk(kt),initp%soil_fracliq(kt))
      soil_water_uint = tl2uint8(initp%soil_tempk(kt),initp%soil_fracliq(kt))
      !------------------------------------------------------------------------------------!

      !----- Initialise the total shedding. -----------------------------------------------!
      leaf_wshed_tot  = 0.d0 
      leaf_qwshed_tot = 0.d0
      leaf_dwshed_tot = 0.d0
      leaf_dew_tot    = 0.d0 
      leaf_qdew_tot   = 0.d0
      leaf_boil_tot   = 0.d0 
      leaf_qboil_tot  = 0.d0
      wood_wshed_tot  = 0.d0 
      wood_qwshed_tot = 0.d0
      wood_dwshed_tot = 0.d0
      wood_dew_tot    = 0.d0 
      wood_qdew_tot   = 0.d0
      wood_boil_tot   = 0.d0 
      wood_qboil_tot  = 0.d0

      !----- Looping over cohorts ---------------------------------------------------------!
      cohortloop: do ico=1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !   Now we find the maximum leaf and wood water possible.                         !
         !---------------------------------------------------------------------------------!
         if (initp%leaf_resolvable(ico) .or. initp%wood_resolvable(ico)) then
            rk4min_leaf_water     = rk4min_veg_lwater * initp%lai(ico)
            rk4min_wood_water     = rk4min_veg_lwater * initp%wai(ico)
            rk4min_leaf_water_im2 = rk4aux(ibuff)%rk4min_leaf_water_im2(ico)
            rk4max_leaf_water_im2 = rk4aux(ibuff)%rk4max_leaf_water_im2(ico)
            rk4min_wood_water_im2 = rk4aux(ibuff)%rk4min_wood_water_im2(ico)
            rk4max_wood_water_im2 = rk4aux(ibuff)%rk4max_wood_water_im2(ico)
            min_leaf_water        = rk4leaf_drywhc    * initp%lai(ico)
            max_leaf_water        = rk4leaf_maxwhc    * initp%lai(ico)
            min_wood_water        = rk4leaf_drywhc    * initp%wai(ico)
            max_wood_water        = rk4leaf_maxwhc    * initp%wai(ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    In case water is to be removed or added, we will need to update the       !
            ! vegetation internal energy.  We want to preserve the temperature, though,    !
            ! because exchanges happen as mass flux or latent heat flux (fast dew or       !
            ! boiling).                                                                    !
            !------------------------------------------------------------------------------!
            !----- Leaf. ------------------------------------------------------------------!
            call uextcm2tl8( initp%leaf_energy   (ico)                                     &
                           , initp%leaf_water    (ico)                                     &
                           + initp%leaf_water_im2(ico)                                     &
                           , initp%leaf_hcap     (ico)                                     &
                           , initp%leaf_temp     (ico)                                     &
                           , initp%leaf_fliq     (ico) )
            leaf_water_hint    = tq2enthalpy8(initp%leaf_temp(ico),1.d0,.true.)
            leaf_water_uint    = tl2uint8(initp%leaf_temp(ico),initp%leaf_fliq(ico))
            leaf_water_zint    =          initp%leaf_fliq(ico)   * wdnsi8                  &
                               + ( 1.d0 - initp%leaf_fliq(ico) ) * fdnsi8
            old_leaf_energy    = initp%leaf_energy   (ico)
            old_leaf_water     = initp%leaf_water    (ico)
            old_leaf_water_im2 = initp%leaf_water_im2(ico)
            old_leaf_temp      = initp%leaf_temp     (ico)
            old_leaf_fliq      = initp%leaf_fliq     (ico)
            !----- Wood. ------------------------------------------------------------------!
            call uextcm2tl8( initp%wood_energy   (ico)                                     &
                           , initp%wood_water    (ico)                                     &
                           + initp%wood_water_im2(ico)                                     &
                           , initp%wood_hcap     (ico)                                     &
                           , initp%wood_temp     (ico)                                     &
                           , initp%wood_fliq     (ico) )
            wood_water_hint    = tq2enthalpy8(initp%wood_temp(ico),1.d0,.true.)
            wood_water_uint    = tl2uint8(initp%wood_temp(ico),initp%wood_fliq(ico))
            wood_water_zint    =          initp%wood_fliq(ico)   * wdnsi8                  &
                               + ( 1.d0 - initp%wood_fliq(ico) ) * fdnsi8
            old_wood_energy    = initp%wood_energy   (ico)
            old_wood_water     = initp%wood_water    (ico)
            old_wood_water_im2 = initp%wood_water_im2(ico)
            old_wood_temp      = initp%wood_temp     (ico)
            old_wood_fliq      = initp%wood_fliq     (ico)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Check whether we can solve leaves in this cohort...                          !
         !---------------------------------------------------------------------------------!
         if (initp%leaf_resolvable(ico)) then

            if (initp%leaf_water(ico) > max_leaf_water) then

               !---------------------------------------------------------------------------!
               !    Too much water over these leaves, we shall shed the excess to the      !
               ! ground.                                                                   !
               !---------------------------------------------------------------------------!
               leaf_wshed      = initp%leaf_water(ico) - max_leaf_water
               leaf_qwshed     = leaf_wshed * leaf_water_uint
               leaf_dwshed     = leaf_wshed * leaf_water_zint

               !----- Add the contribution of this cohort to the total shedding. ----------!
               leaf_wshed_tot  = leaf_wshed_tot  + leaf_wshed
               leaf_qwshed_tot = leaf_qwshed_tot + leaf_qwshed
               leaf_dwshed_tot = leaf_dwshed_tot + leaf_dwshed

               !----- Update water mass and energy. ---------------------------------------!
               initp%leaf_water (ico) = initp%leaf_water (ico) - leaf_wshed
               initp%veg_water  (ico) = initp%veg_water  (ico) - leaf_wshed
               initp%leaf_energy(ico) = initp%leaf_energy(ico) - leaf_qwshed
               initp%veg_energy (ico) = initp%veg_energy (ico) - leaf_qwshed

               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_wshed_lg(ico) = initp%avg_wshed_lg(ico) + leaf_wshed  ! * hdidi
               end if
               if (print_detailed) then
                  initp%cfx_qwshed  (ico) = initp%cfx_qwshed  (ico) + leaf_qwshed ! * hdidi
               end if
               !---------------------------------------------------------------------------!
    

            elseif (initp%leaf_water(ico) < min_leaf_water) then
               !---------------------------------------------------------------------------!
               !    If leaf_water is tiny and positive, exchange moisture with the air by  !
               ! donating the total amount as "boiling" (fast evaporation or sublimation). !
               ! In case the total is tiny but negative, exchange moisture with the air,   !
               ! "stealing" moisture as fast "dew/frost" condensation.                     !
               !---------------------------------------------------------------------------!
               leaf_boil  = max(0.d0,  initp%leaf_water(ico))
               leaf_dew   = max(0.d0,- initp%leaf_water(ico))
               leaf_qboil = leaf_boil * leaf_water_hint
               leaf_qdew  = leaf_dew  * leaf_water_hint
               !---------------------------------------------------------------------------!


               !----- Add the contribution of this cohort to the total boiling. -----------!
               leaf_boil_tot  = leaf_boil_tot  + leaf_boil
               leaf_dew_tot   = leaf_dew_tot   + leaf_dew
               leaf_qboil_tot = leaf_qboil_tot + leaf_qboil
               leaf_qdew_tot  = leaf_qdew_tot  + leaf_qdew
               !---------------------------------------------------------------------------!


               !----- Update cohort state variables. --------------------------------------!
               initp%leaf_water (ico)  = 0.d0
               initp%veg_water  (ico)  = initp%veg_water(ico)    + leaf_dew  - leaf_boil
               initp%leaf_energy(ico)  = initp%leaf_energy(ico)  + leaf_qdew - leaf_qboil
               initp%veg_energy (ico)  = initp%veg_energy(ico)   + leaf_qdew - leaf_qboil
               !---------------------------------------------------------------------------!


               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_vapor_lc(ico) = initp%avg_vapor_lc(ico)                        &
                                          + (leaf_boil  - leaf_dew ) ! * hdidi
               end if
               if (print_detailed) then
                  initp%cfx_qwflxlc (ico) = initp%cfx_qwflxlc(ico)                         &
                                          + (leaf_qboil - leaf_qdew) ! * hdidi
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !       Check whether leaf internal water is bounded.                          !
            !------------------------------------------------------------------------------!
            if (initp%leaf_water_im2(ico) > rk4max_leaf_water_im2) then
               !---------------------------------------------------------------------------!
               !     Leaves have too much water.  If possible, send water back to wood.    !
               ! If wood is also saturated, leaves expel the excess water as shedding.     !
               !---------------------------------------------------------------------------!
               leaf_excess     = rk4max_leaf_water_im2 - initp%leaf_water_im2(ico)
               wood_demand     = min(leaf_excess,max(0.d0                                  &
                                    ,rk4max_wood_water_im2-initp%wood_water_im2(ico)))
               leaf_wshed      = leaf_excess - wood_demand
               wood_qdemand    = wood_demand * leaf_water_uint
               leaf_qwshed     = leaf_wshed  * leaf_water_uint
               leaf_dwshed     = leaf_wshed  * leaf_water_zint
               leaf_qexcess    = wood_qdemand + leaf_qwshed
               !---------------------------------------------------------------------------!

               !----- Add the contribution of this cohort to the total shedding. ----------!
               leaf_wshed_tot  = leaf_wshed_tot  + leaf_wshed
               leaf_qwshed_tot = leaf_qwshed_tot + leaf_qwshed
               leaf_dwshed_tot = leaf_dwshed_tot + leaf_dwshed
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Exchange water and internal energy.                                   !
               !---------------------------------------------------------------------------!
               initp%leaf_water_im2(ico) = initp%leaf_water_im2(ico) - leaf_excess
               initp%wood_water_im2(ico) = initp%wood_water_im2(ico) + wood_demand
               initp%veg_water_im2 (ico) = initp%veg_water_im2 (ico) - leaf_wshed
               initp%leaf_energy   (ico) = initp%leaf_energy   (ico) - leaf_qexcess
               initp%wood_energy   (ico) = initp%wood_energy   (ico) + wood_qdemand
               initp%veg_energy    (ico) = initp%veg_energy    (ico) - leaf_qwshed
               !---------------------------------------------------------------------------!


               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_wflux_wl(ico) = initp%avg_wflux_wl(ico) - wood_demand !*hdidi
                  initp%avg_wshed_lg(ico) = initp%avg_wshed_lg(ico) + leaf_wshed  !*hdidi
               end if
               if (print_detailed) then
                  initp%cfx_qwflux_wl(ico) = initp%cfx_qwflux_wl(ico) - wood_qdemand !*hdidi
                  initp%cfx_qwshed   (ico) = initp%cfx_qwshed   (ico) + leaf_qwshed  !*hdidi
               end if
               !---------------------------------------------------------------------------!


            elseif (initp%leaf_water_im2(ico) < rk4min_leaf_water_im2) then
               !---------------------------------------------------------------------------!
               !     Leaves have too little water.  If possible, pull it from wood. If     !
               ! wood is also dry, leaves steal water from the canopy air space as fast    !
               ! "dew/frost" condensation.                                                 !
               !---------------------------------------------------------------------------!
               leaf_demand     = rk4min_leaf_water_im2 - initp%leaf_water_im2(ico)
               wood_excess     = min(leaf_demand,max(0.d0                                  &
                                    ,initp%wood_water_im2(ico)-rk4min_wood_water))
               leaf_dew        = leaf_demand - wood_excess
               wood_qexcess    = wood_excess * wood_water_uint
               leaf_qdew       = leaf_dew    * leaf_water_hint
               leaf_qdemand    = wood_qexcess + leaf_qdew
               !---------------------------------------------------------------------------!


               !----- Add the contribution of this cohort to the total boiling. -----------!
               leaf_dew_tot   = leaf_dew_tot   + leaf_dew
               leaf_qdew_tot  = leaf_qdew_tot  + leaf_qdew
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Exchange water and internal energy.                                   !
               !---------------------------------------------------------------------------!
               initp%leaf_water_im2(ico) = initp%leaf_water_im2(ico) + leaf_demand
               initp%wood_water_im2(ico) = initp%wood_water_im2(ico) - wood_excess
               initp%veg_water_im2 (ico) = initp%veg_water_im2 (ico) + leaf_dew
               initp%leaf_energy   (ico) = initp%leaf_energy   (ico) + leaf_qdemand
               initp%wood_energy   (ico) = initp%wood_energy   (ico) - wood_qexcess
               initp%veg_energy    (ico) = initp%veg_energy    (ico) + leaf_qdew
               !---------------------------------------------------------------------------!


               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_wflux_wl (ico) = initp%avg_wflux_wl (ico) + wood_excess  !*hdidi
                  initp%avg_transp   (ico) = initp%avg_transp   (ico) - leaf_dew     !*hdidi
               end if
               if (print_detailed) then
                  initp%cfx_qwflux_wl(ico) = initp%cfx_qwflux_wl(ico) + wood_qexcess !*hdidi
                  initp%cfx_qtransp  (ico) = initp%cfx_qtransp  (ico) - leaf_qdew    !*hdidi
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!

         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Check whether we can solve wood in this cohort...                            !
         !---------------------------------------------------------------------------------!
         if (initp%wood_resolvable(ico)) then

            if (initp%wood_water(ico) > max_wood_water) then

               !---------------------------------------------------------------------------!
               !    Too much water over the wood, we shall shed the excess to the ground.  !
               !---------------------------------------------------------------------------!
               wood_wshed  = initp%wood_water(ico) - max_wood_water
               wood_qwshed = wood_wshed                                                    &
                           * tl2uint8(initp%wood_temp(ico),initp%wood_fliq(ico))
               wood_dwshed = wood_wshed * ( initp%wood_fliq(ico) * wdnsi8                  &
                                          + (1.d0-initp%wood_fliq(ico)) * fdnsi8)
               !---------------------------------------------------------------------------!


               !----- Add the contribution of this cohort to the total shedding. ----------!
               wood_wshed_tot  = wood_wshed_tot  + wood_wshed
               wood_qwshed_tot = wood_qwshed_tot + wood_qwshed
               wood_dwshed_tot = wood_dwshed_tot + wood_dwshed
               !---------------------------------------------------------------------------!



               !----- Update water mass and energy. ---------------------------------------!
               initp%wood_water (ico) = initp%wood_water (ico) - wood_wshed
               initp%veg_water  (ico) = initp%veg_water  (ico) - wood_wshed
               initp%wood_energy(ico) = initp%wood_energy(ico) - wood_qwshed
               initp%veg_energy (ico) = initp%veg_energy (ico) - wood_qwshed
               !---------------------------------------------------------------------------!



               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_wshed_wg(ico) = initp%avg_wshed_wg(ico) + wood_wshed  ! * hdidi
               end if
               if (print_detailed) then
                  initp%cfx_qwshed  (ico) = initp%cfx_qwshed  (ico) + wood_qwshed ! * hdidi
               end if
               !---------------------------------------------------------------------------!


            elseif (initp%wood_water(ico) < min_wood_water) then
               !---------------------------------------------------------------------------!
               !    If wood_water is tiny and positive, exchange moisture with the air by  !
               ! donating the total amount as "boiling" (fast evaporation or sublimation). !
               ! In case the total is tiny but negative, exchange moisture with the air,   !
               ! "stealing" moisture as fast "dew/frost" condensation.                     !
               !---------------------------------------------------------------------------!
               wood_boil  = max(0.d0,  initp%wood_water(ico))
               wood_dew   = max(0.d0,- initp%wood_water(ico))
               wood_qboil = wood_boil * tq2enthalpy8(initp%wood_temp(ico),1.d0,.true.)
               wood_qdew  = wood_dew  * tq2enthalpy8(initp%wood_temp(ico),1.d0,.true.)
               !---------------------------------------------------------------------------!


               !----- Add the contribution of this cohort to the total boiling. -----------!
               wood_boil_tot  = wood_boil_tot  + wood_boil
               wood_dew_tot   = wood_dew_tot   + wood_dew
               wood_qboil_tot = wood_qboil_tot + wood_qboil
               wood_qdew_tot  = wood_qdew_tot  + wood_qdew
               !---------------------------------------------------------------------------!

               !----- Update cohort state variables. --------------------------------------!
               initp%wood_water (ico) = 0.d0
               initp%veg_water  (ico) = initp%veg_water  (ico)  + wood_dew  - wood_boil
               initp%wood_energy(ico) = initp%wood_energy(ico)  + wood_qdew - wood_qboil
               initp%veg_energy (ico) = initp%veg_energy (ico)  + wood_qdew - wood_qboil
               !---------------------------------------------------------------------------!

               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_vapor_wc(ico) = initp%avg_vapor_wc(ico)                        &
                                          + (wood_boil  - wood_dew ) ! * hdidi
               end if
               if (print_detailed) then
                  initp%cfx_qwflxwc (ico) = initp%cfx_qwflxwc (ico)                        &
                                          + (wood_qboil - wood_qdew) ! * hdidi
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !       Check whether wood internal water is bounded.                          !
            !------------------------------------------------------------------------------!
            if (initp%wood_water_im2(ico) > rk4max_wood_water_im2) then
               !---------------------------------------------------------------------------!
               !     Wood has too much water.  If possible, send water to leaves.  If      !
               ! leaves are also saturated, wood expels the excess water to the top soil   !
               ! layer.  In case this is not sufficient (unlikely as this excess is tiny), !
               ! the excess water goes to shedding.                                        !
               !---------------------------------------------------------------------------!
               !----- First guess. --------------------------------------------------------!
               wood_excess = rk4max_wood_water_im2 - initp%wood_water_im2(ico)
               leaf_demand = rk4max_leaf_water_im2 - initp%leaf_water_im2(ico)
               soil_demand = (soil8(nstop)%soilbp-initp%soil_water(kt)) * dslz8(kt) * wdns8
               !----- Bounded guess. ------------------------------------------------------!
               leaf_demand = max(0.d0,min(wood_excess            ,leaf_demand))
               soil_demand = max(0.d0,min(wood_excess-leaf_demand,soil_demand))
               wood_wshed  = max(0.d0,wood_excess - leaf_demand - soil_demand)
               !----- Energy associated with transfers. -----------------------------------!
               leaf_qdemand = leaf_demand * wood_water_uint
               soil_qdemand = soil_demand * wood_water_uint
               wood_qwshed  = wood_wshed  * wood_water_uint
               wood_dwshed  = wood_wshed  * wood_water_zint
               wood_qexcess = leaf_qdemand + soil_qdemand + wood_qwshed
               !---------------------------------------------------------------------------!



               !----- Add the contribution of this cohort to the total shedding. ----------!
               wood_wshed_tot  = wood_wshed_tot  + wood_wshed
               wood_qwshed_tot = wood_qwshed_tot + wood_qwshed
               wood_dwshed_tot = wood_dwshed_tot + wood_dwshed
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Exchange water and internal energy.                                   !
               !---------------------------------------------------------------------------!
               initp%wood_water_im2(ico) = initp%wood_water_im2 (ico) - wood_excess
               initp%leaf_water_im2(ico) = initp%leaf_water_im2 (ico) + leaf_demand
               initp%veg_water_im2 (ico) = initp%veg_water_im2  (ico) - soil_demand        &
                                                                      - wood_wshed
               initp%soil_water    ( kt) = initp%soil_water     ( kt) + soil_demand        &
                                                                      * dslzi8(kt) * wdnsi8
               initp%wood_energy   (ico) = initp%wood_energy    (ico) - wood_qexcess
               initp%leaf_energy   (ico) = initp%leaf_energy    (ico) + leaf_qdemand
               initp%veg_energy    (ico) = initp%veg_energy     (ico) - soil_qdemand       &
                                                                      - wood_qwshed
               initp%soil_energy   ( kt) = initp%soil_energy    ( kt) + soil_qdemand       &
                                                                      * dslzi8(kt)
               !---------------------------------------------------------------------------!



               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_wflux_wl          (ico) = initp%avg_wflux_wl          (ico)    &
                                                    + leaf_demand
                  initp%avg_wflux_gw          (ico) = initp%avg_wflux_wl          (ico)    &
                                                    - soil_demand
                  initp%avg_wshed_lg          (ico) = initp%avg_wshed_lg          (ico)    &
                                                    + leaf_wshed
                  initp%avg_wflux_gw_layer (kt,ico) = initp%avg_wflux_gw_layer (kt,ico)    &
                                                    - soil_demand
               end if
               if (print_detailed) then
                  initp%cfx_qwflux_wl         (ico) = initp%cfx_qwflux_wl         (ico)    &
                                                    + leaf_qdemand
                  initp%cfx_qwflux_gw         (ico) = initp%cfx_qwflux_gw         (ico)    &
                                                    - soil_qdemand
                  initp%cfx_qwshed            (ico) = initp%cfx_qwshed            (ico)    &
                                                    + leaf_qwshed 
                  initp%cfx_qwflux_gw_layer(kt,ico) = initp%cfx_qwflux_gw_layer(kt,ico)    &
                                                    - soil_qdemand
               end if
               !---------------------------------------------------------------------------!


            elseif (initp%wood_water_im2(ico) < rk4min_wood_water_im2) then
               !---------------------------------------------------------------------------!
               !     Leaves have too little water.  If possible, pull it from wood. If     !
               ! leaf is also dry, take water from the top soil layer.  As a last          !
               ! resource, extract from the canopy air space (very unlikely to occur).                                                !
               !---------------------------------------------------------------------------!
               !----- First guess. --------------------------------------------------------!
               wood_demand = rk4min_wood_water_im2 - initp%wood_water_im2(ico)
               leaf_excess = initp%leaf_water_im2(ico) - rk4min_leaf_water_im2
               soil_excess = (initp%soil_water(kt)-soil8(nstop)%soilcp) * dslz8(kt) * wdns8
               !----- Bounded guess. ------------------------------------------------------!
               leaf_excess = max(0.d0,min(wood_demand            ,leaf_excess))
               soil_excess = max(0.d0,min(wood_demand-leaf_excess,soil_excess))
               wood_dew    = max(0.d0,wood_demand - leaf_excess - soil_excess)
               !----- Energy associated with transfers. -----------------------------------!
               leaf_qexcess = leaf_excess * leaf_water_uint
               soil_qexcess = soil_excess * soil_water_uint
               wood_qdew    = wood_dew    * wood_water_hint
               wood_qdemand = leaf_qexcess + soil_qexcess + wood_qdew
               !---------------------------------------------------------------------------!


               !----- Add the contribution of this cohort to the total boiling. -----------!
               wood_dew_tot   = wood_dew_tot   + wood_dew
               wood_qdew_tot  = wood_qdew_tot  + wood_qdew
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Exchange water and internal energy.                                   !
               !---------------------------------------------------------------------------!
               initp%wood_water_im2(ico) = initp%wood_water_im2(ico) + wood_demand
               initp%leaf_water_im2(ico) = initp%leaf_water_im2(ico) - leaf_excess
               initp%veg_water_im2 (ico) = initp%veg_water_im2 (ico) + soil_excess         &
                                                                     + wood_dew
               initp%soil_water    ( kt) = initp%soil_water    ( kt) - soil_excess         &
                                                                     * dslzi8(kt) * wdnsi8
               initp%wood_energy   (ico) = initp%wood_energy   (ico) + wood_qdemand
               initp%leaf_energy   (ico) = initp%leaf_energy   (ico) - leaf_qexcess
               initp%veg_energy    (ico) = initp%veg_energy    (ico) + soil_qexcess        &
                                                                     + wood_qdew
               initp%soil_energy   ( kt) = initp%soil_energy   ( kt) - soil_qexcess        &
                                                                     * dslzi8(kt)
               !---------------------------------------------------------------------------!


               !----- Update fluxes if needed be. -----------------------------------------!
               if (fast_diagnostics) then
                  initp%avg_wflux_wl          (ico) = initp%avg_wflux_wl          (ico)    &
                                                    - leaf_excess
                  initp%avg_wflux_gw          (ico) = initp%avg_wflux_wl          (ico)    &
                                                    + soil_excess
                  initp%avg_transp            (ico) = initp%avg_transp            (ico)    &
                                                    - wood_dew
                  initp%avg_wflux_gw_layer (kt,ico) = initp%avg_wflux_gw_layer (kt,ico)    &
                                                    + soil_excess
               end if
               if (print_detailed) then
                  initp%cfx_qwflux_wl         (ico) = initp%cfx_qwflux_wl         (ico)    &
                                                    - leaf_qexcess
                  initp%cfx_qwflux_gw         (ico) = initp%cfx_qwflux_gw         (ico)    &
                                                    + soil_qexcess
                  initp%cfx_qtransp           (ico) = initp%cfx_qtransp           (ico)    &
                                                    - leaf_qdew
                  initp%cfx_qwflux_gw_layer(kt,ico) = initp%cfx_qwflux_gw_layer(kt,ico)    &
                                                    + soil_qexcess
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!

         end if
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    The water that fell from the leaves and branches must go somewhere...  Here we  !
      ! decide which place is the most suitable.  In case there is already a temporary     !
      ! surface water layer we can add the water there, otherwise we dump it into the      !
      ! virtual layer, which may or may not become a temporary surface water layer.        !
      !------------------------------------------------------------------------------------!
      ksn = initp%nlev_sfcwater
      select case(initp%flag_sfcwater)
      case (0)
         !------ No temporary water, shed the water into the virtual layer. ---------------!
         initp%virtual_water  = initp%virtual_water  + leaf_wshed_tot  + wood_wshed_tot
         initp%virtual_energy = initp%virtual_energy + leaf_qwshed_tot + wood_qwshed_tot
         initp%virtual_depth  = initp%virtual_depth  + leaf_dwshed_tot + wood_dwshed_tot
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !     There is a temporary water, shed the excess water to the temporary surface  !
         ! water layer.                                                                    !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_mass(ksn)   = initp%sfcwater_mass(ksn)                             &
                                    + leaf_wshed_tot  + wood_wshed_tot
         initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn)                           &
                                    + leaf_qwshed_tot + wood_qwshed_tot
         initp%sfcwater_depth(ksn)  = initp%sfcwater_depth(ksn)                            &
                                    + leaf_dwshed_tot + wood_dwshed_tot
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!



      !----- Update the canopy air specific humidity and enthalpy. ------------------------!
      initp%can_shv      = initp%can_shv                                                   &
                         + ( leaf_boil_tot  + wood_boil_tot                                &
                           - leaf_dew_tot   - wood_dew_tot )                               &
                         * initp%wcapcani
      initp%can_enthalpy = initp%can_enthalpy                                              &
                         + ( leaf_qboil_tot + wood_qboil_tot                               &
                           - leaf_qdew_tot  - wood_qdew_tot )                              &
                         * initp%hcapcani
      !------------------------------------------------------------------------------------!

      return
   end subroutine adjust_veg_properties
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine print_errmax(errmax,yerr,yscal,cpatch,y)
      use rk4_coms              , only : rk4patchtype       & ! Structure
                                       , ibranch_thermo     & ! intent(in)
                                       , rk4site            & ! intent(in)
                                       , checkbudget        ! ! intent(in)
      use ed_state_vars         , only : patchtype          ! ! Structure
      use grid_coms             , only : nzg                ! ! intent(in)
      use physiology_coms       , only : plant_hydro_scheme ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target       :: yerr,yscal,y
      type(patchtype)    , target       :: cpatch
      real(kind=8)       , intent(out)  :: errmax
      !----- Local variables --------------------------------------------------------------!
      integer                           :: ico
      integer                           :: k
      logical                           :: troublemaker
      !----- Constants --------------------------------------------------------------------!
      character(len=28)  , parameter    :: onefmt = '(a16,1x,3(es12.4,1x),11x,l1)'
      character(len=34)  , parameter    :: lyrfmt = '(a16,1x,i6,1x,3(es12.4,1x),11x,l1)'
      character(len=34)  , parameter    :: cohfmt = '(a16,1x,i6,1x,6(es12.4,1x),11x,l1)'
      !------------------------------------------------------------------------------------!


      write(unit=*,fmt='(80a)'    ) ('=',k=1,80)
      write(unit=*,fmt='(a)'      ) '  ..... PRINTING MAXIMUM ERROR INFORMATION: .....'
      write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
      write(unit=*,fmt='(a)'      ) 
      write(unit=*,fmt='(a)'      ) ' Patch level variables, single layer:'
      write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
      write(unit=*,fmt='(5(a,1x))')  'Name            ','   Max.Error','   Abs.Error'&
                                   &,'       Scale','Problem(T|F)'

      errmax       = max(0.d0,abs(yerr%can_enthalpy/yscal%can_enthalpy))
      troublemaker = large_error(yerr%can_enthalpy,yscal%can_enthalpy)
      write(unit=*,fmt=onefmt) 'CAN_ENTHALPY:',errmax,yerr%can_enthalpy,yscal%can_enthalpy &
                                              ,troublemaker

      errmax       = max(errmax,abs(yerr%can_shv/yscal%can_shv))
      troublemaker = large_error(yerr%can_shv,yscal%can_shv)
      write(unit=*,fmt=onefmt) 'CAN_SHV:',errmax,yerr%can_shv,yscal%can_shv,troublemaker

      errmax = max(errmax,abs(yerr%can_co2/yscal%can_co2))
      troublemaker = large_error(yerr%can_co2,yscal%can_co2)
      write(unit=*,fmt=onefmt) 'CAN_CO2:',errmax,yerr%can_co2,yscal%can_co2,troublemaker

      errmax = max(errmax,abs(yerr%virtual_energy/yscal%virtual_energy))
      troublemaker = large_error(yerr%virtual_energy,yscal%virtual_energy)
      write(unit=*,fmt=onefmt) 'VIRTUAL_ENERGY:',errmax,yerr%virtual_energy                &
                                                ,yscal%virtual_energy,troublemaker

      errmax = max(errmax,abs(yerr%virtual_water/yscal%virtual_water))
      troublemaker = large_error(yerr%virtual_water,yscal%virtual_water)
      write(unit=*,fmt=onefmt) 'VIRTUAL_WATER:',errmax,yerr%virtual_water                  &
                                               ,yscal%virtual_water,troublemaker

      write(unit=*,fmt='(80a)') ('-',k=1,80)
      write(unit=*,fmt='(a)'  ) 
      write(unit=*,fmt='(80a)') ('-',k=1,80)
      write(unit=*,fmt='(a)'      ) ' Patch level variables, soil layers:'
      write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'            &
                                   &,'   Abs.Error','       Scale','Problem(T|F)'

      do k=rk4site%lsl,nzg
         errmax = max(errmax,abs(yerr%soil_water(k)/yscal%soil_water(k)))
         troublemaker = large_error(yerr%soil_water(k),yscal%soil_water(k))
         write(unit=*,fmt=lyrfmt) 'SOIL_WATER:',k,errmax,yerr%soil_water(k)                &
                                               ,yscal%soil_water(k),troublemaker

         errmax       = max(errmax,abs(yerr%soil_energy(k)/yscal%soil_energy(k)))
         troublemaker = large_error(yerr%soil_energy(k),yscal%soil_energy(k))
         write(unit=*,fmt=lyrfmt) 'SOIL_ENERGY:',k,errmax,yerr%soil_energy(k)              &
                                                ,yscal%soil_energy(k),troublemaker
      enddo

      if (yerr%nlev_sfcwater > 0) then
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'  ) 
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'      ) ' Patch level variables, water/snow layers:'
         write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'         &
                                   &,'   Abs.Error','       Scale','Problem(T|F)'
         do k=1,yerr%nlev_sfcwater
            errmax = max(errmax,abs(yerr%sfcwater_energy(k)/yscal%sfcwater_energy(k)))
            troublemaker = large_error(yerr%sfcwater_energy(k),yscal%sfcwater_energy(k))
            write(unit=*,fmt=lyrfmt) 'SFCWATER_ENERGY:',k,errmax,yerr%sfcwater_energy(k)   &
                                                       ,yscal%sfcwater_energy(k)           &
                                                       ,troublemaker

            errmax = max(errmax,abs(yerr%sfcwater_mass(k)/yscal%sfcwater_mass(k)))
            troublemaker = large_error(yerr%sfcwater_mass(k),yscal%sfcwater_mass(k))
            write(unit=*,fmt=lyrfmt) 'SFCWATER_MASS:',k,errmax,yerr%sfcwater_mass(k)       &
                                                     ,yscal%sfcwater_mass(k),troublemaker
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Choose what to print based on wood thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (1)
         !----- Combined case, print vegetation only. -------------------------------------!
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'  ) 
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'      ) ' Veg-level variables (only the resolvable ones):'
         write(unit=*,fmt='(9(a,1x))')         'Name            ','   PFT','         LAI'  &
                                            ,'         WAI','         TAI','   Max.Error'  &
                                            ,'   Abs.Error','       Scale','Problem(T|F)'
         do ico = 1,cpatch%ncohorts
            if (y%veg_resolvable(ico)) then
               errmax       = max(errmax,abs(yerr%veg_water(ico)/yscal%veg_water(ico)))
               troublemaker = large_error(yerr%veg_water(ico),yscal%veg_water(ico))
               write(unit=*,fmt=cohfmt) 'VEG_WATER:',cpatch%pft(ico),y%lai(ico),y%wai(ico) &
                                                     ,y%tai(ico),errmax                    &
                                                     ,yerr%veg_water(ico)                  &
                                                     ,yscal%veg_water(ico),troublemaker
                    

               errmax       = max(errmax,abs(yerr%veg_energy(ico)/yscal%veg_energy(ico)))
               troublemaker = large_error(yerr%veg_energy(ico),yscal%veg_energy(ico))
               write(unit=*,fmt=cohfmt) 'VEG_ENERGY:',cpatch%pft(ico),cpatch%lai(ico)      &
                                                      ,y%wai(ico),y%tai(ico)               &
                                                      ,errmax,yerr%veg_energy(ico)         &
                                                      ,yscal%veg_energy(ico)               &
                                                      ,troublemaker

            end if
         end do
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !     Leaf-only or two separate pools, print leaves if they are resolved.         !
         !---------------------------------------------------------------------------------!
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'  ) 
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'      ) ' Leaf-level variables (only the resolvable ones):'
         write(unit=*,fmt='(9(a,1x))')         'Name            ','   PFT','         LAI'  &
                                            ,'         WAI','         TAI','   Max.Error'  &
                                            ,'   Abs.Error','       Scale','Problem(T|F)'
         do ico = 1,cpatch%ncohorts
            if (y%leaf_resolvable(ico)) then
               errmax = max(errmax,abs(yerr%leaf_water(ico)/yscal%leaf_water(ico)))
               troublemaker = large_error(yerr%leaf_water(ico),yscal%leaf_water(ico))
               write(unit=*,fmt=cohfmt) 'LEAF_WATER:',cpatch%pft(ico),y%lai(ico)           &
                                                     ,y%wai(ico),y%tai(ico),errmax         &
                                                     ,yerr%leaf_water(ico)                 &
                                                     ,yscal%leaf_water(ico),troublemaker

               errmax = max(errmax,abs(yerr%leaf_energy(ico)/yscal%leaf_energy(ico)))
               troublemaker = large_error(yerr%leaf_energy(ico),yscal%leaf_energy(ico))
               write(unit=*,fmt=cohfmt) 'LEAF_ENERGY:',cpatch%pft(ico),cpatch%lai(ico)     &
                                                      ,y%wai(ico),y%tai(ico)               &
                                                      ,errmax,yerr%leaf_energy(ico)        &
                                                      ,yscal%leaf_energy(ico)              &
                                                      ,troublemaker

            end if
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Branchwood, this is solved only if ibranch_thermo is 2.                      !
         !---------------------------------------------------------------------------------!
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'  ) 
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'      ) ' Wood-level variables (only the resolvable ones):'
         write(unit=*,fmt='(9(a,1x))')         'Name            ','   PFT','         LAI'  &
                                            ,'         WAI','         TAI','   Max.Error'  &
                                            ,'   Abs.Error','       Scale','Problem(T|F)'
         do ico = 1,cpatch%ncohorts
            if (y%wood_resolvable(ico)) then
               errmax = max(errmax,abs(yerr%wood_water(ico)/yscal%wood_water(ico)))
               troublemaker = large_error(yerr%wood_water(ico),yscal%wood_water(ico))
               write(unit=*,fmt=cohfmt) 'WOOD_WATER:',cpatch%pft(ico),y%lai(ico)           &
                                                     ,y%wai(ico),y%tai(ico),errmax         &
                                                     ,yerr%wood_water(ico)                 &
                                                     ,yscal%wood_water(ico),troublemaker

               errmax = max(errmax,abs(yerr%wood_energy(ico)/yscal%wood_energy(ico)))
               troublemaker = large_error(yerr%wood_energy(ico),yscal%wood_energy(ico))
               write(unit=*,fmt=cohfmt) 'WOOD_ENERGY:',cpatch%pft(ico),cpatch%lai(ico)     &
                                                      ,y%wai(ico),y%tai(ico)               &
                                                      ,errmax,yerr%wood_energy(ico)        &
                                                      ,yscal%wood_energy(ico),troublemaker
            end if
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Plant internal water.  Print these only if plant hydraulics is enabled.       !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Skip error calculation -- plant hydraulics is disabled. -------------------!
         continue
         !---------------------------------------------------------------------------------!
      case default
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'  ) 
         write(unit=*,fmt='(80a)') ('-',k=1,80)
         write(unit=*,fmt='(a)'      ) ' Wood/Leaf WATER_IM2  (only the resolvable ones):'
         write(unit=*,fmt='(9(a,1x))')         'Name            ','   PFT','         LAI'  &
                                            ,'         WAI','         TAI','   Max.Error'  &
                                            ,'   Abs.Error','       Scale','Problem(T|F)'
         do ico = 1,cpatch%ncohorts
            if (y%leaf_resolvable(ico)) then
               errmax       = max( errmax                                                  &
                                 , abs(yerr%leaf_water_im2(ico)/yscal%leaf_water_im2(ico)))
               troublemaker = large_error( yerr%leaf_water_im2 (ico)                       &
                                         , yscal%leaf_water_im2(ico) )
               write(unit=*,fmt=cohfmt) 'LEAF_WATER_IM2:'                                  &
                                       ,cpatch%pft(ico),y%lai(ico),y%wai(ico)              &
                                       ,y%tai(ico),errmax                                  &
                                       ,yerr%leaf_water_im2(ico)                           &
                                       ,yscal%leaf_water_im2(ico),troublemaker
            end if

            if (y%wood_resolvable(ico)) then
               errmax       = max( errmax                                                  &
                                 , abs(yerr%wood_water_im2(ico)/yscal%wood_water_im2(ico)))
               troublemaker = large_error( yerr%wood_water_im2 (ico)                       &
                                         , yscal%wood_water_im2(ico) )
               write(unit=*,fmt=cohfmt) 'WOOD_WATER_IM2:'                                  &
                                       ,cpatch%pft(ico),y%lai(ico),y%wai(ico)              &
                                       ,y%tai(ico),errmax                                  &
                                       ,yerr%wood_water_im2(ico)                           &
                                       ,yscal%wood_water_im2(ico),troublemaker
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Here we just need to make sure the user is checking mass, otherwise these      !
      ! variables will not be computed.  If this turns out to be essential, we will make   !
      ! this permanent and not dependent on checkbudget.  The only one that is not checked !
      ! is the runoff, because it is computed only after a step is accepted.               !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
         write(unit=*,fmt='(a)'      ) 
         write(unit=*,fmt='(a)'      ) ' Budget variables, single layer:'
         write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
         write(unit=*,fmt='(5(a,1x))')  'Name            ','   Max.Error','   Abs.Error'   &
                                      &,'       Scale','Problem(T|F)'
         errmax = max(errmax                                                               &
                     ,abs(yerr%co2budget_loss2atm/yscal%co2budget_loss2atm))
         troublemaker = large_error(yerr%co2budget_loss2atm                                &
                                   ,yscal%co2budget_loss2atm)
         write(unit=*,fmt=onefmt) 'CO2LOSS2ATM:',errmax,yerr%co2budget_loss2atm            &
                                 ,yscal%co2budget_loss2atm,troublemaker

         errmax = max(errmax                                                               &
                     ,abs(yerr%ebudget_netrad/yscal%ebudget_netrad))
         troublemaker = large_error(yerr%ebudget_netrad                                    &
                                   ,yscal%ebudget_netrad)
         write(unit=*,fmt=onefmt) 'ENNETRAD:',errmax,yerr%ebudget_netrad                   &
                                 ,yscal%ebudget_netrad,troublemaker

         errmax = max(errmax                                                               &
                     ,abs(yerr%ebudget_loss2atm/yscal%ebudget_loss2atm))
         troublemaker = large_error(yerr%ebudget_loss2atm                                  &
                                   ,yscal%ebudget_loss2atm)
         write(unit=*,fmt=onefmt) 'ENLOSS2ATM:',errmax,yerr%ebudget_loss2atm               &
                                 ,yscal%ebudget_loss2atm,troublemaker

         errmax = max(errmax                                                               &
                     ,abs(yerr%wbudget_loss2atm/yscal%wbudget_loss2atm))
         troublemaker = large_error(yerr%wbudget_loss2atm                                  &
                                   ,yscal%wbudget_loss2atm)
         write(unit=*,fmt=onefmt) 'H2OLOSS2ATM:',errmax,yerr%wbudget_loss2atm              &
                                 ,yscal%wbudget_loss2atm,troublemaker

         errmax = max(errmax,abs( yerr%ebudget_loss2drainage                               &
                                / yscal%ebudget_loss2drainage))
         troublemaker = large_error(yerr%ebudget_loss2drainage                             &
                                   ,yscal%ebudget_loss2drainage)
         write(unit=*,fmt=onefmt) 'ENDRAINAGE:',errmax                                     &
                                 ,yerr%ebudget_loss2drainage                               &
                                 ,yscal%ebudget_loss2drainage,troublemaker

         errmax = max(errmax,abs( yerr%wbudget_loss2drainage                               &
                                / yscal%wbudget_loss2drainage))
         troublemaker = large_error(yerr%wbudget_loss2drainage                             &
                                   ,yscal%wbudget_loss2drainage)
         write(unit=*,fmt=onefmt) 'H2ODRAINAGE:',errmax                                    &
                                 ,yerr%wbudget_loss2drainage                               &
                                 ,yscal%wbudget_loss2drainage,troublemaker

         errmax = max(errmax                                                               &
                     ,abs(yerr%co2budget_storage/yscal%co2budget_storage))
         troublemaker = large_error(yerr%co2budget_storage                                 &
                                   ,yscal%co2budget_storage)
         write(unit=*,fmt=onefmt) 'CO2STORAGE:',errmax,yerr%co2budget_storage              &
                                 ,yscal%co2budget_storage,troublemaker

         errmax = max(errmax                                                               &
                     ,abs(yerr%ebudget_storage/yscal%ebudget_storage))
         troublemaker = large_error(yerr%ebudget_storage                                   &
                                   ,yscal%ebudget_storage)
         write(unit=*,fmt=onefmt) 'ENSTORAGE:',errmax,yerr%ebudget_storage                 &
                                 ,yscal%ebudget_storage,troublemaker

         errmax = max(errmax                                                               &
                     ,abs(yerr%wbudget_storage/yscal%wbudget_storage))
         troublemaker = large_error(yerr%wbudget_storage                                   &
                                   ,yscal%wbudget_storage)
         write(unit=*,fmt=onefmt) 'H2OSTORAGE:',errmax,yerr%wbudget_storage                &
                                 ,yscal%wbudget_storage,troublemaker
      end if

      write(unit=*,fmt='(a)'  ) 
      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(a)'  ) 

      return
   end subroutine print_errmax
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine prints the patch and cohort information when the model falls       !
   ! apart...                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine print_csiteipa(csite, ipa)
      use rk4_coms              , only : rk4site       ! ! intent(in)
      use ed_state_vars         , only : sitetype      & ! structure
                                       , patchtype     ! ! structure
      use ed_misc_coms          , only : current_time  ! ! intent(in)
      use grid_coms             , only : nzg           ! ! intent(in)
      use ed_max_dims           , only : n_pft         ! ! intent(in)
      use consts_coms           , only : day_sec       & ! intent(in)
                                       , umol_2_kgC    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: ipa
      !----- Local variable ---------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: k
      real                         :: leaf_growth_resp
      real                         :: root_growth_resp
      real                         :: sapa_growth_resp
      real                         :: sapb_growth_resp
      real                         :: barka_growth_resp
      real                         :: barkb_growth_resp
      real                         :: leaf_storage_resp
      real                         :: root_storage_resp
      real                         :: sapa_storage_resp
      real                         :: sapb_storage_resp
      real                         :: barka_storage_resp
      real                         :: barkb_storage_resp
      real                         :: pss_lai
      real                         :: pss_wai
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)


      !----- Find the total patch LAI and WAI. --------------------------------------------!
      pss_lai = 0.0
      pss_wai = 0.0
      do ico=1,cpatch%ncohorts
         pss_lai = pss_lai + cpatch%lai(ico)
         pss_wai = pss_wai + cpatch%wai(ico)
      end do
      !------------------------------------------------------------------------------------!


      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(80a)') ('=',k=1,80)

      write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (csite) ||||'

      write(unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,3(i2.2,a))')                              &
            'Time:',current_time%month,'/',current_time%date,'/',current_time%year         &
                   ,current_time%hour,':',current_time%min,':',current_time%sec,' UTC'
      write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
      write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
    
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a)'  ) 'Leaf information (only the resolvable ones shown): '
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),12(a12,1x))')                                           &
            '    PFT','KRDEPTH','      NPLANT','         LAI','         DBH'               &
                               ,'      BDEADA','      BDEADB','       BLEAF'               &
                               ,' LEAF_ENERGY','  LEAF_WATER','LEAF_H2O_IM2'               &
                               ,'   LEAF_HCAP','   LEAF_TEMP','   LEAF_FLIQ'
      do ico = 1,cpatch%ncohorts
         if (cpatch%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),12(es12.4,1x))')                                   &
                  cpatch%pft(ico),cpatch%krdepth(ico)                                      &
                 ,cpatch%nplant(ico),cpatch%lai(ico),cpatch%dbh(ico),cpatch%bdeada(ico)    &
                 ,cpatch%bdeadb(ico),cpatch%bleaf(ico),cpatch%leaf_energy(ico)             &
                 ,cpatch%leaf_water(ico),cpatch%leaf_water_im2(ico),cpatch%leaf_hcap(ico)  &
                 ,cpatch%leaf_temp(ico),cpatch%leaf_fliq(ico)
         end if
      end do
      write (unit=*,fmt='(2(a7,1x),6(a12,1x))')                                            &
            '    PFT','KRDEPTH','         LAI','     FS_OPEN','         FSW'               &
                               ,'         FSN','         GPP','   LEAF_RESP'
      do ico = 1,cpatch%ncohorts
         if (cpatch%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),6(es12.4,1x))')                                    &
                  cpatch%pft(ico),cpatch%krdepth(ico)                                      &
                 ,cpatch%lai(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)      &
                 ,cpatch%gpp(ico),cpatch%leaf_respiration(ico)
         end if
      end do
      write (unit=*,fmt='(2(a7,1x),3(a12,1x),12(a17,1x))')                                 &
            '    PFT','KRDEPTH','         LAI','   ROOT_RESP','   STEM_RESP'               &
           ,'  LEAF_STORE_RESP','  ROOT_STORE_RESP','  SAPA_STORE_RESP'                    &
           ,'  SAPB_STORE_RESP',' BARKA_STORE_RESP',' BARKB_STORE_RESP'                    &
           ,' LEAF_GROWTH_RESP',' ROOT_GROWTH_RESP',' SAPA_GROWTH_RESP'                    &
           ,' SAPB_GROWTH_RESP','BARKA_GROWTH_RESP','BARKB_GROWTH_RESP'
      do ico = 1,cpatch%ncohorts
         if (cpatch%leaf_resolvable(ico)) then
            leaf_growth_resp   = cpatch%leaf_growth_resp  (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            root_growth_resp   = cpatch%root_growth_resp  (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            sapa_growth_resp   = cpatch%sapa_growth_resp  (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            sapb_growth_resp   = cpatch%sapb_growth_resp  (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            barka_growth_resp  = cpatch%barka_growth_resp (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            barkb_growth_resp  = cpatch%barkb_growth_resp (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            leaf_storage_resp  = cpatch%leaf_storage_resp (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            root_storage_resp  = cpatch%root_storage_resp (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            sapa_storage_resp  = cpatch%sapa_storage_resp (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            sapb_storage_resp  = cpatch%sapb_storage_resp (ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            barka_storage_resp = cpatch%barka_storage_resp(ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)
            barkb_storage_resp = cpatch%barkb_storage_resp(ico) * cpatch%nplant(ico)       &
                               / (day_sec * umol_2_kgC)

            write(unit=*,fmt='(2(i7,1x),3(es12.5,1x),12(es17.5,1x))')                      &
                    cpatch%pft(ico)  ,  cpatch%krdepth(ico)                                &
                 ,  cpatch%lai(ico)  ,  cpatch%root_respiration(ico)                       &
                 ,  cpatch%stem_respiration(ico)                                           &
                 ,  leaf_storage_resp,  root_storage_resp,  sapa_storage_resp              &
                 ,  sapb_storage_resp, barka_storage_resp, barkb_storage_resp              &
                 ,   leaf_growth_resp,   root_growth_resp,   sapa_growth_resp              &
                 ,   sapb_growth_resp,  barka_growth_resp,  barkb_growth_resp
         end if
      end do
      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a)'  ) 'Wood information (only the resolvable ones shown): '
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),12(a12,1x))')                                           &
            '    PFT','KRDEPTH','      NPLANT','         WAI','         DBH'               &
                               ,'      BDEADA','      BDEADB','   BSAPWOODA'               &
                               ,'   BSAPWOODB','      BBARKA','      BBARKB'               &
                               ,' WOOD_ENERGY','   WOOD_TEMP','  WOOD_WATER'
      do ico = 1,cpatch%ncohorts
         if (cpatch%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),12(es12.4,1x))')                                   &
                  cpatch%pft(ico),cpatch%krdepth(ico)                                      &
                 ,cpatch%nplant(ico),cpatch%wai(ico),cpatch%dbh(ico),cpatch%bdeada(ico)    &
                 ,cpatch%bdeadb(ico),cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico)           &
                 ,cpatch%bbarka(ico),cpatch%bbarkb(ico),cpatch%wood_energy(ico)            &
                 ,cpatch%wood_temp(ico),cpatch%wood_water(ico)
         end if
      end do
      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(14(a12,1x))')  '   DIST_TYPE','         AGE','        AREA'      &
                                        ,'          RH','      FGC_RH','      FSC_RH'      &
                                        ,'     STGC_RH','     STSC_RH','      MSC_RH'      &
                                        ,'      SSC_RH','      PSC_RH','AVGDAILY_TMP'      &
                                        ,'     SUM_CHD','     SUM_DGD'
      write (unit=*,fmt='(i12,1x,13(es12.4,1x))')  csite%dist_type(ipa),csite%age(ipa)     &
            ,csite%area(ipa),csite%rh(ipa),csite%fgc_rh(ipa),csite%fsc_rh(ipa)             &
            ,csite%stgc_rh(ipa),csite%stsc_rh(ipa),csite%msc_rh(ipa),csite%ssc_rh(ipa)     &
            ,csite%psc_rh(ipa),csite%avg_daily_temp(ipa),csite%sum_chd(ipa)                &
            ,csite%sum_dgd(ipa)

      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(9(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','VEG_DISPLACE'       &
                                       ,'   PATCH_LAI','   PATCH_WAI','        HTRY'       &
                                       ,'    CAN_RHOS','    CAN_DMOL','   CAN_DEPTH'
      write (unit=*,fmt='(9(es12.4,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)       &
                                         ,csite%veg_displace(ipa),pss_lai,pss_wai          &
                                         ,csite%htry(ipa),csite%can_rhos(ipa)              &
                                         ,csite%can_dmol(ipa),csite%can_depth(ipa)

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(7(a12,1x))')  '   CAN_THEIV','    CAN_TEMP','     CAN_SHV'       &
                                       ,'    CAN_PRSS','     CAN_CO2','   CAN_VPDEF'       &
                                       ,'       GGNET'
      write (unit=*,fmt='(7(es12.4,1x))')  csite%can_theiv (ipa),csite%can_temp  (ipa)     &
                                         , csite%can_shv   (ipa),csite%can_prss  (ipa)     &
                                         , csite%can_co2   (ipa),csite%can_vpdef (ipa)     &
                                         , csite%ggnet     (ipa)

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(10(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'      &
                                        ,'       TSTAR','        ZETA','     RI_BULK'      &
                                        ,'     RLONG_G','    RSHORT_G','       PAR_G'      &
                                        ,'     RLONG_S'
      write (unit=*,fmt='(10(es12.4,1x))') csite%ustar(ipa),csite%qstar(ipa)               &
                                          ,csite%cstar(ipa),csite%tstar(ipa)               &
                                          ,csite%zeta(ipa),csite%ribulk(ipa)               &
                                          ,csite%rlong_g(ipa),csite%rshort_g(ipa)          &
                                          ,csite%par_g(ipa),csite%rlong_s(ipa)

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(a5,1x,a12)') '  PFT','       REPRO'
      do k=1,n_pft
         write (unit=*,fmt='(i5,1x,es12.4)') k,csite%repro(k,ipa)
      end do

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'       &
                                      &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
      do k=rk4site%lsl,nzg
         write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,rk4site%ntext_soil(k)             &
               ,csite%soil_energy(k,ipa),csite%soil_tempk(k,ipa),csite%soil_water(k,ipa)   &
               ,csite%soil_fracliq(k,ipa)
      end do
      
      if (csite%nlev_sfcwater(ipa) >= 1) then
         write (unit=*,fmt='(80a)') ('-',k=1,80)
         write (unit=*,fmt='(a5,1x,7(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'    &
                                          ,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'    &
                                          ,'    RSHORT_S','       PAR_S'
         do k=1,csite%nlev_sfcwater(ipa)
            write (unit=*,fmt='(i5,1x,7(es12.4,1x))') k,csite%sfcwater_energy(k,ipa)       &
                  ,csite%sfcwater_tempk(k,ipa),csite%sfcwater_mass(k,ipa)                  &
                  ,csite%sfcwater_fracliq(k,ipa),csite%sfcwater_depth(k,ipa)               &
                  ,csite%rshort_s(k,ipa),csite%par_s(k,ipa)
         end do
      end if

      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(a)'  ) ' '
      return
   end subroutine print_csiteipa
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine is similar to print_csite, except that it also prints the          !
   ! outcome of the Runge-Kutta integrator.                                                !
   !---------------------------------------------------------------------------------------!
   subroutine print_rk4patch(y,csite,ipa)
      use rk4_coms              , only : rk4patchtype          & ! structure
                                       , rk4site               ! ! intent(in)
      use ed_state_vars         , only : sitetype              & ! structure
                                       , patchtype             ! ! structure
      use grid_coms             , only : nzg                   ! ! intent(in)
      use ed_misc_coms          , only : current_time          ! ! intent(in)
      use consts_coms           , only : pio1808               ! ! intent(in)
      use therm_lib8            , only : thetaeiv8             & ! function
                                       , vpdefil8              ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: y
      type(sitetype)     , target     :: csite
      integer            , intent(in) :: ipa
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer    :: cpatch
      integer                         :: k
      integer                         :: ico
      real(kind=8)                    :: y_can_rvap
      real(kind=8)                    :: y_can_theiv
      real(kind=8)                    :: y_can_vpdef
      real(kind=4)                    :: pss_lai
      real(kind=4)                    :: pss_wai
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)



      !----- Find the total patch LAI and WAI. --------------------------------------------!
      pss_lai = 0.0
      pss_wai = 0.0
      do ico=1,cpatch%ncohorts
         pss_lai = pss_lai + cpatch%lai(ico)
         pss_wai = pss_wai + cpatch%wai(ico)
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the ice-vapour equivalent potential temperature and vapour pressure      !
      ! deficit (output only).                                                             !
      !------------------------------------------------------------------------------------!
      y_can_rvap  = y%can_shv / (1.d0 - y%can_shv)
      y_can_theiv = thetaeiv8(y%can_theta,y%can_prss,y%can_temp,y_can_rvap,y_can_rvap)
      y_can_vpdef = vpdefil8 (y%can_prss,y%can_temp,y%can_shv,.true.)
      !------------------------------------------------------------------------------------!

      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(80a)') ('=',k=1,80)

      write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (rk4patch) ||||'

      write(unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,3(i2.2,a))')                              &
            'Time:',current_time%month,'/',current_time%date,'/',current_time%year         &
                   ,current_time%hour,':',current_time%min,':',current_time%sec,' UTC'
      write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
      write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(80a)')         ('-',k=1,80)
      write (unit=*,fmt='(a)')           ' ATMOSPHERIC CONDITIONS: '
      write (unit=*,fmt='(a,1x,es12.4)') ' Longitude                  : ',rk4site%lon
      write (unit=*,fmt='(a,1x,es12.4)') ' Latitude                   : ',rk4site%lat
      write (unit=*,fmt='(a,1x,es12.4)') ' Air temperature (Ref. hgt.): ',rk4site%atm_tmp
      write (unit=*,fmt='(a,1x,es12.4)') ' Air potential temp.        : ',rk4site%atm_theta
      write (unit=*,fmt='(a,1x,es12.4)') ' Air theta_Eiv              : ',rk4site%atm_theiv
      write (unit=*,fmt='(a,1x,es12.4)') ' Air vapour pres. deficit   : ',rk4site%atm_vpdef
      write (unit=*,fmt='(a,1x,es12.4)') ' H2Ov mixing ratio          : ',rk4site%atm_shv
      write (unit=*,fmt='(a,1x,es12.4)') ' CO2  mixing ratio          : ',rk4site%atm_co2
      write (unit=*,fmt='(a,1x,es12.4)') ' Pressure                   : ',rk4site%atm_prss
      write (unit=*,fmt='(a,1x,es12.4)') ' Exner function             : ',rk4site%atm_exner
      write (unit=*,fmt='(a,1x,es12.4)') ' Prescribed u*              : ',rk4site%atm_ustar
      write (unit=*,fmt='(a,1x,es12.4)') ' Height                     : ',rk4site%geoht
      write (unit=*,fmt='(a,1x,es12.4)') ' Precip. mass  flux         : ',rk4site%pcpg
      write (unit=*,fmt='(a,1x,es12.4)') ' Precip. heat  flux         : ',rk4site%qpcpg
      write (unit=*,fmt='(a,1x,es12.4)') ' Precip. depth flux         : ',rk4site%dpcpg
      write (unit=*,fmt='(a,1x,es12.4)') ' Downward SW radiation      : ',rk4site%rshort
      write (unit=*,fmt='(a,1x,es12.4)') ' Downward LW radiation      : ',rk4site%rlong
      write (unit=*,fmt='(a,1x,es12.4)') ' Zenith angle (deg)         : '                  &
                                                              ,acos(rk4site%cosz) / pio1808

      write (unit=*,fmt='(80a)') ('=',k=1,80)
      write (unit=*,fmt='(a)'  ) 'Leaf information (only those resolvable are shown): '
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),9(a12,1x))')                                            &
            '    PFT','KRDEPTH','      NPLANT','      HEIGHT','         DBH'               &
                               ,'      BDEADA','      BDEADB','       BLEAF'               &
                               ,'     FS_OPEN','         FSW','         FSN'
      do ico = 1,cpatch%ncohorts
         if (cpatch%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),9(es12.4,1x))')                                    &
                  cpatch%pft(ico),cpatch%krdepth(ico)                                      &
                 ,cpatch%nplant(ico),cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdeada(ico)   &
                 ,cpatch%bdeadb(ico),cpatch%bleaf (ico),cpatch%fs_open(ico)                &
                 ,cpatch%fsw(ico),cpatch%fsn(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),5(a12,1x))')                                            &
         '    PFT','KRDEPTH','         LAI','         GPP','   LEAF_RESP','   ROOT_RESP'   &
        ,'   STEM_RESP'
      do ico = 1,cpatch%ncohorts
         if (cpatch%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),5(es12.4,1x))')                                    &
               cpatch%pft(ico), cpatch%krdepth(ico)                                        &
              ,y%lai(ico),y%gpp(ico),y%leaf_resp(ico),y%root_resp(ico),y%stem_resp(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),10(a12,1x))')                                           &
         '    PFT','KRDEPTH','         LAI','         WAI','         TAI',' LEAF_ENERGY'   &
             ,'  LEAF_WATER','   LEAF_HCAP','LEAF_H2O_IM2','   LEAF_TEMP','   LEAF_FLIQ'   &
             ,'    LINT_SHV'
      do ico = 1,cpatch%ncohorts
         if (y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),10(es12.4,1x))')                                   &
                   cpatch%pft(ico), cpatch%krdepth(ico)                                    &
                  ,y%lai(ico),y%wai(ico),y%tai(ico),y%leaf_energy(ico),y%leaf_water(ico)   &
                  ,y%leaf_water_im2(ico),y%leaf_hcap(ico),y%leaf_temp(ico)                 &
                  ,y%leaf_fliq(ico),y%lint_shv(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                            &
                 '    PFT','KRDEPTH','         LAI','      HEIGHT','   LEAF_TEMP'          &
                     ,'    VEG_WIND','  LEAF_REYNO','LEAF_GRASHOF',' LEAF_NUFORC'          &
                     ,' LEAF_NUFREE'
      do ico = 1,cpatch%ncohorts
         if (y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))')                                    &
                   cpatch%pft(ico),cpatch%krdepth(ico)                                     &
                  ,y%lai(ico),cpatch%hite(ico),y%leaf_temp(ico),y%veg_wind(ico)            &
                  ,y%leaf_reynolds(ico),y%leaf_grashof(ico),y%leaf_nussforc(ico)           &
                  ,y%leaf_nussfree(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),7(a12,1x))')                                            &
                 '    PFT','KRDEPTH','         LAI','      HEIGHT','    LEAF_GBH'          &
                     ,'    LEAF_GBW','  GSW_CLOSED','    GSW_OPEN','     FS_OPEN'
      do ico = 1,cpatch%ncohorts
         if (y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),7(es12.4,1x))')                                    &
                   cpatch%pft(ico),cpatch%krdepth(ico)                                     &
                  ,y%lai(ico),cpatch%hite(ico),y%leaf_gbh(ico),y%leaf_gbw(ico)             &
                  ,y%gsw_closed(ico),y%gsw_open(ico),cpatch%fs_open(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),6(a12,1x))')                                            &
                 '    PFT','KRDEPTH','         LAI','      HEIGHT','    RSHORT_L'          &
                     ,'     RLONG_L','  PAR_L_BEAM','  PAR_L_DIFF'
      do ico = 1,cpatch%ncohorts
         if (y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),6(es12.4,1x))')                                    &
                   cpatch%pft(ico),cpatch%krdepth(ico)                                     &
                  ,y%lai(ico),cpatch%hite(ico),y%rshort_l(ico),y%rlong_l(ico)              &
                  ,cpatch%par_l_beam(ico),cpatch%par_l_diffuse(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('=',k=1,80)
      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(80a)') ('=',k=1,80)
      write (unit=*,fmt='(a)'  ) 'Wood information (only those resolvable are shown): '
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),9(a12,1x))')                                            &
            '    PFT','KRDEPTH','      NPLANT','      HEIGHT','         DBH'               &
                ,'      BDEADA','      BDEADB','   BSAPWOODA','   BSAPWOODB'               &
                ,'      BBARKA','      BBARKB'
      do ico = 1,cpatch%ncohorts
         if (cpatch%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),9(es12.4,1x))')                                    &
                  cpatch%pft(ico),cpatch%krdepth(ico)                                      &
                 ,cpatch%nplant(ico),cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdeada(ico)   &
                 ,cpatch%bdeadb(ico),cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico)           &
                 ,cpatch%bbarka(ico),cpatch%bbarkb(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),9(a12,1x))')                                            &
         '    PFT','KRDEPTH','         LAI','         WAI','         TAI',' WOOD_ENERGY'   &
             ,'  WOOD_WATER','WOOD_H2O_IM2','   WOOD_HCAP','   WOOD_TEMP','   WOOD_FLIQ'
      do ico = 1,cpatch%ncohorts
         if (y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),9(es12.4,1x))')                                    &
                   cpatch%pft(ico), cpatch%krdepth(ico)                                    &
                  ,y%lai(ico),y%wai(ico),y%tai(ico),y%wood_energy(ico),y%wood_water(ico)   &
                  ,y%wood_water_im2(ico),y%wood_hcap(ico),y%wood_temp(ico),y%wood_fliq(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                            &
                 '    PFT','KRDEPTH','         WAI','      HEIGHT','   LEAF_TEMP'          &
                     ,'    VEG_WIND','  WOOD_REYNO','WOOD_GRASHOF',' WOOD_NUFORC'          &
                     ,' WOOD_NUFREE'
      do ico = 1,cpatch%ncohorts
         if (y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))')                                    &
                   cpatch%pft(ico),cpatch%krdepth(ico)                                     &
                  ,y%wai(ico),cpatch%hite(ico),y%wood_temp(ico),y%veg_wind(ico)            &
                  ,y%wood_reynolds(ico),y%wood_grashof(ico),y%wood_nussforc(ico)           &
                  ,y%wood_nussfree(ico)
         end if
      end do
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(2(a7,1x),6(a12,1x))')                                            &
                 '    PFT','KRDEPTH','         WAI','      HEIGHT','    WOOD_GBH'          &
                     ,'    WOOD_GBW','    RSHORT_W','     RLONG_W'
      do ico = 1,cpatch%ncohorts
         if (y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i7,1x),6(es12.4,1x))')                                    &
                   cpatch%pft(ico),cpatch%krdepth(ico)                                     &
                  ,y%wai(ico),cpatch%hite(ico),y%wood_gbh(ico),y%wood_gbw(ico)             &
                  ,y%rshort_w(ico),y%rlong_w(ico) 
         end if
      end do
      write (unit=*,fmt='(80a)') ('=',k=1,80)
      write (unit=*,fmt='(a)'  ) ' '
      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(10(a12,1x))')   '  VEG_HEIGHT','   VEG_ROUGH','VEG_DISPLACE'     &
                                         ,'   PATCH_LAI','   PATCH_WAI','   CAN_DEPTH'     &
                                         ,'    CAN_RHOS','    CAN_DMOL','    CAN_PRSS'     &
                                         ,'       GGNET'

      write (unit=*,fmt='(10(es12.4,1x))') y%veg_height,y%veg_rough,y%veg_displace         &
                                          ,pss_lai,pss_wai,y%can_depth,y%can_rhos          &
                                          ,y%can_dmol,y%can_prss,y%ggnet
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(10(a12,1x))') '     CAN_CO2','   CAN_THEIV','   CAN_THETA'       &
                                       ,'    CAN_TEMP','     CAN_SHV','     CAN_SSH'       &
                                       ,'    CAN_RVAP','   CAN_VPDEF','     CAN_RHV'       &
                                       ,'CAN_ENTHALPY'

      write (unit=*,fmt='(10(es12.4,1x))')  y%can_co2      , y_can_theiv    , y%can_theta  &
                                          , y%can_temp     , y%can_shv      , y%can_ssh    &
                                          , y_can_rvap     , y_can_vpdef    , y%can_rhv    &
                                          , y%can_enthalpy

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'       &
                                       ,'       TSTAR','       ESTAR','        ZETA'       &
                                       ,'     RI_BULK'
      write (unit=*,fmt='(7(es12.4,1x))') y%ustar,y%qstar,y%cstar,y%tstar,y%estar,y%zeta   &
                                         ,y%ribulk

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(9(a12,1x))')  '      FGC_RH','      FSC_RH','     STGC_RH'       &
                                       ,'     STSC_RH','      MSC_RH','      SSC_RH'       &
                                       ,'      PSC_RH','STORAGE_RESP',' GROWTH_RESP'
      write (unit=*,fmt='(9(es12.4,1x))') y%fgc_rh,y%fsc_rh,y%stgc_rh,y%stsc_rh,y%msc_rh   &
                                         ,y%ssc_rh,y%psc_rh,y%commit_storage_resp          &
                                         ,y%commit_growth_resp

      write (unit=*,fmt='(80a)') ('-',k=1,80)


      write (unit=*,fmt='(5(a12,1x))')  '   FLAG_SFCW',' VIRT_ENERGY','  VIRT_WATER'       &
                                       ,'VIRTUAL_TEMP','VIRTUAL_FLIQ'
      write (unit=*,fmt='(i12,1x,4(es12.4,1x))') y%flag_sfcwater,y%virtual_energy          &
                                                ,y%virtual_water,y%virtual_tempk           &
                                                ,y%virtual_fracliq
      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(4(a12,1x))')    '  GROUND_SHV','  GROUND_SSH',' GROUND_TEMP'     &
                                         ,' GROUND_FLIQ'
      write (unit=*,fmt='(4(es12.4,1x))') y%ground_shv, y%ground_ssh, y%ground_temp        &
                                         ,y%ground_fliq

      write (unit=*,fmt='(80a)') ('-',k=1,80)

      write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'       &
                                      &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
      do k=rk4site%lsl,nzg
         write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,rk4site%ntext_soil(k)             &
               ,y%soil_energy(k),y%soil_tempk(k),y%soil_water(k),y%soil_fracliq(k)
      end do
      
      if (csite%nlev_sfcwater(ipa) >= 1) then
         write (unit=*,fmt='(80a)') ('-',k=1,80)
         write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'    &
                                         &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'
         do k=1,csite%nlev_sfcwater(ipa)
            write (unit=*,fmt='(i5,1x,5(es12.4,1x))') k,y%sfcwater_energy(k)               &
                  ,y%sfcwater_tempk(k),y%sfcwater_mass(k),y%sfcwater_fracliq(k)            &
                  ,y%sfcwater_depth(k)
         end do
      end if

      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(80a)') ('=',k=1,80)
      write(unit=*,fmt='(a)'  ) ' '

      !----- Printing the corresponding patch information (with some redundancy) ----------!
      call print_csiteipa(csite, ipa)

      call fatal_error('IFLAG1 problem. The model didn''t converge!','print_rk4patch'&
                    &,'rk4_integ_utils.f90')
      return
   end subroutine print_rk4patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine prints the full state of a given patch, for full debugging       !
   ! purposes.  This will create one file for each patch.  This sub-routine will not print !
   ! the temperature of each cohort, instead it will just compute the average.             !
   !---------------------------------------------------------------------------------------!
   subroutine print_rk4_state(initp,fluxp,csite,ipa,isi,elapsed,hdid)
      use consts_coms  , only : t3ple8        ! ! intent(in)
      use ed_max_dims  , only : str_len       ! ! intent(in)
      use ed_misc_coms , only : current_time  ! ! intent(in)
      use ed_state_vars, only : sitetype      & ! structure
                              , patchtype     ! ! structure
      use grid_coms    , only : nzg           ! ! intent(in)
      use rk4_coms     , only : rk4patchtype  & ! structure
                              , rk4site       & ! intent(in)
                              , detail_pref   ! ! intent(in)
      use therm_lib8   , only : uextcm2tl8    & ! sub-routine
                              , thetaeiv8     & ! function
                              , vpdefil8      ! ! function
      use soil_coms    , only : soil8         ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(rk4patchtype)    , target     :: initp
      type(rk4patchtype)    , target     :: fluxp
      type(sitetype)        , target     :: csite
      integer               , intent(in) :: ipa
      integer               , intent(in) :: isi
      real(kind=8)          , intent(in) :: elapsed
      real(kind=8)          , intent(in) :: hdid
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)       , pointer    :: cpatch
      character(len=str_len)             :: detail_fout
      integer                            :: nsoil
      integer                            :: ico
      integer                            :: leaf_resolve
      integer                            :: wood_resolve
      integer                            :: is_small
      logical                            :: isthere
      real(kind=8)                       :: sum_leaf_energy
      real(kind=8)                       :: sum_leaf_water
      real(kind=8)                       :: sum_leaf_water_im2
      real(kind=8)                       :: sum_leaf_hcap
      real(kind=8)                       :: sum_wood_energy
      real(kind=8)                       :: sum_wood_water_im2
      real(kind=8)                       :: sum_wood_water
      real(kind=8)                       :: sum_wood_hcap
      real(kind=8)                       :: sum_gpp
      real(kind=8)                       :: sum_lai
      real(kind=8)                       :: sum_wai
      real(kind=8)                       :: sum_plresp
      real(kind=8)                       :: qintercepted
      real(kind=8)                       :: avg_leaf_temp
      real(kind=8)                       :: avg_leaf_fliq
      real(kind=8)                       :: avg_wood_temp
      real(kind=8)                       :: avg_wood_fliq
      real(kind=8)                       :: par_b_beam
      real(kind=8)                       :: par_b_diff
      real(kind=8)                       :: nir_b_beam
      real(kind=8)                       :: nir_b_diff
      real(kind=8)                       :: elapsec
      real(kind=8)                       :: can_rvap
      real(kind=8)                       :: can_theiv
      real(kind=8)                       :: can_vpdef
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter :: phfmt='(91(a,1x))'
      character(len=48), parameter ::                                                      &
                                   pbfmt='(3(i13,1x),4(es13.6,1x),3(i13,1x),81(es13.6,1x))'
      character(len=10), parameter :: chfmt='(58(a,1x))'
      character(len=48), parameter ::                                                      &
                                   cbfmt='(3(i13,1x),2(es13.6,1x),4(i13,1x),49(es13.6,1x))'
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Old files are now deleted by sub-routine initialize_misc_stepvars.  This output !
      ! is extremely large, so think twice before turning it for multiple polygons, and    !
      ! for long-term simulations.                                                         !
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     First we loop over all cohorts and add the vegetation energy and water.        !
      !------------------------------------------------------------------------------------!
      sum_leaf_energy     = 0.d0
      sum_leaf_water      = 0.d0
      sum_leaf_water_im2  = 0.d0
      sum_leaf_hcap       = 0.d0
      sum_wood_energy     = 0.d0
      sum_wood_water      = 0.d0
      sum_wood_water_im2  = 0.d0
      sum_wood_hcap       = 0.d0
      sum_gpp             = 0.d0
      sum_plresp          = initp%commit_storage_resp + initp%commit_growth_resp
      sum_lai             = 0.d0
      sum_wai             = 0.d0
      cpatch => csite%patch(ipa)
      do ico=1,cpatch%ncohorts
         if (initp%leaf_resolvable(ico)) then
            !----- Integrate vegetation properties using m2gnd rather than plant. ---------!
            sum_leaf_energy    = sum_leaf_energy    + initp%leaf_energy   (ico)
            sum_leaf_water     = sum_leaf_water     + initp%leaf_water    (ico)
            sum_leaf_water_im2 = sum_leaf_water_im2 + initp%leaf_water_im2(ico)
            sum_leaf_hcap      = sum_leaf_hcap      + initp%leaf_hcap     (ico)
            sum_gpp            = sum_gpp            + initp%gpp           (ico)
            sum_plresp         = sum_plresp         + initp%leaf_resp     (ico)            &
                                                    + initp%root_resp     (ico)            &
                                                    + initp%stem_resp     (ico)
         end if
         if (initp%wood_resolvable(ico)) then
            !----- Integrate vegetation properties using m2gnd rather than plant. ---------!
            sum_wood_energy    = sum_wood_energy    + initp%wood_energy   (ico)
            sum_wood_water     = sum_wood_water     + initp%wood_water    (ico)
            sum_wood_water_im2 = sum_wood_water_im2 + initp%wood_water_im2(ico)
            sum_wood_hcap      = sum_wood_hcap      + initp%wood_hcap     (ico)
         end if

         !----- TAI.  We integrate all cohorts, including those that we skip. -------------!
         sum_lai = sum_lai + initp%lai(ico)
         sum_wai = sum_wai + initp%wai(ico)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the ice-vapour equivalent potential temperature of the canopy air space.  !
      !------------------------------------------------------------------------------------!
      can_rvap   = initp%can_shv / (1.d0 - initp%can_shv)
      can_theiv  = thetaeiv8( initp%can_theta , initp%can_prss  , initp%can_temp           &
                            , can_rvap        , can_rvap        )
      can_vpdef  = vpdefil8 ( initp%can_prss  , initp%can_temp  , initp%can_shv , .true. )
      !------------------------------------------------------------------------------------!

      par_b_beam = dble(csite%par_b_beam   (ipa))
      par_b_diff = dble(csite%par_b_diffuse(ipa))
      nir_b_beam = dble(csite%nir_b_beam   (ipa))
      nir_b_diff = dble(csite%nir_b_diffuse(ipa))


      !------------------------------------------------------------------------------------!
      !     Then we find the average leaf temperature.  If none of the cohorts were        !
      ! solved, or if there is no vegetation, we assign the canopy air temperature.        !
      !------------------------------------------------------------------------------------!
      if (sum_leaf_energy == 0.d0) then
         avg_leaf_temp = initp%can_temp
         if (initp%can_temp == t3ple8) then
            avg_leaf_fliq = 5.d-1
         elseif (initp%can_temp > t3ple8) then
            avg_leaf_fliq = 1.d0
         else
            avg_leaf_fliq = 0.d0
         end if
      else
         call uextcm2tl8(sum_leaf_energy,sum_leaf_water+sum_leaf_water_im2,sum_leaf_hcap   &
                        ,avg_leaf_temp,avg_leaf_fliq)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Then we find the average wood temperature.  If none of the cohorts were        !
      ! solved, or if there is no vegetation, we assign the canopy air temperature.        !
      !------------------------------------------------------------------------------------!
      if (sum_wood_energy == 0.d0) then
         avg_wood_temp = initp%can_temp
         if (initp%can_temp == t3ple8) then
            avg_wood_fliq = 5.d-1
         elseif (initp%can_temp > t3ple8) then
            avg_wood_fliq = 1.d0
         else
            avg_wood_fliq = 0.d0
         end if
      else
         call uextcm2tl8(sum_wood_energy,sum_wood_water+sum_wood_water_im2,sum_wood_hcap   &
                        ,avg_wood_temp,avg_wood_fliq)
      end if
      !------------------------------------------------------------------------------------!


      !----- Compute the hour as elapsed seconds since midnight. --------------------------!
      elapsec = dble(current_time%time) + elapsed
      !------------------------------------------------------------------------------------!


      !----- Find the soil type of the top layer. -----------------------------------------!
      nsoil   = rk4site%ntext_soil(nzg)
      !------------------------------------------------------------------------------------!

      !----- Create the file name. --------------------------------------------------------!
      write (detail_fout,fmt='(2a,2(i4.4,a))')                                             &
                                    trim(detail_pref),'prk4_site_',isi,'_patch_',ipa,'.txt'
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Check whether the file exists or not.  In case it doesn't, create it and add    !
      ! the header.                                                                        !
      !------------------------------------------------------------------------------------!
      inquire(file=trim(detail_fout),exist=isthere)
      if (.not. isthere) then
         open  (unit=83,file=trim(detail_fout),status='replace',action='write')
         write (unit=83,fmt=phfmt)  '         YEAR' , '        MONTH', '          DAY'     &
                                  , '         TIME' , '         HDID', '          LAI'     &
                                  , '          WAI' , '          KSN', 'FLAG.SFCWATER'     &
                                  , '  FLAG.WFLXGC' , '     ATM.PRSS', '     ATM.TEMP'     &
                                  , '      ATM.SHV' , '      ATM.CO2', '     ATM.VELS'     &
                                  , '    ATM.PRATE' , '   ATM.HEIGHT', '     ATM.RHOS'     &
                                  , '   ATM.RELHUM' , '    ATM.THETA', '    ATM.THEIV'     &
                                  , '    ATM.VPDEF' , '   MET.RSHORT', '    MET.RLONG'     &
                                  , '     CAN.PRSS' , '     CAN.TEMP', '      CAN.SHV'     &
                                  , '      CAN.CO2' , '    CAN.DEPTH', '     CAN.RHOS'     &
                                  , '     CAN.DMOL', '   CAN.RELHUM' , '    CAN.THETA'     &
                                  , '    CAN.THEIV', ' CAN.ENTHALPY' , '     SFC.TEMP'     &
                                  , '      SFC.SHV', '    LEAF.TEMP' , '   LEAF.WATER'     &
                                  , '    WOOD.TEMP', '   WOOD.WATER' , '       GGBARE'     &
                                  , '        GGVEG', '        GGNET' , '      OPENCAN'     &
                                  , '    SOIL.TEMP', '   SOIL.WATER' , '       SOILCP'     &
                                  , '       SOILWP', '       SOILFC' , '       SOILBP'     &
                                  , '        USTAR', '        TSTAR' , '        QSTAR'     &
                                  , '        CSTAR', '         ZETA' , '      RI.BULK'     &
                                  , '   GND.RSHORT', '      GND.PAR' , '    GND.RLONG'     &
                                  , '       WFLXLC', '       WFLXWC' , '       WFLXGC'     &
                                  , '       WFLXAC', '       TRANSP' , '        WSHED'     &
                                  , '    INTERCEPT', '  THROUGHFALL' , '       HFLXGC'     &
                                  , '       HFLXLC', '       HFLXWC' , '       HFLXAC'     &
                                  , '       CFLXAC', '       CFLXST' , '       FGC.RH'     &
                                  , '       FSC.RH', '      STGC.RH' , '      STSC.RH'     &
                                  , '       MSC.RH', '       SSC.RH' , '       PSC.RH'     &
                                  , '          GPP', '       PLRESP' , ' PAR.BEAM.TOP'     &
                                  , ' PAR.DIFF.TOP', ' NIR.BEAM.TOP' , ' NIR.DIFF.TOP'     &
                                  , ' PAR.BEAM.BOT', ' PAR.DIFF.BOT' , ' NIR.BEAM.BOT'     &
                                  , ' NIR.DIFF.BOT'
                                  
         close (unit=83,status='keep')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Re-open the file at the last line, and include the current status.             !
      !------------------------------------------------------------------------------------!
      open (unit=83,file=trim(detail_fout),status='old',action='write',position='append')
      write(unit=83,fmt=pbfmt)                                                                &
                     current_time%year     , current_time%month    , current_time%date     &
                   , elapsec               , hdid                  , sum_lai               &
                   , sum_wai               , initp%nlev_sfcwater   , initp%flag_sfcwater   &
                   , initp%flag_wflxgc     , rk4site%atm_prss      , rk4site%atm_tmp       &
                   , rk4site%atm_shv       , rk4site%atm_co2       , initp%vels            &
                   , rk4site%pcpg          , rk4site%geoht         , rk4site%atm_rhos      &
                   , rk4site%atm_rhv       , rk4site%atm_theta     , rk4site%atm_theiv     &
                   , rk4site%atm_vpdef     , rk4site%rshort        , rk4site%rlong         &
                   , initp%can_prss        , initp%can_temp        , initp%can_shv         &
                   , initp%can_co2         , initp%can_depth       , initp%can_rhos        &
                   , initp%can_dmol        , initp%can_rhv         , initp%can_theta       &
                   , can_theiv             , initp%can_enthalpy    , initp%ground_temp     &
                   , initp%ground_shv      , avg_leaf_temp         , sum_leaf_water        &
                   , avg_wood_temp         , sum_wood_water        , initp%ggbare          &
                   , initp%ggveg           , initp%ggnet           , initp%opencan_frac    &
                   , initp%soil_tempk(nzg) , initp%soil_water(nzg) , soil8(nsoil)%soilcp   &
                   , soil8(nsoil)%soilwp   , soil8(nsoil)%sfldcap  , soil8(nsoil)%soilbp   &
                   , initp%ustar           , initp%tstar           , initp%qstar           &
                   , initp%cstar           , initp%zeta            , initp%ribulk          &
                   , fluxp%flx_rshort_gnd  , fluxp%flx_par_gnd     , fluxp%flx_rlong_gnd   &
                   , fluxp%flx_vapor_lc    , fluxp%flx_vapor_wc    , fluxp%flx_vapor_gc    &
                   , fluxp%flx_vapor_ac    , fluxp%flx_transp      , fluxp%flx_wshed_vg    &
                   , fluxp%flx_intercepted , fluxp%flx_throughfall , fluxp%flx_sensible_gc &
                   , fluxp%flx_sensible_lc , fluxp%flx_sensible_wc , fluxp%flx_sensible_ac &
                   , fluxp%flx_carbon_ac   , fluxp%flx_carbon_st   , initp%fgc_rh          &
                   , initp%fsc_rh          , initp%stgc_rh         , initp%stsc_rh         &
                   , initp%msc_rh          , initp%ssc_rh          , initp%psc_rh          &
                   , sum_gpp               , sum_plresp            , rk4site%par_beam      &
                   , rk4site%par_diffuse   , rk4site%nir_beam      , rk4site%nir_diffuse   &
                   , par_b_beam            , par_b_diff            , nir_b_beam            &
                   , nir_b_diff
                   
      close(unit=83,status='keep')
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Now it is time to check and output the cohort-level files.                     !
      !------------------------------------------------------------------------------------!
      do ico=1, cpatch%ncohorts

         !----- Find the integer version of logical flags. --------------------------------!
         if (initp%leaf_resolvable(ico)) then
            leaf_resolve = 1
         else
            leaf_resolve = 0
         end if
         if (initp%wood_resolvable(ico)) then
            wood_resolve = 1
         else
            wood_resolve = 0
         end if
         if (initp%is_small(ico)) then
            is_small = 1
         else
            is_small = 0
         end if
         !---------------------------------------------------------------------------------!

         !----- Copy intercepted water. ---------------------------------------------------!
         qintercepted = fluxp%cfx_qintercepted(ico)
         !---------------------------------------------------------------------------------!

         !----- Create the file name. -----------------------------------------------------!
         write (detail_fout,fmt='(2a,3(i4.4,a))')                                          &
                     trim(detail_pref),'crk4_site_',isi,'_patch_',ipa,'_cohort_',ico,'.txt'
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Check whether the file exists or not.  In case it doesn't, create it and add !
         ! the header.                                                                     !
         !---------------------------------------------------------------------------------!
         inquire(file=trim(detail_fout),exist=isthere)
         if (.not. isthere) then
            open  (unit=84,file=trim(detail_fout),status='replace',action='write')
            write (unit=84,fmt=chfmt)                                                      &
                            '             YEAR', '            MONTH', '              DAY'  &
                          , '             TIME', '             HDID', '              PFT'  &
                          , '     LEAF_RESOLVE', '     WOOD_RESOLVE', '         IS_SMALL'  &
                          , '           NPLANT', '           HEIGHT', '              LAI'  &
                          , '              WAI', '       CROWN_AREA', '      LEAF_ENERGY'  &
                          , '       LEAF_WATER', '   LEAF_WATER_IM2', '        LEAF_HCAP'  &
                          , '        LEAF_TEMP', '        LEAF_FLIQ', '      WOOD_ENERGY'  &
                          , '       WOOD_WATER', '   WOOD_WATER_IM2', '        WOOD_HCAP'  &
                          , '        WOOD_TEMP', '        WOOD_FLIQ', '         VEG_WIND'  &
                          , '          FS_OPEN', '       LEAF_REYNO', '     LEAF_GRASHOF'  &
                          , '      LEAF_NUFREE', '      LEAF_NUFORC', '       WOOD_REYNO'  &
                          , '     WOOD_GRASHOF', '      WOOD_NUFREE', '      WOOD_NUFORC'  &
                          , '         LINT_SHV', '         LEAF_GBH', '         LEAF_GBW'  &
                          , '         WOOD_GBH', '         WOOD_GBW', '         GSW_OPEN'  &
                          , '         GSW_CLOS', '              GPP', '        LEAF_RESP'  &
                          , '        ROOT_RESP', '        STEM_RESP', '         RSHORT_L'  &
                          , '          RLONG_L', '         RSHORT_W', '          RLONG_W'  &
                          , '           HFLXLC', '           HFLXWC', '          QWFLXLC'  &
                          , '          QWFLXWC', '           QWSHED', '          QTRANSP'  &
                          , '     QINTERCEPTED'
            close (unit=84,status='keep')
         end if
         !---------------------------------------------------------------------------------!
         



         !---------------------------------------------------------------------------------!
         !     Re-open the file at the last line, and include the current status.          !
         !---------------------------------------------------------------------------------!
         open (unit=84,file=trim(detail_fout),status='old',action='write',position='append')
         write(unit=84,fmt=cbfmt)                                                          &
                        current_time%year             , current_time%month                 &
                      , current_time%date             , elapsec                            &
                      , hdid                          , cpatch%pft(ico)                    &
                      , leaf_resolve                  , wood_resolve                       &
                      , is_small                      , initp%nplant(ico)                  &
                      , cpatch%hite(ico)              , initp%lai(ico)                     &
                      , initp%wai(ico)                , initp%crown_area(ico)              &
                      , initp%leaf_energy(ico)        , initp%leaf_water(ico)              &
                      , initp%leaf_water_im2(ico)     , initp%leaf_hcap(ico)               &
                      , initp%leaf_temp(ico)          , initp%leaf_fliq(ico)               &
                      , initp%wood_energy(ico)        , initp%wood_water(ico)              &
                      , initp%wood_water_im2(ico)     , initp%wood_hcap(ico)               &
                      , initp%wood_temp(ico)          , initp%wood_fliq(ico)               &
                      , initp%veg_wind(ico)           , initp%fs_open(ico)                 &
                      , initp%leaf_reynolds(ico)      , initp%leaf_grashof(ico)            &
                      , initp%leaf_nussfree(ico)      , initp%leaf_nussforc(ico)           &
                      , initp%wood_reynolds(ico)      , initp%wood_grashof(ico)            &
                      , initp%wood_nussfree(ico)      , initp%wood_nussforc(ico)           &
                      , initp%lint_shv(ico)           , initp%leaf_gbh(ico)                &
                      , initp%leaf_gbw(ico)           , initp%wood_gbh(ico)                &
                      , initp%wood_gbw(ico)           , initp%gsw_open(ico)                &
                      , initp%gsw_closed(ico)         , initp%gpp(ico)                     &
                      , initp%leaf_resp(ico)          , initp%root_resp(ico)               &
                      , initp%stem_resp(ico)          , initp%rshort_l(ico)                &
                      , initp%rlong_l(ico)            , initp%rshort_w(ico)                &
                      , initp%rlong_w(ico)            , fluxp%cfx_hflxlc(ico)              &
                      , fluxp%cfx_hflxwc(ico)         , fluxp%cfx_qwflxlc(ico)             &
                      , fluxp%cfx_qwflxwc(ico)        , fluxp%cfx_qwshed(ico)              &
                      , fluxp%cfx_qtransp(ico)        , qintercepted
         close(unit=84,status='keep')
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
      return
   end subroutine print_rk4_state
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks whether the leaf and wood properties are consistent.  Any !
   ! update on long-term dynamics must ensure that updates on internal energy and heat     !
   ! capacity results in the same temperature.  This check is only needed in case of       !
   ! updates in the long-term dynamics.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine sanity_check_veg_energy(csite,ipa)
      use ed_state_vars          , only : sitetype             & ! structure
                                        , patchtype            ! ! structure
      use consts_coms            , only : tiny_num             ! ! intent(in)
      use therm_lib              , only : uextcm2tl            ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)            , target      :: csite
      integer                   , intent(in)  :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)           , pointer     :: cpatch
      integer                                 :: ico
      real                                    :: test_leaf_temp
      real                                    :: test_leaf_fliq
      real                                    :: test_wood_temp
      real                                    :: test_wood_fliq
      logical                                 :: fine_leaf_temp
      logical                                 :: fine_leaf_fliq
      logical                                 :: fine_wood_temp
      logical                                 :: fine_wood_fliq
      integer                                 :: n
      integer                                 :: nproblem
      !----- Local constants. -------------------------------------------------------------!
      character(len=13)         , parameter   :: efmt       = '(a,1x,es12.5)'
      character(len=9)          , parameter   :: ifmt       = '(a,1x,i5)'
      character(len=9)          , parameter   :: lfmt       = '(a,1x,l1)'
      real                      , parameter   :: fine_toler = 0.01
      !------------------------------------------------------------------------------------!


      !----- Current patch. ---------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Check each cohort.  Print all cohorts that may have problem before crashing.  !
      !------------------------------------------------------------------------------------!
      nproblem = 0
      do ico=1,cpatch%ncohorts
         !----- Check leaf thermodynamics. ------------------------------------------------!
         if (cpatch%leaf_resolvable(ico)) then
            call uextcm2tl(cpatch%leaf_energy(ico)                                         &
                          ,cpatch%leaf_water(ico)+cpatch%leaf_water_im2(ico)               &
                          ,cpatch%leaf_hcap(ico),test_leaf_temp,test_leaf_fliq)
            fine_leaf_temp = abs(cpatch%leaf_temp(ico) - test_leaf_temp) <= fine_toler
            fine_leaf_fliq = abs(cpatch%leaf_fliq(ico) - test_leaf_fliq) <= fine_toler     &
                             .or. cpatch%leaf_water(ico) <= tiny_num
         else
            fine_leaf_temp = .true.
            fine_leaf_fliq = .true.
         end if
         !---------------------------------------------------------------------------------!


         !----- Check wood thermodynamics. ------------------------------------------------!
         if (cpatch%wood_resolvable(ico)) then
            call uextcm2tl(cpatch%wood_energy(ico)                                         &
                          ,cpatch%wood_water(ico)+cpatch%wood_water_im2(ico)               &
                          ,cpatch%wood_hcap(ico),test_wood_temp,test_wood_fliq)
            fine_wood_temp = abs(cpatch%wood_temp(ico) - test_wood_temp) <= fine_toler
            fine_wood_fliq = abs(cpatch%wood_fliq(ico) - test_wood_fliq) <= fine_toler     &
                             .or. cpatch%wood_water(ico) <= tiny_num
         else
            fine_wood_temp = .true.
            fine_wood_fliq = .true.
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Print information if anything is wrong.                                     !
         !---------------------------------------------------------------------------------!
         if ( .not. ( fine_leaf_temp .and. fine_leaf_fliq .and.                            &
                      fine_wood_temp .and. fine_wood_fliq ) ) then
            nproblem = nproblem + 1


            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(92a)') ('=',n=1,92)
            write (unit=*,fmt='(92a)') ('=',n=1,92)
            write (unit=*,fmt='(a)'  ) ' Energy/temperature inconsistency detected!!'
            write (unit=*,fmt='(92a)') ('-',n=1,92)
            write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa
            write (unit=*,fmt=ifmt   ) ' + ILU              =',csite%dist_type   (ipa)
            write (unit=*,fmt=efmt   ) ' + AGE              =',csite%age         (ipa)

            write (unit=*,fmt='(a)'  ) ' '
            write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
            write (unit=*,fmt=ifmt   ) ' + PFT              =',cpatch%pft        (ico)
            write (unit=*,fmt=efmt   ) ' + DBH              =',cpatch%dbh        (ico)
            write (unit=*,fmt=efmt   ) ' + HEIGHT           =',cpatch%hite       (ico)

            write (unit=*,fmt='(a)'  ) ' '
            write (unit=*,fmt=lfmt   ) ' + LEAF_RESOLVABLE  =',cpatch%leaf_resolvable(ico)
            write (unit=*,fmt=efmt   ) ' + LAI              =',cpatch%lai            (ico)
            write (unit=*,fmt=efmt   ) ' + ELONGF           =',cpatch%elongf         (ico)
            write (unit=*,fmt=efmt   ) ' + LEAF_ENERGY      =',cpatch%leaf_energy    (ico)
            write (unit=*,fmt=efmt   ) ' + LEAF_WATER       =',cpatch%leaf_water     (ico)
            write (unit=*,fmt=efmt   ) ' + LEAF_WATER_IM2   =',cpatch%leaf_water_im2 (ico)
            write (unit=*,fmt=efmt   ) ' + LEAF_HCAP        =',cpatch%leaf_hcap      (ico)
            write (unit=*,fmt=lfmt   ) ' + FINE_LEAF_TEMP   =',fine_leaf_temp
            write (unit=*,fmt=efmt   ) ' + LEAF_TEMP_MEMORY =',cpatch%leaf_temp      (ico)
            write (unit=*,fmt=efmt   ) ' + LEAF_TEMP_TEST   =',test_leaf_temp
            write (unit=*,fmt=lfmt   ) ' + FINE_LEAF_FLIQ   =',fine_leaf_fliq
            write (unit=*,fmt=efmt   ) ' + LEAF_FLIQ_MEMORY =',cpatch%leaf_fliq      (ico)
            write (unit=*,fmt=efmt   ) ' + LEAF_FLIQ_TEST   =',test_leaf_fliq


            write (unit=*,fmt='(a)'  ) ' '
            write (unit=*,fmt=lfmt   ) ' + WOOD_RESOLVABLE  =',cpatch%wood_resolvable(ico)
            write (unit=*,fmt=efmt   ) ' + WAI              =',cpatch%wai            (ico)
            write (unit=*,fmt=efmt   ) ' + WOOD_ENERGY      =',cpatch%wood_energy    (ico)
            write (unit=*,fmt=efmt   ) ' + WOOD_WATER       =',cpatch%wood_water     (ico)
            write (unit=*,fmt=efmt   ) ' + WOOD_WATER_IM2   =',cpatch%wood_water_im2 (ico)
            write (unit=*,fmt=efmt   ) ' + WOOD_HCAP        =',cpatch%wood_hcap      (ico)
            write (unit=*,fmt=lfmt   ) ' + FINE_WOOD_TEMP   =',fine_wood_temp
            write (unit=*,fmt=efmt   ) ' + WOOD_TEMP_MEMORY =',cpatch%wood_temp      (ico)
            write (unit=*,fmt=efmt   ) ' + WOOD_TEMP_TEST   =',test_wood_temp
            write (unit=*,fmt=lfmt   ) ' + FINE_WOOD_FLIQ   =',fine_wood_fliq
            write (unit=*,fmt=efmt   ) ' + WOOD_FLIQ_MEMORY =',cpatch%wood_fliq      (ico)
            write (unit=*,fmt=efmt   ) ' + WOOD_FLIQ_TEST   =',test_wood_fliq
            write (unit=*,fmt='(92a)') ('=',n=1,92)
            write (unit=*,fmt='(92a)') ('=',n=1,92)
            write (unit=*,fmt='(a)') ' '
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!

      !----- Stop in case there is a problem. ---------------------------------------------!
      if (nproblem > 0) then
         call fatal_error('Long-term dynamics is not updating leaf/wood energy correctly!' &
                         ,'sanity_check_veg_energy','rk4_misc.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine sanity_check_veg_energy
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine sets ED2 acceptable range for several variables (only those whose  !
   ! range dynamically varies).                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine find_derived_thbounds(ibuff,cpatch,can_theta,can_temp,can_shv,can_prss       &
                                   ,can_depth)
      use grid_coms      , only : nzg                ! ! intent(in)
      use rk4_coms       , only : rk4aux             & ! structure
                                , rk4site            & ! structure
                                , rk4eps             & ! intent(in)
                                , rk4max_can_temp    & ! intent(in)
                                , rk4min_can_temp    & ! intent(in)
                                , rk4min_can_shv     & ! intent(in)
                                , rk4max_can_shv     & ! intent(in)
                                , rk4min_can_co2     & ! intent(in)
                                , rk4max_can_co2     & ! intent(in)
                                , print_thbnd        & ! intent(in)
                                , thbnds_fout        ! ! intent(in)
      use consts_coms    , only : hr_sec             & ! intent(in)
                                , min_sec            & ! intent(in)
                                , huge_num8          ! ! intent(in)
      use ed_state_vars  , only : patchtype          ! ! structure
      use therm_lib8     , only : press2exner8       & ! function
                                , extemp2theta8      & ! function
                                , tq2enthalpy8       & ! function
                                , thetaeiv8          & ! function
                                , thetaeivs8         & ! function
                                , idealdenssh8       & ! function
                                , idealdmolsh8       & ! function
                                , reducedpress8      ! ! function
      use soil_coms      , only : soil8              ! ! intent(in)
      use ed_misc_coms   , only : current_time       ! ! intent(in)
      use pft_coms       , only : leaf_psi_min       & ! intent(in)
                                , wood_psi_min       & ! intent(in)
                                , small_psi_min      ! ! intent(in)
      use physiology_coms, only : plant_hydro_scheme ! ! intent(in)
      use plant_hydro    , only : om_buff            & ! intent(in)
                                , psi2twe            ! ! subroutine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)             , target     :: cpatch
      integer                     , intent(in) :: ibuff
      real(kind=8)                , intent(in) :: can_theta
      real(kind=8)                , intent(in) :: can_temp
      real(kind=8)                , intent(in) :: can_shv
      real(kind=8)                , intent(in) :: can_prss
      real(kind=8)                , intent(in) :: can_depth
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                             :: can_exner_try
      real(kind=8)                             :: can_theta_try
      real(kind=8)                             :: can_enthalpy_try
      real(kind=4)                             :: leaf_water_im2_4
      real(kind=4)                             :: wood_water_im2_4
      integer                                  :: k
      integer                                  :: hour
      integer                                  :: minute
      integer                                  :: second
      integer                                  :: nsoil
      integer                                  :: ico
      integer                                  :: ipft
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     File header is now written by initialize_misc_stepvars.                        !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the bounds for pressure.  To avoid the pressure range to be too relaxed,  !
      ! switch one of the dependent variable a time, and use the current values for the    !
      ! other.  In addition, we force pressure to be bounded between the reduced pressure  !
      ! in case the reference height was different by the order of 10%.                    !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4aux(ibuff)%rk4min_can_prss =                                                      &
            reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv               &
                         ,9.d-1*rk4site%geoht,can_theta,can_shv,can_depth)
      rk4aux(ibuff)%rk4max_can_prss =                                                      &
            reducedpress8(rk4site%atm_prss,rk4site%atm_theta,rk4site%atm_shv               &
                         ,1.1d0*rk4site%geoht,can_theta,can_shv,can_depth)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the bounds for potential temperature.  To avoid the pressure range to be  !
      ! too relaxed, switch one of the dependent variable a time, and use the current      !
      ! values for the other.                                                              !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4aux(ibuff)%rk4min_can_theta   =  huge(1.d0)
      rk4aux(ibuff)%rk4max_can_theta   = -huge(1.d0)
      !----- 2. Minimum temperature. ------------------------------------------------------!
      can_exner_try      = press2exner8(can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,rk4min_can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !----- 3. Maximum temperature. ------------------------------------------------------!
      can_exner_try      = press2exner8(can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,rk4max_can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !----- 4. Minimum pressure. ---------------------------------------------------------!
      can_exner_try      = press2exner8(rk4aux(ibuff)%rk4min_can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !----- 5. Maximum pressure. ---------------------------------------------------------!
      can_exner_try      = press2exner8(rk4aux(ibuff)%rk4max_can_prss)
      can_theta_try      = extemp2theta8(can_exner_try,can_temp)
      rk4aux(ibuff)%rk4min_can_theta   = min(rk4aux(ibuff)%rk4min_can_theta,can_theta_try)
      rk4aux(ibuff)%rk4max_can_theta   = max(rk4aux(ibuff)%rk4max_can_theta,can_theta_try)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Minimum and maximum enthalpy.                                                 !
      !------------------------------------------------------------------------------------!
      !----- 1. Initial value, the most extreme one. --------------------------------------!
      rk4aux(ibuff)%rk4min_can_enthalpy =   huge(1.d0)
      rk4aux(ibuff)%rk4max_can_enthalpy = - huge(1.d0)
      !----- 2. Minimum temperature. ------------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(rk4min_can_temp,can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy =                                                  &
                    min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy =                                                  &
                    max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !----- 3. Maximum temperature. ------------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(rk4max_can_temp,can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy =                                                  &
                    min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy =                                                  &
                    max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !----- 4. Minimum specific humidity. ------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(can_temp,rk4min_can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy =                                                  &
                    min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy =                                                  &
                    max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !----- 5. Maximum specific humidity. ------------------------------------------------!
      can_enthalpy_try    = tq2enthalpy8(can_temp,rk4max_can_shv,.true.)
      rk4aux(ibuff)%rk4min_can_enthalpy =                                                  &
                    min(rk4aux(ibuff)%rk4min_can_enthalpy,can_enthalpy_try)
      rk4aux(ibuff)%rk4max_can_enthalpy =                                                  &
                    max(rk4aux(ibuff)%rk4max_can_enthalpy,can_enthalpy_try)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Print boundaries for scalar quantities.                                        !
      !------------------------------------------------------------------------------------!
      if (print_thbnd) then
         hour   = floor(nint(current_time%time) / hr_sec)
         minute = floor((nint(current_time%time) - hour * hr_sec) / min_sec)
         second = floor( mod( nint(current_time%time) - hour * hr_sec - minute * min_sec   &
                            , min_sec                                                   ) )

         open (unit=39,file=trim(thbnds_fout),status='old',action='write'                  &
                      ,position='append')
         write (unit=39,fmt='(6(i12,1x),12(es12.5,1x))')                                   &
                                 current_time%year, current_time%month,  current_time%date &
                              ,               hour,             minute,             second &
                              ,    rk4min_can_temp,    rk4max_can_temp,     rk4min_can_shv &
                              ,     rk4max_can_shv,     rk4min_can_co2,     rk4max_can_co2 &
                              , rk4aux(ibuff)%rk4min_can_theta                             &
                              , rk4aux(ibuff)%rk4max_can_theta                             &
                              , rk4aux(ibuff)%rk4min_can_prss                              &
                              , rk4aux(ibuff)%rk4max_can_prss                              &
                              , rk4aux(ibuff)%rk4min_can_enthalpy                          &
                              , rk4aux(ibuff)%rk4max_can_enthalpy
         close (unit=39,status='keep')
      end if
      !------------------------------------------------------------------------------------!


      !------ Boundaries for plant internal water content. --------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Plant hydraulics is disabled.  Accept everything. -------------------------!
         do ico=1,cpatch%ncohorts
            rk4aux(ibuff)%rk4min_leaf_water_im2(ico) = - huge_num8
            rk4aux(ibuff)%rk4max_leaf_water_im2(ico) = + huge_num8
            rk4aux(ibuff)%rk4min_wood_water_im2(ico) = - huge_num8
            rk4aux(ibuff)%rk4max_wood_water_im2(ico) = + huge_num8
            rk4aux(ibuff)%rk4min_veg_water_im2 (ico) = - huge_num8
            rk4aux(ibuff)%rk4max_veg_water_im2 (ico) = + huge_num8
         end do
         !---------------------------------------------------------------------------------!
      case default
         do ico=1,cpatch%ncohorts
            ipft    = cpatch%pft(ico)



            !------------------------------------------------------------------------------!
            !    By default, when plant hydraulics is active, wood is resolvable whenever  !
            ! leaf is resolvable.                                                          !
            !------------------------------------------------------------------------------!
            if (cpatch%leaf_resolvable(ico)) then

               !----- Find the lower limits. ----------------------------------------------!
               if (cpatch%is_small(ico)) then
                  !----- Small plant.  Use the small-cohort limits. -----------------------!
                  call psi2twe(small_psi_min(ipft),small_psi_min(ipft),ipft                &
                              ,cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bsapwooda(ico)  &
                              ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico),cpatch%bdeadb(ico) &
                              ,cpatch%broot(ico),cpatch%dbh(ico),leaf_water_im2_4          &
                              ,wood_water_im2_4)
                  !------------------------------------------------------------------------!
               else
                  !----- Large plant.  Use the leaf and wood limits. ----------------------!
                  call psi2twe(leaf_psi_min(ipft),wood_psi_min(ipft),ipft                  &
                              ,cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bsapwooda(ico)  &
                              ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico),cpatch%bdeadb(ico) &
                              ,cpatch%broot(ico),cpatch%dbh(ico),leaf_water_im2_4          &
                              ,wood_water_im2_4)
                  !------------------------------------------------------------------------!
               end if
               rk4aux(ibuff)%rk4min_leaf_water_im2(ico) = dble(leaf_water_im2_4)
               rk4aux(ibuff)%rk4min_wood_water_im2(ico) = dble(wood_water_im2_4)
               rk4aux(ibuff)%rk4min_veg_water_im2 (ico) =                                  &
                                                  rk4aux(ibuff)%rk4min_leaf_water_im2(ico) &
                                                + rk4aux(ibuff)%rk4min_wood_water_im2(ico)
               !---------------------------------------------------------------------------!



               !----- Find the upper limits. ----------------------------------------------!
               call psi2twe(0.0,0.0,ipft,cpatch%nplant(ico),cpatch%bleaf(ico)              &
                           ,cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico),cpatch%bdeada(ico) &
                           ,cpatch%bdeadb(ico),cpatch%broot(ico),cpatch%dbh(ico)           &
                           ,leaf_water_im2_4,wood_water_im2_4)
               rk4aux(ibuff)%rk4max_leaf_water_im2(ico) = dble(leaf_water_im2_4)
               rk4aux(ibuff)%rk4max_wood_water_im2(ico) = dble(wood_water_im2_4)
               rk4aux(ibuff)%rk4max_veg_water_im2 (ico) =                                  &
                                                  rk4aux(ibuff)%rk4max_leaf_water_im2(ico) &
                                                + rk4aux(ibuff)%rk4max_wood_water_im2(ico)
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !      Accept everything.                                                   !
               !---------------------------------------------------------------------------!
               rk4aux(ibuff)%rk4min_leaf_water_im2(ico) = - huge_num8
               rk4aux(ibuff)%rk4max_leaf_water_im2(ico) = + huge_num8
               rk4aux(ibuff)%rk4min_wood_water_im2(ico) = - huge_num8
               rk4aux(ibuff)%rk4max_wood_water_im2(ico) = + huge_num8
               rk4aux(ibuff)%rk4min_veg_water_im2 (ico) = - huge_num8
               rk4aux(ibuff)%rk4max_veg_water_im2 (ico) = + huge_num8
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Limits for the soil moisture. ------------------------------------------------!
      do k = rk4site%lsl, nzg
         nsoil = rk4site%ntext_soil(k)
         rk4aux(ibuff)%rk4min_soil_water(k) = soil8(nsoil)%soilcp * (1.d0 - rk4eps)
         rk4aux(ibuff)%rk4max_soil_water(k) = soil8(nsoil)%soilbp * (1.d0 + rk4eps)
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine find_derived_thbounds
   !=======================================================================================!
   !=======================================================================================!









   !=======================================================================================!
   !=======================================================================================!
   !    This function simply checks whether the relative error is large or not.            !
   !---------------------------------------------------------------------------------------!
   logical function large_error(err,scal)
      use rk4_coms , only : rk4eps ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: err  ! Absolute error
      real(kind=8), intent(in) :: scal ! Characteristic scale
      !------------------------------------------------------------------------------------!
      if(scal > 0.d0) then
         large_error = abs(err/scal)/rk4eps > 1.d0
      else
         large_error = .false.
      end if
      return
   end function large_error
   !=======================================================================================!
   !=======================================================================================!
end module rk4_misc
!==========================================================================================!
!==========================================================================================!
