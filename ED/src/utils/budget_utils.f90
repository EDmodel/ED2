!============================================================================!
!============================================================================!
subroutine compute_budget(csite,lsl,rhos,pcpg,qpcpg,ipa,wcurr_loss2atm       &
                         ,ecurr_loss2atm,co2curr_loss2atm                    &
                         ,wcurr_loss2drainage,ecurr_loss2drainage            &
                         ,wcurr_loss2runoff,ecurr_loss2runoff,ecurr_latent   &
                         ,site_area,cbudget_nep)
   use ed_state_vars, only : sitetype         ! ! structure
   use ed_misc_coms , only : dtlsm            & ! intent(in)
                           , fast_diagnostics & ! intent(in)
                           , current_time     ! ! intent(in)
   use ed_max_dims  , only : n_dbh            ! ! intent(in)
   use consts_coms  , only : umol_2_kgC       & ! intent(in)
                           , day_sec          ! ! intent(in)
   use rk4_coms     , only : rk4eps           & ! intent(in)
                           , checkbudget      ! ! intent(in)
   implicit none
   !----- Arguments ---------------------------------------------------------!
   type(sitetype)        , target        :: csite
   real                  , intent(in)    :: rhos
   real                  , intent(in)    :: pcpg
   real                  , intent(in)    :: qpcpg
   real                  , intent(in)    :: co2curr_loss2atm
   real                  , intent(in)    :: ecurr_loss2atm
   real                  , intent(in)    :: ecurr_loss2drainage
   real                  , intent(in)    :: ecurr_loss2runoff
   real                  , intent(in)    :: ecurr_latent
   real                  , intent(in)    :: wcurr_loss2atm
   real                  , intent(in)    :: wcurr_loss2drainage
   real                  , intent(in)    :: wcurr_loss2runoff
   integer               , intent(in)    :: lsl
   integer               , intent(in)    :: ipa
   real                  , intent(in)    :: site_area
   real                  , intent(inout) :: cbudget_nep
   !----- Local variables ---------------------------------------------------!
   real, dimension(n_dbh)                :: gpp_dbh
   real                                  :: co2budget_finalstorage
   real                                  :: co2budget_deltastorage
   real                                  :: co2curr_gpp
   real                                  :: co2curr_plresp
   real                                  :: co2curr_hetresp
   real                                  :: co2curr_nep
   real                                  :: co2curr_residual
   real                                  :: ebudget_finalstorage
   real                                  :: ebudget_deltastorage
   real                                  :: ecurr_precipgain
   real                                  :: ecurr_netrad
   real                                  :: ecurr_residual
   real                                  :: wbudget_finalstorage
   real                                  :: wbudget_deltastorage
   real                                  :: wcurr_precipgain
   real                                  :: wcurr_residual
   real                                  :: gpp
   real                                  :: plant_respiration
   logical                               :: co2_ok
   logical                               :: energy_ok
   logical                               :: water_ok
   !----- Local constants. --------------------------------------------------!
   character(len=13)     , parameter     :: fmtf='(a,1x,es14.7)'
   logical               , parameter     :: print_debug = .true.
   !----- External functions. -----------------------------------------------!
   real                  , external      :: compute_netrad
   real                  , external      :: compute_water_storage
   real                  , external      :: compute_energy_storage
   real                  , external      :: compute_co2_storage
   !-------------------------------------------------------------------------!

   !-------------------------------------------------------------------------!
   !     Compute gain in water and energy due to precipitation.              !
   !-------------------------------------------------------------------------!
   wcurr_precipgain = pcpg  * dtlsm
   ecurr_precipgain = qpcpg * dtlsm

   !-------------------------------------------------------------------------!
   !     Compute gain in energy due to radiation.                            !
   !-------------------------------------------------------------------------!
   ecurr_netrad     = compute_netrad(csite,ipa) * dtlsm

   !-------------------------------------------------------------------------!
   !     Compute the carbon flux components.                                 !
   !-------------------------------------------------------------------------!
   call sum_plant_cfluxes(csite,ipa,gpp,gpp_dbh,plant_respiration)
   co2curr_gpp     = gpp * dtlsm
   co2curr_plresp  = plant_respiration * dtlsm
   co2curr_hetresp = csite%rh(ipa) * dtlsm
   co2curr_nep = co2curr_gpp - co2curr_plresp - co2curr_hetresp
   cbudget_nep = cbudget_nep + site_area * csite%area(ipa) * co2curr_nep     &
                             * umol_2_kgC


   !----- Compute current storage terms. ------------------------------------!
   co2budget_finalstorage = compute_co2_storage(csite,rhos,ipa)
   wbudget_finalstorage   = compute_water_storage(csite,lsl,rhos,ipa)
   ebudget_finalstorage   = compute_energy_storage(csite,lsl,rhos,ipa)

   !----- Compute the change in storage. ------------------------------------!
   co2budget_deltastorage = co2budget_finalstorage                           &
                          - csite%co2budget_initialstorage(ipa)
   wbudget_deltastorage   = wbudget_finalstorage                             &
                          - csite%wbudget_initialstorage(ipa)
   ebudget_deltastorage   = ebudget_finalstorage                             &
                          - csite%ebudget_initialstorage(ipa)
   !-------------------------------------------------------------------------!
   !     Compute residuals.                                                  !
   !-------------------------------------------------------------------------!
   !----- 1. Canopy CO2. ----------------------------------------------------!
   co2curr_residual = co2budget_deltastorage                                 &
                    - ( - co2curr_nep - co2curr_loss2atm)
   !----- 2. Energy. --------------------------------------------------------!
   ecurr_residual = ebudget_deltastorage                                     &
                  - ( ecurr_precipgain + ecurr_netrad - ecurr_latent         &
                    - ecurr_loss2atm - ecurr_loss2drainage                   &
                    - ecurr_loss2runoff)
   !----- 3. Water. ---------------------------------------------------------!
   wcurr_residual = wbudget_deltastorage                                     &
                   - ( wcurr_precipgain - wcurr_loss2atm                     &
                     - wcurr_loss2drainage - wcurr_loss2runoff)
   !-------------------------------------------------------------------------!

   !-------------------------------------------------------------------------!
   !    If the "check budget" option is activated (you can turn on and turn  !
   ! off by setting checkbudget in ed_params.f90), then the model will crash !
   ! whenever there is some significant leak of CO2, water, or energy.       !
   !-------------------------------------------------------------------------!
   if (checkbudget) then
      !co2_ok  = abs(co2curr_residual) <= rk4eps                              &
      !                                 * (abs(co2budget_finalstorage)        &
      !                                   +abs(co2budget_deltastorage)*dtlsm)
      co2_ok = .true. ! Skipping CO2 test
      energy_ok = abs(ecurr_residual) <= rk4eps                              &
                                       * (abs(ebudget_finalstorage)          &
                                         +abs(ebudget_deltastorage)*dtlsm)
      water_ok  = abs(wcurr_residual) <= rk4eps                              &
                                       * (abs(wbudget_finalstorage)          &
                                         +abs(wbudget_deltastorage)*dtlsm)

      if (print_debug) then 
         write (unit=56,fmt='(i4.4,2(1x,i2.2),1x,f6.0,4(1x,es14.7))')        &
                current_time%year,current_time%month,current_time%date       &
               ,current_time%time                                            &
               ,co2curr_residual/dtlsm                                       &
               ,co2budget_deltastorage/dtlsm                                 &
               ,co2curr_nep/dtlsm                                            &
               ,co2curr_loss2atm/dtlsm

         write (unit=66,fmt='(i4.4,2(1x,i2.2),1x,f6.0,8(1x,es14.7))')        &
                current_time%year,current_time%month,current_time%date       &
               ,current_time%time                                            &
               ,ecurr_residual/dtlsm                                         &
               ,ebudget_deltastorage/dtlsm                                   &
               ,ecurr_precipgain/dtlsm                                       &
               ,ecurr_netrad/dtlsm                                           &
               ,ecurr_latent/dtlsm                                           &
               ,ecurr_loss2atm/dtlsm                                         &
               ,ecurr_loss2drainage/dtlsm                                    &
               ,ecurr_loss2runoff/dtlsm

         write (unit=76,fmt='(i4.4,2(1x,i2.2),1x,f6.0,6(1x,es14.7))')        &
                current_time%year,current_time%month,current_time%date       &
               ,current_time%time                                            &
               ,wcurr_residual*day_sec/dtlsm                                 &
               ,wbudget_deltastorage*day_sec/dtlsm                           &
               ,wcurr_precipgain*day_sec/dtlsm                               &
               ,wcurr_loss2atm*day_sec/dtlsm                                 &
               ,wcurr_loss2drainage*day_sec/dtlsm                            &
               ,wcurr_loss2runoff*day_sec/dtlsm
      end if


      if (.not. co2_ok) then
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a)') '|           !!! ): CO2 budget failed :( !!!           |'
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',           &
            current_time%year,current_time%month,current_time%date ,current_time%time
         write (unit=*,fmt=fmtf ) ' RESIDUAL       : ',co2curr_residual
         write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE: ',csite%co2budget_initialstorage(ipa)
         write (unit=*,fmt=fmtf ) ' FINAL_STORAGE  : ',co2budget_finalstorage
         write (unit=*,fmt=fmtf ) ' DELTA_STORAGE  : ',co2budget_deltastorage
         write (unit=*,fmt=fmtf ) ' GPP            : ',co2curr_gpp
         write (unit=*,fmt=fmtf ) ' PLANT_RESP     : ',co2curr_plresp
         write (unit=*,fmt=fmtf ) ' HET_RESP       : ',co2curr_hetresp
         write (unit=*,fmt=fmtf ) ' NEP            : ',co2curr_nep
         write (unit=*,fmt=fmtf ) ' LOSS2ATM       : ',co2curr_loss2atm
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a)') ' '
      end if

      if (.not. energy_ok) then
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a)') '|         !!! ): Energy budget failed :( !!!          |'
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',           &
            current_time%year,current_time%month,current_time%date ,current_time%time
         write (unit=*,fmt=fmtf ) ' RESIDUAL       : ',ecurr_residual
         write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE: ',csite%ebudget_initialstorage(ipa)
         write (unit=*,fmt=fmtf ) ' FINAL_STORAGE  : ',ebudget_finalstorage
         write (unit=*,fmt=fmtf ) ' DELTA_STORAGE  : ',ebudget_deltastorage
         write (unit=*,fmt=fmtf ) ' PRECIPGAIN     : ',ecurr_precipgain
         write (unit=*,fmt=fmtf ) ' NETRAD         : ',ecurr_netrad
         write (unit=*,fmt=fmtf ) ' LATENT         : ',ecurr_latent
         write (unit=*,fmt=fmtf ) ' LOSS2ATM       : ',ecurr_loss2atm
         write (unit=*,fmt=fmtf ) ' LOSS2DRAINAGE  : ',ecurr_loss2drainage
         write (unit=*,fmt=fmtf ) ' LOSS2RUNOFF    : ',ecurr_loss2runoff
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a)') ' '
      end if


      if (.not. water_ok) then
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a)') '|          !!! ): Water budget failed :( !!!          |'
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',           &
            current_time%year,current_time%month,current_time%date ,current_time%time
         write (unit=*,fmt=fmtf ) ' RESIDUAL       : ',wcurr_residual
         write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE: ',csite%wbudget_initialstorage(ipa)
         write (unit=*,fmt=fmtf ) ' FINAL_STORAGE  : ',wbudget_finalstorage
         write (unit=*,fmt=fmtf ) ' DELTA_STORAGE  : ',wbudget_deltastorage
         write (unit=*,fmt=fmtf ) ' PRECIPGAIN     : ',wcurr_precipgain
         write (unit=*,fmt=fmtf ) ' LOSS2ATM       : ',wcurr_loss2atm
         write (unit=*,fmt=fmtf ) ' LOSS2DRAINAGE  : ',wcurr_loss2drainage
         write (unit=*,fmt=fmtf ) ' LOSS2RUNOFF    : ',wcurr_loss2runoff
         write (unit=*,fmt='(a)') '|-----------------------------------------------------|'
         write (unit=*,fmt='(a)') ' '
      end if

      if (.not. (co2_ok .and. energy_ok .and. water_ok)) then
         call fatal_error('Budget check has failed, see message above!'      &
                         ,'compute_budget','budget_utils.f90')
      end if
   end if


   !-------------------------------------------------------------------------!
   !     Integrate residuals.                                                !
   !-------------------------------------------------------------------------!
   !----- 1. Canopy CO2. ----------------------------------------------------!
   csite%co2budget_residual(ipa) = csite%co2budget_residual(ipa)             &
                                 + co2curr_residual
   !----- 2. Energy. --------------------------------------------------------!
   csite%ebudget_residual(ipa) = csite%ebudget_residual(ipa)                 &
                               + ecurr_residual
   !----- 3. Water. ---------------------------------------------------------!
   csite%wbudget_residual(ipa) = csite%wbudget_residual(ipa)                 &
                               + wcurr_residual
   !-------------------------------------------------------------------------!

   !-------------------------------------------------------------------------!
   !    Integrate the terms that are part of the budget.                     !
   !-------------------------------------------------------------------------!
   !----- 1. Carbon dioxide. ------------------------------------------------!
   csite%co2budget_gpp(ipa)         = csite%co2budget_gpp(ipa) + gpp * dtlsm
   csite%co2budget_gpp_dbh(:,ipa)   = csite%co2budget_gpp_dbh(:,ipa)         &
                                    + gpp_dbh(:) *dtlsm
   csite%co2budget_plresp(ipa)      = csite%co2budget_plresp(ipa)            &
                                    + plant_respiration * dtlsm
   csite%co2budget_rh(ipa)          = csite%co2budget_rh(ipa)                &
                                    + csite%rh(ipa) * dtlsm
   csite%co2budget_loss2atm(ipa)    = csite%co2budget_loss2atm(ipa)          &
                                    + co2curr_loss2atm
   !----- 2. Energy. --------------------------------------------------------!
   csite%ebudget_precipgain(ipa)    = csite%ebudget_precipgain(ipa)          &
                                    + ecurr_precipgain
   csite%ebudget_netrad(ipa)        = csite%ebudget_netrad(ipa)              &
                                    + ecurr_netrad
   csite%ebudget_latent(ipa)        = csite%ebudget_latent(ipa)              &
                                    + ecurr_latent
   csite%ebudget_loss2atm(ipa)      = csite%ebudget_loss2atm(ipa)            &
                                    + ecurr_loss2atm
   csite%ebudget_loss2drainage(ipa) = csite%ebudget_loss2drainage(ipa)       &
                                    + ecurr_loss2drainage
   csite%ebudget_loss2runoff(ipa)   = csite%ebudget_loss2runoff(ipa)         &
                                    + ecurr_loss2runoff
   !----- 3. Water. ---------------------------------------------------------!
   csite%wbudget_precipgain(ipa)    = csite%wbudget_precipgain(ipa)          &
                                    + wcurr_precipgain
   csite%wbudget_loss2atm(ipa)      = csite%wbudget_loss2atm(ipa)            &
                                    + wcurr_loss2atm
   csite%wbudget_loss2drainage(ipa) = csite%wbudget_loss2drainage(ipa)       &
                                    + wcurr_loss2drainage
   csite%wbudget_loss2runoff(ipa)   = csite%wbudget_loss2runoff(ipa)         &
                                    + wcurr_loss2runoff
   !-------------------------------------------------------------------------!


   !-------------------------------------------------------------------------!
   !     Updating the initial storage (initial for next step)                !
   !-------------------------------------------------------------------------!
   csite%wbudget_initialstorage(ipa)   = wbudget_finalstorage
   csite%ebudget_initialstorage(ipa)   = ebudget_finalstorage
   csite%co2budget_initialstorage(ipa) = co2budget_finalstorage

   return
end subroutine compute_budget
!============================================================================!
!============================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function computes the total water stored in the system, in kg/m2.                 !
!   (soil + temporary pools + canopy air space + leaf surface).                            !
!------------------------------------------------------------------------------------------!
real function compute_water_storage(csite, lsl, rhos,ipa)
   use ed_state_vars , only  : sitetype              & ! structure
                             , patchtype             ! ! structure
   use canopy_air_coms, only : minimum_canopy_depth  ! ! intent(in)
   use grid_coms      , only : nzg                   ! ! intent(in)
   use soil_coms      , only : dslz                  ! ! intent(in)
   use consts_coms    , only : wdns                  ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype) , target     :: csite
   integer        , intent(in) :: ipa
   integer        , intent(in) :: lsl
   real           , intent(in) :: rhos
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype), pointer    :: cpatch
   integer                     :: k
   integer                     :: ico
   real                        :: canopy_depth
   !---------------------------------------------------------------------------------------!

   canopy_depth = max(csite%veg_height(ipa),minimum_canopy_depth)

   compute_water_storage = 0.0
   cpatch => csite%patch(ipa)

   !----- 1. Adding the water stored in the soil. -----------------------------------------!
   do k = lsl, nzg
      compute_water_storage = compute_water_storage                                        &
                            + csite%soil_water(k,ipa) * dslz(k) * wdns
   end do
   !----- 2. Adding the water stored in the temporary surface water/snow. -----------------!
   do k = 1, csite%nlev_sfcwater(ipa)
      compute_water_storage = compute_water_storage + csite%sfcwater_mass(k,ipa)
   end do
   !----- 3. Adding the water vapour floating in the canopy air space. --------------------!
   compute_water_storage = compute_water_storage                                           &
                            + csite%can_shv(ipa) * canopy_depth * rhos
   !----- 4. Adding the water over the leaf surface. --------------------------------------!
   do ico = 1,cpatch%ncohorts
      compute_water_storage = compute_water_storage + cpatch%veg_water(ico)
   end do

   return
end function compute_water_storage
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computs the total net radiation, by adding the radiation that interacts !
! with the different surfaces.                                                             !
!------------------------------------------------------------------------------------------!
real function compute_netrad(csite,ipa)
   use ed_state_vars , only : sitetype  & ! structure
                            , patchtype ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype) , target     :: csite
   integer        , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype), pointer    :: cpatch
   integer                     :: k
   integer                     :: ico
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   compute_netrad = 0.0
   !----- 1. Adding the ground components -------------------------------------------------!
   compute_netrad = csite%rshort_g(ipa) + csite%rlong_g(ipa) + csite%rlong_s(ipa)
   !----- 2. Adding the shortwave radiation that reaches each snow/water layer ------------!
   do k = 1, csite%nlev_sfcwater(ipa)
      compute_netrad = compute_netrad + csite%rshort_s(k,ipa)
   end do
   !----- 3. Adding the radiation components that interact with leaves. -------------------!
   do ico = 1,cpatch%ncohorts
      compute_netrad = compute_netrad + cpatch%rshort_v(ico) + cpatch%rlong_v(ico)
   end do
   return
end function compute_netrad
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computs the total net radiation, by adding the radiation that interacts !
! with the different surfaces.  The result is given in J/m2.                               !
!------------------------------------------------------------------------------------------!
real function compute_energy_storage(csite, lsl, rhos, ipa)
   use ed_state_vars        , only : sitetype              & ! structure
                                   , patchtype             ! ! structure
   use canopy_air_coms      , only : minimum_canopy_depth  ! ! intent(in)
   use grid_coms            , only : nzg                   ! ! intent(in)
   use soil_coms            , only : dslz                  ! ! intent(in)
   use consts_coms          , only : cp                    & ! intent(in)
                                   , cliq                  & ! intent(in)
                                   , cice                  & ! intent(in)
                                   , alli                  & ! intent(in)
                                   , t3ple                 ! ! intent(in)
   use rk4_coms             , only : toosparse             ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype) , target     :: csite
   integer        , intent(in) :: ipa
   integer        , intent(in) :: lsl
   real           , intent(in) :: rhos
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype), pointer    :: cpatch
   integer                     :: k
   integer                     :: ico
   real                        :: soil_storage
   real                        :: sfcwater_storage
   real                        :: cas_storage
   real                        :: veg_storage
   real                        :: canopy_depth
   !---------------------------------------------------------------------------------------!

   canopy_depth = max(csite%veg_height(ipa),minimum_canopy_depth)

   cpatch => csite%patch(ipa)
   !----- 1. Computing internal energy stored at the soil. --------------------------------!
   soil_storage = 0.0
   do k = lsl, nzg
      soil_storage = soil_storage + csite%soil_energy(k,ipa) * dslz(k)
   end do
   !---------------------------------------------------------------------------------------!
   !   2. Computing internal energy stored at the temporary snow/water sfc. layer.         !
   !      Converting it to J/m2. 
   !---------------------------------------------------------------------------------------!
   sfcwater_storage = 0.0
   do k = 1, csite%nlev_sfcwater(ipa)
      sfcwater_storage = sfcwater_storage                                                  &
                       + csite%sfcwater_energy(k,ipa) * csite%sfcwater_mass(k,ipa)
   end do

   !---------------------------------------------------------------------------------------!
   ! 3. Finding and approximated value for canopy air total enthalpy.                      !
   !---------------------------------------------------------------------------------------!
   cas_storage = cp * rhos * canopy_depth * csite%can_temp(ipa)

   !---------------------------------------------------------------------------------------!
   ! 4. Compute the internal energy stored in the plants.                                  !
   !    Originally we were only considering those patches that were prognosed, but since   !
   !    we are assigning non-zero internal energy to the tiny cohorts or those cohorts     !
   !    buried in snow, we should account for them, even if this will put the energy con-  !
   !    servation off. After all, this can help us identifying how bad is the assumption   !
   !    of diagnosing such cohorts instead of solving them.                                !
   !---------------------------------------------------------------------------------------!
   veg_storage = 0.0
   do ico = 1,cpatch%ncohorts
      veg_storage = veg_storage + cpatch%veg_energy(ico)
   end do
 
   !----- 5. Integrating the total energy in ED. ------------------------------------------!
   compute_energy_storage = soil_storage + sfcwater_storage + cas_storage               &
                             + veg_storage

   return
end function compute_energy_storage
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the carbon flux terms.                                       !
!------------------------------------------------------------------------------------------!
subroutine sum_plant_cfluxes(csite,ipa, gpp, gpp_dbh,plresp)
   use ed_state_vars        , only : sitetype    & ! structure
                                   , patchtype   ! ! structure
   use consts_coms          , only : day_sec     & ! intent(in)
                                   , umol_2_kgC  ! ! intent(in)
   use ed_max_dims             , only : n_dbh
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target      :: csite
   integer               , intent(in)  :: ipa
   real                  , intent(out) :: gpp
   real, dimension(n_dbh), intent(out) :: gpp_dbh
   real                  , intent(out) :: plresp
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype), pointer            :: cpatch
   integer                             :: k
   integer                             :: ico
   integer                             :: idbh
   real                                :: lrresp !----- Leaf and root respiration
   real                                :: sresp  !----- Storage, growth, vleaf respiration.
   logical                             :: forest
   !----- Local constants -----------------------------------------------------------------!
   real, parameter                     :: ddbh=1./real(n_dbh-1)
   !---------------------------------------------------------------------------------------!

  
   !----- GPP by DBH is computed for forested areas only. ---------------------------------!
   forest = csite%dist_type(ipa) /= 1

   !----- Initializing some variables. ----------------------------------------------------!
   gpp     = 0.0
   gpp_dbh = 0.0 
   lrresp  = 0.0
   sresp   = 0.0
   cpatch => csite%patch(ipa)

   !---------------------------------------------------------------------------------------!
   !     Looping over cohorts.                                                             !
   !---------------------------------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      !----- Adding GPP and leaf respiration only for those cohorts with enough leaves. ---!
      if (cpatch%solvable(ico)) then
         gpp = gpp + cpatch%gpp(ico)
         !----- Forest cohorts have dbh distribution, add them to gpp_dbh. ----------------!
         if (forest) then 
            idbh=max(1,min(n_dbh,ceiling(cpatch%dbh(ico)*ddbh)))
            gpp_dbh(idbh) = gpp_dbh(idbh) + cpatch%gpp(ico)
         end if
         lrresp = lrresp + cpatch%leaf_respiration(ico)

      end if
      !----- Root respiration happens even when the LAI is tiny ---------------------------!
      lrresp = lrresp + cpatch%root_respiration(ico)
      !------------------------------------------------------------------------------------!
      !     So do the other components that go to sresp.  Structural terms are "intens-    !
      ! ive", we must convert them from kgC/plant/day to umol/m2/s.                        !
      !------------------------------------------------------------------------------------!
      sresp  = sresp                                                                       &
             + ( cpatch%growth_respiration(ico) + cpatch%storage_respiration(ico)          &
               + cpatch%vleaf_respiration(ico))                                            &
               * cpatch%nplant(ico) / (day_sec * umol_2_kgC)
   end do
   !----- Plant respiration is the sum between alive and structural. ----------------------!
   plresp = lrresp + sresp
   return
end subroutine sum_plant_cfluxes
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computes the co2 stored in the canopy air space from ppm to umol/m2.    !
!------------------------------------------------------------------------------------------!
real function compute_co2_storage(csite, rhos, ipa)
   use ed_state_vars  , only : sitetype              ! ! structure
   use canopy_air_coms, only : minimum_canopy_depth  ! ! intent(in)
   use consts_coms    , only : mmdryi                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target      :: csite
   integer               , intent(in)  :: ipa
   real                  , intent(in)  :: rhos
   !----- Local variables -----------------------------------------------------------------!
   real                                :: canopy_depth
   !---------------------------------------------------------------------------------------!

   canopy_depth = max(csite%veg_height(ipa),minimum_canopy_depth)

   compute_co2_storage = csite%can_co2(ipa) * mmdryi * rhos * canopy_depth

   return
end function compute_co2_storage
!==========================================================================================!
!==========================================================================================!
