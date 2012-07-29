!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the soil respiration terms (root and heterotrophic).        !
!------------------------------------------------------------------------------------------!
subroutine soil_respiration(csite,ipa,mzg,ntext_soil)

   use ed_state_vars, only : sitetype                 & ! structure
                           , patchtype                ! ! structure
   use soil_coms    , only : soil                     & ! intent(in)
                           , dslz                     & ! intent(in)
                           , slz                      ! ! intent(in)
   use decomp_coms  , only : k_rh_active              ! ! intent(in)
   use consts_coms  , only : wdns                     ! ! intent(in)
   use therm_lib    , only : uextcm2tl                ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)                , target     :: csite
   integer                       , intent(in) :: ipa
   integer                       , intent(in) :: mzg
   integer       , dimension(mzg), intent(in) :: ntext_soil
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype)               , pointer    :: cpatch
   integer                                    :: ico
   integer                                    :: ipft
   integer                                    :: k
   integer                                    :: kroot
   integer                                    :: nsoil
   real                                       :: Lc
   real                                       :: rel_soil_moist
   real                                       :: sum_soil_energy
   real                                       :: sum_soil_water
   real                                       :: sum_soil_hcap
   real                                       :: sum_soil_slmsts
   real                                       :: sum_soil_soilcp
   real                                       :: avg_soil_temp
   real                                       :: avg_soil_fliq
   !----- External functions. -------------------------------------------------------------!
   real                          , external   :: het_resp_weight
   real                          , external   :: root_resp_norm
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find the current root respiration.  This is done cohort by cohort because the    !
   ! parameters may be different depending on the PFT, and also because each layer has a   !
   ! different temperature.                                                                !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   do ico = 1,cpatch%ncohorts
      ipft  = cpatch%pft(ico)
      kroot = cpatch%krdepth(ico)

      !------------------------------------------------------------------------------------!
      !    Add "intensive" contribution of each layer, assuming that the roots are equally !
      ! spread throughout the entire depth.                                                !
      !------------------------------------------------------------------------------------!
      cpatch%root_respiration(ico) = 0.0
      do k = kroot,mzg
         cpatch%root_respiration(ico) = cpatch%root_respiration(ico)                       &
                                      + root_resp_norm(ipft,csite%soil_tempk(k,ipa))       &
                                      * dslz(k)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Now we make the value in umol/m2/s, by dividing by the total depth and        !
      ! multiplying by the total root biomass.  The minus sign is because slz is negative. !
      !------------------------------------------------------------------------------------!
      cpatch%root_respiration(ico) = - cpatch%root_respiration(ico) * cpatch%broot(ico)    &
                                     * cpatch%nplant(ico) / slz(kroot)
      !------------------------------------------------------------------------------------!


      !----- Add this time step to the mean and daily mean root respiration. --------------!
      cpatch%mean_root_resp(ico)   = cpatch%mean_root_resp(ico)                            &
                                   + cpatch%root_respiration(ico)
      cpatch%today_root_resp(ico)  = cpatch%today_root_resp(ico)                           &
                                   + cpatch%root_respiration(ico)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Integrate the soil extensive properties, plus the minimum and maximum possible    !
   ! soil water content of the active layer.                                               !
   !---------------------------------------------------------------------------------------!
   sum_soil_energy = 0.0
   sum_soil_hcap   = 0.0
   sum_soil_water  = 0.0
   sum_soil_slmsts = 0.0
   sum_soil_soilcp = 0.0
   do k = k_rh_active,mzg
      nsoil = ntext_soil(k)

      !------------------------------------------------------------------------------------!
      !    Convert the units so energy is in J/m2, heat capacity in J/m2/K, and water in   !
      ! kg/m2.                                                                             !
      !------------------------------------------------------------------------------------!
      sum_soil_energy = sum_soil_energy + csite%soil_energy(k,ipa)        * dslz(k)
      sum_soil_hcap   = sum_soil_hcap   + soil(nsoil)%slcpd               * dslz(k)
      sum_soil_water  = sum_soil_water  + csite%soil_water (k,ipa) * wdns * dslz(k)
      sum_soil_slmsts = sum_soil_slmsts + soil(nsoil)%slmsts       * wdns * dslz(k)
      sum_soil_soilcp = sum_soil_soilcp + soil(nsoil)%soilcp       * wdns * dslz(k)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!



   !----- Find the average temperature and the relative soil moisture. --------------------!
   call uextcm2tl(sum_soil_energy,sum_soil_water,sum_soil_hcap,avg_soil_temp,avg_soil_fliq)
   rel_soil_moist = min( 1.0, max(0.0, ( sum_soil_water  - sum_soil_soilcp )               &
                                     / ( sum_soil_slmsts - sum_soil_soilcp ) ) )
   !---------------------------------------------------------------------------------------!



   !----- Compute soil/temperature modulation of heterotrophic respiration. ---------------!
   csite%A_decomp(ipa) = het_resp_weight(avg_soil_temp,rel_soil_moist)
   !---------------------------------------------------------------------------------------!



   !----- Compute nitrogen immobilization factor. -----------------------------------------!
   call resp_f_decomp(csite,ipa, Lc)
   !---------------------------------------------------------------------------------------!



   !----- Compute heterotrophic respiration. ----------------------------------------------!
   call resp_rh(csite,ipa, Lc)
   !---------------------------------------------------------------------------------------!



   !----- Update averaged variables. ------------------------------------------------------!
   csite%today_A_decomp (ipa) = csite%today_A_decomp (ipa) + csite%A_decomp(ipa)
   csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa)                                 &
                              + csite%A_decomp       (ipa) * csite%f_decomp(ipa)
   csite%mean_rh        (ipa) = csite%mean_rh        (ipa) + csite%rh      (ipa)
   csite%mean_cwd_rh    (ipa) = csite%mean_cwd_rh    (ipa) + csite%cwd_rh  (ipa)
   !---------------------------------------------------------------------------------------!

   return
end subroutine soil_respiration
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function determines the normalised root respiration (umol/kg_fine_root/s) at a  !
! given soil layer.                                                                        !
!------------------------------------------------------------------------------------------!
real function root_resp_norm(ipft,soil_temp)
   use pft_coms       , only : root_respiration_factor  & ! intent(in)
                             , rrf_low_temp             & ! intent(in)
                             , rrf_high_temp            & ! intent(in)
                             , rrf_decay_e              & ! intent(in)
                             , rrf_hor                  & ! intent(in)
                             , rrf_q10                  ! ! intent(in)
   use farq_leuning   , only : arrhenius                & ! function
                             , collatz                  ! ! function
   use soil_coms      , only : soil8                    ! ! intent(in)
   use rk4_coms       , only : tiny_offset              ! ! intent(in)
   use physiology_coms, only : iphysiol                 ! ! intent(in)
   use consts_coms    , only : lnexp_min8               & ! intent(in)
                             , lnexp_max8               & ! intent(in)
                             , t008                     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in) :: ipft
   real(kind=4), intent(in) :: soil_temp
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: soil_temp8
   real(kind=8)             :: rrf08
   real(kind=8)             :: rrf_low_temp8
   real(kind=8)             :: rrf_high_temp8
   real(kind=8)             :: rrf_decay_e8
   real(kind=8)             :: rrf_hor8
   real(kind=8)             :: rrf_q108
   real(kind=8)             :: lnexplow
   real(kind=8)             :: lnexphigh
   real(kind=8)             :: tlow_fun
   real(kind=8)             :: thigh_fun
   real(kind=8)             :: rrf8
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)             :: sngloff
   !---------------------------------------------------------------------------------------!

   !----- Copy some variables to double precision temporaries. ----------------------------!
   soil_temp8      = dble(soil_temp                    )
   rrf08           = dble(root_respiration_factor(ipft))
   rrf_low_temp8   = dble(rrf_low_temp           (ipft)) + t008
   rrf_high_temp8  = dble(rrf_high_temp          (ipft)) + t008
   rrf_decay_e8    = dble(rrf_decay_e            (ipft))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Compute the functions that will control the Rrf function for low and high temper-  !
   ! ature.  In order to avoid floating point exceptions, we check whether the temperature !
   ! will make the exponential too large or too small.                                     !
   !---------------------------------------------------------------------------------------!
   !----- Low temperature. ----------------------------------------------------------------!
   lnexplow  = rrf_decay_e8 * (rrf_low_temp8  - soil_temp8)
   lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
   tlow_fun  = 1.d0 +  exp(lnexplow)
   !----- High temperature. ---------------------------------------------------------------!
   lnexphigh = rrf_decay_e8 * (soil_temp8 - rrf_high_temp8)
   lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
   thigh_fun = 1.d0 + exp(lnexphigh)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Decide which functional form to use based on the physiology.  This is just to make !
   ! it look similar to the leaf respiration respiration.                                  !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,1)
      rrf_hor8 = dble(rrf_hor(ipft))
      rrf8     = arrhenius(soil_temp8,rrf08,rrf_hor8) / (tlow_fun * thigh_fun)
   case (2,3)
      rrf_q108 = dble(rrf_q10(ipft))
      rrf8     = collatz(soil_temp8,rrf08,rrf_q108)   / (tlow_fun * thigh_fun)
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert result to single precision.                                               !
   !---------------------------------------------------------------------------------------!
   root_resp_norm = sngloff(rrf8,tiny_offset)
   !---------------------------------------------------------------------------------------!


   return
end function root_resp_norm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function computes the heterotrophic respiration limitation factor, which        !
! includes limitations due to temperature and soil moisture.                               !
!------------------------------------------------------------------------------------------!
real function het_resp_weight(soil_tempk,rel_soil_moist)

   use decomp_coms, only : resp_temperature_increase  & ! intent(in)
                         , resp_opt_water             & ! intent(in)
                         , resp_water_below_opt       & ! intent(in)
                         , resp_water_above_opt       & ! intent(in)
                         , decomp_scheme              & ! intent(in)
                         , rh_lloyd_1                 & ! intent(in)
                         , rh_lloyd_2                 & ! intent(in)
                         , rh_lloyd_3                 & ! intent(in)
                         , rh_decay_low               & ! intent(in)
                         , rh_decay_high              & ! intent(in)
                         , rh_low_temp                & ! intent(in)
                         , rh_high_temp               & ! intent(in)
                         , rh_decay_dry               & ! intent(in)
                         , rh_decay_wet               & ! intent(in)
                         , rh_dry_smoist              & ! intent(in)
                         , rh_wet_smoist              ! ! intent(in)
   use consts_coms, only : lnexp_min                  & ! intent(in)
                         , lnexp_max                  ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in) :: soil_tempk
   real(kind=4), intent(in) :: rel_soil_moist
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4)             :: temperature_limitation
   real(kind=4)             :: water_limitation
   real(kind=4)             :: lnexplloyd
   real(kind=4)             :: lnexplow
   real(kind=4)             :: lnexphigh
   real(kind=4)             :: tlow_fun
   real(kind=4)             :: thigh_fun
   real(kind=4)             :: lnexpdry
   real(kind=4)             :: lnexpwet
   real(kind=4)             :: smdry_fun
   real(kind=4)             :: smwet_fun
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the temperature dependence.                                                  !
   !---------------------------------------------------------------------------------------!
   select case(decomp_scheme)
   case (0)
      !----- Use original ED-2.1 exponential temperature dependence. ----------------------!
      temperature_limitation = min( 1.0                                                    &
                                  , exp( resp_temperature_increase * (soil_tempk-318.15)))
      !------------------------------------------------------------------------------------!
   case (1) 
      !----- Use Lloyd and Taylor (1994) temperature dependence. --------------------------!
      lnexplloyd             = rh_lloyd_1 * ( rh_lloyd_2 - 1. / (soil_tempk - rh_lloyd_3))
      lnexplloyd             = max(lnexp_min,min(lnexp_max,lnexplloyd))
      temperature_limitation = min( 1.0, resp_temperature_increase * exp(lnexplloyd) )
      !------------------------------------------------------------------------------------!
   case (2)
      !------------------------------------------------------------------------------------!
      !      Similar to the original ED-1.0 formulation, which is based on the CENTURY     !
      ! model.  The change in the functional form is to avoid power of negative numbers,   !
      ! but the coefficients were tuned to give a similar curve.                           !
      !------------------------------------------------------------------------------------!
      !----- Low temperature limitation. --------------------------------------------------!
      lnexplow               = rh_decay_low * (rh_low_temp - soil_tempk)
      lnexplow               = max(lnexp_min,min(lnexp_max,lnexplow))
      tlow_fun               = 1.0 + exp(lnexplow)
      !----- High temperature limitation. -------------------------------------------------!
      lnexphigh              = rh_decay_high * (soil_tempk - rh_high_temp)
      lnexphigh              = max(lnexp_min,min(lnexp_max,lnexphigh))
      thigh_fun              = 1.0 + exp(lnexphigh)
      !----- Temperature limitation is a combination of both. -----------------------------!
      temperature_limitation = 1.0 / (tlow_fun * thigh_fun )
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the limitation due to soil moisture.                                         !
   !---------------------------------------------------------------------------------------!
   select case (decomp_scheme)
   case (0,1)
      !----- ED-2.1 default, also used when decomp_scheme is 1. ---------------------------!
      if (rel_soil_moist <= resp_opt_water)then
         water_limitation = exp( (rel_soil_moist - resp_opt_water) * resp_water_below_opt)
      else
         water_limitation = exp( (resp_opt_water - rel_soil_moist) * resp_water_above_opt)
      end if
      !------------------------------------------------------------------------------------!
   case (2)
      !----- Dry soil limitation. ---------------------------------------------------------!
      lnexpdry         = rh_decay_dry * (rh_dry_smoist - rel_soil_moist)
      lnexpdry         = max(lnexp_min,min(lnexp_max,lnexpdry))
      smdry_fun        = 1.0 + exp(lnexpdry)
      !----- Wet soil limitation. ---------------------------------------------------------!
      lnexpwet         = rh_decay_wet * (rel_soil_moist - rh_wet_smoist)
      lnexpwet         = max(lnexp_min,min(lnexp_max,lnexpwet))
      smwet_fun        = 1.0 + exp(lnexpwet)
      !----- Soil moisture limitation is a combination of both. ---------------------------!
      water_limitation = 1.0 / (smdry_fun * smwet_fun)
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !----- Compute the weight, which is just the combination of both. ----------------------!
   het_resp_weight = temperature_limitation * water_limitation
   !---------------------------------------------------------------------------------------!

   return
end function het_resp_weight
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the Nitrogen immobilization factor.                         !
!------------------------------------------------------------------------------------------!
subroutine resp_f_decomp(csite,ipa,Lc)

   use ed_state_vars, only : sitetype               ! ! structure
   use decomp_coms  , only : r_stsc                 & ! intent(in)
                           , N_immobil_supply_scale & ! intent(in)
                           , K1                     & ! intent(in)
                           , n_decomp_lim           ! ! intent(in)
   use pft_coms     , only : c2n_structural         & ! intent(in)
                           , c2n_slow               ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target      :: csite
   integer       , intent(in)  :: ipa
   real          , intent(out) :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                        :: N_immobilization_demand
   !---------------------------------------------------------------------------------------!

 
   if (csite%structural_soil_C(ipa) > 0.0) then
      if (csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa)) then
         Lc = 0.049787 ! = exp(-3.0)
      else
         Lc = exp(-3.0 * csite%structural_soil_L(ipa)/csite%structural_soil_C(ipa))
      end if
   else
      Lc=0.0
   end if
   
   if (n_decomp_lim == 1) then
      N_immobilization_demand = csite%A_decomp(ipa) * Lc * K1                              &
                              * csite%structural_soil_C(ipa)                               &
                              * ((1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)
      
      csite%f_decomp(ipa)     = N_immobil_supply_scale * csite%mineralized_soil_N(ipa)     &
                              / ( N_immobilization_demand                                  &
                                + N_immobil_supply_scale  * csite%mineralized_soil_N(ipa))
   else
      !----- Option for no plant N limitation. --------------------------------------------!
      csite%f_decomp(ipa)     = 1.0
   end if

   return
end subroutine resp_f_decomp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the heterotrophic respiration.                              !
!------------------------------------------------------------------------------------------!
subroutine resp_rh(csite,ipa,Lc)

   use ed_state_vars, only : sitetype       ! ! structure
   use consts_coms  , only : kgCday_2_umols ! ! intent(in)
   use decomp_coms  , only : K1             & ! intent(in)
                           , K2             & ! intent(in)
                           , K3             & ! intent(in)
                           , r_fsc          & ! intent(in)
                           , r_ssc          & ! intent(in)
                           , r_stsc         & ! intent(in)
                           , cwd_frac       ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype), target       :: csite
   integer       , intent(in)   :: ipa
   real          , intent(in)   :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                         :: fast_C_loss
   real                         :: structural_C_loss
   real                         :: slow_C_loss
   !---------------------------------------------------------------------------------------!



   !----- The following variables have units of [umol_CO2/m2/s]. --------------------------!
   fast_C_loss       = kgCday_2_umols * csite%A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
   structural_C_loss = kgCday_2_umols * csite%A_decomp(ipa) * Lc * K1                      &
                     * csite%structural_soil_C(ipa)* csite%f_decomp(ipa)
   slow_C_loss       = kgCday_2_umols * csite%A_decomp(ipa) * K3 * csite%slow_soil_C(ipa)
   !---------------------------------------------------------------------------------------!

   !----- Find the heterotrophic respiration and the fraction due to CWD. -----------------!
   csite%rh(ipa)     = r_fsc * fast_C_loss + r_stsc * structural_C_loss                    &
                     + r_ssc * slow_C_loss
   csite%cwd_rh(ipa) = cwd_frac * (r_stsc * structural_C_loss + r_ssc * slow_C_loss)
   !---------------------------------------------------------------------------------------!

   return
end subroutine resp_rh
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will update the soil carbon and nitrogen pools.                      !
!------------------------------------------------------------------------------------------!
subroutine update_C_and_N_pools(cgrid)
   
   use ed_state_vars, only : edtype          & ! structure
                           , polygontype     & ! structure
                           , sitetype        ! ! structure
   use decomp_coms  , only : K1              & ! intent(in)
                           , K2              & ! intent(in)
                           , K3              & ! intent(in)
                           , r_stsc          ! ! intent(in)
   use pft_coms     , only : c2n_slow        & ! intent(in)
                           , c2n_structural  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly
   type(sitetype)   , pointer  :: csite
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   real                        :: Lc
   real                        :: fast_C_loss
   real                        :: fast_N_loss
   real                        :: structural_C_loss
   real                        :: structural_L_loss
   real                        :: slow_C_input
   real                        :: slow_C_loss
   !---------------------------------------------------------------------------------------!

   polygonloop: do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         
         csite => cpoly%site(isi)

         patchloop: do ipa = 1,csite%npatches

            if (csite%structural_soil_C(ipa) > 0.0) then
               if (csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa)) then
                  Lc = 0.049787 ! = exp(-3.0)
               else
                  Lc = exp( -3.0 * csite%structural_soil_L(ipa)                            &
                          /  csite%structural_soil_C(ipa))
               end if
            else
               Lc=0.0
            end if
      
            !----- Fast pools. ------------------------------------------------------------!
            fast_C_loss = csite%today_A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
            fast_N_loss = csite%today_A_decomp(ipa) * K2 * csite%fast_soil_N(ipa)

            !----- Structural pools. ------------------------------------------------------!
            structural_C_loss = csite%today_Af_decomp(ipa) * Lc * K1                       &
                              * csite%structural_soil_C(ipa)
            structural_L_loss = csite%today_Af_decomp(ipa) * Lc * K1                       &
                              * csite%structural_soil_L(ipa)

            !----- Slow pools. ------------------------------------------------------------!
            slow_C_input = (1.0 - r_stsc) * structural_C_loss
            slow_C_loss  = csite%today_A_decomp(ipa) * K3 * csite%slow_soil_C(ipa)
            
            !----- Mineralized pool. ------------------------------------------------------!
            csite%mineralized_N_input = fast_N_loss + slow_C_loss / c2n_slow
            csite%mineralized_N_loss  = csite%total_plant_nitrogen_uptake(ipa)             &
                                      + csite%today_Af_decomp(ipa) * Lc * K1               &
                                      * csite%structural_soil_C(ipa)                       &
                                      * ( (1.0 - r_stsc) / c2n_slow - 1.0 / c2n_structural)


            !------------------------------------------------------------------------------!
            !      All carbon fluxes have units kgC/m2/day, and we are updating on the     !
            ! daily time step.                                                             !
            !------------------------------------------------------------------------------!
            csite%fast_soil_C(ipa)       = csite%fast_soil_C(ipa) + csite%fsc_in(ipa)      &
                                         - fast_C_loss
            csite%structural_soil_C(ipa) = csite%structural_soil_C(ipa)                    &
                                         + csite%ssc_in(ipa) - structural_C_loss
            csite%structural_soil_L(ipa) = csite%structural_soil_L(ipa)                    &
                                         + csite%ssl_in(ipa) - structural_L_loss
            csite%slow_soil_C(ipa)       = csite%slow_soil_C(ipa) + slow_C_input           &
                                         - slow_C_loss
            
            !------------------------------------------------------------------------------!
            !      All nitrogen fluxes have units kgN/m2/day, and we are updating on the   !
            ! daily time step.                                                             !
            !------------------------------------------------------------------------------!
            csite%fast_soil_N(ipa)        = csite%fast_soil_N(ipa) + csite%fsn_in(ipa)     &
                                          - fast_N_loss
            csite%mineralized_soil_N(ipa) = csite%mineralized_soil_N(ipa)                  &
                                          + csite%mineralized_N_input(ipa)                 &
                                          - csite%mineralized_N_loss(ipa)

            !------------------------------------------------------------------------------!
            !      Force all pools to be either zero or positive.                          !
            !------------------------------------------------------------------------------!
            csite%fast_soil_C(ipa)        = max(0.0,csite%fast_soil_C(ipa))
            csite%structural_soil_C(ipa)  = max(0.0,csite%structural_soil_C(ipa))
            csite%structural_soil_L(ipa)  = max(0.0,csite%structural_soil_L(ipa))
            csite%slow_soil_C(ipa)        = max(0.0,csite%slow_soil_C(ipa))
            csite%fast_soil_N(ipa)        = max(0.0,csite%fast_soil_N(ipa))
            csite%mineralized_soil_N(ipa) = max(0.0,csite%mineralized_soil_N(ipa))
            
         end do patchloop
      end do siteloop
   end do polygonloop
   
   return
end subroutine update_C_and_N_pools
!==========================================================================================!
!==========================================================================================!
