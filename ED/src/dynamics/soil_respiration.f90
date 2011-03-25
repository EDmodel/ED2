!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the soil respiration terms (root and heterotrophic).        !
!------------------------------------------------------------------------------------------!
subroutine soil_respiration(csite,ipa,mzg,ntext_soil)

   use ed_state_vars, only : sitetype                 & ! structure
                           , patchtype                ! ! structure
   use soil_coms    , only : soil                     ! ! intent(in)
   use pft_coms     , only : root_respiration_factor  ! ! intent(in)

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
   real                                       :: r_resp_temp_fac
   real                                       :: Lc
   real                                       :: r_resp
   !----- External functions. -------------------------------------------------------------!
   real                          , external   :: resp_weight
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      This is the temperature dependence of root respiration.  Same for all cohorts.   !
   !---------------------------------------------------------------------------------------!
   r_resp_temp_fac = 1.0                                                                   &
                   / (1.0 + exp(0.4 * ( 278.15 - csite%soil_tempk(mzg,ipa) ) ) )           &
                   / (1.0 + exp(0.4 * ( csite%soil_tempk(mzg,ipa) - 318.15 ) ) )           &
                   * exp( 10.41 - 3000.0/csite%soil_tempk(mzg,ipa) )

   cpatch => csite%patch(ipa)
   do ico = 1,cpatch%ncohorts
      ipft = cpatch%pft(ico)
      r_resp = root_respiration_factor(ipft) * r_resp_temp_fac * cpatch%broot(ico)        &
             * cpatch%nplant(ico)
      cpatch%root_respiration(ico) = r_resp
      cpatch%mean_root_resp(ico)   = cpatch%mean_root_resp(ico)  + r_resp
      cpatch%today_root_resp(ico)  = cpatch%today_root_resp(ico) + r_resp
   end do

   !----- Compute soil/temperature modulation of heterotrophic respiration. ---------------!
   csite%A_decomp(ipa) = resp_weight(csite%soil_tempk(mzg,ipa),csite%soil_water(mzg,ipa)   &
                                    ,soil(ntext_soil(mzg))%slmsts)

   !----- Compute nitrogen immobilization factor. -----------------------------------------!
   call resp_f_decomp(csite,ipa, Lc)

   !----- Compute heterotrophic respiration. ----------------------------------------------!
   call resp_rh(csite,ipa, Lc)

   !----- Update averaged variables. ------------------------------------------------------!
   csite%today_A_decomp(ipa)  = csite%today_A_decomp(ipa)  + csite%A_decomp(ipa)
   csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa)                                 &
                              + csite%A_decomp(ipa) * csite%f_decomp(ipa)
   csite%mean_rh(ipa)         = csite%mean_rh(ipa) + csite%rh(ipa)

   return
end subroutine soil_respiration
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function computes the respiration limitation factor, which includes limitations !
! due to temperature and moisture.                                                         !
!------------------------------------------------------------------------------------------!
real function resp_weight(soil_tempk,soil_water,slmsts)

   use decomp_coms, only : resp_temperature_increase  & ! intent(in)
                         , resp_opt_water             & ! intent(in)
                         , resp_water_below_opt       & ! intent(in)
                         , resp_water_above_opt       & ! intent(in)
                         , LloydTaylor                ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in) :: soil_tempk
   real, intent(in) :: soil_water
   real, intent(in) :: slmsts
   !----- Local variables. ----------------------------------------------------------------!
   real             :: temperature_limitation
   real             :: water_limitation
   real             :: rel_soil_moist
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the temperature dependence.                                                  !
   !---------------------------------------------------------------------------------------!
   if (LloydTaylor) then 
      !----- Use Lloyd and Taylor (1994) temperature dependence. --------------------------!
      temperature_limitation = min( 1.0                                                    &
                                  , resp_temperature_increase                              &
                                  * exp(308.56 * (1./56.02 - 1./(soil_tempk-227.15)) ) )
   else 
      !----- Use original exponential temperature dependence. -----------------------------!
      temperature_limitation = min( 1.0                                                    &
                                  , exp( resp_temperature_increase * (soil_tempk-318.15)))
   end if


   !---------------------------------------------------------------------------------------!
   !     Find the relative soil moisture, then the moisture dependence.                    !
   !---------------------------------------------------------------------------------------!
   rel_soil_moist = soil_water/slmsts
   if (rel_soil_moist <= resp_opt_water)then
      water_limitation = exp( (rel_soil_moist - resp_opt_water) * resp_water_below_opt)
   else
      water_limitation = exp( (resp_opt_water - rel_soil_moist) * resp_water_above_opt)
   end if
   
   !----- Compute the weight, which is just the combination of both. ----------------------!
   resp_weight = temperature_limitation * water_limitation
      
   return
end function resp_weight
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
   type(sitetype), target      :: csite
   integer       , intent(in)  :: ipa
   real          , intent(in)  :: Lc
   !----- Local variables. ----------------------------------------------------------------!
   real                        :: fast_C_loss
   real                        :: structural_C_loss
   real                        :: slow_C_loss
   !---------------------------------------------------------------------------------------!



   !----- The following variables have units of [kgC/m2/day]. -----------------------------!
   fast_C_loss       = csite%A_decomp(ipa) * K2 * csite%fast_soil_C(ipa)
   structural_C_loss = csite%A_decomp(ipa) * Lc * K1 * csite%structural_soil_C(ipa)        &
                     * csite%f_decomp(ipa)
   slow_C_loss       = csite%A_decomp(ipa) * K3 * csite%slow_soil_C(ipa)

   !----- The following variables have units of [umol_CO2/m2/s]. --------------------------!
   csite%rh(ipa)     = kgCday_2_umols * ( r_fsc * fast_C_loss + r_stsc * structural_C_loss &
                                        + r_ssc * slow_C_loss)                             
   csite%cwd_rh(ipa) = kgCday_2_umols * (r_stsc * structural_C_loss + r_ssc * slow_C_loss) &
                     * cwd_frac

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
