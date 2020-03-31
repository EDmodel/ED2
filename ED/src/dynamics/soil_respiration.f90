module soil_respiration
  contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the soil respiration terms (root and heterotrophic).     !
   !---------------------------------------------------------------------------------------!
   subroutine soil_respiration_driver(csite,ipa,mzg,ntext_soil)

      use ed_state_vars, only : sitetype                 & ! structure
                              , patchtype                ! ! structure
      use soil_coms    , only : soil                     & ! intent(in)
                              , dslz                     & ! intent(in)
                              , slz                      ! ! intent(in)
      use decomp_coms  , only : decomp_scheme            & ! intent(in)
                              , k_rh_active              ! ! intent(in)
      use consts_coms  , only : wdns                     & ! intent(in)
                              , umols_2_kgCyr            ! ! intent(in)
      use therm_lib    , only : uextcm2tl                ! ! function
      use ed_misc_coms , only : dtlsm                    & ! intent(in)
                              , dtlsm_o_frqsum           ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                , target     :: csite
      integer                       , intent(in) :: ipa
      integer                       , intent(in) :: mzg
      integer       , dimension(mzg), intent(in) :: ntext_soil
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)               , pointer    :: cpatch
      integer                                    :: ico
      integer                                    :: ipft
      integer                                    :: k
      integer                                    :: kroot
      integer                                    :: nsoil
      real                                       :: lyr_soil_oxygen
      real                                       :: lyr_soil_moist
      real                                       :: rel_soil_oxygen
      real                                       :: rel_soil_moist
      real                                       :: sum_soil_energy
      real                                       :: sum_soil_water
      real                                       :: sum_soil_hcap
      real                                       :: sum_soil_slmsts
      real                                       :: sum_soil_soilcp
      real                                       :: avg_soil_temp
      real                                       :: avg_soil_fliq
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the current root respiration.  This is done cohort by cohort because the !
      ! parameters may be different depending on the PFT, and also because each layer has  !
      ! a different temperature.                                                           !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         ipft  = cpatch%pft(ico)
         kroot = cpatch%krdepth(ico)

         !---------------------------------------------------------------------------------!
         !    Add "intensive" contribution of each layer, assuming that the roots are      !
         ! equally spread throughout the entire depth.                                     !
         !---------------------------------------------------------------------------------!
         cpatch%root_respiration(ico) = 0.0
         do k = kroot,mzg
            cpatch%root_respiration(ico) = cpatch%root_respiration(ico)                    &
                                         + root_resp_norm(ipft,csite%soil_tempk(k,ipa))    &
                                         * dslz(k)
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Now we make the value in umol/m2/s, by dividing by the total depth and     !
         ! multiplying by the total root biomass.  The minus sign is because slz is        !
         ! negative.                                                                       !
         !---------------------------------------------------------------------------------!
         cpatch%root_respiration(ico) = - cpatch%root_respiration(ico) * cpatch%broot(ico) &
                                        * cpatch%nplant(ico) / slz(kroot)
         !---------------------------------------------------------------------------------!


         !----- Add this time step to the daily mean root respiration. --------------------!
         cpatch%today_root_resp(ico)  = cpatch%today_root_resp(ico)                        &
                                      + cpatch%root_respiration(ico) * dtlsm
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     The following is for output only, we switch the units to kgC/plant/yr.      !
         !---------------------------------------------------------------------------------!
         cpatch%fmean_root_resp(ico)   = cpatch%fmean_root_resp (ico)                      &
                                       + cpatch%root_respiration(ico)                      &
                                       * umols_2_kgCyr * dtlsm_o_frqsum                    &
                                       / cpatch%nplant          (ico)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!




      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !     Find the scaling factor for decomposition of surface materials, using the top  !
      ! soil layer to define their temperature and moisture.                               !
      !------------------------------------------------------------------------------------!
      k                   = mzg
      nsoil               = ntext_soil(k)

      !------------------------------------------------------------------------------------!
      !   Relative moisture may be defined in terms of moisture or potential, depending on !
      ! the decomposition module.                                                          !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (5)
         !---------------------------------------------------------------------------------!
         !    Use (logarithmic) matric potential to define relative moisture, following    !
         ! K13.                                                                            !
         !                                                                                 !
         !   Koven CD, Riley WJ, Subin ZM, Tang JY, Torn MS, Collins WD, Bonan GB,         !
         !      Lawrence DM, Swenson SC. 2013. The effect of vertically resolved soil      !
         !      biogeochemistry and alternate soil C and N models on C dynamics of CLM4.   !
         !      Biogeosciences, 10:7109-7131. doi:10.5194/bg-10-7109-2013.                 !
         !---------------------------------------------------------------------------------!
         rel_soil_moist = min(1.0, max( 0.0                                                &
                                      , log(csite%soil_mstpot(k,ipa)/soil(nsoil)%slpotcp)  &
                                      / log(soil(nsoil)%slpotfc     /soil(nsoil)%slpotcp) ))
         !---------------------------------------------------------------------------------!
      case default
         !------ Use soil moisture. -------------------------------------------------------!
         rel_soil_moist = min(1.0, max( 0.0                                                &
                                      , (csite%soil_water(k,ipa) - soil(nsoil)%soilcp)     &
                                      / (soil(nsoil)%slmsts      - soil(nsoil)%soilcp) ) )
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      When soil moisture is excessively high (i.e., above field capacity), oxygen   !
      ! may become a limiting factor for decomposition.  Most approaches account for this  !
      ! as a single term in the the moisture function, but decomp_scheme=5 treats these    !
      ! as separate limitations.                                                           !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (5)
         !---------------------------------------------------------------------------------!
         !    Use (logarithmic) matric potential to define relative moisture, following    !
         ! K13.                                                                            !
         !                                                                                 !
         !   Koven CD, Riley WJ, Subin ZM, Tang JY, Torn MS, Collins WD, Bonan GB,         !
         !      Lawrence DM, Swenson SC. 2013. The effect of vertically resolved soil      !
         !      biogeochemistry and alternate soil C and N models on C dynamics of CLM4.   !
         !      Biogeosciences, 10:7109-7131. doi:10.5194/bg-10-7109-2013.                 !
         !---------------------------------------------------------------------------------!
         rel_soil_oxygen = min(1.0, max( 0.0                                               &
                                       , log(csite%soil_mstpot(k,ipa)/soil(nsoil)%slpots)  &
                                       / log(soil(nsoil)%slpotfc     /soil(nsoil)%slpots) ))
         !---------------------------------------------------------------------------------!
      case default
         !------ Set as the complement of rel_soil_moist. ---------------------------------!
         rel_soil_oxygen = 1.0 - rel_soil_moist
         !---------------------------------------------------------------------------------!
      end select
      csite%A_decomp(ipa) = het_resp_weight(csite%soil_tempk(k,ipa),rel_soil_moist         &
                                           ,rel_soil_oxygen)
      !------------------------------------------------------------------------------------!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




      !------------------------------------------------------------------------------------!
      !     Integrate the soil extensive properties, plus the minimum and maximum possible !
      ! soil water content of the active layer.                                            !
      !------------------------------------------------------------------------------------!
      sum_soil_energy = 0.0
      sum_soil_hcap   = 0.0
      sum_soil_water  = 0.0
      sum_soil_slmsts = 0.0
      sum_soil_soilcp = 0.0
      rel_soil_moist  = 0.0
      rel_soil_oxygen = 0.0
      do k = k_rh_active,mzg
         nsoil = ntext_soil(k)

         !---------------------------------------------------------------------------------!
         !    Convert the units so energy is in J/m2, heat capacity in J/m2/K, and water   !
         ! in kg/m2.                                                                       !
         !---------------------------------------------------------------------------------!
         sum_soil_energy = sum_soil_energy + csite%soil_energy(k,ipa)        * dslz(k)
         sum_soil_hcap   = sum_soil_hcap   + soil(nsoil)%slcpd               * dslz(k)
         sum_soil_water  = sum_soil_water  + csite%soil_water (k,ipa) * wdns * dslz(k)
         sum_soil_slmsts = sum_soil_slmsts + soil(nsoil)%slmsts       * wdns * dslz(k)
         sum_soil_soilcp = sum_soil_soilcp + soil(nsoil)%soilcp       * wdns * dslz(k)
         !---------------------------------------------------------------------------------!



         !----- Find the relative soil moisture and soil "oxygen". ------------------------!
         select case (decomp_scheme)
         case (5)
            !------ Compute relative value for layer. -------------------------------------!
            lyr_soil_moist  = log(csite%soil_mstpot(k,ipa)/soil(nsoil)%slpotcp)            &
                            / log(soil(nsoil)%slpotfc     /soil(nsoil)%slpotcp)
            lyr_soil_oxygen = log(csite%soil_mstpot(k,ipa)/soil(nsoil)%slpots )            &
                            / log(soil(nsoil)%slpotfc     /soil(nsoil)%slpots )
            lyr_soil_moist  = max(0.,min(1.,lyr_soil_moist ))
            lyr_soil_oxygen = max(0.,min(1.,lyr_soil_oxygen))
            rel_soil_moist  = rel_soil_moist  + lyr_soil_moist  * dslz(k)
            rel_soil_oxygen = rel_soil_oxygen + lyr_soil_oxygen * dslz(k)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !----- Find the average temperature and the relative soil moisture. -----------------!
      call uextcm2tl(sum_soil_energy,sum_soil_water,sum_soil_hcap                          &
                    ,avg_soil_temp,avg_soil_fliq)
      select case (decomp_scheme)
      case (5)
         !------ Normalise relative soil moisture. ----------------------------------------!
         rel_soil_moist  = - rel_soil_moist  / slz(k_rh_active)
         rel_soil_oxygen = - rel_soil_oxygen / slz(k_rh_active)
         !---------------------------------------------------------------------------------!
      case default
         !------ Relative soil moisture based on total water content. ---------------------!
         rel_soil_moist  = min( 1.0, max(0.0, ( sum_soil_water  - sum_soil_soilcp )        &
                                            / ( sum_soil_slmsts - sum_soil_soilcp ) ) )
         rel_soil_oxygen = 1.0 - rel_soil_moist
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Compute soil/temperature modulation of soil heterotrophic respiration. -------!
      csite%B_decomp(ipa) = het_resp_weight(avg_soil_temp,rel_soil_moist,rel_soil_oxygen)
      !------------------------------------------------------------------------------------!



      !----- Compute nitrogen immobilization factor. --------------------------------------!
      call resp_f_decomp(csite,ipa)
      !------------------------------------------------------------------------------------!



      !----- Compute heterotrophic respiration. -------------------------------------------!
      nsoil = ntext_soil(mzg)
      call resp_rh(csite,ipa,nsoil)
      !------------------------------------------------------------------------------------!



      !----- The output is converted to kgC/m2/yr. ----------------------------------------!
      csite%fmean_rh     (ipa) = csite%fmean_rh     (ipa)                                  &
                               + csite%rh           (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_fgc_rh (ipa) = csite%fmean_fgc_rh (ipa)                                  &
                               + csite%fgc_rh       (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_fsc_rh (ipa) = csite%fmean_fsc_rh (ipa)                                  &
                               + csite%fsc_rh       (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_stgc_rh(ipa) = csite%fmean_stgc_rh(ipa)                                  &
                               + csite%stgc_rh      (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_stsc_rh(ipa) = csite%fmean_stsc_rh(ipa)                                  &
                               + csite%stsc_rh      (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_msc_rh (ipa) = csite%fmean_msc_rh (ipa)                                  &
                               + csite%msc_rh       (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_ssc_rh (ipa) = csite%fmean_ssc_rh (ipa)                                  &
                               + csite%ssc_rh       (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      csite%fmean_psc_rh (ipa) = csite%fmean_psc_rh (ipa)                                  &
                               + csite%psc_rh       (ipa) * umols_2_kgCyr * dtlsm_o_frqsum
      !------------------------------------------------------------------------------------!

      return
   end subroutine soil_respiration_driver
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the Nitrogen immobilization factor.                      !
   !---------------------------------------------------------------------------------------!
   subroutine resp_f_decomp(csite,ipa)

      use ed_state_vars, only : sitetype               ! ! structure
      use decomp_coms  , only : e_lignin               & ! intent(in)
                              , r_stsc_o               & ! intent(in)
                              , r_stsc_l               & ! intent(in)
                              , N_immobil_supply_scale & ! intent(in)
                              , decay_rate_stsc        & ! intent(in)
                              , n_decomp_lim           & ! intent(in)
                              , c2n_structural         & ! intent(in)
                              , c2n_slow               ! ! intent(in)
      use consts_coms  , only : lnexp_min              & ! intent(in)
                              , lnexp_max              ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target      :: csite
      integer       , intent(in)  :: ipa
      !----- Local variables. -------------------------------------------------------------!
      real                        :: N_immobilization_demand
      real                        :: ln_Lc
      real                        :: fg_lignin
      real                        :: fs_lignin
      real                        :: er_stgc
      real                        :: er_stsc
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the penalty factor for structural composition due to lignin content.     !
      !------------------------------------------------------------------------------------!
      if (csite%structural_grnd_C(ipa) > 0.0) then
         fg_lignin = csite%structural_grnd_L(ipa) / csite%structural_grnd_C(ipa)
      else
         fg_lignin = 1.0
      end if
      ln_Lc = max(lnexp_min,min(lnexp_max,-e_lignin * fg_lignin))
      csite%Lg_decomp(ipa) = exp(ln_Lc)
      if (csite%structural_soil_C(ipa) > 0.0) then
         fs_lignin = csite%structural_soil_L(ipa) / csite%structural_soil_C(ipa)
      else
         fs_lignin = 1.0
      end if
      ln_Lc = max(lnexp_min,min(lnexp_max,-e_lignin * fs_lignin))
      csite%Ls_decomp(ipa) = exp(ln_Lc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the respiration rate (this may be double counting the penalty).          !
      !------------------------------------------------------------------------------------!
      er_stgc = (1. - fg_lignin) * r_stsc_o + fg_lignin * r_stsc_l
      er_stsc = (1. - fs_lignin) * r_stsc_o + fs_lignin * r_stsc_l
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the penalty factor due to N limitation.                                  !
      !------------------------------------------------------------------------------------!
      select case (n_decomp_lim)
      case (1)
         N_immobilization_demand = csite%A_decomp(ipa) * decay_rate_stsc                   &
                                 * ( (1.0 - er_stgc) / c2n_slow - 1.0 / c2n_structural   ) &
                                 * ( csite%Lg_decomp(ipa) * csite%structural_grnd_C(ipa) ) &
                                 + csite%B_decomp(ipa) * decay_rate_stsc                   &
                                 * ( (1.0 - er_stsc) / c2n_slow - 1.0 / c2n_structural   ) &
                                 * ( csite%Ls_decomp(ipa) * csite%structural_soil_C(ipa) )
         csite%f_decomp(ipa)     = N_immobil_supply_scale * csite%mineralized_soil_N(ipa)  &
                                 / ( N_immobilization_demand                               &
                                   + N_immobil_supply_scale                                &
                                   * csite%mineralized_soil_N(ipa))
      case default
         !----- Option for no plant N limitation. -----------------------------------------!
         csite%f_decomp(ipa)     = 1.0
      end select
      !------------------------------------------------------------------------------------!




      !----- Update averaged variables. ---------------------------------------------------!
      csite%today_A_decomp (ipa) = csite%today_A_decomp (ipa) + csite%A_decomp(ipa)
      csite%today_B_decomp (ipa) = csite%today_B_decomp (ipa) + csite%B_decomp(ipa)
      csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa)                              &
                                 + csite%A_decomp       (ipa) * csite%f_decomp(ipa)
      csite%today_Bf_decomp(ipa) = csite%today_Bf_decomp(ipa)                              &
                                 + csite%B_decomp       (ipa) * csite%f_decomp(ipa)
      !------------------------------------------------------------------------------------!

      return
   end subroutine resp_f_decomp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the heterotrophic respiration.                           !
   !---------------------------------------------------------------------------------------!
   subroutine resp_rh(csite,ipa,ntext)

      use ed_state_vars, only : sitetype        ! ! structure
      use ed_misc_coms , only : dtlsm           & ! intent(in)
                              , current_time    ! ! intent(in)
      use consts_coms  , only : kgCday_2_umols  & ! intent(in)
                              , umols_2_kgCyr   ! ! intent(in)
      use soil_coms    , only : soil            ! ! look-up table
      use decomp_coms  , only : decomp_scheme   & ! intent(in)
                              , decay_rate_fsc  & ! intent(in)
                              , decay_rate_stsc & ! intent(in)
                              , decay_rate_msc  & ! intent(in)
                              , decay_rate_ssc  & ! intent(in)
                              , decay_rate_psc  & ! intent(in)
                              , r_fsc           & ! intent(in)
                              , r_stsc_l        & ! intent(in)
                              , r_stsc_o        & ! intent(in)
                              , r_msc_int       & ! intent(in)
                              , r_msc_slp       & ! intent(in)
                              , r_ssc           & ! intent(in)
                              , r_psc           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)   , target     :: csite
      integer          , intent(in) :: ipa
      integer          , intent(in) :: ntext
      !----- Local variables. -------------------------------------------------------------!
      real                          :: fg_C_loss
      real                          :: fs_C_loss
      real                          :: fg_N_loss
      real                          :: fs_N_loss
      real                          :: stg_C_loss
      real                          :: sts_C_loss
      real                          :: stg_L_loss
      real                          :: sts_L_loss
      real                          :: stg_N_loss
      real                          :: sts_N_loss
      real                          :: ms_C_loss
      real                          :: ss_C_loss
      real                          :: ps_C_loss
      real                          :: fg_C_resp
      real                          :: fs_C_resp
      real                          :: stg_C_resp
      real                          :: sts_C_resp
      real                          :: ms_C_resp
      real                          :: ss_C_resp
      real                          :: ps_C_resp
      real                          :: tot_C_resp
      real                          :: er_fsc
      real                          :: er_stgc
      real                          :: er_stsc
      real                          :: er_msc
      real                          :: er_ssc
      real                          :: er_psc
      real                          :: fg_lignin
      real                          :: fs_lignin
      real                          :: AfL_decomp
      real                          :: BfL_decomp
      !----- Local constants. -------------------------------------------------------------!
      logical          , parameter  :: print_debug = .false.
      character(len=12), parameter  :: rhetfile    = 'het_resp.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      logical          , save       :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (print_debug) then
            open (unit=84,file=rhetfile,status='replace',action='write')
            write (unit=84,fmt='(31(a,1x))')                                               &
                     '  YEAR',      ' MONTH',      '   DAY',      '  HOUR',      '   MIN'  &
              ,      '   IPA','       C_FGC','       C_FSC','      C_STGC','      C_STSC'  &
              ,'       C_MSC','       C_SSC','       C_PSC','      F_FAST','   FG_STRUCT'  &
              ,'   FS_STRUCT','   F_MICROBE','      F_SLOW','      RH_FGC','      RH_FSC'  &
              ,'     RH_STGC','     RH_STSC','      RH_MSC','      RH_SSC','      RH_PSC'  &
              ,'    RH_TOTAL','    A_DECOMP','    B_DECOMP','    F_DECOMP','   LG_DECOMP'  &
              ,'   LS_DECOMP'
            close (unit=84,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!


      !----- Decay rate for structural carbon. --------------------------------------------!
      AfL_decomp = csite%A_decomp(ipa) * csite%f_decomp(ipa) * csite%Lg_decomp(ipa)
      BfL_decomp = csite%B_decomp(ipa) * csite%f_decomp(ipa) * csite%Ls_decomp(ipa)
      !------------------------------------------------------------------------------------!


      !----- The following variables have units of [kgC/m2/day]. --------------------------!
      fg_C_loss  = csite%A_decomp(ipa) * decay_rate_fsc  * csite%fast_grnd_C       (ipa)
      fs_C_loss  = csite%B_decomp(ipa) * decay_rate_fsc  * csite%fast_soil_C       (ipa)
      fg_N_loss  = csite%A_decomp(ipa) * decay_rate_fsc  * csite%fast_grnd_N       (ipa)
      fs_N_loss  = csite%B_decomp(ipa) * decay_rate_fsc  * csite%fast_soil_N       (ipa)
      stg_C_loss = AfL_decomp          * decay_rate_stsc * csite%structural_grnd_C (ipa)
      sts_C_loss = BfL_decomp          * decay_rate_stsc * csite%structural_soil_C (ipa)
      stg_L_loss = AfL_decomp          * decay_rate_stsc * csite%structural_grnd_L (ipa)
      sts_L_loss = BfL_decomp          * decay_rate_stsc * csite%structural_soil_L (ipa)
      stg_N_loss = AfL_decomp          * decay_rate_stsc * csite%structural_grnd_N (ipa)
      sts_N_loss = BfL_decomp          * decay_rate_stsc * csite%structural_soil_N (ipa)
      ms_C_loss  = csite%B_decomp(ipa) * decay_rate_msc  * csite%microbial_soil_C  (ipa)
      ss_C_loss  = csite%B_decomp(ipa) * decay_rate_ssc  * csite%slow_soil_C       (ipa)
      ps_C_loss  = csite%B_decomp(ipa) * decay_rate_psc  * csite%passive_soil_C    (ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the penalty factor for structural composition due to lignin content.     !
      !------------------------------------------------------------------------------------!
      if (csite%structural_grnd_C(ipa) > 0.0) then
         fg_lignin = csite%structural_grnd_L(ipa) / csite%structural_grnd_C(ipa)
      else
         fg_lignin = 1.0
      end if
      if (csite%structural_soil_C(ipa) > 0.0) then
         fs_lignin = csite%structural_soil_L(ipa) / csite%structural_soil_C(ipa)
      else
         fs_lignin = 1.0
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the respiration rate (this may be double counting the penalty).          !
      !------------------------------------------------------------------------------------!
      er_stgc = (1. - fg_lignin) * r_stsc_o + fg_lignin * r_stsc_l
      er_stsc = (1. - fs_lignin) * r_stsc_o + fs_lignin * r_stsc_l
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide the respiration factors based on the method.                           !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (5)
         !---------------------------------------------------------------------------------!
         !     Five-pool CENTURY model.  We solve transfer between soil pools, so they     !
         ! can all have sub-unity fractions that go to heterotrophic respiration.          !
         !---------------------------------------------------------------------------------!
         er_fsc  = r_fsc
         er_msc  = r_msc_int + r_msc_slp * soil(ntext)%xsand
         er_ssc  = r_ssc
         er_psc  = r_psc
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !    Three-pool CENTURY model (ED-2.0 style).  Slow pool must respire 100% of the !
         ! decomposition.                                                                  !
         !---------------------------------------------------------------------------------!
         er_fsc  = r_fsc
         er_msc  = 1.0    ! Dummy value as microbial pool doesn't exist in this scheme.
         er_ssc  = 1.0    ! This has to be 1.0 because it is the only outlet for SSC.
         er_psc  = 1.0    ! Dummy value as passive pool doesn't exist in this scheme.
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Find the heterotrophic respiration (total and components). -------------------!
      csite%fgc_rh (ipa) = er_fsc  * kgCday_2_umols * fg_C_loss
      csite%fsc_rh (ipa) = er_fsc  * kgCday_2_umols * fs_C_loss
      csite%stgc_rh(ipa) = er_stgc * kgCday_2_umols * stg_C_loss
      csite%stsc_rh(ipa) = er_stsc * kgCday_2_umols * sts_C_loss
      csite%msc_rh (ipa) = er_msc  * kgCday_2_umols * ms_C_loss
      csite%ssc_rh (ipa) = er_ssc  * kgCday_2_umols * ss_C_loss
      csite%psc_rh (ipa) = er_psc  * kgCday_2_umols * ps_C_loss
      csite%rh     (ipa) = csite%fgc_rh (ipa) + csite%fsc_rh (ipa) + csite%stgc_rh(ipa)    &
                         + csite%stsc_rh(ipa) + csite%msc_rh (ipa) + csite%ssc_rh (ipa)    &
                         + csite%psc_rh (ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Save respiration and total carbon loss.  Today_rh is used for carbon           !
      ! conservation verification only, but the daily losses are needed to conserve carbon !
      ! even when patch dynamics occurs (today_A_decomp and related variables may not      !
      ! conserve carbon in case fusion or fission occurs as they are relative rates).      !
      !------------------------------------------------------------------------------------!
      csite%today_rh        (ipa) = csite%today_rh        (ipa) + csite%rh(ipa) * dtlsm
      csite%today_fg_C_loss (ipa) = csite%today_fg_C_loss (ipa) + fg_C_loss     * dtlsm
      csite%today_fs_C_loss (ipa) = csite%today_fs_C_loss (ipa) + fs_C_loss     * dtlsm
      csite%today_fg_N_loss (ipa) = csite%today_fg_N_loss (ipa) + fg_N_loss     * dtlsm
      csite%today_fs_N_loss (ipa) = csite%today_fs_N_loss (ipa) + fs_N_loss     * dtlsm
      csite%today_stg_C_loss(ipa) = csite%today_stg_C_loss(ipa) + stg_C_loss    * dtlsm
      csite%today_sts_C_loss(ipa) = csite%today_sts_C_loss(ipa) + sts_C_loss    * dtlsm
      csite%today_stg_L_loss(ipa) = csite%today_stg_L_loss(ipa) + stg_L_loss    * dtlsm
      csite%today_sts_L_loss(ipa) = csite%today_sts_L_loss(ipa) + sts_L_loss    * dtlsm
      csite%today_stg_N_loss(ipa) = csite%today_stg_N_loss(ipa) + stg_N_loss    * dtlsm
      csite%today_sts_N_loss(ipa) = csite%today_sts_N_loss(ipa) + sts_N_loss    * dtlsm
      csite%today_ms_C_loss (ipa) = csite%today_ms_C_loss (ipa) + ms_C_loss     * dtlsm
      csite%today_ss_C_loss (ipa) = csite%today_ss_C_loss (ipa) + ss_C_loss     * dtlsm
      csite%today_ps_C_loss (ipa) = csite%today_ps_C_loss (ipa) + ps_C_loss     * dtlsm
      !------------------------------------------------------------------------------------!



      !----- Write debugging information. -------------------------------------------------!
      if (print_debug) then
         !----- Convert units for output. -------------------------------------------------!
         fg_C_resp  = umols_2_kgCyr * csite%fgc_rh (ipa)
         fs_C_resp  = umols_2_kgCyr * csite%fsc_rh (ipa)
         stg_C_resp = umols_2_kgCyr * csite%stgc_rh(ipa)
         sts_C_resp = umols_2_kgCyr * csite%stsc_rh(ipa)
         ms_C_resp  = umols_2_kgCyr * csite%msc_rh (ipa)
         ss_C_resp  = umols_2_kgCyr * csite%ssc_rh (ipa)
         ps_C_resp  = umols_2_kgCyr * csite%psc_rh (ipa)
         tot_C_resp = umols_2_kgCyr * csite%rh     (ipa)
         !---------------------------------------------------------------------------------!



         !----- Append step to the output file. -------------------------------------------!
         open (unit=84,file=rhetfile,status='old',position='append',action='write')
         write(unit=84,fmt='(6(i6,1x),25(f12.6,1x))')                                      &
                        current_time%year,current_time%month,current_time%date             &
                       ,current_time%hour,current_time%min,ipa                             &
                       ,csite%fast_grnd_C(ipa),csite%fast_soil_C(ipa)                      &
                       ,csite%structural_grnd_C(ipa),csite%structural_soil_C(ipa)          &
                       ,csite%microbial_soil_C(ipa),csite%slow_soil_C(ipa)                 &
                       ,csite%passive_soil_C(ipa),er_fsc,er_stgc,er_stsc,er_msc,er_ssc     &
                       ,fg_C_resp,fs_C_resp,stg_C_resp,sts_C_resp,ms_C_resp,ss_C_resp      &
                       ,ps_C_resp,tot_C_resp,csite%A_decomp(ipa),csite%B_decomp(ipa)       &
                       ,csite%f_decomp(ipa),csite%Lg_decomp(ipa),csite%Ls_decomp(ipa)
         close (unit=84,status='keep')
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine resp_rh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will update the soil carbon and nitrogen pools.                   !
   !---------------------------------------------------------------------------------------!
   subroutine update_C_and_N_pools(cgrid)
      use ed_state_vars, only : edtype          & ! structure
                              , polygontype     & ! structure
                              , sitetype        ! ! structure
      use grid_coms    , only : nzg             ! ! intent(in)
      use soil_coms    , only : soil            ! ! look-up table
      use decomp_coms  , only : decomp_scheme   & ! intent(in)
                              , r_fsc           & ! intent(in)
                              , r_stsc_o        & ! intent(in)
                              , r_stsc_l        & ! intent(in)
                              , r_msc_int       & ! intent(in)
                              , r_msc_slp       & ! intent(in)
                              , r_ssc           & ! intent(in)
                              , r_psc           & ! intent(in)
                              , fx_msc_psc_int  & ! intent(in)
                              , fx_msc_psc_slp  & ! intent(in)
                              , fx_ssc_psc_int  & ! intent(in)
                              , fx_ssc_psc_slp  & ! intent(in)
                              , c2n_slow        & ! intent(in)
                              , c2n_structural  ! ! intent(in)
      use consts_coms  , only : umol_2_kgC      & ! intent(in)
                              , day_sec         ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target   :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer  :: cpoly
      type(sitetype)   , pointer  :: csite
      integer                     :: ipy
      integer                     :: isi
      integer                     :: ipa
      integer                     :: nsoil
      real                        :: fg_C_loss
      real                        :: fs_C_loss
      real                        :: fg_N_loss
      real                        :: fs_N_loss
      real                        :: stg_C_loss
      real                        :: sts_C_loss
      real                        :: stg_L_loss
      real                        :: sts_L_loss
      real                        :: stg_N_loss
      real                        :: sts_N_loss
      real                        :: stg_N_demand
      real                        :: sts_N_demand
      real                        :: ms_C_input
      real                        :: ms_C_loss
      real                        :: ss_C_input
      real                        :: ss_C_loss
      real                        :: ps_C_input
      real                        :: ps_C_loss
      real                        :: er_fsc
      real                        :: er_stgc
      real                        :: er_stsc
      real                        :: er_msc
      real                        :: er_ssc
      real                        :: er_psc
      real                        :: fg_lignin
      real                        :: fs_lignin
      real                        :: fl_stg
      real                        :: fl_sts
      real                        :: ex_fgc_msc
      real                        :: ex_fgc_ssc
      real                        :: ex_fgc_psc
      real                        :: ex_fsc_msc
      real                        :: ex_fsc_ssc
      real                        :: ex_fsc_psc
      real                        :: ex_stgc_msc
      real                        :: ex_stgc_ssc
      real                        :: ex_stgc_psc
      real                        :: ex_stsc_msc
      real                        :: ex_stsc_ssc
      real                        :: ex_stsc_psc
      real                        :: ex_msc_ssc
      real                        :: ex_msc_psc
      real                        :: ex_ssc_msc
      real                        :: ex_ssc_psc
      real                        :: ex_psc_msc
      real                        :: ex_psc_ssc
      real                        :: fast_grnd_C_in
      real                        :: fast_soil_C_in
      real                        :: structural_grnd_C_in
      real                        :: structural_soil_C_in
      real                        :: microbial_soil_C_in
      real                        :: slow_soil_C_in
      real                        :: passive_soil_C_in
      !------------------------------------------------------------------------------------!

      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         nsoil = cpoly%ntext_soil(nzg,isi)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !      Find the penalty factor for structural composition due to lignin        !
            ! content.                                                                     !
            !------------------------------------------------------------------------------!
            if (csite%structural_grnd_C(ipa) > 0.0) then
               fg_lignin = csite%structural_grnd_L(ipa) / csite%structural_grnd_C(ipa)
            else
               fg_lignin = 1.0
            end if
            if (csite%structural_soil_C(ipa) > 0.0) then
               fs_lignin = csite%structural_soil_L(ipa) / csite%structural_soil_C(ipa)
            else
               fs_lignin = 1.0
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Find the respiration rate (this may be double counting the penalty).    !
            !------------------------------------------------------------------------------!
            er_stgc = (1. - fg_lignin) * r_stsc_o + fg_lignin * r_stsc_l
            er_stsc = (1. - fs_lignin) * r_stsc_o + fs_lignin * r_stsc_l
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Decide the respiration factors based on the method.                     !
            !------------------------------------------------------------------------------!
            select case (decomp_scheme)
            case (5)
               !---------------------------------------------------------------------------!
               !     Five-pool CENTURY model.  We solve transfer between soil pools, so    !
               ! they can all have sub-unity fractions that go to heterotrophic            !
               ! respiration.                                                              !
               !---------------------------------------------------------------------------!
               er_fsc  = r_fsc
               er_msc  = r_msc_int + r_msc_slp * soil(nsoil)%xsand
               er_ssc  = r_ssc
               er_psc  = r_psc
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Find transfer rates between soil pools.                               !
               !---------------------------------------------------------------------------!
               ex_fgc_msc = (1.0 - er_fsc)
               ex_fgc_ssc = 0.0
               ex_fgc_psc = 0.0

               ex_fsc_msc = (1.0 - er_fsc)
               ex_fsc_ssc = 0.0
               ex_fsc_psc = 0.0

               ex_msc_psc = min( 1.0 - er_msc                                              &
                               , fx_msc_psc_int + fx_msc_psc_slp * soil(nsoil)%xclay )
               ex_msc_ssc = 1.0 - er_msc - ex_msc_psc

               ex_ssc_psc = min( 1.0 - er_ssc                                              &
                               , fx_ssc_psc_int + fx_ssc_psc_slp * soil(nsoil)%xclay )
               ex_ssc_msc = 1.0 - er_ssc - ex_ssc_psc

               ex_psc_ssc = 1.0 - er_psc
               ex_psc_msc = 0.0
               !---------------------------------------------------------------------------!
            case default
               !---------------------------------------------------------------------------!
               !    Three-pool CENTURY model (ED-2.0 style).  Slow pool must respire 100%  !
               ! of the decomposition.                                                     !
               !---------------------------------------------------------------------------!
               er_fsc  = r_fsc
               er_msc  = 1.0    ! Microbial pool doesn't exist in this scheme.
               er_ssc  = 1.0    ! Respiration is the only outlet for SSC
               er_psc  = 1.0    ! Passive pool doesn't exist in this scheme.
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     No transfer between soil pools.                                       !
               !---------------------------------------------------------------------------!
               ex_fgc_msc  = 0.0
               ex_fgc_ssc  = (1.0 - er_fsc )
               ex_fgc_psc  = 0.0

               ex_fsc_msc  = 0.0
               ex_fsc_ssc  = (1.0 - er_fsc )
               ex_fsc_psc  = 0.0

               ex_stgc_msc = 0.0
               ex_stgc_ssc = (1.0 - er_stgc)
               ex_stgc_psc = 0.0

               ex_stsc_msc = 0.0
               ex_stsc_ssc = (1.0 - er_stsc)
               ex_stsc_psc = 0.0

               ex_msc_ssc  = 0.0
               ex_msc_psc  = 0.0

               ex_ssc_msc  = 0.0
               ex_ssc_psc  = 0.0

               ex_psc_msc  = 0.0
               ex_psc_ssc  = 0.0
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !       Find the transition rates between all pools.                           !
            !------------------------------------------------------------------------------!
            patchloop: do ipa = 1,csite%npatches



               !----- Save soil carbon state prior to the soil carbon dynamics. -----------!
               fast_grnd_C_in       = csite%fast_grnd_C      (ipa)
               fast_soil_C_in       = csite%fast_soil_C      (ipa)
               structural_grnd_C_in = csite%structural_grnd_C(ipa)
               structural_soil_C_in = csite%structural_soil_C(ipa)
               microbial_soil_C_in  = csite%microbial_soil_C (ipa)
               slow_soil_C_in       = csite%slow_soil_C      (ipa)
               passive_soil_C_in    = csite%passive_soil_C   (ipa)
               !---------------------------------------------------------------------------!



               !----- Fast pools. ---------------------------------------------------------!
               fg_C_loss  = csite%today_fg_C_loss(ipa)
               fs_C_loss  = csite%today_fs_C_loss(ipa)
               fg_N_loss  = csite%today_fg_N_loss(ipa)
               fs_N_loss  = csite%today_fs_N_loss(ipa)
               !---------------------------------------------------------------------------!



               !----- Structural pools. ---------------------------------------------------!
               stg_C_loss   = csite%today_stg_C_loss(ipa)
               sts_C_loss   = csite%today_sts_C_loss(ipa)
               stg_L_loss   = csite%today_stg_L_loss(ipa)
               sts_L_loss   = csite%today_sts_L_loss(ipa)
               stg_N_loss   = csite%today_stg_N_loss(ipa)
               sts_N_loss   = csite%today_sts_N_loss(ipa)
               !----- Demand to maintain stoichiometry (see Moorcroft et al. 2001). -------!
               stg_N_demand = stg_C_loss                                                   &
                            * ( (1.0 - er_stgc) * (1./c2n_slow) - (1./c2n_structural) )
               sts_N_demand = sts_C_loss                                                   &
                            * ( (1.0 - er_stsc) * (1./c2n_slow) - (1./c2n_structural) )
               !---------------------------------------------------------------------------!


               !----- Microbial pool (carbon only). ---------------------------------------!
               ms_C_loss  = csite%today_ms_C_loss(ipa)
               !---------------------------------------------------------------------------!


               !----- Slow pool (carbon only). --------------------------------------------!
               ss_C_loss  = csite%today_ss_C_loss(ipa)
               !---------------------------------------------------------------------------!


               !----- Passive pool (carbon only). -----------------------------------------!
               ps_C_loss  = csite%today_ps_C_loss(ipa)
               !---------------------------------------------------------------------------!




               !----- Mineralized pool. ---------------------------------------------------!
               csite%mineralized_N_input(ipa) = fg_N_loss  + fs_N_loss                     &
                                              + stg_N_loss + sts_N_loss
               csite%mineralized_N_loss (ipa) = csite%total_plant_nitrogen_uptake(ipa)     &
                                              + stg_N_demand + sts_N_demand
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Inputs for the microbial and humified (slow) soil pools depend on the  !
               ! decomposition scheme.                                                     !
               !---------------------------------------------------------------------------!
               select case (decomp_scheme)
               case (5)
                  !------------------------------------------------------------------------!
                  !    Find ratio of decayed structural carbon that goes to microbial and  !
                  ! humified (slow) carbon.                                                !
                  !------------------------------------------------------------------------!
                  if (csite%structural_grnd_C(ipa) > 0.0) then
                     fl_stg = csite%structural_grnd_L(ipa) / csite%structural_grnd_C(ipa)
                  else
                     fl_stg = 1.0
                  end if
                  if (csite%structural_soil_C(ipa) > 0.0) then
                     fl_sts = csite%structural_soil_L(ipa) / csite%structural_soil_C(ipa)
                  else
                     fl_sts = 1.0
                  end if
                  ex_stgc_msc = (1.0 - er_stgc) * (1.0 - fl_stg)
                  ex_stgc_ssc = (1.0 - er_stsc) *        fl_stg
                  ex_stgc_psc = 0.0

                  ex_stsc_msc = (1.0 - er_stsc) * (1.0 - fl_sts)
                  ex_stsc_ssc = (1.0 - er_stsc) *        fl_sts
                  ex_stsc_psc = 0.0
                  !------------------------------------------------------------------------!
               end select
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Microbial pool.  All metabolic litter that is not respired becomes    !
               ! microbial, and passive soil carbon cannot become microbial.               !
               !---------------------------------------------------------------------------!
               ms_C_input = ex_fgc_msc     * fg_C_loss                                     &
                          + ex_fsc_msc     * fs_C_loss                                     &
                          + ex_stgc_msc    * stg_C_loss                                    &
                          + ex_stsc_msc    * sts_C_loss                                    &
                          + ex_ssc_msc     * ss_C_loss                                     &
                          + ex_psc_msc     * ps_C_loss
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Humified (slow) pool.  Lignified structural carbon that doesn't       !
               ! become CO2 becomes humus.  Some passive carbon may revert back to         !
               ! humified, but this should be a small fraction given the decay rate.       !
               !---------------------------------------------------------------------------!
               ss_C_input = ex_fgc_ssc  * fg_C_loss                                        &
                          + ex_fsc_ssc  * fs_C_loss                                        &
                          + ex_stgc_ssc * stg_C_loss                                       &
                          + ex_stsc_ssc * sts_C_loss                                       &
                          + ex_msc_ssc  * ms_C_loss                                        &
                          + ex_psc_ssc  * ps_C_loss
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Passive pool.  Only microbial and humified soil can contribute to     !
               ! the passive soil pool.                                                    !
               !---------------------------------------------------------------------------!
               ps_C_input = ex_fgc_psc  * fg_C_loss                                        &
                          + ex_fsc_psc  * fs_C_loss                                        &
                          + ex_stgc_psc * stg_C_loss                                       &
                          + ex_stsc_psc * sts_C_loss                                       &
                          + ex_msc_psc  * ms_C_loss                                        &
                          + ex_ssc_psc  * ss_C_loss
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      All carbon fluxes have units kgC/m2/day, and we are updating on the  !
               ! daily time step.                                                          !
               !---------------------------------------------------------------------------!
               csite%fast_grnd_C(ipa)       = csite%fast_grnd_C      (ipa)                 &
                                            + csite%fgc_in           (ipa) - fg_C_loss
               csite%fast_soil_C(ipa)       = csite%fast_soil_C      (ipa)                 &
                                            + csite%fsc_in           (ipa) - fs_C_loss
               csite%structural_grnd_C(ipa) = csite%structural_grnd_C(ipa)                 &
                                            + csite%stgc_in          (ipa) - stg_C_loss
               csite%structural_soil_C(ipa) = csite%structural_soil_C(ipa)                 &
                                            + csite%stsc_in          (ipa) - sts_C_loss
               csite%structural_grnd_L(ipa) = csite%structural_grnd_L(ipa)                 &
                                            + csite%stgl_in          (ipa) - stg_L_loss
               csite%structural_soil_L(ipa) = csite%structural_soil_L(ipa)                 &
                                            + csite%stsl_in          (ipa) - sts_L_loss
               csite%microbial_soil_C(ipa)  = csite%microbial_soil_C (ipa)                 &
                                            + ms_C_input                   - ms_C_loss
               csite%slow_soil_C(ipa)       = csite%slow_soil_C      (ipa)                 &
                                            + ss_C_input                   - ss_C_loss
               csite%passive_soil_C(ipa)    = csite%passive_soil_C   (ipa)                 &
                                            + ps_C_input                   - ps_C_loss
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      All nitrogen fluxes have units kgN/m2/day, and we are updating on    !
               ! the daily time step.  Currently the microbial pool does not have          !
               ! nitrogen, but this should be changed in the future.                       !
               !---------------------------------------------------------------------------!
               csite%fast_grnd_N       (ipa) = csite%fast_grnd_N(ipa) + csite%fgn_in(ipa)  &
                                             - fg_N_loss
               csite%fast_soil_N       (ipa) = csite%fast_soil_N(ipa) + csite%fsn_in(ipa)  &
                                             - fs_N_loss
               csite%structural_grnd_N (ipa) = csite%structural_grnd_N(ipa)                &
                                             + csite%stgn_in          (ipa)                &
                                             + stg_N_demand - stg_N_loss
               csite%structural_soil_N (ipa) = csite%structural_soil_N(ipa)                &
                                             + csite%stsn_in          (ipa)                &
                                             + sts_N_demand - sts_N_loss
               csite%mineralized_soil_N(ipa) = csite%mineralized_soil_N (ipa)              &
                                             + csite%mineralized_N_input(ipa)              &
                                             - csite%mineralized_N_loss (ipa)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    As the last step, we make sure that there is no soil C being smuggled  !
               ! or evaded from this patch.  We also check that all pools are non-negative !
               ! (tiny negative values will be fixed, but not large numbers).              !
               !                                                                           !
               ! Disclaimer: this check is currently looking only at soil C, in the future !
               ! we should also include nitrogen checks.                                   !
               !---------------------------------------------------------------------------!
               call check_budget_soilc(csite,ipa,fast_grnd_C_in,fast_soil_C_in             &
                                      ,structural_grnd_C_in,structural_soil_C_in           &
                                      ,microbial_soil_C_in,slow_soil_C_in                  &
                                      ,passive_soil_C_in,er_fsc,er_stgc,er_stsc,er_msc     &
                                      ,er_ssc,er_psc,ex_fgc_msc,ex_fgc_ssc,ex_fgc_psc      &
                                      ,ex_fsc_msc,ex_fsc_ssc,ex_fsc_psc,ex_stgc_msc        &
                                      ,ex_stgc_ssc,ex_stgc_psc,ex_stsc_msc,ex_stsc_ssc     &
                                      ,ex_stsc_psc,ex_msc_ssc,ex_msc_psc,ex_ssc_msc        &
                                      ,ex_ssc_psc,ex_psc_msc,ex_psc_ssc,fg_C_loss          &
                                      ,fs_C_loss,stg_C_loss,sts_C_loss,ms_C_input          &
                                      ,ms_C_loss,ss_C_input,ss_C_loss,ps_C_input,ps_C_loss)
               !---------------------------------------------------------------------------!

            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polygonloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_C_and_N_pools
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine resets the necromass inputs.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine zero_litter_inputs(cgrid)
      use ed_state_vars, only : edtype          & ! structure
                              , polygontype     & ! structure
                              , sitetype        ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target   :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer  :: cpoly
      type(sitetype)   , pointer  :: csite
      integer                     :: ipy
      integer                     :: isi
      integer                     :: ipa
      !------------------------------------------------------------------------------------!

      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)


            !------------------------------------------------------------------------------!
            !       Find the transition rates between all pools.                           !
            !------------------------------------------------------------------------------!
            patchloop: do ipa = 1,csite%npatches
               !---------------------------------------------------------------------------!
               !     This used to be reset in phenology, however, patch and cohort fusion  !
               ! and termination may occur in between the calls and carbon may go          !
               ! unaccounted.  This is the safest place to reset, because the litter       !
               ! inputs have just been sent to the litter pools.                           !
               !---------------------------------------------------------------------------!
               csite%fgc_in (ipa) = 0.0
               csite%fsc_in (ipa) = 0.0
               csite%fgn_in (ipa) = 0.0
               csite%fsn_in (ipa) = 0.0
               csite%stgc_in(ipa) = 0.0
               csite%stsc_in(ipa) = 0.0
               csite%stgl_in(ipa) = 0.0
               csite%stsl_in(ipa) = 0.0
               csite%stgn_in(ipa) = 0.0
               csite%stsn_in(ipa) = 0.0
               !---------------------------------------------------------------------------!

            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polygonloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine zero_litter_inputs
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function determines the normalised root respiration (umol/kgC_fine_root/s)   !
   ! at a given soil layer.                                                                !
   !---------------------------------------------------------------------------------------!
   real function root_resp_norm(ipft,soil_temp)
      use pft_coms       , only : root_respiration_factor  & ! intent(in)
                                , rrf_low_temp             & ! intent(in)
                                , rrf_high_temp            & ! intent(in)
                                , rrf_decay_elow           & ! intent(in)
                                , rrf_decay_ehigh          & ! intent(in)
                                , rrf_hor                  & ! intent(in)
                                , rrf_q10                  ! ! intent(in)
      use farq_leuning   , only : arrhenius                & ! function
                                , collatz                  ! ! function
      use rk4_coms       , only : tiny_offset              ! ! intent(in)
      use physiology_coms, only : iphysiol                 ! ! intent(in)
      use consts_coms    , only : lnexp_min8               & ! intent(in)
                                , lnexp_max8               & ! intent(in)
                                , t008                     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer     , intent(in) :: ipft
      real(kind=4), intent(in) :: soil_temp
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: soil_temp8
      real(kind=8)             :: rrf08
      real(kind=8)             :: rrf_low_temp8
      real(kind=8)             :: rrf_high_temp8
      real(kind=8)             :: rrf_decay_elow8
      real(kind=8)             :: rrf_decay_ehigh8
      real(kind=8)             :: rrf_hor8
      real(kind=8)             :: rrf_q108
      real(kind=8)             :: lnexplow
      real(kind=8)             :: lnexphigh
      real(kind=8)             :: tlow_fun
      real(kind=8)             :: thigh_fun
      real(kind=8)             :: rrf8
      !----- External functions. ----------------------------------------------------------!
      real(kind=4), external   :: sngloff
      !------------------------------------------------------------------------------------!

      !----- Copy some variables to double precision temporaries. -------------------------!
      soil_temp8       = dble(soil_temp                    )
      rrf08            = dble(root_respiration_factor(ipft))
      rrf_low_temp8    = dble(rrf_low_temp           (ipft)) + t008
      rrf_high_temp8   = dble(rrf_high_temp          (ipft)) + t008
      rrf_decay_elow8  = dble(rrf_decay_elow         (ipft))
      rrf_decay_ehigh8 = dble(rrf_decay_ehigh        (ipft))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute the functions that will control the Rrf function for low and high       !
      ! temperature.  In order to avoid floating point exceptions, we check whether the    !
      ! temperature will make the exponential too large or too small.                      !
      !------------------------------------------------------------------------------------!
      !----- Low temperature. -------------------------------------------------------------!
      lnexplow  = rrf_decay_elow8  * (rrf_low_temp8  - soil_temp8)
      lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
      tlow_fun  = 1.d0 +  exp(lnexplow)
      !----- High temperature. ------------------------------------------------------------!
      lnexphigh = rrf_decay_ehigh8 * (soil_temp8 - rrf_high_temp8)
      lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
      thigh_fun = 1.d0 + exp(lnexphigh)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide which functional form to use based on the physiology.  This is just to   !
      ! make it look similar to the leaf respiration respiration.                          !
      !------------------------------------------------------------------------------------!
      select case (iphysiol)
      case (0,1)
         rrf_hor8 = dble(rrf_hor(ipft))
         rrf8     = arrhenius(soil_temp8,rrf08,rrf_hor8) / (tlow_fun * thigh_fun)
      case (2,3)
         rrf_q108 = dble(rrf_q10(ipft))
         rrf8     = collatz(soil_temp8,rrf08,rrf_q108)   / (tlow_fun * thigh_fun)
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Convert result to single precision.                                            !
      !------------------------------------------------------------------------------------!
      root_resp_norm = sngloff(rrf8,tiny_offset)
      !------------------------------------------------------------------------------------!


      return
   end function root_resp_norm
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the heterotrophic respiration limitation factor, which     !
   ! includes limitations due to temperature and soil moisture.                            !
   !---------------------------------------------------------------------------------------!
   real function het_resp_weight(soil_tempk,rel_soil_moist,rel_soil_oxygen)

      use decomp_coms , only : resp_temperature_increase  & ! intent(in)
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
                             , rh_wet_smoist              & ! intent(in)
                             , rh_moyano12_a0             & ! intent(in)
                             , rh_moyano12_a1             & ! intent(in)
                             , rh_moyano12_a2             & ! intent(in)
                             , rh08                       & ! intent(in)
                             , rh_q108                    & ! intent(in)
                             , rh_p_smoist                & ! intent(in)
                             , rh_p_oxygen                ! ! intent(in)
      use farq_leuning, only : collatz                    ! ! function
      use consts_coms , only : lnexp_min                  & ! intent(in)
                             , lnexp_max                  & ! intent(in)
                             , almost_zero                & ! intent(in)
                             , almost_one                 ! ! intent(in)
      use rk4_coms    , only : tiny_offset                ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: soil_tempk
      real(kind=4), intent(in) :: rel_soil_moist
      real(kind=4), intent(in) :: rel_soil_oxygen
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: temperature_scale
      real(kind=4)             :: water_limitation
      real(kind=4)             :: oxygen_limitation
      real(kind=4)             :: lnexplloyd
      real(kind=4)             :: lnexplow
      real(kind=4)             :: lnexphigh
      real(kind=4)             :: tlow_fun
      real(kind=4)             :: thigh_fun
      real(kind=4)             :: lnexpdry
      real(kind=4)             :: lnexpwet
      real(kind=4)             :: smdry_fun
      real(kind=4)             :: smwet_fun
      real(kind=4)             :: rel_smoist_bnd
      real(kind=4)             :: rel_oxygen_bnd
      real(kind=8)             :: soil_tempk8
      real(kind=8)             :: temperature_scale8
      !----- External functions. ----------------------------------------------------------!
      real(kind=4), external   :: sngloff
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     In case moisture and oxygen functions are close to zero or one, simplify the   !
      ! values to avoid degenerate floating point operations.                              !
      !------------------------------------------------------------------------------------!
      !----- Soil moisture. ---------------------------------------------------------------!
      if (rel_soil_moist < almost_zero) then
         rel_smoist_bnd = 0.0
      else if (rel_soil_moist > almost_one) then
         rel_smoist_bnd = 1.0
      else
         rel_smoist_bnd = rel_soil_moist
      end if
      !----- Soil oxygen. -----------------------------------------------------------------!
      if (rel_soil_oxygen < almost_zero) then
         rel_oxygen_bnd = 0.0
      else if (rel_soil_oxygen > almost_one) then
         rel_oxygen_bnd = 1.0
      else
         rel_oxygen_bnd = rel_soil_oxygen
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the temperature dependence.                                               !
      !------------------------------------------------------------------------------------!
      select case(decomp_scheme)
      case (0,3)
         !----- Use original ED-2.1 exponential temperature dependence. -------------------!
         temperature_scale = min( 1.0                                                      &
                                , exp(resp_temperature_increase * (soil_tempk-318.15)))
         !---------------------------------------------------------------------------------!
      case (1,4)
         !----- Use Lloyd and Taylor (1994) temperature dependence. -----------------------!
         lnexplloyd        = rh_lloyd_1 * ( rh_lloyd_2 - 1./(soil_tempk - rh_lloyd_3))
         lnexplloyd        = max(lnexp_min,min(lnexp_max,lnexplloyd))
         temperature_scale = min( 1.0, resp_temperature_increase * exp(lnexplloyd) )
         !---------------------------------------------------------------------------------!
      case (2)
         !---------------------------------------------------------------------------------!
         !      Similar to the original ED-1.0 formulation, which is based on the CENTURY  !
         ! model.  The change in the functional form is to avoid power of negative         !
         ! numbers, but the coefficients were tuned to give a similar curve.               !
         !---------------------------------------------------------------------------------!
         !----- Low temperature limitation. -----------------------------------------------!
         lnexplow          = rh_decay_low * (rh_low_temp - soil_tempk)
         lnexplow          = max(lnexp_min,min(lnexp_max,lnexplow))
         tlow_fun          = 1.0 + exp(lnexplow)
         !----- High temperature limitation. ----------------------------------------------!
         lnexphigh         = rh_decay_high * (soil_tempk - rh_high_temp)
         lnexphigh         = max(lnexp_min,min(lnexp_max,lnexphigh))
         thigh_fun         = 1.0 + exp(lnexphigh)
         !----- Temperature scale is a combination of both. -------------------------------!
         temperature_scale = 1.0 / (tlow_fun * thigh_fun )
         !---------------------------------------------------------------------------------!
      case (5)
         !---------------------------------------------------------------------------------!
         !     Heterotrophic respiration is modulated by a Q10 function, based on K13.     !
         !                                                                                 !
         !  Koven CD, Riley WJ, Subin ZM, Tang JY, Torn MS, Collins WD, Bonan GB, Lawrence !
         !     DM, Swenson SC. 2013. The effect of vertically resolved soil biogeo-        !
         !     chemistry and alternate soil C and N models on C dynamics of CLM4.          !
         !     Biogeosciences, 10: 7109-7131. doi:10.5194/bg-10-7109-2013.                 !
         !---------------------------------------------------------------------------------!
         soil_tempk8        = dble(soil_tempk)
         temperature_scale8 = collatz(soil_tempk8,rh08,rh_q108)
         temperature_scale  = sngloff(temperature_scale8,tiny_offset)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the limitation due to (low) soil moisture.                                !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (0,1)
         !----- ED-2.1 default, also used when decomp_scheme is 1. ------------------------!
         if (rel_soil_moist <= resp_opt_water)then
            water_limitation = exp((rel_smoist_bnd - resp_opt_water) * resp_water_below_opt)
         else
            water_limitation = 1.0
         end if
         !---------------------------------------------------------------------------------!
      case (2)
         !----- Dry soil limitation. ------------------------------------------------------!
         lnexpdry         = rh_decay_dry * (rh_dry_smoist - rel_smoist_bnd)
         lnexpdry         = max(lnexp_min,min(lnexp_max,lnexpdry))
         smdry_fun        = 1.0 + exp(lnexpdry)
         water_limitation = 1.0 / smdry_fun
         !---------------------------------------------------------------------------------!
      case (3,4)
         !---------------------------------------------------------------------------------!
         !   From Jaclyn Matthes: Empirical equation from meta-analysis in M12.            !
         !                                                                                 !
         !  Moyano FE, Vasilyeva N, Bouckaert L, Cook F, Craine J, Curiel Yuste J, Don A,  !
         !     Epron D, Formanek P, Franzluebbers A et al. 2012. The moisture response of  !
         !     soil heterotrophic respiration: interaction with soil properties.           !
         !     Biogeosciences, 9: 1173-1182. doi:10.5194/bg-9-1173-2012 (M12).             !
         !---------------------------------------------------------------------------------!
         water_limitation = rh_moyano12_a0 + rh_moyano12_a1 * rel_smoist_bnd               &
                          + rh_moyano12_a2 + rel_smoist_bnd * rel_smoist_bnd
         !---------------------------------------------------------------------------------!
      case (5)
         !---------------------------------------------------------------------------------!
         !     Moisture function based on K13, but with extra power factor based on some   !
         ! local optimisation.                                                             !
         !                                                                                 !
         !  Koven CD, Riley WJ, Subin ZM, Tang JY, Torn MS, Collins WD, Bonan GB, Lawrence !
         !     DM, Swenson SC. 2013. The effect of vertically resolved soil biogeo-        !
         !     chemistry and alternate soil C and N models on C dynamics of CLM4.          !
         !     Biogeosciences, 10: 7109-7131. doi:10.5194/bg-10-7109-2013 (K13).           !
         !---------------------------------------------------------------------------------!
         water_limitation = rel_smoist_bnd ** rh_p_smoist
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the limitation due to low oxygen (aka too much moisture).                 !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (0,1)
         !----- ED-2.1 default, also used when decomp_scheme is 1. ------------------------!
         if (rel_soil_moist <= resp_opt_water)then
            oxygen_limitation = 1.0
         else
            lnexpwet          = (resp_opt_water - rel_smoist_bnd) * resp_water_above_opt
            lnexpwet          = max(lnexp_min,min(lnexp_max,lnexpwet))
            oxygen_limitation = exp(lnexpwet)
         end if
         !---------------------------------------------------------------------------------!
      case (2)
         !----- Wet soil limitation. ------------------------------------------------------!
         lnexpwet          = rh_decay_wet * (rel_smoist_bnd - rh_wet_smoist)
         lnexpwet          = max(lnexp_min,min(lnexp_max,lnexpwet))
         smwet_fun         = 1.0 + exp(lnexpwet)
         oxygen_limitation = 1.0 / smwet_fun
         !---------------------------------------------------------------------------------!
      case (3,4)
         !------ Not needed, set to one. --------------------------------------------------!
         water_limitation = rh_moyano12_a0 + rh_moyano12_a1 * rel_smoist_bnd               &
                          + rh_moyano12_a2 + rel_soil_moist * rel_smoist_bnd
         !---------------------------------------------------------------------------------!
      case (5)
         !---------------------------------------------------------------------------------!
         !     Oxygen function is similar to the dry limitation (decaying as oxygen        !
         ! becomes limiting as soils become saturated.                                     !
         !---------------------------------------------------------------------------------!
         oxygen_limitation = rel_oxygen_bnd ** rh_p_oxygen
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Compute the weight, which is just the combination of both. -------------------!
      het_resp_weight = temperature_scale * water_limitation * oxygen_limitation
      !------------------------------------------------------------------------------------!

      return
   end function het_resp_weight
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks whether carbon is conserved when soil carbon pools are    !
   ! updated.                                                                              !
   !=======================================================================================!
   !=======================================================================================!
   subroutine check_budget_soilc(csite,ipa,fast_grnd_C_in,fast_soil_C_in                   &
                                ,structural_grnd_C_in,structural_soil_C_in                 &
                                ,microbial_soil_C_in,slow_soil_C_in,passive_soil_C_in      &
                                ,er_fsc,er_stgc,er_stsc,er_msc,er_ssc,er_psc,ex_fgc_msc    &
                                ,ex_fgc_ssc,ex_fgc_psc,ex_fsc_msc,ex_fsc_ssc,ex_fsc_psc    &
                                ,ex_stgc_msc,ex_stgc_ssc,ex_stgc_psc,ex_stsc_msc           &
                                ,ex_stsc_ssc,ex_stsc_psc,ex_msc_ssc,ex_msc_psc,ex_ssc_msc  &
                                ,ex_ssc_psc,ex_psc_msc,ex_psc_ssc,fg_C_loss,fs_C_loss      &
                                ,stg_C_loss,sts_C_loss,ms_C_input,ms_C_loss,ss_C_input     &
                                ,ss_C_loss,ps_C_input,ps_C_loss)

      use ed_state_vars, only : sitetype          ! ! structure
      use ed_misc_coms , only : current_time      ! ! intent(in)
      use budget_utils , only : tol_carbon_budget ! ! intent(in)
      use consts_coms  , only : umol_2_kgC        & ! intent(in)
                              , onethird          & ! intent(in)
                              , day_sec           ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)   , target     :: csite
      integer          , intent(in) :: ipa
      real             , intent(in) :: fast_grnd_C_in
      real             , intent(in) :: fast_soil_C_in
      real             , intent(in) :: structural_grnd_C_in
      real             , intent(in) :: structural_soil_C_in
      real             , intent(in) :: microbial_soil_C_in
      real             , intent(in) :: slow_soil_C_in
      real             , intent(in) :: passive_soil_C_in
      real             , intent(in) :: er_fsc
      real             , intent(in) :: er_stgc
      real             , intent(in) :: er_stsc
      real             , intent(in) :: er_msc
      real             , intent(in) :: er_ssc
      real             , intent(in) :: er_psc
      real             , intent(in) :: ex_fgc_msc
      real             , intent(in) :: ex_fgc_ssc
      real             , intent(in) :: ex_fgc_psc
      real             , intent(in) :: ex_fsc_msc
      real             , intent(in) :: ex_fsc_ssc
      real             , intent(in) :: ex_fsc_psc
      real             , intent(in) :: ex_stgc_msc
      real             , intent(in) :: ex_stgc_ssc
      real             , intent(in) :: ex_stgc_psc
      real             , intent(in) :: ex_stsc_msc
      real             , intent(in) :: ex_stsc_ssc
      real             , intent(in) :: ex_stsc_psc
      real             , intent(in) :: ex_msc_ssc
      real             , intent(in) :: ex_msc_psc
      real             , intent(in) :: ex_ssc_msc
      real             , intent(in) :: ex_ssc_psc
      real             , intent(in) :: ex_psc_msc
      real             , intent(in) :: ex_psc_ssc
      real             , intent(in) :: fg_C_loss
      real             , intent(in) :: fs_C_loss
      real             , intent(in) :: stg_C_loss
      real             , intent(in) :: sts_C_loss
      real             , intent(in) :: ms_C_input
      real             , intent(in) :: ms_C_loss
      real             , intent(in) :: ss_C_input
      real             , intent(in) :: ss_C_loss
      real             , intent(in) :: ps_C_input
      real             , intent(in) :: ps_C_loss
       !----- Local variables. -------------------------------------------------------------!
      real                          :: soil_C_initial
      real                          :: soil_C_final
      real                          :: net_C_input
      real                          :: net_C_loss
      real                          :: delta_soil_C
      real                          :: resid_soil_C
      real                          :: toler_soil_C
      real                          :: fgc_ok_min
      real                          :: fsc_ok_min
      real                          :: stgc_ok_min
      real                          :: stsc_ok_min
      real                          :: msc_ok_min
      real                          :: ssc_ok_min
      real                          :: psc_ok_min
      real                          :: carbon_miss
      real                          :: today_rh
      real                          :: fracloss_rh
      real                          :: nettrans_rh
      real                          :: toler_rh
      logical                       :: neg_soil_C
      logical, dimension(3)         :: rh_violation
      logical                       :: soil_C_violation
      !----- Local constants. -------------------------------------------------------------!
      real             , parameter  :: soil_C_min = 0.1
      real             , parameter  :: rh_min     = 3.e-4
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter  :: fmti='(a,1x,i14)'
      character(len= 9), parameter  :: fmtl='(a,1x,l1)'
      character(len=13), parameter  :: fmtf='(a,1x,es14.7)'
      character(len=27), parameter  :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!


      !----- Find the minimum acceptable soil C stock. ------------------------------------!
      fgc_ok_min  = - tol_carbon_budget * max( soil_C_min, fast_grnd_C_in       )
      fsc_ok_min  = - tol_carbon_budget * max( soil_C_min, fast_soil_C_in       )
      stgc_ok_min = - tol_carbon_budget * max( soil_C_min, structural_grnd_C_in )
      stsc_ok_min = - tol_carbon_budget * max( soil_C_min, structural_soil_C_in )
      msc_ok_min  = - tol_carbon_budget * max( soil_C_min, microbial_soil_C_in  )
      ssc_ok_min  = - tol_carbon_budget * max( soil_C_min, slow_soil_C_in       )
      psc_ok_min  = - tol_carbon_budget * max( soil_C_min, passive_soil_C_in    )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Check every soil carbon pool, to make sure they have sensible numbers.  Tiny    !
      ! negative stocks will be tolerated, but don't fix anything in case any of the pools !
      ! are too negative.                                                                  !
      !------------------------------------------------------------------------------------!
      neg_soil_C = csite%fast_grnd_C       (ipa) < fgc_ok_min   .or.                       &
                   csite%fast_soil_C       (ipa) < fsc_ok_min   .or.                       &
                   csite%structural_grnd_C (ipa) < stgc_ok_min  .or.                       &
                   csite%structural_soil_C (ipa) < stsc_ok_min  .or.                       &
                   csite%microbial_soil_C  (ipa) < msc_ok_min   .or.                       &
                   csite%slow_soil_C       (ipa) < ssc_ok_min   .or.                       &
                   csite%passive_soil_C    (ipa) < psc_ok_min
      !------------------------------------------------------------------------------------!


      if (.not. neg_soil_C) then
         !----- Account for any potential violation of carbon stocks. ---------------------!
         carbon_miss = - min(csite%fast_grnd_C       (ipa),0.0)                            &
                       - min(csite%fast_soil_C       (ipa),0.0)                            &
                       - min(csite%structural_grnd_C (ipa),0.0)                            &
                       - min(csite%structural_soil_C (ipa),0.0)                            &
                       - min(csite%microbial_soil_C  (ipa),0.0)                            &
                       - min(csite%slow_soil_C       (ipa),0.0)                            &
                       - min(csite%passive_soil_C    (ipa),0.0)
         !---------------------------------------------------------------------------------!



         !------ Force all pools to be either zero or positive. ---------------------------!
         csite%fast_grnd_C       (ipa) = max(0.0,csite%fast_grnd_C       (ipa))
         csite%fast_soil_C       (ipa) = max(0.0,csite%fast_soil_C       (ipa))
         csite%structural_grnd_C (ipa) = max(0.0,csite%structural_grnd_C (ipa))
         csite%structural_soil_C (ipa) = max(0.0,csite%structural_soil_C (ipa))
         csite%structural_grnd_L (ipa) = max(0.0,csite%structural_grnd_L (ipa))
         csite%structural_soil_L (ipa) = max(0.0,csite%structural_soil_L (ipa))
         csite%microbial_soil_C  (ipa) = max(0.0,csite%microbial_soil_C  (ipa))
         csite%slow_soil_C       (ipa) = max(0.0,csite%slow_soil_C       (ipa))
         csite%passive_soil_C    (ipa) = max(0.0,csite%passive_soil_C    (ipa))
         csite%fast_grnd_N       (ipa) = max(0.0,csite%fast_grnd_N       (ipa))
         csite%fast_soil_N       (ipa) = max(0.0,csite%fast_soil_N       (ipa))
         csite%mineralized_soil_N(ipa) = max(0.0,csite%mineralized_soil_N(ipa))
         !---------------------------------------------------------------------------------!
      else
         !----- Set missing carbon to zero so the code works with debugging. --------------!
         carbon_miss = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the net input (litter flux and mortality) and the total output             !
      ! (heterotrophic respiration).                                                       !
      !------------------------------------------------------------------------------------!
      net_C_input  = csite%fgc_in (ipa)  + csite%fsc_in (ipa)                              &
                   + csite%stgc_in(ipa)  + csite%stsc_in(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the initial and final soil carbon stocks.                                  !
      !------------------------------------------------------------------------------------!
      soil_C_initial = fast_grnd_C_in       + fast_soil_C_in       + structural_grnd_C_in  &
                     + structural_soil_C_in + microbial_soil_C_in  + slow_soil_C_in        &
                     + passive_soil_C_in
      soil_C_final   = csite%fast_grnd_C       (ipa) + csite%fast_soil_C       (ipa)       &
                     + csite%structural_grnd_C (ipa) + csite%structural_soil_C (ipa)       &
                     + csite%microbial_soil_C  (ipa) + csite%slow_soil_C       (ipa)       &
                     + csite%passive_soil_C    (ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Compare the heterotrophic respiration with two metrics: the residual from       !
      ! transition between pools, and the fraction of carbon loss that is supposed to be   !
      ! heterotrophic respiration.  These three should be all the same.                    !
      !------------------------------------------------------------------------------------!
      today_rh        = csite%today_rh(ipa) * umol_2_kgC * day_sec
      fracloss_rh     = er_fsc  * ( fg_C_loss  + fs_C_loss  )                              &
                      + er_stgc * stg_C_loss                                               &
                      + er_stsc * sts_C_loss                                               &
                      + er_msc  * ms_C_loss                                                &
                      + er_ssc  * ss_C_loss                                                &
                      + er_psc  * ps_C_loss
      nettrans_rh     = fg_C_loss + fs_C_loss + stg_C_loss + sts_C_loss + ms_C_loss        &
                      + ss_C_loss + ps_C_loss - ms_C_input - ss_C_input - ps_C_input
      net_C_loss      = onethird * (today_rh + fracloss_rh + nettrans_rh)
      toler_rh        = max(rh_min,net_C_loss) * tol_carbon_budget
      rh_violation(1) = abs(today_rh    - fracloss_rh) > toler_rh
      rh_violation(2) = abs(today_rh    - nettrans_rh) > toler_rh
      rh_violation(3) = abs(fracloss_rh - nettrans_rh) > toler_rh
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the change and the residual (unexplained change).                         !
      !------------------------------------------------------------------------------------!
      delta_soil_C     = net_C_input - net_C_loss
      resid_soil_C     = soil_C_final - soil_C_initial - delta_soil_C - carbon_miss
      toler_soil_C     = tol_carbon_budget * max(soil_C_initial,soil_C_min)
      soil_C_violation = abs(resid_soil_C) > toler_soil_C
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     In case we identify carbon conservation issues, print information on screen    !
      ! and stop the model.                                                                !
      !------------------------------------------------------------------------------------!
      if ( neg_soil_C .or. soil_C_violation .or. any(rh_violation)) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|       !!!   Soil Carbon budget failed   !!!        |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt=fmti )  ' PATCH               : ',ipa
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' FAST_GRND_C_IN        : ',fast_grnd_C_in
         write(unit=*,fmt=fmtf )  ' FAST_SOIL_C_IN        : ',fast_soil_C_in
         write(unit=*,fmt=fmtf )  ' STRUCTURAL_GRND_C_IN  : ',structural_grnd_C_in
         write(unit=*,fmt=fmtf )  ' STRUCTURAL_SOIL_C_IN  : ',structural_soil_C_in
         write(unit=*,fmt=fmtf )  ' MICROBIAL_SOIL_C_IN   : ',microbial_soil_C_in
         write(unit=*,fmt=fmtf )  ' SLOW_SOIL_C_IN        : ',slow_soil_C_in
         write(unit=*,fmt=fmtf )  ' PASSIVE_SOIL_C_IN     : ',passive_soil_C_in
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' FAST_GRND_C_FN        : ',csite%fast_grnd_C      (ipa)
         write(unit=*,fmt=fmtf )  ' FAST_SOIL_C_FN        : ',csite%fast_soil_C      (ipa)
         write(unit=*,fmt=fmtf )  ' STRUCTURAL_GRND_C_FN  : ',csite%structural_grnd_C(ipa)
         write(unit=*,fmt=fmtf )  ' STRUCTURAL_SOIL_C_FN  : ',csite%structural_soil_C(ipa)
         write(unit=*,fmt=fmtf )  ' MICROBIAL_SOIL_C_FN   : ',csite%microbial_soil_C (ipa)
         write(unit=*,fmt=fmtf )  ' SLOW_SOIL_C_FN        : ',csite%slow_soil_C      (ipa)
         write(unit=*,fmt=fmtf )  ' PASSIVE_SOIL_C_FN     : ',csite%passive_soil_C   (ipa)
         write(unit=*,fmt=fmtl )  ' NEGATIVE SOIL C       : ',neg_soil_C
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' FGC_INPUT             : ',csite%fgc_in(ipa)
         write(unit=*,fmt=fmtf )  ' FSC_INPUT             : ',csite%fsc_in(ipa)
         write(unit=*,fmt=fmtf )  ' STGC_INPUT            : ',csite%stgc_in(ipa)
         write(unit=*,fmt=fmtf )  ' STSC_INPUT            : ',csite%stsc_in(ipa)
         write(unit=*,fmt=fmtf )  ' MSC_INPUT             : ',ms_C_input
         write(unit=*,fmt=fmtf )  ' SSC_INPUT             : ',ss_C_input
         write(unit=*,fmt=fmtf )  ' PSC_INPUT             : ',ps_C_input
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' FGC_LOSS              : ',fg_C_loss
         write(unit=*,fmt=fmtf )  ' FSC_LOSS              : ',fs_C_loss
         write(unit=*,fmt=fmtf )  ' STGC_LOSS             : ',stg_C_loss
         write(unit=*,fmt=fmtf )  ' STSC_LOSS             : ',sts_C_loss
         write(unit=*,fmt=fmtf )  ' MSC_LOSS              : ',ms_C_loss
         write(unit=*,fmt=fmtf )  ' SSC_LOSS              : ',ss_C_loss
         write(unit=*,fmt=fmtf )  ' PSC_LOSS              : ',ps_C_loss
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' FRAC_FGC_CO2          : ',er_fsc
         write(unit=*,fmt=fmtf )  ' FRAC_FGC_MSC          : ',ex_fgc_msc
         write(unit=*,fmt=fmtf )  ' FRAC_FGC_SSC          : ',ex_fgc_ssc
         write(unit=*,fmt=fmtf )  ' FRAC_FGC_PSC          : ',ex_fgc_psc
         write(unit=*,fmt='(a)')  ' '
         write(unit=*,fmt=fmtf )  ' FRAC_FSC_CO2          : ',er_fsc
         write(unit=*,fmt=fmtf )  ' FRAC_FSC_MSC          : ',ex_fsc_msc
         write(unit=*,fmt=fmtf )  ' FRAC_FSC_SSC          : ',ex_fsc_ssc
         write(unit=*,fmt=fmtf )  ' FRAC_FSC_PSC          : ',ex_fsc_psc
         write(unit=*,fmt='(a)')  ' '
         write(unit=*,fmt=fmtf )  ' FRAC_STGC_CO2         : ',er_stgc
         write(unit=*,fmt=fmtf )  ' FRAC_STGC_MSC         : ',ex_stgc_msc
         write(unit=*,fmt=fmtf )  ' FRAC_STGC_SSC         : ',ex_stgc_ssc
         write(unit=*,fmt=fmtf )  ' FRAC_STGC_PSC         : ',ex_stgc_psc
         write(unit=*,fmt='(a)')  ' '
         write(unit=*,fmt=fmtf )  ' FRAC_STSC_CO2         : ',er_stsc
         write(unit=*,fmt=fmtf )  ' FRAC_STSC_MSC         : ',ex_stsc_msc
         write(unit=*,fmt=fmtf )  ' FRAC_STSC_SSC         : ',ex_stsc_ssc
         write(unit=*,fmt=fmtf )  ' FRAC_STSC_PSC         : ',ex_stsc_psc
         write(unit=*,fmt='(a)')  ' '
         write(unit=*,fmt=fmtf )  ' FRAC_MSC_CO2          : ',er_msc
         write(unit=*,fmt=fmtf )  ' FRAC_MSC_SSC          : ',ex_msc_ssc
         write(unit=*,fmt=fmtf )  ' FRAC_MSC_PSC          : ',ex_msc_psc
         write(unit=*,fmt='(a)')  ' '
         write(unit=*,fmt=fmtf )  ' FRAC_SSC_CO2          : ',er_ssc
         write(unit=*,fmt=fmtf )  ' FRAC_SSC_MSC          : ',ex_ssc_msc
         write(unit=*,fmt=fmtf )  ' FRAC_SSC_PSC          : ',ex_ssc_psc
         write(unit=*,fmt='(a)')  ' '
         write(unit=*,fmt=fmtf )  ' FRAC_PSC_CO2          : ',er_psc
         write(unit=*,fmt=fmtf )  ' FRAC_PSC_MSC          : ',ex_psc_msc
         write(unit=*,fmt=fmtf )  ' FRAC_PSC_SSC          : ',ex_psc_ssc
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' TODAY_RH (T)          : ',today_rh
         write(unit=*,fmt=fmtf )  ' FRACLOSS_RH (F)       : ',fracloss_rh
         write(unit=*,fmt=fmtf )  ' NETTRANS_RH (N)       : ',nettrans_rh
         write(unit=*,fmt=fmtf )  ' TOLERANCE_RH          : ',toler_rh
         write(unit=*,fmt=fmtl )  ' T-F VIOLATION         : ',rh_violation(1)
         write(unit=*,fmt=fmtl )  ' T-N VIOLATION         : ',rh_violation(2)
         write(unit=*,fmt=fmtl )  ' F-N VIOLATION         : ',rh_violation(3)
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' SOIL_C_INITIAL        : ',soil_C_initial
         write(unit=*,fmt=fmtf )  ' SOIL_C_FINAL          : ',soil_C_final
         write(unit=*,fmt=fmtf )  ' NET_C_INPUT           : ',net_C_input
         write(unit=*,fmt=fmtf )  ' NET_C_LOSS            : ',net_C_loss
         write(unit=*,fmt=fmtf )  ' DELTA_SOIL_C          : ',delta_soil_C
         write(unit=*,fmt=fmtf )  ' CARBON_MISS           : ',carbon_miss
         write(unit=*,fmt=fmtf )  ' RESIDUAL_SOIL_C       : ',resid_soil_C
         write(unit=*,fmt=fmtf )  ' TOLERANCE_SOIL_C      : ',toler_soil_C
         write(unit=*,fmt=fmtl )  ' SOIL CARBON VIOLATION : ',soil_C_violation
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  ' '


         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_budget_soilc','soil_respiration.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_budget_soilc
   !=======================================================================================!
   !=======================================================================================!
end module soil_respiration
!==========================================================================================!
!==========================================================================================!
