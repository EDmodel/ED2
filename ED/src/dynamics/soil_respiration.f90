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
      use decomp_coms  , only : k_rh_active              ! ! intent(in)
      use consts_coms  , only : wdns                     & ! intent(in)
                              , umols_2_kgCyr            ! ! intent(in)
      use therm_lib    , only : uextcm2tl                ! ! function
      use ed_misc_coms , only : dtlsm_o_frqsum           ! ! intent(in)
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
      real                                       :: Lc
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
                                      + cpatch%root_respiration(ico)
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




      !------------------------------------------------------------------------------------!
      !     Find the scaling factor for decomposition of surface materials, using the top  !
      ! soil layer to define their temperature and moisture.                               !
      !------------------------------------------------------------------------------------!
      k                   = mzg
      nsoil               = ntext_soil(k)
      rel_soil_moist      = min(1.0, max(0.0                                               &
                                        , (csite%soil_water(k,ipa) - soil(nsoil)%soilcp)   &
                                        / (soil(nsoil)%slmsts      - soil(nsoil)%soilcp) ) )
      csite%A_decomp(ipa) = het_resp_weight(csite%soil_tempk(k,ipa),rel_soil_moist)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Integrate the soil extensive properties, plus the minimum and maximum possible !
      ! soil water content of the active layer.                                            !
      !------------------------------------------------------------------------------------!
      sum_soil_energy = 0.0
      sum_soil_hcap   = 0.0
      sum_soil_water  = 0.0
      sum_soil_slmsts = 0.0
      sum_soil_soilcp = 0.0
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
      end do
      !------------------------------------------------------------------------------------!



      !----- Find the average temperature and the relative soil moisture. -----------------!
      call uextcm2tl(sum_soil_energy,sum_soil_water,sum_soil_hcap                          &
                    ,avg_soil_temp,avg_soil_fliq)
      rel_soil_moist = min( 1.0, max(0.0, ( sum_soil_water  - sum_soil_soilcp )            &
                                        / ( sum_soil_slmsts - sum_soil_soilcp ) ) )
      !------------------------------------------------------------------------------------!



      !----- Compute soil/temperature modulation of soil heterotrophic respiration. -------!
      csite%B_decomp(ipa) = het_resp_weight(avg_soil_temp,rel_soil_moist)
      !------------------------------------------------------------------------------------!



      !----- Compute nitrogen immobilization factor. --------------------------------------!
      nsoil = ntext_soil(mzg)
      call resp_f_decomp(csite,ipa, Lc,nsoil)
      !------------------------------------------------------------------------------------!



      !----- Compute heterotrophic respiration. -------------------------------------------!
      call resp_rh(csite,ipa, Lc,nsoil)
      !------------------------------------------------------------------------------------!



      !----- Update averaged variables. ---------------------------------------------------!
      csite%today_A_decomp (ipa) = csite%today_A_decomp (ipa) + csite%A_decomp(ipa)
      csite%today_B_decomp (ipa) = csite%today_B_decomp (ipa) + csite%B_decomp(ipa)
      csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa)                              &
                                 + csite%A_decomp       (ipa) * csite%f_decomp(ipa)
      csite%today_Bf_decomp(ipa) = csite%today_Bf_decomp(ipa)                              &
                                 + csite%B_decomp       (ipa) * csite%f_decomp(ipa)
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
      !------------------------------------------------------------------------------------!

      return
   end subroutine soil_respiration_driver
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the Nitrogen immobilization factor.                      !
   !---------------------------------------------------------------------------------------!
   subroutine resp_f_decomp(csite,ipa,Lc,ntext)

      use ed_state_vars, only : sitetype               ! ! structure
      use decomp_coms  , only : decomp_scheme          & ! intent(in)
                              , r_stsc                 & ! intent(in)
                              , N_immobil_supply_scale & ! intent(in)
                              , decay_rate_stsc        & ! intent(in)
                              , n_decomp_lim           & ! intent(in)
                              , xrothc_a               & ! intent(in)
                              , xrothc_b               & ! intent(in)
                              , xrothc_c               & ! intent(in)
                              , xrothc_d               & ! intent(in)
                              , c2n_structural         & ! intent(in)
                              , c2n_slow               ! ! intent(in)
      use soil_coms    , only : soil                   ! ! look-up table

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target      :: csite
      integer       , intent(in)  :: ipa
      real          , intent(out) :: Lc
      integer       , intent(in)  :: ntext
      !----- Local variables. -------------------------------------------------------------!
      real                        :: N_immobilization_demand
      real                        :: xrothc
      real                        :: er_stsc
      !------------------------------------------------------------------------------------!

    
      if (csite%structural_soil_C(ipa) > 0.0) then
         if (csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa)) then
            Lc = 0.049787 ! = exp(-3.0)
         else
            Lc = exp(-3.0 * csite%structural_soil_L(ipa)/csite%structural_soil_C(ipa))
         end if
      else
         Lc=0.0
      end if


      !------------------------------------------------------------------------------------!
      !      Decide the respiration factors based on the method.                           !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (2)
         xrothc  = xrothc_a * (xrothc_b + xrothc_c * exp( xrothc_d * soil(ntext)%xclay) )
         er_stsc = xrothc / (1.0 + xrothc)
      case default
         er_stsc = r_stsc
      end select
      !------------------------------------------------------------------------------------!



      if (n_decomp_lim == 1) then
         N_immobilization_demand = csite%A_decomp(ipa) * Lc * decay_rate_stsc              &
                                 * csite%structural_soil_C(ipa)                            &
                                 * ((1.0 - er_stsc) / c2n_slow - 1.0 / c2n_structural)
         
         csite%f_decomp(ipa)     = N_immobil_supply_scale * csite%mineralized_soil_N(ipa)  &
                                 / ( N_immobilization_demand                               &
                                   + N_immobil_supply_scale                                &
                                   * csite%mineralized_soil_N(ipa))
      else
         !----- Option for no plant N limitation. -----------------------------------------!
         csite%f_decomp(ipa)     = 1.0
      end if

      return
   end subroutine resp_f_decomp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the heterotrophic respiration.                           !
   !---------------------------------------------------------------------------------------!
   subroutine resp_rh(csite,ipa,Lc,ntext)

      use ed_state_vars, only : sitetype        ! ! structure
      use ed_misc_coms , only : current_time    ! ! intent(in)
      use consts_coms  , only : kgCday_2_umols  & ! intent(in)
                              , umols_2_kgCyr   ! ! intent(in)
      use soil_coms    , only : soil            ! ! look-up table
      use decomp_coms  , only : decomp_scheme   & ! intent(in)
                              , decay_rate_stsc & ! intent(in)
                              , decay_rate_msc  & ! intent(in)
                              , decay_rate_fsc  & ! intent(in)
                              , decay_rate_ssc  & ! intent(in)
                              , r_fsc           & ! intent(in)
                              , r_msc           & ! intent(in)
                              , r_stsc          & ! intent(in)
                              , xrothc_a        & ! intent(in)
                              , xrothc_b        & ! intent(in)
                              , xrothc_c        & ! intent(in)
                              , xrothc_d        ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)   , target     :: csite
      integer          , intent(in) :: ipa
      real             , intent(in) :: Lc
      integer          , intent(in) :: ntext
      !----- Local variables. -------------------------------------------------------------!
      real                          :: fgc_C_loss
      real                          :: fsc_C_loss
      real                          :: stgc_C_loss
      real                          :: stsc_C_loss
      real                          :: msc_C_loss
      real                          :: ssc_C_loss
      real                          :: tot_C_loss
      real                          :: er_fsc
      real                          :: er_ssc
      real                          :: er_msc
      real                          :: er_stsc
      real                          :: xrothc
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
            write (unit=84,fmt='(27(a,1x))')                                               &
                     '  YEAR',      ' MONTH',      '   DAY',      '  HOUR',      '   MIN'  &
              ,      '   IPA','       C_FGC','       C_FSC','      C_STGC','      C_STSC'  &
              ,'       C_MSC','       C_SSC','      F_FAST','    F_STRUCT','   F_MICROBE'  &
              ,'      F_SLOW','      RH_FGC','      RH_FSC','     RH_STGC','     RH_STSC'  &
              ,'      RH_MSC','      RH_SSC','    RH_TOTAL','    A_DECOMP','    B_DECOMP'  &
              ,'    F_DECOMP','          LC'
            close (unit=84,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!



      !----- The following variables have units of [umol_CO2/m2/s]. -----------------------!
      fgc_C_loss        = kgCday_2_umols  * csite%A_decomp(ipa)                            &
                        * decay_rate_fsc  * csite%fast_grnd_C(ipa)
      fsc_C_loss        = kgCday_2_umols  * csite%B_decomp(ipa)                            &
                        * decay_rate_fsc  * csite%fast_soil_C(ipa)
      stgc_C_loss       = kgCday_2_umols  * csite%A_decomp(ipa) * csite%f_decomp(ipa) * Lc &
                        * decay_rate_stsc * csite%structural_grnd_C(ipa)
      stsc_C_loss       = kgCday_2_umols  * csite%B_decomp(ipa) * csite%f_decomp(ipa) * Lc &
                        * decay_rate_stsc * csite%structural_soil_C(ipa)
      msc_C_loss        = kgCday_2_umols * csite%B_decomp(ipa)                             &
                        * decay_rate_msc * csite%microbial_soil_C(ipa)
      ssc_C_loss        = kgCday_2_umols * csite%B_decomp(ipa)                             &
                        * decay_rate_ssc * csite%slow_soil_C(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide the respiration factors based on the method.                           !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (2)
         xrothc  = xrothc_a * (xrothc_b + xrothc_c * exp( xrothc_d * soil(ntext)%xclay) )
         er_fsc  = xrothc / (1.0 + xrothc)
         er_msc  = er_fsc
         er_stsc = er_fsc
         er_ssc  = er_fsc
      case default
         er_fsc  = r_fsc
         er_msc  = r_msc
         er_stsc = r_stsc
         er_ssc  = 1.0    ! This has to be 1.0 because it is the only outlet for SSC.
      end select
      !------------------------------------------------------------------------------------!



      !----- Find the heterotrophic respiration (total and components). -------------------!
      csite%fgc_rh (ipa) = er_fsc  * fgc_C_loss
      csite%fsc_rh (ipa) = er_fsc  * fsc_C_loss
      csite%stgc_rh(ipa) = er_stsc * stgc_C_loss
      csite%stsc_rh(ipa) = er_stsc * stsc_C_loss
      csite%msc_rh (ipa) = er_msc  * msc_C_loss
      csite%ssc_rh (ipa) = er_ssc  * ssc_C_loss
      csite%rh     (ipa) = csite%fgc_rh (ipa) + csite%fsc_rh (ipa) + csite%stgc_rh(ipa)    &
                         + csite%stsc_rh(ipa) + csite%msc_rh (ipa) + csite%ssc_rh (ipa)
      !------------------------------------------------------------------------------------!



      !----- Write debugging information. -------------------------------------------------!
      if (print_debug) then
         !----- Convert units for output. -------------------------------------------------!
         fgc_C_loss  = umols_2_kgCyr * csite%fgc_rh (ipa)
         fsc_C_loss  = umols_2_kgCyr * csite%fsc_rh (ipa)
         stgc_C_loss = umols_2_kgCyr * csite%stgc_rh(ipa)
         stsc_C_loss = umols_2_kgCyr * csite%stsc_rh(ipa)
         msc_C_loss  = umols_2_kgCyr * csite%msc_rh (ipa)
         ssc_C_loss  = umols_2_kgCyr * csite%ssc_rh (ipa)
         tot_C_loss  = umols_2_kgCyr * csite%rh     (ipa)
         !---------------------------------------------------------------------------------!



         !----- Append step to the output file. -------------------------------------------!
         open (unit=84,file=rhetfile,status='old',position='append',action='write')
         write (unit=84,fmt='(6(i6,1x),21(f12.6,1x))')                                     &
                        current_time%year,current_time%month,current_time%date             &
                       ,current_time%hour,current_time%min,ipa                             &
                       ,csite%fast_grnd_C(ipa),csite%fast_soil_C(ipa)                      &
                       ,csite%structural_grnd_C(ipa),csite%structural_soil_C(ipa)          &
                       ,csite%microbial_soil_C(ipa),csite%slow_soil_C(ipa)                 &
                       ,er_fsc,er_stsc,er_msc,er_ssc,fgc_C_loss,fsc_C_loss,stgc_C_loss     &
                       ,stsc_C_loss,msc_C_loss,ssc_C_loss,tot_C_loss,csite%A_decomp(ipa)   &
                       ,csite%B_decomp(ipa),csite%f_decomp(ipa),Lc
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
                              , decay_rate_fsc  & ! intent(in)
                              , decay_rate_stsc & ! intent(in)
                              , decay_rate_msc  & ! intent(in)
                              , decay_rate_ssc  & ! intent(in)
                              , r_fsc           & ! intent(in)
                              , r_msc           & ! intent(in)
                              , r_stsc          & ! intent(in)
                              , xrothc_a        & ! intent(in)
                              , xrothc_b        & ! intent(in)
                              , xrothc_c        & ! intent(in)
                              , xrothc_d        & ! intent(in)
                              , fx_msc          & ! intent(in)
                              , c2n_slow        & ! intent(in)
                              , c2n_structural  ! ! intent(in)
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
      real                        :: Lc
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
      real                        :: er_stsc
      real                        :: ei_msc
      real                        :: ei_ssc
      real                        :: xrothc
      !------------------------------------------------------------------------------------!

      polygonloop: do ipy = 1,cgrid%npolygons

         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            
            !------------------------------------------------------------------------------!
            !      Define the scheme-dependent, heterotrophic respiration and transfer     !
            ! factors.                                                                     !
            !------------------------------------------------------------------------------!
            nsoil = cpoly%ntext_soil(nzg,isi)
            select case (decomp_scheme)
            case (2)
               !---------------------------------------------------------------------------!
               !   RothC model, following Sierra et al. (2012).  Loss of decayed necro-    !
               ! mass as CO2 depends on soil moisture.                                     !
               !---------------------------------------------------------------------------!
               xrothc  = xrothc_a                                                          &
                       * (xrothc_b + xrothc_c * exp( xrothc_d * soil(nsoil)%xclay) )
               er_stsc = xrothc / (1.0 + xrothc)
               ei_msc  = fx_msc / (1.0 + xrothc)
               ei_ssc  = (1.0 - fx_msc) / (1.0 + xrothc)
               !---------------------------------------------------------------------------!
            case default
               er_stsc = r_stsc
               !------ Not used. ----------------------------------------------------------!
               ei_msc  = 0.0
               ei_ssc  = 1.0
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!


            csite => cpoly%site(isi)

            patchloop: do ipa = 1,csite%npatches

               if (csite%structural_soil_C(ipa) > 0.0) then
                  if (csite%structural_soil_L(ipa) == csite%structural_soil_C(ipa)) then
                     Lc = 0.049787 ! = exp(-3.0)
                  else
                     Lc = exp( -3.0 * csite%structural_soil_L(ipa)                         &
                             /  csite%structural_soil_C(ipa))
                  end if
               else
                  Lc=0.0
               end if
         
               !----- Fast pools. ---------------------------------------------------------!
               fg_C_loss  = csite%today_A_decomp(ipa)         * decay_rate_fsc             &
                          * csite%fast_grnd_C(ipa)
               fs_C_loss  = csite%today_B_decomp(ipa)         * decay_rate_fsc             &
                          * csite%fast_soil_C(ipa)
               fg_N_loss  = csite%today_A_decomp(ipa)         * decay_rate_fsc             &
                          * csite%fast_grnd_N(ipa)
               fs_N_loss  = csite%today_B_decomp(ipa)         * decay_rate_fsc             &
                          * csite%fast_soil_N(ipa)

               !----- Structural pools. ---------------------------------------------------!
               stg_C_loss   = csite%today_Af_decomp(ipa)   * Lc * decay_rate_stsc          &
                            * csite%structural_grnd_C(ipa)
               sts_C_loss   = csite%today_Bf_decomp(ipa)   * Lc * decay_rate_stsc          &
                            * csite%structural_soil_C(ipa)
               stg_L_loss   = csite%today_Af_decomp(ipa)   * Lc * decay_rate_stsc          &
                            * csite%structural_grnd_L(ipa)
               sts_L_loss   = csite%today_Bf_decomp(ipa)   * Lc * decay_rate_stsc          &
                            * csite%structural_soil_L(ipa)
               stg_N_loss   = csite%today_Af_decomp(ipa)   * Lc * decay_rate_stsc          &
                            * csite%structural_grnd_N(ipa)
               sts_N_loss   = csite%today_Bf_decomp(ipa)   * Lc * decay_rate_stsc          &
                            * csite%structural_soil_N(ipa)
               !----- Demand to maintain stoichiometry (see Moorcroft et al. 2001). -------!
               stg_N_demand = stg_C_loss                                                   &
                            * ( (1.0 - er_stsc) * (1./c2n_slow) - (1./c2n_structural) )
               sts_N_demand = sts_C_loss                                                   &
                            * ( (1.0 - er_stsc) * (1./c2n_slow) - (1./c2n_structural) )
               !---------------------------------------------------------------------------!


               !----- Microbial pool (carbon only). ---------------------------------------!
               ms_C_loss  = csite%today_B_decomp(ipa)         * decay_rate_msc             &
                          * csite%microbial_soil_C(ipa)
               !---------------------------------------------------------------------------!


               !----- Slow pool (carbon only). --------------------------------------------!
               ss_C_loss  = csite%today_B_decomp(ipa)         * decay_rate_ssc             &
                          * csite%slow_soil_C(ipa)
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
               case (2)
                  !------------------------------------------------------------------------!
                  !     Microbial pool.  Note that part of the loss term is added back as  !
                  ! input, so we represent the (1-alpha_3,3) term in Sierra et al. (2012)  !
                  ! equation 21.                                                           !
                  !------------------------------------------------------------------------!
                  ms_C_input = ei_msc * ( fg_C_loss  + fs_C_loss  + stg_C_loss             &
                                        + sts_C_loss + ms_C_loss  + ss_C_loss  )
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Humified (slow) pool.  Note that part of the loss term is added    !
                  ! back as input, so we represent the (1-alpha_4,4) term in Sierra et al. !
                  ! (2012) equation 21.                                                    !
                  !------------------------------------------------------------------------!
                  ss_C_input = ei_ssc * ( fg_C_loss  + fs_C_loss  + stg_C_loss             &
                                        + sts_C_loss + ms_C_loss  + ss_C_loss  )
                  !------------------------------------------------------------------------!
               case default
                  !----- Microbial pool is inactive. --------------------------------------!
                  ms_C_input = 0.0
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Slow carbon receives everything that decayed and was not lost as   !
                  ! CO2 (heterotrophic respiration).                                       !
                  !------------------------------------------------------------------------!
                  ss_C_input = (1.0 - r_fsc ) * ( fg_C_loss  + fs_C_loss  )                &
                             + (1.0 - r_stsc) * ( stg_C_loss + sts_C_loss )                &
                             + (1.0 - r_msc ) * ms_C_loss
                  !------------------------------------------------------------------------!
               end select
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
               !      Force all pools to be either zero or positive.                       !
               !---------------------------------------------------------------------------!
               csite%fast_grnd_C       (ipa) = max(0.0,csite%fast_grnd_C       (ipa))
               csite%fast_soil_C       (ipa) = max(0.0,csite%fast_soil_C       (ipa))
               csite%structural_grnd_C (ipa) = max(0.0,csite%structural_grnd_C (ipa))
               csite%structural_soil_C (ipa) = max(0.0,csite%structural_soil_C (ipa))
               csite%structural_grnd_L (ipa) = max(0.0,csite%structural_grnd_L (ipa))
               csite%structural_soil_L (ipa) = max(0.0,csite%structural_soil_L (ipa))
               csite%microbial_soil_C  (ipa) = max(0.0,csite%microbial_soil_C  (ipa))
               csite%slow_soil_C       (ipa) = max(0.0,csite%slow_soil_C       (ipa))
               csite%fast_grnd_N       (ipa) = max(0.0,csite%fast_grnd_N       (ipa))
               csite%fast_soil_N       (ipa) = max(0.0,csite%fast_soil_N       (ipa))
               csite%mineralized_soil_N(ipa) = max(0.0,csite%mineralized_soil_N(ipa))
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
      real(kind=4)             :: sngloff
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
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: soil_tempk
      real(kind=4), intent(in) :: rel_soil_moist
      !----- Local variables. -------------------------------------------------------------!
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
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the temperature dependence.                                               !
      !------------------------------------------------------------------------------------!
      select case(decomp_scheme)
      case (0)
         !----- Use original ED-2.1 exponential temperature dependence. -------------------!
         temperature_limitation = min( 1.0                                                 &
                                     , exp(resp_temperature_increase * (soil_tempk-318.15)))
         !---------------------------------------------------------------------------------!
      case (1) 
         !----- Use Lloyd and Taylor (1994) temperature dependence. -----------------------!
         lnexplloyd             = rh_lloyd_1 * ( rh_lloyd_2 - 1./(soil_tempk - rh_lloyd_3))
         lnexplloyd             = max(lnexp_min,min(lnexp_max,lnexplloyd))
         temperature_limitation = min( 1.0, resp_temperature_increase * exp(lnexplloyd) )
         !---------------------------------------------------------------------------------!
      case (2)
         !---------------------------------------------------------------------------------!
         !      Similar to the original ED-1.0 formulation, which is based on the CENTURY  !
         ! model.  The change in the functional form is to avoid power of negative         !
         ! numbers, but the coefficients were tuned to give a similar curve.               !
         !---------------------------------------------------------------------------------!
         !----- Low temperature limitation. -----------------------------------------------!
         lnexplow               = rh_decay_low * (rh_low_temp - soil_tempk)
         lnexplow               = max(lnexp_min,min(lnexp_max,lnexplow))
         tlow_fun               = 1.0 + exp(lnexplow)
         !----- High temperature limitation. ----------------------------------------------!
         lnexphigh              = rh_decay_high * (soil_tempk - rh_high_temp)
         lnexphigh              = max(lnexp_min,min(lnexp_max,lnexphigh))
         thigh_fun              = 1.0 + exp(lnexphigh)
         !----- Temperature limitation is a combination of both. --------------------------!
         temperature_limitation = 1.0 / (tlow_fun * thigh_fun )
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the limitation due to soil moisture.                                      !
      !------------------------------------------------------------------------------------!
      select case (decomp_scheme)
      case (0,1)
         !----- ED-2.1 default, also used when decomp_scheme is 1. ------------------------!
         if (rel_soil_moist <= resp_opt_water)then
            water_limitation = exp((rel_soil_moist - resp_opt_water) * resp_water_below_opt)
         else
            water_limitation = exp((resp_opt_water - rel_soil_moist) * resp_water_above_opt)
         end if
         !---------------------------------------------------------------------------------!
      case (2)
         !----- Dry soil limitation. ------------------------------------------------------!
         lnexpdry         = rh_decay_dry * (rh_dry_smoist - rel_soil_moist)
         lnexpdry         = max(lnexp_min,min(lnexp_max,lnexpdry))
         smdry_fun        = 1.0 + exp(lnexpdry)
         !----- Wet soil limitation. ------------------------------------------------------!
         lnexpwet         = rh_decay_wet * (rel_soil_moist - rh_wet_smoist)
         lnexpwet         = max(lnexp_min,min(lnexp_max,lnexpwet))
         smwet_fun        = 1.0 + exp(lnexpwet)
         !----- Soil moisture limitation is a combination of both. ------------------------!
         water_limitation = 1.0 / (smdry_fun * smwet_fun)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Compute the weight, which is just the combination of both. -------------------!
      het_resp_weight = temperature_limitation * water_limitation
      !------------------------------------------------------------------------------------!

      return
   end function het_resp_weight
   !=======================================================================================!
   !=======================================================================================!
end module soil_respiration
!==========================================================================================!
!==========================================================================================!
