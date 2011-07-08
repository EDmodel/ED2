!==========================================================================================!
!==========================================================================================!
!                           |----------------------------------|                           !
!                           |** FREQUENT AVERAGE SUBROUTINES **|                           !
!                           |----------------------------------|                           !
!==========================================================================================!
!==========================================================================================!
!     This subroutine increments the time averaged polygon met-forcing variables.  These   !
! will be normalized by the output period to give time averages of each quantity.  The     !
! polygon level variables are derived from the weighted spatial average from the site      !
! level quantities.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine int_met_avg(cgrid)
   use ed_state_vars , only : edtype      & ! structure
                            , polygontype & ! structure
                            , sitetype    & ! structure
                            , patchtype   ! ! structure
   use ed_misc_coms  , only : dtlsm       & ! intent(in)
                            , frqsum      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy,isi,ipa,ico
   real                        :: frqsumi,tfact
   real                        :: polygon_area_i
   !---------------------------------------------------------------------------------------!

   !----- Some aliases. -------------------------------------------------------------------!
   frqsumi = 1.0 / frqsum
   tfact = dtlsm * frqsumi

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      polygon_area_i = 1. / sum(cpoly%area)

      do isi = 1,cpoly%nsites
         !---- Site-level averages.  ------------------------------------------------------!
         cpoly%avg_atm_tmp(isi)        = cpoly%avg_atm_tmp(isi)                            &
                                       + cpoly%met(isi)%atm_tmp  * tfact
         cpoly%avg_atm_shv(isi)        = cpoly%avg_atm_shv(isi)                            &
                                       + cpoly%met(isi)%atm_shv  * tfact
         cpoly%avg_atm_prss(isi)       = cpoly%avg_atm_prss(isi)                           &
                                       + cpoly%met(isi)%prss     * tfact

         !----- Now the polygon-level averages. -------------------------------------------!
         cgrid%avg_nir_beam(ipy)       = cgrid%avg_nir_beam(ipy)                           &
                                       + cpoly%met(isi)%nir_beam * cpoly%area(isi)         &
                                       * tfact * polygon_area_i

         cgrid%avg_nir_diffuse(ipy)    = cgrid%avg_nir_diffuse(ipy)                        &
                                       + cpoly%met(isi)%nir_diffuse * cpoly%area(isi)      &
                                       * tfact * polygon_area_i

         cgrid%avg_par_beam(ipy)       = cgrid%avg_par_beam(ipy)                           &
                                       + cpoly%met(isi)%par_beam * cpoly%area(isi)         &
                                       * tfact * polygon_area_i

         cgrid%avg_par_diffuse(ipy)    = cgrid%avg_par_diffuse(ipy)                        &
                                       + cpoly%met(isi)%par_diffuse * cpoly%area(isi)      &
                                       * tfact * polygon_area_i

         cgrid%avg_atm_tmp(ipy)        = cgrid%avg_atm_tmp(ipy)                            &
                                       + cpoly%met(isi)%atm_tmp * cpoly%area(isi)          &
                                       * tfact * polygon_area_i

         cgrid%avg_atm_shv(ipy)        = cgrid%avg_atm_shv(ipy)                            &
                                       + cpoly%met(isi)%atm_shv * cpoly%area(isi)          &
                                       * tfact * polygon_area_i

         cgrid%avg_rshort(ipy)         = cgrid%avg_rshort(ipy)                             &
                                       + cpoly%met(isi)%rshort * cpoly%area(isi)           &
                                       * tfact * polygon_area_i

         cgrid%avg_rshort_diffuse(ipy) = cgrid%avg_rshort_diffuse(ipy)                     &
                                       + cpoly%met(isi)%rshort_diffuse * cpoly%area(isi)   &
                                       * tfact * polygon_area_i

         cgrid%avg_rlong(ipy)          = cgrid%avg_rlong(ipy)                              &
                                       + cpoly%met(isi)%rlong * cpoly%area(isi)            &
                                       * tfact * polygon_area_i

         cgrid%avg_pcpg(ipy)           = cgrid%avg_pcpg(ipy)                               &
                                       + cpoly%met(isi)%pcpg * cpoly%area(isi)             &
                                       * tfact * polygon_area_i

         cgrid%avg_qpcpg(ipy)          = cgrid%avg_qpcpg(ipy)                              &
                                       + cpoly%met(isi)%qpcpg * cpoly%area(isi)            &
                                       * tfact * polygon_area_i

         cgrid%avg_dpcpg(ipy)          = cgrid%avg_dpcpg(ipy)                              &
                                       + cpoly%met(isi)%dpcpg * cpoly%area(isi)            &
                                       * tfact * polygon_area_i

         cgrid%avg_vels(ipy)           = cgrid%avg_vels(ipy)                               &
                                       + cpoly%met(isi)%vels * cpoly%area(isi)             &
                                       * tfact * polygon_area_i

         cgrid%avg_atm_prss(ipy)       = cgrid%avg_atm_prss(ipy)                           &
                                       + cpoly%met(isi)%prss * cpoly%area(isi)             &
                                       * tfact * polygon_area_i

         cgrid%avg_exner(ipy)          = cgrid%avg_exner(ipy)                              &
                                       + cpoly%met(isi)%exner * cpoly%area(isi)            &
                                       * tfact * polygon_area_i

         cgrid%avg_geoht(ipy)          = cgrid%avg_geoht(ipy)                              &
                                       + cpoly%met(isi)%geoht * cpoly%area(isi)            &
                                       * tfact * polygon_area_i

         cgrid%avg_atm_co2(ipy)        = cgrid%avg_atm_co2(ipy)                            &
                                       + cpoly%met(isi)%atm_co2 * cpoly%area(isi)          &
                                       * tfact * polygon_area_i

      end do
   end do
   return
end subroutine int_met_avg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!       The following sub-routine scales several variables that are integrated during one  !
! output step (frqsum) to actual rates.                                                    !
!------------------------------------------------------------------------------------------!
subroutine normalize_averaged_vars(cgrid,frqsum,dtlsm)
   use grid_coms    , only : nzg           ! ! intent(in)
   use ed_misc_coms , only : radfrq        & ! intent(in)
                           , current_time  ! ! intent(in)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   & ! structure
                           , sitetype      & ! structure
                           , patchtype     ! ! structure

   implicit none
   !----- Arguments.  ---------------------------------------------------------------------!
   type(edtype)      , target     :: cgrid
   real              , intent(in) :: frqsum
   real              , intent(in) :: dtlsm
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer    :: cpoly
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   real                           :: tfact
   real                           :: frqsumi
   integer                        :: k
   !---------------------------------------------------------------------------------------!
   

   !---------------------------------------------------------------------------------------!
   !     Find some useful conversion factors.                                              !
   ! 1. FRQSUMI = inverse of the elapsed time between two analyses (or one day).  This     !
   !              should be used by variables that are fluxes but are currently holding    !
   !              the integral over time.                                                  !
   ! 2. TFACT   = number of times steps since last analysis.  This should be used by       !
   !              variables that were added every time step, but not really integrated     !
   !              over time.                                                               !
   !---------------------------------------------------------------------------------------!
   frqsumi = 1.0 / frqsum
   tfact   = dtlsm * frqsumi
   !---------------------------------------------------------------------------------------!

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)

            !------------------------------------------------------------------------------!
            !     The following variables are integrated over time, we must divide them by !
            ! the total time.                                                              !
            !------------------------------------------------------------------------------!
            csite%avg_rshort_gnd     (ipa) = csite%avg_rshort_gnd     (ipa) * frqsumi
            csite%avg_rlong_gnd      (ipa) = csite%avg_rlong_gnd      (ipa) * frqsumi
            csite%avg_rlongup        (ipa) = csite%avg_rlongup        (ipa) * frqsumi
            csite%avg_albedo         (ipa) = csite%avg_albedo         (ipa) * frqsumi
            csite%avg_albedo_beam    (ipa) = csite%avg_albedo_beam    (ipa) * frqsumi
            csite%avg_albedo_diffuse (ipa) = csite%avg_albedo_diffuse (ipa) * frqsumi
            csite%avg_rlong_albedo   (ipa) = csite%avg_rlong_albedo   (ipa) * frqsumi
            csite%aux                (ipa) = csite%aux                (ipa) * frqsumi
            csite%avg_vapor_lc       (ipa) = csite%avg_vapor_lc       (ipa) * frqsumi
            csite%avg_vapor_wc       (ipa) = csite%avg_vapor_wc       (ipa) * frqsumi
            csite%avg_dew_cg         (ipa) = csite%avg_dew_cg         (ipa) * frqsumi
            csite%avg_vapor_gc       (ipa) = csite%avg_vapor_gc       (ipa) * frqsumi
            csite%avg_wshed_vg       (ipa) = csite%avg_wshed_vg       (ipa) * frqsumi
            csite%avg_intercepted    (ipa) = csite%avg_intercepted    (ipa) * frqsumi
            csite%avg_throughfall    (ipa) = csite%avg_throughfall    (ipa) * frqsumi
            csite%avg_vapor_ac       (ipa) = csite%avg_vapor_ac       (ipa) * frqsumi
            csite%avg_transp         (ipa) = csite%avg_transp         (ipa) * frqsumi
            csite%avg_evap           (ipa) = csite%avg_evap           (ipa) * frqsumi
            csite%avg_runoff         (ipa) = csite%avg_runoff         (ipa) * frqsumi
            csite%avg_drainage       (ipa) = csite%avg_drainage       (ipa) * frqsumi
            csite%avg_sensible_lc    (ipa) = csite%avg_sensible_lc    (ipa) * frqsumi
            csite%avg_sensible_wc    (ipa) = csite%avg_sensible_wc    (ipa) * frqsumi
            csite%avg_qwshed_vg      (ipa) = csite%avg_qwshed_vg      (ipa) * frqsumi
            csite%avg_qintercepted   (ipa) = csite%avg_qintercepted   (ipa) * frqsumi
            csite%avg_qthroughfall   (ipa) = csite%avg_qthroughfall   (ipa) * frqsumi
            csite%avg_sensible_gc    (ipa) = csite%avg_sensible_gc    (ipa) * frqsumi
            csite%avg_sensible_ac    (ipa) = csite%avg_sensible_ac    (ipa) * frqsumi
            csite%avg_carbon_ac      (ipa) = csite%avg_carbon_ac      (ipa) * frqsumi
            csite%avg_runoff_heat    (ipa) = csite%avg_runoff_heat    (ipa) * frqsumi
            csite%avg_drainage_heat  (ipa) = csite%avg_drainage_heat  (ipa) * frqsumi
            csite%avg_rk4step        (ipa) = csite%avg_rk4step        (ipa) * frqsumi

         
            do k=cpoly%lsl(isi),nzg
               csite%avg_sensible_gg(k,ipa) = csite%avg_sensible_gg(k,ipa) * frqsumi
               csite%avg_smoist_gg(k,ipa)   = csite%avg_smoist_gg(k,ipa)   * frqsumi
               csite%avg_transloss(k,ipa)   = csite%avg_transloss(k,ipa)   * frqsumi
               csite%aux_s(k,ipa)           = csite%aux_s(k,ipa)           * frqsumi
            end do
            !------------------------------------------------------------------------------!
            
            !------------------------------------------------------------------------------!
            !     Available water is added every dtlsm, so we normalise using tfact.       !
            !------------------------------------------------------------------------------!
            csite%avg_available_water(ipa) = csite%avg_available_water(ipa)       * tfact
            
            do ico=1,cpatch%ncohorts
               !---------------------------------------------------------------------------!
               !      The carbon fluxes were updated every time step, but they don't have  !
               ! integral units, so we normalise using tfact.  Their units will become     !
               ! µmol/m²/s.                                                                !
               !---------------------------------------------------------------------------!
               cpatch%mean_leaf_resp(ico)    = cpatch%mean_leaf_resp(ico)    * tfact
               cpatch%mean_root_resp(ico)    = cpatch%mean_root_resp(ico)    * tfact
               cpatch%mean_gpp(ico)          = cpatch%mean_gpp(ico)          * tfact
               cpatch%mean_storage_resp(ico) = cpatch%mean_storage_resp(ico) * tfact
               cpatch%mean_growth_resp(ico)  = cpatch%mean_growth_resp(ico)  * tfact
               cpatch%mean_vleaf_resp(ico)   = cpatch%mean_vleaf_resp(ico)   * tfact
            end do

            !------------------------------------------------------------------------------!
            !      Likewise, heterotrophic respiration was updated every time step, but it !
            !  doesn't have integral values, so we normalise using tfact.  Its units will  !
            ! become µmol/m²/s.                                                            !
            !------------------------------------------------------------------------------!
            csite%mean_rh(ipa)               = csite%mean_rh(ipa)            * tfact
            
            !------------------------------------------------------------------------------!
            !      Budget variables.  They contain integral values, so we must divide by   !
            ! the elapsed time to get them in flux units.                                  !
            !------------------------------------------------------------------------------!
            csite%co2budget_gpp(ipa)         = csite%co2budget_gpp(ipa)         * frqsumi
            csite%co2budget_gpp_dbh(:,ipa)   = csite%co2budget_gpp_dbh(:,ipa)   * frqsumi
            csite%co2budget_plresp(ipa)      = csite%co2budget_plresp(ipa)      * frqsumi
            csite%co2budget_rh(ipa)          = csite%co2budget_rh(ipa)          * frqsumi
            csite%co2budget_loss2atm(ipa)    = csite%co2budget_loss2atm(ipa)    * frqsumi
            csite%co2budget_denseffect(ipa)  = csite%co2budget_denseffect(ipa)  * frqsumi
            csite%co2budget_residual(ipa)    = csite%co2budget_residual(ipa)    * frqsumi
            csite%ebudget_precipgain(ipa)    = csite%ebudget_precipgain(ipa)    * frqsumi
            csite%ebudget_netrad(ipa)        = csite%ebudget_netrad(ipa)        * frqsumi
            csite%ebudget_denseffect(ipa)    = csite%ebudget_denseffect(ipa)    * frqsumi
            csite%ebudget_loss2atm(ipa)      = csite%ebudget_loss2atm(ipa)      * frqsumi
            csite%ebudget_loss2drainage(ipa) = csite%ebudget_loss2drainage(ipa) * frqsumi
            csite%ebudget_loss2runoff(ipa)   = csite%ebudget_loss2runoff(ipa)   * frqsumi
            csite%ebudget_residual(ipa)      = csite%ebudget_residual(ipa)      * frqsumi
            csite%wbudget_precipgain(ipa)    = csite%wbudget_precipgain(ipa)    * frqsumi
            csite%wbudget_loss2atm(ipa)      = csite%wbudget_loss2atm(ipa)      * frqsumi
            csite%wbudget_loss2drainage(ipa) = csite%wbudget_loss2drainage(ipa) * frqsumi
            csite%wbudget_loss2runoff(ipa)   = csite%wbudget_loss2runoff(ipa)   * frqsumi
            csite%wbudget_denseffect(ipa)    = csite%wbudget_denseffect(ipa)    * frqsumi
            csite%wbudget_residual(ipa)      = csite%wbudget_residual(ipa)      * frqsumi
            !------------------------------------------------------------------------------!
         end do
      end do
   end do
   
   return
end subroutine normalize_averaged_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine reset_averaged_vars(cgrid)

   use ed_state_vars, only : edtype      & ! structure
                           , polygontype & ! structure
                           , sitetype    & ! structure
                           , patchtype   ! ! structure
   
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy
   integer                    :: isi
   integer                    :: ipa
   integer                    :: ico
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy = 1,cgrid%npolygons

      cgrid%cbudget_nep          (ipy) = 0.0
      cgrid%avg_nir_beam         (ipy) = 0.0
      cgrid%avg_nir_diffuse      (ipy) = 0.0
      cgrid%avg_par_beam         (ipy) = 0.0
      cgrid%avg_par_diffuse      (ipy) = 0.0
      cgrid%avg_atm_tmp          (ipy) = 0.0
      cgrid%avg_atm_shv          (ipy) = 0.0
      cgrid%avg_rshort           (ipy) = 0.0
      cgrid%avg_rshort_diffuse   (ipy) = 0.0
      cgrid%avg_rlong            (ipy) = 0.0
      cgrid%avg_pcpg             (ipy) = 0.0
      cgrid%avg_qpcpg            (ipy) = 0.0
      cgrid%avg_dpcpg            (ipy) = 0.0
      cgrid%avg_vels             (ipy) = 0.0
      cgrid%avg_atm_prss         (ipy) = 0.0
      cgrid%avg_exner            (ipy) = 0.0
      cgrid%avg_geoht            (ipy) = 0.0
      cgrid%avg_atm_co2          (ipy) = 0.0
      cgrid%avg_albedo           (ipy) = 0.0
      cgrid%avg_albedo_beam      (ipy) = 0.0
      cgrid%avg_albedo_diffuse   (ipy) = 0.0
      cgrid%avg_rlong_albedo     (ipy) = 0.0
      cgrid%avg_rlongup          (ipy) = 0.0

      cgrid%avg_leaf_energy      (ipy) = 0.0
      cgrid%avg_leaf_temp        (ipy) = 0.0
      cgrid%avg_leaf_hcap        (ipy) = 0.0
      cgrid%avg_leaf_fliq        (ipy) = 0.0
      cgrid%avg_leaf_water       (ipy) = 0.0

      cgrid%avg_wood_energy      (ipy) = 0.0
      cgrid%avg_wood_temp        (ipy) = 0.0
      cgrid%avg_wood_hcap        (ipy) = 0.0
      cgrid%avg_wood_fliq        (ipy) = 0.0
      cgrid%avg_wood_water       (ipy) = 0.0

      cgrid%avg_can_temp         (ipy) = 0.0
      cgrid%avg_can_shv          (ipy) = 0.0
      cgrid%avg_can_co2          (ipy) = 0.0
      cgrid%avg_can_rhos         (ipy) = 0.0
      cgrid%avg_can_prss         (ipy) = 0.0
      cgrid%avg_can_theta        (ipy) = 0.0
      cgrid%avg_can_theiv        (ipy) = 0.0
      cgrid%avg_can_depth        (ipy) = 0.0

      cgrid%avg_drainage         (ipy) = 0.0
      cgrid%avg_evap             (ipy) = 0.0
      cgrid%avg_transp           (ipy) = 0.0
      cgrid%avg_soil_temp      (:,ipy) = 0.0
      cgrid%avg_soil_water     (:,ipy) = 0.0
      cgrid%avg_soil_mstpot    (:,ipy) = 0.0
      cgrid%avg_soil_energy    (:,ipy) = 0.0
      cgrid%avg_soil_fracliq   (:,ipy) = 0.0
      cgrid%avg_soil_rootfrac  (:,ipy) = 0.0

      cgrid%avg_vapor_lc         (ipy) = 0.0
      cgrid%avg_vapor_wc         (ipy) = 0.0
      cgrid%avg_dew_cg           (ipy) = 0.0
      cgrid%avg_vapor_gc         (ipy) = 0.0
      cgrid%avg_wshed_vg         (ipy) = 0.0
      cgrid%avg_intercepted      (ipy) = 0.0
      cgrid%avg_throughfall      (ipy) = 0.0
      cgrid%avg_vapor_ac         (ipy) = 0.0
      cgrid%avg_transp           (ipy) = 0.0
      cgrid%avg_evap             (ipy) = 0.0
      cgrid%avg_runoff           (ipy) = 0.0
      cgrid%avg_drainage         (ipy) = 0.0
      cgrid%avg_drainage_heat    (ipy) = 0.0
      cgrid%aux                  (ipy) = 0.0
      cgrid%avg_carbon_ac        (ipy) = 0.0
      cgrid%avg_sensible_lc      (ipy) = 0.0
      cgrid%avg_sensible_wc      (ipy) = 0.0
      cgrid%avg_qwshed_vg        (ipy) = 0.0
      cgrid%avg_qintercepted     (ipy) = 0.0
      cgrid%avg_qthroughfall     (ipy) = 0.0
      cgrid%avg_sensible_gc      (ipy) = 0.0
      cgrid%avg_sensible_ac      (ipy) = 0.0
      cgrid%avg_runoff_heat      (ipy) = 0.0 

      cgrid%aux_s              (:,ipy) = 0.0
      cgrid%avg_smoist_gg      (:,ipy) = 0.0
      cgrid%avg_transloss      (:,ipy) = 0.0
      cgrid%avg_sensible_gg    (:,ipy) = 0.0

      cgrid%avg_soil_wetness     (ipy) = 0.0
      cgrid%avg_skin_temp        (ipy) = 0.0
      cgrid%avg_available_water  (ipy) = 0.0

      cgrid%avg_lai_ebalvars (:,:,ipy) = 0.0

      cgrid%avg_gpp              (ipy) = 0.0
      cgrid%avg_nppleaf          (ipy) = 0.0
      cgrid%avg_nppfroot         (ipy) = 0.0
      cgrid%avg_nppsapwood       (ipy) = 0.0
      cgrid%avg_nppcroot         (ipy) = 0.0
      cgrid%avg_nppseeds         (ipy) = 0.0
      cgrid%avg_nppwood          (ipy) = 0.0
      cgrid%avg_nppdaily         (ipy) = 0.0
      cgrid%avg_leaf_resp        (ipy) = 0.0
      cgrid%avg_root_resp        (ipy) = 0.0
      cgrid%avg_growth_resp      (ipy) = 0.0
      cgrid%avg_storage_resp     (ipy) = 0.0
      cgrid%avg_vleaf_resp       (ipy) = 0.0
      cgrid%avg_plant_resp       (ipy) = 0.0
      cgrid%avg_growth_resp      (ipy) = 0.0
      cgrid%avg_storage_resp     (ipy) = 0.0
      cgrid%avg_vleaf_resp       (ipy) = 0.0
      cgrid%avg_htroph_resp      (ipy) = 0.0 
      cgrid%avg_leaf_drop        (ipy) = 0.0
      cgrid%avg_leaf_maintenance (ipy) = 0.0
      cgrid%avg_root_maintenance (ipy) = 0.0 

      cgrid%avg_sfcw_depth       (ipy) = 0.0
      cgrid%avg_sfcw_energy      (ipy) = 0.0
      cgrid%avg_sfcw_mass        (ipy) = 0.0
      cgrid%avg_sfcw_tempk       (ipy) = 0.0
      cgrid%avg_sfcw_fracliq     (ipy) = 0.0
      cgrid%avg_bdead            (ipy) = 0.0
      cgrid%avg_balive           (ipy) = 0.0
      cgrid%avg_bleaf            (ipy) = 0.0
      cgrid%avg_broot            (ipy) = 0.0
      cgrid%avg_bsapwood         (ipy) = 0.0
      cgrid%avg_bstorage         (ipy) = 0.0 
      cgrid%avg_bseeds           (ipy) = 0.0
      cgrid%avg_fsc              (ipy) = 0.0
      cgrid%avg_ssc              (ipy) = 0.0 
      cgrid%avg_stsc             (ipy) = 0.0
      cgrid%avg_fsn              (ipy) = 0.0
      cgrid%avg_msn              (ipy) = 0.0 


      cpoly => cgrid%polygon(ipy)
      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         cpoly%avg_atm_tmp(isi)          = 0.0
         cpoly%avg_atm_shv(isi)          = 0.0
         cpoly%avg_atm_prss(isi)         = 0.0

         cpoly%avg_soil_temp(:,isi)      = 0.0
         cpoly%avg_soil_water(:,isi)     = 0.0
         cpoly%avg_soil_mstpot(:,isi)    = 0.0
         cpoly%avg_soil_energy(:,isi)    = 0.0
         cpoly%avg_soil_fracliq(:,isi)   = 0.0
         cpoly%avg_soil_rootfrac(:,isi)  = 0.0



         patchloop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)

            !----------------------------------------------------------------!
            ! Zeroing CO2 budget variables.                                  !
            !----------------------------------------------------------------!
            csite%co2budget_gpp(ipa)            = 0.0
            csite%co2budget_gpp_dbh(:,ipa)      = 0.0
            csite%co2budget_rh(ipa)             = 0.0
            csite%co2budget_plresp(ipa)         = 0.0
            csite%co2budget_residual(ipa)       = 0.0
            csite%co2budget_loss2atm(ipa)       = 0.0
            csite%co2budget_denseffect(ipa)     = 0.0

            !----------------------------------------------------------------!
            ! Zeroing water budget variables.                                !
            !----------------------------------------------------------------!
            csite%wbudget_precipgain(ipa)       = 0.0
            csite%wbudget_loss2atm(ipa)         = 0.0
            csite%wbudget_loss2runoff(ipa)      = 0.0
            csite%wbudget_loss2drainage(ipa)    = 0.0
            csite%wbudget_denseffect(ipa)       = 0.0
            csite%wbudget_residual(ipa)         = 0.0


            !----------------------------------------------------------------!
            ! Zeroing energy budget variables.                               !
            !----------------------------------------------------------------!
            csite%ebudget_precipgain(ipa)       = 0.0
            csite%ebudget_netrad(ipa)           = 0.0
            csite%ebudget_loss2atm(ipa)         = 0.0
            csite%ebudget_loss2runoff(ipa)      = 0.0
            csite%ebudget_loss2drainage(ipa)    = 0.0
            csite%ebudget_denseffect(ipa)       = 0.0
            csite%ebudget_residual(ipa)         = 0.0
            !----------------------------------------------------------------!

            csite%avg_carbon_ac(ipa)        = 0.0
            csite%avg_vapor_lc(ipa)         = 0.0
            csite%avg_vapor_wc(ipa)         = 0.0
            csite%avg_dew_cg(ipa)           = 0.0
            csite%avg_vapor_gc(ipa)         = 0.0
            csite%avg_wshed_vg(ipa)         = 0.0
            csite%avg_intercepted(ipa)      = 0.0
            csite%avg_throughfall(ipa)      = 0.0
            csite%avg_vapor_ac(ipa)         = 0.0
            csite%avg_transp(ipa)           = 0.0
            csite%avg_evap(ipa)             = 0.0
            csite%avg_rshort_gnd(ipa)       = 0.0
            csite%avg_rlong_gnd(ipa)        = 0.0
            csite%avg_rlongup         (ipa) = 0.0
            csite%avg_albedo          (ipa) = 0.0
            csite%avg_albedo_beam     (ipa) = 0.0
            csite%avg_albedo_diffuse  (ipa) = 0.0
            csite%avg_rlong_albedo    (ipa) = 0.0
            csite%avg_smoist_gg(:,ipa)      = 0.0
            csite%avg_transloss(:,ipa)      = 0.0
            csite%avg_runoff(ipa)           = 0.0
            csite%avg_runoff_heat(ipa)      = 0.0
            csite%avg_drainage(ipa)         = 0.0
            csite%avg_drainage_heat(ipa)    = 0.0
            csite%avg_sensible_lc(ipa)      = 0.0
            csite%avg_sensible_wc(ipa)      = 0.0
            csite%avg_qwshed_vg(ipa)        = 0.0
            csite%avg_qintercepted(ipa)     = 0.0
            csite%avg_qthroughfall(ipa)     = 0.0
            csite%avg_sensible_gc(ipa)      = 0.0
            csite%avg_sensible_ac(ipa)      = 0.0
            csite%avg_sensible_gg(:,ipa)    = 0.0
            csite%avg_runoff_heat(ipa)      = 0.0
            csite%avg_rk4step(ipa)          = 0.0
            csite%avg_available_water(ipa)  = 0.0
            csite%aux(ipa)                  = 0.0
            csite%aux_s(:,ipa)              = 0.0
            csite%mean_rh(ipa)              = 0.0
         
            cohortloop: do ico=1,cpatch%ncohorts
               cpatch%leaf_respiration(ico)      = 0.0
               cpatch%root_respiration(ico)      = 0.0
               cpatch%gpp(ico)                   = 0.0
               cpatch%mean_leaf_resp(ico)        = 0.0
               cpatch%mean_root_resp(ico)        = 0.0
               cpatch%mean_gpp(ico)              = 0.0
               cpatch%mean_storage_resp(ico)     = 0.0 
               cpatch%mean_growth_resp(ico)      = 0.0
               cpatch%mean_vleaf_resp(ico)       = 0.0
            end do cohortloop
         end do patchloop
      end do siteloop
   end do polyloop


   return
end subroutine reset_averaged_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!                             |-------------------------------|                            !
!                             |** DAILY AVERAGE SUBROUTINES **|                            !
!                             |-------------------------------|                            !
!==========================================================================================!
!==========================================================================================!
!    This subroutine integrates the daily averages for state vars, plus the mean diurnal   !
! cycle if the user wants that as well.  This is called after each time step in case at    !
! least one of the daily mean, monthly mean, or the mean diurnal is sought.                !
!------------------------------------------------------------------------------------------!
subroutine integrate_ed_daily_output_state(cgrid)
   use ed_state_vars        , only : edtype              & ! structure
                                   , polygontype         & ! structure
                                   , sitetype            & ! structure
                                   , patchtype           ! ! structure
   use grid_coms            , only : nzg                 & ! intent(in)
                                   , time                ! ! intent(in)
   use ed_max_dims          , only : n_dbh               & ! intent(in)
                                   , n_pft               & ! intent(in)
                                   , n_dist_types        & ! intent(in)
                                   , n_mort              ! ! intent(in)
   use ed_misc_coms         , only : dtlsm               & ! intent(in)
                                   , iqoutput            & ! intent(in)
                                   , frqfast             & ! intent(in)
                                   , ndcycle             & ! intent(in)
                                   , current_time        ! ! intent(in)
   use pft_coms             , only : sla                 ! ! intent(in)
   use canopy_radiation_coms, only : rshort_twilight_min & ! intent(in)
                                   , cosz_min            ! ! intent(in)
   use consts_coms          , only : day_sec             ! ! intent(in)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   integer                     :: ico
   integer                     :: it
   real                        :: poly_area_i
   real                        :: site_area_i
   real                        :: forest_site
   real                        :: forest_site_i
   real                        :: forest_poly
   real                        :: poly_lai
   real                        :: site_lai
   real                        :: patch_lai
   real                        :: patch_lai_i
   real                        :: poly_lma
   real                        :: site_lma
   real                        :: patch_lma
   real                        :: sss_can_theta
   real                        :: sss_can_theiv
   real                        :: sss_can_shv
   real                        :: sss_can_co2
   real                        :: sss_can_prss
   real                        :: sss_gnd_temp
   real                        :: sss_gnd_shv
   real                        :: sss_leaf_water
   real                        :: sss_leaf_energy
   real                        :: sss_leaf_hcap
   real                        :: sss_wood_water
   real                        :: sss_wood_energy
   real                        :: sss_wood_hcap
   real                        :: pss_leaf_water
   real                        :: pss_leaf_energy
   real                        :: pss_leaf_hcap
   real                        :: pss_wood_water
   real                        :: pss_wood_energy
   real                        :: pss_wood_hcap
   real                        :: rshort_tot
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the index corresponding to this time of the day for the mean diurnal cycle   !
   ! averages.                                                                             !
   !---------------------------------------------------------------------------------------!
   it = ceiling(mod(current_time%time,day_sec)/frqfast)
   if (it == 0) it = nint(day_sec/frqfast)
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      !----- Initialize auxiliary variables to add sitetype variables. --------------------!
      site_lma         = 0.
      site_lai         = 0.
      sss_leaf_energy   = 0.
      sss_leaf_water    = 0.
      sss_leaf_hcap     = 0.
      sss_wood_energy   = 0.
      sss_wood_water    = 0.
      sss_wood_hcap     = 0.
      sss_can_theta    = 0.
      sss_can_theiv    = 0.
      sss_can_shv      = 0.
      sss_can_co2      = 0.
      sss_can_prss     = 0.
      sss_gnd_temp     = 0.
      sss_gnd_shv      = 0.
      forest_poly      = 0.
      !------------------------------------------------------------------------------------!

      siteloop: do isi=1, cpoly%nsites
         !---------------------------------------------------------------------------------!
         !     Determine if there was any beam radiation, and compute the total (I think   !
         ! rshort is always the total, but not sure.  This is the same test done in the    !
         ! radiation driver.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpoly%cosaoi(isi) <= cosz_min) then
            rshort_tot = cpoly%met(isi)%rshort_diffuse
         else
            rshort_tot = cpoly%met(isi)%rshort
         end if
         !----- Include this time step only if it is still day time. ----------------------!
         if (rshort_tot > rshort_twilight_min) then
            cpoly%daylight(isi) = cpoly%daylight(isi) + dtlsm
         end if
         !---------------------------------------------------------------------------------!

         csite => cpoly%site(isi)
         
         !----- Inverse of total site area (sum of all patches' area). --------------------!
         site_area_i=1./sum(csite%area)

         !----- Forest areas. -------------------------------------------------------------!
         forest_site           = sum(csite%area,csite%dist_type /= 1)
         if (forest_site > 1.0e-6) then
            forest_site_i         = 1./forest_site
         else
            forest_site_i         = 0.0
         end if
         forest_poly           = forest_poly + forest_site

         !----- Initialize auxiliary variables to add patchtype variables. ----------------!
         patch_lai         = 0.
         patch_lma         = 0.
         pss_leaf_energy   = 0.
         pss_leaf_water    = 0.
         pss_leaf_hcap     = 0.
         pss_wood_energy   = 0.
         pss_wood_water    = 0.
         pss_wood_hcap     = 0.

         !----- Loop over the patches to normalize the sum of all cohorts. ----------------!
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)

            if (cpatch%ncohorts > 0) then
               pss_leaf_energy = pss_leaf_energy + sum(cpatch%leaf_energy)*csite%area(ipa)
               pss_leaf_water  = pss_leaf_water  + sum(cpatch%leaf_water )*csite%area(ipa)
               pss_leaf_hcap   = pss_leaf_hcap   + sum(cpatch%leaf_hcap  )*csite%area(ipa)
               pss_wood_energy = pss_wood_energy + sum(cpatch%wood_energy)*csite%area(ipa)
               pss_wood_water  = pss_wood_water  + sum(cpatch%wood_water )*csite%area(ipa)
               pss_wood_hcap   = pss_wood_hcap   + sum(cpatch%wood_hcap  )*csite%area(ipa)
               patch_lai_i = 1./max(tiny(1.),sum(cpatch%lai,cpatch%leaf_resolvable))
               patch_lai = patch_lai                                                       &
                         + csite%area(ipa) * sum(cpatch%lai,cpatch%leaf_resolvable)
               patch_lma = patch_lma                                                       &
                         + csite%area(ipa) * sum(cpatch%lai/sla(cpatch%pft)                &
                                                ,cpatch%leaf_resolvable)
            end if

            do ico=1,cpatch%ncohorts
               if (cpatch%leaf_resolvable(ico)) then
                  cpatch%dmean_par_l       (ico) = cpatch%dmean_par_l       (ico)          &
                                                 + cpatch%par_l             (ico)
                  cpatch%dmean_par_l_beam  (ico) = cpatch%dmean_par_l_beam  (ico)          &
                                                 + cpatch%par_l_beam        (ico)
                  cpatch%dmean_par_l_diff  (ico) = cpatch%dmean_par_l_diff  (ico)          &
                                                 + cpatch%par_l_diffuse     (ico)

                  !------------------------------------------------------------------------!
                  !    Integrate the photosynthesis-related variables and the light level  !
                  ! only if it is day time.                                                !
                  !------------------------------------------------------------------------!
                  if (rshort_tot > rshort_twilight_min) then
                     cpatch%dmean_fs_open     (ico) = cpatch%dmean_fs_open     (ico)       &
                                                    + cpatch%fs_open           (ico)
                     cpatch%dmean_fsw         (ico) = cpatch%dmean_fsw         (ico)       &
                                                    + cpatch%fsw               (ico)
                     cpatch%dmean_fsn         (ico) = cpatch%dmean_fsn         (ico)       &
                                                    + cpatch%fsn               (ico)
                     cpatch%dmean_psi_open    (ico) = cpatch%dmean_psi_open    (ico)       &
                                                    + cpatch%psi_open          (ico)
                     cpatch%dmean_psi_closed  (ico) = cpatch%dmean_psi_closed  (ico)       &
                                                    + cpatch%psi_closed        (ico)
                     cpatch%dmean_water_supply(ico) = cpatch%dmean_water_supply(ico)       &
                                                    + cpatch%water_supply      (ico)
                     cpatch%dmean_light_level(ico)  = cpatch%dmean_light_level (ico)       &
                                                    + cpatch%light_level       (ico)
                     cpatch%dmean_light_level_beam(ico) =                                  &
                                               cpatch%dmean_light_level_beam(ico)          &
                                             + cpatch%light_level_beam(ico)
                     cpatch%dmean_light_level_diff(ico) =                                  &
                                               cpatch%dmean_light_level_diff(ico)          &
                                             + cpatch%light_level_diff(ico)
                     cpatch%dmean_beamext_level(ico)    = cpatch%dmean_beamext_level(ico)  &
                                                        + cpatch%beamext_level(ico)
                     cpatch%dmean_diffext_level(ico)    = cpatch%dmean_diffext_level(ico)  &
                                                        + cpatch%diffext_level(ico)
                     cpatch%dmean_lambda_light(ico)     = cpatch%dmean_lambda_light(ico)   &
                                                        + cpatch%lambda_light(ico)
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Integrate the mean diurnal cycle only if the mean diurnal cycle    !
                  ! output is sought.  This is also done for night time as the length of   !
                  ! the day may vary.                                                      !
                  !------------------------------------------------------------------------!
                  if (iqoutput > 0) then
                     cpatch%qmean_par_l       (it,ico) = cpatch%qmean_par_l       (it,ico) &
                                                       + cpatch%par_l                (ico)
                     cpatch%qmean_par_l_beam  (it,ico) = cpatch%qmean_par_l_beam  (it,ico) &
                                                       + cpatch%par_l_beam           (ico)
                     cpatch%qmean_par_l_diff  (it,ico) = cpatch%qmean_par_l_diff  (it,ico) &
                                                       + cpatch%par_l_diffuse        (ico)
                     cpatch%qmean_fs_open     (it,ico) = cpatch%qmean_fs_open     (it,ico) &
                                                       + cpatch%fs_open              (ico)
                     cpatch%qmean_fsw         (it,ico) = cpatch%qmean_fsw         (it,ico) &
                                                       + cpatch%fsw                  (ico)
                     cpatch%qmean_fsn         (it,ico) = cpatch%qmean_fsn         (it,ico) &
                                                       + cpatch%fsn                  (ico)
                     cpatch%qmean_psi_open    (it,ico) = cpatch%qmean_psi_open    (it,ico) &
                                                       + cpatch%psi_open             (ico)
                     cpatch%qmean_psi_closed  (it,ico) = cpatch%qmean_psi_closed  (it,ico) &
                                                       + cpatch%psi_closed           (ico)
                     cpatch%qmean_water_supply(it,ico) = cpatch%qmean_water_supply(it,ico) &
                                                       + cpatch%water_supply         (ico)
                  end if
                  !------------------------------------------------------------------------!

               end if
            end do

            if (rshort_tot > rshort_twilight_min) then
               csite%dmean_lambda_light(ipa) = csite%dmean_lambda_light(ipa)               &
                                             + csite%lambda_light(ipa)
            end if
         end do patchloop

         !---------------------------------------------------------------------------------!
         !    Variables already average at the sitetype level, just add them to polygon-   !
         ! type level.                                                                     !
         !---------------------------------------------------------------------------------!

         site_lai        = site_lai        + patch_lai * cpoly%area(isi)
         site_lma        = site_lma        + patch_lma * cpoly%area(isi)
         sss_leaf_energy = sss_leaf_energy + (pss_leaf_energy*site_area_i)*cpoly%area(isi)
         sss_leaf_water  = sss_leaf_water  + (pss_leaf_water *site_area_i)*cpoly%area(isi)
         sss_leaf_hcap   = sss_leaf_hcap   + (pss_leaf_hcap  *site_area_i)*cpoly%area(isi)
         sss_wood_energy = sss_wood_energy + (pss_wood_energy*site_area_i)*cpoly%area(isi)
         sss_wood_water  = sss_wood_water  + (pss_wood_water *site_area_i)*cpoly%area(isi)
         sss_wood_hcap   = sss_wood_hcap   + (pss_wood_hcap  *site_area_i)*cpoly%area(isi)


         sss_can_theta  = sss_can_theta + cpoly%area(isi)                                  &
                                        * ( sum(csite%can_theta   * csite%area)            &
                                          * site_area_i )
         sss_can_theiv  = sss_can_theiv + cpoly%area(isi)                                  &
                                        * ( sum(csite%can_theiv   * csite%area)            &
                                          * site_area_i )
         sss_can_shv    = sss_can_shv   + cpoly%area(isi)                                  &
                                        * ( sum(csite%can_shv     * csite%area)            &
                                          * site_area_i )
         sss_can_co2    = sss_can_co2   + cpoly%area(isi)                                  &
                                        * ( sum(csite%can_co2     * csite%area)            &
                                          * site_area_i )
         sss_can_prss   = sss_can_prss  + cpoly%area(isi)                                  &
                                        * ( sum(csite%can_prss    * csite%area)            &
                                          * site_area_i )
         sss_gnd_temp   = sss_gnd_temp  + cpoly%area(isi)                                  &
                                        * ( sum(csite%ground_temp * csite%area)            &
                                          * site_area_i )
         sss_gnd_shv    = sss_gnd_shv   + cpoly%area(isi)                                  &
                                        * ( sum(csite%ground_shv  * csite%area)            &
                                          * site_area_i )
      end do siteloop
      
      !------------------------------------------------------------------------------------!
      !    Variables already averaged at the polygontype level, just add them to edtype    !
      ! level.                                                                             !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_leaf_energy(ipy)  = cgrid%dmean_leaf_energy(ipy)                         &
                                    + sss_leaf_energy * poly_area_i
      cgrid%dmean_leaf_water(ipy)   = cgrid%dmean_leaf_water(ipy)                          &
                                    + sss_leaf_water * poly_area_i
      cgrid%dmean_leaf_hcap(ipy)    = cgrid%dmean_leaf_hcap(ipy)                           &
                                    + sss_leaf_hcap  * poly_area_i
      cgrid%dmean_wood_energy(ipy)  = cgrid%dmean_wood_energy(ipy)                         &
                                    + sss_wood_energy * poly_area_i
      cgrid%dmean_wood_water(ipy)   = cgrid%dmean_wood_water(ipy)                          &
                                    + sss_wood_water * poly_area_i
      cgrid%dmean_wood_hcap(ipy)    = cgrid%dmean_wood_hcap(ipy)                           &
                                    + sss_wood_hcap  * poly_area_i
      cgrid%dmean_can_theta(ipy)    = cgrid%dmean_can_theta(ipy)                           &
                                    + sss_can_theta * poly_area_i
      cgrid%dmean_can_theiv(ipy)    = cgrid%dmean_can_theiv(ipy)                           &
                                    + sss_can_theiv * poly_area_i
      cgrid%dmean_can_shv(ipy)      = cgrid%dmean_can_shv(ipy)                             &
                                    + sss_can_shv   * poly_area_i
      cgrid%dmean_can_co2(ipy)      = cgrid%dmean_can_co2(ipy)                             &
                                    + sss_can_co2   * poly_area_i
      cgrid%dmean_can_prss(ipy)     = cgrid%dmean_can_prss(ipy)                            &
                                    + sss_can_prss  * poly_area_i
      cgrid%dmean_gnd_temp(ipy)     = cgrid%dmean_gnd_temp(ipy)                            &
                                    + sss_gnd_temp  * poly_area_i
      cgrid%dmean_gnd_shv (ipy)     = cgrid%dmean_gnd_shv (ipy)                            &
                                    + sss_gnd_shv   * poly_area_i

      !------------------------------------------------------------------------------------!
      !    Meteorological forcing variables:  even though we read in at polygon level,     !
      ! here we run the site-level weighted average, because this is what ED really uses.  !
      !------------------------------------------------------------------------------------!
      do isi=1,cpoly%nsites
         cgrid%dmean_atm_temp(ipy)      = cgrid%dmean_atm_temp(ipy)                        &
                                        + cpoly%met(isi)%atm_tmp                           &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_rshort(ipy)        = cgrid%dmean_rshort(ipy)                          &
                                        + ( cpoly%met(isi)%nir_beam                        &
                                          + cpoly%met(isi)%nir_diffuse                     &
                                          + cpoly%met(isi)%par_beam                        &
                                          + cpoly%met(isi)%par_diffuse )                   &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_rshort_diff(ipy)   = cgrid%dmean_rshort_diff(ipy)                     &
                                        + ( cpoly%met(isi)%nir_diffuse                     &
                                          + cpoly%met(isi)%par_diffuse )                   &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_rlong(ipy)         = cgrid%dmean_rlong(ipy)                           &
                                        + cpoly%met(isi)%rlong                             &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_atm_shv (ipy)      = cgrid%dmean_atm_shv (ipy)                        &
                                        + cpoly%met(isi)%atm_shv                           &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_atm_co2 (ipy)      = cgrid%dmean_atm_co2 (ipy)                        &
                                        + cpoly%met(isi)%atm_co2                           &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_atm_prss(ipy)      = cgrid%dmean_atm_prss(ipy)                        &
                                        + cpoly%met(isi)%prss                              &
                                        * cpoly%area(isi) * poly_area_i
         cgrid%dmean_atm_vels(ipy)      = cgrid%dmean_atm_vels(ipy)                        &
                                        + cpoly%met(isi)%vels                              &
                                        * cpoly%area(isi) * poly_area_i
      end do
      !------------------------------------------------------------------------------------!

    
      if(site_lai > tiny(1.))then
         cgrid%avg_lma(ipy) = site_lma/site_lai
      end if

      !------------------------------------------------------------------------------------!
      !     If we are integrating the mean diurnal cycle, sum them too.                    !
      !------------------------------------------------------------------------------------!
      if (iqoutput > 0) then
         cgrid%qmean_leaf_energy   (it,ipy) = cgrid%qmean_leaf_energy             (it,ipy) &
                                            + sss_leaf_energy * poly_area_i
         cgrid%qmean_leaf_water    (it,ipy) = cgrid%qmean_leaf_water              (it,ipy) &
                                            + sss_leaf_water * poly_area_i
         cgrid%qmean_leaf_hcap     (it,ipy) = cgrid%qmean_leaf_hcap               (it,ipy) &
                                            + sss_leaf_hcap  * poly_area_i
         cgrid%qmean_wood_energy   (it,ipy) = cgrid%qmean_wood_energy             (it,ipy) &
                                            + sss_wood_energy * poly_area_i
         cgrid%qmean_wood_water    (it,ipy) = cgrid%qmean_wood_water              (it,ipy) &
                                            + sss_wood_water * poly_area_i
         cgrid%qmean_wood_hcap     (it,ipy) = cgrid%qmean_wood_hcap               (it,ipy) &
                                            + sss_wood_hcap  * poly_area_i
         cgrid%qmean_can_theta     (it,ipy) = cgrid%qmean_can_theta               (it,ipy) &
                                            + sss_can_theta * poly_area_i
         cgrid%qmean_can_theiv     (it,ipy) = cgrid%qmean_can_theiv               (it,ipy) &
                                            + sss_can_theiv * poly_area_i
         cgrid%qmean_can_shv       (it,ipy) = cgrid%qmean_can_shv                 (it,ipy) &
                                            + sss_can_shv   * poly_area_i
         cgrid%qmean_can_co2       (it,ipy) = cgrid%qmean_can_co2                 (it,ipy) &
                                            + sss_can_co2   * poly_area_i
         cgrid%qmean_can_prss      (it,ipy) = cgrid%qmean_can_prss                (it,ipy) &
                                            + sss_can_prss  * poly_area_i
         cgrid%qmean_gnd_temp      (it,ipy) = cgrid%qmean_gnd_temp                (it,ipy) &
                                            + sss_gnd_temp  * poly_area_i
         cgrid%qmean_gnd_shv       (it,ipy) = cgrid%qmean_gnd_shv                 (it,ipy) &
                                            + sss_gnd_shv   * poly_area_i
         !---------------------------------------------------------------------------------!
         !    Meteorological forcing variables:  even though we read in at polygon level,  !
         ! here we run the site-level weighted average, because this is what ED really     !
         ! uses.                                                                           !
         !---------------------------------------------------------------------------------!
         do isi=1,cpoly%nsites
            cgrid%qmean_atm_temp      (it,ipy) = cgrid%qmean_atm_temp             (it,ipy) &
                                               + cpoly%met(isi)%atm_tmp                    &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_rshort        (it,ipy) = cgrid%qmean_rshort               (it,ipy) &
                                               + ( cpoly%met(isi)%nir_beam                 &
                                                 + cpoly%met(isi)%nir_diffuse              &
                                                 + cpoly%met(isi)%par_beam                 &
                                                 + cpoly%met(isi)%par_diffuse )            &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_rshort_diff   (it,ipy) = cgrid%qmean_rshort_diff          (it,ipy) &
                                               + ( cpoly%met(isi)%nir_diffuse              &
                                                 + cpoly%met(isi)%par_diffuse )            &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_rlong         (it,ipy) = cgrid%qmean_rlong                (it,ipy) &
                                               + cpoly%met(isi)%rlong                      &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_atm_shv       (it,ipy) = cgrid%qmean_atm_shv              (it,ipy) &
                                               + cpoly%met(isi)%atm_shv                    &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_atm_co2       (it,ipy) = cgrid%qmean_atm_co2              (it,ipy) &
                                               + cpoly%met(isi)%atm_co2                    &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_atm_prss      (it,ipy) = cgrid%qmean_atm_prss             (it,ipy) &
                                               + cpoly%met(isi)%prss                       &
                                               * cpoly%area(isi) * poly_area_i
            cgrid%qmean_atm_vels      (it,ipy) = cgrid%qmean_atm_vels             (it,ipy) &
                                               + cpoly%met(isi)%vels                       &
                                               * cpoly%area(isi) * poly_area_i
         end do
      end if
   end do polyloop
      
   return
end subroutine integrate_ed_daily_output_state
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine integrates the daily average.  This is called at the analysis time in !
! case at least one of daily means, monthly means, and mean diurnal cycle output is        !
! sought.  We take advantage from the previously averaged variables.                       !
!------------------------------------------------------------------------------------------!
subroutine integrate_ed_daily_output_flux(cgrid)
   use ed_state_vars        , only : edtype              & ! structure
                                   , polygontype         & ! structure
                                   , sitetype            & ! structure
                                   , patchtype           ! ! structure
   use grid_coms            , only : nzg                 & ! intent(in)
                                   , time                ! ! intent(in)
   use ed_max_dims          , only : n_dbh               & ! intent(in)
                                   , n_pft               & ! intent(in)
                                   , n_dist_types        ! ! intent(in)
   use ed_misc_coms         , only : dtlsm               & ! intent(in)
                                   , ddbhi               & ! intent(in)
                                   , iqoutput            & ! intent(in)
                                   , frqfast             & ! intent(in)
                                   , ndcycle             & ! intent(in)
                                   , current_time        ! ! intent(in)
   use pft_coms             , only : c2n_leaf            ! ! intent(in)
   use consts_coms          , only : day_sec             & ! intent(in)
                                   , umols_2_kgCyr       ! ! intent(in)
   use canopy_radiation_coms, only : rshort_twilight_min & ! intent(in)
                                   , cosz_min            ! ! intent(in)
   implicit none
   
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)      , target    :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer   :: cpoly
   type(sitetype)    , pointer   :: csite
   type(patchtype)   , pointer   :: cpatch
   integer                       :: ipy
   integer                       :: isi
   integer                       :: ipa
   integer                       :: ico
   integer                       :: k
   integer                       :: ilu
   integer                       :: idbh
   integer                       :: it
   real                          :: poly_area_i
   real                          :: site_area_i
   real                          :: forest_site
   real                          :: forest_site_i
   real                          :: forest_poly
   real                          :: sitesum_gpp
   real, dimension(n_dbh)        :: sitesum_gpp_dbh
   real, dimension(n_dist_types) :: sitesum_gpp_lu
   real, dimension(n_dist_types) :: sitesum_rh_lu
   real, dimension(n_dist_types) :: sitesum_nep_lu
   real                          :: sitesum_rh
   real                          :: sitesum_evap
   real                          :: sitesum_transp
   real                          :: sitesum_Nuptake
   real                          :: sitesum_sensible_tot
   real                          :: sitesum_co2_residual
   real                          :: sitesum_energy_residual
   real                          :: sitesum_water_residual
   real                          :: sitesum_root_litter
   real                          :: sitesum_leaf_litter
   real                          :: sitesum_root_litterN
   real                          :: sitesum_leaf_litterN
   real                          :: sitesum_Nmin_input
   real                          :: sitesum_Nmin_loss
   real                          :: patchsum_leaf_litter
   real                          :: patchsum_root_litter
   real                          :: patchsum_leaf_litterN
   real                          :: patchsum_root_litterN
   real                          :: rshort_tot
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the index corresponding to this time of the day for the mean diurnal cycle   !
   ! averages.                                                                             !
   !---------------------------------------------------------------------------------------!
   it = ceiling(mod(current_time%time,day_sec)/frqfast)
   if (it == 0) it = nint(day_sec/frqfast)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   !
   !---------------------------------------------------------------------------------------!
   !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.        !
   ! This works in off-line runs, but it causes memory leaks (and crashes) in the coupled  !
   ! runs over the ocean, where cgrid%npolygons can be 0 if one of the sub-domains falls   !
   ! entirely over the ocean.  Thanks!                                                     !
   !---------------------------------------------------------------------------------------!
   ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the safe
   !                 !      way of initialising the variable.
   !---------------------------------------------------------------------------------------!
   polyloop: do ipy=1,cgrid%npolygons

      !------------------------------------------------------------------------------------!
      !     This is the right and safe place to initialise polygon-level (cgrid) vari-     !
      ! ables, so in case npolygons is zero this will not cause memory leaks.  I know,     !
      ! this never happens in off-line runs, but it is quite common in coupled runs...     !
      ! Whenever one of the nodes receives a sub-domain where all the points are over the  !
      ! ocean, ED will not assign any polygon in that sub-domain, which means that that    !
      ! node will have 0 polygons, and the variables cannot be allocated.  If you try to   !
      ! access the polygon level variable outside the loop, then the model crashes due to  !
      ! segmentation violation (a bad thing), whereas by putting the variables here both   !
      ! the off-line model and the coupled runs will work, because this loop will be skip- !
      ! ped when there is no polygon.                                                      !
      !------------------------------------------------------------------------------------!
      cgrid%Nbiomass_uptake  (ipy) = 0.
      cgrid%Cleaf_litter_flux(ipy) = 0.
      cgrid%Croot_litter_flux(ipy) = 0.
      cgrid%Nleaf_litter_flux(ipy) = 0.
      cgrid%Nroot_litter_flux(ipy) = 0.
      cgrid%Ngross_min       (ipy) = 0.
      cgrid%Nnet_min         (ipy) = 0.
      !------------------------------------------------------------------------------------!


      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      !---- Initialize auxiliary variables to add sitetype variables. ---------------------!
      sitesum_leaf_litter     = 0.
      sitesum_root_litter     = 0.
      sitesum_leaf_litterN    = 0.
      sitesum_root_litterN    = 0.
      sitesum_Nuptake         = 0.
      sitesum_Nmin_loss       = 0.
      sitesum_Nmin_input      = 0.
      sitesum_gpp_dbh         = 0.
      sitesum_evap            = 0.
      sitesum_transp          = 0.
      sitesum_co2_residual    = 0.
      sitesum_water_residual  = 0.
      sitesum_energy_residual = 0.
      forest_poly             = 0.

      siteloop: do isi=1, cpoly%nsites
         !---------------------------------------------------------------------------------!
         !     Determine if there was any beam radiation, and compute the total (I think   !
         ! rshort is always the total, but not sure.  This is the same test done in the    !
         ! radiation driver.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpoly%cosaoi(isi) <= cosz_min) then
            rshort_tot = cpoly%met(isi)%rshort_diffuse
         else
            rshort_tot = cpoly%met(isi)%rshort
         end if
         !---------------------------------------------------------------------------------!


         csite => cpoly%site(isi)
         
         !----- Inverse of total site area (sum of all patches' area). --------------------!
         site_area_i=1./sum(csite%area)

         !----- Forest areas. -------------------------------------------------------------!
         forest_site           = sum(csite%area,csite%dist_type /= 1)
         if (forest_site > 1.0e-6) then
            forest_site_i         = 1./forest_site
         else
            forest_site_i         = 0.0
         end if
         forest_poly           = forest_poly + forest_site

         !----- Initialize auxiliary variables to add patchtype variables. ----------------!
         patchsum_leaf_litter  = 0.
         patchsum_root_litter  = 0.
         patchsum_leaf_litterN = 0.
         patchsum_root_litterN = 0.

         !----- Looping through the patches to normalize the sum of all cohorts. ----------!
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)
            if (cpatch%ncohorts > 0) then


               !---------------------------------------------------------------------------!
               !     Integrate the mean diurnal cycle only if the mean diurnal cycle       !
               ! output is sought.  This is also done for night time as the length of      !
               ! the day may vary.                                                         !
               !---------------------------------------------------------------------------!
               if (iqoutput > 0) then
                  do ico=1,cpatch%ncohorts
                     if (cpatch%leaf_resolvable(ico)) then
                        cpatch%qmean_gpp      (it,ico) = cpatch%qmean_gpp         (it,ico) &
                                                       + cpatch%mean_gpp             (ico) &
                                                       * umols_2_kgCyr                     &
                                                       / cpatch%nplant(ico)
                        cpatch%qmean_leaf_resp(it,ico) = cpatch%qmean_leaf_resp   (it,ico) &
                                                       + cpatch%mean_leaf_resp       (ico) &
                                                       * umols_2_kgCyr                     &
                                                       / cpatch%nplant(ico)
                        cpatch%qmean_root_resp(it,ico) = cpatch%qmean_root_resp   (it,ico) &
                                                       + cpatch%mean_root_resp       (ico) &
                                                       * umols_2_kgCyr                     &
                                                       / cpatch%nplant(ico)
                     end if
                  end do
               end if
               patchsum_root_litter  = patchsum_root_litter                                &
                                     + sum(cpatch%root_maintenance) * csite%area(ipa) 
               patchsum_leaf_litter  = patchsum_leaf_litter                                &
                                     + sum(cpatch%leaf_maintenance) * csite%area(ipa) 
               patchsum_root_litterN = patchsum_root_litterN                               &
                                     + sum(cpatch%root_maintenance*c2n_leaf(cpatch%pft))   &
                                     * csite%area(ipa) 
               patchsum_leaf_litterN = patchsum_leaf_litterN                               &
                                     + sum(cpatch%leaf_maintenance*c2n_leaf(cpatch%pft))   & 
                                     * csite%area(ipa)
             end if
            csite%dmean_co2_residual      (ipa) = csite%dmean_co2_residual        (ipa)    &
                                                + csite%co2budget_residual        (ipa)
            csite%dmean_energy_residual   (ipa) = csite%dmean_energy_residual     (ipa)    &
                                                + csite%ebudget_residual          (ipa)
            csite%dmean_water_residual    (ipa) = csite%dmean_water_residual      (ipa)    &
                                                + csite%wbudget_residual          (ipa)
            csite%dmean_rh                (ipa) = csite%dmean_rh                  (ipa)    &
                                                + csite%co2budget_rh              (ipa)    &
                                                * umols_2_kgCyr
            csite%dmean_rk4step           (ipa) = csite%dmean_rk4step             (ipa)    &
                                                + csite%avg_rk4step               (ipa)
            if (iqoutput > 0) then
               csite%qmean_rh            (it,ipa)  = csite%qmean_rh            (it,ipa)    &
                                                   + csite%co2budget_rh           (ipa)    &
                                                   * umols_2_kgCyr
               csite%qmean_albedo        (it,ipa)  = csite%qmean_albedo        (it,ipa)    &
                                                   + csite%avg_albedo             (ipa)
               csite%qmean_albedo_beam   (it,ipa)  = csite%qmean_albedo_beam   (it,ipa)    &
                                                   + csite%avg_albedo_beam        (ipa)
               csite%qmean_albedo_diffuse(it,ipa)  = csite%qmean_albedo_diffuse(it,ipa)    &
                                                   + csite%avg_albedo_diffuse     (ipa)
            end if
            
            
            if (rshort_tot > rshort_twilight_min) then
               csite%dmean_albedo        (ipa)  = csite%dmean_albedo        (ipa)          &
                                                + csite%avg_albedo          (ipa)
               csite%dmean_albedo_beam   (ipa)  = csite%dmean_albedo_beam   (ipa)          &
                                                + csite%avg_albedo_beam     (ipa)
               csite%dmean_albedo_diffuse(ipa)  = csite%dmean_albedo_diffuse(ipa)          &
                                                + csite%avg_albedo_diffuse  (ipa)
            end if
         end do patchloop
         
         !---------------------------------------------------------------------------------!
         !      Variables already averaged at the sitetype level, add them to polygontype  !
         ! level.                                                                          !
         !---------------------------------------------------------------------------------!         
         sitesum_evap            = sitesum_evap                                            &
                                 + (sum(csite%avg_evap * csite%area) * site_area_i)        &
                                 * cpoly%area(isi)
         sitesum_transp          = sitesum_transp                                          &
                                 + (sum(csite%avg_transp * csite%area) * site_area_i)      &
                                 * cpoly%area(isi)
         
         sitesum_co2_residual    = sitesum_co2_residual                                    &
                                 + (sum(csite%co2budget_residual*csite%area)*site_area_i)  &
                                 * cpoly%area(isi)
         sitesum_water_residual  = sitesum_water_residual                                  &
                                 + (sum(csite%wbudget_residual*csite%area)*site_area_i)    &
                                 * cpoly%area(isi)
         sitesum_energy_residual = sitesum_energy_residual                                 &
                                 + (sum(csite%ebudget_residual*csite%area)*site_area_i)    &
                                 * cpoly%area(isi)

         cpoly%dmean_co2_residual(isi)    = cpoly%dmean_co2_residual(isi)                  &
                                          + sum(csite%co2budget_residual*csite%area)       &
                                          * site_area_i
         cpoly%dmean_energy_residual(isi) = cpoly%dmean_energy_residual(isi)               &
                                          + sum(csite%ebudget_residual  *csite%area)       &
                                          * site_area_i
         cpoly%dmean_water_residual(isi)  = cpoly%dmean_water_residual(isi)                &
                                          + sum(csite%wbudget_residual  *csite%area)       &
                                          * site_area_i

         sitesum_root_litter    = sitesum_root_litter                                      &
                                + (patchsum_root_litter * site_area_i) * cpoly%area(isi)
         sitesum_leaf_litter    = sitesum_leaf_litter                                      &
                                + (patchsum_leaf_litter * site_area_i) * cpoly%area(isi)

         sitesum_root_litterN    = sitesum_root_litterN                                    &
                                 + (patchsum_root_litterN * site_area_i) * cpoly%area(isi)
         sitesum_leaf_litterN    = sitesum_leaf_litterN                                    &
                                 + (patchsum_leaf_litterN * site_area_i) * cpoly%area(isi)
         sitesum_Nmin_loss       = sitesum_Nmin_loss                                       &
                                 + sum(csite%mineralized_N_loss*csite%area) * site_area_i  &
                                 * cpoly%area(isi)
         sitesum_Nmin_input      = sitesum_Nmin_input                                      &
                                 + sum(csite%mineralized_N_input*csite%area)* site_area_i  &
                                 * cpoly%area(isi)
         sitesum_Nuptake         = sitesum_Nuptake                                         &
                                 + (sum(csite%total_plant_nitrogen_uptake * csite%area)    &
                                   *site_area_i ) * cpoly%area(isi)

         dbhloop: do idbh=1,n_dbh
            sitesum_gpp_dbh(idbh) = sitesum_gpp_dbh(idbh)                                  &
                                  + ( sum(csite%co2budget_gpp_dbh(idbh,:) * csite%area)    &
                                    * site_area_i) * cpoly%area(isi)
         end do dbhloop

      end do siteloop

      cgrid%Croot_litter_flux(ipy)     = cgrid%Croot_litter_flux(ipy)                      &
                                       + sitesum_root_litter * poly_area_i
      cgrid%Cleaf_litter_flux(ipy)     = cgrid%Cleaf_litter_flux(ipy)                      &
                                       + sitesum_leaf_litter * poly_area_i

      cgrid%Nroot_litter_flux(ipy)     = cgrid%Nroot_litter_flux(ipy)                      &
                                       + sitesum_root_litterN * poly_area_i
      cgrid%Nleaf_litter_flux(ipy)     = cgrid%Nleaf_litter_flux(ipy)                      &
                                       + sitesum_leaf_litterN    * poly_area_i

      cgrid%Ngross_min(ipy)            = cgrid%Ngross_min(ipy)                             &
                                       + sitesum_Nmin_input * poly_area_i
      cgrid%Nnet_min(ipy)              = cgrid%Nnet_min(ipy)                               &
                                       + (sitesum_Nmin_input-sitesum_Nmin_loss)            &
                                       * poly_area_i
      cgrid%Nbiomass_uptake(ipy)       = cgrid%Nbiomass_uptake(ipy)                        &
                                       + sitesum_Nuptake*poly_area_i


      cgrid%dmean_co2_residual(ipy)    = cgrid%dmean_co2_residual(ipy)                     &
                                       + sitesum_co2_residual * poly_area_i * umols_2_kgCyr
      cgrid%dmean_energy_residual(ipy) = cgrid%dmean_energy_residual(ipy)                  &
                                       + sitesum_energy_residual * poly_area_i
      cgrid%dmean_water_residual(ipy)  = cgrid%dmean_water_residual(ipy)                   &
                                       + sitesum_water_residual * poly_area_i
      
      do idbh=1,n_dbh
         cgrid%dmean_gpp_dbh(idbh,ipy) = cgrid%dmean_gpp_dbh(idbh,ipy)                     &
                                       + sitesum_gpp_dbh(idbh)* poly_area_i
      end do

      !----- These variables are already averaged at gridtype, just add them up. ----------!
      do k=1,nzg
         cgrid%dmean_soil_temp(k,ipy)   = cgrid%dmean_soil_temp(k,ipy)                     &
                                        + cgrid%avg_soil_temp(k,ipy) 
         cgrid%dmean_soil_water(k,ipy)  = cgrid%dmean_soil_water(k,ipy)                    &
                                        + cgrid%avg_soil_water(k,ipy)
         cgrid%dmean_soil_mstpot(k,ipy) = cgrid%dmean_soil_mstpot(k,ipy)                   &
                                        + cgrid%avg_soil_mstpot(k,ipy)
         cgrid%dmean_transloss(k,ipy)   = cgrid%dmean_transloss(k,ipy)                     &
                                        + cgrid%avg_transloss(k,ipy)
       end do
      cgrid%dmean_sensible_lc(ipy) = cgrid%dmean_sensible_lc(ipy)                          &
                                   + cgrid%avg_sensible_lc(ipy) 
      cgrid%dmean_sensible_wc(ipy) = cgrid%dmean_sensible_wc(ipy)                          &
                                   + cgrid%avg_sensible_wc(ipy) 
      cgrid%dmean_sensible_gc(ipy) = cgrid%dmean_sensible_gc(ipy)                          &
                                   + cgrid%avg_sensible_gc(ipy)
      cgrid%dmean_sensible_ac(ipy) = cgrid%dmean_sensible_ac(ipy)                          &
                                   + cgrid%avg_sensible_ac(ipy)

      cgrid%dmean_rshort_gnd  (ipy) = cgrid%dmean_rshort_gnd   (ipy)                       &
                                    + cgrid%avg_rshort_gnd     (ipy)
      cgrid%dmean_rlong_gnd   (ipy) = cgrid%dmean_rlong_gnd    (ipy)                       &
                                    + cgrid%avg_rlong_gnd      (ipy)
      cgrid%dmean_rlongup     (ipy) = cgrid%dmean_rlongup      (ipy)                       &
                                    + cgrid%avg_rlongup        (ipy)
      cgrid%dmean_rlong_albedo(ipy) = cgrid%dmean_rlong_albedo (ipy)                       &
                                    + cgrid%avg_rlong_albedo   (ipy)
      !------ Integrate the NEE with the conventional NEE sign (< 0 = uptake). ------------!
      cgrid%dmean_nee       (ipy) = cgrid%dmean_nee         (ipy)                          &
                                  - cgrid%avg_carbon_ac     (ipy) * umols_2_kgCyr
      cgrid%dmean_leaf_resp (ipy) = cgrid%dmean_leaf_resp   (ipy)                          &
                                  + cgrid%avg_leaf_resp     (ipy) * umols_2_kgCyr
      cgrid%dmean_root_resp (ipy) = cgrid%dmean_root_resp   (ipy)                          &
                                  + cgrid%avg_root_resp     (ipy) * umols_2_kgCyr
      cgrid%dmean_plresp    (ipy) = cgrid%dmean_plresp      (ipy)                          &
                                  + ( cgrid%avg_leaf_resp   (ipy)                            &
                                    + cgrid%avg_root_resp   (ipy)                          &
                                    + cgrid%avg_growth_resp (ipy)                          &
                                    + cgrid%avg_storage_resp(ipy)                          &
                                    + cgrid%avg_vleaf_resp  (ipy) ) * umols_2_kgCyr
      cgrid%dmean_nep       (ipy) = cgrid%dmean_nep         (ipy)                          &
                                  + ( cgrid%avg_gpp         (ipy)                            &
                                    - cgrid%avg_leaf_resp   (ipy)                          &
                                    - cgrid%avg_root_resp   (ipy)                          &
                                    - cgrid%avg_growth_resp (ipy)                          &
                                    - cgrid%avg_storage_resp(ipy)                          &
                                    - cgrid%avg_vleaf_resp  (ipy)                          &
                                    - cgrid%avg_htroph_resp (ipy) ) * umols_2_kgCyr
      cgrid%dmean_pcpg      (ipy) = cgrid%dmean_pcpg      (ipy) + cgrid%avg_pcpg      (ipy)
      cgrid%dmean_evap      (ipy) = cgrid%dmean_evap      (ipy) + cgrid%avg_evap      (ipy)
      cgrid%dmean_transp    (ipy) = cgrid%dmean_transp    (ipy) + cgrid%avg_transp    (ipy)
      cgrid%dmean_runoff    (ipy) = cgrid%dmean_runoff    (ipy) + cgrid%avg_runoff    (ipy)
      cgrid%dmean_drainage  (ipy) = cgrid%dmean_drainage  (ipy) + cgrid%avg_drainage  (ipy)
      cgrid%dmean_vapor_lc  (ipy) = cgrid%dmean_vapor_lc  (ipy) + cgrid%avg_vapor_lc  (ipy)
      cgrid%dmean_vapor_wc  (ipy) = cgrid%dmean_vapor_wc  (ipy) + cgrid%avg_vapor_wc  (ipy)
      cgrid%dmean_vapor_gc  (ipy) = cgrid%dmean_vapor_gc  (ipy) + cgrid%avg_vapor_gc  (ipy)
      cgrid%dmean_vapor_ac  (ipy) = cgrid%dmean_vapor_ac  (ipy) + cgrid%avg_vapor_ac  (ipy)
      
      !------------------------------------------------------------------------------------!
      !     Integrate the mean diurnal cycle in case the mean diurnal cycle is sought.     !
      !------------------------------------------------------------------------------------!
      if (iqoutput > 0) then
         !------ Use the local site sum to integrate the following variables. -------------!
         cgrid%qmean_gpp            (it,ipy) = cgrid%qmean_gpp                   (it,ipy)  &
                                             + cgrid%avg_gpp                        (ipy)  &
                                             * umols_2_kgCyr
         cgrid%qmean_leaf_resp      (it,ipy) = cgrid%qmean_leaf_resp             (it,ipy)  &
                                             + cgrid%avg_leaf_resp                  (ipy)  &
                                             * umols_2_kgCyr
         cgrid%qmean_root_resp      (it,ipy) = cgrid%qmean_root_resp             (it,ipy)  &
                                             + cgrid%avg_root_resp                  (ipy)  &
                                             * umols_2_kgCyr
         cgrid%qmean_plresp         (it,ipy) = cgrid%qmean_plresp                (it,ipy)  &
                                             + cgrid%avg_plant_resp                 (ipy)  &
                                             * umols_2_kgCyr
         cgrid%qmean_nep            (it,ipy) = cgrid%qmean_nep                   (it,ipy)  &
                                             + ( cgrid%avg_gpp                      (ipy)  &
                                               - cgrid%avg_leaf_resp                (ipy)  &
                                               - cgrid%avg_root_resp                (ipy)  &
                                               - cgrid%avg_growth_resp              (ipy)  &
                                               - cgrid%avg_storage_resp             (ipy)  &
                                               - cgrid%avg_vleaf_resp               (ipy)  &
                                               - cgrid%avg_htroph_resp              (ipy)) &
                                             * umols_2_kgCyr
         cgrid%qmean_rh             (it,ipy) = cgrid%qmean_rh                    (it,ipy)  &
                                             + cgrid%avg_htroph_resp                (ipy)  &
                                             * umols_2_kgCyr
         !------ Integrate the NEE with the conventional NEE sign (< 0 = uptake). ---------!
         cgrid%qmean_nee            (it,ipy) = cgrid%qmean_nee                    (it,ipy) &
                                             - cgrid%avg_carbon_ac                   (ipy) &
                                             * umols_2_kgCyr
         !----- Variables that were previously integrated to this time step. --------------!
         cgrid%qmean_sensible_lc    (it,ipy) = cgrid%qmean_sensible_lc            (it,ipy) &
                                             + cgrid%avg_sensible_lc                 (ipy)
         cgrid%qmean_sensible_wc    (it,ipy) = cgrid%qmean_sensible_wc            (it,ipy) &
                                             + cgrid%avg_sensible_wc                 (ipy)
         cgrid%qmean_sensible_gc    (it,ipy) = cgrid%qmean_sensible_gc            (it,ipy) &
                                             + cgrid%avg_sensible_gc                 (ipy)
         cgrid%qmean_sensible_ac    (it,ipy) = cgrid%qmean_sensible_ac            (it,ipy) &
                                             + cgrid%avg_sensible_ac                 (ipy)
         cgrid%qmean_rshort_gnd     (it,ipy) = cgrid%qmean_rshort_gnd             (it,ipy) &
                                             + cgrid%avg_rshort_gnd                  (ipy)
         cgrid%qmean_rlong_gnd      (it,ipy) = cgrid%qmean_rlong_gnd              (it,ipy) &
                                             + cgrid%avg_rlong_gnd                   (ipy)
         cgrid%qmean_rlongup        (it,ipy) = cgrid%qmean_rlongup                (it,ipy) &
                                             + cgrid%avg_rlongup                     (ipy)
         cgrid%qmean_rlong_albedo   (it,ipy) = cgrid%qmean_rlong_albedo           (it,ipy) &
                                             + cgrid%avg_rlong_albedo                (ipy)
         cgrid%qmean_albedo         (it,ipy) = cgrid%qmean_albedo                 (it,ipy) &
                                             + cgrid%avg_albedo                      (ipy)
         cgrid%qmean_albedo_beam    (it,ipy) = cgrid%qmean_albedo_beam            (it,ipy) &
                                             + cgrid%avg_albedo_beam                 (ipy)
         cgrid%qmean_albedo_diffuse (it,ipy) = cgrid%qmean_albedo_diffuse         (it,ipy) &
                                             + cgrid%avg_albedo_diffuse              (ipy)
         cgrid%qmean_pcpg           (it,ipy) = cgrid%qmean_pcpg                   (it,ipy) &
                                             + cgrid%avg_pcpg                        (ipy)
         cgrid%qmean_evap           (it,ipy) = cgrid%qmean_evap                   (it,ipy) &
                                             + cgrid%avg_evap                        (ipy)
         cgrid%qmean_transp         (it,ipy) = cgrid%qmean_transp                 (it,ipy) &
                                             + cgrid%avg_transp                      (ipy)
         cgrid%qmean_runoff         (it,ipy) = cgrid%qmean_runoff                 (it,ipy) &
                                             + cgrid%avg_runoff                      (ipy)
         cgrid%qmean_drainage       (it,ipy) = cgrid%qmean_drainage               (it,ipy) &
                                             + cgrid%avg_drainage                    (ipy)
         cgrid%qmean_vapor_lc       (it,ipy) = cgrid%qmean_vapor_lc               (it,ipy) &
                                             + cgrid%avg_vapor_lc                    (ipy)
         cgrid%qmean_vapor_wc       (it,ipy) = cgrid%qmean_vapor_wc               (it,ipy) &
                                             + cgrid%avg_vapor_wc                    (ipy)
         cgrid%qmean_vapor_gc       (it,ipy) = cgrid%qmean_vapor_gc               (it,ipy) &
                                             + cgrid%avg_vapor_gc                    (ipy)
         cgrid%qmean_vapor_ac       (it,ipy) = cgrid%qmean_vapor_ac               (it,ipy) &
                                             + cgrid%avg_vapor_ac                    (ipy)
         !----- These variables are already averaged at gridtype, just add them up. -------!
         do k=1,nzg
            cgrid%qmean_soil_temp (k,it,ipy)  = cgrid%qmean_soil_temp           (k,it,ipy) &
                                              + cgrid%avg_soil_temp                (k,ipy)
            cgrid%qmean_soil_water(k,it,ipy)  = cgrid%qmean_soil_water          (k,it,ipy) &
                                              + cgrid%avg_soil_water               (k,ipy)
            cgrid%qmean_soil_mstpot(k,it,ipy) = cgrid%qmean_soil_mstpot         (k,it,ipy) &
                                              + cgrid%avg_soil_mstpot              (k,ipy)
         end do

         !----- Integrate the mean sum of squares. ----------------------------------------!
         cgrid%qmsqu_gpp            (it,ipy)  = cgrid%qmsqu_gpp                 (it,ipy)   &
                                              + cgrid%avg_gpp                      (ipy)   &
                                              * cgrid%avg_gpp                      (ipy)   &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_leaf_resp      (it,ipy)  = cgrid%qmsqu_leaf_resp            (it,ipy)  &
                                              + cgrid%avg_leaf_resp                 (ipy)  &
                                              * cgrid%avg_leaf_resp                 (ipy)  &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_root_resp      (it,ipy)  = cgrid%qmsqu_root_resp            (it,ipy)  &
                                              + cgrid%avg_root_resp                 (ipy)  &
                                              * cgrid%avg_root_resp                 (ipy)  &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_plresp         (it,ipy)  = cgrid%qmsqu_plresp               (it,ipy)  &
                                              + cgrid%avg_plant_resp                (ipy)  &
                                              * cgrid%avg_plant_resp                (ipy)  &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_nep            (it,ipy)  = cgrid%qmsqu_nep                  (it,ipy)  &
                                              + ( cgrid%avg_gpp                     (ipy)  &
                                                - cgrid%avg_leaf_resp               (ipy)  &
                                                - cgrid%avg_root_resp               (ipy)  &
                                                - cgrid%avg_growth_resp             (ipy)  &
                                                - cgrid%avg_storage_resp            (ipy)  &
                                                - cgrid%avg_vleaf_resp              (ipy)  &
                                                - cgrid%avg_htroph_resp             (ipy)) &
                                              * ( cgrid%avg_gpp                     (ipy)  &
                                                - cgrid%avg_leaf_resp               (ipy)  &
                                                - cgrid%avg_root_resp               (ipy)  &
                                                - cgrid%avg_growth_resp             (ipy)  &
                                                - cgrid%avg_storage_resp            (ipy)  &
                                                - cgrid%avg_vleaf_resp              (ipy)  &
                                                - cgrid%avg_htroph_resp             (ipy)) &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_rh             (it,ipy)  = cgrid%qmsqu_rh                   (it,ipy)  &
                                              + cgrid%avg_htroph_resp               (ipy)  &
                                              * cgrid%avg_htroph_resp               (ipy)  &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_nee            (it,ipy)  = cgrid%qmsqu_nee                 (it,ipy)   &
                                              + cgrid%avg_carbon_ac                (ipy)   &
                                              * cgrid%avg_carbon_ac                (ipy)   &
                                              * umols_2_kgCyr * umols_2_kgCyr

         cgrid%qmsqu_sensible_ac    (it,ipy)  = cgrid%qmsqu_sensible_ac         (it,ipy)   &
                                              + cgrid%avg_sensible_ac              (ipy)   &
                                              * cgrid%avg_sensible_ac              (ipy)

         cgrid%qmsqu_sensible_lc    (it,ipy)  = cgrid%qmsqu_sensible_lc         (it,ipy)   &
                                              + cgrid%avg_sensible_lc              (ipy)   &
                                              * cgrid%avg_sensible_lc              (ipy)

         cgrid%qmsqu_sensible_wc    (it,ipy)  = cgrid%qmsqu_sensible_wc         (it,ipy)   &
                                              + cgrid%avg_sensible_wc              (ipy)   &
                                              * cgrid%avg_sensible_wc              (ipy)

         cgrid%qmsqu_sensible_gc    (it,ipy)  = cgrid%qmsqu_sensible_gc         (it,ipy)   &
                                              + cgrid%avg_sensible_gc              (ipy)   &
                                              * cgrid%avg_sensible_gc              (ipy)

         cgrid%qmsqu_evap           (it,ipy)  = cgrid%qmsqu_evap                (it,ipy)   &
                                              + cgrid%avg_evap                     (ipy)   &
                                              * cgrid%avg_evap                     (ipy)

         cgrid%qmsqu_transp         (it,ipy)  = cgrid%qmsqu_transp              (it,ipy)   &
                                              + cgrid%avg_transp                   (ipy)   &
                                              * cgrid%avg_transp                   (ipy)

         cgrid%qmsqu_vapor_ac    (it,ipy)     = cgrid%qmsqu_vapor_ac            (it,ipy)   &
                                              + cgrid%avg_vapor_ac                 (ipy)   &
                                              * cgrid%avg_vapor_ac                 (ipy)

         cgrid%qmsqu_vapor_lc    (it,ipy)     = cgrid%qmsqu_vapor_lc            (it,ipy)   &
                                              + cgrid%avg_vapor_lc                 (ipy)   &
                                              * cgrid%avg_vapor_lc                 (ipy)

         cgrid%qmsqu_vapor_wc    (it,ipy)     = cgrid%qmsqu_vapor_wc            (it,ipy)   &
                                              + cgrid%avg_vapor_wc                 (ipy)   &
                                              * cgrid%avg_vapor_wc                 (ipy)

         cgrid%qmsqu_vapor_gc    (it,ipy)     = cgrid%qmsqu_vapor_gc            (it,ipy)   &
                                              + cgrid%avg_vapor_gc                 (ipy)   &
                                              * cgrid%avg_vapor_gc                 (ipy)
      end if
      !------------------------------------------------------------------------------------!
   end do polyloop

   return
end subroutine integrate_ed_daily_output_flux
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will scale the daily averages of GPP and some respiration variables  !
! to normal units.  These variables are not for output, so they are done separatedly.      !
! There are also some output variables here, because these depend on the average of the    !
! gpp, and leaf and root respiration and would need to be calculated again otherwise.      !
! Some of the 5-D arrays are also integrated here for the same reason.                     !
!------------------------------------------------------------------------------------------!
subroutine normalize_ed_daily_vars(cgrid,timefac1)
   use ed_state_vars , only : edtype        & ! structure
                            , polygontype   & ! structure
                            , sitetype      & ! structure
                            , patchtype     ! ! structure
   use ed_max_dims   , only : n_pft         & ! intent(in)
                            , n_age         & ! intent(in)
                            , n_dist_types  & ! intent(in)
                            , n_dbh         ! ! intent(in)
   use ed_misc_coms  , only : imoutput      & ! intent(in)
                            , idoutput      & ! intent(in)
                            , iqoutput      & ! intent(in)
                            , ddbhi         & ! intent(in)
                            , dagei         ! ! intent(in)
   use consts_coms   , only : umols_2_kgCyr ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)                                  , target     :: cgrid
   real                                          , intent(in) :: timefac1
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)                             , pointer    :: cpoly
   type(sitetype)                                , pointer    :: csite
   type(patchtype)                               , pointer    :: cpatch
   integer                                                    :: ipy
   integer                                                    :: isi
   integer                                                    :: ipa
   integer                                                    :: ico
   integer                                                    :: ipft
   integer                                                    :: ilu
   integer                                                    :: idbh
   integer                                                    :: iage
   real                                                       :: pss_gpp
   real                                                       :: pss_leaf_resp
   real                                                       :: pss_root_resp
   real                                                       :: sss_gpp
   real                                                       :: sss_leaf_resp
   real                                                       :: sss_root_resp
   real                                                       :: poly_area_i
   real                                                       :: site_area_i
   !----- Locally saved variables. --------------------------------------------------------!
   logical           , save       :: first_time = .true.
   logical           , save       :: save_daily
   logical           , save       :: save_monthly
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      first_time   = .false.
      save_daily   = imoutput > 0 .or. idoutput > 0 .or. iqoutput > 0
      save_monthly = imoutput > 0 .or. iqoutput > 0
   end if

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      !----- This part is done only if arrays are sought. ---------------------------------!
      if (save_daily) then
         poly_area_i             = 1./sum(cpoly%area)
         sss_gpp                    = 0.
         sss_leaf_resp              = 0.
         sss_root_resp              = 0.
      end if
      
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         

         if (save_daily) then
            site_area_i               = 1./ sum(csite%area)
            pss_gpp                   = 0.
            pss_leaf_resp             = 0.
            pss_root_resp             = 0.
         end if
         
         patchloop: do ipa=1,csite%npatches

            csite%today_A_decomp (ipa) = csite%today_A_decomp(ipa)  * timefac1
            csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa) * timefac1

            !----- Copy the decomposition terms to the daily mean if they are sought. -----!
            if (save_daily) then
               csite%dmean_A_decomp(ipa)  = csite%today_A_decomp(ipa)
               csite%dmean_Af_decomp(ipa) = csite%today_Af_decomp(ipa)
               !----- Integrate the monthly mean. -----------------------------------------!
               if (save_monthly) then
                  csite%mmean_A_decomp(ipa)  = csite%mmean_A_decomp(ipa)                   &
                                             + csite%dmean_A_decomp(ipa)
                  csite%mmean_Af_decomp(ipa) = csite%mmean_Af_decomp(ipa)                  &
                                             + csite%dmean_Af_decomp(ipa)
               end if
            end if

            cpatch => csite%patch(ipa)
            
            !----- Included a loop so it won't crash with empty cohorts... ----------------!
            cohortloop: do ico=1,cpatch%ncohorts
               cpatch%today_gpp(ico)       = cpatch%today_gpp(ico)       * timefac1
               cpatch%today_gpp_pot(ico)   = cpatch%today_gpp_pot(ico)   * timefac1
               cpatch%today_gpp_max(ico)   = cpatch%today_gpp_max(ico)   * timefac1
               cpatch%today_leaf_resp(ico) = cpatch%today_leaf_resp(ico) * timefac1
               cpatch%today_root_resp(ico) = cpatch%today_root_resp(ico) * timefac1

               !---------------------------------------------------------------------------!
               !    We now update the daily means of GPP, and leaf and root respiration,   !
               ! and we convert them to kgC/plant/yr.                                      !
               !---------------------------------------------------------------------------!
               if (save_daily) then
                  cpatch%dmean_gpp(ico)       = cpatch%today_gpp(ico)                      &
                                              * umols_2_kgCyr / cpatch%nplant(ico)
                  cpatch%dmean_leaf_resp(ico) = cpatch%today_leaf_resp(ico)                &
                                              * umols_2_kgCyr / cpatch%nplant(ico)
                  cpatch%dmean_root_resp(ico) = cpatch%today_root_resp(ico)                &
                                              * umols_2_kgCyr / cpatch%nplant(ico)
                  pss_gpp                     = pss_gpp                                    &
                                              + cpatch%today_gpp(ico)                      &
                                              * csite%area(ipa)                            &
                                              * umols_2_kgCyr                           
                  pss_leaf_resp               = pss_leaf_resp                              &
                                              + cpatch%today_leaf_resp(ico)                &
                                              * csite%area(ipa)                            &
                                              * umols_2_kgCyr
                  pss_root_resp               = pss_root_resp                              &
                                              + cpatch%today_root_resp(ico)                &
                                              * csite%area(ipa)                            &
                                              * umols_2_kgCyr
               end if

               !---------------------------------------------------------------------------!
               !    We update the following monthly means here because these dmean vari-   !
               ! ables will be discarded before integrate_ed_monthly_output_vars is        !
               ! called.                                                                   !
               !---------------------------------------------------------------------------!
               if (save_monthly) then 
                  cpatch%mmean_gpp(ico)           = cpatch%mmean_gpp(ico)                  &
                                                  + cpatch%dmean_gpp(ico)
                  cpatch%mmean_leaf_resp(ico)     = cpatch%mmean_leaf_resp(ico)            &
                                                  + cpatch%dmean_leaf_resp(ico)
                  cpatch%mmean_root_resp(ico)     = cpatch%mmean_root_resp(ico)            &
                                                  + cpatch%dmean_root_resp(ico)
               end if
            end do cohortloop
         end do patchloop
         if (save_daily) then
            sss_gpp       = sss_gpp       + pss_gpp       * site_area_i * cpoly%area(isi)
            sss_leaf_resp = sss_leaf_resp + pss_leaf_resp * site_area_i * cpoly%area(isi)
            sss_root_resp = sss_root_resp + pss_root_resp * site_area_i * cpoly%area(isi)
         end if
         
      end do siteloop

      if (save_daily) then
         cgrid%dmean_gpp(ipy)           = sss_gpp          * poly_area_i
         cgrid%dmean_leaf_resp(ipy)     = sss_leaf_resp    * poly_area_i
         cgrid%dmean_root_resp(ipy)     = sss_root_resp    * poly_area_i
      end if
      
      if (save_monthly) then
         cgrid%mmean_gpp(ipy)       = cgrid%mmean_gpp(ipy)                                 &
                                    + cgrid%dmean_gpp(ipy)
         cgrid%mmean_leaf_resp(ipy) = cgrid%mmean_leaf_resp(ipy)                           &
                                    + cgrid%dmean_leaf_resp(ipy)
         cgrid%mmean_root_resp(ipy) = cgrid%mmean_root_resp(ipy)                           &
                                    + cgrid%dmean_root_resp(ipy)
      end if
   end do polyloop
   
   return
end subroutine normalize_ed_daily_vars
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This subroutine will scale the daily NPP allocation terms                            !
!------------------------------------------------------------------------------------------!
subroutine normalize_ed_dailyNPP_vars(cgrid)
   use ed_state_vars , only : edtype        & ! structure
                            , polygontype   & ! structure
                            , sitetype      & ! structure
                            , patchtype     ! ! structure
   use ed_misc_coms  , only : imoutput      & ! intent(in)
                            , idoutput      & ! intent(in)
                            , iqoutput      ! ! intent(in)
   use consts_coms   , only : yr_day ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)                                  , target     :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)                             , pointer    :: cpoly
   type(sitetype)                                , pointer    :: csite
   type(patchtype)                               , pointer    :: cpatch
   integer                                                    :: ipy
   integer                                                    :: isi
   integer                                                    :: ipa
   integer                                                    :: ico
   real                                                       :: pss_nppleaf    
   real                                                       :: pss_nppfroot   
   real                                                       :: pss_nppsapwood 
   real                                                       :: pss_nppcroot   
   real                                                       :: pss_nppseeds
   real                                                       :: pss_nppwood    
   real                                                       :: pss_nppdaily 
   real                                                       :: sss_nppleaf    
   real                                                       :: sss_nppfroot   
   real                                                       :: sss_nppsapwood 
   real                                                       :: sss_nppcroot   
   real                                                       :: sss_nppseeds
   real                                                       :: sss_nppwood    
   real                                                       :: sss_nppdaily   
   real                                                       :: poly_area_i
   real                                                       :: site_area_i
   !----- Locally saved variables. --------------------------------------------------------!
   logical           , save       :: first_time = .true.
   logical           , save       :: save_daily
   logical           , save       :: save_monthly
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      first_time   = .false.
      save_daily   = imoutput > 0 .or. idoutput > 0 .or. iqoutput > 0
      save_monthly = imoutput > 0 .or. iqoutput > 0
   end if

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      !----- This part is done only if arrays are sought. ---------------------------------!
      if (save_daily) then
         poly_area_i             = 1./sum(cpoly%area)
         sss_nppleaf                = 0.
         sss_nppfroot               = 0.
         sss_nppsapwood             = 0.
         sss_nppcroot               = 0.
         sss_nppseeds               = 0.
         sss_nppwood                = 0.
         sss_nppdaily               = 0.
      end if
      
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         

         if (save_daily) then
            site_area_i               = 1./ sum(csite%area)
            pss_nppleaf               = 0.
            pss_nppfroot              = 0.
            pss_nppsapwood            = 0.
            pss_nppcroot              = 0.
            pss_nppseeds              = 0.
            pss_nppwood               = 0.
            pss_nppdaily              = 0.
         end if
         
         patchloop: do ipa=1,csite%npatches

            cpatch => csite%patch(ipa)
            
            !----- Included a loop so it won't crash with empty cohorts... ----------------!
            cohortloop: do ico=1,cpatch%ncohorts
 
               !---------------------------------------------------------------------------!
               !    We now update the daily means of NPP allocation terms                  !
               ! and we convert them to kgC/plant/yr                                      !
               !---------------------------------------------------------------------------!
               if (save_daily) then                                              
                  cpatch%dmean_nppleaf(ico)   = cpatch%today_nppleaf(ico)                  &
                                              * yr_day / cpatch%nplant(ico)
                  cpatch%dmean_nppfroot(ico)  = cpatch%today_nppfroot(ico)                 &
                                              * yr_day / cpatch%nplant(ico)                                              
                  cpatch%dmean_nppsapwood(ico)= cpatch%today_nppsapwood(ico)               &
                                              * yr_day / cpatch%nplant(ico)                                            
                  cpatch%dmean_nppcroot(ico)  = cpatch%today_nppcroot(ico)                 &
                                              * yr_day / cpatch%nplant(ico)
                  cpatch%dmean_nppseeds(ico)  = cpatch%today_nppseeds(ico)                 &
                                              * yr_day / cpatch%nplant(ico)
                  cpatch%dmean_nppwood(ico)   = cpatch%today_nppwood(ico)                  &
                                              * yr_day / cpatch%nplant(ico)
                  cpatch%dmean_nppdaily(ico)  = cpatch%today_nppdaily(ico)                 &
                                              * yr_day / cpatch%nplant(ico)
                  pss_nppleaf                 = pss_nppleaf                                &
                                              + cpatch%today_nppleaf(ico)                  &
                                              * yr_day * csite%area(ipa)                            
                  pss_nppfroot                = pss_nppfroot                               &
                                              + cpatch%today_nppfroot(ico)                 &
                                              * yr_day * csite%area(ipa)                            
                  pss_nppsapwood              = pss_nppsapwood                             &
                                              + cpatch%today_nppsapwood(ico)               &
                                              * yr_day * csite%area(ipa)                            
                  pss_nppcroot                = pss_nppcroot                               &
                                              + cpatch%today_nppcroot(ico)                 &
                                              * yr_day * csite%area(ipa)                    
                  pss_nppseeds                = pss_nppseeds                               &
                                              + cpatch%today_nppseeds(ico)                 &
                                              * yr_day * csite%area(ipa) 
                  pss_nppwood                 = pss_nppwood                                &
                                              + cpatch%today_nppwood(ico)                  &
                                              * yr_day * csite%area(ipa)                            
                  pss_nppdaily                = pss_nppdaily                               &
                                              + cpatch%today_nppdaily(ico)                 &
                                              * yr_day * csite%area(ipa)                            
               end if

               !---------------------------------------------------------------------------!
               !    We update the following monthly means here because these dmean vari-   !
               ! ables will be discarded before integrate_ed_monthly_output_vars is        !
               ! called.                                                                   !
               !---------------------------------------------------------------------------!
               if (save_monthly) then 
                  cpatch%mmean_nppleaf(ico)       = cpatch%mmean_nppleaf(ico)              &
                                                  + cpatch%dmean_nppleaf(ico)
                  cpatch%mmean_nppfroot(ico)      = cpatch%mmean_nppfroot(ico)             &
                                                  + cpatch%dmean_nppfroot(ico)
                  cpatch%mmean_nppsapwood(ico)    = cpatch%mmean_nppsapwood(ico)           &
                                                  + cpatch%dmean_nppsapwood(ico)
                  cpatch%mmean_nppcroot(ico)      = cpatch%mmean_nppcroot(ico)             &
                                                  + cpatch%dmean_nppcroot(ico)
                  cpatch%mmean_nppseeds(ico)      = cpatch%mmean_nppseeds(ico)             &
                                                  + cpatch%dmean_nppseeds(ico)
                  cpatch%mmean_nppwood(ico)       = cpatch%mmean_nppwood(ico)              &
                                                  + cpatch%dmean_nppwood(ico)
                  cpatch%mmean_nppdaily(ico)      = cpatch%mmean_nppdaily(ico)             &
                                                  + cpatch%dmean_nppdaily(ico)
               end if
            end do cohortloop
         end do patchloop
         if (save_daily) then
            sss_nppleaf   = sss_nppleaf   + pss_nppleaf   * site_area_i * cpoly%area(isi)
            sss_nppfroot  = sss_nppfroot  + pss_nppfroot  * site_area_i * cpoly%area(isi)
            sss_nppsapwood = sss_nppsapwood  + pss_nppsapwood  * site_area_i               &
                              * cpoly%area(isi)
            sss_nppcroot  = sss_nppcroot  + pss_nppcroot  * site_area_i * cpoly%area(isi)
            sss_nppseeds  = sss_nppseeds  + pss_nppseeds  * site_area_i * cpoly%area(isi)
            sss_nppwood   = sss_nppwood   + pss_nppwood   * site_area_i * cpoly%area(isi)
            sss_nppdaily  = sss_nppdaily  + pss_nppdaily  * site_area_i * cpoly%area(isi)
         end if
         
      end do siteloop

      if (save_daily) then
         cgrid%dmean_nppleaf(ipy)       = sss_nppleaf      * poly_area_i
         cgrid%dmean_nppfroot(ipy)      = sss_nppfroot     * poly_area_i
         cgrid%dmean_nppsapwood(ipy)    = sss_nppsapwood   * poly_area_i
         cgrid%dmean_nppcroot(ipy)      = sss_nppcroot     * poly_area_i
         cgrid%dmean_nppseeds(ipy)      = sss_nppseeds     * poly_area_i
         cgrid%dmean_nppwood(ipy)       = sss_nppwood      * poly_area_i
         cgrid%dmean_nppdaily(ipy)      = sss_nppdaily     * poly_area_i
      end if
      
      if (save_monthly) then
         cgrid%mmean_nppleaf(ipy)   = cgrid%mmean_nppleaf(ipy)                             &
                                    + cgrid%dmean_nppleaf(ipy)
         cgrid%mmean_nppfroot(ipy)  = cgrid%mmean_nppfroot(ipy)                            &
                                    + cgrid%dmean_nppfroot(ipy)
         cgrid%mmean_nppsapwood(ipy)= cgrid%mmean_nppsapwood(ipy)                          &
                                    + cgrid%dmean_nppsapwood(ipy)
         cgrid%mmean_nppcroot(ipy)  = cgrid%mmean_nppcroot(ipy)                            &
                                    + cgrid%dmean_nppcroot(ipy)
         cgrid%mmean_nppseeds(ipy)  = cgrid%mmean_nppseeds(ipy)                            &
                                    + cgrid%dmean_nppseeds(ipy)
         cgrid%mmean_nppwood(ipy)   = cgrid%mmean_nppwood(ipy)                             &
                                    + cgrid%dmean_nppwood(ipy)
         cgrid%mmean_nppdaily(ipy)  = cgrid%mmean_nppdaily(ipy)                            &
                                    + cgrid%dmean_nppdaily(ipy)
      end if
   end do polyloop
   
   return
end subroutine normalize_ed_dailyNPP_vars
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
!    This subroutine normalize the sum before writing the daily analysis. It also computes !
! some of the variables that didn't need to be computed every time step, like LAI.         !
!------------------------------------------------------------------------------------------!
subroutine normalize_ed_daily_output_vars(cgrid)
   use ed_state_vars        , only : edtype        & ! structure
                                   , polygontype   & ! structure
                                   , sitetype      & ! structure
                                   , patchtype     ! ! structure
   use grid_coms            , only : nzg           ! ! intent(in)
   use ed_max_dims          , only : n_pft         & ! intent(in)
                                   , n_dbh         & ! intent(in)
                                   , n_age         & ! intent(in)
                                   , n_dist_types  ! ! intent(in)
   use consts_coms          , only : cpi           & ! intent(in)
                                   , alvl          & ! intent(in)
                                   , day_sec       & ! intent(in)
                                   , umols_2_kgCyr & ! intent(in)
                                   , yr_day        & ! intent(in)
                                   , p00i          & ! intent(in)
                                   , rocp          ! ! intent(in)
   use ed_misc_coms         , only : dtlsm         & ! intent(in)
                                   , frqsum        & ! intent(in)
                                   , ddbhi         & ! intent(in)
                                   , dagei         ! ! intent(in)
   use pft_coms             , only : init_density  ! ! intent(in)
   use therm_lib            , only : qwtk          & ! subroutine
                                   , idealdenssh   ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)                       , target     :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)                  , pointer    :: cpoly
   type(sitetype)                     , pointer    :: csite
   type(patchtype)                    , pointer    :: cpatch
   integer                                         :: ipy
   integer                                         :: isi
   integer                                         :: ipa
   integer                                         :: ico
   integer                                         :: ipft
   integer                                         :: ilu
   integer                                         :: idbh
   integer                                         :: iage
   integer                                         :: k
   logical                                         :: forest
   logical                                         :: any_resolvable
   real                                            :: poly_area_i
   real                                            :: site_area_i
   real                                            :: forest_area_i
   real                                            :: patch_laiall_i
   real                                            :: pss_fsw         , sss_fsw
   real                                            :: pss_fsn         , sss_fsn
   real                                            :: pss_fs_open     , sss_fs_open
   real                                            :: pss_storage_resp, sss_storage_resp
   real                                            :: pss_vleaf_resp  , sss_vleaf_resp
   real                                            :: pss_growth_resp , sss_growth_resp
   real                                            :: sss_rh
   real                                            :: sss_albedo
   real                                            :: sss_albedo_beam
   real                                            :: sss_albedo_diffuse
   real                                            :: veg_fliq
   real                                            :: dtlsm_o_daylight
   real                                            :: frqsum_o_daylight
   !----- Locally saved variables. --------------------------------------------------------!
   logical                            , save       :: find_factors    = .true.
   real                               , save       :: dtlsm_o_daysec  = 1.e34
   real                               , save       :: frqsum_o_daysec = 1.e34
   !---------------------------------------------------------------------------------------!


   !----- Compute the normalization factors. This is done once. ---------------------------!
   if (find_factors) then
      dtlsm_o_daysec  = dtlsm/day_sec
      frqsum_o_daysec = frqsum/day_sec
      find_factors    = .false.
   end if

   !----- Reset area indices. -------------------------------------------------------------!
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      cgrid%lai_pft            (:,ipy) = 0.
      cgrid%wpa_pft            (:,ipy) = 0.
      cgrid%wai_pft            (:,ipy) = 0.
      do isi=1,cpoly%nsites
         cpoly%lai_pft (:,isi)  = 0.
         cpoly%wpa_pft (:,isi)  = 0.
         cpoly%wai_pft (:,isi)  = 0.
      end do
   end do


   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)


      !------------------------------------------------------------------------------------!
      !    State variables, updated every time step, so these are normalized by            !
      ! dtlsm/day_sec.                                                                     !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_leaf_energy(ipy)  = cgrid%dmean_leaf_energy(ipy)  * dtlsm_o_daysec
      cgrid%dmean_leaf_hcap(ipy)    = cgrid%dmean_leaf_hcap(ipy)    * dtlsm_o_daysec
      cgrid%dmean_leaf_water(ipy)   = cgrid%dmean_leaf_water(ipy)   * dtlsm_o_daysec
      cgrid%dmean_wood_energy(ipy)  = cgrid%dmean_wood_energy(ipy)  * dtlsm_o_daysec
      cgrid%dmean_wood_hcap(ipy)    = cgrid%dmean_wood_hcap(ipy)    * dtlsm_o_daysec
      cgrid%dmean_wood_water(ipy)   = cgrid%dmean_wood_water(ipy)   * dtlsm_o_daysec
      cgrid%dmean_can_theta(ipy)    = cgrid%dmean_can_theta(ipy)    * dtlsm_o_daysec
      cgrid%dmean_can_theiv(ipy)    = cgrid%dmean_can_theiv(ipy)    * dtlsm_o_daysec
      cgrid%dmean_can_shv(ipy)      = cgrid%dmean_can_shv(ipy)      * dtlsm_o_daysec
      cgrid%dmean_can_co2(ipy)      = cgrid%dmean_can_co2(ipy)      * dtlsm_o_daysec
      cgrid%dmean_can_prss(ipy)     = cgrid%dmean_can_prss(ipy)     * dtlsm_o_daysec
      cgrid%dmean_gnd_temp(ipy)     = cgrid%dmean_gnd_temp(ipy)     * dtlsm_o_daysec
      cgrid%dmean_gnd_shv (ipy)     = cgrid%dmean_gnd_shv (ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_temp(ipy)     = cgrid%dmean_atm_temp(ipy)     * dtlsm_o_daysec
      cgrid%dmean_rshort(ipy)       = cgrid%dmean_rshort(ipy)       * dtlsm_o_daysec
      cgrid%dmean_rshort_diff(ipy)  = cgrid%dmean_rshort_diff(ipy)  * dtlsm_o_daysec
      cgrid%dmean_rlong(ipy)        = cgrid%dmean_rlong(ipy)        * dtlsm_o_daysec
      cgrid%dmean_atm_shv(ipy)      = cgrid%dmean_atm_shv(ipy)      * dtlsm_o_daysec
      cgrid%dmean_atm_co2(ipy)      = cgrid%dmean_atm_co2(ipy)      * dtlsm_o_daysec
      cgrid%dmean_atm_prss(ipy)     = cgrid%dmean_atm_prss(ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_vels(ipy)     = cgrid%dmean_atm_vels(ipy)     * dtlsm_o_daysec

      !------------------------------------------------------------------------------------!
      !     Find the canopy variables that are not conserved when pressure changes.        !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_can_temp(ipy)     = cgrid%dmean_can_theta(ipy)                           &
                                    * (p00i * cgrid%dmean_can_prss(ipy)) ** rocp
      cgrid%dmean_can_rhos(ipy)     = idealdenssh (cgrid%dmean_can_prss(ipy)               &
                                                  ,cgrid%dmean_can_temp(ipy)               &
                                                  ,cgrid%dmean_can_shv (ipy) )
      !------------------------------------------------------------------------------------!




      !----- Find the leaf temperature, only when the mean heat capacity is non-zero. -----!
      if (cgrid%dmean_leaf_hcap(ipy) > 0.) then
         call qwtk(cgrid%dmean_leaf_energy(ipy),cgrid%dmean_leaf_water(ipy)                &
                  ,cgrid%dmean_leaf_hcap(ipy),cgrid%dmean_leaf_temp(ipy),veg_fliq)
      else
         cgrid%dmean_leaf_temp(ipy) = cgrid%dmean_gnd_temp(ipy)
      end if
      !------------------------------------------------------------------------------------!




      !----- Find the leaf temperature, only when the mean heat capacity is non-zero. -----!
      if (cgrid%dmean_wood_hcap(ipy) > 0.) then
         call qwtk(cgrid%dmean_wood_energy(ipy),cgrid%dmean_wood_water(ipy)                &
                  ,cgrid%dmean_wood_hcap(ipy),cgrid%dmean_wood_temp(ipy),veg_fliq)
      else
         cgrid%dmean_wood_temp(ipy) = cgrid%dmean_gnd_temp(ipy)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     State variables, updated every frqsum, so these are normalized by              !
      ! frqsum/day_sec.                                                                    !
      !------------------------------------------------------------------------------------!
      do k=1,nzg
         cgrid%dmean_soil_temp  (k,ipy) = cgrid%dmean_soil_temp  (k,ipy) * frqsum_o_daysec
         cgrid%dmean_soil_water (k,ipy) = cgrid%dmean_soil_water (k,ipy) * frqsum_o_daysec
         cgrid%dmean_soil_mstpot(k,ipy) = cgrid%dmean_soil_mstpot(k,ipy) * frqsum_o_daysec
         cgrid%dmean_transloss  (k,ipy) = cgrid%dmean_transloss  (k,ipy) * frqsum_o_daysec
      end do
      !----- Precipitation and runoff. ----------------------------------------------------!
      cgrid%dmean_pcpg     (ipy)  = cgrid%dmean_pcpg     (ipy) * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_runoff   (ipy)  = cgrid%dmean_runoff   (ipy) * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_drainage (ipy)  = cgrid%dmean_drainage (ipy) * frqsum_o_daysec ! kg/m2/s

      !----- Vapor flux. ------------------------------------------------------------------!
      cgrid%dmean_vapor_lc(ipy)   = cgrid%dmean_vapor_lc(ipy)  * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_vapor_wc(ipy)   = cgrid%dmean_vapor_wc(ipy)  * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_vapor_gc(ipy)   = cgrid%dmean_vapor_gc(ipy)  * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_vapor_ac(ipy)   = cgrid%dmean_vapor_ac(ipy)  * frqsum_o_daysec ! kg/m2/s


      !------------------------------------------------------------------------------------!
      !     Flux variables, updated every frqsum, so these are normalized by               !
      ! frqsum/day_sec.                                                                    !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_evap        (ipy)  = cgrid%dmean_evap        (ipy)  * frqsum_o_daysec
      cgrid%dmean_transp      (ipy)  = cgrid%dmean_transp      (ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_lc (ipy)  = cgrid%dmean_sensible_lc (ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_wc (ipy)  = cgrid%dmean_sensible_wc (ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_gc (ipy)  = cgrid%dmean_sensible_gc (ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_ac (ipy)  = cgrid%dmean_sensible_ac (ipy)  * frqsum_o_daysec
      cgrid%dmean_rshort_gnd  (ipy)  = cgrid%dmean_rshort_gnd  (ipy)  * frqsum_o_daysec
      cgrid%dmean_rlong_gnd   (ipy)  = cgrid%dmean_rlong_gnd   (ipy)  * frqsum_o_daysec
      cgrid%dmean_rlongup     (ipy)  = cgrid%dmean_rlongup     (ipy)  * frqsum_o_daysec
      cgrid%dmean_rlong_albedo(ipy)  = cgrid%dmean_rlong_albedo(ipy)  * frqsum_o_daysec

      !------------------------------------------------------------------------------------!
      !      Carbon flux variables should be total flux integrated over the day at this    !
      ! point in umol/m2/s.  We just multiply by one year in seconds and convert to kgC,   !
      ! so the units will be kgC/m2/yr.                                                    !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_nee        (ipy)  = cgrid%dmean_nee      (ipy) * frqsum_o_daysec
      cgrid%dmean_plresp     (ipy)  = cgrid%dmean_plresp   (ipy) * frqsum_o_daysec
      cgrid%dmean_nep        (ipy)  = cgrid%dmean_nep      (ipy) * frqsum_o_daysec
      cgrid%dmean_gpp_dbh  (:,ipy)  = cgrid%dmean_gpp_dbh(:,ipy) * frqsum_o_daysec

      cgrid%dmean_co2_residual   (ipy) = cgrid%dmean_co2_residual   (ipy)                  &
                                       * frqsum_o_daysec
      cgrid%dmean_energy_residual(ipy) = cgrid%dmean_energy_residual(ipy) * frqsum_o_daysec
      cgrid%dmean_water_residual (ipy) = cgrid%dmean_water_residual (ipy) * frqsum_o_daysec


      poly_area_i = 1./sum(cpoly%area)
      sss_growth_resp     = 0.
      sss_storage_resp    = 0.
      sss_vleaf_resp      = 0.
      sss_rh              = 0.
      sss_albedo          = 0.
      sss_albedo_beam     = 0.
      sss_albedo_diffuse  = 0.
      sss_fsn             = 0.
      sss_fsw             = 0.
      sss_fs_open         = 0.

      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !     Find the average day length.                                                !
         !---------------------------------------------------------------------------------!
         dtlsm_o_daylight  = dtlsm  / cpoly%daylight(isi)
         frqsum_o_daylight = frqsum / cpoly%daylight(isi)
         !---------------------------------------------------------------------------------!

         cpoly%dmean_co2_residual(isi)    = cpoly%dmean_co2_residual(isi)                  &
                                          * umols_2_kgCyr * frqsum_o_daysec
         cpoly%dmean_energy_residual(isi) = cpoly%dmean_energy_residual(isi)               &
                                          * frqsum_o_daysec
         cpoly%dmean_water_residual(isi)  = cpoly%dmean_water_residual(isi)                &
                                          * frqsum_o_daysec
         
         site_area_i = 1./sum(csite%area)

         !---------------------------------------------------------------------------------!
         !     Finding the total "forest" area.  By forest we mean the fraction of land    !
         ! that is not agriculture, even if the area is not a forest.                      !
         !---------------------------------------------------------------------------------!
         forest_area_i = sum(csite%area,csite%dist_type /= 1)
         if (forest_area_i > 1.e-6) then
            forest_area_i = 1. / forest_area_i
         else
            forest_area_i = 0. ! Tiny forest area, we will neglect it in this site. 
         end if
         !---------------------------------------------------------------------------------!

         !----- Initialize auxiliary variables to add patchtype variables. ----------------!
         pss_fsn          = 0.
         pss_fsw          = 0.
         pss_fs_open      = 0.
         pss_growth_resp  = 0.
         pss_storage_resp = 0.
         pss_vleaf_resp   = 0.

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            
            
            any_resolvable = .false.
            if (cpatch%ncohorts > 0) then
               any_resolvable = any(cpatch%leaf_resolvable(1:cpatch%ncohorts))
            end if


            cohortloop: do ico=1,cpatch%ncohorts

               !---------------------------------------------------------------------------!
               !     These variables must be scaled.  They are updated every time step.    !
               !---------------------------------------------------------------------------!
               cpatch%dmean_par_l       (ico) = cpatch%dmean_par_l       (ico)             &
                                              * dtlsm_o_daysec
               cpatch%dmean_par_l_beam  (ico) = cpatch%dmean_par_l_beam  (ico)             &
                                              * dtlsm_o_daysec
               cpatch%dmean_par_l_diff  (ico) = cpatch%dmean_par_l_diff  (ico)             &
                                              * dtlsm_o_daysec
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The light level, the fraction of open stomates, and the water demand  !
               ! and supply variables are averaged over the length of day light only.  We  !
               ! find this variable only if there is any day light (this is to avoid       !
               ! problems with polar nights).                                              !
               !---------------------------------------------------------------------------!
               if (cpoly%daylight(isi) >= dtlsm) then
                  cpatch%dmean_fs_open         (ico) = cpatch%dmean_fs_open         (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_fsw             (ico) = cpatch%dmean_fsw             (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_fsn             (ico) = cpatch%dmean_fsn             (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_psi_open        (ico) = cpatch%dmean_psi_open        (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_psi_closed      (ico) = cpatch%dmean_psi_closed      (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_water_supply    (ico) = cpatch%dmean_water_supply    (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_light_level     (ico) = cpatch%dmean_light_level     (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_light_level_beam(ico) = cpatch%dmean_light_level_beam(ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_light_level_diff(ico) = cpatch%dmean_light_level_diff(ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_beamext_level   (ico) = cpatch%dmean_beamext_level   (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_diffext_level   (ico) = cpatch%dmean_diffext_level   (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_lambda_light    (ico) = cpatch%dmean_lambda_light    (ico)  &
                                                     * dtlsm_o_daylight
               else
                  cpatch%dmean_fs_open         (ico) = 0.
                  cpatch%dmean_fsw             (ico) = 0.
                  cpatch%dmean_fsn             (ico) = 0.
                  cpatch%dmean_psi_open        (ico) = 0.
                  cpatch%dmean_psi_closed      (ico) = 0.
                  cpatch%dmean_water_supply    (ico) = 0.
                  cpatch%dmean_light_level     (ico) = 0.
                  cpatch%dmean_light_level_beam(ico) = 0.
                  cpatch%dmean_light_level_diff(ico) = 0.
                  cpatch%dmean_beamext_level   (ico) = 0.
                  cpatch%dmean_diffext_level   (ico) = 0.
                  cpatch%dmean_lambda_light    (ico) = 0.
               end if
            end do cohortloop


            !------------------------------------------------------------------------------!
            !     Integrate the fraction of open stomata 
            !------------------------------------------------------------------------------!
            if (any_resolvable) then
               patch_laiall_i = 1./max(tiny(1.),sum(cpatch%lai,cpatch%leaf_resolvable))
               pss_fsn     = pss_fsn + csite%area(ipa)                                     &
                           * ( sum(cpatch%dmean_fsn * cpatch%lai,cpatch%leaf_resolvable)   &
                             * patch_laiall_i)
               pss_fsw     = pss_fsw + csite%area(ipa)                                     &
                           * ( sum(cpatch%dmean_fsw * cpatch%lai,cpatch%leaf_resolvable)   &
                             * patch_laiall_i)
               pss_fs_open = pss_fs_open + csite%area(ipa)                                 &
                           * ( sum(cpatch%dmean_fs_open * cpatch%lai                       &
                                  ,cpatch%leaf_resolvable)                                 &
                             * patch_laiall_i)
            end if

            !------------------------------------------------------------------------------!
            !     "Forest" here means non-agricultural patch, it may be a naturally occur- !
            ! ring open canopy biome.                                                      !
            !------------------------------------------------------------------------------!
            forest = csite%dist_type(ipa) /= 1

            !----- CO2 residual is now in kgC/m2/yr!!! ------------------------------------!
            csite%dmean_co2_residual(ipa)    = csite%dmean_co2_residual(ipa)               &
                                             * umols_2_kgCyr * frqsum_o_daysec
            csite%dmean_energy_residual(ipa) = csite%dmean_energy_residual(ipa)            &
                                             * frqsum_o_daysec
            csite%dmean_water_residual(ipa)  = csite%dmean_water_residual(ipa)             &
                                             * frqsum_o_daysec
            csite%dmean_rk4step(ipa)         = csite%dmean_rk4step(ipa)                    &
                                             * frqsum_o_daysec
            !------------------------------------------------------------------------------!
            !     The light level is averaged over the length of day light only.  We find  !
            ! this variable only if there is any day light (this is to avoid problems with !
            ! polar nights).                                                               !
            !------------------------------------------------------------------------------!
            if (cpoly%daylight(isi) >= dtlsm) then
               csite%dmean_lambda_light(ipa)    = csite%dmean_lambda_light(ipa)            &
                                                * dtlsm / cpoly%daylight(isi)
            else
               csite%dmean_lambda_light(ipa)    = 0.0
            end if
            !------------------------------------------------------------------------------!
            !     Heterotrophic respiration is currently the integral over a day, given    !
            ! in µmol(CO2)/m²/s, so we multiply by the number of seconds in a year and     !
            ! convert to kgC, so the final units will be kgC/m2/yr.                        !
            !------------------------------------------------------------------------------!
            csite%dmean_rh(ipa)              = csite%dmean_rh(ipa) * frqsum_o_daysec
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Albedo is updated over daylight hours only.                               !
            !------------------------------------------------------------------------------!
            csite%dmean_albedo(ipa)              = csite%dmean_rh(ipa) * frqsum_o_daylight
            csite%dmean_albedo(ipa)              = csite%dmean_rh(ipa) * frqsum_o_daylight
            csite%dmean_albedo(ipa)              = csite%dmean_rh(ipa) * frqsum_o_daylight
            !------------------------------------------------------------------------------!

            if (cpatch%ncohorts > 0) then
               pss_growth_resp  = pss_growth_resp + csite%area(ipa)                        &
                                * sum(cpatch%growth_respiration  * cpatch%nplant)          &
                                * yr_day
               pss_storage_resp = pss_storage_resp + csite%area(ipa)                       &
                                * sum(cpatch%storage_respiration * cpatch%nplant)          &
                                * yr_day
               pss_vleaf_resp   = pss_vleaf_resp   + csite%area(ipa)                       &
                                * sum(cpatch%vleaf_respiration   * cpatch%nplant)          &
                                * yr_day
               do ipft=1,n_pft
                  cpoly%lai_pft(ipft,isi)  = cpoly%lai_pft(ipft,isi)                       &
                                           + sum(cpatch%lai,cpatch%pft == ipft)            &
                                           * csite%area(ipa) * site_area_i
                  cpoly%wpa_pft(ipft,isi)  = cpoly%wpa_pft(ipft,isi)                       &
                                           + sum(cpatch%wpa,cpatch%pft == ipft)            &
                                           * csite%area(ipa) * site_area_i
                  cpoly%wai_pft(ipft,isi)  = cpoly%wai_pft(ipft,isi)                       &
                                           + sum(cpatch%wai,cpatch%pft == ipft)            &
                                           * csite%area(ipa) * site_area_i
               end do

            end if

         end do patchloop

         !----- Add this patch to the site sum. -------------------------------------------!

         sss_fsn          = sss_fsn          + (pss_fsn          * site_area_i)            &
                                             * cpoly%area(isi)
         sss_fsw          = sss_fsw          + (pss_fsw          * site_area_i)            &
                                             * cpoly%area(isi)
         sss_fs_open      = sss_fs_open      + (pss_fs_open      * site_area_i)            &
                                             * cpoly%area(isi)
         sss_growth_resp  = sss_growth_resp  + (pss_growth_resp  * site_area_i)            &
                                             * cpoly%area(isi)
         sss_storage_resp = sss_storage_resp + (pss_storage_resp * site_area_i)            &
                                             * cpoly%area(isi)
         sss_vleaf_resp   = sss_vleaf_resp   + (pss_vleaf_resp   * site_area_i)            &
                                             * cpoly%area(isi)
         sss_rh              = sss_rh                                                      &
                             + (sum(csite%dmean_rh * csite%area) * site_area_i)            &
                             * cpoly%area(isi)
         sss_albedo          = sss_albedo                                                  &
                             + (sum(csite%dmean_albedo * csite%area) * site_area_i)        &
                             * cpoly%area(isi)
         sss_albedo_beam     = sss_albedo_beam                                             &
                             + (sum(csite%dmean_albedo_beam * csite%area) * site_area_i)   &
                             * cpoly%area(isi)
         sss_albedo_diffuse  = sss_albedo_diffuse                                          &
                             + (sum(csite%dmean_albedo_diffuse * csite%area)*site_area_i)  &
                             * cpoly%area(isi)
      end do siteloop
      
      !------------------------------------------------------------------------------------!
      !     Find the area indices per PFT class.                                           !
      !------------------------------------------------------------------------------------!
      do ipft=1,n_pft
         cgrid%lai_pft(ipft,ipy)  = cgrid%lai_pft(ipft,ipy)                                &
                                  + sum(cpoly%lai_pft(ipft,:)*cpoly%area) * poly_area_i
         cgrid%wpa_pft(ipft,ipy)  = cgrid%wpa_pft(ipft,ipy)                                &
                                  + sum(cpoly%wpa_pft(ipft,:)*cpoly%area) * poly_area_i
         cgrid%wai_pft(ipft,ipy)  = cgrid%wai_pft(ipft,ipy)                                &
                                  + sum(cpoly%wai_pft(ipft,:)*cpoly%area) * poly_area_i
      end do
      !------------------------------------------------------------------------------------!

      cgrid%dmean_fsn(ipy)     = cgrid%dmean_fsn(ipy)     + sss_fsn     * poly_area_i
      cgrid%dmean_fsw(ipy)     = cgrid%dmean_fsw(ipy)     + sss_fsw     * poly_area_i
      cgrid%dmean_fs_open(ipy) = cgrid%dmean_fs_open(ipy) + sss_fs_open * poly_area_i
      
      cgrid%dmean_rh            (ipy)  = cgrid%dmean_rh            (ipy)                   &
                                       + sss_rh             * poly_area_i
      cgrid%dmean_albedo        (ipy)  = cgrid%dmean_albedo        (ipy)                   &
                                       + sss_albedo         * poly_area_i
      cgrid%dmean_albedo_beam   (ipy)  = cgrid%dmean_albedo_beam   (ipy)                   &
                                       + sss_albedo_beam    * poly_area_i
      cgrid%dmean_albedo_diffuse(ipy)  = cgrid%dmean_albedo_diffuse(ipy)                   &
                                       + sss_albedo_diffuse * poly_area_i
      cgrid%dmean_growth_resp(ipy)  = cgrid%dmean_growth_resp(ipy)                         &
                                    + sss_growth_resp  * poly_area_i
      cgrid%dmean_storage_resp(ipy) = cgrid%dmean_storage_resp(ipy)                        &
                                    + sss_storage_resp * poly_area_i
      cgrid%dmean_vleaf_resp(ipy)   = cgrid%dmean_vleaf_resp(ipy)                          &
                                    + sss_vleaf_resp   * poly_area_i

   end do polyloop

   return
 end subroutine normalize_ed_daily_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zero_ed_daily_vars(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine resets the daily_averages for variables actually used in the          !
! integration.                                                                             !
!------------------------------------------------------------------------------------------!
   use ed_state_vars        , only: edtype,polygontype,sitetype,patchtype
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy
   integer                    :: isi
   integer                    :: ipa
   integer                    :: ico
   !---------------------------------------------------------------------------------------!
   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
            
      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            
            !----- Reset variables stored in sitetype. ------------------------------------!
            csite%today_A_decomp(ipa)  = 0.0
            csite%today_Af_decomp(ipa) = 0.0

            !----- Reset variables stored in patchtype. -----------------------------------!
            do ico = 1, cpatch%ncohorts
               cpatch%today_gpp      (ico) = 0.0
               cpatch%today_nppleaf  (ico) = 0.0
               cpatch%today_nppfroot (ico) = 0.0
               cpatch%today_nppsapwood (ico) = 0.0
               cpatch%today_nppcroot (ico) = 0.0
               cpatch%today_nppseeds (ico) = 0.0
               cpatch%today_nppwood  (ico) = 0.0
               cpatch%today_nppdaily (ico) = 0.0
               cpatch%today_gpp_pot  (ico) = 0.0
               cpatch%today_gpp_max  (ico) = 0.0
               cpatch%today_leaf_resp(ico) = 0.0
               cpatch%today_root_resp(ico) = 0.0
            end do
         end do
      end do
   end do
   return
end subroutine zero_ed_daily_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine resets the daily_averages once the daily average was written and used !
! to compute the monthly mean (in case the latter was requested).                          !
!------------------------------------------------------------------------------------------!
subroutine zero_ed_daily_output_vars(cgrid)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   & ! structure
                           , sitetype      & ! structure
                           , patchtype     ! ! structure
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy
   integer                    :: isi
   integer                    :: ipa
   integer                    :: ico
   !---------------------------------------------------------------------------------------!
   

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      
      !----- Variables stored in edtype. --------------------------------------------------!
      cgrid%dmean_pcpg           (ipy) = 0.
      cgrid%dmean_drainage       (ipy) = 0.
      cgrid%dmean_runoff         (ipy) = 0.
      cgrid%dmean_vapor_lc       (ipy) = 0.
      cgrid%dmean_vapor_wc       (ipy) = 0.
      cgrid%dmean_vapor_gc       (ipy) = 0.
      cgrid%dmean_vapor_ac       (ipy) = 0.
      cgrid%dmean_gpp            (ipy) = 0.
      cgrid%dmean_nppleaf        (ipy) = 0.
      cgrid%dmean_nppfroot       (ipy) = 0.
      cgrid%dmean_nppsapwood     (ipy) = 0.
      cgrid%dmean_nppcroot       (ipy) = 0.
      cgrid%dmean_nppseeds       (ipy) = 0.
      cgrid%dmean_nppwood        (ipy) = 0.
      cgrid%dmean_nppdaily       (ipy) = 0.
      cgrid%dmean_evap           (ipy) = 0.
      cgrid%dmean_transp         (ipy) = 0.
      cgrid%dmean_sensible_lc    (ipy) = 0.
      cgrid%dmean_sensible_wc    (ipy) = 0.
      cgrid%dmean_sensible_gc    (ipy) = 0.
      cgrid%dmean_sensible_ac    (ipy) = 0.
      cgrid%dmean_rshort_gnd     (ipy) = 0.
      cgrid%dmean_rlong_gnd      (ipy) = 0.
      cgrid%dmean_rlongup        (ipy) = 0.
      cgrid%dmean_rlong_albedo   (ipy) = 0.
      cgrid%dmean_albedo         (ipy) = 0.
      cgrid%dmean_albedo_beam    (ipy) = 0.
      cgrid%dmean_albedo_diffuse (ipy) = 0.

      cgrid%dmean_nee            (ipy) = 0.
      cgrid%dmean_plresp         (ipy) = 0.
      cgrid%dmean_rh             (ipy) = 0.
      cgrid%dmean_leaf_resp      (ipy) = 0.
      cgrid%dmean_root_resp      (ipy) = 0.
      cgrid%dmean_growth_resp    (ipy) = 0.
      cgrid%dmean_storage_resp   (ipy) = 0.
      cgrid%dmean_vleaf_resp     (ipy) = 0.
      cgrid%dmean_nep            (ipy) = 0.
      cgrid%dmean_soil_temp    (:,ipy) = 0.
      cgrid%dmean_soil_water   (:,ipy) = 0.
      cgrid%dmean_soil_mstpot  (:,ipy) = 0.
      cgrid%dmean_transloss    (:,ipy) = 0.
      cgrid%dmean_gpp_dbh      (:,ipy) = 0.
      cgrid%dmean_fs_open        (ipy) = 0.
      cgrid%dmean_fsw            (ipy) = 0.
      cgrid%dmean_fsn            (ipy) = 0.
      cgrid%dmean_leaf_energy    (ipy) = 0.
      cgrid%dmean_leaf_hcap      (ipy) = 0.
      cgrid%dmean_leaf_water     (ipy) = 0.
      cgrid%dmean_leaf_temp      (ipy) = 0.
      cgrid%dmean_wood_energy    (ipy) = 0.
      cgrid%dmean_wood_hcap      (ipy) = 0.
      cgrid%dmean_wood_water     (ipy) = 0.
      cgrid%dmean_wood_temp      (ipy) = 0.
      cgrid%dmean_can_temp       (ipy) = 0.
      cgrid%dmean_can_shv        (ipy) = 0.
      cgrid%dmean_can_co2        (ipy) = 0.
      cgrid%dmean_can_rhos       (ipy) = 0.
      cgrid%dmean_can_prss       (ipy) = 0.
      cgrid%dmean_can_theta      (ipy) = 0.
      cgrid%dmean_can_theiv      (ipy) = 0.
      cgrid%dmean_gnd_temp       (ipy) = 0.
      cgrid%dmean_gnd_shv        (ipy) = 0.
      cgrid%dmean_atm_temp       (ipy) = 0.
      cgrid%dmean_rshort         (ipy) = 0.
      cgrid%dmean_rshort_diff    (ipy) = 0.
      cgrid%dmean_rlong          (ipy) = 0.
      cgrid%dmean_atm_shv        (ipy) = 0.
      cgrid%dmean_atm_co2        (ipy) = 0.
      cgrid%dmean_atm_prss       (ipy) = 0.
      cgrid%dmean_atm_vels       (ipy) = 0.
      cgrid%lai_pft            (:,ipy) = 0.
      cgrid%wpa_pft            (:,ipy) = 0.
      cgrid%wai_pft            (:,ipy) = 0.
      cgrid%dmean_co2_residual   (ipy) = 0.
      cgrid%dmean_energy_residual(ipy) = 0.
      cgrid%dmean_water_residual (ipy) = 0.

      !----- Reset variables stored in polygontype. ---------------------------------------!
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)

         cpoly%lai_pft              (:,isi) = 0.
         cpoly%wpa_pft              (:,isi) = 0.
         cpoly%wai_pft              (:,isi) = 0.
         cpoly%dmean_co2_residual     (isi) = 0.
         cpoly%dmean_energy_residual  (isi) = 0.
         cpoly%dmean_water_residual   (isi) = 0.
         cpoly%daylight               (isi) = 0.

         do ipa=1,csite%npatches
            csite%dmean_co2_residual   (ipa) = 0.
            csite%dmean_energy_residual(ipa) = 0.
            csite%dmean_water_residual (ipa) = 0.
            csite%dmean_rh             (ipa) = 0.
            csite%dmean_rk4step        (ipa) = 0.
            csite%dmean_lambda_light   (ipa) = 0.
            csite%dmean_A_decomp       (ipa) = 0.
            csite%dmean_Af_decomp      (ipa) = 0.
            csite%dmean_albedo         (ipa) = 0.
            csite%dmean_albedo_beam    (ipa) = 0.
            csite%dmean_albedo_diffuse (ipa) = 0.

            cpatch => csite%patch(ipa)
            do ico=1, cpatch%ncohorts
               cpatch%dmean_gpp(ico)              = 0.
               cpatch%dmean_nppleaf(ico)          = 0.
               cpatch%dmean_nppfroot(ico)         = 0.
               cpatch%dmean_nppsapwood(ico)       = 0.
               cpatch%dmean_nppcroot(ico)         = 0.
               cpatch%dmean_nppseeds(ico)         = 0.
               cpatch%dmean_nppwood(ico)          = 0.
               cpatch%dmean_nppdaily(ico)         = 0.
               cpatch%dmean_leaf_resp(ico)        = 0.
               cpatch%dmean_root_resp(ico)        = 0.
               cpatch%dmean_par_l(ico)            = 0.
               cpatch%dmean_par_l_beam(ico)       = 0.
               cpatch%dmean_par_l_diff(ico)       = 0.
               cpatch%dmean_fs_open(ico)          = 0.
               cpatch%dmean_fsw(ico)              = 0.
               cpatch%dmean_fsn(ico)              = 0.
               cpatch%dmean_psi_open(ico)         = 0.
               cpatch%dmean_psi_closed(ico)       = 0.
               cpatch%dmean_water_supply(ico)     = 0.
               cpatch%dmean_light_level(ico)      = 0.
               cpatch%dmean_light_level_beam(ico) = 0.
               cpatch%dmean_light_level_diff(ico) = 0.
               cpatch%dmean_beamext_level(ico)    = 0.
               cpatch%dmean_diffext_level(ico)    = 0.
               cpatch%dmean_lambda_light(ico)     = 0.
            end do
         end do
      end do
   end do

   return
end subroutine zero_ed_daily_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!                            |---------------------------------|                           !
!                            |** MONTHLY AVERAGE SUBROUTINES **|                           !
!                            |---------------------------------|                           !
!==========================================================================================!
!==========================================================================================!
!    This subroutine integrates the monthly average. This is called after the daily means  !
! were integrated and normalized.                                                          !
!------------------------------------------------------------------------------------------!
subroutine integrate_ed_monthly_output_vars(cgrid)
   use ed_state_vars, only : edtype        & ! structure
                           , polygontype   & ! structure
                           , sitetype      & ! structure
                           , patchtype     ! ! structure
   use ed_max_dims  , only : n_dbh         & ! intent(in)
                           , n_pft         & ! intent(in) 
                           , n_dist_types  & ! intent(in)
                           , n_mort        ! ! intent(in)
   use consts_coms ,  only : yr_day
   implicit none
   !----- Argument. -----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   integer                     :: ico
   integer                     :: imt
   !---------------------------------------------------------------------------------------!
   
   poly_loop: do ipy=1,cgrid%npolygons
      !------------------------------------------------------------------------------------!
      ! First the mean variables that can be computed from the daily averages              !
      !------------------------------------------------------------------------------------!
      cgrid%mmean_fs_open (ipy) = cgrid%mmean_fs_open (ipy) + cgrid%dmean_fs_open (ipy)
      cgrid%mmean_fsw     (ipy) = cgrid%mmean_fsw     (ipy) + cgrid%dmean_fsw     (ipy)
      cgrid%mmean_fsn     (ipy) = cgrid%mmean_fsn     (ipy) + cgrid%dmean_fsn     (ipy)
      cgrid%mmean_evap    (ipy) = cgrid%mmean_evap    (ipy) + cgrid%dmean_evap    (ipy)
      cgrid%mmean_transp  (ipy) = cgrid%mmean_transp  (ipy) + cgrid%dmean_transp  (ipy)

      cgrid%mmean_vapor_ac      (ipy) = cgrid%mmean_vapor_ac      (ipy)                    &
                                      + cgrid%dmean_vapor_ac      (ipy)
      cgrid%mmean_vapor_gc      (ipy) = cgrid%mmean_vapor_gc      (ipy)                    &
                                      + cgrid%dmean_vapor_gc      (ipy)
      cgrid%mmean_vapor_lc      (ipy) = cgrid%mmean_vapor_lc      (ipy)                    &
                                      + cgrid%dmean_vapor_lc      (ipy)
      cgrid%mmean_vapor_wc      (ipy) = cgrid%mmean_vapor_wc      (ipy)                    &
                                      + cgrid%dmean_vapor_wc      (ipy)
      cgrid%mmean_sensible_ac   (ipy) = cgrid%mmean_sensible_ac   (ipy)                    &
                                      + cgrid%dmean_sensible_ac   (ipy)
      cgrid%mmean_sensible_gc   (ipy) = cgrid%mmean_sensible_gc   (ipy)                    &
                                      + cgrid%dmean_sensible_gc   (ipy)
      cgrid%mmean_sensible_lc   (ipy) = cgrid%mmean_sensible_lc   (ipy)                    &
                                      + cgrid%dmean_sensible_lc   (ipy)
      cgrid%mmean_sensible_wc   (ipy) = cgrid%mmean_sensible_wc   (ipy)                    &
                                      + cgrid%dmean_sensible_wc   (ipy)
      cgrid%mmean_rshort_gnd    (ipy) = cgrid%mmean_rshort_gnd    (ipy)                    &
                                      + cgrid%dmean_rshort_gnd    (ipy)
      cgrid%mmean_rlong_gnd     (ipy) = cgrid%mmean_rlong_gnd     (ipy)                    &
                                      + cgrid%dmean_rlong_gnd     (ipy)
      cgrid%mmean_rlongup       (ipy) = cgrid%mmean_rlongup       (ipy)                    &
                                      + cgrid%dmean_rlongup       (ipy)
      cgrid%mmean_rlong_albedo  (ipy) = cgrid%mmean_rlong_albedo  (ipy)                    &
                                      + cgrid%dmean_rlong_albedo  (ipy)
      cgrid%mmean_albedo        (ipy) = cgrid%mmean_albedo        (ipy)                    &
                                      + cgrid%dmean_albedo        (ipy)
      cgrid%mmean_albedo_beam   (ipy) = cgrid%mmean_albedo_beam   (ipy)                    &
                                      + cgrid%dmean_albedo_beam   (ipy)
      cgrid%mmean_albedo_diffuse(ipy) = cgrid%mmean_albedo_diffuse(ipy)                    &
                                      + cgrid%dmean_albedo_diffuse(ipy)
      cgrid%mmean_nee           (ipy) = cgrid%mmean_nee           (ipy)                    &
                                      + cgrid%dmean_nee           (ipy)
      cgrid%mmean_nep           (ipy) = cgrid%mmean_nep           (ipy)                    &
                                      + cgrid%dmean_nep           (ipy)
      cgrid%mmean_plresp        (ipy) = cgrid%mmean_plresp        (ipy)                    &
                                      + cgrid%dmean_plresp        (ipy)
      cgrid%mmean_rh            (ipy) = cgrid%mmean_rh            (ipy)                    &
                                      + cgrid%dmean_rh            (ipy)
      cgrid%mmean_growth_resp   (ipy) = cgrid%mmean_growth_resp   (ipy)                    &
                                      + cgrid%dmean_growth_resp   (ipy)
      cgrid%mmean_storage_resp  (ipy) = cgrid%mmean_storage_resp  (ipy)                    &
                                      + cgrid%dmean_storage_resp  (ipy)
      cgrid%mmean_vleaf_resp    (ipy) = cgrid%mmean_vleaf_resp    (ipy)                    &
                                      + cgrid%dmean_vleaf_resp    (ipy)

      cgrid%mmean_soil_temp   (:,ipy) = cgrid%mmean_soil_temp   (:,ipy)                    &
                                      + cgrid%dmean_soil_temp   (:,ipy)
      cgrid%mmean_soil_water  (:,ipy) = cgrid%mmean_soil_water  (:,ipy)                    &
                                      + cgrid%dmean_soil_water  (:,ipy)
      cgrid%mmean_soil_mstpot (:,ipy) = cgrid%mmean_soil_mstpot (:,ipy)                    &
                                      + cgrid%dmean_soil_mstpot (:,ipy)
      cgrid%mmean_transloss   (:,ipy) = cgrid%mmean_transloss   (:,ipy)                    &
                                      + cgrid%dmean_transloss   (:,ipy)
      cgrid%mmean_gpp_dbh     (:,ipy) = cgrid%mmean_gpp_dbh     (:,ipy)                    &
                                      + cgrid%dmean_gpp_dbh     (:,ipy)

      cgrid%mmean_lai_pft     (:,ipy) = cgrid%mmean_lai_pft     (:,ipy)                    &
                                      + cgrid%lai_pft           (:,ipy)
      cgrid%mmean_wpa_pft     (:,ipy) = cgrid%mmean_wpa_pft     (:,ipy)                    &
                                      + cgrid%wpa_pft           (:,ipy)
      cgrid%mmean_wai_pft     (:,ipy) = cgrid%mmean_wai_pft     (:,ipy)                    &
                                      + cgrid%wai_pft           (:,ipy)

      cgrid%mmean_leaf_energy   (ipy) = cgrid%mmean_leaf_energy   (ipy)                    &
                                      + cgrid%dmean_leaf_energy   (ipy)
      cgrid%mmean_leaf_hcap     (ipy) = cgrid%mmean_leaf_hcap     (ipy)                    &
                                      + cgrid%dmean_leaf_hcap     (ipy)
      cgrid%mmean_leaf_water    (ipy) = cgrid%mmean_leaf_water    (ipy)                    & 
                                      + cgrid%dmean_leaf_water    (ipy)
      cgrid%mmean_wood_energy   (ipy) = cgrid%mmean_wood_energy   (ipy)                    &
                                      + cgrid%dmean_wood_energy   (ipy)
      cgrid%mmean_wood_hcap     (ipy) = cgrid%mmean_wood_hcap     (ipy)                    &
                                      + cgrid%dmean_wood_hcap     (ipy)
      cgrid%mmean_wood_water    (ipy) = cgrid%mmean_wood_water    (ipy)                    & 
                                      + cgrid%dmean_wood_water    (ipy)
      cgrid%mmean_can_theta     (ipy) = cgrid%mmean_can_theta     (ipy)                    &
                                      + cgrid%dmean_can_theta     (ipy)
      cgrid%mmean_can_theiv     (ipy) = cgrid%mmean_can_theiv     (ipy)                    &
                                      + cgrid%dmean_can_theiv     (ipy)
      cgrid%mmean_can_shv       (ipy) = cgrid%mmean_can_shv       (ipy)                    &
                                      + cgrid%dmean_can_shv       (ipy)
      cgrid%mmean_can_co2       (ipy) = cgrid%mmean_can_co2       (ipy)                    &
                                      + cgrid%dmean_can_co2       (ipy)
      cgrid%mmean_can_prss      (ipy) = cgrid%mmean_can_prss      (ipy)                    &
                                      + cgrid%dmean_can_prss      (ipy)
      cgrid%mmean_gnd_temp      (ipy) = cgrid%mmean_gnd_temp      (ipy)                    &
                                      + cgrid%dmean_gnd_temp      (ipy)
      cgrid%mmean_gnd_shv       (ipy) = cgrid%mmean_gnd_shv       (ipy)                    &
                                      + cgrid%dmean_gnd_shv       (ipy)
      cgrid%mmean_atm_temp      (ipy) = cgrid%mmean_atm_temp      (ipy)                    &
                                      + cgrid%dmean_atm_temp      (ipy)
      cgrid%mmean_rshort        (ipy) = cgrid%mmean_rshort        (ipy)                    &
                                      + cgrid%dmean_rshort        (ipy)
      cgrid%mmean_rshort_diff   (ipy) = cgrid%mmean_rshort_diff   (ipy)                    &
                                      + cgrid%dmean_rshort_diff   (ipy)
      cgrid%mmean_rlong         (ipy) = cgrid%mmean_rlong         (ipy)                    &
                                      + cgrid%dmean_rlong         (ipy)

      cgrid%mmean_atm_shv       (ipy) = cgrid%mmean_atm_shv       (ipy)                    &
                                      + cgrid%dmean_atm_shv       (ipy)
      cgrid%mmean_atm_co2       (ipy) = cgrid%mmean_atm_co2       (ipy)                    &
                                      + cgrid%dmean_atm_co2       (ipy)
      cgrid%mmean_atm_prss      (ipy) = cgrid%mmean_atm_prss      (ipy)                    &
                                      + cgrid%dmean_atm_prss      (ipy)
      cgrid%mmean_atm_vels      (ipy) = cgrid%mmean_atm_vels      (ipy)                    &
                                      + cgrid%dmean_atm_vels      (ipy)
      cgrid%mmean_pcpg          (ipy) = cgrid%mmean_pcpg          (ipy)                    &
                                      + cgrid%dmean_pcpg          (ipy)
      cgrid%mmean_runoff        (ipy) = cgrid%mmean_runoff        (ipy)                    &
                                      + cgrid%dmean_runoff        (ipy)
      cgrid%mmean_drainage      (ipy) = cgrid%mmean_drainage      (ipy)                    &
                                      + cgrid%dmean_drainage      (ipy)

      cgrid%mmean_co2_residual   (ipy) = cgrid%mmean_co2_residual   (ipy)                  &
                                       + cgrid%dmean_co2_residual   (ipy)
      cgrid%mmean_energy_residual(ipy) = cgrid%mmean_energy_residual(ipy)                  &
                                       + cgrid%dmean_energy_residual(ipy)
      cgrid%mmean_water_residual (ipy) = cgrid%mmean_water_residual (ipy)                  &
                                       + cgrid%dmean_water_residual (ipy)

      !------------------------------------------------------------------------------------!
      !    These are going to be the monthly mean of the sum of the squares.  We save this !
      ! instead of the standard deviation because we may want to aggregate several months  !
      ! and it is a lot simpler to do that with the sum of the squares.                    !
      !------------------------------------------------------------------------------------!
      cgrid%mmsqu_gpp        (ipy) = cgrid%mmsqu_gpp        (ipy)                          &
                                   + cgrid%dmean_gpp        (ipy)                          &
                                   * cgrid%dmean_gpp        (ipy)
      cgrid%mmsqu_leaf_resp  (ipy) = cgrid%mmsqu_leaf_resp  (ipy)                          &
                                   + cgrid%dmean_leaf_resp  (ipy)                          &
                                   * cgrid%dmean_leaf_resp  (ipy)
      cgrid%mmsqu_root_resp  (ipy) = cgrid%mmsqu_root_resp  (ipy)                          &
                                   + cgrid%dmean_root_resp  (ipy)                          &
                                   * cgrid%dmean_root_resp  (ipy)
      cgrid%mmsqu_plresp     (ipy) = cgrid%mmsqu_plresp     (ipy)                          &
                                   + cgrid%dmean_plresp     (ipy)                          &
                                   * cgrid%dmean_plresp     (ipy)
      cgrid%mmsqu_nee        (ipy) = cgrid%mmsqu_nee        (ipy)                          &
                                   + cgrid%dmean_nee        (ipy)                          &
                                   * cgrid%dmean_nee        (ipy)
      cgrid%mmsqu_nep        (ipy) = cgrid%mmsqu_nep        (ipy)                          &
                                   + cgrid%dmean_nep        (ipy)                          &
                                   * cgrid%dmean_nep        (ipy)
      cgrid%mmsqu_rh         (ipy) = cgrid%mmsqu_rh         (ipy)                          &
                                   + cgrid%dmean_rh         (ipy)                          &
                                   * cgrid%dmean_rh         (ipy)
      cgrid%mmsqu_sensible_ac(ipy) = cgrid%mmsqu_sensible_ac(ipy)                          &
                                   + cgrid%dmean_sensible_ac(ipy)                          &
                                   * cgrid%dmean_sensible_ac(ipy)
      cgrid%mmsqu_sensible_lc(ipy) = cgrid%mmsqu_sensible_lc(ipy)                          &
                                   + cgrid%dmean_sensible_lc(ipy)                          &
                                   * cgrid%dmean_sensible_lc(ipy)
      cgrid%mmsqu_sensible_wc(ipy) = cgrid%mmsqu_sensible_wc(ipy)                          &
                                   + cgrid%dmean_sensible_wc(ipy)                          &
                                   * cgrid%dmean_sensible_wc(ipy)
      cgrid%mmsqu_sensible_gc(ipy) = cgrid%mmsqu_sensible_gc(ipy)                          &
                                   + cgrid%dmean_sensible_gc(ipy)                          &
                                   * cgrid%dmean_sensible_gc(ipy)
      cgrid%mmsqu_evap       (ipy) = cgrid%mmsqu_evap       (ipy)                          &
                                   + cgrid%dmean_evap       (ipy)                          &
                                   * cgrid%dmean_evap       (ipy)
      cgrid%mmsqu_transp     (ipy) = cgrid%mmsqu_transp     (ipy)                          &
                                   + cgrid%dmean_transp     (ipy)                          &
                                   * cgrid%dmean_transp     (ipy)
      cgrid%mmsqu_vapor_ac(ipy)    = cgrid%mmsqu_vapor_ac   (ipy)                          &
                                   + cgrid%dmean_vapor_ac   (ipy)                          &
                                   * cgrid%dmean_vapor_ac   (ipy)
      cgrid%mmsqu_vapor_lc(ipy)    = cgrid%mmsqu_vapor_lc   (ipy)                          &
                                   + cgrid%dmean_vapor_lc   (ipy)                          &
                                   * cgrid%dmean_vapor_lc   (ipy)
      cgrid%mmsqu_vapor_wc(ipy)    = cgrid%mmsqu_vapor_wc   (ipy)                          &
                                   + cgrid%dmean_vapor_wc   (ipy)                          &
                                   * cgrid%dmean_vapor_wc   (ipy)
      cgrid%mmsqu_vapor_gc(ipy)    = cgrid%mmsqu_vapor_gc   (ipy)                          &
                                   + cgrid%dmean_vapor_gc   (ipy)                          &
                                   * cgrid%dmean_vapor_gc   (ipy)

      cpoly => cgrid%polygon(ipy)
      site_loop: do isi=1,cpoly%nsites
         cpoly%mmean_co2_residual(isi)    = cpoly%mmean_co2_residual(isi)                  &
                                          + cpoly%dmean_co2_residual(isi)
         cpoly%mmean_energy_residual(isi) = cpoly%mmean_energy_residual(isi)               &
                                          + cpoly%dmean_energy_residual(isi)
         cpoly%mmean_water_residual(isi)  = cpoly%mmean_water_residual(isi)                &
                                          + cpoly%dmean_water_residual(isi)

         csite => cpoly%site(isi)
         patch_loop: do ipa=1,csite%npatches
            csite%mmean_co2_residual(ipa)    = csite%mmean_co2_residual(ipa)               &
                                             + csite%dmean_co2_residual(ipa)
            csite%mmean_energy_residual(ipa) = csite%mmean_energy_residual(ipa)            &
                                             + csite%dmean_energy_residual(ipa)
            csite%mmean_water_residual(ipa)  = csite%mmean_water_residual(ipa)             &
                                             + csite%dmean_water_residual(ipa)

            csite%mmean_rh            (ipa)  = csite%mmean_rh            (ipa)             &
                                             + csite%dmean_rh            (ipa)
            csite%mmean_albedo        (ipa)  = csite%mmean_albedo        (ipa)             &
                                             + csite%dmean_albedo        (ipa)
            csite%mmean_albedo_beam   (ipa)  = csite%mmean_albedo_beam   (ipa)             &
                                             + csite%dmean_albedo_beam   (ipa)
            csite%mmean_albedo_diffuse(ipa)  = csite%mmean_albedo_diffuse(ipa)             &
                                             + csite%dmean_albedo_diffuse(ipa)

            csite%mmean_rk4step(ipa)         = csite%mmean_rk4step(ipa)                    &
                                             + csite%dmean_rk4step(ipa)

            csite%mmean_lambda_light(ipa)    = csite%mmean_lambda_light(ipa)               &
                                             + csite%dmean_lambda_light(ipa)

            cpatch => csite%patch(ipa)
            cohort_loop: do ico=1,cpatch%ncohorts
               cpatch%mmean_fs_open     (ico) = cpatch%mmean_fs_open     (ico)             &
                                              + cpatch%dmean_fs_open     (ico)
               cpatch%mmean_fsw         (ico) = cpatch%mmean_fsw         (ico)             &
                                              + cpatch%dmean_fsw         (ico)
               cpatch%mmean_fsn         (ico) = cpatch%mmean_fsn         (ico)             &
                                              + cpatch%dmean_fsn         (ico)
               cpatch%mmean_psi_open    (ico) = cpatch%mmean_psi_open    (ico)             &
                                              + cpatch%dmean_psi_open    (ico)
               cpatch%mmean_psi_closed  (ico) = cpatch%mmean_psi_closed  (ico)             &
                                              + cpatch%dmean_psi_closed  (ico)
               cpatch%mmean_water_supply(ico) = cpatch%mmean_water_supply(ico)             &
                                              + cpatch%dmean_water_supply(ico)
               cpatch%mmean_par_l       (ico) = cpatch%mmean_par_l       (ico)             &
                                              + cpatch%dmean_par_l       (ico)
               cpatch%mmean_par_l_beam  (ico) = cpatch%mmean_par_l_beam  (ico)             &
                                              + cpatch%dmean_par_l_beam  (ico)
               cpatch%mmean_par_l_diff  (ico) = cpatch%mmean_par_l_diff  (ico)             &
                                              + cpatch%dmean_par_l_diff  (ico)

               !---------------------------------------------------------------------------!
               !     The following variables are all converted to kgC/plant/yr.            !
               !---------------------------------------------------------------------------!
               cpatch%mmean_leaf_maintenance(ico) = cpatch%mmean_leaf_maintenance(ico)     &
                                                  + cpatch%leaf_maintenance(ico)           &
                                                  * yr_day
               cpatch%mmean_root_maintenance(ico) = cpatch%mmean_root_maintenance(ico)     &
                                                  + cpatch%root_maintenance(ico)           &
                                                  * yr_day
               cpatch%mmean_leaf_drop (ico)       = cpatch%leaf_drop(ico)                  &
                                                  + cpatch%leaf_drop(ico)                  &
                                                  * yr_day
               cpatch%mmean_growth_resp(ico)      = cpatch%mmean_growth_resp(ico)          &
                                                  + cpatch%growth_respiration(ico)         &
                                                  * yr_day
               cpatch%mmean_storage_resp(ico)     = cpatch%mmean_storage_resp(ico)         &
                                                  + cpatch%storage_respiration(ico)        &
                                                  * yr_day
               cpatch%mmean_vleaf_resp(ico)       = cpatch%mmean_vleaf_resp(ico)           &
                                                  + cpatch%vleaf_respiration(ico)          &
                                                  * yr_day
               !---------------------------------------------------------------------------!
               !    Light level, a simple average now.  We currently ignore that different !
               ! days have different day light lenghts.                                    !
               !---------------------------------------------------------------------------!
               cpatch%mmean_light_level(ico)      = cpatch%mmean_light_level(ico)          &
                                                  + cpatch%dmean_light_level(ico)
               cpatch%mmean_light_level_beam(ico) = cpatch%mmean_light_level_beam(ico)     &
                                                  + cpatch%dmean_light_level_beam(ico)
               cpatch%mmean_light_level_diff(ico) = cpatch%mmean_light_level_diff(ico)     &
                                                  + cpatch%dmean_light_level_diff(ico)
               cpatch%mmean_beamext_level(ico)    = cpatch%mmean_beamext_level(ico)        &
                                                  + cpatch%dmean_beamext_level(ico)
               cpatch%mmean_diffext_level(ico)    = cpatch%mmean_diffext_level(ico)        &
                                                  + cpatch%dmean_diffext_level(ico)
               cpatch%mmean_lambda_light(ico)     = cpatch%mmean_lambda_light(ico)         &
                                                  + cpatch%dmean_lambda_light(ico)

               !----- Mortality rates. ----------------------------------------------------!
               do imt=1,n_mort
                  cpatch%mmean_mort_rate(imt,ico) = cpatch%mmean_mort_rate(imt,ico)        &
                                                  + cpatch%mort_rate(imt,ico)
               end do

            end do cohort_loop

         end do patch_loop
      end do site_loop
   end do poly_loop

   return
end subroutine integrate_ed_monthly_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine normalize the sum before writing the mobthly analysis. It also        !
! computes some of the variables that didn't need to be computed every day, like AGB.      !
!------------------------------------------------------------------------------------------!
subroutine normalize_ed_monthly_output_vars(cgrid)
   use ed_state_vars        , only : edtype        & ! structure
                                   , polygontype   & ! structure
                                   , sitetype      & ! structure
                                   , patchtype     ! ! structure
   use ed_misc_coms         , only : current_time  & ! intent(in)
                                   , simtime       & ! intent(in)
                                   , ddbhi         & ! intent(in)
                                   , dagei         & ! intent(in)
                                   , iqoutput      & ! intent(in)
                                   , dtlsm         & ! intent(in)
                                   , frqfast       & ! intent(in)
                                   , ndcycle       ! ! intent(in)
   use ed_max_dims          , only : n_pft         & ! intent(in)
                                   , n_dbh         & ! intent(in)
                                   , n_age         & ! intent(in)
                                   , n_dist_types  & ! intent(in)
                                   , n_mort        ! ! intent(in)
   use consts_coms          , only : p00i          & ! intent(in)
                                   , rocp          & ! intent(in)
                                   , pio4          & ! intent(in)
                                   , umol_2_kgC    & ! intent(in)
                                   , umols_2_kgCyr & ! intent(in)
                                   , day_sec       & ! intent(in)
                                   , yr_day        ! ! intent(in)
   use pft_coms             , only : init_density  ! ! intent(in)
   use therm_lib            , only : idealdenssh   & ! function
                                   , qwtk          ! ! function
   use allometry            , only : ed_biomass    ! ! function

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)                          , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)                     , pointer :: cpoly
   type(sitetype)                        , pointer :: csite
   type(patchtype)                       , pointer :: cpatch
   type(simtime)                                   :: lastmonth
   real                                            :: ndaysi
   real                                            :: poly_area_i
   real                                            :: forest_area_i
   real                                            :: site_area_i
   real, dimension(n_pft)                          :: pss_bseeds_pft
   real, dimension(n_pft)                          :: sss_bseeds_pft
   real, dimension(n_pft,n_dbh)                    :: pss_bseeds
   real, dimension(n_pft,n_dbh)                    :: pss_pldens
   real                                            :: sss_rh
   real                                            :: sss_fsw
   real                                            :: sss_fsn
   real                                            :: sss_fs_open
   real                                            :: pss_fsw
   real                                            :: pss_fsn
   real                                            :: pss_fs_open
   real                                            :: patch_laiall_i
   integer                                         :: lmon
   integer                                         :: ipy
   integer                                         :: isi
   integer                                         :: ipa
   integer                                         :: ico
   integer                                         :: imt
   integer                                         :: t
   integer                                         :: ipft
   integer                                         :: idbh
   integer                                         :: iage
   integer                                         :: ilu
   integer                                         :: jlu
   logical                                         :: any_resolvable
   logical                                         :: forest
   real                                            :: veg_fliq
   real                                            :: cohort_seeds
   !----- Locally saved variables. --------------------------------------------------------!
   logical                            , save       :: find_factors    = .true.
   real                               , save       :: dtlsm_o_frqfast = 1.e34
   !---------------------------------------------------------------------------------------!



   !----- Compute the normalisation factors.  This is done once. --------------------------!
   if (find_factors) then
      dtlsm_o_frqfast = dtlsm  / frqfast
      find_factors    = .false.
   end if
   !---------------------------------------------------------------------------------------!


  
   !---------------------------------------------------------------------------------------!
   !     Finding the inverse of number of days used for this monthly integral.             !
   !---------------------------------------------------------------------------------------!
   call lastmonthdate(current_time,lastmonth,ndaysi)
   lmon = lastmonth%month

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      !------------------------------------------------------------------------------------!
      !      First normalize the variables previously defined.                             !
      !------------------------------------------------------------------------------------!
      cgrid%mmean_fs_open        (ipy) = cgrid%mmean_fs_open        (ipy) * ndaysi
      cgrid%mmean_fsw            (ipy) = cgrid%mmean_fsw            (ipy) * ndaysi
      cgrid%mmean_fsn            (ipy) = cgrid%mmean_fsn            (ipy) * ndaysi
      cgrid%mmean_gpp            (ipy) = cgrid%mmean_gpp            (ipy) * ndaysi
      cgrid%mmean_nppleaf        (ipy) = cgrid%mmean_nppleaf        (ipy) * ndaysi
      cgrid%mmean_nppfroot       (ipy) = cgrid%mmean_nppfroot       (ipy) * ndaysi
      cgrid%mmean_nppsapwood     (ipy) = cgrid%mmean_nppsapwood     (ipy) * ndaysi
      cgrid%mmean_nppcroot       (ipy) = cgrid%mmean_nppcroot       (ipy) * ndaysi
      cgrid%mmean_nppseeds       (ipy) = cgrid%mmean_nppseeds       (ipy) * ndaysi
      cgrid%mmean_nppwood        (ipy) = cgrid%mmean_nppwood        (ipy) * ndaysi
      cgrid%mmean_nppdaily       (ipy) = cgrid%mmean_nppdaily       (ipy) * ndaysi
      cgrid%mmean_evap           (ipy) = cgrid%mmean_evap           (ipy) * ndaysi
      cgrid%mmean_transp         (ipy) = cgrid%mmean_transp         (ipy) * ndaysi
      cgrid%mmean_vapor_ac       (ipy) = cgrid%mmean_vapor_ac       (ipy) * ndaysi
      cgrid%mmean_vapor_gc       (ipy) = cgrid%mmean_vapor_gc       (ipy) * ndaysi
      cgrid%mmean_vapor_lc       (ipy) = cgrid%mmean_vapor_lc       (ipy) * ndaysi
      cgrid%mmean_vapor_wc       (ipy) = cgrid%mmean_vapor_wc       (ipy) * ndaysi
      cgrid%mmean_sensible_ac    (ipy) = cgrid%mmean_sensible_ac    (ipy) * ndaysi
      cgrid%mmean_sensible_gc    (ipy) = cgrid%mmean_sensible_gc    (ipy) * ndaysi
      cgrid%mmean_sensible_lc    (ipy) = cgrid%mmean_sensible_lc    (ipy) * ndaysi
      cgrid%mmean_sensible_wc    (ipy) = cgrid%mmean_sensible_wc    (ipy) * ndaysi
      cgrid%mmean_rshort_gnd     (ipy) = cgrid%mmean_rshort_gnd     (ipy) * ndaysi
      cgrid%mmean_rlong_gnd      (ipy) = cgrid%mmean_rlong_gnd      (ipy) * ndaysi
      cgrid%mmean_rlongup        (ipy) = cgrid%mmean_rlongup        (ipy) * ndaysi
      cgrid%mmean_rlong_albedo   (ipy) = cgrid%mmean_rlong_albedo   (ipy) * ndaysi
      cgrid%mmean_albedo         (ipy) = cgrid%mmean_albedo         (ipy) * ndaysi
      cgrid%mmean_albedo_beam    (ipy) = cgrid%mmean_albedo_beam    (ipy) * ndaysi
      cgrid%mmean_albedo_diffuse (ipy) = cgrid%mmean_albedo_diffuse (ipy) * ndaysi
      cgrid%mmean_nee            (ipy) = cgrid%mmean_nee            (ipy) * ndaysi
      cgrid%mmean_nep            (ipy) = cgrid%mmean_nep            (ipy) * ndaysi
      cgrid%mmean_plresp         (ipy) = cgrid%mmean_plresp         (ipy) * ndaysi
      cgrid%mmean_rh             (ipy) = cgrid%mmean_rh             (ipy) * ndaysi
      cgrid%mmean_leaf_resp      (ipy) = cgrid%mmean_leaf_resp      (ipy) * ndaysi
      cgrid%mmean_root_resp      (ipy) = cgrid%mmean_root_resp      (ipy) * ndaysi
      cgrid%mmean_growth_resp    (ipy) = cgrid%mmean_growth_resp    (ipy) * ndaysi
      cgrid%mmean_storage_resp   (ipy) = cgrid%mmean_storage_resp   (ipy) * ndaysi
      cgrid%mmean_vleaf_resp     (ipy) = cgrid%mmean_vleaf_resp     (ipy) * ndaysi
      cgrid%mmean_soil_temp    (:,ipy) = cgrid%mmean_soil_temp    (:,ipy) * ndaysi
      cgrid%mmean_soil_water   (:,ipy) = cgrid%mmean_soil_water   (:,ipy) * ndaysi
      cgrid%mmean_soil_mstpot  (:,ipy) = cgrid%mmean_soil_mstpot  (:,ipy) * ndaysi
      cgrid%mmean_transloss    (:,ipy) = cgrid%mmean_transloss    (:,ipy) * ndaysi
      cgrid%mmean_gpp_dbh      (:,ipy) = cgrid%mmean_gpp_dbh      (:,ipy) * ndaysi
      cgrid%mmean_leaf_energy    (ipy) = cgrid%mmean_leaf_energy    (ipy) * ndaysi
      cgrid%mmean_leaf_hcap      (ipy) = cgrid%mmean_leaf_hcap      (ipy) * ndaysi
      cgrid%mmean_leaf_water     (ipy) = cgrid%mmean_leaf_water     (ipy) * ndaysi
      cgrid%mmean_wood_energy    (ipy) = cgrid%mmean_wood_energy    (ipy) * ndaysi
      cgrid%mmean_wood_hcap      (ipy) = cgrid%mmean_wood_hcap      (ipy) * ndaysi
      cgrid%mmean_wood_water     (ipy) = cgrid%mmean_wood_water     (ipy) * ndaysi
      cgrid%mmean_can_theta      (ipy) = cgrid%mmean_can_theta      (ipy) * ndaysi
      cgrid%mmean_can_theiv      (ipy) = cgrid%mmean_can_theiv      (ipy) * ndaysi
      cgrid%mmean_can_shv        (ipy) = cgrid%mmean_can_shv        (ipy) * ndaysi
      cgrid%mmean_can_co2        (ipy) = cgrid%mmean_can_co2        (ipy) * ndaysi
      cgrid%mmean_can_prss       (ipy) = cgrid%mmean_can_prss       (ipy) * ndaysi
      cgrid%mmean_gnd_temp       (ipy) = cgrid%mmean_gnd_temp       (ipy) * ndaysi
      cgrid%mmean_gnd_shv        (ipy) = cgrid%mmean_gnd_shv        (ipy) * ndaysi
      cgrid%mmean_atm_temp       (ipy) = cgrid%mmean_atm_temp       (ipy) * ndaysi
      cgrid%mmean_rshort         (ipy) = cgrid%mmean_rshort         (ipy) * ndaysi
      cgrid%mmean_rshort_diff    (ipy) = cgrid%mmean_rshort_diff    (ipy) * ndaysi
      cgrid%mmean_rlong          (ipy) = cgrid%mmean_rlong          (ipy) * ndaysi
      cgrid%mmean_atm_shv        (ipy) = cgrid%mmean_atm_shv        (ipy) * ndaysi
      cgrid%mmean_atm_co2        (ipy) = cgrid%mmean_atm_co2        (ipy) * ndaysi
      cgrid%mmean_atm_prss       (ipy) = cgrid%mmean_atm_prss       (ipy) * ndaysi
      cgrid%mmean_atm_vels       (ipy) = cgrid%mmean_atm_vels       (ipy) * ndaysi
      cgrid%mmean_pcpg           (ipy) = cgrid%mmean_pcpg           (ipy) * ndaysi
      cgrid%mmean_runoff         (ipy) = cgrid%mmean_runoff         (ipy) * ndaysi
      cgrid%mmean_drainage       (ipy) = cgrid%mmean_drainage       (ipy) * ndaysi
      cgrid%mmean_lai_pft      (:,ipy) = cgrid%mmean_lai_pft      (:,ipy) * ndaysi
      cgrid%mmean_wpa_pft      (:,ipy) = cgrid%mmean_wpa_pft      (:,ipy) * ndaysi
      cgrid%mmean_wai_pft      (:,ipy) = cgrid%mmean_wai_pft      (:,ipy) * ndaysi

      cgrid%mmean_co2_residual(ipy)    = cgrid%mmean_co2_residual(ipy)    * ndaysi
      cgrid%mmean_energy_residual(ipy) = cgrid%mmean_energy_residual(ipy) * ndaysi
      cgrid%mmean_water_residual(ipy)  = cgrid%mmean_water_residual(ipy)  * ndaysi

      !------------------------------------------------------------------------------------!
      !   Here we normalise the sum of squares, so it is becomes the mean.  In case you    !
      ! want to find the standard deviation, you may want to use the following equation:   !
      ! X    = variable.                                                                   !
      ! Xm   = monthly mean of the variable.                                               !
      ! Xmsq = monthly mean of the sum of the squares of this variable.                    !
      !                                                                                    !
      !              __________________          ________________                          !
      !             / SUM_i[X_i - Xm]²          / (Xmsq  - Xm²)                            !
      ! sigma = \  /  ----------------   =  \  / ---------------                           !
      !          \/       N - 1              \/      1 - 1/N                               !
      !                                                                                    !
      !------------------------------------------------------------------------------------!
      cgrid%mmsqu_gpp         (ipy) = cgrid%mmsqu_gpp         (ipy) * ndaysi
      cgrid%mmsqu_leaf_resp   (ipy) = cgrid%mmsqu_leaf_resp   (ipy) * ndaysi
      cgrid%mmsqu_root_resp   (ipy) = cgrid%mmsqu_root_resp   (ipy) * ndaysi
      cgrid%mmsqu_plresp      (ipy) = cgrid%mmsqu_plresp      (ipy) * ndaysi
      cgrid%mmsqu_nee         (ipy) = cgrid%mmsqu_nee         (ipy) * ndaysi
      cgrid%mmsqu_nep         (ipy) = cgrid%mmsqu_nep         (ipy) * ndaysi
      cgrid%mmsqu_rh          (ipy) = cgrid%mmsqu_rh          (ipy) * ndaysi
      cgrid%mmsqu_sensible_ac (ipy) = cgrid%mmsqu_sensible_ac (ipy) * ndaysi
      cgrid%mmsqu_sensible_lc (ipy) = cgrid%mmsqu_sensible_lc (ipy) * ndaysi
      cgrid%mmsqu_sensible_wc (ipy) = cgrid%mmsqu_sensible_wc (ipy) * ndaysi
      cgrid%mmsqu_sensible_gc (ipy) = cgrid%mmsqu_sensible_gc (ipy) * ndaysi
      cgrid%mmsqu_evap        (ipy) = cgrid%mmsqu_evap        (ipy) * ndaysi
      cgrid%mmsqu_transp      (ipy) = cgrid%mmsqu_transp      (ipy) * ndaysi
      cgrid%mmsqu_vapor_ac    (ipy) = cgrid%mmsqu_vapor_ac    (ipy) * ndaysi
      cgrid%mmsqu_vapor_lc    (ipy) = cgrid%mmsqu_vapor_lc    (ipy) * ndaysi
      cgrid%mmsqu_vapor_wc    (ipy) = cgrid%mmsqu_vapor_wc    (ipy) * ndaysi
      cgrid%mmsqu_vapor_gc    (ipy) = cgrid%mmsqu_vapor_gc    (ipy) * ndaysi
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Mean canopy air properties.                                                  !
      !------------------------------------------------------------------------------------!
      cgrid%mmean_can_temp    (ipy) = cgrid%mmean_can_theta(ipy)                           &
                                    * (p00i * cgrid%mmean_can_prss(ipy)) ** rocp
      cgrid%mmean_can_rhos    (ipy) = idealdenssh (cgrid%mmean_can_prss(ipy)               &
                                                  ,cgrid%mmean_can_temp(ipy)               &
                                                  ,cgrid%mmean_can_shv (ipy) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the average leaf and stem temperatures.  These are calculated only when  !
      ! there is some heat storage.                                                        !
      !------------------------------------------------------------------------------------!
      if (cgrid%mmean_leaf_hcap(ipy) > 0.) then
         call qwtk(cgrid%mmean_leaf_energy(ipy),cgrid%mmean_leaf_water(ipy)                &
                  ,cgrid%mmean_leaf_hcap(ipy),cgrid%mmean_leaf_temp(ipy),veg_fliq)
      else
         cgrid%mmean_leaf_temp(ipy) = cgrid%mmean_can_temp(ipy)
      end if
      if (cgrid%mmean_wood_hcap(ipy) > 0.) then
         call qwtk(cgrid%mmean_wood_energy(ipy),cgrid%mmean_wood_water(ipy)                &
                  ,cgrid%mmean_wood_hcap(ipy),cgrid%mmean_wood_temp(ipy),veg_fliq)
      else
         cgrid%mmean_wood_temp(ipy) = cgrid%mmean_can_temp(ipy)
      end if
      !------------------------------------------------------------------------------------!



      !---- Find AGB and basal area per PFT -----------------------------------------------!
      poly_area_i = 1./sum(cpoly%area)

      do ipft = 1,n_pft
        do idbh =1,n_dbh
          cgrid%agb_pft(ipft,ipy) = cgrid%agb_pft(ipft,ipy)                                &
                                  + sum(cpoly%agb(ipft,idbh,:)*cpoly%area)*poly_area_i
          cgrid%ba_pft(ipft,ipy)  = cgrid%ba_pft(ipft,ipy)                                 &
                                  + sum(cpoly%basal_area(ipft,idbh,:)*cpoly%area)          &
                                  * poly_area_i
        end do
      end do

      !----- Finding disturbance rates per source and target land use types. --------------!
      do ilu = 1,n_dist_types
         do jlu = 1,n_dist_types
          cgrid%disturbance_rates(ilu,jlu,ipy) = cgrid%disturbance_rates(ilu,jlu,ipy)      &
                                               + sum( cpoly%disturbance_rates(ilu,jlu,:)   &
                                                    * cpoly%area)                          &
                                               * poly_area_i
         end do
      end do


      !------------------------------------------------------------------------------------!
      !    Find a few other variables that are either updated every month, or that         !
      ! depend on site-/patch-/cohort- level variables.                                    !
      !------------------------------------------------------------------------------------!
      !----- Flush the PFT, LU, AGE, and Size (DBH) variables to zero. --------------------!
      cgrid%bseeds_pft (  :,ipy) = 0.
      cgrid%bseeds     (:,:,ipy) = 0.
      cgrid%pldens     (:,:,ipy) = 0.
      sss_bseeds_pft         (:) = 0.

      !----- Looping over all sites. ------------------------------------------------------!
      siteloop: do isi = 1, cpoly%nsites
         csite => cpoly%site(isi)

         !----- Finding the polygon-level monthly mean for residuals. ---------------------!
         cpoly%mmean_co2_residual(isi)    = cpoly%mmean_co2_residual(isi)    * ndaysi
         cpoly%mmean_energy_residual(isi) = cpoly%mmean_energy_residual(isi) * ndaysi
         cpoly%mmean_water_residual(isi)  = cpoly%mmean_water_residual(isi)  * ndaysi


         site_area_i = 1./sum(csite%area)
         !---------------------------------------------------------------------------------!
         !     Finding the total "forest" area.  By forest we mean the fraction of land    !
         ! that is not agriculture, even if the area is not a forest.                      !
         !---------------------------------------------------------------------------------!
         forest_area_i = sum(csite%area,csite%dist_type /= 1)
         if (forest_area_i > 1.e-6) then
            forest_area_i = 1. / forest_area_i
         else
            forest_area_i = 0. ! Tiny forest area, we will neglect it in this site. 
         end if
         !---------------------------------------------------------------------------------!
         
         !----- Flushing all site-level variables to zero before integrating site. --------!
         cpoly%pldens       (:,:,isi) = 0.
         cpoly%bseeds       (:,:,isi) = 0.

         !----- Flushing all patch-level variables to zero before integrating patch. ------!
         pss_bseeds_pft       (:) = 0.
         pss_bseeds         (:,:) = 0.
         pss_pldens         (:,:) = 0.

         patchloop: do ipa=1,csite%npatches

            !------------------------------------------------------------------------------!
            !    Residual of fast-scale budgets.  We hope that this is tiny...             !
            !------------------------------------------------------------------------------!
            csite%mmean_co2_residual(ipa)    = csite%mmean_co2_residual(ipa)    * ndaysi
            csite%mmean_energy_residual(ipa) = csite%mmean_energy_residual(ipa) * ndaysi
            csite%mmean_water_residual(ipa)  = csite%mmean_water_residual(ipa)  * ndaysi
            csite%mmean_rh(ipa)              = csite%mmean_rh(ipa)              * ndaysi
            csite%mmean_rk4step(ipa)         = csite%mmean_rk4step(ipa)         * ndaysi
            csite%mmean_lambda_light(ipa)    = csite%mmean_lambda_light(ipa)    * ndaysi
            csite%mmean_A_decomp(ipa)        = csite%mmean_A_decomp(ipa)        * ndaysi
            csite%mmean_Af_decomp(ipa)       = csite%mmean_Af_decomp(ipa)       * ndaysi
            csite%mmean_albedo(ipa)          = csite%mmean_albedo(ipa)          * ndaysi
            csite%mmean_albedo_beam(ipa)     = csite%mmean_albedo_beam(ipa)     * ndaysi
            csite%mmean_albedo_diffuse(ipa)  = csite%mmean_albedo_diffuse(ipa)  * ndaysi

            !------------------------------------------------------------------------------!
            !     Determining whether this is an agricultural patch or not.  Age and size  !
            ! distribution is done only for primary and secondary vegetation.              !
            !------------------------------------------------------------------------------!
            forest = csite%dist_type(ipa) /= 1

            cpatch => csite%patch(ipa)
            cohortloop: do ico = 1, cpatch%ncohorts
               !----- Find the carbon fluxes. ---------------------------------------------!
               cpatch%mmean_gpp         (ico) = cpatch%mmean_gpp         (ico) * ndaysi
               cpatch%mmean_nppleaf     (ico) = cpatch%mmean_nppleaf     (ico) * ndaysi
               cpatch%mmean_nppfroot    (ico) = cpatch%mmean_nppfroot    (ico) * ndaysi
               cpatch%mmean_nppsapwood  (ico) = cpatch%mmean_nppsapwood  (ico) * ndaysi
               cpatch%mmean_nppcroot    (ico) = cpatch%mmean_nppcroot    (ico) * ndaysi
               cpatch%mmean_nppseeds    (ico) = cpatch%mmean_nppseeds    (ico) * ndaysi
               cpatch%mmean_nppwood     (ico) = cpatch%mmean_nppwood     (ico) * ndaysi
               cpatch%mmean_nppdaily    (ico) = cpatch%mmean_nppdaily    (ico) * ndaysi
               cpatch%mmean_leaf_resp   (ico) = cpatch%mmean_leaf_resp   (ico) * ndaysi
               cpatch%mmean_root_resp   (ico) = cpatch%mmean_root_resp   (ico) * ndaysi
               cpatch%mmean_growth_resp (ico) = cpatch%mmean_growth_resp (ico) * ndaysi
               cpatch%mmean_storage_resp(ico) = cpatch%mmean_storage_resp(ico) * ndaysi
               cpatch%mmean_vleaf_resp  (ico) = cpatch%mmean_vleaf_resp  (ico) * ndaysi
               cpatch%mmean_fsw         (ico) = cpatch%mmean_fsw         (ico) * ndaysi
               cpatch%mmean_fsn         (ico) = cpatch%mmean_fsn         (ico) * ndaysi
               cpatch%mmean_fs_open     (ico) = cpatch%mmean_fs_open     (ico) * ndaysi
               cpatch%mmean_psi_open    (ico) = cpatch%mmean_psi_open    (ico) * ndaysi
               cpatch%mmean_psi_closed  (ico) = cpatch%mmean_psi_closed  (ico) * ndaysi
               cpatch%mmean_water_supply(ico) = cpatch%mmean_water_supply(ico) * ndaysi
               cpatch%mmean_par_l       (ico) = cpatch%mmean_par_l       (ico) * ndaysi
               cpatch%mmean_par_l_beam  (ico) = cpatch%mmean_par_l_beam  (ico) * ndaysi
               cpatch%mmean_par_l_diff  (ico) = cpatch%mmean_par_l_diff  (ico) * ndaysi
               cpatch%mmean_leaf_maintenance (ico) = cpatch%mmean_leaf_maintenance(ico)    &
                                                   * ndaysi
               cpatch%mmean_root_maintenance (ico) = cpatch%mmean_root_maintenance(ico)    &
                                                   * ndaysi
               cpatch%mmean_leaf_drop   (ico) = cpatch%mmean_leaf_drop   (ico) * ndaysi
               !----- Mean carbon balance is re-scaled so it will be in kgC/plant/yr. -----!
               cpatch%mmean_cb          (ico) = cpatch%mmean_cb(ico) * ndaysi * yr_day

               !----- Find the mortality rates. -------------------------------------------!
               do imt=1,n_mort
                  cpatch%mmean_mort_rate(imt,ico) = cpatch%mmean_mort_rate(imt,ico)*ndaysi
               end do

               !----- Finding the light level, ignoring changes in day time length... -----!
               cpatch%mmean_light_level (ico)      = cpatch%mmean_light_level(ico)         &
                                                   * ndaysi
               cpatch%mmean_light_level_beam (ico) = cpatch%mmean_light_level_beam(ico)    &
                                                   * ndaysi
               cpatch%mmean_light_level_diff (ico) = cpatch%mmean_light_level_diff(ico)    &
                                                   * ndaysi
               cpatch%mmean_beamext_level (ico)    = cpatch%mmean_beamext_level(ico)       &
                                                   * ndaysi
               cpatch%mmean_diffext_level (ico)    = cpatch%mmean_diffext_level(ico)       &
                                                   * ndaysi
               cpatch%mmean_lambda_light(ico)      = cpatch%mmean_lambda_light(ico)        &
                                                   * ndaysi

               !----- Define to which PFT this cohort belongs. ----------------------------!
               ipft = cpatch%pft(ico)

               !----- Computing the total seed mass of this cohort. -----------------------!
               cohort_seeds   = cpatch%nplant(ico) * cpatch%bseeds(ico)
               
               pss_bseeds_pft(ipft) = pss_bseeds_pft(ipft)                                 &
                                         + cohort_seeds * csite%area(ipa)

               if (forest) then
                  !----- Define to which size (DBH) class this cohort belongs. ------------!
                  idbh = max(1,min(n_dbh,ceiling(cpatch%dbh(ico)*ddbhi)))

                  !----- Increment the plant density. -------------------------------------!
                  pss_pldens(ipft,idbh) = pss_pldens(ipft,idbh)                            &
                                             + cpatch%nplant(ico) * csite%area(ipa)
                  pss_bseeds(ipft,idbh) = pss_bseeds(ipft,idbh)                            &
                                             + cohort_seeds * csite%area(ipa)
               end if

            end do cohortloop
         end do patchloop
         !---------------------------------------------------------------------------------!
         !     We now increment the site-level variables.                                  !
         !---------------------------------------------------------------------------------!
         !----- PFT classes. --------------------------------------------------------------!
         do ipft = 1,n_pft
            sss_bseeds_pft(ipft) = sss_bseeds_pft(ipft)                                    &
                                 + pss_bseeds_pft(ipft) * site_area_i                      &
                                 * cpoly%area(isi)
         end do
         !----- Size (DBH) classes. -------------------------------------------------------!
         do idbh = 1,n_dbh
            do ipft=1,n_pft
               cpoly%pldens(ipft,idbh,isi) = cpoly%pldens(ipft,idbh,isi)                   &
                                           + pss_pldens(ipft,idbh) * forest_area_i
               cpoly%bseeds(ipft,idbh,isi) = cpoly%bseeds(ipft,idbh,isi)                   &
                                           + pss_bseeds(ipft,idbh) * forest_area_i
            end do
         end do
         !---------------------------------------------------------------------------------!
      end do siteloop

      !------------------------------------------------------------------------------------!
      !   Incrementing the polygon-level variables.                                        !
      !------------------------------------------------------------------------------------!
      !----- PFT classes. -----------------------------------------------------------------!
      do ipft = 1,n_pft
         cgrid%bseeds_pft(ipft,ipy) = cgrid%bseeds_pft(ipft,ipy)                           &
                                    + sss_bseeds_pft(ipft) * poly_area_i
      end do
      !----- Size (DBH) classes. ----------------------------------------------------------!
      do isi=1,cpoly%nsites
         do idbh = 1,n_dbh
            do ipft=1,n_pft
               cgrid%pldens(ipft,idbh,ipy) = cgrid%pldens(ipft,idbh,ipy)                   &
                                           + cpoly%pldens(ipft,idbh,isi)                   &
                                           * cpoly%area(isi) * poly_area_i
               cgrid%bseeds(ipft,idbh,ipy) = cgrid%bseeds(ipft,idbh,ipy)                   &
                                           + cpoly%bseeds(ipft,idbh,isi)                   &
                                           * cpoly%area(isi) * poly_area_i
            end do
         end do
      end do


      !------------------------------------------------------------------------------------!
      !      The mean diurnal cycle is normalised here.                                    !
      !------------------------------------------------------------------------------------!
      if (iqoutput > 0) then
         do t=1,ndcycle

            !----- Initialise site sums (auxiliary variables). ----------------------------!
            poly_area_i = 1./sum(cpoly%area)
            sss_rh           = 0.
            sss_fsn          = 0.
            sss_fsw          = 0.
            sss_fs_open      = 0.
            !------------------------------------------------------------------------------!

            do isi=1,cpoly%nsites
               csite => cpoly%site(isi)

               !----- Initialise patch sums (auxiliary variables). ------------------------!
               pss_fsn          = 0.
               pss_fsw          = 0.
               pss_fs_open      = 0.
               !---------------------------------------------------------------------------!

               do ipa = 1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !----- Find whether there is at least one cohort that is solved. --------!
                  any_resolvable = .false.
                  if (cpatch%ncohorts > 0) then
                     any_resolvable = any(cpatch%leaf_resolvable(1:cpatch%ncohorts))
                  end if
                  !------------------------------------------------------------------------!


                  do ico=1,cpatch%ncohorts
                     !----- Convert GPP and plant respiration to kgC/plant/year. ----------!
                     cpatch%qmean_gpp         (t,ico) = cpatch%qmean_gpp           (t,ico) &
                                                      * ndaysi
                     cpatch%qmean_leaf_resp   (t,ico) = cpatch%qmean_leaf_resp     (t,ico) &
                                                      * ndaysi
                     cpatch%qmean_root_resp   (t,ico) = cpatch%qmean_root_resp     (t,ico) &
                                                      * ndaysi
                     cpatch%qmean_par_l       (t,ico) = cpatch%qmean_par_l         (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_par_l_beam  (t,ico) = cpatch%qmean_par_l_beam    (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_par_l_diff  (t,ico) = cpatch%qmean_par_l_diff    (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_fs_open     (t,ico) = cpatch%qmean_fs_open       (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_fsw         (t,ico) = cpatch%qmean_fsw           (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_fsn         (t,ico) = cpatch%qmean_fsn           (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_psi_open    (t,ico) = cpatch%qmean_psi_open      (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_psi_closed  (t,ico) = cpatch%qmean_psi_closed    (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                     cpatch%qmean_water_supply(t,ico) = cpatch%qmean_water_supply  (t,ico) &
                                                      * ndaysi * dtlsm_o_frqfast
                  end do

                  !------------------------------------------------------------------------!
                  !     Integrate the fraction of open stomata.                            !
                  !------------------------------------------------------------------------!
                  if (any_resolvable) then
                     patch_laiall_i = 1.0                                                  &
                                    / max(tiny(1.),sum(cpatch%lai,cpatch%leaf_resolvable))
                     pss_fsn        = pss_fsn + csite%area(ipa)                            &
                                    * (sum( cpatch%qmean_fsn(t,:) * cpatch%lai             &
                                          , cpatch%leaf_resolvable) * patch_laiall_i)
                     pss_fsw        = pss_fsw + csite%area(ipa)                            &
                                    * (sum( cpatch%qmean_fsw(t,:) * cpatch%lai             &
                                          , cpatch%leaf_resolvable) * patch_laiall_i)
                     pss_fs_open    = pss_fs_open + csite%area(ipa)                        &
                                    * (sum( cpatch%qmean_fs_open(t,:) * cpatch%lai         &
                                          , cpatch%leaf_resolvable) * patch_laiall_i)
                  end if
                  !------------------------------------------------------------------------!

                  csite%qmean_rh            (t,ipa) = csite%qmean_rh            (t,ipa)    &
                                                    * ndaysi
                  csite%qmean_albedo        (t,ipa) = csite%qmean_albedo        (t,ipa)    &
                                                    * ndaysi
                  csite%qmean_albedo_beam   (t,ipa) = csite%qmean_albedo_beam   (t,ipa)    &
                                                    * ndaysi
                  csite%qmean_albedo_diffuse(t,ipa) = csite%qmean_albedo_diffuse(t,ipa)    &
                                                    * ndaysi
               end do

               !----- Add this patch to the site sum. -------------------------------------!
               sss_fsn          = sss_fsn          + (pss_fsn          * site_area_i)      &
                                                   * cpoly%area(isi)
               sss_fsw          = sss_fsw          + (pss_fsw          * site_area_i)      &
                                                   * cpoly%area(isi)
               sss_fs_open      = sss_fs_open      + (pss_fs_open      * site_area_i)      &
                                                   * cpoly%area(isi)
               sss_rh           = sss_rh + ( sum(csite%qmean_rh(t,:) * csite%area)         &
                                           * site_area_i) * cpoly%area(isi)
               !---------------------------------------------------------------------------!

            end do

            cgrid%qmean_fsn           (t,ipy) = cgrid%qmean_fsn           (t,ipy)          &
                                              + sss_fsn     * poly_area_i
            cgrid%qmean_fsw           (t,ipy) = cgrid%qmean_fsw           (t,ipy)          &
                                              + sss_fsw     * poly_area_i
            cgrid%qmean_fs_open       (t,ipy) = cgrid%qmean_fs_open       (t,ipy)          &
                                              + sss_fs_open * poly_area_i
            cgrid%qmean_rh            (t,ipy) = cgrid%qmean_rh            (t,ipy)          &
                                              + sss_rh      * poly_area_i      

            cgrid%qmean_leaf_energy   (t,ipy) = cgrid%qmean_leaf_energy    (t,ipy)         &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_leaf_water    (t,ipy) = cgrid%qmean_leaf_water     (t,ipy)         &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_leaf_hcap     (t,ipy) = cgrid%qmean_leaf_hcap      (t,ipy)         &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_wood_energy   (t,ipy) = cgrid%qmean_wood_energy    (t,ipy)         &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_wood_water    (t,ipy) = cgrid%qmean_wood_water     (t,ipy)         &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_wood_hcap     (t,ipy) = cgrid%qmean_wood_hcap      (t,ipy)         &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_can_theta     (t,ipy) = cgrid%qmean_can_theta     (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_can_theiv     (t,ipy) = cgrid%qmean_can_theiv     (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_can_shv       (t,ipy) = cgrid%qmean_can_shv       (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_can_co2       (t,ipy) = cgrid%qmean_can_co2       (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_can_prss      (t,ipy) = cgrid%qmean_can_prss      (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_gnd_temp      (t,ipy) = cgrid%qmean_gnd_temp      (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_gnd_shv       (t,ipy) = cgrid%qmean_gnd_shv       (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_atm_temp      (t,ipy) = cgrid%qmean_atm_temp      (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_rshort        (t,ipy) = cgrid%qmean_rshort        (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_rshort_diff   (t,ipy) = cgrid%qmean_rshort_diff   (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_rlong         (t,ipy) = cgrid%qmean_rlong         (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_atm_shv       (t,ipy) = cgrid%qmean_atm_shv       (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_atm_co2       (t,ipy) = cgrid%qmean_atm_co2       (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_atm_prss      (t,ipy) = cgrid%qmean_atm_prss      (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_atm_vels      (t,ipy) = cgrid%qmean_atm_vels      (t,ipy)          &
                                              * ndaysi * dtlsm_o_frqfast
            cgrid%qmean_gpp           (t,ipy) = cgrid%qmean_gpp           (t,ipy)  * ndaysi
            cgrid%qmean_leaf_resp     (t,ipy) = cgrid%qmean_leaf_resp     (t,ipy)  * ndaysi
            cgrid%qmean_root_resp     (t,ipy) = cgrid%qmean_root_resp     (t,ipy)  * ndaysi
            cgrid%qmean_plresp        (t,ipy) = cgrid%qmean_plresp        (t,ipy)  * ndaysi
            cgrid%qmean_nep           (t,ipy) = cgrid%qmean_nep           (t,ipy)  * ndaysi
            cgrid%qmean_rh            (t,ipy) = cgrid%qmean_rh            (t,ipy)  * ndaysi
            cgrid%qmean_sensible_lc   (t,ipy) = cgrid%qmean_sensible_lc   (t,ipy)  * ndaysi
            cgrid%qmean_sensible_wc   (t,ipy) = cgrid%qmean_sensible_wc   (t,ipy)  * ndaysi
            cgrid%qmean_sensible_gc   (t,ipy) = cgrid%qmean_sensible_gc   (t,ipy)  * ndaysi
            cgrid%qmean_sensible_ac   (t,ipy) = cgrid%qmean_sensible_ac   (t,ipy)  * ndaysi
            cgrid%qmean_rshort_gnd    (t,ipy) = cgrid%qmean_rshort_gnd    (t,ipy)  * ndaysi
            cgrid%qmean_rlong_gnd     (t,ipy) = cgrid%qmean_rlong_gnd     (t,ipy)  * ndaysi
            cgrid%qmean_rlongup       (t,ipy) = cgrid%qmean_rlongup       (t,ipy)  * ndaysi
            cgrid%qmean_rlong_albedo  (t,ipy) = cgrid%qmean_rlong_albedo  (t,ipy)  * ndaysi
            cgrid%qmean_albedo        (t,ipy) = cgrid%qmean_albedo        (t,ipy)  * ndaysi
            cgrid%qmean_albedo_beam   (t,ipy) = cgrid%qmean_albedo_beam   (t,ipy)  * ndaysi
            cgrid%qmean_albedo_diffuse(t,ipy) = cgrid%qmean_albedo_diffuse(t,ipy)  * ndaysi
            cgrid%qmean_nee           (t,ipy) = cgrid%qmean_nee           (t,ipy)  * ndaysi
            cgrid%qmean_pcpg          (t,ipy) = cgrid%qmean_pcpg          (t,ipy)  * ndaysi
            cgrid%qmean_evap          (t,ipy) = cgrid%qmean_evap          (t,ipy)  * ndaysi
            cgrid%qmean_transp        (t,ipy) = cgrid%qmean_transp        (t,ipy)  * ndaysi
            cgrid%qmean_runoff        (t,ipy) = cgrid%qmean_runoff        (t,ipy)  * ndaysi
            cgrid%qmean_drainage      (t,ipy) = cgrid%qmean_drainage      (t,ipy)  * ndaysi
            cgrid%qmean_vapor_lc      (t,ipy) = cgrid%qmean_vapor_lc      (t,ipy)  * ndaysi
            cgrid%qmean_vapor_wc      (t,ipy) = cgrid%qmean_vapor_wc      (t,ipy)  * ndaysi
            cgrid%qmean_vapor_gc      (t,ipy) = cgrid%qmean_vapor_gc      (t,ipy)  * ndaysi
            cgrid%qmean_vapor_ac      (t,ipy) = cgrid%qmean_vapor_ac      (t,ipy)  * ndaysi
            cgrid%qmean_soil_temp   (:,t,ipy) = cgrid%qmean_soil_temp   (:,t,ipy)  * ndaysi
            cgrid%qmean_soil_water  (:,t,ipy) = cgrid%qmean_soil_water  (:,t,ipy)  * ndaysi
            cgrid%qmean_soil_mstpot (:,t,ipy) = cgrid%qmean_soil_mstpot (:,t,ipy)  * ndaysi
            cgrid%qmsqu_gpp           (t,ipy) = cgrid%qmsqu_gpp           (t,ipy)  * ndaysi
            cgrid%qmsqu_leaf_resp     (t,ipy) = cgrid%qmsqu_leaf_resp     (t,ipy)  * ndaysi
            cgrid%qmsqu_root_resp     (t,ipy) = cgrid%qmsqu_root_resp     (t,ipy)  * ndaysi
            cgrid%qmsqu_plresp        (t,ipy) = cgrid%qmsqu_plresp        (t,ipy)  * ndaysi
            cgrid%qmsqu_nee           (t,ipy) = cgrid%qmsqu_nee           (t,ipy)  * ndaysi
            cgrid%qmsqu_nep           (t,ipy) = cgrid%qmsqu_nep           (t,ipy)  * ndaysi
            cgrid%qmsqu_rh            (t,ipy) = cgrid%qmsqu_rh            (t,ipy)  * ndaysi
            cgrid%qmsqu_sensible_ac   (t,ipy) = cgrid%qmsqu_sensible_ac   (t,ipy)  * ndaysi
            cgrid%qmsqu_sensible_lc   (t,ipy) = cgrid%qmsqu_sensible_lc   (t,ipy)  * ndaysi
            cgrid%qmsqu_sensible_wc   (t,ipy) = cgrid%qmsqu_sensible_wc   (t,ipy)  * ndaysi
            cgrid%qmsqu_sensible_gc   (t,ipy) = cgrid%qmsqu_sensible_gc   (t,ipy)  * ndaysi
            cgrid%qmsqu_evap          (t,ipy) = cgrid%qmsqu_evap          (t,ipy)  * ndaysi
            cgrid%qmsqu_transp        (t,ipy) = cgrid%qmsqu_transp        (t,ipy)  * ndaysi
            cgrid%qmsqu_vapor_ac      (t,ipy) = cgrid%qmsqu_vapor_ac      (t,ipy)  * ndaysi
            cgrid%qmsqu_vapor_lc      (t,ipy) = cgrid%qmsqu_vapor_lc      (t,ipy)  * ndaysi
            cgrid%qmsqu_vapor_wc      (t,ipy) = cgrid%qmsqu_vapor_wc      (t,ipy)  * ndaysi
            cgrid%qmsqu_vapor_gc      (t,ipy) = cgrid%qmsqu_vapor_gc      (t,ipy)  * ndaysi

            !------------------------------------------------------------------------------!
            !     Find the derived average propertiesof the canopy air space.              !
            !------------------------------------------------------------------------------!
            cgrid%qmean_can_temp (t,ipy) = cgrid%qmean_can_theta(t,ipy)                    &
                                         *  (p00i * cgrid%qmean_can_prss(t,ipy)) ** rocp
            cgrid%qmean_can_rhos (t,ipy) = idealdenssh (cgrid%qmean_can_prss(t,ipy)        &
                                                       ,cgrid%qmean_can_temp(t,ipy)        &
                                                       ,cgrid%qmean_can_shv (t,ipy))


            !------------------------------------------------------------------------------!
            !     Find the mean leaf and stem temperatures, provided that the mean heat    !
            ! capacity is not zero.  In case the heat capacity is zero, use the mean       !
            ! canopy air space temperature.                                                !
            !------------------------------------------------------------------------------!
            if (cgrid%qmean_leaf_hcap(t,ipy) > 0.) then
               call qwtk(cgrid%qmean_leaf_energy(t,ipy),cgrid%qmean_leaf_water(t,ipy)      &
                        ,cgrid%qmean_leaf_hcap(t,ipy),cgrid%qmean_leaf_temp(t,ipy)         &
                        ,veg_fliq)
            else
               cgrid%qmean_leaf_temp(t,ipy) = cgrid%qmean_can_temp (t,ipy)
            end if
            if (cgrid%qmean_wood_hcap(t,ipy) > 0.) then
               call qwtk(cgrid%qmean_wood_energy(t,ipy),cgrid%qmean_wood_water(t,ipy)      &
                        ,cgrid%qmean_wood_hcap(t,ipy),cgrid%qmean_wood_temp(t,ipy)         &
                        ,veg_fliq)
            else
               cgrid%qmean_wood_temp(t,ipy) = cgrid%qmean_can_temp (t,ipy)
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

   end do polyloop

   return
end subroutine normalize_ed_monthly_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zero_ed_monthly_output_vars(cgrid)
   use ed_state_vars , only : edtype        & ! structure
                            , polygontype   & ! structure
                            , sitetype      & ! structure
                            , patchtype     ! ! structure
   use ed_misc_coms  , only : iqoutput      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy
   integer                    :: isi
   integer                    :: ipa
   integer                    :: ico
   !---------------------------------------------------------------------------------------!

   !----- The loop is necessary for coupled runs (when npolygons may be 0) ----------------!
   do ipy=1,cgrid%npolygons
      cgrid%mmean_fs_open            (ipy) = 0.
      cgrid%mmean_fsw                (ipy) = 0.
      cgrid%mmean_fsn                (ipy) = 0.
      cgrid%mmean_gpp                (ipy) = 0.
      cgrid%mmean_nppleaf            (ipy) = 0.
      cgrid%mmean_nppfroot           (ipy) = 0.
      cgrid%mmean_nppsapwood         (ipy) = 0.
      cgrid%mmean_nppcroot           (ipy) = 0.
      cgrid%mmean_nppseeds           (ipy) = 0.
      cgrid%mmean_nppwood            (ipy) = 0.
      cgrid%mmean_nppdaily           (ipy) = 0.
      cgrid%mmean_evap               (ipy) = 0.
      cgrid%mmean_transp             (ipy) = 0.
      cgrid%mmean_vapor_ac           (ipy) = 0.
      cgrid%mmean_vapor_gc           (ipy) = 0.
      cgrid%mmean_vapor_lc           (ipy) = 0.
      cgrid%mmean_vapor_wc           (ipy) = 0.
      cgrid%mmean_sensible_ac        (ipy) = 0.
      cgrid%mmean_sensible_gc        (ipy) = 0.
      cgrid%mmean_sensible_lc        (ipy) = 0.
      cgrid%mmean_sensible_wc        (ipy) = 0.
      cgrid%mmean_rshort_gnd         (ipy) = 0.
      cgrid%mmean_rlong_gnd          (ipy) = 0.
      cgrid%mmean_rlongup            (ipy) = 0.
      cgrid%mmean_rlong_albedo       (ipy) = 0.
      cgrid%mmean_albedo             (ipy) = 0.
      cgrid%mmean_albedo_beam        (ipy) = 0.
      cgrid%mmean_albedo_diffuse     (ipy) = 0.
      cgrid%mmean_nee                (ipy) = 0.
      cgrid%mmean_nep                (ipy) = 0.
      cgrid%mmean_plresp             (ipy) = 0.
      cgrid%mmean_rh                 (ipy) = 0.
      cgrid%mmean_leaf_resp          (ipy) = 0.
      cgrid%mmean_root_resp          (ipy) = 0.
      cgrid%mmean_growth_resp        (ipy) = 0.
      cgrid%mmean_storage_resp       (ipy) = 0.
      cgrid%mmean_vleaf_resp         (ipy) = 0.
      cgrid%mmean_soil_temp        (:,ipy) = 0.
      cgrid%mmean_soil_water       (:,ipy) = 0.
      cgrid%mmean_soil_mstpot      (:,ipy) = 0.
      cgrid%mmean_transloss        (:,ipy) = 0.
      cgrid%mmean_gpp_dbh          (:,ipy) = 0.
      cgrid%mmean_leaf_energy        (ipy) = 0.
      cgrid%mmean_leaf_hcap          (ipy) = 0.
      cgrid%mmean_leaf_water         (ipy) = 0.
      cgrid%mmean_leaf_temp          (ipy) = 0.
      cgrid%mmean_wood_energy        (ipy) = 0.
      cgrid%mmean_wood_hcap          (ipy) = 0.
      cgrid%mmean_wood_water         (ipy) = 0.
      cgrid%mmean_wood_temp          (ipy) = 0.
      cgrid%mmean_can_theta          (ipy) = 0.
      cgrid%mmean_can_theiv          (ipy) = 0.
      cgrid%mmean_can_prss           (ipy) = 0.
      cgrid%mmean_can_temp           (ipy) = 0.
      cgrid%mmean_can_shv            (ipy) = 0.
      cgrid%mmean_can_co2            (ipy) = 0.
      cgrid%mmean_can_rhos           (ipy) = 0.
      cgrid%mmean_gnd_temp           (ipy) = 0.
      cgrid%mmean_gnd_shv            (ipy) = 0.
      cgrid%mmean_atm_temp           (ipy) = 0.
      cgrid%mmean_rshort             (ipy) = 0.
      cgrid%mmean_rshort_diff        (ipy) = 0.
      cgrid%mmean_rlong              (ipy) = 0.
      cgrid%mmean_atm_shv            (ipy) = 0.
      cgrid%mmean_atm_co2            (ipy) = 0.
      cgrid%mmean_atm_prss           (ipy) = 0.
      cgrid%mmean_atm_vels           (ipy) = 0.
      cgrid%mmean_pcpg               (ipy) = 0.
      cgrid%mmean_runoff             (ipy) = 0.
      cgrid%mmean_drainage           (ipy) = 0.
      cgrid%mmean_lai_pft          (:,ipy) = 0.
      cgrid%mmean_wpa_pft          (:,ipy) = 0.
      cgrid%mmean_wai_pft          (:,ipy) = 0.
      cgrid%agb_pft                (:,ipy) = 0.
      cgrid%ba_pft                 (:,ipy) = 0.
      cgrid%mmsqu_gpp                (ipy) = 0.
      cgrid%mmsqu_leaf_resp          (ipy) = 0.
      cgrid%mmsqu_root_resp          (ipy) = 0.
      cgrid%mmsqu_plresp             (ipy) = 0.
      cgrid%mmsqu_nee                (ipy) = 0.
      cgrid%mmsqu_nep                (ipy) = 0.
      cgrid%mmsqu_rh                 (ipy) = 0.
      cgrid%mmsqu_sensible_ac        (ipy) = 0.
      cgrid%mmsqu_sensible_lc        (ipy) = 0.
      cgrid%mmsqu_sensible_wc        (ipy) = 0.
      cgrid%mmsqu_sensible_gc        (ipy) = 0.
      cgrid%mmsqu_evap               (ipy) = 0.
      cgrid%mmsqu_transp             (ipy) = 0.
      cgrid%mmsqu_vapor_ac           (ipy) = 0.
      cgrid%mmsqu_vapor_lc           (ipy) = 0.
      cgrid%mmsqu_vapor_wc           (ipy) = 0.
      cgrid%mmsqu_vapor_gc           (ipy) = 0.
      cgrid%disturbance_rates    (:,:,ipy) = 0.

      cgrid%mmean_co2_residual       (ipy) = 0.
      cgrid%mmean_energy_residual    (ipy) = 0.
      cgrid%mmean_water_residual     (ipy) = 0.

      cpoly => cgrid%polygon(ipy)
      do isi = 1, cpoly%nsites

         cpoly%mmean_co2_residual    (isi) = 0.
         cpoly%mmean_energy_residual (isi) = 0.
         cpoly%mmean_water_residual  (isi) = 0.

         csite => cpoly%site(isi)
         do ipa=1,csite%npatches
            csite%mmean_co2_residual      (ipa) = 0.
            csite%mmean_energy_residual   (ipa) = 0.
            csite%mmean_water_residual    (ipa) = 0.
            csite%mmean_rh                (ipa) = 0.
            csite%mmean_rk4step           (ipa) = 0.
            csite%mmean_lambda_light      (ipa) = 0.
            csite%mmean_A_decomp          (ipa) = 0.
            csite%mmean_Af_decomp         (ipa) = 0.
            csite%mmean_albedo            (ipa) = 0.
            csite%mmean_albedo_beam       (ipa) = 0.
            csite%mmean_albedo_diffuse    (ipa) = 0.

            cpatch=> csite%patch(ipa)
            do ico=1,cpatch%ncohorts
               cpatch%mmean_par_l             (ico) = 0.
               cpatch%mmean_par_l_beam        (ico) = 0.
               cpatch%mmean_par_l_diff        (ico) = 0.
               cpatch%mmean_fs_open           (ico) = 0.
               cpatch%mmean_fsw               (ico) = 0.
               cpatch%mmean_fsn               (ico) = 0.
               cpatch%mmean_psi_open          (ico) = 0.
               cpatch%mmean_psi_closed        (ico) = 0.
               cpatch%mmean_water_supply      (ico) = 0.
               cpatch%mmean_leaf_maintenance  (ico) = 0.
               cpatch%mmean_root_maintenance  (ico) = 0.
               cpatch%mmean_leaf_drop         (ico) = 0.
               cpatch%mmean_gpp               (ico) = 0.
               cpatch%mmean_nppleaf           (ico) = 0.
               cpatch%mmean_nppfroot          (ico) = 0.
               cpatch%mmean_nppsapwood        (ico) = 0.
               cpatch%mmean_nppcroot          (ico) = 0.
               cpatch%mmean_nppseeds          (ico) = 0.
               cpatch%mmean_nppwood           (ico) = 0.
               cpatch%mmean_nppdaily          (ico) = 0.
               cpatch%mmean_leaf_resp         (ico) = 0.
               cpatch%mmean_root_resp         (ico) = 0.
               cpatch%mmean_growth_resp       (ico) = 0.
               cpatch%mmean_storage_resp      (ico) = 0.
               cpatch%mmean_vleaf_resp        (ico) = 0.
               cpatch%mmean_light_level       (ico) = 0.
               cpatch%mmean_light_level_beam  (ico) = 0.
               cpatch%mmean_light_level_diff  (ico) = 0.
               cpatch%mmean_beamext_level     (ico) = 0.
               cpatch%mmean_diffext_level     (ico) = 0.
               cpatch%mmean_lambda_light      (ico) = 0.
               cpatch%mmean_mort_rate       (:,ico) = 0.
            end do
         end do
      end do


      !------------------------------------------------------------------------------------!
      !      The mean diurnal cycle is flushed here.                                       !
      !------------------------------------------------------------------------------------!
      if (iqoutput > 0) then
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)
               do ico=1,cpatch%ncohorts
                  !----- Convert GPP and plant respiration to kgC/plant/year. -------------!
                  cpatch%qmean_gpp         (:,ico) = 0.0
                  cpatch%qmean_leaf_resp   (:,ico) = 0.0
                  cpatch%qmean_root_resp   (:,ico) = 0.0
                  cpatch%qmean_par_l       (:,ico) = 0.0
                  cpatch%qmean_par_l_beam  (:,ico) = 0.0
                  cpatch%qmean_par_l_diff  (:,ico) = 0.0
                  cpatch%qmean_fs_open     (:,ico) = 0.0
                  cpatch%qmean_fsw         (:,ico) = 0.0
                  cpatch%qmean_fsn         (:,ico) = 0.0
                  cpatch%qmean_psi_open    (:,ico) = 0.0
                  cpatch%qmean_psi_closed  (:,ico) = 0.0
                  cpatch%qmean_water_supply(:,ico) = 0.0
               end do

               csite%qmean_rh              (:,ipa) = 0.0
               csite%qmean_albedo          (:,ipa) = 0.0
               csite%qmean_albedo_beam     (:,ipa) = 0.0
               csite%qmean_albedo_diffuse  (:,ipa) = 0.0
            end do
         end do

         cgrid%qmean_fs_open       (:,ipy) = 0.0
         cgrid%qmean_fsw           (:,ipy) = 0.0
         cgrid%qmean_fsn           (:,ipy) = 0.0
         cgrid%qmean_leaf_energy   (:,ipy) = 0.0
         cgrid%qmean_leaf_water    (:,ipy) = 0.0
         cgrid%qmean_leaf_hcap     (:,ipy) = 0.0
         cgrid%qmean_wood_energy   (:,ipy) = 0.0
         cgrid%qmean_wood_water    (:,ipy) = 0.0
         cgrid%qmean_wood_hcap     (:,ipy) = 0.0
         cgrid%qmean_can_theta     (:,ipy) = 0.0
         cgrid%qmean_can_theiv     (:,ipy) = 0.0
         cgrid%qmean_can_shv       (:,ipy) = 0.0
         cgrid%qmean_can_co2       (:,ipy) = 0.0
         cgrid%qmean_can_prss      (:,ipy) = 0.0
         cgrid%qmean_gnd_temp      (:,ipy) = 0.0
         cgrid%qmean_gnd_shv       (:,ipy) = 0.0
         cgrid%qmean_atm_temp      (:,ipy) = 0.0
         cgrid%qmean_rshort        (:,ipy) = 0.0
         cgrid%qmean_rshort_diff   (:,ipy) = 0.0
         cgrid%qmean_rlong         (:,ipy) = 0.0
         cgrid%qmean_atm_shv       (:,ipy) = 0.0
         cgrid%qmean_atm_co2       (:,ipy) = 0.0
         cgrid%qmean_atm_prss      (:,ipy) = 0.0
         cgrid%qmean_atm_vels      (:,ipy) = 0.0
         cgrid%qmean_gpp           (:,ipy) = 0.0
         cgrid%qmean_leaf_resp     (:,ipy) = 0.0
         cgrid%qmean_root_resp     (:,ipy) = 0.0
         cgrid%qmean_plresp        (:,ipy) = 0.0
         cgrid%qmean_nep           (:,ipy) = 0.0
         cgrid%qmean_rh            (:,ipy) = 0.0
         cgrid%qmean_sensible_lc   (:,ipy) = 0.0
         cgrid%qmean_sensible_wc   (:,ipy) = 0.0
         cgrid%qmean_sensible_gc   (:,ipy) = 0.0
         cgrid%qmean_sensible_ac   (:,ipy) = 0.0
         cgrid%qmean_rshort_gnd    (:,ipy) = 0.0
         cgrid%qmean_rlong_gnd     (:,ipy) = 0.0
         cgrid%qmean_rlongup       (:,ipy) = 0.0
         cgrid%qmean_rlong_albedo  (:,ipy) = 0.0
         cgrid%qmean_albedo        (:,ipy) = 0.0
         cgrid%qmean_albedo_beam   (:,ipy) = 0.0
         cgrid%qmean_albedo_diffuse(:,ipy) = 0.0
         cgrid%qmean_nee           (:,ipy) = 0.0
         cgrid%qmean_pcpg          (:,ipy) = 0.0
         cgrid%qmean_evap          (:,ipy) = 0.0
         cgrid%qmean_transp        (:,ipy) = 0.0
         cgrid%qmean_runoff        (:,ipy) = 0.0
         cgrid%qmean_drainage      (:,ipy) = 0.0
         cgrid%qmean_vapor_lc      (:,ipy) = 0.0
         cgrid%qmean_vapor_wc      (:,ipy) = 0.0
         cgrid%qmean_vapor_gc      (:,ipy) = 0.0
         cgrid%qmean_vapor_ac      (:,ipy) = 0.0
         cgrid%qmean_soil_temp   (:,:,ipy) = 0.0
         cgrid%qmean_soil_water  (:,:,ipy) = 0.0
         cgrid%qmean_soil_mstpot (:,:,ipy) = 0.0
         cgrid%qmsqu_gpp           (:,ipy) = 0.0
         cgrid%qmsqu_leaf_resp     (:,ipy) = 0.0
         cgrid%qmsqu_root_resp     (:,ipy) = 0.0
         cgrid%qmsqu_plresp        (:,ipy) = 0.0
         cgrid%qmsqu_nee           (:,ipy) = 0.0
         cgrid%qmsqu_nep           (:,ipy) = 0.0
         cgrid%qmsqu_rh            (:,ipy) = 0.0
         cgrid%qmsqu_sensible_ac   (:,ipy) = 0.0
         cgrid%qmsqu_sensible_lc   (:,ipy) = 0.0
         cgrid%qmsqu_sensible_wc   (:,ipy) = 0.0
         cgrid%qmsqu_sensible_gc   (:,ipy) = 0.0
         cgrid%qmsqu_evap          (:,ipy) = 0.0
         cgrid%qmsqu_transp        (:,ipy) = 0.0
         cgrid%qmsqu_vapor_ac      (:,ipy) = 0.0
         cgrid%qmsqu_vapor_lc      (:,ipy) = 0.0
         cgrid%qmsqu_vapor_wc      (:,ipy) = 0.0
         cgrid%qmsqu_vapor_gc      (:,ipy) = 0.0
      end if
      !------------------------------------------------------------------------------------!
      
   end do

   return
end subroutine zero_ed_monthly_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!                             |--------------------------------|                           !
!                             |** YEARLY AVERAGE SUBROUTINES **|                           !
!                             |--------------------------------|                           !
!==========================================================================================!
!==========================================================================================!
subroutine update_ed_yearly_vars(cgrid)

   use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
   use ed_max_dims, only: n_pft, n_dbh
   use consts_coms, only: pi1
   use allometry, only: ed_biomass
  
   implicit none

   type(edtype),target       :: cgrid
   type(polygontype),pointer :: cpoly
   type(sitetype),pointer    :: csite
   type(patchtype),pointer   :: cpatch
   integer :: ipy,isi,ipa,ico

   ! All agb's are in tC/ha/y; all basal areas are in m2/ha/y.
  
   do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)
      
      cgrid%total_basal_area(ipy) = 0.0
      cgrid%total_basal_area_growth(ipy) = 0.0
      cgrid%total_basal_area_mort(ipy) = 0.0
      cgrid%total_basal_area_recruit(ipy) = 0.0
      cgrid%total_agb(ipy) = 0.0
      cgrid%total_agb_growth(ipy) = 0.0
      cgrid%total_agb_mort(ipy) = 0.0
      cgrid%total_agb_recruit(ipy) = 0.0
      
      ! Loop over sites
      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)
         
         ! Do growth, mortality, harvesting.
         cgrid%total_agb(ipy) = cgrid%total_agb(ipy) + sum(cpoly%agb(:,:,isi)) * cpoly%area(isi)
         
         cgrid%total_basal_area(ipy) = cgrid%total_basal_area(ipy) +  &
              sum(cpoly%basal_area(:,:,isi)) * cpoly%area(isi)
         
         cgrid%total_agb_growth(ipy) = cgrid%total_agb_growth(ipy) +  &
              sum(cpoly%agb_growth(:,:,isi)) * cpoly%area(isi)

         cgrid%total_agb_mort(ipy) = cgrid%total_agb_mort(ipy) +  &
              sum(cpoly%agb_mort(1:n_pft, 1:n_dbh,:)) * cpoly%area(isi)

         cgrid%total_basal_area_growth(ipy) = cgrid%total_basal_area_growth(ipy) +  &
              sum(cpoly%basal_area_growth(1:n_pft,2:n_dbh,isi)) * cpoly%area(isi)

         cgrid%total_basal_area_mort(ipy) = cgrid%total_basal_area_mort(ipy) +  &
              sum(cpoly%basal_area_mort(1:n_pft,2:n_dbh,isi)) * cpoly%area(isi)

         !     cgrid%total_agb_cut =   &
         !          sum(cgrid%cs(1)%agb_cut(1:n_pft, 1:n_dbh)) * 10.0

         ! Loop over cohorts to get recruitment. 
         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
        
            ! Loop over cohorts
            do ico = 1,cpatch%ncohorts

               if(cpatch%new_recruit_flag(ico) == 1)then
                  cgrid%total_agb_recruit(ipy) = cgrid%total_agb_recruit(ipy) +   &
                       cpatch%agb(ico) * cpatch%nplant(ico) * csite%area(ipa) * &
                       cpoly%area(isi)
                  cgrid%total_basal_area_recruit(ipy) =   &
                       cgrid%total_basal_area_recruit(ipy) +   &
                       cpatch%basarea(ico) * cpatch%nplant(ico) *&
                       csite%area(ipa) * cpoly%area(isi)
                  cpatch%new_recruit_flag(ico) = 0
               endif
               cpatch%first_census(ico) = 1

            enddo
            
         enddo

      enddo


   enddo

   return
end subroutine update_ed_yearly_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zero_ed_yearly_vars(cgrid)

   use ed_max_dims, only: n_pft, n_dbh
   use ed_state_vars,only:edtype,polygontype

   implicit none
   integer :: ipy
   type(edtype),target       :: cgrid
   type(polygontype),pointer :: cpoly

   do ipy = 1,cgrid%npolygons
      
      cpoly => cgrid%polygon(ipy)
      
      cpoly%agb_growth        = 0.0
      cpoly%agb_mort          = 0.0
      cpoly%agb_cut           = 0.0
!      cpoly%agb_recruit       = 0.0
      cpoly%basal_area_growth = 0.0
      cpoly%basal_area_mort   = 0.0
      cpoly%basal_area_cut    = 0.0
!      cpoly%basal_area_recruit= 0.0
      
   enddo

   return
end subroutine zero_ed_yearly_vars
!==========================================================================================!
!==========================================================================================!
