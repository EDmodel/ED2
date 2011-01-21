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

         cgrid%avg_albedt(ipy)         = cgrid%avg_albedt(ipy)                             &
                                       + 0.5 * ( cpoly%albedo_beam(isi)                    &
                                               + cpoly%albedo_diffuse(isi) )               &
                                       * cpoly%area(isi) * tfact * polygon_area_i

         cgrid%avg_rlongup(ipy)        = cgrid%avg_rlongup(ipy)                            &
                                       + cpoly%rlongup(isi) * cpoly%area(isi)              &
                                       * tfact * polygon_area_i

      end do
   end do
   return
end subroutine int_met_avg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine normalize_averaged_vars(cgrid,frqsum,dtlsm)

   use grid_coms, only: nzg
   use ed_misc_coms, only: radfrq, current_time
   use ed_state_vars,only:edtype,polygontype,sitetype,patchtype

   

   implicit none
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy,isi,ipa,ico
   real                        :: frqsum,dtlsm,tfact,frqsumi
   integer                     :: k

   !-------------------------------------------------------------------------!
   !  The following procedure normalizes these variables to their true rate. !
   !  They have been integrated over the period frqanl, and must be divided  !
   !  by that period to get their value per unit time.                       !
   !-------------------------------------------------------------------------!

   
   frqsumi = 1.0 / frqsum
   tfact = dtlsm * frqsumi
   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            csite%avg_netrad(ipa)       = csite%avg_netrad(ipa)        * frqsumi
            csite%aux(ipa)              = csite%aux(ipa)               * frqsumi
            csite%avg_vapor_vc(ipa)     = csite%avg_vapor_vc(ipa)      * frqsumi
            csite%avg_dew_cg(ipa)       = csite%avg_dew_cg(ipa)        * frqsumi
            csite%avg_vapor_gc(ipa)     = csite%avg_vapor_gc(ipa)      * frqsumi
            csite%avg_wshed_vg(ipa)     = csite%avg_wshed_vg(ipa)      * frqsumi
            csite%avg_intercepted(ipa)  = csite%avg_intercepted(ipa)   * frqsumi
            csite%avg_throughfall(ipa)  = csite%avg_throughfall(ipa)   * frqsumi
            csite%avg_vapor_ac(ipa)     = csite%avg_vapor_ac(ipa)      * frqsumi
            csite%avg_transp(ipa)       = csite%avg_transp(ipa)        * frqsumi
            csite%avg_evap(ipa)         = csite%avg_evap(ipa)          * frqsumi
            csite%avg_runoff(ipa)       = csite%avg_runoff(ipa)        * frqsumi
            csite%avg_drainage(ipa)     = csite%avg_drainage(ipa)      * frqsumi
            csite%avg_sensible_vc(ipa)  = csite%avg_sensible_vc(ipa)   * frqsumi
            csite%avg_qwshed_vg(ipa)    = csite%avg_qwshed_vg(ipa)     * frqsumi
            csite%avg_qintercepted(ipa) = csite%avg_qintercepted(ipa)  * frqsumi
            csite%avg_qthroughfall(ipa) = csite%avg_qthroughfall(ipa)  * frqsumi
            csite%avg_sensible_gc(ipa)  = csite%avg_sensible_gc(ipa)   * frqsumi
            csite%avg_sensible_ac(ipa)  = csite%avg_sensible_ac(ipa)   * frqsumi
            csite%avg_carbon_ac(ipa)    = csite%avg_carbon_ac(ipa)     * frqsumi
            csite%avg_runoff_heat(ipa)  = csite%avg_runoff_heat(ipa)   * frqsumi
            csite%avg_drainage_heat(ipa)= csite%avg_drainage_heat(ipa) * frqsumi
            csite%avg_rk4step(ipa)      = csite%avg_rk4step(ipa)       * frqsumi

         
            do k=cpoly%lsl(isi),nzg
               csite%avg_sensible_gg(k,ipa) = csite%avg_sensible_gg(k,ipa) * frqsumi
               csite%avg_smoist_gg(k,ipa)   = csite%avg_smoist_gg(k,ipa)   * frqsumi
               csite%avg_smoist_gc(k,ipa)   = csite%avg_smoist_gc(k,ipa)   * frqsumi
               csite%aux_s(k,ipa)           = csite%aux_s(k,ipa)           * frqsumi
            end do
            
            !     Available water is added every dtlsm, so we normalise using tfact...
            csite%avg_available_water(ipa) = csite%avg_available_water(ipa)       * tfact
            
            do ico=1,cpatch%ncohorts
               !-- Normalization of cohort level state variables
               !-- leaf and root respiration are state variables, ie
               !-- used to update other quanities in the code. They
               !-- are updated in the photosynthesis driver.  They are
               !-- only normalized briefly for the output driver
               !-- to achieve a rate, then they are zero'd in
               !-- reset_averaged_vars.
               
               ! The following variables were assumed to be integrated
               !cpatch%leaf_respiration = cpatch%leaf_respiration * dtlsm
               !cpatch%root_respiration = cpatch%root_respiration * dtlsm

               ! For IO - they should be replaced by these guys
               ! units: [umol/m2/s * steps]    -> [umol/m2/s]
               cpatch%mean_leaf_resp(ico)    = cpatch%mean_leaf_resp(ico)    * tfact
               cpatch%mean_root_resp(ico)    = cpatch%mean_root_resp(ico)    * tfact
               cpatch%mean_gpp(ico)          = cpatch%mean_gpp(ico)          * tfact
               cpatch%mean_storage_resp(ico) = cpatch%mean_storage_resp(ico) * tfact
               cpatch%mean_growth_resp(ico)  = cpatch%mean_growth_resp(ico)  * tfact
               cpatch%mean_vleaf_resp(ico)   = cpatch%mean_vleaf_resp(ico)   * tfact

            end do
            
            !------ Budget variables. -----------------------------------------------------!
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

   use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
   
   implicit none
   type(edtype),target :: cgrid
   type(polygontype),pointer :: cpoly
   type(sitetype),pointer :: csite
   type(patchtype),pointer :: cpatch
   integer :: ipy,isi,ipa,ico

   do ipy = 1,cgrid%npolygons

      ! Should this be here as well?
      cgrid%cbudget_nep(ipy)       = 0.0

      ! Reset the meteorological diagnostic
   
      cgrid%avg_nir_beam(ipy)       = 0.0
      cgrid%avg_nir_diffuse(ipy)    = 0.0
      cgrid%avg_par_beam(ipy)       = 0.0
      cgrid%avg_par_diffuse(ipy)    = 0.0
      cgrid%avg_atm_tmp(ipy)        = 0.0
      cgrid%avg_atm_shv(ipy)        = 0.0
      cgrid%avg_rshort(ipy)         = 0.0
      cgrid%avg_rshort_diffuse(ipy) = 0.0
      cgrid%avg_rlong(ipy)          = 0.0
      cgrid%avg_pcpg(ipy)           = 0.0
      cgrid%avg_qpcpg(ipy)          = 0.0
      cgrid%avg_dpcpg(ipy)          = 0.0
      cgrid%avg_vels(ipy)           = 0.0
      cgrid%avg_atm_prss(ipy)       = 0.0
      cgrid%avg_exner(ipy)          = 0.0
      cgrid%avg_geoht(ipy)          = 0.0
      cgrid%avg_atm_co2(ipy)        = 0.0
      cgrid%avg_albedt(ipy)         = 0.0
      cgrid%avg_rlongup(ipy)        = 0.0
   
      !
      cgrid%avg_drainage(ipy)       = 0.0
      cgrid%avg_evap(ipy)           = 0.0
      cgrid%avg_transp(ipy)         = 0.0
      cgrid%avg_soil_temp(:,ipy)    = 0.0
      cgrid%avg_soil_water(:,ipy)   = 0.0
      cgrid%avg_soil_energy(:,ipy)  = 0.0
      cgrid%avg_soil_fracliq(:,ipy) = 0.0


      cpoly => cgrid%polygon(ipy)
      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         cpoly%avg_atm_tmp(isi)          = 0.0
         cpoly%avg_atm_shv(isi)          = 0.0
         cpoly%avg_atm_prss(isi)         = 0.0

         cpoly%avg_soil_temp(:,isi)      = 0.0
         cpoly%avg_soil_water(:,isi)     = 0.0
         cpoly%avg_soil_energy(:,isi)    = 0.0

         do ipa = 1,csite%npatches
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
            csite%avg_vapor_vc(ipa)         = 0.0
            csite%avg_dew_cg(ipa)           = 0.0
            csite%avg_vapor_gc(ipa)         = 0.0
            csite%avg_wshed_vg(ipa)         = 0.0
            csite%avg_intercepted(ipa)      = 0.0
            csite%avg_throughfall(ipa)      = 0.0
            csite%avg_vapor_ac(ipa)         = 0.0
            csite%avg_transp(ipa)           = 0.0
            csite%avg_evap(ipa)             = 0.0
            csite%avg_netrad(ipa)           = 0.0
            csite%avg_smoist_gg(:,ipa)      = 0.0
            csite%avg_smoist_gc(:,ipa)      = 0.0
            csite%avg_runoff(ipa)           = 0.0
            csite%avg_runoff_heat(ipa)      = 0.0
            csite%avg_drainage(ipa)         = 0.0
            csite%avg_drainage_heat(ipa)    = 0.0
            csite%avg_sensible_vc(ipa)      = 0.0
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
         
            do ico=1,cpatch%ncohorts
               cpatch%leaf_respiration(ico)      = 0.0
               cpatch%root_respiration(ico)      = 0.0
               cpatch%gpp(ico)                   = 0.0
               cpatch%mean_leaf_resp(ico)        = 0.0
               cpatch%mean_root_resp(ico)        = 0.0
               cpatch%mean_gpp(ico)              = 0.0
               cpatch%mean_storage_resp(ico)     = 0.0 
               cpatch%mean_growth_resp(ico)      = 0.0
               cpatch%mean_vleaf_resp(ico)       = 0.0
            end do
         end do


      end do

      cpoly%avg_soil_temp(:,:)      = 0.0
      cpoly%avg_soil_water(:,:)     = 0.0
      cpoly%avg_soil_energy(:,:)    = 0.0
      cpoly%avg_soil_fracliq(:,:)   = 0.0

   end do


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
!    This subroutine integrates the daily averages for state vars. This is called after    !
! each time step in case daily or monthly averages were requested.                         !
!------------------------------------------------------------------------------------------!
subroutine integrate_ed_daily_output_state(cgrid)
   use ed_state_vars        , only : edtype              & ! structure
                                   , polygontype         & ! structure
                                   , sitetype            & ! structure
                                   , patchtype           ! ! structure
   use grid_coms            , only : nzg                 ! ! intent(in)
   use ed_max_dims          , only : n_dbh               & ! intent(in)
                                   , n_pft               & ! intent(in)
                                   , n_dist_types        & ! intent(in)
                                   , n_mort              ! ! intent(in)
   use ed_misc_coms         , only : dtlsm               ! ! intent(in)
   use pft_coms             , only : sla                 ! ! intent(in)
   use canopy_radiation_coms, only : rshort_twilight_min ! ! intent(in)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy, isi, ipa, ico, imt
   real                        :: poly_area_i,site_area_i
   real                        :: forest_site,forest_site_i, forest_poly
   real                        :: poly_lai,site_lai,patch_lai,patch_lai_i
   real                        :: poly_lma,site_lma,patch_lma
   real                        :: sss_can_theta, sss_can_theiv, sss_can_shv
   real                        :: sss_can_co2, sss_can_prss
   real                        :: pss_veg_water, pss_veg_energy, pss_veg_hcap
   real                        :: sss_veg_water, sss_veg_energy, sss_veg_hcap
   real                        :: rshort_tot
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      !----- Initialize auxiliary variables to add sitetype variables. --------------------!
      site_lma         = 0.
      site_lai         = 0.
      sss_veg_energy   = 0.
      sss_veg_water    = 0.
      sss_veg_hcap     = 0.
      sss_can_theta    = 0.
      sss_can_theiv    = 0.
      sss_can_shv      = 0.
      sss_can_co2      = 0.
      sss_can_prss     = 0.
      forest_poly      = 0.

      
      siteloop: do isi=1, cpoly%nsites


         !---------------------------------------------------------------------------------!
         !     Determine if there was any beam radiation, and compute the total (I think   !
         ! rshort is always the total, but not sure.  This is the same test done in the    !
         ! radiation driver.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpoly%cosaoi(isi) <= 0.0) then
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
         patch_lai        = 0.
         patch_lma        = 0.
         pss_veg_energy   = 0.
         pss_veg_water    = 0.
         pss_veg_hcap     = 0.

         !----- Looping through the patches to normalize the sum of all cohorts. ----------!
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)

            if (cpatch%ncohorts > 0) then
               pss_veg_energy = pss_veg_energy + sum(cpatch%veg_energy) * csite%area(ipa)
               pss_veg_water  = pss_veg_water  + sum(cpatch%veg_water ) * csite%area(ipa)
               pss_veg_hcap   = pss_veg_hcap   + sum(cpatch%hcapveg   ) * csite%area(ipa)
               patch_lai_i = 1./max(tiny(1.),sum(cpatch%lai,cpatch%solvable))
               patch_lai = patch_lai + csite%area(ipa)*sum(cpatch%lai,cpatch%solvable)
               patch_lma = patch_lma + csite%area(ipa)*sum(cpatch%lai/sla(cpatch%pft)      &
                                                          ,cpatch%solvable)
            end if

            do ico=1,cpatch%ncohorts
               if (cpatch%solvable(ico)) then
                  cpatch%dmean_par_v       (ico) = cpatch%dmean_par_v       (ico)          &
                                                 + cpatch%par_v             (ico)
                  cpatch%dmean_par_v_beam  (ico) = cpatch%dmean_par_v_beam  (ico)          &
                                                 + cpatch%par_v_beam        (ico)
                  cpatch%dmean_par_v_diff  (ico) = cpatch%dmean_par_v_diff  (ico)          &
                                                 + cpatch%par_v_diffuse     (ico)

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
                     cpatch%dmean_norm_par_beam(ico)    = cpatch%dmean_norm_par_beam(ico)  &
                                                        + cpatch%norm_par_beam(ico)
                     cpatch%dmean_norm_par_diff(ico)    = cpatch%dmean_norm_par_diff(ico)  &
                                                        + cpatch%norm_par_diff(ico)
                     cpatch%dmean_lambda_light(ico)     = cpatch%dmean_lambda_light(ico)   &
                                                        + cpatch%lambda_light(ico)
                  end if

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

         site_lai         = site_lai       + patch_lai * cpoly%area(isi)
         site_lma         = site_lma       + patch_lma * cpoly%area(isi)
         sss_veg_energy   = sss_veg_energy + (pss_veg_energy*site_area_i) * cpoly%area(isi)
         sss_veg_water    = sss_veg_water  + (pss_veg_water *site_area_i) * cpoly%area(isi)
         sss_veg_hcap     = sss_veg_hcap   + (pss_veg_hcap  *site_area_i) * cpoly%area(isi)


         sss_can_theta  = sss_can_theta                                                    &
                        + cpoly%area(isi) * (sum(csite%can_theta*csite%area) * site_area_i)
         sss_can_theiv  = sss_can_theiv                                                    &
                        + cpoly%area(isi) * (sum(csite%can_theiv*csite%area) * site_area_i)
         sss_can_shv    = sss_can_shv                                                      &
                        + cpoly%area(isi) * (sum(csite%can_shv  *csite%area) * site_area_i)
         sss_can_co2    = sss_can_co2                                                      &
                        + cpoly%area(isi) * (sum(csite%can_co2  *csite%area) * site_area_i)
         sss_can_prss   = sss_can_prss                                                     &
                        + cpoly%area(isi) * (sum(csite%can_prss *csite%area) * site_area_i)
      end do siteloop
      
      !------------------------------------------------------------------------------------!
      !    Variables already averaged at the polygontype level, just add them to edtype    !
      ! level.                                                                             !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_veg_energy(ipy)   = cgrid%dmean_veg_energy(ipy)                          &
                                    + sss_veg_energy * poly_area_i
      cgrid%dmean_veg_water(ipy)    = cgrid%dmean_veg_water(ipy)                           &
                                    + sss_veg_water * poly_area_i
      cgrid%dmean_veg_hcap(ipy)     = cgrid%dmean_veg_hcap(ipy)                            &
                                    + sss_veg_hcap  * poly_area_i
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

      !------------------------------------------------------------------------------------!
      !    Variables already at edtype level, simple integration only.                     !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_atm_temp(ipy) = cgrid%dmean_atm_temp(ipy) + cgrid%met(ipy)%atm_tmp
      cgrid%dmean_rshort(ipy)   = cgrid%dmean_rshort(ipy)   &
           + cgrid%met(ipy)%nir_beam    &
           + cgrid%met(ipy)%nir_diffuse &
           + cgrid%met(ipy)%par_beam    &
           + cgrid%met(ipy)%par_diffuse

      cgrid%dmean_rlong(ipy)    = cgrid%dmean_rlong(ipy)    + cgrid%met(ipy)%rlong
      cgrid%dmean_atm_shv (ipy) = cgrid%dmean_atm_shv (ipy) + cgrid%met(ipy)%atm_shv
      cgrid%dmean_atm_prss(ipy) = cgrid%dmean_atm_prss(ipy) + cgrid%met(ipy)%prss
      cgrid%dmean_atm_vels(ipy) = cgrid%dmean_atm_vels(ipy) + cgrid%met(ipy)%vels
    
      if(site_lai > tiny(1.))then
         cgrid%avg_lma(ipy) = site_lma/site_lai
      end if

   end do polyloop
      
   return
end subroutine integrate_ed_daily_output_state
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine integrate_ed_daily_output_flux(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine integrates the daily average. This is called at the analysis time in  !
! case daily or monthly averages were requested. We take advantage from the previously     !
! averaged variables.                                                                      !
!------------------------------------------------------------------------------------------!
   use ed_state_vars        , only : edtype       & ! structure
                                   , polygontype  & ! structure
                                   , sitetype     & ! structure
                                   , patchtype    ! ! structure
   use grid_coms            , only : nzg          ! ! intent(in)
   use ed_max_dims          , only : n_dbh        & ! intent(in)
                                   , n_pft        & ! intent(in)
                                   , n_dist_types ! ! intent(in)
   use ed_misc_coms         , only : dtlsm        & ! intent(in)
                                   , ddbhi        ! ! intent(in)
   use pft_coms             , only : c2n_leaf
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
   real                          :: poly_area_i
   real                          :: site_area_i
   real                          :: forest_site
   real                          :: forest_site_i
   real                          :: forest_poly
   real                          :: sitesum_leaf_resp
   real                          :: sitesum_root_resp 
   real                          :: sitesum_plresp
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
   real                          :: patchsum_leaf_resp
   real                          :: patchsum_root_resp
   real                          :: patchsum_leaf_litter
   real                          :: patchsum_root_litter
   real                          :: patchsum_leaf_litterN
   real                          :: patchsum_root_litterN

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
      sitesum_rh              = 0.
      sitesum_leaf_resp       = 0.
      sitesum_root_resp       = 0.
      sitesum_plresp          = 0.
      sitesum_gpp             = 0.
      sitesum_gpp_dbh         = 0.
      sitesum_evap            = 0.
      sitesum_transp          = 0.
      sitesum_co2_residual    = 0.
      sitesum_water_residual  = 0.
      sitesum_energy_residual = 0.
      forest_poly             = 0.

      siteloop: do isi=1, cpoly%nsites
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
         patchsum_leaf_resp    = 0.
         patchsum_root_resp    = 0.
         patchsum_leaf_litter  = 0.
         patchsum_root_litter  = 0.
         patchsum_leaf_litterN = 0.
         patchsum_root_litterN = 0.

         !----- Looping through the patches to normalize the sum of all cohorts. ----------!
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)
            if (cpatch%ncohorts > 0) then
               patchsum_leaf_resp    = patchsum_leaf_resp                                  &
                                     + sum(cpatch%mean_leaf_resp, cpatch%solvable)         &
                                     * csite%area(ipa) 
               patchsum_root_resp    = patchsum_root_resp                                  &
                                     + sum(cpatch%mean_root_resp)   * csite%area(ipa)  
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
            csite%dmean_co2_residual(ipa)    = csite%dmean_co2_residual(ipa)               &
                                             + csite%co2budget_residual(ipa)
            csite%dmean_energy_residual(ipa) = csite%dmean_energy_residual(ipa)            &
                                             + csite%ebudget_residual(ipa)
            csite%dmean_water_residual(ipa)  = csite%dmean_water_residual(ipa)             &
                                             + csite%wbudget_residual(ipa)
            csite%dmean_rh(ipa)              = csite%dmean_rh(ipa)                         &
                                             + csite%co2budget_rh(ipa)
            csite%dmean_rk4step(ipa)         = csite%dmean_rk4step(ipa)                    &
                                             + csite%avg_rk4step(ipa)
         end do patchloop
         
         !---------------------------------------------------------------------------------!
         !      Variables already averaged at the sitetype level, add them to polygontype  !
         ! level.                                                                          !
         !---------------------------------------------------------------------------------!
         sitesum_leaf_resp       = sitesum_leaf_resp                                       &
                                 + (patchsum_leaf_resp * site_area_i) * cpoly%area(isi)
         sitesum_root_resp       = sitesum_root_resp                                       &
                                 + (patchsum_root_resp * site_area_i) * cpoly%area(isi)

         sitesum_gpp             = sitesum_gpp                                             &
                                 + (sum(csite%co2budget_gpp * csite%area) * site_area_i)   &
                                 * cpoly%area(isi)
         sitesum_rh              = sitesum_rh                                              &
                                 + (sum(csite%co2budget_rh    *csite%area) * site_area_i)  &
                                 * cpoly%area(isi)
         sitesum_plresp          = sitesum_plresp                                          &
                                 + (sum(csite%co2budget_plresp * csite%area)* site_area_i) &
                                 * cpoly%area(isi)
         
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

      cgrid%dmean_leaf_resp(ipy)       = cgrid%dmean_leaf_resp(ipy)                        &
                                       + sitesum_leaf_resp * poly_area_i
      cgrid%dmean_root_resp(ipy)       = cgrid%dmean_root_resp(ipy)                        &
                                       + sitesum_root_resp * poly_area_i

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

      cgrid%dmean_plresp(ipy)          = cgrid%dmean_plresp(ipy)                           &
                                       + sitesum_plresp * poly_area_i
      cgrid%dmean_nep(ipy)             = cgrid%dmean_nep(ipy)                              &
                                       + (sitesum_gpp-sitesum_rh-sitesum_plresp)           &
                                       * poly_area_i

      cgrid%dmean_co2_residual(ipy)    = cgrid%dmean_co2_residual(ipy)                     &
                                       + sitesum_co2_residual * poly_area_i
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
      end do
      cgrid%dmean_sensible_vc(ipy) = cgrid%dmean_sensible_vc(ipy)                          &
                                   + cgrid%avg_sensible_vc(ipy) 
      cgrid%dmean_sensible_gc(ipy) = cgrid%dmean_sensible_gc(ipy)                          &
                                   + cgrid%avg_sensible_gc(ipy)
      cgrid%dmean_sensible_ac(ipy) = cgrid%dmean_sensible_ac(ipy)                          &
                                   + cgrid%avg_sensible_ac(ipy)

      !------ Integrate the NEE with the conventional NEE sign (< 0 = uptake). ------------!
      cgrid%dmean_nee     (ipy) = cgrid%dmean_nee(ipy)      - cgrid%avg_carbon_ac(ipy)

      cgrid%dmean_pcpg    (ipy) = cgrid%dmean_pcpg    (ipy) + cgrid%avg_pcpg    (ipy)
      cgrid%dmean_evap    (ipy) = cgrid%dmean_evap    (ipy) + cgrid%avg_evap    (ipy)
      cgrid%dmean_transp  (ipy) = cgrid%dmean_transp  (ipy) + cgrid%avg_transp  (ipy)
      cgrid%dmean_runoff  (ipy) = cgrid%dmean_runoff  (ipy) + cgrid%avg_runoff  (ipy)
      cgrid%dmean_drainage(ipy) = cgrid%dmean_drainage(ipy) + cgrid%avg_drainage(ipy)
      cgrid%dmean_vapor_vc(ipy) = cgrid%dmean_vapor_vc(ipy) + cgrid%avg_vapor_vc(ipy)
      cgrid%dmean_vapor_gc(ipy) = cgrid%dmean_vapor_gc(ipy) + cgrid%avg_vapor_gc(ipy)
      cgrid%dmean_vapor_ac(ipy) = cgrid%dmean_vapor_ac(ipy) + cgrid%avg_vapor_ac(ipy)
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
      save_daily   = imoutput > 0 .or. idoutput > 0
      save_monthly = imoutput > 0
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
         cgrid%dmean_gpp(ipy)       = sss_gpp       * poly_area_i
         cgrid%dmean_leaf_resp(ipy) = sss_leaf_resp * poly_area_i
         cgrid%dmean_root_resp(ipy) = sss_root_resp * poly_area_i
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
   logical                                         :: any_solvable
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
   real                                            :: veg_fliq
   real                                            :: dtlsm_o_daylight
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
      cgrid%dmean_veg_energy(ipy)   = cgrid%dmean_veg_energy(ipy)   * dtlsm_o_daysec
      cgrid%dmean_veg_hcap(ipy)     = cgrid%dmean_veg_hcap(ipy)     * dtlsm_o_daysec
      cgrid%dmean_veg_water(ipy)    = cgrid%dmean_veg_water(ipy)    * dtlsm_o_daysec
      cgrid%dmean_can_theta(ipy)    = cgrid%dmean_can_theta(ipy)    * dtlsm_o_daysec
      cgrid%dmean_can_theiv(ipy)    = cgrid%dmean_can_theiv(ipy)    * dtlsm_o_daysec
      cgrid%dmean_can_shv(ipy)      = cgrid%dmean_can_shv(ipy)      * dtlsm_o_daysec
      cgrid%dmean_can_co2(ipy)      = cgrid%dmean_can_co2(ipy)      * dtlsm_o_daysec
      cgrid%dmean_can_prss(ipy)     = cgrid%dmean_can_prss(ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_temp(ipy)     = cgrid%dmean_atm_temp(ipy)     * dtlsm_o_daysec
      cgrid%dmean_rshort(ipy)       = cgrid%dmean_rshort(ipy)       * dtlsm_o_daysec
      cgrid%dmean_rlong(ipy)        = cgrid%dmean_rlong(ipy)        * dtlsm_o_daysec
      cgrid%dmean_atm_shv(ipy)      = cgrid%dmean_atm_shv(ipy)      * dtlsm_o_daysec
      cgrid%dmean_atm_prss(ipy)     = cgrid%dmean_atm_prss(ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_vels(ipy)     = cgrid%dmean_atm_vels(ipy)     * dtlsm_o_daysec

      !------------------------------------------------------------------------------------!
      !     Finding the canopy variables that are not conserved when pressure changes.     !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_can_temp(ipy)     = cgrid%dmean_can_theta(ipy)                           &
                                    * (p00i * cgrid%dmean_can_prss(ipy)) ** rocp
      cgrid%dmean_can_rhos(ipy)     = idealdenssh (cgrid%dmean_can_prss(ipy)               &
                                                  ,cgrid%dmean_can_temp(ipy)               &
                                                  ,cgrid%dmean_can_shv (ipy) )

      !----- Finding vegetation temperature -----------------------------------------------!
      call qwtk(cgrid%dmean_veg_energy(ipy),cgrid%dmean_veg_water(ipy)                     &
               ,cgrid%dmean_veg_hcap(ipy),cgrid%dmean_veg_temp(ipy),veg_fliq)
      
      !------------------------------------------------------------------------------------!
      !     State variables, updated every frqsum, so these are normalized by              !
      ! frqsum/day_sec.                                                                    !
      !------------------------------------------------------------------------------------!
      do k=1,nzg
         cgrid%dmean_soil_temp (k,ipy) = cgrid%dmean_soil_temp (k,ipy) * frqsum_o_daysec
         cgrid%dmean_soil_water(k,ipy) = cgrid%dmean_soil_water(k,ipy) * frqsum_o_daysec
      end do
      !----- Precipitation and runoff. ----------------------------------------------------!
      cgrid%dmean_pcpg     (ipy)  = cgrid%dmean_pcpg     (ipy) * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_runoff   (ipy)  = cgrid%dmean_runoff   (ipy) * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_drainage (ipy)  = cgrid%dmean_drainage (ipy) * frqsum_o_daysec ! kg/m2/s

      !----- Vapor flux. ------------------------------------------------------------------!
      cgrid%dmean_vapor_vc(ipy)   = cgrid%dmean_vapor_vc(ipy)  * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_vapor_gc(ipy)   = cgrid%dmean_vapor_gc(ipy)  * frqsum_o_daysec ! kg/m2/s
      cgrid%dmean_vapor_ac(ipy)   = cgrid%dmean_vapor_ac(ipy)  * frqsum_o_daysec ! kg/m2/s


      !------------------------------------------------------------------------------------!
      !     Flux variables, updated every frqsum, so these are normalized by               !
      ! frqsum/day_sec.                                                                    !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_evap       (ipy)  = cgrid%dmean_evap       (ipy)  * frqsum_o_daysec
      cgrid%dmean_transp     (ipy)  = cgrid%dmean_transp     (ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_vc(ipy)  = cgrid%dmean_sensible_vc(ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_gc(ipy)  = cgrid%dmean_sensible_gc(ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_ac(ipy)  = cgrid%dmean_sensible_ac(ipy)  * frqsum_o_daysec

      !------------------------------------------------------------------------------------!
      !      Carbon flux variables should be total flux integrated over the day at this    !
      ! point in umol/m2/s.  We just multiply by one year in seconds and convert to kgC,   !
      ! so the units will be kgC/m2/yr.                                                    !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_nee        (ipy)  = cgrid%dmean_nee(ipy)                                 &
                                    * umols_2_kgCyr * frqsum_o_daysec
      cgrid%dmean_plresp     (ipy)  = cgrid%dmean_plresp    (ipy)                          &
                                    * umols_2_kgCyr * frqsum_o_daysec
      cgrid%dmean_nep        (ipy)  = cgrid%dmean_nep       (ipy)                          &
                                    * umols_2_kgCyr * frqsum_o_daysec
      cgrid%dmean_gpp_dbh  (:,ipy)  = cgrid%dmean_gpp_dbh (:,ipy)                          &
                                    * umols_2_kgCyr * frqsum_o_daysec

      cgrid%dmean_co2_residual   (ipy) = cgrid%dmean_co2_residual   (ipy)                  &
                                       * umols_2_kgCyr * frqsum_o_daysec
      cgrid%dmean_energy_residual(ipy) = cgrid%dmean_energy_residual(ipy) * frqsum_o_daysec
      cgrid%dmean_water_residual (ipy) = cgrid%dmean_water_residual (ipy) * frqsum_o_daysec


      poly_area_i = 1./sum(cpoly%area)
      sss_growth_resp  = 0.
      sss_storage_resp = 0.
      sss_vleaf_resp   = 0.
      sss_rh           = 0.
      sss_fsn          = 0.
      sss_fsw          = 0.
      sss_fs_open      = 0.

      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !     Find the average day length.                                                !
         !---------------------------------------------------------------------------------!
         dtlsm_o_daylight = dtlsm / cpoly%daylight(isi)
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
            
            
            any_solvable = .false.
            if (cpatch%ncohorts > 0) then
               any_solvable = any(cpatch%solvable(1:cpatch%ncohorts))
            end if


            cohortloop: do ico=1,cpatch%ncohorts

               !---------------------------------------------------------------------------!
               !     These variables must be scaled.  They are updated every time step.    !
               !---------------------------------------------------------------------------!
               cpatch%dmean_par_v       (ico) = cpatch%dmean_par_v       (ico)             &
                                              * dtlsm_o_daysec
               cpatch%dmean_par_v_beam  (ico) = cpatch%dmean_par_v_beam  (ico)             &
                                              * dtlsm_o_daysec
               cpatch%dmean_par_v_diff  (ico) = cpatch%dmean_par_v_diff  (ico)             &
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
                  cpatch%dmean_norm_par_beam   (ico) = cpatch%dmean_norm_par_beam   (ico)  &
                                                     * dtlsm_o_daylight
                  cpatch%dmean_norm_par_diff   (ico) = cpatch%dmean_norm_par_diff   (ico)  &
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
                  cpatch%dmean_norm_par_beam   (ico) = 0.
                  cpatch%dmean_norm_par_diff   (ico) = 0.
                  cpatch%dmean_lambda_light    (ico) = 0.
               end if
            end do cohortloop

            
            
            !------------------------------------------------------------------------------!
            !     Integrate the fraction of open stomata 
            !------------------------------------------------------------------------------!
            if (any_solvable) then
               patch_laiall_i = 1./max(tiny(1.),sum(cpatch%lai,cpatch%solvable))
               pss_fsn     = pss_fsn + csite%area(ipa)                                     &
                           * ( sum(cpatch%dmean_fsn * cpatch%lai,cpatch%solvable)          &
                             * patch_laiall_i)
               pss_fsw     = pss_fsw + csite%area(ipa)                                     &
                           * ( sum(cpatch%dmean_fsw * cpatch%lai,cpatch%solvable)          &
                             * patch_laiall_i)
               pss_fs_open = pss_fs_open + csite%area(ipa)                                 &
                           * ( sum(cpatch%dmean_fs_open * cpatch%lai,cpatch%solvable)      &
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
            csite%dmean_rh(ipa)              = csite%dmean_rh(ipa) * umols_2_kgCyr         &
                                             * frqsum_o_daysec

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

         !----- Adding this patch to the site sum. ----------------------------------------!

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
         sss_rh           = sss_rh  + (sum(csite%dmean_rh * csite%area) * site_area_i)     &
                                    * cpoly%area(isi)
      end do siteloop
      
      !------------------------------------------------------------------------------------!
      !     Finding the area indices per PFT class.                                        !
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
      
      cgrid%dmean_rh(ipy)      = cgrid%dmean_rh(ipy)      + sss_rh      * poly_area_i      
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
      cgrid%dmean_vapor_vc       (ipy) = 0.
      cgrid%dmean_vapor_gc       (ipy) = 0.
      cgrid%dmean_vapor_ac       (ipy) = 0.
      cgrid%dmean_gpp            (ipy) = 0.
      cgrid%dmean_evap           (ipy) = 0.
      cgrid%dmean_transp         (ipy) = 0.
      cgrid%dmean_sensible_vc    (ipy) = 0.
      cgrid%dmean_sensible_gc    (ipy) = 0.
      cgrid%dmean_sensible_ac    (ipy) = 0.
 
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
      cgrid%dmean_gpp_dbh      (:,ipy) = 0.
      cgrid%dmean_fs_open        (ipy) = 0.
      cgrid%dmean_fsw            (ipy) = 0.
      cgrid%dmean_fsn            (ipy) = 0.
      cgrid%dmean_veg_energy     (ipy) = 0.
      cgrid%dmean_veg_hcap       (ipy) = 0.
      cgrid%dmean_veg_water      (ipy) = 0.
      cgrid%dmean_veg_temp       (ipy) = 0.
      cgrid%dmean_can_temp       (ipy) = 0.
      cgrid%dmean_can_shv        (ipy) = 0.
      cgrid%dmean_can_co2        (ipy) = 0.
      cgrid%dmean_can_rhos       (ipy) = 0.
      cgrid%dmean_can_prss       (ipy) = 0.
      cgrid%dmean_can_theta      (ipy) = 0.
      cgrid%dmean_can_theiv      (ipy) = 0.
      cgrid%dmean_atm_temp       (ipy) = 0.
      cgrid%dmean_rshort         (ipy) = 0.
      cgrid%dmean_rlong          (ipy) = 0.
      cgrid%dmean_atm_shv        (ipy) = 0.
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

            cpatch => csite%patch(ipa)
            do ico=1, cpatch%ncohorts
               cpatch%dmean_gpp(ico)              = 0.
               cpatch%dmean_leaf_resp(ico)        = 0.
               cpatch%dmean_root_resp(ico)        = 0.
               cpatch%dmean_par_v(ico)            = 0.
               cpatch%dmean_par_v_beam(ico)       = 0.
               cpatch%dmean_par_v_diff(ico)       = 0.
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
               cpatch%dmean_norm_par_beam   (ico) = 0.
               cpatch%dmean_norm_par_diff   (ico) = 0.
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
      cgrid%mmean_vapor_vc      (ipy) = cgrid%mmean_vapor_vc      (ipy)                    &
                                      + cgrid%dmean_vapor_vc      (ipy)
      cgrid%mmean_sensible_ac   (ipy) = cgrid%mmean_sensible_ac   (ipy)                    &
                                      + cgrid%dmean_sensible_ac   (ipy)
      cgrid%mmean_sensible_gc   (ipy) = cgrid%mmean_sensible_gc   (ipy)                    &
                                      + cgrid%dmean_sensible_gc   (ipy)
      cgrid%mmean_sensible_vc   (ipy) = cgrid%mmean_sensible_vc   (ipy)                    &
                                      + cgrid%dmean_sensible_vc   (ipy)
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

      cgrid%mmean_gpp_dbh     (:,ipy) = cgrid%mmean_gpp_dbh     (:,ipy)                    &
                                      + cgrid%dmean_gpp_dbh     (:,ipy)

      cgrid%mmean_lai_pft     (:,ipy) = cgrid%mmean_lai_pft     (:,ipy)                    &
                                      + cgrid%lai_pft           (:,ipy)
      cgrid%mmean_wpa_pft     (:,ipy) = cgrid%mmean_wpa_pft     (:,ipy)                    &
                                      + cgrid%wpa_pft           (:,ipy)
      cgrid%mmean_wai_pft     (:,ipy) = cgrid%mmean_wai_pft     (:,ipy)                    &
                                      + cgrid%wai_pft           (:,ipy)

      cgrid%mmean_veg_energy    (ipy) = cgrid%mmean_veg_energy    (ipy)                    &
                                      + cgrid%dmean_veg_energy    (ipy)
      cgrid%mmean_veg_hcap      (ipy) = cgrid%mmean_veg_hcap      (ipy)                    &
                                      + cgrid%dmean_veg_hcap      (ipy)
      cgrid%mmean_veg_water     (ipy) = cgrid%mmean_veg_water     (ipy)                    & 
                                      + cgrid%dmean_veg_water     (ipy)
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
      cgrid%mmean_atm_temp      (ipy) = cgrid%mmean_atm_temp      (ipy)                    &
                                      + cgrid%dmean_atm_temp      (ipy)
      cgrid%mmean_rshort        (ipy) = cgrid%mmean_rshort        (ipy)                    &
                                      + cgrid%dmean_rshort        (ipy)
      cgrid%mmean_rlong         (ipy) = cgrid%mmean_rlong         (ipy)                    &
                                      + cgrid%dmean_rlong        (ipy)

      cgrid%mmean_atm_shv       (ipy) = cgrid%mmean_atm_shv       (ipy)                    &
                                      + cgrid%dmean_atm_shv       (ipy)
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
      !    During the integration stage we keep the sum of squares, it will be converted   !
      ! to standard deviation right before the monthly output.                             !
      !------------------------------------------------------------------------------------!
      cgrid%stdev_gpp     (ipy) = cgrid%stdev_gpp     (ipy)                                &
                                + cgrid%dmean_gpp     (ipy)    ** 2
      cgrid%stdev_evap    (ipy) = cgrid%stdev_evap    (ipy)                                &
                                + cgrid%dmean_evap    (ipy)    ** 2
      cgrid%stdev_transp  (ipy) = cgrid%stdev_transp  (ipy)                                &
                                + cgrid%dmean_transp  (ipy)    ** 2
      cgrid%stdev_nep     (ipy) = cgrid%stdev_nep     (ipy)                                &
                                + cgrid%dmean_nep     (ipy)    ** 2
      cgrid%stdev_rh      (ipy) = cgrid%stdev_rh      (ipy)                                &
                                + cgrid%dmean_rh      (ipy)    ** 2
      cgrid%stdev_sensible(ipy) = cgrid%stdev_sensible(ipy)                                &
                                + ( cgrid%dmean_sensible_gc(ipy)                           &
                                  + cgrid%dmean_sensible_vc(ipy)) ** 2

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

            csite%mmean_rh(ipa)              = csite%mmean_rh(ipa) + csite%dmean_rh(ipa)

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
               cpatch%mmean_par_v       (ico) = cpatch%mmean_par_v       (ico)             &
                                              + cpatch%dmean_par_v       (ico)
               cpatch%mmean_par_v_beam  (ico) = cpatch%mmean_par_v_beam  (ico)             &
                                              + cpatch%dmean_par_v_beam  (ico)
               cpatch%mmean_par_v_diff  (ico) = cpatch%mmean_par_v_diff  (ico)             &
                                              + cpatch%dmean_par_v_diff  (ico)

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
               cpatch%mmean_norm_par_beam   (ico) = cpatch%mmean_norm_par_beam(ico)        &
                                                  + cpatch%dmean_norm_par_beam(ico)
               cpatch%mmean_norm_par_diff   (ico) = cpatch%mmean_norm_par_diff(ico)        &
                                                  + cpatch%dmean_norm_par_diff(ico)
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
                                   , dagei         ! ! intent(in)
   use ed_max_dims          , only : n_pft         & ! intent(in)
                                   , n_dbh         & ! intent(in)
                                   , n_age         & ! intent(in)
                                   , n_dist_types  & ! intent(in)
                                   , n_mort        ! ! intent(in)
   use consts_coms          , only : p00i          & ! intent(in)
                                   , rocp          & ! intent(in)
                                   , pio4          & ! intent(in)
                                   , umol_2_kgC    & ! intent(in)
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
   real, dimension(n_pft)                          :: pss_bseeds_pft, sss_bseeds_pft
   real, dimension(n_pft,n_dbh)                    :: pss_bseeds    , pss_pldens
   integer                                         :: lmon
   integer                                         :: ipy
   integer                                         :: isi
   integer                                         :: ipa
   integer                                         :: ico
   integer                                         :: imt
   integer                                         :: ipft
   integer                                         :: idbh
   integer                                         :: iage
   integer                                         :: ilu
   integer                                         :: jlu
   logical                                         :: forest
   real                                            :: srnonm1
   real                                            :: veg_fliq
   real                                            :: cohort_seeds
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
      cgrid%mmean_evap           (ipy) = cgrid%mmean_evap           (ipy) * ndaysi
      cgrid%mmean_transp         (ipy) = cgrid%mmean_transp         (ipy) * ndaysi
      cgrid%mmean_vapor_ac       (ipy) = cgrid%mmean_vapor_ac       (ipy) * ndaysi
      cgrid%mmean_vapor_gc       (ipy) = cgrid%mmean_vapor_gc       (ipy) * ndaysi
      cgrid%mmean_vapor_vc       (ipy) = cgrid%mmean_vapor_vc       (ipy) * ndaysi
      cgrid%mmean_sensible_ac    (ipy) = cgrid%mmean_sensible_ac    (ipy) * ndaysi
      cgrid%mmean_sensible_gc    (ipy) = cgrid%mmean_sensible_gc    (ipy) * ndaysi
      cgrid%mmean_sensible_vc    (ipy) = cgrid%mmean_sensible_vc    (ipy) * ndaysi
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
      cgrid%mmean_gpp_dbh      (:,ipy) = cgrid%mmean_gpp_dbh      (:,ipy) * ndaysi
      cgrid%mmean_veg_energy     (ipy) = cgrid%mmean_veg_energy     (ipy) * ndaysi
      cgrid%mmean_veg_hcap       (ipy) = cgrid%mmean_veg_hcap       (ipy) * ndaysi
      cgrid%mmean_veg_water      (ipy) = cgrid%mmean_veg_water      (ipy) * ndaysi
      cgrid%mmean_can_theta      (ipy) = cgrid%mmean_can_theta      (ipy) * ndaysi
      cgrid%mmean_can_theiv      (ipy) = cgrid%mmean_can_theiv      (ipy) * ndaysi
      cgrid%mmean_can_shv        (ipy) = cgrid%mmean_can_shv        (ipy) * ndaysi
      cgrid%mmean_can_co2        (ipy) = cgrid%mmean_can_co2        (ipy) * ndaysi
      cgrid%mmean_can_prss       (ipy) = cgrid%mmean_can_prss       (ipy) * ndaysi
      cgrid%mmean_atm_temp       (ipy) = cgrid%mmean_atm_temp       (ipy) * ndaysi
      cgrid%mmean_rshort         (ipy) = cgrid%mmean_rshort         (ipy) * ndaysi
      cgrid%mmean_rlong          (ipy) = cgrid%mmean_rlong          (ipy) * ndaysi
      cgrid%mmean_atm_shv        (ipy) = cgrid%mmean_atm_shv        (ipy) * ndaysi
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
      !   Here we convert the sum of squares into standard deviation. The standard devi-   !
      ! ation can be written in two different ways, and we will use the latter because it  !
      ! doesn't require previous knowledge of the mean.                                    !
      !              __________________          ____________________________________      !
      !             / SUM_i[X_i - Xm]²          /  / SUM_i[X_i²]        \      1           !
      ! sigma = \  /  ----------------   =  \  /  |  -----------  - Xm²  | ---------       !
      !          \/       N - 1              \/    \      N             /   1 - 1/N        !
      !                                                                                    !
      ! srnonm1 is the square root of 1 / (1 - 1/N)                                        !
      !------------------------------------------------------------------------------------!
      srnonm1 = sqrt(1./(1.0-ndaysi))
      !------------------------------------------------------------------------------------!
      cgrid%stdev_gpp     (ipy) = srnonm1 * sqrt( cgrid%stdev_gpp     (ipy) * ndaysi       &
                                                - cgrid%mmean_gpp     (ipy) ** 2)
      cgrid%stdev_evap    (ipy) = srnonm1 * sqrt( cgrid%stdev_evap    (ipy) * ndaysi       &
                                                - cgrid%mmean_evap    (ipy) ** 2)
      cgrid%stdev_transp  (ipy) = srnonm1 * sqrt( cgrid%stdev_transp  (ipy) * ndaysi       &
                                                - cgrid%mmean_transp  (ipy) ** 2)
      cgrid%stdev_nep     (ipy) = srnonm1 * sqrt( cgrid%stdev_nep     (ipy) * ndaysi       &
                                                - cgrid%mmean_nep     (ipy) ** 2)
      cgrid%stdev_rh      (ipy) = srnonm1 * sqrt( cgrid%stdev_rh      (ipy) * ndaysi       &
                                                - cgrid%mmean_rh      (ipy) ** 2)
      cgrid%stdev_sensible(ipy) = srnonm1 * sqrt( cgrid%stdev_sensible(ipy) * ndaysi       &
                                                - ( cgrid%mmean_sensible_vc(ipy)           &
                                                  + cgrid%mmean_sensible_gc(ipy))** 2)
  
      !---- Finding the derived average properties from vegetation and canopy air space. --!
      call qwtk(cgrid%mmean_veg_energy(ipy),cgrid%mmean_veg_water(ipy)                     &
               ,cgrid%mmean_veg_hcap(ipy),cgrid%mmean_veg_temp(ipy),veg_fliq)

      cgrid%mmean_can_temp    (ipy) = cgrid%mmean_can_theta(ipy)                           &
                                    * (p00i * cgrid%mmean_can_prss(ipy)) ** rocp
      cgrid%mmean_can_rhos    (ipy) = idealdenssh (cgrid%mmean_can_prss(ipy)               &
                                                  ,cgrid%mmean_can_temp(ipy)               &
                                                  ,cgrid%mmean_can_shv (ipy) )

      !---- Finding AGB and basal area per PFT --------------------------------------------!
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
      !    Finding a few other variables that are either updated every month, or that      !
      ! depend on site-/patch-/cohort- level variables.                                    !
      !------------------------------------------------------------------------------------!
      !----- Flushing the PFT, LU, AGE, and Size (DBH) variables to zero. -----------------!
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

            !------------------------------------------------------------------------------!
            !     Determining whether this is an agricultural patch or not.  Age and size  !
            ! distribution is done only for primary and secondary vegetation.              !
            !------------------------------------------------------------------------------!
            forest = csite%dist_type(ipa) /= 1

            cpatch => csite%patch(ipa)
            cohortloop: do ico = 1, cpatch%ncohorts
               !----- Finding the carbon fluxes. ------------------------------------------!
               cpatch%mmean_gpp         (ico) = cpatch%mmean_gpp         (ico) * ndaysi
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
               cpatch%mmean_par_v       (ico) = cpatch%mmean_par_v       (ico) * ndaysi
               cpatch%mmean_par_v_beam  (ico) = cpatch%mmean_par_v_beam  (ico) * ndaysi
               cpatch%mmean_par_v_diff  (ico) = cpatch%mmean_par_v_diff  (ico) * ndaysi
               cpatch%mmean_leaf_maintenance (ico) = cpatch%mmean_leaf_maintenance(ico)    &
                                                   * ndaysi
               cpatch%mmean_root_maintenance (ico) = cpatch%mmean_root_maintenance(ico)    &
                                                   * ndaysi
               cpatch%mmean_leaf_drop   (ico) = cpatch%mmean_leaf_drop   (ico) * ndaysi
               !----- Mean carbon balance is re-scaled so it will be in kgC/plant/yr. -----!
               cpatch%mmean_cb          (ico) = cpatch%cb(lmon,ico) * ndaysi * yr_day

               !----- Finding the mortality rates. ----------------------------------------!
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
               cpatch%mmean_norm_par_beam(ico)     = cpatch%mmean_norm_par_beam(ico)       &
                                                   * ndaysi
               cpatch%mmean_norm_par_diff(ico)     = cpatch%mmean_norm_par_diff(ico)       &
                                                   * ndaysi
               cpatch%mmean_lambda_light(ico)      = cpatch%mmean_lambda_light(ico)        &
                                                   * ndaysi

               !----- Defining to which PFT this cohort belongs. --------------------------!
               ipft = cpatch%pft(ico)

               !----- Computing the total seed mass of this cohort. -----------------------!
               cohort_seeds   = cpatch%nplant(ico) * cpatch%bseeds(ico)
               
               pss_bseeds_pft(ipft) = pss_bseeds_pft(ipft)                                 &
                                         + cohort_seeds * csite%area(ipa)

               if (forest) then
                  !----- Defining to which size (DBH) class this cohort belongs. ----------!
                  idbh = max(1,min(n_dbh,ceiling(cpatch%dbh(ico)*ddbhi)))

                  !----- Incrementing the plant density. ----------------------------------!
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
      cgrid%mmean_evap               (ipy) = 0.
      cgrid%mmean_transp             (ipy) = 0.
      cgrid%mmean_vapor_ac           (ipy) = 0.
      cgrid%mmean_vapor_gc           (ipy) = 0.
      cgrid%mmean_vapor_vc           (ipy) = 0.
      cgrid%mmean_sensible_ac        (ipy) = 0.
      cgrid%mmean_sensible_gc        (ipy) = 0.
      cgrid%mmean_sensible_vc        (ipy) = 0.
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
      cgrid%mmean_gpp_dbh          (:,ipy) = 0.
      cgrid%mmean_veg_energy         (ipy) = 0.
      cgrid%mmean_veg_hcap           (ipy) = 0.
      cgrid%mmean_veg_water          (ipy) = 0.
      cgrid%mmean_veg_temp           (ipy) = 0.
      cgrid%mmean_can_theta          (ipy) = 0.
      cgrid%mmean_can_theiv          (ipy) = 0.
      cgrid%mmean_can_prss           (ipy) = 0.
      cgrid%mmean_can_temp           (ipy) = 0.
      cgrid%mmean_can_shv            (ipy) = 0.
      cgrid%mmean_can_co2            (ipy) = 0.
      cgrid%mmean_can_rhos           (ipy) = 0.
      cgrid%mmean_atm_temp           (ipy) = 0.
      cgrid%mmean_rshort         (ipy) = 0.
      cgrid%mmean_rlong          (ipy) = 0.
      cgrid%mmean_atm_shv            (ipy) = 0.
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
      cgrid%stdev_gpp                (ipy) = 0.
      cgrid%stdev_evap               (ipy) = 0.
      cgrid%stdev_transp             (ipy) = 0.
      cgrid%stdev_sensible           (ipy) = 0.
      cgrid%stdev_nep                (ipy) = 0.
      cgrid%stdev_rh                 (ipy) = 0.
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

            cpatch=> csite%patch(ipa)
            do ico=1,cpatch%ncohorts
               cpatch%mmean_par_v             (ico) = 0.
               cpatch%mmean_par_v_beam        (ico) = 0.
               cpatch%mmean_par_v_diff        (ico) = 0.
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
               cpatch%mmean_norm_par_beam     (ico) = 0.
               cpatch%mmean_norm_par_diff     (ico) = 0.
               cpatch%mmean_lambda_light      (ico) = 0.
               cpatch%mmean_mort_rate       (:,ico) = 0.
            end do
         end do
      end do
      
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
