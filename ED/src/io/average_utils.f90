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
         !---- Site-level averages.  Pressure is particularly important for enthalpy. -----!
         cpoly%avg_atm_tmp(isi)        = cpoly%avg_atm_tmp(isi)                            &
                                       + cpoly%met(isi)%atm_tmp  * tfact
         cpoly%avg_atm_shv(isi)        = cpoly%avg_atm_shv(isi)                            &
                                       + cpoly%met(isi)%atm_shv  * tfact
         cpoly%avg_prss(isi)           = cpoly%avg_prss(isi)                               &
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

         cgrid%avg_prss(ipy)           = cgrid%avg_prss(ipy)                               &
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
   use ed_misc_coms, only: radfrq
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
            csite%avg_vapor_ac(ipa)     = csite%avg_vapor_ac(ipa)      * frqsumi
            csite%avg_transp(ipa)       = csite%avg_transp(ipa)        * frqsumi
            csite%avg_evap(ipa)         = csite%avg_evap(ipa)          * frqsumi
            csite%avg_runoff(ipa)       = csite%avg_runoff(ipa)        * frqsumi
            csite%avg_drainage(ipa)     = csite%avg_drainage(ipa)      * frqsumi
            csite%avg_sensible_vc(ipa)  = csite%avg_sensible_vc(ipa)   * frqsumi
            csite%avg_qwshed_vg(ipa)    = csite%avg_qwshed_vg(ipa)     * frqsumi
            csite%avg_sensible_gc(ipa)  = csite%avg_sensible_gc(ipa)   * frqsumi
            csite%avg_sensible_ac(ipa)  = csite%avg_sensible_ac(ipa)   * frqsumi
            csite%avg_carbon_ac(ipa)    = csite%avg_carbon_ac(ipa)     * frqsumi
            csite%avg_runoff_heat(ipa)  = csite%avg_runoff_heat(ipa)   * frqsumi
         
            do k=cpoly%lsl(isi),nzg
               csite%avg_sensible_gg(k,ipa) = csite%avg_sensible_gg(k,ipa) * frqsumi
               csite%avg_smoist_gg(k,ipa)   = csite%avg_smoist_gg(k,ipa)   * frqsumi
               csite%avg_smoist_gc(k,ipa)   = csite%avg_smoist_gc(k,ipa)   * frqsumi
               csite%aux_s(k,ipa)           = csite%aux_s(k,ipa)           * frqsumi
            end do
            
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
               ! units: [umol/m2/s * steps] -> [umol/m2/s] 
               cpatch%mean_leaf_resp(ico) = cpatch%mean_leaf_resp(ico) * tfact
               cpatch%mean_root_resp(ico) = cpatch%mean_root_resp(ico) * tfact
               cpatch%mean_gpp(ico)       = cpatch%mean_gpp(ico)       * tfact
               
            end do
            
            !------ Budget variables. -----------------------------------------------------!
            csite%co2budget_gpp(ipa)         = csite%co2budget_gpp(ipa)         * frqsumi
            csite%co2budget_gpp_dbh(:,ipa)   = csite%co2budget_gpp_dbh(:,ipa)   * frqsumi
            csite%co2budget_plresp(ipa)      = csite%co2budget_plresp(ipa)      * frqsumi
            csite%co2budget_rh(ipa)          = csite%co2budget_rh(ipa)          * frqsumi
            csite%co2budget_loss2atm(ipa)    = csite%co2budget_loss2atm(ipa)    * frqsumi
            csite%co2budget_residual(ipa)    = csite%co2budget_residual(ipa)    * frqsumi
            csite%ebudget_precipgain(ipa)    = csite%ebudget_precipgain(ipa)    * frqsumi
            csite%ebudget_netrad(ipa)        = csite%ebudget_netrad(ipa)        * frqsumi
            csite%ebudget_loss2atm(ipa)      = csite%ebudget_loss2atm(ipa)      * frqsumi
            csite%ebudget_loss2drainage(ipa) = csite%ebudget_loss2drainage(ipa) * frqsumi
            csite%ebudget_loss2runoff(ipa)   = csite%ebudget_loss2runoff(ipa)   * frqsumi
            csite%ebudget_residual(ipa)      = csite%ebudget_residual(ipa)      * frqsumi
            csite%wbudget_precipgain(ipa)    = csite%wbudget_precipgain(ipa)    * frqsumi
            csite%wbudget_loss2atm(ipa)      = csite%wbudget_loss2atm(ipa)      * frqsumi
            csite%wbudget_loss2drainage(ipa) = csite%wbudget_loss2drainage(ipa) * frqsumi
            csite%wbudget_loss2runoff(ipa)   = csite%wbudget_loss2runoff(ipa)   * frqsumi
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
      cgrid%avg_prss(ipy)           = 0.0
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
         cpoly%avg_prss(isi)             = 0.0

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

            !----------------------------------------------------------------!
            ! Zeroing water budget variables.                                !
            !----------------------------------------------------------------!
            csite%wbudget_precipgain(ipa)       = 0.0
            csite%wbudget_loss2atm(ipa)         = 0.0
            csite%wbudget_loss2runoff(ipa)      = 0.0
            csite%wbudget_loss2drainage(ipa)    = 0.0
            csite%wbudget_residual(ipa)         = 0.0


            !----------------------------------------------------------------!
            ! Zeroing energy budget variables.                               !
            !----------------------------------------------------------------!
            csite%ebudget_precipgain(ipa)       = 0.0
            csite%ebudget_netrad(ipa)           = 0.0
            csite%ebudget_loss2atm(ipa)         = 0.0
            csite%ebudget_loss2runoff(ipa)      = 0.0
            csite%ebudget_loss2drainage(ipa)    = 0.0
            csite%ebudget_residual(ipa)         = 0.0
            !----------------------------------------------------------------!

            csite%avg_carbon_ac(ipa)        = 0.0
            csite%avg_vapor_vc(ipa)         = 0.0
            csite%avg_dew_cg(ipa)           = 0.0
            csite%avg_vapor_gc(ipa)         = 0.0
            csite%avg_wshed_vg(ipa)         = 0.0
            csite%avg_vapor_ac(ipa)         = 0.0
            csite%avg_transp(ipa)           = 0.0
            csite%avg_evap(ipa)             = 0.0
            csite%avg_netrad(ipa)           = 0.0
            csite%avg_smoist_gg(:,ipa)      = 0.0
            csite%avg_smoist_gc(:,ipa)      = 0.0
            csite%avg_runoff(ipa)           = 0.0
            csite%avg_runoff_heat(ipa)      = 0.0
            csite%avg_drainage(ipa)         = 0.0
            csite%avg_sensible_vc(ipa)      = 0.0
            csite%avg_qwshed_vg(ipa)        = 0.0
            csite%avg_sensible_gc(ipa)      = 0.0
            csite%avg_sensible_ac(ipa)      = 0.0
            csite%avg_sensible_gg(:,ipa)    = 0.0
            csite%avg_runoff_heat(ipa)      = 0.0
            csite%aux(ipa)                  = 0.0
            csite%aux_s(:,ipa)              = 0.0
         
            do ico=1,cpatch%ncohorts
               cpatch%leaf_respiration(ico)      = 0.0
               cpatch%root_respiration(ico)      = 0.0
               cpatch%gpp(ico)                   = 0.0
               cpatch%mean_leaf_resp(ico)        = 0.0
               cpatch%mean_root_resp(ico)        = 0.0
               cpatch%mean_gpp(ico)              = 0.0
            end do
         end do


      enddo

      cpoly%avg_soil_temp(:,:)      = 0.0
      cpoly%avg_soil_water(:,:)     = 0.0
      cpoly%avg_soil_energy(:,:)    = 0.0
      cpoly%avg_soil_fracliq(:,:)   = 0.0

   enddo


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
   use ed_state_vars        , only : edtype        & ! structure
                                   , polygontype   & ! structure
                                   , sitetype      & ! structure
                                   , patchtype     ! ! structure
   use grid_coms            , only : nzg           ! ! intent(in)
   use ed_max_dims             , only : n_dbh         & ! intent(in)
                                   , n_pft         & ! intent(in)
                                   , n_dist_types  ! ! structure
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy, isi, ipa
   logical                     :: any_solvable
   real                        :: poly_area_i,site_area_i, patch_lai_i
   real                        :: forest_site,forest_site_i, forest_poly
   real                        :: sss_fsn, sss_fsw, pss_fsn, pss_fsw
   real                        :: sss_can_enthalpy, sss_can_shv, sss_can_co2, sss_can_rhos
   real                        :: pss_veg_water, pss_veg_energy, pss_veg_hcap
   real                        :: sss_veg_water, sss_veg_energy, sss_veg_hcap
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      !----- Initialize auxiliary variables to add sitetype variables. --------------------!
      sss_fsn          = 0.
      sss_fsw          = 0.
      sss_veg_energy   = 0.
      sss_veg_water    = 0.
      sss_veg_hcap     = 0.
      sss_can_enthalpy = 0.
      sss_can_shv      = 0.
      sss_can_co2      = 0.
      sss_can_rhos     = 0.
      forest_poly      = 0.
      
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
         pss_fsn          = 0.
         pss_fsw          = 0.
         
         pss_veg_energy        = 0.
         pss_veg_water         = 0.
         pss_veg_hcap          = 0.

         !----- Looping through the patches to normalize the sum of all cohorts. ----------!
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)

            any_solvable = .false.
            if (cpatch%ncohorts > 0) then
               any_solvable = any(cpatch%solvable(1:cpatch%ncohorts))

               pss_veg_energy = pss_veg_energy + sum(cpatch%veg_energy) * csite%area(ipa)
               pss_veg_water  = pss_veg_water  + sum(cpatch%veg_water ) * csite%area(ipa)
               pss_veg_hcap   = pss_veg_hcap   + sum(cpatch%hcapveg   ) * csite%area(ipa)
            end if

            if (any_solvable) then
               patch_lai_i = 1./max(tiny(1.),sum(cpatch%lai,cpatch%solvable))

               pss_fsn = pss_fsn + csite%area(ipa)                                         &
                       * (sum(cpatch%fsn * cpatch%lai,cpatch%solvable) * patch_lai_i)
               pss_fsw = pss_fsw + csite%area(ipa)                                         &
                       * (sum(cpatch%fsw * cpatch%lai,cpatch%solvable) * patch_lai_i)
            end if
            
         end do patchloop

         !---------------------------------------------------------------------------------!
         !    Variables already average at the sitetype level, just add them to polygon-   !
         ! type level.                                                                     !
         !---------------------------------------------------------------------------------!
         sss_fsn          = sss_fsn        + (pss_fsn       *site_area_i) * cpoly%area(isi)
         sss_fsw          = sss_fsw        + (pss_fsw       *site_area_i) * cpoly%area(isi)
         sss_veg_energy   = sss_veg_energy + (pss_veg_energy*site_area_i) * cpoly%area(isi)
         sss_veg_water    = sss_veg_water  + (pss_veg_water *site_area_i) * cpoly%area(isi)
         sss_veg_hcap     = sss_veg_hcap   + (pss_veg_hcap  *site_area_i) * cpoly%area(isi)
         sss_can_enthalpy = sss_can_enthalpy                                               &
                          + cpoly%area(isi) * ( sum(csite%can_enthalpy * csite%area)       &
                                              * site_area_i)
         sss_can_shv    = sss_can_shv                                                      &
                        + cpoly%area(isi) * (sum(csite%can_shv * csite%area)  * site_area_i)
         sss_can_co2    = sss_can_co2                                                      &
                        + cpoly%area(isi) * (sum(csite%can_co2 * csite%area)  * site_area_i)
         sss_can_rhos   = sss_can_rhos                                                     &
                        + cpoly%area(isi) * (sum(csite%can_rhos * csite%area) * site_area_i)
      end do siteloop
      
      !------------------------------------------------------------------------------------!
      !    Variables already averaged at the polygontype level, just add them to edtype    !
      ! level.                                                                             !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_fsn(ipy)          = cgrid%dmean_fsn(ipy) + sss_fsn * poly_area_i
      cgrid%dmean_fsw(ipy)          = cgrid%dmean_fsw(ipy) + sss_fsw * poly_area_i
      cgrid%dmean_veg_energy(ipy)   = cgrid%dmean_veg_energy(ipy)                          &
                                    + sss_veg_energy * poly_area_i
      cgrid%dmean_veg_water(ipy)    = cgrid%dmean_veg_water(ipy)                           &
                                    + sss_veg_water * poly_area_i
      cgrid%dmean_veg_hcap(ipy)     = cgrid%dmean_veg_hcap(ipy)                            &
                                    + sss_veg_hcap  * poly_area_i
      cgrid%dmean_can_enthalpy(ipy) = cgrid%dmean_can_enthalpy(ipy)                        &
                                    + sss_can_enthalpy  * poly_area_i
      cgrid%dmean_can_shv(ipy)      = cgrid%dmean_can_shv(ipy)                             &
                                    + sss_can_shv   * poly_area_i
      cgrid%dmean_can_co2(ipy)      = cgrid%dmean_can_co2(ipy)                             &
                                    + sss_can_co2   * poly_area_i
      cgrid%dmean_can_rhos(ipy)     = cgrid%dmean_can_rhos(ipy)                            &
                                    + sss_can_rhos  * poly_area_i

      !------------------------------------------------------------------------------------!
      !    Variables already at edtype level, simple integration only.                     !
      !------------------------------------------------------------------------------------!
      cgrid%dmean_atm_temp(ipy) = cgrid%dmean_atm_temp(ipy) + cgrid%met(ipy)%atm_tmp
      cgrid%dmean_atm_shv(ipy)  = cgrid%dmean_atm_shv(ipy)  + cgrid%met(ipy)%atm_shv
      cgrid%dmean_atm_prss(ipy) = cgrid%dmean_atm_prss(ipy) + cgrid%met(ipy)%prss
      cgrid%dmean_atm_vels(ipy) = cgrid%dmean_atm_vels(ipy) + cgrid%met(ipy)%vels


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
   use ed_state_vars        , only : edtype,polygontype,sitetype,patchtype
   use grid_coms            , only : nzg
   use ed_max_dims             , only : n_dbh, n_pft, n_dist_types
   use ed_misc_coms            , only : dtlsm
   implicit none
   
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy, isi, ipa, k,lu,dbh
   real                        :: poly_area_i,site_area_i
   real                        :: forest_site,forest_site_i, forest_poly
   real, parameter             :: ddbh=1./real(n_dbh-1)
   ! These are auxiliary variables for averaging sitetype variables
   !----------------------------------------------------------------------------------------
   real :: sitesum_leaf_resp, sitesum_root_resp 
   real :: sitesum_plresp,sitesum_gpp,sitesum_rh
   real :: sitesum_evap,sitesum_transp,sitesum_sensible_tot
   real :: sitesum_co2_residual,sitesum_energy_residual,sitesum_water_residual
   real, dimension(n_dist_types) :: sitesum_gpp_lu, sitesum_rh_lu,sitesum_nep_lu
   real, dimension(n_dbh)        :: sitesum_gpp_dbh
   ! These are auxiliary variables for averaging patchtype variables
   !----------------------------------------------------------------------------------------
   real :: patchsum_leaf_resp   , patchsum_root_resp
   
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      ! Initialize auxiliary variables to add sitetype variables
      sitesum_rh              = 0.
      sitesum_leaf_resp       = 0.
      sitesum_root_resp       = 0.
      sitesum_plresp          = 0.
      sitesum_gpp             = 0.
      sitesum_rh              = 0.
      sitesum_gpp_lu          = 0.
      sitesum_nep_lu          = 0.
      sitesum_rh_lu           = 0.
      sitesum_gpp_dbh         = 0.
      sitesum_evap            = 0.
      sitesum_transp          = 0.
      sitesum_sensible_tot    = 0.
      sitesum_co2_residual    = 0.
      sitesum_water_residual  = 0.
      sitesum_energy_residual = 0.
      
      forest_poly          = 0.
      
      siteloop: do isi=1, cpoly%nsites
         csite => cpoly%site(isi)
         
         ! Inverse of total site area (sum of all patches' area)
         site_area_i=1./sum(csite%area)
         ! Forest areas
         forest_site           = sum(csite%area,csite%dist_type /= 1)
         if (forest_site > 1.0e-6) then
            forest_site_i         = 1./forest_site
         else
            forest_site_i         = 0.0
         end if
         forest_poly           = forest_poly + forest_site

         ! Initialize auxiliary variables to add patchtype variables
         patchsum_leaf_resp    = 0.
         patchsum_root_resp    = 0.

         ! Looping through the patches to normalize the sum of all cohorts.
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)
            if (cpatch%ncohorts > 0) then
               patchsum_leaf_resp = patchsum_leaf_resp + sum(cpatch%mean_leaf_resp, cpatch%solvable)  * csite%area(ipa) 
               patchsum_root_resp = patchsum_root_resp + sum(cpatch%mean_root_resp                 )  * csite%area(ipa) 
            end if
            csite%dmean_co2_residual(ipa)    = csite%dmean_co2_residual(ipa)    + csite%co2budget_residual(ipa)
            csite%dmean_energy_residual(ipa) = csite%dmean_energy_residual(ipa) + csite%ebudget_residual(ipa)
            csite%dmean_water_residual(ipa)  = csite%dmean_water_residual(ipa)  + csite%wbudget_residual(ipa)
         end do patchloop
         
         ! Variables already average at the sitetype level, just add them to polygontype level
         sitesum_leaf_resp    = sitesum_leaf_resp    + (patchsum_leaf_resp      *site_area_i) * cpoly%area(isi)
         sitesum_root_resp    = sitesum_root_resp    + (patchsum_root_resp      *site_area_i) * cpoly%area(isi)

         sitesum_rh           = sitesum_rh           + (sum(csite%co2budget_rh    *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_gpp          = sitesum_gpp          + (sum(csite%co2budget_gpp   *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_plresp       = sitesum_plresp       + (sum(csite%co2budget_plresp*csite%area) * site_area_i) * cpoly%area(isi)
         
         sitesum_evap         = sitesum_evap         + (sum(csite%avg_evap        *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_transp       = sitesum_transp       + (sum(csite%avg_transp      *csite%area) * site_area_i) * cpoly%area(isi)
         
         sitesum_co2_residual    = sitesum_co2_residual    + (sum(csite%co2budget_residual*csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_water_residual  = sitesum_water_residual  + (sum(csite%wbudget_residual  *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_energy_residual = sitesum_energy_residual + (sum(csite%ebudget_residual  *csite%area) * site_area_i) * cpoly%area(isi)

         cpoly%dmean_co2_residual(isi)    = cpoly%dmean_co2_residual(isi)    + sum(csite%co2budget_residual*csite%area) * site_area_i
         cpoly%dmean_energy_residual(isi) = cpoly%dmean_energy_residual(isi) + sum(csite%ebudget_residual  *csite%area) * site_area_i
         cpoly%dmean_water_residual(isi)  = cpoly%dmean_water_residual(isi)  + sum(csite%wbudget_residual  *csite%area) * site_area_i

         luloop: do lu=1,n_dist_types
            sitesum_rh_lu(lu)    = sitesum_rh_lu(lu) + &
                                   (sum(csite%co2budget_rh*csite%area,csite%dist_type == lu)   * site_area_i) * cpoly%area(isi)
            sitesum_gpp_lu(lu)   = sitesum_gpp_lu(lu) + &
                                   (sum(csite%co2budget_gpp*csite%area,csite%dist_type == lu)  * site_area_i) * cpoly%area(isi)
            sitesum_nep_lu(lu)   = sitesum_nep_lu(lu) + &
                                   (sum((csite%co2budget_gpp-csite%co2budget_rh-csite%co2budget_plresp)*csite%area &
                                        ,csite%dist_type == lu) * site_area_i) * cpoly%area(isi)
         end do luloop
         
         dbhloop: do dbh=1,n_dbh
            sitesum_gpp_dbh(dbh) = sitesum_gpp_dbh(dbh) + &
                 (sum(csite%co2budget_gpp_dbh(dbh,:)*csite%area) * site_area_i) * cpoly%area(isi)
         end do dbhloop

      end do siteloop
      
      cgrid%dmean_leaf_resp(ipy)    = cgrid%dmean_leaf_resp(ipy)    + sitesum_leaf_resp    * poly_area_i
      cgrid%dmean_root_resp(ipy)    = cgrid%dmean_root_resp(ipy)    + sitesum_root_resp    * poly_area_i
      
      cgrid%dmean_rh(ipy)           = cgrid%dmean_rh(ipy)           + sitesum_rh           * poly_area_i
      cgrid%dmean_gpp(ipy)          = cgrid%dmean_gpp(ipy)          + sitesum_gpp          * poly_area_i
      cgrid%dmean_plresp(ipy)       = cgrid%dmean_plresp(ipy)       + sitesum_plresp       * poly_area_i
      cgrid%dmean_nep(ipy)          = cgrid%dmean_nep(ipy)          + (sitesum_gpp-sitesum_rh-sitesum_plresp)     &
                                                                                           * poly_area_i

      cgrid%dmean_co2_residual(ipy)    = cgrid%dmean_co2_residual(ipy)    + sitesum_co2_residual    * poly_area_i
      cgrid%dmean_energy_residual(ipy) = cgrid%dmean_energy_residual(ipy) + sitesum_energy_residual  * poly_area_i
      cgrid%dmean_water_residual(ipy)  = cgrid%dmean_water_residual(ipy)  + sitesum_water_residual * poly_area_i

      do lu=1,n_dist_types
         cgrid%dmean_gpp_lu(lu,ipy) = cgrid%dmean_gpp_lu(lu,ipy)    + sitesum_gpp_lu(lu)   * poly_area_i
         cgrid%dmean_nep_lu(lu,ipy) = cgrid%dmean_nep_lu(lu,ipy)    + sitesum_nep_lu(lu)   * poly_area_i
         cgrid%dmean_rh_lu(lu,ipy)  = cgrid%dmean_rh_lu(lu,ipy)     + sitesum_rh_lu(lu)    * poly_area_i
      end do
      
      do dbh=1,n_dbh
         cgrid%dmean_gpp_dbh(dbh,ipy) = cgrid%dmean_gpp_dbh(dbh,ipy) + sitesum_gpp_dbh(dbh)* poly_area_i
      end do

      !These variables are already averaged at gridtype, just add them up
      do k=1,nzg
         cgrid%dmean_soil_temp(k,ipy)   = cgrid%dmean_soil_temp(k,ipy)   + cgrid%avg_soil_temp(k,ipy) 
         cgrid%dmean_soil_water(k,ipy)  = cgrid%dmean_soil_water(k,ipy)  + cgrid%avg_soil_water(k,ipy)
      end do
      cgrid%dmean_evap(ipy)         = cgrid%dmean_evap(ipy)         + cgrid%avg_evap(ipy)     
      cgrid%dmean_transp(ipy)       = cgrid%dmean_transp(ipy)       + cgrid%avg_transp(ipy)   
      cgrid%dmean_sensible_vc(ipy)  = cgrid%dmean_sensible_vc(ipy)  + cgrid%avg_sensible_vc(ipy) 
      cgrid%dmean_sensible_gc(ipy)  = cgrid%dmean_sensible_gc(ipy)  + cgrid%avg_sensible_gc(ipy) 
      cgrid%dmean_sensible_ac(ipy)  = cgrid%dmean_sensible_ac(ipy)  + cgrid%avg_sensible_ac(ipy) 

      cgrid%dmean_pcpg(ipy)         = cgrid%dmean_pcpg(ipy)         + cgrid%avg_pcpg(ipy)     
      cgrid%dmean_runoff(ipy)       = cgrid%dmean_runoff(ipy)       + cgrid%avg_runoff(ipy)
      cgrid%dmean_drainage(ipy)     = cgrid%dmean_drainage(ipy)     + cgrid%avg_drainage(ipy)
      cgrid%dmean_vapor_vc(ipy)     = cgrid%dmean_vapor_vc(ipy)     + cgrid%avg_vapor_vc(ipy) 
      cgrid%dmean_vapor_gc(ipy)     = cgrid%dmean_vapor_gc(ipy)     + cgrid%avg_vapor_gc(ipy) 
      cgrid%dmean_vapor_ac(ipy)     = cgrid%dmean_vapor_ac(ipy)     + cgrid%avg_vapor_ac(ipy) 


   end do polyloop

   return
end subroutine integrate_ed_daily_output_flux
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine normalize_ed_daily_vars(cgrid,timefac1)
   use ed_state_vars , only : edtype,polygontype,sitetype,patchtype
   use ed_max_dims      , only : n_pft
   implicit none
   real, intent(in)         :: timefac1 ! Daily sum          => daily average
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   
   integer                  :: ipy,isi,ipa,ico
   real, dimension(n_pft)   :: patchsum_lai
   
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         
         csite%dmean_A_decomp  = csite%dmean_A_decomp  * timefac1
         csite%dmean_Af_decomp = csite%dmean_Af_decomp * timefac1
         
         patchsum_lai = 0.
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            !----- Included a loop so it won't crash with empty cohorts... ----------------!
            cohortloop: do ico=1,cpatch%ncohorts
               cpatch%dmean_gpp(ico)       = cpatch%dmean_gpp(ico)       * timefac1
               cpatch%dmean_gpp_pot(ico)   = cpatch%dmean_gpp_pot(ico)   * timefac1
               cpatch%dmean_gpp_max(ico)   = cpatch%dmean_gpp_max(ico)   * timefac1
               cpatch%dmean_leaf_resp(ico) = cpatch%dmean_leaf_resp(ico) * timefac1
               cpatch%dmean_root_resp(ico) = cpatch%dmean_root_resp(ico) * timefac1
            end do cohortloop
         end do patchloop
      end do siteloop
   
   end do polyloop
   
   return
end subroutine normalize_ed_daily_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine normalize_ed_daily_output_vars(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine normalize the sum before writing the daily analysis. It also computes !
! some of the variables that didn't need to be computed every time step, like LAI.         !
!------------------------------------------------------------------------------------------!
   use grid_coms     , only : nzg
   use ed_state_vars , only : edtype,polygontype,sitetype,patchtype
   use ed_max_dims      , only : n_pft,n_dist_types
   use consts_coms   , only : cpi, alvl,day_sec,umol_2_kgC
   use ed_misc_coms     , only : dtlsm,frqsum
   use therm_lib     , only : qwtk,hpq2temp
   implicit none
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                  :: ipy,isi,ipa,ipft,ilu,k
   real                     :: polygon_area_i, site_area_i
   real :: sitesum_storage_resp, sitesum_vleaf_resp, sitesum_growth_resp
   real :: patchsum_storage_resp, patchsum_vleaf_resp , patchsum_growth_resp
   real :: veg_fliq
   
   logical           , save :: find_factors=.true.
   real              , save :: dtlsm_o_daysec=1.E34, frqsum_o_daysec=1.E34

   !!! RESET area indices
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      cgrid%lai_pft            (:,ipy) = 0.
      cgrid%lai_lu             (:,ipy) = 0.
      cgrid%wpa_pft            (:,ipy) = 0.
      cgrid%wpa_lu             (:,ipy) = 0.
      cgrid%wai_pft            (:,ipy) = 0.
      cgrid%wai_lu             (:,ipy) = 0.
      do isi=1,cpoly%nsites
         cpoly%lai_pft (:,isi)  = 0.
         cpoly%lai_lu  (:,isi)  = 0.
         cpoly%wpa_pft (:,isi)  = 0.
         cpoly%wpa_lu  (:,isi)  = 0.
         cpoly%wai_pft (:,isi)  = 0.
         cpoly%wai_lu  (:,isi)  = 0.
      end do
   end do

   
   ! Computing the normalization factors. This is done once.
   if (find_factors) then
      dtlsm_o_daysec  = dtlsm/day_sec
      frqsum_o_daysec = frqsum/day_sec
      find_factors    = .false.
   end if
   

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      ! State variables, updated every time step, so these are normalized by dtlsm/day_sec
      cgrid%dmean_fsw(ipy)          = cgrid%dmean_fsw(ipy)          * dtlsm_o_daysec
      cgrid%dmean_fsn(ipy)          = cgrid%dmean_fsn(ipy)          * dtlsm_o_daysec
      cgrid%dmean_veg_energy(ipy)   = cgrid%dmean_veg_energy(ipy)   * dtlsm_o_daysec
      cgrid%dmean_veg_hcap(ipy)     = cgrid%dmean_veg_hcap(ipy)     * dtlsm_o_daysec
      cgrid%dmean_veg_water(ipy)    = cgrid%dmean_veg_water(ipy)    * dtlsm_o_daysec
      cgrid%dmean_can_enthalpy(ipy) = cgrid%dmean_can_enthalpy(ipy) * dtlsm_o_daysec
      cgrid%dmean_can_shv(ipy)      = cgrid%dmean_can_shv(ipy)      * dtlsm_o_daysec
      cgrid%dmean_can_co2(ipy)      = cgrid%dmean_can_co2(ipy)      * dtlsm_o_daysec
      cgrid%dmean_can_rhos(ipy)     = cgrid%dmean_can_rhos(ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_temp(ipy)     = cgrid%dmean_atm_temp(ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_shv(ipy)      = cgrid%dmean_atm_shv(ipy)      * dtlsm_o_daysec
      cgrid%dmean_atm_prss(ipy)     = cgrid%dmean_atm_prss(ipy)     * dtlsm_o_daysec
      cgrid%dmean_atm_vels(ipy)     = cgrid%dmean_atm_vels(ipy)     * dtlsm_o_daysec

      !----- Finding canopy temperature. --------------------------------------------------!
      cgrid%dmean_can_temp(ipy)     = hpq2temp(cgrid%dmean_can_enthalpy(ipy)               &
                                              ,cgrid%dmean_atm_prss(ipy)                   &
                                              ,cgrid%dmean_can_shv(ipy))

      !----- Finding vegetation temperature -----------------------------------------------!
      call qwtk(cgrid%dmean_veg_energy(ipy),cgrid%dmean_veg_water(ipy)                     &
               ,cgrid%dmean_veg_hcap(ipy),cgrid%dmean_veg_temp(ipy),veg_fliq)
      
      ! State variables, updated every frqsum, so these are normalized by frqsum/day_sec
      do k=1,nzg
         cgrid%dmean_soil_temp (k,ipy) = cgrid%dmean_soil_temp (k,ipy)   * frqsum_o_daysec
         cgrid%dmean_soil_water(k,ipy) = cgrid%dmean_soil_water(k,ipy)   * frqsum_o_daysec
      end do
      ! Precip and runoff
      cgrid%dmean_pcpg       (ipy)  = cgrid%dmean_pcpg       (ipy)  * frqsum_o_daysec ! kg/m2/sec
      cgrid%dmean_runoff     (ipy)  = cgrid%dmean_runoff     (ipy)  * frqsum_o_daysec ! kg/m2/sec
      cgrid%dmean_drainage   (ipy)  = cgrid%dmean_drainage   (ipy)  * frqsum_o_daysec ! kg/m2/sec

      ! vapor flux
      cgrid%dmean_vapor_vc(ipy)  = cgrid%dmean_vapor_vc(ipy)  * frqsum_o_daysec
      cgrid%dmean_vapor_gc(ipy)  = cgrid%dmean_vapor_gc(ipy)  * frqsum_o_daysec
      cgrid%dmean_vapor_ac(ipy)  = cgrid%dmean_vapor_ac(ipy)  * frqsum_o_daysec


      ! Flux variables, updated every frqsum, so these are normalized by frqsum/day_sec
      cgrid%dmean_evap       (ipy)  = cgrid%dmean_evap       (ipy)  * frqsum_o_daysec
      cgrid%dmean_transp     (ipy)  = cgrid%dmean_transp     (ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_vc(ipy)  = cgrid%dmean_sensible_vc(ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_gc(ipy)  = cgrid%dmean_sensible_gc(ipy)  * frqsum_o_daysec
      cgrid%dmean_sensible_ac(ipy)  = cgrid%dmean_sensible_ac(ipy)  * frqsum_o_daysec


      ! Carbon flux variables should be total flux integrated over the day at this point, [µmol/m²/s] -> [kgC/m²/day]
      cgrid%dmean_leaf_resp  (ipy)  = cgrid%dmean_leaf_resp (ipy)   * umol_2_kgC * day_sec
      cgrid%dmean_root_resp  (ipy)  = cgrid%dmean_root_resp (ipy)   * umol_2_kgC * day_sec

      ! Carbon flux variables should be total flux integrated over the day at this point, [µmol/m²/day] -> [kgC/m²/day]
      cgrid%dmean_rh         (ipy)  = cgrid%dmean_rh        (ipy)   * umol_2_kgC
      cgrid%dmean_gpp        (ipy)  = cgrid%dmean_gpp       (ipy)   * umol_2_kgC
      cgrid%dmean_plresp     (ipy)  = cgrid%dmean_plresp    (ipy)   * umol_2_kgC
      cgrid%dmean_nep        (ipy)  = cgrid%dmean_nep       (ipy)   * umol_2_kgC
   
      cgrid%dmean_gpp_lu     (:,ipy)  = cgrid%dmean_gpp_lu    (:,ipy)   * umol_2_kgC
      cgrid%dmean_rh_lu      (:,ipy)  = cgrid%dmean_rh_lu     (:,ipy)   * umol_2_kgC
      cgrid%dmean_nep_lu     (:,ipy)  = cgrid%dmean_nep_lu    (:,ipy)   * umol_2_kgC
      cgrid%dmean_gpp_dbh    (:,ipy)  = cgrid%dmean_gpp_dbh   (:,ipy)   * umol_2_kgC

      cgrid%dmean_co2_residual(ipy)    = cgrid%dmean_co2_residual(ipy)    * frqsum_o_daysec
      cgrid%dmean_energy_residual(ipy) = cgrid%dmean_energy_residual(ipy) * frqsum_o_daysec
      cgrid%dmean_water_residual(ipy)  = cgrid%dmean_water_residual(ipy)  * frqsum_o_daysec


      polygon_area_i = 1./sum(cpoly%area)

      sitesum_growth_resp  = 0.
      sitesum_storage_resp = 0.
      sitesum_vleaf_resp   = 0.

      
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)

         cpoly%dmean_co2_residual(isi)    = cpoly%dmean_co2_residual(isi)    * frqsum_o_daysec
         cpoly%dmean_energy_residual(isi) = cpoly%dmean_energy_residual(isi) * frqsum_o_daysec
         cpoly%dmean_water_residual(isi)  = cpoly%dmean_water_residual(isi)  * frqsum_o_daysec
         

         site_area_i = 1./sum(csite%area)
         
         patchsum_growth_resp  = 0.
         patchsum_storage_resp = 0.
         patchsum_vleaf_resp   = 0.

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            csite%dmean_co2_residual(ipa)    = csite%dmean_co2_residual(ipa)    * frqsum_o_daysec
            csite%dmean_energy_residual(ipa) = csite%dmean_energy_residual(ipa) * frqsum_o_daysec
            csite%dmean_water_residual(ipa)  = csite%dmean_water_residual(ipa)  * frqsum_o_daysec
            
            if (cpatch%ncohorts > 0) then
               ! These are in kg/plant/day, converting to kg/m²/day
               patchsum_growth_resp  = patchsum_growth_resp + csite%area(ipa)              &
                                     * sum(cpatch%growth_respiration  * cpatch%nplant)
               patchsum_storage_resp = patchsum_storage_resp + csite%area(ipa)             &
                                     * sum(cpatch%storage_respiration * cpatch%nplant)
               patchsum_vleaf_resp   = patchsum_vleaf_resp   + csite%area(ipa)             &
                                     * sum(cpatch%vleaf_respiration   * cpatch%nplant)
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

            ilu = csite%dist_type(ipa)
            cpoly%lai_lu(ilu,isi)  = cpoly%lai_lu(ilu,isi)                                 &
                                   + csite%lai(ipa) * csite%area(ipa) * site_area_i
            cpoly%wpa_lu(ilu,isi)  = cpoly%wpa_lu(ilu,isi)                                 &
                                   + csite%wpa(ipa) * csite%area(ipa) * site_area_i
            cpoly%wai_lu(ilu,isi)  = cpoly%wai_lu(ilu,isi)                                 &
                                   + csite%wai(ipa) * csite%area(ipa) * site_area_i

         end do patchloop

         sitesum_growth_resp  = sitesum_growth_resp  + (patchsum_growth_resp    *site_area_i) * cpoly%area(isi)
         sitesum_storage_resp = sitesum_storage_resp + (patchsum_storage_resp   *site_area_i) * cpoly%area(isi)
         sitesum_vleaf_resp   = sitesum_vleaf_resp   + (patchsum_vleaf_resp     *site_area_i) * cpoly%area(isi)


      end do siteloop
      
      do ipft=1,n_pft
         cgrid%lai_pft(ipft,ipy)  = cgrid%lai_pft(ipft,ipy)                                &
                                  + sum(cpoly%lai_pft(ipft,:)*cpoly%area) * polygon_area_i
         cgrid%wpa_pft(ipft,ipy)  = cgrid%wpa_pft(ipft,ipy)                                &
                                  + sum(cpoly%wpa_pft(ipft,:)*cpoly%area) * polygon_area_i
         cgrid%wai_pft(ipft,ipy)  = cgrid%wai_pft(ipft,ipy)                                &
                                  + sum(cpoly%wai_pft(ipft,:)*cpoly%area) * polygon_area_i
      end do
      
      do ilu=1,n_dist_types
         cgrid%lai_lu(ilu,ipy)  = cgrid%lai_lu(ilu,ipy)                                    &
                                + sum(cpoly%lai_lu(ilu,:)*cpoly%area) * polygon_area_i
         cgrid%wpa_lu(ilu,ipy)  = cgrid%wpa_lu(ilu,ipy)                                    &
                                + sum(cpoly%wpa_lu(ilu,:)*cpoly%area) * polygon_area_i
         cgrid%wai_lu(ilu,ipy)  = cgrid%wai_lu(ilu,ipy)                                    &
                                + sum(cpoly%wai_lu(ilu,:)*cpoly%area) * polygon_area_i
      end do
      
      cgrid%dmean_growth_resp(ipy)  = cgrid%dmean_growth_resp(ipy)  + sitesum_growth_resp  * polygon_area_i
      cgrid%dmean_storage_resp(ipy) = cgrid%dmean_storage_resp(ipy) + sitesum_storage_resp * polygon_area_i
      cgrid%dmean_vleaf_resp(ipy)   = cgrid%dmean_vleaf_resp(ipy)   + sitesum_vleaf_resp   * polygon_area_i

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
   type(edtype)     , target  :: cgrid
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy,isi,ipa,ico
   
   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
            
      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)


         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            
            !--------------------------------------!
            ! Reset variables stored in sitetype   !
            !--------------------------------------!
            csite%dmean_A_decomp(ipa)  = 0.0
            csite%dmean_Af_decomp(ipa) = 0.0

            !-------------------------------------!
            ! Reset variables stored in patchtype !
            !-------------------------------------!
            do ico = 1, cpatch%ncohorts
               cpatch%dmean_gpp      (ico) = 0.0
               cpatch%dmean_gpp_pot  (ico) = 0.0
               cpatch%dmean_gpp_max  (ico) = 0.0
               cpatch%dmean_leaf_resp(ico) = 0.0
               cpatch%dmean_root_resp(ico) = 0.0
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
subroutine zero_ed_daily_output_vars(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine resets the daily_averages once the daily average was written and used !
! to compute the monthly mean (in case the latter was requested).                          !
!------------------------------------------------------------------------------------------!
   use ed_state_vars        , only: edtype,polygontype,sitetype,patchtype
   implicit none
   type(edtype)     , target  :: cgrid
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   integer                    :: ipy,isi,ipa
   

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      
      !-------------------------------------!
      ! Variables stored in edtype          !
      !-------------------------------------!
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
 
      cgrid%dmean_plresp         (ipy) = 0.
      cgrid%dmean_rh             (ipy) = 0.
      cgrid%dmean_leaf_resp      (ipy) = 0.
      cgrid%dmean_root_resp      (ipy) = 0.
      cgrid%dmean_growth_resp    (ipy) = 0.
      cgrid%dmean_storage_resp   (ipy) = 0.
      cgrid%dmean_vleaf_resp     (ipy) = 0.
      cgrid%dmean_nep            (ipy) = 0.
      cgrid%dmean_fsw            (ipy) = 0.
      cgrid%dmean_fsn            (ipy) = 0.
      cgrid%dmean_soil_temp    (:,ipy) = 0.
      cgrid%dmean_soil_water   (:,ipy) = 0.
      cgrid%dmean_gpp_lu       (:,ipy) = 0.
      cgrid%dmean_rh_lu        (:,ipy) = 0.
      cgrid%dmean_nep_lu       (:,ipy) = 0.
      cgrid%dmean_gpp_dbh      (:,ipy) = 0.
      cgrid%dmean_fsw            (ipy) = 0.
      cgrid%dmean_fsn            (ipy) = 0.
      cgrid%dmean_veg_energy     (ipy) = 0.
      cgrid%dmean_veg_hcap       (ipy) = 0.
      cgrid%dmean_veg_water      (ipy) = 0.
      cgrid%dmean_veg_temp       (ipy) = 0.
      cgrid%dmean_can_enthalpy   (ipy) = 0.
      cgrid%dmean_can_temp       (ipy) = 0.
      cgrid%dmean_can_shv        (ipy) = 0.
      cgrid%dmean_can_co2        (ipy) = 0.
      cgrid%dmean_can_rhos       (ipy) = 0.
      cgrid%dmean_atm_temp       (ipy) = 0.
      cgrid%dmean_atm_shv        (ipy) = 0.
      cgrid%dmean_atm_prss       (ipy) = 0.
      cgrid%dmean_atm_vels       (ipy) = 0.
      cgrid%lai_pft            (:,ipy) = 0.
      cgrid%lai_lu             (:,ipy) = 0.
      cgrid%wpa_pft            (:,ipy) = 0.
      cgrid%wpa_lu             (:,ipy) = 0.
      cgrid%wai_pft            (:,ipy) = 0.
      cgrid%wai_lu             (:,ipy) = 0.
      cgrid%dmean_co2_residual   (ipy) = 0.
      cgrid%dmean_energy_residual(ipy) = 0.
      cgrid%dmean_water_residual (ipy) = 0.

      !-----------------------------------------!
      ! Reset variables stored in polygontype   !
      !-----------------------------------------!
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)

         cpoly%lai_pft           (:,isi)  = 0.
         cpoly%lai_lu            (:,isi)  = 0.
         cpoly%wpa_pft           (:,isi)  = 0.
         cpoly%wpa_lu            (:,isi)  = 0.
         cpoly%wai_pft           (:,isi)  = 0.
         cpoly%wai_lu            (:,isi)  = 0.
         cpoly%dmean_co2_residual   (isi) = 0.
         cpoly%dmean_energy_residual(isi) = 0.
         cpoly%dmean_water_residual (isi) = 0.

         do ipa=1,csite%npatches
            csite%dmean_co2_residual   (ipa) = 0.
            csite%dmean_energy_residual(ipa) = 0.
            csite%dmean_water_residual (ipa) = 0.       
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
subroutine integrate_ed_monthly_output_vars(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine integrates the monthly average. This is called after the daily means  !
! were integrated and normalized.                                                          !
!------------------------------------------------------------------------------------------!
   use ed_state_vars, only : edtype       & ! structure
                           , polygontype  & ! structure
                           , sitetype     ! ! structure
   use ed_max_dims  , only : n_dbh        & ! intent(in)
                           , n_pft        & ! intent(in) 
                           , n_dist_types ! ! intent(in)
   implicit none
   !----- Argument. -----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   !---------------------------------------------------------------------------------------!
   
   do ipy=1,cgrid%npolygons
      !------------------------------------------------------------------------------------!
      ! First the mean variables that can be computed from the daily averages              !
      !------------------------------------------------------------------------------------!
      cgrid%mmean_gpp           (ipy) = cgrid%mmean_gpp           (ipy) +  cgrid%dmean_gpp           (ipy)
      cgrid%mmean_evap          (ipy) = cgrid%mmean_evap          (ipy) +  cgrid%dmean_evap          (ipy)
      cgrid%mmean_transp        (ipy) = cgrid%mmean_transp        (ipy) +  cgrid%dmean_transp        (ipy)

      cgrid%mmean_sensible_ac   (ipy) = cgrid%mmean_sensible_ac   (ipy) +  cgrid%dmean_sensible_ac   (ipy)
      cgrid%mmean_sensible_gc   (ipy) = cgrid%mmean_sensible_gc   (ipy) +  cgrid%dmean_sensible_gc   (ipy)
      cgrid%mmean_sensible_vc   (ipy) = cgrid%mmean_sensible_vc   (ipy) +  cgrid%dmean_sensible_vc   (ipy)
      cgrid%mmean_nep           (ipy) = cgrid%mmean_nep           (ipy) +  cgrid%dmean_nep           (ipy)
      cgrid%mmean_plresp        (ipy) = cgrid%mmean_plresp        (ipy) +  cgrid%dmean_plresp        (ipy)
      cgrid%mmean_rh            (ipy) = cgrid%mmean_rh            (ipy) +  cgrid%dmean_rh            (ipy)
      cgrid%mmean_leaf_resp     (ipy) = cgrid%mmean_leaf_resp     (ipy) +  cgrid%dmean_leaf_resp     (ipy)
      cgrid%mmean_root_resp     (ipy) = cgrid%mmean_root_resp     (ipy) +  cgrid%dmean_root_resp     (ipy)
      cgrid%mmean_growth_resp   (ipy) = cgrid%mmean_growth_resp   (ipy) +  cgrid%dmean_growth_resp   (ipy)
      cgrid%mmean_storage_resp  (ipy) = cgrid%mmean_storage_resp  (ipy) +  cgrid%dmean_storage_resp  (ipy)
      cgrid%mmean_vleaf_resp    (ipy) = cgrid%mmean_vleaf_resp    (ipy) +  cgrid%dmean_vleaf_resp    (ipy)

      cgrid%mmean_soil_temp   (:,ipy) = cgrid%mmean_soil_temp   (:,ipy) +  cgrid%dmean_soil_temp   (:,ipy)
      cgrid%mmean_soil_water  (:,ipy) = cgrid%mmean_soil_water  (:,ipy) +  cgrid%dmean_soil_water  (:,ipy)
      cgrid%mmean_gpp_lu      (:,ipy) = cgrid%mmean_gpp_lu      (:,ipy) +  cgrid%dmean_gpp_lu      (:,ipy)
      cgrid%mmean_rh_lu       (:,ipy) = cgrid%mmean_rh_lu       (:,ipy) +  cgrid%dmean_rh_lu       (:,ipy)
      cgrid%mmean_nep_lu      (:,ipy) = cgrid%mmean_nep_lu      (:,ipy) +  cgrid%dmean_nep_lu      (:,ipy)
      cgrid%mmean_gpp_dbh     (:,ipy) = cgrid%mmean_gpp_dbh     (:,ipy) +  cgrid%dmean_gpp_dbh     (:,ipy)
      cgrid%mmean_lai_pft     (:,ipy) = cgrid%mmean_lai_pft     (:,ipy) +  cgrid%lai_pft           (:,ipy)
      cgrid%mmean_lai_lu      (:,ipy) = cgrid%mmean_lai_lu      (:,ipy) +  cgrid%lai_lu            (:,ipy)
      cgrid%mmean_wpa_pft     (:,ipy) = cgrid%mmean_wpa_pft     (:,ipy) +  cgrid%wpa_pft           (:,ipy)
      cgrid%mmean_wpa_lu      (:,ipy) = cgrid%mmean_wpa_lu      (:,ipy) +  cgrid%wpa_lu            (:,ipy)
      cgrid%mmean_wai_pft     (:,ipy) = cgrid%mmean_wai_pft     (:,ipy) +  cgrid%wai_pft           (:,ipy)
      cgrid%mmean_wai_lu      (:,ipy) = cgrid%mmean_wai_lu      (:,ipy) +  cgrid%wai_lu            (:,ipy)

      cgrid%mmean_veg_energy    (ipy) = cgrid%mmean_veg_energy    (ipy) + cgrid%dmean_veg_energy     (ipy)
      cgrid%mmean_veg_hcap      (ipy) = cgrid%mmean_veg_hcap      (ipy) + cgrid%dmean_veg_hcap       (ipy)
      cgrid%mmean_veg_water     (ipy) = cgrid%mmean_veg_water     (ipy) + cgrid%dmean_veg_water      (ipy)
      cgrid%mmean_veg_temp      (ipy) = cgrid%mmean_veg_temp      (ipy) + cgrid%dmean_veg_temp       (ipy)
      cgrid%mmean_can_enthalpy  (ipy) = cgrid%mmean_can_enthalpy  (ipy) + cgrid%dmean_can_enthalpy   (ipy)
      cgrid%mmean_can_temp      (ipy) = cgrid%mmean_can_temp      (ipy) + cgrid%dmean_can_temp       (ipy)
      cgrid%mmean_can_shv       (ipy) = cgrid%mmean_can_shv       (ipy) + cgrid%dmean_can_shv        (ipy)
      cgrid%mmean_can_co2       (ipy) = cgrid%mmean_can_co2       (ipy) + cgrid%dmean_can_co2        (ipy)
      cgrid%mmean_can_rhos      (ipy) = cgrid%mmean_can_rhos      (ipy) + cgrid%dmean_can_rhos       (ipy)
      cgrid%mmean_atm_temp      (ipy) = cgrid%mmean_atm_temp      (ipy) + cgrid%dmean_atm_temp       (ipy)
      cgrid%mmean_atm_shv       (ipy) = cgrid%mmean_atm_shv       (ipy) + cgrid%dmean_atm_shv        (ipy)
      cgrid%mmean_atm_prss      (ipy) = cgrid%mmean_atm_prss      (ipy) + cgrid%dmean_atm_prss       (ipy)
      cgrid%mmean_atm_vels      (ipy) = cgrid%mmean_atm_vels      (ipy) + cgrid%dmean_atm_vels       (ipy)
      cgrid%mmean_pcpg          (ipy) = cgrid%mmean_pcpg          (ipy) + cgrid%dmean_pcpg           (ipy)

      cgrid%mmean_co2_residual(ipy)    = cgrid%mmean_co2_residual(ipy)    + cgrid%dmean_co2_residual(ipy)   
      cgrid%mmean_energy_residual(ipy) = cgrid%mmean_energy_residual(ipy) + cgrid%dmean_energy_residual(ipy)
      cgrid%mmean_water_residual(ipy)  = cgrid%mmean_water_residual(ipy)  + cgrid%dmean_water_residual(ipy) 

      !------------------------------------------------------------------------------------!
      !    During the integration stage we keep the sum of squares, it will be converted   !
      ! to standard deviation right before the monthly output.                             !
      !------------------------------------------------------------------------------------!
      cgrid%stdev_gpp     (ipy) = cgrid%stdev_gpp     (ipy) +  cgrid%dmean_gpp     (ipy) ** 2
      cgrid%stdev_evap    (ipy) = cgrid%stdev_evap    (ipy) +  cgrid%dmean_evap    (ipy) ** 2
      cgrid%stdev_transp  (ipy) = cgrid%stdev_transp  (ipy) +  cgrid%dmean_transp  (ipy) ** 2
      cgrid%stdev_sensible(ipy) = cgrid%stdev_sensible(ipy) +  cgrid%dmean_sensible_ac(ipy) ** 2
      cgrid%stdev_nep     (ipy) = cgrid%stdev_nep     (ipy) +  cgrid%dmean_nep     (ipy) ** 2
      cgrid%stdev_rh      (ipy) = cgrid%stdev_rh      (ipy) +  cgrid%dmean_rh      (ipy) ** 2
      
      cpoly => cgrid%polygon(ipy)
      do isi=1,cpoly%nsites
         cpoly%mmean_co2_residual(isi)    = cpoly%mmean_co2_residual(isi)    + cpoly%dmean_co2_residual(isi)   
         cpoly%mmean_energy_residual(isi) = cpoly%mmean_energy_residual(isi) + cpoly%dmean_energy_residual(isi)
         cpoly%mmean_water_residual(isi)  = cpoly%mmean_water_residual(isi)  + cpoly%dmean_water_residual(isi) 

         csite => cpoly%site(isi)
         do ipa=1,csite%npatches
            csite%mmean_co2_residual(ipa)    = csite%mmean_co2_residual(ipa)    + csite%dmean_co2_residual(ipa)   
            csite%mmean_energy_residual(ipa) = csite%mmean_energy_residual(ipa) + csite%dmean_energy_residual(ipa)
            csite%mmean_water_residual(ipa)  = csite%mmean_water_residual(ipa)  + csite%dmean_water_residual(ipa) 
         end do

      end do
      
   end do

   return
end subroutine integrate_ed_monthly_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine normalize_ed_monthly_output_vars(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine normalize the sum before writing the mobthly analysis. It also        !
! computes some of the variables that didn't need to be computed every day, like AGB.      !
!------------------------------------------------------------------------------------------!
   use ed_state_vars, only: edtype        & ! structure
                          , polygontype   & ! structure
                          , sitetype      & ! structure
                          , patchtype     ! ! structure
   use ed_misc_coms , only: current_time  & ! intent(in)
                          , simtime       ! ! intent(in)
   use ed_max_dims  , only: n_pft         & ! intent(in)
                          , n_dbh         & ! intent(in)
                          , n_dist_types  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   type(simtime)               :: lastmonth
   real                        :: ndaysi
   real                        :: polygon_area_i
   real                        :: site_area_i
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   integer                     :: ico
   integer                     :: ipft
   integer                     :: dbh
   integer                     :: ilu
   integer                     :: jlu
   real                        :: srnonm1
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   !     Finding the inverse of number of days used for this monthly integral.             !
   !---------------------------------------------------------------------------------------!
   call lastmonthdate(current_time,lastmonth,ndaysi)

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      !------------------------------------------------------------------------------------!
      !      First normalize the variables previously defined.                             !
      !------------------------------------------------------------------------------------!
      cgrid%mmean_gpp            (ipy) = cgrid%mmean_gpp            (ipy) * ndaysi
      cgrid%mmean_evap           (ipy) = cgrid%mmean_evap           (ipy) * ndaysi
      cgrid%mmean_transp         (ipy) = cgrid%mmean_transp         (ipy) * ndaysi
      cgrid%mmean_sensible_ac    (ipy) = cgrid%mmean_sensible_ac    (ipy) * ndaysi
      cgrid%mmean_sensible_gc    (ipy) = cgrid%mmean_sensible_gc    (ipy) * ndaysi
      cgrid%mmean_sensible_vc    (ipy) = cgrid%mmean_sensible_vc    (ipy) * ndaysi
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
      cgrid%mmean_gpp_lu       (:,ipy) = cgrid%mmean_gpp_lu       (:,ipy) * ndaysi
      cgrid%mmean_rh_lu        (:,ipy) = cgrid%mmean_rh_lu        (:,ipy) * ndaysi
      cgrid%mmean_nep_lu       (:,ipy) = cgrid%mmean_nep_lu       (:,ipy) * ndaysi
      cgrid%mmean_gpp_dbh      (:,ipy) = cgrid%mmean_gpp_dbh      (:,ipy) * ndaysi
      cgrid%mmean_veg_energy     (ipy) = cgrid%mmean_veg_energy     (ipy) * ndaysi
      cgrid%mmean_veg_hcap       (ipy) = cgrid%mmean_veg_hcap       (ipy) * ndaysi
      cgrid%mmean_veg_water      (ipy) = cgrid%mmean_veg_water      (ipy) * ndaysi
      cgrid%mmean_veg_temp       (ipy) = cgrid%mmean_veg_temp       (ipy) * ndaysi
      cgrid%mmean_can_enthalpy   (ipy) = cgrid%mmean_can_enthalpy   (ipy) * ndaysi
      cgrid%mmean_can_temp       (ipy) = cgrid%mmean_can_temp       (ipy) * ndaysi
      cgrid%mmean_can_shv        (ipy) = cgrid%mmean_can_shv        (ipy) * ndaysi
      cgrid%mmean_can_co2        (ipy) = cgrid%mmean_can_co2        (ipy) * ndaysi
      cgrid%mmean_can_rhos       (ipy) = cgrid%mmean_can_rhos       (ipy) * ndaysi
      cgrid%mmean_atm_temp       (ipy) = cgrid%mmean_atm_temp       (ipy) * ndaysi
      cgrid%mmean_atm_shv        (ipy) = cgrid%mmean_atm_shv        (ipy) * ndaysi
      cgrid%mmean_atm_prss       (ipy) = cgrid%mmean_atm_prss       (ipy) * ndaysi
      cgrid%mmean_atm_vels       (ipy) = cgrid%mmean_atm_vels       (ipy) * ndaysi
      cgrid%mmean_pcpg           (ipy) = cgrid%mmean_pcpg           (ipy) * ndaysi
      cgrid%mmean_lai_pft      (:,ipy) = cgrid%mmean_lai_pft      (:,ipy) * ndaysi
      cgrid%mmean_lai_lu       (:,ipy) = cgrid%mmean_lai_lu       (:,ipy) * ndaysi
      cgrid%mmean_wpa_pft      (:,ipy) = cgrid%mmean_wpa_pft      (:,ipy) * ndaysi
      cgrid%mmean_wpa_lu       (:,ipy) = cgrid%mmean_wpa_lu       (:,ipy) * ndaysi
      cgrid%mmean_wai_pft      (:,ipy) = cgrid%mmean_wai_pft      (:,ipy) * ndaysi
      cgrid%mmean_wai_lu       (:,ipy) = cgrid%mmean_wai_lu       (:,ipy) * ndaysi

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
      cgrid%stdev_gpp     (ipy) = srnonm1 * sqrt(cgrid%stdev_gpp     (ipy) * ndaysi - cgrid%mmean_gpp     (ipy) ** 2)
      cgrid%stdev_evap    (ipy) = srnonm1 * sqrt(cgrid%stdev_evap    (ipy) * ndaysi - cgrid%mmean_evap    (ipy) ** 2)
      cgrid%stdev_transp  (ipy) = srnonm1 * sqrt(cgrid%stdev_transp  (ipy) * ndaysi - cgrid%mmean_transp  (ipy) ** 2)
      cgrid%stdev_sensible(ipy) = srnonm1 * sqrt(cgrid%stdev_sensible(ipy) * ndaysi - cgrid%mmean_sensible_ac(ipy) ** 2)
      cgrid%stdev_nep     (ipy) = srnonm1 * sqrt(cgrid%stdev_nep     (ipy) * ndaysi - cgrid%mmean_nep     (ipy) ** 2)
      cgrid%stdev_rh      (ipy) = srnonm1 * sqrt(cgrid%stdev_rh      (ipy) * ndaysi - cgrid%mmean_rh      (ipy) ** 2)
  
      !---- Finding AGB and basal area per PFT --------------------------------------------!
      polygon_area_i = 1./sum(cpoly%area)
      do ipft = 1,n_pft
        do dbh =1,n_dbh
          cgrid%agb_pft(ipft,ipy) = cgrid%agb_pft(ipft,ipy)                                &
                                  + sum(cpoly%agb(ipft,dbh,:)*cpoly%area)*polygon_area_i
          cgrid%ba_pft(ipft,ipy)  = cgrid%ba_pft(ipft,ipy)                                 &
                                  + sum(cpoly%basal_area(ipft,dbh,:)*cpoly%area)           &
                                  * polygon_area_i
        end do
      end do

      !----- Finding AGB per land use type ------------------------------------------------!
      do ilu = 1,n_dist_types
          cgrid%agb_lu(ilu,ipy) = cgrid%agb_lu(ilu,ipy)                                    &
                                + sum(cpoly%agb_lu(ilu,:)*cpoly%area)*polygon_area_i
      end do

      !----- Finding disturbance rates per source and target land use types. --------------!
      do ilu = 1,n_dist_types
         do jlu = 1,n_dist_types
          cgrid%disturbance_rates(ilu,jlu,ipy) = cgrid%disturbance_rates(ilu,jlu,ipy)      &
                                               + sum( cpoly%disturbance_rates(ilu,jlu,:)   &
                                                    * cpoly%area)                          &
                                               * polygon_area_i
         end do
      end do

      !----- Finding the fractional area covered by each land use and each PFT ------------!
      cgrid%area_pft(:,ipy) = 0.
      cgrid%area_lu(:,ipy)  = 0.
      siteloop: do isi = 1, cpoly%nsites
         csite => cpoly%site(isi)

         cpoly%mmean_co2_residual(isi)    = cpoly%mmean_co2_residual(isi)    * ndaysi
         cpoly%mmean_energy_residual(isi) = cpoly%mmean_energy_residual(isi) * ndaysi
         cpoly%mmean_water_residual(isi)  = cpoly%mmean_water_residual(isi)  * ndaysi


         site_area_i = 1./sum(csite%area)
         patchloop: do ipa=1,csite%npatches

            csite%mmean_co2_residual(ipa)    = csite%mmean_co2_residual(ipa)    * ndaysi
            csite%mmean_energy_residual(ipa) = csite%mmean_energy_residual(ipa) * ndaysi
            csite%mmean_water_residual(ipa)  = csite%mmean_water_residual(ipa)  * ndaysi

            ilu = csite%dist_type(ipa)
            cgrid%area_lu(ilu,ipy) = cgrid%area_lu(ilu,ipy)                                &
                                   + csite%area(ipa) * cpoly%area(isi)                     &
                                   * site_area_i * polygon_area_i
            cpatch => csite%patch(ipa)
            if (cpatch%ncohorts > 0) then
               pftloop: do ipft=1,n_pft
                  if (any(cpatch%pft(:) == ipft)) then
                     cgrid%area_pft(ipft,ipy) = cgrid%area_pft(ipft,ipy)                   &
                                              + csite%area(ipa) * cpoly%area(isi)          &
                                              * site_area_i * polygon_area_i
                  end if
               end do pftloop
            end if
         end do patchloop
      end do siteloop

   end do polyloop

   return
end subroutine normalize_ed_monthly_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zero_ed_monthly_output_vars(cgrid)
   use ed_state_vars,only:edtype,polygontype,sitetype

   implicit none
   type(edtype)     , target  :: cgrid
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   integer                    :: ipy,isi,ipa

   !----- The loop is necessary for coupled runs (when npolygons may be 0) ----------------!
   do ipy=1,cgrid%npolygons
      cgrid%mmean_gpp            (ipy) = 0.
      cgrid%mmean_evap           (ipy) = 0.
      cgrid%mmean_transp         (ipy) = 0.
      cgrid%mmean_sensible_ac    (ipy) = 0.
      cgrid%mmean_sensible_gc    (ipy) = 0.
      cgrid%mmean_sensible_vc    (ipy) = 0.
      cgrid%mmean_nep            (ipy) = 0.
      cgrid%mmean_plresp         (ipy) = 0.
      cgrid%mmean_rh             (ipy) = 0.
      cgrid%mmean_leaf_resp      (ipy) = 0.
      cgrid%mmean_root_resp      (ipy) = 0.
      cgrid%mmean_growth_resp    (ipy) = 0.
      cgrid%mmean_storage_resp   (ipy) = 0.
      cgrid%mmean_vleaf_resp     (ipy) = 0.
      cgrid%mmean_soil_temp    (:,ipy) = 0.
      cgrid%mmean_soil_water   (:,ipy) = 0.
      cgrid%mmean_gpp_lu       (:,ipy) = 0.
      cgrid%mmean_rh_lu        (:,ipy) = 0.
      cgrid%mmean_nep_lu       (:,ipy) = 0.
      cgrid%mmean_gpp_dbh      (:,ipy) = 0.
      cgrid%mmean_veg_energy     (ipy) = 0.
      cgrid%mmean_veg_hcap       (ipy) = 0.
      cgrid%mmean_veg_water      (ipy) = 0.
      cgrid%mmean_veg_temp       (ipy) = 0.
      cgrid%mmean_can_enthalpy   (ipy) = 0.
      cgrid%mmean_can_temp       (ipy) = 0.
      cgrid%mmean_can_shv        (ipy) = 0.
      cgrid%mmean_can_co2        (ipy) = 0.
      cgrid%mmean_can_rhos       (ipy) = 0.
      cgrid%mmean_atm_temp       (ipy) = 0.
      cgrid%mmean_atm_shv        (ipy) = 0.
      cgrid%mmean_atm_prss       (ipy) = 0.
      cgrid%mmean_atm_vels       (ipy) = 0.
      cgrid%mmean_pcpg           (ipy) = 0.
      cgrid%mmean_lai_pft      (:,ipy) = 0.
      cgrid%mmean_lai_lu       (:,ipy) = 0.
      cgrid%mmean_wpa_pft      (:,ipy) = 0.
      cgrid%mmean_wpa_lu       (:,ipy) = 0.
      cgrid%mmean_wai_pft      (:,ipy) = 0.
      cgrid%mmean_wai_lu       (:,ipy) = 0.
      cgrid%agb_pft            (:,ipy) = 0.
      cgrid%ba_pft             (:,ipy) = 0.
      cgrid%agb_lu             (:,ipy) = 0.
      cgrid%stdev_gpp            (ipy) = 0.
      cgrid%stdev_evap           (ipy) = 0.
      cgrid%stdev_transp         (ipy) = 0.
      cgrid%stdev_sensible       (ipy) = 0.
      cgrid%stdev_nep            (ipy) = 0.
      cgrid%stdev_rh             (ipy) = 0.
      cgrid%disturbance_rates(:,:,ipy) = 0.

      cgrid%mmean_co2_residual(ipy)    = 0.
      cgrid%mmean_energy_residual(ipy) = 0.
      cgrid%mmean_water_residual(ipy)  = 0.
      
      cpoly => cgrid%polygon(ipy)
      do isi = 1, cpoly%nsites

         cpoly%mmean_co2_residual(isi)    = 0.
         cpoly%mmean_energy_residual(isi) = 0.
         cpoly%mmean_water_residual(isi)  = 0.

         csite => cpoly%site(isi)
         do ipa=1,csite%npatches
            csite%mmean_co2_residual(ipa)    = 0.
            csite%mmean_energy_residual(ipa) = 0.
            csite%mmean_water_residual(ipa)  = 0.
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
                       ed_biomass(cpatch%bdead(ico), cpatch%balive(ico),   &
                       cpatch%bleaf(ico), cpatch%pft(ico), cpatch%hite(ico),   &
                       cpatch%bstorage(ico)) * csite%area(ipa) * &
                       cpoly%area(isi) * 10.0 * cpatch%nplant(ico)
                  cgrid%total_basal_area_recruit(ipy) =   &
                       cgrid%total_basal_area_recruit(ipy) +   &
                       pi1 * 0.25 * cpatch%dbh(ico)**2 * &
                       cpatch%nplant(ico) * csite%area(ipa) * cpoly%area(isi)
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
