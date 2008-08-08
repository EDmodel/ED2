!==========================================================================================!
!==========================================================================================!
!                           |----------------------------------|                           !
!                           |** FREQUENT AVERAGE SUBROUTINES **|                           !
!                           |----------------------------------|                           !
!==========================================================================================!
!==========================================================================================!
subroutine normalize_averaged_vars_ar(cgrid,frqsum,dtlsm)

   use grid_coms, only: nzg
   use misc_coms, only: radfrq
   use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
   use ed_misc_coms,only:diag_veg_heating
   

   implicit none
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy,isi,ipa,ico
   real                        :: frqsum,dtlsm,tfact,frqsumi
   real                        :: gpp_sum
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

         csite%aux(:)              = csite%aux(:)               * frqsumi
         csite%avg_vapor_vc(:)     = csite%avg_vapor_vc(:)      * frqsumi
         csite%avg_dew_cg(:)       = csite%avg_dew_cg(:)        * frqsumi
         csite%avg_vapor_gc(:)     = csite%avg_vapor_gc(:)      * frqsumi
         csite%avg_wshed_vg(:)     = csite%avg_wshed_vg(:)      * frqsumi
         csite%avg_vapor_ac(:)     = csite%avg_vapor_ac(:)      * frqsumi
         csite%avg_transp(:)       = csite%avg_transp(:)        * frqsumi
         csite%avg_evap(:)         = csite%avg_evap(:)          * frqsumi
         csite%avg_runoff(:)       = csite%avg_runoff(:)        * tfact ! kg (water) / m2 / s
         csite%avg_sensible_vc(:)  = csite%avg_sensible_vc(:)   * frqsumi
         csite%avg_sensible_2cas(:)= csite%avg_sensible_2cas(:) * frqsumi
         csite%avg_qwshed_vg(:)    = csite%avg_qwshed_vg(:)     * frqsumi
         csite%avg_sensible_gc(:)  = csite%avg_sensible_gc(:)   * frqsumi
         csite%avg_sensible_ac(:)  = csite%avg_sensible_ac(:)   * frqsumi
         csite%avg_sensible_tot(:) = csite%avg_sensible_tot(:)  * frqsumi
         csite%avg_carbon_ac(:)    = csite%avg_carbon_ac(:)     * frqsumi
         csite%avg_runoff_heat(:)  = csite%avg_runoff_heat(:)   * tfact
         ! csite%avg_heatstor_veg(:) = csite%avg_heatstor_veg(:)  * tfact ! CHECK THIS?
         
         do k=cpoly%lsl(isi),nzg
            csite%avg_sensible_gg(k,:) = csite%avg_sensible_gg(k,:) * frqsumi
            csite%avg_smoist_gg(k,:)   = csite%avg_smoist_gg(k,:)   * frqsumi
            csite%avg_smoist_gc(k,:)   = csite%avg_smoist_gc(k,:)   * frqsumi
            csite%aux_s(k,:)           = csite%aux_s(k,:)           * frqsumi
         end do

         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            
            if (cpatch%ncohorts>0) then
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
               cpatch%mean_leaf_resp = cpatch%mean_leaf_resp * tfact
               cpatch%mean_root_resp = cpatch%mean_root_resp * tfact
               cpatch%mean_gpp       = cpatch%mean_gpp       * tfact
               
               !  Normalize the vegetation heating/cooling rates
               !  
               if(diag_veg_heating) then
                  cpatch%co_srad_h = cpatch%co_srad_h * frqsumi
                  cpatch%co_lrad_h = cpatch%co_lrad_h * frqsumi
                  cpatch%co_sens_h = cpatch%co_sens_h * frqsumi
                  cpatch%co_evap_h = cpatch%co_evap_h * frqsumi
                  cpatch%co_liqr_h = cpatch%co_liqr_h * frqsumi
               endif
            end if
         end do
      end do
   end do
   
   return
end subroutine normalize_averaged_vars_ar
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
      cpoly => cgrid%polygon(ipy)
      do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !----------------------------------------------------------------!
         ! Zeroing CO2 budget variables                                   !
         !----------------------------------------------------------------!
         csite%co2budget_gpp      = 0.0
         csite%co2budget_gpp_dbh  = 0.0
         csite%co2budget_rh       = 0.0
         csite%co2budget_plresp   = 0.0
         csite%wbudget_precipgain = 0.0
         csite%ebudget_precipgain = 0.0
         csite%ebudget_netrad     = 0.0
         !----------------------------------------------------------------!

         csite%avg_carbon_ac    = 0.0
         csite%avg_vapor_vc     = 0.0
         csite%avg_dew_cg       = 0.0
         csite%avg_vapor_gc     = 0.0
         csite%avg_wshed_vg     = 0.0
         csite%avg_vapor_ac     = 0.0
         csite%avg_transp       = 0.0
         csite%avg_evap         = 0.0
         csite%avg_smoist_gg    = 0.0
         csite%avg_smoist_gc    = 0.0
         csite%avg_runoff       = 0.0
         csite%avg_sensible_vc  = 0.0
         csite%avg_sensible_2cas= 0.0
         csite%avg_qwshed_vg    = 0.0
         csite%avg_sensible_gc  = 0.0
         csite%avg_sensible_ac  = 0.0
         csite%avg_sensible_tot = 0.0
         csite%avg_sensible_gg  = 0.0
         csite%avg_runoff_heat  = 0.0
         csite%aux              = 0.0
         csite%aux_s            = 0.0
         
         csite%avg_heatstor_veg = 0.0  !SHOULD THIS BE ZERO'D ALSO?

         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            if (cpatch%ncohorts>0) then
               cpatch%leaf_respiration    = 0.0
               cpatch%root_respiration    = 0.0
               !Still checking if these should be zerod, should not be necessary
               !They are not integrated variables...
               !cpatch%mean_leaf_resp      = 0.0
               !cpatch%mean_root_resp      = 0.0
               cpatch%gpp                 = 0.0

               cpatch%co_srad_h = 0.0
               cpatch%co_lrad_h = 0.0
               cpatch%co_sens_h = 0.0
               cpatch%co_evap_h = 0.0
               cpatch%co_liqr_h = 0.0

            endif
         end do


      enddo

      cpoly%avg_soil_temp(:,:)      = 0.0
      cpoly%avg_soil_water(:,:)     = 0.0
   enddo

   ! Should this be here as well?
   cgrid%cbudget_nep       = 0.0

   ! Reset the meteorological diagnostic
   
   cgrid%avg_nir_beam(:)       = 0.0
   cgrid%avg_nir_diffuse(:)    = 0.0
   cgrid%avg_par_beam(:)       = 0.0
   cgrid%avg_par_diffuse(:)    = 0.0
   cgrid%avg_atm_tmp(:)        = 0.0
   cgrid%avg_atm_shv(:)        = 0.0
   cgrid%avg_rhos(:)           = 0.0
   cgrid%avg_theta(:)          = 0.0
   cgrid%avg_rshort(:)         = 0.0
   cgrid%avg_rshort_diffuse(:) = 0.0
   cgrid%avg_rlong(:)          = 0.0
   cgrid%avg_pcpg(:)           = 0.0
   cgrid%avg_qpcpg(:)          = 0.0
   cgrid%avg_dpcpg(:)          = 0.0
   cgrid%avg_vels(:)           = 0.0
   cgrid%avg_prss(:)           = 0.0
   cgrid%avg_exner(:)          = 0.0
   cgrid%avg_geoht(:)          = 0.0
   cgrid%avg_atm_co2(:)        = 0.0
   cgrid%avg_albedt(:)         = 0.0
   cgrid%avg_rlongup(:)        = 0.0
   
   !
   cgrid%avg_evap(:)           = 0.0
   cgrid%avg_transp(:)         = 0.0
   cgrid%avg_sensible_tot(:)   = 0.0
   cgrid%avg_soil_temp(:,:)    = 0.0
   cgrid%avg_soil_water(:,:)   = 0.0

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
subroutine integrate_ed_daily_output_state(cgrid)
!------------------------------------------------------------------------------------------!
!    This subroutine integrates the daily averages for state vars. This is called after    !
! each time step in case daily or monthly averages were requested.                         !
!------------------------------------------------------------------------------------------!
   use ed_state_vars        , only: edtype,polygontype,sitetype,patchtype
   use grid_coms            , only : nzg
   use max_dims             , only : n_dbh, n_pft, n_dist_types
   use canopy_radiation_coms, only: lai_min
   implicit none
   
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy, isi, ipa, ico, k,lu,dbh
   real                        :: poly_area_i,site_area_i, patch_lai_i
   real                        :: forest_site,forest_site_i, forest_poly
   ! These are auxiliary variables for averaging sitetype variables
   !----------------------------------------------------------------------------------------
   real :: sitesum_fsn, sitesum_fsw
   ! These are auxiliary variables for averaging patchtype variables
   !----------------------------------------------------------------------------------------
   real :: patchsum_fsn, patchsum_fsw
   real, dimension(n_dbh)        :: patchsum_gpp_dbh

         
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      ! Initialize auxiliary variables to add sitetype variables
      sitesum_fsn          = 0.
      sitesum_fsw          = 0.
      
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
         patchsum_fsn          = 0.
         patchsum_fsw          = 0.

         ! Looping through the patches to normalize the sum of all cohorts.
         patchloop: do ipa=1, csite%npatches
            cpatch => csite%patch(ipa)
            
            patch_lai_i = 1./max(tiny(1.),sum(cpatch%lai,cpatch%lai > lai_min))
            
            patchsum_fsn = patchsum_fsn + (sum(cpatch%fsn * cpatch%lai,cpatch%lai > lai_min) * patch_lai_i) * csite%area(ipa)
            patchsum_fsw = patchsum_fsw + (sum(cpatch%fsw * cpatch%lai,cpatch%lai > lai_min) * patch_lai_i) * csite%area(ipa)

         end do patchloop
         
         ! Variables already average at the sitetype level, just add them to polygontype level
         sitesum_fsn          = sitesum_fsn          + (patchsum_fsn            *site_area_i) * cpoly%area(isi)
         sitesum_fsw          = sitesum_fsw          + (patchsum_fsw            *site_area_i) * cpoly%area(isi)
         
      end do siteloop
      
      cgrid%dmean_fsn(ipy)          = cgrid%dmean_fsn(ipy)          + sitesum_fsn          * poly_area_i
      cgrid%dmean_fsw(ipy)          = cgrid%dmean_fsw(ipy)          + sitesum_fsw          * poly_area_i
      

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
   use max_dims             , only : n_dbh, n_pft, n_dist_types
   use canopy_radiation_coms, only : lai_min
   use misc_coms            , only : dtlsm
   implicit none
   
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                     :: ipy, isi, ipa, ico, k,lu,dbh
   real                        :: poly_area_i,site_area_i
   real                        :: forest_site,forest_site_i, forest_poly
   real, parameter             :: ddbh=1./(n_dbh-1)
   ! These are auxiliary variables for averaging sitetype variables
   !----------------------------------------------------------------------------------------
   real :: sitesum_leaf_resp, sitesum_root_resp 
   real :: sitesum_plresp,sitesum_gpp,sitesum_rh
   real :: sitesum_evap,sitesum_transp,sitesum_sensible_tot
   real, dimension(n_dist_types) :: sitesum_gpp_lu, sitesum_rh_lu,sitesum_nep_lu
   real, dimension(n_dbh)        :: sitesum_gpp_dbh
   ! These are auxiliary variables for averaging patchtype variables
   !----------------------------------------------------------------------------------------
   real :: patchsum_leaf_resp   , patchsum_root_resp
   
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      poly_area_i=1./sum(cpoly%area)

      
      ! Initialize auxiliary variables to add sitetype variables
      sitesum_rh           = 0.
      sitesum_leaf_resp    = 0.
      sitesum_root_resp    = 0.
      sitesum_plresp       = 0.
      sitesum_gpp          = 0.
      sitesum_rh           = 0.
      sitesum_gpp_lu       = 0.
      sitesum_nep_lu       = 0.
      sitesum_rh_lu        = 0.
      sitesum_gpp_dbh      = 0.
      sitesum_evap         = 0.
      sitesum_transp       = 0.
      sitesum_sensible_tot = 0.
      
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

            patchsum_leaf_resp = patchsum_leaf_resp + sum(cpatch%mean_leaf_resp, cpatch%lai > lai_min)  * csite%area(ipa) 
            patchsum_root_resp = patchsum_root_resp + sum(cpatch%mean_root_resp                      )  * csite%area(ipa) 

         end do patchloop
         
         ! Variables already average at the sitetype level, just add them to polygontype level
         sitesum_leaf_resp    = sitesum_leaf_resp    + (patchsum_leaf_resp      *site_area_i) * cpoly%area(isi)
         sitesum_root_resp    = sitesum_root_resp    + (patchsum_root_resp      *site_area_i) * cpoly%area(isi)

         sitesum_rh           = sitesum_rh           + (sum(csite%co2budget_rh    *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_gpp          = sitesum_gpp          + (sum(csite%co2budget_gpp   *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_plresp       = sitesum_plresp       + (sum(csite%co2budget_plresp*csite%area) * site_area_i) * cpoly%area(isi)
         
         sitesum_evap         = sitesum_evap         + (sum(csite%avg_evap        *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_transp       = sitesum_transp       + (sum(csite%avg_transp      *csite%area) * site_area_i) * cpoly%area(isi)
         sitesum_sensible_tot = sitesum_sensible_tot + (sum(csite%avg_sensible_tot*csite%area) * site_area_i) * cpoly%area(isi)
         
         
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
      do lu=1,n_dist_types
         cgrid%dmean_gpp_lu(lu,ipy) = cgrid%dmean_gpp_lu(lu,ipy)    + sitesum_gpp_lu(lu)   * poly_area_i
         cgrid%dmean_nep_lu(lu,ipy) = cgrid%dmean_nep_lu(lu,ipy)    + sitesum_nep_lu(lu)   * poly_area_i
         cgrid%dmean_rh_lu(lu,ipy)  = cgrid%dmean_rh_lu(lu,ipy)     + sitesum_rh_lu(lu)    * poly_area_i
      end do
      
      do dbh=1,n_dbh
         cgrid%dmean_gpp_dbh(dbh,ipy) = cgrid%dmean_gpp_dbh(dbh,ipy) + sitesum_gpp_dbh(dbh)* poly_area_i
      end do

   end do polyloop

   !These variables are already averaged at gridtype, just add them up
   cgrid%dmean_soil_temp(:,:)  = cgrid%dmean_soil_temp(:,:)  + cgrid%avg_soil_temp(:,:) 
   cgrid%dmean_soil_water(:,:) = cgrid%dmean_soil_water(:,:) + cgrid%avg_soil_water(:,:)
   
   cgrid%dmean_evap(:)         = cgrid%dmean_evap(:)         + cgrid%avg_evap(:)     
   cgrid%dmean_transp(:)       = cgrid%dmean_transp(:)       + cgrid%avg_transp(:)   
   cgrid%dmean_sensible_vc(:)  = cgrid%dmean_sensible_vc(:)  + cgrid%avg_sensible_vc(:) 
   cgrid%dmean_sensible_gc(:)  = cgrid%dmean_sensible_gc(:)  + cgrid%avg_sensible_gc(:) 
   cgrid%dmean_sensible_ac(:)  = cgrid%dmean_sensible_ac(:)  + cgrid%avg_sensible_ac(:) 
   cgrid%dmean_sensible(:)     = cgrid%dmean_sensible(:)     + cgrid%avg_sensible_tot(:) 

   return
end subroutine integrate_ed_daily_output_flux
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine normalize_ed_daily_vars(cgrid,timefac1)
   use ed_state_vars , only : edtype,polygontype,sitetype,patchtype
   use max_dims      , only : n_pft
   implicit none
   real, intent(in)         :: timefac1 ! Daily sum          => daily average
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   
   integer                  :: ipy,isi,ipa,ipft
   real                     :: polygon_area_i, site_area_i
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
            
            cpatch%dmean_gpp       = cpatch%dmean_gpp       * timefac1
            cpatch%dmean_gpp_pot   = cpatch%dmean_gpp_pot   * timefac1
            cpatch%dmean_gpp_max   = cpatch%dmean_gpp_max   * timefac1
            cpatch%dmean_leaf_resp = cpatch%dmean_leaf_resp * timefac1
            cpatch%dmean_root_resp = cpatch%dmean_root_resp * timefac1
         
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
   use ed_state_vars , only : edtype,polygontype,sitetype,patchtype
   use max_dims      , only : n_pft
   use consts_coms   , only : alvl,day_sec,umol_2_kgC
   use misc_coms     , only : dtlsm,frqsum
   implicit none
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   type(patchtype)   , pointer :: cpatch
   integer                  :: ipy,isi,ipa,ipft
   real                     :: polygon_area_i, site_area_i
   real :: sitesum_storage_resp, sitesum_vleaf_resp, sitesum_growth_resp
   real :: patchsum_storage_resp, patchsum_vleaf_resp , patchsum_growth_resp
   
   logical           , save :: find_factors=.true.
   real              , save :: dtlsm_o_daysec=1.E34, frqsum_o_daysec=1.E34
   
   ! Computing the normalization factors. This is done once.
   if (find_factors) then
      dtlsm_o_daysec  = dtlsm/day_sec
      frqsum_o_daysec = frqsum/day_sec
      find_factors    = .false.
   end if
   
   ! State variables, updated every time step, so these are normalized by dtlsm/day_sec
   cgrid%dmean_fsw               = cgrid%dmean_fsw               * dtlsm_o_daysec
   cgrid%dmean_fsn               = cgrid%dmean_fsn               * dtlsm_o_daysec
   
   ! State variables, updated every frqsum, so these are normalized by frqsum/day_sec
   cgrid%dmean_soil_temp         = cgrid%dmean_soil_temp         * frqsum_o_daysec
   cgrid%dmean_soil_water        = cgrid%dmean_soil_water        * frqsum_o_daysec

   ! Flux variables, updated every frqsum, so these are normalized by frqsum/day_sec
   cgrid%dmean_evap              = cgrid%dmean_evap              * frqsum_o_daysec
   cgrid%dmean_transp            = cgrid%dmean_transp            * frqsum_o_daysec
   cgrid%dmean_sensible_vc       = cgrid%dmean_sensible_vc       * frqsum_o_daysec
   cgrid%dmean_sensible_gc       = cgrid%dmean_sensible_gc       * frqsum_o_daysec
   cgrid%dmean_sensible_ac       = cgrid%dmean_sensible_ac       * frqsum_o_daysec
   cgrid%dmean_sensible          = cgrid%dmean_sensible          * frqsum_o_daysec

   ! Carbon flux variables should be total flux integrated over the day at this point, [µmol/m²/s] -> [kgC/m²/day]
   cgrid%dmean_leaf_resp         = cgrid%dmean_leaf_resp         * umol_2_kgC * day_sec
   cgrid%dmean_root_resp         = cgrid%dmean_root_resp         * umol_2_kgC * day_sec

   ! Carbon flux variables should be total flux integrated over the day at this point, [µmol/m²/day] -> [kgC/m²/day]
   cgrid%dmean_rh                = cgrid%dmean_rh                * umol_2_kgC
   cgrid%dmean_gpp               = cgrid%dmean_gpp               * umol_2_kgC
   cgrid%dmean_plresp            = cgrid%dmean_plresp            * umol_2_kgC
   cgrid%dmean_nep               = cgrid%dmean_nep               * umol_2_kgC
   
   cgrid%dmean_gpp_lu            = cgrid%dmean_gpp_lu            * umol_2_kgC
   cgrid%dmean_rh_lu             = cgrid%dmean_rh_lu             * umol_2_kgC
   cgrid%dmean_nep_lu            = cgrid%dmean_nep_lu            * umol_2_kgC
   cgrid%dmean_gpp_dbh           = cgrid%dmean_gpp_dbh           * umol_2_kgC

   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      polygon_area_i = 1./sum(cpoly%area)

      sitesum_growth_resp  = 0.
      sitesum_storage_resp = 0.
      sitesum_vleaf_resp   = 0.
      
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         
         site_area_i = 1./sum(csite%area)
         
         patchsum_growth_resp  = 0.
         patchsum_storage_resp = 0.
         patchsum_vleaf_resp   = 0.

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            
            ! These are in kg/plant/day, converting to kg/m²/day
            patchsum_growth_resp  = patchsum_growth_resp  + sum(cpatch%growth_respiration  * cpatch%nplant) * csite%area(ipa)
            patchsum_storage_resp = patchsum_storage_resp + sum(cpatch%storage_respiration * cpatch%nplant) * csite%area(ipa)
            patchsum_vleaf_resp   = patchsum_vleaf_resp   + sum(cpatch%vleaf_respiration   * cpatch%nplant) * csite%area(ipa)

            do ipft=1,n_pft
               cpoly%lai_pft(ipft,isi) = cpoly%lai_pft(ipft,isi) + &
                    sum(cpatch%lai,cpatch%pft == ipft) * csite%area(ipa) * site_area_i
            end do
         end do patchloop

         sitesum_growth_resp  = sitesum_growth_resp  + (patchsum_growth_resp    *site_area_i) * cpoly%area(isi)
         sitesum_storage_resp = sitesum_storage_resp + (patchsum_storage_resp   *site_area_i) * cpoly%area(isi)
         sitesum_vleaf_resp   = sitesum_vleaf_resp   + (patchsum_vleaf_resp     *site_area_i) * cpoly%area(isi)

      end do siteloop
      
      do ipft=1,n_pft
         cgrid%lai_pft(ipft,ipy) = cgrid%lai_pft(ipft,ipy) + sum(cpoly%lai_pft(ipft,:)*cpoly%area) * polygon_area_i
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

         !--------------------------------------!
         ! Reset variables stored in sitetype   !
         !--------------------------------------!
         csite%dmean_A_decomp  = 0.0
         csite%dmean_Af_decomp = 0.0

         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            
            !-------------------------------------!
            ! Reset variables stored in patchtype !
            !-------------------------------------!
            cpatch%dmean_gpp       = 0.0
            cpatch%dmean_gpp_pot   = 0.0
            cpatch%dmean_gpp_max   = 0.0
            cpatch%dmean_leaf_resp = 0.0
            cpatch%dmean_root_resp = 0.0
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
   type(patchtype)  , pointer :: cpatch
   integer                    :: ipy,isi,ipa,ico
   
   !--------------------------------!
   ! Variables stored in edtype     !
   !--------------------------------!
   cgrid%dmean_gpp               = 0.
   cgrid%dmean_evap              = 0.
   cgrid%dmean_transp            = 0.
   cgrid%dmean_sensible_vc       = 0.
   cgrid%dmean_sensible_gc       = 0.
   cgrid%dmean_sensible_ac       = 0.
   cgrid%dmean_sensible          = 0.
   cgrid%dmean_plresp            = 0.
   cgrid%dmean_rh                = 0.
   cgrid%dmean_leaf_resp         = 0.
   cgrid%dmean_root_resp         = 0.
   cgrid%dmean_growth_resp       = 0.
   cgrid%dmean_storage_resp      = 0.
   cgrid%dmean_vleaf_resp        = 0.
   cgrid%dmean_nep               = 0.
   cgrid%dmean_soil_temp         = 0.
   cgrid%dmean_soil_water        = 0.
   cgrid%dmean_fsw               = 0.
   cgrid%dmean_fsn               = 0.
   cgrid%dmean_gpp_lu            = 0.
   cgrid%dmean_rh_lu             = 0.
   cgrid%dmean_nep_lu            = 0.
   cgrid%dmean_gpp_dbh           = 0.
   cgrid%lai_pft                 = 0.

   do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      
      !-----------------------------------------!
      ! Reset variables stored in polygontype   !
      !-----------------------------------------!
      cpoly%lai_pft = 0.
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
   use ed_state_vars, only : edtype
   use      max_dims, only : n_dbh,n_pft, n_dist_types
   implicit none
   
   type(edtype)      , target  :: cgrid
   
   !------------------------------------------------------------------------------!
   ! First the mean variables that can be computed from the daily averages        !
   !------------------------------------------------------------------------------!
   cgrid%mmean_gpp          = cgrid%mmean_gpp           +  cgrid%dmean_gpp
   cgrid%mmean_evap         = cgrid%mmean_evap          +  cgrid%dmean_evap
   cgrid%mmean_transp       = cgrid%mmean_transp        +  cgrid%dmean_transp
   cgrid%mmean_sensible     = cgrid%mmean_sensible      +  cgrid%dmean_sensible
   cgrid%mmean_sensible_ac  = cgrid%mmean_sensible_ac   +  cgrid%dmean_sensible_ac
   cgrid%mmean_sensible_gc  = cgrid%mmean_sensible_gc   +  cgrid%dmean_sensible_gc
   cgrid%mmean_sensible_vc  = cgrid%mmean_sensible_vc   +  cgrid%dmean_sensible_vc
   cgrid%mmean_nep          = cgrid%mmean_nep           +  cgrid%dmean_nep
   cgrid%mmean_soil_temp    = cgrid%mmean_soil_temp     +  cgrid%dmean_soil_temp
   cgrid%mmean_soil_water   = cgrid%mmean_soil_water    +  cgrid%dmean_soil_water
   cgrid%mmean_plresp       = cgrid%mmean_plresp        +  cgrid%dmean_plresp
   cgrid%mmean_rh           = cgrid%mmean_rh            +  cgrid%dmean_rh
   cgrid%mmean_leaf_resp    = cgrid%mmean_leaf_resp     +  cgrid%dmean_leaf_resp
   cgrid%mmean_root_resp    = cgrid%mmean_root_resp     +  cgrid%dmean_root_resp
   cgrid%mmean_growth_resp  = cgrid%mmean_growth_resp   +  cgrid%dmean_growth_resp
   cgrid%mmean_storage_resp = cgrid%mmean_storage_resp  +  cgrid%dmean_storage_resp
   cgrid%mmean_vleaf_resp   = cgrid%mmean_vleaf_resp    +  cgrid%dmean_vleaf_resp
   cgrid%mmean_gpp_lu       = cgrid%mmean_gpp_lu        +  cgrid%dmean_gpp_lu
   cgrid%mmean_rh_lu        = cgrid%mmean_rh_lu         +  cgrid%dmean_rh_lu
   cgrid%mmean_nep_lu       = cgrid%mmean_nep_lu        +  cgrid%dmean_nep_lu
   cgrid%mmean_gpp_dbh      = cgrid%mmean_gpp_dbh       +  cgrid%dmean_gpp_dbh
   cgrid%mmean_lai_pft      = cgrid%mmean_lai_pft       +  cgrid%lai_pft

   !------------------------------------------------------------------------------------!
   !    During the integration stage we keep the sum of squares, it will be converted   !
   ! to standard deviation right before the monthly output.                             !
   !------------------------------------------------------------------------------------!
   cgrid%stdev_gpp          = cgrid%stdev_gpp           +  cgrid%dmean_gpp           ** 2
   cgrid%stdev_evap         = cgrid%stdev_evap          +  cgrid%dmean_evap          ** 2
   cgrid%stdev_transp       = cgrid%stdev_transp        +  cgrid%dmean_transp        ** 2
   cgrid%stdev_sensible     = cgrid%stdev_sensible      +  cgrid%dmean_sensible      ** 2
   cgrid%stdev_nep          = cgrid%stdev_nep           +  cgrid%dmean_nep           ** 2
   cgrid%stdev_rh           = cgrid%stdev_rh            +  cgrid%dmean_rh            ** 2

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
   use ed_state_vars, only: edtype,polygontype
   use     misc_coms, only: current_time,simtime
   use      max_dims, only: n_pft,n_dbh
   implicit none
  
   type(edtype)      , target  :: cgrid
   type(polygontype) , pointer :: cpoly
   
   type(simtime)               :: lastmonth
   real                        :: ndaysi,polygon_area_i
   integer                     :: ipy,ipft,dbh
   real                        :: srnonm1
   !----------------------------------------------------------------------!
   ! Finding the inverse of number of days used for this monthly integral !
   !----------------------------------------------------------------------!
   call lastmonthdate(current_time,lastmonth,ndaysi)
   
   !-----------------------------------------------------------!
   ! First normalize the variables previously defined          !
   !-----------------------------------------------------------!
   cgrid%mmean_gpp          = cgrid%mmean_gpp           * ndaysi
   cgrid%mmean_evap         = cgrid%mmean_evap          * ndaysi
   cgrid%mmean_transp       = cgrid%mmean_transp        * ndaysi
   cgrid%mmean_sensible     = cgrid%mmean_sensible      * ndaysi
   cgrid%mmean_sensible_ac  = cgrid%mmean_sensible_ac   * ndaysi
   cgrid%mmean_sensible_gc  = cgrid%mmean_sensible_gc   * ndaysi
   cgrid%mmean_sensible_vc  = cgrid%mmean_sensible_vc   * ndaysi
   cgrid%mmean_nep          = cgrid%mmean_nep           * ndaysi
   cgrid%mmean_soil_temp    = cgrid%mmean_soil_temp     * ndaysi
   cgrid%mmean_soil_water   = cgrid%mmean_soil_water    * ndaysi
   cgrid%mmean_plresp       = cgrid%mmean_plresp        * ndaysi
   cgrid%mmean_rh           = cgrid%mmean_rh            * ndaysi
   cgrid%mmean_leaf_resp    = cgrid%mmean_leaf_resp     * ndaysi
   cgrid%mmean_root_resp    = cgrid%mmean_root_resp     * ndaysi
   cgrid%mmean_growth_resp  = cgrid%mmean_growth_resp   * ndaysi
   cgrid%mmean_storage_resp = cgrid%mmean_storage_resp  * ndaysi
   cgrid%mmean_vleaf_resp   = cgrid%mmean_vleaf_resp    * ndaysi
   cgrid%mmean_gpp_lu       = cgrid%mmean_gpp_lu        * ndaysi
   cgrid%mmean_rh_lu        = cgrid%mmean_rh_lu         * ndaysi
   cgrid%mmean_nep_lu       = cgrid%mmean_nep_lu        * ndaysi
   cgrid%mmean_gpp_dbh      = cgrid%mmean_gpp_dbh       * ndaysi
   cgrid%mmean_lai_pft      = cgrid%mmean_lai_pft       * ndaysi
  
!------------------------------------------------------------------------------------------!
!   Here we convert the sum of squares into standard deviation. The standard deviation can !
! be written in two different ways, and we will use the latter because it doesn't require  !
! previous knowledge of the mean.                                                          !
!              __________________          ____________________________________            !
!             / SUM_i[X_i - Xm]²          /  / SUM_i[X_i²]        \      1                 !
! sigma = \  /  ----------------   =  \  /  |  -----------  - Xm²  | ---------             !
!          \/       N - 1              \/    \      N             /   1 - 1/N              !
!                                                                                          !
! srnonm1 is the square root of 1 / (1 - 1/N)                                              !
!------------------------------------------------------------------------------------------!
   srnonm1 = sqrt(1./(1-ndaysi))
!------------------------------------------------------------------------------------------!
   cgrid%stdev_gpp      = srnonm1 * sqrt(cgrid%stdev_gpp      * ndaysi - cgrid%mmean_gpp      ** 2)
   cgrid%stdev_evap     = srnonm1 * sqrt(cgrid%stdev_evap     * ndaysi - cgrid%mmean_evap     ** 2)
   cgrid%stdev_transp   = srnonm1 * sqrt(cgrid%stdev_transp   * ndaysi - cgrid%mmean_transp   ** 2)
   cgrid%stdev_sensible = srnonm1 * sqrt(cgrid%stdev_sensible * ndaysi - cgrid%mmean_sensible ** 2)
   cgrid%stdev_nep      = srnonm1 * sqrt(cgrid%stdev_nep      * ndaysi - cgrid%mmean_nep      ** 2)
   cgrid%stdev_rh       = srnonm1 * sqrt(cgrid%stdev_rh       * ndaysi - cgrid%mmean_rh       ** 2)
  
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      
      polygon_area_i = 1./sum(cpoly%area)
      do ipft = 1,n_pft
        do dbh =1,n_dbh
          cgrid%agb_pft(ipft,ipy) = cgrid%agb_pft(ipft,ipy) &
               + sum(cpoly%agb(ipft,dbh,:)*cpoly%area)*polygon_area_i
          cgrid%ba_pft(ipft,ipy)  = cgrid%ba_pft(ipft,ipy)  &
               + sum(cpoly%basal_area(ipft,dbh,:)*cpoly%area)*polygon_area_i
        end do
      end do
   end do

   return
end subroutine normalize_ed_monthly_output_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zero_ed_monthly_output_vars(cgrid)
   use ed_state_vars,only:edtype

   implicit none
   type(edtype),target       :: cgrid

   cgrid%mmean_gpp          = 0.
   cgrid%mmean_evap         = 0.
   cgrid%mmean_transp       = 0.
   cgrid%mmean_sensible     = 0.
   cgrid%mmean_sensible_ac  = 0.
   cgrid%mmean_sensible_gc  = 0.
   cgrid%mmean_sensible_vc  = 0.
   cgrid%mmean_nep          = 0.
   cgrid%mmean_soil_temp    = 0.
   cgrid%mmean_soil_water   = 0.
   cgrid%mmean_plresp       = 0.
   cgrid%mmean_rh           = 0.
   cgrid%mmean_leaf_resp    = 0.
   cgrid%mmean_root_resp    = 0.
   cgrid%mmean_growth_resp  = 0.
   cgrid%mmean_storage_resp = 0.
   cgrid%mmean_vleaf_resp   = 0.
   cgrid%mmean_gpp_lu       = 0.
   cgrid%mmean_rh_lu        = 0.
   cgrid%mmean_nep_lu       = 0.
   cgrid%mmean_gpp_dbh      = 0.
   cgrid%mmean_lai_pft      = 0.
   cgrid%agb_pft            = 0.
   cgrid%ba_pft             = 0.
   cgrid%stdev_gpp          = 0.
   cgrid%stdev_evap         = 0.
   cgrid%stdev_transp       = 0.
   cgrid%stdev_sensible     = 0.
   cgrid%stdev_nep          = 0.
   cgrid%stdev_rh           = 0.
 
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
subroutine update_ed_yearly_vars_ar(cgrid)

   use ed_state_vars,only:edtype,polygontype,sitetype,patchtype
   use max_dims, only: n_pft, n_dbh
   use consts_coms, only: pi1
  
   implicit none

   type(edtype),target       :: cgrid
   type(polygontype),pointer :: cpoly
   type(sitetype),pointer    :: csite
   type(patchtype),pointer   :: cpatch
   integer :: ipy,isi,ipa,ico
   real, external :: ed_biomass

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
end subroutine update_ed_yearly_vars_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine zero_ed_yearly_vars_ar(cgrid)

   use max_dims, only: n_pft, n_dbh
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
end subroutine zero_ed_yearly_vars_ar
!==========================================================================================!
!==========================================================================================!
