!==========================================================================================!
!==========================================================================================!
!     This module contains the wrappers for the Runge-Kutta integration scheme.            !
!==========================================================================================!
!==========================================================================================!
module rk4_driver_ar

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      Main driver of short-time scale dynamics of the Runge-Kutta integrator for the   !
   ! land surface model.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_timestep_ar(cgrid,ifm,integration_buff)
      use ed_state_vars  , only : integration_vars_ar  & ! structure
                                , edtype               & ! structure
                                , polygontype          & ! structure
                                , sitetype             & ! structure
                                , patchtype            ! ! structure
      use grid_coms      , only : nzg                  ! ! intent(in)
      use max_dims       , only : n_dbh                ! ! intent(in)
      use misc_coms      , only : dtlsm                ! ! intent(in)
      use consts_coms    , only : umol_2_kgC           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target      :: cgrid
      type(integration_vars_ar) , target      :: integration_buff
      integer                   , intent (in) :: ifm
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
      type(patchtype)           , pointer     :: cpatch
      integer                                 :: ipy,isi,ipa
      integer, dimension(nzg)                 :: ed_ktrans
      real   , dimension(n_dbh)               :: gpp_dbh
      real                                    :: sum_lai_rbi
      real                                    :: gpp
      real                                    :: plant_respiration
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: compute_netrad_ar
      !------------------------------------------------------------------------------------!
      
      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Get velocity for aerodynamic resistance. ----------------------------!
               if (csite%can_temp(ipa) < cpoly%met(isi)%atm_tmp) then
                  cpoly%met(isi)%vels = cpoly%met(isi)%vels_stab
               else
                  cpoly%met(isi)%vels = cpoly%met(isi)%vels_unstab
               end if

               !----- Get photosynthesis, stomatal conductance, and transpiration. --------!
               call canopy_photosynthesis_ar(csite,ipa,cpoly%met(isi)%vels                 &
                          ,cpoly%met(isi)%rhos,cpoly%met(isi)%prss,ed_ktrans               &
                          ,csite%ntext_soil(:,ipa),csite%soil_water(:,ipa)                 &
                          ,csite%soil_fracliq(:,ipa),cpoly%lsl(isi),sum_lai_rbi            &
                          ,cpoly%leaf_aging_factor(:,isi),cpoly%green_leaf_factor(:,isi) )

               !----- Compute root and heterotrophic respiration. -------------------------!
               call soil_respiration_ar(csite,ipa)

               csite%wbudget_precipgain(ipa) = csite%wbudget_precipgain(ipa)               &
                                             + cpoly%met(isi)%pcpg * dtlsm
               csite%ebudget_precipgain(ipa) = csite%ebudget_precipgain(ipa)               &
                                             + cpoly%met(isi)%qpcpg * dtlsm
               csite%ebudget_netrad(ipa)     = csite%ebudget_netrad(ipa)                   &
                                             + compute_netrad_ar(csite,ipa) * dtlsm

               !----- Compute the carbon flux components. ---------------------------------!
               call sum_plant_cfluxes_ar(csite,ipa,gpp,gpp_dbh,plant_respiration)
               csite%co2budget_gpp(ipa)       = csite%co2budget_gpp(ipa) + gpp * dtlsm
               csite%co2budget_gpp_dbh(:,ipa) = csite%co2budget_gpp_dbh(:,ipa) + gpp_dbh(:) *dtlsm
               csite%co2budget_plresp(ipa)    = csite%co2budget_plresp(ipa)                &
                                              + plant_respiration * dtlsm
               csite%co2budget_rh(ipa)        = csite%co2budget_rh(ipa)                    &
                                              + csite%rh(ipa) * dtlsm
               cgrid%cbudget_nep(ipy)         = cgrid%cbudget_nep(ipy)                     &
                                              + cpoly%area(isi) * csite%area(ipa) * dtlsm  &
                                              * (gpp - plant_respiration - csite%rh(ipa))  &
                                              * umol_2_kgC

               !---------------------------------------------------------------------------!
               !    This is the driver for the integration process...                      !
               !---------------------------------------------------------------------------!
               call integrate_patch_ar(csite,ipa,isi,ipy,ifm,integration_buff              &
                                      ,cpoly%met(isi)%rhos,cpoly%met(isi)%vels             &
                                      ,cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv       &
                                      ,cpoly%met(isi)%atm_co2,cpoly%met(isi)%geoht         &
                                      ,cpoly%met(isi)%exner,cpoly%met(isi)%pcpg            &
                                      ,cpoly%met(isi)%qpcpg,cpoly%met(isi)%dpcpg           &
                                      ,cpoly%met(isi)%prss,cpoly%met(isi)%atm_shv          &
                                      ,cpoly%met(isi)%geoht,cpoly%lsl(isi),cgrid%lon(ipy)  &
                                      ,cgrid%lat(ipy))

               !---------------------------------------------------------------------------!
               !    Update the minimum monthly temperature, based on canopy temperature.
               !---------------------------------------------------------------------------!
               if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                  cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
               end if

            end do patchloop
         end do siteloop
      end do polygonloop

      return
   end subroutine rk4_timestep_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_ar(csite,ipa,isi,ipy,ifm,integration_buff,rhos,vels,atm_tmp  &
                                ,rv,atm_co2,zoff,exner,pcpg,qpcpg,dpcpg,prss,atm_shv,geoht &
                                ,lsl,lon,lat)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            & ! structure
                                 , integration_vars_ar  & ! structure
                                 , rk4patchtype         ! ! structure
      use misc_coms       , only : dtlsm                ! ! intent(in)
      use soil_coms       , only : soil_rough           & ! intent(in)
                                 , snow_rough           ! ! intent(in)
      use canopy_air_coms , only : exar                 ! ! intent(in)
      use consts_coms     , only : vonk                 & ! intent(in)
                                 , cp                   ! ! intent(in)
      use rk4_coms        , only : tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4                & ! intent(inout)
                                 , dtrk4i               ! ! intent(inout)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(integration_vars_ar), target     :: integration_buff
      type(sitetype)           , target     :: csite
      integer                  , intent(in) :: ifm
      integer                  , intent(in) :: ipy
      integer                  , intent(in) :: isi
      integer                  , intent(in) :: ipa
      integer                  , intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      real                     , intent(in) :: rhos
      real                     , intent(in) :: vels
      real                     , intent(in) :: atm_tmp
      real                     , intent(in) :: atm_shv
      real                     , intent(in) :: rv
      real                     , intent(in) :: atm_co2
      real                     , intent(in) :: zoff
      real                     , intent(in) :: exner
      real                     , intent(in) :: pcpg
      real                     , intent(in) :: qpcpg
      real                     , intent(in) :: dpcpg
      real                     , intent(in) :: prss
      real                     , intent(in) :: geoht
      real                     , intent(in) :: lon
      real                     , intent(in) :: lat
      !----- Local variables --------------------------------------------------------------!
      type(rk4patchtype)       , pointer    :: initp
      real                                  :: hbeg
      real                                  :: factv
      real                                  :: aux
      real                                  :: zveget
      real                                  :: zdisp
      !----- Locally saved variable -------------------------------------------------------!
      logical, save :: first_time=.true.
      !------------------------------------------------------------------------------------!



      !----- Assigning some constants which will remain the same throughout the run. ------!
      if (first_time) then
         first_time = .false.
         tbeg   = 0.0
         tend   = dtlsm
         dtrk4  = tend - tbeg
         dtrk4i = 1./dtrk4
      end if

      !------------------------------------------------------------------------------------!
      !    Making sure that all buffers are flushed to zero.                               !
      !------------------------------------------------------------------------------------!
      call zero_rk4_patch(integration_buff%initp)
      call zero_rk4_patch(integration_buff%yscal)
      call zero_rk4_patch(integration_buff%y)
      call zero_rk4_patch(integration_buff%dydx)
      call zero_rk4_patch(integration_buff%yerr)
      call zero_rk4_patch(integration_buff%ytemp)
      call zero_rk4_patch(integration_buff%ak2)
      call zero_rk4_patch(integration_buff%ak3)
      call zero_rk4_patch(integration_buff%ak4)
      call zero_rk4_patch(integration_buff%ak5)
      call zero_rk4_patch(integration_buff%ak6)
      call zero_rk4_patch(integration_buff%ak7)
      call zero_rk4_cohort(integration_buff%initp)
      call zero_rk4_cohort(integration_buff%yscal)
      call zero_rk4_cohort(integration_buff%y)
      call zero_rk4_cohort(integration_buff%dydx)
      call zero_rk4_cohort(integration_buff%yerr)
      call zero_rk4_cohort(integration_buff%ytemp)
      call zero_rk4_cohort(integration_buff%ak2)
      call zero_rk4_cohort(integration_buff%ak3)
      call zero_rk4_cohort(integration_buff%ak4)
      call zero_rk4_cohort(integration_buff%ak5)
      call zero_rk4_cohort(integration_buff%ak6)
      call zero_rk4_cohort(integration_buff%ak7)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Set up the integration patch.                                                  !
      !------------------------------------------------------------------------------------!
      initp => integration_buff%initp
      call copy_patch_init_ar(csite,ipa, initp, rhos, lsl)

      !------------------------------------------------------------------------------------!
      !      Initial step size.  Experience has shown that giving this too large a value   !
      ! causes the integrator to fail (e.g., soil layers become supersaturated).           !
      !------------------------------------------------------------------------------------!
      hbeg = csite%htry(ipa)


      !------------------------------------------------------------------------------------!
      !      This is Bob Walko's recommended way of calculating the resistance.  Note that !
      ! temperature, not potential temperature, is input here.                             !
      !------------------------------------------------------------------------------------!
      initp%rough = max(soil_rough, csite%veg_rough(ipa)) * (1.0 - csite%snowfac(ipa))     &
                  + snow_rough
      zveget      = csite%veg_height(ipa) * (1.0 - csite%snowfac(ipa))
      zdisp       = 0.63 * zveget

      !------------------------------------------------------------------------------------!
      !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm, !
      ! so they must be zeroed at each call.                                               !
      !------------------------------------------------------------------------------------!
      initp%upwp = 0.
      initp%tpwp = 0.
      initp%rpwp = 0.
      initp%wpwp = 0.

      !----- Finding the characteristic scales (a.k.a. starts). ---------------------------!
      call ed_stars(atm_tmp,atm_shv,atm_co2,initp%can_temp,geoht,vels,initp%rough          &
                   ,initp%ustar,initp%rstar,initp%tstar,initp%cstar,initp%can_shv          &
                   ,initp%can_co2)
      !----- Converting T-star to theta-star. ---------------------------------------------!
      initp%tstar = initp%tstar * cp / exner  

      !----- Finding the aerodynamic resistance due to vegetation. ------------------------!
      factv        = 1.0 / (vonk * initp%ustar)
      aux          = exp(exar * (1. - (zdisp + initp%rough) / zveget))
      initp%rasveg = factv * zveget / (exar * (zveget - zdisp)) * (exp(exar) - aux)

      !----- Go into the ODE integrator. --------------------------------------------------!
      call odeint_ar(hbeg,csite,ipa,isi,ipy,ifm,integration_buff,rhos,vels,atm_tmp,atm_shv &
                    ,atm_co2,geoht,exner, pcpg, qpcpg, dpcpg, prss, lsl)

      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = rhos*initp%upwp/dtlsm
      initp%tpwp = rhos*initp%tpwp/dtlsm
      initp%rpwp = rhos*initp%rpwp/dtlsm
      initp%wpwp = rhos*initp%wpwp/dtlsm
      
      
      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp_ar(tend-tbeg,initp,csite,ipa,isi,ipy,lsl,atm_tmp,atm_shv,atm_co2   &
                          ,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg,lon,lat)

      return
   end subroutine integrate_patch_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp_ar(hdid,initp,csite,ipa,isi,ipy,lsl,atm_tmp,atm_shv,atm_co2     &
                             ,prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg,lon,lat)
      use ed_state_vars        , only : sitetype          & ! structure
                                      , patchtype         & ! structure
                                      , rk4patchtype      & ! structure
                                      , edgrid_g          ! ! structure
      use consts_coms          , only : day_sec           & ! intent(in)
                                      , t3ple             ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics  ! ! intent(in)
      use soil_coms            , only : soil              & ! intent(in)
                                      , slz               & ! intent(in)
                                      , min_sfcwater_mass ! ! intent(in) 
      use ed_misc_coms         , only : fast_diagnostics  ! ! intent(in)
      use grid_coms            , only : nzg               & ! intent(in)
                                      , nzs               ! ! intent(in)
      use canopy_radiation_coms, only : veg_temp_min      ! ! intent(in)
      use rk4_coms             , only : rk4max_veg_temp   ! ! intent(in)
      use therm_lib            , only : qwtk              ! ! subroutine
      use ed_therm_lib         , only : calc_hcapveg      & ! function
                                      , ed_grndvap        ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target     :: initp
      type(sitetype)    , target     :: csite
      integer           , intent(in) :: ipa
      integer           , intent(in) :: ipy
      integer           , intent(in) :: isi
      integer           , intent(in) :: lsl
      real              , intent(in) :: atm_tmp
      real              , intent(in) :: atm_shv
      real              , intent(in) :: atm_co2
      real              , intent(in) :: prss,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg
      real              , intent(in) :: lon,lat
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer    :: cpatch
      integer                        :: ico
      integer                        :: k
      integer                        :: kclosest
      integer                        :: ksn
      integer                        :: nsoil
      integer                        :: nlsw1
      real                           :: hdid
      real                           :: veg_fliq
      real                           :: available_water
      real                           :: surface_temp
      real                           :: surface_fliq
      !----- Local contants ---------------------------------------------------------------!
      real, parameter                :: tendays_sec=10.*day_sec
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.                                                   !
      !------------------------------------------------------------------------------------!
      csite%can_temp(ipa) = initp%can_temp
      csite%can_shv(ipa)  = initp%can_shv
      csite%can_co2(ipa)  = initp%can_co2

      csite%ustar(ipa)    = initp%ustar
      csite%tstar(ipa)    = initp%tstar
      csite%rstar(ipa)    = initp%rstar
      csite%cstar(ipa)    = initp%cstar

      csite%upwp(ipa)     = initp%upwp
      csite%wpwp(ipa)     = initp%wpwp
      csite%tpwp(ipa)     = initp%tpwp
      csite%rpwp(ipa)     = initp%rpwp

      !------------------------------------------------------------------------------------!
      !    These variables are fast scale fluxes, and they may not be allocated, so just   !
      ! check this before copying.                                                         !
      !------------------------------------------------------------------------------------!
      if(fast_diagnostics) then
         csite%wbudget_loss2atm(ipa)     = initp%wbudget_loss2atm
         csite%ebudget_loss2atm(ipa)     = initp%ebudget_loss2atm
         csite%co2budget_loss2atm(ipa)   = initp%co2budget_loss2atm
         csite%ebudget_latent(ipa)       = initp%ebudget_latent
         csite%avg_vapor_vc(ipa)         = initp%avg_vapor_vc
         csite%avg_dew_cg(ipa)           = initp%avg_dew_cg
         csite%avg_vapor_gc(ipa)         = initp%avg_vapor_gc
         csite%avg_wshed_vg(ipa)         = initp%avg_wshed_vg
         csite%avg_vapor_ac(ipa)         = initp%avg_vapor_ac
         csite%avg_transp(ipa)           = initp%avg_transp
         csite%avg_evap(ipa)             = initp%avg_evap
         csite%avg_netrad(ipa)           = initp%avg_netrad
         csite%aux(ipa)                  = initp%aux
         csite%avg_sensible_vc(ipa)      = initp%avg_sensible_vc
         csite%avg_sensible_2cas(ipa)    = initp%avg_sensible_2cas
         csite%avg_qwshed_vg(ipa)        = initp%avg_qwshed_vg
         csite%avg_sensible_gc(ipa)      = initp%avg_sensible_gc
         csite%avg_sensible_ac(ipa)      = initp%avg_sensible_ac
         csite%avg_sensible_tot(ipa)     = initp%avg_sensible_tot
         csite%avg_carbon_ac(ipa)        = initp%avg_carbon_ac
         do k = lsl, nzg
            csite%avg_sensible_gg(k,ipa) = initp%avg_sensible_gg(k)
            csite%avg_smoist_gg(k,ipa)   = initp%avg_smoist_gg(k)
            csite%avg_smoist_gc(k,ipa)   = initp%avg_smoist_gc(k)
            csite%aux_s(k,ipa)           = initp%aux_s(k)
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following is not a pure diagnostic, it is used for phenology and mortality !
      ! functions, preserve this variable and its dependencies in all contexts.            !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! paw_avg10d - 10-day average of plant available water.                              !
      !     I don't think this is currently used, but it may be turned on for drought-     !
      ! -related  phenology.                                                               !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         available_water = 0.0
         do k = cpatch%krdepth(ico), nzg - 1
            nsoil = csite%ntext_soil(k,ipa)
            available_water = available_water                                              &
                            + real( (initp%soil_water(k) - dble(soil(nsoil)%soilcp))       &
                                  * dble(slz(k+1)-slz(k))                                  &
                                  / dble(soil(nsoil)%slmsts - soil(nsoil)%soilcp) )
         end do
         nsoil = csite%ntext_soil(nzg,ipa)
         available_water = available_water                                                 &
                         + real( (initp%soil_water(nzg) - dble(soil(nsoil)%soilcp))        &
                               * dble(-1.0*slz(nzg))                                       &
                               / dble(soil(nsoil)%slmsts -soil(nsoil)%soilcp)) 
         available_water = available_water / (-1.0*slz(cpatch%krdepth(ico)))

         cpatch%paw_avg10d(ico) = cpatch%paw_avg10d(ico)*(1.0-hdid/tendays_sec)            &
                                + available_water*hdid/tendays_sec
      end do

      
      do k = lsl, nzg
         csite%soil_water(k,ipa)   = initp%soil_water(k)
         csite%soil_energy(k,ipa)  = initp%soil_energy(k)
         csite%soil_tempk(k,ipa)   = initp%soil_tempk(k)
         csite%soil_fracliq(k,ipa) = initp%soil_fracliq(k)
      end do
      

      !------------------------------------------------------------------------------------!
      !    Surface water energy is computed in J/m² inside the integrator. Converting it   !
      ! back to J/kg in the layers that surface water/snow still exists.                   !
      !------------------------------------------------------------------------------------!
      csite%nlev_sfcwater(ipa) = initp%nlev_sfcwater
      do k = 1, csite%nlev_sfcwater(ipa)
         csite%sfcwater_depth(k,ipa)   = initp%sfcwater_depth(k)
         csite%sfcwater_mass(k,ipa)    = initp%sfcwater_mass(k)
         csite%sfcwater_tempk(k,ipa)   = initp%sfcwater_tempk(k)
         csite%sfcwater_fracliq(k,ipa) = initp%sfcwater_fracliq(k)
         csite%sfcwater_energy(k,ipa)  = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
      end do
      !------------------------------------------------------------------------------------!
      !    For the layers that no longer exist, assign zeroes for prognostic variables,    !
      ! and something for temperature and liquid fraction (just to avoid singularities,    !
      ! and funny numbers in the output, but these values are meaningless and should never !
      ! be used).                                                                          !
      !------------------------------------------------------------------------------------!
      do k = csite%nlev_sfcwater(ipa)+1,nzs
         csite%sfcwater_energy(k,ipa)  = 0.
         csite%sfcwater_mass(k,ipa)    = 0.
         csite%sfcwater_depth(k,ipa)   = 0.
         if (k == 1) then
            csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
            csite%sfcwater_tempk(k,ipa)   = csite%soil_tempk(nzg,ipa)
         else
            csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
            csite%sfcwater_tempk(k,ipa)   = csite%sfcwater_tempk(k-1,ipa)
         end if
      end do
      !------------------------------------------------------------------------------------!



      
      !------------------------------------------------------------------------------------!
      !     Cohort variables.  Here we must check whether the cohort was really solved or  !
      ! it was skipped after being flagged as "unsafe".  Here the reason why it was flag-  !
      ! ged as such matters.                                                               !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         if (initp%solvable(ico)) then
            !------------------------------------------------------------------------------!
            !    The cohort was solved, update internal energy and water, and re-calculate !
            ! temperature.  Note that energy may need to be scaled back.                   !
            !------------------------------------------------------------------------------!
            cpatch%veg_water(ico)  = initp%veg_water(ico)
            cpatch%veg_energy(ico) = initp%veg_energy(ico)                                 &
                                   + (cpatch%hcapveg(ico)-initp%hcapveg(ico))              &
                                   * initp%veg_temp(ico)
            call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)     &
                     ,cpatch%veg_temp(ico),veg_fliq)
         elseif (cpatch%hite(ico) <=  csite%total_snow_depth(ipa)) then
            !------------------------------------------------------------------------------!
            !    For plants buried in snow, fix the leaf temperature to the snow temper-   !
            ! ature of the layer that is the closest to the leaves.                        !
            !------------------------------------------------------------------------------!
            kclosest = 1
            do k = csite%nlev_sfcwater(ipa), 1, -1
               if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
            end do
            cpatch%veg_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
            cpatch%veg_water(ico)  = 0.
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
         else
            !------------------------------------------------------------------------------!
            !     For plants with minimal foliage or very sparse patches, fix the leaf     !
            ! temperature to the canopy air space and force veg_water to be zero.          !
            !------------------------------------------------------------------------------!
            cpatch%veg_temp(ico)   = csite%can_temp(ipa)
            cpatch%veg_water(ico)  = 0. 
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
         end if

         !---------------------------------------------------------------------------------!
         !     Final sanity check.  This should be removed soon, since it should never     !
         ! happen (well, if this still happens, then it's a bug, and we should remove the  !
         ! bug first...).                                                                  !
         !---------------------------------------------------------------------------------!
         if (cpatch%veg_temp(ico) < veg_temp_min .or.                                      &
             cpatch%veg_temp(ico) > rk4max_veg_temp   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL VEG_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',lat
            write (unit=*,fmt='(a,1x,i6)')     ' + POLYGON:      ',ipy
            write (unit=*,fmt='(a,1x,i6)')     ' + SITE:         ',isi
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,es12.5)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.5)') '   - ENERGY:     ',initp%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',initp%veg_water(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',initp%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',initp%veg_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',initp%hcapveg(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.5)') '   - ENERGY:     ',cpatch%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',cpatch%veg_water(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',cpatch%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',veg_fliq
            write (unit=*,fmt='(a,1x,es12.5)') '   - WATER:      ',cpatch%hcapveg(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_patch_ar(initp, csite,ipa, lsl,atm_tmp,atm_shv,atm_co2,prss           &
                               ,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do


      ksn   = csite%nlev_sfcwater(ipa)
      nsoil = csite%ntext_soil(nzg,ipa)
      nlsw1 = max(1, ksn)
      call ed_grndvap(ksn,nsoil,csite%soil_water(nzg,ipa),csite%soil_energy(nzg,ipa)       &
                     ,csite%sfcwater_energy(nlsw1,ipa),rhos,csite%can_shv(ipa)             &
                     ,csite%ground_shv(ipa),csite%surface_ssh(ipa),surface_temp            &
                     ,surface_fliq)
      return
   end subroutine initp2modelp_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Currently not in use.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_atm_fluxes_ar(csite,cpoly,ipa,isi)
      use ed_state_vars , only : polygontype  & ! structure
                               , sitetype     & ! structure
                               , patchtype    ! ! structure
      use consts_coms   , only : cpi          ! ! intent(in)
      implicit none
    
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype) , target     :: cpoly
      type(sitetype)    , target     :: csite
      integer           , intent(in) :: ipa
      integer           , intent(in) :: isi
      !------------------------------------------------------------------------------------!

      call fatal_error('Decide how to set vels in canopy_atm_fluxes.'                      &
                     &,'canopy_atm_fluxes_ar','rk4_driver.f90')

      !----- Calculate turbulent fluxes between atmosphere and canopy. --------------------!
      !    pis = cpoly%pi0 * cpi
      !    thetacan = pss%can_temp / pis
      !    if(thetacan.lt.cpoly%theta)then
      !       cpoly%vels = cpoly%vels_stab
      !    else
      !       cpoly%vels = cpoly%vels_unstab
      !    endif

      return
   end subroutine canopy_atm_fluxes_ar
   !=======================================================================================!
   !=======================================================================================! 
end module rk4_driver_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This function computes the total water stored in the system, in kg/m2.                 !
!   (soil + temporary pools + canopy air space + leaf surface).                            !
!------------------------------------------------------------------------------------------!
real function compute_water_storage_ar(csite, lsl, rhos,ipa)
   use ed_state_vars , only : sitetype   & ! structure
                            , patchtype  ! ! structure
   use grid_coms     , only : nzg        ! ! intent(in)
   use soil_coms     , only : dslz       ! ! intent(in)
   use consts_coms   , only : wdns       ! ! intent(in)
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
   !---------------------------------------------------------------------------------------!

   compute_water_storage_ar = 0.0
   cpatch => csite%patch(ipa)

   !----- 1. Adding the water stored in the soil. -----------------------------------------!
   do k = lsl, nzg
      compute_water_storage_ar = compute_water_storage_ar                                  &
                               + real(csite%soil_water(k,ipa)) * dslz(k) * wdns
   end do
   !----- 2. Adding the water stored in the temporary surface water/snow. -----------------!
   do k = 1, csite%nlev_sfcwater(ipa)
      compute_water_storage_ar = compute_water_storage_ar + csite%sfcwater_mass(k,ipa)
   end do
   !----- 3. Adding the water vapour floating in the canopy air space. --------------------!
   compute_water_storage_ar = compute_water_storage_ar                                     &
                            + csite%can_shv(ipa) * csite%veg_height(ipa) * rhos
   !----- 4. Adding the water over the leaf surface. --------------------------------------!
   do ico = 1,cpatch%ncohorts
      compute_water_storage_ar = compute_water_storage_ar + cpatch%veg_water(ico)
   end do

   return
end function compute_water_storage_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computs the total net radiation, by adding the radiation that interacts !
! with the different surfaces.                                                             !
!------------------------------------------------------------------------------------------!
real function compute_netrad_ar(csite,ipa)
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

   compute_netrad_ar = 0.0
   !----- 1. Adding the ground components -------------------------------------------------!
   compute_netrad_ar = csite%rshort_g(ipa) + csite%rlong_g(ipa) + csite%rlong_s(ipa)
   !----- 2. Adding the shortwave radiation that reaches each snow/water layer ------------!
   do k = 1, csite%nlev_sfcwater(ipa)
      compute_netrad_ar = compute_netrad_ar + csite%rshort_s(k,ipa)
   end do
   !----- 3. Adding the radiation components that interact with leaves. -------------------!
   do ico = 1,cpatch%ncohorts
      compute_netrad_ar = compute_netrad_ar + cpatch%rshort_v(ico) + cpatch%rlong_v(ico)
   end do
   return
end function compute_netrad_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computs the total net radiation, by adding the radiation that interacts !
! with the different surfaces.  The result is given in J/m2.                               !
!------------------------------------------------------------------------------------------!
real function compute_energy_storage_ar(csite, lsl, rhos, ipa)
   use ed_state_vars        , only : sitetype   & ! structure
                                   , patchtype  ! ! structure
   use grid_coms            , only : nzg        ! ! intent(in)
   use soil_coms            , only : dslz       ! ! intent(in)
   use consts_coms          , only : cp         & ! intent(in)
                                   , cliq       & ! intent(in)
                                   , cice       & ! intent(in)
                                   , alli       & ! intent(in)
                                   , t3ple      ! ! intent(in)
   use rk4_coms             , only : toosparse  ! ! intent(in)
   use canopy_radiation_coms, only : lai_min    ! ! intent(in)
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
   !---------------------------------------------------------------------------------------!

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
   cas_storage = cp * rhos * csite%veg_height(ipa) * csite%can_temp(ipa)

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
      !if(csite%lai(ipa)   > lai_min                     .and.                             &
      !   cpatch%hite(ico) > csite%total_snow_depth(ipa) .and.                             &
      !   (.not. toosparse)                                   ) then
      !   veg_storage = veg_storage + cpatch%veg_energy(ico)
      !end if
   end do
 
   !----- 5. Integrating the total energy in ED. ------------------------------------------!
   compute_energy_storage_ar = soil_storage + sfcwater_storage + cas_storage               &
                             + veg_storage

   return
end function compute_energy_storage_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the carbon flux terms.                                       !
!------------------------------------------------------------------------------------------!
subroutine sum_plant_cfluxes_ar(csite,ipa, gpp, gpp_dbh,plresp)
   use ed_state_vars        , only : sitetype    & ! structure
                                   , patchtype   ! ! structure
   use consts_coms          , only : day_sec     & ! intent(in)
                                   , umol_2_kgC  ! ! intent(in)
   use canopy_radiation_coms, only : lai_min     ! ! intent(in)
   use max_dims             , only : n_dbh
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
      if (cpatch%lai(ico) > lai_min) then
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
      ! ive", we must convert them from umol/plant/day to kgC/m2/s.                        !
      !------------------------------------------------------------------------------------!
      sresp  = sresp                                                                       &
             + ( cpatch%growth_respiration(ico) + cpatch%storage_respiration(ico)          &
               + cpatch%vleaf_respiration(ico))                                            &
               * cpatch%nplant(ico) / (day_sec * umol_2_kgC)
   end do
   !----- Plant respiration is the sum between alive and structural. ----------------------!
   plresp = lrresp + sresp
   return
end subroutine sum_plant_cfluxes_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function computes the co2 stored in the canopy air space from ppm to kgC/m2.     !
!------------------------------------------------------------------------------------------!
real function compute_co2_storage_ar(csite, rhos, ipa)
   use ed_state_vars, only : sitetype ! ! structure
   use consts_coms  , only : mmdryi   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target      :: csite
   integer               , intent(in)  :: ipa
   real                  , intent(in)  :: rhos
   !---------------------------------------------------------------------------------------!

   compute_co2_storage_ar = csite%can_co2(ipa) * mmdryi * rhos * csite%veg_height(ipa)

   return
end function compute_co2_storage_ar
!==========================================================================================!
!==========================================================================================!
