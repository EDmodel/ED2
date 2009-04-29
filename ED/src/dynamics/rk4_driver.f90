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
   subroutine rk4_timestep_ar(cgrid,ifm)
      use rk4_coms               , only : integration_vars_ar  & ! structure
                                        , integration_buff     & ! structure
                                        , rk4patchtype         & ! structure
                                        , rk4met               ! ! intent(out)
      use ed_state_vars          , only : edtype               & ! structure
                                        , polygontype          & ! structure
                                        , sitetype             & ! structure
                                        , patchtype            ! ! structure
      use grid_coms              , only : nzg                  ! ! intent(in)
      use max_dims               , only : n_dbh                ! ! intent(in)
      use misc_coms              , only : dtlsm                ! ! intent(in)
      use consts_coms            , only : umol_2_kgC           ! ! intent(in)
      use canopy_struct_dynamics , only : canopy_turbulence ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target      :: cgrid
      integer                   , intent (in) :: ifm
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
      type(patchtype)           , pointer     :: cpatch
      integer                                 :: ipy,isi,ipa
      integer, dimension(nzg)                 :: ed_ktrans
      real                                    :: time_py_start,time_py_spent
      real   , dimension(n_dbh)               :: gpp_dbh
      real                                    :: sum_lai_rbi
      real                                    :: gpp
      real                                    :: plant_respiration
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: compute_netrad_ar
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!
      
      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         time_py_start = walltime(0.) 

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

               !---------------------------------------------------------------------------!
               !    Copy the meteorological variables to the rk4met structure.             !
               !---------------------------------------------------------------------------!
               call copy_met_2_rk4met(cpoly%met(isi)%rhos,cpoly%met(isi)%vels              &
                                     ,cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv        &
                                     ,cpoly%met(isi)%atm_co2,cpoly%met(isi)%geoht          &
                                     ,cpoly%met(isi)%exner,cpoly%met(isi)%pcpg             &
                                     ,cpoly%met(isi)%qpcpg,cpoly%met(isi)%dpcpg            &
                                     ,cpoly%met(isi)%prss,cpoly%met(isi)%geoht             &
                                     ,cpoly%lsl(isi),cgrid%lon(ipy),cgrid%lat(ipy))


               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init_ar(csite,ipa,integration_buff%initp)



               !---------------------------------------------------------------------------!
               !     Calculate the canopy geometry, and the scalar transport coefficients. !
               !---------------------------------------------------------------------------!
               call canopy_turbulence(csite,integration_buff%initp,isi,ipa,.true.)



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
               csite%co2budget_gpp_dbh(:,ipa) = csite%co2budget_gpp_dbh(:,ipa)             & 
                                              + gpp_dbh(:) *dtlsm
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
               call integrate_patch_ar(csite,ipa,isi,ipy,ifm,integration_buff)

               !---------------------------------------------------------------------------!
               !    Update the minimum monthly temperature, based on canopy temperature.
               !---------------------------------------------------------------------------!
               if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                  cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
               end if

            end do patchloop
         end do siteloop

         time_py_spent = walltime(time_py_start)
         cgrid%walltime_py(ipy) = cgrid%walltime_py(ipy) + dble(time_py_spent)

      end do polygonloop

      return
   end subroutine rk4_timestep_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_ar(csite,ipa,isi,ipy,ifm,integration_buff)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use misc_coms       , only : dtlsm                ! ! intent(in)
      use soil_coms       , only : soil_rough           & ! intent(in)
                                 , snow_rough           ! ! intent(in)
      use canopy_air_coms , only : exar8                ! ! intent(in)
      use consts_coms     , only : vonk8                & ! intent(in)
                                 , cp8                  & ! intent(in)
                                 , cpi8                 ! ! intent(in)
      use rk4_coms        , only : integration_vars_ar  & ! structure
                                 , rk4patchtype         & ! structure
                                 , rk4met               & ! intent(inout)
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4                & ! intent(inout)
                                 , dtrk4i               & ! intent(inout)
                                 , ibranch_thermo       & ! intent(in)
                                 , effarea_water        & ! intent(out)
                                 , effarea_heat         ! ! intent(out)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(integration_vars_ar), target     :: integration_buff
      type(sitetype)           , target     :: csite
      integer                  , intent(in) :: ifm
      integer                  , intent(in) :: ipy
      integer                  , intent(in) :: isi
      integer                  , intent(in) :: ipa
      !----- Local variables --------------------------------------------------------------!
      type(rk4patchtype)       , pointer    :: initp
      real(kind=8)                          :: hbeg
      !----- Locally saved variable -------------------------------------------------------!
      logical                  , save       :: first_time=.true.
      !------------------------------------------------------------------------------------!

      !----- Making an alias. -------------------------------------------------------------!
      initp => integration_buff%initp

      !----- Assigning some constants which will remain the same throughout the run. ------!
      if (first_time) then
         first_time = .false.
         tbeg   = 0.d0
         tend   = dble(dtlsm)
         dtrk4  = tend - tbeg
         dtrk4i = 1.d0/dtrk4
         
         !---------------------------------------------------------------------------------!
         !    The area factor for heat and water exchange between canopy and vegetation is !
         ! applied only on LAI, and it depends on how we are considering the branches and  !
         ! twigs.  If their area isn't explicitly defined, we add a 0.2 factor to the      !
         ! area because BAI will be 0.  Otherwise, we don't add anything to the LAI, and   !
         ! let BAI to do the job.                                                          !
         !---------------------------------------------------------------------------------!
         select case (ibranch_thermo)
         case (0)
            effarea_water = 1.2d0
            effarea_heat  = 2.2d0
         case (1)
            effarea_water = 1.0d0
            effarea_heat  = 2.0d0
         end select
      end if

      !------------------------------------------------------------------------------------!
      !      Initial step size.  Experience has shown that giving this too large a value   !
      ! causes the integrator to fail (e.g., soil layers become supersaturated).           !
      !------------------------------------------------------------------------------------!
      hbeg = dble(csite%htry(ipa))


      !------------------------------------------------------------------------------------!
      !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm, !
      ! so they must be zeroed at each call.                                               !
      !------------------------------------------------------------------------------------!
      initp%upwp = 0.d0
      initp%tpwp = 0.d0
      initp%qpwp = 0.d0
      initp%wpwp = 0.d0

      !----- Go into the ODE integrator. --------------------------------------------------!
      call odeint_ar(hbeg,csite,ipa,isi,ipy,ifm,integration_buff)

      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = rk4met%rhos*initp%upwp * dtrk4i
      initp%tpwp = rk4met%rhos*initp%tpwp * dtrk4i
      initp%qpwp = rk4met%rhos*initp%qpwp * dtrk4i
      initp%wpwp = rk4met%rhos*initp%wpwp * dtrk4i
      
      
      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp_ar(tend-tbeg,initp,csite,ipa,isi,ipy)

      return
   end subroutine integrate_patch_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp_ar(hdid,initp,csite,ipa,isi,ipy)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4met               & ! intent(in)
                                      , rk4min_veg_temp      & ! intent(in)
                                      , rk4max_veg_temp      & ! intent(in)
                                      , tiny_offset          ! ! intent(in) 
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            & ! structure
                                      , edgrid_g             ! ! structure
      use consts_coms          , only : day_sec              & ! intent(in)
                                      , t3ple8               ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , slz8                 ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     ! ! intent(in)
      use grid_coms            , only : nzg                  & ! intent(in)
                                      , nzs                  ! ! intent(in)
      use therm_lib            , only : qwtk                 ! ! subroutine
      use ed_therm_lib         , only : ed_grndvap           ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target     :: initp
      type(sitetype)    , target     :: csite
      real(kind=8)      , intent(in) :: hdid
      integer           , intent(in) :: ipa
      integer           , intent(in) :: ipy
      integer           , intent(in) :: isi
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer    :: cpatch
      integer                        :: ico
      integer                        :: k
      integer                        :: kclosest
      integer                        :: ksn
      integer                        :: nsoil
      integer                        :: nlsw1
      real(kind=8)                   :: available_water
      real(kind=8)                   :: tmp_energy
      real                           :: surface_temp
      real                           :: surface_fliq
      !----- Local contants ---------------------------------------------------------------!
      real        , parameter        :: tendays_sec=10.*day_sec
      !----- External function ------------------------------------------------------------!
      real        , external         :: sngloff
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.  All floating point variables are converted back  !
      ! to single precision.                                                               !
      !------------------------------------------------------------------------------------!
      csite%can_temp(ipa) = sngloff(initp%can_temp,tiny_offset)
      csite%can_shv(ipa)  = sngloff(initp%can_shv ,tiny_offset)
      csite%can_co2(ipa)  = sngloff(initp%can_co2 ,tiny_offset)

      csite%ustar(ipa)    = sngloff(initp%ustar   ,tiny_offset)
      csite%tstar(ipa)    = sngloff(initp%tstar   ,tiny_offset)
      csite%qstar(ipa)    = sngloff(initp%qstar   ,tiny_offset)
      csite%cstar(ipa)    = sngloff(initp%cstar   ,tiny_offset)

      csite%upwp(ipa)     = sngloff(initp%upwp    ,tiny_offset)
      csite%wpwp(ipa)     = sngloff(initp%wpwp    ,tiny_offset)
      csite%tpwp(ipa)     = sngloff(initp%tpwp    ,tiny_offset)
      csite%qpwp(ipa)     = sngloff(initp%qpwp    ,tiny_offset)

      !------------------------------------------------------------------------------------!
      !    These variables are fast scale fluxes, and they may not be allocated, so just   !
      ! check this before copying.                                                         !
      !------------------------------------------------------------------------------------!
      if(fast_diagnostics) then
         csite%wbudget_loss2atm(ipa)     =sngloff(initp%wbudget_loss2atm  ,tiny_offset)
         csite%ebudget_loss2atm(ipa)     =sngloff(initp%ebudget_loss2atm  ,tiny_offset)
         csite%co2budget_loss2atm(ipa)   =sngloff(initp%co2budget_loss2atm,tiny_offset)
         csite%ebudget_latent(ipa)       =sngloff(initp%ebudget_latent    ,tiny_offset)
         csite%avg_vapor_vc(ipa)         =sngloff(initp%avg_vapor_vc      ,tiny_offset)
         csite%avg_dew_cg(ipa)           =sngloff(initp%avg_dew_cg        ,tiny_offset)
         csite%avg_vapor_gc(ipa)         =sngloff(initp%avg_vapor_gc      ,tiny_offset)
         csite%avg_wshed_vg(ipa)         =sngloff(initp%avg_wshed_vg      ,tiny_offset)
         csite%avg_vapor_ac(ipa)         =sngloff(initp%avg_vapor_ac      ,tiny_offset)
         csite%avg_transp(ipa)           =sngloff(initp%avg_transp        ,tiny_offset)
         csite%avg_evap(ipa)             =sngloff(initp%avg_evap          ,tiny_offset)
         csite%avg_drainage(ipa)         =sngloff(initp%avg_drainage      ,tiny_offset)
         csite%avg_netrad(ipa)           =sngloff(initp%avg_netrad        ,tiny_offset)
         csite%avg_sensible_vc(ipa)      =sngloff(initp%avg_sensible_vc   ,tiny_offset)
         csite%avg_sensible_2cas(ipa)    =sngloff(initp%avg_sensible_2cas ,tiny_offset)
         csite%avg_qwshed_vg(ipa)        =sngloff(initp%avg_qwshed_vg     ,tiny_offset)
         csite%avg_sensible_gc(ipa)      =sngloff(initp%avg_sensible_gc   ,tiny_offset)
         csite%avg_sensible_ac(ipa)      =sngloff(initp%avg_sensible_ac   ,tiny_offset)
         csite%avg_sensible_tot(ipa)     =sngloff(initp%avg_sensible_tot  ,tiny_offset)
         csite%avg_carbon_ac(ipa)        =sngloff(initp%avg_carbon_ac     ,tiny_offset)
         do k = rk4met%lsl, nzg
            csite%avg_sensible_gg(k,ipa) =sngloff(initp%avg_sensible_gg(k),tiny_offset)
            csite%avg_smoist_gg(k,ipa)   =sngloff(initp%avg_smoist_gg(k)  ,tiny_offset)
            csite%avg_smoist_gc(k,ipa)   =sngloff(initp%avg_smoist_gc(k)  ,tiny_offset)
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
      ! paw_avg - 10-day average of plant available water.                              !
      !     I don't think this is currently used, but it may be turned on for drought-     !
      ! -related  phenology.                                                               !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         available_water = 0.d0
         do k = cpatch%krdepth(ico), nzg - 1
            nsoil = csite%ntext_soil(k,ipa)
            available_water = available_water                                              &
                            + (initp%soil_water(k) - soil8(nsoil)%soilcp)                  &
                            * (slz8(k+1)-slz8(k))                                          &
                            / (soil8(nsoil)%slmsts - soil8(nsoil)%soilcp)
         end do
         nsoil = csite%ntext_soil(nzg,ipa)
         available_water = available_water                                                 &
                         + (initp%soil_water(nzg) - soil8(nsoil)%soilcp)                   &
                         * (-1.d0*slz8(nzg))                                               &
                         / (soil8(nsoil)%slmsts -soil8(nsoil)%soilcp) 
         available_water = available_water / (-1.d0*slz8(cpatch%krdepth(ico)))


         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)      &
                                + sngl(available_water)*sngl(hdid)/tendays_sec
      end do

      
      do k = rk4met%lsl, nzg
         csite%soil_water(k,ipa)   = sngloff(initp%soil_water(k)  ,tiny_offset)
         csite%soil_energy(k,ipa)  = sngloff(initp%soil_energy(k) ,tiny_offset)
         csite%soil_tempk(k,ipa)   = sngloff(initp%soil_tempk(k)  ,tiny_offset)
         csite%soil_fracliq(k,ipa) = sngloff(initp%soil_fracliq(k),tiny_offset)
      end do
      

      !------------------------------------------------------------------------------------!
      !    Surface water energy is computed in J/m² inside the integrator. Converting it   !
      ! back to J/kg in the layers that surface water/snow still exists.                   !
      !------------------------------------------------------------------------------------!
      csite%nlev_sfcwater(ipa) = initp%nlev_sfcwater
      do k = 1, csite%nlev_sfcwater(ipa)
         csite%sfcwater_depth(k,ipa)   = sngloff(initp%sfcwater_depth(k)   ,tiny_offset)
         csite%sfcwater_mass(k,ipa)    = sngloff(initp%sfcwater_mass(k)    ,tiny_offset)
         csite%sfcwater_tempk(k,ipa)   = sngloff(initp%sfcwater_tempk(k)   ,tiny_offset)
         csite%sfcwater_fracliq(k,ipa) = sngloff(initp%sfcwater_fracliq(k) ,tiny_offset)
         tmp_energy                    = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
         csite%sfcwater_energy(k,ipa)  = sngloff(tmp_energy                ,tiny_offset)
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
            cpatch%veg_water(ico)  = sngloff(initp%veg_water(ico),tiny_offset)
            tmp_energy             = initp%veg_energy(ico)                                 &
                                   + (dble(cpatch%hcapveg(ico))-initp%hcapveg(ico))        &
                                   * initp%veg_temp(ico)

            cpatch%veg_energy(ico) = sngloff(tmp_energy,tiny_offset)
            call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)     &
                     ,cpatch%veg_temp(ico),cpatch%veg_fliq(ico))
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
            cpatch%veg_fliq(ico)   = 0.
            cpatch%veg_water(ico)  = 0.
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
         else
            !------------------------------------------------------------------------------!
            !     For plants with minimal foliage or very sparse patches, fix the leaf     !
            ! temperature to the canopy air space and force veg_water to be zero.          !
            !------------------------------------------------------------------------------!
            cpatch%veg_temp(ico)   = csite%can_temp(ipa)
            cpatch%veg_fliq(ico)   = 0.
            cpatch%veg_water(ico)  = 0. 
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
         end if

         !---------------------------------------------------------------------------------!
         !     Final sanity check.  This should be removed soon, since it should never     !
         ! happen (well, if this still happens, then it's a bug, and we should remove the  !
         ! bug first...).                                                                  !
         !---------------------------------------------------------------------------------!
         if (cpatch%veg_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%veg_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL VEG_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4met%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4met%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + POLYGON:      ',ipy
            write (unit=*,fmt='(a,1x,i6)')     ' + SITE:         ',isi
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',initp%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',initp%veg_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',initp%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',initp%veg_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',initp%hcapveg(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',cpatch%veg_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',cpatch%veg_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',cpatch%veg_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',cpatch%veg_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',cpatch%hcapveg(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch_ar(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do


      ksn   = csite%nlev_sfcwater(ipa)
      nsoil = csite%ntext_soil(nzg,ipa)
      nlsw1 = max(1, ksn)
      call ed_grndvap(ksn,nsoil,csite%soil_water(nzg,ipa),csite%soil_energy(nzg,ipa)       &
                     ,csite%sfcwater_energy(nlsw1,ipa),sngl(rk4met%rhos)                   &
                     ,csite%can_shv(ipa),csite%ground_shv(ipa),csite%surface_ssh(ipa)      &
                     ,surface_temp,surface_fliq)
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
