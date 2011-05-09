!==========================================================================================!
!==========================================================================================!
!     This module contains the wrappers for the Runge-Kutta integration scheme.            !
!==========================================================================================!
!==========================================================================================!
module rk4_driver

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      Main driver of short-time scale dynamics of the Runge-Kutta integrator           !
   !      for the land surface model.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_timestep(cgrid,ifm)
      use rk4_coms               , only : integration_vars     & ! structure
                                        , rk4patchtype         & ! structure
                                        , zero_rk4_patch       & ! subroutine
                                        , zero_rk4_cohort      & ! subroutine
                                        , integration_buff     & ! intent(out)
                                        , rk4site              ! ! intent(out)
      use ed_state_vars          , only : edtype               & ! structure
                                        , polygontype          & ! structure
                                        , sitetype             & ! structure
                                        , patchtype            ! ! structure
      use met_driver_coms        , only : met_driv_state       ! ! structure
      use grid_coms              , only : nzg                  & ! intent(in)
                                        , nzs                  ! ! intent(in)
      use canopy_struct_dynamics , only : canopy_turbulence8   ! ! subroutine
      use ed_misc_coms           , only : current_time         ! ! intent(in)
      implicit none

      !----------- Use MPI timing calls, need declarations --------------------------------!
      include 'mpif.h'
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target      :: cgrid
      integer                   , intent (in) :: ifm
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
      type(patchtype)           , pointer     :: cpatch
      type(met_driv_state)      , pointer     :: cmet
      integer                                 :: ipy
      integer                                 :: isi
      integer                                 :: ipa
      integer                                 :: iun
      integer, dimension(:)     , allocatable :: ed_ktrans
      integer                                 :: nsteps
      real                                    :: wcurr_loss2atm
      real                                    :: ecurr_loss2atm
      real                                    :: co2curr_loss2atm
      real                                    :: wcurr_loss2drainage
      real                                    :: ecurr_loss2drainage
      real                                    :: wcurr_loss2runoff
      real                                    :: ecurr_loss2runoff
      real                                    :: old_can_theiv
      real                                    :: old_can_shv
      real                                    :: old_can_co2
      real                                    :: old_can_rhos
      real                                    :: old_can_temp
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!


      !----- Allocate the auxiliary variables. --------------------------------------------!
      allocate(ed_ktrans(nzg))

      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               
               !----- Reset all buffers to zero, as a safety measure. ---------------------!
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

               !----- Get velocity for aerodynamic resistance. ----------------------------!
               if (csite%can_theta(ipa) < cmet%atm_theta) then
                  cmet%vels = cmet%vels_stab
               else
                  cmet%vels = cmet%vels_unstab
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Update roughness and canopy depth.                                     !
               !---------------------------------------------------------------------------!
               call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs                        &
                                             ,cpoly%ntext_soil(:,isi))
               call update_patch_derived_props(csite,cpoly%lsl(isi),cmet%prss,ipa)
               !---------------------------------------------------------------------------!


               !----- Save the previous thermodynamic state. ------------------------------!
               old_can_theiv    = csite%can_theiv(ipa)
               old_can_shv      = csite%can_shv(ipa)
               old_can_co2      = csite%can_co2(ipa)
               old_can_rhos     = csite%can_rhos(ipa)
               old_can_temp     = csite%can_temp(ipa)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Copy the meteorological variables to the rk4site structure.            !
               !---------------------------------------------------------------------------!
               call copy_met_2_rk4site(nzg,cmet%vels,cmet%atm_theiv,cmet%atm_theta         &
                                      ,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2,cmet%geoht   &
                                      ,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg          &
                                      ,cmet%prss,cmet%rshort,cmet%rlong,cmet%geoht         &
                                      ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)              &
                                      ,cpoly%green_leaf_factor(:,isi)                      &
                                      ,cgrid%lon(ipy),cgrid%lat(ipy))

               !----- Compute current storage terms. --------------------------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)

               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init(csite,ipa,integration_buff%initp)



               !----- Get photosynthesis, stomatal conductance, and transpiration. --------!
               call canopy_photosynthesis(csite,cmet,nzg,ipa,ed_ktrans,cpoly%lsl(isi)      &
                                         ,cpoly%ntext_soil(:,isi)                          &
                                         ,cpoly%leaf_aging_factor(:,isi)                   &
                                         ,cpoly%green_leaf_factor(:,isi))

               !----- Compute root and heterotrophic respiration. -------------------------!
               call soil_respiration(csite,ipa,nzg,cpoly%ntext_soil(:,isi))

               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init_carbon(csite,ipa,integration_buff%initp)

               !---------------------------------------------------------------------------!
               !    This is the driver for the integration process...                      !
               !---------------------------------------------------------------------------!
               call integrate_patch_rk4(csite,integration_buff%initp,ipa,wcurr_loss2atm    &
                                       ,ecurr_loss2atm,co2curr_loss2atm                    &
                                       ,wcurr_loss2drainage,ecurr_loss2drainage            &
                                       ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps)

               !----- Add the number of steps into the step counter. ----------------------!
               cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)

               !---------------------------------------------------------------------------!
               !    Update the minimum monthly temperature, based on canopy temperature.   !
               !---------------------------------------------------------------------------!
               if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                  cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
               end if
               
               !---------------------------------------------------------------------------!
               !     Compute the residuals.                                                !
               !---------------------------------------------------------------------------!
               call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa           &
                                  ,wcurr_loss2atm,ecurr_loss2atm,co2curr_loss2atm          &
                                  ,wcurr_loss2drainage,ecurr_loss2drainage                 &
                                  ,wcurr_loss2runoff,ecurr_loss2runoff,cpoly%area(isi)     &
                                  ,cgrid%cbudget_nep(ipy),old_can_theiv,old_can_shv        &
                                  ,old_can_co2,old_can_rhos,old_can_temp)

            end do patchloop
         end do siteloop

      end do polygonloop

      !----- De-allocate scratch variables. -----------------------------------------------!
      deallocate (ed_ktrans)

      return
   end subroutine rk4_timestep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_rk4(csite,initp,ipa,wcurr_loss2atm,ecurr_loss2atm            &
                                 ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage &
                                 ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use ed_misc_coms    , only : dtlsm                ! ! intent(in)
      use soil_coms       , only : soil_rough           & ! intent(in)
                                 , snow_rough           ! ! intent(in)
      use canopy_air_coms , only : exar8                ! ! intent(in)
      use consts_coms     , only : vonk8                & ! intent(in)
                                 , cp8                  & ! intent(in)
                                 , cpi8                 ! ! intent(in)
      use rk4_coms        , only : integration_vars     & ! structure
                                 , rk4patchtype         & ! structure
                                 , rk4site              & ! intent(inout)
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4                & ! intent(inout)
                                 , dtrk4i               ! ! intent(inout)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      type(rk4patchtype)    , target      :: initp
      integer               , intent(in)  :: ipa
      real                  , intent(out) :: wcurr_loss2atm
      real                  , intent(out) :: ecurr_loss2atm
      real                  , intent(out) :: co2curr_loss2atm
      real                  , intent(out) :: wcurr_loss2drainage
      real                  , intent(out) :: ecurr_loss2drainage
      real                  , intent(out) :: wcurr_loss2runoff
      real                  , intent(out) :: ecurr_loss2runoff
      integer               , intent(out) :: nsteps
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                          :: hbeg
      !----- Locally saved variable -------------------------------------------------------!
      logical                  , save       :: first_time=.true.
      !------------------------------------------------------------------------------------!

      !----- Assign some constants which will remain the same throughout the run. ---------!
      if (first_time) then
         first_time = .false.
         tbeg   = 0.d0
         tend   = dble(dtlsm)
         dtrk4  = tend - tbeg
         dtrk4i = 1.d0/dtrk4
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
      initp%cpwp = 0.d0
      initp%wpwp = 0.d0

      !----- Go into the ODE integrator. --------------------------------------------------!

      call odeint(hbeg,csite,ipa,nsteps)

      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
      initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
      initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
      initp%cpwp = initp%can_rhos * initp%cpwp * dtrk4i
      initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i
      
      
      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp(tend-tbeg,initp,csite,ipa,wcurr_loss2atm,ecurr_loss2atm            &
                       ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage           &
                       ,wcurr_loss2runoff,ecurr_loss2runoff)

      return
   end subroutine integrate_patch_rk4
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp(hdid,initp,csite,ipa,wbudget_loss2atm,ebudget_loss2atm          &
                          ,co2budget_loss2atm,wbudget_loss2drainage,ebudget_loss2drainage  &
                          ,wbudget_loss2runoff,ebudget_loss2runoff)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4site              & ! intent(in)
                                      , rk4min_veg_temp      & ! intent(in)
                                      , rk4max_veg_temp      & ! intent(in)
                                      , tiny_offset          & ! intent(in) 
                                      , checkbudget          ! ! intent(in)
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            & ! structure
                                      , edgrid_g             ! ! structure
      use consts_coms          , only : day_sec              & ! intent(in)
                                      , t3ple8               ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     & ! intent(in)
                                      , dtlsm                ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , slz8                 ! ! intent(in)
      use grid_coms            , only : nzg                  & ! intent(in)
                                      , nzs                  ! ! intent(in)
      use therm_lib            , only : qwtk                 & ! subroutine
                                      , rslif                ! ! function
      use canopy_air_coms      , only : i_blyr_condct        ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target      :: initp
      type(sitetype)    , target      :: csite
      real(kind=8)      , intent(in)  :: hdid
      integer           , intent(in)  :: ipa
      real              , intent(out) :: wbudget_loss2atm
      real              , intent(out) :: ebudget_loss2atm
      real              , intent(out) :: co2budget_loss2atm
      real              , intent(out) :: wbudget_loss2drainage
      real              , intent(out) :: ebudget_loss2drainage
      real              , intent(out) :: wbudget_loss2runoff
      real              , intent(out) :: ebudget_loss2runoff
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      integer                         :: mould
      integer                         :: ico
      integer                         :: k
      integer                         :: kclosest
      integer                         :: ksn
      integer                         :: nsoil
      integer                         :: nlsw1
      real(kind=8)                    :: tmp_energy
      real(kind=8)                    :: available_water
      !----- Local contants ---------------------------------------------------------------!
      real        , parameter         :: tendays_sec=10.*day_sec
      !----- External function ------------------------------------------------------------!
      real        , external          :: sngloff
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.  All floating point variables are converted back  !
      ! to single precision.                                                               !
      !------------------------------------------------------------------------------------!
      csite%can_theiv(ipa)        = sngloff(initp%can_theiv       ,tiny_offset)
      csite%can_theta(ipa)        = sngloff(initp%can_theta       ,tiny_offset)
      csite%can_prss(ipa)         = sngloff(initp%can_prss        ,tiny_offset)
      csite%can_temp(ipa)         = sngloff(initp%can_temp        ,tiny_offset)
      csite%can_shv(ipa)          = sngloff(initp%can_shv         ,tiny_offset)
      csite%can_co2(ipa)          = sngloff(initp%can_co2         ,tiny_offset)
      csite%can_rhos(ipa)         = sngloff(initp%can_rhos        ,tiny_offset)
      csite%can_depth(ipa)        = sngloff(initp%can_depth       ,tiny_offset)
      csite%veg_displace(ipa)     = sngloff(initp%veg_displace    ,tiny_offset)
      csite%rough(ipa)            = sngloff(initp%rough           ,tiny_offset)
      csite%snowfac(ipa)          = sngloff(initp%snowfac         ,tiny_offset)
      csite%total_sfcw_depth(ipa) = sngloff(initp%total_sfcw_depth,tiny_offset)

      csite%ggbare(ipa)           = sngloff(initp%ggbare          ,tiny_offset)
      csite%ggveg (ipa)           = sngloff(initp%ggveg           ,tiny_offset)
      csite%ggnet (ipa)           = sngloff(initp%ggnet           ,tiny_offset)

      csite%ustar (ipa)           = sngloff(initp%ustar           ,tiny_offset)
      csite%tstar (ipa)           = sngloff(initp%tstar           ,tiny_offset)
      csite%qstar (ipa)           = sngloff(initp%qstar           ,tiny_offset)
      csite%cstar (ipa)           = sngloff(initp%cstar           ,tiny_offset)

      csite%zeta  (ipa)           = sngloff(initp%zeta            ,tiny_offset)
      csite%ribulk(ipa)           = sngloff(initp%ribulk          ,tiny_offset)

      csite%upwp  (ipa)           = sngloff(initp%upwp            ,tiny_offset)
      csite%wpwp  (ipa)           = sngloff(initp%wpwp            ,tiny_offset)
      csite%tpwp  (ipa)           = sngloff(initp%tpwp            ,tiny_offset)
      csite%qpwp  (ipa)           = sngloff(initp%qpwp            ,tiny_offset)
      csite%cpwp  (ipa)           = sngloff(initp%cpwp            ,tiny_offset)

      !------------------------------------------------------------------------------------!
      !    These variables are fast scale fluxes, and they may not be allocated, so just   !
      ! check this before copying.                                                         !
      !------------------------------------------------------------------------------------!
      if (fast_diagnostics) then

         csite%avg_vapor_vc(ipa)         =sngloff(initp%avg_vapor_vc      ,tiny_offset)
         csite%avg_dew_cg(ipa)           =sngloff(initp%avg_dew_cg        ,tiny_offset)
         csite%avg_vapor_gc(ipa)         =sngloff(initp%avg_vapor_gc      ,tiny_offset)
         csite%avg_wshed_vg(ipa)         =sngloff(initp%avg_wshed_vg      ,tiny_offset)
         csite%avg_intercepted(ipa)      =sngloff(initp%avg_intercepted   ,tiny_offset)
         csite%avg_throughfall(ipa)      =sngloff(initp%avg_throughfall   ,tiny_offset)
         csite%avg_vapor_ac(ipa)         =sngloff(initp%avg_vapor_ac      ,tiny_offset)
         csite%avg_transp(ipa)           =sngloff(initp%avg_transp        ,tiny_offset)
         csite%avg_evap(ipa)             =sngloff(initp%avg_evap          ,tiny_offset)
         csite%avg_drainage(ipa)         =sngloff(initp%avg_drainage      ,tiny_offset)
         csite%avg_drainage_heat(ipa)    =sngloff(initp%avg_drainage_heat ,tiny_offset)
         csite%avg_netrad(ipa)           =sngloff(initp%avg_netrad        ,tiny_offset)
         csite%avg_sensible_vc(ipa)      =sngloff(initp%avg_sensible_vc   ,tiny_offset)
         csite%avg_qwshed_vg(ipa)        =sngloff(initp%avg_qwshed_vg     ,tiny_offset)
         csite%avg_qintercepted(ipa)     =sngloff(initp%avg_qintercepted  ,tiny_offset)
         csite%avg_qthroughfall(ipa)     =sngloff(initp%avg_qthroughfall  ,tiny_offset)
         csite%avg_sensible_gc(ipa)      =sngloff(initp%avg_sensible_gc   ,tiny_offset)
         csite%avg_sensible_ac(ipa)      =sngloff(initp%avg_sensible_ac   ,tiny_offset)
         csite%avg_carbon_ac(ipa)        =sngloff(initp%avg_carbon_ac     ,tiny_offset)
         do k = rk4site%lsl, nzg
            csite%avg_sensible_gg(k,ipa) =sngloff(initp%avg_sensible_gg(k)   ,tiny_offset)
            csite%avg_smoist_gg(k,ipa)   =sngloff(initp%avg_smoist_gg(k)     ,tiny_offset)
            csite%avg_transloss(k,ipa)   =sngloff(initp%avg_transloss(k)     ,tiny_offset)
         end do
      end if

      if(checkbudget) then
         co2budget_loss2atm    = sngloff(initp%co2budget_loss2atm   ,tiny_offset)
         ebudget_loss2atm      = sngloff(initp%ebudget_loss2atm     ,tiny_offset)
         ebudget_loss2drainage = sngloff(initp%ebudget_loss2drainage,tiny_offset)
         ebudget_loss2runoff   = sngloff(initp%ebudget_loss2runoff  ,tiny_offset)
         wbudget_loss2atm      = sngloff(initp%wbudget_loss2atm     ,tiny_offset)
         wbudget_loss2drainage = sngloff(initp%wbudget_loss2drainage,tiny_offset)
         wbudget_loss2runoff   = sngloff(initp%wbudget_loss2runoff  ,tiny_offset)
      else
         co2budget_loss2atm             = 0.
         ebudget_loss2atm               = 0.
         ebudget_loss2drainage          = 0.
         ebudget_loss2runoff            = 0.
         wbudget_loss2atm               = 0.
         wbudget_loss2drainage          = 0.
         wbudget_loss2runoff            = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The following is not a pure diagnostic, it is used for phenology and mortality !
      ! functions, preserve this variable and its dependencies in all contexts.            !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! paw_avg - 10-day average of plant available water.                                 !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         available_water = 0.d0
         do k = cpatch%krdepth(ico), nzg - 1
            nsoil = rk4site%ntext_soil(k)
            available_water = available_water                                              &
                            + max(0.d0,(initp%soil_water(k) - soil8(nsoil)%soilwp))        &
                            * (slz8(k+1)-slz8(k))                                          &
                            / (soil8(nsoil)%slmsts - soil8(nsoil)%soilwp)
         end do
         nsoil = rk4site%ntext_soil(nzg)
         available_water = available_water                                                 &
                         + max(0.d0,(initp%soil_water(nzg) - soil8(nsoil)%soilwp))         &
                         * (-1.d0*slz8(nzg))                                               &
                         / (soil8(nsoil)%slmsts -soil8(nsoil)%soilwp) 
         available_water = available_water / (-1.d0*slz8(cpatch%krdepth(ico)))


         cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)            &
                             + sngl(available_water)*sngl(hdid)/tendays_sec
      end do

      
      do k = rk4site%lsl, nzg
         csite%soil_water(k,ipa)   = sngloff(initp%soil_water(k)  ,tiny_offset)
         csite%soil_energy(k,ipa)  = sngloff(initp%soil_energy(k) ,tiny_offset)
         csite%soil_tempk(k,ipa)   = sngloff(initp%soil_tempk(k)  ,tiny_offset)
         csite%soil_fracliq(k,ipa) = sngloff(initp%soil_fracliq(k),tiny_offset)
      end do
      

      !------------------------------------------------------------------------------------!
      !    Surface water energy is computed in J/m² inside the integrator. Convert it back !
      ! to J/kg in the layers that surface water/snow still exists.                        !
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

         if (initp%resolvable(ico)) then
            select case (i_blyr_condct)
            case (-1)
               !---------------------------------------------------------------------------!
               !    The cohort was solved, update internal energy and water, and re-       !
               ! calculate temperature.  Note that energy may need to be scaled back.      !
               !---------------------------------------------------------------------------!
               cpatch%veg_water(ico)  = sngloff(initp%veg_water(ico),tiny_offset)
               tmp_energy             = initp%veg_energy(ico)                              &
                                      + (dble(cpatch%hcapveg(ico))-initp%hcapveg(ico))     &
                                      * initp%veg_temp(ico)

               cpatch%veg_energy(ico) = sngloff(tmp_energy,tiny_offset)
               call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)  &
                        ,cpatch%veg_temp(ico),cpatch%veg_fliq(ico))

            case default
               !---------------------------------------------------------------------------!
               !    The cohort was solved, update water and internal energy, and re-       !
               ! calculate the temperature and leaf intercellular specific humidity.  The  !
               ! vegetation dry heat capacity is constant within one time step, so it does !
               ! not need to be updated.                                                   !
               !---------------------------------------------------------------------------!
               cpatch%veg_water(ico)  = sngloff(initp%veg_water(ico) , tiny_offset)
               cpatch%veg_energy(ico) = sngloff(initp%veg_energy(ico), tiny_offset)
               call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)  &
                        ,cpatch%veg_temp(ico),cpatch%veg_fliq(ico))

            end select

            !------------------------------------------------------------------------------!
            !     The intercellular specific humidity is always assumed to be at           !
            ! saturation for a given temperature.  Find the saturation mixing ratio, then  !
            ! convert it to specific humidity.                                             !
            !------------------------------------------------------------------------------!
            cpatch%lint_shv(ico) = rslif(csite%can_prss(ipa),cpatch%veg_temp(ico))
            cpatch%lint_shv(ico) = cpatch%lint_shv(ico) / (1. + cpatch%lint_shv(ico))
            !----- Convert the wind. ------------------------------------------------------!
            cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Copy the conductances.                                                   !
            !------------------------------------------------------------------------------!
            cpatch%gbh       (ico) = sngloff(initp%gbh       (ico), tiny_offset)
            cpatch%gbw       (ico) = sngloff(initp%gbw       (ico), tiny_offset)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Divide the values of water demand by the time step to obtain the average !
            ! value over the past DTLSM period.                                            !
            !------------------------------------------------------------------------------!
            cpatch%psi_open  (ico) = sngloff(initp%psi_open  (ico),tiny_offset) / hdid
            cpatch%psi_closed(ico) = sngloff(initp%psi_closed(ico),tiny_offset) / hdid

         elseif (cpatch%hite(ico) <=  csite%total_sfcw_depth(ipa)) then
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
            !------------------------------------------------------------------------------!
            !     The intercellular specific humidity is always assumed to be at           !
            ! saturation for a given temperature.  Find the saturation mixing ratio, then  !
            ! convert it to specific humidity.                                             !
            !------------------------------------------------------------------------------!
            cpatch%lint_shv(ico) = rslif(csite%can_prss(ipa),cpatch%veg_temp(ico))
            cpatch%lint_shv(ico) = cpatch%lint_shv(ico) / (1. + cpatch%lint_shv(ico))
            !----- Copy the meteorological wind to here. ----------------------------------!
            cpatch%veg_wind(ico) = sngloff(rk4site%vels, tiny_offset)
            !----- Make water demand 0. ---------------------------------------------------!
            cpatch%psi_open  (ico) = 0.0
            cpatch%psi_closed(ico) = 0.0

         else
            !------------------------------------------------------------------------------!
            !     For plants with minimal foliage or very sparse patches, fix the leaf     !
            ! temperature to the canopy air space and force veg_water to be zero.          !
            !------------------------------------------------------------------------------!
            cpatch%veg_temp(ico)   = csite%can_temp(ipa)
            cpatch%veg_fliq(ico)   = 0.
            cpatch%veg_water(ico)  = 0. 
            cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
            !------------------------------------------------------------------------------!
            !     The intercellular specific humidity is always assumed to be at           !
            ! saturation for a given temperature.  Find the saturation mixing ratio, then  !
            ! convert it to specific humidity.                                             !
            !------------------------------------------------------------------------------!
            cpatch%lint_shv(ico) = rslif(csite%can_prss(ipa),cpatch%veg_temp(ico))
            cpatch%lint_shv(ico) = cpatch%lint_shv(ico) / (1. + cpatch%lint_shv(ico))
            !----- Copy the meteorological wind to here. ----------------------------------!
            cpatch%veg_wind(ico) = sngloff(rk4site%vels, tiny_offset)
            !----- Make water demand 0. ---------------------------------------------------!
            cpatch%psi_open  (ico) = 0.0
            cpatch%psi_closed(ico) = 0.0
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
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4site%lat
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
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do

      !------ Copy the ground variables to the output. ------------------------------------!
      csite%ground_shv (ipa) = sngloff(initp%ground_shv , tiny_offset)
      csite%ground_ssh (ipa) = sngloff(initp%ground_ssh , tiny_offset)
      csite%ground_temp(ipa) = sngloff(initp%ground_temp, tiny_offset)
      csite%ground_fliq(ipa) = sngloff(initp%ground_fliq, tiny_offset)
      return
   end subroutine initp2modelp
   !=======================================================================================!
   !=======================================================================================! 
end module rk4_driver

!==========================================================================================!
!==========================================================================================!
