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
      use ed_misc_coms           , only : current_time         & ! intent(in)
                                        , dtlsm                ! ! intent(in)
      use therm_lib              , only : tq2enthalpy          ! ! function
      use budget_utils           , only : update_budget        & ! function
                                        , compute_budget       ! ! function
      !$ use omp_lib
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

      type(rk4patchtype)       , pointer      :: initp
      type(rk4patchtype)       , pointer      :: yscal
      type(rk4patchtype)       , pointer      :: y
      type(rk4patchtype)       , pointer      :: dydx
      type(rk4patchtype)       , pointer      :: yerr
      type(rk4patchtype)       , pointer      :: ytemp
      type(rk4patchtype)       , pointer      :: ak2
      type(rk4patchtype)       , pointer      :: ak3
      type(rk4patchtype)       , pointer      :: ak4
      type(rk4patchtype)       , pointer      :: ak5
      type(rk4patchtype)       , pointer      :: ak6
      type(rk4patchtype)       , pointer      :: ak7
      integer                                 :: ipy
      integer                                 :: isi
      integer                                 :: ipa
      integer                                 :: iun
      integer                                 :: nsteps
      integer                                 :: imon
      real                                    :: wcurr_loss2atm
      real                                    :: ecurr_netrad
      real                                    :: ecurr_loss2atm
      real                                    :: co2curr_loss2atm
      real                                    :: wcurr_loss2drainage
      real                                    :: ecurr_loss2drainage
      real                                    :: wcurr_loss2runoff
      real                                    :: ecurr_loss2runoff
      real                                    :: ecurr_prsseffect
      real                                    :: old_can_enthalpy
      real                                    :: old_can_shv
      real                                    :: old_can_co2
      real                                    :: old_can_rhos
      real                                    :: old_can_temp
      real                                    :: old_can_prss
      real                                    :: old_can_depth
      real                                    :: patch_vels
      integer                                 :: ibuff
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!

      


      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)

            !------------------------------------------------------------------------------!
            !     Update the monthly rainfall.                                             !
            !------------------------------------------------------------------------------!
            imon                             = current_time%month
            cpoly%avg_monthly_pcpg(imon,isi) = cpoly%avg_monthly_pcpg(imon,isi)            &
                                             + cmet%pcpg * dtlsm
            !------------------------------------------------------------------------------!

            !---------------------------------------------------------------------------!
            !    Copy the meteorological variables to the rk4site structure.            !
            !---------------------------------------------------------------------------!
            call copy_met_2_rk4site(nzg,cmet%atm_ustar,cmet%atm_theiv,cmet%atm_vpdef,   &
                  cmet%atm_theta,        &
                  cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2,cmet%geoht,   &
                  cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg,          &
                  cmet%prss,cmet%rshort,cmet%rlong,cmet%par_beam,      &
                  cmet%par_diffuse,cmet%nir_beam,cmet%nir_diffuse,     &
                  cmet%geoht,cpoly%lsl(isi),cpoly%ntext_soil(:,isi),   &
                  cpoly%green_leaf_factor(:,isi),cgrid%lon(ipy),       &
                  cgrid%lat(ipy),cgrid%cosz(ipy))
            !---------------------------------------------------------------------------!

            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(      &
            !$OMP initp,yscal,y,dydx,yerr,ytemp,ak2,ak3,    &
            !$OMP ak4,ak5,ak6,ak7,cpatch,patch_vels,        &
            !$OMP old_can_co2,old_can_rhos,old_can_temp,    &
            !$OMP old_can_prss,old_can_enthalpy,            &
            !$OMP old_can_shv,ecurr_netrad,                 &
            !$OMP wcurr_loss2atm,ecurr_loss2atm,            &
            !$OMP co2curr_loss2atm,                         &
            !$OMP wcurr_loss2drainage,ecurr_loss2drainage,  &
            !$OMP wcurr_loss2runoff,ecurr_loss2runoff,nsteps )

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               ibuff = 1
               !$ ibuff = OMP_get_thread_num()+1

               initp => integration_buff(ibuff)%initp
               yscal => integration_buff(ibuff)%yscal
               y     => integration_buff(ibuff)%y
               dydx  => integration_buff(ibuff)%dydx
               yerr  => integration_buff(ibuff)%yerr
               ytemp => integration_buff(ibuff)%ytemp
               ak2   => integration_buff(ibuff)%ak2
               ak3   => integration_buff(ibuff)%ak3
               ak4   => integration_buff(ibuff)%ak4
               ak5   => integration_buff(ibuff)%ak5
               ak6   => integration_buff(ibuff)%ak6
               ak7   => integration_buff(ibuff)%ak7



               !----- Reset all buffers to zero, as a safety measure. ---------------------!
               call zero_rk4_patch(initp)
               call zero_rk4_patch(yscal)
               call zero_rk4_patch(y)
               call zero_rk4_patch(dydx)
               call zero_rk4_patch(yerr)
               call zero_rk4_patch(ytemp)
               call zero_rk4_patch(ak2)
               call zero_rk4_patch(ak3)
               call zero_rk4_patch(ak4)
               call zero_rk4_patch(ak5)
               call zero_rk4_patch(ak6)
               call zero_rk4_patch(ak7)
               call zero_rk4_cohort(initp)
               call zero_rk4_cohort(yscal)
               call zero_rk4_cohort(y)
               call zero_rk4_cohort(dydx)
               call zero_rk4_cohort(yerr)
               call zero_rk4_cohort(ytemp)
               call zero_rk4_cohort(ak2)
               call zero_rk4_cohort(ak3)
               call zero_rk4_cohort(ak4)
               call zero_rk4_cohort(ak5)
               call zero_rk4_cohort(ak6)
               call zero_rk4_cohort(ak7)

               !----- Get velocity for aerodynamic resistance. ----------------------------!
               if (csite%can_theta(ipa) < cmet%atm_theta) then
                  patch_vels = cmet%vels_stab
               else
                  patch_vels = cmet%vels_unstab
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
               old_can_shv      = csite%can_shv  (ipa)
               old_can_co2      = csite%can_co2  (ipa)
               old_can_rhos     = csite%can_rhos (ipa)
               old_can_temp     = csite%can_temp (ipa)
               old_can_prss     = csite%can_prss (ipa)
               old_can_enthalpy = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
               !---------------------------------------------------------------------------!

               !----- Compute current storage terms. --------------------------------------!
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init(csite,ipa,initp,patch_vels)
               !---------------------------------------------------------------------------!


               !----- Get photosynthesis, stomatal conductance, and transpiration. --------!
               call canopy_photosynthesis(csite,cmet,nzg,ipa,cpoly%lsl(isi)                &
                                         ,cpoly%ntext_soil(:,isi)                          &
                                         ,cpoly%leaf_aging_factor(:,isi)                   &
                                         ,cpoly%green_leaf_factor(:,isi))
               !---------------------------------------------------------------------------!


               !----- Compute root and heterotrophic respiration. -------------------------!
               call soil_respiration(csite,ipa,nzg,cpoly%ntext_soil(:,isi))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Set up the integration patch.                                         !
               !---------------------------------------------------------------------------!
               call copy_patch_init_carbon(csite,ipa,initp)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    This is the driver for the integration process...                      !
               !---------------------------------------------------------------------------!
               call integrate_patch_rk4(csite,initp,ipa                                    &
                                       ,cpoly%nighttime(isi),wcurr_loss2atm                &
                                       ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm       &
                                       ,wcurr_loss2drainage,ecurr_loss2drainage            &
                                       ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps)
               !---------------------------------------------------------------------------!


               !----- Add the number of steps into the step counter. ----------------------!
               !----- workload accumulation is order-independent, so this can stay shared 
               cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Update the minimum monthly temperature, based on canopy temperature.   !
               !---------------------------------------------------------------------------!
               if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                  cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Compute the residuals.                                                !
               !---------------------------------------------------------------------------!
               call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa           &
                                  ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm              &
                                  ,co2curr_loss2atm,wcurr_loss2drainage                    &
                                  ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff &
                                  ,cpoly%area(isi),cgrid%cbudget_nep(ipy),old_can_enthalpy &
                                  ,old_can_shv,old_can_co2,old_can_rhos,old_can_temp       &
                                  ,old_can_prss)
               !---------------------------------------------------------------------------!
            end do patchloop
            !$OMP END PARALLEL DO

            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polygonloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine rk4_timestep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_rk4(csite,initp,ipa,nighttime,wcurr_loss2atm,ecurr_netrad    &
                                 ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage      &
                                 ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff  &
                                 ,nsteps)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use ed_misc_coms    , only : dtlsm                ! ! intent(in)
      use soil_coms       , only : soil_rough           & ! intent(in)
                                 , snow_rough           ! ! intent(in)
      use canopy_air_coms , only : exar8                ! ! intent(in)
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
      logical               , intent(in)  :: nighttime
      real                  , intent(out) :: wcurr_loss2atm
      real                  , intent(out) :: ecurr_netrad
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
      call initp2modelp(tend-tbeg,initp,csite,ipa,nighttime,wcurr_loss2atm,ecurr_netrad    &
                       ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage                &
                       ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff)

      return
   end subroutine integrate_patch_rk4
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will copy the variables from the integration buffer to the state  !
   ! patch and cohorts.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine initp2modelp(hdid,initp,csite,ipa,nighttime,wbudget_loss2atm,ebudget_netrad  &
                          ,ebudget_loss2atm,co2budget_loss2atm,wbudget_loss2drainage       &
                          ,ebudget_loss2drainage,wbudget_loss2runoff,ebudget_loss2runoff)
      use rk4_coms             , only : rk4patchtype         & ! structure
                                      , rk4site              & ! intent(in)
                                      , rk4min_veg_temp      & ! intent(in)
                                      , rk4max_veg_temp      & ! intent(in)
                                      , tiny_offset          & ! intent(in) 
                                      , checkbudget          & ! intent(in)
                                      , ibranch_thermo       ! ! intent(in)
      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            & ! structure
                                      , edgrid_g             ! ! structure
      use consts_coms          , only : day_sec              & ! intent(in)
                                      , t3ple                & ! intent(in)
                                      , t3ple8               & ! intent(in)
                                      , wdns8                ! ! intent(in)
      use ed_misc_coms         , only : fast_diagnostics     & ! intent(in)
                                      , writing_long         & ! intent(in)
                                      , dtlsm                & ! intent(in)
                                      , frqsum               ! ! intent(in)
      use soil_coms            , only : soil8                & ! intent(in)
                                      , dslz8                & ! intent(in)
                                      , slz8                 & ! intent(in)
                                      , slzt8                ! ! intent(in)
      use grid_coms            , only : nzg                  & ! intent(in)
                                      , nzs                  ! ! intent(in)
      use therm_lib            , only : thetaeiv             & ! subroutine
                                      , vpdefil              & ! subroutine
                                      , uextcm2tl            & ! subroutine
                                      , cmtl2uext            & ! subroutine
                                      , qslif                ! ! function
      use phenology_coms       , only : spot_phen            ! ! intent(in)
      use allometry            , only : h2crownbh            ! ! function
      use disturb_coms         , only : include_fire         & ! intent(in)
                                      , k_fire_first         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target      :: initp
      type(sitetype)    , target      :: csite
      real(kind=8)      , intent(in)  :: hdid
      integer           , intent(in)  :: ipa
      logical           , intent(in)  :: nighttime
      real              , intent(out) :: wbudget_loss2atm
      real              , intent(out) :: ebudget_netrad
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
      integer                         :: ipft
      integer                         :: k
      integer                         :: ka
      integer                         :: kroot
      integer                         :: kclosest
      integer                         :: ksn
      integer                         :: nsoil
      integer                         :: nlsw1
      real(kind=8)                    :: tmp_energy
      real(kind=8)                    :: available_water
      real(kind=8)                    :: gnd_water
      real(kind=8)                    :: psiplusz
      real(kind=8)                    :: mcheight
      real(kind=4)                    :: step_waterdef
      real(kind=4)                    :: can_rvap
      !----- Local contants ---------------------------------------------------------------!
      real        , parameter         :: tendays_sec    = 10. * day_sec
      real        , parameter         :: thirtydays_sec = 30. * day_sec
      !----- External function ------------------------------------------------------------!
      real        , external          :: sngloff
      !----- Locally saved variables. -----------------------------------------------------!
      real        , save              :: dtlsm_o_frqsum
      logical     , save              :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- Assign the constant scaling factor. ------------------------------------------!
      if (first_time) then
         first_time     = .false.
         dtlsm_o_frqsum = dtlsm / frqsum
      end if
      !------------------------------------------------------------------------------------!


      !----- Alias for the cohorts. -------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Most variables require just a simple copy.  More comments will be made next to !
      ! those in which this is not true.  All floating point variables are converted back  !
      ! to single precision.                                                               !
      !------------------------------------------------------------------------------------!
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

      !------------------------------------------------------------------------------------!
      !    Find the ice-vapour equivalent potential temperature.  This is done outside the !
      ! integrator because it is an iterative method and currently we are not using it as  !
      ! a prognostic variable.                                                             !
      !------------------------------------------------------------------------------------!
      can_rvap                    = csite%can_shv(ipa) / ( 1.0 - csite%can_shv(ipa))
      csite%can_theiv(ipa)        = thetaeiv(csite%can_theta (ipa), csite%can_prss(ipa)    &
                                            ,csite%can_temp  (ipa), can_rvap               &
                                            ,can_rvap             )
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Find the vapour pressure deficit, which is diagnostic only.                     !
      !------------------------------------------------------------------------------------!
      csite%can_vpdef(ipa)        = vpdefil(csite%can_prss(ipa),csite%can_temp(ipa)        &
                                           ,csite%can_shv(ipa) ,.true.)
      !------------------------------------------------------------------------------------!



      !------ Copy the ground variables to the output. ------------------------------------!
      csite%ground_shv (ipa) = sngloff(initp%ground_shv , tiny_offset)
      csite%ground_ssh (ipa) = sngloff(initp%ground_ssh , tiny_offset)
      csite%ground_temp(ipa) = sngloff(initp%ground_temp, tiny_offset)
      csite%ground_fliq(ipa) = sngloff(initp%ground_fliq, tiny_offset)
      !------------------------------------------------------------------------------------!



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
         csite%fmean_vapor_gc        (ipa) = sngloff(initp%avg_vapor_gc       ,tiny_offset)
         csite%fmean_throughfall     (ipa) = sngloff(initp%avg_throughfall    ,tiny_offset)
         csite%fmean_vapor_ac        (ipa) = sngloff(initp%avg_vapor_ac       ,tiny_offset)
         csite%fmean_drainage        (ipa) = sngloff(initp%avg_drainage       ,tiny_offset)
         csite%fmean_qdrainage       (ipa) = sngloff(initp%avg_qdrainage      ,tiny_offset)
         csite%fmean_qthroughfall    (ipa) = sngloff(initp%avg_qthroughfall   ,tiny_offset)
         csite%fmean_sensible_gc     (ipa) = sngloff(initp%avg_sensible_gc    ,tiny_offset)
         csite%fmean_sensible_ac     (ipa) = sngloff(initp%avg_sensible_ac    ,tiny_offset)
         csite%fmean_carbon_ac       (ipa) = sngloff(initp%avg_carbon_ac      ,tiny_offset)
         csite%fmean_carbon_st       (ipa) = sngloff(initp%avg_carbon_st      ,tiny_offset)
         csite%fmean_ustar           (ipa) = sngloff(initp%avg_ustar          ,tiny_offset)
         csite%fmean_tstar           (ipa) = sngloff(initp%avg_tstar          ,tiny_offset)
         csite%fmean_qstar           (ipa) = sngloff(initp%avg_qstar          ,tiny_offset)
         csite%fmean_cstar           (ipa) = sngloff(initp%avg_cstar          ,tiny_offset)
         do k = rk4site%lsl, nzg
            csite%fmean_sensible_gg(k,ipa) = sngloff(initp%avg_sensible_gg(k) ,tiny_offset)
            csite%fmean_smoist_gg  (k,ipa) = sngloff(initp%avg_smoist_gg  (k) ,tiny_offset)
            csite%fmean_transloss  (k,ipa) = sngloff(initp%avg_transloss  (k) ,tiny_offset)
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Cohort-level variables.                                                     !
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts
            cpatch%fmean_sensible_lc   (ico) = sngloff(initp%avg_sensible_lc    (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_sensible_wc   (ico) = sngloff(initp%avg_sensible_wc    (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_vapor_lc      (ico) = sngloff(initp%avg_vapor_lc       (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_vapor_wc      (ico) = sngloff(initp%avg_vapor_wc       (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_transp        (ico) = sngloff(initp%avg_transp         (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_intercepted_al(ico) = sngloff(initp%avg_intercepted_al (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_intercepted_aw(ico) = sngloff(initp%avg_intercepted_aw (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_wshed_lg      (ico) = sngloff(initp%avg_wshed_lg       (ico)      &
                                                     ,tiny_offset)
            cpatch%fmean_wshed_wg      (ico) = sngloff(initp%avg_wshed_wg       (ico)      &
                                                     ,tiny_offset)
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      if(checkbudget) then
         co2budget_loss2atm    = sngloff(initp%co2budget_loss2atm   ,tiny_offset)
         ebudget_netrad        = sngloff(initp%ebudget_netrad       ,tiny_offset)
         ebudget_loss2atm      = sngloff(initp%ebudget_loss2atm     ,tiny_offset)
         ebudget_loss2drainage = sngloff(initp%ebudget_loss2drainage,tiny_offset)
         ebudget_loss2runoff   = sngloff(initp%ebudget_loss2runoff  ,tiny_offset)
         wbudget_loss2atm      = sngloff(initp%wbudget_loss2atm     ,tiny_offset)
         wbudget_loss2drainage = sngloff(initp%wbudget_loss2drainage,tiny_offset)
         wbudget_loss2runoff   = sngloff(initp%wbudget_loss2runoff  ,tiny_offset)
      else
         co2budget_loss2atm             = 0.
         ebudget_netrad                 = 0.
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
      !     Update the water deficit.  This is done as a 30-day running average.           !
      !------------------------------------------------------------------------------------!
      step_waterdef                   = sngloff(initp%water_deficit,tiny_offset)
      csite%avg_monthly_waterdef(ipa) = csite%avg_monthly_waterdef(ipa) + step_waterdef
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     This variable is the monthly mean ground water that will be used to control    !
      ! fire disturbance.                                                                  !
      !------------------------------------------------------------------------------------!
      gnd_water = 0.d0
      !----- Add temporary surface water. -------------------------------------------------!
      do k=1,initp%nlev_sfcwater
         gnd_water = gnd_water + initp%sfcwater_mass(k)
      end do
      !----- Find the bottommost layer to consider. ---------------------------------------!
      select case(include_fire)
      case (1)
         ka = rk4site%lsl
      case default
         ka = k_fire_first
      end select
      !----- Add soil moisture. -----------------------------------------------------------!
      do k=ka,nzg
         gnd_water = gnd_water + initp%soil_water(k) * dslz8(k) * wdns8
      end do
      !----- Add to the monthly mean. -----------------------------------------------------!
      csite%avg_monthly_gndwater(ipa) = csite%avg_monthly_gndwater(ipa)                    &
                                      + sngloff(gnd_water,tiny_offset)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! paw_avg - 10-day average of relative plant available water.  The relative value    !
      !           depends on whether the user wants to define phenology based on soil      !
      !           moisture or soil potential.                                              !
      !------------------------------------------------------------------------------------!
      if (spot_phen) then
         do ico = 1,cpatch%ncohorts
            ipft  = cpatch%pft(ico)
            kroot = cpatch%krdepth(ico)

            available_water = 0.d0
            do k = kroot, nzg
               nsoil            = rk4site%ntext_soil(k)
               mcheight         = 5.d-1 * ( dble(cpatch%hite(ico))                         &
                                          + dble(h2crownbh(cpatch%hite(ico),ipft)) )
               psiplusz         = slzt8(k) - mcheight                                      &
                                + soil8(nsoil)%slpots                                      &
                                / (initp%soil_water(k) / soil8(nsoil)%slmsts)              &
                                ** soil8(nsoil)%slbs
               available_water  = available_water                                          &
                                + max(0.d0,(psiplusz - soil8(nsoil)%slpotwp)) * dslz8(k)   &
                                / (soil8(nsoil)%slpotld - soil8(nsoil)%slpotwp)
            end do
            available_water     = available_water / abs(slz8(kroot))
            cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)         &
                                + sngl(available_water)*sngl(hdid)/tendays_sec
         end do
      else
         do ico = 1,cpatch%ncohorts
            available_water = 0.d0
            kroot           = cpatch%krdepth(ico)
            do k = kroot, nzg
               nsoil            = rk4site%ntext_soil(k)
               available_water  = available_water                                          &
                                + max(0.d0,(initp%soil_water(k)   - soil8(nsoil)%soilwp))  &
                                * dslz8(k) / (soil8(nsoil)%soilld - soil8(nsoil)%soilwp)
            end do
            available_water     = available_water / abs(slz8(kroot))
            cpatch%paw_avg(ico) = cpatch%paw_avg(ico)*(1.0-sngl(hdid)/tendays_sec)         &
                                + sngl(available_water)*sngl(hdid)/tendays_sec
         end do
      end if
      
      do k = rk4site%lsl, nzg
         csite%soil_water  (k,ipa) = sngloff(initp%soil_water  (k),tiny_offset)
         csite%soil_mstpot (k,ipa) = sngloff(initp%soil_mstpot (k),tiny_offset)
         csite%soil_energy (k,ipa) = sngloff(initp%soil_energy (k),tiny_offset)
         csite%soil_tempk  (k,ipa) = sngloff(initp%soil_tempk  (k),tiny_offset)
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
            csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk  (nzg,ipa)
         else
            csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
            csite%sfcwater_tempk  (k,ipa) = csite%sfcwater_tempk  (k-1,ipa)
         end if
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Cohort variables.  Here we must check whether the cohort was really solved or  !
      ! it was skipped after being flagged as "unsafe".  In case the cohort was skipped,   !
      ! we must check whether it was because it was too small or because it was buried in  !
      ! snow.                                                                              !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         select case (ibranch_thermo)
         case (1)
            !------------------------------------------------------------------------------!
            !  VEGETATION -- Leaf and branchwood were solved together, so they must remain !
            !                in thermal equilibrium.                                       !
            !------------------------------------------------------------------------------!
            if (initp%veg_resolvable(ico)) then

               !---------------------------------------------------------------------------!
               !     Copy vegetation wind.                                                 !
               !---------------------------------------------------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    LEAVES.  It is always safe to copy internal energy and standing water, !
               !             but we must check whether leaves were truly resolved or not   !
               !             before copying the other variables.                           !
               !---------------------------------------------------------------------------!
               cpatch%leaf_water (ico) = sngloff(initp%leaf_water (ico) , tiny_offset)
               cpatch%leaf_energy(ico) = sngloff(initp%leaf_energy(ico) , tiny_offset)


               if (initp%leaf_resolvable(ico)) then
                  !------------------------------------------------------------------------!
                  !    Leaves were solved, find the temperature and liquid fraction from   !
                  ! internal energy.                                                       !
                  !------------------------------------------------------------------------!
                  call uextcm2tl(cpatch%leaf_energy(ico),cpatch%leaf_water(ico)            &
                                ,cpatch%leaf_hcap(ico),cpatch%leaf_temp(ico)               &
                                ,cpatch%leaf_fliq(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The intercellular specific humidity is always assumed to be at     !
                  ! saturation for a given temperature.  Find the saturation mixing ratio, !
                  ! then convert it to specific humidity.                                  !
                  !------------------------------------------------------------------------!
                  cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the leaf-level vapour pressure deficit using canopy pressure  !
                  ! and humitdity, but leaf temperature.                                   !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                  &
                                                  , cpatch%leaf_temp(ico)                  &
                                                  , csite%can_shv   (ipa), .true.)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Copy the conductances.                                             !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico), tiny_offset)
                  cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico), tiny_offset)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Divide the values of water demand by the time step to obtain the   !
                  ! average value over the past hdid period.                               !
                  !------------------------------------------------------------------------!
                  cpatch%psi_open  (ico) = sngloff(initp%psi_open  (ico),tiny_offset)      &
                                         / sngl(hdid)
                  cpatch%psi_closed(ico) = sngloff(initp%psi_closed(ico),tiny_offset)      &
                                         / sngl(hdid)
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !    We solved leaf and branchwood together, the combined pool was re-   !
                  ! solvable but leaves weren't.  We copy the leaf temperature and liquid  !
                  ! fraction from the integrator, so they remain in thermal equilibrium    !
                  ! with branchwood.                                                       !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_temp(ico) = sngloff(initp%leaf_temp(ico) , tiny_offset)
                  cpatch%leaf_fliq(ico) = sngloff(initp%leaf_fliq(ico) , tiny_offset)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     The intercellular specific humidity is always assumed to be at     !
                  ! saturation for a given temperature.  Find the saturation mixing ratio, !
                  ! then convert it to specific humidity.                                  !
                  !------------------------------------------------------------------------!
                  cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the leaf-level vapour pressure deficit using canopy pressure  !
                  ! and humitdity, but leaf temperature.                                   !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                  &
                                                  , cpatch%leaf_temp(ico)                  &
                                                  , csite%can_shv   (ipa), .true.)
                  !------------------------------------------------------------------------!


                  !----- Set water demand and conductances to zero. -----------------------!
                  cpatch%psi_open  (ico) = 0.0
                  cpatch%psi_closed(ico) = 0.0
                  cpatch%leaf_gbh  (ico) = 0.0
                  cpatch%leaf_gbw  (ico) = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    BRANCHES.  It is always safe to copy internal energy and standing      !
               !               water,  but we must check whether branches were truly       !
               !               resolved or not before copying the other variables.         !
               !---------------------------------------------------------------------------!
               cpatch%wood_water (ico) = sngloff(initp%wood_water (ico) , tiny_offset)
               cpatch%wood_energy(ico) = sngloff(initp%wood_energy(ico) , tiny_offset)
               if (initp%wood_resolvable(ico)) then
                  !------------------------------------------------------------------------!
                  !    Branches were solved, find the temperature and liquid fraction from !
                  ! internal energy.                                                       !
                  !------------------------------------------------------------------------!
                  call uextcm2tl(cpatch%wood_energy(ico),cpatch%wood_water(ico)            &
                                ,cpatch%wood_hcap(ico),cpatch%wood_temp(ico)               &
                                ,cpatch%wood_fliq(ico))
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Copy the conductances.                                             !
                  !------------------------------------------------------------------------!
                  cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico), tiny_offset)
                  cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico), tiny_offset)
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !    We solved leaf and branchwood together, the combined pool was re-   !
                  ! solvable but leaves weren't.  We copy the leaf temperature and liquid  !
                  ! fraction from the integrator, so they remain in thermal equilibrium    !
                  ! with branchwood.                                                       !
                  !------------------------------------------------------------------------!
                  cpatch%wood_temp(ico) = sngloff(initp%wood_temp(ico) , tiny_offset)
                  cpatch%wood_fliq(ico) = sngloff(initp%wood_fliq(ico) , tiny_offset)
                  !------------------------------------------------------------------------!


                  !----- Set the conductances to zero. ------------------------------------!
                  cpatch%wood_gbh(ico) = 0.0
                  cpatch%wood_gbw(ico) = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            elseif (cpatch%hite(ico) <=  csite%total_sfcw_depth(ipa)) then
               !---------------------------------------------------------------------------!
               !    For plants buried in snow, fix the leaf and branch temperatures to the !
               ! snow temperature of the layer that is the closest to the cohort top.      !
               !---------------------------------------------------------------------------!
               kclosest = 1
               do k = csite%nlev_sfcwater(ipa), 1, -1
                  if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
               end do
               !---------------------------------------------------------------------------!


               cpatch%leaf_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
               cpatch%wood_temp(ico)   = cpatch%leaf_temp(ico)

               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico)   = 0.5
                  cpatch%wood_fliq(ico)   = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico)   = 1.0
                  cpatch%wood_fliq(ico)   = 1.0
               else
                  cpatch%leaf_fliq(ico)   = 0.0
                  cpatch%wood_fliq(ico)   = 0.0
               end if
               cpatch%leaf_water(ico)  = 0.
               cpatch%wood_water(ico)  = 0.

               !---------------------------------------------------------------------------!
               !     Find the internal energy diagnostically...                            !
               !---------------------------------------------------------------------------!
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap (ico)                 &
                                                  , cpatch%leaf_water(ico)                 &
                                                  , cpatch%leaf_temp (ico)                 &
                                                  , cpatch%leaf_fliq (ico)                 )
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap (ico)                 &
                                                  , cpatch%wood_water(ico)                 &
                                                  , cpatch%wood_temp (ico)                 &
                                                  , cpatch%wood_fliq (ico)                 )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!

               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !----- Set water demand and conductances to zero. --------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               cpatch%leaf_gbh  (ico) = 0.0
               cpatch%leaf_gbw  (ico) = 0.0
               cpatch%wood_gbh  (ico) = 0.0
               cpatch%wood_gbw  (ico) = 0.0
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !     For plants with minimal foliage or very sparse patches, fix the leaf  !
               ! and branch temperatures to the canopy air space and force leaf and branch !
               ! intercepted water to be zero.                                             !
               !---------------------------------------------------------------------------!
               cpatch%leaf_temp(ico) = csite%can_temp(ipa)
               cpatch%wood_temp(ico) = cpatch%leaf_temp(ico)

               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico) = 0.5
                  cpatch%wood_fliq(ico) = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico) = 1.0
                  cpatch%wood_fliq(ico) = 1.0
               else
                  cpatch%leaf_fliq(ico) = 0.0
                  cpatch%wood_fliq(ico) = 0.0
               end if
               cpatch%leaf_water(ico)   = 0.
               cpatch%wood_water(ico)   = 0.
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the internal energy diagnostically...                            !
               !---------------------------------------------------------------------------!
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap (ico)                 &
                                                  , cpatch%leaf_water(ico)                 &
                                                  , cpatch%leaf_temp (ico)                 &
                                                  , cpatch%leaf_fliq (ico)                 )
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap (ico)                 &
                                                  , cpatch%wood_water(ico)                 &
                                                  , cpatch%wood_temp (ico)                 &
                                                  , cpatch%wood_fliq (ico)                 )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!

               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !----- Set water demand and conductances to zero. --------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               cpatch%leaf_gbh  (ico) = 0.0
               cpatch%leaf_gbw  (ico) = 0.0
               cpatch%wood_gbh  (ico) = 0.0
               cpatch%wood_gbw  (ico) = 0.0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         case (0,2)
            !------------------------------------------------------------------------------!
            !  VEGETATION -- Leaf and branchwood were solved separately, so they are       !
            !                analysed independently.                                       !
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !  LEAVES                                                                      !
            !------------------------------------------------------------------------------!
            if (initp%leaf_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !     Leaves were solved, update water and internal energy, and re-         !
               ! calculate the temperature and leaf intercellular specific humidity.  The  !
               ! vegetation dry heat capacity is constant within one time step, so it      !
               ! doesn't need to be updated.                                               !
               !---------------------------------------------------------------------------!
               cpatch%leaf_water(ico)  = sngloff(initp%leaf_water(ico) , tiny_offset)
               cpatch%leaf_energy(ico) = sngloff(initp%leaf_energy(ico), tiny_offset)
               call uextcm2tl(cpatch%leaf_energy(ico),cpatch%leaf_water(ico)               &
                             ,cpatch%leaf_hcap(ico),cpatch%leaf_temp(ico)                  &
                             ,cpatch%leaf_fliq(ico))

               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!



               !----- Convert the wind. ---------------------------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Copy the conductances.                                                !
               !---------------------------------------------------------------------------!
               cpatch%leaf_gbh(ico) = sngloff(initp%leaf_gbh(ico), tiny_offset)
               cpatch%leaf_gbw(ico) = sngloff(initp%leaf_gbw(ico), tiny_offset)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Divide the values of water demand by the time step to obtain the      !
               ! average value over the past hdid period.                                  !
               !---------------------------------------------------------------------------!
               cpatch%psi_open  (ico) = sngloff(initp%psi_open  (ico),tiny_offset)         &
                                      / sngl(hdid)
               cpatch%psi_closed(ico) = sngloff(initp%psi_closed(ico),tiny_offset)         &
                                      / sngl(hdid)
               !---------------------------------------------------------------------------!

            elseif (cpatch%hite(ico) <=  csite%total_sfcw_depth(ipa)) then
               !---------------------------------------------------------------------------!
               !    For plants buried in snow, fix the leaf temperature to the snow        !
               ! temperature of the layer that is the closest to the leaves.               !
               !---------------------------------------------------------------------------!
               kclosest = 1
               do k = csite%nlev_sfcwater(ipa), 1, -1
                  if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
               end do
               cpatch%leaf_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico)   = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico)   = 1.0
               else
                  cpatch%leaf_fliq(ico)   = 0.0
               end if
               cpatch%leaf_water(ico)  = 0.
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap (ico)                 &
                                                  , cpatch%leaf_water(ico)                 &
                                                  , cpatch%leaf_temp (ico)                 &
                                                  , cpatch%leaf_fliq (ico)                 )
               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !----- Set water demand and conductances to zero. --------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               cpatch%leaf_gbh  (ico) = 0.0
               cpatch%leaf_gbw  (ico) = 0.0
               !---------------------------------------------------------------------------!

            else
               !---------------------------------------------------------------------------!
               !     For plants with minimal foliage or very sparse patches, fix the leaf  !
               ! temperature to the canopy air space and force leaf_water to be zero.      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_temp(ico)   = csite%can_temp(ipa)
               if (cpatch%leaf_temp(ico) == t3ple) then
                  cpatch%leaf_fliq(ico)   = 0.5
               elseif (cpatch%leaf_temp(ico) > t3ple) then
                  cpatch%leaf_fliq(ico)   = 1.0
               else
                  cpatch%leaf_fliq(ico)   = 0.0
               end if
               cpatch%leaf_water(ico)  = 0. 
               cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap (ico)                 &
                                                  , cpatch%leaf_water(ico)                 &
                                                  , cpatch%leaf_temp (ico)                 &
                                                  , cpatch%leaf_fliq (ico)                 )
               !---------------------------------------------------------------------------!
               !     The intercellular specific humidity is always assumed to be at        !
               ! saturation for a given temperature.  Find the saturation mixing ratio,    !
               ! then convert it to specific humidity.                                     !
               !---------------------------------------------------------------------------!
               cpatch%lint_shv  (ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the leaf-level vapour pressure deficit using canopy pressure     !
               ! and humitdity, but leaf temperature.                                      !
               !---------------------------------------------------------------------------!
               cpatch%leaf_vpdef(ico) = vpdefil( csite%can_prss  (ipa)                     &
                                               , cpatch%leaf_temp(ico)                     &
                                               , csite%can_shv   (ipa), .true.)
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind  (ico) = sngloff(initp%vels, tiny_offset)
               !----- Set water demand and conductances to zero. --------------------------!
               cpatch%psi_open  (ico) = 0.0
               cpatch%psi_closed(ico) = 0.0
               cpatch%leaf_gbh  (ico) = 0.0
               cpatch%leaf_gbw  (ico) = 0.0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!





            !------------------------------------------------------------------------------!
            !  WOOD                                                                        !
            !------------------------------------------------------------------------------!
            if (initp%wood_resolvable(ico)) then
               !---------------------------------------------------------------------------!
               !     Wood was solved, update water and internal energy, and recalculate    !
               ! the temperature.  The wood dry heat capacity is constant within one time  !
               ! step, so it doesn't need to be updated.                                   !
               !---------------------------------------------------------------------------!
               cpatch%wood_water(ico)  = sngloff(initp%wood_water(ico) , tiny_offset)
               cpatch%wood_energy(ico) = sngloff(initp%wood_energy(ico), tiny_offset)
               call uextcm2tl(cpatch%wood_energy(ico),cpatch%wood_water(ico)               &
                             ,cpatch%wood_hcap(ico),cpatch%wood_temp(ico)                  &
                             ,cpatch%wood_fliq(ico))

               !----- Convert the wind. ---------------------------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%veg_wind(ico),tiny_offset)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Copy the conductances.                                                !
               !---------------------------------------------------------------------------!
               cpatch%wood_gbh(ico) = sngloff(initp%wood_gbh(ico), tiny_offset)
               cpatch%wood_gbw(ico) = sngloff(initp%wood_gbw(ico), tiny_offset)
               !---------------------------------------------------------------------------!

            elseif (cpatch%hite(ico) <=  csite%total_sfcw_depth(ipa)) then
               !---------------------------------------------------------------------------!
               !    For plants buried in snow, fix the wood temperature to the snow        !
               ! temperature of the layer that is the closest to the branches.             !
               !---------------------------------------------------------------------------!
               kclosest = 1
               do k = csite%nlev_sfcwater(ipa), 1, -1
                  if (sum(csite%sfcwater_depth(1:k,ipa)) > cpatch%hite(ico)) kclosest = k
               end do
               cpatch%wood_temp(ico)   = csite%sfcwater_tempk(kclosest,ipa)
               if (cpatch%wood_temp(ico) == t3ple) then
                  cpatch%wood_fliq(ico)   = 0.5
               elseif (cpatch%wood_temp(ico) > t3ple) then
                  cpatch%wood_fliq(ico)   = 1.0
               else
                  cpatch%wood_fliq(ico)   = 0.0
               end if
               cpatch%wood_water(ico)  = 0.
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap (ico)                 &
                                                  , cpatch%wood_water(ico)                 &
                                                  , cpatch%wood_temp (ico)                 &
                                                  , cpatch%wood_fliq (ico)                 )
               !---------------------------------------------------------------------------!

               !----- Copy the meteorological wind to here. -------------------------------!
               cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               !---------------------------------------------------------------------------!


               !----- Set the conductances to zero. ---------------------------------------!
               cpatch%wood_gbh(ico) = 0.0
               cpatch%wood_gbw(ico) = 0.0
               !---------------------------------------------------------------------------!

            else
               !---------------------------------------------------------------------------!
               !     For very sparse patches of for when wood thermodynamics is off, fix   !
               ! the wood temperature to the canopy air space and force wood_water to be   !
               ! zero.                                                                     !
               !---------------------------------------------------------------------------!
               cpatch%wood_temp(ico)   = csite%can_temp(ipa)
               if (cpatch%wood_temp(ico) == t3ple) then
                  cpatch%wood_fliq(ico)   = 0.5
               elseif (cpatch%wood_temp(ico) > t3ple) then
                  cpatch%wood_fliq(ico)   = 1.0
               else
                  cpatch%wood_fliq(ico)   = 0.0
               end if
               cpatch%wood_water(ico)  = 0.
               cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap (ico)                 &
                                                  , cpatch%wood_water(ico)                 &
                                                  , cpatch%wood_temp (ico)                 &
                                                  , cpatch%wood_fliq (ico)                 )


               !----- Set the conductances to zero. ---------------------------------------!
               cpatch%wood_gbh(ico) = 0.0
               cpatch%wood_gbw(ico) = 0.0
               !---------------------------------------------------------------------------!


               !----- Copy the meteorological wind to here. -------------------------------!
               if (.not. cpatch%leaf_resolvable(ico)) then
                  cpatch%veg_wind(ico) = sngloff(initp%vels, tiny_offset)
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Final sanity check.  This should be removed soon, since it should never     !
         ! happen (well, if this still happens, then it's a bug, and we should remove the  !
         ! bug first...).                                                                  !
         !---------------------------------------------------------------------------------!
         if (cpatch%leaf_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%leaf_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL LEAF_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4site%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    '
            write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',initp%leaf_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',initp%leaf_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',initp%leaf_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',initp%leaf_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',initp%leaf_hcap(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',cpatch%leaf_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',cpatch%leaf_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',cpatch%leaf_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',cpatch%leaf_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',cpatch%leaf_hcap(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         if (cpatch%wood_temp(ico) < sngl(rk4min_veg_temp) .or.                             &
             cpatch%wood_temp(ico) > sngl(rk4max_veg_temp)   ) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'FINAL WOOD_TEMP IS WRONG IN INITP2MODELP'
            write (unit=*,fmt='(80a)')         ('-',k=1,80)
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:    ',rk4site%lon
            write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:     ',rk4site%lat
            write (unit=*,fmt='(a,1x,i6)')     ' + PATCH:        ',ipa
            write (unit=*,fmt='(a,1x,i6)')     ' + COHORT:       ',ico
            write (unit=*,fmt='(a)')           ' + PATCH AGE:    '
            write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:        ',csite%age(ipa)
            write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:  ',csite%dist_type(ipa)
            write (unit=*,fmt='(a)')           ' + BUFFER_COHORT (initp):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',initp%wood_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',initp%wood_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',initp%wood_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',initp%wood_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',initp%wood_hcap(ico)
            write (unit=*,fmt='(a)')           ' + STATE_COHORT (cpatch):'
            write (unit=*,fmt='(a,1x,es12.4)') '   - ENERGY:     ',cpatch%wood_energy(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - WATER:      ',cpatch%wood_water(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - TEMPERATURE:',cpatch%wood_temp(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - FRACLIQ:    ',cpatch%wood_fliq(ico)
            write (unit=*,fmt='(a,1x,es12.4)') '   - HEAT_CAP:   ',cpatch%wood_hcap(ico)
            write (unit=*,fmt='(80a)') ('-',k=1,80)
            call print_rk4patch(initp, csite,ipa)
            call fatal_error('extreme vegetation temperature','initp2modelp'               &
                            &,'rk4_driver.f90')
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Integrate the average state variables.  Notice that many variables (e.g.,     !
      ! temperature, density, and soil matric potential) are NOT integrated here: instead, !
      ! we find the averaged value after we normalise the average of the prognostic        !
      ! variables.  Same thing for aggregated variables at the patch level.                !
      !------------------------------------------------------------------------------------!
      csite%fmean_can_theiv  (ipa) = csite%fmean_can_theiv  (ipa)                          &
                                   + csite%can_theiv        (ipa) * dtlsm_o_frqsum
      csite%fmean_can_theta  (ipa) = csite%fmean_can_theta  (ipa)                          &
                                   + csite%can_theta        (ipa) * dtlsm_o_frqsum
      csite%fmean_can_vpdef  (ipa) = csite%fmean_can_vpdef  (ipa)                          &
                                   + csite%can_vpdef        (ipa) * dtlsm_o_frqsum
      csite%fmean_can_shv    (ipa) = csite%fmean_can_shv    (ipa)                          &
                                   + csite%can_shv          (ipa) * dtlsm_o_frqsum
      csite%fmean_can_co2    (ipa) = csite%fmean_can_co2    (ipa)                          &
                                   + csite%can_co2          (ipa) * dtlsm_o_frqsum
      csite%fmean_can_prss   (ipa) = csite%fmean_can_prss   (ipa)                          &
                                   + csite%can_prss         (ipa) * dtlsm_o_frqsum
      csite%fmean_gnd_temp   (ipa) = csite%fmean_gnd_temp   (ipa)                          &
                                   + csite%ground_temp      (ipa) * dtlsm_o_frqsum
      csite%fmean_gnd_shv    (ipa) = csite%fmean_gnd_shv    (ipa)                          &
                                   + csite%ground_shv       (ipa) * dtlsm_o_frqsum
      csite%fmean_can_ggnd   (ipa) = csite%fmean_can_ggnd   (ipa)                          &
                                   + csite%ggnet            (ipa) * dtlsm_o_frqsum
      !------------------------------------------------------------------------------------!
      !       Snow/pounding layers.  We keep track of the total, not individual layers.    !
      ! Energy will be integrated as an extensive variable, we will convert it by the      !
      ! output time only.                                                                  !
      !------------------------------------------------------------------------------------!
      do k=1,csite%nlev_sfcwater(ipa)
         csite%fmean_sfcw_depth (ipa) = csite%fmean_sfcw_depth   (ipa)                     &
                                      + csite%sfcwater_depth   (k,ipa) * dtlsm_o_frqsum
         csite%fmean_sfcw_energy(ipa) = csite%fmean_sfcw_energy  (ipa)                     &
                                      + csite%sfcwater_energy  (k,ipa)                     &
                                      * csite%sfcwater_mass    (k,ipa) * dtlsm_o_frqsum
         csite%fmean_sfcw_mass  (ipa) = csite%fmean_sfcw_mass    (ipa)                     &
                                      + csite%sfcwater_mass    (k,ipa) * dtlsm_o_frqsum
      end do
      !------ Cohort-level variables. -----------------------------------------------------!
      do ico=1,cpatch%ncohorts
         cpatch%fmean_leaf_energy(ico) = cpatch%fmean_leaf_energy(ico)                     &
                                       + cpatch%leaf_energy      (ico) * dtlsm_o_frqsum
         cpatch%fmean_leaf_water (ico) = cpatch%fmean_leaf_water (ico)                     &
                                       + cpatch%leaf_water       (ico) * dtlsm_o_frqsum
         cpatch%fmean_leaf_hcap  (ico) = cpatch%fmean_leaf_hcap  (ico)                     &
                                       + cpatch%leaf_hcap        (ico) * dtlsm_o_frqsum
         cpatch%fmean_leaf_vpdef (ico) = cpatch%fmean_leaf_vpdef (ico)                     &
                                       + cpatch%leaf_vpdef       (ico) * dtlsm_o_frqsum
         cpatch%fmean_wood_energy(ico) = cpatch%fmean_wood_energy(ico)                     &
                                       + cpatch%wood_energy      (ico) * dtlsm_o_frqsum
         cpatch%fmean_wood_water (ico) = cpatch%fmean_wood_water (ico)                     &
                                       + cpatch%wood_water       (ico) * dtlsm_o_frqsum
         cpatch%fmean_wood_hcap  (ico) = cpatch%fmean_wood_hcap  (ico)                     &
                                       + cpatch%wood_hcap        (ico) * dtlsm_o_frqsum
         cpatch%fmean_leaf_gsw   (ico) = cpatch%fmean_leaf_gsw   (ico)                     &
                                       + cpatch%leaf_gsw         (ico) * dtlsm_o_frqsum
         cpatch%fmean_leaf_gbw   (ico) = cpatch%fmean_leaf_gbw   (ico)                     &
                                       + cpatch%leaf_gbw         (ico) * dtlsm_o_frqsum
         cpatch%fmean_wood_gbw   (ico) = cpatch%fmean_wood_gbw   (ico)                     &
                                       + cpatch%wood_gbw         (ico) * dtlsm_o_frqsum
         cpatch%fmean_psi_open   (ico) = cpatch%fmean_psi_open   (ico)                     &
                                       + cpatch%psi_open         (ico) * dtlsm_o_frqsum
         cpatch%fmean_psi_closed (ico) = cpatch%fmean_psi_closed (ico)                     &
                                       + cpatch%psi_closed       (ico) * dtlsm_o_frqsum
         cpatch%fmean_fs_open    (ico) = cpatch%fmean_fs_open    (ico)                     &
                                       + cpatch%fs_open          (ico) * dtlsm_o_frqsum
         cpatch%fmean_fsw        (ico) = cpatch%fmean_fsw        (ico)                     &
                                       + cpatch%fsw              (ico) * dtlsm_o_frqsum
         cpatch%fmean_fsn        (ico) = cpatch%fmean_fsn        (ico)                     &
                                       + cpatch%fsn              (ico) * dtlsm_o_frqsum
         !---------------------------------------------------------------------------------!
         !     The penalty factor for water and nitrogen are meaningful only during the    !
         ! day.  For the daily means we must add only when it is daytime, so we integrate  !
         ! them here too.                                                                  !
         !---------------------------------------------------------------------------------!
         if (.not. nighttime .and. writing_long) then
            cpatch%dmean_fs_open    (ico) = cpatch%dmean_fs_open    (ico)                  &
                                          + cpatch%fs_open          (ico)                  &
                                          * dtlsm
            cpatch%dmean_fsw        (ico) = cpatch%dmean_fsw        (ico)                  &
                                          + cpatch%fsw              (ico)                  &
                                          * dtlsm
            cpatch%dmean_fsn        (ico) = cpatch%dmean_fsn        (ico)                  &
                                          + cpatch%fsn              (ico)                  &
                                          * dtlsm
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------ Soil variables. -------------------------------------------------------------!
      do k = rk4site%lsl, nzg
         csite%fmean_soil_energy(k,ipa) = csite%fmean_soil_energy(k,ipa)                   &
                                        + csite%soil_energy      (k,ipa)                   &
                                        * dtlsm_o_frqsum
         csite%fmean_soil_mstpot(k,ipa) = csite%fmean_soil_mstpot(k,ipa)                   &
                                        + csite%soil_mstpot      (k,ipa)                   &
                                        * dtlsm_o_frqsum
         csite%fmean_soil_water (k,ipa) = csite%fmean_soil_water (k,ipa)                   &
                                        + csite%soil_water       (k,ipa)                   &
                                        * dtlsm_o_frqsum
      end do
      !------------------------------------------------------------------------------------!
     return
   end subroutine initp2modelp
   !=======================================================================================!
   !=======================================================================================! 
end module rk4_driver

!==========================================================================================!
!==========================================================================================!
