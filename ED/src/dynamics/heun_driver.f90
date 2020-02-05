!==========================================================================================!
!==========================================================================================!
module heun_driver
  contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is the main driver for the Heun's (RK2) integration scheme.       !
   !---------------------------------------------------------------------------------------!
   subroutine heun_timestep(cgrid)
      use rk4_coms               , only : integration_vars           & ! structure
                                        , rk4patchtype               & ! structure
                                        , zero_rk4_patch             & ! subroutine
                                        , zero_rk4_cohort            & ! subroutine
                                        , integration_buff           ! ! intent(out)
      use ed_para_coms           , only : nthreads                   ! ! intent(in)
      use ed_state_vars          , only : edtype                     & ! structure
                                        , polygontype                & ! structure
                                        , sitetype                   & ! structure
                                        , patchtype                  ! ! structure
      use met_driver_coms        , only : met_driv_state             ! ! structure
      use grid_coms              , only : nzg                        ! ! intent(in)
      use ed_misc_coms           , only : current_time               & ! intent(in)
                                        , dtlsm                      ! ! intent(in)
      use ed_max_dims            , only : n_dbh                      ! ! intent(in)
      use budget_utils           , only : update_cbudget_committed   & ! function
                                        , compute_budget             ! ! function
      use soil_respiration       , only : soil_respiration_driver    ! ! function
      use photosyn_driv          , only : canopy_photosynthesis      ! ! function
      use update_derived_utils   , only : update_patch_derived_props ! ! subroutine
      use rk4_integ_utils        , only : copy_met_2_rk4site         ! ! subroutine
      use rk4_misc               , only : sanity_check_veg_energy    ! ! sub-routine
      use rk4_copy_patch         , only : copy_rk4patch_init         ! ! sub-routine
      use plant_hydro            , only : plant_hydro_driver         ! ! subroutine
      use therm_lib              , only : tq2enthalpy                ! ! function
      !$ use omp_lib

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)             , target       :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)        , pointer      :: cpoly
      type(sitetype)           , pointer      :: csite
      type(met_driv_state)     , pointer      :: cmet
      type(rk4patchtype)       , pointer      :: initp
      type(rk4patchtype)       , pointer      :: ytemp
      type(rk4patchtype)       , pointer      :: yerr
      type(rk4patchtype)       , pointer      :: yscal
      type(rk4patchtype)       , pointer      :: dydx
      type(rk4patchtype)       , pointer      :: y
      integer                                 :: ipy
      integer                                 :: isi
      integer                                 :: ipa
      integer                                 :: imon
      integer                                 :: nsteps
      real                                    :: wcurr_loss2atm
      real                                    :: ecurr_netrad
      real                                    :: ecurr_loss2atm
      real                                    :: co2curr_loss2atm
      real                                    :: wcurr_loss2drainage
      real                                    :: ecurr_loss2drainage
      real                                    :: wcurr_loss2runoff
      real                                    :: ecurr_loss2runoff
      real                                    :: co2curr_denseffect
      real                                    :: ecurr_denseffect
      real                                    :: wcurr_denseffect
      real                                    :: ecurr_prsseffect
      real                                    :: old_can_prss
      real                                    :: old_can_enthalpy
      real                                    :: old_can_temp
      real                                    :: old_can_shv
      real                                    :: old_can_co2
      real                                    :: old_can_rhos
      real                                    :: old_can_dmol
      real                                    :: patch_vels
      real                                    :: rshort_tot
      integer                                 :: ibuff
      integer                                 :: npa_thread
      integer                                 :: ita
      !----- Local constants. -------------------------------------------------------------!
      logical                   , parameter   :: test_energy_sanity = .false.
      !------------------------------------------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)


            !----- Find the number of patches per thread. ---------------------------------!
            npa_thread = ceiling(real(csite%npatches) / real(nthreads))
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Update the monthly rainfall.                                             !
            !------------------------------------------------------------------------------!
            imon                             = current_time%month
            cpoly%avg_monthly_pcpg(imon,isi) = cpoly%avg_monthly_pcpg(imon,isi)            &
                                             + cmet%pcpg * dtlsm
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Copy the meteorological variables to the rk4site structure.               !
            !------------------------------------------------------------------------------!
            call copy_met_2_rk4site(nzg,cmet%atm_ustar,cmet%atm_theiv,cmet%atm_vpdef       &
                                   ,cmet%atm_theta,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2  &
                                   ,cmet%geoht,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg  &
                                   ,cmet%prss,cmet%rshort,cmet%rlong,cmet%par_beam         &
                                   ,cmet%par_diffuse,cmet%nir_beam,cmet%nir_diffuse        &
                                   ,cmet%geoht,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)      &
                                   ,cpoly%green_leaf_factor(:,isi),cgrid%lon(ipy)          &
                                   ,cgrid%lat(ipy),cgrid%cosz(ipy))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !  MLO - Changed the parallel do loop to account for cases in which the number !
            !        of threads is less than the number of patches.                        !
            !------------------------------------------------------------------------------!
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(                                     &
            !$OMP  ipa,ita,initp,ytemp,yerr,yscal,dydx,y,patch_vels,old_can_prss           &
            !$OMP ,old_can_enthalpy,old_can_temp,old_can_shv,old_can_co2,old_can_rhos      &
            !$OMP ,old_can_dmol,ecurr_netrad,wcurr_loss2atm,ecurr_loss2atm                 &
            !$OMP ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage                &
            !$OMP ,wcurr_loss2runoff,ecurr_loss2runoff,co2curr_denseffect,ecurr_denseffect &
            !$OMP ,wcurr_denseffect,ecurr_prsseffect,rshort_tot,nsteps)
            threadloop: do ibuff=1,nthreads
               !------ Update pointers. ---------------------------------------------------!
               initp => integration_buff(ibuff)%initp
               ytemp => integration_buff(ibuff)%ytemp
               yerr  => integration_buff(ibuff)%yerr
               yscal => integration_buff(ibuff)%yscal
               dydx  => integration_buff(ibuff)%dydx
               y     => integration_buff(ibuff)%y
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Loop through tasks.  We don't assign contiguous blocks of patches to  !
               ! each thread because patches are sorted by age and older patches have more !
               ! cohorts and are likely to be slower.                                      !
               !---------------------------------------------------------------------------!
               taskloop: do ita=1,npa_thread
                  !------------------------------------------------------------------------!
                  !     Find out which patch to solve.  In case the number of patches      !
                  ! is not a perfect multiple of number of threads, some patch numbers     !
                  ! will exceed csite%npatches in the last iteration, in which we can      !
                  ! terminate the loop.                                                    !
                  !------------------------------------------------------------------------!
                  ipa = ibuff + (ita - 1) * nthreads
                  if (ipa > csite%npatches) exit taskloop
                  !------------------------------------------------------------------------!

                  !----- Reset all buffers to zero, as a safety measure. ------------------!
                  call zero_rk4_patch(initp)
                  call zero_rk4_patch(ytemp)
                  call zero_rk4_patch(yerr)
                  call zero_rk4_patch(yscal)
                  call zero_rk4_patch(dydx)
                  call zero_rk4_patch(y)
                  call zero_rk4_cohort(initp)
                  call zero_rk4_cohort(ytemp)
                  call zero_rk4_cohort(yerr)
                  call zero_rk4_cohort(yscal)
                  call zero_rk4_cohort(dydx)
                  call zero_rk4_cohort(y)

                  !----- Get velocity for aerodynamic resistance. -------------------------!
                  if (csite%can_theta(ipa) < cmet%atm_theta) then
                     patch_vels = cmet%vels_stab
                  else
                     patch_vels = cmet%vels_unstab
                  end if
                  !------------------------------------------------------------------------!



                  !----- Save the previous thermodynamic state. ---------------------------!
                  old_can_prss     = csite%can_prss(ipa)
                  old_can_enthalpy = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa)    &
                                                ,.true.)
                  old_can_temp     = csite%can_temp(ipa)
                  old_can_shv      = csite%can_shv (ipa)
                  old_can_co2      = csite%can_co2 (ipa)
                  old_can_rhos     = csite%can_rhos(ipa)
                  old_can_dmol     = csite%can_dmol(ipa)
                  !------------------------------------------------------------------------!



                  !----- Find incoming radiation used by the radiation driver. ------------!
                  if (cpoly%nighttime(isi)) then
                     rshort_tot = 0.0
                  else
                     rshort_tot = cmet%par_beam * csite%fbeam(ipa) + cmet%par_diffuse      &
                                + cmet%nir_beam * csite%fbeam(ipa) + cmet%nir_diffuse
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Test whether temperature and energy are reasonable.               !
                  !------------------------------------------------------------------------!
                  if (test_energy_sanity) then
                     call sanity_check_veg_energy(csite,ipa)
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Get plant water flow driven by plant hydraulics.  This must be     !
                  ! placed before canopy_photosynthesis, because plant_hydro_driver needs  !
                  ! fs_open from the previous timestep.                                    !
                  !------------------------------------------------------------------------!
                  call plant_hydro_driver(csite,ipa,cpoly%ntext_soil(:,isi))
                  !------------------------------------------------------------------------!



                  !----- Get photosynthesis, stomatal conductance, and transpiration. -----!
                  call canopy_photosynthesis(csite,cmet,nzg,ipa,ibuff                      &
                                            ,cpoly%ntext_soil(:,isi)                       &
                                            ,cpoly%leaf_aging_factor(:,isi)                &
                                            ,cpoly%green_leaf_factor(:,isi))
                  !------------------------------------------------------------------------!



                  !----- Compute root and heterotrophic respiration. ----------------------!
                  call soil_respiration_driver(csite,ipa,nzg,cpoly%ntext_soil(:,isi))
                  !------------------------------------------------------------------------!



                  !----- Update the committed carbon change pool. -------------------------!
                  call update_cbudget_committed(csite,ipa)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Set up the integration patch.                                      !
                  !------------------------------------------------------------------------!
                  call copy_rk4patch_init(csite,ipa,ibuff,initp,patch_vels                 &
                                         ,old_can_enthalpy,old_can_rhos,old_can_dmol       &
                                         ,ecurr_prsseffect)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     This is the step in which the derivatives are computed, we use a   !
                  ! structure that is very similar to the Runge-Kutta, though a simpler    !
                  ! one.                                                                   !
                  !------------------------------------------------------------------------!
                  call integrate_patch_heun(csite,ipa,isi,ibuff,cpoly%nighttime(isi)       &
                                           ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm     &
                                           ,co2curr_loss2atm,wcurr_loss2drainage           &
                                           ,ecurr_loss2drainage,wcurr_loss2runoff          &
                                           ,ecurr_loss2runoff,co2curr_denseffect           &
                                           ,ecurr_denseffect,wcurr_denseffect,nsteps)
                  !------------------------------------------------------------------------!



                  !----- Add the number of steps into the step counter. -------------------!
                  cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !   Update the minimum monthly temperature, based on canopy temperature. !
                  !------------------------------------------------------------------------!
                  if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                     cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Update roughness and canopy depth.  This should be done after the   !
                  ! integration.                                                           !
                  !------------------------------------------------------------------------!
                  call update_patch_derived_props(csite,ipa,.false.)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Compute the residuals.                                             !
                  !------------------------------------------------------------------------!
                  call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg            &
                                     ,rshort_tot,cmet%rlong,ipa,wcurr_loss2atm             &
                                     ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm         &
                                     ,wcurr_loss2drainage,ecurr_loss2drainage              &
                                     ,wcurr_loss2runoff,ecurr_loss2runoff                  &
                                     ,co2curr_denseffect,ecurr_denseffect,wcurr_denseffect &
                                     ,ecurr_prsseffect,cpoly%area(isi)                     &
                                     ,cgrid%cbudget_nep(ipy),old_can_prss                  &
                                     ,old_can_enthalpy,old_can_temp,old_can_shv            &
                                     ,old_can_co2,old_can_rhos,old_can_dmol)
                  !------------------------------------------------------------------------!
               end do taskloop
               !---------------------------------------------------------------------------!
            end do threadloop
            !$OMP END PARALLEL DO
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine heun_timestep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process using the Heun method.  Notice !
   ! that most of the Heun method utilises the subroutines from Runge-Kutta.               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_heun(csite,ipa,isi,ibuff,nighttime,wcurr_loss2atm            &
                                  ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm            &
                                  ,wcurr_loss2drainage,ecurr_loss2drainage                 &
                                  ,wcurr_loss2runoff,ecurr_loss2runoff,co2curr_denseffect  &
                                  ,ecurr_denseffect,wcurr_denseffect,nsteps)
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use rk4_coms        , only : integration_vars     & ! structure
                                 , integration_buff     & ! structure
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4i               ! ! intent(inout)
      use rk4_copy_patch  , only : initp2modelp         ! ! subroutine

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      integer               , intent(in)  :: ipa
      integer               , intent(in)  :: isi
      integer               , intent(in)  :: ibuff
      logical               , intent(in)  :: nighttime
      real                  , intent(out) :: wcurr_loss2atm
      real                  , intent(out) :: ecurr_netrad
      real                  , intent(out) :: ecurr_loss2atm
      real                  , intent(out) :: co2curr_loss2atm
      real                  , intent(out) :: wcurr_loss2drainage
      real                  , intent(out) :: ecurr_loss2drainage
      real                  , intent(out) :: wcurr_loss2runoff
      real                  , intent(out) :: ecurr_loss2runoff
      real                  , intent(out) :: co2curr_denseffect
      real                  , intent(out) :: ecurr_denseffect
      real                  , intent(out) :: wcurr_denseffect
      integer               , intent(out) :: nsteps
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                        :: hbeg
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Initial step size.  Experience has shown that giving this too large a value   !
      ! causes the integrator to fail (e.g., soil layers become supersaturated).           !
      !------------------------------------------------------------------------------------!
      hbeg = dble(csite%htry(ipa))

      !------------------------------------------------------------------------------------!
      !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm, !
      ! so they must be zeroed at each call.                                               !
      !------------------------------------------------------------------------------------!
      integration_buff(ibuff)%initp%upwp = 0.d0
      integration_buff(ibuff)%initp%tpwp = 0.d0
      integration_buff(ibuff)%initp%qpwp = 0.d0
      integration_buff(ibuff)%initp%cpwp = 0.d0
      integration_buff(ibuff)%initp%wpwp = 0.d0

      !----- Go into the ODE integrator using Euler. --------------------------------------!
      call heun_integ(hbeg,csite,ipa,isi,ibuff,nsteps)

      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      integration_buff(ibuff)%initp%upwp = integration_buff(ibuff)%initp%can_rhos          &
                                         * integration_buff(ibuff)%initp%upwp * dtrk4i
      integration_buff(ibuff)%initp%tpwp = integration_buff(ibuff)%initp%can_rhos          &
                                         * integration_buff(ibuff)%initp%tpwp * dtrk4i
      integration_buff(ibuff)%initp%qpwp = integration_buff(ibuff)%initp%can_rhos          &
                                         * integration_buff(ibuff)%initp%qpwp * dtrk4i
      integration_buff(ibuff)%initp%cpwp = integration_buff(ibuff)%initp%can_dmol          &
                                         * integration_buff(ibuff)%initp%cpwp * dtrk4i
      integration_buff(ibuff)%initp%wpwp = integration_buff(ibuff)%initp%can_rhos          &
                                         * integration_buff(ibuff)%initp%wpwp * dtrk4i

      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp(tend-tbeg,integration_buff(ibuff)%initp,csite,ipa,nighttime        &
                       ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm        &
                       ,wcurr_loss2drainage,ecurr_loss2drainage,wcurr_loss2runoff          &
                       ,ecurr_loss2runoff,co2curr_denseffect,ecurr_denseffect              &
                       ,wcurr_denseffect)
      !------------------------------------------------------------------------------------!

      return
   end subroutine integrate_patch_heun
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! Subroutine heun_integ                                                                 !
   !                                                                                       !
   !     This subroutine will drive the integration of several ODEs that drive the fast-   !
   ! -scale state variables using the Heun's method (a 2nd order Runge-Kutta).             !
   !---------------------------------------------------------------------------------------!
   subroutine heun_integ(h1,csite,ipa,isi,ibuff,nsteps)
      use ed_state_vars  , only : sitetype                  & ! structure
                                , patchtype                 & ! structure
                                , polygontype               ! ! structure
      use rk4_coms       , only : integration_vars          & ! structure
                                , rk4patchtype              & ! structure
                                , integration_buff          & ! structure
                                , rk4site                   & ! intent(in)
                                , print_detailed            & ! intent(in)
                                , maxstp                    & ! intent(in)
                                , tbeg                      & ! intent(in)
                                , tend                      & ! intent(in)
                                , dtrk4                     & ! intent(in)
                                , dtrk4i                    & ! intent(in)
                                , tiny_offset               & ! intent(in)
                                , checkbudget               & ! intent(in)
                                , zero_rk4_patch            & ! subroutine
                                , zero_rk4_cohort           & ! subroutine
                                , hmin                      & ! intent(in)
                                , rk4eps                    & ! intent(in)
                                , rk4epsi                   & ! intent(in)
                                , safety                    & ! intent(in)
                                , pgrow                     & ! intent(in)
                                , pshrnk                    & ! intent(in)
                                , norm_rk4_fluxes           & ! sub-routine
                                , reset_rk4_fluxes          ! ! sub-routine
      use rk4_derivs     , only : leaf_derivs               ! ! subroutine
      use rk4_copy_patch , only : copy_rk4_patch            ! ! sub-routine
      use rk4_integ_utils, only : get_yscal                 & ! sub-routine
                                , get_errmax                & ! sub-routine
                                , inc_rk4_patch             & ! sub-routine
                                , rk4_sanity_check          & ! sub-routine
                                , print_sanity_check        ! ! sub-routine
      use rk4_misc       , only : print_rk4patch            & ! sub-routine
                                , print_errmax              & ! sub-routine
                                , adjust_veg_properties     & ! sub-routine
                                , adjust_topsoil_properties & ! sub-routine
                                , adjust_sfcw_properties    & ! sub-routine
                                , update_diagnostic_vars    & ! sub-routine
                                , update_density_vars       & ! sub-routine
                                , print_rk4_state           ! ! sub-routine
      use ed_misc_coms   , only : fast_diagnostics          & ! intent(in)
                                , dtlsm                     ! ! intent(in)
      use grid_coms      , only : nzg                       & ! intent(in)
                                , nzs                       & ! intent(in)
                                , time                      ! ! intent(in)
      use soil_coms      , only : runoff_time_i             & ! intent(in)
                                , simplerunoff              ! ! intent(in)
      use consts_coms    , only : t3ple8                    & ! intent(in)
                                , wdnsi8                    ! ! intent(in)
      use therm_lib8     , only : tl2uint8                  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)            , target      :: csite         !< Current site
      integer                   , intent(in)  :: ipa           !< Current patch ID
      integer                   , intent(in)  :: isi           !< Current site ID
      integer                   , intent(in)  :: ibuff         !< Current thread
      real(kind=8)              , intent(in)  :: h1            !< First guess of delta-t
      integer                   , intent(out) :: nsteps        !< Number of steps taken.
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)           , pointer     :: cpatch        !< Current patch
      type(rk4patchtype)        , pointer     :: initp         !< Current RK4 patch
      type(rk4patchtype)        , pointer     :: ylast         !< Last RK4 patch
      logical                                 :: reject_step   !< Flag: reject the step
      logical                                 :: reject_result !< Flag: reject the result
      logical                                 :: minstep       !< Minimum time step reached
      logical                                 :: stuck         !< Flag: tiny step
      logical                                 :: test_reject   !< Reject the test
      logical                                 :: gapstep       !< This a "gap-filler" step
      integer                                 :: i             !< Step counter
      integer                                 :: k             !< Format counter
      integer                                 :: ksn           !< # of snow/water layers
      real(kind=8)                            :: x             !< Elapsed time
      real(kind=8)                            :: xnew          !< Elapsed time + h
      real(kind=8)                            :: newh          !< New time step suggested
      real(kind=8)                            :: oldh          !< Old time step
      real(kind=8)                            :: h             !< Current delta-t attempt
      real(kind=8)                            :: fgrow         !< Delta-t increase factor
      real(kind=8)                            :: hgoal         !< Delta-t ignoring overstep
      real(kind=8)                            :: hnext         !< Next delta-t
      real(kind=8)                            :: hdid          !< delta-t that worked (???)
      real(kind=8)                            :: qwfree        !< Free water int. energy
      real(kind=8)                            :: wfreeb        !< Free water
      real(kind=8)                            :: errmax        !< Max. error of this step
      real(kind=8)                            :: elaptime      !< Absolute elapsed time
      !----- External function. -----------------------------------------------------------!
      real                      , external    :: sngloff
      !------------------------------------------------------------------------------------!


      !----- Use some aliases for simplicity. ---------------------------------------------!
      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !     Create temporary patches.                                                      !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(integration_buff(ibuff)%initp,integration_buff(ibuff)%y,cpatch)

      !------------------------------------------------------------------------------------!
      ! Set initial time and stepsize.                                                     !
      !------------------------------------------------------------------------------------!
      x = tbeg
      h = h1
      if (dtrk4 < 0.d0) h = -h1

      !----- Define total elapsed time. ---------------------------------------------------!
      elaptime = time + x

      !------------------------------------------------------------------------------------!
      ! Begin timestep loop                                                                !
      !------------------------------------------------------------------------------------!
      timesteploop: do i=1,maxstp


         !----- Get initial derivatives ---------------------------------------------------!
         call leaf_derivs(integration_buff(ibuff)%y,integration_buff(ibuff)%dydx,csite,ipa &
                         ,ibuff,h,.false.)

         !----- Get scalings used to determine stability ----------------------------------!
         call get_yscal(integration_buff(ibuff)%y,integration_buff(ibuff)%dydx,h           &
                       ,integration_buff(ibuff)%yscal,cpatch)

         !----- Be sure not to overstep ---------------------------------------------------!
         hgoal   = h
         gapstep = (x+h-tend)*(x+h-tbeg) > 0.d0
         if (gapstep) h = tend -x
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Here we will perform the Heun's integration using the time step.  As in     !
         ! Runge-Kutta, we also check whether the integration is going well and if needed  !
         ! we shrink the intermediate time steps.                                          !
         !---------------------------------------------------------------------------------!
         reject_step =  .false.
         hstep:   do

            !------------------------------------------------------------------------------!
            ! 1. Try a step of varying size.                                               !
            !------------------------------------------------------------------------------!
            call heun_stepper(h,csite,ipa,ibuff,reject_step,reject_result)

            !------------------------------------------------------------------------------!
            !     Here we check the error of this step.  Three outcomes are possible:      !
            ! 1.  The updated values make no sense.  Reject step, assign a large error and !
            !     try again with a smaller time step;                                      !
            ! 2.  The updated values are reasonable, but the error is large.  Reject step  !
            !     and try again with a smaller time step;                                  !
            ! 3.  The updated values are reasonable, and the error is small.  Accept step  !
            !     and try again with a larger time step.                                   !
            !------------------------------------------------------------------------------!
            if (reject_step .or. reject_result) then
               !---------------------------------------------------------------------------!
               !    If step was already rejected, that means the step had finished pre-    !
               ! maturely, so we assign a standard large error (10.0).                     !
               !---------------------------------------------------------------------------!
               errmax = 1.d1
            else
               call get_errmax(errmax,integration_buff(ibuff)%yerr                         &
                              ,integration_buff(ibuff)%yscal,csite%patch(ipa)              &
                              ,integration_buff(ibuff)%ytemp)
               !----- Scale the error based on the prescribed tolerance. ------------------!
               errmax = errmax * rk4epsi
            end if


            !------------------------------------------------------------------------------!
            ! 3. If that error was large, then calculate a new step size to try.  There    !
            !    are two types of new tries.  If step failed to be minimally reasonable    !
            !    (rejected) we have assigned a standard large error (10.0).  Otherwise a   !
            !    new step is calculated based on the size of that error.  Hopefully, those !
            !    new steps should be less than the previous h.  If the error was small,    !
            !    i.e. less than rk4eps, then we are done with this step, and we can move   !
            !    forward time: x = x + h                                                   !
            !------------------------------------------------------------------------------!
            if (errmax > 1.d0) then
               !---------------------------------------------------------------------------!
               !      We are about to shrink the time step, so this can't be considered    !
               ! gap step.                                                                 !
               !---------------------------------------------------------------------------!
               gapstep = .false.
               !---------------------------------------------------------------------------!


               !----- Define new step and check whether it is sufficiently long or not. ---!
               oldh    = h
               newh    = safety * h * errmax**pshrnk
               minstep = (newh == h) .or. newh < hmin
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Define next time, and check whether it is long enough to change time.  !
               !---------------------------------------------------------------------------!
               h       = max(1.d-1*h, newh)
               hgoal   = h
               xnew    = x + h
               stuck   = xnew == x
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               ! 3a. Here is the moment of truth... If we reached a tiny step and yet the  !
               !     model didn't converge, then we print various values to inform the     !
               !     user and abort the run.  Please, don't hate the messenger.            !
               !---------------------------------------------------------------------------!
               if (minstep .or. stuck) then

                  write (unit=*,fmt='(80a)')         ('=',k=1,80)
                  write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN HEUN_INTEG'
                  write (unit=*,fmt='(80a)')         ('-',k=1,80)
                  write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:     ',rk4site%lon
                  write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:      ',rk4site%lat
                  write (unit=*,fmt='(a)')           ' + PATCH INFO:    '
                  write (unit=*,fmt='(a,1x,i6)')     '   - NUMBER:      ',ipa
                  write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:         ',csite%age(ipa)
                  write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:   '                  &
                                                                      ,csite%dist_type(ipa)
                  write (unit=*,fmt='(a,1x,l1)')     ' + REJECT_STEP:   ',reject_step
                  write (unit=*,fmt='(a,1x,l1)')     ' + REJECT_RESULT: ',reject_result
                  write (unit=*,fmt='(a,1x,l1)')     ' + MINSTEP:       ',minstep
                  write (unit=*,fmt='(a,1x,l1)')     ' + STUCK:         ',stuck
                  write (unit=*,fmt='(a,1x,es12.4)') ' + ERRMAX:        ',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' + X:             ',x
                  write (unit=*,fmt='(a,1x,es12.4)') ' + H:             ',h
                  write (unit=*,fmt='(a,1x,es12.4)') ' + OLDH:          ',oldh
                  write (unit=*,fmt='(a,1x,es12.4)') ' + NEWH:          ',newh
                  write (unit=*,fmt='(a,1x,es12.4)') ' + SAFETY:        ',safety
                  write (unit=*,fmt='(80a)') ('-',k=1,80)
                  if (reject_step) then
                     write (unit=*,fmt='(a)') '   Likely to be a rejected step problem.'
                     write (unit=*,fmt='(80a)') ('=',k=1,80)
                  else
                     write (unit=*,fmt='(a)') '   Likely to be an errmax problem.'
                     write (unit=*,fmt='(80a)') ('=',k=1,80)
                  end if

                  if (reject_result) then
                     !----- Run the LSM sanity check but this time we force the print. ----!
                     call rk4_sanity_check(ibuff,integration_buff(ibuff)%ytemp,test_reject &
                                          ,csite,ipa,integration_buff(ibuff)%dydx,h,.true.)
                     call print_sanity_check(integration_buff(ibuff)%y,csite,ipa)
                  elseif (reject_step) then
                     call rk4_sanity_check(ibuff,integration_buff(ibuff)%ak3,test_reject   &
                                          ,csite,ipa,integration_buff(ibuff)%dydx,h,.true.)
                     call print_sanity_check(integration_buff(ibuff)%y,csite,ipa)
                  else
                     call print_errmax(errmax,integration_buff(ibuff)%yerr                 &
                                      ,integration_buff(ibuff)%yscal,csite%patch(ipa)      &
                                      ,integration_buff(ibuff)%y)
                     write (unit=*,fmt='(80a)') ('=',k=1,80)
                     write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax
                     write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax*rk4eps
                     write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                     write (unit=*,fmt='(80a)') ('=',k=1,80)
                  end if
                  call print_rk4patch(integration_buff(ibuff)%y, csite,ipa)
               end if

            else
               !---------------------------------------------------------------------------!
               ! 3b.  Great, it worked, so now we can advance to the next step.  We just   !
               !      need  to make some minor adjustments.                                !
               !---------------------------------------------------------------------------!
               !----- i.   Update leaf properties to avoid negative water. ----------------!
               call adjust_veg_properties(integration_buff(ibuff)%ytemp,h,csite,ipa,ibuff)
               !----- ii.  Update of top soil properties to avoid off-bounds moisture. ----!
               call adjust_topsoil_properties(integration_buff(ibuff)%ytemp,h)
               !----- iii.  Make snow layers stable and positively defined. ---------------!
               call adjust_sfcw_properties(nzg,nzs,integration_buff(ibuff)%ytemp,h         &
                                          ,csite,ipa)
               !----- iv. Update the diagnostic variables. --------------------------------!
               call update_diagnostic_vars(integration_buff(ibuff)%ytemp,csite,ipa,ibuff)
               !----- v.  Update the density variables. -----------------------------------!
               call update_density_vars(integration_buff(ibuff)%ytemp                      &
                                       ,integration_buff(ibuff)%y     )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               ! 3c. Set up h for the next time.  And here we can relax h for the next     !
               !     step, and try something faster, unless this is a "gap step" (shorter  !
               !     time step just to close the full thermodynamic time step).            !
               !---------------------------------------------------------------------------!
               if (gapstep) then
                  hnext = hgoal
               else
                  fgrow = min(5.d0,max(safety*errmax**pgrow,1.d0))
                  hnext = max(2.d0*hmin, min(dble(dtlsm), fgrow * hgoal))
               end if
               !---------------------------------------------------------------------------!



               !------ 3d. Normalise the fluxes if the user wants detailed debugging. -----!
               if (print_detailed) then
                  call norm_rk4_fluxes(integration_buff(ibuff)%ytemp,h)
                  call print_rk4_state(integration_buff(ibuff)%y                           &
                                      ,integration_buff(ibuff)%ytemp,csite,ipa,isi,x,h)
               end if

               !----- 3e. Copy the temporary structure to the intermediate state. ---------!
               call copy_rk4_patch(integration_buff(ibuff)%ytemp,integration_buff(ibuff)%y &
                                  ,csite%patch(ipa))

               !---------------------------------------------------------------------------!
               !    3f. Flush step-by-step fluxes to zero if the user wants detailed       !
               !        debugging.                                                         !
               !---------------------------------------------------------------------------!
               if (print_detailed) then
                  call reset_rk4_fluxes(integration_buff(ibuff)%y)
               end if

               !----- 3g. Update time. ----------------------------------------------------!
               x        = x + h
               hdid     = h
               h        = hnext
               elaptime = elaptime + h

               exit hstep
            end if
            !------------------------------------------------------------------------------!
         end do hstep
         !---------------------------------------------------------------------------------!



         !----- If the integration reached the next step, make some final adjustments -----!
         if ((x-tend)*dtrk4 >= 0.d0) then

            !----- Use pointers to simplify code. -----------------------------------------!
            initp => integration_buff(ibuff)%initp
            ylast => integration_buff(ibuff)%y
            !------------------------------------------------------------------------------!


            !------ Copy the temporary patch to the next intermediate step ----------------!
            call copy_rk4_patch(ylast,initp,cpatch)
            !------------------------------------------------------------------------------!

            !----- Number of temporary surface water layers. ------------------------------!
            ksn = initp%nlev_sfcwater
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !   Make temporary surface liquid water disappear.  This will not happen       !
            ! immediately, but liquid water will decay with the time scale defined by      !
            ! runoff_time scale. If the time scale is too tiny, then it will be forced to  !
            ! be hdid (no reason to be faster than that).                                  !
            !------------------------------------------------------------------------------!
            if ( simplerunoff .and. ksn >= 1) then
               !---------------------------------------------------------------------------!
               !    Test mass and liquid fraction inside the if block to avoid bound       !
               ! check errors.                                                             !
               !---------------------------------------------------------------------------!
               if ( initp%sfcwater_mass   (ksn) > 0.d0  .and.                              &
                    initp%sfcwater_fracliq(ksn) > 1.d-1        ) then

                  !----- Find the amount of water to be lost as runoff. -------------------!
                  wfreeb = min(1.d0, dtrk4 * runoff_time_i) * initp%sfcwater_mass(ksn)     &
                         * (initp%sfcwater_fracliq(ksn) - 1.d-1) / 9.d-1
                  qwfree = wfreeb * tl2uint8(initp%sfcwater_tempk(ksn),1.d0)
                  !------------------------------------------------------------------------!



                  !----- Update the state variables. --------------------------------------!
                  initp%sfcwater_mass  (ksn) = initp%sfcwater_mass  (ksn) - wfreeb
                  initp%sfcwater_depth (ksn) = initp%sfcwater_depth (ksn) - wfreeb * wdnsi8
                  initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn) - qwfree
                  !------------------------------------------------------------------------!



                  !----- Compute runoff for output ----------------------------------------!
                  if (fast_diagnostics) then
                     !---------------------------------------------------------------------!
                     !      There is no need to divide wfreeb and qwfree  by time step,    !
                     ! which will be done in subroutine normalize_averaged_vars.           !
                     !---------------------------------------------------------------------!
                     csite%runoff       (ipa) = csite%runoff(ipa)                          &
                                              + sngloff(wfreeb * dtrk4i,tiny_offset)
                     csite%fmean_runoff (ipa) = csite%fmean_runoff(ipa)                    &
                                              + sngloff(wfreeb,tiny_offset)
                     csite%fmean_qrunoff(ipa) = csite%fmean_qrunoff(ipa)                   &
                                              + sngloff(qwfree,tiny_offset)
                  end if
                  if (checkbudget) then
                     !---------------------------------------------------------------------!
                     !      To make sure that the previous values of wbudget_loss2runoff   !
                     ! and ebudget_loss2runoff are accumulated to the next time step.      !
                     !---------------------------------------------------------------------!
                     initp%wbudget_loss2runoff = initp%wbudget_loss2runoff + wfreeb
                     initp%ebudget_loss2runoff = initp%ebudget_loss2runoff + qwfree
                     initp%wbudget_storage     = initp%wbudget_storage     - wfreeb
                     initp%ebudget_storage     = initp%ebudget_storage     - qwfree
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!


                  !----- Make sure all diagnostic variables are consistent and bounded. ---!
                  call adjust_sfcw_properties(nzg,nzs,initp,dtrk4,csite,ipa)
                  call update_diagnostic_vars(initp,csite,ipa,ibuff)
                  call update_density_vars   (initp,ylast)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!



            !------ Update the substep for next time and leave ----------------------------!
            csite%hprev(ipa) = csite%htry(ipa)
            csite%htry(ipa)  = sngl(hnext)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Update the average time step.  The square of DTLSM (tend-tbeg) is needed !
            ! because we will divide this by the time between t0 and t0+frqsum.            !
            !------------------------------------------------------------------------------!
            csite%fmean_rk4step(ipa) = csite%fmean_rk4step(ipa)                            &
                                     + sngl((tend-tbeg)*(tend-tbeg))/real(i)
            nsteps = i
            !------------------------------------------------------------------------------!
            return
         end if
         !---------------------------------------------------------------------------------!

         !----- Use hnext as the next substep ---------------------------------------------!
         h = hnext
         !---------------------------------------------------------------------------------!
      end do timesteploop
      !------------------------------------------------------------------------------------!


      !----- If it reached this point, that is really bad news... -------------------------!
      write (unit=*,fmt='(a)') ' ==> Too many steps in routine heun_integ'
      call print_rk4patch(integration_buff(ibuff)%ytemp, csite,ipa)
      !------------------------------------------------------------------------------------!

      return
   end subroutine heun_integ
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the integration step of the Heun's method.                !
   ! Here we use several buffers, so a quick guide of what each one means...               !
   !                                                                                       !
   ! 1. integration_buff(ibuff)%y      => This is the state y(t)                           !
   ! 2. integration_buff(ibuff)%ytemp  => This is the state y(t+h)                         !
   ! 3. integration_buff(ibuff)%yerr   => This is the error e(t+h)                         !
   ! 4. integration_buff(ibuff)%dydx   => This is the Runge-Kutta term K1                  !
   ! 5. integration_buff(ibuff)%ak2    => This is the Runge-Kutta term K2                  !
   ! 6. integration_buff(ibuff)%ak3    => This is the next step using Euler:               !
   !                                      ye(t+h) = y(t) + K1 * h                          !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine heun_stepper(h,csite,ipa,ibuff,reject_step,reject_result)
      use rk4_coms       , only : integration_buff       & ! structure
                                , integration_vars       & ! structure
                                , zero_rk4_patch         & ! subroutine
                                , zero_rk4_cohort        & ! subroutine
                                , print_diags            & ! intent(in)
                                , heun_a2                & ! intent(in)
                                , heun_b21               & ! intent(in)
                                , heun_c1                & ! intent(in)
                                , heun_c2                & ! intent(in)
                                , heun_dc1               & ! intent(in)
                                , heun_dc2               ! ! intent(in)
      use ed_state_vars  , only : sitetype               & ! structure
                                , patchtype              ! ! structure
      use rk4_copy_patch , only : copy_rk4_patch         ! ! sub-routine
      use rk4_derivs     , only : leaf_derivs            ! ! sub-routine
      use rk4_integ_utils, only : inc_rk4_patch          & ! sub-routine
                                , rk4_sanity_check       ! ! sub-routine
      use rk4_misc       , only : update_diagnostic_vars & ! sub-routine
                                , print_rk4patch         ! ! sub-routine
                                

      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)    , target      :: csite
      integer           , intent(in)  :: ipa
      integer           , intent(in)  :: ibuff
      logical           , intent(out) :: reject_step
      logical           , intent(out) :: reject_result
      real(kind=8)      , intent(in)  :: h
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Start and assume that nothing went wrong up to this point... If we find any    !
      ! seriously bad step, quit and reduce the time step without even bothering to try    !
      ! further.                                                                           !
      !------------------------------------------------------------------------------------!
      reject_step   = .false.
      reject_result = .false.
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!



      !----- yeuler is the temporary array with the Euler step with no correction. --------!
      call copy_rk4_patch(integration_buff(ibuff)%y,integration_buff(ibuff)%ak3,cpatch)
      call inc_rk4_patch (integration_buff(ibuff)%ak3,integration_buff(ibuff)%dydx         &
                         ,heun_b21*h, cpatch)
      call update_diagnostic_vars(integration_buff(ibuff)%ak3      ,csite,ipa,ibuff)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Check to see if the Euler result makes sense.  Since we will compute the        !
      ! derivative correction using it, the Euler step must be bounded.  If not, reject    !
      ! the step and try a smaller step size.                                              !
      !------------------------------------------------------------------------------------!
      call rk4_sanity_check(ibuff,integration_buff(ibuff)%ak3,reject_step,csite,ipa        &
                           ,integration_buff(ibuff)%dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute the second term (correction) of the derivative, using the Euler's      !
      ! predicted state.                                                                   !
      !------------------------------------------------------------------------------------!
      call leaf_derivs(integration_buff(ibuff)%ak3,integration_buff(ibuff)%ak2, csite,ipa  &
                      ,ibuff,h,.false.)
      !------------------------------------------------------------------------------------!



      !----- We now combine both derivatives and find the potential next step. ------------!
      call copy_rk4_patch(integration_buff(ibuff)%y,integration_buff(ibuff)%ytemp, cpatch)
      call inc_rk4_patch(integration_buff(ibuff)%ytemp,integration_buff(ibuff)%dydx        &
                        ,heun_c1*h, cpatch)
      call inc_rk4_patch(integration_buff(ibuff)%ytemp,integration_buff(ibuff)%ak2         &
                        ,heun_c2*h, cpatch)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update the diagnostic properties and make final adjustments.  This time we    !
      ! will run the full adjustment, to make sure that the step will be rejected          !
      ! especially if there are issues with the top soil properties.                       !
      !------------------------------------------------------------------------------------!
      call update_diagnostic_vars(integration_buff(ibuff)%ytemp,csite,ipa,ibuff)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Check to see if this attempt of advancing one time step makes sense.  If not,   !
      ! reject the result and try a smaller step size.                                     !
      !------------------------------------------------------------------------------------!
      call rk4_sanity_check(ibuff,integration_buff(ibuff)%ytemp,reject_result,csite,ipa    &
                           ,integration_buff(ibuff)%ak2,h,print_diags)
      if(reject_result)return
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute the estimate of the error associated with the step.                    !
      !------------------------------------------------------------------------------------!
      call zero_rk4_patch (integration_buff(ibuff)%yerr)
      call zero_rk4_cohort(integration_buff(ibuff)%yerr)
      call inc_rk4_patch(integration_buff(ibuff)%yerr,integration_buff(ibuff)%dydx         &
                        ,heun_dc1*h,cpatch)
      call inc_rk4_patch(integration_buff(ibuff)%yerr,integration_buff(ibuff)%ak2          &
                        ,heun_dc2*h,cpatch)
      !------------------------------------------------------------------------------------!

      return
   end subroutine heun_stepper
   !=======================================================================================!
   !=======================================================================================!

end module heun_driver
!==========================================================================================!
!==========================================================================================!
