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
   subroutine rk4_timestep(cgrid)
      use rk4_coms               , only : integration_vars           & ! structure
                                        , rk4patchtype               & ! structure
                                        , zero_rk4_patch             & ! subroutine
                                        , zero_rk4_cohort            & ! subroutine
                                        , integration_buff           ! ! intent(out)
      use ed_para_coms           , only : nthreads                   ! ! intent(in)
      use ed_state_vars          , only : edtype                     & ! structure
                                        , polygontype                & ! structure
                                        , sitetype                   ! ! structure
      use met_driver_coms        , only : met_driv_state             ! ! structure
      use grid_coms              , only : nzg                        ! ! intent(in)
      use ed_misc_coms           , only : current_time               & ! intent(in)
                                        , dtlsm                      ! ! intent(in)
      use budget_utils           , only : update_cbudget_committed   & ! function
                                        , compute_budget             ! ! function
      use soil_respiration       , only : soil_respiration_driver    ! ! sub-routine
      use stem_resp_driv         , only : stem_respiration           ! ! function
      use photosyn_driv          , only : canopy_photosynthesis      ! ! sub-routine
      use rk4_misc               , only : sanity_check_veg_energy    ! ! sub-routine
      use rk4_copy_patch         , only : copy_rk4patch_init         ! ! sub-routine
      use rk4_integ_utils        , only : copy_met_2_rk4site         ! ! sub-routine
      use update_derived_utils   , only : update_patch_derived_props ! ! sub-routine
      use plant_hydro            , only : plant_hydro_driver         ! ! sub-routine
      use therm_lib              , only : tq2enthalpy                ! ! function
      !$ use omp_lib
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(edtype)              , target      :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
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
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!



      polygonloop: do ipy = 1,cgrid%npolygons
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
            !$OMP  ipa,ita,initp,yscal,y,dydx,yerr,ytemp,ak2,ak3,ak4,ak5,ak6,ak7           &
            !$OMP ,patch_vels,old_can_prss,old_can_enthalpy,old_can_temp,old_can_shv       &
            !$OMP ,old_can_co2,old_can_rhos,old_can_dmol,ecurr_netrad,wcurr_loss2atm       &
            !$OMP ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage &
            !$OMP ,wcurr_loss2runoff,ecurr_loss2runoff,co2curr_denseffect,ecurr_denseffect &
            !$OMP ,wcurr_denseffect,ecurr_prsseffect,rshort_tot,nsteps)
            threadloop: do ibuff=1,nthreads
               !------ Update pointers. ---------------------------------------------------!
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
                  !------------------------------------------------------------------------!

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

                  !----- Compute stem respiration. ----------------------------------------!
                  call stem_respiration(csite,ipa)
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
                  !    This is the driver for the integration process...                   !
                  !------------------------------------------------------------------------!
                  call integrate_patch_rk4(csite,initp,ipa,isi,ibuff                       &
                                          ,cpoly%nighttime(isi),wcurr_loss2atm             &
                                          ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm    &
                                          ,wcurr_loss2drainage,ecurr_loss2drainage         &
                                          ,wcurr_loss2runoff,ecurr_loss2runoff             &
                                          ,co2curr_denseffect,ecurr_denseffect             &
                                          ,wcurr_denseffect,nsteps)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Add the number of steps into the step counter. Workload            !
                  ! accumulation is order-independent, so this can stay shared.            !
                  !------------------------------------------------------------------------!
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
   subroutine integrate_patch_rk4(csite,initp,ipa,isi,ibuff,nighttime,wcurr_loss2atm       &
                                 ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm             &
                                 ,wcurr_loss2drainage,ecurr_loss2drainage                  &
                                 ,wcurr_loss2runoff,ecurr_loss2runoff,co2curr_denseffect   &
                                 ,ecurr_denseffect,wcurr_denseffect,nsteps)
      use rk4_integ_utils , only : odeint               ! ! sub-routine
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use rk4_coms        , only : integration_vars     & ! structure
                                 , rk4patchtype         & ! structure
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4i               ! ! intent(inout)
      use rk4_copy_patch  , only : initp2modelp         ! ! sub-routine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      type(rk4patchtype)    , target      :: initp
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
      !------------------------------------------------------------------------------------!


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
      call odeint(csite,ipa,isi,ibuff,nsteps)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
      initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
      initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
      initp%cpwp = initp%can_dmol * initp%cpwp * dtrk4i
      initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i


      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp(tend-tbeg,initp,csite,ipa,nighttime,wcurr_loss2atm,ecurr_netrad    &
                       ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage                &
                       ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff            &
                       ,co2curr_denseffect,ecurr_denseffect,wcurr_denseffect)
      !------------------------------------------------------------------------------------!


      return
   end subroutine integrate_patch_rk4
   !=======================================================================================!
   !=======================================================================================!
end module rk4_driver

!==========================================================================================!
!==========================================================================================!
