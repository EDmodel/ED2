module hybrid_driver
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is the main driver for the Forward/Backward (FB)                  !
   !     Euler integration scheme.                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine hybrid_timestep(cgrid)
     use rk4_coms              , only : integration_vars           & ! structure
                                      , rk4patchtype               & ! structure
                                      , zero_rk4_patch             & ! subroutine
                                      , zero_rk4_cohort            & ! subroutine
                                      , zero_bdf2_patch            &
                                      , integration_buff           & ! intent(out)
                                      , bdf2patchtype              &
                                      , tbeg                       &
                                      , tend                       &
                                      , dtrk4i                     !
     use ed_para_coms          , only : nthreads                   ! ! intent(in)
     use rk4_driver            , only : initp2modelp
     use ed_state_vars         , only : edtype                     & ! structure
                                      , polygontype                & ! structure
                                      , sitetype                   ! ! structure
     use met_driver_coms       , only : met_driv_state             ! ! structure
     use grid_coms             , only : nzg                        & ! intent(in)
                                      , nzs                        ! ! intent(in)
     use ed_misc_coms          , only : current_time               & ! intent(in)
                                      , dtlsm                      ! ! intent(in)
     use budget_utils          , only : update_cbudget_committed   & ! function
                                      , compute_budget             ! ! function
     use soil_respiration      , only : soil_respiration_driver    ! ! function
     use photosyn_driv         , only : canopy_photosynthesis      ! ! function
     use update_derived_utils  , only : update_patch_thermo_props  & ! subroutine
                                      , update_patch_derived_props ! ! subroutine
     use rk4_integ_utils       , only : copy_met_2_rk4site         & ! subroutine
                                      , rk4_sanity_check           ! ! subroutine
     use rk4_misc              , only : copy_patch_init            & ! subroutine
                                      , sanity_check_veg_energy    ! ! subroutine
     use plant_hydro           , only : plant_hydro_driver         ! ! subroutine
     use rk4_copy_patch        , only : copy_rk4_patch             ! ! subroutine
     use therm_lib             , only : tq2enthalpy                ! ! function

     !$  use omp_lib

     implicit none
     !----- Arguments ---------------------------------------------------------------------!
     type(edtype)             , target      :: cgrid
     !----- Local variables ---------------------------------------------------------------!
     type(polygontype)        , pointer     :: cpoly
     type(sitetype)           , pointer     :: csite
     type(met_driv_state)     , pointer     :: cmet
     type(rk4patchtype)       , pointer     :: initp
     type(rk4patchtype)       , pointer     :: dinitp
     type(rk4patchtype)       , pointer     :: ytemp
     type(bdf2patchtype)      , pointer     :: yprev
     integer                                :: ipy
     integer                                :: isi
     integer                                :: ipa
     integer                                :: imon
     integer                                :: nsteps
     real                                   :: patch_vels
     real                                   :: wcurr_loss2atm
     real                                   :: ecurr_loss2atm
     real                                   :: co2curr_loss2atm
     real                                   :: ecurr_netrad
     real                                   :: wcurr_loss2drainage
     real                                   :: ecurr_loss2drainage
     real                                   :: wcurr_loss2runoff
     real                                   :: ecurr_loss2runoff
     real                                   :: old_can_prss
     real                                   :: old_can_enthalpy
     real                                   :: old_can_temp
     real                                   :: old_can_shv
     real                                   :: old_can_co2
     real                                   :: old_can_rhos
     real                                   :: old_can_dmol
     real                                   :: mid_can_rhos
     real                                   :: mid_can_dmol
     real                                   :: wtime0
     real(kind=8)                           :: hbeg
     integer                                :: ibuff
     integer                                :: npa_thread
     integer                                :: ita
     !----- Local constants. -----------------------------------------------!
     logical                  , parameter   :: test_energy_sanity = .false.
     !----- External functions. --------------------------------------------!
     real, external                         :: walltime
     !----------------------------------------------------------------------!

     polyloop: do ipy = 1,cgrid%npolygons
        cpoly => cgrid%polygon(ipy)

        wtime0=walltime(0.)

        siteloop: do isi = 1,cpoly%nsites
           csite => cpoly%site(isi)
           cmet  => cpoly%met(isi)

           !----- Find the number of patches per thread. ---------------------------------!
           npa_thread = ceiling(real(csite%npatches) / real(nthreads))
           !------------------------------------------------------------------------------!

           !---------------------------------------------------------------------!
           !     Update the monthly rainfall.                                    !
           !---------------------------------------------------------------------!
           imon                             = current_time%month
           cpoly%avg_monthly_pcpg(imon,isi) = cpoly%avg_monthly_pcpg(imon,isi)   &
                                            + cmet%pcpg * dtlsm
           !---------------------------------------------------------------------!

           call copy_met_2_rk4site(nzg,cmet%atm_ustar,cmet%atm_theiv         &
                ,cmet%atm_vpdef      &
                ,cmet%atm_theta,cmet%atm_tmp,cmet%atm_shv   &
                ,cmet%atm_co2,cmet%geoht,cmet%exner         &
                ,cmet%pcpg,cmet%qpcpg,cmet%dpcpg,cmet%prss  &
                ,cmet%rshort,cmet%rlong,cmet%par_beam       &
                ,cmet%par_diffuse,cmet%nir_beam             &
                ,cmet%nir_diffuse,cmet%geoht                &
                ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)     &
                ,cpoly%green_leaf_factor(:,isi)             &
                ,cgrid%lon(ipy),cgrid%lat(ipy)              &
                ,cgrid%cosz(ipy))




           !------------------------------------------------------------------------------!
           !  MLO - Changed the parallel do loop to account for cases in which the number !
           !        of threads is less than the number of patches.                        !
           !------------------------------------------------------------------------------!
           !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(          &
           !$OMP  ita,ipa,initp,ytemp,dinitp,yprev,hbeg,nsteps &
           !$OMP ,patch_vels,old_can_prss,old_can_enthalpy     &
           !$OMP ,old_can_temp,old_can_shv,old_can_co2         &
           !$OMP ,old_can_rhos,old_can_dmol,mid_can_rhos       &
           !$OMP ,mid_can_dmol,wcurr_loss2atm,ecurr_netrad     &
           !$OMP ,ecurr_loss2atm,co2curr_loss2atm              &
           !$OMP ,wcurr_loss2drainage,ecurr_loss2drainage      &
           !$OMP ,wcurr_loss2runoff,ecurr_loss2runoff)
           threadloop: do ibuff=1,nthreads
              initp => integration_buff(ibuff)%initp
              ytemp => integration_buff(ibuff)%ytemp
              dinitp => integration_buff(ibuff)%dinitp
              yprev  => integration_buff(ibuff)%yprev



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

                 !----- Reset all buffers to zero, as a safety measure. ------------!
                 call zero_rk4_patch(initp)
                 call zero_rk4_patch(ytemp)
                 call zero_rk4_patch(dinitp)
                 call zero_bdf2_patch(yprev)
                 call zero_rk4_cohort(initp)
                 call zero_rk4_cohort(ytemp)
                 call zero_rk4_cohort(dinitp)

                 !----- Get velocity for aerodynamic resistance. -------------------!
                 if (csite%can_theta(ipa) < cmet%atm_theta) then
                    patch_vels = cmet%vels_stab
                    cmet%vels  = cmet%vels_stab
                 else
                    patch_vels = cmet%vels_unstab
                    cmet%vels  = cmet%vels_unstab
                 end if

                 !------------------------------------------------------------------!


                 !----- Save the previous thermodynamic state. ---------------------!
                 old_can_prss     = csite%can_prss(ipa)
                 old_can_enthalpy = tq2enthalpy(csite%can_temp(ipa)                 &
                                               ,csite%can_shv (ipa),.true.)
                 old_can_temp     = csite%can_temp(ipa)
                 old_can_shv      = csite%can_shv (ipa)
                 old_can_co2      = csite%can_co2 (ipa)
                 old_can_rhos     = csite%can_rhos(ipa)
                 old_can_dmol     = csite%can_dmol(ipa)
                 !------------------------------------------------------------------!



                 !------------------------------------------------------------------!
                 !      Test whether temperature and energy are reasonable.         !
                 !------------------------------------------------------------------!
                 if (test_energy_sanity) then
                    call sanity_check_veg_energy(csite,ipa)
                 end if
                 !------------------------------------------------------------------!

                 !------------------------------------------------------------------!
                 !     Set up the buffer for the previous step's leaf temperature   !
                 !------------------------------------------------------------------!
                 call copy_bdf2_prev(csite,ipa,yprev)

                 !------------------------------------------------------------------!
                 !   Get plant water flow driven by plant hydraulics.  This must be !
                 ! placed before canopy_photosynthesis because plant_hydro_driver   !
                 ! needs fs_open from previous timestep.                            !
                 !------------------------------------------------------------------!
                 call plant_hydro_driver(csite,ipa,cpoly%ntext_soil(:,isi))
                 !------------------------------------------------------------------!


                 !----- Get photosynthesis, stomatal conductance,
                 !                                    and transpiration. -----------!
                 call canopy_photosynthesis(csite,cmet,nzg,ipa,ibuff,               &
                      cpoly%ntext_soil(:,isi),cpoly%leaf_aging_factor(:,isi),       &
                      cpoly%green_leaf_factor(:,isi))

                 !----- Compute root and heterotrophic respiration. ----------------!
                 call soil_respiration_driver(csite,ipa,nzg,cpoly%ntext_soil(:,isi))


                 !----- Update the committed carbon change pool. -------------------!
                 call update_cbudget_committed(csite,ipa)
                 !------------------------------------------------------------------!

                 !------------------------------------------------------------------!
                 !     Set up the integration patch.                                !
                 !------------------------------------------------------------------!
                 call copy_patch_init(csite,ipa,ibuff,initp,patch_vels,mid_can_rhos &
                                     ,mid_can_dmol)
                 !------------------------------------------------------------------!

                 !------------------------------------------------------------------!
                 !  Perform the forward and backward step.  It is possible this will!
                 !  be done over a series of sub-steps.  1)derivs,2)forward,3)back  !
                 !  4) check stability and error 5) repeat as shorter or continue   !
                 !------------------------------------------------------------------!
   !               call integrate_patch_hybrid(csite,                                 &
   !                    integration_buff(ibuff)%yprev,integration_buff(ibuff)%initp,             &
   !                    integration_buff(ibuff)%dinitp,integration_buff(ibuff)%ytemp,            &
   !                    ipa,wcurr_loss2atm,                                           &
   !                    ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage,          &
   !                    ecurr_loss2drainage,wcurr_loss2runoff,                        &
   !                    ecurr_loss2runoff,ecurr_netrad,nsteps)




                  !--------------------------------------------------------------------------!
                  ! Initial step size.  Experience has shown that giving this too large a    !
                  ! value causes the integrator to fail (e.g., soil layers become            !
                  ! supersaturated).                                                         !
                  !--------------------------------------------------------------------------!
                  hbeg = dble(csite%htry(ipa))

                  !--------------------------------------------------------------------------!
                  ! Zero the canopy-atmosphere flux values.  These values are updated        !
                  ! every dtlsm, so they must be zeroed at each call.                        !
                  !--------------------------------------------------------------------------!
                  initp%upwp = 0.d0
                  initp%tpwp = 0.d0
                  initp%qpwp = 0.d0
                  initp%cpwp = 0.d0
                  initp%wpwp = 0.d0

                  !----- Go into the ODE integrator using Euler. ----------------------------!

                  call hybrid_integ(hbeg,csite,yprev,initp,dinitp,                         &
                       ytemp,ipa,isi,ibuff,nsteps)

                  !--------------------------------------------------------------------------!
                  !  Normalize canopy-atmosphere flux values.  These values are updated ever !
                  ! dtlsm, so they must be normalized every time.                            !
                  !--------------------------------------------------------------------------!
                  initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
                  initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
                  initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
                  initp%cpwp = initp%can_dmol * initp%cpwp * dtrk4i
                  initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i

                  !--------------------------------------------------------------------------!
                  ! Move the state variables from the integrated patch to the model patch.   !
                  !--------------------------------------------------------------------------!
                  call initp2modelp(tend-tbeg,initp,csite,ipa,cpoly%nighttime(isi)   &
                                   ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm       &
                                   ,co2curr_loss2atm,wcurr_loss2drainage             &
                                   ,ecurr_loss2drainage,wcurr_loss2runoff            &
                                   ,ecurr_loss2runoff)


                  !----- Add the number of steps into the step counter. -------------!
                  cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)

                  !------------------------------------------------------------------!
                  !    Update the minimum monthly temperature,                       !
                  !    based on canopy temperature.                                  !
                  !------------------------------------------------------------------!
                  if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                     cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
                  end if
                 !------------------------------------------------------------------!



                 !------------------------------------------------------------------!
                 !    Update roughness and canopy depth.  This should be done after !
                 ! the integration.                                                 !
                 !------------------------------------------------------------------!
                 call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs,&
                      cpoly%ntext_soil(:,isi))
                 call update_patch_derived_props(csite,ipa,.false.)
                 !------------------------------------------------------------------!

                  !------------------------------------------------------------------!
                  !     Compute the residuals.                                       !
                  !------------------------------------------------------------------!

                  call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa  &
                       ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm                   &
                       ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage     &
                       ,wcurr_loss2runoff,ecurr_loss2runoff,cpoly%area(isi)          &
                       ,cgrid%cbudget_nep(ipy),old_can_prss,old_can_enthalpy         &
                       ,old_can_temp,old_can_shv,old_can_co2,old_can_rhos            &
                       ,old_can_dmol,mid_can_rhos,mid_can_dmol)
               end do taskloop
               !---------------------------------------------------------------------------!
            end do threadloop
            !$OMP END PARALLEL DO
            !------------------------------------------------------------------------------!



            !-------------------------------------------------------------------!
         end do siteloop
         !----------------------------------------------------------------------!

         cgrid%walltime_py(ipy) = cgrid%walltime_py(ipy)+walltime(wtime0)

      end do polyloop

      return
    end subroutine hybrid_timestep
    !============================================================================!
    !============================================================================!


    !============================================================================!
    !============================================================================!
    !  This subroutine will drive the integration of several ODEs that drive     !
    !  the fast-scale state variables.                                           !
    !----------------------------------------------------------------------------!
    subroutine hybrid_integ(h1,csite,yprev,initp,dinitp,ytemp,ipa,isi,ibuff      &
                           ,nsteps)
      use ed_state_vars  , only : sitetype                  & ! structure
                                , patchtype                 ! ! structure
      use rk4_coms       , only : integration_vars          & ! structure
                                , rk4patchtype              & ! structure
                                , bdf2patchtype             & ! structure
                                , rk4site                   & ! intent(in)
                                , print_diags               & ! intent(in)
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
                                , safety                    & ! intent(in)
                                , pgrow                     & ! intent(in)
                                , pshrnk                    & ! intent(in)
                                , print_detailed            & ! intent(in)
                                , norm_rk4_fluxes           & ! sub-routine
                                , reset_rk4_fluxes          ! ! sub-routine
      use rk4_integ_utils, only : rk4_sanity_check          & ! subroutine
                                , inc_rk4_patch             & ! subroutine
                                , print_sanity_check        ! ! subroutine
      use rk4_derivs     , only : leaf_derivs               ! ! subroutine
      use rk4_copy_patch , only : copy_rk4_patch            ! ! sub-routine
      use rk4_misc       , only : print_rk4patch            & ! sub-routine
                                , adjust_veg_properties     & ! sub-routine
                                , adjust_topsoil_properties & ! sub-routine
                                , adjust_sfcw_properties    & ! sub-routine
                                , update_diagnostic_vars    & ! sub-routine
                                , print_rk4_state           ! ! sub-routine
      use ed_misc_coms   , only : fast_diagnostics          & ! intent(in)
                                , dtlsm                     ! ! intent(in)
      use grid_coms      , only : nzg                       & ! intent(in)
                                , nzs                       & ! intent(in)
                                , time                      ! ! intent(in)
      use soil_coms      , only : runoff_time_i             & ! intent(in)
                                , simplerunoff              ! ! intent(in)
      use consts_coms    , only : cliq8                     & ! intent(in)
                                , t3ple8                    & ! intent(in)
                                , tsupercool_liq8               & ! intent(in)
                                , wdnsi8                    ! ! intent(in)
      use bdf2_solver    , only : bdf2_integ                ! ! sub-routine
      implicit none
      !----- Arguments ----------------------------------------------------------!
      type(sitetype)            , target      :: csite   ! Current site

      type(rk4patchtype)        , target      :: initp   ! Current integ. patch
      type(rk4patchtype)        , target      :: dinitp  ! Integration derivative
      type(rk4patchtype)        , target      :: ytemp   ! Patch at n+1
      type(bdf2patchtype)       , target      :: yprev   ! Patch at n-1

      integer                   , intent(in)  :: ipa     ! Current patch ID
      integer                   , intent(in)  :: isi     ! Current site ID
      integer                   , intent(in)  :: ibuff   ! Multithread ID
      real(kind=8)              , intent(in)  :: h1      ! First guess of delta-t
      integer                   , intent(out) :: nsteps  ! Number of steps taken.
      !----- Local variables ----------------------------------------------------!
      type(patchtype)           , pointer     :: cpatch           ! Current patch
      logical                                 :: restart_step
      logical                                 :: reject_step
      logical                                 :: minstep
      logical                                 :: stuck  
      logical                                 :: test_reject 
      logical                                 :: gapstep
      integer                                 :: i         
      integer                                 :: k            ! Format counter
      integer                                 :: ksn          ! # of snow/water
                                                              ! layers
      real(kind=8)                            :: x            ! Elapsed time
      real(kind=8)                            :: xnew         ! Elapsed time + h
      real(kind=8)                            :: newh         ! New time step
                                                              ! suggested
      real(kind=8)                            :: oldh         ! Old time step
      real(kind=8)                            :: h            ! Current delta-t
                                                              ! attempt
      real(kind=8)                            :: htrunc
      real(kind=8)                            :: fgrow        ! Delta-t increase factor
      real(kind=8)                            :: hgoal        ! Delta-t ignoring overstep
      real(kind=8)                            :: hnext        ! Next delta-t
      real(kind=8)                            :: qwfree       ! Free water
                                                              ! internal energy
      real(kind=8)                            :: wfreeb       ! Free water
      real(kind=8)                            :: errmax       ! Maximum error
                                                              ! of this step
      real(kind=8)                            :: elaptime     ! Absolute elapsed
                                                              ! time.

      !----- External function. -------------------------------------------------!
      real                      , external    :: sngloff
      !--------------------------------------------------------------------------!

      !----- Use some aliases for simplicity. -----------------------------------!
      cpatch => csite%patch(ipa)

      !--------------------------------------------------------------------------!
      ! Set initial time and stepsize.                                           !
      !--------------------------------------------------------------------------!
      x = tbeg
      h = h1
      if (dtrk4 < 0.d0) h = -h1

      !----- Define total elapsed time. -----------------------------------------!
      elaptime = time + x

      !--------------------------------------------------------------------------!
      ! Begin timestep loop                                                      !
      !--------------------------------------------------------------------------!
      timesteploop: do i=1,maxstp

         !----- Be sure not to overstep -----------------------------------------!
         hgoal   = h
         gapstep = (x+h-tend)*(x+h-tbeg) > 0.d0
         if (gapstep) h = tend - x

         reject_step =  .false.
         hstep:   do

            call leaf_derivs(initp,dinitp,csite,ipa,ibuff,h,.true.)

            !---------------------------------------------------------------------!
            ! Very simple analysis of derivative.  ie try to reduce drastic
            ! changes in key state variables.
            !---------------------------------------------------------------------!
            call fb_dy_step_trunc(initp,restart_step,dinitp,h,htrunc)

            if (restart_step) then

               oldh    = h
               newh    = htrunc
               minstep = (newh == h) .or. newh < hmin

               if(minstep)then

                  call fail_whale()
                  write(*,*) "hybrid euler truncation converged"
                  print*,htrunc,h
                  stop
               end if


               !------------------------------------------------------------------------------!
               !      We are about to shrink the time step, so this can't be considered gap   !
               ! step.                                                                        !
               !------------------------------------------------------------------------------!
               gapstep = .false.
               !------------------------------------------------------------------------------!



               !----- Define new step and check whether it is sufficiently long or not. ------!
               h       = max(1.d-1*h, newh)
               hgoal   = h
               xnew    = x + h
               stuck   = xnew == x
               !------------------------------------------------------------------------------!


               cycle hstep
            end if


            !---------------------------------------------------------------------------------!
            !  MLO -> RGK.                                                                    !
            !    I did not see any reason to have separate sub-routines for having separate   !
            ! sub-routines for copying and integrating the hybrid step.  I replaced them with !
            ! the rk4 ones, to reduce the number of points to change the code (I already      !
            ! noticed a few additions missing in the hybrid step counterparts).               !
            !---------------------------------------------------------------------------------!
            !---- Copy patch to the temporary structure. -------------------------------------!
            call copy_rk4_patch(initp,ytemp,cpatch)
            !---- Integrate the forward step. ------------------------------------------------!
            call inc_rk4_patch(ytemp,dinitp,h,cpatch)
            !---------------------------------------------------------------------------------!

            !--------------------------------------------------------------------!
            !   Integrate the implicit/backwards step                            !
            !--------------------------------------------------------------------!
            call bdf2_integ(ibuff,cpatch,yprev,initp,ytemp,dinitp,h, &
                 dble(csite%hprev(ipa)))

            !----- Perform a sanity check on canopy,leaf and wood stuff ---------!
            call rk4_sanity_check(ibuff,ytemp,reject_step,csite,ipa,dinitp,h,print_diags)

            !---------------------------------------------------------------------------------!
            !     Here we check the error of this step.  Three outcomes are possible:         !
            ! 1.  The updated values make no sense.  Reject step, assign a large error and    !
            !     try again with a smaller time step;                                         !
            ! 2.  The updated values are reasonable, but the error is large.  Reject step and !
            !     try again with a smaller time step;                                         !
            ! 3.  The updated values are reasonable, and the error is small.  Accept step and !
            !     try again with a larger time step.                                          !
            !---------------------------------------------------------------------------------!
            if (reject_step) then
               !------------------------------------------------------------------------------!
               !    If step was already rejected, that means the step had finished premature- !
               ! ly, so we assign a standard large error (10.0).                              !
               !------------------------------------------------------------------------------!
               errmax = 1.d1

            else

               errmax = 1.d-1

            end if

            !---------------------------------------------------------------------------------!
            ! 3. If the step failed, then calculate a new shorter step size to try.           !
            !---------------------------------------------------------------------------------!

            if (reject_step) then
               !------------------------------------------------------------------------------!
               !      We are about to shrink the time step, so this can't be considered gap   !
               ! step.                                                                        !
               !------------------------------------------------------------------------------!
               gapstep = .false.
               !------------------------------------------------------------------------------!


               !----- Define new step and check whether it is sufficiently long or not. ------!
               oldh    = h
               newh    = safety * h * errmax**pshrnk
               minstep = (newh == h) .or. newh < hmin
               !------------------------------------------------------------------------------!


               !----- Define next time, and check whether it is long enough to change time. --!
               h       = max(1.d-1*h, newh)
               hgoal   = h
               xnew    = x + h
               stuck   = xnew == x
               !------------------------------------------------------------------------------!



               !------------------------------------------------------------------------------!
               ! 3a. Here is the moment of truth... If we reached a tiny step and yet the     !
               !     model didn't converge, then we print various values to inform the user   !
               !     and abort the run.  Please, don't hate the messenger.                    !
               !------------------------------------------------------------------------------!
               if (minstep .or. stuck ) then 

                  write (unit=*,fmt='(80a)')         ('=',k=1,80)
                  write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN EULER_INT'
                  write (unit=*,fmt='(80a)')         ('-',k=1,80)
                  write (unit=*,fmt='(a,1x,f11.6)')   ' + LONGITUDE:     ',rk4site%lon
                  write (unit=*,fmt='(a,1x,f11.6)')   ' + LATITUDE:      ',rk4site%lat
                  write (unit=*,fmt='(a)')           ' + PATCH INFO:    '
                  write (unit=*,fmt='(a,1x,i6)')     '   - NUMBER:      ',ipa
                  write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:         ',csite%age(ipa)
                  write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:   ',csite%dist_type(ipa)
                  write (unit=*,fmt='(a,1x,l1)')     ' + MINSTEP:       ',minstep
                  write (unit=*,fmt='(a,1x,l1)')     ' + STUCK:         ',stuck
                  write (unit=*,fmt='(a,1x,es12.4)') ' + ERRMAX:        ',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' + X:             ',x
                  write (unit=*,fmt='(a,1x,es12.4)') ' + H:             ',h
                  write (unit=*,fmt='(a,1x,es12.4)') ' + OLDH:          ',oldh
                  write (unit=*,fmt='(a,1x,es12.4)') ' + NEWH:          ',newh
                  write (unit=*,fmt='(a,1x,es12.4)') ' + SAFETY:        ',safety
                  write (unit=*,fmt='(80a)') ('-',k=1,80)
                  write (unit=*,fmt='(a)') '   Likely to be a rejected step problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)

                  call rk4_sanity_check(ibuff,ytemp,test_reject,csite,ipa,dinitp,h,.true.)
                  call print_sanity_check(ytemp,csite,ipa)
                  call print_rk4patch(ytemp, csite,ipa)
               end if

            else
               !------------------------------------------------------------------------------!
               ! 3b.  Great, it worked, so now we can advance to the next step.  We just need !
               !      to do some minor adjustments before...                                  !
               !------------------------------------------------------------------------------!
               call adjust_veg_properties(ytemp,h,csite,ipa,ibuff)

               !----- ii.  Final update of top soil properties to avoid off-bounds moisture. -!
               call adjust_topsoil_properties(ytemp,h,ibuff)

               !----- ii. Make temporary surface water stable and positively defined. --------!
               call adjust_sfcw_properties(nzg,nzs,ytemp,h,ibuff)

               !----- iii.  Update the diagnostic variables. ---------------------------------!
               call update_diagnostic_vars(ytemp, csite,ipa,ibuff)
               !------------------------------------------------------------------------------!

               !------------------------------------------------------------------------------!
               ! 3c. Set up h for the next time.  And here we can relax h for the next step,  !
               !    and try something faster, unless this is a "gap step" (shorter time step  !
               !    just to close the full thermodynamic time step).                          !
               !------------------------------------------------------------------------------!
               if (gapstep) then
                  hnext = hgoal
               else
                  fgrow = min(5.d0,1.d0+sqrt(2.d0),safety*errmax**pgrow)
                  hnext = max(2.d0*hmin, min(dble(dtlsm), fgrow * hgoal))
               end if
               !------------------------------------------------------------------------------!


               !------ 3d. Normalise the fluxes if the user wants detailed debugging. --------!
               if (print_detailed) then
                  call norm_rk4_fluxes(ytemp,h)
                  call print_rk4_state(initp,ytemp,csite,ipa,isi,x,h)
               end if

               !----- 3e. Copy the temporary structure to the intermediate state. ------------!
               call copy_initp2prev(initp,yprev,csite%patch(ipa))

               call copy_rk4_patch(ytemp, initp,csite%patch(ipa))

               !------------------------------------------------------------------------------!
               !    3f. Flush step-by-step fluxes to zero if the user wants detailed          !
               !        debugging.                                                            !
               !------------------------------------------------------------------------------!
               if (print_detailed) then
                  call reset_rk4_fluxes(initp)
               end if

               !----- 3g. Update time. -------------------------------------------------------!
               csite%hprev(ipa) = sngl(h)
               x        = x + h
               h        = hnext
               elaptime = elaptime + h

               exit hstep
            end if
         end do hstep


         !----- If the integration reached the next step, make some final adjustments --------!
         if((x-tend)*dtrk4 >= 0.d0)then

            ksn = initp%nlev_sfcwater

            !---------------------------------------------------------------------------------!
            !   Make temporary surface liquid water disappear.  This will not happen          !
            ! immediately, but liquid water will decay with the time scale defined by         !
            ! runoff_time scale. If the time scale is too tiny, then it will be forced to be  !
            ! hdid (no reason to be faster than that).                                        !
            !---------------------------------------------------------------------------------!
            if (simplerunoff .and. ksn >= 1) then

               if (initp%sfcwater_mass(ksn)    > 0.d0           .and.                         &
                   initp%sfcwater_fracliq(ksn) > 1.d-1        ) then

                  wfreeb = min(1.d0,dtrk4*runoff_time_i) * initp%sfcwater_mass(ksn)           &
                         * (initp%sfcwater_fracliq(ksn) - 1.d-1) / 9.d-1

                  qwfree = wfreeb * cliq8 * (initp%sfcwater_tempk(ksn) - tsupercool_liq8 )

                  initp%sfcwater_mass(ksn) = initp%sfcwater_mass(ksn)   - wfreeb

                  initp%sfcwater_depth(ksn) = initp%sfcwater_depth(ksn) - wfreeb * wdnsi8

                  !----- Recompute the energy removing runoff --------------------------------!
                  initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn) - qwfree

                  call adjust_sfcw_properties(nzg,nzs,initp,dtrk4,ibuff)
                  call update_diagnostic_vars(initp,csite,ipa,ibuff)

                  !----- Compute runoff for output -------------------------------------------!
                  if (fast_diagnostics) then
                     csite%runoff       (ipa) = csite%runoff(ipa)                             &
                                              + sngloff(wfreeb * dtrk4i,tiny_offset)
                     csite%fmean_runoff (ipa) = csite%fmean_runoff(ipa)                       &
                                              + sngloff(wfreeb * dtrk4i,tiny_offset)
                     csite%fmean_qrunoff(ipa) = csite%fmean_qrunoff(ipa)                      &
                                              + sngloff(qwfree * dtrk4i,tiny_offset)
                  end if
                  if (checkbudget) then
                     initp%wbudget_loss2runoff = initp%wbudget_loss2runoff + wfreeb
                     initp%ebudget_loss2runoff = initp%ebudget_loss2runoff + qwfree
                     initp%wbudget_storage     = initp%wbudget_storage - wfreeb
                     initp%ebudget_storage     = initp%ebudget_storage - qwfree
                  end if
               end if
            end if

            !------ Update the substep for next time and leave -------------------------------!

            csite%htry(ipa)  = sngl(hnext)

            call copy_prev2patch(yprev,csite,ipa)

            !---------------------------------------------------------------------------------!
            !     Update the average time step.  The square of DTLSM (tend-tbeg) is needed    !
            ! because we will divide this by the time between t0 and t0+frqsum.               !
            !---------------------------------------------------------------------------------!
            csite%fmean_rk4step(ipa) = csite%fmean_rk4step(ipa)                               &
                 + sngl((tend-tbeg)*(tend-tbeg))/real(i)
            nsteps = i
            return
         end if

         !----- Use hnext as the next substep ------------------------------------------------!
         h = hnext
      end do timesteploop

      !----- If it reached this point, that is really bad news... ----------------------------!
      write (unit=*,fmt='(a)') ' ==> Too many steps in routine euler_integ'
      call print_rk4patch(ytemp, csite,ipa)

      return
    end subroutine hybrid_integ


    !=========================================================================================!
    !=========================================================================================!



    !=============================================================!

    subroutine copy_initp2prev(initp,yprev,cpatch)

      use rk4_coms             , only : rk4patchtype,bdf2patchtype
      use ed_state_vars        , only : patchtype

      implicit none

      type(rk4patchtype), target     :: initp      ! Main memory
      type(bdf2patchtype), target    :: yprev      ! Buffer memory
      type(patchtype),target         :: cpatch
      integer                        :: ico

      yprev%can_temp = initp%can_temp

      do ico=1,cpatch%ncohorts
         yprev%leaf_temp(ico) = initp%leaf_temp(ico)
         yprev%wood_temp(ico) = initp%wood_temp(ico)
      end do

      return
    end subroutine copy_initp2prev

    subroutine copy_prev2patch(yprev,csite,ipa)

      use rk4_coms             , only : bdf2patchtype  & ! structure
                                      , tiny_offset    ! ! intent(in)
      use ed_state_vars        , only : patchtype,sitetype

      implicit none

      type(bdf2patchtype), target    :: yprev      ! Buffer memory
      type(patchtype),pointer         :: cpatch
      type(sitetype),target          :: csite
      integer                        :: ico,ipa
      real, external :: sngloff

      cpatch => csite%patch(ipa)
      csite%can_temp_pv(ipa) = sngloff(yprev%can_temp,tiny_offset)

      do ico=1,cpatch%ncohorts
         cpatch%leaf_temp_pv(ico) = sngloff(yprev%leaf_temp(ico),tiny_offset)
         cpatch%wood_temp_pv(ico) = sngloff(yprev%wood_temp(ico),tiny_offset)
      end do

      return
    end subroutine copy_prev2patch



    !=============================================================!

    subroutine copy_bdf2_prev(csite,ipa,yprev)

      use ed_state_vars        , only : patchtype,sitetype   ! ! structure
      use rk4_coms             , only : bdf2patchtype          ! structure

      implicit none

      type(sitetype)    , target     :: csite
      type(patchtype)   , pointer    :: cpatch     ! Main memory
      type(bdf2patchtype), target     :: yprev      ! Buffer memory
      integer                        :: ipa
      integer                        :: ico

      cpatch => csite%patch(ipa)
      yprev%can_temp = csite%can_temp_pv(ipa)

      do ico=1,cpatch%ncohorts
         yprev%leaf_temp(ico) = cpatch%leaf_temp_pv(ico)
         yprev%wood_temp(ico) = cpatch%wood_temp_pv(ico)
      end do

      return
    end subroutine copy_bdf2_prev

   !=============================================================!



    ! ========================================================================= !

    subroutine fb_dy_step_trunc(y,restart_step,dydx,h,hmin)

      use rk4_coms               , only : rk4patchtype          & ! structure
                                        , integration_vars      & ! structure
                                        , integ_err             & ! intent(inout)
                                        , record_err            ! ! intent(in)
      use therm_lib8             , only : eslif8
      use consts_coms            , only : ep8

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target      :: y
      type(rk4patchtype) , target      :: dydx
      real(kind=8)       , intent(in)  :: h
      real(kind=8)       , intent(out) :: hmin
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                     :: max_dshv_can
      logical                          :: restart_step
      !------------------------------------------------------------------------------------!

      restart_step = .false.

      ! ---------------- Maximum step change in canopy CO2 (PPM) --------------------------!

   !!   max_dco2_can = 20.d0
   !!   hmin = max_dco2_can/abs(dydx%can_co2)

   !!   if ( h > max_dco2_can/abs(dydx%can_co2) .and. record_err) &
   !!        integ_err(6,1) = integ_err(6,1) + 1_8

      ! ---------------- Maximum change in canopy Relative Humidity ----------------------!

      max_dshv_can = 0.15d0
      hmin = (max_dshv_can*y%can_shv)/max(abs(dydx%can_shv),1.0d-10)

      if ( h > max_dshv_can/max(abs(dydx%can_shv),1.0d-10) .and. record_err) then
         integ_err(3,1) = integ_err(3,1) + 1_8
      end if

      ! ---------------- Maximum change in the first soil layer --------------------------!

   !!   max_dwater_soil = 0.5d0  ! Maximum change in relative soil moisture

   !!   do k=rk4site%lsl,nzg

   !!      hmin_tmp = max_dwater_soil/(abs(dydx%soil_water(k))/soil8(rk4site%ntext_soil(k))%slmsts)
   !!      hmin = min(hmin,hmin_tmp)

   !!      if ( h > hmin_tmp .and. record_err) &
   !!           integ_err(osow+k,1) = integ_err(osow+k,1) + 1_8

   !!   end do

      if (hmin < 0.99999*h) then
         restart_step = .true.
      else
         restart_step = .false.
      end if

      return
    end subroutine fb_dy_step_trunc

    !===========================================================
end module hybrid_driver

