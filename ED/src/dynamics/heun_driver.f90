!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main driver for the Heun's (RK2) integration scheme.          !
!------------------------------------------------------------------------------------------!
subroutine heun_timestep(cgrid)
   use rk4_coms              , only : integration_vars   & ! structure
                                    , rk4patchtype       & ! structure
                                    , zero_rk4_patch     & ! subroutine
                                    , zero_rk4_cohort    & ! subroutine
                                    , integration_buff   & ! intent(out)
                                    , rk4site            ! ! intent(out)
   use ed_state_vars         , only : edtype             & ! structure
                                    , polygontype        & ! structure
                                    , sitetype           & ! structure
                                    , patchtype          ! ! structure
   use met_driver_coms       , only : met_driv_state     ! ! structure
   use grid_coms             , only : nzg                & ! intent(in)
                                    , nzs                ! ! intent(in)
   use ed_misc_coms          , only : current_time       & ! intent(in)
                                    , dtlsm              ! ! intent(in)
   use ed_max_dims           , only : n_dbh              ! ! intent(in)
   use soil_coms             , only : soil_rough         & ! intent(in)
                                    , snow_rough         ! ! intent(in)
   use therm_lib             , only : tq2enthalpy        ! ! function
   use budget_utils          , only : update_budget      & ! function
                                    , compute_budget     ! ! function
   ! OMP use omp_lib

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)             , target      :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype)        , pointer     :: cpoly
   type(sitetype)           , pointer     :: csite
   type(patchtype)          , pointer     :: cpatch
   type(met_driv_state)     , pointer     :: cmet
   integer                                :: ipy
   integer                                :: isi
   integer                                :: ipa
   integer                                :: imon
   integer                                :: nsteps
   real                                   :: wcurr_loss2atm
   real                                   :: ecurr_netrad
   real                                   :: ecurr_loss2atm
   real                                   :: co2curr_loss2atm
   real                                   :: wcurr_loss2drainage
   real                                   :: ecurr_loss2drainage
   real                                   :: wcurr_loss2runoff
   real                                   :: ecurr_loss2runoff
   real                                   :: old_can_enthalpy
   real                                   :: old_can_shv
   real                                   :: old_can_co2
   real                                   :: old_can_rhos
   real                                   :: old_can_temp
   real                                   :: old_can_prss
   real                                   :: patch_vels
   integer                                :: ibuff
   !----- Local constants. ----------------------------------------------------------------!
   logical                   , parameter  :: test_energy_sanity = .false.
   !---------------------------------------------------------------------------------------!

   ibuff = 1

   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)
         cmet  => cpoly%met(isi)

         !---------------------------------------------------------------------------------!
         !     Update the monthly rainfall.                                                !
         !---------------------------------------------------------------------------------!
         imon                             = current_time%month
         cpoly%avg_monthly_pcpg(imon,isi) = cpoly%avg_monthly_pcpg(imon,isi)               &
                                          + cmet%pcpg * dtlsm
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Copy the meteorological variables to the rk4site structure.                  !
         !---------------------------------------------------------------------------------!

         call copy_met_2_rk4site(nzg,cmet%atm_ustar,cmet%atm_theiv,cmet%atm_vpdef          &
                                ,cmet%atm_theta,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2     &
                                ,cmet%geoht,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg     &
                                ,cmet%prss,cmet%rshort,cmet%rlong,cmet%par_beam            &
                                ,cmet%par_diffuse,cmet%nir_beam,cmet%nir_diffuse           &
                                ,cmet%geoht,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)         &
                                ,cpoly%green_leaf_factor(:,isi),cgrid%lon(ipy)             &
                                ,cgrid%lat(ipy),cgrid%cosz(ipy))
         !---------------------------------------------------------------------------------!


         patchloop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)


            !----- Reset all buffers to zero, as a safety measure. ------------------------!
            call zero_rk4_patch(integration_buff(ibuff)%initp)
            call zero_rk4_patch(integration_buff(ibuff)%ytemp)
            call zero_rk4_patch(integration_buff(ibuff)%yerr)
            call zero_rk4_patch(integration_buff(ibuff)%yscal)
            call zero_rk4_patch(integration_buff(ibuff)%dydx)
            call zero_rk4_cohort(integration_buff(ibuff)%initp)
            call zero_rk4_cohort(integration_buff(ibuff)%ytemp)
            call zero_rk4_cohort(integration_buff(ibuff)%yerr)
            call zero_rk4_cohort(integration_buff(ibuff)%yscal)
            call zero_rk4_cohort(integration_buff(ibuff)%dydx)

            !----- Get velocity for aerodynamic resistance. -------------------------------!
            if (csite%can_theta(ipa) < cmet%atm_theta) then
!               cmet%vels = cmet%vels_stab
               patch_vels = cmet%vels_stab
            else
!               cmet%vels = cmet%vels_unstab
               patch_vels = cmet%vels_stab
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Update roughness and canopy depth.                                        !
            !------------------------------------------------------------------------------!
            call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs,cpoly%ntext_soil(:,isi))
            call update_patch_derived_props(csite,ipa)
            !------------------------------------------------------------------------------!



            !----- Save the previous thermodynamic state. ---------------------------------!
            old_can_shv      = csite%can_shv  (ipa)
            old_can_co2      = csite%can_co2  (ipa)
            old_can_rhos     = csite%can_rhos (ipa)
            old_can_temp     = csite%can_temp (ipa)
            old_can_prss     = csite%can_prss (ipa)
            old_can_enthalpy = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
            !------------------------------------------------------------------------------!

            !----- Compute current storage terms. -----------------------------------------!
            call update_budget(csite,cpoly%lsl(isi),ipa)
            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !      Test whether temperature and energy are reasonable.                     !
            !------------------------------------------------------------------------------!
            if (test_energy_sanity) then
               call sanity_check_veg_energy(csite,ipa)
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Set up the integration patch.                                            !
            !------------------------------------------------------------------------------!
            call copy_patch_init(csite,ipa,integration_buff(ibuff)%initp,patch_vels)
            !------------------------------------------------------------------------------!

            !----- Get photosynthesis, stomatal conductance, and transpiration. -----------!
            call canopy_photosynthesis(csite,cmet,nzg,ipa                                  &
                                      ,cpoly%ntext_soil(:,isi)                             &
                                      ,cpoly%leaf_aging_factor(:,isi)                      &
                                      ,cpoly%green_leaf_factor(:,isi))
            !------------------------------------------------------------------------------!



            !----- Compute root and heterotrophic respiration. ----------------------------!
            call soil_respiration(csite,ipa,nzg,cpoly%ntext_soil(:,isi))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Set up the remaining, carbon-dependent variables to the buffer.          !
            !------------------------------------------------------------------------------!
            call copy_patch_init_carbon(csite,ipa,integration_buff(ibuff)%initp)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     This is the step in which the derivatives are computed, we use a         !
            ! structure that is very similar to the Runge-Kutta, though a simpler one.     !
            !------------------------------------------------------------------------------!
            call integrate_patch_heun(csite,ipa,isi,cpoly%nighttime(isi),wcurr_loss2atm    &
                                     ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm         &
                                     ,wcurr_loss2drainage,ecurr_loss2drainage              &
                                     ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps)
            !------------------------------------------------------------------------------!



            !----- Add the number of steps into the step counter. -------------------------!
            cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Update the minimum monthly temperature, based on canopy temperature.      !
            !------------------------------------------------------------------------------!
            if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
               cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Compute the residuals.                                                   !
            !------------------------------------------------------------------------------!
            call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa              &
                               ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm                 &
                               ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage   &
                               ,wcurr_loss2runoff,ecurr_loss2runoff,cpoly%area(isi)        &
                               ,cgrid%cbudget_nep(ipy),old_can_enthalpy,old_can_shv        &
                               ,old_can_co2,old_can_rhos,old_can_prss)
            !------------------------------------------------------------------------------!
         end do patchloop
         !---------------------------------------------------------------------------------!
      end do siteloop
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine heun_timestep
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the integration process using the Heun method.  Notice    !
! that most of the Heun method utilises the subroutines from Runge-Kutta.                  !
!------------------------------------------------------------------------------------------!
subroutine integrate_patch_heun(csite,ipa,isi,nighttime,wcurr_loss2atm,ecurr_netrad        &
                               ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage        &
                               ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff    &
                               ,nsteps)
   use ed_state_vars   , only : sitetype             & ! structure
                              , patchtype            ! ! structure
   use ed_misc_coms    , only : dtlsm                ! ! intent(in)
   use soil_coms       , only : soil_rough           & ! intent(in)
                              , snow_rough           ! ! intent(in)
   use canopy_air_coms , only : exar8                ! ! intent(in)
   use rk4_coms        , only : integration_vars     & ! structure
                              , integration_buff     & ! structure
                              , rk4site              & ! intent(inout)
                              , zero_rk4_patch       & ! subroutine
                              , zero_rk4_cohort      & ! subroutine
                              , tbeg                 & ! intent(inout)
                              , tend                 & ! intent(inout)
                              , dtrk4                & ! intent(inout)
                              , dtrk4i               ! ! intent(inout)
   use rk4_driver      , only : initp2modelp         ! ! subroutine
   ! OMP use omp_lib

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target      :: csite
   integer               , intent(in)  :: ipa
   integer               , intent(in)  :: isi
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
   integer                             :: ibuff
   !----- Local variables -----------------------------------------------------------------!
   real(kind=8)                        :: hbeg
   !---------------------------------------------------------------------------------------!

   ibuff = 1

   !---------------------------------------------------------------------------------------!
   !      Initial step size.  Experience has shown that giving this too large a value      !
   ! causes the integrator to fail (e.g., soil layers become supersaturated).              !
   !---------------------------------------------------------------------------------------!
   hbeg = dble(csite%htry(ipa))

   !---------------------------------------------------------------------------------------!
   !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm,    !
   ! so they must be zeroed at each call.                                                  !
   !---------------------------------------------------------------------------------------!
   integration_buff(ibuff)%initp%upwp = 0.d0
   integration_buff(ibuff)%initp%tpwp = 0.d0
   integration_buff(ibuff)%initp%qpwp = 0.d0
   integration_buff(ibuff)%initp%cpwp = 0.d0
   integration_buff(ibuff)%initp%wpwp = 0.d0

   !----- Go into the ODE integrator using Euler. -----------------------------------------!
   call heun_integ(hbeg,csite,ipa,isi,nsteps)

   !---------------------------------------------------------------------------------------!
   !      Normalize canopy-atmosphere flux values.  These values are updated every         !
   ! dtlsm, so they must be normalized every time.                                         !
   !---------------------------------------------------------------------------------------!
   integration_buff(ibuff)%initp%upwp = integration_buff(ibuff)%initp%can_rhos             &
                                      * integration_buff(ibuff)%initp%upwp * dtrk4i
   integration_buff(ibuff)%initp%tpwp = integration_buff(ibuff)%initp%can_rhos             &
                                      * integration_buff(ibuff)%initp%tpwp * dtrk4i
   integration_buff(ibuff)%initp%qpwp = integration_buff(ibuff)%initp%can_rhos             &
                                      * integration_buff(ibuff)%initp%qpwp * dtrk4i
   integration_buff(ibuff)%initp%cpwp = integration_buff(ibuff)%initp%can_rhos             &
                                      * integration_buff(ibuff)%initp%cpwp * dtrk4i
   integration_buff(ibuff)%initp%wpwp = integration_buff(ibuff)%initp%can_rhos             &
                                      * integration_buff(ibuff)%initp%wpwp * dtrk4i

   !---------------------------------------------------------------------------------------!
   ! Move the state variables from the integrated patch to the model patch.                !
   !---------------------------------------------------------------------------------------!
   call initp2modelp(tend-tbeg,integration_buff(ibuff)%initp,csite,ipa,nighttime           &
                    ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm           &
                    ,wcurr_loss2drainage,ecurr_loss2drainage,wcurr_loss2runoff             &
                    ,ecurr_loss2runoff)

   return
end subroutine integrate_patch_heun
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine heun_integ                                                                    !
!                                                                                          !
!     This subroutine will drive the integration of several ODEs that drive the fast-scale !
! state variables using the Heun's method (a 2nd order Runge-Kutta).                       !
!------------------------------------------------------------------------------------------!
subroutine heun_integ(h1,csite,ipa,isi,nsteps)
   use ed_state_vars  , only : sitetype               & ! structure
                             , patchtype              & ! structure
                             , polygontype            ! ! structure
   use rk4_coms       , only : integration_vars       & ! structure
                             , integration_buff       & ! structure
                             , rk4site                & ! intent(in)
                             , print_diags            & ! intent(in)
                             , print_detailed         & ! intent(in)
                             , maxstp                 & ! intent(in)
                             , tbeg                   & ! intent(in)
                             , tend                   & ! intent(in)
                             , dtrk4                  & ! intent(in)
                             , dtrk4i                 & ! intent(in)
                             , tiny_offset            & ! intent(in)
                             , checkbudget            & ! intent(in)
                             , zero_rk4_patch         & ! subroutine
                             , zero_rk4_cohort        & ! subroutine
                             , hmin                   & ! intent(in)
                             , rk4eps                 & ! intent(in)
                             , rk4epsi                & ! intent(in)
                             , safety                 & ! intent(in)
                             , pgrow                  & ! intent(in)
                             , pshrnk                 & ! intent(in)
                             , errcon                 & ! intent(in)
                             , norm_rk4_fluxes        & ! sub-routine
                             , reset_rk4_fluxes       ! ! sub-routine
   use rk4_stepper    , only : rk4_sanity_check       & ! subroutine
                             , print_sanity_check     ! ! subroutine
   use ed_misc_coms   , only : fast_diagnostics       ! ! intent(in)
   use hydrology_coms , only : useRUNOFF              ! ! intent(in)
   use grid_coms      , only : nzg                    & ! intent(in)
                             , nzs                    & ! intent(in)
                             , time                   ! ! intent(in)
   use soil_coms      , only : dslz8                  & ! intent(in)
                             , runoff_time            & ! intent(in)
                             , runoff_time_i          & ! intent(in)
                             , simplerunoff           ! ! intent(in)
   use consts_coms    , only : t3ple8                 & ! intent(in)
                             , wdnsi8                 ! ! intent(in)
   use therm_lib8     , only : tl2uint8               ! ! intent(in)
   ! OMP use omp_lib

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite            ! Current site
   integer                   , intent(in)  :: ipa              ! Current patch ID
   integer                   , intent(in)  :: isi              ! Current site ID
   real(kind=8)              , intent(in)  :: h1               ! First guess of delta-t
   integer                   , intent(out) :: nsteps           ! Number of steps taken.
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch           ! Current patch
   logical                                 :: reject_step      ! Should I reject the step?
   logical                                 :: reject_result    ! Should I reject the result?
   logical                                 :: minstep          ! Minimum time step reached
   logical                                 :: stuck            ! Tiny step, it won't advance
   logical                                 :: test_reject      ! Reject the test
   integer                                 :: i                ! Step counter
   integer                                 :: k                ! Format counter
   integer                                 :: ksn              ! # of snow/water layers
   real(kind=8)                            :: x                ! Elapsed time
   real(kind=8)                            :: xnew             ! Elapsed time + h
   real(kind=8)                            :: newh             ! New time step suggested
   real(kind=8)                            :: oldh             ! Old time step
   real(kind=8)                            :: h                ! Current delta-t attempt
   real(kind=8)                            :: hnext            ! Next delta-t
   real(kind=8)                            :: hdid             ! delta-t that worked (???)
   real(kind=8)                            :: qwfree           ! Free water internal energy
   real(kind=8)                            :: wfreeb           ! Free water
   real(kind=8)                            :: errmax           ! Maximum error of this step
   real(kind=8)                            :: elaptime         ! Absolute elapsed time.
   integer                                 :: ibuff
   !----- External function. --------------------------------------------------------------!
   real                      , external    :: sngloff
   !---------------------------------------------------------------------------------------!

   ibuff = 1

   !----- Use some aliases for simplicity. ------------------------------------------------!
   cpatch => csite%patch(ipa)

   !---------------------------------------------------------------------------------------!
   !     Create temporary patches.                                                         !
   !---------------------------------------------------------------------------------------!
   call copy_rk4_patch(integration_buff(ibuff)%initp,integration_buff(ibuff)%y,cpatch)

   !---------------------------------------------------------------------------------------!
   ! Set initial time and stepsize.                                                        !
   !---------------------------------------------------------------------------------------!
   x = tbeg
   h = h1
   if (dtrk4 < 0.d0) h = -h1

   !----- Define total elapsed time. ------------------------------------------------------!
   elaptime = time + x

   !---------------------------------------------------------------------------------------!
   ! Begin timestep loop                                                                   !
   !---------------------------------------------------------------------------------------!
   timesteploop: do i=1,maxstp


      !----- Get initial derivatives ------------------------------------------------------!
      call leaf_derivs(integration_buff(ibuff)%y,integration_buff(ibuff)%dydx,csite,ipa,h  &
                      ,.false.)

      !----- Get scalings used to determine stability -------------------------------------!
      call get_yscal(integration_buff(ibuff)%y,integration_buff(ibuff)%dydx,h              &
                    ,integration_buff(ibuff)%yscal,cpatch)

      !----- Be sure not to overstep ------------------------------------------------------!
      if((x+h-tend)*(x+h-tbeg) > 0.d0) h=tend-x

      !------------------------------------------------------------------------------------!
      !     Here we will perform the Heun's integration using the time step.  As in Runge- !
      ! Kutta, we also check whether the integration is going well and if needed we shrink !
      ! the intermediate time steps.                                                       !
      !------------------------------------------------------------------------------------!
      reject_step =  .false.
      hstep:   do

         !---------------------------------------------------------------------------------!
         ! 1. Try a step of varying size.                                                  !
         !---------------------------------------------------------------------------------!
         call heun_stepper(x,h,csite,ipa,reject_step,reject_result)

         !---------------------------------------------------------------------------------!
         !     Here we check the error of this step.  Three outcomes are possible:         !
         ! 1.  The updated values make no sense.  Reject step, assign a large error and    !
         !     try again with a smaller time step;                                         !
         ! 2.  The updated values are reasonable, but the error is large.  Reject step and !
         !     try again with a smaller time step;                                         !
         ! 3.  The updated values are reasonable, and the error is small.  Accept step and !
         !     try again with a larger time step.                                          !
         !---------------------------------------------------------------------------------!
         if (reject_step .or. reject_result) then
            !------------------------------------------------------------------------------!
            !    If step was already rejected, that means the step had finished premature- !
            ! ly, so we assign a standard large error (10.0).                              !
            !------------------------------------------------------------------------------!
            errmax = 1.d1
         else
            call get_errmax(errmax,integration_buff(ibuff)%yerr                            &
                           ,integration_buff(ibuff)%yscal,csite%patch(ipa)                 &
                           ,integration_buff(ibuff)%initp,integration_buff(ibuff)%ytemp)
            !----- Scale the error based on the prescribed tolerance. ---------------------!
            errmax = errmax * rk4epsi
         end if


         !---------------------------------------------------------------------------------!
         ! 3. If that error was large, then calculate a new step size to try.  There are   !
         !    two types of new tries.  If step failed to be minimally reasonable (reject-  !
         !    ed) we have assigned a standard large error (10.0).  Otherwise a new step is !
         !    calculated based on the size of that error.  Hopefully, those new steps      !
         !    should be less than the previous h.  If the error was small, i.e. less than  !
         !    rk4eps, then we are done with this step, and we can move forward             !
         !    time: x = x + h                                                              !
         !---------------------------------------------------------------------------------!
         if (errmax > 1.d0) then
            !----- Defining new step and checking if it can be. ---------------------------!
            oldh    = h
            newh    = safety * h * errmax**pshrnk
            minstep = (newh == h) .or. newh < hmin

            !----- Defining next time, and checking if it really added something. ---------!
            h       = max(1.d-1*h, newh)
            xnew    = x + h
            stuck   = xnew == x

            !------------------------------------------------------------------------------!
            ! 3a. Here is the moment of truth... If we reached a tiny step and yet the     !
            !     model didn't converge, then we print various values to inform the user   !
            !     and abort the run.  Please, don't hate the messenger.                    !
            !------------------------------------------------------------------------------!
            if (minstep .or. stuck) then

               write (unit=*,fmt='(80a)')         ('=',k=1,80)
               write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN HEUN_INTEG'
               write (unit=*,fmt='(80a)')         ('-',k=1,80)
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:     ',rk4site%lon
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:      ',rk4site%lat
               write (unit=*,fmt='(a)')           ' + PATCH INFO:    '
               write (unit=*,fmt='(a,1x,i6)')     '   - NUMBER:      ',ipa
               write (unit=*,fmt='(a,1x,es12.4)') '   - AGE:         ',csite%age(ipa)
               write (unit=*,fmt='(a,1x,i6)')     '   - DIST_TYPE:   ',csite%dist_type(ipa)
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
                  !----- Run the LSM sanity check but this time we force the print. -------!
                  call rk4_sanity_check(integration_buff(ibuff)%ytemp,test_reject,csite    &
                                       ,ipa,integration_buff(ibuff)%dydx,h,.true.)
                  call print_sanity_check(integration_buff(ibuff)%y,csite,ipa)
               elseif (reject_step) then
                  call rk4_sanity_check(integration_buff(ibuff)%ak3,test_reject,csite,ipa  &
                                       ,integration_buff(ibuff)%dydx,h,.true.)
                  call print_sanity_check(integration_buff(ibuff)%y,csite,ipa)
               else
                  call print_errmax(errmax,integration_buff(ibuff)%yerr                    &
                                   ,integration_buff(ibuff)%yscal,csite%patch(ipa)         &
                                   ,integration_buff(ibuff)%y,integration_buff(ibuff)%ytemp)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax*rk4eps
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_rk4patch(integration_buff(ibuff)%y, csite,ipa)
            end if

         else
            !------------------------------------------------------------------------------!
            ! 3b.  Great, it worked, so now we can advance to the next step.  We just need !
            !      to do some minor adjustments before...                                  !
            !------------------------------------------------------------------------------!
            !----- i.   Final update of leaf properties to avoid negative water. ----------!
            call adjust_veg_properties(integration_buff(ibuff)%ytemp,h,csite,ipa)
            !----- ii.  Final update of top soil properties to avoid off-bounds moisture. -!
            call adjust_topsoil_properties(integration_buff(ibuff)%ytemp,h,csite,ipa)
            !----- iii.  Make snow layers stable and positively defined. ------------------!
            call adjust_sfcw_properties(nzg,nzs,integration_buff(ibuff)%ytemp,h,csite,ipa)
            !----- iv. Update the diagnostic variables. -----------------------------------!
            call update_diagnostic_vars(integration_buff(ibuff)%ytemp,csite,ipa)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            ! 3c. Set up h for the next time.  And here we can relax h for the next step,  !
            !    and try something faster.                                                 !
            !------------------------------------------------------------------------------!
            if (errmax > errcon) then
               hnext = safety * h * errmax**pgrow
            else
               hnext = 5.d0 * h
            endif
            hnext = max(2.d0*hmin,hnext)

            !------ 3d. Normalise the fluxes if the user wants detailed debugging. --------!
            if (print_detailed) then
               call norm_rk4_fluxes(integration_buff(ibuff)%ytemp,h)
               call print_rk4_state(integration_buff(ibuff)%y                              &
                                   ,integration_buff(ibuff)%ytemp,csite,ipa,isi,x,h)
            end if

            !----- 3e. Copy the temporary structure to the intermediate state. ------------!
            call copy_rk4_patch(integration_buff(ibuff)%ytemp,integration_buff(ibuff)%y    &
                               ,csite%patch(ipa))

            !------------------------------------------------------------------------------!
            !    3f. Flush step-by-step fluxes to zero if the user wants detailed          !
            !        debugging.                                                            !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               call reset_rk4_fluxes(integration_buff(ibuff)%y)
            end if

            !----- 3g. Update time. -------------------------------------------------------!
            x        = x + h
            hdid     = h
            h        = hnext
            elaptime = elaptime + h

            exit hstep
         end if
      end do hstep

      !----- If the integration reached the next step, make some final adjustments --------!
      if((x-tend)*dtrk4 >= 0.d0)then

         ksn = integration_buff(ibuff)%y%nlev_sfcwater

         !---------------------------------------------------------------------------------!
         !   Make temporary surface liquid water disappear.  This will not happen          !
         ! immediately, but liquid water will decay with the time scale defined by         !
         ! runoff_time scale. If the time scale is too tiny, then it will be forced to be  !
         ! hdid (no reason to be faster than that).                                        !
         !---------------------------------------------------------------------------------!
         if (simplerunoff .and. ksn >= 1) then

            if (integration_buff(ibuff)%y%sfcwater_mass(ksn)    > 0.d0  .and.              &
                integration_buff(ibuff)%y%sfcwater_fracliq(ksn) > 1.d-1 ) then

               wfreeb = min(1.d0,dtrk4*runoff_time_i)                                      &
                      * integration_buff(ibuff)%y%sfcwater_mass(ksn)                       &
                      * (integration_buff(ibuff)%y%sfcwater_fracliq(ksn) - 1.d-1) / 9.d-1
               qwfree = wfreeb                                                             &
                      * tl2uint8(integration_buff(ibuff)%y%sfcwater_tempk(ksn),1.d0)

               integration_buff(ibuff)%y%sfcwater_mass(ksn) =                              &
                                  integration_buff(ibuff)%y%sfcwater_mass(ksn)  - wfreeb
               integration_buff(ibuff)%y%sfcwater_depth(ksn) =                             &
                                  integration_buff(ibuff)%y%sfcwater_depth(ksn)            &
                                - wfreeb * wdnsi8
               !----- Recompute the energy removing runoff --------------------------------!
               integration_buff(ibuff)%y%sfcwater_energy(ksn) =                            &
                                  integration_buff(ibuff)%y%sfcwater_energy(ksn) - qwfree

               call adjust_sfcw_properties(nzg,nzs,integration_buff(ibuff)%y,dtrk4,csite   &
                                          ,ipa)
               call update_diagnostic_vars(integration_buff(ibuff)%y,csite,ipa)

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
                  integration_buff(ibuff)%y%wbudget_loss2runoff = wfreeb                   &
                        + integration_buff(ibuff)%y%wbudget_loss2runoff
                  integration_buff(ibuff)%y%ebudget_loss2runoff = qwfree                   &
                        + integration_buff(ibuff)%y%ebudget_loss2runoff
                  integration_buff(ibuff)%y%wbudget_storage     =                          &
                                     integration_buff(ibuff)%y%wbudget_storage - wfreeb
                  integration_buff(ibuff)%y%ebudget_storage    =                           &
                                     integration_buff(ibuff)%y%ebudget_storage - qwfree
               end if
            end if
         end if

         !------ Copy the temporary patch to the next intermediate step -------------------!
         call copy_rk4_patch(integration_buff(ibuff)%y,integration_buff(ibuff)%initp,cpatch)

         !------ Update the substep for next time and leave -------------------------------!
         csite%htry(ipa) = sngl(hnext)

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
   write (unit=*,fmt='(a)') ' ==> Too many steps in routine heun_integ'
   call print_rk4patch(integration_buff(ibuff)%ytemp, csite,ipa)

   return
end subroutine heun_integ
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the integration step of the Heun's method.                   !
! Here we use several buffers, so a quick guide of what each one means...                  !
!                                                                                          !
! 1. integration_buff(ibuff)%y      => This is the state y(t)                              !
! 2. integration_buff(ibuff)%ytemp  => This is the state y(t+h)                            !
! 3. integration_buff(ibuff)%yerr   => This is the error e(t+h)                            !
! 4. integration_buff(ibuff)%dydx   => This is the Runge-Kutta term K1                     !
! 5. integration_buff(ibuff)%ak2    => This is the Runge-Kutta term K2                     !
! 6. integration_buff(ibuff)%ak3    => This is the next step using Euler:                  !
!                                      ye(t+h) = y(t) + K1 * h                             !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine heun_stepper(x,h,csite,ipa,reject_step,reject_result)

   use rk4_coms      , only : integration_buff    & ! structure
                            , integration_vars    & ! structure
                            , zero_rk4_patch      & ! subroutine
                            , zero_rk4_cohort     & ! subroutine
                            , rk4site             & ! intent(in)
                            , print_diags         & ! intent(in)
                            , heun_a2             & ! intent(in)
                            , heun_b21            & ! intent(in)
                            , heun_c1             & ! intent(in)
                            , heun_c2             & ! intent(in)
                            , heun_dc1            & ! intent(in)
                            , heun_dc2            ! ! intent(in)
   use rk4_stepper   , only : rk4_sanity_check    ! ! subroutine
   use ed_state_vars , only : sitetype            & ! structure
                            , patchtype           ! ! structure
   use grid_coms     , only : nzg                 & ! intent(in)
                            , nzs                 ! ! structure
   ! OMP use omp_lib

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)    , target      :: csite
   integer           , intent(in)  :: ipa
   logical           , intent(out) :: reject_step
   logical           , intent(out) :: reject_result
   real(kind=8)      , intent(in)  :: x
   real(kind=8)      , intent(in)  :: h
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)   , pointer     :: cpatch
   integer                         :: ibuff
   !---------------------------------------------------------------------------------------!

   ibuff = 1

   !---------------------------------------------------------------------------------------!
   !     Start and assume that nothing went wrong up to this point... If we find any       !
   ! seriously bad step, quit and reduce the time step without even bothering to try       !
   ! further.                                                                              !
   !---------------------------------------------------------------------------------------!
   reject_step   = .false.
   reject_result = .false.
   cpatch => csite%patch(ipa)
   !---------------------------------------------------------------------------------------!



   !----- yeuler is the temporary array with the Euler step with no correction. -----------!
   call copy_rk4_patch(integration_buff(ibuff)%y,integration_buff(ibuff)%ak3,cpatch)
   call inc_rk4_patch (integration_buff(ibuff)%ak3,integration_buff(ibuff)%dydx            &
                      ,heun_b21*h, cpatch)
   call update_diagnostic_vars(integration_buff(ibuff)%ak3      ,csite,ipa)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Check to see if the Euler result makes sense.  Since we will compute the           !
   ! derivative correction using it, the Euler step must be bounded.  If not, reject the   !
   ! step and try a smaller step size.                                                     !
   !---------------------------------------------------------------------------------------!
   call rk4_sanity_check(integration_buff(ibuff)%ak3,reject_step,csite,ipa                 &
                        ,integration_buff(ibuff)%dydx,h,print_diags)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute the second term (correction) of the derivative, using the Euler's         !
   ! predicted state.                                                                      !
   !---------------------------------------------------------------------------------------!
   call leaf_derivs(integration_buff(ibuff)%ak3,integration_buff(ibuff)%ak2, csite,ipa,h   &
                   ,.false.)
   !---------------------------------------------------------------------------------------!



   !----- We now combine both derivatives and find the potential next step. ---------------!
   call copy_rk4_patch(integration_buff(ibuff)%y,integration_buff(ibuff)%ytemp, cpatch)
   call inc_rk4_patch(integration_buff(ibuff)%ytemp,integration_buff(ibuff)%dydx           &
                     ,heun_c1*h, cpatch)
   call inc_rk4_patch(integration_buff(ibuff)%ytemp,integration_buff(ibuff)%ak2            &
                     ,heun_c2*h, cpatch)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Update the diagnostic properties and make final adjustments.  This time we will  !
   ! run the full adjustment, to make sure that the step will be rejected especially if    !
   ! there are issues with the top soil properties.                                        !
   !---------------------------------------------------------------------------------------!
   call update_diagnostic_vars   (integration_buff(ibuff)%ytemp        ,csite,ipa)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check to see if this attempt of advancing one time step makes sense.  If not,      !
   ! reject the result and try a smaller step size.                                        !
   !---------------------------------------------------------------------------------------!
   call rk4_sanity_check(integration_buff(ibuff)%ytemp,reject_result,csite,ipa             &
                        ,integration_buff(ibuff)%ak2,h,print_diags)
   if(reject_result)return
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute the estimate of the error associated with the step.                       !
   !---------------------------------------------------------------------------------------!
   call zero_rk4_patch (integration_buff(ibuff)%yerr)
   call zero_rk4_cohort(integration_buff(ibuff)%yerr)
   call inc_rk4_patch(integration_buff(ibuff)%yerr,integration_buff(ibuff)%dydx            &
                     ,heun_dc1*h,cpatch)
   call inc_rk4_patch(integration_buff(ibuff)%yerr,integration_buff(ibuff)%ak2             &
                     ,heun_dc2*h,cpatch)
   !---------------------------------------------------------------------------------------!

   return
end subroutine heun_stepper
!==========================================================================================!
!==========================================================================================!
