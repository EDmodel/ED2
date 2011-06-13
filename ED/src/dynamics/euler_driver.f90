!==========================================================================================!
!==========================================================================================!
!     This subroutine is the main driver for the Euler integration scheme.                 !
!------------------------------------------------------------------------------------------!
subroutine euler_timestep(cgrid)
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
   use ed_misc_coms          , only : dtlsm              ! ! intent(in)
   use ed_max_dims           , only : n_dbh              ! ! intent(in)
   use soil_coms             , only : soil_rough         & ! intent(in)
                                    , snow_rough         ! ! intent(in)
   use consts_coms           , only : cp                 & ! intent(in)
                                    , mmdryi             & ! intent(in)
                                    , day_sec            & ! intent(in)
                                    , umol_2_kgC         ! ! intent(in)
   use canopy_struct_dynamics, only : canopy_turbulence8 ! ! subroutine

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
   integer                                :: ico
   integer                                :: nsteps
   real                                   :: thetaatm
   real                                   :: thetacan
   real                                   :: rasveg
   real                                   :: storage_decay
   real                                   :: leaf_flux
   real                                   :: veg_tai
   real                                   :: wcurr_loss2atm
   real                                   :: ecurr_loss2atm
   real                                   :: co2curr_loss2atm
   real                                   :: wcurr_loss2drainage
   real                                   :: ecurr_loss2drainage
   real                                   :: wcurr_loss2runoff
   real                                   :: ecurr_loss2runoff
   real                                   :: old_can_theiv
   real                                   :: old_can_shv
   real                                   :: old_can_co2
   real                                   :: old_can_rhos
   real                                   :: old_can_temp
   real                                   :: fm
   !----- External functions. -------------------------------------------------------------!
   real, external                         :: compute_netrad
   !---------------------------------------------------------------------------------------!

   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)
         cmet  => cpoly%met(isi)

         patchloop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)


            !----- Reset all buffers to zero, as a safety measure. ------------------------!
            call zero_rk4_patch(integration_buff%initp)
            call zero_rk4_patch(integration_buff%dinitp)
            call zero_rk4_patch(integration_buff%ytemp)
            call zero_rk4_patch(integration_buff%yerr)
            call zero_rk4_patch(integration_buff%yscal)
            call zero_rk4_patch(integration_buff%dydx)
            call zero_rk4_cohort(integration_buff%initp)
            call zero_rk4_cohort(integration_buff%dinitp)
            call zero_rk4_cohort(integration_buff%ytemp)
            call zero_rk4_cohort(integration_buff%yerr)
            call zero_rk4_cohort(integration_buff%yscal)
            call zero_rk4_cohort(integration_buff%dydx)

            !----- Get velocity for aerodynamic resistance. -------------------------------!
            if (csite%can_theta(ipa) < cmet%atm_theta) then
               cmet%vels = cmet%vels_stab
            else
               cmet%vels = cmet%vels_unstab
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Update roughness and canopy depth.                                        !
            !------------------------------------------------------------------------------!
            call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs,cpoly%ntext_soil(:,isi))
            call update_patch_derived_props(csite,cpoly%lsl(isi),cmet%prss,ipa)
            !------------------------------------------------------------------------------!



            !----- Save the previous thermodynamic state. ---------------------------------!
            old_can_theiv    = csite%can_theiv(ipa)
            old_can_shv      = csite%can_shv(ipa)
            old_can_co2      = csite%can_co2(ipa)
            old_can_rhos     = csite%can_rhos(ipa)
            old_can_temp     = csite%can_temp(ipa)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Copy the meteorological variables to the rk4site structure.               !
            !------------------------------------------------------------------------------!
            call copy_met_2_rk4site(nzg,cmet%vels,cmet%atm_theiv,cmet%atm_theta            &
                                   ,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2,cmet%geoht      &
                                   ,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg,cmet%prss   &
                                   ,cmet%rshort,cmet%rlong,cmet%geoht,cpoly%lsl(isi)       &
                                   ,cpoly%ntext_soil(:,isi)                                &
                                   ,cpoly%green_leaf_factor(:,isi)                         &
                                   ,cgrid%lon(ipy),cgrid%lat(ipy),cgrid%cosz(ipy))

            !----- Compute current storage terms. -----------------------------------------!
            call update_budget(csite,cpoly%lsl(isi),ipa,ipa)


            !------------------------------------------------------------------------------!
            !     Set up the integration patch.                                            !
            !------------------------------------------------------------------------------!
            call copy_patch_init(csite,ipa,integration_buff%initp)

            !----- Get photosynthesis, stomatal conductance, and transpiration. -----------!
            call canopy_photosynthesis(csite,cmet,nzg,ipa,cpoly%lsl(isi)                   &
                                      ,cpoly%ntext_soil(:,isi)                             &
                                      ,cpoly%leaf_aging_factor(:,isi)                      &
                                      ,cpoly%green_leaf_factor(:,isi))

            !----- Compute root and heterotrophic respiration. ----------------------------!
            call soil_respiration(csite,ipa,nzg,cpoly%ntext_soil(:,isi))

            !------------------------------------------------------------------------------!
            !     Set up the remaining, carbon-dependent variables to the buffer.          !
            !------------------------------------------------------------------------------!
            call copy_patch_init_carbon(csite,ipa,integration_buff%initp)


            !------------------------------------------------------------------------------!
            !     This is the step in which the derivatives are computed, we a structure   !
            ! that is very similar to the Runge-Kutta, though a simpler one.               !
            !------------------------------------------------------------------------------!
            call integrate_patch_euler(csite,integration_buff%initp                        &
                                      ,integration_buff%dinitp,integration_buff%ytemp      &
                                      ,integration_buff%yscal,integration_buff%yerr        &
                                      ,integration_buff%dydx,ipa,wcurr_loss2atm            &
                                      ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage &
                                      ,ecurr_loss2drainage,wcurr_loss2runoff               &
                                      ,ecurr_loss2runoff,nsteps)

            !----- Add the number of steps into the step counter. -------------------------!
            cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)

            !------------------------------------------------------------------------------!
            !    Update the minimum monthly temperature, based on canopy temperature.      !
            !------------------------------------------------------------------------------!
            if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
               cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
            end if
               
            !------------------------------------------------------------------------------!
            !     Compute the residuals.                                                   !
            !------------------------------------------------------------------------------!
            call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa              &
                               ,wcurr_loss2atm,ecurr_loss2atm,co2curr_loss2atm             &
                               ,wcurr_loss2drainage,ecurr_loss2drainage,wcurr_loss2runoff  &
                               ,ecurr_loss2runoff,cpoly%area(isi),cgrid%cbudget_nep(ipy)   &
                               ,old_can_theiv,old_can_shv,old_can_co2,old_can_rhos         &
                               ,old_can_temp)
         end do patchloop
      end do siteloop
   end do polyloop

   return
end subroutine euler_timestep
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the integration process using the Euler method.  Notice   !
! that most of the Euler method utilises the subroutines from Runge-Kutta.                 !
!------------------------------------------------------------------------------------------!
subroutine integrate_patch_euler(csite,initp,dinitp,ytemp,yscal,yerr,dydx,ipa              &
                                ,wcurr_loss2atm,ecurr_loss2atm,co2curr_loss2atm            &
                                ,wcurr_loss2drainage,ecurr_loss2drainage,wcurr_loss2runoff &
                                ,ecurr_loss2runoff,nsteps)
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
   use rk4_driver      , only : initp2modelp         ! ! subroutine

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target      :: csite
   type(rk4patchtype)    , target      :: initp
   type(rk4patchtype)    , target      :: dinitp
   type(rk4patchtype)    , target      :: ytemp
   type(rk4patchtype)    , target      :: yscal
   type(rk4patchtype)    , target      :: yerr
   type(rk4patchtype)    , target      :: dydx
   integer               , intent(in)  :: ipa
   real                  , intent(out) :: wcurr_loss2atm
   real                  , intent(out) :: ecurr_loss2atm
   real                  , intent(out) :: co2curr_loss2atm
   real                  , intent(out) :: wcurr_loss2drainage
   real                  , intent(out) :: ecurr_loss2drainage
   real                  , intent(out) :: wcurr_loss2runoff
   real                  , intent(out) :: ecurr_loss2runoff
   integer               , intent(out) :: nsteps
   !----- Local variables -----------------------------------------------------------------!
   real(kind=8)                        :: hbeg
   !----- Locally saved variable ----------------------------------------------------------!
   logical                  , save     :: first_time=.true.
   !---------------------------------------------------------------------------------------!

   !----- Assigning some constants which will remain the same throughout the run. ---------!
   if (first_time) then
      first_time = .false.
      tbeg   = 0.d0
      tend   = dble(dtlsm)
      dtrk4  = tend - tbeg
      dtrk4i = 1.d0/dtrk4
   end if

   !---------------------------------------------------------------------------------------!
   !      Initial step size.  Experience has shown that giving this too large a value      !
   ! causes the integrator to fail (e.g., soil layers become supersaturated).              !
   !---------------------------------------------------------------------------------------!
   hbeg = dble(csite%htry(ipa))

   !---------------------------------------------------------------------------------------!
   !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm,    !
   ! so they must be zeroed at each call.                                                  !
   !---------------------------------------------------------------------------------------!
   initp%upwp = 0.d0
   initp%tpwp = 0.d0
   initp%qpwp = 0.d0
   initp%cpwp = 0.d0
   initp%wpwp = 0.d0

   !----- Go into the ODE integrator using Euler. -----------------------------------------!
   call euler_integ(hbeg,csite,initp,dinitp,ytemp,yscal,yerr,dydx,ipa,nsteps)

   !---------------------------------------------------------------------------------------!
   !      Normalize canopy-atmosphere flux values.  These values are updated every         !
   ! dtlsm, so they must be normalized every time.                                         !
   !---------------------------------------------------------------------------------------!
   initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
   initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
   initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
   initp%cpwp = initp%can_rhos * initp%cpwp * dtrk4i
   initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i
      
   !---------------------------------------------------------------------------------------!
   ! Move the state variables from the integrated patch to the model patch.                !
   !---------------------------------------------------------------------------------------!
   call initp2modelp(tend-tbeg,initp,csite,ipa,wcurr_loss2atm,ecurr_loss2atm               &
                    ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage              &
                    ,wcurr_loss2runoff,ecurr_loss2runoff)

   return
end subroutine integrate_patch_euler
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Subroutine odeint                                                                        !
!                                                                                          !
!     This subroutine will drive the integration of several ODEs that drive the fast-scale !
! state variables.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine euler_integ(h1,csite,initp,dinitp,ytemp,yscal,yerr,dydx,ipa,nsteps)
   use ed_state_vars  , only : sitetype               & ! structure
                             , patchtype              & ! structure
                             , polygontype            ! ! structure
   use rk4_coms       , only : integration_vars       & ! structure
                             , rk4patchtype           & ! structure
                             , rk4site                & ! intent(in)
                             , print_diags            & ! intent(in)
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
                             , print_detailed         & ! intent(in)
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
                             , runoff_time            ! ! intent(in)
   use consts_coms    , only : cliq8                  & ! intent(in)
                             , t3ple8                 & ! intent(in)
                             , tsupercool8            & ! intent(in)
                             , wdnsi8                 ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite            ! Current site
   type(rk4patchtype)        , target      :: initp            ! Current integ. patch
   type(rk4patchtype)        , target      :: dinitp           ! Integration derivative
   type(rk4patchtype)        , target      :: ytemp            ! Temporary integ. patch
   type(rk4patchtype)        , target      :: yscal            ! Scale for error analysis
   type(rk4patchtype)        , target      :: yerr             ! Patch integration error
   type(rk4patchtype)        , target      :: dydx             ! Patch integration error
   integer                   , intent(in)  :: ipa              ! Current patch ID
   real(kind=8)              , intent(in)  :: h1               ! First guess of delta-t
   integer                   , intent(out) :: nsteps           ! Number of steps taken.
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch           ! Current patch
   logical                                 :: reject_step      ! Should I reject the step?
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
   !----- Saved variables -----------------------------------------------------------------!
   logical                   , save        :: first_time=.true.
   logical                   , save        :: simplerunoff
   real(kind=8)              , save        :: runoff_time_i
   !----- External function. --------------------------------------------------------------!
   real                      , external    :: sngloff
   !---------------------------------------------------------------------------------------!
   
   !----- Checking whether we will use runoff or not, and saving this check to save time. -!
   if (first_time) then
      simplerunoff = useRUNOFF == 0 .and. runoff_time /= 0.
      if (runoff_time /= 0.) then
         runoff_time_i = 1.d0/dble(runoff_time)
      else 
         runoff_time_i = 0.d0
      end if
      first_time   = .false.
   end if

   !----- Use some aliases for simplicity. ------------------------------------------------!
   cpatch => csite%patch(ipa)

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
      call leaf_derivs(initp,dinitp,csite,ipa)

      !----- Get scalings used to determine stability -------------------------------------!
      call get_yscal(initp,dinitp,h,yscal,cpatch)

      !----- Be sure not to overstep ------------------------------------------------------!
      if((x+h-tend)*(x+h-tbeg) > 0.d0) h=tend-x

      !------------------------------------------------------------------------------------!
      !     Here we will perform the Euler integration using the time step.  As in Runge-  !
      ! Kutta, we also check whether the integration is going well and if needed we shrink !
      ! the intermediate time steps.  However, we cannot estimate the errors as in Runge-  !
      ! Kutta or Euler, so steps are shrunk only when it fails the sanity check.           !
      !------------------------------------------------------------------------------------!
      reject_step =  .false.
      hstep:   do

         !----- Copy patch to the temporary structure. ------------------------------------!
         call copy_rk4_patch(initp, ytemp, cpatch)

         !---------------------------------------------------------------------------------!
         !    Integrate, then update and correct diagnostic variables to avoid overshoot-  !
         ! ing, provided that the overshooting is small.                                   !
         !---------------------------------------------------------------------------------!
         call inc_rk4_patch(ytemp,dinitp,h,cpatch)
         call update_diagnostic_vars(ytemp,csite,ipa)

         !----- Perform a sanity check. ---------------------------------------------------!
         call rk4_sanity_check(ytemp,reject_step,csite,ipa,dinitp,h,print_diags)

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

            !----- Take the derivative of the upcoming step. ------------------------------!
            call leaf_derivs(ytemp,dydx,csite,ipa)
         end if


         !---------------------------------------------------------------------------------!
         ! 3. If the step failed, then calculate a new shorter step size to try.           !
         !---------------------------------------------------------------------------------!
         if (reject_step) then
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
               write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN EULER_INT'
               write (unit=*,fmt='(80a)')         ('-',k=1,80)
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LONGITUDE:     ',rk4site%lon
               write (unit=*,fmt='(a,1x,f9.4)')   ' + LATITUDE:      ',rk4site%lat
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

               call rk4_sanity_check(ytemp,test_reject,csite,ipa,dinitp,h,.true.)
               call print_sanity_check(ytemp,csite,ipa)
               call print_rk4patch(ytemp, csite,ipa)
            end if

         else
            !------------------------------------------------------------------------------!
            ! 3b.  Great, it worked, so now we can advance to the next step.  We just need !
            !      to do some minor adjustments before...                                  !
            !------------------------------------------------------------------------------!
            !----- i.   Final update of leaf properties to avoid negative water. ----------!
            call adjust_veg_properties(ytemp,h,csite,ipa)
            !----- ii.  Final update of top soil properties to avoid off-bounds moisture. -!
            call adjust_topsoil_properties(ytemp,h,csite,ipa)
            !----- ii. Make temporary surface water stable and positively defined. --------!
            call adjust_sfcw_properties(nzg,nzs,ytemp,h,csite,ipa)
            !----- iii.  Update the diagnostic variables. ---------------------------------!
            call update_diagnostic_vars(ytemp, csite,ipa)
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
               call norm_rk4_fluxes(ytemp,h)
               call print_rk4_state(initp,ytemp,csite,ipa,x,h)
            end if

            !----- 3e. Copy the temporary structure to the intermediate state. ------------!
            call copy_rk4_patch(ytemp, initp,csite%patch(ipa))
            call copy_rk4_patch(dydx ,dinitp,csite%patch(ipa))

            !------------------------------------------------------------------------------!
            !    3f. Flush step-by-step fluxes to zero if the user wants detailed          !
            !        debugging.                                                            !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               call reset_rk4_fluxes(initp)
            end if

            !----- 3g. Update time. -------------------------------------------------------!
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
         
            if (initp%sfcwater_mass(ksn)    > 0.d0   .and.                                 &
                initp%sfcwater_fracliq(ksn) > 1.d-1        ) then

               wfreeb = min(1.d0,dtrk4*runoff_time_i) * initp%sfcwater_mass(ksn)           &
                      * (initp%sfcwater_fracliq(ksn) - 1.d-1) / 9.d-1
               qwfree = wfreeb * cliq8 * (initp%sfcwater_tempk(ksn) - tsupercool8 )

               initp%sfcwater_mass(ksn) = initp%sfcwater_mass(ksn)   - wfreeb
               initp%sfcwater_depth(ksn) = initp%sfcwater_depth(ksn) - wfreeb * wdnsi8
               !----- Recompute the energy removing runoff --------------------------------!
               initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn) - qwfree

               call adjust_sfcw_properties(nzg,nzs,initp,dtrk4,csite,ipa)
               call update_diagnostic_vars(initp,csite,ipa)

               !----- Compute runoff for output -------------------------------------------!
               if (fast_diagnostics) then
                  csite%runoff(ipa) = csite%runoff(ipa)                                    &
                                    + sngloff(wfreeb * dtrk4i,tiny_offset)
                  csite%avg_runoff(ipa) = csite%avg_runoff(ipa)                            &
                                        + sngloff(wfreeb * dtrk4i,tiny_offset)
                  csite%avg_runoff_heat(ipa) = csite%avg_runoff_heat(ipa)                  &
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
         csite%htry(ipa) = sngl(hnext)

         !---------------------------------------------------------------------------------!
         !     Update the average time step.  The square of DTLSM (tend-tbeg) is needed    !
         ! because we will divide this by the time between t0 and t0+frqsum.               !
         !---------------------------------------------------------------------------------!
         csite%avg_rk4step(ipa) = csite%avg_rk4step(ipa)                                   &
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
end subroutine euler_integ
!==========================================================================================!
!==========================================================================================!



