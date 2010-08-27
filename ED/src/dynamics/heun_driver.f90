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
   use grid_coms             , only : nzg                ! ! intent(in)
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
   type(edtype)             , target    :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype)        , pointer   :: cpoly
   type(sitetype)           , pointer   :: csite
   type(patchtype)          , pointer   :: cpatch
   type(met_driv_state)     , pointer   :: cmet
   integer                              :: ipy
   integer                              :: isi
   integer                              :: ipa
   integer                              :: ico
   integer                              :: nsteps
   integer, dimension(nzg)              :: ed_ktrans
   real                                 :: thetaatm
   real                                 :: thetacan
   real                                 :: rasveg
   real                                 :: storage_decay
   real                                 :: leaf_flux
   real                                 :: veg_tai
   real                                 :: sum_lai_rbi
   real                                 :: wcurr_loss2atm
   real                                 :: ecurr_loss2atm
   real                                 :: co2curr_loss2atm
   real                                 :: wcurr_loss2drainage
   real                                 :: ecurr_loss2drainage
   real                                 :: wcurr_loss2runoff
   real                                 :: ecurr_loss2runoff
   real                                 :: old_can_enthalpy
   real                                 :: old_can_shv
   real                                 :: old_can_co2
   real                                 :: old_can_rhos
   real                                 :: fm
   !----- External functions. -------------------------------------------------------------!
   real, external                       :: compute_netrad
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
            call zero_rk4_patch(integration_buff%ytemp)
            call zero_rk4_patch(integration_buff%yerr)
            call zero_rk4_patch(integration_buff%yscal)
            call zero_rk4_patch(integration_buff%dydx)
            call zero_rk4_cohort(integration_buff%initp)
            call zero_rk4_cohort(integration_buff%ytemp)
            call zero_rk4_cohort(integration_buff%yerr)
            call zero_rk4_cohort(integration_buff%yscal)
            call zero_rk4_cohort(integration_buff%dydx)

            !----- Save the previous thermodynamic state. ---------------------------------!
            old_can_enthalpy = csite%can_enthalpy(ipa)
            old_can_shv      = csite%can_shv(ipa)
            old_can_co2      = csite%can_co2(ipa)
            old_can_rhos     = csite%can_rhos(ipa)

            !----- Get velocity for aerodynamic resistance. -------------------------------!
            if (csite%can_theta(ipa) < cmet%atm_theta) then
               cmet%vels = cmet%vels_stab
            else
               cmet%vels = cmet%vels_unstab
            end if

            !------------------------------------------------------------------------------!
            !    Copy the meteorological variables to the rk4site structure.               !
            !------------------------------------------------------------------------------!
            call copy_met_2_rk4site(cmet%vels,cmet%atm_enthalpy,cmet%atm_theta             &
                                   ,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2,cmet%geoht      &
                                   ,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg,cmet%prss   &
                                   ,cmet%geoht,cpoly%lsl(isi)                              &
                                   ,cpoly%green_leaf_factor(:,isi)                         &
                                   ,cgrid%lon(ipy),cgrid%lat(ipy))

            !----- Compute current storage terms. -----------------------------------------!
            call update_budget(csite,cpoly%lsl(isi),ipa,ipa)


            !------------------------------------------------------------------------------!
            !     Set up the integration patch.                                            !
            !------------------------------------------------------------------------------!
            call copy_patch_init(csite,ipa,integration_buff%initp)

            !------------------------------------------------------------------------------!
            !     Here we compute canopy turbulence-related variables, such as the rough-  !
            ! ness scale, the characteristic scales (stars) and canopy resistance and      !
            ! capacities.                                                                  !
            !------------------------------------------------------------------------------!
            call canopy_turbulence8(csite,integration_buff%initp,ipa,.true.)

            !----- Get photosynthesis, stomatal conductance, and transpiration. -----------!
            call canopy_photosynthesis(csite,ipa,cmet%vels,cmet%atm_tmp,cmet%prss          &
                                      ,ed_ktrans,csite%ntext_soil(:,ipa)                   &
                                      ,csite%soil_water(:,ipa),csite%soil_fracliq(:,ipa)   &
                                      ,cpoly%lsl(isi),sum_lai_rbi                          &
                                      ,cpoly%leaf_aging_factor(:,isi)                      &
                                      ,cpoly%green_leaf_factor(:,isi))

            !----- Compute root and heterotrophic respiration. ----------------------------!
            call soil_respiration(csite,ipa)

            !------------------------------------------------------------------------------!
            !     Set up the remaining, carbon-dependent variables to the buffer.          !
            !------------------------------------------------------------------------------!
            call copy_patch_init_carbon(csite,ipa,integration_buff%initp)


            !------------------------------------------------------------------------------!
            !     This is the step in which the derivatives are computed, we a structure   !
            ! that is very similar to the Runge-Kutta, though a simpler one.               !
            !------------------------------------------------------------------------------!
            call integrate_patch_heun(csite,integration_buff%initp,integration_buff%y      &
                                     ,integration_buff%ytemp,integration_buff%yscal        &
                                     ,integration_buff%yerr,integration_buff%dydx          &
                                     ,integration_buff%ak2,integration_buff%ak3,ipa        &
                                     ,wcurr_loss2atm,ecurr_loss2atm,co2curr_loss2atm       &
                                     ,wcurr_loss2drainage,ecurr_loss2drainage              &
                                     ,wcurr_loss2runoff,ecurr_loss2runoff,nsteps)

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
                               ,old_can_enthalpy,old_can_shv,old_can_co2,old_can_rhos)
         end do patchloop
      end do siteloop
   end do polyloop

   return
end subroutine heun_timestep
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the integration process using the Heun method.  Notice    !
! that most of the Heun method utilises the subroutines from Runge-Kutta.                  !
!------------------------------------------------------------------------------------------!
subroutine integrate_patch_heun(csite,initp,y0,ytry,yscal,yerr,ak1,ak2,yeuler,ipa          &
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
                              , dtrk4i               & ! intent(inout)
                              , ibranch_thermo       & ! intent(in)
                              , effarea_water        & ! intent(out)
                              , effarea_heat         ! ! intent(out)
   use rk4_driver      , only : initp2modelp         ! ! subroutine

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target      :: csite
   type(rk4patchtype)    , target      :: initp
   type(rk4patchtype)    , target      :: y0
   type(rk4patchtype)    , target      :: ytry
   type(rk4patchtype)    , target      :: yscal
   type(rk4patchtype)    , target      :: yerr
   type(rk4patchtype)    , target      :: ak1
   type(rk4patchtype)    , target      :: ak2
   type(rk4patchtype)    , target      :: yeuler
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
      
      !------------------------------------------------------------------------------------!
      !    The area factor for heat and water exchange between canopy and vegetation is    !
      ! applied only on LAI, and it depends on how we are considering the branches and     !
      ! twigs.  If their area isn't explicitly defined, we add a 0.2 factor to the area    !
      ! because WPA will be 0.  Otherwise, we don't add anything to the LAI, and let WPA   !
      ! to do the job.                                                                     !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (0)
         effarea_water = 1.2d0
         effarea_heat  = 2.2d0
      case default
         effarea_water = 1.0d0
         effarea_heat  = 2.0d0
      end select
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
   call heun_integ(hbeg,csite,initp,y0,ytry,yscal,yerr,ak1,ak2,yeuler,ipa,nsteps)

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
subroutine heun_integ(h1,csite,initp,y0,ytry,yscal,yerr,ak1,ak2,yeuler,ipa,nsteps)
   use ed_state_vars  , only : sitetype               & ! structure
                             , patchtype              & ! structure
                             , polygontype            ! ! structure
   use rk4_coms       , only : integration_vars       & ! structure
                             , rk4patchtype           & ! structure
                             , rk4site                & ! intent(in)
                             , rk4min_sfcwater_mass   & ! intent(in)
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
                             , errcon                 ! ! intent(in)
   use rk4_stepper    , only : rk4_sanity_check       & ! subroutine
                             , print_sanity_check     ! ! subroutine
   use ed_misc_coms   , only : fast_diagnostics       ! ! intent(in)
   use hydrology_coms , only : useRUNOFF              ! ! intent(in)
   use grid_coms      , only : nzg                    & ! intent(in)
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
   type(rk4patchtype)        , target      :: y0               ! Temporary integ. patch
   type(rk4patchtype)        , target      :: ytry             ! Temporary integ. patch
   type(rk4patchtype)        , target      :: yscal            ! Scale for error analysis
   type(rk4patchtype)        , target      :: yerr             ! Patch integration error
   type(rk4patchtype)        , target      :: ak1              ! Euler derivative
   type(rk4patchtype)        , target      :: ak2              ! Correction derivative
   type(rk4patchtype)        , target      :: yeuler           ! Euler integration
   integer                   , intent(in)  :: ipa              ! Current patch ID
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
   !----- Saved variables -----------------------------------------------------------------!
   logical                   , save        :: first_time=.true.
   logical                   , save        :: simplerunoff
   real(kind=8)              , save        :: runoff_time_i
   !----- External function. --------------------------------------------------------------!
   real                      , external    :: sngloff
   !----- Local constant. -----------------------------------------------------------------!
   logical                   , parameter   :: print_intstep = .false.
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

      if (print_intstep) then
         write (unit=66,fmt='(27(a,1x))')                                                  &
                       '  ELAPSED_TIME','     TIME_STEP','       ATM_SHV','       CAN_SHV' &
                      ,'    GROUND_SHV','     ATM_THETA','     CAN_THETA','      ATM_TEMP' &
                      ,'      CAN_TEMP','  ATM_ENTHALPY','  CAN_ENTHALPY','      ATM_PRSS' &
                      ,'      CAN_PRSS','      ATM_RHOS','      CAN_RHOS','      ATM_VELS' &
                      ,'     CAN_DEPTH','          RAIN','         USTAR','         TSTAR' &
                      ,'         ESTAR','         QSTAR','          ZETA','        RIBULK' &
                      ,'     SFCW_MASS','    SFCW_DEPTH','    SOIL_WATER'
      end if
   end if

   !----- Use some aliases for simplicity. ------------------------------------------------!
   cpatch => csite%patch(ipa)

   !---------------------------------------------------------------------------------------!
   !    If top snow layer is too thin for computational stability, have it evolve in       !
   ! thermal equilibrium with top soil layer.                                              !
   !---------------------------------------------------------------------------------------!
   call redistribute_snow(initp, csite,ipa)
   call update_diagnostic_vars(initp,csite,ipa)

   !---------------------------------------------------------------------------------------!
   !     Create temporary patches.                                                         !
   !---------------------------------------------------------------------------------------!
   call copy_rk4_patch(initp,y0,cpatch)

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
      call leaf_derivs(y0,ak1,csite,ipa)

      !----- Get scalings used to determine stability -------------------------------------!
      call get_yscal(y0,ak1,h,yscal,cpatch)

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
         call heun_stepper(y0,ytry,yerr,ak1,ak2,yeuler,x,h,csite,ipa,reject_step           &
                          ,reject_result)

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
            call get_errmax(errmax,yerr,yscal,csite%patch(ipa),initp,ytry)
            !----- Scale the error based on the prescribed tolerance. ---------------------!
            errmax = errmax * rk4epsi
         end if


         !---------------------------------------------------------------------------------!
         ! 3. If that error was large, then calculate a new step size to try.  There are   !
         !    two types of new tries.  If step failed to be minimally reasonable (reject-  !
         !    ed) we have assigned a standard large error (10.0).  Otherwise a new step is !
         !    calculated based on the size of that error.  Hopefully, those new steps      !
         !    should be less than the previous h.  If the error was small, i.e. less then  !
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
                  call rk4_sanity_check(ytry,test_reject,csite,ipa,ak1,h,.true.)
                  call print_sanity_check(y0,csite,ipa)
               elseif (reject_step) then
                  call rk4_sanity_check(yeuler,test_reject,csite,ipa,ak1,h,.true.)
                  call print_sanity_check(y0,csite,ipa)
               else
                  call print_errmax(errmax,yerr,yscal,csite%patch(ipa),y0,ytry)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax*rk4eps
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_rk4patch(y0, csite,ipa)
            end if

         else
            !------------------------------------------------------------------------------!
            ! 3b.  Great, it worked, so now we can advance to the next step.  We just need !
            !      to do some minor adjustments before...                                  !
            !------------------------------------------------------------------------------!
            !----- i.   Final update of leaf properties to avoid negative water. ----------!
            call adjust_veg_properties(ytry,h,csite,ipa)
            !----- ii.  Final update of top soil properties to avoid off-bounds moisture. -!
            call adjust_topsoil_properties(ytry,h,csite,ipa)
            !----- iii. Make snow layers stable and positively defined. -------------------!
            call redistribute_snow(ytry, csite,ipa)
            !----- iv.  Update the diagnostic variables. ----------------------------------!
            call update_diagnostic_vars(ytry, csite,ipa)
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

            !----- 3d. Copying the temporary structure to the intermediate state. ---------!
            call copy_rk4_patch(ytry, y0,csite%patch(ipa))

            !----- 3e. Print the output if the user wants it. -----------------------------!
            if (print_intstep) then
               write (unit=66,fmt='(27(es14.7,1x))')                                       &
                     elaptime,h,rk4site%atm_shv,y0%can_shv,y0%ground_shv,rk4site%atm_theta &
                    ,y0%can_theta,rk4site%atm_tmp,y0%can_temp,rk4site%atm_enthalpy         &
                    ,y0%can_enthalpy,rk4site%atm_prss,y0%can_prss,rk4site%rhos,y0%can_rhos &
                    ,rk4site%vels,y0%can_depth,rk4site%pcpg,y0%ustar,y0%tstar,y0%estar     &
                    ,y0%qstar,y0%zeta,y0%ribulk,y0%sfcwater_mass(1),y0%sfcwater_depth(1)   &
                    ,y0%soil_water(nzg)
            end if

            !----- 3f. Updating time. -----------------------------------------------------!
            x        = x + h
            h        = hnext
            elaptime = elaptime + h

            exit hstep
         end if
      end do hstep

      !----- If the integration reached the next step, make some final adjustments --------!
      if((x-tend)*dtrk4 >= 0.d0)then

         ksn = y0%nlev_sfcwater

         !---------------------------------------------------------------------------------!
         !   Make temporary surface liquid water disappear.  This will not happen          !
         ! immediately, but liquid water will decay with the time scale defined by         !
         ! runoff_time scale. If the time scale is too tiny, then it will be forced to be  !
         ! hdid (no reason to be faster than that).                                        !
         !---------------------------------------------------------------------------------!
         if (simplerunoff .and. ksn >= 1) then
         
            if (y0%sfcwater_mass(ksn) > 0.d0 .and. y0%sfcwater_fracliq(ksn) > 1.d-1 ) then

               wfreeb = min(1.d0,dtrk4*runoff_time_i) * y0%sfcwater_mass(ksn)              &
                      * (y0%sfcwater_fracliq(ksn) - 1.d-1) / 9.d-1
               qwfree = wfreeb * cliq8 * (y0%sfcwater_tempk(ksn) - tsupercool8 )

               y0%sfcwater_mass(ksn)  = y0%sfcwater_mass(ksn)   - wfreeb
               y0%sfcwater_depth(ksn) = y0%sfcwater_depth(ksn) - wfreeb * wdnsi8
               !----- Recompute the energy removing runoff --------------------------------!
               y0%sfcwater_energy(ksn) = y0%sfcwater_energy(ksn) - qwfree

               call redistribute_snow(y0,csite,ipa)
               call update_diagnostic_vars(y0,csite,ipa)

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
                  y0%wbudget_loss2runoff = wfreeb
                  y0%ebudget_loss2runoff = qwfree
                  y0%wbudget_storage     = y0%wbudget_storage - wfreeb
                  y0%ebudget_storage     = y0%ebudget_storage - qwfree
               end if

            else
               csite%runoff(ipa)          = 0.0
               csite%avg_runoff(ipa)      = 0.0
               csite%avg_runoff_heat(ipa) = 0.0
               y0%wbudget_loss2runoff  = 0.d0
               y0%ebudget_loss2runoff  = 0.d0
            end if
         else
            csite%runoff(ipa)          = 0.0
            csite%avg_runoff(ipa)      = 0.0
            csite%avg_runoff_heat(ipa) = 0.0
            y0%wbudget_loss2runoff  = 0.d0
            y0%ebudget_loss2runoff  = 0.d0
         end if

         !------ Copy the temporary patch to the next intermediate step -------------------!
         call copy_rk4_patch(y0,initp,cpatch)

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
   write (unit=*,fmt='(a)') ' ==> Too many steps in routine heun_integ'
   call print_rk4patch(ytry, csite,ipa)

   return
end subroutine heun_integ
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine heun_stepper(y0,ytry,yerr,ak1,ak2,yeuler,x,h,csite,ipa,reject_step,reject_result)

   use rk4_coms      , only : rk4patchtype        & ! structure
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
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype), target      :: y0
   type(rk4patchtype), target      :: ytry
   type(rk4patchtype), target      :: yerr
   type(rk4patchtype), target      :: ak1
   type(rk4patchtype), target      :: ak2
   type(rk4patchtype), target      :: yeuler
   type(sitetype)    , target      :: csite
   integer           , intent(in)  :: ipa
   logical           , intent(out) :: reject_step
   logical           , intent(out) :: reject_result
   real(kind=8)      , intent(in)  :: x
   real(kind=8)      , intent(in)  :: h
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)   , pointer     :: cpatch
   real(kind=8)                    :: combh
   !---------------------------------------------------------------------------------------!


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
   call copy_rk4_patch(y0,yeuler,cpatch)
   call inc_rk4_patch (yeuler,ak1, heun_b21*h, cpatch)
   combh = heun_b21*h
   call adjust_veg_properties    (yeuler,combh,csite,ipa)
   call adjust_topsoil_properties(yeuler,combh,csite,ipa)
   call redistribute_snow        (yeuler      ,csite,ipa)
   call update_diagnostic_vars   (yeuler      ,csite,ipa)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check to see if the Euler result makes sense.  Since we will compute the           !
   ! derivative correction using it, the Euler step must be bounded.  If not, reject the   !
   ! step and try a smaller step size.                                                     !
   !---------------------------------------------------------------------------------------!
   call rk4_sanity_check(yeuler, reject_step, csite, ipa,ak1,h,print_diags)
   if (reject_step) return
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute the second term (correction) of the derivative, using the Euler's         !
   ! predicted state.                                                                      !
   !---------------------------------------------------------------------------------------!
   call leaf_derivs(yeuler, ak2, csite,ipa)
   !---------------------------------------------------------------------------------------!



   !----- We now combine both derivatives and update y0. ----------------------------------!
   call copy_rk4_patch(y0, ytry, cpatch)
   call inc_rk4_patch(ytry,  ak1, heun_c1*h, cpatch)
   call inc_rk4_patch(ytry,  ak2, heun_c2*h, cpatch)
   combh = (heun_c1+heun_c2) * h ! Which should be h
   call adjust_veg_properties(ytry,combh,csite,ipa)
   call adjust_topsoil_properties(ytry,combh,csite,ipa)
   call redistribute_snow(ytry, csite,ipa)
   call update_diagnostic_vars(ytry, csite,ipa)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Check to see if this attempt of advancing one time step makes sense.  If not,      !
   ! reject the result and try a smaller step size.                                        !
   !---------------------------------------------------------------------------------------!
   call rk4_sanity_check(ytry, reject_result, csite,ipa,ak2,h,print_diags)
   if(reject_result)return
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Compute the estimate of the error associated with the step.                       !
   !---------------------------------------------------------------------------------------!
   call zero_rk4_patch (yerr)
   call zero_rk4_cohort(yerr)
   call inc_rk4_patch(yerr, ak1,  heun_dc1*h, cpatch)
   call inc_rk4_patch(yerr, ak2,  heun_dc2*h, cpatch)
   !---------------------------------------------------------------------------------------!

   return
end subroutine heun_stepper
!==========================================================================================!
!==========================================================================================!
