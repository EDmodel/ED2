module hybrid_driver
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is the main driver for the Forward/Backward (FB)                  !
   !     Euler integration scheme.                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine hybrid_timestep(cgrid)
     use rk4_coms              , only : integration_vars   & ! structure
                                      , rk4patchtype       & ! structure
                                      , zero_rk4_patch     & ! subroutine
                                      , zero_rk4_cohort    & ! subroutine
                                      , zero_bdf2_patch    &
                                      , integration_buff   & ! intent(out)
                                      , bdf2patchtype      &
                                      , tbeg               &
                                      , tend               &
                                      , dtrk4i             !
     use ed_para_coms          , only : nthreads           ! ! intent(in)
     use rk4_driver            , only : initp2modelp
     use ed_state_vars         , only : edtype             & ! structure
                                      , polygontype        & ! structure
                                      , sitetype           ! ! structure
     use met_driver_coms       , only : met_driv_state     ! ! structure
     use grid_coms             , only : nzg                & ! intent(in)
                                      , nzs                ! ! intent(in)
     use ed_misc_coms          , only : current_time       & ! intent(in)
                                      , dtlsm              ! ! intent(in)
     use therm_lib             , only : tq2enthalpy        ! ! function
     use budget_utils          , only : update_budget      & ! function
                                      , compute_budget     ! ! function
     use soil_respiration      , only : soil_respiration_driver    ! ! function
     use photosyn_driv         , only : canopy_photosynthesis      ! ! function
     use update_derived_utils  , only : update_patch_thermo_props  & ! subroutine
                                      , update_patch_derived_props ! ! subroutine
     use rk4_integ_utils       , only : copy_met_2_rk4site         ! ! subroutine
     use rk4_misc              , only : copy_patch_init            & ! subroutine
                                      , copy_patch_init_carbon     & ! subroutine
                                      , sanity_check_veg_energy    ! ! subroutine

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
     real                                   :: old_can_shv
     real                                   :: old_can_co2
     real                                   :: old_can_rhos
     real                                   :: old_can_temp
     real                                   :: old_can_prss
     real                                   :: old_can_enthalpy
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
           !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE( &
           !$OMP ita,ipa,initp,ytemp,dinitp,yprev,hbeg,nsteps,&
           !$OMP patch_vels,old_can_shv,                    &
           !$OMP old_can_co2,old_can_rhos,old_can_temp,     &
           !$OMP old_can_prss,old_can_enthalpy,             &
           !$OMP wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm,&
           !$OMP co2curr_loss2atm,wcurr_loss2drainage,      &
           !$OMP ecurr_loss2drainage,wcurr_loss2runoff,     &
           !$OMP ecurr_loss2runoff)
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

                 !------------------------------------------------------------------!
                 !    Update roughness and canopy depth.                            !
                 !------------------------------------------------------------------!
                 call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs,&
                      cpoly%ntext_soil(:,isi))
                 call update_patch_derived_props(csite,ipa)
                 !------------------------------------------------------------------!

                 !----- Save the previous thermodynamic state. ---------------------!
                 old_can_shv      = csite%can_shv(ipa)
                 old_can_co2      = csite%can_co2(ipa)
                 old_can_rhos     = csite%can_rhos(ipa)
                 old_can_temp     = csite%can_temp(ipa)
                 old_can_prss     = csite%can_prss (ipa)
                 old_can_enthalpy = tq2enthalpy(csite%can_temp(ipa)                 &
                                               ,csite%can_shv(ipa),.true.)
                 !------------------------------------------------------------------!


                 !----- Compute current storage terms. -----------------------------!
                 call update_budget(csite,cpoly%lsl(isi),ipa)
                 !------------------------------------------------------------------!



                 !------------------------------------------------------------------!
                 !      Test whether temperature and energy are reasonable.         !
                 !------------------------------------------------------------------!
                 if (test_energy_sanity) then
                    call sanity_check_veg_energy(csite,ipa)
                 end if
                 !------------------------------------------------------------------!



                 !------------------------------------------------------------------!
                 !     Set up the integration patch.                                !
                 !------------------------------------------------------------------!
                 call copy_patch_init(csite,ipa,ibuff,initp,patch_vels)

                 !------------------------------------------------------------------!
                 !     Set up the buffer for the previous step's leaf temperature   !
                 !------------------------------------------------------------------!
                 call copy_bdf2_prev(csite,ipa,yprev)

                 !----- Get photosynthesis, stomatal conductance,
                 !                                    and transpiration. -----------!
                 call canopy_photosynthesis(csite,cmet,nzg,ipa,ibuff,               &
                      cpoly%ntext_soil(:,isi),cpoly%leaf_aging_factor(:,isi),       &
                      cpoly%green_leaf_factor(:,isi))

                  !----- Compute root and heterotrophic respiration. ----------------!
                 call soil_respiration_driver(csite,ipa,nzg,cpoly%ntext_soil(:,isi))

                 !------------------------------------------------------------------!
                 ! Set up the remaining, carbon-dependent variables to the buffer.  !
                 !------------------------------------------------------------------!
                 call copy_patch_init_carbon(csite,ipa,initp)

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
                  initp%cpwp = initp%can_rhos * initp%cpwp * dtrk4i
                  initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i

                  !--------------------------------------------------------------------------!
                  ! Move the state variables from the integrated patch to the model patch.   !
                  !--------------------------------------------------------------------------!
                  call initp2modelp(tend-tbeg,initp,csite,ipa,cpoly%nighttime(isi)           &
                                   ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm               &
                                   ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage &
                                   ,wcurr_loss2runoff,ecurr_loss2runoff)


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
                  !     Compute the residuals.                                       !
                  !------------------------------------------------------------------!

                  call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg,ipa       &
                       ,wcurr_loss2atm,ecurr_netrad,ecurr_loss2atm                 &
                       ,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage   &
                       ,wcurr_loss2runoff,ecurr_loss2runoff,cpoly%area(isi)        &
                       ,cgrid%cbudget_nep(ipy),old_can_enthalpy,old_can_shv        &
                       ,old_can_co2,old_can_rhos,old_can_prss)
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

            !--------------------------------------------------------------------!
            !   Copy patch to the temporary structure                            !
            !   Note that this routine also calculates the size of the matrix    !
            !   used in the implicit step.                                       !
            !   nsolve = 1 + n_leaf_cohorts + n_wood_cohorts                     !
            !--------------------------------------------------------------------!
            call copy_fb_patch(initp,ytemp,cpatch)

            !--------------------------------------------------------------------!
            !   Integrate the forward step                                       !
            !--------------------------------------------------------------------!
            call inc_fwd_patch(ytemp,dinitp,h,cpatch)

            !--------------------------------------------------------------------!
            !   Integrate the implicit/backwards step                            !
            !--------------------------------------------------------------------!
            call bdf2_integ(ibuff,cpatch,yprev,initp,ytemp,dinitp,h, &
                 dble(csite%hprev(ipa)))

            !----- Perform a sanity check on canopy,leaf and wood stuff ---------!
            call fb_sanity_check(ytemp,reject_step,csite,ipa,dinitp,h,print_diags)

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

                  call fb_sanity_check(ytemp,test_reject,csite,ipa,dinitp,h,.true.)
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


    subroutine copy_fb_patch(sourcep, targetp, cpatch)

     use rk4_coms      , only : rk4site           & ! intent(in)
                               , rk4patchtype      & ! structure
                               , checkbudget       & ! intent(in)
                               , print_detailed    ! ! intent(in)
      use ed_state_vars , only : sitetype          & ! structure
                               , patchtype         ! ! structure
      use grid_coms     , only : nzg               & ! intent(in)
                               , nzs               ! ! intent(in)
      use ed_max_dims   , only : n_pft             ! ! intent(in)
      use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)

      implicit none
      !----- Arguments -----------------------------------------------------------------------!
      type(rk4patchtype) , target     :: sourcep
      type(rk4patchtype) , target     :: targetp
      type(patchtype)    , target     :: cpatch
      integer                         :: nsolve
      !----- Local variable ------------------------------------------------------------------!
      integer                         :: k
      !---------------------------------------------------------------------------------------!

      targetp%can_enthalpy     = sourcep%can_enthalpy
      targetp%can_theta        = sourcep%can_theta
      targetp%can_temp         = sourcep%can_temp
      targetp%can_shv          = sourcep%can_shv
      targetp%can_co2          = sourcep%can_co2
      targetp%can_rhos         = sourcep%can_rhos
      targetp%can_prss         = sourcep%can_prss
      targetp%can_exner        = sourcep%can_exner
      targetp%can_cp           = sourcep%can_cp
      targetp%can_depth        = sourcep%can_depth
      targetp%can_rhv          = sourcep%can_rhv
      targetp%can_ssh          = sourcep%can_ssh
      targetp%veg_height       = sourcep%veg_height
      targetp%veg_displace     = sourcep%veg_displace
      targetp%veg_rough        = sourcep%veg_rough
      targetp%opencan_frac     = sourcep%opencan_frac
      targetp%total_sfcw_depth = sourcep%total_sfcw_depth
      targetp%snowfac          = sourcep%snowfac

      targetp%atm_enthalpy     = sourcep%atm_enthalpy
      targetp%vels             = sourcep%vels


      targetp%ggbare           = sourcep%ggbare
      targetp%ggveg            = sourcep%ggveg
      targetp%ggnet            = sourcep%ggnet
      targetp%ggsoil           = sourcep%ggsoil

      targetp%flag_wflxgc      = sourcep%flag_wflxgc

      targetp%virtual_water    = sourcep%virtual_water
      targetp%virtual_energy   = sourcep%virtual_energy
      targetp%virtual_depth    = sourcep%virtual_depth
      targetp%virtual_tempk    = sourcep%virtual_tempk
      targetp%virtual_fracliq  = sourcep%virtual_fracliq

      targetp%rough            = sourcep%rough

      targetp%upwp             = sourcep%upwp
      targetp%wpwp             = sourcep%wpwp
      targetp%tpwp             = sourcep%tpwp
      targetp%qpwp             = sourcep%qpwp
      targetp%cpwp             = sourcep%cpwp

      targetp%ground_shv       = sourcep%ground_shv
      targetp%ground_ssh       = sourcep%ground_ssh
      targetp%ground_temp      = sourcep%ground_temp
      targetp%ground_fliq      = sourcep%ground_fliq

      targetp%ustar            = sourcep%ustar
      targetp%cstar            = sourcep%cstar
      targetp%tstar            = sourcep%tstar
      targetp%estar            = sourcep%estar
      targetp%qstar            = sourcep%qstar
      targetp%zeta             = sourcep%zeta
      targetp%ribulk           = sourcep%ribulk
      targetp%rasveg           = sourcep%rasveg

      targetp%fgc_rh           = sourcep%fgc_rh
      targetp%fsc_rh           = sourcep%fsc_rh
      targetp%stgc_rh          = sourcep%stgc_rh
      targetp%stsc_rh          = sourcep%stsc_rh
      targetp%msc_rh           = sourcep%msc_rh
      targetp%ssc_rh           = sourcep%ssc_rh

      targetp%water_deficit    = sourcep%water_deficit


      do k=rk4site%lsl,nzg
         targetp%soil_water            (k) = sourcep%soil_water            (k)
         targetp%soil_mstpot           (k) = sourcep%soil_mstpot           (k)
         targetp%soil_energy           (k) = sourcep%soil_energy           (k)
         targetp%soil_tempk            (k) = sourcep%soil_tempk            (k)
         targetp%soil_fracliq          (k) = sourcep%soil_fracliq          (k)
      end do

      targetp%nlev_sfcwater   = sourcep%nlev_sfcwater
      targetp%flag_sfcwater   = sourcep%flag_sfcwater

      do k=1,nzs
         targetp%sfcwater_mass     (k) = sourcep%sfcwater_mass     (k)
         targetp%sfcwater_energy   (k) = sourcep%sfcwater_energy   (k)
         targetp%sfcwater_depth    (k) = sourcep%sfcwater_depth    (k)
         targetp%sfcwater_tempk    (k) = sourcep%sfcwater_tempk    (k)
         targetp%sfcwater_fracliq  (k) = sourcep%sfcwater_fracliq  (k)
      end do

      nsolve=1

      do k=1,cpatch%ncohorts
         targetp%leaf_resolvable   (k) = sourcep%leaf_resolvable   (k)
         targetp%wood_resolvable   (k) = sourcep%wood_resolvable   (k)

         targetp%leaf_energy       (k) = sourcep%leaf_energy       (k)
         targetp%leaf_water        (k) = sourcep%leaf_water        (k)
         targetp%leaf_temp         (k) = sourcep%leaf_temp         (k)
         targetp%leaf_fliq         (k) = sourcep%leaf_fliq         (k)
         targetp%leaf_hcap         (k) = sourcep%leaf_hcap         (k)
         targetp%leaf_reynolds     (k) = sourcep%leaf_reynolds     (k)
         targetp%leaf_grashof      (k) = sourcep%leaf_grashof      (k)
         targetp%leaf_nussfree     (k) = sourcep%leaf_nussfree     (k)
         targetp%leaf_nussforc     (k) = sourcep%leaf_nussforc     (k)
         targetp%leaf_gbh          (k) = sourcep%leaf_gbh          (k)
         targetp%leaf_gbw          (k) = sourcep%leaf_gbw          (k)
         targetp%rshort_l          (k) = sourcep%rshort_l          (k)
         targetp%rlong_l           (k) = sourcep%rlong_l           (k)

         targetp%wood_resolvable   (k) = sourcep%wood_resolvable   (k)

         targetp%wood_energy       (k) = sourcep%wood_energy       (k)
         targetp%wood_water        (k) = sourcep%wood_water        (k)
         targetp%wood_temp         (k) = sourcep%wood_temp         (k)
         targetp%wood_fliq         (k) = sourcep%wood_fliq         (k)
         targetp%wood_hcap         (k) = sourcep%wood_hcap         (k)
         targetp%wood_reynolds     (k) = sourcep%wood_reynolds     (k)
         targetp%wood_grashof      (k) = sourcep%wood_grashof      (k)
         targetp%wood_nussfree     (k) = sourcep%wood_nussfree     (k)
         targetp%wood_nussforc     (k) = sourcep%wood_nussforc     (k)
         targetp%wood_gbh          (k) = sourcep%wood_gbh          (k)
         targetp%wood_gbw          (k) = sourcep%wood_gbw          (k)
         targetp%rshort_w          (k) = sourcep%rshort_w          (k)
         targetp%rlong_w           (k) = sourcep%rlong_w           (k)

         targetp%veg_resolvable    (k) = sourcep%veg_resolvable    (k)
         targetp%veg_energy        (k) = sourcep%veg_energy        (k)
         targetp%veg_water         (k) = sourcep%veg_water         (k)
         targetp%veg_hcap          (k) = sourcep%veg_hcap          (k)

         targetp%veg_wind          (k) = sourcep%veg_wind          (k)
         targetp%lint_shv          (k) = sourcep%lint_shv          (k)
         targetp%nplant            (k) = sourcep%nplant            (k)
         targetp%lai               (k) = sourcep%lai               (k)
         targetp%wai               (k) = sourcep%wai               (k)
         targetp%tai               (k) = sourcep%tai               (k)
         targetp%crown_area        (k) = sourcep%crown_area        (k)
         targetp%elongf            (k) = sourcep%elongf            (k)
         targetp%gsw_open          (k) = sourcep%gsw_open          (k)
         targetp%gsw_closed        (k) = sourcep%gsw_closed        (k)
         targetp%psi_open          (k) = sourcep%psi_open          (k)
         targetp%psi_closed        (k) = sourcep%psi_closed        (k)
         targetp%fs_open           (k) = sourcep%fs_open           (k)
         targetp%gpp               (k) = sourcep%gpp               (k)
         targetp%leaf_resp         (k) = sourcep%leaf_resp         (k)
         targetp%root_resp         (k) = sourcep%root_resp         (k)
         targetp%leaf_growth_resp  (k) = sourcep%leaf_growth_resp  (k)
         targetp%root_growth_resp  (k) = sourcep%root_growth_resp  (k)
         targetp%sapa_growth_resp  (k) = sourcep%sapa_growth_resp  (k)
         targetp%sapb_growth_resp  (k) = sourcep%sapb_growth_resp  (k)
         targetp%barka_growth_resp (k) = sourcep%barka_growth_resp (k)
         targetp%barkb_growth_resp (k) = sourcep%barkb_growth_resp (k)
         targetp%leaf_storage_resp (k) = sourcep%leaf_storage_resp (k)
         targetp%root_storage_resp (k) = sourcep%root_storage_resp (k)
         targetp%sapa_storage_resp (k) = sourcep%sapa_storage_resp (k)
         targetp%sapb_storage_resp (k) = sourcep%sapb_storage_resp (k)
         targetp%barka_storage_resp(k) = sourcep%barka_storage_resp(k)
         targetp%barkb_storage_resp(k) = sourcep%barkb_storage_resp(k)
      end do

      if (checkbudget) then
         targetp%co2budget_storage      = sourcep%co2budget_storage
         targetp%co2budget_loss2atm     = sourcep%co2budget_loss2atm
         targetp%ebudget_netrad         = sourcep%ebudget_netrad
         targetp%ebudget_loss2atm       = sourcep%ebudget_loss2atm
         targetp%ebudget_loss2drainage  = sourcep%ebudget_loss2drainage
         targetp%ebudget_loss2runoff    = sourcep%ebudget_loss2runoff
         targetp%wbudget_loss2atm       = sourcep%wbudget_loss2atm
         targetp%wbudget_loss2drainage  = sourcep%wbudget_loss2drainage
         targetp%wbudget_loss2runoff    = sourcep%wbudget_loss2runoff
         targetp%ebudget_storage        = sourcep%ebudget_storage
         targetp%wbudget_storage        = sourcep%wbudget_storage
      end if
      if (fast_diagnostics) then
         targetp%avg_ustar              = sourcep%avg_ustar
         targetp%avg_tstar              = sourcep%avg_tstar
         targetp%avg_qstar              = sourcep%avg_qstar
         targetp%avg_cstar              = sourcep%avg_cstar
         targetp%avg_carbon_ac          = sourcep%avg_carbon_ac
         targetp%avg_carbon_st          = sourcep%avg_carbon_st
         targetp%avg_vapor_gc           = sourcep%avg_vapor_gc
         targetp%avg_throughfall        = sourcep%avg_throughfall
         targetp%avg_vapor_ac           = sourcep%avg_vapor_ac
         targetp%avg_qthroughfall       = sourcep%avg_qthroughfall
         targetp%avg_sensible_gc        = sourcep%avg_sensible_gc
         targetp%avg_sensible_ac        = sourcep%avg_sensible_ac
         targetp%avg_drainage           = sourcep%avg_drainage
         targetp%avg_qdrainage          = sourcep%avg_qdrainage


         do k=rk4site%lsl,nzg
            targetp%avg_sensible_gg(k) = sourcep%avg_sensible_gg(k)
            targetp%avg_smoist_gg(k)   = sourcep%avg_smoist_gg(k)
            targetp%avg_transloss(k)   = sourcep%avg_transloss(k)
         end do


         do k=1,cpatch%ncohorts
            targetp%avg_sensible_lc    (k) = sourcep%avg_sensible_lc   (k)
            targetp%avg_sensible_wc    (k) = sourcep%avg_sensible_wc   (k)
            targetp%avg_vapor_lc       (k) = sourcep%avg_vapor_lc      (k)
            targetp%avg_vapor_wc       (k) = sourcep%avg_vapor_wc      (k)
            targetp%avg_transp         (k) = sourcep%avg_transp        (k)
            targetp%avg_intercepted_al (k) = sourcep%avg_intercepted_al(k)
            targetp%avg_intercepted_aw (k) = sourcep%avg_intercepted_aw(k)
            targetp%avg_wshed_lg       (k) = sourcep%avg_wshed_lg      (k)
            targetp%avg_wshed_wg       (k) = sourcep%avg_wshed_wg      (k)
         end do
      end if

      if (print_detailed) then
         targetp%flx_carbon_ac          = sourcep%flx_carbon_ac
         targetp%flx_carbon_st          = sourcep%flx_carbon_st
         targetp%flx_vapor_lc           = sourcep%flx_vapor_lc
         targetp%flx_vapor_wc           = sourcep%flx_vapor_wc
         targetp%flx_vapor_gc           = sourcep%flx_vapor_gc
         targetp%flx_wshed_vg           = sourcep%flx_wshed_vg
         targetp%flx_intercepted        = sourcep%flx_intercepted
         targetp%flx_throughfall        = sourcep%flx_throughfall
         targetp%flx_vapor_ac           = sourcep%flx_vapor_ac
         targetp%flx_transp             = sourcep%flx_transp
         targetp%flx_rshort_gnd         = sourcep%flx_rshort_gnd
         targetp%flx_par_gnd            = sourcep%flx_par_gnd
         targetp%flx_rlong_gnd          = sourcep%flx_rlong_gnd
         targetp%flx_sensible_lc        = sourcep%flx_sensible_lc
         targetp%flx_sensible_wc        = sourcep%flx_sensible_wc
         targetp%flx_qwshed_vg          = sourcep%flx_qwshed_vg
         targetp%flx_qintercepted       = sourcep%flx_qintercepted
         targetp%flx_qthroughfall       = sourcep%flx_qthroughfall
         targetp%flx_sensible_gc        = sourcep%flx_sensible_gc
         targetp%flx_sensible_ac        = sourcep%flx_sensible_ac
         targetp%flx_drainage           = sourcep%flx_drainage
         targetp%flx_qdrainage          = sourcep%flx_qdrainage

         do k=rk4site%lsl,nzg
            targetp%flx_sensible_gg(k) = sourcep%flx_sensible_gg(k)
            targetp%flx_smoist_gg(k)   = sourcep%flx_smoist_gg(k)
            targetp%flx_transloss(k)   = sourcep%flx_transloss(k)
         end do

         do k=1,cpatch%ncohorts
            targetp%cfx_hflxlc      (k) = sourcep%cfx_hflxlc      (k)
            targetp%cfx_hflxwc      (k) = sourcep%cfx_hflxwc      (k)
            targetp%cfx_qwflxlc     (k) = sourcep%cfx_qwflxlc     (k)
            targetp%cfx_qwflxwc     (k) = sourcep%cfx_qwflxwc     (k)
            targetp%cfx_qwshed      (k) = sourcep%cfx_qwshed      (k)
            targetp%cfx_qtransp     (k) = sourcep%cfx_qtransp     (k)
            targetp%cfx_qintercepted(k) = sourcep%cfx_qintercepted(k)
         end do
      end if

    end subroutine copy_fb_patch


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


    subroutine inc_fwd_patch(rkp, inc, fac, cpatch)
      use ed_state_vars , only : sitetype           & ! structure
                               , patchtype          ! ! structure
      use rk4_coms      , only : rk4patchtype       & ! structure
                               , rk4site            & ! intent(in)
                               , checkbudget        & ! intent(in)
                               , print_detailed     ! ! intent(in)
      use grid_coms     , only : nzg                ! ! intent(in)
      use ed_misc_coms  , only : fast_diagnostics   ! ! intent(in)

      implicit none

      !----- Arguments -----------------------------------------------------------------------!
      type(rk4patchtype) , target     :: rkp    ! Temporary patch with previous state
      type(rk4patchtype) , target     :: inc    ! Temporary patch with its derivatives
      type(patchtype)    , target     :: cpatch ! Current patch (for characteristics)
      real(kind=8)       , intent(in) :: fac    ! Increment factor
      !----- Local variables -----------------------------------------------------------------!
      integer                         :: ico    ! Cohort ID
      integer                         :: k      ! Counter
      !---------------------------------------------------------------------------------------!




      rkp%can_enthalpy = rkp%can_enthalpy + fac * inc%can_enthalpy
      rkp%can_shv      = rkp%can_shv      + fac * inc%can_shv
      rkp%can_co2      = rkp%can_co2      + fac * inc%can_co2

      do k=rk4site%lsl,nzg
         rkp%soil_water(k)       = rkp%soil_water(k)  + fac * inc%soil_water(k)
         rkp%soil_energy(k)      = rkp%soil_energy(k) + fac * inc%soil_energy(k)
      end do

      do k=1,rkp%nlev_sfcwater
         rkp%sfcwater_mass(k)   = rkp%sfcwater_mass(k)   + fac * inc%sfcwater_mass(k)
         rkp%sfcwater_energy(k) = rkp%sfcwater_energy(k) + fac * inc%sfcwater_energy(k)
         rkp%sfcwater_depth(k)  = rkp%sfcwater_depth(k)  + fac * inc%sfcwater_depth(k)
      end do

      rkp%virtual_energy  = rkp%virtual_energy  + fac * inc%virtual_energy
      rkp%virtual_water   = rkp%virtual_water   + fac * inc%virtual_water
      rkp%virtual_depth   = rkp%virtual_depth   + fac * inc%virtual_depth


      rkp%upwp = rkp%upwp + fac * inc%upwp
      rkp%wpwp = rkp%wpwp + fac * inc%wpwp
      rkp%tpwp = rkp%tpwp + fac * inc%tpwp
      rkp%qpwp = rkp%qpwp + fac * inc%qpwp
      rkp%cpwp = rkp%cpwp + fac * inc%cpwp

      rkp%water_deficit   = rkp%water_deficit      + fac * inc%water_deficit

      do ico = 1,cpatch%ncohorts
         rkp%leaf_water (ico) = rkp%leaf_water (ico) + fac * inc%leaf_water (ico)
         rkp%leaf_energy(ico) = rkp%leaf_energy(ico) + fac * inc%leaf_energy(ico)
         rkp%wood_water (ico) = rkp%wood_water (ico) + fac * inc%wood_water (ico)
         rkp%wood_energy(ico) = rkp%wood_energy(ico) + fac * inc%wood_energy(ico)
         rkp%veg_water (ico)  = rkp%veg_water  (ico) + fac * inc%veg_water  (ico)
         rkp%veg_energy(ico)  = rkp%veg_energy (ico) + fac * inc%veg_energy (ico)

         rkp%psi_open  (ico) = rkp%psi_open  (ico) + fac * inc%psi_open  (ico)
         rkp%psi_closed(ico) = rkp%psi_closed(ico) + fac * inc%psi_closed(ico)
      end do

      if (checkbudget) then

         rkp%co2budget_storage      = rkp%co2budget_storage     + fac * inc%co2budget_storage
         rkp%co2budget_loss2atm     = rkp%co2budget_loss2atm    + fac * inc%co2budget_loss2atm

         rkp%wbudget_storage       = rkp%wbudget_storage       + fac * inc%wbudget_storage
         rkp%wbudget_loss2atm      = rkp%wbudget_loss2atm      + fac * inc%wbudget_loss2atm
         rkp%wbudget_loss2drainage = rkp%wbudget_loss2drainage                                &
                                   + fac * inc%wbudget_loss2drainage

         rkp%ebudget_storage       = rkp%ebudget_storage       + fac * inc%ebudget_storage
         rkp%ebudget_netrad        = rkp%ebudget_netrad        + fac * inc%ebudget_netrad
         rkp%ebudget_loss2atm      = rkp%ebudget_loss2atm      + fac * inc%ebudget_loss2atm
         rkp%ebudget_loss2drainage = rkp%ebudget_loss2drainage                                &
                                   + fac * inc%ebudget_loss2drainage
      end if
      if (fast_diagnostics) then
         rkp%avg_ustar          = rkp%avg_ustar          + fac * inc%avg_ustar
         rkp%avg_tstar          = rkp%avg_tstar          + fac * inc%avg_tstar
         rkp%avg_qstar          = rkp%avg_qstar          + fac * inc%avg_qstar
         rkp%avg_cstar          = rkp%avg_cstar          + fac * inc%avg_cstar


         rkp%avg_carbon_ac      = rkp%avg_carbon_ac      + fac * inc%avg_carbon_ac
         rkp%avg_carbon_st      = rkp%avg_carbon_st      + fac * inc%avg_carbon_st

         rkp%avg_vapor_gc       = rkp%avg_vapor_gc       + fac * inc%avg_vapor_gc
         rkp%avg_throughfall    = rkp%avg_throughfall    + fac * inc%avg_throughfall
         rkp%avg_vapor_ac       = rkp%avg_vapor_ac       + fac * inc%avg_vapor_ac
         rkp%avg_drainage       = rkp%avg_drainage       + fac * inc%avg_drainage
         rkp%avg_qdrainage      = rkp%avg_qdrainage      + fac * inc%avg_qdrainage
         rkp%avg_qthroughfall   = rkp%avg_qthroughfall   + fac * inc%avg_qthroughfall
         rkp%avg_sensible_gc    = rkp%avg_sensible_gc    + fac * inc%avg_sensible_gc
         rkp%avg_sensible_ac    = rkp%avg_sensible_ac    + fac * inc%avg_sensible_ac


         do k=rk4site%lsl,nzg
            rkp%avg_sensible_gg(k)  = rkp%avg_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
            rkp%avg_smoist_gg(k)    = rkp%avg_smoist_gg(k)    + fac * inc%avg_smoist_gg(k)
            rkp%avg_transloss(k)    = rkp%avg_transloss(k)    + fac * inc%avg_transloss(k)
         end do


         do k=1,cpatch%ncohorts
            rkp%avg_sensible_lc   (k) =       rkp%avg_sensible_lc   (k)                       &
                                      + fac * inc%avg_sensible_lc   (k)
            rkp%avg_sensible_wc   (k) =       rkp%avg_sensible_wc   (k)                       &
                                      + fac * inc%avg_sensible_wc   (k)
            rkp%avg_vapor_lc      (k) =       rkp%avg_vapor_lc      (k)                       &
                                      + fac * inc%avg_vapor_lc      (k)
            rkp%avg_vapor_wc      (k) =       rkp%avg_vapor_wc      (k)                       &
                                      + fac * inc%avg_vapor_wc      (k)
            rkp%avg_transp        (k) =       rkp%avg_transp        (k)                       &
                                      + fac * inc%avg_transp        (k)
            rkp%avg_intercepted_al(k) =       rkp%avg_intercepted_al(k)                       &
                                      + fac * inc%avg_intercepted_al(k)
            rkp%avg_intercepted_aw(k) =       rkp%avg_intercepted_aw(k)                       &
                                      + fac * inc%avg_intercepted_aw(k)
            rkp%avg_wshed_lg      (k) =       rkp%avg_wshed_lg      (k)                       &
                                      + fac * inc%avg_wshed_lg      (k)
            rkp%avg_wshed_wg      (k) =       rkp%avg_wshed_wg      (k)                       &
                                      + fac * inc%avg_wshed_wg      (k)
         end do

      end if

      !---------------------------------------------------------------------------------------!
      !    Increment the instantaneous fluxes.  The derivative term should be the same as the !
      ! the full fluxes, the only difference is that these variables are normalised and       !
      ! re-set after each time step.                                                          !
      !---------------------------------------------------------------------------------------!
      if (print_detailed) then
         rkp%flx_carbon_ac      = rkp%flx_carbon_ac      + fac * inc%avg_carbon_ac
         rkp%flx_carbon_st      = rkp%flx_carbon_st      + fac * inc%avg_carbon_st

         rkp%flx_vapor_gc       = rkp%flx_vapor_gc       + fac * inc%avg_vapor_gc
         rkp%flx_throughfall    = rkp%flx_throughfall    + fac * inc%avg_throughfall
         rkp%flx_vapor_ac       = rkp%flx_vapor_ac       + fac * inc%avg_vapor_ac
         rkp%flx_drainage       = rkp%flx_drainage       + fac * inc%avg_drainage
         rkp%flx_qdrainage      = rkp%flx_qdrainage      + fac * inc%avg_qdrainage
         rkp%flx_qthroughfall   = rkp%flx_qthroughfall   + fac * inc%avg_qthroughfall
         rkp%flx_sensible_gc    = rkp%flx_sensible_gc    + fac * inc%avg_sensible_gc
         rkp%flx_sensible_ac    = rkp%flx_sensible_ac    + fac * inc%avg_sensible_ac

         do k=rk4site%lsl,nzg
            rkp%flx_sensible_gg(k)  = rkp%flx_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
            rkp%flx_smoist_gg(k)    = rkp%flx_smoist_gg(k)    + fac * inc%avg_smoist_gg(k)
            rkp%flx_transloss(k)    = rkp%flx_transloss(k)    + fac * inc%avg_transloss(k)
         end do

         do ico = 1,cpatch%ncohorts
            rkp%flx_vapor_lc           =         rkp%flx_vapor_lc                             &
                                       + fac *   inc%avg_vapor_lc       (ico)
            rkp%flx_vapor_wc           =         rkp%flx_vapor_wc                             &
                                       + fac *   inc%avg_vapor_wc       (ico)
            rkp%flx_wshed_vg           =         rkp%flx_wshed_vg                             &
                                       + fac * ( inc%avg_wshed_lg       (ico)                 &
                                               + inc%avg_wshed_wg       (ico) )
            rkp%flx_intercepted        =       rkp%flx_intercepted                            &
                                       + fac * ( inc%avg_intercepted_al (ico)                 &
                                               + inc%avg_intercepted_aw (ico) )
            rkp%flx_sensible_lc        =         rkp%flx_sensible_lc                          &
                                       + fac *   inc%avg_sensible_lc    (ico)
            rkp%flx_sensible_wc        =         rkp%flx_sensible_wc                          &
                                       + fac *   inc%avg_sensible_wc    (ico)
            rkp%flx_qwshed_vg          =         rkp%flx_qwshed_vg                            &
                                       + fac *   inc%cfx_qwshed         (ico)
            rkp%flx_qintercepted       =         rkp%flx_qintercepted                         &
                                       + fac *   inc%cfx_qintercepted   (ico)
            rkp%cfx_hflxlc      (ico)  =         rkp%cfx_hflxlc         (ico)                 &
                                       + fac *   inc%cfx_hflxlc         (ico)
            rkp%cfx_hflxwc      (ico)  =         rkp%cfx_hflxwc         (ico)                 &
                                       + fac *   inc%cfx_hflxwc         (ico)
            rkp%cfx_qwflxlc     (ico)  =         rkp%cfx_qwflxlc        (ico)                 &
                                       + fac *   inc%cfx_qwflxlc        (ico)
            rkp%cfx_qwflxwc     (ico)  =         rkp%cfx_qwflxwc        (ico)                 &
                                       + fac *   inc%cfx_qwflxwc        (ico)
            rkp%cfx_qwshed      (ico)  =         rkp%cfx_qwshed         (ico)                 &
                                       + fac *   inc%cfx_qwshed         (ico)
            rkp%cfx_qtransp     (ico)  =         rkp%cfx_qtransp        (ico)                 &
                                       + fac *   inc%cfx_qtransp        (ico)
            rkp%cfx_qintercepted(ico)  =         rkp%cfx_qintercepted   (ico)                 &
                                       + fac *   inc%cfx_qintercepted   (ico)
         end do
      end if

      !---------------------------------------------------------------------------------------!

      return
    end subroutine inc_fwd_patch

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

    subroutine fb_sanity_check(y,reject_step, csite,ipa,dydx,h, &
         print_problems)
      use rk4_coms              , only : rk4patchtype          & ! structure
                                       , integration_vars      & ! structure
                                       , rk4site               & ! intent(in)
                                       , rk4aux                & ! intent(in)
                                       , rk4min_can_shv        & ! intent(in)
                                       , rk4min_can_temp       & ! intent(in)
                                       , rk4max_can_temp       & ! intent(in)
                                       , rk4min_can_co2        & ! intent(in)
                                       , rk4max_can_co2        & ! intent(in)
                                       , rk4max_veg_temp       & ! intent(in)
                                       , rk4min_veg_temp       & ! intent(in)
                                       , rk4min_veg_lwater     & ! intent(in)
                                       , rk4min_sfcw_temp      & ! intent(in)
                                       , rk4max_sfcw_temp      & ! intent(in)
                                       , rk4max_soil_temp      & ! intent(in)
                                       , rk4min_soil_temp      & ! intent(in)
                                       , rk4min_sfcw_mass      & ! intent(in)
                                       , rk4min_virt_water     & ! intent(in)
                                       , rk4water_stab_thresh  & ! intent(in)
                                       , integ_err             & ! intent(inout)
                                       , record_err            & ! intent(in)
                                       , osow                  & ! intent(in)
                                       , osoe                  & ! intent(in)
                                       , oswe                  & ! intent(in)
                                       , oswm                  ! ! intent(in)
      use ed_state_vars         , only : sitetype              & ! structure
                                       , patchtype             ! ! structure
      use grid_coms             , only : nzg                   ! ! intent(in)
      use therm_lib8            , only : eslif8
      use consts_coms           , only : ep8
      !$ use omp_lib

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target      :: y
      type(rk4patchtype) , target      :: dydx
      type(sitetype)     , target      :: csite
      logical            , intent(in)  :: print_problems
      logical            , intent(out) :: reject_step
      real(kind=8)       , intent(in)  :: h
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer     :: cpatch
      integer                          :: k
      integer                          :: ksn
      real(kind=8)                     :: rk4min_leaf_water
      real(kind=8)                     :: rk4min_wood_water
      real(kind=8)                     :: fbmax_can_shv
      integer                          :: ipa
      integer                          :: ico
      logical                          :: cflag7
      logical                          :: cflag8
      logical                          :: cflag9
      logical                          :: cflag10
      integer                          :: ibuff

      ibuff = 1
      !$ ibuff = OMP_get_thread_num()+1
      !------------------------------------------------------------------------------------!

      !----- Be optimistic and start assuming that things are fine. -----------------------!
      reject_step = .false.
      !------------------------------------------------------------------------------------!

      fbmax_can_shv = ep8*eslif8(320.d0)/y%can_prss


      if ( y%can_shv > fbmax_can_shv .or. y%can_shv < rk4min_can_shv  ) then
         reject_step = .true.
         if(record_err) integ_err(3,2) = integ_err(3,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Canopy air sp. humidity is off-track...'
               write(unit=*,fmt='(a)')           '-------------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP )/Dt:',dydx%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' '
            elseif (.not. record_err) then
               return
            end if
         end if
         !------------------------------------------------------------------------------------!



         !------------------------------------------------------------------------------------!
         !   Check whether the canopy air temperature is off.                                 !
         !------------------------------------------------------------------------------------!
         if (y%can_temp > rk4max_can_temp .or. y%can_temp < rk4min_can_temp) then
            reject_step = .true.
            if(record_err) integ_err(4,2) = integ_err(4,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' + Canopy air temperature is off-track...'
               write(unit=*,fmt='(a)')           '-------------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP )/Dt:',dydx%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' '
            elseif (.not. record_err) then
               return
            end if
         end if
         !------------------------------------------------------------------------------------!



         !------------------------------------------------------------------------------------!
         !   Check whether the canopy air pressure is off.                                    !
         !------------------------------------------------------------------------------------!
      if (y%can_prss > rk4aux(ibuff)%rk4max_can_prss .or. y%can_prss < rk4aux(ibuff)%rk4min_can_prss) then
            reject_step = .true.
            if(record_err) integ_err(5,2) = integ_err(5,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' + Canopy air pressure is off-track...'
               write(unit=*,fmt='(a)')           '-------------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP )/Dt:',dydx%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' '
            elseif (.not. record_err) then
               return
            end if
         end if
         !------------------------------------------------------------------------------------!

         !------------------------------------------------------------------------------------!
         !   Check whether the canopy air co2 is off.            !
         !------------------------------------------------------------------------------------!
         if (y%can_co2 > rk4max_can_co2 .or. y%can_co2 < rk4min_can_co2) then
            reject_step = .true.
            if(record_err) integ_err(6,2) = integ_err(6,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' + Canopy air CO2 is off-track...'
               write(unit=*,fmt='(a)')           '-------------------------------------------'
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
               write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_TEMP )/Dt:',dydx%can_theta
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
               write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
               write(unit=*,fmt='(a)')           '==========================================='
               write(unit=*,fmt='(a)')           ' '
               stop
            elseif (.not. record_err) then
               return
            end if
         end if
         !------------------------------------------------------------------------------------!



         !------------------------------------------------------------------------------------!
         !     Check leaf properties, but only for those cohorts with sufficient LAI.         !
         !------------------------------------------------------------------------------------!
         cpatch => csite%patch(ipa)
         cflag7 = .false.
         cflag8 = .false.
         leafloop: do ico = 1,cpatch%ncohorts
            if (.not. y%leaf_resolvable(ico)) cycle leafloop

            !----- Find the minimum leaf surface water. --------------------------------------!
            rk4min_leaf_water = rk4min_veg_lwater * y%lai(ico)

            !----- Check leaf surface water. -------------------------------------------------!
            if (y%leaf_water(ico) < rk4min_leaf_water) then
               reject_step = .true.
               if(record_err) cflag7 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a)')           ' + Leaf surface water is off-track...'
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_HCAP:     ',y%leaf_hcap(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_TEMP:     ',y%leaf_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_FRACLIQ:  ',y%leaf_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_ENERGY:   ',y%leaf_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER:    ',y%leaf_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' MIN_LEAF_WATER:',rk4min_leaf_water
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBH:      ',y%leaf_gbh(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBW:      ',y%leaf_gbw(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_REYNOLDS: ',y%leaf_reynolds(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GRASHOF:  ',y%leaf_grashof(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFREE:   ',y%leaf_nussfree(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFORC:   ',y%leaf_nussforc(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_EN)/Dt: ',dydx%leaf_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_WAT)/Dt:',dydx%leaf_water(ico)
                  write(unit=*,fmt='(a)')           '========================================'
               elseif (.not. record_err) then
                  return
               end if
            end if

            !----- Check leaf temperature. ---------------------------------------------------!
            if (y%leaf_temp(ico) > rk4max_veg_temp .or.                                       &
                y%leaf_temp(ico) < rk4min_veg_temp      ) then
               reject_step = .true.
               if(record_err) cflag8 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a)')           ' + Leaf temperature is off-track...'
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_HCAP:     ',y%leaf_hcap(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_TEMP:     ',y%leaf_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_FRACLIQ:  ',y%leaf_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_ENERGY:   ',y%leaf_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER:    ',y%leaf_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' MIN_LEAF_WATER:',rk4min_leaf_water
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBH:      ',y%leaf_gbh(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBW:      ',y%leaf_gbw(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_REYNOLDS: ',y%leaf_reynolds(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GRASHOF:  ',y%leaf_grashof(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFREE:   ',y%leaf_nussfree(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFORC:   ',y%leaf_nussforc(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_EN)/Dt: ',dydx%leaf_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_WAT)/Dt:',dydx%leaf_water(ico)
                  write(unit=*,fmt='(a)')           '========================================'
               elseif (.not. record_err) then
                  return
               end if
            end if
         end do leafloop
         if(record_err .and. cflag7) integ_err(7,2) = integ_err(7,2) + 1_8
         if(record_err .and. cflag8) integ_err(8,2) = integ_err(8,2) + 1_8
         !------------------------------------------------------------------------------------!




         !------------------------------------------------------------------------------------!
         !     Check wood properties, but only for those cohorts with sufficient LAI.         !
         !------------------------------------------------------------------------------------!
         cpatch => csite%patch(ipa)
         cflag9  = .false.
         cflag10 = .false.
         woodloop: do ico = 1,cpatch%ncohorts
            if (.not. y%wood_resolvable(ico)) cycle woodloop

            !----- Find the minimum wood surface water. --------------------------------------!
            rk4min_wood_water = rk4min_veg_lwater * y%wai(ico)

            !----- Check wood surface water. -------------------------------------------------!
            if (y%wood_water(ico) < rk4min_wood_water) then
               reject_step = .true.
               if(record_err) cflag9 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a)')           ' + Wood surface water is off-track...'
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_HCAP:     ',y%wood_hcap(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_TEMP:     ',y%wood_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_FRACLIQ:  ',y%wood_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_ENERGY:   ',y%wood_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER:    ',y%wood_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' MIN_WOOD_WATER:',rk4min_wood_water
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBH:      ',y%wood_gbh(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBW:      ',y%wood_gbw(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_REYNOLDS: ',y%wood_reynolds(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GRASHOF:  ',y%wood_grashof(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFREE:   ',y%wood_nussfree(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFORC:   ',y%wood_nussforc(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_EN)/Dt: ',dydx%wood_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_WAT)/Dt:',dydx%wood_water(ico)
                  write(unit=*,fmt='(a)')           '========================================'
               elseif (.not. record_err) then
                  return
               end if
            end if

            !----- Check wood temperature. ---------------------------------------------------!
            if (y%wood_temp(ico) > rk4max_veg_temp .or.                                       &
                y%wood_temp(ico) < rk4min_veg_temp      ) then
               reject_step = .true.
               if(record_err) cflag10 = .true.
               if (print_problems) then
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a)')           ' + Wood temperature is off-track...'
                  write(unit=*,fmt='(a)')           '========================================'
                  write(unit=*,fmt='(a,1x,i6)')     ' PFT:           ',cpatch%pft(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' HEIGHT:        ',cpatch%hite(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LAI:           ',y%lai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WAI:           ',y%wai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' TAI:           ',y%tai(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' NPLANT:        ',y%nplant(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' CROWN_AREA:    ',y%crown_area(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_HCAP:     ',y%wood_hcap(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_TEMP:     ',y%wood_temp(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_FRACLIQ:  ',y%wood_fliq(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_ENERGY:   ',y%wood_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER:    ',y%wood_water(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' MIN_WOOD_WATER:',rk4min_wood_water
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBH:      ',y%wood_gbh(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBW:      ',y%wood_gbw(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_REYNOLDS: ',y%wood_reynolds(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GRASHOF:  ',y%wood_grashof(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFREE:   ',y%wood_nussfree(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFORC:   ',y%wood_nussforc(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_EN)/Dt: ',dydx%wood_energy(ico)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_WAT)/Dt:',dydx%wood_water(ico)
                  write(unit=*,fmt='(a)')           '========================================'
               elseif (.not. record_err) then
                  return
               end if
            end if
         end do woodloop
         if(record_err .and. cflag9 ) integ_err( 9,2) = integ_err( 9,2) + 1_8
         if(record_err .and. cflag10) integ_err(10,2) = integ_err(10,2) + 1_8
         !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check the water mass of the virtual pool.  The energy is checked only when     !
      ! there is enough mass.                                                              !
      !------------------------------------------------------------------------------------!
      if (y%virtual_water < rk4min_virt_water) then
         reject_step = .true.
         if(record_err) integ_err(12,2) = integ_err(12,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Virtual layer mass is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_ENERGY:   ',y%virtual_energy
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_WATER:    ',y%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_DEPTH:    ',y%virtual_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_TEMPK:    ',y%virtual_tempk
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_FLIQ :    ',y%virtual_fracliq
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_WATER)/Dt: ',dydx%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_ENERGY)/Dt:',dydx%virtual_energy
            write(unit=*,fmt='(a)')           '==========================================='
         elseif (.not. record_err) then
            return
         end if
      elseif (y%virtual_water > 5.d-1 * rk4water_stab_thresh .and.                         &
           (y%virtual_tempk < rk4min_sfcw_temp .or. y%virtual_tempk > rk4max_sfcw_temp)) &
           then
         reject_step = .true.
         if(record_err) integ_err(11,2) = integ_err(11,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + Virtual layer temp. is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_ENERGY:   ',y%virtual_energy
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_WATER:    ',y%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_DEPTH:    ',y%virtual_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_TEMPK:    ',y%virtual_tempk
            write(unit=*,fmt='(a,1x,es12.4)') ' VIRTUAL_FLIQ :    ',y%virtual_fracliq
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_WATER)/Dt: ',dydx%virtual_water
            write(unit=*,fmt='(a,1x,es12.4)') ' D(VIRT_ENERGY)/Dt:',dydx%virtual_energy
            write(unit=*,fmt='(a)')           '==========================================='
         elseif (.not. record_err) then
            return
         end if
         return
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Checking whether the soil layers have decent moisture and temperatures.         !
      !------------------------------------------------------------------------------------!
      do k=rk4site%lsl,nzg
         !----- Soil moisture -------------------------------------------------------------!
         if (y%soil_water(k)< rk4aux(ibuff)%rk4min_soil_water(k) .or.                                 &
              y%soil_water(k)> rk4aux(ibuff)%rk4max_soil_water(k) ) then
            reject_step = .true.
            if(record_err) integ_err(osow+k,2) = integ_err(osow+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Soil layer water is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,f12.4)') ' H:           ',h
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_MSTPOT: ',y%soil_mstpot(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SOIL_E)/Dt:',dydx%soil_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SOIL_M)/Dt:',dydx%soil_water(k)
               if (k == nzg .and. y%nlev_sfcwater > 0) then
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(1)
               end if
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Soil temperature ----------------------------------------------------------!
         if (y%soil_tempk(k) > rk4max_soil_temp .or. y%soil_tempk(k) < rk4min_soil_temp )  &
              then
            reject_step = .true.
            if(record_err) integ_err(osoe+k,2) = integ_err(osoe+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Soil layer temp is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,f12.4)') ' H:           ',h
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SOIL_MSTPOT: ',y%soil_mstpot(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SOIL_E)/Dt:',dydx%soil_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SOIL_M)/Dt:',dydx%soil_water(k)
               if (k == nzg .and. y%nlev_sfcwater > 0) then
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(1)
               end if
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Check whether the temporary snow/water layer(s) has(ve) reasonable values.      !
      !------------------------------------------------------------------------------------!
      ksn = y%nlev_sfcwater

      do k=1, ksn
         !----- Temperature ---------------------------------------------------------------!
         if (y%sfcwater_tempk(k) < rk4min_sfcw_temp .or.                                   &
              y%sfcwater_tempk(k) > rk4max_sfcw_temp      ) then
            reject_step = .true.
            if(record_err) integ_err(oswe+ksn,2) = integ_err(oswe+ksn,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Snow/pond temperature is off...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' This layer:    ',k
               write(unit=*,fmt='(a,1x,i6)')     ' # of layers:   ',y%nlev_sfcwater
               write(unit=*,fmt='(a,1x,i6)')     ' Stability flag:',y%flag_sfcwater
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_TEMP:     ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_ENERGY:   ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_MASS:     ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_DEPTH:    ',y%sfcwater_depth(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_E)/Dt:  ',dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_M)/Dt:  ',dydx%sfcwater_mass(k)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if

         !----- Mass ----------------------------------------------------------------------!
         if (y%sfcwater_mass(k) < rk4min_sfcw_mass) then
            reject_step = .true.
            if(record_err) integ_err(oswm+ksn,2) = integ_err(oswm+ksn,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Snow/pond mass is off...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' This layer:    ',k
               write(unit=*,fmt='(a,1x,i6)')     ' # of layers:   ',y%nlev_sfcwater
               write(unit=*,fmt='(a,1x,i6)')     ' Stability flag:',y%flag_sfcwater
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_TEMP:     ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_ENERGY:   ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_MASS:     ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' SFCW_DEPTH:    ',y%sfcwater_depth(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_E)/Dt:  ',dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,f12.4)') ' D(SFCW_M)/Dt:  ',dydx%sfcwater_mass(k)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      !------------------------------------------------------------------------------------!

      if (reject_step .and. print_problems) then

         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('=',k=1,78)
         write(unit=*,fmt='(a,1x,f12.4)') ' TIMESTEP:          ',h
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           '         ---- SANITY CHECK BOUNDS ----'
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 1. CANOPY AIR SPACE: '
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(6(a,1x))')     '   MIN_THEIV','   MAX_THEIV','     MIN_SHV'    &
              ,'     MAX_SHV','     MIN_RHV','     MAX_RHV'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(4(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_THETA'    &
              ,'   MAX_THETA'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(4(a,1x))')     '    MIN_PRSS','    MAX_PRSS','     MIN_CO2'    &
              ,'     MAX_CO2'
         write(unit=*,fmt='(4(f12.5,1x))') rk4aux(ibuff)%rk4min_can_prss ,rk4aux(ibuff)%rk4max_can_prss              &
              ,rk4min_can_co2  ,rk4max_can_co2
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 2. LEAF PROPERTIES: '
         write(unit=*,fmt='(3(a,1x))')     '    MIN_TEMP','    MAX_TEMP','  MIN_LWATER'
         write(unit=*,fmt='(3(f12.5,1x))') rk4min_veg_temp ,rk4max_veg_temp               &
              ,rk4min_veg_lwater
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 3. SURFACE WATER / VIRTUAL POOL PROPERTIES: '
         write(unit=*,fmt='(3(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_WMASS'
         write(unit=*,fmt='(3(f12.5,1x))') rk4min_sfcw_temp ,rk4max_sfcw_temp             &
              ,rk4min_sfcw_mass
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 4. SOIL (TEXTURE CLASS AT TOP LAYER): '
         write(unit=*,fmt='(4(a,1x))')     '   MIN_WATER','   MAX_WATER','    MIN_TEMP'    &
              ,'    MAX_TEMP'
         write(unit=*,fmt='(4(f12.5,1x))') rk4aux(ibuff)%rk4min_soil_water(nzg),rk4aux(ibuff)%rk4max_soil_water(nzg)  &
              ,rk4min_soil_temp      ,rk4max_soil_temp
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('=',k=1,78)
         write(unit=*,fmt='(a)')           ' '
      end if

      return
    end subroutine fb_sanity_check


 
end module hybrid_driver

