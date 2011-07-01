!==========================================================================================!
!==========================================================================================!
! Subroutine odeint                                                                        !
!                                                                                          !
!     This subroutine will drive the integration of several ODEs that drive the fast-scale !
! state variables.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine odeint(h1,csite,ipa,nsteps)

   use ed_state_vars  , only : sitetype               & ! structure
                             , patchtype              & ! structure
                             , polygontype
   use rk4_coms       , only : integration_vars       & ! structure
                             , integration_buff       & ! intent(inout)
                             , rk4site                & ! intent(in)
                             , rk4tiny_sfcw_mass      & ! intent(in)
                             , maxstp                 & ! intent(in)
                             , tbeg                   & ! intent(in)
                             , tend                   & ! intent(in)
                             , dtrk4                  & ! intent(in)
                             , dtrk4i                 & ! intent(in)
                             , tiny_offset            & ! intent(in)
                             , checkbudget            & ! intent(in)
                             , print_detailed         & ! intent(in)
                             , norm_rk4_fluxes        ! ! sub-routine
   use rk4_stepper    , only : rkqs                   ! ! subroutine
   use ed_misc_coms   , only : fast_diagnostics       ! ! intent(in)
   use hydrology_coms , only : useRUNOFF              ! ! intent(in)
   use grid_coms      , only : nzg                    & ! intent(in)
                             , nzs                    ! ! intent(in)
   use soil_coms      , only : dslz8                  & ! intent(in)
                             , runoff_time            ! ! intent(in)
   use consts_coms    , only : cliq8                  & ! intent(in)
                             , t3ple8                 & ! intent(in)
                             , tsupercool8            & ! intent(in)
                             , wdnsi8                 ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite            ! Current site
   integer                   , intent(in)  :: ipa              ! Current patch ID
   real(kind=8)              , intent(in)  :: h1               ! First guess of delta-t
   integer                   , intent(out) :: nsteps           ! Number of steps taken.
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch           ! Current patch
   integer                                 :: i                ! Step counter
   integer                                 :: ksn              ! # of snow/water layers
   real(kind=8)                            :: x                ! Elapsed time
   real(kind=8)                            :: h                ! Current delta-t attempt
   real(kind=8)                            :: hnext            ! Next delta-t
   real(kind=8)                            :: hdid             ! delta-t that worked (???)
   real(kind=8)                            :: qwfree           ! Free water internal energy
   real(kind=8)                            :: wfreeb           ! Free water 
   !----- Saved variables -----------------------------------------------------------------!
   logical                   , save        :: first_time=.true.
   logical                   , save        :: simplerunoff
   real(kind=8)              , save        :: runoff_time_i
   !----- External function. --------------------------------------------------------------!
   real                      , external    :: sngloff
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Check whether we will use runoff or not, and saving this check to save time.      !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      simplerunoff = useRUNOFF == 0 .and. runoff_time /= 0.
      if (runoff_time /= 0.) then
         runoff_time_i = 1.d0/dble(runoff_time)
      else 
         runoff_time_i = 0.d0
      end if
      first_time   = .false.
   end if
   !---------------------------------------------------------------------------------------!


   cpatch => csite%patch(ipa)

   !---------------------------------------------------------------------------------------!
   !     Copy the initial patch to the one we use for integration.                         !
   !---------------------------------------------------------------------------------------!
   call copy_rk4_patch(integration_buff%initp, integration_buff%y,cpatch)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! Set initial time and stepsize.                                                        !
   !---------------------------------------------------------------------------------------!
   x = tbeg
   h = h1
   if (dtrk4 < 0.d0) h = -h1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Begin timestep loop                                                                   !
   !---------------------------------------------------------------------------------------!
   timesteploop: do i=1,maxstp

      !----- Get initial derivatives ------------------------------------------------------!
      call leaf_derivs(integration_buff%y,integration_buff%dydx,csite,ipa)

      !----- Get scalings used to determine stability -------------------------------------!
      call get_yscal(integration_buff%y, integration_buff%dydx,h,integration_buff%yscal    &
                       ,cpatch)

      !----- Be sure not to overstep ------------------------------------------------------!
      if((x+h-tend)*(x+h-tbeg) > 0.d0) h=tend-x

      !----- Take the step ----------------------------------------------------------------!
      call rkqs(x,h,hdid,hnext,csite,ipa)

      !----- If the integration reached the next step, make some final adjustments --------!
      if((x-tend)*dtrk4 >= 0.d0)then

         ksn = integration_buff%y%nlev_sfcwater

         !---------------------------------------------------------------------------------!
         !   Make temporary surface liquid water disappear.  This will not happen          !
         ! immediately, but liquid water will decay with the time scale defined by         !
         ! runoff_time scale. If the time scale is too tiny, then it will be forced to be  !
         ! hdid (no reason to be faster than that).                                        !
         !---------------------------------------------------------------------------------!
         if (simplerunoff .and. ksn >= 1) then
            if (integration_buff%y%sfcwater_mass(ksn)    > 0.d0   .and.                    &
                integration_buff%y%sfcwater_fracliq(ksn) > 1.d-1) then

               wfreeb = min(1.d0, dtrk4 * runoff_time_i)                                   &
                      * integration_buff%y%sfcwater_mass(ksn)                              &
                      * (integration_buff%y%sfcwater_fracliq(ksn) - 1.d-1) / 9.d-1

               qwfree = wfreeb                                                             &
                      * cliq8 * (integration_buff%y%sfcwater_tempk(ksn) - tsupercool8 )

               integration_buff%y%sfcwater_mass(ksn) =                                     &
                                   integration_buff%y%sfcwater_mass(ksn)                   &
                                 - wfreeb

               integration_buff%y%sfcwater_depth(ksn) =                                    &
                                   integration_buff%y%sfcwater_depth(ksn)                  &
                                 - wfreeb*wdnsi8

               !----- Recompute the energy removing runoff --------------------------------!
               integration_buff%y%sfcwater_energy(ksn) =                                   &
                                     integration_buff%y%sfcwater_energy(ksn) - qwfree

               call adjust_sfcw_properties(nzg,nzs,integration_buff%y,dtrk4,csite,ipa)
               call update_diagnostic_vars(integration_buff%y,csite,ipa)

               !----- Compute runoff for output -------------------------------------------!
               if (fast_diagnostics) then
                  !------------------------------------------------------------------------!
                  !      There is no need to divide wfreeb and qwfree  by time step, which !
                  ! will be done in subroutine normalize_averaged_vars.                    !
                  !------------------------------------------------------------------------!
                  csite%runoff(ipa)          = csite%runoff(ipa)                           &
                                             + sngloff(wfreeb,tiny_offset)
                  csite%avg_runoff(ipa)      = csite%avg_runoff(ipa)                       &
                                             + sngloff(wfreeb,tiny_offset)
                  csite%avg_runoff_heat(ipa) = csite%avg_runoff_heat(ipa)                  &
                                             + sngloff(qwfree,tiny_offset)
               end if
               if (checkbudget) then
                  !------------------------------------------------------------------------!
                  !      To make sure that the previous values of wbudget_loss2runoff and  !
                  ! ebudget_loss2runoff are accumulated to the next time step.             !
                  !------------------------------------------------------------------------!
                  integration_buff%y%wbudget_loss2runoff = wfreeb                          &
                                    + integration_buff%y%wbudget_loss2runoff
                  integration_buff%y%ebudget_loss2runoff = qwfree                          &
                                    + integration_buff%y%ebudget_loss2runoff
                  integration_buff%y%wbudget_storage =                                     &
                                      integration_buff%y%wbudget_storage - wfreeb
                  integration_buff%y%ebudget_storage =                                     &
                                      integration_buff%y%ebudget_storage - qwfree
               end if
            end if
         end if

         !------ Copy the temporary patch to the next intermediate step -------------------!
         call copy_rk4_patch(integration_buff%y,integration_buff%initp, cpatch)
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
   write (unit=*,fmt='(a)') ' ==> Too many steps in routine odeint'
   call print_rk4patch(integration_buff%y, csite,ipa)

   return
end subroutine odeint
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies the meteorological variables to the Runge-Kutta buffer.  This  !
! is to ensure all variables are in double precision, so consistent with the buffer vari-  !
! ables.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine copy_met_2_rk4site(mzg,vels,atm_theiv,atm_theta,atm_tmp,atm_shv,atm_co2,zoff    &
                             ,exner,pcpg,qpcpg,dpcpg,prss,rshort,rlong,geoht,lsl           &
                             ,ntext_soil,green_leaf_factor,lon,lat,cosz)
   use ed_max_dims    , only : n_pft         ! ! intent(in)
   use rk4_coms       , only : rk4site       ! ! structure
   use canopy_air_coms, only : ubmin8        ! ! intent(in)
   use therm_lib8     , only : rehuil8       & ! function
                             , idealdenssh8  ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in) :: mzg
   integer                  , intent(in) :: lsl
   real                     , intent(in) :: vels
   real                     , intent(in) :: atm_theiv
   real                     , intent(in) :: atm_theta
   real                     , intent(in) :: atm_tmp
   real                     , intent(in) :: atm_shv
   real                     , intent(in) :: atm_co2
   real                     , intent(in) :: zoff
   real                     , intent(in) :: exner
   real                     , intent(in) :: pcpg
   real                     , intent(in) :: qpcpg
   real                     , intent(in) :: dpcpg
   real                     , intent(in) :: prss
   real                     , intent(in) :: rshort
   real                     , intent(in) :: rlong
   real                     , intent(in) :: geoht
   integer, dimension(mzg)  , intent(in) :: ntext_soil
   real   , dimension(n_pft), intent(in) :: green_leaf_factor
   real                     , intent(in) :: lon
   real                     , intent(in) :: lat
   real                     , intent(in) :: cosz
   !----- Local variables. ----------------------------------------------------------------!
   integer             :: ipft
   !---------------------------------------------------------------------------------------!

   
   !----- Copy the integer variables. -----------------------------------------------------!
   rk4site%lsl               = lsl
   rk4site%ntext_soil(:)     = 0
   rk4site%ntext_soil(1:mzg) = ntext_soil(1:mzg)

   !----- Convert to double precision. ----------------------------------------------------!
   rk4site%atm_theiv             = dble(atm_theiv           )
   rk4site%atm_theta             = dble(atm_theta           )
   rk4site%atm_tmp               = dble(atm_tmp             )
   rk4site%atm_shv               = dble(atm_shv             )
   rk4site%atm_co2               = dble(atm_co2             )
   rk4site%zoff                  = dble(zoff                )
   rk4site%atm_exner             = dble(exner               )
   rk4site%pcpg                  = dble(pcpg                )
   rk4site%qpcpg                 = dble(qpcpg               )
   rk4site%dpcpg                 = dble(dpcpg               )
   rk4site%atm_prss              = dble(prss                )
   rk4site%rshort                = dble(rshort              )
   rk4site%rlong                 = dble(rlong               )
   rk4site%geoht                 = dble(geoht               )
   rk4site%lon                   = dble(lon                 )
   rk4site%lat                   = dble(lat                 )
   rk4site%cosz                  = dble(cosz                )
   rk4site%green_leaf_factor(:)  = dble(green_leaf_factor(:))
   !---------------------------------------------------------------------------------------!


   !----- Find the other variables that require a little math. ----------------------------!
   rk4site%vels                  = max(ubmin8,dble(vels))
   rk4site%atm_lntheta           = log(rk4site%atm_theta)
   rk4site%atm_rvap              = rk4site%atm_shv / (1.d0 - rk4site%atm_shv)
   rk4site%atm_rhv               = rehuil8(rk4site%atm_prss,rk4site%atm_tmp                &
                                          ,rk4site%atm_rvap)
   rk4site%atm_rhos              = idealdenssh8(rk4site%atm_prss,rk4site%atm_tmp           &
                                               ,rk4site%atm_shv)
   !---------------------------------------------------------------------------------------!


   return
end subroutine copy_met_2_rk4site
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutines increment the derivative into the previous guess to create the new   !
! guess.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine inc_rk4_patch(rkp, inc, fac, cpatch)
   use ed_state_vars , only : sitetype           & ! structure
                            , patchtype          ! ! structure
   use rk4_coms      , only : rk4patchtype       & ! structure
                            , rk4site            & ! intent(in)
                            , checkbudget        & ! intent(in)
                            , print_detailed     ! ! intent(in)
   use grid_coms     , only : nzg                & ! intent(in)
                            , nzs                ! ! intent(in)
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

   rkp%can_lntheta  = rkp%can_lntheta  + fac * inc%can_lntheta
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

  
   do ico = 1,cpatch%ncohorts
      rkp%leaf_water (ico) = rkp%leaf_water (ico) + fac * inc%leaf_water (ico)
      rkp%leaf_energy(ico) = rkp%leaf_energy(ico) + fac * inc%leaf_energy(ico)
      rkp%wood_water (ico) = rkp%wood_water (ico) + fac * inc%wood_water (ico)
      rkp%wood_energy(ico) = rkp%wood_energy(ico) + fac * inc%wood_energy(ico)

      rkp%psi_open  (ico) = rkp%psi_open  (ico) + fac * inc%psi_open  (ico)
      rkp%psi_closed(ico) = rkp%psi_closed(ico) + fac * inc%psi_closed(ico)
   end do

   if(checkbudget) then

     rkp%co2budget_storage      = rkp%co2budget_storage     + fac * inc%co2budget_storage
     rkp%co2budget_loss2atm     = rkp%co2budget_loss2atm    + fac * inc%co2budget_loss2atm

      rkp%wbudget_storage       = rkp%wbudget_storage       + fac * inc%wbudget_storage
      rkp%wbudget_loss2atm      = rkp%wbudget_loss2atm      + fac * inc%wbudget_loss2atm
      rkp%wbudget_loss2drainage = rkp%wbudget_loss2drainage                                &
                                + fac * inc%wbudget_loss2drainage

      rkp%ebudget_storage       = rkp%ebudget_storage       + fac * inc%ebudget_storage
      rkp%ebudget_loss2atm      = rkp%ebudget_loss2atm      + fac * inc%ebudget_loss2atm
      rkp%ebudget_loss2drainage = rkp%ebudget_loss2drainage                                &
                                + fac * inc%ebudget_loss2drainage
   end if
   if (fast_diagnostics) then
      rkp%avg_carbon_ac      = rkp%avg_carbon_ac      + fac * inc%avg_carbon_ac
      
      rkp%avg_vapor_lc       = rkp%avg_vapor_lc       + fac * inc%avg_vapor_lc
      rkp%avg_vapor_wc       = rkp%avg_vapor_wc       + fac * inc%avg_vapor_wc
      rkp%avg_dew_cg         = rkp%avg_dew_cg         + fac * inc%avg_dew_cg
      rkp%avg_vapor_gc       = rkp%avg_vapor_gc       + fac * inc%avg_vapor_gc
      rkp%avg_wshed_vg       = rkp%avg_wshed_vg       + fac * inc%avg_wshed_vg
      rkp%avg_intercepted    = rkp%avg_intercepted    + fac * inc%avg_intercepted
      rkp%avg_throughfall    = rkp%avg_throughfall    + fac * inc%avg_throughfall
      rkp%avg_vapor_ac       = rkp%avg_vapor_ac       + fac * inc%avg_vapor_ac
      rkp%avg_transp         = rkp%avg_transp         + fac * inc%avg_transp
      rkp%avg_evap           = rkp%avg_evap           + fac * inc%avg_evap
      rkp%avg_drainage       = rkp%avg_drainage       + fac * inc%avg_drainage
      rkp%avg_drainage_heat  = rkp%avg_drainage_heat  + fac * inc%avg_drainage_heat
      rkp%avg_rshort_gnd     = rkp%avg_rshort_gnd     + fac * inc%avg_rshort_gnd
      rkp%avg_rlong_gnd      = rkp%avg_rlong_gnd      + fac * inc%avg_rlong_gnd
      rkp%avg_sensible_lc    = rkp%avg_sensible_lc    + fac * inc%avg_sensible_lc
      rkp%avg_sensible_wc    = rkp%avg_sensible_wc    + fac * inc%avg_sensible_wc
      rkp%avg_qwshed_vg      = rkp%avg_qwshed_vg      + fac * inc%avg_qwshed_vg
      rkp%avg_qintercepted   = rkp%avg_qintercepted   + fac * inc%avg_qintercepted
      rkp%avg_qthroughfall   = rkp%avg_qthroughfall   + fac * inc%avg_qthroughfall
      rkp%avg_sensible_gc    = rkp%avg_sensible_gc    + fac * inc%avg_sensible_gc
      rkp%avg_sensible_ac    = rkp%avg_sensible_ac    + fac * inc%avg_sensible_ac

      do k=rk4site%lsl,nzg
         rkp%avg_sensible_gg(k)  = rkp%avg_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
         rkp%avg_smoist_gg(k)    = rkp%avg_smoist_gg(k)    + fac * inc%avg_smoist_gg(k)  
         rkp%avg_transloss(k)    = rkp%avg_transloss(k)    + fac * inc%avg_transloss(k)  
      end do

   end if

   !---------------------------------------------------------------------------------------!
   !    Increment the instantaneous fluxes.  The derivative term should be the same as the !
   ! the full fluxes, the only difference is that these variables are normalised and       !
   ! re-set after each time step.                                                          !
   !---------------------------------------------------------------------------------------!
   if (print_detailed) then
      rkp%flx_carbon_ac      = rkp%flx_carbon_ac      + fac * inc%avg_carbon_ac
      
      rkp%flx_vapor_lc       = rkp%flx_vapor_lc       + fac * inc%avg_vapor_lc
      rkp%flx_vapor_wc       = rkp%flx_vapor_wc       + fac * inc%avg_vapor_wc
      rkp%flx_dew_cg         = rkp%flx_dew_cg         + fac * inc%avg_dew_cg
      rkp%flx_vapor_gc       = rkp%flx_vapor_gc       + fac * inc%avg_vapor_gc
      rkp%flx_wshed_vg       = rkp%flx_wshed_vg       + fac * inc%avg_wshed_vg
      rkp%flx_intercepted    = rkp%flx_intercepted    + fac * inc%avg_intercepted
      rkp%flx_throughfall    = rkp%flx_throughfall    + fac * inc%avg_throughfall
      rkp%flx_vapor_ac       = rkp%flx_vapor_ac       + fac * inc%avg_vapor_ac
      rkp%flx_transp         = rkp%flx_transp         + fac * inc%avg_transp
      rkp%flx_evap           = rkp%flx_evap           + fac * inc%avg_evap
      rkp%flx_drainage       = rkp%flx_drainage       + fac * inc%avg_drainage
      rkp%flx_drainage_heat  = rkp%flx_drainage_heat  + fac * inc%avg_drainage_heat
      rkp%flx_rshort_gnd     = rkp%flx_rshort_gnd     + fac * inc%avg_rshort_gnd
      rkp%flx_rlong_gnd      = rkp%flx_rlong_gnd      + fac * inc%avg_rlong_gnd
      rkp%flx_sensible_lc    = rkp%flx_sensible_lc    + fac * inc%avg_sensible_lc
      rkp%flx_sensible_wc    = rkp%flx_sensible_wc    + fac * inc%avg_sensible_wc
      rkp%flx_qwshed_vg      = rkp%flx_qwshed_vg      + fac * inc%avg_qwshed_vg
      rkp%flx_qintercepted   = rkp%flx_qintercepted   + fac * inc%avg_qintercepted
      rkp%flx_qthroughfall   = rkp%flx_qthroughfall   + fac * inc%avg_qthroughfall
      rkp%flx_sensible_gc    = rkp%flx_sensible_gc    + fac * inc%avg_sensible_gc
      rkp%flx_sensible_ac    = rkp%flx_sensible_ac    + fac * inc%avg_sensible_ac

      do k=rk4site%lsl,nzg
         rkp%flx_sensible_gg(k)  = rkp%flx_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
         rkp%flx_smoist_gg(k)    = rkp%flx_smoist_gg(k)    + fac * inc%avg_smoist_gg(k)  
         rkp%flx_transloss(k)    = rkp%flx_transloss(k)    + fac * inc%avg_transloss(k)  
      end do

      do ico = 1,cpatch%ncohorts
         rkp%cfx_hflxlc      (ico)  =       rkp%cfx_hflxlc      (ico)                      &
                                    + fac * inc%cfx_hflxlc      (ico)
         rkp%cfx_hflxwc      (ico)  =       rkp%cfx_hflxwc      (ico)                      &
                                    + fac * inc%cfx_hflxwc      (ico)
         rkp%cfx_qwflxlc     (ico)  =       rkp%cfx_qwflxlc     (ico)                      &
                                    + fac * inc%cfx_qwflxlc     (ico)
         rkp%cfx_qwflxwc     (ico)  =       rkp%cfx_qwflxwc     (ico)                      &
                                    + fac * inc%cfx_qwflxwc     (ico)
         rkp%cfx_qwshed      (ico)  =       rkp%cfx_qwshed      (ico)                      &
                                    + fac * inc%cfx_qwshed      (ico)
         rkp%cfx_qtransp     (ico)  =       rkp%cfx_qtransp     (ico)                      &
                                    + fac * inc%cfx_qtransp     (ico)
         rkp%cfx_qintercepted(ico)  =       rkp%cfx_qintercepted(ico)                      &
                                    + fac * inc%cfx_qintercepted(ico)
      end do

   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine inc_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the error scale for the integrated variables, which will be     !
! later used to define the relative error.                                                 !
!------------------------------------------------------------------------------------------!
subroutine get_yscal(y,dy,htry,yscal,cpatch)
   use ed_state_vars        , only : patchtype             ! ! structure
   use rk4_coms             , only : rk4patchtype          & ! structure
                                   , rk4site               & ! intent(in)
                                   , tiny_offset           & ! intent(in)
                                   , huge_offset           & ! intent(in)
                                   , rk4water_stab_thresh  & ! intent(in)
                                   , rk4tiny_sfcw_mass     & ! intent(in)
                                   , rk4leaf_drywhc        & ! intent(in)
                                   , checkbudget           ! ! intent(in)
   use grid_coms            , only : nzg                   & ! intent(in)
                                   , nzs                   ! ! intent(in)
   use consts_coms          , only : cliq8                 & ! intent(in)
                                   , qliqt38               & ! intent(in)
                                   , wdnsi8                ! ! intent(in)
   use soil_coms            , only : isoilbc               & ! intent(in)
                                   , dslzi8                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype), target     :: y                     ! Struct. with the guesses
   type(rk4patchtype), target     :: dy                    ! Struct. with their derivatives
   type(rk4patchtype), target     :: yscal                 ! Struct. with their scales
   type(patchtype)   , target     :: cpatch                ! Current patch
   real(kind=8)      , intent(in) :: htry                  ! Time-step we are trying
   !----- Local variables -----------------------------------------------------------------!
   real(kind=8)                   :: meanscale_sfcw_mass   ! Average Sfc. water mass scale
   real(kind=8)                   :: meanscale_sfcw_energy ! Average Sfc. water en. scale
   real(kind=8)                   :: meanscale_sfcw_depth  ! Average Sfc. water depth scale
   integer                        :: k                     ! Counter
   integer                        :: ico                   ! Current cohort ID
   !---------------------------------------------------------------------------------------!

   yscal%can_lntheta =  abs(y%can_lntheta) + abs(dy%can_lntheta  * htry)
   yscal%can_shv     =  abs(y%can_shv    ) + abs(dy%can_shv      * htry)
   yscal%can_co2     =  abs(y%can_co2    ) + abs(dy%can_co2      * htry)

   !---------------------------------------------------------------------------------------!
   !     We don't solve pressure prognostically, so the scale cannot be computed based on  !
   ! the derivative.  Also, pressure is a variable with a large absolute variable and tiny !
   ! variation, so the scale is given by the scale of the variability...                   !
   !---------------------------------------------------------------------------------------!
   yscal%can_prss    = abs(rk4site%atm_prss - y%can_prss)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     We determine the scale for all layers but the top soil, which will be done        !
   ! differently depending on the status of the temporary surface water.                   !
   !---------------------------------------------------------------------------------------!
   do k=1,rk4site%lsl-1
      yscal%soil_water (k) = huge_offset
      yscal%soil_energy(k) = huge_offset
   end do
   do k=rk4site%lsl,nzg-1
      yscal%soil_water (k) = abs(y%soil_water(k) ) + abs(dy%soil_water(k)  * htry)
      yscal%soil_energy(k) = abs(y%soil_energy(k)) + abs(dy%soil_energy(k) * htry)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Temporary surface layers require a special approach. The number of layers may     !
   ! vary during the integration process, so we must make sure that all layers initially   !
   ! have a scale.  Also, if the total mass is small, we must be more tolerant to avoid    !
   ! overestimating the error because of the small size.                                   !
   !---------------------------------------------------------------------------------------!
   select case(y%flag_sfcwater)
   case(0)
      !------------------------------------------------------------------------------------!
      !     No temporary surface water, we can skip the check as no layer will be created  !
      ! until the upcoming step.                                                           !
      !------------------------------------------------------------------------------------!
      do k=1,nzs
         yscal%sfcwater_mass(k)   = huge_offset
         yscal%sfcwater_energy(k) = huge_offset
         yscal%sfcwater_depth(k)  = huge_offset
      end do
      !------------------------------------------------------------------------------------!


      !------ Soil scale won't be affected by the temporary surface water. ----------------!
      yscal%soil_water (nzg) = abs(y%soil_water (nzg)) + abs(dy%soil_water (nzg) * htry)
      yscal%soil_energy(nzg) = abs(y%soil_energy(nzg)) + abs(dy%soil_energy(nzg) * htry)
      !------------------------------------------------------------------------------------!

   case(1)
      !------------------------------------------------------------------------------------!
      !    Low stability threshold, there can't be more than one layer, and the energy     !
      ! will be solved together with the top soil layer.  Therefore, we skip energy check, !
      ! but attribute the scale for mass and depth.                                        !
      !------------------------------------------------------------------------------------!
      yscal%sfcwater_mass  (1) = abs(y%sfcwater_mass   (1)       )                         &
                               + abs(dy%sfcwater_mass  (1) * htry)
      yscal%sfcwater_depth (1) = abs(y%sfcwater_depth  (1)       )                         &
                               + abs(dy%sfcwater_energy(1) * htry)
      yscal%sfcwater_energy(1) = huge_offset
      do k=2,nzs
         yscal%sfcwater_mass  (k) = huge_offset
         yscal%sfcwater_energy(k) = huge_offset
         yscal%sfcwater_depth (k) = huge_offset
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Soil mass won't be affected by the temporary surface water, but soil energy    !
      ! will include the energy associated with the temporary surface water energy.  In    !
      ! reality the derivative of the temporary surface water will be zero, but we add for !
      ! code consistency.                                                                  !
      !------------------------------------------------------------------------------------!
      yscal%soil_water (nzg) = abs(y%soil_water (nzg)) + abs(dy%soil_water (nzg) * htry)
      yscal%soil_energy(nzg) = abs(y%soil_energy(nzg)) + abs(dy%soil_energy(nzg) * htry)   &
                             + dslzi8(nzg) * ( abs(y%sfcwater_energy(1))                   &
                                             + abs(dy%sfcwater_energy(1) * htry) )
      !------------------------------------------------------------------------------------!

   case (2)
      !----- Computationally stable layer. ------------------------------------------------!
      meanscale_sfcw_mass   = 0.d0
      meanscale_sfcw_energy = 0.d0
      meanscale_sfcw_depth  = 0.d0
      do k=1,y%nlev_sfcwater
         yscal%sfcwater_mass  (k) = abs(y%sfcwater_mass   (k)       )                      &
                                  + abs(dy%sfcwater_mass  (k) * htry)
         yscal%sfcwater_energy(k) = abs(y%sfcwater_energy (k)       )                      &
                                  + abs(dy%sfcwater_energy(k) * htry)
         yscal%sfcwater_depth (k) = abs(y%sfcwater_depth  (k)       )                      &
                                  + abs(dy%sfcwater_energy(k) * htry)
      end do
      do k=y%nlev_sfcwater+1,nzs
         yscal%sfcwater_mass  (k) = huge_offset
         yscal%sfcwater_energy(k) = huge_offset
         yscal%sfcwater_depth (k) = huge_offset
      end do


      !------ Soil scale won't be affected by the temporary surface water. ----------------!
      yscal%soil_water (nzg) = abs(y%soil_water (nzg)) + abs(dy%soil_water (nzg) * htry)
      yscal%soil_energy(nzg) = abs(y%soil_energy(nzg)) + abs(dy%soil_energy(nzg) * htry)
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !----- Scale for the virtual water pools -----------------------------------------------!
   if (abs(y%virtual_water) > 1.d-2*rk4water_stab_thresh) then
      yscal%virtual_water   = abs(y%virtual_water)  + abs(dy%virtual_water*htry)
      yscal%virtual_energy  = abs(y%virtual_energy) + abs(dy%virtual_energy*htry)
   elseif (abs(y%virtual_water) > rk4tiny_sfcw_mass) then
      yscal%virtual_water   = 1.d-2*rk4water_stab_thresh
      yscal%virtual_energy  = (yscal%virtual_water / abs(y%virtual_water))                 &
                            * (abs(y%virtual_energy) + abs(dy%virtual_energy*htry))
   else
      yscal%virtual_water   = huge_offset
      yscal%virtual_energy  = huge_offset
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Scale for leaf water and energy. In case the plants have few or no leaves, or the  !
   ! plant is buried in snow, we assign huge values for typical scale, thus preventing     !
   ! unecessary small steps.                                                               !
   !    Also, if the cohort has almost no water, make the scale less strict.               !
   !---------------------------------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      yscal%leaf_resolvable(ico) = y%leaf_resolvable(ico)
      if (y%leaf_resolvable(ico)) then
         yscal%leaf_energy(ico) = abs( y%leaf_energy(ico))                                 &
                                     + abs(dy%leaf_energy(ico) * htry)
         yscal%leaf_temp(ico)   = abs( y%leaf_temp(ico))
         yscal%leaf_water(ico)  = max( abs(y%leaf_water(ico))                              &
                                     + abs(dy%leaf_water(ico)  * htry)                     &
                                    , rk4leaf_drywhc * y%lai(ico))
      else
         yscal%leaf_water(ico)  = huge_offset
         yscal%leaf_energy(ico) = huge_offset
         yscal%leaf_temp(ico)   = huge_offset
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Scale for wood water and energy. In case the user doesn't want to solve branch     !
   ! thermodynamics, the wood area index is too small, or the plant is buried in snow, we  !
   ! assign huge values for typical scale, thus preventing unecessary small steps.         !
   !    Also, if the cohort has almost no water, make the scale less strict.               !
   !---------------------------------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      yscal%wood_resolvable(ico) = y%wood_resolvable(ico)
      if (y%wood_resolvable(ico)) then
         yscal%wood_energy(ico) = abs( y%wood_energy(ico))                                 &
                                     + abs(dy%wood_energy(ico) * htry)
         yscal%wood_temp(ico)   = abs( y%wood_temp(ico))
         yscal%wood_water(ico)  = max( abs(y%wood_water(ico))                              &
                                     + abs(dy%wood_water(ico)  * htry)                     &
                                    , rk4leaf_drywhc * y%wai(ico))
      else
         yscal%wood_water(ico)  = huge_offset
         yscal%wood_energy(ico) = huge_offset
         yscal%wood_temp(ico)   = huge_offset
      end if
   end do
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Here we just need to make sure the user is checking mass, otherwise  these vari-  !
   ! ables will not be computed at all.  If this turns out to be essential, we will make   !
   ! this permanent and not dependent on checkbudget.  The only one that is not checked is !
   ! the runoff, because it is computed after a step was accepted.                         !
   !---------------------------------------------------------------------------------------!
   if (checkbudget) then
      !------------------------------------------------------------------------------------!
      !    If this is the very first time step, or if we are misfortuned, we may have a    !
      ! situation in which the derivative is numerically zero, and making the check will   !
      ! become impossible because the scale would be ridiculously tiny, so we skip the     !
      ! check this time and hope that everything will be alright next step.                !
      !------------------------------------------------------------------------------------!
      if (abs(y%co2budget_loss2atm)  < tiny_offset .and.                                   &
          abs(dy%co2budget_loss2atm) < tiny_offset) then
         yscal%co2budget_loss2atm = 1.d-1
      else 
         yscal%co2budget_loss2atm = abs(y%co2budget_loss2atm)                              &
                                  + abs(dy%co2budget_loss2atm*htry)
         yscal%co2budget_loss2atm = max(yscal%co2budget_loss2atm,1.d-1)
      end if

      if (abs(y%ebudget_loss2atm)  < tiny_offset .and.                                     &
          abs(dy%ebudget_loss2atm) < tiny_offset) then
         yscal%ebudget_loss2atm = 1.d0
      else 
         yscal%ebudget_loss2atm = abs(y%ebudget_loss2atm)                                  &
                                + abs(dy%ebudget_loss2atm*htry)
         yscal%ebudget_loss2atm = max(yscal%ebudget_loss2atm,1.d0)
      end if

      if (abs(y%wbudget_loss2atm)  < tiny_offset .and.                                     &
          abs(dy%wbudget_loss2atm) < tiny_offset) then
         yscal%wbudget_loss2atm      = 1.d-6
      else 
         yscal%wbudget_loss2atm = abs(y%wbudget_loss2atm)                                  &
                                + abs(dy%wbudget_loss2atm*htry)
         yscal%wbudget_loss2atm = max(yscal%wbudget_loss2atm,1.d-6)
      end if

      if (abs(y%ebudget_storage)  < tiny_offset .and.                                      &
          abs(dy%ebudget_storage) < tiny_offset) then
         yscal%ebudget_storage = huge_offset
      else 
         yscal%ebudget_storage = abs(y%ebudget_storage)                                    &
                               + abs(dy%ebudget_storage*htry)
      end if

      if (abs(y%co2budget_storage)  < tiny_offset .and.                                    &
          abs(dy%co2budget_storage) < tiny_offset) then
         yscal%co2budget_storage = huge_offset
      else 
         yscal%co2budget_storage = abs(y%co2budget_storage)                                &
                                 + abs(dy%co2budget_storage*htry)
      end if

      if (abs(y%wbudget_storage)  < tiny_offset .and.                                      &
          abs(dy%wbudget_storage) < tiny_offset) then
         yscal%wbudget_storage      = huge_offset
      else 
         yscal%wbudget_storage = abs(y%wbudget_storage)                                    &
                               + abs(dy%wbudget_storage*htry)
      end if

      !------------------------------------------------------------------------------------!
      !     Drainage terms will be checked only if the boundary condition is free drain-   !
      ! age.                                                                               !
      !------------------------------------------------------------------------------------!
      if (isoilbc == 0 .or. (abs(y%ebudget_loss2drainage)  < tiny_offset .and.             &
                             abs(dy%ebudget_loss2drainage) < tiny_offset)      ) then
         yscal%ebudget_loss2drainage = huge_offset
      else 
         yscal%ebudget_loss2drainage = abs(y%ebudget_loss2drainage)                        &
                                     + abs(dy%ebudget_loss2drainage*htry)
      end if
      if (isoilbc == 0 .or. (abs(y%wbudget_loss2drainage)  < tiny_offset .and.             &
                             abs(dy%wbudget_loss2drainage) < tiny_offset)      ) then
         yscal%wbudget_loss2drainage = huge_offset
      else 
         yscal%wbudget_loss2drainage = abs(y%wbudget_loss2drainage)                        &
                                     + abs(dy%wbudget_loss2drainage*htry)
      end if

   else 
      yscal%co2budget_storage       = huge_offset
      yscal%co2budget_loss2atm      = huge_offset
      yscal%ebudget_loss2atm        = huge_offset
      yscal%wbudget_loss2atm        = huge_offset
      yscal%ebudget_storage         = huge_offset
      yscal%wbudget_storage         = huge_offset
      yscal%ebudget_loss2drainage   = huge_offset
      yscal%wbudget_loss2drainage   = huge_offset
   end if

   return
end subroutine get_yscal
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine loops through the integrating variables, seeking for the largest      !
! error.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine get_errmax(errmax,yerr,yscal,cpatch,y,ytemp)

   use rk4_coms              , only : rk4patchtype       & ! structure
                                    , rk4eps             & ! intent(in)
                                    , rk4site            & ! intent(in)
                                    , checkbudget        & ! intent(in)
                                    , integ_err          & ! intent(inout)
                                    , record_err         & ! intent(in)
                                    , osow               & ! intent(in)
                                    , osoe               & ! intent(in)
                                    , oswe               & ! intent(in)
                                    , oswm               ! ! intent(in)
   use ed_state_vars         , only : patchtype          ! ! structure
   use grid_coms             , only : nzg                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target      :: yerr             ! Error structure
   type(rk4patchtype) , target      :: yscal            ! Scale structure
   type(rk4patchtype) , target      :: y                ! Structure with previous value
   type(rk4patchtype) , target      :: ytemp            ! Structure with attempted values
   type(patchtype)    , target      :: cpatch           ! Current patch
   real(kind=8)       , intent(out) :: errmax           ! Maximum error
   !----- Local variables -----------------------------------------------------------------!
   integer                          :: ico              ! Current cohort ID
   real(kind=8)                     :: errh2o           ! Scratch error variable
   real(kind=8)                     :: errene           ! Scratch error variable
   real(kind=8)                     :: err              ! Scratch error variable
   real(kind=8)                     :: errh2oMAX        ! Scratch error variable
   real(kind=8)                     :: erreneMAX        ! Scratch error variable
   real(kind=8)                     :: scal_err_prss    ! Scaling factor for CAS pressure
   integer                          :: k                ! Counter
   !---------------------------------------------------------------------------------------!

   !----- Initialise the error with an optimistic number... -------------------------------!
   errmax = 0.d0



   !---------------------------------------------------------------------------------------!
   !    We now check each variable error, and keep track of the worst guess, which will    !
   ! be our worst guess in the end.                                                        !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Get the canopy air space errors.  Only the ice-vapour equivalent potential        !
   ! temperature, water vapour mixing ratio and carbon dioxide mixing ratio are accounted. !
   ! Temperature and density will be also checked for sanity.                              !
   !---------------------------------------------------------------------------------------!
   err    = abs(yerr%can_lntheta/yscal%can_lntheta)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(1,1) = integ_err(1,1) + 1_8 

   err    = abs(yerr%can_shv/yscal%can_shv)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(2,1) = integ_err(2,1) + 1_8 

   err    = abs(yerr%can_co2/yscal%can_co2)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(6,1) = integ_err(6,1) + 1_8 

   !---------------------------------------------------------------------------------------!
   !      Pressure is only a semi-prognostic variable, which means that we don't really    !
   ! have a good error estimate, but we try to get a sense of the error.                   !
   !---------------------------------------------------------------------------------------!
   err    = abs(yerr%can_prss/yscal%can_prss)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(5,1) = integ_err(5,1) + 1_8 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Get the worst error only amongst the cohorts in which leaf properties were        !
   ! computed.                                                                             !
   !---------------------------------------------------------------------------------------!
   errh2oMAX  = 0.d0
   erreneMAX  = 0.d0
   do ico = 1,cpatch%ncohorts
      if (yscal%leaf_resolvable(ico)) then
         errh2o     = abs(yerr%leaf_water (ico) / yscal%leaf_water (ico))
         errene     = abs(yerr%leaf_energy(ico) / yscal%leaf_energy(ico))
         errmax     = max(errmax,errh2o,errene)
         errh2oMAX  = max(errh2oMAX ,errh2o )
         erreneMAX  = max(erreneMAX ,errene )
      end if
   end do
   if(cpatch%ncohorts > 0 .and. record_err) then
      if (errh2oMAX  > rk4eps) integ_err(7,1) = integ_err(7,1) + 1_8
      if (erreneMAX  > rk4eps) integ_err(8,1) = integ_err(8,1) + 1_8
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Get the worst error only amongst the cohorts in which wood properties were        !
   ! computed.                                                                             !
   !---------------------------------------------------------------------------------------!
   errh2oMAX  = 0.d0
   erreneMAX  = 0.d0
   do ico = 1,cpatch%ncohorts
      if (yscal%wood_resolvable(ico)) then
         errh2o     = abs(yerr%wood_water (ico) / yscal%wood_water (ico))
         errene     = abs(yerr%wood_energy(ico) / yscal%wood_energy(ico))
         errmax     = max(errmax,errh2o,errene)
         errh2oMAX  = max(errh2oMAX ,errh2o )
         erreneMAX  = max(erreneMAX ,errene )
      end if
   end do
   if(cpatch%ncohorts > 0 .and. record_err) then
      if (errh2oMAX  > rk4eps) integ_err( 9,1) = integ_err( 9,1) + 1_8
      if (erreneMAX  > rk4eps) integ_err(10,1) = integ_err(10,1) + 1_8
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Virtual pool.                                                                     !
   !---------------------------------------------------------------------------------------!
   err    = abs(yerr%virtual_energy/yscal%virtual_energy)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(11,1) = integ_err(11,1) + 1_8

   err    = abs(yerr%virtual_water/yscal%virtual_water)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(12,1) = integ_err(12,1) + 1_8      
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Here we just need to make sure the user is checking mass, otherwise these         !
   ! variables will not be computed at all.  The only one that is not checked is the       !
   ! runoff, because it is computed only after a step was accepted.                        !
   !---------------------------------------------------------------------------------------!
   if (checkbudget) then
      err    = abs(yerr%co2budget_storage/yscal%co2budget_storage)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(13,1) = integ_err(13,1) + 1_8

      err    = abs(yerr%co2budget_loss2atm/yscal%co2budget_loss2atm)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(14,1) = integ_err(14,1) + 1_8

      err    = abs(yerr%ebudget_loss2atm/yscal%ebudget_loss2atm)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(15,1) = integ_err(15,1) + 1_8

      err    = abs(yerr%wbudget_loss2atm/yscal%wbudget_loss2atm)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(16,1) = integ_err(16,1) + 1_8

      err    = abs(yerr%ebudget_loss2drainage/yscal%ebudget_loss2drainage)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(17,1) = integ_err(17,1) + 1_8

      err    = abs(yerr%wbudget_loss2drainage/yscal%wbudget_loss2drainage)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(18,1) = integ_err(18,1) + 1_8

      err    = abs(yerr%ebudget_storage/yscal%ebudget_storage)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(19,1) = integ_err(19,1) + 1_8

      err    = abs(yerr%wbudget_storage/yscal%wbudget_storage)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(20,1) = integ_err(20,1) + 1_8
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Soil properties.  Notice that the index depends on the number of layers.          !
   !---------------------------------------------------------------------------------------!
   do k=rk4site%lsl,nzg
      err    = abs(yerr%soil_water(k)/yscal%soil_water(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(osow+k,1) = integ_err(osow+k,1) + 1_8 
   end do

   do k=rk4site%lsl,nzg
      err    = abs(yerr%soil_energy(k)/yscal%soil_energy(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(osoe+k,1) = integ_err(osoe+k,1) + 1_8
   end do
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Surface water/snow properties.  Notice that the index also depends on the number  !
   ! of layers.                                                                            !
   !---------------------------------------------------------------------------------------!
   do k=1,ytemp%nlev_sfcwater
      err = abs(yerr%sfcwater_energy(k) / yscal%sfcwater_energy(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(oswe+k,1) = integ_err(oswe+k,1) + 1_8
   end do

   do k=1,ytemp%nlev_sfcwater
      err    = abs(yerr%sfcwater_mass(k) / yscal%sfcwater_mass(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(oswm+k,1) = integ_err(oswm+k,1) + 1_8
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine get_errmax
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies the values to different buffers inside the RK4 integration     !
! scheme.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine copy_rk4_patch(sourcep, targetp, cpatch)

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
   !----- Local variable ------------------------------------------------------------------!
   integer                         :: k
   !---------------------------------------------------------------------------------------!

   targetp%can_theiv        = sourcep%can_theiv
   targetp%can_lntheta      = sourcep%can_lntheta
   targetp%can_theta        = sourcep%can_theta
   targetp%can_temp         = sourcep%can_temp
   targetp%can_shv          = sourcep%can_shv
   targetp%can_co2          = sourcep%can_co2
   targetp%can_rhos         = sourcep%can_rhos
   targetp%can_prss         = sourcep%can_prss
   targetp%can_exner        = sourcep%can_exner
   targetp%can_depth        = sourcep%can_depth
   targetp%can_rvap         = sourcep%can_rvap
   targetp%can_rhv          = sourcep%can_rhv
   targetp%can_ssh          = sourcep%can_ssh
   targetp%veg_height       = sourcep%veg_height
   targetp%veg_displace     = sourcep%veg_displace
   targetp%veg_rough        = sourcep%veg_rough
   targetp%opencan_frac     = sourcep%opencan_frac
   targetp%total_sfcw_depth = sourcep%total_sfcw_depth
   targetp%snowfac          = sourcep%snowfac

   targetp%ggbare           = sourcep%ggbare
   targetp%ggveg            = sourcep%ggveg
   targetp%ggnet            = sourcep%ggnet

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

   targetp%cwd_rh           = sourcep%cwd_rh
   targetp%rh               = sourcep%rh

   do k=rk4site%lsl,nzg      
      targetp%soil_water            (k) = sourcep%soil_water            (k)
      targetp%soil_energy           (k) = sourcep%soil_energy           (k)
      targetp%soil_tempk            (k) = sourcep%soil_tempk            (k)
      targetp%soil_fracliq          (k) = sourcep%soil_fracliq          (k)
   end do

   targetp%nlev_sfcwater   = sourcep%nlev_sfcwater
   targetp%flag_sfcwater   = sourcep%flag_sfcwater

   do k=1,nzs
      targetp%sfcwater_mass   (k) = sourcep%sfcwater_mass   (k)
      targetp%sfcwater_energy (k) = sourcep%sfcwater_energy (k)
      targetp%sfcwater_depth  (k) = sourcep%sfcwater_depth  (k)
      targetp%sfcwater_tempk  (k) = sourcep%sfcwater_tempk  (k)
      targetp%sfcwater_fracliq(k) = sourcep%sfcwater_fracliq(k)
   end do

   do k=1,cpatch%ncohorts
      targetp%leaf_resolvable (k) = sourcep%leaf_resolvable (k)
      targetp%leaf_energy     (k) = sourcep%leaf_energy     (k)
      targetp%leaf_water      (k) = sourcep%leaf_water      (k)
      targetp%leaf_temp       (k) = sourcep%leaf_temp       (k)
      targetp%leaf_fliq       (k) = sourcep%leaf_fliq       (k)
      targetp%leaf_hcap       (k) = sourcep%leaf_hcap       (k)
      targetp%leaf_reynolds   (k) = sourcep%leaf_reynolds   (k)
      targetp%leaf_grashof    (k) = sourcep%leaf_grashof    (k)
      targetp%leaf_nussfree   (k) = sourcep%leaf_nussfree   (k)
      targetp%leaf_nussforc   (k) = sourcep%leaf_nussforc   (k)
      targetp%leaf_gbh        (k) = sourcep%leaf_gbh        (k)
      targetp%leaf_gbw        (k) = sourcep%leaf_gbw        (k)
      targetp%rshort_l        (k) = sourcep%rshort_l        (k)
      targetp%rlong_l         (k) = sourcep%rlong_l         (k)

      targetp%wood_resolvable (k) = sourcep%wood_resolvable (k)
      targetp%wood_energy     (k) = sourcep%wood_energy     (k)
      targetp%wood_water      (k) = sourcep%wood_water      (k)
      targetp%wood_temp       (k) = sourcep%wood_temp       (k)
      targetp%wood_fliq       (k) = sourcep%wood_fliq       (k)
      targetp%wood_hcap       (k) = sourcep%wood_hcap       (k)
      targetp%wood_reynolds   (k) = sourcep%wood_reynolds   (k)
      targetp%wood_grashof    (k) = sourcep%wood_grashof    (k)
      targetp%wood_nussfree   (k) = sourcep%wood_nussfree   (k)
      targetp%wood_nussforc   (k) = sourcep%wood_nussforc   (k)
      targetp%wood_gbh        (k) = sourcep%wood_gbh        (k)
      targetp%wood_gbw        (k) = sourcep%wood_gbw        (k)
      targetp%rshort_w        (k) = sourcep%rshort_w        (k)
      targetp%rlong_w         (k) = sourcep%rlong_w         (k)

      targetp%veg_wind        (k) = sourcep%veg_wind        (k)
      targetp%lint_shv        (k) = sourcep%lint_shv        (k)
      targetp%nplant          (k) = sourcep%nplant          (k)
      targetp%lai             (k) = sourcep%lai             (k)
      targetp%wai             (k) = sourcep%wai             (k)
      targetp%wpa             (k) = sourcep%wpa             (k)
      targetp%tai             (k) = sourcep%tai             (k)
      targetp%crown_area      (k) = sourcep%crown_area      (k)
      targetp%elongf          (k) = sourcep%elongf          (k)
      targetp%gsw_open        (k) = sourcep%gsw_open        (k)
      targetp%gsw_closed      (k) = sourcep%gsw_closed      (k)
      targetp%psi_open        (k) = sourcep%psi_open        (k)
      targetp%psi_closed      (k) = sourcep%psi_closed      (k)
      targetp%fs_open         (k) = sourcep%fs_open         (k)
      targetp%gpp             (k) = sourcep%gpp             (k)
      targetp%leaf_resp       (k) = sourcep%leaf_resp       (k)
      targetp%root_resp       (k) = sourcep%root_resp       (k)
      targetp%growth_resp     (k) = sourcep%growth_resp     (k)
      targetp%storage_resp    (k) = sourcep%storage_resp    (k)
      targetp%vleaf_resp      (k) = sourcep%vleaf_resp      (k)
   end do

   if (checkbudget) then
      targetp%co2budget_storage      = sourcep%co2budget_storage
      targetp%co2budget_loss2atm     = sourcep%co2budget_loss2atm
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
      targetp%avg_carbon_ac          = sourcep%avg_carbon_ac
      targetp%avg_vapor_lc           = sourcep%avg_vapor_lc
      targetp%avg_vapor_wc           = sourcep%avg_vapor_wc
      targetp%avg_dew_cg             = sourcep%avg_dew_cg
      targetp%avg_vapor_gc           = sourcep%avg_vapor_gc
      targetp%avg_wshed_vg           = sourcep%avg_wshed_vg
      targetp%avg_intercepted        = sourcep%avg_intercepted
      targetp%avg_throughfall        = sourcep%avg_throughfall
      targetp%avg_vapor_ac           = sourcep%avg_vapor_ac
      targetp%avg_transp             = sourcep%avg_transp
      targetp%avg_evap               = sourcep%avg_evap
      targetp%avg_rshort_gnd         = sourcep%avg_rshort_gnd
      targetp%avg_rlong_gnd          = sourcep%avg_rlong_gnd
      targetp%avg_sensible_lc        = sourcep%avg_sensible_lc
      targetp%avg_sensible_wc        = sourcep%avg_sensible_wc
      targetp%avg_qwshed_vg          = sourcep%avg_qwshed_vg
      targetp%avg_qintercepted       = sourcep%avg_qintercepted
      targetp%avg_qthroughfall       = sourcep%avg_qthroughfall
      targetp%avg_sensible_gc        = sourcep%avg_sensible_gc
      targetp%avg_sensible_ac        = sourcep%avg_sensible_ac
      targetp%avg_drainage           = sourcep%avg_drainage
      targetp%avg_drainage_heat      = sourcep%avg_drainage_heat

      do k=rk4site%lsl,nzg
         targetp%avg_sensible_gg(k) = sourcep%avg_sensible_gg(k)
         targetp%avg_smoist_gg(k)   = sourcep%avg_smoist_gg(k)
         targetp%avg_transloss(k)   = sourcep%avg_transloss(k)
      end do
   end if

   if (print_detailed) then
      targetp%flx_carbon_ac          = sourcep%flx_carbon_ac
      targetp%flx_vapor_lc           = sourcep%flx_vapor_lc
      targetp%flx_vapor_wc           = sourcep%flx_vapor_wc
      targetp%flx_dew_cg             = sourcep%flx_dew_cg
      targetp%flx_vapor_gc           = sourcep%flx_vapor_gc
      targetp%flx_wshed_vg           = sourcep%flx_wshed_vg
      targetp%flx_intercepted        = sourcep%flx_intercepted
      targetp%flx_throughfall        = sourcep%flx_throughfall
      targetp%flx_vapor_ac           = sourcep%flx_vapor_ac
      targetp%flx_transp             = sourcep%flx_transp
      targetp%flx_evap               = sourcep%flx_evap
      targetp%flx_rshort_gnd         = sourcep%flx_rshort_gnd
      targetp%flx_rlong_gnd          = sourcep%flx_rlong_gnd
      targetp%flx_sensible_lc        = sourcep%flx_sensible_lc
      targetp%flx_sensible_wc        = sourcep%flx_sensible_wc
      targetp%flx_qwshed_vg          = sourcep%flx_qwshed_vg
      targetp%flx_qintercepted       = sourcep%flx_qintercepted
      targetp%flx_qthroughfall       = sourcep%flx_qthroughfall
      targetp%flx_sensible_gc        = sourcep%flx_sensible_gc
      targetp%flx_sensible_ac        = sourcep%flx_sensible_ac
      targetp%flx_drainage           = sourcep%flx_drainage
      targetp%flx_drainage_heat      = sourcep%flx_drainage_heat

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



   return
end subroutine copy_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will perform the allocation for the Runge-Kutta integrator structure, !
! and initialize it as well.                                                               !
!------------------------------------------------------------------------------------------!
subroutine initialize_rk4patches(init)
   use soil_coms     , only : nzg                   & ! intent(in)
                            , nzs                   ! ! intent(in)
   use ed_state_vars , only : edgrid_g              & ! intent(inout)
                            , edtype                & ! structure
                            , polygontype           & ! structure
                            , sitetype              & ! structure
                            , patchtype             ! ! structure
   use rk4_coms      , only : integration_buff      & ! structure
                            , deallocate_rk4_coh    & ! structure
                            , allocate_rk4_patch    & ! structure
                            , allocate_rk4_coh      & ! structure
                            , allocate_rk4_aux      ! ! structure
   use ed_misc_coms  , only : integration_scheme    ! ! intent(in)
   use grid_coms     , only : ngrids                ! ! intent(in)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   logical           , intent(in) :: init
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)      , pointer    :: cgrid
   type(polygontype) , pointer    :: cpoly
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   integer                        :: maxcohort
   integer                        :: igr
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   !---------------------------------------------------------------------------------------!

   if (init) then
      !------------------------------------------------------------------------------------!
      !     If this is initialization, make sure soil and sfcwater arrays are allocated.   !
      !------------------------------------------------------------------------------------!
      allocate(integration_buff%initp )
      allocate(integration_buff%yscal )
      allocate(integration_buff%y     )
      allocate(integration_buff%dydx  )
      allocate(integration_buff%yerr  )
      allocate(integration_buff%ytemp )

      call allocate_rk4_patch(integration_buff%initp )
      call allocate_rk4_patch(integration_buff%yscal )
      call allocate_rk4_patch(integration_buff%y     )
      call allocate_rk4_patch(integration_buff%dydx  )
      call allocate_rk4_patch(integration_buff%yerr  )
      call allocate_rk4_patch(integration_buff%ytemp )

      !------ Allocate and initialise the auxiliary structure. ----------------------------!
      call allocate_rk4_aux(nzg,nzs)


      !------------------------------------------------------------------------------------!
      !     The following structures are allocated/deallocated depending on the            !
      ! integration method.                                                                !
      !------------------------------------------------------------------------------------!
      select case(integration_scheme) 
      case (0) !----- Euler. --------------------------------------------------------------!
         allocate(integration_buff%dinitp)

         call allocate_rk4_patch(integration_buff%dinitp)

      case (1) !----- Runge-Kutta. --------------------------------------------------------!
         allocate(integration_buff%ak2)
         allocate(integration_buff%ak3)
         allocate(integration_buff%ak4)
         allocate(integration_buff%ak5)
         allocate(integration_buff%ak6)
         allocate(integration_buff%ak7)

         call allocate_rk4_patch(integration_buff%ak2)
         call allocate_rk4_patch(integration_buff%ak3)
         call allocate_rk4_patch(integration_buff%ak4)
         call allocate_rk4_patch(integration_buff%ak5)
         call allocate_rk4_patch(integration_buff%ak6)
         call allocate_rk4_patch(integration_buff%ak7)

      case (2) !----- Heun's. -------------------------------------------------------------!
         allocate(integration_buff%ak2)
         allocate(integration_buff%ak3)

         call allocate_rk4_patch(integration_buff%ak2)
         call allocate_rk4_patch(integration_buff%ak3)
      end select
      !------------------------------------------------------------------------------------!
   else
      !------------------------------------------------------------------------------------!
      !    If this is not initialization, deallocate cohort memory from integration        !
      ! patches.                                                                           !
      !------------------------------------------------------------------------------------!
      call deallocate_rk4_coh(integration_buff%initp )
      call deallocate_rk4_coh(integration_buff%yscal )
      call deallocate_rk4_coh(integration_buff%y     )
      call deallocate_rk4_coh(integration_buff%dydx  )
      call deallocate_rk4_coh(integration_buff%yerr  )
      call deallocate_rk4_coh(integration_buff%ytemp )

      !------------------------------------------------------------------------------------!
      !     The following structures are allocated/deallocated depending on the            !
      ! integration method.                                                                !
      !------------------------------------------------------------------------------------!
      select case(integration_scheme) 
      case (0) !----- Euler. --------------------------------------------------------------!
         call deallocate_rk4_coh(integration_buff%dinitp)
      case (1) !----- Runge-Kutta. --------------------------------------------------------!
         call deallocate_rk4_coh(integration_buff%ak2)
         call deallocate_rk4_coh(integration_buff%ak3)
         call deallocate_rk4_coh(integration_buff%ak4)
         call deallocate_rk4_coh(integration_buff%ak5)
         call deallocate_rk4_coh(integration_buff%ak6)
         call deallocate_rk4_coh(integration_buff%ak7)
      case (2) !----- Heun's. -------------------------------------------------------------!
         call deallocate_rk4_coh(integration_buff%ak2)
         call deallocate_rk4_coh(integration_buff%ak3)
      end select
      !------------------------------------------------------------------------------------!
   end if

   !----- Find maximum number of cohorts amongst all patches ------------------------------!
   maxcohort = 1
   do igr = 1,ngrids
      cgrid => edgrid_g(igr)
      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)
               maxcohort = max(maxcohort,cpatch%ncohorts)
            end do
         end do
      end do
   end do
   ! write (unit=*,fmt='(a,1x,i5)') 'Maxcohort = ',maxcohort

   !----- Create new memory in each of the integration patches. ---------------------------!
   call allocate_rk4_coh(maxcohort,integration_buff%initp )
   call allocate_rk4_coh(maxcohort,integration_buff%yscal )
   call allocate_rk4_coh(maxcohort,integration_buff%y     )
   call allocate_rk4_coh(maxcohort,integration_buff%dydx  )
   call allocate_rk4_coh(maxcohort,integration_buff%yerr  )
   call allocate_rk4_coh(maxcohort,integration_buff%ytemp )

   !---------------------------------------------------------------------------------------!
   !     The following structures are allocated/deallocated depending on the integration   !
   ! method.                                                                               !
   !---------------------------------------------------------------------------------------!
   select case(integration_scheme) 
   case (0) !----- Euler. -----------------------------------------------------------------!
      call allocate_rk4_coh(maxcohort,integration_buff%dinitp)
   case (1) !----- Runge-Kutta. -----------------------------------------------------------!
      call allocate_rk4_coh(maxcohort,integration_buff%ak2   )
      call allocate_rk4_coh(maxcohort,integration_buff%ak3   )
      call allocate_rk4_coh(maxcohort,integration_buff%ak4   )
      call allocate_rk4_coh(maxcohort,integration_buff%ak5   )
      call allocate_rk4_coh(maxcohort,integration_buff%ak6   )
      call allocate_rk4_coh(maxcohort,integration_buff%ak7   )
   case (2) !----- Heun's. ----------------------------------------------------------------!
      call allocate_rk4_coh(maxcohort,integration_buff%ak2   )
      call allocate_rk4_coh(maxcohort,integration_buff%ak3   )
   end select

   return
end subroutine initialize_rk4patches
!==========================================================================================!
!==========================================================================================!

