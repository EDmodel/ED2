module rk4_integ_utils
   contains

   !=======================================================================================!
   !=======================================================================================!
   ! Subroutine odeint                                                                     !
   !                                                                                       !
   !     This subroutine will drive the integration of several ODEs that drive the fast-   !
   ! -scale state variables.                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine odeint(csite,ipa,isi,ibuff,nsteps)

      use ed_state_vars  , only : sitetype               & ! structure
                                , patchtype              & ! structure
                                , polygontype            ! ! structure
      use rk4_coms       , only : integration_vars       & ! structure
                                , rk4patchtype           & ! structure
                                , integration_buff       & ! intent(inout)
                                , maxstp                 & ! intent(in)
                                , tbeg                   & ! intent(in)
                                , tend                   & ! intent(in)
                                , dtrk4                  & ! intent(in)
                                , dtrk4i                 & ! intent(in)
                                , tiny_offset            & ! intent(in)
                                , checkbudget            & ! intent(in)
                                , norm_rk4_fluxes        ! ! sub-routine
      use rk4_copy_patch , only : copy_rk4_patch         ! ! sub-routine
      use rk4_derivs     , only : leaf_derivs            ! ! sub-routine
      use rk4_misc       , only : adjust_sfcw_properties & ! sub-routine
                                , update_diagnostic_vars & ! sub-routine
                                , update_density_vars    & ! sub-routine
                                , print_rk4patch         ! ! sub-routine
      use ed_misc_coms   , only : fast_diagnostics       ! ! intent(in)
      use grid_coms      , only : nzg                    & ! intent(in)
                                , nzs                    ! ! intent(in)
      use soil_coms      , only : runoff_time_i          & ! intent(in)
                                , simplerunoff           ! ! intent(in)
      use consts_coms    , only : wdnsi8                 ! ! intent(in)
      use therm_lib8     , only : tl2uint8               ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)            , target      :: csite   !< Current site
      integer                   , intent(in)  :: ipa     !< Current patch ID
      integer                   , intent(in)  :: isi     !< Current site ID
      integer                   , intent(in)  :: ibuff   !< Current thread ID
      integer                   , intent(out) :: nsteps  !< Number of steps taken.
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)           , pointer     :: cpatch  !< Current patch
      type(rk4patchtype)        , pointer     :: initp   !< Current RK4 patch
      type(rk4patchtype)        , pointer     :: ylast   !< Last RK4 patch
      integer                                 :: i       !< Step counter
      integer                                 :: ksn     !< # of snow/water layers
      real(kind=8)                            :: x       !< Elapsed time
      real(kind=8)                            :: hgoal   !< Current delta-t attempt
      real(kind=8)                            :: hbeg    !< Initial delta-t
      real(kind=8)                            :: h       !< Current delta-t attempt
      real(kind=8)                            :: hnext   !< Next delta-t
      real(kind=8)                            :: hdid    !< delta-t that worked (???)
      real(kind=8)                            :: qwfree  !< Free water internal energy
      real(kind=8)                            :: wfreeb  !< Free water 
      !----- External function. -----------------------------------------------------------!
      real                      , external    :: sngloff
      !------------------------------------------------------------------------------------!



      cpatch => csite%patch(ipa)


      !------------------------------------------------------------------------------------!
      !      Initial step size.  Experience has shown that giving this too large a value   !
      ! causes the integrator to fail (e.g., soil layers become supersaturated).           !
      !------------------------------------------------------------------------------------!
      hbeg = dble(csite%htry(ipa))

      !------------------------------------------------------------------------------------!
      !     Copy the initial patch to the one we use for integration.                      !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(integration_buff(ibuff)%initp, integration_buff(ibuff)%y,cpatch)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Set initial time and stepsize.                                                     !
      !------------------------------------------------------------------------------------!
      x = tbeg
      h = hbeg
      if (dtrk4 < 0.d0) h = -hbeg
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! Begin timestep loop                                                                !
      !------------------------------------------------------------------------------------!
      timesteploop: do i=1,maxstp

         !----- Get initial derivatives ---------------------------------------------------!
         call leaf_derivs(integration_buff(ibuff)%y,integration_buff(ibuff)%dydx,csite,ipa &
                         ,ibuff,h,.false.)

         !----- Get scalings used to determine stability ----------------------------------!
         call get_yscal(integration_buff(ibuff)%y, integration_buff(ibuff)%dydx,h          &
                       ,integration_buff(ibuff)%yscal,cpatch)

         !---------------------------------------------------------------------------------!
         !     Be sure not to overstep, but save the intended step to control the next     !
         ! step and thus avoid unnecessary shrinking.                                      !
         !---------------------------------------------------------------------------------!
         hgoal = h
         if ((x+h-tend)*(x+h-tbeg) > 0.d0) h = tend - x
         !---------------------------------------------------------------------------------!


         !----- Take the step -------------------------------------------------------------!
         call rkqs(x,h,hgoal,hdid,hnext,csite,ipa,isi,ibuff)
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
      write (unit=*,fmt='(a)') ' ==> Too many steps in routine odeint'
      call print_rk4patch(integration_buff(ibuff)%y, csite,ipa)
      !------------------------------------------------------------------------------------!

      return
   end subroutine odeint
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine copies the meteorological variables to the Runge-Kutta buffer.     !
   ! This is to ensure all variables are in double precision, so consistent with the       !
   ! buffer variables.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine copy_met_2_rk4site(mzg,atm_ustar,atm_theiv,atm_vpdef,atm_theta,atm_tmp       &
                                ,atm_shv,atm_co2,zoff,exner,pcpg,qpcpg,dpcpg,prss,rshort   &
                                ,rlong,par_beam,par_diffuse,nir_beam,nir_diffuse,geoht     &
                                ,lsl,ntext_soil,green_leaf_factor,lon,lat,cosz)
      use ed_max_dims    , only : n_pft         ! ! intent(in)
      use rk4_coms       , only : rk4site       ! ! structure
      use canopy_air_coms, only : ustmin8       ! ! intent(in)
      use therm_lib8     , only : rehuil8       & ! function
                                , reducedpress8 & ! function
                                , tq2enthalpy8  & ! function
                                , press2exner8  & ! function
                                , extheta2temp8 & ! function
                                , idealdenssh8  ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer                  , intent(in) :: mzg
      integer                  , intent(in) :: lsl
      real                     , intent(in) :: atm_ustar
   !   real                     , intent(in) :: vels
      real                     , intent(in) :: atm_theiv
      real                     , intent(in) :: atm_vpdef
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
      real                     , intent(in) :: par_beam
      real                     , intent(in) :: par_diffuse
      real                     , intent(in) :: nir_beam
      real                     , intent(in) :: nir_diffuse
      real                     , intent(in) :: geoht
      integer, dimension(mzg)  , intent(in) :: ntext_soil
      real   , dimension(n_pft), intent(in) :: green_leaf_factor
      real                     , intent(in) :: lon
      real                     , intent(in) :: lat
      real                     , intent(in) :: cosz
      !------------------------------------------------------------------------------------!

      
      !----- Copy the integer variables. --------------------------------------------------!
      rk4site%lsl               = lsl
      rk4site%ntext_soil(:)     = 0
      rk4site%ntext_soil(1:mzg) = ntext_soil(1:mzg)

      !----- Convert to double precision. -------------------------------------------------!
      rk4site%atm_theiv             = dble(atm_theiv           )
      rk4site%atm_vpdef             = dble(atm_vpdef           )
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
      rk4site%par_beam              = dble(par_beam            )
      rk4site%par_diffuse           = dble(par_diffuse         )
      rk4site%nir_beam              = dble(nir_beam            )
      rk4site%nir_diffuse           = dble(nir_diffuse         )
      rk4site%geoht                 = dble(geoht               )
      rk4site%lon                   = dble(lon                 )
      rk4site%lat                   = dble(lat                 )
      rk4site%cosz                  = dble(cosz                )
      rk4site%green_leaf_factor(:)  = dble(green_leaf_factor(:))
      !------------------------------------------------------------------------------------!



      !----- Find the other variables that require a little math. -------------------------!
      rk4site%atm_ustar = max(ustmin8,dble(atm_ustar))
      rk4site%atm_rhv   = rehuil8(rk4site%atm_prss,rk4site%atm_tmp,rk4site%atm_shv,.true.)
      rk4site%atm_rhos  = idealdenssh8(rk4site%atm_prss,rk4site%atm_tmp,rk4site%atm_shv)
      !------------------------------------------------------------------------------------!


      return
   end subroutine copy_met_2_rk4site
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutines increment the derivative into the previous guess to create the    !
   ! new guess.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine inc_rk4_patch(rkp, inc, fac, cpatch)
      use ed_state_vars , only : sitetype           & ! structure
                               , patchtype          ! ! structure
      use rk4_coms      , only : rk4patchtype       & ! structure
                               , rk4site            & ! intent(in)
                               , checkbudget        & ! intent(in)
                               , print_detailed     ! ! intent(in)
      use grid_coms     , only : nzg                ! ! intent(in)
      use ed_misc_coms  , only : fast_diagnostics   ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: rkp    ! Temporary patch with previous state
      type(rk4patchtype) , target     :: inc    ! Temporary patch with its derivatives
      type(patchtype)    , target     :: cpatch ! Current patch (for characteristics)
      real(kind=8)       , intent(in) :: fac    ! Increment factor
      !----- Local variables --------------------------------------------------------------!
      integer                         :: ico    ! Cohort ID
      integer                         :: k      ! Counter
      !------------------------------------------------------------------------------------!

      rkp%can_enthalpy  = rkp%can_enthalpy + fac * inc%can_enthalpy
      rkp%can_shv       = rkp%can_shv      + fac * inc%can_shv
      rkp%can_co2       = rkp%can_co2      + fac * inc%can_co2

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

      rkp%water_deficit   = rkp%water_deficit   + fac * inc%water_deficit

      rkp%upwp = rkp%upwp + fac * inc%upwp
      rkp%wpwp = rkp%wpwp + fac * inc%wpwp
      rkp%tpwp = rkp%tpwp + fac * inc%tpwp
      rkp%qpwp = rkp%qpwp + fac * inc%qpwp
      rkp%cpwp = rkp%cpwp + fac * inc%cpwp

      do ico = 1,cpatch%ncohorts
         rkp%leaf_water    (ico) = rkp%leaf_water    (ico) + fac * inc%leaf_water    (ico)
         rkp%leaf_water_im2(ico) = rkp%leaf_water_im2(ico) + fac * inc%leaf_water_im2(ico)
         rkp%leaf_energy   (ico) = rkp%leaf_energy   (ico) + fac * inc%leaf_energy   (ico)
         rkp%wood_water    (ico) = rkp%wood_water    (ico) + fac * inc%wood_water    (ico)
         rkp%wood_water_im2(ico) = rkp%wood_water_im2(ico) + fac * inc%wood_water_im2(ico)
         rkp%wood_energy   (ico) = rkp%wood_energy   (ico) + fac * inc%wood_energy   (ico)
         rkp%veg_water     (ico) = rkp%veg_water     (ico) + fac * inc%veg_water     (ico)
         rkp%veg_water_im2 (ico) = rkp%veg_water_im2 (ico) + fac * inc%veg_water_im2 (ico)
         rkp%veg_energy    (ico) = rkp%veg_energy    (ico) + fac * inc%veg_energy    (ico)


         rkp%psi_open  (ico) = rkp%psi_open  (ico) + fac * inc%psi_open  (ico)
         rkp%psi_closed(ico) = rkp%psi_closed(ico) + fac * inc%psi_closed(ico)
      end do

      if (checkbudget) then

         rkp%co2budget_storage      = rkp%co2budget_storage                                &
                                    + fac * inc%co2budget_storage
         rkp%co2budget_loss2atm     = rkp%co2budget_loss2atm                               &
                                    + fac * inc%co2budget_loss2atm

         rkp%wbudget_storage       = rkp%wbudget_storage       + fac * inc%wbudget_storage
         rkp%wbudget_loss2atm      = rkp%wbudget_loss2atm      + fac * inc%wbudget_loss2atm
         rkp%wbudget_loss2drainage = rkp%wbudget_loss2drainage                             &
                                   + fac * inc%wbudget_loss2drainage

         rkp%ebudget_storage       = rkp%ebudget_storage       + fac * inc%ebudget_storage
         rkp%ebudget_netrad        = rkp%ebudget_netrad        + fac * inc%ebudget_netrad
         rkp%ebudget_loss2atm      = rkp%ebudget_loss2atm      + fac * inc%ebudget_loss2atm
         rkp%ebudget_loss2drainage = rkp%ebudget_loss2drainage                             &
                                   + fac * inc%ebudget_loss2drainage
      end if
      if (fast_diagnostics) then
         rkp%avg_ustar          = rkp%avg_ustar          + fac * inc%avg_ustar
         rkp%avg_tstar          = rkp%avg_tstar          + fac * inc%avg_tstar
         rkp%avg_qstar          = rkp%avg_qstar          + fac * inc%avg_qstar
         rkp%avg_cstar          = rkp%avg_cstar          + fac * inc%avg_cstar


         rkp%avg_carbon_ac      = rkp%avg_carbon_ac      + fac * inc%avg_carbon_ac
         rkp%avg_carbon_st      = rkp%avg_carbon_st      + fac * inc%avg_carbon_st

         rkp%avg_throughfall    = rkp%avg_throughfall    + fac * inc%avg_throughfall
         rkp%avg_vapor_ac       = rkp%avg_vapor_ac       + fac * inc%avg_vapor_ac
         rkp%avg_vapor_gc       = rkp%avg_vapor_gc       + fac * inc%avg_vapor_gc
         rkp%avg_drainage       = rkp%avg_drainage       + fac * inc%avg_drainage
         rkp%avg_qdrainage      = rkp%avg_qdrainage      + fac * inc%avg_qdrainage
         rkp%avg_qthroughfall   = rkp%avg_qthroughfall   + fac * inc%avg_qthroughfall
         rkp%avg_sensible_gc    = rkp%avg_sensible_gc    + fac * inc%avg_sensible_gc
         rkp%avg_sensible_ac    = rkp%avg_sensible_ac    + fac * inc%avg_sensible_ac

         do k=rk4site%lsl,nzg
            rkp%avg_sensible_gg(k)  = rkp%avg_sensible_gg(k)                               &
                                    + fac * inc%avg_sensible_gg(k)
            rkp%avg_smoist_gg(k)    = rkp%avg_smoist_gg(k) + fac * inc%avg_smoist_gg(k)  
            rkp%avg_transloss(k)    = rkp%avg_transloss(k) + fac * inc%avg_transloss(k)  
         end do


         do ico=1,cpatch%ncohorts
            rkp%avg_sensible_lc        (ico) =       rkp%avg_sensible_lc     (ico)         &
                                             + fac * inc%avg_sensible_lc     (ico)
            rkp%avg_sensible_wc        (ico) =       rkp%avg_sensible_wc     (ico)         &
                                             + fac * inc%avg_sensible_wc     (ico)
            rkp%avg_vapor_lc           (ico) =       rkp%avg_vapor_lc        (ico)         &
                                             + fac * inc%avg_vapor_lc        (ico)
            rkp%avg_vapor_wc           (ico) =       rkp%avg_vapor_wc        (ico)         &
                                             + fac * inc%avg_vapor_wc        (ico)
            rkp%avg_transp             (ico) =       rkp%avg_transp          (ico)         &
                                             + fac * inc%avg_transp          (ico)
            rkp%avg_intercepted_al     (ico) =       rkp%avg_intercepted_al  (ico)         &
                                             + fac * inc%avg_intercepted_al  (ico)
            rkp%avg_intercepted_aw     (ico) =       rkp%avg_intercepted_aw  (ico)         &
                                             + fac * inc%avg_intercepted_aw  (ico)
            rkp%avg_wshed_lg           (ico) =       rkp%avg_wshed_lg        (ico)         &
                                             + fac * inc%avg_wshed_lg        (ico)
            rkp%avg_wshed_wg           (ico) =       rkp%avg_wshed_wg        (ico)         &
                                             + fac * inc%avg_wshed_wg        (ico)
            rkp%avg_wflux_wl           (ico) =       rkp%avg_wflux_wl        (ico)         &
                                             + fac * inc%avg_wflux_wl        (ico)
            rkp%avg_wflux_gw           (ico) =       rkp%avg_wflux_gw        (ico)         &
                                             + fac * inc%avg_wflux_gw        (ico)
            do k=rk4site%lsl,nzg
               rkp%avg_wflux_gw_layer(k,ico) =       rkp%avg_wflux_gw_layer(k,ico)         &
                                             + fac * inc%avg_wflux_gw_layer(k,ico)
            end do

         end do

      end if

      !------------------------------------------------------------------------------------!
      !    Increment the instantaneous fluxes.  The derivative term should be the same as  !
      ! the full fluxes, the only difference is that these variables are normalised and    !
      ! re-set after each time step.                                                       !
      !------------------------------------------------------------------------------------!
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
            rkp%flx_sensible_gg(k) = rkp%flx_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
            rkp%flx_smoist_gg  (k) = rkp%flx_smoist_gg  (k)  + fac * inc%avg_smoist_gg  (k)
            rkp%flx_transloss  (k) = rkp%flx_transloss  (k)  + fac * inc%avg_transloss  (k)
         end do

         do ico = 1,cpatch%ncohorts
            rkp%flx_vapor_lc                  =         rkp%flx_vapor_lc                   &
                                              + fac *   inc%avg_vapor_lc         (ico)
            rkp%flx_vapor_wc                  =         rkp%flx_vapor_wc                   &
                                              + fac *   inc%avg_vapor_wc         (ico)
            rkp%flx_wshed_vg                  =         rkp%flx_wshed_vg                   &
                                              + fac * ( inc%avg_wshed_lg         (ico)     &
                                                      + inc%avg_wshed_wg         (ico) )
            rkp%flx_wflux_wl                  =         rkp%flx_wflux_wl                   &
                                              + fac *   inc%avg_wflux_wl         (ico)
            rkp%flx_wflux_gw                  =         rkp%flx_wflux_gw                   &
                                              + fac *   inc%avg_wflux_gw         (ico)
            rkp%flx_intercepted               =       rkp%flx_intercepted                  &
                                              + fac * ( inc%avg_intercepted_al   (ico)     &
                                                      + inc%avg_intercepted_aw   (ico) )
            rkp%flx_sensible_lc               =         rkp%flx_sensible_lc                &
                                              + fac *   inc%avg_sensible_lc      (ico)
            rkp%flx_sensible_wc               =         rkp%flx_sensible_wc                &
                                              + fac *   inc%avg_sensible_wc      (ico)
            rkp%flx_qwshed_vg                 =         rkp%flx_qwshed_vg                  &
                                              + fac *   inc%cfx_qwshed           (ico)
            rkp%flx_qintercepted              =         rkp%flx_qintercepted               &
                                              + fac *   inc%cfx_qintercepted     (ico)
            rkp%cfx_hflxlc              (ico) =         rkp%cfx_hflxlc           (ico)     &
                                              + fac *   inc%cfx_hflxlc           (ico)
            rkp%cfx_hflxwc              (ico) =         rkp%cfx_hflxwc           (ico)     &
                                              + fac *   inc%cfx_hflxwc           (ico)
            rkp%cfx_qwflxlc             (ico) =         rkp%cfx_qwflxlc          (ico)     &
                                              + fac *   inc%cfx_qwflxlc          (ico)
            rkp%cfx_qwflxwc             (ico) =         rkp%cfx_qwflxwc          (ico)     &
                                              + fac *   inc%cfx_qwflxwc          (ico)
            rkp%cfx_qwshed              (ico) =         rkp%cfx_qwshed           (ico)     &
                                              + fac *   inc%cfx_qwshed           (ico)
            rkp%cfx_qtransp             (ico) =         rkp%cfx_qtransp          (ico)     &
                                              + fac *   inc%cfx_qtransp          (ico)
            rkp%cfx_qintercepted        (ico) =         rkp%cfx_qintercepted     (ico)     &
                                              + fac *   inc%cfx_qintercepted     (ico)
            rkp%cfx_qwflux_wl           (ico) =         rkp%cfx_qwflux_wl        (ico)     &
                                              + fac *   inc%cfx_qwflux_wl        (ico)
            rkp%cfx_qwflux_gw           (ico) =         rkp%cfx_qwflux_gw        (ico)     &
                                              + fac *   inc%cfx_qwflux_gw        (ico)
            do k=rk4site%lsl,nzg
               rkp%cfx_qwflux_gw_layer(k,ico) =         rkp%cfx_qwflux_gw_layer(k,ico)     &
                                              + fac *   inc%cfx_qwflux_gw_layer(k,ico)
            end do
         end do

      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine inc_rk4_patch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine finds the error scale for the integrated variables, which will be  !
   ! later used to define the relative error.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine get_yscal(y,dy,htry,yscal,cpatch)
      use ed_state_vars        , only : patchtype             ! ! structure
      use rk4_coms             , only : rk4patchtype          & ! structure
                                      , rk4site               & ! intent(in)
                                      , ibranch_thermo        & ! intent(in)
                                      , tiny_offset           & ! intent(in)
                                      , huge_offset           & ! intent(in)
                                      , rk4water_stab_thresh  & ! intent(in)
                                      , rk4tiny_sfcw_mass     & ! intent(in)
                                      , rk4leaf_drywhc        & ! intent(in)
                                      , checkbudget           ! ! intent(in)
      use grid_coms            , only : nzg                   & ! intent(in)
                                      , nzs                   ! ! intent(in)
      use consts_coms          , only : wdnsi8                ! ! intent(in)
      use soil_coms            , only : isoilbc               & ! intent(in)
                                      , dslzi8                ! ! intent(in)
      use physiology_coms      , only : plant_hydro_scheme    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype), target     :: y      ! Struct. with the guesses
      type(rk4patchtype), target     :: dy     ! Struct. with their derivatives
      type(rk4patchtype), target     :: yscal  ! Struct. with their scales
      type(patchtype)   , target     :: cpatch ! Current patch
      real(kind=8)      , intent(in) :: htry   ! Time-step we are trying
      !----- Local variables --------------------------------------------------------------!
      real(kind=8) :: meanscale_sfcw_mass   ! Average Sfc. water mass scale
      real(kind=8) :: meanscale_sfcw_energy ! Average Sfc. water en. scale
      real(kind=8) :: meanscale_sfcw_depth  ! Average Sfc. water depth scale
      integer      :: k                     ! Counter
      integer      :: ico                   ! Current cohort ID
      !------------------------------------------------------------------------------------!

      yscal%can_enthalpy =  abs(y%can_enthalpy) + abs(dy%can_enthalpy * htry)
      yscal%can_shv      =  abs(y%can_shv     ) + abs(dy%can_shv      * htry)
      yscal%can_co2      =  abs(y%can_co2     ) + abs(dy%can_co2      * htry)

      !------------------------------------------------------------------------------------!
      !     We don't solve pressure prognostically, so the scale cannot be computed based  !
      ! on the derivative.  Also, pressure is a variable with a large absolute variable    !
      ! and tiny variation, so the scale is given by the scale of the variability...       !
      !------------------------------------------------------------------------------------!
      yscal%can_prss    = abs(rk4site%atm_prss - y%can_prss)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     We determine the scale for all layers but the top soil, which will be done     !
      ! differently depending on the status of the temporary surface water.                !
      !------------------------------------------------------------------------------------!
      do k=1,rk4site%lsl-1
         yscal%soil_water (k) = huge_offset
         yscal%soil_energy(k) = huge_offset
      end do
      do k=rk4site%lsl,nzg-1
         yscal%soil_water (k) = abs(y%soil_water(k) ) + abs(dy%soil_water(k)  * htry)
         yscal%soil_energy(k) = abs(y%soil_energy(k)) + abs(dy%soil_energy(k) * htry)
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Temporary surface layers require a special approach. The number of layers may  !
      ! vary during the integration process, so we must make sure that all layers          !
      ! initially have a scale.  Also, if the total mass is small, we must be more         !
      ! tolerant to avoid overestimating the error because of the small size.              !
      !------------------------------------------------------------------------------------!
      select case(y%flag_sfcwater)
      case(0)
         !---------------------------------------------------------------------------------!
         !     No temporary surface water, we can skip the check as no layer will be       !
         ! created until the upcoming step.                                                !
         !---------------------------------------------------------------------------------!
         do k=1,nzs
            yscal%sfcwater_mass(k)   = huge_offset
            yscal%sfcwater_energy(k) = huge_offset
            yscal%sfcwater_depth(k)  = huge_offset
         end do
         !---------------------------------------------------------------------------------!


         !------ Soil scale won't be affected by the temporary surface water. -------------!
         yscal%soil_water (nzg) = abs(y%soil_water (nzg)) + abs(dy%soil_water (nzg) * htry)
         yscal%soil_energy(nzg) = abs(y%soil_energy(nzg)) + abs(dy%soil_energy(nzg) * htry)
         !---------------------------------------------------------------------------------!

      case(1)
         !---------------------------------------------------------------------------------!
         !    Low stability threshold, there can't be more than one layer, and the energy  !
         ! will be solved together with the top soil layer.  Therefore, we skip energy     !
         ! check, but attribute the scale for mass and depth.                              !
         !---------------------------------------------------------------------------------!
         yscal%sfcwater_mass  (1) = abs(y%sfcwater_mass   (1)       )                      &
                                  + abs(dy%sfcwater_mass  (1) * htry)
         yscal%sfcwater_depth (1) = abs(y%sfcwater_depth  (1)       )                      &
                                  + abs(dy%sfcwater_depth (1) * htry)
         yscal%sfcwater_energy(1) = huge_offset
         do k=2,nzs
            yscal%sfcwater_mass  (k) = huge_offset
            yscal%sfcwater_energy(k) = huge_offset
            yscal%sfcwater_depth (k) = huge_offset
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Soil mass won't be affected by the temporary surface water, but soil energy !
         ! will include the energy associated with the temporary surface water energy.  In !
         ! reality the derivative of the temporary surface water will be zero, but we add  !
         ! for code consistency.                                                           !
         !---------------------------------------------------------------------------------!
         yscal%soil_water (nzg) = abs(y%soil_water (nzg)) + abs(dy%soil_water (nzg)*htry)
         yscal%soil_energy(nzg) = abs(y%soil_energy(nzg)) + abs(dy%soil_energy(nzg)*htry)  &
                                + dslzi8(nzg) * ( abs(y%sfcwater_energy(1))                &
                                                + abs(dy%sfcwater_energy(1) * htry) )
         !---------------------------------------------------------------------------------!

      case (2)
         !----- Computationally stable layer. ---------------------------------------------!
         meanscale_sfcw_mass   = 0.d0
         meanscale_sfcw_energy = 0.d0
         meanscale_sfcw_depth  = 0.d0
         do k=1,y%nlev_sfcwater
            yscal%sfcwater_mass  (k) = abs(y%sfcwater_mass   (k)       )                   &
                                     + abs(dy%sfcwater_mass  (k) * htry)
            yscal%sfcwater_energy(k) = abs(y%sfcwater_energy (k)       )                   &
                                     + abs(dy%sfcwater_energy(k) * htry)
            yscal%sfcwater_depth (k) = abs(y%sfcwater_depth  (k)       )                   &
                                     + abs(dy%sfcwater_depth (k) * htry)
         end do
         do k=y%nlev_sfcwater+1,nzs
            yscal%sfcwater_mass  (k) = huge_offset
            yscal%sfcwater_energy(k) = huge_offset
            yscal%sfcwater_depth (k) = huge_offset
         end do


         !------ Soil scale won't be affected by the temporary surface water. -------------!
         yscal%soil_water (nzg) = abs(y%soil_water (nzg)) + abs(dy%soil_water (nzg) * htry)
         yscal%soil_energy(nzg) = abs(y%soil_energy(nzg)) + abs(dy%soil_energy(nzg) * htry)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !----- Scale for the virtual water pools --------------------------------------------!
      if (abs(y%virtual_water) > 1.d-2*rk4water_stab_thresh) then
         yscal%virtual_water   = abs(y%virtual_water)  + abs(dy%virtual_water*htry)
         yscal%virtual_energy  = abs(y%virtual_energy) + abs(dy%virtual_energy*htry)
      elseif (abs(y%virtual_water) > rk4tiny_sfcw_mass) then
         yscal%virtual_water   = 1.d-2*rk4water_stab_thresh
         yscal%virtual_energy  = (yscal%virtual_water / abs(y%virtual_water))              &
                               * (abs(y%virtual_energy) + abs(dy%virtual_energy*htry))
      else
         yscal%virtual_water   = huge_offset
         yscal%virtual_energy  = huge_offset
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Scale for leaf, wood, and vegetation water and energy.  In case the plants have !
      ! few or no leaves, or the plant is buried in snow, we assign huge values for        !
      ! typical scale, thus preventing unecessary small steps.                             !
      !    Also, if the cohort has almost no water, make the scale less strict.            !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (1)
         !----- Combined leaf+wood solution. ----------------------------------------------!
         do ico=1,cpatch%ncohorts
            !----- Copy the logical tests. ------------------------------------------------!
            yscal%leaf_resolvable(ico) = y%leaf_resolvable(ico)
            yscal%wood_resolvable(ico) = y%wood_resolvable(ico)
            yscal%veg_resolvable (ico) = y%veg_resolvable (ico)
            yscal%is_small       (ico) = y%is_small       (ico)
            !------------------------------------------------------------------------------!



            !----- Find the scale only if we must solve veg. ------------------------------!
            if (y%veg_resolvable(ico)) then
               yscal%veg_energy(ico) = abs( y%veg_energy(ico))                             &
                                          + abs(dy%veg_energy(ico) * htry)
               yscal%veg_water(ico)  = max( abs(y%veg_water(ico))                          &
                                          + abs(dy%veg_water(ico)  * htry)                 &
                                          , rk4leaf_drywhc * y%tai(ico))  
            else
               yscal%veg_water(ico)  = huge_offset
               yscal%veg_energy(ico) = huge_offset
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Check plant hydraulics variables.                                       !
            !------------------------------------------------------------------------------!
            select case (plant_hydro_scheme)
            case (0)
               !---------------------------------------------------------------------------!
               !    Plant hydraulics is not enabled, make changes in internal water always !
               ! acceptable.                                                               !
               !---------------------------------------------------------------------------!
               yscal%veg_water_im2 (ico)    = huge_offset
               !---------------------------------------------------------------------------!
            case default
               !---------------------------------------------------------------------------!
               !    Simulation running with plant hydraulics.  Calculate the scale         !
               ! similarly to leaf/wood energy.                                            !
               !---------------------------------------------------------------------------!
               if (yscal%veg_resolvable(ico)) then
                  yscal%veg_water_im2(ico) = abs( y%veg_water_im2(ico))                    &
                                           + abs(dy%veg_water_im2(ico) * htry )
               else
                  yscal%veg_water_im2(ico) = huge_offset
               end if
               !---------------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !   No need to scale wood and leaf correctly, make it always acceptable.       !
            !------------------------------------------------------------------------------!
            yscal%leaf_water    (ico) = huge_offset
            yscal%leaf_water_im2(ico) = huge_offset
            yscal%leaf_energy   (ico) = huge_offset
            yscal%leaf_temp     (ico) = huge_offset
            yscal%wood_water    (ico) = huge_offset
            yscal%wood_water_im2(ico) = huge_offset
            yscal%wood_energy   (ico) = huge_offset
            yscal%wood_temp     (ico) = huge_offset
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!

      case (0,2)
         !---------------------------------------------------------------------------------!
         !      Either we are solving leaves only, or leaves and branches are treated as   !
         ! independent pools.                                                              !
         !---------------------------------------------------------------------------------!
         do ico = 1,cpatch%ncohorts
            !----- Copy the logical tests. ------------------------------------------------!
            yscal%leaf_resolvable(ico) = y%leaf_resolvable(ico)
            yscal%wood_resolvable(ico) = y%wood_resolvable(ico)
            yscal%veg_resolvable (ico) = y%veg_resolvable (ico)
            yscal%is_small       (ico) = y%is_small       (ico)
            !------------------------------------------------------------------------------!


            !----- Find the scale only if we must solve leaves. ---------------------------!
            if (y%leaf_resolvable(ico)) then
               yscal%leaf_energy(ico) = abs( y%leaf_energy(ico))                           &
                                           + abs(dy%leaf_energy(ico) * htry)
               yscal%leaf_temp(ico)   = abs( y%leaf_temp(ico))
               yscal%leaf_water(ico)  = max( abs(y%leaf_water(ico))                        &
                                           + abs(dy%leaf_water(ico)  * htry)               &
                                           , rk4leaf_drywhc * y%lai(ico))
            else
               yscal%leaf_water(ico)  = huge_offset
               yscal%leaf_energy(ico) = huge_offset
               yscal%leaf_temp(ico)   = huge_offset
            end if
            !------------------------------------------------------------------------------!


            !----- Find the scale only if we must solve wood. -----------------------------!
            if (y%wood_resolvable(ico)) then
               yscal%wood_energy(ico) = abs( y%wood_energy(ico))                           &
                                           + abs(dy%wood_energy(ico) * htry)
               yscal%wood_temp(ico)   = abs( y%wood_temp(ico))
               yscal%wood_water(ico)  = max( abs(y%wood_water(ico))                        &
                                           + abs(dy%wood_water(ico)  * htry)               &
                                           , rk4leaf_drywhc * y%wai(ico))
            else
               yscal%wood_water (ico) = huge_offset
               yscal%wood_energy(ico) = huge_offset
               yscal%wood_temp  (ico) = huge_offset
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Check plant hydraulics variables.                                       !
            !------------------------------------------------------------------------------!
            select case (plant_hydro_scheme)
            case (0)
               !---------------------------------------------------------------------------!
               !    Plant hydraulics is not enabled, make changes in internal water always !
               ! acceptable.                                                               !
               !---------------------------------------------------------------------------!
               yscal%leaf_water_im2(ico)    = huge_offset
               yscal%wood_water_im2(ico)    = huge_offset
               !---------------------------------------------------------------------------!
            case default
               !---------------------------------------------------------------------------!
               !    Simulation running with plant hydraulics.  Calculate the scale         !
               ! similarly to leaf/wood energy.                                            !
               !---------------------------------------------------------------------------!
               if (yscal%leaf_resolvable(ico)) then
                  yscal%leaf_water_im2(ico) = abs( y%leaf_water_im2(ico))                  &
                                            + abs(dy%leaf_water_im2(ico) * htry )
               else
                   yscal%leaf_water_im2(ico) = huge_offset
               end if
               if (yscal%wood_resolvable(ico)) then
                  yscal%wood_water_im2(ico) = abs( y%wood_water_im2(ico))                  &
                                            + abs(dy%wood_water_im2(ico) * htry)
               else
                  yscal%wood_water_im2(ico) = huge_offset
               end if
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !----- No need to scale veg correctly, let's make it always acceptable. -------!
            yscal%veg_water    (ico) = huge_offset
            yscal%veg_water_im2(ico) = huge_offset
            yscal%veg_energy   (ico) = huge_offset
            !------------------------------------------------------------------------------!
         end do
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Here we just need to make sure the user is checking mass, otherwise  these     !
      ! variables will not be computed at all.  If this turns out to be essential, we will !
      ! make this permanent and not dependent on checkbudget.  The only one that is not    !
      ! checked is the runoff, because it is computed after a step was accepted.           !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         !---------------------------------------------------------------------------------!
         !    If this is the very first time step, or if we are misfortuned, we may have a !
         ! situation in which the derivative is numerically zero, and making the check     !
         ! will become impossible because the scale would be ridiculously tiny, so we skip !
         ! the check this time and hope that everything will be alright next step.         !
         !---------------------------------------------------------------------------------!
         if (abs(y%co2budget_loss2atm)  < tiny_offset .and.                                &
             abs(dy%co2budget_loss2atm) < tiny_offset) then
            yscal%co2budget_loss2atm = 1.d-1
         else 
            yscal%co2budget_loss2atm = abs(y%co2budget_loss2atm)                           &
                                     + abs(dy%co2budget_loss2atm*htry)
            yscal%co2budget_loss2atm = max(yscal%co2budget_loss2atm,1.d-1)
         end if

         if (abs(y%ebudget_netrad)  < tiny_offset .and.                                    &
             abs(dy%ebudget_netrad) < tiny_offset) then
            yscal%ebudget_netrad  = 1.d0
         else 
            yscal%ebudget_netrad  = abs(y%ebudget_netrad)                                  &
                                  + abs(dy%ebudget_netrad*htry)
            yscal%ebudget_netrad = max(yscal%ebudget_netrad,1.d0)
         end if

         if (abs(y%ebudget_loss2atm)  < tiny_offset .and.                                  &
             abs(dy%ebudget_loss2atm) < tiny_offset) then
            yscal%ebudget_loss2atm = 1.d0
         else 
            yscal%ebudget_loss2atm = abs(y%ebudget_loss2atm)                               &
                                   + abs(dy%ebudget_loss2atm*htry)
            yscal%ebudget_loss2atm = max(yscal%ebudget_loss2atm,1.d0)
         end if

         if (abs(y%wbudget_loss2atm)  < tiny_offset .and.                                  &
             abs(dy%wbudget_loss2atm) < tiny_offset) then
            yscal%wbudget_loss2atm      = 1.d-6
         else 
            yscal%wbudget_loss2atm = abs(y%wbudget_loss2atm)                               &
                                   + abs(dy%wbudget_loss2atm*htry)
            yscal%wbudget_loss2atm = max(yscal%wbudget_loss2atm,1.d-6)
         end if

         if (abs(y%ebudget_storage)  < tiny_offset .and.                                   &
             abs(dy%ebudget_storage) < tiny_offset) then
            yscal%ebudget_storage = huge_offset
         else 
            yscal%ebudget_storage = abs(y%ebudget_storage)                                 &
                                  + abs(dy%ebudget_storage*htry)
         end if

         if (abs(y%co2budget_storage)  < tiny_offset .and.                                 &
             abs(dy%co2budget_storage) < tiny_offset) then
            yscal%co2budget_storage = huge_offset
         else 
            yscal%co2budget_storage = abs(y%co2budget_storage)                             &
                                    + abs(dy%co2budget_storage*htry)
         end if

         if (abs(y%wbudget_storage)  < tiny_offset .and.                                   &
             abs(dy%wbudget_storage) < tiny_offset) then
            yscal%wbudget_storage      = huge_offset
         else 
            yscal%wbudget_storage = abs(y%wbudget_storage)                                 &
                                  + abs(dy%wbudget_storage*htry)
         end if

         !---------------------------------------------------------------------------------!
         !     Drainage terms will be checked only if the boundary condition is free       !
         ! drainage.                                                                       !
         !---------------------------------------------------------------------------------!
         if (isoilbc == 0 .or. (abs(y%ebudget_loss2drainage)  < tiny_offset .and.          &
                                abs(dy%ebudget_loss2drainage) < tiny_offset)      ) then
            yscal%ebudget_loss2drainage = huge_offset
         else 
            yscal%ebudget_loss2drainage = abs(y%ebudget_loss2drainage)                     &
                                        + abs(dy%ebudget_loss2drainage*htry)
         end if
         if (isoilbc == 0 .or. (abs(y%wbudget_loss2drainage)  < tiny_offset .and.          &
                                abs(dy%wbudget_loss2drainage) < tiny_offset)      ) then
            yscal%wbudget_loss2drainage = huge_offset
         else 
            yscal%wbudget_loss2drainage = abs(y%wbudget_loss2drainage)                     &
                                        + abs(dy%wbudget_loss2drainage*htry)
         end if

      else 
         yscal%co2budget_storage       = huge_offset
         yscal%co2budget_loss2atm      = huge_offset
         yscal%ebudget_netrad          = huge_offset
         yscal%ebudget_loss2atm        = huge_offset
         yscal%wbudget_loss2atm        = huge_offset
         yscal%ebudget_storage         = huge_offset
         yscal%wbudget_storage         = huge_offset
         yscal%ebudget_loss2drainage   = huge_offset
         yscal%wbudget_loss2drainage   = huge_offset
      end if

      return
   end subroutine get_yscal
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine loops through the integrating variables, seeking for the largest   !
   ! error.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine get_errmax(errmax,yerr,yscal,cpatch,ytemp)

      use rk4_coms              , only : rk4patchtype       & ! structure
                                       , ibranch_thermo     & ! intent(in)
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
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target      :: yerr      ! Error structure
      type(rk4patchtype) , target      :: yscal     ! Scale structure
      type(rk4patchtype) , target      :: ytemp     ! Structure with attempted values
      type(patchtype)    , target      :: cpatch    ! Current patch
      real(kind=8)       , intent(out) :: errmax    ! Maximum error
      !----- Local variables --------------------------------------------------------------!
      integer                          :: ico       ! Current cohort ID
      real(kind=8)                     :: errh2o    ! Scratch error variable
      real(kind=8)                     :: errene    ! Scratch error variable
      real(kind=8)                     :: err       ! Scratch error variable
      real(kind=8)                     :: errh2oMAX ! Scratch error variable
      real(kind=8)                     :: erreneMAX ! Scratch error variable
      integer                          :: k         ! Counter
      !------------------------------------------------------------------------------------!

      !----- Initialise the error with an optimistic number... ----------------------------!
      errmax = 0.d0



      !------------------------------------------------------------------------------------!
      !    We now check each variable error, and keep track of the worst guess, which will !
      ! be our worst guess in the end.                                                     !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Get the canopy air space errors.  Only the ice-vapour equivalent potential     !
      ! temperature, water vapour mixing ratio and carbon dioxide mixing ratio are         !
      ! accounted.  Temperature and density will be also checked for sanity.               !
      !------------------------------------------------------------------------------------!
      err    = abs(yerr%can_enthalpy/yscal%can_enthalpy)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(1,1) = integ_err(1,1) + 1_8 

      err    = abs(yerr%can_shv/yscal%can_shv)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(3,1) = integ_err(3,1) + 1_8 

      err    = abs(yerr%can_co2/yscal%can_co2)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(5,1) = integ_err(5,1) + 1_8 
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Get the worst error only amongst the cohorts in which leaf or wood properties  !
      ! were computed.                                                                     !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (1)
         !---------------------------------------------------------------------------------!
         !     The combined leaf+branch pool is being solved. Check the veg variables      !
         ! only, but add the error to both leaf and wood.                                  !
         !---------------------------------------------------------------------------------!
         errh2oMAX  = 0.d0
         erreneMAX  = 0.d0
         do ico = 1,cpatch%ncohorts
            if (yscal%veg_resolvable(ico)) then
               errh2o     = abs(yerr%veg_water (ico) / yscal%veg_water (ico))
               errene     = abs(yerr%veg_energy(ico) / yscal%veg_energy(ico))
               errmax     = max(errmax,errh2o,errene)
               errh2oMAX  = max(errh2oMAX ,errh2o )
               erreneMAX  = max(erreneMAX ,errene )
            end if
         end do
         if(cpatch%ncohorts > 0 .and. record_err) then
            if (errh2oMAX  > rk4eps) then
               integ_err( 6,1) = integ_err( 6,1) + 1_8
               integ_err( 9,1) = integ_err( 9,1) + 1_8
            end if
            if (erreneMAX  > rk4eps) then
               integ_err( 7,1) = integ_err( 7,1) + 1_8
               integ_err(10,1) = integ_err(10,1) + 1_8
            end if
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Combined leaf+wood internal water pool                                      !
         !---------------------------------------------------------------------------------!
         ! leaf
         errh2oMAX  = 0.d0
         do ico = 1,cpatch%ncohorts
            if (yscal%veg_resolvable(ico)) then
               errh2o     = abs(yerr%veg_water_im2 (ico) / yscal%veg_water_im2 (ico))
               errmax     = max(errmax,errh2o)
               errh2oMAX  = max(errh2oMAX ,errh2o )
            end if
         end do
         if (cpatch%ncohorts > 0 .and. record_err .and. errh2oMAX  > rk4eps) then
            integ_err( 8,1) = integ_err( 8,1) + 1_8
            integ_err(11,1) = integ_err(11,1) + 1_8
         end if
         !---------------------------------------------------------------------------------!

      case default
         !---------------------------------------------------------------------------------!
         !     Either we are solving leaves only, or both leaf and branch pools are being  !
         ! solved independently. Check both leaf and wood variables.                       !
         !---------------------------------------------------------------------------------!
         !----- Leaves. -------------------------------------------------------------------!
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
            if (errh2oMAX  > rk4eps) integ_err(6,1) = integ_err(6,1) + 1_8
            if (erreneMAX  > rk4eps) integ_err(7,1) = integ_err(7,1) + 1_8
         end if
         !----- Wood. ---------------------------------------------------------------------!
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
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Leaf/wood internal water pool                                               !
         !---------------------------------------------------------------------------------!
         ! leaf
         errh2oMAX  = 0.d0
         do ico = 1,cpatch%ncohorts
            if (yscal%leaf_resolvable(ico)) then
               errh2o     = abs(yerr%leaf_water_im2 (ico) / yscal%leaf_water_im2 (ico))
               errmax     = max(errmax,errh2o)
               errh2oMAX  = max(errh2oMAX ,errh2o )
            end if
         end do
         if(cpatch%ncohorts > 0 .and. record_err) then
            if (errh2oMAX  > rk4eps) integ_err(8,1) = integ_err(8,1) + 1_8
         end if
         ! wood
         errh2oMAX  = 0.d0
         do ico = 1,cpatch%ncohorts
            if (yscal%wood_resolvable(ico)) then
               errh2o     = abs(yerr%wood_water_im2 (ico) / yscal%wood_water_im2 (ico))
               errmax     = max(errmax,errh2o)
               errh2oMAX  = max(errh2oMAX ,errh2o )
            end if
         end do
         if(cpatch%ncohorts > 0 .and. record_err) then
            if (errh2oMAX  > rk4eps) integ_err(11,1) = integ_err(11,1) + 1_8
         end if
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Virtual pool.                                                                  !
      !------------------------------------------------------------------------------------!
      err    = abs(yerr%virtual_energy/yscal%virtual_energy)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(12,1) = integ_err(12,1) + 1_8

      err    = abs(yerr%virtual_water/yscal%virtual_water)
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(13,1) = integ_err(13,1) + 1_8      
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Here we just need to make sure the user is checking mass, otherwise these      !
      ! variables will not be computed at all.  The only one that is not checked is the    !
      ! runoff, because it is computed only after a step was accepted.                     !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         err    = abs(yerr%co2budget_storage/yscal%co2budget_storage)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(14,1) = integ_err(14,1) + 1_8

         err    = abs(yerr%co2budget_loss2atm/yscal%co2budget_loss2atm)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(15,1) = integ_err(15,1) + 1_8

         err    = abs(yerr%ebudget_netrad/yscal%ebudget_netrad)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(16,1) = integ_err(16,1) + 1_8

         err    = abs(yerr%ebudget_loss2atm/yscal%ebudget_loss2atm)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(17,1) = integ_err(17,1) + 1_8

         err    = abs(yerr%wbudget_loss2atm/yscal%wbudget_loss2atm)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(18,1) = integ_err(18,1) + 1_8

         err    = abs(yerr%ebudget_loss2drainage/yscal%ebudget_loss2drainage)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(19,1) = integ_err(19,1) + 1_8

         err    = abs(yerr%wbudget_loss2drainage/yscal%wbudget_loss2drainage)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(20,1) = integ_err(20,1) + 1_8

         err    = abs(yerr%ebudget_storage/yscal%ebudget_storage)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(21,1) = integ_err(21,1) + 1_8

         err    = abs(yerr%wbudget_storage/yscal%wbudget_storage)
         errmax = max(errmax,err)
         if(record_err .and. err > rk4eps) integ_err(22,1) = integ_err(22,1) + 1_8
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Soil properties.  Notice that the index depends on the number of layers.       !
      !------------------------------------------------------------------------------------!
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
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Surface water/snow properties.  Notice that the index also depends on the      !
      ! number of layers.                                                                  !
      !------------------------------------------------------------------------------------!
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
      !------------------------------------------------------------------------------------!

      return
   end subroutine get_errmax
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will perform the allocation for the Runge-Kutta integrator         !
   ! structure, and initialize it as well.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine initialize_rk4patches(init)
      use soil_coms             , only : nzg                   & ! intent(in)
                                       , nzs                   ! ! intent(in)
      use ed_state_vars         , only : edgrid_g              & ! intent(inout)
                                       , edtype                & ! structure
                                       , polygontype           & ! structure
                                       , sitetype              & ! structure
                                       , patchtype             ! ! structure
      use rk4_coms              , only : integration_buff      & ! structure
                                       , rk4aux                & ! structure
                                       , deallocate_rk4_coh    & ! sub-routine
                                       , deallocate_rk4_aux    & ! sub-routine
                                       , allocate_rk4_patch    & ! sub-routine
                                       , allocate_rk4_coh      & ! sub-routine
                                       , allocate_rk4_aux      & ! sub-routine
                                       , allocate_bdf2_patch   & ! sub-routine
                                       , deallocate_bdf2_patch ! ! sub-routine
      use ed_para_coms          , only : nthreads              ! ! intent(in)
      use canopy_layer_coms     , only : tai_lyr_max           ! ! intent(in)
      use ed_misc_coms          , only : integration_scheme    & ! intent(in)
                                       , ibigleaf              ! ! intent(in)
      use grid_coms             , only : ngrids                ! ! intent(in)
      use c34constants          , only : thispft               & ! intent(out)
                                       , met                   & ! intent(out)
                                       , aparms                & ! intent(out)
                                       , stopen                & ! intent(out)
                                       , stclosed              & ! intent(out)
                                       , rubpsatlim            & ! intent(out)
                                       , thirdlim              & ! intent(out)
                                       , lightlim              ! ! intent(out)
      use canopy_radiation_coms , only : radscr                & ! intent(out)
                                       , alloc_radscratch      & ! sub-routine
                                       , dealloc_radscratch    & ! sub-routine
                                       , nullify_radscratch    ! ! sub-routine

      implicit none
      !----- Argument ---------------------------------------------------------------------!
      logical           , intent(in) :: init
      !----- Local variables --------------------------------------------------------------!
      type(edtype)      , pointer    :: cgrid
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      integer                        :: maxcohort
      integer                        :: cohort_count
      integer                        :: igr
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ibuff
      !------------------------------------------------------------------------------------!

      ! With openmp, we need to initialize as many buffers as there are threads


      if (init) then

         !---------------------------------------------------------------------------------!
         ! Initialize the photosynthesis arrays.
         !---------------------------------------------------------------------------------!
         allocate(thispft   (nthreads))
         allocate(met       (nthreads))
         allocate(aparms    (nthreads))
         allocate(stopen    (nthreads))
         allocate(stclosed  (nthreads))
         allocate(rubpsatlim(nthreads))
         allocate(thirdlim  (nthreads))
         allocate(lightlim  (nthreads))

         !---------------------------------------------------------------------------------!
         ! Initialize radiation scratch space                                              !
         !---------------------------------------------------------------------------------!
         allocate(radscr(nthreads))
         do ibuff=1,nthreads
            call nullify_radscratch(radscr(ibuff))
         end do

         !---------------------------------------------------------------------------------!
         !     If this initialization, make sure soil and sfcwater arrays are allocated.   !
         !---------------------------------------------------------------------------------!

         allocate(integration_buff(nthreads))
         allocate(rk4aux(nthreads))

         select case (integration_scheme)
         case (3)

            do ibuff=1,nthreads
               allocate(integration_buff(ibuff)%initp)
               allocate(integration_buff(ibuff)%ytemp)
               call allocate_rk4_patch(integration_buff(ibuff)%initp )
               call allocate_rk4_patch(integration_buff(ibuff)%ytemp )
            end do
         
         case default
            
            do ibuff=1,nthreads
               allocate(integration_buff(ibuff)%initp )
               allocate(integration_buff(ibuff)%yscal )
               allocate(integration_buff(ibuff)%y     )
               allocate(integration_buff(ibuff)%dydx  )
               allocate(integration_buff(ibuff)%yerr  )
               allocate(integration_buff(ibuff)%ytemp )
               call allocate_rk4_patch(integration_buff(ibuff)%initp )
               call allocate_rk4_patch(integration_buff(ibuff)%yscal )
               call allocate_rk4_patch(integration_buff(ibuff)%y     )
               call allocate_rk4_patch(integration_buff(ibuff)%dydx  )
               call allocate_rk4_patch(integration_buff(ibuff)%yerr  )
               call allocate_rk4_patch(integration_buff(ibuff)%ytemp )
            end do

         end select


         !---------------------------------------------------------------------------------!
         !     The following structures are allocated/deallocated depending on the         !
         ! integration method.                                                             !
         !---------------------------------------------------------------------------------!
         select case(integration_scheme) 
         case (0) !----- Euler. -----------------------------------------------------------!

            do ibuff=1,nthreads
               allocate(integration_buff(ibuff)%dinitp)
               call allocate_rk4_patch(integration_buff(ibuff)%dinitp)
            end do

         case (1) !----- Runge-Kutta. -----------------------------------------------------!

            do ibuff=1,nthreads
               allocate(integration_buff(ibuff)%ak2)
               allocate(integration_buff(ibuff)%ak3)
               allocate(integration_buff(ibuff)%ak4)
               allocate(integration_buff(ibuff)%ak5)
               allocate(integration_buff(ibuff)%ak6)
               allocate(integration_buff(ibuff)%ak7)
               call allocate_rk4_patch(integration_buff(ibuff)%ak2)
               call allocate_rk4_patch(integration_buff(ibuff)%ak3)
               call allocate_rk4_patch(integration_buff(ibuff)%ak4)
               call allocate_rk4_patch(integration_buff(ibuff)%ak5)
               call allocate_rk4_patch(integration_buff(ibuff)%ak6)
               call allocate_rk4_patch(integration_buff(ibuff)%ak7)
            end do
            
         case (2) !----- Heun's. ----------------------------------------------------------!

            do ibuff=1,nthreads
               allocate(integration_buff(ibuff)%ak2)
               allocate(integration_buff(ibuff)%ak3)
               call allocate_rk4_patch(integration_buff(ibuff)%ak2)
               call allocate_rk4_patch(integration_buff(ibuff)%ak3)
            end do

         case (3) !----- Hybrid (forward Euler/BDF2)---------------------------------------!
            
            do ibuff=1,nthreads
               allocate(integration_buff(ibuff)%dinitp)
               call allocate_rk4_patch(integration_buff(ibuff)%dinitp)
               allocate(integration_buff(ibuff)%yprev)
            end do
               
         end select
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    If this is not initialization, deallocate cohort memory from integration     !
         ! patches.                                                                        !
         !---------------------------------------------------------------------------------!

         select case (integration_scheme)
         case (3)

            do ibuff=1,nthreads
               call deallocate_rk4_coh(integration_buff(ibuff)%initp )
               call deallocate_rk4_coh(integration_buff(ibuff)%ytemp )
            end do
         case default
            do ibuff=1,nthreads
               call deallocate_rk4_coh(integration_buff(ibuff)%initp )
               call deallocate_rk4_coh(integration_buff(ibuff)%yscal )
               call deallocate_rk4_coh(integration_buff(ibuff)%y     )
               call deallocate_rk4_coh(integration_buff(ibuff)%dydx  )
               call deallocate_rk4_coh(integration_buff(ibuff)%yerr  )
               call deallocate_rk4_coh(integration_buff(ibuff)%ytemp )
            end do
         end select

         !------ De-allocate the auxiliary structure. -------------------------------------!
         do ibuff=1,nthreads
            call deallocate_rk4_aux(rk4aux(ibuff))
         end do

         !---------------------------------------------------------------------------------!
         !     The following structures are allocated/deallocated depending on the         !
         ! integration method.                                                             !
         !---------------------------------------------------------------------------------!
         select case(integration_scheme) 
         case (0) !----- Euler. -----------------------------------------------------------!
            do ibuff=1,nthreads
               call deallocate_rk4_coh(integration_buff(ibuff)%dinitp)
            end do
         case (1) !----- Runge-Kutta. -----------------------------------------------------!
            do ibuff=1,nthreads
               call deallocate_rk4_coh(integration_buff(ibuff)%ak2)
               call deallocate_rk4_coh(integration_buff(ibuff)%ak3)
               call deallocate_rk4_coh(integration_buff(ibuff)%ak4)
               call deallocate_rk4_coh(integration_buff(ibuff)%ak5)
               call deallocate_rk4_coh(integration_buff(ibuff)%ak6)
               call deallocate_rk4_coh(integration_buff(ibuff)%ak7)
            end do
         case (2) !----- Heun's. ----------------------------------------------------------!
            do ibuff=1,nthreads
               call deallocate_rk4_coh(integration_buff(ibuff)%ak2)
               call deallocate_rk4_coh(integration_buff(ibuff)%ak3)
            end do
         case (3) !----- Hybrid -----------------------------------------------------------!
            do ibuff=1,nthreads
               call deallocate_rk4_coh(integration_buff(ibuff)%dinitp)
               call deallocate_bdf2_patch(integration_buff(ibuff)%yprev)
            end do
         end select

         !---------------------------------------------------------------------------------!
      end if

      !----- Find maximum number of cohorts amongst all patches ---------------------------!
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

      !----- Create new memory in each of the integration patches. ------------------------!
      select case (integration_scheme)
      case (3)
         do ibuff=1,nthreads
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%initp )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ytemp )
         end do
      case default
         do ibuff=1,nthreads
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%initp )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%yscal )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%y     )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%dydx  )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%yerr  )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ytemp )
         end do
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Allocate/deallocate structures depending on the integration method.            !
      !------------------------------------------------------------------------------------!
      select case(integration_scheme) 
      case (0) !----- Euler. --------------------------------------------------------------!
         do ibuff=1,nthreads
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%dinitp)
         end do
      case (1) !----- Runge-Kutta. --------------------------------------------------------!
         do ibuff=1,nthreads
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak2   )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak3   )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak4   )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak5   )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak6   )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak7   )
         end do
      case (2) !----- Heun's. -------------------------------------------------------------!
         do ibuff=1,nthreads
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak2   )
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%ak3   )
         end do
      case (3) !----- Hybrid --------------------------------------------------------------!
         do ibuff=1,nthreads
            call allocate_rk4_coh(maxcohort,integration_buff(ibuff)%dinitp)
            call allocate_bdf2_patch(integration_buff(ibuff)%yprev,maxcohort)
         end do
      end select
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! Initialize radiation scratch space                                                 !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
      case (1)
         !---- Big leaf.  Use the maximum LAI. --------------------------------------------!
         maxcohort = 0
         do igr = 1,ngrids
            cgrid => edgrid_g(igr)
            do ipy = 1,cgrid%npolygons
               cpoly => cgrid%polygon(ipy)
               do isi = 1,cpoly%nsites
                  csite => cpoly%site(isi)
                  do ipa = 1,csite%npatches
                     cpatch => csite%patch(ipa)
                     cohort_count = ceiling( (cpatch%lai(1) + cpatch%wai(1)) / tai_lyr_max )
                     maxcohort = max(maxcohort, cohort_count)
                  end do
               end do
            end do
         end do
      end select
    
      do ibuff=1,nthreads
         call dealloc_radscratch(radscr(ibuff))
      end do
      do ibuff=1,nthreads
         call alloc_radscratch(radscr(ibuff),maxcohort)
      end do

      !------ Allocate and initialise the auxiliary structure. ----------------------------!
      do ibuff=1,nthreads
         call allocate_rk4_aux(rk4aux(ibuff),nzg,nzs,maxcohort)
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine initialize_rk4patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine initialize the multiple time step-related variables.            !
   !                                                                                       !
   !  MLO.  RGK, I changed a couple of things in this sub-routine.  I didn't see any       !
   !        reason to define the runoff and time variables inside the select case (all     !
   !        schemes can access variables from rk4_coms).  Also, I moved the detailed       !
   !        output outside the case selection, but check whether it should print detailed  !
   !        output.  You are right about it not working for multiple polygons, and         !
   !        ed_opspec.F90 already checks this.  It would overwrite in the case of multiple !
   !        sites per run, but I fixed this.                                               !
   !---------------------------------------------------------------------------------------!
   subroutine initialize_misc_stepvars()
      use rk4_coms      , only  : tbeg                & ! intent(inout)
                                , tend                & ! intent(inout)
                                , dtrk4               & ! intent(inout)
                                , dtrk4i              & ! intent(inout)
                                , print_detailed      & ! intent(in)
                                , detail_pref         & ! intent(in)
                                , print_thbnd         & ! intent(in)
                                , thbnds_fout         ! ! intent(in)
      use ed_misc_coms   , only : dtlsm               ! ! intent(in)
      use ed_max_dims    , only : str_len             ! ! intent(in)
      use soil_coms      , only : runoff_time         & ! intent(in)
                                , runoff_time_i       & ! intent(out)
                                , simplerunoff        ! ! intent(out)
      use hydrology_coms , only : useRUNOFF           ! ! intent(in)
      use ed_state_vars,   only : edgrid_g            ! ! structure
      use grid_coms,       only : ngrids              ! ! intent(in)

      implicit none
      !------ Local variables. ------------------------------------------------------------!
      integer                :: igr
      integer                :: ipy
      integer                :: isi
      integer                :: ipa
      integer                :: ico
      character(len=str_len) :: detail_fout
      logical                :: isthere
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Check whether printing detailed output.  Despite the loop, this type of output !
      ! works for single polygon runs.                                                     !
      !------------------------------------------------------------------------------------!
      if (print_detailed) then
         do igr=1,ngrids
            do ipy=1,edgrid_g(igr)%npolygons
               do isi=1,edgrid_g(igr)%polygon(ipy)%nsites
                  do ipa = 1,edgrid_g(igr)%polygon(ipy)%site(isi)%npatches
                     !---------------------------------------------------------------------!
                     ! Patch level files.                                                  !
                     !---------------------------------------------------------------------!
                     write (detail_fout,fmt='(2a,2(i4.4,a))')                              &
                           trim(detail_pref),'prk4_site_',isi,'_patch_',ipa,'.txt'
                     
                     inquire(file=trim(detail_fout),exist=isthere)
                     if (isthere) then
                        !---- Open the file to delete when closing. -----------------------!
                        open (unit=83,file=trim(detail_fout),status='old',action='write')
                        close(unit=83,status='delete')
                     end if
                     !---------------------------------------------------------------------!
                     !---------------------------------------------------------------------!
                     ! Cohort level files.                                                 !
                     !---------------------------------------------------------------------!
                     do ico = 1,edgrid_g(igr)%polygon(ipy)%site(isi)%patch(ipa)%ncohorts
                        write (detail_fout,fmt='(2a,3(i4.4,a))')                           &
                           trim(detail_pref),'crk4_site_',isi,'_patch_',ipa,'_cohort_',ico &
                                            ,'.txt'
                        inquire(file=trim(detail_fout),exist=isthere)
                        if (isthere) then
                           !---- Open the file to delete when closing. --------------------!
                           open (unit=84,file=trim(detail_fout),status='old',action='write')
                           close(unit=84,status='delete')
                        end if
                     end do
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Write header for thermodynamic boundaries in case this output is active.       !
      !------------------------------------------------------------------------------------!
      if (print_thbnd) then
         open (unit=39,file=trim(thbnds_fout),status='replace',action='write')
         write(unit=39,fmt='(18(a,1x))')  '        YEAR','       MONTH','         DAY'     &
                                         ,'        HOUR','        MINU','        SECO'     &
                                         ,'    MIN_TEMP','    MAX_TEMP','     MIN_SHV'     &
                                         ,'     MAX_SHV','     MIN_CO2','     MAX_CO2'     &
                                         ,'   MIN_THETA','   MAX_THETA','    MIN_PRSS'     &
                                         ,'    MAX_PRSS','MIN_ENTHALPY','MAX_ENTHALPY'
         close(unit=39,status='keep')
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Check whether we will use runoff or not, and saving this check to save time.   !
      !------------------------------------------------------------------------------------!
      simplerunoff = useRUNOFF == 0 .and. runoff_time /= 0.
      if (runoff_time /= 0.) then
         runoff_time_i = 1.0 / runoff_time
      else 
         runoff_time_i = 0.0
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Set time-step related variables, which shall remain constant throughout the     !
      ! entire simulation.                                                                 !
      !------------------------------------------------------------------------------------!
      tbeg   = 0.d0
      tend   = dble(dtlsm)
      dtrk4  = tend - tbeg
      dtrk4i = 1.d0/dtrk4
      !------------------------------------------------------------------------------------!

      return
   end subroutine initialize_misc_stepvars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine is the main Runge-Kutta step driver.                                !
   !---------------------------------------------------------------------------------------!
   subroutine rkqs(x,htry,hgoal,hdid,hnext,csite,ipa,isi,ibuff)

      use rk4_coms      , only : rk4patchtype              & ! structure
                               , integration_buff          & ! intent(inout)
                               , rk4site                   & ! intent(in)
                               , hmin                      & ! intent(in)
                               , rk4eps                    & ! intent(in)
                               , rk4epsi                   & ! intent(in)
                               , safety                    & ! intent(in)
                               , pgrow                     & ! intent(in)
                               , pshrnk                    & ! intent(in)
                               , print_detailed            & ! intent(in)
                               , norm_rk4_fluxes           & ! intent(in)
                               , reset_rk4_fluxes          ! ! intent(in)
      use rk4_copy_patch, only : copy_rk4_patch            ! ! sub-routine
      use rk4_misc      , only : print_errmax              & ! sub-routine
                               , print_rk4patch            & ! sub-routine
                               , adjust_veg_properties     & ! sub-routine
                               , adjust_topsoil_properties & ! sub-routine
                               , adjust_sfcw_properties    & ! sub-routine
                               , update_diagnostic_vars    & ! sub-routine
                               , update_density_vars       & ! sub-routine
                               , print_rk4_state           ! ! sub-routine
      use ed_state_vars , only : sitetype                  & ! structure
                               , patchtype                 ! ! structure
      use grid_coms     , only : nzg                       & ! intent(in)
                               , nzs                       ! ! intent(in)
      use ed_misc_coms  , only : dtlsm                     ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer                  , intent(in)    :: ipa
      integer                  , intent(in)    :: isi
      integer                  , intent(in)    :: ibuff
      type(sitetype)           , target        :: csite
      real(kind=8)             , intent(in)    :: htry
      real(kind=8)             , intent(inout) :: hgoal
      real(kind=8)             , intent(inout) :: x
      real(kind=8)             , intent(out)   :: hdid
      real(kind=8)             , intent(out)   :: hnext
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                             :: h,errmax,xnew,newh,oldh
      real(kind=8)                             :: fgrow
      logical                                  :: reject_step,reject_result
      logical                                  :: minstep,stuck,test_reject
      logical                                  :: gapstep
      integer                                  :: k
      !------------------------------------------------------------------------------------!

      gapstep     =  htry < hgoal
      h           =  htry
      reject_step =  .false.
      hstep:   do

         !---------------------------------------------------------------------------------!
         ! 1. Try a step of varying size.                                                  !
         !---------------------------------------------------------------------------------!
         call rkck(ibuff,integration_buff(ibuff)%y,integration_buff(ibuff)%dydx            &
                  ,integration_buff(ibuff)%ytemp,integration_buff(ibuff)%yerr              &
                  ,integration_buff(ibuff)%ak2,integration_buff(ibuff)%ak3                 &
                  ,integration_buff(ibuff)%ak4,integration_buff(ibuff)%ak5                 &
                  ,integration_buff(ibuff)%ak6,integration_buff(ibuff)%ak7,h,csite,ipa     &
                  ,reject_step,reject_result)

         !---------------------------------------------------------------------------------!
         ! 2. Check to see how accurate the step was.  Errors were calculated by integrat- !
         !    ing the derivative of that last step.                                        !
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
                           ,integration_buff(ibuff)%ytemp)
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
            if (minstep .or. stuck) then

               write (unit=*,fmt='(80a)')         ('=',k=1,80)
               write (unit=*,fmt='(a)')           '   STEPSIZE UNDERFLOW IN RKQS'
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
               if (reject_step .or. reject_result) then
                  write (unit=*,fmt='(a)') '   Likely to be a rejected step problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               else
                  write (unit=*,fmt='(a)') '   Likely to be an errmax problem.'
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if

               if (reject_result) then
                  !----- Run the LSM sanity check but this time we force the print. -------!
                  call rk4_sanity_check(ibuff,integration_buff(ibuff)%ytemp,test_reject    &
                                       ,csite,ipa,integration_buff(ibuff)%dydx,h,.true.)
                  call print_sanity_check(integration_buff(ibuff)%y,csite,ipa)
               elseif (reject_step) then
                  call rk4_sanity_check(ibuff,integration_buff(ibuff)%ak7,test_reject      &
                                       ,csite,ipa,integration_buff(ibuff)%dydx,h,.true.)
                  call print_sanity_check(integration_buff(ibuff)%y,csite,ipa)
               else
                  call print_errmax(errmax,integration_buff(ibuff)%yerr                    &
                                   ,integration_buff(ibuff)%yscal,csite%patch(ipa)         &
                                   ,integration_buff(ibuff)%y)
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Rel. errmax:',errmax
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Raw errmax: ',errmax*rk4eps
                  write (unit=*,fmt='(a,1x,es12.4)') ' - Epsilon:',rk4eps
                  write (unit=*,fmt='(80a)') ('=',k=1,80)
               end if
               call print_rk4patch(integration_buff(ibuff)%y, csite,ipa)
            endif
         
         else

            !------------------------------------------------------------------------------!
            ! 3b.  Great, it worked, so now we can advance to the next step.  We just need !
            !      to do some minor adjustments before...                                  !
            !------------------------------------------------------------------------------!
            !----- i.   Final update of leaf properties to avoid negative water. ----------!
            call adjust_veg_properties(integration_buff(ibuff)%ytemp,h,csite,ipa,ibuff)
            !----- ii.  Final update of top soil properties to avoid off-bounds moisture. -!
            call adjust_topsoil_properties(integration_buff(ibuff)%ytemp,h)
            !----- iii. Make temporary surface water stable and positively defined. -------!
            call adjust_sfcw_properties(nzg,nzs,integration_buff(ibuff)%ytemp,h,csite,ipa)
            !----- iv.  Update the diagnostic variables. ----------------------------------!
            call update_diagnostic_vars(integration_buff(ibuff)%ytemp,csite,ipa,ibuff)
            !----- v.  Update density. ----------------------------------------------------!
            call update_density_vars(integration_buff(ibuff)%ytemp                         &
                                    ,integration_buff(ibuff)%y    )
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            ! 3c. Set up h for the next time.  And here we can relax h for the next step,  !
            !    and try something faster, unless this is a "gap step" (shorter time step  !
            !    just to close the full thermodynamic time step).                          !
            !------------------------------------------------------------------------------!
            if (gapstep) then
               hnext = hgoal
            else
               fgrow = min(5.d0,max(safety*errmax**pgrow,1.d0))
               hnext = max(2.d0*hmin, min(dble(dtlsm), fgrow * hgoal))
            end if
            !------------------------------------------------------------------------------!


            !------ 3d. Normalise the fluxes if the user wants detailed debugging. --------!
            if (print_detailed) then
               call norm_rk4_fluxes(integration_buff(ibuff)%ytemp,h)
               call print_rk4_state(integration_buff(ibuff)%y                              &
                                   ,integration_buff(ibuff)%ytemp,csite,ipa,isi,x,h)
            end if
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    3e. Copy the temporary structure to the intermediate state.               !
            !------------------------------------------------------------------------------!
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
            x    = x + h
            hdid = h

            exit hstep
         end if
      end do hstep

      return
   end subroutine rkqs
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will update the variables and perform the actual time stepping.    !
   !---------------------------------------------------------------------------------------!
   subroutine rkck(ibuff,y,dydx,yout,yerr,ak2,ak3,ak4,ak5,ak6,ak7,h,csite,ipa              &
                  ,reject_step,reject_result)

      use rk4_coms       , only : rk4patchtype           & ! structure
                                , integration_vars       & ! structure
                                , print_diags            & ! intent(in)
                                , rk4_a2                 & ! intent(in)
                                , rk4_a3                 & ! intent(in)
                                , rk4_a4                 & ! intent(in)
                                , rk4_a5                 & ! intent(in)
                                , rk4_a6                 & ! intent(in)
                                , rk4_b21                & ! intent(in)
                                , rk4_b31                & ! intent(in)
                                , rk4_b32                & ! intent(in)
                                , rk4_b41                & ! intent(in)
                                , rk4_b42                & ! intent(in)
                                , rk4_b43                & ! intent(in)
                                , rk4_b51                & ! intent(in)
                                , rk4_b52                & ! intent(in)
                                , rk4_b53                & ! intent(in)
                                , rk4_b54                & ! intent(in)
                                , rk4_b61                & ! intent(in)
                                , rk4_b62                & ! intent(in)
                                , rk4_b63                & ! intent(in)
                                , rk4_b64                & ! intent(in)
                                , rk4_b65                & ! intent(in)
                                , rk4_c1                 & ! intent(in)
                                , rk4_c3                 & ! intent(in)
                                , rk4_c4                 & ! intent(in)
                                , rk4_c6                 & ! intent(in)
                                , rk4_dc5                & ! intent(in)
                                , rk4_dc1                & ! intent(in)
                                , rk4_dc3                & ! intent(in)
                                , rk4_dc4                & ! intent(in)
                                , rk4_dc6                & ! intent(in)
                                , zero_rk4_patch         & ! intent(in)
                                , zero_rk4_cohort        ! ! intent(in)
      use rk4_copy_patch , only : copy_rk4_patch         ! ! sub-routine
      use rk4_derivs     , only : leaf_derivs            ! ! sub-routine
      use rk4_misc       , only : update_diagnostic_vars & ! sub-routine
                                , print_rk4patch         ! ! sub-routine
      use ed_state_vars  , only : sitetype               & ! structure
                                , patchtype              ! ! structure
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      integer           , intent(in)  :: ibuff
      integer           , intent(in)  :: ipa
      real(kind=8)      , intent(in)  :: h
      type(rk4patchtype), target      :: y,dydx,yout,yerr
      type(rk4patchtype), target      :: ak2,ak3,ak4,ak5,ak6,ak7
      type(sitetype)    , target      :: csite
      logical           , intent(out) :: reject_step
      logical           , intent(out) :: reject_result
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)   , pointer     :: cpatch
      real(kind=8)                    :: combh
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



      !------------------------------------------------------------------------------------!
      !     For each stage we checked the sanity, so we avoid using non-sense values to    !
      ! advance.  Also, we estimate the derivative of pressure after each stage, and in    !
      ! the end we estimate the full step derivative as the weighted average for each of   !
      ! these partial steps taken.                                                         !
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Second stage (the first was the Euler, whose derivative is already computed    !
      ! and saved in dydx.                                                                 !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b21*h, cpatch)
      combh = rk4_b21*h
      call update_diagnostic_vars(ak7, csite,ipa,ibuff)
      call rk4_sanity_check(ibuff,ak7, reject_step, csite, ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak2, csite, ipa,ibuff,h,.false.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Third stage.                                                                    !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b31*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b32*h, cpatch)
      combh = (rk4_b31+rk4_b32)*h
      call update_diagnostic_vars(ak7, csite,ipa,ibuff)
      call rk4_sanity_check(ibuff,ak7,reject_step,csite,ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak3, csite,ipa,ibuff,h,.false.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fourth stage.                                                                   !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b41*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b42*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, rk4_b43*h, cpatch)
      combh = (rk4_b41+rk4_b42+rk4_b43)*h
      call update_diagnostic_vars(ak7, csite,ipa,ibuff)
      call rk4_sanity_check(ibuff,ak7, reject_step, csite,ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak4, csite, ipa,ibuff,h,.false.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fifth stage.                                                                    !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b51*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b52*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, rk4_b53*h, cpatch)
      call inc_rk4_patch(ak7,  ak4, rk4_b54*h, cpatch)
      combh = (rk4_b51+rk4_b52+rk4_b53+rk4_b54)*h
      call update_diagnostic_vars(ak7, csite,ipa,ibuff)
      call rk4_sanity_check(ibuff,ak7,reject_step,csite,ipa,dydx,h,print_diags)
      if (reject_step) return
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak5, csite, ipa,ibuff,h,.false.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Sixth stage.                                                                    !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, ak7, cpatch)
      call inc_rk4_patch(ak7, dydx, rk4_b61*h, cpatch)
      call inc_rk4_patch(ak7,  ak2, rk4_b62*h, cpatch)
      call inc_rk4_patch(ak7,  ak3, rk4_b63*h, cpatch)
      call inc_rk4_patch(ak7,  ak4, rk4_b64*h, cpatch)
      call inc_rk4_patch(ak7,  ak5, rk4_b65*h, cpatch)
      combh = (rk4_b61+rk4_b62+rk4_b63+rk4_b64+rk4_b65)*h
      call update_diagnostic_vars(ak7, csite,ipa,ibuff)
      call rk4_sanity_check(ibuff,ak7, reject_step, csite,ipa,dydx,h,print_diags)
      if(reject_step)return
      !------------------------------------------------------------------------------------!


      !------ Get the new derivative evaluation. ------------------------------------------!
      call leaf_derivs(ak7, ak6, csite,ipa,ibuff,h,.false.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Aggregate all derivatives to make the new guess.                                !
      !------------------------------------------------------------------------------------!
      call copy_rk4_patch(y, yout, cpatch)
      call inc_rk4_patch(yout, dydx, rk4_c1*h, cpatch)
      call inc_rk4_patch(yout,  ak3, rk4_c3*h, cpatch)
      call inc_rk4_patch(yout,  ak4, rk4_c4*h, cpatch)
      call inc_rk4_patch(yout,  ak6, rk4_c6*h, cpatch)
      combh = (rk4_c1+rk4_c3+rk4_c4+rk4_c6)*h
      call update_diagnostic_vars   (yout, csite,ipa,ibuff)
      call rk4_sanity_check(ibuff,yout, reject_result, csite,ipa,dydx,h,print_diags)
      !------------------------------------------------------------------------------------!
      if(reject_result)return
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Estimate the error for this step.                                               !
      !------------------------------------------------------------------------------------!
      call zero_rk4_patch (yerr)
      call zero_rk4_cohort(yerr)
      call inc_rk4_patch(yerr, dydx, rk4_dc1*h, cpatch)
      call inc_rk4_patch(yerr, ak3,  rk4_dc3*h, cpatch)
      call inc_rk4_patch(yerr, ak4,  rk4_dc4*h, cpatch)
      call inc_rk4_patch(yerr, ak5,  rk4_dc5*h, cpatch)
      call inc_rk4_patch(yerr, ak6,  rk4_dc6*h, cpatch)
      !------------------------------------------------------------------------------------!

      return
   end subroutine rkck
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will check for potentially serious problems.  Note that the upper  !
   ! and lower bound are defined in rk4_coms.f90, so if you need to change any limit for   !
   ! some reason, you can adjust there.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_sanity_check(ibuff,y,reject_step, csite,ipa,dydx,h,print_problems)
      use rk4_coms              , only : rk4patchtype          & ! structure
                                       , integration_vars      & ! structure
                                       , rk4site               & ! structure
                                       , rk4aux                & ! structure
                                       , rk4eps                & ! intent(in)
                                       , rk4max_can_shv        & ! intent(in)
                                       , rk4min_can_shv        & ! intent(in)
                                       , rk4min_can_rhv        & ! intent(in)
                                       , rk4max_can_rhv        & ! intent(in)
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

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer            , intent(in)  :: ibuff
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
      real(kind=8)                     :: rk4min_leaf_water_im2
      real(kind=8)                     :: rk4max_leaf_water_im2
      real(kind=8)                     :: rk4min_wood_water_im2
      real(kind=8)                     :: rk4max_wood_water_im2
      integer                          :: ipa
      integer                          :: ico
      integer                          :: ipft
      logical                          :: cflag6
      logical                          :: cflag7
      logical                          :: cflag8
      logical                          :: cflag9
      logical                          :: cflag10
      logical                          :: cflag11
      !------------------------------------------------------------------------------------!


      !----- Be optimistic and start assuming that things are fine. -----------------------!
      reject_step = .false.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   Check whether the canopy air specific enthalpy is off.                           !
      !------------------------------------------------------------------------------------!
      if (y%can_enthalpy > rk4aux(ibuff)%rk4max_can_enthalpy   .or.                        &
          y%can_enthalpy < rk4aux(ibuff)%rk4min_can_enthalpy        ) then
         reject_step = .true.
         if(record_err) integ_err(1,2) = integ_err(1,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + CAS specific enthalpy is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_ENTHALPY:      ',y%can_enthalpy
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DMOL:          ',y%can_dmol
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_ENTHALPY)/Dt:',dydx%can_enthalpy
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
      !   Check whether the canopy air potential temperature is off.                       !
      !------------------------------------------------------------------------------------!
      if (  y%can_theta > rk4aux(ibuff)%rk4max_can_theta .or.                              &
            y%can_theta < rk4aux(ibuff)%rk4min_can_theta ) then
         reject_step = .true.
         if(record_err) integ_err(2,2) = integ_err(2,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + CAS potential temp. is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_ENTHALPY:      ',y%can_enthalpy
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DMOL:          ',y%can_dmol
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_ENTHALPY)/Dt:',dydx%can_enthalpy
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
      !   Check whether the canopy air equivalent potential temperature is off.            !
      !------------------------------------------------------------------------------------!
      if ( y%can_shv > rk4max_can_shv  .or. y%can_shv < rk4min_can_shv ) then
         reject_step = .true.
         if(record_err) integ_err(3,2) = integ_err(3,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + CAS specific humidity is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_ENTHALPY:      ',y%can_enthalpy
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DMOL:          ',y%can_dmol
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_ENTHALPY)/Dt:',dydx%can_enthalpy
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
            write(unit=*,fmt='(a)')           ' + CAS temperature is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_ENTHALPY:      ',y%can_enthalpy
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DMOL:          ',y%can_dmol
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_ENTHALPY)/Dt:',dydx%can_enthalpy
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
      !   Check whether the canopy air CO2 molar count is off.                             !
      !------------------------------------------------------------------------------------!
      if ( y%can_co2 > rk4max_can_co2              .or.                                    &
           y%can_co2 < rk4min_can_co2                   ) then
         reject_step = .true.
         if(record_err) integ_err(5,2) = integ_err(5,2) + 1_8
         if (print_problems) then
            write(unit=*,fmt='(a)')           '==========================================='
            write(unit=*,fmt='(a)')           ' + CAS CO2 mixing ratio is off-track...'
            write(unit=*,fmt='(a)')           '-------------------------------------------'
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_ENTHALPY:      ',y%can_enthalpy
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DMOL:          ',y%can_dmol
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
            write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
            write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_ENTHALPY)/Dt:',dydx%can_enthalpy
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
      !     Check leaf properties, but only for those cohorts with sufficient LAI.         !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cflag6  = .false.
      cflag7  = .false.
      cflag8  = .false.
      leafloop: do ico = 1,cpatch%ncohorts
         if (.not. y%leaf_resolvable(ico)) cycle leafloop

         !----- Find the boundaries for leaf water (surface and internal). ----------------!
         rk4min_leaf_water     = rk4min_veg_lwater * y%lai(ico)
         rk4min_leaf_water_im2 = rk4aux(ibuff)%rk4min_leaf_water_im2(ico) * (1.d0 - rk4eps)
         rk4max_leaf_water_im2 = rk4aux(ibuff)%rk4max_leaf_water_im2(ico) * (1.d0 + rk4eps)
         !---------------------------------------------------------------------------------!



         !----- Check leaf surface water. -------------------------------------------------!
         if (y%leaf_water(ico) < rk4min_leaf_water) then
            reject_step = .true.
            if(record_err) cflag6 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Leaf surface water is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' ICO:           ',ico
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
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER_IM2:',y%leaf_water_im2(ico)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_IM2)/Dt:',dydx%leaf_water_im2(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
         !---------------------------------------------------------------------------------!


         !----- Check leaf temperature. ---------------------------------------------------!
         if (y%leaf_temp(ico) > rk4max_veg_temp .or.                                       &
             y%leaf_temp(ico) < rk4min_veg_temp      ) then
            reject_step = .true.
            if(record_err) cflag7 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Leaf temperature is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' ICO:           ',ico
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
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER_IM2:',y%leaf_water_im2(ico)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_IM2)/Dt:',dydx%leaf_water_im2(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
         !---------------------------------------------------------------------------------!



         !----- Check leaf internal water. ------------------------------------------------!
         if (y%leaf_water_im2(ico) > rk4max_leaf_water_im2 .or.                            &
             y%leaf_water_im2(ico) < rk4min_leaf_water_im2      ) then
            reject_step = .true.
            if(record_err) cflag8 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Leaf internal water is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' ICO:           ',ico
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
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_WATER_IM2:',y%leaf_water_im2(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' MIN_LEAFW_IM2: ',rk4min_leaf_water_im2
               write(unit=*,fmt='(a,1x,es12.4)') ' MAX_LEAFW_IM2: ',rk4max_leaf_water_im2
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBH:      ',y%leaf_gbh(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GBW:      ',y%leaf_gbw(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_REYNOLDS: ',y%leaf_reynolds(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_GRASHOF:  ',y%leaf_grashof(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFREE:   ',y%leaf_nussfree(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LEAF_NUFORC:   ',y%leaf_nussforc(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_EN)/Dt: ',dydx%leaf_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_WAT)/Dt:',dydx%leaf_water(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(LEAF_IM2)/Dt:',dydx%leaf_water_im2(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
         !---------------------------------------------------------------------------------!
      end do leafloop
      if(record_err .and. cflag6) integ_err( 6,2) = integ_err( 6,2) + 1_8
      if(record_err .and. cflag7) integ_err( 7,2) = integ_err( 7,2) + 1_8
      if(record_err .and. cflag8) integ_err( 8,2) = integ_err( 8,2) + 1_8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check wood properties, but only for those cohorts with sufficient LAI.         !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cflag9  = .false.
      cflag10 = .false.
      cflag11 = .false.
      woodloop: do ico = 1,cpatch%ncohorts
         if (.not. y%wood_resolvable(ico)) cycle woodloop

         !----- Find the boundaries for wood water (surface and internal). ----------------!
         rk4min_wood_water = rk4min_veg_lwater * y%wai(ico)
         rk4min_wood_water_im2 = rk4aux(ibuff)%rk4min_wood_water_im2(ico) * (1.d0 - rk4eps)
         rk4max_wood_water_im2 = rk4aux(ibuff)%rk4max_wood_water_im2(ico) * (1.d0 + rk4eps)
         !---------------------------------------------------------------------------------!


         !----- Check wood surface water. -------------------------------------------------!
         if (y%wood_water(ico) < rk4min_wood_water) then
            reject_step = .true.
            if(record_err) cflag9 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Wood surface water is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' ICO:           ',ico
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
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER_IM2:',y%wood_water_im2(ico)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_IM2)/Dt:',dydx%wood_water_im2(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
         !---------------------------------------------------------------------------------!



         !----- Check wood temperature. ---------------------------------------------------!
         if (y%wood_temp(ico) > rk4max_veg_temp .or.                                       &
             y%wood_temp(ico) < rk4min_veg_temp      ) then
            reject_step = .true.
            if(record_err) cflag10 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Wood temperature is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' ICO:           ',ico
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
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER_IM2:',y%wood_water_im2(ico)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_IM2)/Dt:',dydx%wood_water_im2(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
         !---------------------------------------------------------------------------------!


         !----- Check wood internal water. ------------------------------------------------!
         if (y%wood_water_im2(ico) > rk4max_wood_water_im2 .or.                            &
             y%wood_water_im2(ico) < rk4min_wood_water_im2      ) then
            reject_step = .true.
            if(record_err) cflag11 = .true.
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Wood internal water is off-track...'
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a,1x,i6)')     ' ICO:           ',ico
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
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_WATER_IM2:',y%wood_water_im2(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' VEG_WIND:      ',y%veg_wind(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' LINT_SHV:      ',y%lint_shv(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' MIN_WOODW_IM2: ',rk4min_wood_water_im2
               write(unit=*,fmt='(a,1x,es12.4)') ' MAX_WOODW_IM2: ',rk4max_wood_water_im2
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBH:      ',y%wood_gbh(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GBW:      ',y%wood_gbw(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_REYNOLDS: ',y%wood_reynolds(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_GRASHOF:  ',y%wood_grashof(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFREE:   ',y%wood_nussfree(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' WOOD_NUFORC:   ',y%wood_nussforc(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_EN)/Dt: ',dydx%wood_energy(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_WAT)/Dt:',dydx%wood_water(ico)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(WOOD_IM2)/Dt:',dydx%wood_water_im2(ico)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
         !---------------------------------------------------------------------------------!
      end do woodloop
      if(record_err .and. cflag9 ) integ_err( 9,2) = integ_err( 9,2) + 1_8
      if(record_err .and. cflag10) integ_err(10,2) = integ_err(10,2) + 1_8
      if(record_err .and. cflag11) integ_err(11,2) = integ_err(11,2) + 1_8
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Check the water mass of the virtual pool.  The energy is checked only when     !
      ! there is enough mass.                                                              !
      !------------------------------------------------------------------------------------!
      if (y%virtual_water < rk4min_virt_water) then
         reject_step = .true.
         if(record_err) integ_err(13,2) = integ_err(13,2) + 1_8
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
         if(record_err) integ_err(12,2) = integ_err(12,2) + 1_8
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
         if (y%soil_water(k)< rk4aux(ibuff)%rk4min_soil_water(k) .or.                                    &
             y%soil_water(k)> rk4aux(ibuff)%rk4max_soil_water(k) ) then
            reject_step = .true.
            if(record_err) integ_err(osow+k,2) = integ_err(osow+k,2) + 1_8
            if (print_problems) then
               write(unit=*,fmt='(a)')           '========================================'
               write(unit=*,fmt='(a)')           ' + Soil layer water is off-track...'
               write(unit=*,fmt='(a)')           '----------------------------------------'
               write(unit=*,fmt='(a,1x,i6)')     ' Level:       ',k
               write(unit=*,fmt='(a,1x,es12.4)') ' H:           ',h
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_MSTPOT: ',y%soil_mstpot(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_E)/Dt:',dydx%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_M)/Dt:',dydx%soil_water(k)
               if (k == nzg .and. y%nlev_sfcwater > 0) then
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(1)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' H:           ',h
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_TEMPK:  ',y%soil_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_FLIQ :  ',y%soil_fracliq(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_ENERGY: ',y%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_WATER:  ',y%soil_water(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SOIL_MSTPOT: ',y%soil_mstpot(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_E)/Dt:',dydx%soil_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SOIL_M)/Dt:',dydx%soil_water(k)
               if (k == nzg .and. y%nlev_sfcwater > 0) then
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:   ',y%sfcwater_tempk(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY: ',y%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:   ',y%sfcwater_mass(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:  ',y%sfcwater_depth(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:',dydx%sfcwater_energy(1)
                  write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:',dydx%sfcwater_mass(1)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:     ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY:   ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:     ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:    ',y%sfcwater_depth(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:  ',dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:  ',dydx%sfcwater_mass(k)
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
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_TEMP:     ',y%sfcwater_tempk(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_ENERGY:   ',y%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_MASS:     ',y%sfcwater_mass(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' SFCW_DEPTH:    ',y%sfcwater_depth(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_E)/Dt:  ',dydx%sfcwater_energy(k)
               write(unit=*,fmt='(a,1x,es12.4)') ' D(SFCW_M)/Dt:  ',dydx%sfcwater_mass(k)
               write(unit=*,fmt='(a)')           '========================================'
            elseif (.not. record_err) then
               return
            end if
         end if
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Print bounds to help debugging.                                                !
      !------------------------------------------------------------------------------------!
      if (reject_step .and. print_problems) then
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('=',k=1,78)
         write(unit=*,fmt='(a,1x,es12.4)') ' TIMESTEP:          ',h
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           '         ---- SANITY CHECK BOUNDS ----'
         write(unit=*,fmt='(a)')           ' '

         !----- 1. Canopy air space. ------------------------------------------------------!
         write(unit=*,fmt='(a)')           ' 1. CANOPY AIR SPACE: '
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(4(a,1x))')     '     MIN_SHV','     MAX_SHV','     MIN_RHV'    &
                                          ,'     MAX_RHV'
         write(unit=*,fmt='(4(es12.5,1x))')  rk4min_can_shv  ,rk4max_can_shv               &
                                            ,rk4min_can_rhv  ,rk4max_can_rhv
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(6(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_THETA'    &
                                          ,'   MAX_THETA','MIN_ENTHALPY','MAX_ENTHALPY'
         write(unit=*,fmt='(6(es12.5,1x))') rk4min_can_temp    ,rk4max_can_temp            &
                                           ,rk4aux(ibuff)%rk4min_can_theta                 &
                                           ,rk4aux(ibuff)%rk4max_can_theta                 &
                                           ,rk4aux(ibuff)%rk4min_can_enthalpy              &
                                           ,rk4aux(ibuff)%rk4max_can_enthalpy
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(4(a,1x))')     '    MIN_PRSS','    MAX_PRSS','     MIN_CO2'    &
                                          ,'     MAX_CO2'
         write(unit=*,fmt='(4(es12.5,1x))') rk4aux(ibuff)%rk4min_can_prss                  &
                                           ,rk4aux(ibuff)%rk4max_can_prss                  &
                                           ,rk4min_can_co2  ,rk4max_can_co2
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         !---------------------------------------------------------------------------------!


         !----- 2. Leaves. ----------------------------------------------------------------!
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 2. LEAF PROPERTIES (resolvable cohorts only): '
         write(unit=*,fmt='(7(a,1x))')               '  ICO',       ' IPFT','   MIN_WATER' &
                                             ,' MIN_WAT_IM2',' MAX_WAT_IM2','    MIN_TEMP' &
                                             ,'    MAX_TEMP'
         prob_leaf_loop: do ico = 1,cpatch%ncohorts
            if (.not. (y%leaf_resolvable(ico))) cycle prob_leaf_loop

            !----- Find the cohort-specific boundaries and print them. --------------------!
            ipft                  = cpatch%pft(ico)
            rk4min_leaf_water     = rk4min_veg_lwater * y%lai(ico)
            rk4min_leaf_water_im2 = rk4aux(ibuff)%rk4min_leaf_water_im2(ico)               &
                                  * (1.d0 - rk4eps)
            rk4max_leaf_water_im2 = rk4aux(ibuff)%rk4max_leaf_water_im2(ico)               &
                                  * (1.d0 + rk4eps)

            write(unit=*,fmt='(2(i5,1x),5(es12.5,1x))') ico,ipft,rk4min_leaf_water         &
                                                       ,rk4min_leaf_water_im2              &
                                                       ,rk4max_leaf_water_im2              &
                                                       ,rk4min_veg_temp ,rk4max_veg_temp
            !------------------------------------------------------------------------------!
         end do prob_leaf_loop
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         !---------------------------------------------------------------------------------!


         !----- 3. Wood. ------------------------------------------------------------------!
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 3. WOOD PROPERTIES (resolvable only): '
         write(unit=*,fmt='(6(a,1x))')               '  ICO',       ' IPFT','   MIN_WATER' &
                                             ,' MIN_WAT_IM2',' MAX_WAT_IM2','    MIN_TEMP' &
                                             ,'    MAX_TEMP'
         prob_wood_loop: do ico = 1,cpatch%ncohorts
            if (.not. (y%wood_resolvable(ico))) cycle prob_wood_loop


            !----- Find the cohort-specific boundaries and print them. --------------------!
            ipft                  = cpatch%pft(ico)
            rk4min_wood_water     = rk4min_veg_lwater * y%lai(ico)
            rk4min_wood_water_im2 = rk4aux(ibuff)%rk4min_wood_water_im2(ico)               &
                                  * (1.d0 - rk4eps)
            rk4max_wood_water_im2 = rk4aux(ibuff)%rk4max_wood_water_im2(ico)               &
                                  * (1.d0 + rk4eps)

            write(unit=*,fmt='(2(i5,1x),5(es12.5,1x))') ico,ipft,rk4min_wood_water         &
                                                       ,rk4min_wood_water_im2              &
                                                       ,rk4max_wood_water_im2              &
                                                       ,rk4min_veg_temp ,rk4max_veg_temp
            !------------------------------------------------------------------------------!
         end do prob_wood_loop
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         !---------------------------------------------------------------------------------!


         !----- 4. Temporary surface water (real or virtual). -----------------------------!
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 4. SURFACE WATER / VIRTUAL POOL PROPERTIES: '
         write(unit=*,fmt='(3(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_WMASS'
         write(unit=*,fmt='(3(es12.5,1x))') rk4min_sfcw_temp ,rk4max_sfcw_temp             &
                                           ,rk4min_sfcw_mass
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(78a)')         ('-',k=1,78)
         !---------------------------------------------------------------------------------!


         !----- 5. Soil (currently only the top layer). -----------------------------------!
         write(unit=*,fmt='(a)')           ' '
         write(unit=*,fmt='(a)')           ' 5. SOIL (TEXTURE CLASS AT TOP LAYER): '
         write(unit=*,fmt='(4(a,1x))')     '   MIN_WATER','   MAX_WATER','    MIN_TEMP'    &
                                          ,'    MAX_TEMP'
         write(unit=*,fmt='(4(es12.5,1x))') rk4aux(ibuff)%rk4min_soil_water(nzg)           &
                                           ,rk4aux(ibuff)%rk4max_soil_water(nzg)  &
                                           ,rk4min_soil_temp      ,rk4max_soil_temp
         write(unit=*,fmt='(78a)')         ('=',k=1,78)
         write(unit=*,fmt='(a)')           ' '
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine rk4_sanity_check
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This will print the values whenever the step didn't converge due to crazy values  !
   ! no matter how small the steps were (reject_step=.true.).                              !
   !---------------------------------------------------------------------------------------!
   subroutine print_sanity_check(y, csite, ipa)

      use rk4_coms              , only : rk4patchtype  & ! structure
                                       , rk4site       ! ! intent(in)
      use ed_state_vars         , only : sitetype      & ! structure
                                       , patchtype     ! ! structure
      use grid_coms             , only : nzg           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: y
      type(sitetype)     , target     :: csite
      integer            , intent(in) :: ipa
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)    , pointer    :: cpatch
      integer                         :: ico,k
      !------------------------------------------------------------------------------------!

      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(a,20x,a,20x,a)') '======','SANITY CHECK','======'
      write(unit=*,fmt='(78a)') ('=',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a5,4(1x,a12))') 'LEVEL','  SOIL_TEMPK','SOIL_FRACLIQ'             &
                                                ,'  SOIL_WATER','SOIL_MSTPOT'
      do k=rk4site%lsl,nzg
         write(unit=*,fmt='(i5,4(1x,es12.4))') &
              k, y%soil_tempk(k), y%soil_fracliq(k), y%soil_water(k),y%soil_mstpot(k)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a5,4(1x,a12))') 'LEVEL','  OLD_SOIL_T','OLD_SOIL_FLQ'             &
                                               &,'OLD_SOIL_H2O','OLD_SOIL_POT'
      do k=rk4site%lsl,nzg
         write(unit=*,fmt='(i5,4(1x,es12.4))')                                             &
              k, csite%soil_tempk(k,ipa), csite%soil_fracliq(k,ipa)                        &
               , csite%soil_water(k,ipa), csite%soil_mstpot (k,ipa)
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP=     ',y%can_temp
      write (unit=*,fmt='(a,1x,es12.4)') ' OLD_CAN_TEMP= ',csite%can_temp(ipa)
      write (unit=*,fmt='(a,1x,es12.4)') ' CAN_VAPOR=    ',y%can_shv
      write (unit=*,fmt='(a,1x,es12.4)') ' OLD_CAN_VAP=  ',csite%can_shv(ipa)
      write (unit=*,fmt='(a,1x,i12)')    ' #LEV_SFCH2O=  ',y%nlev_sfcwater
      write (unit=*,fmt='(a,1x,i12)')    ' OLD_#_SFCH2O= ',csite%nlev_sfcwater(ipa)
      if(y%nlev_sfcwater == 1) then
         write(unit=*,fmt='(a,1x,es12.4)') 'SFCWATER_TEMPK=',y%sfcwater_tempk(1)
      end if
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      cpatch => csite%patch(ipa)
      write (unit=*,fmt='(2(a5,1x),7(a12,1x))')                                            &
         '  COH','  PFT','         LAI','         WAI','         TAI',' LEAF_ENERGY'       &
                        ,' OLD_LEAF_EN','   LEAF_TEMP','OLD_LEAF_TMP'
      do ico = 1,cpatch%ncohorts
         if(y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),7(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%tai(ico),y%leaf_energy(ico)     &
                  ,cpatch%leaf_energy(ico),y%leaf_temp(ico),cpatch%leaf_temp(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(2(a5,1x),9(a12,1x))') &
         '  COH','  PFT','         LAI','         WAI','         TAI','  LEAF_WATER'       &
                        ,'OLD_LEAF_H2O','LEAF_H2O_IM2',' OLD_LFW_IM2','   LEAF_HCAP'       &
                        ,'   LEAF_FLIQ'
      do ico = 1,cpatch%ncohorts
         if(y%leaf_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),9(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%tai(ico),y%leaf_water(ico)      &
                  ,cpatch%leaf_water(ico),y%leaf_water_im2(ico),cpatch%leaf_water_im2(ico) &
                  ,cpatch%leaf_hcap(ico),y%leaf_hcap(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a)') ' '

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      cpatch => csite%patch(ipa)
      write (unit=*,fmt='(2(a5,1x),7(a12,1x))')                                            &
         '  COH','  PFT','         LAI','         WAI','         TAI',' WOOD_ENERGY'       &
                        ,' OLD_WOOD_EN','   WOOD_TEMP','OLD_WOOD_TMP'
      do ico = 1,cpatch%ncohorts
         if(y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),7(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%tai(ico),y%wood_energy(ico)     &
                  ,cpatch%wood_energy(ico),y%wood_temp(ico),cpatch%wood_temp(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write (unit=*,fmt='(2(a5,1x),9(a12,1x))') &
         '  COH','  PFT','         LAI','         WAI','         TAI','  WOOD_WATER'       &
                        ,'OLD_WOOD_H2O','WOOD_H2O_IM2',' OLD_WDW_IM2','   WOOD_HCAP'       &
                        ,'   WOOD_FLIQ'
      do ico = 1,cpatch%ncohorts
         if(y%wood_resolvable(ico)) then
            write(unit=*,fmt='(2(i5,1x),9(es12.4,1x))')                                    &
               ico,cpatch%pft(ico),y%lai(ico),y%wai(ico),y%tai(ico),y%wood_water(ico)      &
                  ,cpatch%wood_water(ico),y%wood_water_im2(ico),cpatch%wood_water_im2(ico) &
                  ,cpatch%wood_hcap(ico),y%wood_hcap(ico)
         end if
      end do
      write(unit=*,fmt='(78a)') ('-',k=1,78)
      write(unit=*,fmt='(a)') ' '


      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(78a)') ('=',k=1,78)
      write(unit=*,fmt='(a)') ' '

      return
   end subroutine print_sanity_check
   !=======================================================================================!
   !=======================================================================================!
end module rk4_integ_utils
!==========================================================================================!
!==========================================================================================!
