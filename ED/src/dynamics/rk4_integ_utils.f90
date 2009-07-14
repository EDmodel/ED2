!==========================================================================================!
!==========================================================================================!
! Subroutine odeint                                                                     !
!                                                                                          !
!     This subroutine will drive the integration of several ODEs that drive the fast-scale !
! state variables.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine odeint(h1,csite,ipa,isi,ipy,ifm,integration_buff)

   use ed_state_vars  , only : sitetype               & ! structure
                             , patchtype              ! ! structure
   use rk4_coms       , only : integration_vars    & ! structure
                             , rk4met                 & ! intent(in)
                             , rk4min_sfcwater_mass   & ! intent(in)
                             , maxstp                 & ! intent(in)
                             , tbeg                   & ! intent(in)
                             , tend                   & ! intent(in)
                             , dtrk4                  & ! intent(in)
                             , dtrk4i                 & ! intent(in)
                             , tiny_offset            ! ! intent(in)
   use rk4_stepper , only : rkqs                ! ! subroutine
   use ed_misc_coms   , only : fast_diagnostics       ! ! intent(in)
   use hydrology_coms , only : useRUNOFF              ! ! intent(in)
   use grid_coms      , only : nzg                    ! ! intent(in)
   use soil_coms      , only : dslz8                  & ! intent(in)
                             , runoff_time            ! ! intent(in)
   use consts_coms    , only : cliq8                  & ! intent(in)
                             , t3ple8                 & ! intent(in)
                             , tsupercool8            & ! intent(in)
                             , wdnsi8                 ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(integration_vars) , target      :: integration_buff ! RK4 variables
   type(sitetype)            , target      :: csite            ! Current site
   integer                   , intent(in)  :: ipa              ! Current patch ID
   integer                   , intent(in)  :: isi              ! Current site ID
   integer                   , intent(in)  :: ipy              ! Current polygon ID
   integer                   , intent(in)  :: ifm              ! Current grid ID
   real(kind=8)              , intent(in)  :: h1               ! First guess of delta-t
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

   !---------------------------------------------------------------------------------------!
   !    If top snow layer is too thin for computational stability, have it evolve in       !
   ! thermal equilibrium with top soil layer.                                              !
   !---------------------------------------------------------------------------------------!
   call redistribute_snow(integration_buff%initp, csite,ipa)
   call update_diagnostic_vars(integration_buff%initp,csite,ipa)



   !---------------------------------------------------------------------------------------!
   !     Create temporary patches.                                                         !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   call copy_rk4_patch(integration_buff%initp, integration_buff%y,cpatch)


   !---------------------------------------------------------------------------------------!
   ! Set initial time and stepsize.                                                        !
   !---------------------------------------------------------------------------------------!
   x = tbeg
   h = h1
   if (dtrk4 < 0.d0) h = -h1

   !---------------------------------------------------------------------------------------!
   ! Begin timestep loop                                                                   !
   !---------------------------------------------------------------------------------------!
   timesteploop: do i=1,maxstp

      !----- Get initial derivatives ------------------------------------------------------!
      call leaf_derivs(integration_buff%y,integration_buff%dydx,csite,ipa,isi,ipy)

      !----- Get scalings used to determine stability -------------------------------------!
      call get_yscal(integration_buff%y, integration_buff%dydx,h,integration_buff%yscal &
                       ,cpatch)

      !----- Be sure not to overstep ------------------------------------------------------!
      if((x+h-tend)*(x+h-tbeg) > 0.d0) h=tend-x

      !----- Take the step ----------------------------------------------------------------!
      call rkqs(integration_buff,x,h,hdid,hnext,csite,ipa,isi,ipy,ifm)

      !----- If the integration reached the next step, make some final adjustments --------!
      if((x-tend)*dtrk4 >= 0.d0)then

         csite%wbudget_loss2runoff(ipa) = 0.0
         csite%ebudget_loss2runoff(ipa) = 0.0
         ksn = integration_buff%y%nlev_sfcwater

         !---------------------------------------------------------------------------------!
         !   Make temporary surface liquid water disappear.  This will not happen          !
         ! immediately, but liquid water will decay with the time scale defined by         !
         ! runoff_time scale. If the time scale is too tiny, then it will be forced to be  !
         ! hdid (no reason to be faster than that).                                        !
         !---------------------------------------------------------------------------------!
         if (simplerunoff .and. ksn >= 1) then
         
            if (integration_buff%y%sfcwater_mass(ksn) > 0.d0   .and.                       &
                integration_buff%y%sfcwater_fracliq(ksn) > 1.d-1) then
               wfreeb = min(1.d0,dtrk4*runoff_time_i)                                      &
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

               call redistribute_snow(integration_buff%y,csite,ipa)
               call update_diagnostic_vars(integration_buff%y,csite,ipa)

               !----- Compute runoff for output -------------------------------------------!
               if(fast_diagnostics) then
                  csite%runoff(ipa) = csite%runoff(ipa)                                    &
                                    + sngloff(wfreeb * dtrk4i,tiny_offset)
                  csite%avg_runoff(ipa) = csite%avg_runoff(ipa)                            &
                                        + sngloff(wfreeb * dtrk4i,tiny_offset)
                  csite%avg_runoff_heat(ipa) = csite%avg_runoff_heat(ipa)                  &
                                             + sngloff(qwfree * dtrk4i,tiny_offset)
                  csite%wbudget_loss2runoff(ipa) = sngloff(wfreeb,tiny_offset)
                  csite%ebudget_loss2runoff(ipa) = sngloff(qwfree,tiny_offset)
               end if

            else
               csite%runoff(ipa)              = 0.0
               csite%avg_runoff(ipa)          = 0.0
               csite%avg_runoff_heat(ipa)     = 0.0
               csite%wbudget_loss2runoff(ipa) = 0.0
               csite%ebudget_loss2runoff(ipa) = 0.0
            end if
         else
            csite%runoff(ipa)              = 0.0
            csite%avg_runoff(ipa)          = 0.0
            csite%avg_runoff_heat(ipa)     = 0.0
            csite%wbudget_loss2runoff(ipa) = 0.0
            csite%ebudget_loss2runoff(ipa) = 0.0
         end if

         !------ Copying the temporary patch to the next intermediate step ----------------!
         call copy_rk4_patch(integration_buff%y,integration_buff%initp, cpatch)
         !------ Updating the substep for next time and leave -----------------------------!
         csite%htry(ipa) = sngl(hnext)

         return
      end if
      
      !----- Use hnext as the next substep ------------------------------------------------!
      h = hnext
   end do timesteploop

   !----- If it reached this point, that is really bad news... ----------------------------!
   print*,'Too many steps in routine odeint'
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
subroutine copy_met_2_rk4met(rhos,vels,atm_tmp,atm_shv,atm_co2,zoff,exner,pcpg,qpcpg,dpcpg &
                            ,prss,geoht,lsl,lon,lat)
   use rk4_coms   , only : rk4met ! ! structure
   use consts_coms, only : cpi8   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: lsl
   real   , intent(in) :: rhos
   real   , intent(in) :: vels
   real   , intent(in) :: atm_tmp
   real   , intent(in) :: atm_shv
   real   , intent(in) :: atm_co2
   real   , intent(in) :: zoff
   real   , intent(in) :: exner
   real   , intent(in) :: pcpg
   real   , intent(in) :: qpcpg
   real   , intent(in) :: dpcpg
   real   , intent(in) :: prss
   real   , intent(in) :: geoht
   real   , intent(in) :: lon
   real   , intent(in) :: lat
   !---------------------------------------------------------------------------------------!

   
   rk4met%lsl     = lsl
   !----- Converting to double precision. -------------------------------------------------!
   rk4met%rhos      = dble(rhos   )
   rk4met%vels      = dble(vels   )
   rk4met%atm_tmp   = dble(atm_tmp)
   rk4met%atm_shv   = dble(atm_shv)
   rk4met%atm_co2   = dble(atm_co2)
   rk4met%zoff      = dble(zoff   )
   rk4met%exner     = dble(exner  )
   rk4met%pcpg      = dble(pcpg   )
   rk4met%qpcpg     = dble(qpcpg  )
   rk4met%dpcpg     = dble(dpcpg  )
   rk4met%prss      = dble(prss   )
   rk4met%geoht     = dble(geoht  )
   
   rk4met%atm_theta = cpi8 * rk4met%exner * rk4met%atm_tmp

   rk4met%lon       = dble(lon    )
   rk4met%lat       = dble(lat    )

   return
end subroutine copy_met_2_rk4met
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies that variables that are integrated by the Runge-Kutta solver   !
! to a buffer structure.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine copy_patch_init(sourcesite,ipa,targetp)
   use ed_state_vars        , only : sitetype              & ! structure
                                   , patchtype             ! ! structure
   use grid_coms            , only : nzg                   & ! intent(in)
                                   , nzs                   ! ! intent(in) 
   use ed_misc_coms         , only : fast_diagnostics      ! ! intent(in)
   use consts_coms          , only : cpi8                  ! ! intent(in)
   use rk4_coms             , only : rk4patchtype          & ! structure
                                   , rk4met                & ! structure
                                   , hcapveg_ref           & ! intent(in)
                                   , hcapveg_coh_min       & ! intent(in)
                                   , rk4eps                & ! intent(in)
                                   , min_height            & ! intent(in)
                                   , any_solvable          & ! intent(out)
                                   , zoveg                 & ! intent(out)
                                   , zveg                  & ! intent(out)
                                   , wcapcan               & ! intent(out)
                                   , wcapcani              & ! intent(out)
                                   , hcapcani              & ! intent(out)
                                   , rk4water_stab_thresh  & ! intent(in)
                                   , rk4min_sfcwater_mass  ! ! intent(in)
   use ed_max_dims             , only : n_pft                 ! ! intent(in)
   use canopy_radiation_coms, only : tai_min               ! ! intent(in)
   use therm_lib            , only : qwtk8                 ! ! subroutine
   use allometry            , only : dbh2bl                ! ! function
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)    , target     :: targetp
   type(sitetype)        , target     :: sourcesite
   integer               , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch
   real(kind=8)                       :: hvegpat_min
   real(kind=8)                       :: hcap_scale
   integer                            :: ico
   integer                            :: ipft
   integer                            :: k
   !---------------------------------------------------------------------------------------!


   targetp%can_temp  = dble(sourcesite%can_temp(ipa))
   targetp%can_shv   = dble(sourcesite%can_shv(ipa))
   targetp%can_co2   = dble(sourcesite%can_co2(ipa))

   do k = rk4met%lsl, nzg
      targetp%soil_water(k)   = dble(sourcesite%soil_water(k,ipa))
      targetp%soil_energy(k)  = dble(sourcesite%soil_energy(k,ipa))
      targetp%soil_tempk(k)   = dble(sourcesite%soil_tempk(k,ipa))
      targetp%soil_fracliq(k) = dble(sourcesite%soil_fracliq(k,ipa))
   end do

   do k = 1, nzs
      targetp%sfcwater_mass(k)    = dble(sourcesite%sfcwater_mass(k,ipa))
      targetp%sfcwater_depth(k)   = dble(sourcesite%sfcwater_depth(k,ipa))
      !----- Converting sfcwater_energy to J/m² inside the Runge-Kutta integrator. --------!
      targetp%sfcwater_energy(k)  = dble(sourcesite%sfcwater_energy(k,ipa))                &
                                  * dble(sourcesite%sfcwater_mass(k,ipa))
      targetp%sfcwater_tempk(k)   = dble(sourcesite%sfcwater_tempk(k,ipa))
      targetp%sfcwater_fracliq(k) = dble(sourcesite%sfcwater_fracliq(k,ipa))
   end do


   targetp%ustar = dble(sourcesite%ustar(ipa))
   targetp%cstar = dble(sourcesite%cstar(ipa))
   targetp%tstar = dble(sourcesite%tstar(ipa))
   targetp%qstar = dble(sourcesite%qstar(ipa))


   targetp%upwp = dble(sourcesite%upwp(ipa))
   targetp%wpwp = dble(sourcesite%wpwp(ipa))
   targetp%tpwp = dble(sourcesite%tpwp(ipa))
   targetp%qpwp = dble(sourcesite%qpwp(ipa))

  
   targetp%nlev_sfcwater = sourcesite%nlev_sfcwater(ipa)


   !----- The virtual pools should be always zero, they are temporary entities ------------!
   targetp%virtual_water = 0.0d0
   targetp%virtual_heat  = 0.0d0
   targetp%virtual_depth = 0.0d0

   if (targetp%nlev_sfcwater == 0) then
      targetp%virtual_flag = 2
   else
      if (targetp%sfcwater_mass(1) < rk4min_sfcwater_mass) then
         targetp%virtual_flag = 2
      elseif (targetp%sfcwater_mass(1) < rk4water_stab_thresh) then
         targetp%virtual_flag = 1
      else
         targetp%virtual_flag = 0
      end if
   end if

   !---------------------------------------------------------------------------------------!
   !     Here we find the minimum patch-level leaf heat capacity.  If the total patch leaf !
   ! heat capacity is less than this, we scale the cohorts heat capacity inside the        !
   ! integrator, so it preserves the proportional heat capacity and prevents the pool to   !
   ! be too small.                                                                         !
   !---------------------------------------------------------------------------------------!
   cpatch => sourcesite%patch(ipa)
   sourcesite%hcapveg(ipa) = 0.
   sourcesite%lai(ipa)     = 0.
   sourcesite%wpa(ipa)     = 0.
   sourcesite%wai(ipa)     = 0.
   do ico=1,cpatch%ncohorts
      sourcesite%hcapveg(ipa) = sourcesite%hcapveg(ipa) + cpatch%hcapveg(ico)
      sourcesite%lai(ipa)     = sourcesite%lai(ipa)     + cpatch%lai(ico)
      sourcesite%wpa(ipa)     = sourcesite%wpa(ipa)     + cpatch%wpa(ico)
      sourcesite%wai(ipa)     = sourcesite%wai(ipa)     + cpatch%wai(ico)
   end do
   
   any_solvable = .false.
   do ico=1, cpatch%ncohorts
      !----- Copying the flag that determines whether this cohort is numerically stable. --!
      targetp%solvable(ico) = cpatch%solvable(ico)
      if (targetp%solvable(ico)) any_solvable = .true.
   end do

   if ((sourcesite%lai(ipa)+sourcesite%wai(ipa)) > tai_min) then
      hvegpat_min = hcapveg_ref * max(dble(cpatch%hite(1)),min_height)
      hcap_scale  = max(1.d0,hvegpat_min / sourcesite%hcapveg(ipa))
   else
      hcap_scale  = 1.d0
   end if

   do ico = 1,cpatch%ncohorts
      ipft=cpatch%pft(ico)
    
      !----- Copying the leaf area index and total (leaf+branch+twig) area index. ---------!
      targetp%lai(ico)  = dble(cpatch%lai(ico))
      targetp%wpa(ico)  = dble(cpatch%wpa(ico))
      targetp%tai(ico) = targetp%lai(ico) + dble(cpatch%wai(ico))

      !------------------------------------------------------------------------------------!
      !    If the cohort is too small, we give some extra heat capacity, so the model can  !
      ! run in a stable range inside the integrator.  At the end this extra heat capacity  !
      ! will be removed.                                                                   !
      !------------------------------------------------------------------------------------!
      targetp%hcapveg(ico) = dble(cpatch%hcapveg(ico)) * hcap_scale

      !------------------------------------------------------------------------------------!
      !     Checking whether this is considered a "safe" one or not.  In case it is, we    !
      ! copy water, temperature, and liquid fraction, and scale energy and heat capacity   !
      ! as needed.  Otherwise, just fill with some safe values, but the cohort won't be    !
      ! really solved.                                                                     !
      !------------------------------------------------------------------------------------!
      targetp%veg_water(ico)     = dble(cpatch%veg_water(ico))

      if (targetp%solvable(ico)) then
         targetp%veg_energy(ico)    = dble(cpatch%veg_energy(ico))                         &
                                    + (targetp%hcapveg(ico)-dble(cpatch%hcapveg(ico)))     &
                                    * dble(cpatch%veg_temp(ico))
         call qwtk8(targetp%veg_energy(ico),targetp%veg_water(ico),targetp%hcapveg(ico)    &
                   ,targetp%veg_temp(ico),targetp%veg_fliq(ico))
      else
         targetp%veg_fliq(ico)   = dble(cpatch%veg_fliq(ico))
         targetp%veg_temp(ico)   = dble(cpatch%veg_temp(ico))
         targetp%veg_energy(ico) = targetp%hcapveg(ico) * targetp%veg_temp(ico)
      end if
   end do
   !----- Diagnostics variables -----------------------------------------------------------!
   if(fast_diagnostics) then

      targetp%wbudget_loss2atm   = dble(sourcesite%wbudget_loss2atm(ipa)  )
      targetp%ebudget_loss2atm   = dble(sourcesite%ebudget_loss2atm(ipa)  )
      targetp%co2budget_loss2atm = dble(sourcesite%co2budget_loss2atm(ipa))
      targetp%ebudget_latent     = dble(sourcesite%ebudget_latent(ipa)    )
      targetp%avg_carbon_ac      = dble(sourcesite%avg_carbon_ac(ipa)     )
      targetp%avg_vapor_vc       = dble(sourcesite%avg_vapor_vc(ipa)      )
      targetp%avg_dew_cg         = dble(sourcesite%avg_dew_cg(ipa)        )
      targetp%avg_vapor_gc       = dble(sourcesite%avg_vapor_gc(ipa)      )
      targetp%avg_wshed_vg       = dble(sourcesite%avg_wshed_vg(ipa)      )
      targetp%avg_vapor_ac       = dble(sourcesite%avg_vapor_ac(ipa)      )
      targetp%avg_transp         = dble(sourcesite%avg_transp(ipa)        )
      targetp%avg_evap           = dble(sourcesite%avg_evap(ipa)          )
      targetp%avg_drainage       = dble(sourcesite%avg_drainage(ipa)      )
      targetp%avg_netrad         = dble(sourcesite%avg_netrad(ipa)        )
      targetp%avg_sensible_vc    = dble(sourcesite%avg_sensible_vc(ipa)   )
      targetp%avg_sensible_2cas  = dble(sourcesite%avg_sensible_2cas(ipa) )
      targetp%avg_qwshed_vg      = dble(sourcesite%avg_qwshed_vg(ipa)     )
      targetp%avg_sensible_gc    = dble(sourcesite%avg_sensible_gc(ipa)   )
      targetp%avg_sensible_ac    = dble(sourcesite%avg_sensible_ac(ipa)   )
      targetp%avg_sensible_tot   = dble(sourcesite%avg_sensible_tot(ipa)  )

      do k = rk4met%lsl, nzg
         targetp%avg_sensible_gg(k) = dble(sourcesite%avg_sensible_gg(k,ipa))
         targetp%avg_smoist_gg(k)   = dble(sourcesite%avg_smoist_gg(k,ipa)  )
         targetp%avg_smoist_gc(k)   = dble(sourcesite%avg_smoist_gc(k,ipa)  )
      end do
   end if

   return
end subroutine copy_patch_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutines increment the derivative into the previous guess to create the new   !
! guess.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine inc_rk4_patch(rkp, inc, fac, cpatch)
   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         ! ! structure
   use rk4_coms      , only : rk4patchtype      & ! structure
                            , rk4met            ! ! intent(in)
   use grid_coms     , only : nzg               & ! intent(in)
                            , nzs               ! ! intent(in)
   use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)
  
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


   rkp%can_temp = rkp%can_temp  + fac * inc%can_temp
   rkp%can_shv  = rkp%can_shv   + fac * inc%can_shv
   rkp%can_co2  = rkp%can_co2   + fac * inc%can_co2

   do k=rk4met%lsl,nzg
      rkp%soil_water(k)       = rkp%soil_water(k)  + fac * inc%soil_water(k)
      rkp%soil_energy(k)      = rkp%soil_energy(k) + fac * inc%soil_energy(k)
   end do

   do k=1,rkp%nlev_sfcwater
      rkp%sfcwater_mass(k)   = rkp%sfcwater_mass(k)   + fac * inc%sfcwater_mass(k)
      rkp%sfcwater_energy(k) = rkp%sfcwater_energy(k) + fac * inc%sfcwater_energy(k)
      rkp%sfcwater_depth(k)  = rkp%sfcwater_depth(k)  + fac * inc%sfcwater_depth(k)
   end do

   rkp%virtual_heat  = rkp%virtual_heat  + fac * inc%virtual_heat
   rkp%virtual_water = rkp%virtual_water + fac * inc%virtual_water
   rkp%virtual_depth = rkp%virtual_depth + fac * inc%virtual_depth

  
   rkp%upwp = rkp%upwp + fac * inc%upwp
   rkp%wpwp = rkp%wpwp + fac * inc%wpwp
   rkp%tpwp = rkp%tpwp + fac * inc%tpwp
   rkp%qpwp = rkp%qpwp + fac * inc%qpwp

  
   do ico = 1,cpatch%ncohorts
      rkp%veg_water(ico)     = rkp%veg_water(ico) + fac * inc%veg_water(ico)
      rkp%veg_energy(ico)    = rkp%veg_energy(ico) + fac * inc%veg_energy(ico)
   enddo

   if(fast_diagnostics) then

      rkp%wbudget_loss2atm   = rkp%wbudget_loss2atm   + fac * inc%wbudget_loss2atm
      rkp%ebudget_loss2atm   = rkp%ebudget_loss2atm   + fac * inc%ebudget_loss2atm
      rkp%co2budget_loss2atm = rkp%co2budget_loss2atm + fac * inc%co2budget_loss2atm
      rkp%ebudget_latent     = rkp%ebudget_latent     + fac * inc%ebudget_latent

      rkp%avg_carbon_ac      = rkp%avg_carbon_ac      + fac * inc%avg_carbon_ac
      
      rkp%avg_vapor_vc       = rkp%avg_vapor_vc       + fac * inc%avg_vapor_vc
      rkp%avg_dew_cg         = rkp%avg_dew_cg         + fac * inc%avg_dew_cg
      rkp%avg_vapor_gc       = rkp%avg_vapor_gc       + fac * inc%avg_vapor_gc
      rkp%avg_wshed_vg       = rkp%avg_wshed_vg       + fac * inc%avg_wshed_vg
      rkp%avg_vapor_ac       = rkp%avg_vapor_ac       + fac * inc%avg_vapor_ac
      rkp%avg_transp         = rkp%avg_transp         + fac * inc%avg_transp
      rkp%avg_evap           = rkp%avg_evap           + fac * inc%avg_evap
      rkp%avg_drainage       = rkp%avg_drainage       + fac * inc%avg_drainage
      rkp%avg_netrad         = rkp%avg_netrad         + fac * inc%avg_netrad
      rkp%avg_sensible_vc    = rkp%avg_sensible_vc    + fac * inc%avg_sensible_vc
      rkp%avg_sensible_2cas  = rkp%avg_sensible_2cas  + fac * inc%avg_sensible_2cas
      rkp%avg_qwshed_vg      = rkp%avg_qwshed_vg      + fac * inc%avg_qwshed_vg
      rkp%avg_sensible_gc    = rkp%avg_sensible_gc    + fac * inc%avg_sensible_gc
      rkp%avg_sensible_ac    = rkp%avg_sensible_ac    + fac * inc%avg_sensible_ac
      rkp%avg_sensible_tot   = rkp%avg_sensible_tot   + fac * inc%avg_sensible_tot

      do k=rk4met%lsl,nzg
         rkp%avg_sensible_gg(k)  = rkp%avg_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
         rkp%avg_smoist_gg(k)    = rkp%avg_smoist_gg(k)    + fac * inc%avg_smoist_gg(k)  
         rkp%avg_smoist_gc(k)    = rkp%avg_smoist_gc(k)    + fac * inc%avg_smoist_gc(k)  
      end do

   end if

   return
end subroutine inc_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the error scale for the integrated variables, which will be     !
! later used to define the relative error.                                                 !
!------------------------------------------------------------------------------------------!
subroutine get_yscal(y, dy, htry, yscal, cpatch)
   use ed_state_vars        , only : patchtype            ! ! structure
   use rk4_coms             , only : rk4patchtype         & ! structure
                                   , rk4met               & ! intent(in)
                                   , tiny_offset          & ! intent(in)
                                   , huge_offset          & ! intent(in)
                                   , rk4water_stab_thresh & ! intent(in)
                                   , rk4min_sfcwater_mass & ! intent(in)
                                   , rk4dry_veg_lwater    ! ! intent(in)
   use grid_coms            , only : nzg                  & ! intent(in)
                                   , nzs                  ! ! intent(in)
   use consts_coms          , only : cliq8                & ! intent(in)
                                   , qliqt38              ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype), target     :: y                ! Structure with the guesses
   type(rk4patchtype), target     :: dy               ! Structure with their derivatives
   type(rk4patchtype), target     :: yscal            ! Structure with their scales
   type(patchtype)   , target     :: cpatch           ! Current patch
   real(kind=8)      , intent(in) :: htry             ! Time-step we are trying
   !----- Local variables -----------------------------------------------------------------!
   real(kind=8)                   :: tot_sfcw_mass    ! Total surface water/snow mass
   integer                        :: k                ! Counter
   integer                        :: ico              ! Current cohort ID
   !---------------------------------------------------------------------------------------!

  
   yscal%can_temp = abs(y%can_temp) + abs(dy%can_temp*htry) + tiny_offset
   yscal%can_shv  = abs(y%can_shv)  + abs(dy%can_shv*htry)  + tiny_offset
   yscal%can_co2  = abs(y%can_co2)  + abs(dy%can_co2*htry)  + tiny_offset
  
   yscal%upwp = max(abs(y%upwp) + abs(dy%upwp*htry),1.d0)
   yscal%wpwp = max(abs(y%wpwp) + abs(dy%wpwp*htry),1.d0)


  
   do k=rk4met%lsl,nzg
      yscal%soil_water(k)  = abs(y%soil_water(k))  + abs(dy%soil_water(k)*htry)            &
                           + tiny_offset
      yscal%soil_energy(k) = abs(y%soil_energy(k)) + abs(dy%soil_energy(k)*htry)           &
                           + tiny_offset
   end do

   tot_sfcw_mass = 0.d0
   do k=1,y%nlev_sfcwater
      tot_sfcw_mass = tot_sfcw_mass + y%sfcwater_mass(k)
   end do
   tot_sfcw_mass = abs(tot_sfcw_mass)
   
   if (tot_sfcw_mass > 1.d-2*rk4water_stab_thresh) then
      !----- Computationally stable layer. ------------------------------------------------!
      do k=1,nzs
         yscal%sfcwater_mass(k)   = abs(y%sfcwater_mass(k))                                &
                                  + abs(dy%sfcwater_mass(k)*htry)
         yscal%sfcwater_energy(k) = abs(y%sfcwater_energy(k))                              &
                                  + abs(dy%sfcwater_energy(k)*htry)
         yscal%sfcwater_depth(k)  = abs(y%sfcwater_depth(k))                               &
                                  + abs(dy%sfcwater_depth(k)*htry)
      end do
   else
      !----- Low stability threshold ------------------------------------------------------!
      do k=1,nzs
         if(abs(y%sfcwater_mass(k)) > rk4min_sfcwater_mass)then
            yscal%sfcwater_mass(k) = 1.d-2*rk4water_stab_thresh
            yscal%sfcwater_energy(k) = ( yscal%sfcwater_mass(k) / abs(y%sfcwater_mass(k))) &
                                     * ( abs( y%sfcwater_energy(k))                        &
                                       + abs(dy%sfcwater_energy(k)*htry))
            yscal%sfcwater_depth(k)  = ( yscal%sfcwater_mass(k) / abs(y%sfcwater_mass(k))) &
                                     * abs(y%sfcwater_depth(k))                            &
                                     + abs(dy%sfcwater_depth(k)*htry)
         else
            yscal%sfcwater_mass(k)   = huge_offset
            yscal%sfcwater_energy(k) = huge_offset
            yscal%sfcwater_depth(k)  = huge_offset
         end if
      end do
   end if

   !----- Scale for the virtual water pools -----------------------------------------------!
   if (abs(y%virtual_water) > 1.d-2*rk4water_stab_thresh) then
      yscal%virtual_water = abs(y%virtual_water) + abs(dy%virtual_water*htry)
      yscal%virtual_heat  = abs(y%virtual_heat) + abs(dy%virtual_heat*htry)
   elseif (abs(y%virtual_water) > rk4min_sfcwater_mass) then
      yscal%virtual_water = 1.d-2*rk4water_stab_thresh
      yscal%virtual_heat  = (yscal%virtual_water / abs(y%virtual_water))                   &
                          * (abs(y%virtual_heat) + abs(dy%virtual_heat*htry))
   else
      yscal%virtual_water = huge_offset
      yscal%virtual_heat  = huge_offset
   end if

   !---------------------------------------------------------------------------------------!
   !    Scale for leaf water and energy. In case the plants have few or no leaves, or the  !
   ! plant is buried in snow, we assign huge values for typical scale, thus preventing     !
   ! unecessary small steps.                                                               !
   !    Also, if the cohort is tiny and has almost no water, make the scale less strict.   !
   !---------------------------------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      if (.not. y%solvable(ico)) then
         yscal%solvable(ico)   = .false.
         yscal%veg_water(ico)  = huge_offset
         yscal%veg_energy(ico) = huge_offset
      elseif (y%veg_water(ico) > rk4dry_veg_lwater*y%tai(ico)) then
         yscal%solvable(ico)   = .true.
         yscal%veg_water(ico)  = abs(y%veg_water(ico)) + abs(dy%veg_water(ico)*htry)
         yscal%veg_energy(ico) = abs(y%veg_energy(ico)) + abs(dy%veg_energy(ico)*htry)
      else
         yscal%solvable(ico)   = .true.
         yscal%veg_water(ico)  = rk4dry_veg_lwater*y%tai(ico)
         yscal%veg_energy(ico) = max(yscal%veg_water(ico)*qliqt38                          &
                                    ,abs(y%veg_energy(ico)) + abs(dy%veg_energy(ico)*htry))
      end if
   end do


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

   use rk4_coms              , only : rk4patchtype  & ! structure
                                    , rk4eps        & ! intent(in)
                                    , rk4met        ! ! intent(in)
   use ed_state_vars         , only : patchtype     ! ! structure
   use grid_coms             , only : nzg           ! ! intent(in)
   use ed_misc_coms             , only : integ_err     & ! intent(in)
                                    , record_err    ! ! intent(in)
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
   integer                          :: k                ! Counter
   !---------------------------------------------------------------------------------------!

   !----- Initialize error ----------------------------------------------------------------!
   errmax = 0.d0

   !---------------------------------------------------------------------------------------!
   !    We know check each variable error, and keep track of the worst guess, which will   !
   ! be our worst guess in the end.                                                        !
   !---------------------------------------------------------------------------------------!
   
   err    = abs(yerr%can_temp/yscal%can_temp)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(1,1) = integ_err(1,1) + 1_8 

   err    = abs(yerr%can_shv/yscal%can_shv)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(2,1) = integ_err(2,1) + 1_8 

   err    = abs(yerr%can_co2/yscal%can_co2)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(3,1) = integ_err(3,1) + 1_8 
  
   do k=rk4met%lsl,nzg
      err    = abs(yerr%soil_water(k)/yscal%soil_water(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(3+k,1) = integ_err(3+k,1) + 1_8 
   end do

   do k=rk4met%lsl,nzg
      err    = abs(yerr%soil_energy(k)/yscal%soil_energy(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(15+k,1) = integ_err(15+k,1) + 1_8      
   enddo

   do k=1,ytemp%nlev_sfcwater
      err = abs(yerr%sfcwater_energy(k) / yscal%sfcwater_energy(k))
      errmax = max(errmax,err)
      if(record_err .and. err .gt. rk4eps) integ_err(27+k,1) = integ_err(27+k,1) + 1_8      
   enddo

   do k=1,ytemp%nlev_sfcwater
      err    = abs(yerr%sfcwater_mass(k) / yscal%sfcwater_mass(k))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(32+k,1) = integ_err(32+k,1) + 1_8      
   enddo

   err    = abs(yerr%virtual_heat/yscal%virtual_heat)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(38,1) = integ_err(38,1) + 1_8      

   err    = abs(yerr%virtual_water/yscal%virtual_water)
   errmax = max(errmax,err)
   if(record_err .and. err > rk4eps) integ_err(39,1) = integ_err(39,1) + 1_8      

   !---------------------------------------------------------------------------------------!
   !     Getting the worst error only amongst the cohorts in which leaf properties were    !
   ! computed.                                                                             !
   !---------------------------------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      errh2oMAX = 0.d0
      erreneMAX = 0.d0
      if (yscal%solvable(ico)) then
         errh2o     = abs(yerr%veg_water(ico)/yscal%veg_water(ico))
         errene     = abs(yerr%veg_energy(ico)/yscal%veg_energy(ico))
         errmax     = max(errmax,errh2o,errene)
         errh2oMAX  = max(errh2oMAX,errh2o)
         erreneMAX  = max(erreneMAX,errene)
      end if
   end do
   if(cpatch%ncohorts > 0 .and. record_err) then
      if(errh2oMAX > rk4eps) integ_err(40,1) = integ_err(40,1) + 1_8
      if(erreneMAX > rk4eps) integ_err(41,1) = integ_err(41,1) + 1_8
   end if

   return
end subroutine get_errmax
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine print_errmax(errmax,yerr,yscal,cpatch,y,ytemp)
   use rk4_coms              , only : rk4patchtype  & ! Structure
                                    , rk4eps        & ! intent(in)
                                    , rk4met        ! ! intent(in)
   use ed_state_vars         , only : patchtype     ! ! Structure
   use grid_coms             , only : nzg           & ! intent(in)
                                    , nzs           ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target       :: yerr,yscal,y,ytemp
   type(patchtype)    , target       :: cpatch
   real(kind=8)       , intent(out)  :: errmax
   !----- Local variables -----------------------------------------------------------------!
   integer                           :: ico
   integer                           :: k
   logical                           :: troublemaker
   !----- Constants -----------------------------------------------------------------------!
   character(len=28)  , parameter    :: onefmt = '(a16,1x,3(es12.4,1x),11x,l1)'
   character(len=34)  , parameter    :: lyrfmt = '(a16,1x,i6,1x,3(es12.4,1x),11x,l1)'
   character(len=34)  , parameter    :: cohfmt = '(a16,1x,i6,1x,7(es12.4,1x),11x,l1)'
   !----- Functions -----------------------------------------------------------------------!
   logical            , external     :: large_error
   !---------------------------------------------------------------------------------------!


   write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) '  >>>>> PRINTING MAXIMUM ERROR INFORMATION: '
   write(unit=*,fmt='(a)'      ) 
   write(unit=*,fmt='(a)'      ) ' Patch level variables, single layer:'
   write(unit=*,fmt='(5(a,1x))')  'Name            ','   Max.Error','   Abs.Error'&
                                &,'       Scale','Problem(T|F)'

   errmax       = max(0.0,abs(yerr%can_temp/yscal%can_temp))
   troublemaker = large_error(yerr%can_temp,yscal%can_temp)
   write(unit=*,fmt=onefmt) 'CAN_TEMP:',errmax,yerr%can_temp,yscal%can_temp,troublemaker

   errmax       = max(errmax,abs(yerr%can_shv/yscal%can_shv))
   troublemaker = large_error(yerr%can_shv,yscal%can_shv)
   write(unit=*,fmt=onefmt) 'CAN_SHV:',errmax,yerr%can_shv,yscal%can_shv,troublemaker

   errmax = max(errmax,abs(yerr%can_co2/yscal%can_co2))
   troublemaker = large_error(yerr%can_co2,yscal%can_co2)
   write(unit=*,fmt=onefmt) 'CAN_CO2:',errmax,yerr%can_co2,yscal%can_co2,troublemaker

  
   errmax = max(errmax,abs(yerr%virtual_heat/yscal%virtual_heat))
   troublemaker = large_error(yerr%virtual_heat,yscal%virtual_heat)
   write(unit=*,fmt=onefmt) 'VIRTUAL_HEAT:',errmax,yerr%virtual_heat,yscal%virtual_heat    &
                                           ,troublemaker

   errmax = max(errmax,abs(yerr%virtual_water/yscal%virtual_water))
   troublemaker = large_error(yerr%virtual_water,yscal%virtual_water)
   write(unit=*,fmt=onefmt) 'VIRTUAL_WATER:',errmax,yerr%virtual_water,yscal%virtual_water &
                                            ,troublemaker

   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) ' Patch level variables, soil layers:'
   write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'               &
                                &,'   Abs.Error','       Scale','Problem(T|F)'

   do k=rk4met%lsl,nzg
      errmax = max(errmax,abs(yerr%soil_water(k)/yscal%soil_water(k)))
      troublemaker = large_error(yerr%soil_water(k),yscal%soil_water(k))
      write(unit=*,fmt=lyrfmt) 'SOIL_WATER:',k,errmax,yerr%soil_water(k)                   &
                                            ,yscal%soil_water(k),troublemaker

      errmax       = max(errmax,abs(yerr%soil_energy(k)/yscal%soil_energy(k)))
      troublemaker = large_error(yerr%soil_energy(k),yscal%soil_energy(k))
      write(unit=*,fmt=lyrfmt) 'SOIL_ENERGY:',k,errmax,yerr%soil_energy(k)                 &
                                             ,yscal%soil_energy(k),troublemaker
   enddo

   if (yerr%nlev_sfcwater > 0) then
      write(unit=*,fmt='(a)'  ) 
      write(unit=*,fmt='(80a)') ('-',k=1,80)
      write(unit=*,fmt='(a)'      ) ' Patch level variables, water/snow layers:'
      write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'      &
                                &,'   Abs.Error','       Scale','Problem(T|F)'
      do k=1,yerr%nlev_sfcwater
         errmax       = max(errmax,abs(yerr%sfcwater_energy(k)/yscal%sfcwater_energy(k)))
         troublemaker = large_error(yerr%sfcwater_energy(k),yscal%sfcwater_energy(k))
         write(unit=*,fmt=lyrfmt) 'SFCWATER_ENERGY:',k,errmax,yerr%sfcwater_energy(k)      &
                                                    ,yscal%sfcwater_energy(k),troublemaker

         errmax       = max(errmax,abs(yerr%sfcwater_mass(k)/yscal%sfcwater_mass(k)))
         troublemaker = large_error(yerr%sfcwater_mass(k),yscal%sfcwater_mass(k))
         write(unit=*,fmt=lyrfmt) 'SFCWATER_MASS:',k,errmax,yerr%sfcwater_mass(k)          &
                                                  ,yscal%sfcwater_mass(k),troublemaker
      end do
   end if

   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) ' Cohort_level variables (only the solvable ones):'
   write(unit=*,fmt='(7(a,1x))')  'Name            ','         PFT','         LAI'         &
                                     ,'         WPA','         TAI','   Max.Error'         &
                                     ,'   Abs.Error','       Scale','Problem(T|F)'
   do ico = 1,cpatch%ncohorts
      if (y%solvable(ico)) then
         errmax       = max(errmax,abs(yerr%veg_water(ico)/yscal%veg_water(ico)))
         troublemaker = large_error(yerr%veg_water(ico),yscal%veg_water(ico))
         write(unit=*,fmt=cohfmt) 'VEG_WATER:',cpatch%pft(ico),y%lai(ico),y%wpa(ico)       &
                                              ,y%tai(ico),errmax,yerr%veg_water(ico)       &
                                              ,yscal%veg_water(ico),troublemaker
              

         errmax       = max(errmax,abs(yerr%veg_energy(ico)/yscal%veg_energy(ico)))
         troublemaker = large_error(yerr%veg_energy(ico),yscal%veg_energy(ico))
         write(unit=*,fmt=cohfmt) 'VEG_ENERGY:',cpatch%pft(ico),cpatch%lai(ico),y%wpa(ico) &
                                               ,y%tai(ico),errmax,yerr%veg_energy(ico)     &
                                               ,yscal%veg_energy(ico),troublemaker
      end if
   end do

   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)

   return
end subroutine print_errmax
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This function simply checks whether the relative error is large or not.               !
!------------------------------------------------------------------------------------------!
logical function large_error(err,scal)
   use rk4_coms , only : rk4eps ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real(kind=8), intent(in) :: err  ! Absolute error
   real(kind=8), intent(in) :: scal ! Characteristic scale
   !---------------------------------------------------------------------------------------!
   if(scal > 0.d0) then
      large_error = abs(err/scal)/rk4eps > 1.d0
   else
      large_error = .false.
   end if
   return
end function large_error
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine is called before the sanity check, and updates the diagnostic vari-  !
! ables, namely the temperature and liquid fraction of leaf water, soil layers and         !
! temporary snow/pond layers.                                                                      !
!------------------------------------------------------------------------------------------!
subroutine update_diagnostic_vars(initp, csite,ipa)
   use rk4_coms             , only : rk4met               & ! intent(in)
                                   , rk4min_sfcwater_mass & ! intent(in)
                                   , rk4patchtype         ! ! structure
   use ed_state_vars        , only : sitetype             & ! structure
                                   , patchtype            ! ! structure
   use soil_coms            , only : soil8                ! ! intent(in)
   use grid_coms            , only : nzg                  & ! intent(in)
                                   , nzs                  ! ! intent(in)
   use therm_lib            , only : qwtk8                & ! subroutine
                                   , qtk8                 ! ! subroutine
   use consts_coms          , only : wdns8                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)        , pointer :: cpatch
   integer                          :: ico
   integer                          :: k
   real(kind=8)                     :: soilhcap
   !---------------------------------------------------------------------------------------!


   !----- Updating soil temperature and liquid water fraction. ----------------------------!
   do k = rk4met%lsl, nzg - 1
      soilhcap = soil8(csite%ntext_soil(k,ipa))%slcpd
      call qwtk8(initp%soil_energy(k),initp%soil_water(k)*wdns8,soilhcap                   &
                ,initp%soil_tempk(k),initp%soil_fracliq(k))
   end do

   !---------------------------------------------------------------------------------------!
   !    Updating surface water temperature and liquid water fraction, remembering that in- !
   ! side the RK4 integration, surface water energy is in J/m². The abs is necessary be-   !
   ! cause surface mass may indeed become too negative during the integration process and  !
   ! if it happens, we want the step to be rejected.                                       !
   !---------------------------------------------------------------------------------------!
   do k = 1, nzs
      if(abs(initp%sfcwater_mass(k)) > rk4min_sfcwater_mass)  then
           call qtk8(initp%sfcwater_energy(k)/initp%sfcwater_mass(k)                       &
                    ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
      elseif (k == 1) then
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
         initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
      else
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
         initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
      end if
   end do


   cpatch => csite%patch(ipa)

   !----- Looping over cohorts ------------------------------------------------------------!
   cohortloop: do ico=1,cpatch%ncohorts
      !----- Checking whether this is a prognostic cohort... ------------------------------!
      if (initp%solvable(ico)) then
         !----- Lastly we update leaf temperature and liquid fraction. --------------------!
         call qwtk8(initp%veg_energy(ico),initp%veg_water(ico),initp%hcapveg(ico)          &
                   ,initp%veg_temp(ico),initp%veg_fliq(ico))
      end if

   end do cohortloop

   return
end subroutine update_diagnostic_vars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine performs the following tasks:                                         !
! 1. Check how many layers of temporary water or snow we have, and include the virtual     !
!    pools at the topmost if needed;                                                       !
! 2. Force thermal equilibrium between topmost soil layer and a single snow/water layer    !
!    if the layer is too thin;                                                             !
! 3. Compute the amount of mass each layer has, and redistribute them accordingly.         !
! 4. Percolates excessive liquid water if needed.                                          !
!------------------------------------------------------------------------------------------!
subroutine redistribute_snow(initp,csite,ipa)

   use rk4_coms      , only : rk4patchtype         & ! structure
                            , rk4min_sfcw_mass     & ! intent(in)
                            , rk4min_virt_water    & ! intent(in)
                            , rk4water_stab_thresh & ! intent(in)
                            , rk4min_sfcwater_mass & ! intent(in)
                            , rk4snowmin           ! ! intent(in)
   use ed_state_vars , only : sitetype             & ! structure
                            , patchtype            ! ! structure
   use grid_coms     , only : nzs                  & ! intent(in)
                            , nzg                  ! ! intent(in)
   use soil_coms     , only : soil8                & ! intent(in)
                            , dslz8                & ! intent(in)
                            , dslzi8               & ! intent(in)
                            , thick                & ! intent(in)
                            , thicknet             ! ! intent(in)
   use consts_coms   , only : cice8                & ! intent(in)
                            , cliq8                & ! intent(in)
                            , t3ple8               & ! intent(in)
                            , wdns8                & ! intent(in)
                            , tsupercool8          & ! intent(in)
                            , qliqt38              & ! intent(in)
                            , wdnsi8               ! ! intent(in)
   use therm_lib     , only : qtk8                 & ! subroutine
                            , qwtk8                ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp
   type(sitetype)         , target     :: csite
   integer                , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: kold
   integer                             :: newlayers
   integer                             :: nlayers
   integer                             :: ksn
   integer                             :: ksnnew
   integer                             :: k
   !----- Control variables ---------------------------------------------------------------!
   real(kind=8)                        :: wtold
   real(kind=8)                        :: wtnew
   real(kind=8), dimension(nzs)        :: newsfcw_mass
   real(kind=8), dimension(nzs)        :: newsfcw_energy
   real(kind=8), dimension(nzs)        :: newsfcw_depth
   real(kind=8)                        :: wdiff
   real(kind=8)                        :: totsnow
   real(kind=8)                        :: depthgain
   real(kind=8)                        :: wfree
   real(kind=8)                        :: qwfree
   real(kind=8)                        :: qw
   real(kind=8)                        :: w
   real(kind=8)                        :: wfreeb
   real(kind=8)                        :: depthloss
   real(kind=8)                        :: snden
   real(kind=8)                        :: sndenmin
   real(kind=8)                        :: sndenmax
   real(kind=8)                        :: qwt
   real(kind=8)                        :: wt
   real(kind=8)                        :: soilhcap
   real(kind=8)                        :: free_surface_water_demand
   integer                             :: nsoil
   !----- Constants -----------------------------------------------------------------------!
   logical                , parameter  :: debug = .false.
   !---------------------------------------------------------------------------------------!


   !----- Initializing # of layers alias --------------------------------------------------!
   ksn       = initp%nlev_sfcwater

   if (ksn >= 1) then
      !------------------------------------------------------------------------------------!
      ! 1. There used to exist temporary water/snow layers here.  Check total mass to see  !
      !    whether there is still enough mass.                                             !
      !------------------------------------------------------------------------------------!
      totsnow = sum(initp%sfcwater_mass(1:ksn))
      if (totsnow < rk4min_sfcw_mass) then
         !----- Temporary layer is too negative, break it so the step can be rejected. ----!
         return
      elseif (totsnow <= rk4min_sfcwater_mass) then
         !---------------------------------------------------------------------------------!
         ! 1.a. Too little or negative mass.  Eliminate layers, ensuring that it will  not !
         !      leak mass or energy, by "stealing" them from the top soil  !
         !      layer.                                                                     !
         !---------------------------------------------------------------------------------!
         initp%sfcwater_energy(1) = sum(initp%sfcwater_energy(1:ksn))
         initp%sfcwater_mass(1)   = sum(initp%sfcwater_mass(1:ksn))
         initp%soil_energy(nzg)   = initp%soil_energy(nzg)                                 &
                                  + initp%sfcwater_energy(1) * dslzi8(nzg)
         initp%soil_water(nzg)    = initp%soil_water(nzg)                                  &
                                  + initp%sfcwater_mass(1)   * dslzi8(nzg)
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                     &
                   ,soil8(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)           &
                   ,initp%soil_fracliq(nzg))
         initp%sfcwater_mass      = 0.d0
         initp%sfcwater_energy    = 0.d0
         initp%sfcwater_tempk     = initp%soil_tempk(nzg)
         initp%sfcwater_fracliq   = 0.d0
         initp%sfcwater_depth     = 0.d0       
         initp%nlev_sfcwater      = 0
         ksnnew                   = 0
      else
         !---------------------------------------------------------------------------------!
         ! 1.b.  Still something there, nothing changes at least not for the time being.   !
         !---------------------------------------------------------------------------------!
         ksnnew = ksn
         wfree               = 0.d0
         qwfree              = 0.d0
         depthgain           = 0.d0
      end if
   else
      !------------------------------------------------------------------------------------!
      ! 2.  No temporary layer, dealing with virtual layer.  Check whether the virtual     !
      !     layer would be thick enough to create a pond, otherwise skip the entire thing. !
      !------------------------------------------------------------------------------------!
      if (initp%virtual_water < rk4min_virt_water) then
         !----- Virtual layer is too negative, break it so the step can be rejected. ------!
         return
      elseif (initp%virtual_water <= rk4min_sfcwater_mass) then
         !---------------------------------------------------------------------------------!
         ! 2.a. Too little or negative mass in the virtual layer.  No layer will be creat- !
         !      ed, but before eliminating it, just make sure mass and energy will be      !
         !      conserved.                                                                 !
         !---------------------------------------------------------------------------------!
         ksnnew = 0
         initp%soil_energy(nzg)   = initp%soil_energy(nzg)                                 &
                                  + initp%virtual_heat * dslzi8(nzg)
         initp%soil_water(nzg)    = initp%soil_water(nzg)                                  &
                                  + initp%virtual_water * dslzi8(nzg)
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                     &
                   ,soil8(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)           &
                   ,initp%soil_fracliq(nzg))
         initp%virtual_water      = 0.d0
         initp%virtual_heat       = 0.d0
         initp%virtual_depth      = 0.d0
      else
         !---------------------------------------------------------------------------------!
         ! 2.b. No temporary layer, significant mass to add.  ksnnew will be at least one. !
         !      If there was no layer before, create one.                                  !
         !---------------------------------------------------------------------------------!
         wfree               = initp%virtual_water
         qwfree              = initp%virtual_heat
         depthgain           = initp%virtual_depth
         initp%virtual_water = 0.d0
         initp%virtual_heat  = 0.d0
         initp%virtual_depth = 0.d0
         ksnnew = 1
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! 3. We now update the diagnostic variables, and ensure the layers are stable.  Loop    !
   !    over layers, from top to bottom this time.                                         !
   !---------------------------------------------------------------------------------------!
   totsnow =0.d0
   do k = ksnnew,1,-1

      !----- Update current mass and energy of temporary layer ----------------------------!
      qw = initp%sfcwater_energy(k) + qwfree
      w = initp%sfcwater_mass(k) + wfree

      !------------------------------------------------------------------------------------!
      !    Single layer, and this is a very thin one, which can cause numerical instabili- !
      ! ty. Force a fast heat exchange between this thin layer and the soil topmost level, !
      ! bringing both layers to a thermal equilibrium.                                     !
      !------------------------------------------------------------------------------------!
      if (ksnnew == 1 .and. initp%sfcwater_mass(k) < rk4water_stab_thresh) then

         !---------------------------------------------------------------------------------!
         !     Total internal energy and water of the combined system, in J/m² and kg/m²,  !
         ! respectively.                                                                   !
         !---------------------------------------------------------------------------------!
         qwt = qw + initp%soil_energy(nzg) * dslz8(nzg)
         wt  = w  + initp%soil_water(nzg)  * dslz8(nzg) * wdns8

         !----- Finding the equilibrium temperature and liquid/ice partition. -------------!
         soilhcap = soil8(csite%ntext_soil(nzg,ipa))%slcpd * dslz8(nzg)
         call qwtk8(qwt,wt,soilhcap,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))

         !---------------------------------------------------------------------------------!
         !    Computing internal energy of the temporary layer with the temperature and    !
         ! liquid/ice distribution we just found, for the mass the layer has.              !
         !---------------------------------------------------------------------------------!
         qw = w                                                                            &
            * (initp%sfcwater_fracliq(k) * cliq8 *(initp%sfcwater_tempk(k)-tsupercool8)    &
                  + (1.d0-initp%sfcwater_fracliq(k)) * cice8 * initp%sfcwater_tempk(k))

         !---------------------------------------------------------------------------------!
         !    Set the properties of top soil layer. Since internal energy is an extensive  !
         ! quantity, we can simply take the difference to be the soil internal energy,     !
         ! just remembering that we need to convert it back to J/m³. The other properties  !
         ! can be copied from the surface layer because we assumed phase and temperature   !
         ! equilibrium.                                                                    !
         !---------------------------------------------------------------------------------!
         initp%soil_energy(nzg)  = (qwt - qw) * dslzi8(nzg)
         initp%soil_tempk(nzg)   = initp%sfcwater_tempk(k)
         initp%soil_fracliq(nzg) = initp%sfcwater_fracliq(k)
      else
         !----- Layer is computationally stable, just update the temperature and phase ----!
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                     &
                   ,soil8(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)           &
                   ,initp%soil_fracliq(nzg))
         call qtk8(qw/w,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
      end if


      !------------------------------------------------------------------------------------!
      !    Shed liquid in excess of a 1:9 liquid-to-ice ratio through percolation.  Limit  !
      ! this shed amount (wfreeb) in lowest snow layer to amount top soil layer can hold.  !
      !------------------------------------------------------------------------------------!
      if (w > rk4min_sfcwater_mass) then
         wfreeb = max(0.d0, w * (initp%sfcwater_fracliq(k)-1.d-1)/9.d-1)
      else
         wfreeb = 0.0
      end if

      if (k == 1)then
           !----- Do "greedy" infiltration. -----------------------------------------------!
           nsoil = csite%ntext_soil(nzg,ipa)
           free_surface_water_demand = max(0.d0                                            &
                                          ,soil8(nsoil)%slmsts - initp%soil_water(nzg))    &
                                     * wdns8 * dslz8(nzg)
           wfreeb = min(wfreeb,free_surface_water_demand)
           qwfree = wfreeb * cliq8 * (initp%sfcwater_tempk(k)-tsupercool8)
           !----- Update topmost soil moisture and energy, updating temperature and phase -!
           initp%soil_water(nzg)  = initp%soil_water(nzg)  + wfreeb*wdnsi8*dslzi8(nzg)
           initp%soil_energy(nzg) = initp%soil_energy(nzg) + qwfree * dslzi8(nzg)
           soilhcap = soil8(nsoil)%slcpd
           call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*wdns8                   &
                     ,soilhcap,initp%soil_tempk(nzg),initp%soil_fracliq(nzg))
      else
         !---- Not the first layer, just shed all free water, and compute its energy ------!
         qwfree = wfreeb * cliq8 * (initp%sfcwater_tempk(k)-tsupercool8)
      end if
      depthloss = wfreeb * wdnsi8
      
      !----- Remove water and internal energy losses due to percolation -------------------!
      initp%sfcwater_mass(k)  = w - wfreeb
      initp%sfcwater_depth(k) = initp%sfcwater_depth(k) + depthgain - depthloss
      if(initp%sfcwater_mass(k) > rk4min_sfcwater_mass) then
         initp%sfcwater_energy(k) = qw - qwfree
         call qtk8(initp%sfcwater_energy(k)/initp%sfcwater_mass(k),initp%sfcwater_tempk(k) &
                 ,initp%sfcwater_fracliq(k))
      else
         initp%sfcwater_energy(k) = 0.d0
         initp%sfcwater_mass(k)   = 0.d0
         initp%sfcwater_depth(k)  = 0.d0
         if (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end if

      !----- Integrate total "snow" -------------------------------------------------------!
      totsnow = totsnow + initp%sfcwater_mass(k)

      !----- Calculate density and depth of snow ------------------------------------------!
      snden    = initp%sfcwater_mass(k) / max(1.0d-6,initp%sfcwater_depth(k))
      sndenmax = wdns8
      sndenmin = max(3.d1, 2.d2 * (wfree + wfreeb)                                         &
               / max(rk4min_sfcwater_mass,initp%sfcwater_mass(k)))
      snden    = min(sndenmax, max(sndenmin,snden))
      initp%sfcwater_depth(k) = initp%sfcwater_mass(k) / snden

      !----- Set up input to next layer ---------------------------------------------------!
      wfree = wfreeb
      depthgain = depthloss
   end do

   !---------------------------------------------------------------------------------------!
   ! 4. Re-distribute snow layers to maintain prescribed distribution of mass.             !
   !---------------------------------------------------------------------------------------!
   if (totsnow <= rk4min_sfcwater_mass .or. ksnnew == 0) then
      initp%nlev_sfcwater = 0
      !----- Making sure that the unused layers have zero in everything -------------------!
      do k = 1, nzs
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         if (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do
   else
      !---- Check whether there is enough snow for a new layer. ---------------------------!
      nlayers   = ksnnew
      newlayers = 1
      do k = 1,nzs
         !----- Checking whether we need 
         if (      initp%sfcwater_mass(k)   >  rk4min_sfcwater_mass                        &
             .and. rk4snowmin * thicknet(k) <= totsnow                                     &
             .and. initp%sfcwater_energy(k) <  initp%sfcwater_mass(k)*qliqt38 ) then

            newlayers = newlayers + 1
         end if
      end do
      newlayers = min(newlayers, nzs, nlayers + 1)
      initp%nlev_sfcwater = newlayers
      kold  = 1
      wtnew = 1.d0
      wtold = 1.d0
      do k = 1,newlayers
         newsfcw_mass(k)   = totsnow * thick(k,newlayers)
         newsfcw_energy(k) = 0.d0
         newsfcw_depth(k)  = 0.d0
         !----- Finding new layer properties ----------------------------------------------!
         find_layer: do

            !----- Difference between old and new snow ------------------------------------!
            wdiff = wtnew * newsfcw_mass(k) - wtold * initp%sfcwater_mass(kold)  

            if (wdiff > 0.d0) then
               newsfcw_energy(k) = newsfcw_energy(k) + wtold * initp%sfcwater_energy(kold)
               newsfcw_depth(k)  = newsfcw_depth(k)  + wtold * initp%sfcwater_depth(kold)
               wtnew  = wtnew - wtold * initp%sfcwater_mass(kold) / newsfcw_mass(k)
               kold   = kold + 1
               wtold  = 1.0
               if (kold > nlayers) exit find_layer
            else
               newsfcw_energy(k) = newsfcw_energy(k) + wtnew * newsfcw_mass(k)             &
                                 * initp%sfcwater_energy(kold)                             &
                                 / max(rk4min_sfcwater_mass,initp%sfcwater_mass(kold))
               newsfcw_depth(k)  = newsfcw_depth(k)  + wtnew * newsfcw_mass(k)             &
                                 * initp%sfcwater_depth(kold)                              &
                                 / max(rk4min_sfcwater_mass,initp%sfcwater_mass(kold))
               wtold = wtold - wtnew * newsfcw_mass(k)                                     &
                             / max(rk4min_sfcwater_mass,initp%sfcwater_mass(kold))
               wtnew = 1.
               exit find_layer
            end if
         end do find_layer
      end do

      !----- Updating the water/snow layer properties -------------------------------------!
      do k = 1,newlayers
         initp%sfcwater_mass(k)   = newsfcw_mass(k)
         initp%sfcwater_energy(k) = newsfcw_energy(k)
         initp%sfcwater_depth(k)  = newsfcw_depth(k)
         if (newsfcw_mass(k) > rk4min_sfcwater_mass) then
            call qtk8(initp%sfcwater_energy(k)/initp%sfcwater_mass(k)                      &
                    ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
         elseif (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do

      !----- Making sure that the unused layers have zero in everything -------------------!
      do k = newlayers + 1, nzs
         initp%sfcwater_mass(k)    = 0.d0
         initp%sfcwater_energy(k)  = 0.d0
         initp%sfcwater_depth(k)   = 0.d0
         if (k == 1) then
            initp%sfcwater_tempk(k)   = initp%soil_tempk(nzg)
            initp%sfcwater_fracliq(k) = initp%soil_fracliq(nzg)
         else
            initp%sfcwater_tempk(k)   = initp%sfcwater_tempk(k-1)
            initp%sfcwater_fracliq(k) = initp%sfcwater_fracliq(k-1)
         end if
      end do
   end if

   return
end subroutine redistribute_snow
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will ensure that leaf water is positively defined.  Depending on its  !
! derivative, it can go under zero, in which case we must correct the derivatives rather   !
! than forcing it to be zero.  This guarantees mass conservation.  Likewise, if in the end !
! of the step the leaf water is over the maximum, we remove the excess through shedding.   !
!    After this is checked, we then update the remaining leaf properties, namely the       !
! temperature and liquid water fraction.                                                   !
!------------------------------------------------------------------------------------------!
subroutine adjust_veg_properties(initp,hdid,csite,ipa)
   use rk4_coms             , only : rk4patchtype      & ! structure
                                   , rk4met            & ! intent(in)
                                   , rk4eps            & ! intent(in)
                                   , rk4min_veg_lwater & ! intent(in)
                                   , wcapcani          & ! intent(in)
                                   , rk4dry_veg_lwater & ! intent(in)
                                   , rk4fullveg_lwater ! ! intent(in)
   use ed_state_vars        , only : sitetype          & ! structure
                                   , patchtype         ! ! structure
   use consts_coms          , only : cice8             & ! intent(in)
                                   , cliq8             & ! intent(in)
                                   , alvl8             & ! intent(in)
                                   , alvi8             & ! intent(in)
                                   , t3ple8            & ! intent(in)
                                   , wdns8             & ! intent(in)
                                   , idns8             & ! intent(in)
                                   , tsupercool8       & ! intent(in)
                                   , qliqt38           & ! intent(in)
                                   , wdnsi8            ! ! intent(in)
   use therm_lib            , only : qwtk8             ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp  ! Integration buffer
   type(sitetype)         , target     :: csite  ! Current site
   integer                , intent(in) :: ipa    ! Current patch ID
   real(kind=8)           , intent(in) :: hdid   ! Time step 
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)        , pointer    :: cpatch
   integer                             :: ico
   integer                             :: ksn
   real(kind=8)                        :: rk4min_leaf_water
   real(kind=8)                        :: min_leaf_water
   real(kind=8)                        :: max_leaf_water
   real(kind=8)                        :: veg_wshed
   real(kind=8)                        :: veg_qwshed
   real(kind=8)                        :: veg_dwshed
   real(kind=8)                        :: veg_dew
   real(kind=8)                        :: veg_qdew
   real(kind=8)                        :: hdidi
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)
   
   !----- Inverse of time increment -------------------------------------------------------!
   hdidi = 1.d0 / hdid

   !----- Looping over cohorts ------------------------------------------------------------!
   cohortloop: do ico=1,cpatch%ncohorts
      !----- Checking whether this is a prognostic cohort... ------------------------------!
      if (initp%solvable(ico)) then
         !---------------------------------------------------------------------------------!
         !   Now we find the maximum leaf water possible.                                  !
         !---------------------------------------------------------------------------------!
         rk4min_leaf_water = rk4min_veg_lwater * initp%tai(ico)
         min_leaf_water    = rk4dry_veg_lwater * initp%tai(ico)
         max_leaf_water    = rk4fullveg_lwater * initp%tai(ico)

         !------ Leaf water is too negative, break it so the step can be rejected. --------!
         if (initp%veg_water(ico) < rk4min_leaf_water) then
            return
         !----- Shedding excessive water to the ground ------------------------------------!
         elseif (initp%veg_water(ico) > max_leaf_water) then
            veg_wshed  = (initp%veg_water(ico)-max_leaf_water)
            veg_qwshed = veg_wshed                                                         &
                       * (initp%veg_fliq(ico) * cliq8 * (initp%veg_temp(ico)-tsupercool8)    &
                         + (1.d0-initp%veg_fliq(ico)) * cice8 * initp%veg_temp(ico))
            veg_dwshed = veg_wshed                                                         &
                       / (initp%veg_fliq(ico) * wdns8 + (1.d0-initp%veg_fliq(ico))*idns8)

            !----- Updating water mass and energy. ----------------------------------------!
            initp%veg_water(ico)  = initp%veg_water(ico)  - veg_wshed
            initp%veg_energy(ico) = initp%veg_energy(ico) - veg_qwshed
            
            !----- Updating virtual pool --------------------------------------------------!
            ksn = initp%nlev_sfcwater
            if (ksn > 0) then
               initp%sfcwater_mass(ksn)   = initp%sfcwater_mass(ksn)   + veg_wshed
               initp%sfcwater_energy(ksn) = initp%sfcwater_energy(ksn) + veg_qwshed
               initp%sfcwater_depth(ksn)  = initp%sfcwater_depth(ksn)  + veg_dwshed
            else
               initp%virtual_water   = initp%virtual_water + veg_wshed
               initp%virtual_heat    = initp%virtual_heat  + veg_qwshed
               initp%virtual_depth   = initp%virtual_depth + veg_dwshed
            end if
            !----- Updating output fluxes -------------------------------------------------!
            initp%avg_wshed_vg  = initp%avg_wshed_vg  + veg_wshed  * hdidi
            initp%avg_qwshed_vg = initp%avg_qwshed_vg + veg_qwshed * hdidi

         !---------------------------------------------------------------------------------!
         !    If veg_water is tiny or negative, exchange moisture with the air, "stealing" !
         ! moisture as fast "dew/frost" condensation if it is negative, or "donating" the  !
         ! remaining as "boiling" (fast evaporation).                                      !
         !---------------------------------------------------------------------------------!
         elseif (initp%veg_water(ico) < min_leaf_water) then
            veg_dew = - initp%veg_water(ico)
            if (initp%can_temp >= t3ple8) then
               veg_qdew = veg_dew * alvl8
            else
               veg_qdew = veg_dew * alvi8
            end if

            !----- Updating state variables -----------------------------------------------!
            initp%veg_water(ico)  = 0.d0
            initp%veg_energy(ico) = initp%veg_energy(ico)  + veg_qdew
            initp%can_shv         = initp%can_shv          - veg_dew * wcapcani

            !----- Updating output flux ---------------------------------------------------!
            initp%avg_vapor_vc    = initp%avg_vapor_vc - veg_dew * hdidi
         end if

         !----- Lastly we update leaf temperature and liquid fraction. --------------------!
         call qwtk8(initp%veg_energy(ico),initp%veg_water(ico),initp%hcapveg(ico)           &
                  ,initp%veg_temp(ico),initp%veg_fliq(ico))
      end if

   end do cohortloop

   return
end subroutine adjust_veg_properties
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies the values to different buffers inside the RK4 integration     !
! scheme.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine copy_rk4_patch(sourcep, targetp, cpatch)

   use rk4_coms      , only : rk4met            & ! intent(in)
                            , rk4patchtype      ! ! structure
   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         ! ! structure
   use grid_coms     , only : nzg               & ! intent(in)
                            , nzs               ! ! intent(in)
   use ed_max_dims      , only : n_pft             ! ! intent(in)
   use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: sourcep
   type(rk4patchtype) , target     :: targetp
   type(patchtype)    , target     :: cpatch
   !----- Local variable ------------------------------------------------------------------!
   integer                         :: k
   !---------------------------------------------------------------------------------------!

   targetp%can_temp      = sourcep%can_temp
   targetp%can_shv       = sourcep%can_shv
   targetp%can_co2       = sourcep%can_co2

   targetp%virtual_water = sourcep%virtual_water
   targetp%virtual_heat  = sourcep%virtual_heat
   targetp%virtual_depth = sourcep%virtual_depth

   targetp%rough         = sourcep%rough
 
   targetp%upwp          = sourcep%upwp
   targetp%wpwp          = sourcep%wpwp
   targetp%tpwp          = sourcep%tpwp
   targetp%qpwp          = sourcep%qpwp

   targetp%ground_shv    = sourcep%ground_shv
   targetp%surface_ssh   = sourcep%surface_ssh
   targetp%surface_temp  = sourcep%surface_temp
   targetp%surface_fliq  = sourcep%surface_fliq

   targetp%nlev_sfcwater = sourcep%nlev_sfcwater
   targetp%ustar         = sourcep%ustar
   targetp%cstar         = sourcep%cstar
   targetp%tstar         = sourcep%tstar
   targetp%qstar         = sourcep%qstar
   targetp%virtual_flag  = sourcep%virtual_flag
   targetp%rasveg        = sourcep%rasveg

   do k=rk4met%lsl,nzg
      
      targetp%soil_water(k)             = sourcep%soil_water(k)
      targetp%soil_energy(k)            = sourcep%soil_energy(k)
      targetp%soil_tempk(k)             = sourcep%soil_tempk(k)
      targetp%soil_fracliq(k)           = sourcep%soil_fracliq(k)
      targetp%available_liquid_water(k) = sourcep%available_liquid_water(k)
      targetp%extracted_water(k)        = sourcep%extracted_water(k)
      targetp%psiplusz(k)               = sourcep%psiplusz(k)
      targetp%soilair99(k)              = sourcep%soilair99(k)
      targetp%soilair01(k)              = sourcep%soilair01(k)
      targetp%soil_liq(k)               = sourcep%soil_liq(k)
   end do

   do k=1,nzs
      targetp%sfcwater_mass(k)    = sourcep%sfcwater_mass(k)   
      targetp%sfcwater_energy(k)  = sourcep%sfcwater_energy(k) 
      targetp%sfcwater_depth(k)   = sourcep%sfcwater_depth(k)  
      targetp%sfcwater_tempk(k)   = sourcep%sfcwater_tempk(k)  
      targetp%sfcwater_fracliq(k) = sourcep%sfcwater_fracliq(k)
   end do

   do k=1,cpatch%ncohorts
      targetp%veg_water(k)   = sourcep%veg_water(k)
      targetp%veg_energy(k)  = sourcep%veg_energy(k)
      targetp%veg_temp(k)    = sourcep%veg_temp(k)
      targetp%veg_fliq(k)    = sourcep%veg_fliq(k)
      targetp%hcapveg(k)     = sourcep%hcapveg(k)
      targetp%lai(k)         = sourcep%lai(k)
      targetp%wpa(k)         = sourcep%wpa(k)
      targetp%tai(k)         = sourcep%tai(k)
      targetp%solvable(k)    = sourcep%solvable(k)
   end do

   if (fast_diagnostics) then
      targetp%wbudget_loss2atm   = sourcep%wbudget_loss2atm
      targetp%co2budget_loss2atm = sourcep%co2budget_loss2atm
      targetp%ebudget_loss2atm   = sourcep%ebudget_loss2atm
      targetp%ebudget_latent     = sourcep%ebudget_latent
      targetp%avg_carbon_ac      = sourcep%avg_carbon_ac
      targetp%avg_vapor_vc       = sourcep%avg_vapor_vc
      targetp%avg_dew_cg         = sourcep%avg_dew_cg  
      targetp%avg_vapor_gc       = sourcep%avg_vapor_gc
      targetp%avg_wshed_vg       = sourcep%avg_wshed_vg
      targetp%avg_vapor_ac       = sourcep%avg_vapor_ac
      targetp%avg_transp         = sourcep%avg_transp  
      targetp%avg_evap           = sourcep%avg_evap   
      targetp%avg_netrad         = sourcep%avg_netrad   
      targetp%avg_sensible_vc    = sourcep%avg_sensible_vc  
      targetp%avg_sensible_2cas  = sourcep%avg_sensible_2cas
      targetp%avg_qwshed_vg      = sourcep%avg_qwshed_vg    
      targetp%avg_sensible_gc    = sourcep%avg_sensible_gc  
      targetp%avg_sensible_ac    = sourcep%avg_sensible_ac  
      targetp%avg_sensible_tot   = sourcep%avg_sensible_tot 

      do k=rk4met%lsl,nzg
         targetp%avg_sensible_gg(k) = sourcep%avg_sensible_gg(k)
         targetp%avg_smoist_gg(k)   = sourcep%avg_smoist_gg(k)  
         targetp%avg_smoist_gc(k)   = sourcep%avg_smoist_gc(k)  
      end do
   end if



   return
end subroutine copy_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine prints the patch and cohort information when the model falls apart... !
!------------------------------------------------------------------------------------------!
subroutine print_csiteipa(csite, ipa)
   use rk4_coms              , only : rk4met        ! ! intent(in)
   use ed_state_vars         , only : sitetype      & ! structure
                                    , patchtype     ! ! structure
   use ed_misc_coms             , only : current_time  ! ! intent(in)
   use grid_coms             , only : nzs           & ! intent(in)
                                    , nzg           ! ! intent(in)
   use ed_max_dims              , only : n_pft         ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: ipa
   !----- Local variable ------------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   integer                      :: ico 
   integer                      :: k   
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (csite) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,f12.0,1x,a)')                                &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%time,'UTC'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
 
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Cohort information (only the solvable ones shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),11(a12,1x))')                                              &
         '    PFT','KRDEPTH','      NPLANT','         LAI','         DBH','       BDEAD'   &
                           &,'      BALIVE','  VEG_ENERGY','    VEG_TEMP','   VEG_WATER'   &
                           &,'     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if (cpatch%solvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),11(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico) &
              ,cpatch%nplant(ico),cpatch%lai(ico),cpatch%dbh(ico),cpatch%bdead(ico)        &
              ,cpatch%balive(ico),cpatch%veg_energy(ico),cpatch%veg_temp(ico)              &
              ,cpatch%veg_water(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)
      end if
   end do
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '   DIST_TYPE','         AGE','        AREA'          &
                                   &,'          RH','AVGDAILY_TMP','     SUM_CHD'          &
                                   &,'     SUM_DGD'
   write (unit=*,fmt='(i12,1x,6(es12.4,1x))')  csite%dist_type(ipa),csite%age(ipa)         &
         ,csite%area(ipa),csite%rh(ipa),csite%avg_daily_temp(ipa),csite%sum_chd(ipa)       &
         ,csite%sum_dgd(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','         LAI'          &
                                   &,'        HTRY','     CAN_CO2','    CAN_TEMP'          &
                                   &,'     CAN_SHV'
   write (unit=*,fmt='(7(es12.4,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)          &
         ,csite%lai(ipa),csite%htry(ipa),csite%can_co2(ipa),csite%can_temp(ipa)            &
         ,csite%can_shv(ipa) 

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                   &,'       TSTAR','     RLONG_G','    RSHORT_G'          &
                                   &,'     RLONG_S'
   write (unit=*,fmt='(7(es12.4,1x))') csite%ustar(ipa),csite%qstar(ipa),csite%cstar(ipa)  &
         ,csite%tstar(ipa),csite%rlong_g(ipa),csite%rshort_g(ipa),csite%rlong_s(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,a12)') '  PFT','       REPRO'
   do k=1,n_pft
      write (unit=*,fmt='(i5,1x,es12.4)') k,csite%repro(k,ipa)
   end do

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=rk4met%lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,csite%ntext_soil(k,ipa)              &
            ,csite%soil_energy(k,ipa),csite%soil_tempk(k,ipa),csite%soil_water(k,ipa)      &
            ,csite%soil_fracliq(k,ipa)
   end do
   
   if (csite%nlev_sfcwater(ipa) >= 1) then
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a5,1x,6(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'       &
                                      &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'       &
                                      &,'    RSHORT_S'
      do k=1,csite%nlev_sfcwater(ipa)
         write (unit=*,fmt='(i5,1x,6(es12.4,1x))') k,csite%sfcwater_energy(k,ipa)          &
               ,csite%sfcwater_tempk(k,ipa),csite%sfcwater_mass(k,ipa)                     &
               ,csite%sfcwater_fracliq(k,ipa),csite%sfcwater_depth(k,ipa)                  &
               ,csite%rshort_s(k,ipa)
      end do
   end if

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '
   return
end subroutine print_csiteipa
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is similar to print_csiteipa, except that it also prints the       !
! outcome of the Runge-Kutta integrator.                                                   !
!------------------------------------------------------------------------------------------!
subroutine print_rk4patch(y,csite,ipa)
   use rk4_coms              , only : rk4patchtype         & ! structure
                                    , rk4met               & ! intent(in)
                                    , rk4min_sfcwater_mass ! ! intent(in)
   use ed_state_vars         , only : sitetype             & ! structure
                                    , patchtype            ! ! structure
   use grid_coms             , only : nzg                  & ! intent(in)
                                    , nzs                  ! ! intent(in)
   use ed_misc_coms             , only : current_time         ! ! intent(in)
   use therm_lib             , only : qtk8                 & ! subroutine
                                    , qwtk8                ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: y
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer    :: cpatch
   integer                         :: k
   integer                         :: ico
   real(kind=8)                    :: virtual_temp, virtual_fliq
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (rk4patch) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,f12.0,1x,a)')                                &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%time,'s'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(80a)')         ('-',k=1,80)
   write (unit=*,fmt='(a)')           ' ATMOSPHERIC CONDITIONS: '
   write (unit=*,fmt='(a,1x,es12.4)') ' Air temperature    : ',rk4met%atm_tmp
   write (unit=*,fmt='(a,1x,es12.4)') ' H2Ov mixing ratio  : ',rk4met%atm_shv
   write (unit=*,fmt='(a,1x,es12.4)') ' CO2  mixing ratio  : ',rk4met%atm_co2
   write (unit=*,fmt='(a,1x,es12.4)') ' Pressure           : ',rk4met%prss
   write (unit=*,fmt='(a,1x,es12.4)') ' Exner function     : ',rk4met%exner
   write (unit=*,fmt='(a,1x,es12.4)') ' Air density        : ',rk4met%rhos
   write (unit=*,fmt='(a,1x,es12.4)') ' Wind speed         : ',rk4met%vels
   write (unit=*,fmt='(a,1x,es12.4)') ' Height             : ',rk4met%geoht
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. mass  flux : ',rk4met%pcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. heat  flux : ',rk4met%qpcpg
   write (unit=*,fmt='(a,1x,es12.4)') ' Precip. depth flux : ',rk4met%dpcpg

   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Cohort information (only those solvable are shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','        HITE','         DBH','       BDEAD'   &
                           &,'      BALIVE','     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if (cpatch%solvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),8(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdead(ico)       &
              ,cpatch%balive(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','         LAI','         WPA','         TAI','  VEG_ENERGY'   &
             ,'   VEG_WATER','    VEG_HCAP','    VEG_TEMP','    VEG_FLIQ'
   do ico = 1,cpatch%ncohorts
      if (y%solvable(ico)) then
         write(unit=*,fmt='(2(i7,1x),9(es12.4,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
               ,y%lai(ico),y%wpa(ico),y%tai(ico),y%veg_energy(ico),y%veg_water(ico)        &
               ,y%hcapveg(ico),y%veg_temp(ico),y%veg_fliq(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(6(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','   PATCH_LAI'          &
                                   &,'     CAN_CO2','    CAN_TEMP','     CAN_SHV'
   write (unit=*,fmt='(6(es12.4,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)          &
         ,csite%lai(ipa),y%can_co2,y%can_temp,y%can_shv

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(4(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                   &,'       TSTAR'
   write (unit=*,fmt='(4(es12.4,1x))') y%ustar,y%qstar,y%cstar,y%tstar

   write (unit=*,fmt='(80a)') ('-',k=1,80)
   if (y%virtual_water /= 0.) then
      call qtk8(y%virtual_heat/y%virtual_water,virtual_temp,virtual_fliq)
   else
      virtual_temp = y%soil_tempk(nzg)
      virtual_fliq = y%soil_fracliq(nzg)
   end if


   write (unit=*,fmt='(5(a12,1x))')  'VIRTUAL_FLAG','VIRTUAL_HEAT','  VIRT_WATER'          &
                                   &,'VIRTUAL_TEMP','VIRTUAL_FLIQ'
   write (unit=*,fmt='(i12,1x,4(es12.4,1x))') y%virtual_flag,y%virtual_heat                &
                                             ,y%virtual_water,virtual_temp,virtual_fliq
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(4(a12,1x))')    '  GROUND_SHV',' SURFACE_SSH','SURFACE_TEMP'        &
                                      ,'SURFACE_FLIQ'
   write (unit=*,fmt='(4(es12.4,1x))') y%ground_shv, y%surface_ssh, y%surface_temp         &
                                      ,y%surface_fliq

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=rk4met%lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.4,1x))') k,csite%ntext_soil(k,ipa)              &
            ,y%soil_energy(k),y%soil_tempk(k),y%soil_water(k),y%soil_fracliq(k)
   end do
   
   if (csite%nlev_sfcwater(ipa) >= 1) then
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'       &
                                      &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'
      do k=1,csite%nlev_sfcwater(ipa)
         write (unit=*,fmt='(i5,1x,5(es12.4,1x))') k,y%sfcwater_energy(k)                  &
               ,y%sfcwater_tempk(k),y%sfcwater_mass(k),y%sfcwater_fracliq(k)               &
               ,y%sfcwater_depth(k)
      end do
   end if

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '

   !----- Printing the corresponding patch information (with some redundancy) -------------!
   call print_csiteipa(csite, ipa)

   call fatal_error('IFLAG1 problem. The model didn''t converge!','print_rk4patch'&
                 &,'rk4_integ_utils.f90')
   return
end subroutine print_rk4patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will perform the allocation for the Runge-Kutta integrator structure, !
! and initialize it as well.                                                               !
!------------------------------------------------------------------------------------------!
subroutine initialize_rk4patches(init)

   use ed_state_vars , only : edgrid_g              & ! intent(inout)
                            , edtype                & ! structure
                            , polygontype           & ! structure
                            , sitetype              & ! structure
                            , patchtype             ! ! structure
   use rk4_coms      , only : integration_buff      & ! structure
                            , deallocate_rk4_coh & ! structure
                            , allocate_rk4_patch    & ! structure
                            , allocate_rk4_coh   ! ! structure
   use grid_coms     , only : ngrids                ! ! intent(in)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   integer           , intent(in) :: init
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

   if (init == 0) then
      !------------------------------------------------------------------------------------!
      !    If this is not initialization, deallocate cohort memory from integration        !
      ! patches.                                                                           !
      !------------------------------------------------------------------------------------!
      call deallocate_rk4_coh(integration_buff%initp)
      call deallocate_rk4_coh(integration_buff%yscal)
      call deallocate_rk4_coh(integration_buff%y)
      call deallocate_rk4_coh(integration_buff%dydx)
      call deallocate_rk4_coh(integration_buff%yerr)
      call deallocate_rk4_coh(integration_buff%ytemp)
      call deallocate_rk4_coh(integration_buff%ak2)
      call deallocate_rk4_coh(integration_buff%ak3)
      call deallocate_rk4_coh(integration_buff%ak4)
      call deallocate_rk4_coh(integration_buff%ak5)
      call deallocate_rk4_coh(integration_buff%ak6)
      call deallocate_rk4_coh(integration_buff%ak7)
   else
      !------------------------------------------------------------------------------------!
      !     If this is initialization, make sure soil and sfcwater arrays are allocated.   !
      !------------------------------------------------------------------------------------!
      call allocate_rk4_patch(integration_buff%initp)
      call allocate_rk4_patch(integration_buff%yscal)
      call allocate_rk4_patch(integration_buff%y)
      call allocate_rk4_patch(integration_buff%dydx)
      call allocate_rk4_patch(integration_buff%yerr)
      call allocate_rk4_patch(integration_buff%ytemp)
      call allocate_rk4_patch(integration_buff%ak2)
      call allocate_rk4_patch(integration_buff%ak3)
      call allocate_rk4_patch(integration_buff%ak4)
      call allocate_rk4_patch(integration_buff%ak5)
      call allocate_rk4_patch(integration_buff%ak6)
      call allocate_rk4_patch(integration_buff%ak7)
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
   call allocate_rk4_coh(maxcohort,integration_buff%initp)
   call allocate_rk4_coh(maxcohort,integration_buff%yscal)
   call allocate_rk4_coh(maxcohort,integration_buff%y)
   call allocate_rk4_coh(maxcohort,integration_buff%dydx)
   call allocate_rk4_coh(maxcohort,integration_buff%yerr)
   call allocate_rk4_coh(maxcohort,integration_buff%ytemp)
   call allocate_rk4_coh(maxcohort,integration_buff%ak2)
   call allocate_rk4_coh(maxcohort,integration_buff%ak3)
   call allocate_rk4_coh(maxcohort,integration_buff%ak4)
   call allocate_rk4_coh(maxcohort,integration_buff%ak5)
   call allocate_rk4_coh(maxcohort,integration_buff%ak6)
   call allocate_rk4_coh(maxcohort,integration_buff%ak7)
  
   return
end subroutine initialize_rk4patches
!==========================================================================================!
!==========================================================================================!

