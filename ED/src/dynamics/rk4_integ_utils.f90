!==========================================================================================!
!==========================================================================================!
! Subroutine odeint_ar                                                                     !
!                                                                                          !
!     This subroutine will drive the integration of several ODEs that drive the fast-scale !
! state variables.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine odeint_ar(h1,csite,ipa,isi,ipy,ifm,integration_buff,rhos,vels   &
                    ,atm_tmp,atm_shv,atm_co2,geoht,exner,pcpg,qpcpg,dpcpg,prss,lsl)

   use ed_state_vars  , only : integration_vars_ar & ! structure
                             , sitetype            & ! structure
                             , patchtype           ! ! structure
   use rk4_coms       , only : maxstp              & ! intent(in)
                             , tbeg                & ! intent(in)
                             , tend                & ! intent(in)
                             , dtrk4               & ! intent(in)
                             , dtrk4i              ! ! intent(in)
   use rk4_stepper_ar , only : rkqs_ar             ! ! subroutine
   use ed_misc_coms   , only : fast_diagnostics    ! ! intent(in)
   use hydrology_coms , only : useRUNOFF           ! ! intent(in)
   use grid_coms      , only : nzg                 ! ! intent(in)
   use soil_coms      , only : dslz                & ! intent(in)
                             , min_sfcwater_mass   & ! intent(in)
                             , runoff_time         ! ! intent(in)
   use consts_coms    , only : cliq                & ! intent(in)
                             , t3ple               & ! intent(in)
                             , tsupercool          & ! intent(in)
                             , wdnsi               ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(integration_vars_ar) , target      :: integration_buff ! RK4 variables
   type(sitetype)            , target      :: csite            ! Current site
   integer                   , intent(in)  :: ipa              ! Current patch ID
   integer                   , intent(in)  :: isi              ! Current site ID
   integer                   , intent(in)  :: ipy              ! Current polygon ID
   integer                   , intent(in)  :: ifm              ! Current grid ID
   integer                   , intent(in)  :: lsl              ! Lowest soil point
   real                      , intent(in)  :: rhos             ! Air density
   real                      , intent(in)  :: vels             ! Air wind speed
   real                      , intent(in)  :: atm_tmp          ! Air temperature
   real                      , intent(in)  :: atm_shv          ! Air specific humidity
   real                      , intent(in)  :: atm_co2          ! Air CO2 mixing ratio
   real                      , intent(in)  :: geoht            ! Geopotential height
   real                      , intent(in)  :: exner            ! Exner function
   real                      , intent(in)  :: pcpg             ! Precipitation rate
   real                      , intent(in)  :: qpcpg            ! Precipitation heat rate
   real                      , intent(in)  :: dpcpg            ! Precipitation "depth flux"
   real                      , intent(in)  :: prss             ! Air pressure
   real                                    :: h1               ! First guess of delta-t
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch           ! Current patch
   integer                                 :: i                ! Step counter
   integer                                 :: ksn              ! # of snow/water layers
   real                                    :: x                ! Elapsed time
   real                                    :: h                ! Current delta-t attempt
   real                                    :: hnext            ! Next delta-t
   real                                    :: hdid             ! delta-t that worked (???)
   real                                    :: qwfree           ! Free water internal energy
   real                                    :: wfreeb           ! Free water 
   !----- Saved variables -----------------------------------------------------------------!
   logical, save    :: first_time=.true.
   logical, save    :: simplerunoff
   !---------------------------------------------------------------------------------------!
   
   !----- Checking whether we will use runoff or not, and saving this check to save time. -!
   if (first_time) then
      simplerunoff = useRUNOFF == 0
      first_time   = .false.
   end if


   !---------------------------------------------------------------------------------------!
   !    If top snow layer is too thin for computational stability, have it evolve in       !
   ! thermal equilibrium with top soil layer.                                              !
   !---------------------------------------------------------------------------------------!
   call stabilize_snow_layers_ar(integration_buff%initp, csite,ipa, lsl)



   !---------------------------------------------------------------------------------------!
   !     Create temporary patches.                                                         !
   !---------------------------------------------------------------------------------------!
   cpatch => csite%patch(ipa)
   call copy_rk4_patch_ar(integration_buff%initp, integration_buff%y,cpatch,lsl)


   !---------------------------------------------------------------------------------------!
   ! Set initial time and stepsize.                                                        !
   !---------------------------------------------------------------------------------------!
   x = tbeg
   h = h1
   if (dtrk4 < 0.0) h = -h1

   !---------------------------------------------------------------------------------------!
   ! Begin timestep loop                                                                   !
   !---------------------------------------------------------------------------------------!
   timesteploop: do i=1,maxstp

      !----- Get initial derivatives ------------------------------------------------------!
      call leaf_derivs_ar(integration_buff%y,integration_buff%dydx,csite,ipa,isi,ipy,rhos  &
                         ,prss,pcpg,qpcpg,dpcpg,atm_tmp,exner,geoht,vels,atm_shv,atm_co2   &
                         ,lsl)

      !----- Get scalings used to determine stability -------------------------------------!
      call get_yscal_ar(integration_buff%y, integration_buff%dydx,h,integration_buff%yscal &
                       ,cpatch,csite%total_snow_depth(ipa),lsl)

      !----- Be sure not to overstep ------------------------------------------------------!
      if((x+h-tend)*(x+h-tbeg) > 0.0) h=tend-x

      !----- Take the step ----------------------------------------------------------------!
      call rkqs_ar(integration_buff,x,h,hdid,hnext,csite,ipa,isi,ipy,ifm,rhos,vels,atm_tmp &
                  ,atm_shv,atm_co2,geoht,exner,pcpg,qpcpg,dpcpg,prss,lsl)

      !----- Re-calculate tempks, fracliqs, surface water flags. --------------------------!
      call stabilize_snow_layers_ar(integration_buff%y, csite,ipa, lsl)

      !----- If the integration reached the next step, make some final adjustments --------!
      if((x-tend)*dtrk4 >= 0.0)then

         csite%wbudget_loss2runoff(ipa) = 0.0
         csite%ebudget_loss2runoff(ipa) = 0.0
         ksn = integration_buff%y%nlev_sfcwater

         !---------------------------------------------------------------------------------!
         !   Make temporary surface liquid water disappear.                                !
         !---------------------------------------------------------------------------------!
         if (simplerunoff .and. ksn >= 1) then
         
            if (integration_buff%y%sfcwater_mass(ksn) > 0.0   .and.                        &
                integration_buff%y%sfcwater_fracliq(ksn) > 0.1) then

               wfreeb = integration_buff%y%sfcwater_mass(ksn)                              &
                      * (integration_buff%y%sfcwater_fracliq(ksn) - .1)                    &
                      / 0.9 * min(1.0,runoff_time*hdid) 
               qwfree = wfreeb                                                             &
                      * cliq * (integration_buff%y%sfcwater_tempk(ksn) - tsupercool )

               integration_buff%y%sfcwater_mass(ksn) =                                     &
                                   integration_buff%y%sfcwater_mass(ksn)                   &
                                 - wfreeb

               integration_buff%y%sfcwater_depth(ksn) =                                    &
                                   integration_buff%y%sfcwater_depth(ksn)                  &
                                 - wfreeb*wdnsi

               !----- Recompute the energy removing runoff --------------------------------!
               if (integration_buff%y%sfcwater_mass(ksn) >= min_sfcwater_mass) then
                  integration_buff%y%sfcwater_energy(ksn) =                                &
                                     integration_buff%y%sfcwater_energy(ksn) - qwfree
               else
                  integration_buff%y%sfcwater_energy(ksn) = 0.0
               end if

               call stabilize_snow_layers_ar(integration_buff%y,csite,ipa,lsl)

               !----- Compute runoff for output -------------------------------------------!
               if(fast_diagnostics) then
                  csite%runoff(ipa) = csite%runoff(ipa) + wfreeb * dtrk4i
                  csite%avg_runoff(ipa) = csite%avg_runoff(ipa) + wfreeb
                  csite%avg_runoff_heat(ipa) = csite%avg_runoff_heat(ipa) + qwfree
                  csite%wbudget_loss2runoff(ipa) = wfreeb
                  csite%ebudget_loss2runoff(ipa) = qwfree
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
         call copy_rk4_patch_ar(integration_buff%y,integration_buff%initp, cpatch, lsl)
         !------ Updating the substep for next time and leave -----------------------------!
         csite%htry(ipa) = hnext

         return
      end if
      
      !----- Use hnext as the next substep ------------------------------------------------!
      h = hnext
   end do timesteploop

   !----- If it reached this point, that is really bad news... ----------------------------!
   print*,'Too many steps in routine odeint'
   call print_patch_ar(integration_buff%y, csite,ipa, lsl,atm_tmp,atm_shv,atm_co2,prss     &
                      ,exner,rhos,vels,geoht,pcpg,qpcpg,dpcpg)

   return
end subroutine odeint_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies that variables that are integrated by the Runge-Kutta solver   !
! to a buffer structure.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine copy_patch_init_ar(sourcesite,ipa, targetp, lsl)
   use ed_state_vars , only :  sitetype           & ! structure
                            ,  rk4patchtype       & ! structure
                            ,  patchtype          ! ! structure
   use grid_coms     , only :  nzg                & ! integer(in)
                            ,  nzs                ! ! integer(in) 
   use soil_coms     , only :  water_stab_thresh  & ! integer(in)
                            , min_sfcwater_mass   ! ! integer(in)
   use ed_misc_coms  , only :  fast_diagnostics   ! ! integer(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: targetp
   type(sitetype)     , target     :: sourcesite
   integer            , intent(in) :: lsl
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer    :: cpatch
   integer                         :: ico
   integer                         :: k
   !---------------------------------------------------------------------------------------!



   targetp%can_temp  = sourcesite%can_temp(ipa)
   targetp%can_shv   = sourcesite%can_shv(ipa)
   targetp%can_co2   = sourcesite%can_co2(ipa)

   do k = lsl, nzg
      targetp%soil_water(k)   = sourcesite%soil_water(k,ipa)
      targetp%soil_energy(k)  = sourcesite%soil_energy(k,ipa)
      targetp%soil_tempk(k)   = sourcesite%soil_tempk(k,ipa)
      targetp%soil_fracliq(k) = sourcesite%soil_fracliq(k,ipa)
   end do

   do k = 1, nzs
      targetp%sfcwater_mass(k)    = sourcesite%sfcwater_mass(k,ipa)
      targetp%sfcwater_depth(k)   = sourcesite%sfcwater_depth(k,ipa)
      !----- Converting sfcwater_energy to J/m² inside the Runge-Kutta integrator. --------!
      targetp%sfcwater_energy(k)  = sourcesite%sfcwater_energy(k,ipa)                      &
                                  * sourcesite%sfcwater_mass(k,ipa)
      targetp%sfcwater_tempk(k)   = sourcesite%sfcwater_tempk(k,ipa)
      targetp%sfcwater_fracliq(k) = sourcesite%sfcwater_fracliq(k,ipa)
   end do


   targetp%ustar = sourcesite%ustar(ipa)
   targetp%cstar = sourcesite%cstar(ipa)
   targetp%tstar = sourcesite%tstar(ipa)
   targetp%rstar = sourcesite%rstar(ipa)


   targetp%upwp = sourcesite%upwp(ipa)
   targetp%wpwp = sourcesite%wpwp(ipa)
   targetp%tpwp = sourcesite%tpwp(ipa)
   targetp%rpwp = sourcesite%rpwp(ipa)

  
   targetp%nlev_sfcwater = sourcesite%nlev_sfcwater(ipa)


   !----- The virtual pools should be always zero, they are temporary entities ------------!
   targetp%virtual_water = 0.0
   targetp%virtual_heat  = 0.0
   targetp%virtual_depth = 0.0

   if (targetp%nlev_sfcwater == 0) then
      targetp%virtual_flag = 2
   else
      if (targetp%sfcwater_mass(1) < min_sfcwater_mass) then
         targetp%virtual_flag = 2
      elseif (targetp%sfcwater_mass(1) < water_stab_thresh) then
         targetp%virtual_flag = 1
      else
         targetp%virtual_flag = 0
      end if
   end if

   cpatch => sourcesite%patch(ipa)
   do ico = 1,cpatch%ncohorts
      targetp%veg_water(ico)     = cpatch%veg_water(ico)
      targetp%veg_energy(ico)    = cpatch%veg_energy(ico)
   enddo

   !----- Diagnostics variables -----------------------------------------------------------!
   if(fast_diagnostics) then

      targetp%wbudget_loss2atm   = sourcesite%wbudget_loss2atm(ipa)
      targetp%ebudget_loss2atm   = sourcesite%ebudget_loss2atm(ipa)
      targetp%co2budget_loss2atm = sourcesite%co2budget_loss2atm(ipa)
      targetp%ebudget_latent     = sourcesite%ebudget_latent(ipa)
      targetp%avg_carbon_ac      = sourcesite%avg_carbon_ac(ipa)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! WHY IS THIS COMMENTED OUT? RGK                                             !!!!!
      !   targetp%avg_gpp = sourcesite%avg_gpp(ipa)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      targetp%avg_vapor_vc       = sourcesite%avg_vapor_vc(ipa)
      targetp%avg_dew_cg         = sourcesite%avg_dew_cg(ipa)
      targetp%avg_vapor_gc       = sourcesite%avg_vapor_gc(ipa)
      targetp%avg_wshed_vg       = sourcesite%avg_wshed_vg(ipa)
      targetp%avg_vapor_ac       = sourcesite%avg_vapor_ac(ipa)
      targetp%avg_transp         = sourcesite%avg_transp(ipa)
      targetp%avg_evap           = sourcesite%avg_evap(ipa)
      targetp%avg_netrad         = sourcesite%avg_netrad(ipa)
      targetp%aux                = sourcesite%aux(ipa)
      targetp%avg_sensible_vc    = sourcesite%avg_sensible_vc(ipa)
      targetp%avg_sensible_2cas  = sourcesite%avg_sensible_2cas(ipa)
      targetp%avg_qwshed_vg      = sourcesite%avg_qwshed_vg(ipa)
      targetp%avg_sensible_gc    = sourcesite%avg_sensible_gc(ipa)
      targetp%avg_sensible_ac    = sourcesite%avg_sensible_ac(ipa)
      targetp%avg_sensible_tot   = sourcesite%avg_sensible_tot(ipa)

      do k = lsl, nzg
         targetp%avg_sensible_gg(k) = sourcesite%avg_sensible_gg(k,ipa)
         targetp%avg_smoist_gg(k)   = sourcesite%avg_smoist_gg(k,ipa)
         targetp%avg_smoist_gc(k)   = sourcesite%avg_smoist_gc(k,ipa)
         targetp%aux_s(k)           = sourcesite%aux_s(k,ipa)
      end do
   end if

   return
end subroutine copy_patch_init_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutines increment the derivative into the previous guess to create the new   !
! guess.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine inc_rk4_patch_ar(rkp, inc, fac, cpatch, lsl)
   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         & ! structure
                            , rk4patchtype      ! ! structure
   use grid_coms     , only : nzg               & ! intent(in)
                            , nzs               ! ! intent(in)
   use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)
  
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: rkp    ! Temporary patch with previous state
   type(rk4patchtype) , target     :: inc    ! Temporary patch with its derivatives
   type(patchtype)    , target     :: cpatch ! Current patch (for characteristics)
   real               , intent(in) :: fac    ! Increment factor
   integer            , intent(in) :: lsl    ! Lowest soil level
   !----- Local variables -----------------------------------------------------------------!
   integer                         :: ico    ! Cohort ID
   integer                         :: k      ! Counter
   !---------------------------------------------------------------------------------------!


   rkp%can_temp = rkp%can_temp  + fac * inc%can_temp
   rkp%can_shv  = rkp%can_shv   + fac * inc%can_shv
   rkp%can_co2  = rkp%can_co2   + fac * inc%can_co2

   do k=lsl,nzg
      rkp%soil_water(k)       = rkp%soil_water(k)  + dble(fac) * inc%soil_water(k)
      rkp%soil_energy(k)      = rkp%soil_energy(k) + fac * inc%soil_energy(k)
   end do

   do k=1,rkp%nlev_sfcwater
      rkp%sfcwater_mass(k)   = max(0.0,rkp%sfcwater_mass(k)   + fac * inc%sfcwater_mass(k))
      rkp%sfcwater_energy(k) = rkp%sfcwater_energy(k)         + fac * inc%sfcwater_energy(k)
      rkp%sfcwater_depth(k)  = rkp%sfcwater_depth(k)          + fac * inc%sfcwater_depth(k)
   end do

   rkp%virtual_heat  = rkp%virtual_heat  + fac * inc%virtual_heat
   rkp%virtual_water = rkp%virtual_water + fac * inc%virtual_water
   rkp%virtual_depth = rkp%virtual_depth + fac * inc%virtual_depth

  
   rkp%upwp = rkp%upwp + fac * inc%upwp
   rkp%wpwp = rkp%wpwp + fac * inc%wpwp
   rkp%tpwp = rkp%tpwp + fac * inc%tpwp
   rkp%rpwp = rkp%rpwp + fac * inc%rpwp

  
   do ico = 1,cpatch%ncohorts
      rkp%veg_water(ico)     = max(0.0,rkp%veg_water(ico) + fac * inc%veg_water(ico))
      rkp%veg_energy(ico)    = rkp%veg_energy(ico) + fac * inc%veg_energy(ico)
   enddo

   if(fast_diagnostics) then

      rkp%wbudget_loss2atm   = rkp%wbudget_loss2atm   + fac * inc%wbudget_loss2atm
      rkp%ebudget_loss2atm   = rkp%ebudget_loss2atm   + fac * inc%ebudget_loss2atm
      rkp%co2budget_loss2atm = rkp%co2budget_loss2atm + fac * inc%co2budget_loss2atm
      rkp%ebudget_latent     = rkp%ebudget_latent     + fac * inc%ebudget_latent

      rkp%avg_carbon_ac      = rkp%avg_carbon_ac      + fac * inc%avg_carbon_ac
      rkp%avg_gpp            = rkp%avg_gpp            + fac * inc%avg_gpp
      
      rkp%avg_vapor_vc       = rkp%avg_vapor_vc       + fac * inc%avg_vapor_vc
      rkp%avg_dew_cg         = rkp%avg_dew_cg         + fac * inc%avg_dew_cg  
      rkp%avg_vapor_gc       = rkp%avg_vapor_gc       + fac * inc%avg_vapor_gc
      rkp%avg_wshed_vg       = rkp%avg_wshed_vg       + fac * inc%avg_wshed_vg
      rkp%avg_vapor_ac       = rkp%avg_vapor_ac       + fac * inc%avg_vapor_ac
      rkp%avg_transp         = rkp%avg_transp         + fac * inc%avg_transp  
      rkp%avg_evap           = rkp%avg_evap           + fac * inc%avg_evap  
      rkp%avg_netrad         = rkp%avg_netrad         + fac * inc%avg_netrad      
      rkp%aux                = rkp%aux                + fac * inc%aux
      rkp%avg_sensible_vc    = rkp%avg_sensible_vc    + fac * inc%avg_sensible_vc  
      rkp%avg_sensible_2cas  = rkp%avg_sensible_2cas  + fac * inc%avg_sensible_2cas
      rkp%avg_qwshed_vg      = rkp%avg_qwshed_vg      + fac * inc%avg_qwshed_vg    
      rkp%avg_sensible_gc    = rkp%avg_sensible_gc    + fac * inc%avg_sensible_gc  
      rkp%avg_sensible_ac    = rkp%avg_sensible_ac    + fac * inc%avg_sensible_ac  
      rkp%avg_sensible_tot   = rkp%avg_sensible_tot   + fac * inc%avg_sensible_tot 

      do k=lsl,nzg
         rkp%avg_sensible_gg(k)  = rkp%avg_sensible_gg(k)  + fac * inc%avg_sensible_gg(k)
         rkp%avg_smoist_gg(k)    = rkp%avg_smoist_gg(k)    + fac * inc%avg_smoist_gg(k)  
         rkp%avg_smoist_gc(k)    = rkp%avg_smoist_gc(k)    + fac * inc%avg_smoist_gc(k)  
         rkp%aux_s(k)            = rkp%aux_s(k)            + fac * inc%aux_s(k)
      end do

   end if

   return
end subroutine inc_rk4_patch_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the error scale for the integrated variables, which will be     !
! later used to define the relative error.                                                 !
!------------------------------------------------------------------------------------------!
subroutine get_yscal_ar(y, dy, htry, yscal, cpatch, total_snow_depth, lsl)
   use ed_state_vars        , only : patchtype          & ! subroutine
                                   , rk4patchtype       ! ! subroutine
   use rk4_coms             , only : tiny_offset        ! ! intent(in)
   use grid_coms            , only : nzg                & ! intent(in)
                                   , nzs                ! ! intent(in)
   use soil_coms            , only : min_sfcwater_mass  ! ! intent(in)
   use consts_coms          , only : cliq               & ! intent(in)
                                   , qliqt3             ! ! intent(in)
   use canopy_radiation_coms, only : lai_min            ! ! intent(in)
   use pft_coms             , only : sla                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype), target     :: y                ! Structure with the guesses
   type(rk4patchtype), target     :: dy               ! Structure with their derivatives
   type(rk4patchtype), target     :: yscal            ! Structure with their scales
   type(patchtype)   , target     :: cpatch           ! Current patch
   integer           , intent(in) :: lsl              ! Lowest soil level
   real              , intent(in) :: total_snow_depth ! Snow depth
   real              , intent(in) :: htry             ! Time-step we are trying
   !----- Local variables -----------------------------------------------------------------!
   integer                        :: k      ! Counter
   integer                        :: ico    ! Current cohort ID
   !---------------------------------------------------------------------------------------!

  
   yscal%can_temp = abs(y%can_temp) + abs(dy%can_temp*htry) + tiny_offset
   yscal%can_shv  = abs(y%can_shv)  + abs(dy%can_shv*htry)  + tiny_offset
   yscal%can_co2  = abs(y%can_co2)  + abs(dy%can_co2*htry)  + tiny_offset
  
   yscal%upwp = max(abs(y%upwp) + abs(dy%upwp*htry),1.0)
   yscal%wpwp = max(abs(y%wpwp) + abs(dy%wpwp*htry),1.0)


  
   do k=lsl,nzg
      yscal%soil_water(k)  = abs(y%soil_water(k))  + abs(dy%soil_water(k)*htry)            &
                           + tiny_offset
      yscal%soil_energy(k) = abs(y%soil_energy(k)) + abs(dy%soil_energy(k)*htry)
   end do

   if(y%sfcwater_mass(1) > 0.1 .or. y%nlev_sfcwater > 1)then
      !----- Either frozen or computationally stable layer. -------------------------------!
      do k=1,nzs
         yscal%sfcwater_mass(k)   = abs(y%sfcwater_mass(k))                                &
                                  + abs(dy%sfcwater_mass(k)*htry)     + tiny_offset
         yscal%sfcwater_energy(k) = abs(y%sfcwater_energy(k))                              &
                                  + abs(dy%sfcwater_energy(k)*htry) + tiny_offset
         yscal%sfcwater_depth(k)  = abs(y%sfcwater_depth(k))                               &
                                  + abs(dy%sfcwater_depth(k)*htry)   + tiny_offset
      end do
   else
      !----- Low stability threshold ------------------------------------------------------!
      do k=1,nzs
         yscal%sfcwater_mass(k) = 0.1
         if(y%sfcwater_mass(k) > min_sfcwater_mass)then
            yscal%sfcwater_energy(k) = ( yscal%sfcwater_mass(k) / y%sfcwater_mass(k))      &
                                     * ( abs( y%sfcwater_energy(k))                        &
                                       + abs(dy%sfcwater_energy(k)))
         else
            yscal%sfcwater_energy(k) = 1.0e30
         end if
      end do
   end if

   !----- Scale for the virtual water pools -----------------------------------------------!
   yscal%virtual_water = 0.1
   yscal%virtual_heat  = qliqt3 * yscal%virtual_water

   !---------------------------------------------------------------------------------------!
   !    Scale for leaf water and energy. In case the plants have few or no leaves, or the  !
   ! plant is buried in snow, we assign huge values for typical scale, thus preventing     !
   ! unecessary small steps.                                                               !
   !    Also, if the cohort is tiny and has almost no water, make the scale less strict.   !
   !---------------------------------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      if (cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > total_snow_depth)  &
      then
         yscal%veg_water(ico)  = 0.22
         yscal%veg_energy(ico) = max(abs(y%veg_energy(ico)) + abs(dy%veg_energy(ico)*htry) &
                                    ,yscal%veg_water(ico)*qliqt3)
      else
         yscal%veg_water(ico)  = 1.e30
         yscal%veg_energy(ico) = 1.e30
      end if
   end do


   return
end subroutine get_yscal_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine loops through the integrating variables, seeking for the largest      !
! error.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine get_errmax_ar(errmax,yerr,yscal,cpatch,lsl,y,ytemp)

   use ed_state_vars         , only : patchtype     & ! structure
                                    , rk4patchtype  ! ! structure
   use rk4_coms              , only : rk4eps        ! ! intent(in)
   use grid_coms             , only : nzg           ! ! intent(in)
   use canopy_radiation_coms , only : lai_min       ! ! intent(in)
   use misc_coms             , only : integ_err     & ! intent(in)
                                    , record_err    ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target      :: yerr      ! Error structure
   type(rk4patchtype) , target      :: yscal     ! Scale structure
   type(rk4patchtype) , target      :: y         ! Structure with previous value
   type(rk4patchtype) , target      :: ytemp     ! Structure with attempted values
   type(patchtype)    , target      :: cpatch    ! Current patch
   integer            , intent(in)  :: lsl       ! Lowest soil level
   real               , intent(out) :: errmax    ! Maximum error
   !----- Local variables -----------------------------------------------------------------!
   integer                          :: ico       ! Current cohort ID
   real                             :: errh2o    ! Scratch error variable
   real                             :: errene    ! Scratch error variable
   real                             :: err       ! Scratch error variable
   real                             :: errh2oMAX ! Scratch error variable
   real                             :: erreneMAX ! Scratch error variable
   integer                          :: k         ! Scratch error variable
   !---------------------------------------------------------------------------------------!

   !----- Initialize error ----------------------------------------------------------------!
   errmax = 0.0

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
  
   do k=lsl,nzg
      err    = sngl(abs(yerr%soil_water(k)/yscal%soil_water(k)))
      errmax = max(errmax,err)
      if(record_err .and. err > rk4eps) integ_err(3+k,1) = integ_err(3+k,1) + 1_8 
   end do

   do k=lsl,nzg
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
      errh2oMAX = 0.0
      erreneMAX = 0.0
      if (cpatch%lai(ico) > lai_min) then
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
end subroutine get_errmax_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine print_errmax_ar(errmax,yerr,yscal,cpatch,lsl,y,ytemp)
   use ed_state_vars         , only : patchtype     & ! Structure
                                    , rk4patchtype  ! ! Structure
   use rk4_coms              , only : rk4eps        ! ! intent(in)
   use grid_coms             , only : nzg           & ! intent(in)
                                    , nzs           ! ! intent(in)
   use canopy_radiation_coms , only : lai_min       ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target       :: yerr,yscal,y,ytemp
   type(patchtype)    , target       :: cpatch
   integer            , intent(in)   :: lsl
   real               , intent(out)  :: errmax
   !----- Local variables -----------------------------------------------------------------!
   integer                           :: ico
   integer                           :: k
   real                              :: error_soil_water
   real                              :: scale_soil_water
   logical                           :: troublemaker
   !----- Constants -----------------------------------------------------------------------!
   character(len=28)  , parameter    :: onefmt = '(a16,1x,3(es12.5,1x),11x,l1)'
   character(len=34)  , parameter    :: lyrfmt = '(a16,1x,i6,1x,3(es12.5,1x),11x,l1)'
   character(len=34)  , parameter    :: cohfmt = '(a16,1x,i6,1x,4(es12.5,1x),11x,l1)'
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
   write(unit=*,fmt='(6(a,1x))')  'Name            ',' Level','   Max.Error'         &
                                &,'   Abs.Error','       Scale','Problem(T|F)'

   do k=lsl,nzg
      errmax = sngl(dmax1(dble(errmax),dabs(yerr%soil_water(k)/yscal%soil_water(k))))
      error_soil_water = sngl(yerr%soil_water(k))
      scale_soil_water = sngl(yscal%soil_water(k))
      troublemaker = large_error(error_soil_water,scale_soil_water)
      write(unit=*,fmt=lyrfmt) 'SOIL_WATER:',k,errmax,sngl(yerr%soil_water(k))             &
                                            ,sngl(yscal%soil_water(k)),troublemaker

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
   write(unit=*,fmt='(a)'      ) ' Cohort_level variables (only those with sufficient LAI):'
   write(unit=*,fmt='(7(a,1x))')  'Name            ','         PFT','         LAI'         &
                             &,'   Max.Error','   Abs.Error','       Scale','Problem(T|F)'
   do ico = 1,cpatch%ncohorts
      if (cpatch%lai(ico) > lai_min) then
         errmax       = max(errmax,abs(yerr%veg_water(ico)/yscal%veg_water(ico)))
         troublemaker = large_error(yerr%veg_water(ico),yscal%veg_water(ico))
         write(unit=*,fmt=cohfmt) 'VEG_WATER:',cpatch%pft(ico),cpatch%lai(ico),errmax      &
                                              ,yerr%veg_water(ico),yscal%veg_water(ico)    &
                                              ,troublemaker
              

         errmax       = max(errmax,abs(yerr%veg_energy(ico)/yscal%veg_energy(ico)))
         troublemaker = large_error(yerr%veg_energy(ico),yscal%veg_energy(ico))
         write(unit=*,fmt=cohfmt) 'VEG_ENERGY:',cpatch%pft(ico),cpatch%lai(ico),errmax     &
                                               ,yerr%veg_energy(ico),yscal%veg_energy(ico) &
                                               ,troublemaker
      end if
   end do

   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('-',k=1,80)

   return
end subroutine print_errmax_ar
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
   real, intent(in) :: err  ! Absolute error
   real, intent(in) :: scal ! Characteristic scale
   !---------------------------------------------------------------------------------------!
   if(scal > 0.0) then
      large_error = abs(err/scal)/rk4eps > 1.0
   else
      large_error = .false.
   end if
   return
end function large_error
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine is the driver for snow layer stability check, and it will call the   !
! procedure to redistribute the snow layers should the partition become numerically        !
! dangerous.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine stabilize_snow_layers_ar(initp, csite,ipa, lsl)
   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         & ! structure
                            , rk4patchtype      ! ! structure
   use soil_coms     , only : soil              & ! intent(in)
                            , min_sfcwater_mass ! ! intent(in)
   use grid_coms     , only : nzg               & ! intent(in)
                            , nzs               ! ! intent(in)
   use therm_lib     , only : qwtk8             & ! subroutine
                            , qtk               ! ! subroutine
   use consts_coms   , only : wdns              ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: initp
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: lsl
   integer            , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   integer                         :: k
   real                            :: soilhcap
   !---------------------------------------------------------------------------------------!


   !----- Updating soil temperature and liquid water fraction. ----------------------------!
   do k = lsl, nzg - 1
      soilhcap = soil(csite%ntext_soil(k,ipa))%slcpd
      call qwtk8(initp%soil_energy(k),initp%soil_water(k)*dble(wdns),soilhcap              &
                ,initp%soil_tempk(k),initp%soil_fracliq(k))
   end do

   !---------------------------------------------------------------------------------------!
   !    Updating surface water temperature and liquid water fraction, remembering that in- !
   ! side the RK4 integration, surface water energy is in J/m².                            !
   !---------------------------------------------------------------------------------------!
   do k = 2, nzs
      if(initp%sfcwater_mass(k) > min_sfcwater_mass)  then
           call qtk(initp%sfcwater_energy(k)/initp%sfcwater_mass(k)                        &
                   ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
       else
          initp%sfcwater_energy(k)  = 0.
       end if
   end do
  
   !---------------------------------------------------------------------------------------!
   !    Checking whether the snow layers need to be adjusted based on their changes in     !
   ! mass, forcing excessive water to percolate, and performing a quick thermal adjust if  !
   ! the temporary layer is way too thin.                                                  !
   !---------------------------------------------------------------------------------------!
   call redistribute_snow_ar(initp,csite,ipa)
   
  
   return
end subroutine stabilize_snow_layers_ar
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
subroutine redistribute_snow_ar(initp,csite,ipa)

   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         & ! structure
                            , rk4patchtype      ! ! structure
   use grid_coms     , only : nzs               & ! intent(in)
                            , nzg               ! ! intent(in)
   use soil_coms     , only : soil              & ! intent(in)
                            , water_stab_thresh & ! intent(in)
                            , dslz              & ! intent(in)
                            , dslzi             & ! intent(in)
                            , snowmin           & ! intent(in)
                            , min_sfcwater_mass ! ! intent(in)
   use consts_coms   , only : cice              & ! intent(in)
                            , cliq              & ! intent(in)
                            , t3ple             & ! intent(in)
                            , wdns              & ! intent(in)
                            , tsupercool        & ! intent(in)
                            , qliqt3            & ! intent(in)
                            , wdnsi             ! ! intent(in)
   use therm_lib     , only : qtk               & ! subroutine
                            , qwtk              & ! subroutine
                            , qwtk8             ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)     , target     :: initp
   type(sitetype)         , target     :: csite
   integer                , intent(in) :: ipa
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: kzs
   integer                             :: kold
   integer                             :: newlayers
   integer                             :: nlayers
   integer                             :: ksn
   integer                             :: ksnnew
   integer                             :: k
   real                                :: stretch
   real                                :: thik
   !----- Control variables ---------------------------------------------------------------!
   logical                , save       :: first_call = .true.
   real, dimension(20)    , save       :: thicknet
   real, dimension(20,20) , save       :: thick   
   real                                :: wtold
   real                                :: wtnew
   real, dimension(nzs)                :: newsfcw_mass
   real, dimension(nzs)                :: newsfcw_energy
   real, dimension(nzs)                :: newsfcw_depth
   real                                :: wdiff
   real                                :: totsnow
   real                                :: depthgain
   real                                :: wfree
   real                                :: qwfree
   real                                :: qw
   real                                :: w
   real                                :: wfreeb
   real                                :: depthloss
   real                                :: snden
   real                                :: sndenmin
   real                                :: sndenmax
   real                                :: qwt
   real(kind=8)                        :: wt
   real                                :: soilhcap
   real                                :: free_surface_water_demand
   integer                             :: nsoil
   !----- Constants -----------------------------------------------------------------------!
   logical                , parameter  :: debug = .false.
   !---------------------------------------------------------------------------------------!
   

   !----- Defining some thickness variables only once, at start. --------------------------!
   if (first_call) then
      first_call = .false.

      stretch = 2.0
      do kzs = 1,nzs
         thik = 1.0
         thicknet(kzs) = 0.0
         do k = 1,(kzs+1)/2
            thick(k,kzs)       = thik
            thick(kzs+1-k,kzs) = thik
            thicknet(kzs)      = thicknet(kzs) + 2. * thik
            thik               = thik * stretch
         end do
         if ((kzs+1)/2 .ne. kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
         do k = 1,kzs
            thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
         end do
      end do
   end if

   !----- Initializing total snow and # of layers alias -----------------------------------!
   totsnow = 0.0
   ksn     = initp%nlev_sfcwater
   
   !---------------------------------------------------------------------------------------!
   !    Checking whether we still have water/snow, and depending on this, where to add the !
   ! virtual pool.                                                                         !
   !---------------------------------------------------------------------------------------!
   if (initp%nlev_sfcwater >= 1) then
      !----- Nothing left, no more temporary water/snow layer -----------------------------!
      if (sum(initp%sfcwater_mass) < min_sfcwater_mass) then
         initp%nlev_sfcwater = 0
         ksnnew              = 0
      !----- Still something there, nothing changes at least now --------------------------!
      else
         wfree = 0.0
         qwfree = 0.0
         depthgain = 0.0
         ksnnew = ksn
      end if
   !----- No temporary layer, nothing in the virtual pool, keep it dry --------------------!
   elseif (initp%virtual_water < min_sfcwater_mass) then
      ksnnew = 0
      initp%virtual_water = 0.0
      initp%virtual_heat  = 0.0
      initp%virtual_depth = 0.0
   !----- No temporary layer, significant mass to add, create one layer -------------------!
   else
      wfree = initp%virtual_water
      qwfree = initp%virtual_heat
      depthgain = initp%virtual_depth
      initp%virtual_water = 0.0
      initp%virtual_heat  = 0.0
      initp%virtual_depth = 0.0
      ksnnew = 1
   end if
   !---------------------------------------------------------------------------------------!


   !----- Loop over layers, from top to bottom this time ----------------------------------!
   do k = ksnnew,1,-1

      !----- Update current mass and energy of temporary layer ----------------------------!
      qw = initp%sfcwater_energy(k) + qwfree
      w = initp%sfcwater_mass(k) + wfree

      !------------------------------------------------------------------------------------!
      !    Single layer, and this is a very thin one, which can cause numerical instabili- !
      ! ty. Force a fast heat exchange between this thin layer and the soil topmost level, !
      ! bringing both layers to a thermal equilibrium.                                     !
      !------------------------------------------------------------------------------------!
      if (ksnnew == 1 .and. initp%sfcwater_mass(k) < water_stab_thresh) then

         !---------------------------------------------------------------------------------!
         !     Total internal energy and water of the combined system, in J/m² and kg/m²,  !
         ! respectively.                                                                   !
         !---------------------------------------------------------------------------------!
         qwt = qw      + initp%soil_energy(nzg) * dslz(nzg)
         wt  = dble(w) + initp%soil_water(nzg)  * dble(dslz(nzg)) * dble(wdns)

         !----- Finding the equilibrium temperature and liquid/ice partition. -------------!
         soilhcap = soil(csite%ntext_soil(nzg,ipa))%slcpd * dslz(nzg)
         call qwtk8(qwt,wt,soilhcap,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))

         !---------------------------------------------------------------------------------!
         !    Computing internal energy of the temporary layer with the temperature and    !
         ! liquid/ice distribution we just found, for the mass the layer has.              !
         !---------------------------------------------------------------------------------!
         qw = w * (initp%sfcwater_fracliq(k) * cliq *(initp%sfcwater_tempk(k)-tsupercool)  &
                  + (1.-initp%sfcwater_fracliq(k)) * cice * initp%sfcwater_tempk(k))

         !---------------------------------------------------------------------------------!
         !    Set the properties of top soil layer. Since internal energy is an extensive  !
         ! quantity, we can simply take the difference to be the soil internal energy,     !
         ! just remembering that we need to convert it back to J/m³. The other properties  !
         ! can be copied from the surface layer because we assumed phase and temperature   !
         ! equilibrium.                                                                    !
         !---------------------------------------------------------------------------------!
         initp%soil_energy(nzg) = (qwt - qw) * dslzi(nzg)
         initp%soil_tempk(nzg) = initp%sfcwater_tempk(k)
         initp%soil_fracliq(nzg) = initp%sfcwater_fracliq(k)
      else
         !----- Layer is computationally stable, just update the temperature and phase ----!
         call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*dble(wdns)                &
                   ,soil(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)            &
                   ,initp%soil_fracliq(nzg))
         call qtk(qw/w,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
      end if


      !------------------------------------------------------------------------------------!
      !    Shed liquid in excess of a 1:9 liquid-to-ice ratio through percolation.  Limit  !
      ! this shed amount (wfreeb) in lowest snow layer to amount top soil layer can hold.  !
      !------------------------------------------------------------------------------------!
      if (initp%sfcwater_fracliq(k) == 1.0) then
         wfreeb = w
      else
         wfreeb = max(0.0, w * (initp%sfcwater_fracliq(k)-0.1)/0.9 )
      endif
      if (k == 1)then
           !----- Do "greedy" infiltration. -----------------------------------------------!
           nsoil = csite%ntext_soil(nzg,ipa)
           free_surface_water_demand = sngl(dmax1(dble(0.0)                                &
                                     , dble(soil(nsoil)%slmsts) - initp%soil_water(nzg))   &
                                     * dble(wdns) * dble(dslz(nzg)))
           wfreeb = min(wfreeb,free_surface_water_demand)
           qwfree = wfreeb * cliq * (initp%sfcwater_tempk(k)-tsupercool)
           !----- Update topmost soil moisture and energy, updating temperature and phase -!
           initp%soil_water(nzg)  = initp%soil_water(nzg)  + dble(wfreeb*wdnsi*dslzi(nzg)) 
           initp%soil_energy(nzg) = initp%soil_energy(nzg) + qwfree * dslzi(nzg)
           soilhcap = soil(nsoil)%slcpd
           call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*dble(wdns)              &
                     ,soilhcap,initp%soil_tempk(nzg),initp%soil_fracliq(nzg))
      else
         !---- Not the first layer, just shed all free water, and compute its energy ------!
         qwfree = wfreeb * cliq * (initp%sfcwater_tempk(k)-tsupercool)
      end if
      depthloss = wfreeb * wdnsi
      
      !----- Remove water and internal energy losses due to percolation -------------------!
      initp%sfcwater_mass(k)  = w - wfreeb
      initp%sfcwater_depth(k) = initp%sfcwater_depth(k) + depthgain - depthloss
      if(initp%sfcwater_mass(k) >= min_sfcwater_mass)then
         initp%sfcwater_energy(k) = qw - qwfree
         call qtk(initp%sfcwater_energy(k)/initp%sfcwater_mass(k),initp%sfcwater_tempk(k)  &
                 ,initp%sfcwater_fracliq(k))
      else
         initp%sfcwater_energy(k) = 0.0
      endif

      !----- Integrate total "snow" -------------------------------------------------------!
      totsnow = totsnow + initp%sfcwater_mass(k)

      !----- Calculate density and depth of snow ------------------------------------------!
      snden    = initp%sfcwater_mass(k) / max(1.0e-6,initp%sfcwater_depth(k))
      sndenmax = wdns
      sndenmin = max(30.0, 200.0 * (wfree + wfreeb) / max(1.0e-12,initp%sfcwater_mass(k)))
      snden    = min(sndenmax, max(sndenmin,snden))
      initp%sfcwater_depth(k) = initp%sfcwater_mass(k) / snden

      !----- Set up input to next layer ---------------------------------------------------!
      wfree = wfreeb
      depthgain = depthloss
   enddo

   !----- Re-distribute snow layers to maintain prescribed distribution of mass. ----------!
   if (totsnow < min_sfcwater_mass .or. ksnnew == 0) then
      initp%nlev_sfcwater = 0
      call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg) * dble(wdns)                 &
                ,soil(csite%ntext_soil(nzg,ipa))%slcpd,initp%soil_tempk(nzg)               &
                ,initp%soil_fracliq(nzg))
   else
      !---- Check whether there is enough snow for a new layer. ---------------------------!
      nlayers   = ksnnew
      newlayers = 1
      do k = 1,nzs
         !----- Checking whether we need 
         if (      initp%sfcwater_mass(k)   >= min_sfcwater_mass                           &
             .and. snowmin * thicknet(k)    <= totsnow                                     &
             .and. initp%sfcwater_energy(k) <  initp%sfcwater_mass(k)*qliqt3 ) then

            newlayers = newlayers + 1
         end if
      end do
      newlayers = min(newlayers, nzs, nlayers + 1)
      initp%nlev_sfcwater = newlayers
      kold  = 1
      wtnew = 1.0
      wtold = 1.0
      do k = 1,newlayers
         newsfcw_mass(k)   = totsnow * thick(k,newlayers)
         newsfcw_energy(k) = 0.0
         newsfcw_depth(k)  = 0.0
         !----- Finding new layer properties ----------------------------------------------!
         find_layer: do

            !----- Difference between old and new snow ------------------------------------!
            wdiff = wtnew * newsfcw_mass(k) - wtold * initp%sfcwater_mass(kold)  

            if (wdiff > 0.0) then
               newsfcw_energy(k) = newsfcw_energy(k) + wtold * initp%sfcwater_energy(kold)
               newsfcw_depth(k)  = newsfcw_depth(k)  + wtold * initp%sfcwater_depth(kold)
               wtnew  = wtnew - wtold * initp%sfcwater_mass(kold) / newsfcw_mass(k)
               kold   = kold + 1
               wtold  = 1.0
               if (kold > nlayers) exit find_layer
            else
               newsfcw_energy(k) = newsfcw_energy(k) + wtnew * newsfcw_mass(k)             &
                                 * initp%sfcwater_energy(kold)                             &
                                 / max(min_sfcwater_mass,initp%sfcwater_mass(kold))
               newsfcw_depth(k)  = newsfcw_depth(k)  + wtnew * newsfcw_mass(k)             &
                                 * initp%sfcwater_depth(kold)                              &
                                 / max(min_sfcwater_mass,initp%sfcwater_mass(kold))
               wtold = wtold - wtnew * newsfcw_mass(k)                                     &
                             / max(min_sfcwater_mass,initp%sfcwater_mass(kold))
               wtnew = 1.
               exit find_layer
            end if
         end do find_layer
      end do

      !----- Updating the water/snow layer properties -------------------------------------!
      do k = 1,newlayers
         initp%sfcwater_mass(k)   = newsfcw_mass(k)
         initp%sfcwater_energy(k) = newsfcw_energy(k)
         if (newsfcw_mass(k) >= min_sfcwater_mass) then
            call qtk(initp%sfcwater_energy(k)/initp%sfcwater_mass(k)                       &
                    ,initp%sfcwater_tempk(k),initp%sfcwater_fracliq(k))
         else
            initp%sfcwater_tempk(k) = t3ple
            initp%sfcwater_fracliq(k) = 0.
         end if
         initp%sfcwater_depth(k) = newsfcw_depth(k)
      end do

      !----- Making sure that the unused layers have zero in everything -------------------!
      do k = newlayers + 1, nzs
         initp%sfcwater_mass(k)    = 0.0
         initp%sfcwater_energy(k)  = 0.0
         initp%sfcwater_depth(k)   = 0.0
         initp%sfcwater_tempk(k)   = t3ple
         initp%sfcwater_fracliq(k) = 0.
      end do
   end if

   return
end subroutine redistribute_snow_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine copies the values to different buffers inside the RK4 integration     !
! scheme.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine copy_rk4_patch_ar(sourcep, targetp, cpatch, lsl)

   use ed_state_vars , only : sitetype          & ! structure
                            , patchtype         & ! structure
                            , rk4patchtype      ! ! structure
   use grid_coms     , only : nzg               & ! intent(in)
                            , nzs               ! ! intent(in)
   use max_dims      , only : n_pft             ! ! intent(in)
   use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: sourcep
   type(rk4patchtype) , target     :: targetp
   type(patchtype)    , target     :: cpatch
   integer            , intent(in) :: lsl
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
   targetp%rpwp          = sourcep%rpwp

   targetp%ground_shv    = sourcep%ground_shv
   targetp%surface_ssh   = sourcep%surface_ssh

   targetp%nlev_sfcwater = sourcep%nlev_sfcwater
   targetp%ustar         = sourcep%ustar
   targetp%cstar         = sourcep%cstar
   targetp%tstar         = sourcep%tstar
   targetp%rstar         = sourcep%rstar
   targetp%virtual_flag  = sourcep%virtual_flag
   targetp%rasveg        = sourcep%rasveg
   targetp%root_res_fac  = sourcep%root_res_fac

   do k=lsl,nzg
      
      targetp%soil_water(k)             = sourcep%soil_water(k)
      targetp%soil_energy(k)            = sourcep%soil_energy(k)
      targetp%soil_tempk(k)             = sourcep%soil_tempk(k)
      targetp%soil_fracliq(k)           = sourcep%soil_fracliq(k)
      targetp%available_liquid_water(k) = sourcep%available_liquid_water(k)
      targetp%extracted_water(k)        = sourcep%extracted_water(k)

   end do

   do k=1,nzs
      targetp%sfcwater_mass(k)    = sourcep%sfcwater_mass(k)   
      targetp%sfcwater_energy(k)  = sourcep%sfcwater_energy(k) 
      targetp%sfcwater_depth(k)   = sourcep%sfcwater_depth(k)  
      targetp%sfcwater_tempk(k)   = sourcep%sfcwater_tempk(k)  
      targetp%sfcwater_fracliq(k) = sourcep%sfcwater_fracliq(k)
   end do

   do k = 1, n_pft
      targetp%a_o_max(k) = sourcep%a_o_max(k)
      targetp%a_c_max(k) = sourcep%a_c_max(k)
   end do
   
   do k=1,cpatch%ncohorts
      targetp%veg_water(k)   = sourcep%veg_water(k)
      targetp%veg_energy(k)  = sourcep%veg_energy(k)
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!     WHY IS THIS COMMENTED OUT? IS IT NOT INTEGRATED? REMEMBER TO DOUBLE     !!!!
      !!!!! CHECK AND THEN REMOVE THESE COMMENTS IF SO.                                 !!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    targetp%avg_gpp = sourcep%avg_gpp

      do k=lsl,nzg
         targetp%avg_sensible_gg(k) = sourcep%avg_sensible_gg(k)
         targetp%avg_smoist_gg(k)   = sourcep%avg_smoist_gg(k)  
         targetp%avg_smoist_gc(k)   = sourcep%avg_smoist_gc(k)  
         targetp%aux_s(k)           = sourcep%aux_s(k)
      end do
   end if



   return
end subroutine copy_rk4_patch_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine prints the patch and cohort information when the model falls apart... !
!------------------------------------------------------------------------------------------!
subroutine print_patch_pss_ar(csite, ipa, lsl)
   use ed_state_vars         , only : sitetype      & ! structure
                                    , patchtype     ! ! structure
   use misc_coms             , only : current_time  ! ! intent(in)
   use grid_coms             , only : nzs           & ! intent(in)
                                    , nzg           ! ! intent(in)
   use max_dims              , only : n_pft         ! ! intent(in)
   use canopy_radiation_coms , only : lai_min       ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: lsl
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
   write(unit=*,fmt='(a,1x,es12.5)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
 
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Cohort information (only those with LAI > LAI_MIN shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),11(a12,1x))')                                              &
         '    PFT','KRDEPTH','      NPLANT','         LAI','         DBH','       BDEAD'   &
                           &,'      BALIVE','  VEG_ENERGY','    VEG_TEMP','   VEG_WATER'   &
                           &,'     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if(cpatch%lai(ico) > lai_min)then
         write(unit=*,fmt='(2(i7,1x),11(es12.5,1x))') cpatch%pft(ico), cpatch%krdepth(ico) &
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
   write (unit=*,fmt='(i12,1x,6(es12.5,1x))')  csite%dist_type(ipa),csite%age(ipa)         &
         ,csite%area(ipa),csite%rh(ipa),csite%avg_daily_temp(ipa),csite%sum_chd(ipa)       &
         ,csite%sum_dgd(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','         LAI'          &
                                   &,'        HTRY','     CAN_CO2','    CAN_TEMP'          &
                                   &,'     CAN_SHV'
   write (unit=*,fmt='(7(es12.5,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)          &
         ,csite%lai(ipa),csite%htry(ipa),csite%can_co2(ipa),csite%can_temp(ipa)            &
         ,csite%can_shv(ipa) 

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       RSTAR','       CSTAR'          &
                                   &,'       TSTAR','     RLONG_G','    RSHORT_G'          &
                                   &,'     RLONG_S'
   write (unit=*,fmt='(7(es12.5,1x))') csite%ustar(ipa),csite%rstar(ipa),csite%cstar(ipa)  &
         ,csite%tstar(ipa),csite%rlong_g(ipa),csite%rshort_g(ipa),csite%rlong_s(ipa)

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,a12)') '  PFT','       REPRO'
   do k=1,n_pft
      write (unit=*,fmt='(i5,1x,es12.5)') k,csite%repro(k,ipa)
   end do

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.5,1x))') k,csite%ntext_soil(k,ipa)              &
            ,csite%soil_energy(k,ipa),csite%soil_tempk(k,ipa),csite%soil_water(k,ipa)      &
            ,csite%soil_fracliq(k,ipa)
   end do
   
   if (csite%nlev_sfcwater(ipa) >= 1) then
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a5,1x,6(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'       &
                                      &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'       &
                                      &,'    RSHORT_S'
      do k=1,csite%nlev_sfcwater(ipa)
         write (unit=*,fmt='(i5,1x,6(es12.5,1x))') k,csite%sfcwater_energy(k,ipa)          &
               ,csite%sfcwater_tempk(k,ipa),csite%sfcwater_mass(k,ipa)                     &
               ,csite%sfcwater_fracliq(k,ipa),csite%sfcwater_depth(k,ipa)                  &
               ,csite%rshort_s(k,ipa)
      end do
   end if

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '
   return
end subroutine print_patch_pss_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is similar to print_patch_pss_ar, except that it also prints the      !
! outcome of the Runge-Kutta integrator.                                                   !
!------------------------------------------------------------------------------------------!
subroutine print_patch_ar(y,csite,ipa, lsl,atm_tmp,atm_shv,atm_co2,prss,exner,rhos,vels    &
                         ,geoht,pcpg,qpcpg,dpcpg)
   use ed_state_vars         , only : sitetype          & ! intent(in) 
                                    , patchtype         & ! intent(in) 
                                    , rk4patchtype      ! ! intent(in) 
   use grid_coms             , only : nzg               & ! intent(in) 
                                    , nzs               ! ! intent(in) 
   use canopy_radiation_coms , only : lai_min           ! ! intent(in) 
   use misc_coms             , only : current_time      ! ! intent(in) 
   use therm_lib             , only : qtk               & ! subroutine 
                                    , qwtk              ! ! subroutine
   use soil_coms             , only : min_sfcwater_mass ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) , target     :: y
   type(sitetype)     , target     :: csite
   integer            , intent(in) :: lsl
   integer            , intent(in) :: ipa
   real               , intent(in) :: atm_tmp
   real               , intent(in) :: atm_shv
   real               , intent(in) :: atm_co2
   real               , intent(in) :: prss
   real               , intent(in) :: exner
   real               , intent(in) :: rhos
   real               , intent(in) :: vels
   real               , intent(in) :: geoht
   real               , intent(in) :: pcpg
   real               , intent(in) :: qpcpg
   real               , intent(in) :: dpcpg
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)    , pointer    :: cpatch
   integer                         :: k
   integer                         :: ico
   real                            :: virtual_temp, virtual_fliq, veg_temp, veg_fliq
   !---------------------------------------------------------------------------------------!

   cpatch => csite%patch(ipa)

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (rk4patch) ||||'

   write(unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,f12.0,1x,a)')                                &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%time,'s'
   write(unit=*,fmt='(a,1x,es12.5)') 'Attempted step size:',csite%htry(ipa)
   write (unit=*,fmt='(a,1x,i6)')    'Ncohorts: ',cpatch%ncohorts
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(80a)')         ('-',k=1,80)
   write (unit=*,fmt='(a)')           ' ATMOSPHERIC CONDITIONS: '
   write (unit=*,fmt='(a,1x,es12.5)') ' Air temperature    : ',atm_tmp
   write (unit=*,fmt='(a,1x,es12.5)') ' H2Ov mixing ratio  : ',atm_shv
   write (unit=*,fmt='(a,1x,es12.5)') ' CO2  mixing ratio  : ',atm_co2
   write (unit=*,fmt='(a,1x,es12.5)') ' Pressure           : ',prss
   write (unit=*,fmt='(a,1x,es12.5)') ' Exner function     : ',exner
   write (unit=*,fmt='(a,1x,es12.5)') ' Air density        : ',rhos
   write (unit=*,fmt='(a,1x,es12.5)') ' Wind speed         : ',vels
   write (unit=*,fmt='(a,1x,es12.5)') ' Height             : ',geoht
   write (unit=*,fmt='(a,1x,es12.5)') ' Precip. mass  flux : ',pcpg
   write (unit=*,fmt='(a,1x,es12.5)') ' Precip. heat  flux : ',qpcpg
   write (unit=*,fmt='(a,1x,es12.5)') ' Precip. depth flux : ',dpcpg

   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) 'Cohort information (only those with LAI > LAI_MIN shown): '
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),8(a12,1x))')                                               &
         '    PFT','KRDEPTH','      NPLANT','         LAI','         DBH','       BDEAD'   &
                           &,'      BALIVE','     FS_OPEN','         FSW','         FSN'
   do ico = 1,cpatch%ncohorts
      if (cpatch%lai(ico) > lai_min) then
         write(unit=*,fmt='(2(i7,1x),8(es12.5,1x))') cpatch%pft(ico), cpatch%krdepth(ico)  &
              ,cpatch%nplant(ico),cpatch%lai(ico),cpatch%dbh(ico),cpatch%bdead(ico)        &
              ,cpatch%balive(ico),cpatch%fs_open(ico),cpatch%fsw(ico),cpatch%fsn(ico)
      end if
   end do
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(2(a7,1x),5(a12,1x))')                                               &
         '    PFT','KRDEPTH','  VEG_ENERGY','   VEG_WATER' ,'   VEG_HCAP'                  &
                           &,'    VEG_TEMP','    VEG_FLIQ'
   do ico = 1,cpatch%ncohorts
      if (cpatch%lai(ico) > lai_min) then
         call qwtk(y%veg_energy(ico),y%veg_water(ico),cpatch%hcapveg(ico)                  &
                  ,veg_temp,veg_fliq)
         write(unit=*,fmt='(2(i7,1x),5(es12.5,1x))') cpatch%pft(ico), cpatch%krdepth(ico) &
               ,y%veg_energy(ico),y%veg_water(ico),cpatch%hcapveg(ico),veg_temp,veg_fliq
      end if
   end do
   write (unit=*,fmt='(80a)') ('=',k=1,80)
   write (unit=*,fmt='(a)'  ) ' '
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(6(a12,1x))')  '  VEG_HEIGHT','   VEG_ROUGH','         LAI'          &
                                   &,'     CAN_CO2','    CAN_TEMP','     CAN_SHV'
   write (unit=*,fmt='(6(es12.5,1x))') csite%veg_height(ipa),csite%veg_rough(ipa)          &
         ,csite%lai(ipa),y%can_co2,y%can_temp,y%can_shv

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(4(a12,1x))')  '       USTAR','       RSTAR','       CSTAR'          &
                                   &,'       TSTAR'
   write (unit=*,fmt='(4(es12.5,1x))') y%ustar,y%rstar,y%cstar,y%tstar

   write (unit=*,fmt='(80a)') ('-',k=1,80)
   if (y%virtual_water /= 0.) then
      call qtk(y%virtual_heat/y%virtual_water,virtual_temp,virtual_fliq)
   else
      virtual_temp = y%soil_tempk(nzg)
      virtual_fliq = y%soil_fracliq(nzg)
   end if


   write (unit=*,fmt='(5(a12,1x))')  'VIRTUAL_FLAG','VIRTUAL_HEAT','  VIRT_WATER'          &
                                   &,'VIRTUAL_TEMP','VIRTUAL_FLIQ'
   write (unit=*,fmt='(i12,1x,4(es12.5,1x))') y%virtual_flag,y%virtual_heat                &
                                             ,y%virtual_water,virtual_temp,virtual_fliq
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(2(a12,1x))')    '  GROUND_SHV',' SURFACE_SSH'
   write (unit=*,fmt='(2(es12.5,1x))') y%ground_shv, y%surface_ssh

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZG','  NTEXT_SOIL',' SOIL_ENERGY'          &
                                   &,'  SOIL_TEMPK','  SOIL_WATER','SOIL_FRACLIQ'
   do k=lsl,nzg
      write (unit=*,fmt='(i5,1x,i12,4(es12.5,1x))') k,csite%ntext_soil(k,ipa)              &
            ,y%soil_energy(k),y%soil_tempk(k),y%soil_water(k),y%soil_fracliq(k)
   end do
   
   if (csite%nlev_sfcwater(ipa) >= 1) then
      write (unit=*,fmt='(80a)') ('-',k=1,80)
      write (unit=*,fmt='(a5,1x,5(a12,1x))')   '  KZS',' SFCW_ENERGY','  SFCW_TEMPK'       &
                                      &,'   SFCW_MASS','SFCW_FRACLIQ','  SFCW_DEPTH'
      do k=1,csite%nlev_sfcwater(ipa)
         write (unit=*,fmt='(i5,1x,5(es12.5,1x))') k,y%sfcwater_energy(k)                  &
               ,y%sfcwater_tempk(k),y%sfcwater_mass(k),y%sfcwater_fracliq(k)               &
               ,y%sfcwater_depth(k)
      end do
   end if

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '

   !----- Printing the corresponding patch information (with some redundancy) -------------!
   call print_patch_pss_ar(csite, ipa, lsl)

   call fatal_error('IFLAG1 problem. The model didn''t converge!','print_patch_ar'&
                 &,'rk4_integ_utils.f90')
   return
end subroutine print_patch_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will perform the allocation for the Runge-Kutta integrator structure, !
! and initialize it as well.                                                               !
!------------------------------------------------------------------------------------------!
subroutine initialize_rk4patches_ar(init)

   use ed_state_vars , only : edgrid_g            & ! intent(inout)
                            , edtype              & ! structure
                            , polygontype         & ! structure
                            , sitetype            & ! structure
                            , patchtype           & ! structure
                            , integration_buff_g  ! ! intent(inout)
   use grid_coms     , only : ngrids              ! ! intent(in)
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
      call deallocate_rk4_coh_ar(integration_buff_g%initp)
      call deallocate_rk4_coh_ar(integration_buff_g%yscal)
      call deallocate_rk4_coh_ar(integration_buff_g%y)
      call deallocate_rk4_coh_ar(integration_buff_g%dydx)
      call deallocate_rk4_coh_ar(integration_buff_g%yerr)
      call deallocate_rk4_coh_ar(integration_buff_g%ytemp)
      call deallocate_rk4_coh_ar(integration_buff_g%ak2)
      call deallocate_rk4_coh_ar(integration_buff_g%ak3)
      call deallocate_rk4_coh_ar(integration_buff_g%ak4)
      call deallocate_rk4_coh_ar(integration_buff_g%ak5)
      call deallocate_rk4_coh_ar(integration_buff_g%ak6)
      call deallocate_rk4_coh_ar(integration_buff_g%ak7)
   else
      !------------------------------------------------------------------------------------!
      !     If this is initialization, make sure soil and sfcwater arrays are allocated.   !
      !------------------------------------------------------------------------------------!
      call allocate_rk4_patch(integration_buff_g%initp)
      call allocate_rk4_patch(integration_buff_g%yscal)
      call allocate_rk4_patch(integration_buff_g%y)
      call allocate_rk4_patch(integration_buff_g%dydx)
      call allocate_rk4_patch(integration_buff_g%yerr)
      call allocate_rk4_patch(integration_buff_g%ytemp)
      call allocate_rk4_patch(integration_buff_g%ak2)
      call allocate_rk4_patch(integration_buff_g%ak3)
      call allocate_rk4_patch(integration_buff_g%ak4)
      call allocate_rk4_patch(integration_buff_g%ak5)
      call allocate_rk4_patch(integration_buff_g%ak6)
      call allocate_rk4_patch(integration_buff_g%ak7)
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
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%initp)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%yscal)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%y)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%dydx)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%yerr)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ytemp)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ak2)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ak3)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ak4)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ak5)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ak6)
   call allocate_rk4_coh_ar(maxcohort,integration_buff_g%ak7)
  
   return
end subroutine initialize_rk4patches_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will perform the temporary patch allocation for the RK4 integration.  !
!------------------------------------------------------------------------------------------!
subroutine allocate_rk4_patch(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   use grid_coms     , only : nzg          & ! intent(in)
                            , nzs          ! ! intent(in)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!

   call nullify_rk4_patch(y)

   allocate(y%soil_energy(nzg))
   allocate(y%soil_water(nzg))
   allocate(y%soil_fracliq(nzg))
   allocate(y%soil_tempk(nzg))
   allocate(y%available_liquid_water(nzg))
   allocate(y%extracted_water(nzg))

   allocate(y%sfcwater_energy(nzs))
   allocate(y%sfcwater_mass(nzs))
   allocate(y%sfcwater_depth(nzs))
   allocate(y%sfcwater_fracliq(nzs))
   allocate(y%sfcwater_tempk(nzs))

   !---------------------------------------------------------------------------------------!
   !     Diagnostics - for now we will always allocate the diagnostics, even if they       !
   !                   aren't used.                                                        !
   !---------------------------------------------------------------------------------------!
   allocate(y%avg_smoist_gg(nzg))
   allocate(y%avg_smoist_gc(nzg))
   allocate(y%aux_s(nzg))
   allocate(y%avg_sensible_gg(nzg))

   call zero_rk4_patch(y)

   return
end subroutine allocate_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will nullify all pointers, to make a safe allocation.                !
!------------------------------------------------------------------------------------------!
subroutine nullify_rk4_patch(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!

   nullify(y%soil_energy)
   nullify(y%soil_water)
   nullify(y%soil_fracliq)
   nullify(y%soil_tempk)
   nullify(y%available_liquid_water)
   nullify(y%extracted_water)

   nullify(y%sfcwater_energy)
   nullify(y%sfcwater_mass)
   nullify(y%sfcwater_depth)
   nullify(y%sfcwater_fracliq)
   nullify(y%sfcwater_tempk)

   !---------------------------------------------------------------------------------------!
   !     Diagnostics - for now we will always allocate the diagnostics, even if they       !
   !                   aren't used.                                                        !
   !---------------------------------------------------------------------------------------!
   nullify(y%avg_smoist_gg)
   nullify(y%avg_smoist_gc)
   nullify(y%aux_s)
   nullify(y%avg_sensible_gg)

   return
end subroutine nullify_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Forcing all variables to be zero.                                                     !
!------------------------------------------------------------------------------------------!
subroutine zero_rk4_patch(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   use grid_coms     , only : nzg          & ! intent(in)
                            , nzs          ! ! intent(in)
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!

   y%wbudget_loss2atm               = 0.
   y%ebudget_loss2atm               = 0.
   y%ebudget_latent                 = 0.
   y%co2budget_loss2atm             = 0.
  
   y%can_temp                       = 0.
   y%can_shv                        = 0.
   y%can_co2                        = 0.
  
   y%soil_energy(:)                 = 0.
   y%soil_tempk(:)                  = 0.
   y%soil_fracliq(:)                = 0.
   y%soil_water(:)                  = 0.d0
  
   y%sfcwater_depth(:)              = 0.
   y%sfcwater_mass(:)               = 0.
   y%sfcwater_energy(:)             = 0.
  
   y%virtual_water                  = 0.
   y%virtual_heat                   = 0.
   y%virtual_depth                  = 0.
  
   y%ground_shv                     = 0.
   y%surface_ssh                    = 0.
   y%sfcwater_tempk(:)              = 0.
   y%sfcwater_fracliq(:)            = 0.
   y%nlev_sfcwater                  = 0
   y%net_rough_length               = 0.
  
   y%rough                          = 0.
  
   y%ustar                          = 0.
   y%cstar                          = 0.
   y%tstar                          = 0.
   y%rstar                          = 0.
   y%virtual_flag                   = 0
   y%avg_carbon_ac                  = 0.
  
   y%upwp                           = 0.
   y%wpwp                           = 0.
   y%tpwp                           = 0.
   y%rpwp                           = 0.
  
   y%avg_gpp                        = 0.
  
   y%a_o_max                        = 0.
   y%a_c_max                        = 0.
   y%rasveg                         = 0.
   y%root_res_fac                   = 0.
   y%available_liquid_water(:)      = 0.
   y%extracted_water(:)             = 0.
  

   y%avg_vapor_vc                   = 0.
   y%avg_dew_cg                     = 0.
   y%avg_vapor_gc                   = 0.
   y%avg_wshed_vg                   = 0.
   y%avg_vapor_ac                   = 0.
   y%avg_transp                     = 0.
   y%avg_evap                       = 0.
   y%avg_netrad                     = 0.
   y%avg_smoist_gg                  = 0.
   y%avg_smoist_gc                  = 0.
   y%aux                            = 0.
   y%aux_s                          = 0.
   y%avg_sensible_vc                = 0.
   y%avg_sensible_2cas              = 0.
   y%avg_qwshed_vg                  = 0.
   y%avg_sensible_gc                = 0.
   y%avg_sensible_ac                = 0.
   y%avg_sensible_tot               = 0.
   y%avg_sensible_gg                = 0.
   y%avg_heatstor_veg               = 0.
  

  return
end subroutine zero_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will perform the temporary patch deallocation.                        !
!------------------------------------------------------------------------------------------!
subroutine deallocate_rk4_patch(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   implicit none
   !----- Argument ------------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!

   if (associated(y%soil_energy))             deallocate(y%soil_energy)
   if (associated(y%soil_water))              deallocate(y%soil_water)
   if (associated(y%soil_fracliq))            deallocate(y%soil_fracliq)
   if (associated(y%soil_tempk))              deallocate(y%soil_tempk)
   if (associated(y%available_liquid_water))  deallocate(y%available_liquid_water)
   if (associated(y%extracted_water))         deallocate(y%extracted_water)

   if (associated(y%sfcwater_energy))         deallocate(y%sfcwater_energy)
   if (associated(y%sfcwater_mass))           deallocate(y%sfcwater_mass)
   if (associated(y%sfcwater_depth))          deallocate(y%sfcwater_depth)
   if (associated(y%sfcwater_fracliq))        deallocate(y%sfcwater_fracliq)
   if (associated(y%sfcwater_tempk))          deallocate(y%sfcwater_tempk)
   
   ! Diagnostics
   if (associated(y%avg_smoist_gg))           deallocate(y%avg_smoist_gg)
   if (associated(y%avg_smoist_gc))           deallocate(y%avg_smoist_gc)
   if (associated(y%aux_s))                   deallocate(y%aux_s)
   if (associated(y%avg_sensible_gg))         deallocate(y%avg_sensible_gg)

   return
end subroutine deallocate_rk4_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will allocate the cohorts of the temporary patch.                     !
!------------------------------------------------------------------------------------------!

subroutine allocate_rk4_coh_ar(maxcohort,y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype)              :: y
   integer            , intent(in) :: maxcohort
   !---------------------------------------------------------------------------------------!
   
   call nullify_rk4_cohort(y)

   allocate(y%veg_energy(maxcohort))
   allocate(y%veg_water(maxcohort))

   call zero_rk4_cohort(y)

   return
end subroutine allocate_rk4_coh_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will nullify the cohort pointers for a safe allocation.               !
!------------------------------------------------------------------------------------------!
subroutine nullify_rk4_cohort(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!
       
   nullify(y%veg_energy)
   nullify(y%veg_water)

   return
end subroutine nullify_rk4_cohort
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will initialize the cohort variables with zeroes.                     !
!------------------------------------------------------------------------------------------!
subroutine zero_rk4_cohort(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!

   if(associated(y%veg_energy    ))  y%veg_energy    = 0.
   if(associated(y%veg_water     ))  y%veg_water     = 0.

   return
end subroutine zero_rk4_cohort
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will deallocate the cohorts of the temporary patch.                   !
!------------------------------------------------------------------------------------------!
subroutine deallocate_rk4_coh_ar(y)
   use ed_state_vars , only : rk4patchtype ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(rk4patchtype) :: y
   !---------------------------------------------------------------------------------------!

   if(associated(y%veg_energy))     deallocate(y%veg_energy)
   if(associated(y%veg_water))      deallocate(y%veg_water)

   return
end subroutine deallocate_rk4_coh_ar
!==========================================================================================!
!==========================================================================================!
