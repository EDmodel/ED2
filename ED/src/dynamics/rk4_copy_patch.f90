!==========================================================================================!
!==========================================================================================!
! MODULE: RK4_COPY_PATCH
!
!> \brief   Routine that copies contents of rk4 structures 
!> \details This subroutine copies the content from one RK4 structure to another.  This
!!          has been moved to a separate module to avoid circularities.
!------------------------------------------------------------------------------------------!
module rk4_copy_patch
   contains

   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine copies the values to different buffers inside the RK4 integration  !
   ! scheme.                                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine copy_rk4_patch(sourcep, targetp, cpatch)

      use rk4_coms      , only : rk4patchtype      & ! structure
                               , checkbudget       & ! intent(in)
                               , print_detailed    & ! intent(in)
                               , rk4site
      use ed_state_vars , only : sitetype          & ! structure
                               , patchtype         ! ! structure
      use grid_coms     , only : nzg               & ! intent(in)
                               , nzs               ! ! intent(in)
      use ed_max_dims   , only : n_pft             ! ! intent(in)
      use ed_misc_coms  , only : fast_diagnostics  ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(rk4patchtype) , target     :: sourcep
      type(rk4patchtype) , target     :: targetp
      type(patchtype)    , target     :: cpatch
      !----- Local variable ---------------------------------------------------------------!
      integer                         :: k
      !------------------------------------------------------------------------------------!

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
      targetp%total_sfcw_mass  = sourcep%total_sfcw_mass
      targetp%snowfac          = sourcep%snowfac

   !  These are not incremented
      targetp%vels             = sourcep%vels
      targetp%atm_enthalpy     = sourcep%atm_enthalpy

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

      targetp%cwd_rh           = sourcep%cwd_rh
      targetp%rh               = sourcep%rh

      targetp%water_deficit    = sourcep%water_deficit

      do k=rk4site%lsl,nzg      
         targetp%soil_water            (k) = sourcep%soil_water            (k)
         targetp%soil_energy           (k) = sourcep%soil_energy           (k)
         targetp%soil_mstpot           (k) = sourcep%soil_mstpot           (k)
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

         targetp%veg_resolvable  (k) = sourcep%veg_resolvable  (k)
         targetp%veg_energy      (k) = sourcep%veg_energy      (k)
         targetp%veg_water       (k) = sourcep%veg_water       (k)
         targetp%veg_hcap        (k) = sourcep%veg_hcap        (k)

         targetp%veg_wind        (k) = sourcep%veg_wind        (k)
         targetp%lint_shv        (k) = sourcep%lint_shv        (k)
         targetp%nplant          (k) = sourcep%nplant          (k)
         targetp%lai             (k) = sourcep%lai             (k)
         targetp%wai             (k) = sourcep%wai             (k)
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
         targetp%leaf_growth_resp(k) = sourcep%leaf_growth_resp(k)
         targetp%root_growth_resp(k) = sourcep%root_growth_resp(k)
         targetp%sapa_growth_resp(k) = sourcep%sapa_growth_resp(k)
         targetp%sapb_growth_resp(k) = sourcep%sapb_growth_resp(k)
         targetp%bark_growth_resp(k) = sourcep%bark_growth_resp(k)
         targetp%leaf_storage_resp(k) = sourcep%leaf_storage_resp(k)
         targetp%root_storage_resp(k) = sourcep%root_storage_resp(k)
         targetp%sapa_storage_resp(k) = sourcep%sapa_storage_resp(k)
         targetp%sapb_storage_resp(k) = sourcep%sapb_storage_resp(k)
         targetp%bark_storage_resp(k) = sourcep%bark_storage_resp(k)
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



      return
   end subroutine copy_rk4_patch
   !=======================================================================================!
   !=======================================================================================!
end module rk4_copy_patch
!==========================================================================================!
!==========================================================================================!
