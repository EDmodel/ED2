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
      integer                         :: ico
      integer                         :: k
      !------------------------------------------------------------------------------------!

      targetp%can_enthalpy     = sourcep%can_enthalpy
      targetp%can_theta        = sourcep%can_theta
      targetp%can_temp         = sourcep%can_temp
      targetp%can_shv          = sourcep%can_shv
      targetp%can_co2          = sourcep%can_co2
      targetp%can_rhos         = sourcep%can_rhos
      targetp%can_dmol         = sourcep%can_dmol
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

      targetp%fgc_rh           = sourcep%fgc_rh
      targetp%fsc_rh           = sourcep%fsc_rh
      targetp%stgc_rh          = sourcep%stgc_rh
      targetp%stsc_rh          = sourcep%stsc_rh
      targetp%msc_rh           = sourcep%msc_rh
      targetp%ssc_rh           = sourcep%ssc_rh
      targetp%psc_rh           = sourcep%psc_rh

      targetp%water_deficit    = sourcep%water_deficit

      targetp%commit_storage_resp = sourcep%commit_storage_resp
      targetp%commit_growth_resp  = sourcep%commit_growth_resp

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
         targetp%sfcwater_mass     (k) = sourcep%sfcwater_mass     (k)
         targetp%sfcwater_energy   (k) = sourcep%sfcwater_energy   (k)
         targetp%sfcwater_depth    (k) = sourcep%sfcwater_depth    (k)
         targetp%sfcwater_tempk    (k) = sourcep%sfcwater_tempk    (k)
         targetp%sfcwater_fracliq  (k) = sourcep%sfcwater_fracliq  (k)
      end do

      do ico=1,cpatch%ncohorts
         targetp%leaf_resolvable   (ico) = sourcep%leaf_resolvable   (ico)
         targetp%leaf_energy       (ico) = sourcep%leaf_energy       (ico)
         targetp%leaf_water        (ico) = sourcep%leaf_water        (ico)
         targetp%leaf_water_im2    (ico) = sourcep%leaf_water_im2    (ico)
         targetp%leaf_temp         (ico) = sourcep%leaf_temp         (ico)
         targetp%leaf_fliq         (ico) = sourcep%leaf_fliq         (ico)
         targetp%leaf_hcap         (ico) = sourcep%leaf_hcap         (ico)
         targetp%leaf_reynolds     (ico) = sourcep%leaf_reynolds     (ico)
         targetp%leaf_grashof      (ico) = sourcep%leaf_grashof      (ico)
         targetp%leaf_nussfree     (ico) = sourcep%leaf_nussfree     (ico)
         targetp%leaf_nussforc     (ico) = sourcep%leaf_nussforc     (ico)
         targetp%leaf_gbh          (ico) = sourcep%leaf_gbh          (ico)
         targetp%leaf_gbw          (ico) = sourcep%leaf_gbw          (ico)
         targetp%rshort_l          (ico) = sourcep%rshort_l          (ico)
         targetp%rlong_l           (ico) = sourcep%rlong_l           (ico)

         targetp%wood_resolvable   (ico) = sourcep%wood_resolvable   (ico)
         targetp%wood_energy       (ico) = sourcep%wood_energy       (ico)
         targetp%wood_water        (ico) = sourcep%wood_water        (ico)
         targetp%wood_water_im2    (ico) = sourcep%wood_water_im2    (ico)
         targetp%wood_temp         (ico) = sourcep%wood_temp         (ico)
         targetp%wood_fliq         (ico) = sourcep%wood_fliq         (ico)
         targetp%wood_hcap         (ico) = sourcep%wood_hcap         (ico)
         targetp%wood_reynolds     (ico) = sourcep%wood_reynolds     (ico)
         targetp%wood_grashof      (ico) = sourcep%wood_grashof      (ico)
         targetp%wood_nussfree     (ico) = sourcep%wood_nussfree     (ico)
         targetp%wood_nussforc     (ico) = sourcep%wood_nussforc     (ico)
         targetp%wood_gbh          (ico) = sourcep%wood_gbh          (ico)
         targetp%wood_gbw          (ico) = sourcep%wood_gbw          (ico)
         targetp%rshort_w          (ico) = sourcep%rshort_w          (ico)
         targetp%rlong_w           (ico) = sourcep%rlong_w           (ico)

         targetp%veg_resolvable    (ico) = sourcep%veg_resolvable    (ico)
         targetp%veg_energy        (ico) = sourcep%veg_energy        (ico)
         targetp%veg_water         (ico) = sourcep%veg_water         (ico)
         targetp%veg_water_im2     (ico) = sourcep%veg_water_im2     (ico)
         targetp%veg_hcap          (ico) = sourcep%veg_hcap          (ico)

         targetp%is_small          (ico) = sourcep%is_small          (ico)
         targetp%veg_wind          (ico) = sourcep%veg_wind          (ico)
         targetp%lint_shv          (ico) = sourcep%lint_shv          (ico)
         targetp%nplant            (ico) = sourcep%nplant            (ico)
         targetp%lai               (ico) = sourcep%lai               (ico)
         targetp%wai               (ico) = sourcep%wai               (ico)
         targetp%tai               (ico) = sourcep%tai               (ico)
         targetp%crown_area        (ico) = sourcep%crown_area        (ico)
         targetp%elongf            (ico) = sourcep%elongf            (ico)
         targetp%gsw_open          (ico) = sourcep%gsw_open          (ico)
         targetp%gsw_closed        (ico) = sourcep%gsw_closed        (ico)
         targetp%psi_open          (ico) = sourcep%psi_open          (ico)
         targetp%psi_closed        (ico) = sourcep%psi_closed        (ico)
         targetp%fs_open           (ico) = sourcep%fs_open           (ico)
         targetp%gpp               (ico) = sourcep%gpp               (ico)
         targetp%leaf_resp         (ico) = sourcep%leaf_resp         (ico)
         targetp%root_resp         (ico) = sourcep%root_resp         (ico)
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


         do ico=1,cpatch%ncohorts
            targetp%avg_sensible_lc    (ico) = sourcep%avg_sensible_lc    (ico)
            targetp%avg_sensible_wc    (ico) = sourcep%avg_sensible_wc    (ico)
            targetp%avg_vapor_lc       (ico) = sourcep%avg_vapor_lc       (ico)
            targetp%avg_vapor_wc       (ico) = sourcep%avg_vapor_wc       (ico)
            targetp%avg_transp         (ico) = sourcep%avg_transp         (ico)
            targetp%avg_intercepted_al (ico) = sourcep%avg_intercepted_al (ico)
            targetp%avg_intercepted_aw (ico) = sourcep%avg_intercepted_aw (ico)
            targetp%avg_wshed_lg       (ico) = sourcep%avg_wshed_lg       (ico)
            targetp%avg_wshed_wg       (ico) = sourcep%avg_wshed_wg       (ico)
            targetp%avg_wflux_wl       (ico) = sourcep%avg_wflux_wl       (ico)
            targetp%avg_wflux_gw       (ico) = sourcep%avg_wflux_gw       (ico)
         end do
      end if

      if (print_detailed) then
         targetp%flx_carbon_ac          = sourcep%flx_carbon_ac
         targetp%flx_carbon_st          = sourcep%flx_carbon_st
         targetp%flx_vapor_lc           = sourcep%flx_vapor_lc
         targetp%flx_vapor_wc           = sourcep%flx_vapor_wc
         targetp%flx_vapor_gc           = sourcep%flx_vapor_gc
         targetp%flx_wshed_vg           = sourcep%flx_wshed_vg
         targetp%flx_wflux_wl           = sourcep%flx_wflux_wl
         targetp%flx_wflux_gw           = sourcep%flx_wflux_gw
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

         do ico=1,cpatch%ncohorts
            targetp%cfx_hflxlc      (ico) = sourcep%cfx_hflxlc      (ico)
            targetp%cfx_hflxwc      (ico) = sourcep%cfx_hflxwc      (ico)
            targetp%cfx_qwflxlc     (ico) = sourcep%cfx_qwflxlc     (ico)
            targetp%cfx_qwflxwc     (ico) = sourcep%cfx_qwflxwc     (ico)
            targetp%cfx_qwshed      (ico) = sourcep%cfx_qwshed      (ico)
            targetp%cfx_qtransp     (ico) = sourcep%cfx_qtransp     (ico)
            targetp%cfx_qintercepted(ico) = sourcep%cfx_qintercepted(ico)
            targetp%cfx_qwflux_wl   (ico) = sourcep%cfx_qwflux_wl   (ico)
            targetp%cfx_qwflux_gw   (ico) = sourcep%cfx_qwflux_gw   (ico)
            do k=rk4site%lsl,nzg
               targetp%cfx_qwflux_gw_layer(k,ico) = sourcep%cfx_qwflux_gw_layer(k,ico)
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      return
   end subroutine copy_rk4_patch
   !=======================================================================================!
   !=======================================================================================!
end module rk4_copy_patch
!==========================================================================================!
!==========================================================================================!
