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
      integer                         :: i
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

      do i=1,cpatch%ncohorts
         targetp%leaf_resolvable   (i) = sourcep%leaf_resolvable   (i)
         targetp%leaf_energy       (i) = sourcep%leaf_energy       (i)
         targetp%leaf_water        (i) = sourcep%leaf_water        (i)
         targetp%leaf_water_im2    (i) = sourcep%leaf_water_im2    (i)
         targetp%leaf_temp         (i) = sourcep%leaf_temp         (i)
         targetp%leaf_fliq         (i) = sourcep%leaf_fliq         (i)
         targetp%leaf_hcap         (i) = sourcep%leaf_hcap         (i)
         targetp%leaf_reynolds     (i) = sourcep%leaf_reynolds     (i)
         targetp%leaf_grashof      (i) = sourcep%leaf_grashof      (i)
         targetp%leaf_nussfree     (i) = sourcep%leaf_nussfree     (i)
         targetp%leaf_nussforc     (i) = sourcep%leaf_nussforc     (i)
         targetp%leaf_gbh          (i) = sourcep%leaf_gbh          (i)
         targetp%leaf_gbw          (i) = sourcep%leaf_gbw          (i)
         targetp%rshort_l          (i) = sourcep%rshort_l          (i)
         targetp%rlong_l           (i) = sourcep%rlong_l           (i)

         targetp%wood_resolvable   (i) = sourcep%wood_resolvable   (i)
         targetp%wood_energy       (i) = sourcep%wood_energy       (i)
         targetp%wood_water        (i) = sourcep%wood_water        (i)
         targetp%wood_water_im2    (i) = sourcep%wood_water_im2    (i)
         targetp%wood_temp         (i) = sourcep%wood_temp         (i)
         targetp%wood_fliq         (i) = sourcep%wood_fliq         (i)
         targetp%wood_hcap         (i) = sourcep%wood_hcap         (i)
         targetp%wood_reynolds     (i) = sourcep%wood_reynolds     (i)
         targetp%wood_grashof      (i) = sourcep%wood_grashof      (i)
         targetp%wood_nussfree     (i) = sourcep%wood_nussfree     (i)
         targetp%wood_nussforc     (i) = sourcep%wood_nussforc     (i)
         targetp%wood_gbh          (i) = sourcep%wood_gbh          (i)
         targetp%wood_gbw          (i) = sourcep%wood_gbw          (i)
         targetp%rshort_w          (i) = sourcep%rshort_w          (i)
         targetp%rlong_w           (i) = sourcep%rlong_w           (i)

         targetp%veg_resolvable    (i) = sourcep%veg_resolvable    (i)
         targetp%veg_energy        (i) = sourcep%veg_energy        (i)
         targetp%veg_water         (i) = sourcep%veg_water         (i)
         targetp%veg_water_im2     (i) = sourcep%veg_water_im2     (i)
         targetp%veg_hcap          (i) = sourcep%veg_hcap          (i)

         targetp%veg_wind          (i) = sourcep%veg_wind          (i)
         targetp%lint_shv          (i) = sourcep%lint_shv          (i)
         targetp%nplant            (i) = sourcep%nplant            (i)
         targetp%lai               (i) = sourcep%lai               (i)
         targetp%wai               (i) = sourcep%wai               (i)
         targetp%tai               (i) = sourcep%tai               (i)
         targetp%crown_area        (i) = sourcep%crown_area        (i)
         targetp%elongf            (i) = sourcep%elongf            (i)
         targetp%gsw_open          (i) = sourcep%gsw_open          (i)
         targetp%gsw_closed        (i) = sourcep%gsw_closed        (i)
         targetp%psi_open          (i) = sourcep%psi_open          (i)
         targetp%psi_closed        (i) = sourcep%psi_closed        (i)
         targetp%fs_open           (i) = sourcep%fs_open           (i)
         targetp%gpp               (i) = sourcep%gpp               (i)
         targetp%leaf_resp         (i) = sourcep%leaf_resp         (i)
         targetp%root_resp         (i) = sourcep%root_resp         (i)
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


         do i=1,cpatch%ncohorts
            targetp%avg_sensible_lc    (i) = sourcep%avg_sensible_lc    (i)
            targetp%avg_sensible_wc    (i) = sourcep%avg_sensible_wc    (i)
            targetp%avg_vapor_lc       (i) = sourcep%avg_vapor_lc       (i)
            targetp%avg_vapor_wc       (i) = sourcep%avg_vapor_wc       (i)
            targetp%avg_transp         (i) = sourcep%avg_transp         (i)
            targetp%avg_intercepted_al (i) = sourcep%avg_intercepted_al (i)
            targetp%avg_intercepted_aw (i) = sourcep%avg_intercepted_aw (i)
            targetp%avg_wshed_lg       (i) = sourcep%avg_wshed_lg       (i)
            targetp%avg_wshed_wg       (i) = sourcep%avg_wshed_wg       (i)
            targetp%avg_wflux_wl       (i) = sourcep%avg_wflux_wl       (i)
            targetp%avg_wflux_gw       (i) = sourcep%avg_wflux_gw       (i)
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

         do i=1,cpatch%ncohorts
            targetp%cfx_hflxlc      (i) = sourcep%cfx_hflxlc      (i)
            targetp%cfx_hflxwc      (i) = sourcep%cfx_hflxwc      (i)
            targetp%cfx_qwflxlc     (i) = sourcep%cfx_qwflxlc     (i)
            targetp%cfx_qwflxwc     (i) = sourcep%cfx_qwflxwc     (i)
            targetp%cfx_qwshed      (i) = sourcep%cfx_qwshed      (i)
            targetp%cfx_qtransp     (i) = sourcep%cfx_qtransp     (i)
            targetp%cfx_qintercepted(i) = sourcep%cfx_qintercepted(i)
            targetp%cfx_qwflux_wl   (i) = sourcep%cfx_qwflux_wl   (i)
            targetp%cfx_qwflux_gw   (i) = sourcep%cfx_qwflux_gw   (i)
            do k=rk4site%lsl,nzg
               targetp%cfx_qwflux_gw_layer(k,i) = sourcep%cfx_qwflux_gw_layer(k,i)
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
