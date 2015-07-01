!==========================================================================================!
!==========================================================================================!
!      This module contains functions and routines to evaluate the budgets of water,       !
! enthalpy, and carbon dioxide.                                                            !
!------------------------------------------------------------------------------------------!
module budget_utils

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine simply updates the budget variables.                               !
   !---------------------------------------------------------------------------------------!
   subroutine update_budget(csite,lsl,ipaa,ipaz)
     
      use ed_state_vars, only : sitetype     ! ! structure
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      integer         , intent(in) :: lsl
      integer         , intent(in) :: ipaa
      integer         , intent(in) :: ipaz
      !----- Local variables. -------------------------------------------------------------!
      integer                      :: ipa
      !------------------------------------------------------------------------------------!


      do ipa=ipaa,ipaz
         !---------------------------------------------------------------------------------!
         !      Computing the storage terms for CO2, energy, and water budgets.            !
         !---------------------------------------------------------------------------------!
         csite%co2budget_initialstorage(ipa) = compute_co2_storage(csite,ipa)
         csite%wbudget_initialstorage(ipa)   = compute_water_storage(csite,lsl,ipa)
         csite%ebudget_initialstorage(ipa)   = compute_energy_storage(csite,lsl,ipa)
      end do

      return
   end subroutine update_budget
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine compute_budget(csite,lsl,pcpg,qpcpg,ipa,wcurr_loss2atm,ecurr_netrad          &
                            ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage           &
                            ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff       &
                            ,site_area,cbudget_nep,old_can_enthalpy,old_can_shv            &
                            ,old_can_co2,old_can_rhos,old_can_temp,old_can_prss)
      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use ed_max_dims  , only : str_len            ! ! intent(in)
      use ed_misc_coms , only : dtlsm              & ! intent(in)
                              , fast_diagnostics   & ! intent(in)
                              , current_time       ! ! intent(in)
      use ed_max_dims  , only : n_dbh              ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            & ! intent(in)
                              , rdry               & ! intent(in)
                              , mmdryi             & ! intent(in)
                              , epim1              ! ! intent(in)
      use rk4_coms     , only : rk4eps             & ! intent(in)
                              , print_budget       & ! intent(in)
                              , budget_pref        & ! intent(in)
                              , checkbudget        ! ! intent(in)
      use therm_lib    , only : tq2enthalpy        ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)                          , target        :: csite
      real                                    , intent(inout) :: pcpg
      real                                    , intent(inout) :: qpcpg
      real                                    , intent(inout) :: co2curr_loss2atm
      real                                    , intent(inout) :: ecurr_netrad
      real                                    , intent(inout) :: ecurr_loss2atm
      real                                    , intent(inout) :: ecurr_loss2drainage
      real                                    , intent(inout) :: ecurr_loss2runoff
      real                                    , intent(inout) :: wcurr_loss2atm
      real                                    , intent(inout) :: wcurr_loss2drainage
      real                                    , intent(inout) :: wcurr_loss2runoff
      integer                                 , intent(in)    :: lsl
      integer                                 , intent(in)    :: ipa
      real                                    , intent(in)    :: site_area
      real                                    , intent(inout) :: cbudget_nep
      real                                    , intent(in)    :: old_can_enthalpy
      real                                    , intent(in)    :: old_can_shv
      real                                    , intent(in)    :: old_can_co2
      real                                    , intent(in)    :: old_can_rhos
      real                                    , intent(in)    :: old_can_temp
      real                                    , intent(in)    :: old_can_prss
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)                         , pointer       :: cpatch
      character(len=str_len)                                  :: budget_fout
      real                                                    :: co2budget_finalstorage
      real                                                    :: co2budget_deltastorage
      real                                                    :: co2curr_gpp
      real                                                    :: co2curr_leafresp
      real                                                    :: co2curr_rootresp
      real                                                    :: co2curr_growthresp
      real                                                    :: co2curr_storageresp
      real                                                    :: co2curr_hetresp
      real                                                    :: co2curr_nep
      real                                                    :: co2curr_denseffect
      real                                                    :: co2curr_residual
      real                                                    :: ebudget_finalstorage
      real                                                    :: ebudget_deltastorage
      real                                                    :: ecurr_precipgain
      real                                                    :: ecurr_denseffect
      real                                                    :: ecurr_prsseffect
      real                                                    :: ecurr_residual
      real                                                    :: wbudget_finalstorage
      real                                                    :: wbudget_deltastorage
      real                                                    :: wcurr_precipgain
      real                                                    :: wcurr_denseffect
      real                                                    :: wcurr_residual
      real                                                    :: curr_can_enthalpy
      real                                                    :: gpp
      real                                                    :: leaf_resp
      real                                                    :: root_resp
      real                                                    :: growth_resp
      real                                                    :: storage_resp
      real                                                    :: co2_factor
      real                                                    :: ene_factor
      real                                                    :: h2o_factor
      real                                                    :: patch_lai
      real                                                    :: patch_wai
      integer                                                 :: jpa
      integer                                                 :: ico
      logical                                                 :: isthere
      logical                                                 :: co2_ok
      logical                                                 :: energy_ok
      logical                                                 :: water_ok
      !----- Local constants. -------------------------------------------------------------!
      character(len=13)     , parameter     :: fmtf='(a,1x,es14.7)'
      character(len=10)     , parameter     :: bhfmt='(31(a,1x))'
      character(len=48)     , parameter     :: bbfmt='(3(i13,1x),28(es13.6,1x))'
      !----- Locally saved variables. -----------------------------------------------------!
      logical               , save          :: first_time = .true.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      If this is the first time, we initialise all files with their headers.        !
      !------------------------------------------------------------------------------------!
      if (first_time) then
         do jpa = 1, csite%npatches
            write(budget_fout,fmt='(2a,i4.4,a)') trim(budget_pref),'patch_',jpa,'.txt'
            inquire(file=trim(budget_fout),exist=isthere)
            if (isthere) then
               !---- Open the file to delete when closing. --------------------------------!
               open (unit=86,file=trim(budget_fout),status='old',action='write')
               close(unit=86,status='delete')
            end if
            !------------------------------------------------------------------------------!

            if (print_budget) then
               !---------------------------------------------------------------------------!
               open (unit=86,file=trim(budget_fout),status='replace',action='write')
               write(unit=86,fmt=bhfmt)   '         YEAR' , '        MONTH'                &
                                        , '          DAY' , '         TIME'                &
                                        , '          LAI' , '          WAI'                &
                                        , '       HEIGHT' , '  CO2.STORAGE'                &
                                        , ' CO2.RESIDUAL' , ' CO2.DSTORAGE'                &
                                        , '      CO2.NEP' , ' CO2.DENS.EFF'                &
                                        , ' CO2.LOSS2ATM' , '  ENE.STORAGE'                &
                                        , ' ENE.RESIDUAL' , ' ENE.DSTORAGE'                &
                                        , '   ENE.PRECIP' , '   ENE.NETRAD'                &
                                        , ' ENE.DENS.EFF' , ' ENE.PRSS.EFF'                &
                                        , ' ENE.LOSS2ATM' , ' ENE.DRAINAGE'                &
                                        , '   ENE.RUNOFF' , '  H2O.STORAGE'                &
                                        , ' H2O.RESIDUAL' , ' H2O.DSTORAGE'                &
                                        , '   H2O.PRECIP' , ' H2O.DENS.EFF'                &
                                        , ' H2O.LOSS2ATM' , ' H2O.DRAINAGE'                &
                                        , '   H2O.RUNOFF'
                                        
               close(unit=86,status='keep')
            end if
         end do
         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute gain in water and energy due to precipitation.                         !
      !------------------------------------------------------------------------------------!
      wcurr_precipgain = pcpg  * dtlsm
      ecurr_precipgain = qpcpg * dtlsm
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compute the density and pressure effects.  We seek the conservation of the     !
      ! extensive properties [X/m2], but the canopy air space solves the intensive         !
      ! quantities instead [X/mol or X/kg].  Because density is not constant within the    !
      ! time step, and during the integration we solve the intensive form for the canopy   !
      ! air space, we must subtract the density effect from the residual.  The derivation  !
      ! is shown below.  The storage term is what we aim at, but in reality we solve the   !
      ! equations after the >>>.                                                           !
      !                                                                                    !
      !  dM              d(rho*z*m)                         dm                     d(rho)  !
      ! ---- = I - L -> ------------ = I - L >>> rho * z * ---- = I - L - m * z * -------- !
      !  dt                  dt                             dt                       dt    !
      !                                                                                    !
      ! where M is the extensive propery, I is the input flux, L is the loss flux, z is    !
      ! the canopy air space depth, and rho is the canopy air space density.               !
      !    For the specific case of enthalpy, we also compute the pressure effect between  !
      ! time steps.  We cannot guarantee conservation of enthalpy when we update pressure, !
      ! because of the first law of thermodynamics (the way to address this would be to    !
      ! use equivalent potential temperature, which is enthalpy plus pressure effect).     !
      ! Enthalpy is preserved within one time step, once pressure is updated and remains   !
      ! constant.                                                                          !
      !                                                                                    !
      !   dH         dp                     dh             dp             d(rho)           !
      !  ---- - V * ---- = Q >>> rho * z * ---- = Q + z * ---- - m * z * --------          !
      !   dt         dt                     dt             dt               dt             !
      !                                                                                    !
      ! where p is the canopy air space pressure, Q is the net heat exchange, V is the     !
      ! volume of the canopy air space.                                                    !
      !------------------------------------------------------------------------------------!
      !------ CO2.  Density effect only. --------------------------------------------------!
      co2curr_denseffect  = ddens_dt_effect(old_can_rhos,csite%can_rhos(ipa)               &
                                           ,old_can_co2,csite%can_co2(ipa)                 &
                                           ,csite%can_depth(ipa),mmdryi)
      !------ Water. Density effect only. -------------------------------------------------!
      wcurr_denseffect    = ddens_dt_effect(old_can_rhos,csite%can_rhos(ipa)               &
                                           ,old_can_shv,csite%can_shv(ipa)                 &
                                           ,csite%can_depth(ipa),1.)
      !------ Enthalpy.  Density and pressure effects. ------------------------------------!
      curr_can_enthalpy   = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
      ecurr_denseffect    = ddens_dt_effect(old_can_rhos,csite%can_rhos(ipa)               &
                                           ,old_can_enthalpy, curr_can_enthalpy            &
                                           ,csite%can_depth(ipa),1.0)
      ecurr_prsseffect    = csite%can_depth(ipa) * (csite%can_prss(ipa) - old_can_prss)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compute the carbon flux components.                                            !
      !------------------------------------------------------------------------------------!
      call sum_plant_cfluxes(csite,ipa,gpp,leaf_resp,root_resp,growth_resp                 &
                            ,storage_resp)
      co2curr_gpp         = gpp           * dtlsm
      co2curr_leafresp    = leaf_resp     * dtlsm
      co2curr_rootresp    = root_resp     * dtlsm
      co2curr_growthresp  = growth_resp   * dtlsm
      co2curr_storageresp = storage_resp  * dtlsm
      co2curr_hetresp     = csite%rh(ipa) * dtlsm
      co2curr_nep         = co2curr_gpp - co2curr_leafresp - co2curr_rootresp              &
                          - co2curr_growthresp - co2curr_storageresp - co2curr_hetresp
      cbudget_nep         = cbudget_nep + site_area * csite%area(ipa) * co2curr_nep        &
                                        * umol_2_kgC


      !----- Compute current storage terms. -----------------------------------------------!
      co2budget_finalstorage = compute_co2_storage(csite,ipa)
      wbudget_finalstorage   = compute_water_storage(csite,lsl,ipa)
      ebudget_finalstorage   = compute_energy_storage(csite,lsl,ipa)

      !----- Compute the change in storage. -----------------------------------------------!
      co2budget_deltastorage = co2budget_finalstorage - csite%co2budget_initialstorage(ipa)
      wbudget_deltastorage   = wbudget_finalstorage   - csite%wbudget_initialstorage(ipa)
      ebudget_deltastorage   = ebudget_finalstorage   - csite%ebudget_initialstorage(ipa)
      !------------------------------------------------------------------------------------!
      !     Compute residuals.                                                             !
      !------------------------------------------------------------------------------------!
      !----- 1. Canopy CO2. ---------------------------------------------------------------!
      co2curr_residual = co2budget_deltastorage                                            &
                       - ( - co2curr_nep       - co2curr_loss2atm  )                       &
                       - co2curr_denseffect
      !----- 2. Energy. -------------------------------------------------------------------!
      ecurr_residual   = ebudget_deltastorage                                              &
                       - ( ecurr_precipgain    - ecurr_loss2atm  - ecurr_loss2drainage     &
                         - ecurr_loss2runoff   + ecurr_netrad    + ecurr_prsseffect    )   &
                       - ecurr_denseffect
      !----- 3. Water. --------------------------------------------------------------------!
      wcurr_residual   = wbudget_deltastorage                                              &
                       - ( wcurr_precipgain    - wcurr_loss2atm                            &
                         - wcurr_loss2drainage - wcurr_loss2runoff )                       &
                       - wcurr_denseffect
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Integrate residuals.                                                           !
      !------------------------------------------------------------------------------------!
      !----- 1. Canopy CO2. ---------------------------------------------------------------!
      csite%co2budget_residual(ipa) = csite%co2budget_residual(ipa)  + co2curr_residual
      !----- 2. Energy. -------------------------------------------------------------------!
      csite%ebudget_residual(ipa) = csite%ebudget_residual(ipa)      + ecurr_residual
      !----- 3. Water. --------------------------------------------------------------------!
      csite%wbudget_residual(ipa) = csite%wbudget_residual(ipa)      + wcurr_residual
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Integrate the terms that are part of the budget.                                !
      !------------------------------------------------------------------------------------!
      !----- 1. Carbon dioxide. -----------------------------------------------------------!
      csite%co2budget_gpp(ipa)         = csite%co2budget_gpp(ipa)       + gpp       *dtlsm
      csite%co2budget_plresp(ipa)      = csite%co2budget_plresp(ipa)                       &
                                       + ( leaf_resp + root_resp + growth_resp             &
                                         + storage_resp) * dtlsm
      csite%co2budget_rh(ipa)          = csite%co2budget_rh(ipa)                           &
                                       + csite%rh(ipa) * dtlsm
      csite%co2budget_denseffect(ipa)  = csite%co2budget_denseffect(ipa)                   &
                                       + co2curr_denseffect
      csite%co2budget_loss2atm(ipa)    = csite%co2budget_loss2atm(ipa)                     &
                                       + co2curr_loss2atm
      !----- 2. Energy. -------------------------------------------------------------------!
      csite%ebudget_precipgain(ipa)    = csite%ebudget_precipgain(ipa)   + ecurr_precipgain
      csite%ebudget_netrad(ipa)        = csite%ebudget_netrad(ipa)       + ecurr_netrad
      csite%ebudget_prsseffect(ipa)    = csite%ebudget_prsseffect(ipa)   + ecurr_prsseffect
      csite%ebudget_denseffect(ipa)    = csite%ebudget_denseffect(ipa)   + ecurr_denseffect
      csite%ebudget_loss2atm(ipa)      = csite%ebudget_loss2atm(ipa)     + ecurr_loss2atm
      csite%ebudget_loss2drainage(ipa) = csite%ebudget_loss2drainage(ipa)                  &
                                       + ecurr_loss2drainage
      csite%ebudget_loss2runoff(ipa)   = csite%ebudget_loss2runoff(ipa)                    &
                                       + ecurr_loss2runoff
      !----- 3. Water. --------------------------------------------------------------------!
      csite%wbudget_precipgain(ipa)    = csite%wbudget_precipgain(ipa) + wcurr_precipgain
      csite%wbudget_denseffect(ipa)    = csite%wbudget_denseffect(ipa) + wcurr_denseffect
      csite%wbudget_loss2atm(ipa)      = csite%wbudget_loss2atm(ipa)   + wcurr_loss2atm
      csite%wbudget_loss2drainage(ipa) = csite%wbudget_loss2drainage(ipa)                  &
                                       + wcurr_loss2drainage
      csite%wbudget_loss2runoff(ipa)   = csite%wbudget_loss2runoff(ipa)                    &
                                       + wcurr_loss2runoff
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update density and initial storage for next step.                              !
      !------------------------------------------------------------------------------------!
      csite%wbudget_initialstorage(ipa)   = wbudget_finalstorage
      csite%ebudget_initialstorage(ipa)   = ebudget_finalstorage
      csite%co2budget_initialstorage(ipa) = co2budget_finalstorage


      !------------------------------------------------------------------------------------!
      !    If the "check budget" option is activated (you can turn on and turn off by set- !
      ! ting checkbudget in ed_params.f90), then the model will crash whenever there is    !
      ! some significant leak of CO2, water, or energy.                                    !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         co2_ok  = abs(co2curr_residual) <= rk4eps * ( abs(co2budget_finalstorage)         &
                                                     + abs(co2budget_deltastorage) )
         energy_ok = abs(ecurr_residual) <= rk4eps * ( abs(ebudget_finalstorage)           &
                                                     + abs(ebudget_deltastorage)   )
         water_ok  = abs(wcurr_residual) <= rk4eps * ( abs(wbudget_finalstorage)           &
                                                     + abs(wbudget_deltastorage)   )


         !---------------------------------------------------------------------------------!
         !     Find the patch LAI and WAI for output, only if it is needed.                !
         !---------------------------------------------------------------------------------!
         if (.not. (co2_ok .and. energy_ok .and. water_ok) .or. print_budget) then
            cpatch => csite%patch(ipa)
            patch_lai = 0.0
            patch_wai = 0.0
            do ico=1,cpatch%ncohorts
               patch_lai = patch_lai + cpatch%lai(ico)
               patch_wai = patch_wai + cpatch%wai(ico)
            end do
         end if
         !---------------------------------------------------------------------------------!


         if (.not. co2_ok) then
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') '|           !!! ): CO2 budget failed :( !!!        |'
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',         &
               current_time%year,current_time%month,current_time%date ,current_time%time
            write (unit=*,fmt=fmtf ) ' LAI            : ',patch_lai
            write (unit=*,fmt=fmtf ) ' VEG_HEIGHT     : ',csite%veg_height(ipa)
            write (unit=*,fmt=fmtf ) ' CAN_RHOS       : ',csite%can_rhos(ipa)
            write (unit=*,fmt=fmtf ) ' OLD_CAN_RHOS   : ',old_can_rhos
            write (unit=*,fmt=fmtf ) ' CAN_DEPTH      : ',csite%can_depth(ipa)
            write (unit=*,fmt=fmtf ) ' RESIDUAL       : ',co2curr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE: '                                  &
                                    ,csite%co2budget_initialstorage(ipa)
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE  : ',co2budget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE  : ',co2budget_deltastorage
            write (unit=*,fmt=fmtf ) ' GPP            : ',co2curr_gpp
            write (unit=*,fmt=fmtf ) ' LEAF_RESP      : ',co2curr_leafresp
            write (unit=*,fmt=fmtf ) ' ROOT_RESP      : ',co2curr_rootresp
            write (unit=*,fmt=fmtf ) ' GROWTH_RESP    : ',co2curr_growthresp
            write (unit=*,fmt=fmtf ) ' STORAGE_RESP   : ',co2curr_storageresp
            write (unit=*,fmt=fmtf ) ' HET_RESP       : ',co2curr_hetresp
            write (unit=*,fmt=fmtf ) ' NEP            : ',co2curr_nep
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT : ',co2curr_denseffect
            write (unit=*,fmt=fmtf ) ' LOSS2ATM       : ',co2curr_loss2atm
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') ' '
         end if

         if (.not. energy_ok) then
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') '|         !!! ): Enthalpy budget failed :( !!!     |'
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',         &
               current_time%year,current_time%month,current_time%date ,current_time%time
            write (unit=*,fmt=fmtf ) ' LAI            : ',patch_lai
            write (unit=*,fmt=fmtf ) ' VEG_HEIGHT     : ',csite%veg_height(ipa)
            write (unit=*,fmt=fmtf ) ' CAN_DEPTH      : ',csite%can_depth(ipa)
            write (unit=*,fmt=fmtf ) ' RESIDUAL       : ',ecurr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE: ',csite%ebudget_initialstorage(ipa)
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE  : ',ebudget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE  : ',ebudget_deltastorage
            write (unit=*,fmt=fmtf ) ' PRECIPGAIN     : ',ecurr_precipgain
            write (unit=*,fmt=fmtf ) ' NETRAD         : ',ecurr_netrad
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT : ',ecurr_denseffect
            write (unit=*,fmt=fmtf ) ' PRESSURE_EFFECT: ',ecurr_prsseffect
            write (unit=*,fmt=fmtf ) ' LOSS2ATM       : ',ecurr_loss2atm
            write (unit=*,fmt=fmtf ) ' LOSS2DRAINAGE  : ',ecurr_loss2drainage
            write (unit=*,fmt=fmtf ) ' LOSS2RUNOFF    : ',ecurr_loss2runoff
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') ' '
         end if


         if (.not. water_ok) then
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') '|          !!! ): Water budget failed :( !!!       |'
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',         &
               current_time%year,current_time%month,current_time%date ,current_time%time
            write (unit=*,fmt=fmtf ) ' LAI            : ',patch_lai
            write (unit=*,fmt=fmtf ) ' VEG_HEIGHT     : ',csite%veg_height(ipa)
            write (unit=*,fmt=fmtf ) ' CAN_DEPTH      : ',csite%can_depth(ipa)
            write (unit=*,fmt=fmtf ) ' RESIDUAL       : ',wcurr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE: ',csite%wbudget_initialstorage(ipa)
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE  : ',wbudget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE  : ',wbudget_deltastorage
            write (unit=*,fmt=fmtf ) ' PRECIPGAIN     : ',wcurr_precipgain
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT : ',wcurr_denseffect
            write (unit=*,fmt=fmtf ) ' LOSS2ATM       : ',wcurr_loss2atm
            write (unit=*,fmt=fmtf ) ' LOSS2DRAINAGE  : ',wcurr_loss2drainage
            write (unit=*,fmt=fmtf ) ' LOSS2RUNOFF    : ',wcurr_loss2runoff
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') ' '
         end if






         if (print_budget) then 
            co2_factor =      1. / dtlsm
            ene_factor =      1. / dtlsm
            h2o_factor = day_sec / dtlsm

            !----- Fix the units so the terms are expressed as fluxes. --------------------!
            co2curr_residual       = co2curr_residual       * co2_factor
            co2budget_deltastorage = co2budget_deltastorage * co2_factor
            co2curr_nep            = co2curr_nep            * co2_factor
            co2curr_denseffect     = co2curr_denseffect     * co2_factor
            co2curr_loss2atm       = co2curr_loss2atm       * co2_factor
            ecurr_residual         = ecurr_residual         * ene_factor
            ebudget_deltastorage   = ebudget_deltastorage   * ene_factor
            ecurr_precipgain       = ecurr_precipgain       * ene_factor
            ecurr_netrad           = ecurr_netrad           * ene_factor
            ecurr_denseffect       = ecurr_denseffect       * ene_factor
            ecurr_prsseffect       = ecurr_prsseffect       * ene_factor
            ecurr_loss2atm         = ecurr_loss2atm         * ene_factor
            ecurr_loss2drainage    = ecurr_loss2drainage    * ene_factor
            ecurr_loss2runoff      = ecurr_loss2runoff      * ene_factor
            wcurr_residual         = wcurr_residual         * h2o_factor
            wbudget_deltastorage   = wbudget_deltastorage   * h2o_factor
            wcurr_precipgain       = wcurr_precipgain       * h2o_factor
            wcurr_denseffect       = wcurr_denseffect       * h2o_factor
            wcurr_loss2atm         = wcurr_loss2atm         * h2o_factor
            wcurr_loss2drainage    = wcurr_loss2drainage    * h2o_factor
            wcurr_loss2runoff      = wcurr_loss2runoff      * h2o_factor
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Write the file.                                                         !
            !------------------------------------------------------------------------------!
            write(budget_fout,fmt='(2a,i4.4,a)') trim(budget_pref),'patch_',ipa,'.txt'
            open (unit=86,file=trim(budget_fout),status='old',action='write'               &
                 ,position='append')
            write(unit=86,fmt=bbfmt)                                                       &
                 current_time%year      , current_time%month     , current_time%date       &
               , current_time%time      , patch_lai              , patch_wai               &
               , csite%veg_height(ipa)  , co2budget_finalstorage , co2curr_residual        &
               , co2budget_deltastorage , co2curr_nep            , co2curr_denseffect      &
               , co2curr_loss2atm       , ebudget_finalstorage   , ecurr_residual          &
               , ebudget_deltastorage   , ecurr_precipgain       , ecurr_netrad            &
               , ecurr_denseffect       , ecurr_prsseffect       , ecurr_loss2atm          &
               , ecurr_loss2drainage    , ecurr_loss2runoff      , wbudget_finalstorage    &
               , wcurr_residual         , wbudget_deltastorage   , wcurr_precipgain        &
               , wcurr_denseffect       , wcurr_loss2atm         , wcurr_loss2drainage     &
               , wcurr_loss2runoff
            close(unit=86,status='keep')
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Stop the run in case there is any leak of CO2, enthalpy, or water.         !
         !---------------------------------------------------------------------------------!
         if (.not. (co2_ok .and. energy_ok .and. water_ok)) then
            call fatal_error('Budget check has failed, see message above!'      &
                            ,'compute_budget','budget_utils.f90')
         end if
      end if

      return
   end subroutine compute_budget
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This function computes the total water stored in the system, in kg/m2.              !
   !   (soil + temporary pools + canopy air space + leaf surface).                         !
   !---------------------------------------------------------------------------------------!
   real function compute_water_storage(csite, lsl,ipa)
      use ed_state_vars , only  : sitetype              & ! structure
                                , patchtype             ! ! structure
      use grid_coms      , only : nzg                   ! ! intent(in)
      use soil_coms      , only : dslz                  ! ! intent(in)
      use consts_coms    , only : wdns                  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      integer        , intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      integer                     :: k
      integer                     :: ico
      !------------------------------------------------------------------------------------!


      compute_water_storage = 0.0
      cpatch => csite%patch(ipa)

      !----- 1. Add the water stored in the soil. -----------------------------------------!
      do k = lsl, nzg
         compute_water_storage = compute_water_storage                                     &
                               + csite%soil_water(k,ipa) * dslz(k) * wdns
      end do
      !----- 2. Add the water stored in the temporary surface water/snow. -----------------!
      do k = 1, csite%nlev_sfcwater(ipa)
         compute_water_storage = compute_water_storage + csite%sfcwater_mass(k,ipa)
      end do
      !----- 3. Add the water vapour floating in the canopy air space. --------------------!
      compute_water_storage = compute_water_storage                                        &
                            + csite%can_shv(ipa) * csite%can_depth(ipa)                    &
                            * csite%can_rhos(ipa)
      !----- 4. Add the water on the leaf and wood surfaces. ------------------------------!
      do ico = 1,cpatch%ncohorts
         compute_water_storage = compute_water_storage + cpatch%leaf_water(ico)
         compute_water_storage = compute_water_storage + cpatch%wood_water(ico)
      end do

      return
   end function compute_water_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computs the total net radiation, by adding the radiation that        !
   ! interacts with the different surfaces.                                                !
   !---------------------------------------------------------------------------------------!
   real function compute_netrad(csite,ipa)
      use ed_state_vars , only : sitetype  & ! structure
                               , patchtype ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      integer                     :: k
      integer                     :: ico
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)

      compute_netrad = 0.0
      !----- 1. Add the ground components -------------------------------------------------!
      compute_netrad = csite%rshort_g(ipa) + csite%rlong_g(ipa) + csite%rlong_s(ipa)
      !----- 2. Add the shortwave radiation that reaches each snow/water layer ------------!
      do k = 1, csite%nlev_sfcwater(ipa)
         compute_netrad = compute_netrad + csite%rshort_s(k,ipa)
      end do
      !----- 3. Add the radiation components that is absorbed by leaves and branches. -----!
      do ico = 1,cpatch%ncohorts
         compute_netrad = compute_netrad + cpatch%rshort_l(ico) + cpatch%rlong_l(ico)
         compute_netrad = compute_netrad + cpatch%rshort_w(ico) + cpatch%rlong_w(ico)
      end do
      return
   end function compute_netrad
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computs the total net radiation, by adding the radiation that        !
   ! interacts with the different surfaces.  The result is given in J/m2.                  !
   !---------------------------------------------------------------------------------------!
   real function compute_energy_storage(csite, lsl, ipa)
      use ed_state_vars        , only : sitetype              & ! structure
                                      , patchtype             ! ! structure
      use grid_coms            , only : nzg                   ! ! intent(in)
      use soil_coms            , only : dslz                  ! ! intent(in)
      use therm_lib            , only : tq2enthalpy           ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      integer        , intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      integer                     :: k
      integer                     :: ico
      real                        :: soil_storage
      real                        :: sfcwater_storage
      real                        :: cas_storage
      real                        :: veg_storage
      real                        :: can_enthalpy
      !------------------------------------------------------------------------------------!


      cpatch => csite%patch(ipa)

      !----- 1. Computing internal energy stored at the soil. -----------------------------!
      soil_storage = 0.0
      do k = lsl, nzg
         soil_storage = soil_storage + csite%soil_energy(k,ipa) * dslz(k)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !   2. Computing internal energy stored at the temporary snow/water sfc. layer.      !
      !      Converting it to J/m2. 
      !------------------------------------------------------------------------------------!
      sfcwater_storage = 0.0
      do k = 1, csite%nlev_sfcwater(ipa)
         sfcwater_storage = sfcwater_storage                                               &
                          + csite%sfcwater_energy(k,ipa) * csite%sfcwater_mass(k,ipa)
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 3. Find and value for canopy air total enthalpy .                                  !
      !------------------------------------------------------------------------------------!
      can_enthalpy = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
      cas_storage  = csite%can_rhos(ipa) * csite%can_depth(ipa) * can_enthalpy
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 4. Compute the internal energy stored in the plants.                               !
      !------------------------------------------------------------------------------------!
      veg_storage = 0.0
      do ico = 1,cpatch%ncohorts
         veg_storage = veg_storage + cpatch%leaf_energy(ico) + cpatch%wood_energy(ico)
      end do
      !------------------------------------------------------------------------------------!


      !----- 5. Integrating the total energy in ED. ---------------------------------------!
      compute_energy_storage = soil_storage + sfcwater_storage + cas_storage + veg_storage
      !------------------------------------------------------------------------------------!

      return
   end function compute_energy_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the carbon flux terms.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine sum_plant_cfluxes(csite,ipa, gpp, leaf_resp,root_resp,growth_resp            &
                               ,storage_resp)
      use ed_state_vars        , only : sitetype    & ! structure
                                      , patchtype   ! ! structure
      use consts_coms          , only : day_sec     & ! intent(in)
                                      , umol_2_kgC  ! ! intent(in)
      use ed_max_dims          , only : n_dbh       ! ! intent(in)
      use ed_misc_coms         , only : ddbhi       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      integer               , intent(in)  :: ipa
      real                  , intent(out) :: gpp
      real                  , intent(out) :: leaf_resp
      real                  , intent(out) :: root_resp
      real                  , intent(out) :: growth_resp
      real                  , intent(out) :: storage_resp
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer            :: cpatch
      integer                             :: k
      integer                             :: ico
      integer                             :: idbh
      !------------------------------------------------------------------------------------!

      !----- Initializing some variables. -------------------------------------------------!
      gpp          = 0.0
      leaf_resp    = 0.0
      root_resp    = 0.0
      growth_resp  = 0.0
      storage_resp = 0.0
      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !     Looping over cohorts.                                                          !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         !----- Add GPP and leaf respiration only for those cohorts with enough leaves. ---!
         if (cpatch%leaf_resolvable(ico)) then
            gpp       = gpp       + cpatch%gpp(ico)
            leaf_resp = leaf_resp + cpatch%leaf_respiration(ico)

         end if
         !----- Root respiration happens even when the LAI is tiny ------------------------!
         root_resp = root_resp + cpatch%root_respiration(ico)

         !---------------------------------------------------------------------------------!
         !      Structural terms are "intensive", we must convert them from kgC/plant/day  !
         ! to umol/m2/s.                                                                   !
         !---------------------------------------------------------------------------------!
         growth_resp  = growth_resp                                                        &
                      + cpatch%growth_respiration(ico)  * cpatch%nplant(ico)               &
                      / (day_sec * umol_2_kgC)
         storage_resp = storage_resp                                                       &
                      + cpatch%storage_respiration(ico) * cpatch%nplant(ico)               &
                      / (day_sec * umol_2_kgC)
      end do

      return
   end subroutine sum_plant_cfluxes
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the co2 stored in the canopy air space from ppm to umol/m2. !
   !---------------------------------------------------------------------------------------!
   real function compute_co2_storage(csite,ipa)
      use ed_state_vars  , only : sitetype              ! ! structure
      use consts_coms    , only : mmdryi                ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      integer               , intent(in)  :: ipa
      !------------------------------------------------------------------------------------!

      compute_co2_storage = csite%can_co2(ipa) * mmdryi * csite%can_rhos(ipa)              &
                          * csite%can_depth(ipa)

      return
   end function compute_co2_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the change of the integrated value of a given property in  !
   ! the canopy air space due to change in density.                                        !
   !---------------------------------------------------------------------------------------!
   real function ddens_dt_effect(old_rhos,new_rhos,old_prop,new_prop,can_depth,multi)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real, intent(in) :: old_rhos   ! Density before time integration
      real, intent(in) :: new_rhos   ! Density after time integration
      real, intent(in) :: old_prop   ! Property before time integration
      real, intent(in) :: new_prop   ! Property after time integration
      real, intent(in) :: can_depth  ! Canopy depth
      real, intent(in) :: multi      ! Some scaling constant that may be needed.
      !------------------------------------------------------------------------------------!
      ddens_dt_effect = multi * can_depth * 0.5 * (old_prop + new_prop)                    &
                      * (new_rhos - old_rhos)
      return
   end function ddens_dt_effect
   !=======================================================================================!
   !=======================================================================================!
end module budget_utils
