!==========================================================================================!
!==========================================================================================!
!      This module contains functions and routines to evaluate the budgets of water,       !
! enthalpy, and carbon dioxide.                                                            !
!------------------------------------------------------------------------------------------!
module budget_utils

   implicit none


   !---------------------------------------------------------------------------------------!
   !      This variable has the tolerance for sub-daily checks of energy, water, CO2 and   !
   ! carbon conservation, which are performed every photosynthesis step.                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: tol_subday_budget
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      This variable has the tolerance for long-term carbon checks (phenology, growth   !
   ! of living and structural tissues, reproduction, and decomposition).  This tolerance   !
   ! is normally stricter than the sub-daily because the changes are relatively small.     !
   !---------------------------------------------------------------------------------------!
   real(kind=4) :: tol_carbon_budget
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine initialises all the budget variables.  This is called at the     !
   ! beginning of the simulation, after all variables have been initialised, for all       !
   ! patches.  This subroutine is then called every year after the vegetation dynamics,    !
   ! but only for recent patches.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine ed_init_budget(cgrid,initial)
      use ed_state_vars, only : edtype       & ! structure
                              , polygontype  & ! structure
                              , sitetype     ! ! structure
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      logical          , intent(in) :: initial
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      integer                       :: ipy
      integer                       :: isi
      integer                       :: lsl
      integer                       :: ipa
      !----- Local constants. -------------------------------------------------------------!
      real             , parameter  :: age_min = 1.0 / 24.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop through all polygons, sites, and patches, and reset virtual pool.         !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            lsl   =  cpoly%lsl (isi)
            patchloop: do ipa=1,csite%npatches
               !---------------------------------------------------------------------------!
               !     Reset patch budget variables, but first check whether this is         !
               ! the initialisation step, or if the patch has age zero (just created).     !
               !---------------------------------------------------------------------------!
               if (initial .or. (csite%age(ipa) < age_min)) then
                  !---- Reset patch budget variables. -------------------------------------!
                  call initial_patch_budget(csite,lsl,ipa)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do patchloop
         end do siteloop
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_init_budget
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Call this subroutine whenever all budget terms should be reset.  This should be   !
   ! called only when the patch is new or during the model initialisation.  All budget     !
   ! terms will be set to zero, and the initial storage will be computed based on the      !
   ! carbon stocks.                                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine initial_patch_budget(csite,lsl,ipa)
      use ed_state_vars, only : sitetype ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)   , target     :: csite
      integer          , intent(in) :: lsl
      integer          , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      !------------------------------------------------------------------------------------!



      !----- Set all flux terms to zero. --------------------------------------------------!
      csite%co2budget_loss2atm   (ipa) = 0.0
      csite%co2budget_denseffect (ipa) = 0.0
      csite%co2budget_zcaneffect (ipa) = 0.0
      csite%co2budget_gpp        (ipa) = 0.0
      csite%co2budget_plresp     (ipa) = 0.0
      csite%co2budget_rh         (ipa) = 0.0
      csite%co2budget_residual   (ipa) = 0.0
      csite%cbudget_loss2atm     (ipa) = 0.0
      csite%cbudget_loss2yield   (ipa) = 0.0
      csite%cbudget_seedrain     (ipa) = 0.0
      csite%cbudget_denseffect   (ipa) = 0.0
      csite%cbudget_zcaneffect   (ipa) = 0.0
      csite%cbudget_residual     (ipa) = 0.0
      csite%ebudget_precipgain   (ipa) = 0.0
      csite%ebudget_netrad       (ipa) = 0.0
      csite%ebudget_loss2atm     (ipa) = 0.0
      csite%ebudget_loss2runoff  (ipa) = 0.0
      csite%ebudget_loss2drainage(ipa) = 0.0
      csite%ebudget_denseffect   (ipa) = 0.0
      csite%ebudget_prsseffect   (ipa) = 0.0
      csite%ebudget_hcapeffect   (ipa) = 0.0
      csite%ebudget_wcapeffect   (ipa) = 0.0
      csite%ebudget_zcaneffect   (ipa) = 0.0
      csite%ebudget_residual     (ipa) = 0.0
      csite%wbudget_precipgain   (ipa) = 0.0
      csite%wbudget_loss2atm     (ipa) = 0.0
      csite%wbudget_loss2runoff  (ipa) = 0.0
      csite%wbudget_loss2drainage(ipa) = 0.0
      csite%wbudget_denseffect   (ipa) = 0.0
      csite%wbudget_wcapeffect   (ipa) = 0.0
      csite%wbudget_zcaneffect   (ipa) = 0.0
      csite%wbudget_residual     (ipa) = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise the committed carbon pool.  In case this is a newly disturbed       !
      ! patch, it can skip the test but it must be initialised with the committed          !
      ! emissions.                                                                         !
      !------------------------------------------------------------------------------------!
      call reset_cbudget_committed(csite,ipa,.false.)
      !------------------------------------------------------------------------------------!


      !----- Compute current storage terms. -----------------------------------------------!
      csite%co2budget_initialstorage(ipa) = compute_co2_storage     (csite,ipa)
      csite%cbudget_initialstorage  (ipa) = compute_carbon_storage  (csite,ipa)
      csite%wbudget_initialstorage  (ipa) = compute_water_storage   (csite,lsl,ipa)
      csite%ebudget_initialstorage  (ipa) = compute_enthalpy_storage(csite,lsl,ipa)
      !------------------------------------------------------------------------------------!



      return
   end subroutine initial_patch_budget
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine resets the committed changes in carbon stocks (biomass and         !
   ! necromass).  This pool must be zero at the end of the day when vegetation dynamics is !
   ! on and after the metabolic NPP and heterotrophic respiration are accounted.  In case  !
   ! it is not zero, then stop the run, because carbon is not being conserved.             !
   !---------------------------------------------------------------------------------------!
   subroutine reset_cbudget_committed(csite,ipa,check_budget)
      use ed_state_vars, only : sitetype       & ! structure
                              , patchtype      ! ! structured
      use consts_coms  , only : kgCday_2_umols ! ! intent(in)
      use ed_misc_coms , only : current_time ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)   , target     :: csite
      integer          , intent(in) :: ipa
      logical          , intent(in) :: check_budget
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ico
      real                          :: today_gpp
      real                          :: today_leaf_resp
      real                          :: today_root_resp
      real                          :: today_het_resp
      real                          :: toler_committed
      real                          :: resid_committed
      logical                       :: committed_violation
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter :: fmti='(a,1x,i14)'
      character(len=13), parameter :: fmtf='(a,1x,es14.7)'
      character(len=27), parameter :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!


      !----- Alias for current patch. -----------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Before we flush the committed carbon to zero, we check that we are not         !
      ! violating carbon conservation.  We skip this check in case this is a simulation    !
      ! without dynamic vegetation (we cannot guarantee carbon conservation in this case,  !
      ! as carbon stocks are not updated in response to non-zero NEP.  We also skip the    !
      ! check when the patch is just created.  In both cases, we still initialise the      !
      ! committed carbon, and in case we are not checking conservation of the committed    !
      ! carbon, we must also update the initial storage to avoid crashes in the next sub-  !
      ! daily time step.                                                                   !
      !------------------------------------------------------------------------------------!
      if (check_budget) then
         !----- Make sure we are conserving carbon. ---------------------------------------!
         today_gpp          = 0.
         today_leaf_resp    = 0.
         today_root_resp    = 0.
         do ico=1,cpatch%ncohorts
            today_gpp       = today_gpp       + cpatch%today_gpp      (ico)
            today_leaf_resp = today_leaf_resp + cpatch%today_leaf_resp(ico)
            today_root_resp = today_root_resp + cpatch%today_root_resp(ico)
         end do
         today_gpp           = today_gpp           / kgCday_2_umols
         today_leaf_resp     = today_leaf_resp     / kgCday_2_umols
         today_root_resp     = today_root_resp     / kgCday_2_umols
         today_het_resp      = csite%today_rh(ipa) / kgCday_2_umols
         toler_committed     = tol_subday_budget                                           &
                             * max(today_gpp,today_leaf_resp,today_root_resp,today_het_resp)
         resid_committed     = csite%cbudget_committed(ipa) - today_gpp                    &
                             + today_leaf_resp + today_root_resp + today_het_resp
         committed_violation = abs(resid_committed) > toler_committed
         !---------------------------------------------------------------------------------!


         !----- Stop in case the committed pool is not flushed. ---------------------------!
         if (committed_violation) then
            write(unit=*,fmt='(a)')  '|==================================================|'
            write(unit=*,fmt='(a)')  '|==================================================|'
            write(unit=*,fmt='(a)')  '|    !!!   Committed carbon budget failed   !!!    |'
            write(unit=*,fmt='(a)')  '|--------------------------------------------------|'
            write(unit=*,fmt=fmtt )  ' TIME               : ',current_time%year            &
                                                             ,current_time%month           &
                                                             ,current_time%date            &
                                                             ,current_time%time
            write(unit=*,fmt=fmti )  ' PATCH              : ',ipa
            write(unit=*,fmt=fmti )  ' DIST_TYPE          : ',csite%dist_type(ipa)
            write(unit=*,fmt=fmtf )  ' AGE                : ',csite%age      (ipa)
            write(unit=*,fmt='(a)')  ' -------------------------------------------------- '
            write(unit=*,fmt=fmtf )  ' COMMITTED_CARBON   : ',csite%cbudget_committed(ipa)
            write(unit=*,fmt=fmtf )  ' GPP                : ',today_gpp
            write(unit=*,fmt=fmtf )  ' LEAF RESPIRATION   : ',today_leaf_resp
            write(unit=*,fmt=fmtf )  ' ROOT RESPIRATION   : ',today_root_resp
            write(unit=*,fmt=fmtf )  ' HETEROTROPHIC RESP : ',today_het_resp
            write(unit=*,fmt=fmtf )  ' RESIDUAL           : ',resid_committed
            write(unit=*,fmt=fmtf )  ' TOLERANCE          : ',toler_committed
            write(unit=*,fmt='(a)')  '|==================================================|'
            write(unit=*,fmt='(a)')  '|==================================================|'
            write(unit=*,fmt='(a)')  ' '


            call fatal_error('Budget check has failed, see message above.'                 &
                            ,'reset_cbudget_committed','budget_utils.f90')
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The committed pool exists because carbon fluxes in ecosystem dynamics are not  !
      ! updated every thermodynamics time step, so the carbon lost by the ecosystem        !
      ! (excluding canopy air space) is only transferred to the canopy air space during    !
      ! the following day.  This is the beginning of the day, so we initialise the pool    !
      ! with what is going to be lost as growth and storage respiration.  The committed    !
      ! pool will be updated throughout the day, as the respiration is actually released   !
      ! from the committed pool towards the atmosphere, and gross primary productivity is  !
      ! captures carbon that will only be assimilated by the ecosystem in the day after.   !
      ! This is the initial time, so we ADD respiration terms here because they will be    !
      ! released.                                                                          !
      !                                                                                    !
      !     In case vegetation dynamics is not active, this pool is never reset to account !
      ! for the fact that carbon storage does not change in spite of the NEP not being     !
      ! zero.                                                                              !
      !------------------------------------------------------------------------------------!
      csite%cbudget_committed(ipa) = csite%commit_storage_resp(ipa)                        &
                                   + csite%commit_growth_resp (ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we are not checking the conservation of the committed carbon, we must  !
      ! update the total storage, to avoid false alarms.  We cannot conserve carbon in the !
      ! long term when vegetation dynamics is disabled, because carbon stocks do not       !
      ! change even when NEP is non zero.                                                  !
      !------------------------------------------------------------------------------------!
      if (.not. check_budget) then
         csite%cbudget_initialstorage(ipa) = compute_carbon_storage(csite,ipa)
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine reset_cbudget_committed
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine resets the committed changes in carbon stocks (biomass and         !
   ! necromass).                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine update_cbudget_committed(csite,ipa)
     
      use ed_state_vars, only : sitetype     & ! structure
                              , patchtype    ! ! structure
      use ed_misc_coms , only : dtlsm        & ! intent(in)
                              , current_time ! ! intent(in)
      use consts_coms  , only : umol_2_kgC   & ! intent(in)
                              , day_sec      ! ! intent(in)
      implicit none

      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)   , target    :: csite
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)  , pointer   :: cpatch
      integer                      :: ipa
      integer                      :: jpa
      integer                      :: ico
      real                         :: dtlsm_o_daysec
      real                         :: umol_o_sec_2_kgC
      real                         :: bef_committed
      real                         :: step_delta
      real                         :: step_gpp
      real                         :: step_leaf_resp
      real                         :: step_root_resp
      real                         :: step_growth_resp
      real                         :: step_storage_resp
      real                         :: step_het_resp
      character(len=27)            :: committed_file
      !----- Local parameters, for debugging. ---------------------------------------------!
      logical          , save      :: first_time = .true.
      logical          , parameter :: print_step = .false.
      !------------------------------------------------------------------------------------!


      !----- Alias for current patch. -----------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !----- Alias for time step for growth and storage respiration. ----------------------!
      dtlsm_o_daysec   = dtlsm / day_sec
      umol_o_sec_2_kgC = umol_2_kgC * dtlsm
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Print header with time and variables.                                          !
      !------------------------------------------------------------------------------------!
      if (first_time) then
         first_time = .false.
         if (print_step) then
            do jpa=1,csite%npatches
               write(committed_file,fmt='(a,i4.4,a)') 'check_committed_ipa',jpa,'.txt'
               open(unit=59,file=committed_file,status='replace',action='write')
               write(unit=59,fmt='(13(a,1x))')                                             &
                     '  YEAR',' MONTH','   DAY','        TIME','  BEF_COMMIT'              &
                                ,'  NOW_COMMIT','       DELTA','         GPP'              &
                                ,'   LEAF_RESP','   ROOT_RESP','  STORE_RESP'              &
                                ,'   GROW_RESP','    HET_RESP'
               close(unit=59,status='keep')
            end do
         end if
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Loop through all cohorts, and update the virtual pool based on each cohort's   !
      ! NPP.                                                                               !
      !------------------------------------------------------------------------------------!
      bef_committed     = csite%cbudget_committed(ipa)
      step_gpp          = 0.
      step_leaf_resp    = 0.
      step_root_resp    = 0.
      step_growth_resp  = 0.
      step_storage_resp = 0.
      cohloop: do ico=1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !     Update committed pools based on metabolic NPP.  These variables are scaled  !
         ! by LAI and thus are extensive, no need to multiply by number density.           !
         !---------------------------------------------------------------------------------!
         step_gpp       = step_gpp       + cpatch%gpp(ico)
         step_leaf_resp = step_leaf_resp + cpatch%leaf_respiration(ico)
         step_root_resp = step_root_resp + cpatch%root_respiration(ico)
         !---------------------------------------------------------------------------------!
      end do cohloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Integrate the step.  Notice that metabolic terms are in umol/m2/s and storage  !
      ! and growth respiration are in kgC/m2/day.  We want everything to be in kgC/m2      !
      ! (integrated over dtlsm).                                                           !
      !------------------------------------------------------------------------------------!
      step_gpp          = step_gpp                       * umol_o_sec_2_kgC
      step_leaf_resp    = step_leaf_resp                 * umol_o_sec_2_kgC
      step_root_resp    = step_root_resp                 * umol_o_sec_2_kgC
      step_storage_resp = csite%commit_storage_resp(ipa) * dtlsm_o_daysec
      step_growth_resp  = csite%commit_growth_resp (ipa) * dtlsm_o_daysec
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Integrate committed change in carbon stocks due to heterotrophic respiration.  !
      !------------------------------------------------------------------------------------!
      step_het_resp = csite%rh(ipa) * umol_o_sec_2_kgC
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the committed carbon pool.  Add assimilation and subtract respiration.  !
      !------------------------------------------------------------------------------------!
      step_delta                   = step_gpp          - step_leaf_resp   - step_root_resp &
                                   - step_storage_resp - step_growth_resp - step_het_resp 
      csite%cbudget_committed(ipa) = bef_committed     + step_delta
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Print step in case it's needed.                                                !
      !------------------------------------------------------------------------------------!
      if (print_step) then
         write(committed_file,fmt='(a,i4.4,a)') 'check_committed_ipa',ipa,'.txt'
         open(unit=59,file=committed_file,status='old',position='append',action='write')
         write(unit=59,fmt='(3(i6,1x),f12.1,9(1x,es12.5))')                                &
                 current_time%year,current_time%month,current_time%date,current_time%time  &
                ,bef_committed,csite%cbudget_committed(ipa),step_delta,step_gpp            &
                ,step_leaf_resp,step_root_resp,step_storage_resp,step_growth_resp          &
                ,step_het_resp
         close(unit=59,status='keep')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_cbudget_committed
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine integrates the changes in storage after accounting for sources and !
   ! sinks, and checks whether the model is conserving carbon, carbon dioxide, enthalpy,   !
   ! and water.  By default, ED-2 issues a fatal error in case it detects any violation to !
   ! conservation.  In the future this should be extended to nitrogen pools.               !
   !                                                                                       !
   !    This substoutine now accounts for two effects that affect storage when vegetation  !
   ! dynamics is turned on.  We subtract these effects from the initial storage so we      !
   ! include them in the budget check at the end of the time step.  These additional terms !
   ! are:                                                                                  !
   ! hcapeffect -- change in total enthalpy caused by changes in vegetation heat capacity  !
   !               (growth or mortality).                                                  !
   ! zcaneffect -- change in total storage of all state variables (enthalpy, water,        !
   !               carbon, and CO2) due to changes in canopy air space depth (also linked  !
   !               to growth and mortality).                                               !
   !---------------------------------------------------------------------------------------!
   subroutine compute_budget(csite,lsl,pcpg,qpcpg,ipa,wcurr_loss2atm,ecurr_netrad          &
                            ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage           &
                            ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff       &
                            ,site_area,cbudget_nep,old_can_prss,old_can_enthalpy           &
                            ,old_can_temp,old_can_shv,old_can_co2,old_can_rhos             &
                            ,old_can_dmol,mid_can_rhos,mid_can_dmol)
      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use ed_max_dims  , only : str_len            ! ! intent(in)
      use ed_misc_coms , only : dtlsm              & ! intent(in)
                              , frqsum             & ! intent(in)
                              , current_time       ! ! intent(in)
      use ed_max_dims  , only : n_dbh              ! ! intent(in)
      use consts_coms  , only : umol_2_kgC         & ! intent(in)
                              , day_sec            ! ! intent(in)
      use rk4_coms     , only : print_budget       & ! intent(in)
                              , budget_pref        & ! intent(in)
                              , checkbudget        ! ! intent(in)
      use therm_lib    , only : tq2enthalpy        ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target        :: csite
      real                  , intent(inout) :: pcpg
      real                  , intent(inout) :: qpcpg
      real                  , intent(inout) :: co2curr_loss2atm
      real                  , intent(inout) :: ecurr_netrad
      real                  , intent(inout) :: ecurr_loss2atm
      real                  , intent(inout) :: ecurr_loss2drainage
      real                  , intent(inout) :: ecurr_loss2runoff
      real                  , intent(inout) :: wcurr_loss2atm
      real                  , intent(inout) :: wcurr_loss2drainage
      real                  , intent(inout) :: wcurr_loss2runoff
      integer               , intent(in)    :: lsl
      integer               , intent(in)    :: ipa
      real                  , intent(in)    :: site_area
      real                  , intent(inout) :: cbudget_nep
      real                  , intent(in)    :: old_can_prss
      real                  , intent(in)    :: old_can_enthalpy
      real                  , intent(in)    :: old_can_temp
      real                  , intent(in)    :: old_can_shv
      real                  , intent(in)    :: old_can_co2
      real                  , intent(in)    :: old_can_rhos
      real                  , intent(in)    :: old_can_dmol
      real                  , intent(in)    :: mid_can_rhos
      real                  , intent(in)    :: mid_can_dmol

      !----- Local variables --------------------------------------------------------------!
      type(patchtype)       , pointer       :: cpatch
      character(len=str_len)                :: budget_fout
      real                                  :: co2budget_initialstorage
      real                                  :: co2budget_finalstorage
      real                                  :: co2budget_deltastorage
      real                                  :: co2budget_tolerance
      real                                  :: co2budget_scale
      real                                  :: co2curr_gpp
      real                                  :: co2curr_leafresp
      real                                  :: co2curr_rootresp
      real                                  :: co2curr_storageresp
      real                                  :: co2curr_growthresp
      real                                  :: co2curr_hetresp
      real                                  :: co2curr_nee
      real                                  :: co2curr_nep
      real                                  :: co2curr_denseffect
      real                                  :: co2curr_zcaneffect
      real                                  :: co2curr_residual
      real                                  :: cbudget_initialstorage
      real                                  :: cbudget_committed
      real                                  :: cbudget_finalstorage
      real                                  :: cbudget_deltastorage
      real                                  :: cbudget_tolerance
      real                                  :: cbudget_scale
      real                                  :: ccurr_loss2atm
      real                                  :: ccurr_denseffect
      real                                  :: ccurr_zcaneffect
      real                                  :: ccurr_loss2yield
      real                                  :: ccurr_seedrain
      real                                  :: ccurr_residual
      real                                  :: ebudget_initialstorage
      real                                  :: ebudget_finalstorage
      real                                  :: ebudget_deltastorage
      real                                  :: ebudget_tolerance
      real                                  :: ebudget_scale
      real                                  :: ecurr_precipgain
      real                                  :: ecurr_denseffect
      real                                  :: ecurr_prsseffect
      real                                  :: ecurr_hcapeffect
      real                                  :: ecurr_wcapeffect
      real                                  :: ecurr_zcaneffect
      real                                  :: ecurr_residual
      real                                  :: wbudget_initialstorage
      real                                  :: wbudget_finalstorage
      real                                  :: wbudget_deltastorage
      real                                  :: wbudget_tolerance
      real                                  :: wbudget_scale
      real                                  :: wcurr_precipgain
      real                                  :: wcurr_denseffect
      real                                  :: wcurr_wcapeffect
      real                                  :: wcurr_zcaneffect
      real                                  :: wcurr_residual
      real                                  :: curr_can_enthalpy
      real                                  :: gpp
      real                                  :: leaf_resp
      real                                  :: root_resp
      real                                  :: storage_resp
      real                                  :: growth_resp
      real                                  :: co2_factor
      real                                  :: crb_factor
      real                                  :: ent_factor
      real                                  :: h2o_factor
      real                                  :: patch_lai
      real                                  :: patch_wai
      integer                               :: jpa
      integer                               :: ico
      logical                               :: isthere
      logical                               :: budget_fine
      logical                               :: co2_fine
      logical                               :: carbon_fine
      logical                               :: enthalpy_fine
      logical                               :: water_fine
      !----- Local constants. -------------------------------------------------------------!
      character(len=13)     , parameter     :: fmtf='(a,1x,es14.7)'
      character(len= 9)     , parameter     :: fmti='(a,1x,i7)'
      character(len= 9)     , parameter     :: fmtl='(a,1x,l1)'
      character(len=10)     , parameter     :: bhfmt='(46(a,1x))'
      character(len=48)     , parameter     :: bbfmt='(3(i14,1x),43(es14.7,1x))'
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
               write(unit=86,fmt=bhfmt)   '          YEAR' , '         MONTH'              &
                                        , '           DAY' , '          TIME'              &
                                        , '           LAI' , '           WAI'              &
                                        , '        HEIGHT' , '   CO2.STORAGE'              &
                                        , '  CO2.RESIDUAL' , '  CO2.DSTORAGE'              &
                                        , '       CO2.NEP' , '  CO2.DENS.EFF'              &
                                        , '  CO2.ZCAN.EFF' , '  CO2.LOSS2ATM'              &
                                        , '   CRB.STORAGE' , ' CRB.COMMITTED'              &
                                        , '  CRB.RESIDUAL' , '  CRB.DSTORAGE'              &
                                        , '  CRB.DENS.EFF' , '  CRB.ZCAN.EFF'              &
                                        , '  CRB.SEEDRAIN' , 'CRB.LOSS2YIELD'              &
                                        , '  CRB.LOSS2ATM' , '   ENT.STORAGE'              &
                                        , '  ENT.RESIDUAL' , '  ENT.DSTORAGE'              &
                                        , '    ENT.PRECIP' , '    ENT.NETRAD'              &
                                        , '  ENT.DENS.EFF' , '  ENT.PRSS.EFF'              &
                                        , '  ENT.HCAP.EFF' , '  ENT.WCAP.EFF'              &
                                        , '  ENT.ZCAN.EFF' , '  ENT.LOSS2ATM'              &
                                        , '  ENT.DRAINAGE' , '    ENT.RUNOFF'              &
                                        , '   H2O.STORAGE' , '  H2O.RESIDUAL'              &
                                        , '  H2O.DSTORAGE' , '    H2O.PRECIP'              &
                                        , '  H2O.DENS.EFF' , '  H2O.WCAP.EFF'              &
                                        , '  H2O.ZCAN.EFF' , '  H2O.LOSS2ATM'              &
                                        , '  H2O.DRAINAGE' , '    H2O.RUNOFF'
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
      !     Compute the density effect.  This term corrects for the fact that air expands  !
      ! or contracts as pressure or temperature changes.  Because we impose constant CAS   !
      ! depth, some mass may leak in or out of this somewhat arbitrary layer.  For any     !
      ! extensive [X/m2] property X [X/m2], defined as X = x*z*rho, where x is the         !
      ! intensive property [X/kg or X/mol]], and z is the canopy air space depth, the      !
      ! budget can be written as:                                                          !
      !                                                                                    !
      !    dX                  d(rho)                                                      !
      !   ---- = I - L + x*z ----------                                                    !
      !    dt                    dt                                                        !
      !                                                                                    !
      ! where I is the input flux, L is the loss flux, and rho is the density [kg/m3 or    !
      ! mol/m3].  The density effect is the third term of the right-hand side.             !
      !------------------------------------------------------------------------------------!
      !------ CO2. ------------------------------------------------------------------------!
      co2curr_denseffect  = ddens_dt_effect(csite%can_depth(ipa),old_can_dmol,mid_can_dmol &
                                           ,csite%can_dmol(ipa),old_can_co2                &
                                           ,csite%can_co2(ipa) )
      !------ Carbon.  Derive it from CO2. ------------------------------------------------!
      ccurr_denseffect    = co2curr_denseffect * umol_2_kgC
      !------ Water. ----------------------------------------------------------------------!
      wcurr_denseffect    = ddens_dt_effect(csite%can_depth(ipa),old_can_rhos,mid_can_rhos &
                                           ,csite%can_rhos(ipa),old_can_shv                &
                                           ,csite%can_shv(ipa))
      !------ Enthalpy.  ------------------------------------------------------------------!
      curr_can_enthalpy   = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
      ecurr_denseffect    = ddens_dt_effect(csite%can_depth(ipa),old_can_rhos,mid_can_rhos &
                                           ,csite%can_rhos(ipa),old_can_enthalpy           &
                                           ,curr_can_enthalpy)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    For the specific case of enthalpy, we also compute the pressure effect between  !
      ! time steps.  We cannot guarantee conservation of enthalpy when we update pressure, !
      ! because of the first law of thermodynamics (the way to address this would be to    !
      ! use equivalent potential temperature, which is enthalpy plus pressure effect).     !
      ! Enthalpy is preserved within one time step, once pressure is updated and remains   !
      ! constant.                                                                          !
      !                                                                                    !
      !   dH             dp                                                                !
      !  ---- = Q + z * ----                                                               !
      !   dt             dt                                                                !
      !                                                                                    !
      ! where p is the canopy air space pressure, Q is the net heat exchange, z is the     !
      ! volume per unit area (aka depth) of the canopy air space.                          !
      !------------------------------------------------------------------------------------!
      ecurr_prsseffect    = csite%can_depth(ipa) * (csite%can_prss(ipa) - old_can_prss)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Changes due to the vegetation dynamics.  These are applied instantaneously,    !
      ! but to make their scale consistent with other effects, we report them as fluxes,   !
      ! hence the frqsum term.  Beware that they may be lagged in the history output (-S-  !
      ! files), which shouldn't be used for any research analysis anyway.                  !
      !------------------------------------------------------------------------------------!
      co2curr_zcaneffect = csite%co2budget_zcaneffect(ipa) * frqsum
      ccurr_zcaneffect   = csite%cbudget_zcaneffect  (ipa) * frqsum
      ccurr_loss2yield   = csite%cbudget_loss2yield  (ipa) * frqsum
      ccurr_seedrain     = csite%cbudget_seedrain    (ipa) * frqsum
      wcurr_wcapeffect   = csite%wbudget_wcapeffect  (ipa) * frqsum
      wcurr_zcaneffect   = csite%wbudget_zcaneffect  (ipa) * frqsum
      ecurr_hcapeffect   = csite%ebudget_hcapeffect  (ipa) * frqsum
      ecurr_wcapeffect   = csite%ebudget_wcapeffect  (ipa) * frqsum
      ecurr_zcaneffect   = csite%ebudget_zcaneffect  (ipa) * frqsum
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute the carbon dioxide flux components.                                    !
      !------------------------------------------------------------------------------------!
      call sum_plant_cfluxes(csite,ipa, gpp, leaf_resp,root_resp,storage_resp,growth_resp)
      co2curr_gpp         = gpp           * dtlsm
      co2curr_leafresp    = leaf_resp     * dtlsm
      co2curr_rootresp    = root_resp     * dtlsm
      co2curr_storageresp = storage_resp  * dtlsm
      co2curr_growthresp  = growth_resp   * dtlsm
      co2curr_hetresp     = csite%rh(ipa) * dtlsm
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Compute the NEE contribution to changes CO2 storage.  NEP is used for         !
      ! long-term, stand-scale budget tracking of carbon pools excluding canopy air space, !
      ! and for the CO2-only budget.  For the short-term, patch-scale total carbon budget  !
      ! we use the loss to atmosphere instead, because we also account for carbon storage  !
      ! in the canopy air space.                                                           !
      !------------------------------------------------------------------------------------!
      co2curr_nee    = co2curr_leafresp    + co2curr_rootresp    + co2curr_storageresp     &
                     + co2curr_growthresp  + co2curr_hetresp     - co2curr_gpp
      co2curr_nep    = - co2curr_nee
      cbudget_nep    = cbudget_nep + site_area * csite%area(ipa) * co2curr_nep * umol_2_kgC
      !----- Leverage the CO2 budget loss term to find the carbon equivalent. -------------!
      ccurr_loss2atm = co2curr_loss2atm * umol_2_kgC
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the initial storage.  Because of patch fusion, we must account for the    !
      ! vegetation dynamics internally, and thus the initial storage must exclude recent   !
      ! additions/losses due to vegetation dynamics.  Note that zcaneffect, hcapeffect,    !
      ! and seedrain also appear in the residual calculation, and this is a way to make    !
      ! sure we account for these effects when checking the conservation of the state      !
      ! variables.                                                                         !
      !------------------------------------------------------------------------------------!
      co2budget_initialstorage = csite%co2budget_initialstorage(ipa)
      cbudget_initialstorage   = csite%cbudget_initialstorage  (ipa)
      cbudget_committed        = csite%cbudget_committed       (ipa)
      wbudget_initialstorage   = csite%wbudget_initialstorage  (ipa)
      ebudget_initialstorage   = csite%ebudget_initialstorage  (ipa)
      !------------------------------------------------------------------------------------!



      !----- Compute current storage terms. -----------------------------------------------!
      co2budget_finalstorage = compute_co2_storage     (csite,ipa)
      cbudget_finalstorage   = compute_carbon_storage  (csite,ipa)
      wbudget_finalstorage   = compute_water_storage   (csite,lsl,ipa)
      ebudget_finalstorage   = compute_enthalpy_storage(csite,lsl,ipa)
      !------------------------------------------------------------------------------------!



      !----- Compute the change in storage. -----------------------------------------------!
      co2budget_deltastorage = co2budget_finalstorage - co2budget_initialstorage
      cbudget_deltastorage   = cbudget_finalstorage   - cbudget_initialstorage
      wbudget_deltastorage   = wbudget_finalstorage   - wbudget_initialstorage
      ebudget_deltastorage   = ebudget_finalstorage   - ebudget_initialstorage
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Compute residuals.                                                             !
      !------------------------------------------------------------------------------------!
      !----- 1. Canopy CO2. ---------------------------------------------------------------!
      co2curr_residual = co2budget_deltastorage                                            &
                       - ( co2curr_nee        - co2curr_loss2atm                           &
                         + co2curr_denseffect + co2curr_zcaneffect )
      !----- 2. Carbon. -------------------------------------------------------------------!
      ccurr_residual   = cbudget_deltastorage                                              &
                       - ( ccurr_seedrain   - ccurr_loss2atm   - ccurr_loss2yield          &
                         + ccurr_denseffect + ccurr_zcaneffect )
      !----- 3. Energy. -------------------------------------------------------------------!
      ecurr_residual   = ebudget_deltastorage                                              &
                       - ( ecurr_precipgain    - ecurr_loss2atm    - ecurr_loss2drainage   &
                         - ecurr_loss2runoff   + ecurr_netrad      + ecurr_prsseffect      &
                         + ecurr_denseffect    + ecurr_hcapeffect  + ecurr_wcapeffect      &
                         + ecurr_zcaneffect    )
      !----- 4. Water. --------------------------------------------------------------------!
      wcurr_residual   = wbudget_deltastorage                                              &
                       - ( wcurr_precipgain    - wcurr_loss2atm                            &
                         - wcurr_loss2drainage - wcurr_loss2runoff                         &
                         + wcurr_denseffect    + wcurr_wcapeffect  + wcurr_zcaneffect  )
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Compute scale: we use the combination of the derivatives, as we are more       !
      ! the change may be small when large numbers with opposite signs are combined.                                                !
      !------------------------------------------------------------------------------------!
      !----- 1. Canopy CO2. ---------------------------------------------------------------!
      co2budget_scale  = max( abs(co2budget_finalstorage  )                                &
                            , abs(co2budget_initialstorage) )                              &
                       + max( abs(co2budget_deltastorage  )                                &
                            , abs(co2curr_leafresp        )                                &
                            , abs(co2curr_rootresp        )                                &
                            , abs(co2curr_storageresp     )                                &
                            , abs(co2curr_growthresp      )                                &
                            , abs(co2curr_hetresp         )                                &
                            , abs(co2curr_gpp             )                                &
                            , abs(co2curr_denseffect      )                                &
                            , abs(co2curr_zcaneffect      ) )
      !----- 2. Carbon. -------------------------------------------------------------------!
      cbudget_scale    = max( abs(cbudget_initialstorage  )                                &
                            , abs(cbudget_finalstorage    ) )                              &
                       + max( abs(cbudget_deltastorage    )                                &
                            , abs(ccurr_seedrain          )                                &
                            , abs(ccurr_loss2atm          )                                &
                            , abs(ccurr_loss2yield        )                                &
                            , abs(ccurr_denseffect        )                                &
                            , abs(ccurr_zcaneffect        ) )
      !----- 3. Energy. -------------------------------------------------------------------!
      ebudget_scale    = max( abs(ebudget_initialstorage  )                                &
                            , abs(ebudget_finalstorage    ) )                              &
                       + max( abs(ebudget_deltastorage    )                                &
                            , abs(ecurr_precipgain        )                                &
                            , abs(ecurr_loss2atm          )                                &
                            , abs(ecurr_loss2drainage     )                                &
                            , abs(ecurr_loss2runoff       )                                &
                            , abs(ecurr_netrad            )                                &
                            , abs(ecurr_prsseffect        )                                &
                            , abs(ecurr_denseffect        )                                &
                            , abs(ecurr_hcapeffect        )                                &
                            , abs(ecurr_wcapeffect        )                                &
                            , abs(ecurr_zcaneffect        ) )
      !----- 4. Water. --------------------------------------------------------------------!
      wbudget_scale    = max( abs(wbudget_initialstorage  )                                &
                            , abs(wbudget_finalstorage    ) )                              &
                       + max( abs(wbudget_deltastorage    )                                &
                            , abs(wcurr_precipgain        )                                &
                            , abs(wcurr_loss2atm          )                                &
                            , abs(wcurr_loss2drainage     )                                &
                            , abs(wcurr_loss2runoff       )                                &
                            , abs(wcurr_denseffect        )                                &
                            , abs(wcurr_wcapeffect        )                                &
                            , abs(wcurr_zcaneffect        ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Integrate residuals.                                                           !
      !------------------------------------------------------------------------------------!
      !----- 1. Canopy CO2. ---------------------------------------------------------------!
      csite%co2budget_residual(ipa) = csite%co2budget_residual(ipa)  + co2curr_residual
      !----- 2. Carbon. -------------------------------------------------------------------!
      csite%cbudget_residual  (ipa) = csite%cbudget_residual(ipa)    + ccurr_residual
      !----- 3. Energy. -------------------------------------------------------------------!
      csite%ebudget_residual  (ipa) = csite%ebudget_residual(ipa)    + ecurr_residual
      !----- 4. Water. --------------------------------------------------------------------!
      csite%wbudget_residual  (ipa) = csite%wbudget_residual(ipa)    + wcurr_residual
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Integrate the terms that are part of the budget.                                !
      !------------------------------------------------------------------------------------!
      !----- 1. Carbon dioxide. -----------------------------------------------------------!
      csite%co2budget_gpp        (ipa) = csite%co2budget_gpp(ipa)       + gpp       *dtlsm
      csite%co2budget_plresp     (ipa) = csite%co2budget_plresp(ipa)                       &
                                       + ( leaf_resp          + root_resp                  &
                                         + storage_resp       + growth_resp        )       &
                                       * dtlsm
      csite%co2budget_rh         (ipa) = csite%co2budget_rh(ipa)                           &
                                       + csite%rh(ipa) * dtlsm
      csite%co2budget_denseffect (ipa) = csite%co2budget_denseffect(ipa)                   &
                                       + co2curr_denseffect
      csite%co2budget_loss2atm   (ipa) = csite%co2budget_loss2atm(ipa)                     &
                                       + co2curr_loss2atm
      !----- 2. Carbon. -------------------------------------------------------------------!
      csite%cbudget_denseffect   (ipa) = csite%cbudget_denseffect(ipa)                     &
                                       + ccurr_denseffect
      csite%cbudget_loss2atm     (ipa) = csite%cbudget_loss2atm(ipa)                       &
                                       + ccurr_loss2atm
      !----- 3. Energy. -------------------------------------------------------------------!
      csite%ebudget_precipgain   (ipa) = csite%ebudget_precipgain(ipa)   + ecurr_precipgain
      csite%ebudget_netrad       (ipa) = csite%ebudget_netrad    (ipa)   + ecurr_netrad
      csite%ebudget_prsseffect   (ipa) = csite%ebudget_prsseffect(ipa)   + ecurr_prsseffect
      csite%ebudget_denseffect   (ipa) = csite%ebudget_denseffect(ipa)   + ecurr_denseffect
      csite%ebudget_loss2atm     (ipa) = csite%ebudget_loss2atm  (ipa)   + ecurr_loss2atm
      csite%ebudget_loss2drainage(ipa) = csite%ebudget_loss2drainage(ipa)                  &
                                       + ecurr_loss2drainage
      csite%ebudget_loss2runoff  (ipa) = csite%ebudget_loss2runoff(ipa)                    &
                                       + ecurr_loss2runoff
      !----- 4. Water. --------------------------------------------------------------------!
      csite%wbudget_precipgain   (ipa) = csite%wbudget_precipgain(ipa) + wcurr_precipgain
      csite%wbudget_denseffect   (ipa) = csite%wbudget_denseffect(ipa) + wcurr_denseffect
      csite%wbudget_loss2atm     (ipa) = csite%wbudget_loss2atm  (ipa) + wcurr_loss2atm
      csite%wbudget_loss2drainage(ipa) = csite%wbudget_loss2drainage(ipa)                  &
                                       + wcurr_loss2drainage
      csite%wbudget_loss2runoff  (ipa) = csite%wbudget_loss2runoff(ipa)                    &
                                       + wcurr_loss2runoff
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Update density and initial storage for next step.                              !
      !------------------------------------------------------------------------------------!
      csite%co2budget_initialstorage(ipa) = co2budget_finalstorage
      csite%cbudget_initialstorage  (ipa) = cbudget_finalstorage
      csite%wbudget_initialstorage  (ipa) = wbudget_finalstorage
      csite%ebudget_initialstorage  (ipa) = ebudget_finalstorage
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Flush the vegetation dynamics terms, they only need to be accounted for one     !
      ! time step.                                                                         !
      !------------------------------------------------------------------------------------!
      csite%co2budget_zcaneffect(ipa) = 0.0
      csite%cbudget_zcaneffect  (ipa) = 0.0
      csite%cbudget_loss2yield  (ipa) = 0.0
      csite%cbudget_seedrain    (ipa) = 0.0
      csite%ebudget_hcapeffect  (ipa) = 0.0
      csite%ebudget_wcapeffect  (ipa) = 0.0
      csite%ebudget_zcaneffect  (ipa) = 0.0
      csite%wbudget_wcapeffect  (ipa) = 0.0
      csite%wbudget_zcaneffect  (ipa) = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    If the "check budget" option is activated (you can turn on and turn off by set- !
      ! ting checkbudget in ed_params.f90), then the model will crash whenever there is    !
      ! some significant leak of CO2, water, or energy.                                    !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         !----- Look for violation of conservation in all quantities. ---------------------!
         co2budget_tolerance = tol_subday_budget * co2budget_scale
         cbudget_tolerance   = tol_carbon_budget * cbudget_scale
         ebudget_tolerance   = tol_subday_budget * ebudget_scale
         wbudget_tolerance   = tol_subday_budget * wbudget_scale
         !---------------------------------------------------------------------------------!



         !----- Look for violation of conservation in all quantities. ---------------------!
         co2_fine      = abs(co2curr_residual) <= co2budget_tolerance
         carbon_fine   = abs(ccurr_residual)   <= cbudget_tolerance
         enthalpy_fine = abs(ecurr_residual)   <= ebudget_tolerance
         water_fine    = abs(wcurr_residual)   <= wbudget_tolerance
         budget_fine   = co2_fine .and. carbon_fine .and. enthalpy_fine .and. water_fine
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the patch LAI and WAI for output, only if it is needed.                !
         !---------------------------------------------------------------------------------!
         if ( (.not. budget_fine) .or. print_budget ) then
            cpatch => csite%patch(ipa)
            patch_lai = 0.0
            patch_wai = 0.0
            do ico=1,cpatch%ncohorts
               patch_lai = patch_lai + cpatch%lai(ico)
               patch_wai = patch_wai + cpatch%wai(ico)
            end do
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !       Report all the budget terms in case any of the budget checks have failed. !
         !---------------------------------------------------------------------------------!
         if (.not. budget_fine) then
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') '|       !!!    Sub-daily budget failed    !!!      |'
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : ',         &
               current_time%year,current_time%month,current_time%date ,current_time%time
            write (unit=*,fmt=fmti ) ' IPA               : ',ipa
            write (unit=*,fmt=fmti ) ' DIST_TYPE         : ',csite%dist_type(ipa)
            write (unit=*,fmt=fmtf ) ' AGE               : ',csite%age(ipa)
            write (unit=*,fmt=fmtf ) ' LAI               : ',patch_lai
            write (unit=*,fmt=fmtf ) ' WAI               : ',patch_wai
            write (unit=*,fmt=fmtf ) ' VEG_HEIGHT        : ',csite%veg_height(ipa)
            write (unit=*,fmt=fmtf ) ' CAN_DEPTH         : ',csite%can_depth(ipa)
            write (unit=*,fmt=fmtf ) ' OLD_CAN_PRSS      : ',old_can_prss
            write (unit=*,fmt=fmtf ) ' CAN_PRSS          : ',csite%can_prss(ipa)
            write (unit=*,fmt=fmtf ) ' OLD_CAN_ENTHALPY  : ',old_can_enthalpy
            write (unit=*,fmt=fmtf ) ' CAN_ENTHALPY      : ',curr_can_enthalpy
            write (unit=*,fmt=fmtf ) ' OLD_CAN_TEMP      : ',old_can_temp
            write (unit=*,fmt=fmtf ) ' CAN_TEMP          : ',csite%can_temp(ipa)
            write (unit=*,fmt=fmtf ) ' CAN_SHV           : ',csite%can_shv (ipa)
            write (unit=*,fmt=fmtf ) ' OLD_CAN_SHV       : ',old_can_shv
            write (unit=*,fmt=fmtf ) ' CAN_CO2           : ',csite%can_co2 (ipa)
            write (unit=*,fmt=fmtf ) ' OLD_CAN_CO2       : ',old_can_co2
            write (unit=*,fmt=fmtf ) ' OLD_CAN_RHOS      : ',old_can_rhos
            write (unit=*,fmt=fmtf ) ' MID_CAN_RHOS      : ',mid_can_rhos
            write (unit=*,fmt=fmtf ) ' CAN_RHOS          : ',csite%can_rhos(ipa)
            write (unit=*,fmt=fmtf ) ' OLD_CAN_DMOL      : ',old_can_dmol
            write (unit=*,fmt=fmtf ) ' MID_CAN_DMOL      : ',mid_can_dmol
            write (unit=*,fmt=fmtf ) ' CAN_DMOL          : ',csite%can_dmol(ipa)
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '  Summary'
            write (unit=*,fmt='(a)') ' .................................................. '
            write (unit=*,fmt=fmtl ) ' CO2_FINE          : ',co2_fine
            write (unit=*,fmt=fmtl ) ' CARBON_FINE       : ',carbon_fine
            write (unit=*,fmt=fmtl ) ' ENTHALPY_FINE     : ',enthalpy_fine
            write (unit=*,fmt=fmtl ) ' WATER_FINE        : ',water_fine
            write (unit=*,fmt=fmtf ) ' REL_TOLERANCE     : ',tol_subday_budget
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '  CO2 Budget'
            write (unit=*,fmt='(a)') ' .................................................. '
            write (unit=*,fmt=fmtf ) ' TOLERANCE         : ',co2budget_tolerance
            write (unit=*,fmt=fmtf ) ' RESIDUAL          : ',co2curr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE   : ',co2budget_initialstorage
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE     : ',co2budget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE     : ',co2budget_deltastorage
            write (unit=*,fmt=fmtf ) ' GPP               : ',co2curr_gpp
            write (unit=*,fmt=fmtf ) ' LEAF_RESP         : ',co2curr_leafresp
            write (unit=*,fmt=fmtf ) ' ROOT_RESP         : ',co2curr_rootresp
            write (unit=*,fmt=fmtf ) ' STORAGE_RESP      : ',co2curr_storageresp
            write (unit=*,fmt=fmtf ) ' GROWTH_RESP       : ',co2curr_growthresp
            write (unit=*,fmt=fmtf ) ' HET_RESP          : ',co2curr_hetresp
            write (unit=*,fmt=fmtf ) ' NEP               : ',co2curr_nep
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT    : ',co2curr_denseffect
            write (unit=*,fmt=fmtf ) ' CANDEPTH_EFFECT   : ',co2curr_zcaneffect
            write (unit=*,fmt=fmtf ) ' LOSS2ATM          : ',co2curr_loss2atm
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '  Carbon budget'
            write (unit=*,fmt='(a)') ' .................................................. '
            write (unit=*,fmt=fmtf ) ' TOLERANCE         : ',cbudget_tolerance
            write (unit=*,fmt=fmtf ) ' RESIDUAL          : ',ccurr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE   : ',cbudget_initialstorage
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE     : ',cbudget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE     : ',cbudget_deltastorage
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT    : ',ccurr_denseffect
            write (unit=*,fmt=fmtf ) ' SEEDRAIN          : ',ccurr_seedrain
            write (unit=*,fmt=fmtf ) ' CANDEPTH_EFFECT   : ',ccurr_zcaneffect
            write (unit=*,fmt=fmtf ) ' LOSS2YIELD        : ',ccurr_loss2yield
            write (unit=*,fmt=fmtf ) ' LOSS2ATM          : ',ccurr_loss2atm
            write (unit=*,fmt=fmtf ) ' COMMITTED         : ',cbudget_committed
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '  Enthalpy budget'
            write (unit=*,fmt='(a)') ' .................................................. '
            write (unit=*,fmt=fmtf ) ' TOLERANCE       : ',ebudget_tolerance
            write (unit=*,fmt=fmtf ) ' RESIDUAL        : ',ecurr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE : ',ebudget_initialstorage
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE   : ',ebudget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE   : ',ebudget_deltastorage
            write (unit=*,fmt=fmtf ) ' PRECIPGAIN      : ',ecurr_precipgain
            write (unit=*,fmt=fmtf ) ' NETRAD          : ',ecurr_netrad
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT  : ',ecurr_denseffect
            write (unit=*,fmt=fmtf ) ' PRESSURE_EFFECT : ',ecurr_prsseffect
            write (unit=*,fmt=fmtf ) ' VEG_HCAP_EFFECT : ',ecurr_hcapeffect
            write (unit=*,fmt=fmtf ) ' CAPACITY_EFFECT : ',ecurr_wcapeffect
            write (unit=*,fmt=fmtf ) ' CANDEPTH_EFFECT : ',ecurr_zcaneffect
            write (unit=*,fmt=fmtf ) ' LOSS2ATM        : ',ecurr_loss2atm
            write (unit=*,fmt=fmtf ) ' LOSS2DRAINAGE   : ',ecurr_loss2drainage
            write (unit=*,fmt=fmtf ) ' LOSS2RUNOFF     : ',ecurr_loss2runoff
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '  Water budget'
            write (unit=*,fmt='(a)') ' .................................................. '
            write (unit=*,fmt=fmtf ) ' TOLERANCE       : ',wbudget_tolerance
            write (unit=*,fmt=fmtf ) ' RESIDUAL        : ',wcurr_residual
            write (unit=*,fmt=fmtf ) ' INITIAL_STORAGE : ',wbudget_initialstorage
            write (unit=*,fmt=fmtf ) ' FINAL_STORAGE   : ',wbudget_finalstorage
            write (unit=*,fmt=fmtf ) ' DELTA_STORAGE   : ',wbudget_deltastorage
            write (unit=*,fmt=fmtf ) ' PRECIPGAIN      : ',wcurr_precipgain
            write (unit=*,fmt=fmtf ) ' DENSITY_EFFECT  : ',wcurr_denseffect
            write (unit=*,fmt=fmtf ) ' CAPACITY_EFFECT : ',wcurr_wcapeffect
            write (unit=*,fmt=fmtf ) ' CANDEPTH_EFFECT : ',wcurr_zcaneffect
            write (unit=*,fmt=fmtf ) ' LOSS2ATM        : ',wcurr_loss2atm
            write (unit=*,fmt=fmtf ) ' LOSS2DRAINAGE   : ',wcurr_loss2drainage
            write (unit=*,fmt=fmtf ) ' LOSS2RUNOFF     : ',wcurr_loss2runoff
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') ' '
            write (unit=*,fmt='(a)') '  A note on budget check'
            write (unit=*,fmt='(a)') ' .................................................. '
            write (unit=*,fmt='(a)') '    Budget failure doesn''t necessarily mean a'      &
                                  // ' problem in the RK4 or Euler integrators.  If you'
            write (unit=*,fmt='(a)') ' see NaN in any variable above, and the simulation'  &
                                  // ' time is the first day of the month near 00UTC,'
            write (unit=*,fmt='(a)') ' then it is very likely that some variable has'      &
                                  // ' not been properly initialised in the patch or'
            write (unit=*,fmt='(a)') ' cohort dynamics (e.g. new recruit or new patch).'   &
                                  // '  The best way to spot the error is to compile the'
            write (unit=*,fmt='(a)') ' model with strict debugging options and checks.'
            write (unit=*,fmt='(a)') '     Problems that occur in the middle of the day'   &
                                  // ' and middle of the month more are more likely to'
            write (unit=*,fmt='(a)') ' be true budget violations at sub-daily steps (i.e.' &
                                  // ' thermodynamics, photosynthesis, radiation).'
            write (unit=*,fmt='(a)') '|--------------------------------------------------|'
            write (unit=*,fmt='(a)') ' '
         end if
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         !    Check whether to print a detailed report on all the budget terms.            !
         !---------------------------------------------------------------------------------!
         if (print_budget) then 
            co2_factor =      1. / dtlsm
            crb_factor =      1. / dtlsm
            ent_factor =      1. / dtlsm
            h2o_factor =      1. / dtlsm

            !----- Fix the units so the terms are expressed as fluxes. --------------------!
            co2curr_residual       = co2curr_residual       * co2_factor
            co2budget_deltastorage = co2budget_deltastorage * co2_factor
            co2curr_nep            = co2curr_nep            * co2_factor
            co2curr_denseffect     = co2curr_denseffect     * co2_factor
            co2curr_zcaneffect     = co2curr_zcaneffect     * co2_factor
            co2curr_loss2atm       = co2curr_loss2atm       * co2_factor
            ccurr_residual         = ccurr_residual         * crb_factor
            cbudget_deltastorage   = cbudget_deltastorage   * crb_factor
            ccurr_denseffect       = ccurr_denseffect       * crb_factor
            ccurr_loss2yield       = ccurr_loss2yield       * crb_factor
            ccurr_zcaneffect       = ccurr_zcaneffect       * crb_factor
            ccurr_seedrain         = ccurr_seedrain         * crb_factor
            ccurr_loss2atm         = ccurr_loss2atm         * crb_factor
            ecurr_residual         = ecurr_residual         * ent_factor
            ebudget_deltastorage   = ebudget_deltastorage   * ent_factor
            ecurr_precipgain       = ecurr_precipgain       * ent_factor
            ecurr_netrad           = ecurr_netrad           * ent_factor
            ecurr_denseffect       = ecurr_denseffect       * ent_factor
            ecurr_prsseffect       = ecurr_prsseffect       * ent_factor
            ecurr_hcapeffect       = ecurr_hcapeffect       * ent_factor
            ecurr_wcapeffect       = ecurr_wcapeffect       * ent_factor
            ecurr_zcaneffect       = ecurr_zcaneffect       * ent_factor
            ecurr_loss2atm         = ecurr_loss2atm         * ent_factor
            ecurr_loss2drainage    = ecurr_loss2drainage    * ent_factor
            ecurr_loss2runoff      = ecurr_loss2runoff      * ent_factor
            wcurr_residual         = wcurr_residual         * h2o_factor
            wbudget_deltastorage   = wbudget_deltastorage   * h2o_factor
            wcurr_precipgain       = wcurr_precipgain       * h2o_factor
            wcurr_denseffect       = wcurr_denseffect       * h2o_factor
            wcurr_wcapeffect       = wcurr_wcapeffect       * h2o_factor
            wcurr_zcaneffect       = wcurr_zcaneffect       * h2o_factor
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
               , co2curr_zcaneffect     , co2curr_loss2atm       , cbudget_finalstorage    &
               , cbudget_committed      , ccurr_residual         , cbudget_deltastorage    &
               , ccurr_denseffect       , ccurr_zcaneffect       , ccurr_seedrain          &
               , ccurr_loss2yield       , ccurr_loss2atm         , ebudget_finalstorage    &
               , ecurr_residual         , ebudget_deltastorage   , ecurr_precipgain        &
               , ecurr_netrad           , ecurr_denseffect       , ecurr_prsseffect        &
               , ecurr_hcapeffect       , ecurr_wcapeffect       , ecurr_zcaneffect        &
               , ecurr_loss2atm         , ecurr_loss2drainage    , ecurr_loss2runoff       &
               , wbudget_finalstorage   , wcurr_residual         , wbudget_deltastorage    &
               , wcurr_precipgain       , wcurr_denseffect       , wcurr_wcapeffect        &
               , wcurr_zcaneffect       , wcurr_loss2atm         , wcurr_loss2drainage     &
               , wcurr_loss2runoff
            close(unit=86,status='keep')
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Stop the run in case there is any leak of CO2, enthalpy, or water.         !
         !---------------------------------------------------------------------------------!
         if (.not. budget_fine) then
            call fatal_error('Budget check has failed, see message above!'                 &
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
      !----- 3. Add the water vapour in the canopy air space. -----------------------------!
      compute_water_storage = compute_water_storage                                        &
                            + csite%can_shv(ipa) * csite%can_depth(ipa)                    &
                            * csite%can_rhos(ipa)
      !----- 4. Add the water on the leaf and wood surfaces. ------------------------------!
      do ico = 1,cpatch%ncohorts
         compute_water_storage = compute_water_storage + cpatch%leaf_water    (ico)        &
                                                       + cpatch%leaf_water_im2(ico)        &
                                                       + cpatch%wood_water    (ico)        &
                                                       + cpatch%wood_water_im2(ico)
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

      !----- 1. Add the ground components -------------------------------------------------!
      compute_netrad = csite%rshort_g(ipa) + csite%rlong_g(ipa) + csite%rlong_s(ipa)
      !----- 2. Add the shortwave radiation that reaches each snow/water layer ------------!
      do k = 1, csite%nlev_sfcwater(ipa)
         compute_netrad = compute_netrad + csite%rshort_s(k,ipa)
      end do
      !----- 3. Add the radiation components that is absorbed by leaves and branches. -----!
      do ico = 1,cpatch%ncohorts
         compute_netrad = compute_netrad + cpatch%rshort_l(ico) + cpatch%rlong_l(ico)      &
                                         + cpatch%rshort_w(ico) + cpatch%rlong_w(ico)
      end do
      return
   end function compute_netrad
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the total enthalpy stored in the thermodynamic pools that   !
   ! ED2 keeps track (soil, temporary surface water, vegetation, canopy air space).        !
   ! The result is given in J/m2.                                                          !
   !---------------------------------------------------------------------------------------!
   real function compute_enthalpy_storage(csite, lsl, ipa)
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
      ! 3. Find canopy air specific enthalpy then compute total enthalpy storage.          !
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
      compute_enthalpy_storage = soil_storage + sfcwater_storage + cas_storage + veg_storage
      !------------------------------------------------------------------------------------!

      return
   end function compute_enthalpy_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the carbon flux terms.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine sum_plant_cfluxes(csite,ipa, gpp, leaf_resp,root_resp,storage_resp           &
                               ,growth_resp)
      use ed_state_vars        , only : sitetype        & ! structure
                                      , patchtype       ! ! structure
      use consts_coms          , only : kgCday_2_umols  ! ! intent(in)
      use ed_max_dims          , only : n_dbh           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      integer               , intent(in)  :: ipa
      real                  , intent(out) :: gpp
      real                  , intent(out) :: leaf_resp
      real                  , intent(out) :: root_resp
      real                  , intent(out) :: storage_resp
      real                  , intent(out) :: growth_resp
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)       , pointer     :: cpatch
      integer                             :: ico
      !------------------------------------------------------------------------------------!


      !----- Initialize some variables. ---------------------------------------------------!
      gpp                = 0.0
      leaf_resp          = 0.0
      root_resp          = 0.0
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over cohorts to obtain GPP and metabolic respiration.                     !
      !------------------------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         !----- Add GPP and leaf respiration only for those cohorts with enough leaves. ---!
         if (cpatch%leaf_resolvable(ico)) then
            gpp       = gpp       + cpatch%gpp(ico)
            leaf_resp = leaf_resp + cpatch%leaf_respiration(ico)

         end if
         !----- Root respiration happens even when the LAI is tiny ------------------------!
         root_resp = root_resp + cpatch%root_respiration(ico)
      end do
      !------------------------------------------------------------------------------------!


      !----- Convert storage and growth respiration to umol/m2/s. -------------------------!
      storage_resp = csite%commit_storage_resp(ipa) * kgCday_2_umols
      growth_resp  = csite%commit_growth_resp (ipa) * kgCday_2_umols
      !------------------------------------------------------------------------------------!

      return
   end subroutine sum_plant_cfluxes
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the co2 stored in the canopy air space from ppm to umol/m2. !
   !---------------------------------------------------------------------------------------!
   real function compute_co2_storage(csite,ipa)
      use ed_state_vars, only : sitetype    ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype), target      :: csite
      integer       , intent(in)  :: ipa
      !----- Local variables. -------------------------------------------------------------!
      !------------------------------------------------------------------------------------!

      !----- Storage in umol/m2. ----------------------------------------------------------!
      compute_co2_storage = csite%can_dmol(ipa) * csite%can_depth(ipa) * csite%can_co2(ipa)
      !------------------------------------------------------------------------------------!

      return
   end function compute_co2_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the total carbon stored in the ED-2 system (kgC/m2).        !
   ! This includes carbon stocks from necromass, vegetation, seed bank, canopy air space   !
   ! and the committed changes in C stocks to be updated at the daily time step.           !
   !---------------------------------------------------------------------------------------!
   real function compute_carbon_storage(csite,ipa)
      use ed_state_vars  , only : sitetype              & ! structure
                                , patchtype             ! ! structure
      use ed_max_dims    , only : n_pft                 ! ! intent(in)
      use consts_coms    , only : mmdryi                & ! intent(in)
                                , umol_2_kgC            ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target      :: csite
      integer        , intent(in)  :: ipa
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer     :: cpatch
      integer                      :: ico
      real                         :: necro_storage
      real                         :: repro_storage
      real                         :: cas_storage
      real                         :: veg_storage
      !------------------------------------------------------------------------------------!


      cpatch => csite%patch(ipa)

      !----- 1. Find the total carbon stored in the necromass pools. ----------------------!
      necro_storage = csite%fast_grnd_C      (ipa) + csite%fast_soil_C      (ipa)          &
                    + csite%structural_grnd_C(ipa) + csite%structural_soil_C(ipa)          &
                    + csite%microbial_soil_C (ipa)                                         &
                    + csite%slow_soil_C      (ipa)                                         &
                    + csite%passive_soil_C   (ipa)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   2. Find the total carbon stored in the seed bank.                                !
      !------------------------------------------------------------------------------------!
      repro_storage = sum(csite%repro(:,ipa))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 3. Find carbon storage in the canopy air space based on CO2 mixing ratio.          !
      !------------------------------------------------------------------------------------!
      cas_storage  = compute_co2_storage(csite,ipa) * umol_2_kgC
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 4. Compute the internal energy stored in the plants.                               !
      !------------------------------------------------------------------------------------!
      veg_storage = 0.0
      do ico = 1,cpatch%ncohorts
         veg_storage = veg_storage                                                         &
                     + cpatch%nplant(ico) * ( cpatch%balive(ico) + cpatch%bdeada  (ico)    &
                                            + cpatch%bdeadb(ico) + cpatch%bstorage(ico) )
      end do
      !------------------------------------------------------------------------------------!


      !----- 5. Integrate the total carbon in ED. -----------------------------------------!
      compute_carbon_storage = necro_storage + repro_storage + cas_storage + veg_storage   &
                             + csite%cbudget_committed(ipa)
      !------------------------------------------------------------------------------------!

      return
   end function compute_carbon_storage
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the change of any property in the canopy air space due to  !
   ! changes in density associated with changes in pressure, temperature, and humidity,    !
   ! in order to satisfy the ideal gas law.                                                !
   !---------------------------------------------------------------------------------------!
   real function ddens_dt_effect(depth,old_dens,mid_dens,now_dens,old_prop,now_prop)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real, intent(in) :: depth      ! Depth (volume per unit area)
      real, intent(in) :: old_dens   ! Density before met update
      real, intent(in) :: mid_dens   ! Density after met update
      real, intent(in) :: now_dens   ! Current density
      real, intent(in) :: old_prop   ! Property before met update
      real, intent(in) :: now_prop   ! Property before met update
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     We use the new value of the property because density is updated the property   !
      ! has been updated.                                                                  !
      !------------------------------------------------------------------------------------!
      ddens_dt_effect = depth * ( (mid_dens - old_dens) * old_prop                         &
                                + (now_dens - mid_dens) * now_prop )
      !------------------------------------------------------------------------------------!

      return
   end function ddens_dt_effect
   !=======================================================================================!
   !=======================================================================================!
end module budget_utils
!==========================================================================================!
!==========================================================================================!
