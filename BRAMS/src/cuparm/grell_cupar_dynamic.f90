!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the dynamic control, including the effect of different    !
! cloud sizes interacting, in the case of Grell and Arakawa-Schubert closures.             !
! Whenever applicable, the follow convention is used:                                      !
!       [+] "?_cup"  : These are variables defined at the Grell level (staggered).         !
!       [+] "?d_cld" : These are cloud variables associated with downdrafts.               !
!       [+] "?u_cld" : These are cloud variables associated with updrafts.                 !
!       [+] "x_?"    : Variables modified by an arbitrary mass flux, for ensemble          !
!                      statistics and interaction between clouds.                          !
!       [+] "?0?"    : Whenever a variable like the ones above contains a 0, it means that !
!                      these variables are based on current values. If this is a BRAMS     !
!                      variable, it means the value before applying the tendency,          !
!                      otherwise it is a derived value based on the current large scale    !
!                      These are used only if Grell (1993) or Kain-Fritsch (1990)          !
!                      parameterizations are used.                                         !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine grell_cupar_dynamic(cldd,clds,nclouds,dtime,maxens_cap,maxens_eff,maxens_lsf    &
                              ,maxens_dyn,mgmzp,closure_type,comp_modif_thermo,prec_cld    &
                              ,cld2prec,mynum,i,j)
   use mem_ensemble     , only : ensemble_vars & ! structure
                               , ensemble_e    ! ! intent(inout)

   use mem_scratch_grell, only : &
       co2                & ! intent(in)  - CO2 mixing ratio with forcing.        [    ppm]
      ,co2sur             & ! intent(in)  - surface CO2 mixing ratio              [    ppm]
      ,dzd_cld            & ! intent(in)  - Delta-height for dndraft calculations [      m]
      ,dzu_cld            & ! intent(in)  - Delta-height for updraft calculations [      m]
      ,exner              & ! intent(in)  - Forced Exner funtion                  [ J/kg/K]
      ,exnersur           & ! intent(in)  - Surface Exner function                [ J/kg/K]
      ,mconv              & ! intent(in)  - Integrated moisture convergence       [kg/m²/s]
      ,mkx                & ! intent(in)  - Number of vertical levels             [    ---]
      ,omeg               & ! intent(in)  - Vertical velocity in pressure, dp/dt: [   Pa/s]
      ,p                  & ! intent(in)  - Pressure with forcing                 [     Pa]
      ,psur               & ! intent(in)  - Surface pressure                      [     Pa]
      ,qtot               & ! intent(in)  - Total mixing ratio with forcing.      [  kg/kg]
      ,qtotsur            & ! intent(in)  - surface total mixing ratio            [  kg/kg]
      ,qvap               & ! intent(in)  - Water vapour mix. ratio with forcing  [  kg/kg]
      ,qvapsur            & ! intent(in)  - surface water vapour mixing ratio     [  kg/kg]
      ,qliq               & ! intent(in)  - Liquid water mix. ratio with forcing  [  kg/kg]
      ,qliqsur            & ! intent(in)  - surface liquid water mixing ratio     [  kg/kg]
      ,qice               & ! intent(in)  - Ice mixing ratio with forcing.        [  kg/kg]
      ,qicesur            & ! intent(in)  - surface ice mixing ratio              [  kg/kg]
      ,rho                & ! intent(in)  - Air density                           [  kg/m³]
      ,theiv              & ! intent(in)  - THETA-Eiv with forcing                [      K]
      ,theivsur           & ! intent(in)  - Surface ice-vapour equiv. pot. temp.  [      K]
      ,thil               & ! intent(in)  - THETA-il with forcing                 [      K]
      ,thilsur            & ! intent(in)  - Surface ice-liquid potential temp.    [      K]
      ,t                  & ! intent(in)  - Temperature with forcing.             [      K]
      ,tsur               & ! intent(in)  - Surface temperature                   [      K]
      ,tscal_kf           & ! intent(in)  - Time scale for Kain-Fritsch (1990)    [      s]
      ,z                  & ! intent(in)  - Height                                [      m]
      ,z_cup              & ! intent(in)  - Height at cloud levels                [      m]
      ,cdd                & ! intent(out) - Normalised downdraft detrainment rate [    ---]
      ,cdu                & ! intent(out) - Normalised updraft detrainment rate   [    ---]
      ,mentrd_rate        & ! intent(out) - Normalised downdraft entrainment rate [    ---]
      ,mentru_rate        & ! intent(out) - Normalised updraft entrainment rate   [    ---]
      ,dbyd               & ! intent(out) - Buoyancy associated with downdrafts   [   m/s²]
      ,dbyu               & ! intent(out) - Buoyancy associated with updrafts     [   m/s²]
      ,etad_cld           & ! intent(out) - normalised downdraft mass flux        [    ---]
      ,etau_cld           & ! intent(out) - normalised updraft mass flux          [    ---]
      ,rhod_cld           & ! intent(out) - Downdraft density                     [  kg/m³]
      ,rhou_cld           & ! intent(out) - Updraft density                       [  kg/m³]
      ,qliqd_cld          & ! intent(out) - Liquid water mixing ratio at dndraft  [  kg/kg]
      ,qliqu_cld          & ! intent(out) - Liquid water mixing ratio at updraft  [  kg/kg]
      ,qiced_cld          & ! intent(out) - Ice mixing ratio at downdraft         [  kg/kg]
      ,qiceu_cld          & ! intent(out) - Ice mixing ratio at updraft           [  kg/kg]
      ,x_aad              & ! intent(out) - Updraft work function                 [   J/kg]
      ,x_aau              & ! intent(out) - Updraft work function                 [   J/kg]
      ,x_co2              & ! intent(out) - CO2 mixing ratio                      [    ppm]
      ,x_co2_cup          & ! intent(out) - CO2 mixing ratio                      [    ppm]
      ,x_co2d_cld         & ! intent(out) - CO2 mixing ratio                      [    ppm]
      ,x_co2u_cld         & ! intent(out) - CO2 mixing ratio                      [    ppm]
      ,x_dbyd             & ! intent(out) - Buoyancy acceleration                 [   m/s²]
      ,x_dbyu             & ! intent(out) - Buoyancy acceleration                 [   m/s²]
      ,x_exner_cup        & ! intent(out) - Exner function                        [   J/kg]
      ,x_p_cup            & ! intent(out) - Pressure                              [     Pa]
      ,x_pwav             & ! intent(out) - Integrated condensation               [  kg/kg]
      ,x_pwev             & ! intent(out) - Integrated evaporation                [  kg/kg]
      ,x_pwd_cld          & ! intent(out) - Condensation                          [  kg/kg]
      ,x_pwu_cld          & ! intent(out) - Condensation                          [  kg/kg]
      ,x_qtot             & ! intent(out) - Total mixing ratio                    [  kg/kg]
      ,x_qtot_cup         & ! intent(out) - Total mixing ratio                    [  kg/kg]
      ,x_qtotd_cld        & ! intent(out) - Total water mixing ratio              [  kg/kg]
      ,x_qtotu_cld        & ! intent(out) - Total water mixing ratio              [  kg/kg]
      ,x_qvap             & ! intent(out) - Vapour mixing ratio                   [  kg/kg]
      ,x_qvap_cup         & ! intent(out) - Water vapour mixing ratio             [  kg/kg]
      ,x_qvapd_cld        & ! intent(out) - Vapour mixing ratio                   [  kg/kg]
      ,x_qvapu_cld        & ! intent(out) - Vapour mixing ratio                   [  kg/kg]
      ,x_qliq             & ! intent(out) - Liquid water mixing ratio             [  kg/kg]
      ,x_qliq_cup         & ! intent(out) - Liquid water                          [  kg/kg]
      ,x_qliqd_cld        & ! intent(out) - Liquid water mixing ratio             [  kg/kg]
      ,x_qliqu_cld        & ! intent(out) - Liquid water mixing ratio             [  kg/kg]
      ,x_qice             & ! intent(out) - Ice mixing ratio                      [  kg/kg]
      ,x_qice_cup         & ! intent(out) - Mixing ratio                          [  kg/kg]
      ,x_qiced_cld        & ! intent(out) - Ice mixing ratio                      [  kg/kg]
      ,x_qiceu_cld        & ! intent(out) - Ice mixing ratio                      [  kg/kg]
      ,x_qsat_cup         & ! intent(out) - Saturation mixing ratio               [  kg/kg]
      ,x_qsatd_cld        & ! intent(out) - Sat. mixing ratio                     [  kg/kg]
      ,x_qsatu_cld        & ! intent(out) - Sat. mixing ratio                     [  kg/kg]
      ,x_rho_cup          & ! intent(out) - Density                               [  kg/m³]
      ,x_rhod_cld         & ! intent(out) - Density                               [  kg/m³]
      ,x_rhou_cld         & ! intent(out) - Density                               [  kg/m³]
      ,x_t                & ! intent(out) - Temperature                           [      K]
      ,x_t_cup            & ! intent(out) - Temperature                           [      K]
      ,x_td_cld           & ! intent(out) - Temperature                           [      K]
      ,x_tu_cld           & ! intent(out) - Temperature                           [      K]
      ,x_theiv            & ! intent(out) - Ice-vapour equiv. pot. temp.          [      K]
      ,x_theiv_cup        & ! intent(out) - Ice-vapour equiv. pot. temp.          [      K]
      ,x_theivs_cup       & ! intent(out) - Sat. Ice-vapour equiv. pot. temp.     [      K]
      ,x_theivd_cld       & ! intent(out) - Ice-vapour equiv. pot. temp.          [      K]
      ,x_theivu_cld       & ! intent(out) - Ice-vapour equiv. pot. temp.          [      K]
      ,x_thil             & ! intent(out) - Ice-liquid potential temperature.     [      K]
      ,x_thil_cup         & ! intent(out) - Ice-liquid potential temperature      [      K]
      ,x_thild_cld        & ! intent(out) - Ice-liquid potential temperature      [      K]
      ,x_thilu_cld        & ! intent(out) - Ice-liquid potential temperature      [      K]
      ,zero_scratch_grell ! ! subroutine - Resets scratch variables to zero.
   use rconstants, only: toodry
   use therm_lib , only : thil2tqall
   implicit none
   !---------------------------------------------------------------------------------------!
   ! List of arguments                                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Some dimensions, must come first because some variables depend on them. ---------!
   integer , intent(in)  :: cldd         ! Deepest Grell cloud
   integer , intent(in)  :: clds         ! Shallowest Grell cloud
   integer , intent(in)  :: nclouds      ! # of clouds
   integer , intent(in)  :: maxens_cap   ! Ensemble size, static control
   integer , intent(in)  :: maxens_eff   ! Ensemble size, prec. efficiency
   integer , intent(in)  :: maxens_lsf   ! Ensemble size, large-scale forc.
   integer , intent(in)  :: maxens_dyn   ! Ensemble size, dynamic control
   integer , intent(in)  :: mgmzp        ! Vertical grid size
   integer , intent(in)  :: mynum        ! Node ID
   integer , intent(in)  :: i            ! zonal grid point
   integer , intent(in)  :: j            ! meridional grid point
   !----- Character, containing the closure type(s) for dynamic control -------------------!
   character(len=2), intent(in) :: closure_type! Short name to define the method.
   !----- Logical flags, to bypass uncessary steps ----------------------------------------!
   logical                    , intent(in) :: comp_modif_thermo ! Compute interactions.
   logical, dimension(nclouds), intent(in) :: prec_cld          ! Precipitating cloud
   !----- Miscellaneous variables ---------------------------------------------------------!
   real, intent(in)  :: cld2prec     ! Fraction of cloud water converted to prec. [    ---]
   real, intent(in)  :: dtime        ! Time step.                                 [      s]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Local variables.                                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Counters and level holders ------------------------------------------------------!
   integer                   :: icld          ! Counter for current cloud
   integer                   :: jcld          ! Counter for other cloud
   integer                   :: icap          ! Counter for cap_maxs member
   integer                   :: iedt          ! Counter for prec. eff. member
   integer                   :: imbp          ! Counter for LSF member
   integer                   :: idyn          ! Counter for dynamic control member
   integer                   :: k             ! Level counter
   integer                   :: mboffset      ! Offset for the perturbation of forcing;
   !---------------------------------------------------------------------------------------!
   !    Aliases for some levels/flags, they will be stored in the ensemble structure       !
   ! later.                                                                                !
   !---------------------------------------------------------------------------------------!
   !----- Scalars. ------------------------------------------------------------------------!
   integer :: klod      ! Level in which downdrafts originate
   integer :: kdet      ! Top of downdraft detrainemnt layer
   integer :: kstabi    ! cloud stable layer base
   integer :: kstabm    ! cloud stable layer top
   !----- Combined variables --------------------------------------------------------------!
   logical, dimension(nclouds)            :: comp_dn  ! Downdrafts were computed.
   integer, dimension(nclouds)            :: ierr     ! Cloud error flag        [      ---]
   integer, dimension(nclouds)            :: klou     ! Level of origin of updrafts
   integer, dimension(nclouds)            :: klfc     ! Level of free convection
   integer, dimension(nclouds)            :: klnb     ! Level of newtral buoyancy
   integer, dimension(nclouds)            :: ktop     ! Cloud top
   real   , dimension(nclouds,nclouds)    :: mfke     ! Mass flux kernel        [ J m2/kg2]
   real   , dimension(nclouds)            :: aatot    ! Forced cloud work fctn. [     J/kg]
   real   , dimension(nclouds)            :: aatot0   ! Prev. cloud work fctn.  [     J/kg]
   real   , dimension(nclouds)            :: pwav     ! Int. condensed water    [    kg/kg]
   real   , dimension(nclouds)            :: pwev     ! Int. evaporated water   [    kg/kg]
   real   , dimension(nclouds)            :: prev_dnmf! Prev. dnward mass flux  [  kg/m2/s]
   real   , dimension(nclouds,maxens_dyn) :: upmf_dyn ! Ref. upward mass flux   [  kg/m2/s]
   real   , dimension(nclouds,maxens_dyn) :: dnmf_dyn ! Ref. dnward mass flux   [  kg/m2/s]
   real   , dimension(nclouds,maxens_dyn) :: upmx_dyn ! Max. upward mass flux   [  kg/m2/s]
   real   , dimension(nclouds,maxens_dyn) :: dnmx_dyn ! Max. dnward mass flux   [  kg/m2/s]
   !----- Dummy variables, needed as arguments for subroutines but not really used. -------!
   integer                                :: x_ierr   ! Flag for convection error
   integer                                :: x_klnb   ! Level of neutral buoyancy
   integer                                :: x_ktop   ! Cloud top level
   logical                                :: x_comp_dn! Downdraft flag
   real                                   :: x_qsat
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !    Miscellaneous parameters                                                           !
   !---------------------------------------------------------------------------------------!
   real, dimension(maxens_lsf) :: mbprime ! Arbitrary m. flux to modify environm. [kg/m²/s]
   real, dimension(maxens_lsf) :: one_b   ! 1-b, Krishnamurti et al. (1983)       [    ---]
   real, dimension(nclouds)    :: edt     ! Alias for the downdraft/updraft ratio.
   !------ Printing aux. variables. -------------------------------------------------------!
   logical          , parameter :: printing=.false. ! Printing some debug stuff   [    T|F]
   integer                      :: uni
   character(len=13), parameter :: fmti='(a,1x,i13,1x)'
   character(len=16), parameter :: fmtf='(a,1x,es13.6,1x)'
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Contrary to the static control, now we must consider all clouds together...  This  !
   ! is going to happen in a bunch of nested loops, so we get all permutations.            !
   !---------------------------------------------------------------------------------------!
   icldloop1: do icld = cldd, clds

      stacloop1: do icap = maxens_cap,1,-1

         call zero_scratch_grell(1)

         !---------------------------------------------------------------------------------!
         ! 1. Copying the variables stored at the ensemble structure to temporary arrays.  !
         !---------------------------------------------------------------------------------!
         !----- ierr and klnb have the cloud index because of the dyn. control ensemble. --!
         ierr    (icld) = ensemble_e(icld)%ierr_cap     (icap)
         comp_dn (icld) = ensemble_e(icld)%comp_down_cap(icap)
         klou    (icld) = ensemble_e(icld)%klou_cap     (icap)
         klfc    (icld) = ensemble_e(icld)%klfc_cap     (icap)
         klnb    (icld) = ensemble_e(icld)%klnb_cap     (icap)
         ktop    (icld) = ensemble_e(icld)%ktop_cap     (icap)
         klod           = ensemble_e(icld)%klod_cap     (icap)
         kdet           = ensemble_e(icld)%kdet_cap     (icap)
         kstabi         = ensemble_e(icld)%kstabi_cap   (icap)
         kstabm         = ensemble_e(icld)%kstabm_cap   (icap)

         do k=1,mkx
            cdd(k)         = ensemble_e(icld)%cdd_cap(k,icap)
            cdu(k)         = ensemble_e(icld)%cdu_cap(k,icap)
            mentrd_rate(k) = ensemble_e(icld)%mentrd_rate_cap(k,icap)
            mentru_rate(k) = ensemble_e(icld)%mentru_rate_cap(k,icap)
            dbyd(k)        = ensemble_e(icld)%dbyd_cap(k,icap)
            dbyu(k)        = ensemble_e(icld)%dbyu_cap(k,icap)
            etad_cld(k)    = ensemble_e(icld)%etad_cld_cap(k,icap)
            etau_cld(k)    = ensemble_e(icld)%etau_cld_cap(k,icap)
            rhod_cld(k)    = ensemble_e(icld)%rhod_cld_cap(k,icap)
            rhou_cld(k)    = ensemble_e(icld)%rhou_cld_cap(k,icap)
            qliqd_cld(k)   = ensemble_e(icld)%qliqd_cld_cap(k,icap)
            qliqu_cld(k)   = ensemble_e(icld)%qliqu_cld_cap(k,icap)
            qiced_cld(k)   = ensemble_e(icld)%qiced_cld_cap(k,icap)
            qiceu_cld(k)   = ensemble_e(icld)%qiceu_cld_cap(k,icap)
         end do

         effloop1: do iedt=1,maxens_eff
            !------------------------------------------------------------------------------!
            ! 2. Copying the current epsilon to a shortcut. This is an output variable,    !
            !    however, it will be used as scratch here. At the feedback time the value  !
            !    will be overwritten by the actual output value.                           !
            !------------------------------------------------------------------------------!
            edt(icld) = ensemble_e(icld)%edt_eff(iedt,icap)

            mbprimeloop1: do imbp=maxens_lsf,1,-1

               !---------------------------------------------------------------------------!
               ! 3. Initialise some ensemble-related scratch variables. mbprime is an      !
               !    arbitrary extra mass flux that is applied to the forcing, to get the   !
               !    upward mass flux for all cloud work related parametrisations. This is  !
               !    also a variable that we perturb to get different members for the       !
               !    ensemble.                                                              !
               !---------------------------------------------------------------------------!
               mboffset = (-1)**imbp * imbp /2 !----- This gives 0, 1, -1, 2, -2, ... -----!
               mbprime(imbp) = (4.+real(mboffset))*1.e-3

               !---------------------------------------------------------------------------!
               ! 4. Defining the Kuo (1974) moistening factor b. Following Krishnamurti    !
               !     et al. (1983), this value is around 0.3, with significant oscilla-    !
               !     tions between 0.0 and 0.8. Here we will play it closer to the aver-   !
               !     age, oscillating it between 0.2 and 0.4, equivalent to oscillate      !
               !     (1-b) between 0.6 and 0.8.                                            !
               !---------------------------------------------------------------------------!
               one_b(imbp) = 0.7+0.2*real(mboffset)

               !---------------------------------------------------------------------------!
               ! 5. Applying the change of arbitrary mass flux associated with a cloud of  !
               !    type "jcld" to some thermodynamic variables, and then we compute the   !
               !    cloud work function.                                                   !
               !---------------------------------------------------------------------------!
               jcldloop1: do jcld=cldd,clds
                  call zero_scratch_grell(2)
                  !------ We compute x_aatot only if the j cloud exists... ----------------!
                  if (ensemble_e(jcld)%ierr_cap(icap) == 0) then
                     !---------------------------------------------------------------------! 
                     !     Assign the  original values, then change the conditions from    !
                     ! the surface to the cloud top.                                       ! 
                     !---------------------------------------------------------------------! 
                     do k=1,ensemble_e(jcld)%ktop_cap(icap)
                        x_theiv(k) = theiv(k)                                              &
                                   + mbprime(imbp) * dtime                                 &
                                   * ensemble_e(jcld)%dellatheiv_eff(k,iedt,icap)
                        x_qtot(k)  = max(toodry,qtot(k)                                    &
                                   + mbprime(imbp) * dtime                                 &
                                   * ensemble_e(jcld)%dellaqtot_eff(k,iedt,icap))
                        x_thil(k)  = thil(k)                                               &
                                   + mbprime(imbp) * dtime                                 &
                                   * ensemble_e(jcld)%dellathil_eff(k,iedt,icap)
                        x_co2(k)   = co2(k)                                                &
                                   + mbprime(imbp) * dtime                                 &
                                   * ensemble_e(jcld)%dellaco2_eff(k,iedt,icap)
                     end do
                     do k = ensemble_e(jcld)%ktop_cap(icap)+1,mkx
                        x_theiv(k) = theiv(k)
                        x_qtot(k)  = qtot(k)
                        x_thil(k)  = thil(k)
                        x_co2(k)   = co2(k)
                     end do

                     !---------------------------------------------------------------------!
                     ! 5a. Initialise some variables.                                      !
                     !---------------------------------------------------------------------!
                     x_ierr    = 0
                     x_comp_dn = comp_dn(icld)
                     !----- Cloud work ----------------------------------------------------!
                     x_aad     = 0.
                     x_aau     = 0.
 
                     modif_comp_if: if (comp_modif_thermo .and.                            &
                                        ensemble_e(icld)%ierr_cap(icap) == 0) then
                        !------------------------------------------------------------------!
                        ! 5b. Compute the modified structure, finding the consistent set   !
                        !     and then interpolating them to the cloud levels.             !
                        !------------------------------------------------------------------!
                        do k=1,mkx
                           x_t(k)    = t  (k)
                           call thil2tqall(x_thil(k),exner(k),p(k),x_qtot(k),x_qliq(k)     &
                                          ,x_qice(k),x_t(k),x_qvap(k),x_qsat)
                        end do
                        call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner,x_thil,x_t,x_qtot   &
                                                ,x_qliq,x_qice,x_co2,exnersur,thilsur,tsur &
                                                ,qtotsur,qliqsur,qicesur,co2sur            &
                                                ,x_exner_cup,x_p_cup,x_t_cup,x_thil_cup    &
                                                ,x_qtot_cup,x_qvap_cup,x_qliq_cup          &
                                                ,x_qice_cup,x_qsat_cup,x_co2_cup,x_rho_cup &
                                                ,x_theiv_cup,x_theivs_cup)

                        !------------------------------------------------------------------!
                        ! 5c. Finding the updraft thermodynamics between the updraft       !
                        !     origin and the level of free convection.                     !
                        !------------------------------------------------------------------!
                        call grell_buoy_below_lfc(mkx,mgmzp,klou(icld),klfc(icld)          &
                                                 ,x_exner_cup,x_p_cup,x_theiv_cup          &
                                                 ,x_thil_cup,x_t_cup,x_qtot_cup,x_qvap_cup &
                                                 ,x_qliq_cup,x_qice_cup,x_qsat_cup         &
                                                 ,x_co2_cup,x_rho_cup,x_theivu_cld         &
                                                 ,x_thilu_cld,x_tu_cld,x_qtotu_cld         &
                                                 ,x_qvapu_cld,x_qliqu_cld,x_qiceu_cld      &
                                                 ,x_qsatu_cld,x_co2u_cld,x_rhou_cld,x_dbyu)

                        !------------------------------------------------------------------!
                        ! 5d. Compute the updraft profiles ice-vapour equivalent potential !
                        !     temperature.                                                 !
                        !------------------------------------------------------------------!
                        call grell_theiv_updraft(mkx,mgmzp,klou(icld),klfc(icld),cdu       &
                                                ,mentru_rate,x_theiv,x_theiv_cup,dzu_cld   &
                                                ,x_theivu_cld)

                        !------------------------------------------------------------------!
                        ! 5e. Getting the updraft moisture profile                         !
                        !------------------------------------------------------------------!
                        call grell_most_thermo_updraft(prec_cld(icld),.false.,mkx,mgmzp    &
                                                      ,klfc(icld),ktop(icld),cld2prec,cdu  &
                                                      ,mentru_rate,x_qtot,x_co2,x_p_cup    &
                                                      ,x_exner_cup,x_theiv_cup,x_thil_cup  &
                                                      ,x_t_cup,x_qtot_cup,x_qvap_cup       &
                                                      ,x_qliq_cup,x_qice_cup,x_qsat_cup    &
                                                      ,x_co2_cup,x_rho_cup,x_theivu_cld    &
                                                      ,etau_cld,dzu_cld,x_thilu_cld        &
                                                      ,x_tu_cld,x_qtotu_cld,x_qvapu_cld    &
                                                      ,x_qliqu_cld,x_qiceu_cld,x_qsatu_cld &
                                                      ,x_co2u_cld,x_rhou_cld,x_dbyu        &
                                                      ,x_pwu_cld,x_pwav,x_klnb,x_ktop      &
                                                      ,x_ierr)

                        !------------------------------------------------------------------!
                        ! 5f. Recalculating the updraft cloud work                         !
                        !------------------------------------------------------------------!
                        call grell_cldwork_updraft(mkx,mgmzp,klou(icld),ktop(icld),x_dbyu  &
                                                  ,dzu_cld,etau_cld,x_aau)

                        modif_down: if (comp_dn(icld)) then

                           !---------------------------------------------------------------!
                           ! 5g. Finding moist static energy                               !
                           !---------------------------------------------------------------!
                           call grell_theiv_downdraft(mkx,mgmzp,klod,cdd,mentrd_rate       &
                                                     ,x_theiv,x_theiv_cup,x_theivs_cup     &
                                                     ,dzd_cld,x_theivd_cld)

                           !---------------------------------------------------------------!
                           ! 5h. Moisture properties                                       !
                           !---------------------------------------------------------------!
                           call grell_most_thermo_downdraft(mkx,mgmzp,klod,x_qtot,x_co2    &
                                                           ,mentrd_rate,cdd,x_p_cup        &
                                                           ,x_exner_cup,x_thil_cup,x_t_cup &
                                                           ,x_qtot_cup,x_qvap_cup          &
                                                           ,x_qliq_cup,x_qice_cup          &
                                                           ,x_qsat_cup,x_co2_cup,x_rho_cup &
                                                           ,x_pwav,x_theivd_cld            &
                                                           ,etad_cld,dzd_cld               &
                                                           ,x_thild_cld,x_td_cld           &
                                                           ,x_qtotd_cld,x_qvapd_cld        &
                                                           ,x_qliqd_cld,x_qiced_cld        &
                                                           ,x_qsatd_cld,x_co2d_cld         &
                                                           ,x_rhod_cld,x_dbyd,x_pwd_cld    &
                                                           ,x_pwev,x_ierr)

                           !---------------------------------------------------------------!
                           ! 5i. Computing cloud work function associated with downdrafts. !
                           !---------------------------------------------------------------!
                           call grell_cldwork_downdraft(mkx,mgmzp,klod,x_dbyd,dzd_cld      &
                                                       ,etad_cld,x_aad)
                        end if modif_down
                     end if modif_comp_if

                     !---------------------------------------------------------------------!
                     ! 5j. If the perturbed field produced a cloud, we compute the total   !
                     !     cloud work function for this particular test.  If not, or if    !
                     !     the user is using moisture convergence closure only, this will  !
                     !     be a sum of zeroes.                                             !
                     !---------------------------------------------------------------------!
                     ensemble_e(icld)%x_aatot(jcld,imbp,iedt,icap) = x_aau                 &
                                                                   + edt(icld) * x_aad
                  else
                     !---------------------------------------------------------------------!
                     ! 5k. No cloud here, we assign the cloud work function to be exactly  !
                     !     the cloud work function without the perturbation.               !
                     !     The reference mass flux of cloud jcld is zero.                  !
                     !---------------------------------------------------------------------!
                     ensemble_e(icld)%x_aatot(jcld,imbp,iedt,icap) =                       &
                                                      ensemble_e(icld)%aatot_eff(iedt,icap)
                  end if
               end do jcldloop1
            end do mbprimeloop1
         end do effloop1
      end do stacloop1
   end do icldloop1

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   if (printing) then
      uni = 60 + mynum
      write(unit=uni,fmt='(a)') '---------------------------------------------------------'
      write(unit=uni,fmt=fmti ) ' I        =', i
      write(unit=uni,fmt=fmti ) ' J        =', j
      write(unit=uni,fmt='(a)') ' '
   end if
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!


   !---------------------------------------------------------------------------------------!
   ! 6. Running the dynamic control ensemble, calling it for each combination of cloud     !
   !    efficiency, static control, and large-scale forcing.  Inside this subroutine we    !
   !    will find all members of the dynamic control for all clouds.  The outcome will be  !
   !    a vector of reference updraft mass flux for that member.                           !
   !     mass flux will then have maxens_eff × maxens_lsf × maxens_dyn different values.   !
   !---------------------------------------------------------------------------------------!
   stacloop2: do icap=1,maxens_cap
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      if (printing) then
         write(unit=uni,fmt=fmti ) '  ICAP     =', icap
         write(unit=uni,fmt='(a)') ' '
      end if
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      effloop2: do iedt=1,maxens_eff

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         if (printing) then
            write(unit=uni,fmt=fmti ) '   IEDT     =', iedt
            write(unit=uni,fmt='(a)') ' '
         end if
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         mbprimeloop2: do imbp=1,maxens_lsf
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            if (printing) then
               write(unit=uni,fmt=fmti ) '    IMBP     =', imbp
               write(unit=uni,fmt=fmtf ) '    MBPRIME  =', mbprime(1)
               write(unit=uni,fmt='(a)') ' '
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

            !------------------------------------------------------------------------------!
            ! 6a. Copying some variables to scratch arrays.                                !
            !------------------------------------------------------------------------------!
            mfke   = 0.
            ierr   = 1
            klnb   = 1
            ktop   = 1
            aatot0 = 0.
            aatot  = 0.
            do icld=cldd,clds
               edt(icld)       = ensemble_e(icld)%edt_eff(iedt,icap)
               aatot0(icld)    = ensemble_e(icld)%aatot0_eff(iedt,icap)
               aatot(icld)     = ensemble_e(icld)%aatot_eff(iedt,icap)
               klou(icld)      = ensemble_e(icld)%klou_cap(icap)
               klfc(icld)      = ensemble_e(icld)%klfc_cap(icap)
               klnb(icld)      = ensemble_e(icld)%klnb_cap(icap)
               ktop(icld)      = ensemble_e(icld)%ktop_cap(icap)
               ierr(icld)      = ensemble_e(icld)%ierr_cap(icap)
               pwav(icld)      = ensemble_e(icld)%pwav_cap(icap)
               pwev(icld)      = ensemble_e(icld)%pwev_cap(icap)
               prev_dnmf(icld) = ensemble_e(icld)%prev_dnmf(1)

               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               if (printing) then
                  write(unit=uni,fmt=fmti ) '     ICLD     =', icld
                  write(unit=uni,fmt=fmtf ) '     EDT      =', edt(icld)
                  write(unit=uni,fmt=fmtf ) '     TSCAL_KF =', tscal_kf
                  write(unit=uni,fmt=fmtf ) '     DTIME    =', dtime
                  if (klou(icld) > 0) then
                     write(unit=uni,fmt=fmtf ) '     KLOU     =', z(klou(icld))
                  end if
                  if (klfc(icld) > 0) then
                     write(unit=uni,fmt=fmtf ) '     KLFC     =', z(klfc(icld))
                  end if
                  if (klnb(icld) > 0) then
                     write(unit=uni,fmt=fmtf ) '     KLNB     =', z(klnb(icld))
                  end if
                  if (ktop(icld) > 0) then
                     write(unit=uni,fmt=fmtf ) '     KTOP     =', z(ktop(icld))
                  end if
                  write(unit=uni,fmt=fmti ) '     IERR     =', ierr(icld)
                  write(unit=uni,fmt=fmtf ) '     AATOT0   ='                              &
                                            , ensemble_e(icld)%aatot0_eff(iedt,icap)
                  write(unit=uni,fmt=fmtf ) '     AATOT    ='                              &
                                            , ensemble_e(icld)%aatot_eff(iedt,icap)
                  write(unit=uni,fmt='(a)') ' '
               end if
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!


               !---------------------------------------------------------------------------!
               ! 6b. MFKE is the kernel matrix for this member.                            !
               !---------------------------------------------------------------------------!
               do jcld=cldd,clds
                  mfke(icld,jcld) = ( ensemble_e(icld)%x_aatot(jcld,imbp,iedt,icap)        &
                                    - aatot(icld)) / (mbprime(imbp) * dtime) 

                  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
                  !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
                  if (printing) then
                     write(unit=uni,fmt=fmti ) '      JCLD     =', jcld
                     write(unit=uni,fmt=fmtf ) '      X_AATOT  =',                         &
                                              ensemble_e(icld)%x_aatot(jcld,imbp,iedt,icap)
                     write(unit=uni,fmt=fmtf ) '      MFKE     =', mfke(icld,jcld)
                     write(unit=uni,fmt='(a)') ' '
                  end if
                  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
                  !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               end do
            end do


            !------------------------------------------------------------------------------!
            ! 6c. Solving the reference mass flux for all possible dynamic controls for    !
            !     this member.                                                             !
            !------------------------------------------------------------------------------!
            call grell_dyncontrol_ensemble(nclouds,mgmzp,maxens_dyn,cldd,clds,dtime        &
                                          ,closure_type,comp_dn,tscal_kf,mconv             &
                                          ,omeg,x_p_cup,edt,mbprime(imbp),one_b(imbp)      &
                                          ,aatot0,aatot,ierr,klou,klfc,ktop,mfke,pwav,pwev &
                                          ,prev_dnmf,dnmf_dyn,upmf_dyn,dnmx_dyn,upmx_dyn)
            
            !------------------------------------------------------------------------------!
            ! 6d. Copying back the error flag and mass fluxes to the ensemble structures.  !
            !------------------------------------------------------------------------------!
            do icld=cldd,clds
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               if (printing) then
                  write(unit=uni,fmt=fmti ) '     ICLD     =', icld
                  write(unit=uni,fmt='(a)') ' '
               end if
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

               do idyn=1,maxens_dyn
                  ensemble_e(icld)%dnmf_ens(idyn,imbp,iedt,icap) = dnmf_dyn(icld,idyn)
                  ensemble_e(icld)%upmf_ens(idyn,imbp,iedt,icap) = upmf_dyn(icld,idyn)
                  ensemble_e(icld)%dnmx_ens(idyn,imbp,iedt,icap) = dnmx_dyn(icld,idyn)
                  ensemble_e(icld)%upmx_ens(idyn,imbp,iedt,icap) = upmx_dyn(icld,idyn)
                  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
                  !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
                  if (printing) then
                     write(unit=uni,fmt=fmti ) '      IDYN     =', idyn
                     write(unit=uni,fmt=fmtf ) '      DNMF     =', dnmf_dyn(icld,idyn)
                     write(unit=uni,fmt=fmtf ) '      UPMF     =', upmf_dyn(icld,idyn)
                     write(unit=uni,fmt='(a)') ' '
                  end if
                  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
                  !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

                  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
                  !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
                  !write(unit=60+mynum,fmt='(2(a,1x,f10.4,1x))') '     dnmf   =',ensemble_e(icld)%dnmf_ens(1,imbp,iedt,icap),'upmf   =',ensemble_e(icld)%upmf_ens(1,imbp,iedt,icap)
                  !write(unit=60+mynum,fmt='(a)')                '-------------------------------------------------------------------------------------------------'
                  !write(unit=60+mynum,fmt='(a)'               ) ' '
                  !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
                  !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               end do
            end do
        end do mbprimeloop2
      end do effloop2
   end do stacloop2

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   if (printing) then
      write(unit=uni,fmt='(a)') '---------------------------------------------------------'
      write(unit=uni,fmt='(a)') ' '
   end if
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

   return
end subroutine grell_cupar_dynamic
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the mass flux associated with the dynamic control. This      !
! should be called inside the other loops, so the ensemble variables that should have      !
! the permutation of all ensemble dimensions are fully filled.                             !
!------------------------------------------------------------------------------------------!
subroutine grell_dyncontrol_ensemble(nclouds,mgmzp,maxens_dyn,cldd,clds,dtime,closure_type &
                                    ,comp_down,tscal_kf,mconv,omeg,p_cup,edt,mbprime,one_b &
                                    ,aatot0,aatot,ierr,klou,klfc,ktop,mfke,pwav,pwev       &
                                    ,prev_dnmf,dnmf_dyn,upmf_dyn,dnmx_dyn,upmx_dyn)
   use rconstants  , only : grav           ! ! Gravity acceleration.

   implicit none
   !----- Input array dimensions and boundaries. ------------------------------------------!
   integer                    , intent(in)  :: nclouds      ! # of levels
   integer                    , intent(in)  :: mgmzp        ! # of levels
   integer                    , intent(in)  :: maxens_dyn   ! # of dynamic control members
   integer                    , intent(in)  :: cldd         ! Deepest cloud
   integer                    , intent(in)  :: clds         ! Shallowest cloud
   !----- Input flags ---------------------------------------------------------------------!
   character(len=2)           , intent(in)  :: closure_type ! My dynamic control
   logical, dimension(nclouds), intent(in)  :: comp_down    ! I have downdrafts and rain.
   !----- Input constant properties -------------------------------------------------------!
   real                       , intent(in)  :: dtime        ! Current grid time step
   real                       , intent(in)  :: tscal_kf     ! Kain-Fritsch time scale
   real                       , intent(in)  :: mconv        ! Integ. moisture convergence
   !----- Input large scale variables -----------------------------------------------------!
   real, dimension(mgmzp)     , intent(in)  :: omeg         ! Vertical pressure velocity
   real, dimension(mgmzp)     , intent(in)  :: p_cup        ! Pressure with forcing
   !----- Input ensemble-dependent variables ----------------------------------------------!
   real, dimension(nclouds)   , intent(in)  :: edt          ! Epsilon, I1/I2
   real                       , intent(in)  :: mbprime      ! Arbitrary mass flux
   real                       , intent(in)  :: one_b        ! (1-b) from Krishnamurti '83.
   !----- Input cloud- and ensemble-dependent variables. ----------------------------------!
   real   , dimension(nclouds)        , intent(in) :: aatot0    ! Curr. cloud work function
   real   , dimension(nclouds)        , intent(in) :: aatot     ! Cloud work with forcing
   real   , dimension(nclouds,nclouds), intent(in) :: mfke      ! Mass flux kernel.
   real   , dimension(nclouds)        , intent(in) :: pwav      ! Int. condensed water
   real   , dimension(nclouds)        , intent(in) :: pwev      ! Int. evaporated water
   real   , dimension(nclouds)        , intent(in) :: prev_dnmf ! Int. evaporated water
   integer, dimension(nclouds)        , intent(in) :: ierr      ! Error flag
   integer, dimension(nclouds)        , intent(in) :: klou      ! Updraft origin
   integer, dimension(nclouds)        , intent(in) :: klfc      ! Level of free convection
   integer, dimension(nclouds)        , intent(in) :: ktop      ! Cloud top
   !----- Output variables. ---------------------------------------------------------------!
   real, dimension(nclouds,maxens_dyn), intent(inout) :: dnmf_dyn ! Ref. dndraft mass flux
   real, dimension(nclouds,maxens_dyn), intent(inout) :: upmf_dyn ! Ref. updraft mass flux
   real, dimension(nclouds,maxens_dyn), intent(inout) :: dnmx_dyn ! Max. dndraft mass flux
   real, dimension(nclouds,maxens_dyn), intent(inout) :: upmx_dyn ! Max. updraft mass flux
   !----- Local variables -----------------------------------------------------------------!
   integer                      :: icld      ! Cloud index
   integer                      :: idyn      ! Dynamic control counter
   integer                      :: ksmf      ! Level with the strongest mass flux.
   real                         :: divisor   ! Scratch, with the divisor
   !----- Local constants. ----------------------------------------------------------------!
   real, parameter   :: tinyden=1.e-25 ! Small number for division
   !---------------------------------------------------------------------------------------!
 
   !----- Initialise mass flux to be zero, so if anything goes wrong, they are set. -------!
   do idyn=1,maxens_dyn
      do icld=1,nclouds
         upmf_dyn(icld,idyn) = 0.
         dnmf_dyn(icld,idyn) = 0.
         upmx_dyn(icld,idyn) = 0.
         dnmx_dyn(icld,idyn) = 0.
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !    Solving the dynamic control according to the type of closure we want.              !
   !---------------------------------------------------------------------------------------!
   select case (closure_type)
   !---------------------------------------------------------------------------------------!
   ! 1. Standard Grell (1993), modified quasi-equilibrium buoyant energy.                  !
   !---------------------------------------------------------------------------------------!
   case ('gr')
      call grell_grell_solver(nclouds,cldd,clds,dtime,1.0,aatot0,aatot,mfke,ierr           &
                             ,upmf_dyn(1:nclouds,1),upmx_dyn(1:nclouds,1))

   !---------------------------------------------------------------------------------------!
   ! 2. Standard Arakawa and Schubert (1974), quasi-equilibrium buoyant energy.            !
   !---------------------------------------------------------------------------------------!
   case ('as')
      call grell_arakschu_solver(nclouds,cldd,clds,mgmzp,dtime,p_cup,1,1,ktop,aatot,mfke   &
                                ,ierr,upmf_dyn(1:nclouds,1),upmx_dyn(1:nclouds,1))

   !---------------------------------------------------------------------------------------!
   ! 3. Standard Kain and Fritsch (1990), instability removal.                             !
   !---------------------------------------------------------------------------------------!
   case ('kf')
      call grell_inre_solver(nclouds,cldd,clds,tscal_kf,1.0,aatot0,mfke,ierr               &
                            ,upmf_dyn(1:nclouds,1),upmx_dyn(1:nclouds,1))

   !---------------------------------------------------------------------------------------!
   ! 4. Standard Frank and Cohen (1987), low-level environment mass flux                   !
   !---------------------------------------------------------------------------------------!
   case ('lo')
      loloop: do icld = cldd,clds
         if (ierr(icld) /= 0) cycle loloop
         upmf_dyn(icld,1) = max(0.,-omeg(klou(icld))/grav - prev_dnmf(icld))
         upmx_dyn(icld,1) = upmf_dyn(icld,1)
      end do loloop

   !---------------------------------------------------------------------------------------!
   ! 5. Standard Krishnamurti et al. (1983), moisture convergence.                         !
   !---------------------------------------------------------------------------------------!
   case ('mc')
      mcloop: do icld = cldd, clds
         if (ierr(icld) /= 0) cycle mcloop
         divisor     = sign(max(tinyden,abs(pwav(icld)- edt(icld) * pwev(icld)))           &
                           ,pwav(icld)-edt(icld)*pwev(icld))
         upmf_dyn(icld,1) = max(0.,mconv * 1.4 * one_b / divisor)
         upmx_dyn(icld,1) = upmf_dyn(icld,1)
      end do mcloop

   !---------------------------------------------------------------------------------------!
   ! 6. Ensemble, based on Grell and Dévényi (2002)                                        !
   !    Here for each different dynamic control style, the first element is the standard,  !
   !    like in the previous cases, and the others are perturbations.                      !
   !---------------------------------------------------------------------------------------!
   case ('en','nc')

      !------------------------------------------------------------------------------------!
      ! 6a. Grell (1993), modified quasi-equilibrium buoyant energy.                       !
      !------------------------------------------------------------------------------------!
      call grell_grell_solver(nclouds,cldd,clds,dtime,1.0,aatot0,aatot,mfke,ierr           &
                             ,upmf_dyn(1:nclouds,1),upmx_dyn(1:nclouds,1))
      call grell_grell_solver(nclouds,cldd,clds,dtime,0.9,aatot0,aatot,mfke,ierr           &
                             ,upmf_dyn(1:nclouds,2),upmx_dyn(1:nclouds,2))
      call grell_grell_solver(nclouds,cldd,clds,dtime,1.1,aatot0,aatot,mfke,ierr           &
                             ,upmf_dyn(1:nclouds,3),upmx_dyn(1:nclouds,3))

      !------------------------------------------------------------------------------------!
      ! 6b. Arakawa and Schubert (1974), quasi-equilibrium buoyant energy.                 !
      !------------------------------------------------------------------------------------!
      !------ Computing the upward mass flux ----------------------------------------------!
      call grell_arakschu_solver(nclouds,cldd,clds,mgmzp,dtime,p_cup,1,1,ktop,aatot,mfke   &
                                ,ierr,upmf_dyn(1:nclouds,4),upmx_dyn(1:nclouds,4))
      call grell_arakschu_solver(nclouds,cldd,clds,mgmzp,dtime,p_cup,2,1,ktop,aatot,mfke   &
                                ,ierr,upmf_dyn(1:nclouds,5),upmx_dyn(1:nclouds,5))
      call grell_arakschu_solver(nclouds,cldd,clds,mgmzp,dtime,p_cup,1,2,ktop,aatot,mfke   &
                                ,ierr,upmf_dyn(1:nclouds,6),upmx_dyn(1:nclouds,6))
      call grell_arakschu_solver(nclouds,cldd,clds,mgmzp,dtime,p_cup,2,2,ktop,aatot,mfke   &
                                ,ierr,upmf_dyn(1:nclouds,7),upmx_dyn(1:nclouds,7))

      !------------------------------------------------------------------------------------!
      ! 6c. Kain and Fritsch (1990), instability removal.                                  !
      !------------------------------------------------------------------------------------!
      call grell_inre_solver(nclouds,cldd,clds,tscal_kf,1.0,aatot0,mfke,ierr               &
                            ,upmf_dyn(1:nclouds,8),upmx_dyn(1:nclouds,8))
      call grell_inre_solver(nclouds,cldd,clds,tscal_kf,0.9,aatot0,mfke,ierr               &
                            ,upmf_dyn(1:nclouds,9),upmx_dyn(1:nclouds,9))
      call grell_inre_solver(nclouds,cldd,clds,tscal_kf,1.1,aatot0,mfke,ierr               &
                            ,upmf_dyn(1:nclouds,10),upmx_dyn(1:nclouds,10))

      if (closure_type == 'en') then
         !---------------------------------------------------------------------------------!
         ! 6d. Frank and Cohen (1987), low-level environment mass flux.                    !
         !---------------------------------------------------------------------------------!
         enloloop: do icld=cldd,clds
            if (ierr(icld) /= 0) cycle enloloop
            upmf_dyn(icld,11) = max(0.,-omeg(klou(icld))/grav   - prev_dnmf(icld))
            upmf_dyn(icld,12) = max(0.,-omeg(klfc(icld))/grav - prev_dnmf(icld))
            !----- Picking up the strongest mass flux below the LFC (except klou) ---------!
            ksmf = minloc(omeg(1:(klfc(icld)-1)),dim=1                                     &
                         ,mask=omeg(1:(klfc(icld)-1)) /= omeg(klou(icld)))
            upmf_dyn(icld,13) = max(0.,-omeg(ksmf)/grav - prev_dnmf(icld))
            upmx_dyn(icld,11) = upmf_dyn(icld,11)
            upmx_dyn(icld,12) = upmf_dyn(icld,12)
            upmx_dyn(icld,13) = upmf_dyn(icld,13)
         end do enloloop

         
         !---------------------------------------------------------------------------------!
         ! 6e. Krishnamurti et al. (1983), moisture convergence. Here I am using some      !
         !     different values of (1+eta) (what Grell and Devenyi (2002) called f_emp).   !
         !     I'm oscillating it between (1+0.25) and (1+0.55), which is respectively     !
         !     roughly 0.9 and 1.1 times (1+0.4), their average value                      !
         !---------------------------------------------------------------------------------!
         enmcloop: do icld=cldd,clds
            if (ierr(icld) /= 0) cycle enmcloop
            divisor      = sign(max(tinyden,abs(pwav(icld)- edt(icld) * pwev(icld)))       &
                               ,pwav(icld)-edt(icld)*pwev(icld))
            upmf_dyn(icld,14) = max(0.,mconv * 1.4 * one_b / divisor)
            upmf_dyn(icld,15) = 1.1 * upmf_dyn(icld,14)
            upmf_dyn(icld,16) = 0.9 * upmf_dyn(icld,14)
            upmx_dyn(icld,14) = upmf_dyn(icld,14)
            upmx_dyn(icld,15) = upmf_dyn(icld,15)
            upmx_dyn(icld,16) = upmf_dyn(icld,16)
         end do enmcloop
      end if
   end select

   do idyn=1,maxens_dyn
      do icld=cldd,clds
         dnmf_dyn(icld,idyn) = edt(icld) * upmf_dyn(icld,idyn)
         dnmx_dyn(icld,idyn) = edt(icld) * upmx_dyn(icld,idyn)
      end do
   end do

   return
end subroutine grell_dyncontrol_ensemble
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will solve the mass fluxes of a collection of Grell clouds using     !
! the quasi-equilibrium assumption.  This is used by Grell (1993) dynamic control.         !
!------------------------------------------------------------------------------------------!
subroutine grell_grell_solver(nclouds,cldd,clds,dtime,fac,aatot0,aatot,mfke,ierr,upmf,upmx)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                            , intent(in)    :: cldd    ! Deepest Grell cloud
   integer                            , intent(in)    :: clds    ! Shallowest Grell cloud
   integer                            , intent(in)    :: nclouds ! # of clouds
   real                               , intent(in)    :: dtime   ! Current grid time step
   real                               , intent(in)    :: fac     ! For ensemble variance
   real   , dimension(nclouds)        , intent(in)    :: aatot0  ! Current cloud work fctn.
   real   , dimension(nclouds)        , intent(in)    :: aatot   ! Cloud work with forcing
   real   , dimension(nclouds,nclouds), intent(in)    :: mfke    ! Mass flux kernel.
   !----- Downdraft is inout because I may not compute it ---------------------------------!
   integer, dimension(nclouds)        , intent(in)    :: ierr    ! Error flag
   real   , dimension(nclouds)        , intent(inout) :: upmf    ! Ref. updraft mass flux
   real   , dimension(nclouds)        , intent(inout) :: upmx    ! Max. updraft mass flux
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: icld    ! Cloud index
   integer                              :: jcld    ! Cloud index too
   integer                              :: isol    ! Another cloud index
   integer                              :: jsol    ! Yet another cloud index
   real                                 :: dtimei  ! Scratch, inverse of delta-t
   integer                              :: nsolv   ! # of clouds we are still solving.
   logical                              :: is_sing ! Flag.  The kernel was singular. [T|F]
   !----- Auxilliary variables, containing flags and references. --------------------------!
   logical, dimension(nclouds)          :: okcld   ! This cloud is still ok to use.
   integer, dimension(nclouds)          :: cloud   ! Absolute cloud index.
   !----- Temporary, allocatable arrays. --------------------------------------------------!
   real   , dimension(:,:), allocatable :: kke     ! Subset of mass flux kernel
   real   , dimension(:)  , allocatable :: diagkke ! Diagonal of kke
   real   , dimension(:)  , allocatable :: mfo     ! Minus the large-scale forcing
   real   , dimension(:)  , allocatable :: mb      ! Mass flux
   integer, dimension(:)  , allocatable :: cldidx  ! Cloud index, to copy back to ensemble
   !---------------------------------------------------------------------------------------!


   !----- Initial settings. ---------------------------------------------------------------!
   dtimei = 1. / dtime
   do icld=1,nclouds
      cloud(icld) = icld
   end do

   !----- Discarding non-Grell clouds. ----------------------------------------------------!
   do icld=1,cldd-1
      okcld(icld) = .false.
   end do
   !---------------------------------------------------------------------------------------!
   !    Determine the clouds that may be solved here, and allocate a subset with these     !
   ! clouds only.  We may eliminate some other clouds during this process because they     !
   ! could lead to negative mass fluxes.                                                   !
   !---------------------------------------------------------------------------------------!
   do icld=cldd,clds
      okcld(icld) = ierr(icld) == 0
   end do
   !----- Discarding non-Grell clouds. ----------------------------------------------------!
   do icld=clds+1,nclouds
      okcld(icld) = .false.
   end do

   !----- Assign zero to the fluxes, just in case something prevents the calculation. -----!
   do icld=1,nclouds
      upmf(icld) = 0.
      upmx(icld) = 0.
   end do

   !----- Initialise the cloud arrays. ----------------------------------------------------!
   queq_loop: do
      nsolv = count(okcld)

      !------------------------------------------------------------------------------------!
      !     If no cloud is good, we quit before allocating the arrays.                     !
      !------------------------------------------------------------------------------------!
      if (nsolv == 0) exit queq_loop

      !----- Allocate some vectors we need for solving the linear system. -----------------!
      allocate(kke    (nsolv,nsolv))
      allocate(diagkke      (nsolv))
      allocate(mfo          (nsolv))
      allocate(mb           (nsolv))
      allocate(cldidx       (nsolv))

      !----- Store the actual cloud number of the clouds that exist. ----------------------!
      cldidx = pack(cloud,okcld)

      do isol=1,nsolv
         icld = cldidx(isol)
         !----- MFO is minus the forcing due to large scale. ------------------------------!
         mfo(isol) = ( aatot0(icld) - aatot(icld) ) * dtimei
         !----- KKE is the mass flux kernel. ----------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            kke(isol,jsol) = fac * mfke(icld,jcld)
         end do
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the array of kernels.                                                     !
      !------------------------------------------------------------------------------------!
      call diagon(nsolv,kke,diagkke)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Solve the linear system using a Gaussian elimination method.  The system, as   !
      ! stated in Lord et al. (1982), is : K * Mb = -F.                                    !
      !------------------------------------------------------------------------------------!
      call lisys_solver(nsolv,kke,mfo,mb,is_sing)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------! 
      !     First we check whether the solution is an appropriate one or something went    !
      ! wrong.                                                                             !
      !------------------------------------------------------------------------------------! 
      if (is_sing) then
         !---------------------------------------------------------------------------------!
         !     The matrix is singular or almost singular, so we cannot solve these clouds. !
         ! We can quit this routine after freeing the allocated arrays.                    !
         !---------------------------------------------------------------------------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
         exit queq_loop
      elseif (any(mb < 0.) .or. any(diagkke == 0.)) then
         !---------------------------------------------------------------------------------!
         !     There are some negative members which can't be solved, so we must eliminate !
         ! these and try again.                                                            !
         !---------------------------------------------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)

            if (mb(jsol) < 0. .or. diagkke(jsol) == 0.) then
               okcld(jcld) = .false.
            end if
         end do
         !----- Free memory so it will be ready to be allocated again next time. ----------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     All mass flux terms are zero or positive, we are all set.  Simply copy the  !
         ! mass fluxes to the proper arrays and quit.                                      !
         !---------------------------------------------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            upmf(jcld) = mb(jsol)
            !----- UPMX is the maximum cloud mass flux if it was the only cloud present. --!
            upmx(jcld) = max(0.,mfo(jsol) / kke(jsol,jsol))
         end do
         !----- Free memory before leaving the subroutine. --------------------------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
         exit queq_loop
      end if

   end do queq_loop

   return
end subroutine grell_grell_solver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will solve the mass fluxes of a collection of Grell clouds using     !
! the quasi-equilibrium assumption.  This is what is used by Arakawa-Schubert (1974).  The !
! only difference between this and Grellis that Arakawa-Schubert uses a climatological     !
! value instead of the previous step value.                                                !
!------------------------------------------------------------------------------------------!
subroutine grell_arakschu_solver(nclouds,cldd,clds,mgmzp,dtime,p_cup,clim,whlev,ktop,aatot &
                                ,mfke,ierr,upmf,upmx)
   use grell_coms  , only : pclim  & ! Levels with available climatological aatot
                          , aclim1 & ! Standard cloud work function climatology
                          , aclim2 ! ! Alternative cloud work function climatology
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                            , intent(in)    :: cldd    ! Deepest Grell cloud
   integer                            , intent(in)    :: clds    ! Shallowest Grell cloud
   integer                            , intent(in)    :: nclouds ! # of clouds
   integer                            , intent(in)    :: mgmzp   ! # of vertical levels
   integer                            , intent(in)    :: clim    ! Which climatology to use
   integer                            , intent(in)    :: whlev   ! Which level to use.
   integer, dimension(nclouds)        , intent(in)    :: ktop    ! Cloud top level.
   real                               , intent(in)    :: dtime   ! Current grid time step
   real   , dimension(mgmzp)          , intent(in)    :: p_cup   ! Pressure at cloud levels
   real   , dimension(nclouds)        , intent(in)    :: aatot   ! Cloud work with forcing
   real   , dimension(nclouds,nclouds), intent(in)    :: mfke    ! Mass flux kernel.
   !----- Downdraft is inout because I may not compute it ---------------------------------!
   integer, dimension(nclouds)        , intent(in)    :: ierr    ! Error flag
   real   , dimension(nclouds)        , intent(inout) :: upmf    ! Ref. updraft mass flux
   real   , dimension(nclouds)        , intent(inout) :: upmx    ! Max. updraft mass flux
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: icld    ! Cloud index
   integer                              :: jcld    ! Cloud index too
   integer                              :: isol    ! Another cloud index
   integer                              :: jsol    ! Yet another cloud index
   integer                              :: kclim   ! Climatological level to be used.
   real                                 :: dtimei  ! Scratch, inverse of delta-t
   integer                              :: nsolv   ! # of clouds we are still solving.
   logical                              :: is_sing ! Flag.  The kernel was singular. [T|F]
   !----- Auxilliary variables, containing flags and references. --------------------------!
   logical, dimension(nclouds)          :: okcld   ! This cloud is still ok to use.
   integer, dimension(nclouds)          :: cloud   ! Absolute cloud index.
   !----- Temporary, allocatable arrays. --------------------------------------------------!
   real   , dimension(:,:), allocatable :: kke     ! Subset of mass flux kernel
   real   , dimension(:)  , allocatable :: diagkke ! Diagonal of kke
   real   , dimension(:)  , allocatable :: mfo     ! Minus the large-scale forcing
   real   , dimension(:)  , allocatable :: mb      ! Mass flux
   integer, dimension(:)  , allocatable :: cldidx  ! Cloud index, to copy back to ensemble
   !---------------------------------------------------------------------------------------!


   !----- Initial settings. ---------------------------------------------------------------!
   dtimei = 1. / dtime
   do icld=1,nclouds
      cloud(icld) = icld
   end do

   !----- Discarding non-Grell clouds. ----------------------------------------------------!
   do icld=1,cldd-1
      okcld(icld) = .false.
   end do
   !---------------------------------------------------------------------------------------!
   !    Determine the clouds that may be solved here, and allocate a subset with these     !
   ! clouds only.  We may eliminate some other clouds during this process because they     !
   ! could lead to negative mass fluxes.                                                   !
   !---------------------------------------------------------------------------------------!
   do icld=cldd,clds
      okcld(icld) = ierr(icld) == 0
   end do
   !----- Discarding non-Grell clouds. ----------------------------------------------------!
   do icld=clds+1,nclouds
      okcld(icld) = .false.
   end do

   !----- Assign zero to the fluxes, just in case something prevents the calculation. -----!
   do icld=1,nclouds
      upmf(icld) = 0.
      upmx(icld) = 0.
   end do

   !----- Initialise the cloud arrays. ----------------------------------------------------!
   queq_loop: do
      nsolv = count(okcld)

      !------------------------------------------------------------------------------------!
      !     If no cloud is good, we quit before allocating the arrays.                     !
      !------------------------------------------------------------------------------------!
      if (nsolv == 0) exit queq_loop

      !----- Allocate some vectors we need for solving the linear system. -----------------!
      allocate(kke    (nsolv,nsolv))
      allocate(diagkke      (nsolv))
      allocate(mfo          (nsolv))
      allocate(mb           (nsolv))
      allocate(cldidx       (nsolv))

      !----- Store the actual cloud number of those clouds 
      cldidx = pack(cloud,okcld)

      do isol=1,nsolv
         icld = cldidx(isol)

         !---------------------------------------------------------------------------------!
         !     Finding the reference cloud work function.                                  !
         !---------------------------------------------------------------------------------!
         select case (whlev)
         case (1)
            !----- The closest level. -----------------------------------------------------!
            kclim = minloc(abs(pclim-p_cup(ktop(icld))),dim=1)
         case (2)
            !----- The next-to-closest level. ---------------------------------------------!
            kclim = minloc(abs(pclim-p_cup(ktop(icld))),dim=1)
            kclim = minloc(abs(pclim-p_cup(ktop(icld))),dim=1,mask=(pclim /= pclim(kclim)))
         end select

         !---------------------------------------------------------------------------------!
         !     Finding the large-scale forcing using the level we just found.              !
         !---------------------------------------------------------------------------------!
         select case(clim)
         case (1)
            mfo(isol) = (aclim1(kclim) - aatot(icld)) * dtimei
         case (2)
            mfo(isol) = (aclim2(kclim) - aatot(icld)) * dtimei
         end select

         !----- KKE is the mass flux kernel. ----------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            kke(isol,jsol) = mfke(icld,jcld)
         end do
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the array of kernels.                                                     !
      !------------------------------------------------------------------------------------!
      call diagon(nsolv,kke,diagkke)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Solve the linear system using a Gaussian elimination method.  The system, as   !
      ! stated in Lord et al. (1982), is : K * Mb = -F.                                    !
      !------------------------------------------------------------------------------------!
      call lisys_solver(nsolv,kke,mfo,mb,is_sing)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------! 
      !     First we check whether the solution is an appropriate one or something went    !
      ! wrong.                                                                             !
      !------------------------------------------------------------------------------------! 
      if (is_sing) then
         !---------------------------------------------------------------------------------!
         !     The matrix is singular or almost singular, so we cannot solve these clouds. !
         ! We can quit this routine after freeing the allocated arrays.                    !
         !---------------------------------------------------------------------------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         exit queq_loop
      elseif (any(mb < 0.) .or. any(diagkke == 0.)) then
         !---------------------------------------------------------------------------------!
         !     There are some negative members which can't be solved, so we must eliminate !
         ! these and try again.                                                            !
         !---------------------------------------------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            if (mb(jsol) < 0. .or. diagkke(jsol) == 0. ) then
               okcld(jcld) = .false.
            end if
         end do
         !----- Free memory so it will be ready to be allocated next time. ----------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     All mass flux terms are zero or positive, we are all set.  Simply copy the  !
         ! mass fluxes to the proper arrays and quit.                                      !
         !---------------------------------------------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            upmf(jcld) = mb(jsol)
            !----- UPMX is the maximum cloud mass flux if it was the only cloud present. --!
            upmx(jcld) = max(0.,mfo(jsol) / kke(jsol,jsol))
         end do
         !----- Free memory before leaving the subroutine. --------------------------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
         exit queq_loop
      end if

   end do queq_loop

   return
end subroutine grell_arakschu_solver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will solve the mass fluxes of a collection of Grell clouds using     !
! the instability removal assumption.  This is what is used by Kain-Fristsch (1987) dyna-  !
! mic control.                                                                             !
!------------------------------------------------------------------------------------------!
subroutine grell_inre_solver(nclouds,cldd,clds,tscal,fac,aatot0,mfke,ierr,upmf,upmx)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                            , intent(in)    :: cldd    ! Deepest Grell cloud
   integer                            , intent(in)    :: clds    ! Shallowest Grell cloud
   integer                            , intent(in)    :: nclouds ! # of clouds
   real                               , intent(in)    :: tscal   ! Current grid time step
   real                               , intent(in)    :: fac     ! For ensemble variance
   real   , dimension(nclouds)        , intent(in)    :: aatot0  ! Current cloud work fctn.
   real   , dimension(nclouds,nclouds), intent(in)    :: mfke    ! Mass flux kernel.
   !----- Downdraft is inout because I may not compute it ---------------------------------!
   integer, dimension(nclouds)        , intent(in)    :: ierr    ! Error flag
   real   , dimension(nclouds)        , intent(inout) :: upmf    ! Ref. updraft mass flux
   real   , dimension(nclouds)        , intent(inout) :: upmx    ! Max. updraft mass flux
   !----- Local variables. ----------------------------------------------------------------!
   integer                              :: icld    ! Cloud index
   integer                              :: jcld    ! Cloud index too
   integer                              :: isol    ! Another cloud index
   integer                              :: jsol    ! Yet another cloud index
   real                                 :: tscali  ! Scratch, inverse of (tscal * fac)
   integer                              :: nsolv   ! # of clouds we are still solving.
   logical                              :: is_sing ! Flag.  The kernel was singular. [T|F]
   !----- Auxilliary variables, containing flags and references. --------------------------!
   logical, dimension(nclouds)          :: okcld   ! This cloud is still ok to use.
   integer, dimension(nclouds)          :: cloud   ! Absolute cloud index.
   !----- Temporary, allocatable arrays. --------------------------------------------------!
   real   , dimension(:,:), allocatable :: kke     ! Subset of mass flux kernel
   real   , dimension(:)  , allocatable :: diagkke ! Diagonal of kke
   real   , dimension(:)  , allocatable :: mfo     ! Minus the large-scale forcing
   real   , dimension(:)  , allocatable :: mb      ! Mass flux
   integer, dimension(:)  , allocatable :: cldidx  ! Cloud index, to copy back to ensemble
   !---------------------------------------------------------------------------------------!


   !----- Initial settings. ---------------------------------------------------------------!
   tscali = 1. / (fac * tscal)
   do icld=1,nclouds
      cloud(icld) = icld
   end do

   !----- Discarding non-Grell clouds. ----------------------------------------------------!
   do icld=1,cldd-1
      okcld(icld) = .false.
   end do
   !---------------------------------------------------------------------------------------!
   !    Determine the clouds that may be solved here, and allocate a subset with these     !
   ! clouds only.  We may eliminate some other clouds during this process because they     !
   ! could lead to negative mass fluxes.                                                   !
   !---------------------------------------------------------------------------------------!
   do icld=cldd,clds
      okcld(icld) = ierr(icld) == 0
   end do
   !----- Discarding non-Grell clouds. ----------------------------------------------------!
   do icld=clds+1,nclouds
      okcld(icld) = .false.
   end do

   !----- Assign zero to the fluxes, just in case something prevents the calculation. -----!
   do icld=1,nclouds
      upmf(icld) = 0.
      upmx(icld) = 0.
   end do

   !----- Initialise the cloud arrays. ----------------------------------------------------!
   inre_loop: do
      nsolv = count(okcld)

      !------------------------------------------------------------------------------------!
      !     If no cloud is good, we quit before allocating the arrays.                     !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     If no cloud is good, assign non-zero errors to all clouds, then leave.         !
      !------------------------------------------------------------------------------------!
      if (nsolv == 0) exit inre_loop

      !----- Allocate some vectors we need for solving the linear system. -----------------!
      allocate(kke    (nsolv,nsolv))
      allocate(diagkke      (nsolv))
      allocate(mfo          (nsolv))
      allocate(mb           (nsolv))
      allocate(cldidx       (nsolv))
      !------------------------------------------------------------------------------------!

      !----- Store the actual cloud number of those clouds 
      cldidx = pack(cloud,okcld)

      do isol=1,nsolv
         icld = cldidx(isol)
         !----- MFO is minus the forcing due to large scale. ------------------------------!
         mfo(isol) = - aatot0(cldidx(isol))  * tscali
         !----- KKE is the mass flux kernel. ----------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            kke(isol,jsol) = fac * mfke(cldidx(isol),cldidx(jsol))
         end do
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the array of kernels.                                                     !
      !------------------------------------------------------------------------------------!
      call diagon(nsolv,kke,diagkke)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Solve the linear system using a Gaussian elimination method.  The system, as   !
      ! stated in Lord et al. (1982), is : K * Mb = -F.                                    !
      !------------------------------------------------------------------------------------!
      call lisys_solver(nsolv,kke,mfo,mb,is_sing)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------! 
      !     First we check whether the solution is an appropriate one or something went    !
      ! wrong.                                                                             !
      !------------------------------------------------------------------------------------! 
      if (is_sing) then
         !---------------------------------------------------------------------------------!
         !     The matrix is singular or almost singular, so we cannot solve these clouds. !
         ! We can quit this routine after freeing the allocated arrays.                    !
         !---------------------------------------------------------------------------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         exit inre_loop
      elseif (any(mb < 0.) .or. any(diagkke == 0.)) then
         !---------------------------------------------------------------------------------!
         !     There are some negative members or the diagonal of the kernel has some      !
         ! zeroes, which can't be solved, so we must eliminate these and try again.        !
         !---------------------------------------------------------------------------------!
         do jsol=1,nsolv
            jcld = cldidx(jsol)
            if (mb(jsol) < 0. .or. diagkke(jsol) == 0.) then
               okcld(jcld) = .false.
            end if
         end do
         !----- Free memory so it will be ready to be allocated next time. ----------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     All mass flux terms are zero or positive, we are all set.  Simply copy the  !
         ! mass fluxes to the proper arrays and quit.                                      !
         !---------------------------------------------------------------------------------!
         do jsol=1,nsolv
            jcld       = cldidx(jsol)
            upmf(jcld) = mb(jsol)
            !----- UPMX is the maximum cloud mass flux if it was the only cloud present. --!
            upmx(jcld) = max(0.,mfo(jsol) / kke(jsol,jsol))
         end do
         !----- Free memory before leaving the subroutine. --------------------------------!
         deallocate(kke    )
         deallocate(diagkke)
         deallocate(mfo    )
         deallocate(mb     )
         deallocate(cldidx )
         !---------------------------------------------------------------------------------!
         exit inre_loop
      end if

   end do inre_loop

   return
end subroutine grell_inre_solver
!==========================================================================================!
!==========================================================================================!
