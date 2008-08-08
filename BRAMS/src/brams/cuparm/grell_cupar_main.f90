!==========================================================================================!
! grell_cupar_main.f90                                                                     !
!                                                                                          !
!    This file contains the main Grell's subroutine. This subroutine is the same for both  !
! shallow and deep convection, the only difference is that shallow convection bypasses     !
! everything related to downdrafts.                                                        !
!    Whenever applicable, the follow convention is used:                                   !
!       [+] "?_cup"  : These are variables defined at the Grell level (staggered).         !
!       [+] "?d_cld" : These are cloud variables associated with downdrafts.               !
!       [+] "?u_cld" : These are cloud variables associated with updrafts.                 !
!       [+] "x_?"    : Variables modified by a fudging factor, for ensemble statistics     !
!       [+] "?0?"    : Whenever a variable like the ones above contains a 0, it means that !
!                      these variables are based on current values. If this is a BRAMS     !
!                      variable, it means the value before applying the tendency,          !
!                      otherwise it is a derived value based on the current large scale    !
!                      These are used only if Grell (1993) or Kain-Fritsch (1990)          !
!                      parameterizations are used.                                         !
!                                                                                          !
!------------------------------------------------------------------------------------------!
!==========================================================================================!
!==========================================================================================!
subroutine grell_cupar_main(closure_type,comp_down,comp_noforc_cldwork,comp_modif_thermo   &
                           ,checkmass,iupmethod,maxens_cap,maxens_dyn,maxens_eff           &
                           ,maxens_lsf,mgmzp,cap_maxs,cap_max_increment,dtime,depth_min    &
                           ,depth_max,edtmax,edtmin,inv_ensdim,masstol,max_heat,pmass_left &
                           ,radius,relheight_down,zkbmax,zcutdown,z_detr,edt_eff           &
                           ,dellahe_eff,dellaq_eff,dellaqc_eff,dellat_eff,pw_eff,dnmf_ens  &
                           ,upmf_ens,aad,aau,edt,dnmf,upmf,mynum)
               
   use mem_scratch_grell, only : &
      !----- Variables defined at the initialization process ------------------------------!
       dens_curr          & ! intent(in)  - Current density strenghtening updraft
      ,dzd_cld            & ! intent(in)  - Delta-height for downdraft calculations
      ,dzu_cld            & ! intent(in)  - Delta-height for updraft calculations
      ,he0                & ! intent(in)  - Current Moist static energy
      ,he                 & ! intent(in)  - Moist static energy with forcing
      ,hesur              & ! intent(in)  - Surface moist static energy
      ,hes0               & ! intent(in)  - Current Sat. moist static energy
      ,hes                & ! intent(in)  - Sat. moist static energy with forcing
      ,hessur             & ! intent(in)  - Surface sat. moist static energy
      ,kpbl               & ! intent(in)  - Level of PBL top (Nakanishi/Niino only)
      ,mkx                & ! intent(in)  - Number of vertical levels
      ,mconv              & ! intent(in)  - Integrated moisture convergence
      ,omeg               & ! intent(in)  - Vertical velocity in pressure, dp/dt: [Pa/s]
      ,p0                 & ! intent(in)  - Current pressure
      ,p                  & ! intent(in)  - Pressure with forcing
      ,psur               & ! intent(in)  - Pressure for no convection
      ,prev_dnmf          & ! intent(in)  - Previous call downdraft mass flux
      ,q0                 & ! intent(in)  - Current water vapour mixing ratio
      ,q                  & ! intent(in)  - Water vapour mixing ratio with forcing.
      ,qsur               & ! intent(in)  - surface water mixing ratio
      ,qes0               & ! intent(in)  - Current Saturation mixing ratio
      ,qes                & ! intent(in)  - Saturation mixing ratio with forcing.
      ,qessur             & ! intent(in)  - Surface saturation mixing ratio
      ,rcpg               & ! intent(in)  - Large-scale liquid water mixing ratio 
      ,t0                 & ! intent(in)  - Current Temperature
      ,t                  & ! intent(in)  - Temperature with forcing.
      ,tkeg               & ! intent(in)  - Turbulent Kinetic Energy
      ,tsur               & ! intent(in)  - Surface temperature
      ,tscal_kf           & ! intent(in)  - Time scale for Kain-Fritsch (1990)
      ,upconv             & ! intent(in)  - Flag: is convection happening upstream?
      ,uwind              & ! intent(in)  - Zonal wind
      ,vwind              & ! intent(in)  - Meridional wind
      ,z                  & ! intent(in)  - Height
      ,z_cup              & ! intent(in)  - Height at cloud levels
      !------ Variables defined here ------------------------------------------------------!
      ,cdd                & ! intent(out) - Normalized downdraft detr. function   [   ----]
      ,cdu                & ! intent(out) - Normalized updraft detr. function     [   ----]
      ,etad_cld           & ! intent(out) - Normalized downdraft mass flux        [   ----]
      ,etau_cld           & ! intent(out) - Normalized updraft Mass flux          [   ----]
      ,ierr               & ! intent(out) - Flag for convection error             [   ----]
      ,jmin               & ! intent(out) - Downdraft originating level           [   ----]
      ,k22                & ! intent(out) - Updraft originating level             [   ----]
      ,kbcon              & ! intent(out) - Cloud base                            [   ----]
      ,kdet               & ! intent(out) - Top of downdraft detrainment layer    [   ----]
      ,kstabi             & ! intent(out) - Stable layer bottom                   [   ----]
      ,kstabm             & ! intent(out) - Stable layer top                      [   ----]
      ,ktop               & ! intent(out) - Cloud top                             [   ----]
      ,mentrd_rate        & ! intent(out) - Normalized downdraft entrainment rate [   ----]
      ,mentru_rate        & ! intent(out) - Normalized updraft entrainment rate   [   ----]
      ,outt               & ! intent(out) - Temperature tendency                  [    K/s]
      ,outq               & ! intent(out) - Water vapour mixing ratio tendency    [kg/kg/s]
      ,outqc              & ! intent(out) - Condensed water mixing ratio tendency [kg/kg/s]
      ,p_cup              & ! intent(out) - Pressure with forcing @ cloud levels  [     Pa] 
      ,precip             & ! intent(out) - Precipitation rate                    [kg/m²/s]
      ,q_cup              & ! intent(out) - Mix. ratio w/ forcing @ cloud levels  [  kg/kg] 
      ,qd_cld             & ! intent(out) - Downdraft water vapour mixing ratio   [  kg/kg]
      ,qu_cld             & ! intent(out) - Updraft water vapour mixing ratio     [  kg/kg]
      ,qld_cld            & ! intent(out) - Downdraft liquid water mixing ratio   [  kg/kg]
      ,qlu_cld            & ! intent(out) - Updraft liquid water mixing ratio     [  kg/kg]
      ,rho                & ! intent(in)  - Air density                           [  kg/m³]
      ,t_cup              & ! intent(out) - Temperature w/ forcing @ cloudlevels  [      K] 
      ,td_cld             & ! intent(out) - Downdraft temperature                 [      K]
      ,tu_cld             ! ! intent(out) - Updraft temperature                   [      K]

   use rconstants, only: cpi,alvl,g

   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! List of arguments                                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Character, containing the closure type(s) for dynamic control -------------------!
   character(len=2), intent(in) :: closure_type! Short name to define the method.
   !----- Logical flags, to bypass uncessary steps ----------------------------------------!
   logical , intent(in) :: comp_down           ! I will compute downdrafts.           [T/F]
   logical , intent(in) :: comp_noforc_cldwork ! I will compute no forced cloud work  [T/F]
   logical , intent(in) :: comp_modif_thermo   ! I will compute the x_* variables     [T/F]
   logical , intent(in) :: checkmass           ! I will check mass balance            [T/F]
   !----- Integer, global dimensions ------------------------------------------------------!
   integer , intent(in) :: iupmethod        ! Method to find the updraft originatin level
   integer , intent(in) :: maxens_cap       ! Ensemble size on static control (cap_maxs)
   integer , intent(in) :: maxens_dyn       ! Ensemble size on dynamic control 
   integer , intent(in) :: maxens_eff       ! Ensemble size on precipitation efficiency
   integer , intent(in) :: maxens_lsf       ! Ensemble size on large scale perturbations
   integer , intent(in) :: mgmzp            ! Vertical grid size
   integer , intent(in) :: mynum            ! Node ID, for debugging purposes only

   !----- Real, parameters to define the cloud --------------------------------------------!
   real    , intent(in) :: cap_maxs         ! Maximum depth of capping inversion     [ hPa]
   real    , intent(in) :: cap_max_increment! Extra cap_maxs due to upstream conv.   [ hPa]
   real    , intent(in) :: dtime            ! Time step                              [   s]
   real    , intent(in) :: depth_min        ! Minimum cloud depth to qualify it      [   m]
   real    , intent(in) :: depth_max        ! Maximum cloud depth to qualify it      [   m]
   real    , intent(in) :: edtmax           ! Maximum epsilon (dnmf/upmf)            [ ---]
   real    , intent(in) :: edtmin           ! Minimum epsilon (dnmf/upmf)            [ ---]
   real    , intent(in) :: inv_ensdim       ! Inverse of ensemble dimension size     [ ---]
   real    , intent(in) :: masstol          ! Maximum mass leak allowed to happen    [ ---]
   real    , intent(in) :: max_heat         ! Maximum heating scale                  [K/dy]
   real    , intent(in) :: pmass_left       ! Fraction of mass left at the ground    [ ---]
   real    , intent(in) :: radius           ! Radius, for entrainment rate.          [   m]
   real    , intent(in) :: relheight_down   ! Relative height for downdraft origin   [ ---]
   real    , intent(in) :: zkbmax           ! Top height for updrafts to originate   [   m]
   real    , intent(in) :: zcutdown         ! Top height for downdrafts to originate [   m]
   real    , intent(in) :: z_detr           ! Top height for downdraft detrainment   [   m]

   !----- Ensemble variable, fraction between reference downdraft and updraft -------------!
   real, dimension(maxens_eff,maxens_cap)      , intent(out) :: & 
                           edt_eff          ! Grell's epsilon
   !----- Ensemble variables, changes of the thermodynamic properties per unit of mass. ---!
   real, dimension(mgmzp,maxens_eff,maxens_cap), intent(out) :: &
                           dellahe_eff    & ! Moist static energy
                          ,dellaq_eff     & ! Vapour mixing ratio
                          ,dellaqc_eff    & ! Condensed mixing ratio
                          ,dellat_eff     & ! Temperature
                          ,pw_eff         ! ! Fall-out water
   !----- Ensemble variables, reference mass fluxes, in [kg/m²/s] -------------------------!
   real, dimension(maxens_dyn,maxens_lsf,maxens_eff,maxens_cap), intent(out) :: &
                           dnmf_ens       & ! reference downdraft flux
                          ,upmf_ens       ! reference updraft flux
   !----- Output variables, the reference mass fluxes. Others will be found afterwards. ---!
   real                  , intent(out)   :: aad  ! Downdraft work function        [   J/kg]
   real                  , intent(out)   :: aau  ! Updraft work function          [   J/kg]
   real                  , intent(out)   :: edt  ! dnmf/upmf                      [    ---]
   real                  , intent(out)   :: dnmf ! Reference downward mass flux   [kg/m²/s]
   real                  , intent(out)   :: upmf ! Reference upward mass flux     [kg/m²/s]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Local variables.                                                                   !
   !---------------------------------------------------------------------------------------!
   !----- Counters and level holders ------------------------------------------------------!
   integer                   :: icap          ! Static control (cap_max) counter;
   integer                   :: iedt          ! Precipitation efficiency counter;
   integer                   :: imbp          ! Response variable fudge counter;
   integer                   :: k             ! Generic vertical level counter;
   integer                   :: kbmin         ! Auxiliary variable;
   integer                   :: kbmax         ! Top level that updrafts can generate;
   integer                   :: kzdown        ! Actual maximum origin for downdrafts;
   integer                   :: capoffset     ! Offset for perturbation of cap_max
   integer                   :: mboffset      ! Offset for the perturbation of forcing;

   !---------------------------------------------------------------------------------------!
   !     Environment variables at Grell's staggered grid.                                  !
   !---------------------------------------------------------------------------------------!
   !----- Variables based on values before any forcing is applied -------------------------!
   real   , dimension(mgmzp) :: gamma0_cup    ! Gamma                              [      ]
   real   , dimension(mgmzp) :: he0_cup       ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: hes0_cup      ! Saturation moist static energy     [  J/kg]
   real   , dimension(mgmzp) :: p0_cup        ! Pressure                           [    Pa]
   real   , dimension(mgmzp) :: q0_cup        ! Mixing ratio                       [ kg/kg]
   real   , dimension(mgmzp) :: qes0_cup      ! Saturation mixing ratio            [ kg/kg]
   real   , dimension(mgmzp) :: t0_cup        ! Temperature                        [     K]
   !----- Variables with all time tendencies but this convection. -------------------------!
   real   , dimension(mgmzp) :: gamma_cup     ! Gamma                              [   ---]
   real   , dimension(mgmzp) :: he_cup        ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: hes_cup       ! Saturation moist static energy     [  J/kg]
   !real  , dimension(mgmzp) :: p_cup         ! This is at mem_scratch_grell module
   !real  , dimension(mgmzp) :: q_cup         ! This is at mem_scratch_grell module
   real   , dimension(mgmzp) :: qes_cup       ! Saturation mixing ratio            [ kg/kg]
   !real  , dimension(mgmzp) :: t_cup         ! This is at mem_scratch_grell module
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Downdraft-related variables                                                        !
   !---------------------------------------------------------------------------------------!
   !----- Based on variables without forcing (current) ------------------------------------!
   real   , dimension(mgmzp) :: cd0d          ! Normalized detrainment flux        [   1/m]
   real   , dimension(mgmzp) :: dby0d         ! Downdraft buoyancy term            [  J/kg]
   real   , dimension(mgmzp) :: eta0d_cld     ! Normalized downdraft mass flux     [   ---]
   real   , dimension(mgmzp) :: he0d_cld      ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: pw0d_cld      ! Downdraft evaporation              [ kg/kg]
   real   , dimension(mgmzp) :: q0d_cld       ! Mixing ratio                       [ kg/kg]
   real   , dimension(mgmzp) :: ql0d_cld      ! Condensate mixing ratio            [ kg/kg]
   real   , dimension(mgmzp) :: qtd0_cld      ! Total mixing ratio                 [ kg/kg]
   real                      :: aa0d          ! Downdraft work function            [  J/kg]
   real                      :: pwev0         ! Integrated evaporation             [ kg/kg]
   !----- Based on values with forcing ----------------------------------------------------!
   real   , dimension(mgmzp) :: dbyd          ! Downdraft buoyancy term            [  J/kg]
   real   , dimension(mgmzp) :: hed_cld       ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: pwd_cld       ! Downdraft evaporation              [ kg/kg]
   real   , dimension(mgmzp) :: qtd_cld       ! Total mixing ratio                 [ kg/kg]
   real                      :: pwev          ! Integrated evaporation             [ kg/kg]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Updraft-related variables                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real   , dimension(mgmzp) :: cd0u          ! Normalized detrainment flux        [   1/m]
   real   , dimension(mgmzp) :: dby0u         ! Buoyancy term                      [  J/kg]
   real   , dimension(mgmzp) :: eta0u_cld     ! Normalized updraft mass flux       [   ---]
   real   , dimension(mgmzp) :: he0u_cld      ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: pw0u_cld      ! Condensation                       [ kg/kg]
   real   , dimension(mgmzp) :: q0u_cld       ! Vapour mixing ratio                [ kg/kg]
   real   , dimension(mgmzp) :: ql0u_cld      ! Liquid water mixing ratio          [ kg/kg]
   real   , dimension(mgmzp) :: qt0u_cld      ! Total water mixing ratio           [ kg/kg]
   real                      :: aa0u          ! Updraft work function              [  J/kg]
   real                      :: pwav0         ! Integrated condensation            [ kg/kg]
   !----- Based on values with forcing ----------------------------------------------------!
   real   , dimension(mgmzp) :: dbyu          ! Buoyancy term                      [  J/kg]
   real   , dimension(mgmzp) :: heu_cld       ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: pwu_cld       ! Condensation                       [ kg/kg]
   real   , dimension(mgmzp) :: qtu_cld       ! Total water mixing ratio           [     K]
   real                      :: pwav          ! Integrated condensation            [ kg/kg]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Combined updraft/downdraft                                                         !
   !---------------------------------------------------------------------------------------!
   real                      :: aatot0        ! Current total cloud work function  [  J/kg]
   real                      :: aatot         ! Total cloud work function w/ forc. [  J/kg]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Levels for the non-forced variables                                                !
   !---------------------------------------------------------------------------------------!
   integer :: ierr0     ! Flag for convection error
   integer :: jmin0     ! Downdraft originating level
   integer :: k220      ! Updraft originating level
   integer :: kbcon0    ! Cloud base
   integer :: kdet0     ! Top of downdraft detrainment
   integer :: kstabi0   ! Stable layer bottom
   integer :: kstabm0   ! Stable layer top
   integer :: ktop0     ! Cloud top
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Variables modified by the arbitrary mass flux                                      !
   !---------------------------------------------------------------------------------------!
   !----- Model levels --------------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_he          ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: x_hes         ! Saturation moist static energy     [  J/kg]
   real   , dimension(mgmzp) :: x_p           ! Pressure                           [    Pa]
   real   , dimension(mgmzp) :: x_q           ! Vapour mixing ratio                [ kg/kg]
   real   , dimension(mgmzp) :: x_qes         ! Saturation vapour mixing ratio     [ kg/kg]
   real   , dimension(mgmzp) :: x_ql          ! Condensed phase mixing ratio       [ kg/kg]
   real   , dimension(mgmzp) :: x_t           ! Temperature                        [     K]
   !----- Cloud levels --------------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_gamma_cup   ! Gamma                              [   ---]
   real   , dimension(mgmzp) :: x_he_cup      ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: x_hes_cup     ! Saturation moist static energy     [  J/kg]
   real   , dimension(mgmzp) :: x_p_cup       ! Pressure                           [    Pa]
   real   , dimension(mgmzp) :: x_q_cup       ! Vapour mixing ratio after cumulus  [ kg/kg]
   real   , dimension(mgmzp) :: x_qes_cup     ! Saturation vapour mixing ratio     [ kg/kg]
   real   , dimension(mgmzp) :: x_t_cup       ! Temperature                        [     K]
   !----- Downdraft variables -------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_cdd         ! Normalized detrainment flux        [   1/m]
   real   , dimension(mgmzp) :: x_dbyd        ! Downdraft buoyancy term            [  J/kg]
   real   , dimension(mgmzp) :: x_etad_cld    ! Normalized downdraft mass flux     [   ---]
   real   , dimension(mgmzp) :: x_hed_cld     ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: x_pwd_cld     ! Downdraft evaporation              [ kg/kg]
   real   , dimension(mgmzp) :: x_qd_cld      ! Mixing ratio                       [ kg/kg]
   real   , dimension(mgmzp) :: x_qld_cld     ! Condensate mixing ratio            [ kg/kg]
   real   , dimension(mgmzp) :: x_qtd_cld     ! Total mixing ratio                 [ kg/kg]
   real                      :: x_aad         ! Downdraft work function            [  J/kg]
   real                      :: x_pwev        ! Integrated evaporation             [ kg/kg]
   !----- Updraft variables --------------------------------------------------------------!
   real   , dimension(mgmzp) :: x_cdu         ! Normalized detrainment flux        [   1/m]
   real   , dimension(mgmzp) :: x_dbyu        ! Buoyancy term                      [  J/kg]
   real   , dimension(mgmzp) :: x_etau_cld    ! Normalized updraft mass flux       [   ---]
   real   , dimension(mgmzp) :: x_heu_cld     ! Moist static energy                [  J/kg]
   real   , dimension(mgmzp) :: x_pwu_cld     ! Condensation                       [ kg/kg]
   real   , dimension(mgmzp) :: x_qu_cld      ! Vapour mixing ratio                [ kg/kg]
   real   , dimension(mgmzp) :: x_qlu_cld     ! Liquid water mixing ratio          [ kg/kg]
   real   , dimension(mgmzp) :: x_qtu_cld     ! Total water mixing ratio           [ kg/kg]
   real                      :: x_aau         ! Updraft work function              [  J/kg]
   real                      :: x_pwav        ! Integrated condensation            [ kg/kg]
   !----- Combined variables --------------------------------------------------------------!
   real                      :: x_aatot       ! Combined cloud work function       [  J/kg]
   !----- Other parameters ----------------------------------------------------------------!
   integer :: x_ierr      ! Flag for convection error
   integer :: x_jmin      ! Downdraft originating level
   integer :: x_k22       ! Updraft originating level
   integer :: x_kbcon     ! Cloud base
   integer :: x_kdet      ! Top of downdraft detrainment
   integer :: x_kstabi    ! Stable layer bottom
   integer :: x_kstabm    ! Stable layer top
   integer :: x_ktop      ! Cloud top
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Miscellaneous parameters                                                           !
   !---------------------------------------------------------------------------------------!
   real :: cap_max  ! Actual "Depth" of inversion capping                         [     Pa]
   real :: mbprime  ! Current arbitrary mass flux fudging factor for ensemble     [kg/m²/s]
   real :: one_b    ! (1-b), as defined in Krishnamurti et al. (1983)             [    ---]
   real :: zktop    ! Actual maximum height that downdrafts can originate.        [      m]
   !----- Scratch array -------------------------------------------------------------------!
   real, dimension(mgmzp) :: scrvar  ! Scratch variable
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Initialize reference downdraft and updraft reference.                              !
   !---------------------------------------------------------------------------------------!
   dnmf = 0.
   upmf = 0.


   !---------------------------------------------------------------------------------------!
   !    This is the beginning of the outermost loop, on static control. Currently it is    !
   ! set for cap_maxs sensitivity but any other parameter could be used. Also in case you  !
   ! feel like including more terms, it is just a matter of including more loops.          !
   !                                                                                       !
   !    IMPORTANT: This loop is done from the end to the first element because for all     !
   !               ensemble tests, 1 is the reference, so to avoid storing larger          !
   !               dimension variables, I simply start from the end and when I am done     !
   !               the variables will have the values I need.                              !
   !---------------------------------------------------------------------------------------!
   stacloop: do icap=maxens_cap,1,-1
   
      !------------------------------------------------------------------------------------!
      ! A. Define cap_max for this member. The reference given by the user will be used as !
      !    reference, but the actual cap_max will be lowered in case upstream convection   !
      !    has happened and the user asked to consider upstream convection. Also, cap_max  !
      !    will be perturbed for members other than the first. The value for use in the    !
      !    subroutine will be converted to Pa.                                             !
      !------------------------------------------------------------------------------------!
      capoffset= (-1)**icap * icap / 2 !----- This gives 0,1,-1,2,-2 for icap=1,2,3... ----!
      if (upconv) then
         cap_max = 100. * max(cap_max_increment,cap_maxs + (capoffset-1)*cap_max_increment)
      else
         cap_max = 100. * max(cap_max_increment,cap_maxs + capoffset * cap_max_increment)
      end if
      !------------------------------------------------------------------------------------!
      ! B. This flag states whether convection happened or failed. If ierr=0, then         !
      !    convection  happened, otherwise something kept it to develop. Whenever the      !
      !    parameterization finds a condition to impeach convection, a different number is !
      !    assigned and the static control is cycled. The error table is the following:    !
      !                                                                                    !
      !  1. Too early in the simulation for convection to be called.                       !
      !  2. The level updrafts origin is way too high;                                     !
      !  3. Inversion capping may be too thick so the cloud base is probably out of reach; !
      !  4. Couldn't find a layer with negative buoyancy for downdrafts;                   !
      !  5. Couldn't find a suitable cloud top, it would be above the domain top;          !
      !  6. This cloud would be too thin to fall in this spectral type;                    !
      !  7. This cloud would be too thick to fall in this spectral type;                   !
      !  8. Forced downdraft layer would have positive buoyancy, which doesn't make sense; !
      !  9. Cloud work function associated with updraft is zero.                           !
      ! 10. Reference upward mass flux is zero.                                            !
      !------------------------------------------------------------------------------------!
      ierr   = 0

      !------------------------------------------------------------------------------------!
      ! C. Initialize Entrainment and Detrainment variables.                               !
      !------------------------------------------------------------------------------------!
      mentrd_rate(1:mkx) = .2/radius          ! Entrainment rate associated with dndrafts
      mentru_rate(1:mkx) = mentrd_rate(1:mkx) ! Entrainment rate associated with updrafts
      cdd(1:mkx) = 0.                         ! Detrainment function associated w/ dndrafts
      cdu(1:mkx) = 0.1*mentru_rate(1:mkx)     ! Detrainment function associated w/ updrafts

      !------------------------------------------------------------------------------------!
      ! D. Initialize mass fluxes, since they may never be computed in case convection     !
      !    fails.                                                                          !
      !------------------------------------------------------------------------------------!
      etad_cld(1:mkx) = 0.
      etau_cld(1:mkx) = 0.

      !------------------------------------------------------------------------------------!
      ! E. Initialize Cloud work functions associated with updrafts and downdrafts.        !
      !------------------------------------------------------------------------------------!
      aau    = 0.
      aad    = 0.
      aatot  = 0.
      edt    = 0.
   
      !------------------------------------------------------------------------------------!
      ! F. Initialize all level-related cloud variables. They are all output variables and !
      !    some may never be assigned when convection fails. Except for kstabm, I am set-  !
      !    ting all to a non-sense value to make the point that convection did not happen. !
      !------------------------------------------------------------------------------------!
      jmin   = 0
      k22    = 0
      kbcon  = 0
      kdet   = 0
      ktop   = 0
      kstabi = 0
      kstabm = mkx -2

      !------------------------------------------------------------------------------------!
      ! G. Calculate all thermodynamic properties at the cloud level. The cloud levels are !
      !    staggered in relation to BRAMS model.                                           !
      !------------------------------------------------------------------------------------!
      call grell_thermo_cldlev(mkx,mgmzp,z_cup,t,q,p,tsur,qessur,qsur,hessur,hesur,psur    &
                              ,t_cup,qes_cup,q_cup,hes_cup,he_cup,p_cup,gamma_cup          )


      !------------------------------------------------------------------------------------!
      ! H. Initialize updraft temperature and drafts liquid mixing ratio, in case          !
      !    convection fails. These variables are used in other modules, such as CATT and   !
      !    Harrington.                                                                     !
      !------------------------------------------------------------------------------------!
      td_cld(1:mkx)  = t_cup(1:mkx)
      tu_cld(1:mkx)  = t_cup(1:mkx)
      qd_cld(1:mkx)  = q_cup(1:mkx)
      qu_cld(1:mkx)  = q_cup(1:mkx)
      qld_cld(1:mkx) = 0.
      qlu_cld(1:mkx) = 0.



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! I. Finding some updraft properties, namely the levels in which updrafts originate, !
      !    the cloud base and top, and the updraft-related energy, mass, and moisture      !
      !    properties                                                                      !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 1. Find the topmost level in which updrafts can originate based on this grid.      !
      !------------------------------------------------------------------------------------!
      kbmaxloop: do kbmax=1,mkx
         if (z_cup(kbmax) > zkbmax) exit kbmaxloop
      end do kbmaxloop

      !------------------------------------------------------------------------------------!
      ! 2. Determine a first guess for k22, the level of origin of updrafts. Check if that !
      !    didn't prevent convection to happen.                                            !
      !------------------------------------------------------------------------------------!
      call grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,tkeg,rcpg,he_cup,ierr,k22)
      if (ierr /= 0) cycle stacloop
   
   
      !------------------------------------------------------------------------------------!
      ! 3. Finding the convective cloud base. Two important points here:                   !
      !    a. This call may end up preventing convection, so I must check after the call   !
      !    b. This subroutine may also affect the updraft originating level.               !
      !------------------------------------------------------------------------------------!
      call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,he_cup,hes_cup,p_cup,k22,ierr      &
                                ,kbcon)
      if (ierr /= 0) cycle stacloop

      !------------------------------------------------------------------------------------!
      ! 4. Finding the minimum saturation moist static energy. This will be the bottom of  !
      !    the stable layer.                                                               !
      !------------------------------------------------------------------------------------!
      kstabi=(kbcon - 1) + minloc(hes_cup(kbcon:kstabm),dim=1)

      !------------------------------------------------------------------------------------!
      ! 5. Increasing the detrainment in stable layers provided that there is such layer.  !
      !    this rate increases linearly until a maximum value, currently set to 10 times   !
      !    the entrainment rate.                                                           !
      !------------------------------------------------------------------------------------!
      do k=kstabi,kstabm-1
         cdu(k) = min(10.*mentru_rate(k) , cdu(k-1) + 1.5 * mentru_rate(k))
      end do

      !------------------------------------------------------------------------------------!
      ! 6. Calculate the incloud moist static energy.                                      !
      !------------------------------------------------------------------------------------!
      call grell_he_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,he,he_cup,hes_cup,dzu_cld  &
                           ,heu_cld,dbyu)
                           
      !------------------------------------------------------------------------------------!
      ! 7. Finding the cloud top. Since this may keep convection to happen, I check        !
      !    whether I should move on or break here.                                         !
      !------------------------------------------------------------------------------------!
      call  grell_find_cloud_top(mkx,mgmzp,kbcon,he_cup,heu_cld,dbyu,ierr,ktop)
      if (ierr /= 0) cycle stacloop
      !----- Fixing he0u_cld above cloud top ----------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! 8. Also check whether this cloud qualifies to be in this spectrum size. It need    !
      !    to be thicker than the minimum value provided by the user.                      !
      !------------------------------------------------------------------------------------!
      if (z_cup(ktop)-z_cup(kbcon) < depth_min) then
         ierr = 6
         cycle stacloop
      elseif (z_cup(ktop)-z_cup(kbcon) > depth_max) then
         ierr = 7
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      ! 9. Finding the normalized mass fluxes associated with updrafts. Since we are using !
      !    height-based vertical coordinate, there is no need to find the forced           !
      !    normalized mass flux, they'll be the same, so just copy it afterwards.          !
      !------------------------------------------------------------------------------------!
       call grell_nms_updraft(mkx,mgmzp,k22,kbcon,ktop,mentru_rate,cdu,dzu_cld,etau_cld)

      !------------------------------------------------------------------------------------!
      !10. Finding the moisture profiles associated with updrafts.                         !
      !------------------------------------------------------------------------------------!
      call grell_moist_updraft(comp_down,mkx,mgmzp,k22,kbcon,ktop,radius,q,q_cup,qes_cup   &
                              ,gamma_cup,mentru_rate,cdu,dbyu,dzu_cld,etau_cld,qtu_cld     &
                              ,qlu_cld,qu_cld,pwu_cld,pwav)
      !----- Finding updraft temperature --------------------------------------------------!
      tu_cld(1:ktop) = cpi * (heu_cld(1:ktop) - g * z_cup(1:ktop) - alvl * qu_cld(1:ktop))

      !------------------------------------------------------------------------------------!
      !11. Finding the cloud work function associated with updrafts. If this cloud doesn't !
      !    produce cloud work, break the run, we don't simulate lazy clouds in this model. !
      !------------------------------------------------------------------------------------!
      call grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,t_cup,gamma_cup,dbyu,dzu_cld         &
                                ,etau_cld,aau)
      if (aau == 0.) then
         ierr = 9
         cycle stacloop
      end if

      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! J. Finding the downdraft counterpart of the above properties, namely where the     !
      !    downdrafts detrain all its mass, where they originate, their mass, energy and   !
      !    moisture properties. This should be done only when it is a cloud that has       !
      !    downdrafts. !                                                                   !
      !------------------------------------------------------------------------------------!
      if (comp_down) then
         !---------------------------------------------------------------------------------!
         ! 1. The level in which downdrafts should detrain (almost) all its mass.          !
         !---------------------------------------------------------------------------------!
         kdetloop: do kdet=1,mkx
            if (z_cup(kdet) > z_detr) exit kdetloop
         end do kdetloop

         !---------------------------------------------------------------------------------!
         ! 2. The level in which downdrafts will originate in this cloud. Since this may   !
         !    lead to an impossible cloud, check for success.                              !
         !---------------------------------------------------------------------------------!
         call grell_find_downdraft_origin(mkx,mgmzp,k22,ktop,relheight_down,zcutdown,z_cup &
                                         ,hes_cup,dzd_cld,ierr,kdet,jmin)
         if (ierr /= 0) cycle stacloop

         !---------------------------------------------------------------------------------!
         ! 3. Normalized mass fluxes.                                                      !
         !---------------------------------------------------------------------------------!
         call grell_nms_downdraft(mkx,mgmzp,kdet,jmin,.true.,pmass_left,mentrd_rate,cdd    &
                                 ,z_cup,dzd_cld,etad_cld)

         !---------------------------------------------------------------------------------!
         ! 4. Moist static energy. This may prevent convection to happen, check for        !
         !    success.                                                                     !
         !---------------------------------------------------------------------------------!
         call grell_he_downdraft(mkx,mgmzp,jmin,.true.,cdd,mentrd_rate,he,hes_cup,dzd_cld  &
                                ,ierr,hed_cld,dbyd)
         if (ierr /= 0) cycle stacloop

         !---------------------------------------------------------------------------------!
         ! 5. Moisture properties of downdraft                                             !
         !---------------------------------------------------------------------------------!
         if (comp_noforc_cldwork) then
            call grell_moist_downdraft(mkx,mgmzp,jmin,q0,q0_cup,qes0_cup,gamma0_cup        &
                                      ,mentrd_rate,cdd,dby0d,dzd_cld,eta0d_cld,qtd0_cld    &
                                      ,q0d_cld,ql0d_cld,pw0d_cld,pwev0)
         end if
         call grell_moist_downdraft(mkx,mgmzp,jmin,q,q_cup,qes_cup,gamma_cup,mentrd_rate   &
                                   ,cdd,dbyd,dzd_cld,etad_cld,qtd_cld,qd_cld,qld_cld       &
                                   ,pwd_cld,pwev)
         !----- Finding updraft temperature -----------------------------------------------!
         td_cld(1:jmin) = cpi*(hed_cld(1:jmin) - g*z_cup(1:jmin) - alvl*qd_cld(1:jmin))

         !---------------------------------------------------------------------------------!
         ! 6. Compute the downdraft strength in terms of windshear. Remembering that edt   !
         !    is the fraction between downdraft and updraft.                               !
         !---------------------------------------------------------------------------------!
         call grell_efficiency_ensemble(mkx,mgmzp,maxens_eff,k22,kbcon,ktop,edtmin,edtmax  &
                                       ,pwav,pwev,z_cup,uwind,vwind,edt_eff(1,icap))

         !---------------------------------------------------------------------------------!
         ! 7. Computing cloud work function associated with downdrafts                     !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_downdraft(mkx,mgmzp,jmin,t_cup,gamma_cup,dbyd,dzd_cld     &
                                     ,etad_cld,aad)

      else
         !---------------------------------------------------------------------------------!
         !    Now we assume that if a cloud doesn't have dowdrafts, it doesn't precipitate !
         ! either. In the future we may revisit this assumption. For now, if the cloud is  !
         ! shallow, then we set up some of the above parameters for non existent dowdraft  !
         ! (i.e. no mass flux, downdraft thermodynamic variables coincident to cloud       !
         ! levels etc.). kdet and jmin were already initialized to zero.                   !
         !---------------------------------------------------------------------------------!
         etad_cld     = 0.
         hed_cld      = hes_cup
         dbyd         = 0.
         qd_cld       = qes_cup
         pwd_cld      = 0.
         pwev         = 0.
         edt_eff      = 0.
         aad          = 0.
      end if
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!


      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! K. Finding the non-forced cloud work, which will require computing most            !
      !    subroutines again. This will go and compute even if the tests would prevent     !
      !    the cloud to happen. This part will be skipped if the user is asking for moist- !
      !    ure convergence only.                                                           !
      !------------------------------------------------------------------------------------!
      if (comp_noforc_cldwork) then
         !---------------------------------------------------------------------------------!
         ! i.     Initializing levels                                                      !
         !---------------------------------------------------------------------------------!
         ierr0   = 0
         jmin0   = 0
         k220    = 0
         kbcon0  = 0
         kdet0   = kdet
         ktop0   = 0
         kstabi0 = 0
         kstabm0 = kstabm
      
         !---------------------------------------------------------------------------------!
         ! ii.    Initialize detrainment and cloud work variables.                         !
         !---------------------------------------------------------------------------------!
         cd0d(1:mkx) = 0.                     ! Detrainment function associated w/ dndrafts
         cd0u(1:mkx) = 0.1*mentru_rate(1:mkx) ! Detrainment function associated w/ updrafts
         aa0d   = 0.
         aa0u   = 0.
         aatot0 = 0.

         !---------------------------------------------------------------------------------!
         ! iii.   Calculate all thermodynamic properties at the cloud level and initialize !
         !        draft thermodynamic properties                                           !
         !---------------------------------------------------------------------------------!
         call grell_thermo_cldlev(mkx,mgmzp,z_cup,t0,q0,p0,tsur,qessur,qsur,hessur,hesur   &
                                 ,psur,t0_cup,qes0_cup,q0_cup,hes0_cup,he0_cup,p0_cup      &
                                 ,gamma0_cup)


         !---------------------------------------------------------------------------------!
         ! iv.    Determine a first guess for k22, the cloud base and the bottom of the    !
         !        stable layer                                                             !
         !---------------------------------------------------------------------------------!
         call grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,tkeg,rcpg,he0_cup      &
                                  ,ierr0,k220)
         call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,he0_cup,hes0_cup,p0_cup,k220    &
                                   ,ierr0,kbcon0)

         kstabi0=(kbcon0 - 1) + minloc(hes0_cup(kbcon0:kstabm0),dim=1)

         !---------------------------------------------------------------------------------!
         ! v.     Increasing the detrainment in stable layers.                             !
         !---------------------------------------------------------------------------------!
         do k=kstabi0,kstabm0-1
            cd0u(k) = min(10.*mentru_rate(k) , cd0u(k-1) + 1.5 * mentru_rate(k))
         end do
         !---------------------------------------------------------------------------------!
         ! vi .   Calculate the incloud moist static energy.                               !
         !---------------------------------------------------------------------------------!
         call grell_he_updraft(mkx,mgmzp,k220,kbcon0,cd0u,mentru_rate,he0,he0_cup          &
                              ,hes0_cup,dzu_cld,he0u_cld,dby0u)
         !---------------------------------------------------------------------------------!
         ! vii.   Finding the cloud top.                                                   !
         !---------------------------------------------------------------------------------!
         call  grell_find_cloud_top(mkx,mgmzp,kbcon0,he0_cup,he0u_cld,dby0u,ierr0,ktop0)

         !---------------------------------------------------------------------------------!
         ! viii.  Finding the moisture profile.                                            !
         !---------------------------------------------------------------------------------!
         call grell_nms_updraft(mkx,mgmzp,k220,kbcon0,ktop0,mentru_rate,cd0u,dzu_cld       &
                               ,eta0u_cld)
         !---------------------------------------------------------------------------------!
         ! ix.    Finding the moisture profiles associated with updrafts.                  !
         !---------------------------------------------------------------------------------!
         call grell_moist_updraft(comp_down,mkx,mgmzp,k220,kbcon0,ktop0,radius,q0,q0_cup   &
                                 ,qes0_cup,gamma0_cup,mentru_rate,cd0u,dby0u,dzu_cld       &
                                 ,eta0u_cld,qt0u_cld,ql0u_cld,q0u_cld,pw0u_cld,pwav0)
         !---------------------------------------------------------------------------------!
         ! x.     Finding the cloud work function                                          !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_updraft(mkx,mgmzp,kbcon0,ktop0,t0_cup,gamma0_cup,dby0u         &
                                   ,dzu_cld,eta0u_cld,aa0u)
         
         !----- Downdraft properties ------------------------------------------------------!
         if (comp_down) then
            !------------------------------------------------------------------------------!
            ! xi.    The level in which downdrafts will originate in this cloud.           !
            !------------------------------------------------------------------------------!
            call grell_find_downdraft_origin(mkx,mgmzp,k220,ktop0,relheight_down,zcutdown  &
                                            ,z_cup,hes0_cup,dzd_cld,ierr0,kdet0,jmin0)

            !------------------------------------------------------------------------------!
            ! xiii.  Normalized mass fluxes.                                               !
            !------------------------------------------------------------------------------!
            call grell_nms_downdraft(mkx,mgmzp,kdet0,jmin0,.true.,pmass_left,mentrd_rate   &
                                    ,cd0d,z_cup,dzd_cld,eta0d_cld)
            !------------------------------------------------------------------------------!
            ! xiv.   Moist static energy. This may prevent convection to happen.           !
            !------------------------------------------------------------------------------!
            call grell_he_downdraft(mkx,mgmzp,jmin0,.true.,cd0d,mentrd_rate,he0,hes0_cup  &
                                   ,dzd_cld,ierr0,he0d_cld,dby0d)
            !------------------------------------------------------------------------------!
            ! xv.    Moisture properties of downdraft                                      !
            !------------------------------------------------------------------------------!
            call grell_moist_downdraft(mkx,mgmzp,jmin0,q0,q0_cup,qes0_cup,gamma0_cup       &
                                      ,mentrd_rate,cd0d,dby0d,dzd_cld,eta0d_cld,qtd0_cld   &
                                      ,q0d_cld,ql0d_cld,pw0d_cld,pwev0)

            !------------------------------------------------------------------------------!
            ! xvi.   Downdraft cloud work function.                                        !
            !------------------------------------------------------------------------------!
            call grell_cldwork_downdraft(mkx,mgmzp,jmin0,t0_cup,gamma0_cup,dby0d,dzd_cld   &
                                        ,eta0d_cld,aa0d)

         end if

      end if
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! L. Now we perform the big loop across the precipitation efficiency ensembles.      !
      !------------------------------------------------------------------------------------!
      effloop: do iedt=1,maxens_eff
         
         !---------------------------------------------------------------------------------!
         !    Copying the current epsilon to a shortcut. This is an output variable, how-  !
         ! ever, it will be used as scratch here. At the feedback time the value will be   !
         ! overwritten by the actual output value.                                         !
         !---------------------------------------------------------------------------------!
         edt  = edt_eff(iedt,icap)

         !---------------------------------------------------------------------------------!
         ! 1. Compute the total cloud work function for this edt member                    !
         !---------------------------------------------------------------------------------!
         if (comp_noforc_cldwork) aatot0 = aa0u + edt * aa0d
         aatot  = aau + edt * aad


         !---------------------------------------------------------------------------------!
         ! 2. Compute the changes in moisture and energy at the bottom, considering        !
         !    downdraft effects. If this cloud doesn't have downdrafts, it should be zero  !
         !    so skip it.                                                                  !
         !---------------------------------------------------------------------------------!
         if (comp_down) then
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,he,p_cup,he_cup      &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,hed_cld         &
                                         ,dellahe_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,q,p_cup,q_cup        &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,qd_cld          &
                                         ,dellaq_eff(1,iedt,icap))
         end if
         
         !---------------------------------------------------------------------------------!
         ! 3. Compute the changes in energy at the other levels. In case downdrafts are    !
         !    absent, this will work fine too, since the mass fluxes will all be zero.     !
         !---------------------------------------------------------------------------------!
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,he,p_cup,he_cup,mentrd_rate,mentru_rate,cdd,cdu        &
                                   ,dzd_cld,etad_cld,etau_cld,hed_cld,heu_cld              &
                                   ,dellahe_eff(1,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,q,p_cup,q_cup,mentrd_rate,mentru_rate,cdd,cdu,dzd_cld  &
                                   ,etad_cld,etau_cld,qd_cld,qu_cld                        &
                                   ,dellaq_eff(1,iedt,icap))
         
         !---------------------------------------------------------------------------------!
         ! 4. Compute the cloud water tendency, this is done somewhat differently just     !
         !    because no liquid water entrains the cloud, do we should worry about         !
         !    detrainment terms and save some computation                                  !
         !---------------------------------------------------------------------------------!
         call grell_dellacond_ensemble(mgmzp,kbcon,ktop,p_cup,cdu,dzd_cld,etau_cld,qlu_cld &
                                      ,dellaqc_eff(1,iedt,icap))

         !---------------------------------------------------------------------------------!
         ! 5. Computing the derived values. Temperature tendency is derived from moist     !
         !    static energy and vapour tendencies ,and the precipitation is simply the     !
         !    integral of the fallen-out condensated water                                 !
         !---------------------------------------------------------------------------------!
         do k=1,ktop
            !------ Temperature tendency --------------------------------------------------!
            dellat_eff(k,iedt,icap) = cpi*(dellahe_eff(k,iedt,icap)                        &
                                          -alvl*dellaq_eff(k,iedt,icap))
            !------ Precipitation ---------------------------------------------------------!
            pw_eff(k,iedt,icap) =  pwu_cld(k) + edt * pwd_cld(k)
         end do

         !---------------------------------------------------------------------------------!
         ! 6. This is another big loop, across the large scale force members. This is      !
         !    essentially computing the profiles changed by the cumulus cloud, and then    !
         !    "fudging" the variables a little bit, so we assess how sensitive this cloud  !
         !    is.                                                                          !
         !---------------------------------------------------------------------------------!
         mbprimeloop: do imbp=1,maxens_lsf

            !------------------------------------------------------------------------------!
            ! 6a. Initialize some ensemble-related scratch variables. mbprime is an        !
            !     arbitrary extra mass flux that is applied to the forcing, to get the     !
            !     upward mass flux for all cloud work related parametrisations. This is    !
            !     also a variable that we perturb to get different members for the         !
            !     ensemble                                                                 !
            !------------------------------------------------------------------------------!
            mboffset = (-1)**imbp * imbp /2 !----- This gives 0, 1, -1, 2, -2, ... --------!
            mbprime=(4.+real(mboffset))*1.e-3

            !------------------------------------------------------------------------------!
            ! 6b. Defining the Kuo (1974) moistening factor b. Following Krishnamurti      !
            !     et al. (1983), this value is around 0.3, with significant oscillations   !
            !     between 0.0 and 0.8. Here we will play it closer to the average,         !
            !     oscillating it                                                           !
            !     between 0.2 and 0.4, equivalent to oscillate (1-b) between 0.6 and 0.8.  !
            !------------------------------------------------------------------------------!
            one_b = 0.7+0.2*real(mboffset)

            !------------------------------------------------------------------------------!
            ! 6c. Applying the fudging to some thermodynamic variables                     !
            !------------------------------------------------------------------------------!
         
            do k=1,mkx-1
               x_he(k) = he(k)           + mbprime * dtime * dellahe_eff(k,iedt,icap)
               x_q (k) = max(1.e-8,q(k)  + mbprime * dtime * dellaq_eff(k,iedt,icap))
               x_p (k) = p(k)
            end do
            !----- Boundary condition at mkx, just use environment ------------------------!
            x_he(mkx) = he(mkx)
            x_q (mkx) = max(1.e-8,q(mkx))
            x_p (mkx) = p(mkx)

            !------------------------------------------------------------------------------!
            ! 6d. Initialize  variables that may not be assigned in case this member       !
            !     fails.                                                                   !
            !------------------------------------------------------------------------------!
            !----- Indices ----------------------------------------------------------------!
            x_ierr   = 0
            x_jmin   = 0      ! jmin  
            x_k22    = 0      ! k22   
            x_kbcon  = 0      ! kbcon 
            x_kdet   = kdet   ! kdet  
            x_ktop   = 0      ! ktop  
            x_kstabi = 0      ! kstabi
            x_kstabm = kstabm ! kstabm
            !----- Detrainment fluxes -----------------------------------------------------!
            x_cdd(1:mkx) = 0.                    ! cdd(1:mkx) 
            x_cdu(1:mkx) = 0.1*mentru_rate(1:mkx)! cdu(1:mkx) 
            !----- Cloud work -------------------------------------------------------------!
            x_aad    = 0.
            x_aau    = 0.
            !----- Ensemble variables -----------------------------------------------------!
            x_aatot  = 0.
        
            !------------------------------------------------------------------------------!
            ! OBS: Contrary to before, I won't return or cycle in case convection fails    !
            !      here. But I will need to bypass all following steps to ensure that the  !
            !      cloud work will be zero.                                                !
            !------------------------------------------------------------------------------!
            modif_comp_if: if (comp_modif_thermo) then
               !---------------------------------------------------------------------------!
               ! 6e. Compute the modified structure                                        !
               !---------------------------------------------------------------------------!
               call grell_thermo_modlev(mkx,mgmzp,x_p,x_he,z,x_q,x_t,x_ql,x_qes,x_hes)
               call grell_thermo_cldlev(mkx,mgmzp,z_cup,x_t,x_q,x_p,tsur,qessur,qsur       &
                                       ,hessur,hesur,psur,x_t_cup,x_qes_cup,x_q_cup        &
                                       ,x_hes_cup,x_he_cup,x_p_cup,x_gamma_cup)

               !---------------------------------------------------------------------------!
               ! 6f. Finding the cloud levels. If the modified profile makes the cloud     !
               !     impossible, cycle and go to the next member.                          !
               !---------------------------------------------------------------------------!
               !----- New updraft origin --------------------------------------------------!
               call grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,tkeg,rcpg        &
                                        ,x_he_cup,x_ierr,x_k22)
               !----- New cloud base ------------------------------------------------------!
               call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,x_he_cup,x_hes_cup        &
                                         ,x_p_cup,x_k22,x_ierr,x_kbcon)
               !----- New stable layer bottom ---------------------------------------------!
               x_kstabi = (x_kbcon - 1) + minloc(x_hes_cup(x_kbcon:x_kstabm),dim=1)

               !---------------------------------------------------------------------------!
               ! 6g. Compute the modified detrainment profile                              !
               !---------------------------------------------------------------------------!
               do k=x_kstabi,x_kstabm-1
                  x_cdu(k) = min(10.*mentru_rate(k), x_cdu(k-1) + 1.5 * mentru_rate(k))
               end do

               !---------------------------------------------------------------------------!
               ! 6h. Compute the updraft profiles of mass, energy and moisture for this    !
               !     member.                                                               !
               !---------------------------------------------------------------------------!
               call grell_he_updraft(mkx,mgmzp,x_k22,x_kbcon,x_cdu,mentru_rate,x_he        &
                                    ,x_he_cup,x_hes_cup,dzu_cld,x_heu_cld,x_dbyu)
               
               !---------------------------------------------------------------------------!
               ! 6i. Finding the cloud top and checking the new cloud depth                !
               !---------------------------------------------------------------------------!
               !----- Finding the cloud top. ----------------------------------------------!
               call  grell_find_cloud_top(mkx,mgmzp,x_kbcon,x_he_cup,x_heu_cld,x_dbyu      &
                                         ,x_ierr,x_ktop)

               !---------------------------------------------------------------------------!
               ! 6j. Finding the updraft mass flux                                         !
               !---------------------------------------------------------------------------!
               call grell_nms_updraft(mkx,mgmzp,x_k22,x_kbcon,x_ktop,mentru_rate,x_cdu     &
                                     ,dzu_cld,x_etau_cld)

               !---------------------------------------------------------------------------!
               ! 6k. Getting the updraft moisture profile                                  !
               !---------------------------------------------------------------------------!
               call grell_moist_updraft(comp_down,mkx,mgmzp,x_k22,x_kbcon,x_ktop,radius    &
                                       ,x_q,x_q_cup,x_qes_cup,x_gamma_cup,mentru_rate      &
                                       ,x_cdu,x_dbyu,dzu_cld,x_etau_cld,x_qtu_cld          &
                                       ,x_qlu_cld,x_qu_cld,x_pwu_cld,x_pwav)

               !---------------------------------------------------------------------------!
               ! 6l. Recalculating the updraft cloud work                                  !
               !---------------------------------------------------------------------------!
               call grell_cldwork_updraft(mkx,mgmzp,x_kbcon,x_ktop,x_t_cup,x_gamma_cup     &
                                         ,x_dbyu,dzu_cld,x_etau_cld,x_aau)
         
               modif_down: if (comp_down) then
                  !------------------------------------------------------------------------!
                  ! 6m. Finding the new origin of downdrafts                               !
                  !------------------------------------------------------------------------!
                  call grell_find_downdraft_origin(mkx,mgmzp,x_k22,x_ktop,relheight_down   &
                                                  ,zcutdown,z_cup,x_hes_cup,dzd_cld,x_ierr &
                                                  ,x_kdet,x_jmin)

                  !------------------------------------------------------------------------!
                  ! 6n. Finding the normalized mass flux                                   !
                  !------------------------------------------------------------------------!
                  call grell_nms_downdraft(mkx,mgmzp,x_kdet,x_jmin,.true.,pmass_left       &
                                          ,mentrd_rate,x_cdd,z_cup,dzd_cld,x_etad_cld)
                  
                  !------------------------------------------------------------------------!
                  ! 6o. Finding moist static energy                                        !
                  !------------------------------------------------------------------------!
                  call grell_he_downdraft(mkx,mgmzp,x_jmin,.false.,x_cdd,mentrd_rate,x_he  &
                                         ,x_hes_cup,dzd_cld,x_ierr,x_hed_cld,x_dbyd)

                  !------------------------------------------------------------------------!
                  ! 6p. Moisture properties                                                !
                  !------------------------------------------------------------------------!
                  call grell_moist_downdraft(mkx,mgmzp,x_jmin,x_q,x_q_cup,x_qes_cup        &
                                            ,x_gamma_cup,mentrd_rate,x_cdd,x_dbyd,dzd_cld  &
                                            ,x_etad_cld,x_qtd_cld,x_qd_cld,x_qld_cld       &
                                            ,x_pwd_cld,x_pwev)

                  !------------------------------------------------------------------------!
                  ! 6q. Computing cloud work function associated with downdrafts           !
                  !------------------------------------------------------------------------!
                  call grell_cldwork_downdraft(mkx,mgmzp,x_jmin,x_t_cup,x_gamma_cup,x_dbyd &
                                              ,dzd_cld,x_etad_cld,x_aad)
               end if modif_down
            end if modif_comp_if
            !------------------------------------------------------------------------------!
            ! 6r. If the fudged field produced a cloud, we compute the total cloud work    !
            !     function for this particular test. If not, or if the user is using       !
            !     moisture convergence closure only, this will be a silly 0 + 0.           !
            !------------------------------------------------------------------------------!
            x_aatot   = x_aau + edt * x_aad

            !------------------------------------------------------------------------------!
            ! 6s. Running the dynamic control ensemble, using the different combinations   !
            !     of cloud efficiency and fudged variables. The reference upward mass flux !
            !     will then have maxens_eff × maxens_lsf × maxens_dyn different values.    !
            !------------------------------------------------------------------------------!
            call grell_dyncontrol_ensemble(closure_type,comp_down,mgmzp,maxens_dyn,dtime   &
                                          ,tscal_kf,dens_curr,prev_dnmf,mconv,k22,kbcon    &
                                          ,ktop,omeg,x_p_cup,edt,mbprime,one_b,aatot0      &
                                          ,aatot,x_aatot,pwav,pwev                         &
                                          ,dnmf_ens(1,imbp,iedt,icap)                      &
                                          ,upmf_ens(1,imbp,iedt,icap))
         end do mbprimeloop
      end do effloop
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!
   end do stacloop

   !---------------------------------------------------------------------------------------!
   ! 7. Now that all ensemble members are computed, we will proceed with the feedback      !
   !    to the large scale. If ierr is activated, then I don't need to compute this step.  !
   !---------------------------------------------------------------------------------------!
   if (ierr == 0) then 
      call grell_feedback(comp_down,mgmzp,maxens_cap,maxens_eff,maxens_lsf,maxens_dyn      &
                         ,inv_ensdim,max_heat,ktop,edt_eff,dellahe_eff,dellaq_eff          &
                         ,dellaqc_eff,dellat_eff,pw_eff,dnmf_ens,upmf_ens,ierr,upmf,dnmf   &
                         ,edt,outt,outq,outqc,precip)
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine grell_cupar_main
!==========================================================================================!
!==========================================================================================!
