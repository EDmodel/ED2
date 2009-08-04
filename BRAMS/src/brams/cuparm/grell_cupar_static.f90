!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute all the variables that depend on the static control, and !
! create the arrays of ensemble members that depend on precipitation efficiency and the    !
! cap_maxs variation.  This subroutine works for any cloud, and the only difference is     !
! that if this is a shallow cloud, then it will bypass everything related to downdrafts.   !
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
subroutine grell_cupar_static(comp_down,comp_noforc_cldwork,checkmass,iupmethod,maxens_cap &
                             ,maxens_eff,mgmzp,cap_maxs,cap_max_increment,wnorm_max        &
                             ,wnorm_increment,depth_min,depth_max,edtmax,edtmin,masstol    &
                             ,pmass_left,radius,relheight_down,zkbmax,zcutdown,z_detr,ee   &
                             ,aad,aau,mynum)
   use mem_ensemble     , only : &
       ensemble_vars             ! ! structure - The ensemble scratch structure. ----------!

   use mem_scratch_grell, only : &
      !----- Variables defined at the initialization process ------------------------------!
       co20         & ! intent(in)  - Current CO2 mixing ratio                    [    ppm]
      ,co2          & ! intent(in)  - CO2 mixing ratio with forcing.              [    ppm]
      ,co2sur       & ! intent(in)  - surface CO2 mixing ratio                    [    ppm]
      ,dzd_cld      & ! intent(in)  - Delta-height for downdraft calculations     [      m]
      ,dzu_cld      & ! intent(in)  - Delta-height for updraft calculations       [      m]
      ,exner0       & ! intent(in)  - Current Exner function                      [ J/kg/K]
      ,exner        & ! intent(in)  - Forced Exner funtion                        [ J/kg/K]
      ,exnersur     & ! intent(in)  - Surface Exner function                      [ J/kg/K]
      ,kpbl         & ! intent(in)  - Level of PBL top (Nakanishi/Niino only)     [    ---]
      ,ktpse        & ! intent(in)  - Maximum level allowed for cloud top         [    ---]
      ,mkx          & ! intent(in)  - Number of vertical levels                   [    ---]
      ,p0           & ! intent(in)  - Current pressure                            [     Pa]
      ,p            & ! intent(in)  - Pressure with forcing                       [     Pa]
      ,psur         & ! intent(in)  - Surface pressure                            [     Pa]
      ,qtot0        & ! intent(in)  - Current total mixing ratio                  [  kg/kg]
      ,qtot         & ! intent(in)  - Total mixing ratio with forcing.            [  kg/kg]
      ,qtotsur      & ! intent(in)  - surface total mixing ratio                  [  kg/kg]
      ,qvap0        & ! intent(in)  - Current water vapour mixing ratio           [  kg/kg]
      ,qvap         & ! intent(in)  - Water vapour mixing ratio with forcing.     [  kg/kg]
      ,qvapsur      & ! intent(in)  - surface water vapour mixing ratio           [  kg/kg]
      ,qliq0        & ! intent(in)  - Current liquid water mixing ratio           [  kg/kg]
      ,qliq         & ! intent(in)  - Liquid water mixing ratio with forcing.     [  kg/kg]
      ,qliqsur      & ! intent(in)  - surface liquid water mixing ratio           [  kg/kg]
      ,qice0        & ! intent(in)  - Current ice mixing ratio                    [  kg/kg]
      ,qice         & ! intent(in)  - Ice mixing ratio with forcing.              [  kg/kg]
      ,qicesur      & ! intent(in)  - surface ice mixing ratio                    [  kg/kg]
      ,rho          & ! intent(in)  - Air density                                 [  kg/m設
      ,theiv0       & ! intent(in)  - Current ice-vapour equiv. pot. temp.        [      K]
      ,theiv        & ! intent(in)  - Ice-vapour equiv. pot. temp. with forcing   [      K]
      ,theivsur     & ! intent(in)  - Surface ice-vapour equiv. pot. temp.        [      K]
      ,thil0        & ! intent(in)  - Current ice-liquid potential temp.          [      K]
      ,thil         & ! intent(in)  - Ice-liquid potential temp. with forcing     [      K]
      ,thilsur      & ! intent(in)  - Surface ice-liquid potential temp.          [      K]
      ,t0           & ! intent(in)  - Current ice-liquid potential temp.          [      K]
      ,t            & ! intent(in)  - Temperature with forcing.                   [      K]
      ,tsur         & ! intent(in)  - Surface temperature                         [      K]
      ,tke0         & ! intent(in)  - Turbulent Kinetic Energy                    [   J/kg]
      ,tke          & ! intent(in)  - Turbulent Kinetic Energy with forcing       [   J/kg]
      ,sigw         & ! intent(in)  - Vertical velocity standard deviation        [    m/s]
      ,uwind        & ! intent(in)  - Zonal wind                                  [    m/s]
      ,vwind        & ! intent(in)  - Meridional wind                             [    m/s]
      ,wwind        & ! intent(in)  - Mean vertical velocity                      [    m/s]
      ,z            & ! intent(in)  - Height                                      [      m]
      ,z_cup        ! ! intent(in)  - Height at cloud levels                      [      m]

   implicit none
   
   !---------------------------------------------------------------------------------------!
   ! List of arguments                                                                     !
   !---------------------------------------------------------------------------------------!
   !----- Logical flags, to bypass uncessary steps ----------------------------------------!
   logical , intent(in) :: comp_down           ! I will compute downdrafts.           [T/F]
   logical , intent(in) :: comp_noforc_cldwork ! I will compute no forced cloud work  [T/F]
   logical , intent(in) :: checkmass           ! I will check mass balance            [T/F]
   !----- Integer, global dimensions ------------------------------------------------------!
   integer , intent(in) :: iupmethod        ! Method to find the updraft originatin level
   integer , intent(in) :: maxens_cap       ! Ensemble size on static control (cap_maxs)
   integer , intent(in) :: maxens_eff       ! Ensemble size on precipitation efficiency
   integer , intent(in) :: mgmzp            ! Vertical grid size
   integer , intent(in) :: mynum            ! Node ID, for debugging purposes only

   !----- Real, parameters to define the cloud --------------------------------------------!
   real    , intent(in) :: cap_maxs         ! Maximum depth of capping inversion     [ hPa]
   real    , intent(in) :: cap_max_increment! Extra cap_maxs due to upstream conv.   [ hPa]
   real    , intent(in) :: depth_min        ! Minimum cloud depth to qualify it      [   m]
   real    , intent(in) :: depth_max        ! Maximum cloud depth to qualify it      [   m]
   real    , intent(in) :: edtmax           ! Maximum epsilon (dnmf/upmf)            [ ---]
   real    , intent(in) :: edtmin           ! Minimum epsilon (dnmf/upmf)            [ ---]
   real    , intent(in) :: masstol          ! Maximum mass leak allowed to happen    [ ---]
   real    , intent(in) :: pmass_left       ! Fraction of mass left at the ground    [ ---]
   real    , intent(in) :: radius           ! Radius, for entrainment rate.          [   m]
   real    , intent(in) :: relheight_down   ! Relative height for downdraft origin   [ ---]
   real    , intent(in) :: wnorm_max        ! Maximum normalised w to be considered  [ ---]
   real    , intent(in) :: wnorm_increment  ! Extra wnorm for ensemble               [ ---]
   real    , intent(in) :: zkbmax           ! Top height for updrafts to originate   [   m]
   real    , intent(in) :: zcutdown         ! Top height for downdrafts to originate [   m]
   real    , intent(in) :: z_detr           ! Top height for downdraft detrainment   [   m]

   !----- The ensemble scratch structure. -------------------------------------------------!
   type(ensemble_vars)   , intent(inout) :: ee 

   !----- Output variables, the cloud work functions. Others will be found afterwards. ----!
   real                  , intent(out)   :: aad  ! Downdraft work function        [   J/kg]
   real                  , intent(out)   :: aau  ! Updraft work function          [   J/kg]
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
   !---------------------------------------------------------------------------------------!
   !     Environment variables at Grell's staggered grid.                                  !
   !---------------------------------------------------------------------------------------!
   !----- Variables based on values before any forcing is applied -------------------------!
   real   , dimension(mgmzp) :: co20_cup      ! CO2 mixing ratio                  [    ppm]
   real   , dimension(mgmzp) :: exner0_cup    ! Exner function                    [   J/kg]
   real   , dimension(mgmzp) :: p0_cup        ! Pressure                          [     Pa]
   real   , dimension(mgmzp) :: qtot0_cup     ! Total mixing ratio                [  kg/kg]
   real   , dimension(mgmzp) :: qvap0_cup     ! Water vapour mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qliq0_cup     ! Liquid water                      [  kg/kg]
   real   , dimension(mgmzp) :: qice0_cup     ! Mixing ratio                      [  kg/kg]
   real   , dimension(mgmzp) :: qsat0_cup     ! Saturation mixing ratio           [  kg/kg]
   real   , dimension(mgmzp) :: rho0_cup      ! Density                           [  kg/m設
   real   , dimension(mgmzp) :: t0_cup        ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thil0_cup     ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theiv0_cup    ! Ice-vapour potential temperature  [      K]
   real   , dimension(mgmzp) :: theivs0_cup   ! Sat. Ice-vapour equiv. pot. temp. [      K]
   !----- Variables with all time tendencies but this convection. -------------------------!
   real   , dimension(mgmzp) :: co2_cup       ! CO2 mixing ratio                  [    ppm]
   real   , dimension(mgmzp) :: exner_cup     ! Exner function                    [   J/kg]
   real   , dimension(mgmzp) :: p_cup         ! Pressure                          [     Pa]
   real   , dimension(mgmzp) :: qtot_cup      ! Total mixing ratio                [  kg/kg]
   real   , dimension(mgmzp) :: qvap_cup      ! Water vapour mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qliq_cup      ! Liquid water                      [  kg/kg]
   real   , dimension(mgmzp) :: qice_cup      ! Mixing ratio                      [  kg/kg]
   real   , dimension(mgmzp) :: qsat_cup      ! Saturation mixing ratio           [  kg/kg]
   real   , dimension(mgmzp) :: rho_cup       ! Density                           [  kg/m設
   real   , dimension(mgmzp) :: t_cup         ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thil_cup      ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theiv_cup     ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: theivs_cup    ! Sat. Ice-vapour equiv. pot. temp. [      K]
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Downdraft-related variables                                                        !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real   , dimension(mgmzp) :: co20d_cld     ! CO2 water mixing ratio            [    ppm]
   real   , dimension(mgmzp) :: dby0d         ! Buoyancy acceleration             [   m/s淫
   real   , dimension(mgmzp) :: pw0d_cld      ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtot0d_cld    ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvap0d_cld    ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qliq0d_cld    ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qice0d_cld    ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: qsat0d_cld    ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: rho0d_cld     ! Density                           [  kg/m設
   real   , dimension(mgmzp) :: t0d_cld       ! Temperature                       [      K]
   real   , dimension(mgmzp) :: theiv0d_cld   ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: thil0d_cld    ! Ice-liquid pot. temperature       [      K]
   real                      :: aa0d          ! Updraft work function             [   J/kg]
   real                      :: pwev0         ! Integrated evaporation            [  kg/kg]
   !----- Based on values with forcing ----------------------------------------------------!
   real   , dimension(mgmzp) :: cdd           ! Norm. dndraft detrainment rate    [    ---]
   real   , dimension(mgmzp) :: co2d_cld      ! CO2 water mixing ratio            [    ppm]
   real   , dimension(mgmzp) :: dbyd          ! Buoyancy acceleration             [   m/s淫
   real   , dimension(mgmzp) :: etad_cld      ! normalised downdraft mass flux    [    ---]
   real   , dimension(mgmzp) :: mentrd_rate   ! Norm. dndraft entrainment rate    [    ---]
   real   , dimension(mgmzp) :: pwd_cld       ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtotd_cld     ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvapd_cld     ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qliqd_cld     ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qiced_cld     ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: qsatd_cld     ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: rhod_cld      ! Density                           [  kg/m設
   real   , dimension(mgmzp) :: td_cld        ! Temperature                       [      K]
   real   , dimension(mgmzp) :: theivd_cld    ! Ice-vapour equiv. pot. temp.      [      K]
   real   , dimension(mgmzp) :: thild_cld     ! Ice-liquid pot. temperature       [      K]
   real                      :: pwav          ! Integrated condensation           [  kg/kg]
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Updraft-related variables                                                          !
   !---------------------------------------------------------------------------------------!
   !----- Based on current values (no forcing) --------------------------------------------!
   real   , dimension(mgmzp) :: co20u_cld     ! CO2 water mixing ratio            [    ppm]
   real   , dimension(mgmzp) :: dby0u         ! Buoyancy acceleration             [   m/s淫
   real   , dimension(mgmzp) :: pw0u_cld      ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtot0u_cld    ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvap0u_cld    ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qliq0u_cld    ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qice0u_cld    ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: qsat0u_cld    ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: rho0u_cld     ! Density                           [  kg/m設
   real   , dimension(mgmzp) :: t0u_cld       ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thil0u_cld    ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theiv0u_cld   ! Ice-vapour equiv. pot. temp.      [      K]
   real                      :: aa0u          ! Updraft work function             [   J/kg]
   real                      :: pwav0         ! Integrated condensation           [  kg/kg]
   real                      :: wbuoymin0     ! Minimum buoyant velocity          [    m/s]
   !----- Based on values with forcing ----------------------------------------------------!
   real   , dimension(mgmzp) :: cdu           ! Norm. updraft detrainment rate    [    ---]
   real   , dimension(mgmzp) :: co2u_cld      ! CO2 water mixing ratio            [    ppm]
   real   , dimension(mgmzp) :: dbyu          ! Buoyancy acceleration             [   m/s淫
   real   , dimension(mgmzp) :: etau_cld      ! Normalised updraft mass flux      [    ---]
   real   , dimension(mgmzp) :: mentru_rate   ! Norm. updraft entrainment rate    [    ---]
   real   , dimension(mgmzp) :: pwu_cld       ! Condensation                      [  kg/kg]
   real   , dimension(mgmzp) :: qtotu_cld     ! Total water mixing ratio          [  kg/kg]
   real   , dimension(mgmzp) :: qvapu_cld     ! Vapour mixing ratio               [  kg/kg]
   real   , dimension(mgmzp) :: qliqu_cld     ! Liquid water mixing ratio         [  kg/kg]
   real   , dimension(mgmzp) :: qiceu_cld     ! Ice mixing ratio                  [  kg/kg]
   real   , dimension(mgmzp) :: qsatu_cld     ! Sat. mixing ratio                 [  kg/kg]
   real   , dimension(mgmzp) :: rhou_cld      ! Density                           [  kg/m設
   real   , dimension(mgmzp) :: tu_cld        ! Temperature                       [      K]
   real   , dimension(mgmzp) :: thilu_cld     ! Ice-liquid potential temperature  [      K]
   real   , dimension(mgmzp) :: theivu_cld    ! Ice-vapour equiv. pot. temp.      [      K]
   real                      :: pwev          ! Integrated evaporation            [  kg/kg]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Flag for convection error when using the "0"  variables.                           !
   !---------------------------------------------------------------------------------------!
   integer :: ierr0     ! Flag for convection error
   integer :: ktop0     ! ktop, but it is just a dummy variable.

   !---------------------------------------------------------------------------------------!
   !    Aliases for some levels/flags, they will be stored in the ensemble structure       !
   ! later.                                                                                !
   !---------------------------------------------------------------------------------------!
   !----- Scalars. ------------------------------------------------------------------------!
   integer :: ierr      ! Flag for convection error
   integer :: jmin      ! Level in which downdrafts originate
   integer :: k22       ! Level in which updrafts originate
   integer :: kbcon     ! Level of free convection
   integer :: kdet      ! Top of downdraft detrainemnt layer
   integer :: kstabi    ! cloud stable layer base
   integer :: kstabm    ! cloud stable layer top
   integer :: ktop      ! cloud top
   real    :: wbuoymin  ! Minimum buoyancy velocity
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Miscellaneous parameters                                                           !
   !---------------------------------------------------------------------------------------!
   real :: cap_max    ! Actual "Depth" of inversion capping                       [     Pa]
   real :: wbuoy_max  ! Maximum acceptable buoyant velocity                       [    m/s]
   real :: zktop      ! Actual maximum height that downdrafts can originate.      [      m]
   real :: pmass_kept ! 1.-pmass_left                                             [    ---]
   real :: edt        ! dnmf/upmf                                                 [    ---]
   !----- Scratch array -------------------------------------------------------------------!
   real, dimension(mgmzp) :: scrvar    ! Scratch variable
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !  Setting up the mass that is kept by downdraft                                        !
   !---------------------------------------------------------------------------------------!
   pmass_kept = 1.-pmass_left

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

      do k=1,mkx
         ee%cdd_cap(k,icap)         = 0.
         ee%cdu_cap(k,icap)         = 0.
         ee%mentrd_rate_cap(k,icap) = 0.
         ee%mentru_rate_cap(k,icap) = 0.
         ee%dbyd_cap(k,icap)        = 0.
         ee%dbyu_cap(k,icap)        = 0.
         ee%etad_cld_cap(k,icap)    = 0.
         ee%etau_cld_cap(k,icap)    = 0.
         ee%rhod_cld_cap(k,icap)    = 0.
         ee%rhou_cld_cap(k,icap)    = 0.
         ee%qliqd_cld_cap(k,icap)   = 0.
         ee%qliqu_cld_cap(k,icap)   = 0.
         ee%qiced_cld_cap(k,icap)   = 0.
         ee%qiceu_cld_cap(k,icap)   = 0.
      end do

      capoffset= (-1)**icap * icap / 2 !----- This gives 0,1,-1,2,-2 for icap=1,2,3... ----!
      if (cap_maxs > 0.) then
         !---------------------------------------------------------------------------------!
         ! A1. Define cap_max for this member if the user opted by cap_maxs (i.e.,         !
         !     cap_max > 0). The reference given by the user will be used as               !
         !     reference, but cap_max will be perturbed for members other than the first.  !
         !---------------------------------------------------------------------------------!
         cap_max =  max(cap_max_increment,cap_maxs + capoffset * cap_max_increment)
         wbuoy_max = 0.
      else
         !---------------------------------------------------------------------------------!
         ! A2. If the user opted by the probability, then I use capoffset to perturb the   !
         !     threshold in which I should accept convection. threwnorm is the maximum     !
         !     normalised value of vertical wind that will give a cumulative probability   !
         !     density function (CPDF) greater than what the user asked, assuming that the !
         !     vertical velocity with average equals to wp and standard deviation equals   !
         !     to sigw (from Mellor-Yamada closures).                                      ! 
         !---------------------------------------------------------------------------------!
         cap_max = cap_maxs
         wbuoy_max = wnorm_max + capoffset * wnorm_increment
      end if


      !------------------------------------------------------------------------------------!
      ! B. Initialise ierr to zero, thus allowing convection to happen.  But also setting  !
      !    wbuoymin to an unrealistically large value, so it will not make sense in case   !
      !    convection fails.                                                               !
      !------------------------------------------------------------------------------------!
      ierr     = 0
      wbuoymin = 1.e20

      !------------------------------------------------------------------------------------!
      ! C. Initialize Entrainment and Detrainment variables.                               !
      !------------------------------------------------------------------------------------!
      mentrd_rate = .2/radius             ! Entrain. rate  associated w/ dndrafts
      mentru_rate = mentrd_rate           ! Entrain. rate  associated w/ updrafts
      cdd         = 0.                    ! Detrain. fctn. associated w/ dndrafts
      cdu         = 0.1*mentru_rate       ! Detrain. fctn. associated w/ updrafts

      !------------------------------------------------------------------------------------!
      ! D. Initialize mass fluxes and buoyancy, since they may never be computed in case   !
      !    convection fails.                                                               !
      !------------------------------------------------------------------------------------!
      etad_cld = 0.
      etau_cld = 0.
      dbyd     = 0.
      dbyu     = 0.

      !------------------------------------------------------------------------------------!
      ! E. Initialize Cloud work functions associated with updrafts and downdrafts.        !
      !------------------------------------------------------------------------------------!
      aau       = 0.
      aad       = 0.
      edt       = 0.
   
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
      call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner,thil,t,qtot,qliq,qice,co2,exnersur    &
                              ,thilsur,tsur,qtotsur,qliqsur,qicesur,co2sur,exner_cup,p_cup &
                              ,t_cup,thil_cup,qtot_cup,qvap_cup,qliq_cup,qice_cup,qsat_cup &
                              ,co2_cup,rho_cup,theiv_cup,theivs_cup)

      !------------------------------------------------------------------------------------!
      ! H. Initialise drafts liquid mixing ratio and density.  These will be the           !
      !    convective values in case fails.                                                !
      !------------------------------------------------------------------------------------!
      qliqd_cld = qliq_cup
      qliqu_cld = qliq_cup
      qiced_cld = qice_cup
      qiceu_cld = qice_cup
      rhod_cld  = rho_cup
      rhou_cld  = rho_cup


      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! I. Finding some updraft properties, namely the levels in which updrafts originate, !
      !    the level of free convection and the cloud top, and the updraft-related energy, !
      !    mass, and moisture properties.                                                  !
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
      call grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,wwind,sigw,tke,qice,qliq  &
                               ,theiv_cup,ierr,k22)


      if (ierr /= 0) then
         ee%ierr_cap(icap) = ierr
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      ! 3. Finding the level of free convection. Two important points here:                !
      !    a. This call may end up preventing convection, so I must check after the call   !
      !    b. This subroutine may also affect the updraft originating level.               !
      !------------------------------------------------------------------------------------!
      call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,wbuoy_max,wwind,sigw,exner_cup     &
                               ,p_cup,theiv_cup,thil_cup,t_cup,qtot_cup,qvap_cup,qliq_cup  &
                               ,qice_cup,qsat_cup,co2_cup,rho_cup,dzd_cld,mentru_rate      &
                               ,theivu_cld,thilu_cld,tu_cld,qtotu_cld,qvapu_cld,qliqu_cld  &
                               ,qiceu_cld,qsatu_cld,co2u_cld,rhou_cld,dbyu,k22,ierr,kbcon  &
                               ,wbuoymin)

      if (ierr /= 0) then
         ee%ierr_cap(icap) = ierr
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      ! 4. Finding the minimum saturation thetae_iv. This will be the bottom of the stable !
      !    layer.                                                                          !
      !------------------------------------------------------------------------------------!
      kstabi=(kbcon - 1) + minloc(theivs_cup(kbcon:kstabm),dim=1)

      !------------------------------------------------------------------------------------!
      ! 5. Increasing the detrainment in stable layers provided that there is such layer.  !
      !    this rate increases linearly until a maximum value, currently set to 10 times   !
      !    the entrainment rate.                                                           !
      !------------------------------------------------------------------------------------!
      do k=kstabi,kstabm-1
         cdu(k) = min(10.*mentru_rate(k) , cdu(k-1) + 1.5 * mentru_rate(k))
      end do

      !------------------------------------------------------------------------------------!
      ! 6. Calculate the incloud ice-vapour equivalent potential temperature               !
      !------------------------------------------------------------------------------------!
      call grell_theiv_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,theiv,theiv_cup,dzu_cld &
                              ,theivu_cld)

      !------------------------------------------------------------------------------------!
      ! 7. Finding the normalized mass fluxes associated with updrafts. Since we are using !
      !    height-based vertical coordinate, there is no need to find the forced           !
      !    normalized mass flux, they'll be the same, so just copy it afterwards.          !
      !------------------------------------------------------------------------------------!
      call grell_nms_updraft(mkx,mgmzp,k22,kbcon,ktpse,mentru_rate,cdu,dzu_cld,etau_cld)

      !------------------------------------------------------------------------------------!
      ! 8. Finding the moisture profiles associated with updrafts.                         !
      !------------------------------------------------------------------------------------!
      call grell_most_thermo_updraft(comp_down,.true.,mkx,mgmzp,kbcon,ktpse,cdu            &
                                    ,mentru_rate,qtot,co2,p_cup,exner_cup,theiv_cup        &
                                    ,thil_cup,t_cup,qtot_cup,qvap_cup,qliq_cup,qice_cup    &
                                    ,qsat_cup,co2_cup,rho_cup,theivu_cld,etau_cld,dzu_cld  &
                                    ,thilu_cld,tu_cld,qtotu_cld,qvapu_cld,qliqu_cld        &
                                    ,qiceu_cld,qsatu_cld,co2u_cld,rhou_cld,dbyu,pwu_cld    &
                                    ,pwav,ktop,ierr)

      !------------------------------------------------------------------------------------!
      ! 9. Checking whether we found a cloud top. Since this may keep convection to        !
      !    happen, I check whether I should move on or break here.  Also check whether     !
      !    this cloud qualifies to be in this spectrum size. It need to be thicker than    !
      !    the minimum value provided by the user.                                         !
      !------------------------------------------------------------------------------------!
      if (ierr /= 0) then
         ee%ierr_cap(icap) = ierr
         cycle stacloop
      elseif (z_cup(ktop)-z_cup(k22) < depth_min) then
         ee%ierr_cap(icap) = 6
         cycle stacloop
      elseif (z_cup(ktop)-z_cup(k22) > depth_max) then
         ee%ierr_cap(icap) = 7
         cycle stacloop
      end if

      !------------------------------------------------------------------------------------!
      !10. Finding the cloud work function associated with updrafts. If this cloud doesn't !
      !    produce cloud work, break the run, we don't simulate lazy clouds in this model. !
      !------------------------------------------------------------------------------------!
      call grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,dbyu,dzu_cld,etau_cld,aau)
      if (aau == 0.) then
         ee%ierr_cap(icap) = 10
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
         kdetloop: do kdet = 1,mkx
            if (z_cup(kdet) > z_detr) exit kdetloop
         end do kdetloop
         !---------------------------------------------------------------------------------!
         ! 2. The level in which downdrafts will originate in this cloud. Since this may   !
         !    lead to an impossible cloud, check for success.                              !
         !---------------------------------------------------------------------------------!
         call grell_find_downdraft_origin(mkx,mgmzp,k22,ktop,relheight_down,zcutdown,z_cup &
                                         ,theivs_cup,dzd_cld,ierr,kdet,jmin)
         if (ierr /= 0) then
            ee%ierr_cap(icap) = ierr
            cycle stacloop
         end if

         !---------------------------------------------------------------------------------!
         ! 3. Now that we know kdet, change detrainment rate below this level.             !
         !---------------------------------------------------------------------------------!
         do k=kdet,1,-1
            cdd(k) = mentrd_rate(k)                                                        &
                   + (1. - (pmass_kept*z_cup(k) + pmass_left*z_cup(kdet))                  &
                     / (pmass_kept*z_cup(k+1) + pmass_left*z_cup(kdet)) )/dzd_cld(k)
         end do


         !---------------------------------------------------------------------------------!
         ! 4. Normalized mass fluxes.                                                      !
         !---------------------------------------------------------------------------------!
         call grell_nms_downdraft(mkx,mgmzp,kdet,jmin,mentrd_rate,cdd,z_cup,dzd_cld        &
                                 ,etad_cld)

         !---------------------------------------------------------------------------------!
         ! 5. Moist static energy. This may prevent convection to happen, check for        !
         !    success.                                                                     !
         !---------------------------------------------------------------------------------!
         call grell_theiv_downdraft(mkx,mgmzp,jmin,cdd,mentrd_rate,theiv,theiv_cup         &
                                   ,theivs_cup,dzd_cld,ierr,theivd_cld)
         if (ierr /= 0) then
            ee%ierr_cap(icap) = ierr
            cycle stacloop
         end if

         !---------------------------------------------------------------------------------!
         ! 6. Moisture properties of downdraft                                             !
         !---------------------------------------------------------------------------------!
         call grell_most_thermo_downdraft(mkx,mgmzp,jmin,qtot,co2,mentrd_rate,cdd,p_cup    &
                                         ,exner_cup,thil_cup,t_cup,qtot_cup,qvap_cup       &
                                         ,qliq_cup,qice_cup,qsat_cup,co2_cup,rho_cup,pwav  &
                                         ,theivd_cld,etad_cld,dzd_cld,thild_cld,td_cld     &
                                         ,qtotd_cld,qvapd_cld,qliqd_cld,qiced_cld          &
                                         ,qsatd_cld,co2d_cld,rhod_cld,dbyd,pwd_cld,pwev    &
                                         ,ierr)
         if (ierr /= 0) then
            ee%ierr_cap(icap) = ierr
            cycle stacloop
         end if
         !---------------------------------------------------------------------------------!
         ! 7. Compute the downdraft strength in terms of windshear. Remembering that edt   !
         !    is the fraction between downdraft and updraft.                               !
         !---------------------------------------------------------------------------------!
         call grell_efficiency_ensemble(mkx,mgmzp,maxens_eff,k22,kbcon,ktop,edtmin,edtmax  &
                                       ,pwav,pwev,z_cup,uwind,vwind,dzd_cld                &
                                       ,ee%edt_eff(:,icap))

         !---------------------------------------------------------------------------------!
         ! 8. Checking for water availability and evaporation consistency: we assume that  !
         !    downdraft is always saturated, and it gets the moisture from the rain. How-  !
         !    ever, it cannot require more rain than what is available, so if that would   !
         !    be the case, we don't allow this cloud.                                      !
         !---------------------------------------------------------------------------------!
         do iedt = 1, maxens_eff
            if (pwav + ee%edt_eff(iedt,icap) * pwev < 0. ) then
               !ee%edt_eff(iedt,icap) = - 0.9999 * pwav / pwev
               ee%ierr_cap(icap) = 9
               cycle stacloop
            end if
         end do

         !---------------------------------------------------------------------------------!
         ! 9. Computing cloud work function associated with downdrafts                     !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_downdraft(mkx,mgmzp,jmin,dbyd,dzd_cld,etad_cld,aad)

      else
         !---------------------------------------------------------------------------------!
         !    Now we assume that if a cloud doesn't have dowdrafts, it doesn't precipitate !
         ! either. In the future we may revisit this assumption. For now, if the cloud is  !
         ! shallow, then we set up some of the above parameters for non existent dowdraft  !
         ! (i.e. no mass flux, downdraft thermodynamic variables coincident to cloud       !
         ! levels etc.). kdet and jmin were already initialized to zero.                   !
         !---------------------------------------------------------------------------------!
         etad_cld           = 0.
         dbyd               = 0.
         thild_cld          = thil_cup
         theivd_cld         = theiv_cup
         qtotd_cld          = qtot_cup
         qliqd_cld          = qliq_cup
         qiced_cld          = qice_cup
         co2d_cld           = co2_cup
         pwd_cld            = 0.
         pwev               = 0.
         aad                = 0.
         do iedt=1,maxens_eff
            ee%edt_eff(iedt,icap) = 0.
         end do
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
         ! i.     Initialising error flag                                                  !
         !---------------------------------------------------------------------------------!
         ierr0   = 0
      
         !---------------------------------------------------------------------------------!
         ! ii.    Initialise detrainment and cloud work variables.                         !
         !---------------------------------------------------------------------------------!
         aa0d      = 0.
         aa0u      = 0.

         !---------------------------------------------------------------------------------!
         ! iii.   Calculate all thermodynamic properties at the cloud level and initialize !
         !        draft thermodynamic properties                                           !
         !---------------------------------------------------------------------------------!
         call grell_thermo_cldlev(mkx,mgmzp,z_cup,exner0,thil0,t0,qtot0,qliq0,qice0,co20   &
                                 ,exnersur,thilsur,tsur,qtotsur,qliqsur,qicesur,co2sur     &
                                 ,exner0_cup,p0_cup,t0_cup,thil0_cup,qtot0_cup,qvap0_cup   &
                                 ,qliq0_cup,qice0_cup,qsat0_cup,co20_cup,rho0_cup          &
                                 ,theiv0_cup,theivs0_cup)

         !---------------------------------------------------------------------------------!
         ! iv.    Calculate the thermodynamic properties below the level of free convec-   !
         !        tion.                                                                    !
         !---------------------------------------------------------------------------------!
         call grell_buoy_below_lfc(mkx,mgmzp,k22,kbcon,exner0_cup,p0_cup,theiv0_cup        &
                                  ,thil0_cup,t0_cup,qtot0_cup,qvap0_cup,qliq0_cup          &
                                  ,qice0_cup,qsat0_cup,co20_cup,rho0_cup,theiv0u_cld       &
                                  ,thil0u_cld,t0u_cld,qtot0u_cld,qvap0u_cld,qliq0u_cld     &
                                  ,qice0u_cld,qsat0u_cld,co20u_cld,rho0u_cld,dby0u)

         !---------------------------------------------------------------------------------!
         ! v.     Calculate the incloud ice-vapour equivalent potential temperature        !
         !---------------------------------------------------------------------------------!
         call grell_theiv_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,theiv0,theiv0_cup    &
                                 ,dzu_cld,theiv0u_cld)

         !---------------------------------------------------------------------------------!
         ! vi.    Finding the moisture profiles associated with updrafts.                  !
         !---------------------------------------------------------------------------------!
         call grell_most_thermo_updraft(comp_down,.false.,mkx,mgmzp,kbcon,ktop,cdu         &
                                       ,mentru_rate,qtot0,co20,p0_cup,exner0_cup           &
                                       ,theiv0_cup,thil0_cup,t0_cup,qtot0_cup,qvap0_cup    &
                                       ,qliq0_cup,qice0_cup,qsat0_cup,co20_cup,rho0_cup    &
                                       ,theiv0u_cld,etau_cld,dzu_cld,thil0u_cld,t0u_cld    &
                                       ,qtot0u_cld,qvap0u_cld,qliq0u_cld,qice0u_cld        &
                                       ,qsat0u_cld,co2u_cld,rho0u_cld,dby0u,pw0u_cld,pwav0 &
                                       ,ktop0,ierr0)
         !---------------------------------------------------------------------------------!
         ! vii.   Finding the cloud work function                                          !
         !---------------------------------------------------------------------------------!
         call grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,dby0u,dzu_cld,etau_cld,aa0u)
         
         !----- Downdraft properties ------------------------------------------------------!
         if (comp_down) then

            !------------------------------------------------------------------------------!
            ! viii.  Moist static energy.                                                  !
            !------------------------------------------------------------------------------!
            call grell_theiv_downdraft(mkx,mgmzp,jmin,cdd,mentrd_rate,theiv0,theiv0_cup    &
                                      ,theivs0_cup,dzd_cld,ierr0,theiv0d_cld)
            !------------------------------------------------------------------------------!
            ! ix.    Moisture properties of downdraft                                      !
            !------------------------------------------------------------------------------!
            call grell_most_thermo_downdraft(mkx,mgmzp,jmin,qtot0,co20,mentrd_rate,cdd     &
                                            ,p0_cup,exner0_cup,thil0_cup,t0_cup,qtot0_cup  &
                                            ,qvap0_cup,qliq0_cup,qice0_cup,qsat0_cup       &
                                            ,co20_cup,rho0_cup,pwav0,theiv0d_cld,etad_cld  &
                                            ,dzd_cld,thil0d_cld,t0d_cld,qtot0d_cld         &
                                            ,qvap0d_cld,qliq0d_cld,qice0d_cld,qsat0d_cld   &
                                            ,co2d_cld,rho0d_cld,dby0d,pw0d_cld,pwev0,ierr0)

            !------------------------------------------------------------------------------!
            ! x.     Downdraft cloud work function.                                        !
            !------------------------------------------------------------------------------!
            call grell_cldwork_downdraft(mkx,mgmzp,jmin,dby0d,dzd_cld,etad_cld,aa0d)

         end if

      end if
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      ! L. Now we perform the loop across the precipitation efficiency ensembles to find   !
      !    the normalised tendency of conserved variables due to convection for each       !
      !    ensemble member.                                                                !
      !------------------------------------------------------------------------------------!
      effloop: do iedt=1,maxens_eff
         
         !---------------------------------------------------------------------------------!
         ! 1. Copying the current epsilon to a shortcut. This is an output variable, how-  !
         !    ever, it will be used as scratch here. At the feedback time the value will   !
         !    be overwritten by the actual output value.                                   !
         !---------------------------------------------------------------------------------!
         edt  = ee%edt_eff(iedt,icap)

         !---------------------------------------------------------------------------------!
         ! 2. Compute the total cloud work function for this edt member                    !
         !---------------------------------------------------------------------------------!
         if (comp_noforc_cldwork) ee%aatot0_eff(iedt,icap) = aa0u + edt * aa0d
         ee%aatot_eff(iedt,icap)  = aau + edt * aad


         !---------------------------------------------------------------------------------!
         ! 3. Compute the changes in thetae_iv,theta_il, and total mixing ratio at the     !
         !    bottom, considering downdraft effects. If this cloud doesn't have down-      !
         !    drafts, it should be zero so skip it.                                        !
         !---------------------------------------------------------------------------------!
         if (comp_down) then
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,theiv,p_cup          &
                                         ,theiv_cup,mentrd_rate,cdd,dzd_cld,etad_cld       &
                                         ,theivd_cld,ee%dellatheiv_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,qtot,p_cup,qtot_cup  &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,qtotd_cld       &
                                         ,ee%dellaqtot_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,thil,p_cup,thil_cup  &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,thild_cld       &
                                         ,ee%dellathil_eff(1,iedt,icap))
            call  grell_dellabot_ensemble(mgmzp,checkmass,masstol,edt,co2,p_cup,co2_cup    &
                                         ,mentrd_rate,cdd,dzd_cld,etad_cld,co2d_cld        &
                                         ,ee%dellaco2_eff(1,iedt,icap))
         end if
         
         !---------------------------------------------------------------------------------!
         ! 4. Compute the changes in thetae_iv and total mixing ratio at other levels. In  !
         !    case downdrafts are absent, this will work fine too, since the mass fluxes   !
         !    will be all zero.                                                            !
         !---------------------------------------------------------------------------------!
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,theiv,p_cup,theiv_cup,mentrd_rate,mentru_rate,cdd      &
                                   ,cdu,dzd_cld,etad_cld,etau_cld,theivd_cld,theivu_cld    &
                                   ,ee%dellatheiv_eff(:,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,qtot,p_cup,qtot_cup,mentrd_rate,mentru_rate,cdd,cdu    &
                                   ,dzd_cld,etad_cld,etau_cld,qtotd_cld,qtotu_cld          &
                                   ,ee%dellaqtot_eff(:,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,thil,p_cup,thil_cup,mentrd_rate,mentru_rate,cdd,cdu    &
                                   ,dzd_cld,etad_cld,etau_cld,thild_cld,thilu_cld          &
                                   ,ee%dellathil_eff(:,iedt,icap))
         call grell_dellas_ensemble(mgmzp,checkmass,masstol,edt,kdet,k22,kbcon,jmin,ktop   &
                                   ,co2,p_cup,co2_cup,mentrd_rate,mentru_rate,cdd,cdu      &
                                   ,dzd_cld,etad_cld,etau_cld,co2d_cld,co2u_cld            &
                                   ,ee%dellaco2_eff(:,iedt,icap))


         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write(unit=16,fmt='(a)') '-------------------------------------------------------'
         !write(unit=16,fmt='(2(a,1x,i5,1x))') '   Densities: iedt=',iedt,'icap=',icap
         !write(unit=16,fmt='(4(a,1x,f10.4,1x))') '     k22      =',   z(ee%k22)               &
         !                                            ,'kbcon    =',   z(ee%kbcon)             &
         !                                            ,'jmin     =',   z(ee%jmin)              &
         !                                            ,'ktop     =',   z(ee%ktop)
         !write(unit=16,fmt='(a5,2(3x,5(1x,a11)))')                                         &
         !   'Level','   RHO0_CUP','  RHO0D_CLD','  RHO0U_CLD','      DBY0D','      DBY0U'  &
         !          ,'    RHO_CUP','   RHOD_CLD','   RHOU_CLD','       DBYD','       DBYU'
         !do k=1,mkx
         !   write(unit=16,fmt='(i5,2(3x,5(1x,f11.5)))')                                    &
         !       k,rho0_cup(k),rho0d_cld(k),rho0u_cld(k),dby0d(k),dby0u(k)                  &
         !        ,rho_cup (k),ee%rhod_cld (k),ee%rhou_cld (k),ee%dbyd (k),ee%dbyu (k)
         !end do
         !write(unit=16,fmt='(a)') '-------------------------------------------------------'
         !write(unit=16,fmt='(a)') ' '
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !---------------------------------------------------------------------------------!
         ! 5. Computing the rainfall flux, which will be simply all the fall-out water     !
         !    that wasn't consumed by the downdrafts.                                      !
         !---------------------------------------------------------------------------------!
         do k=1,ktop
            !------ Precipitation ---------------------------------------------------------!
            ee%pw_eff(k,iedt,icap) =  pwu_cld(k) + edt * pwd_cld(k)
         end do
      end do effloop
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!

      !------------------------------------------------------------------------------------!
      ! M.  Before we cycle to the next static control, we copy the variables that depend  !
      !     on the static control.  In case this static control didn't happen and we don't !
      !     reach this point, no problem, the values will keep the original value (zero).  !
      !------------------------------------------------------------------------------------!
      ee%ierr_cap(icap)          = 0
      ee%jmin_cap(icap)          = jmin
      ee%k22_cap(icap)           = k22
      ee%kbcon_cap(icap)         = kbcon
      ee%kdet_cap(icap)          = kdet
      ee%kstabi_cap(icap)        = kstabi
      ee%kstabm_cap(icap)        = kstabm
      ee%ktop_cap(icap)          = ktop

      ee%pwav_cap(icap)          = pwav
      ee%pwev_cap(icap)          = pwev
      ee%wbuoymin_cap(icap)      = wbuoymin
      
      do k=1,mkx
         ee%cdd_cap(k,icap)         = cdd(k)
         ee%cdu_cap(k,icap)         = cdu(k)
         ee%mentrd_rate_cap(k,icap) = mentrd_rate(k)
         ee%mentru_rate_cap(k,icap) = mentru_rate(k)
         ee%dbyd_cap(k,icap)        = dbyd(k)
         ee%dbyu_cap(k,icap)        = dbyu(k)
         ee%etad_cld_cap(k,icap)    = etad_cld(k)
         ee%etau_cld_cap(k,icap)    = etau_cld(k)
         ee%rhod_cld_cap(k,icap)    = rhod_cld(k)
         ee%rhou_cld_cap(k,icap)    = rhou_cld(k)
         ee%qliqd_cld_cap(k,icap)   = qliqd_cld(k)
         ee%qliqu_cld_cap(k,icap)   = qliqu_cld(k)
         ee%qiced_cld_cap(k,icap)   = qiced_cld(k)
         ee%qiceu_cld_cap(k,icap)   = qiceu_cld(k)
      end do

   end do stacloop

   return
end subroutine grell_cupar_static
!==========================================================================================!
!==========================================================================================!
