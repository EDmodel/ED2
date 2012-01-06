!==========================================================================================!
!==========================================================================================!
!     This module contains some memory variables used by the simple lake model.            !
!------------------------------------------------------------------------------------------!
module lake_coms

   use consts_coms, only : grav ! ! intent(in)

   !---------------------------------------------------------------------------------------!
   !    This variable structure has the site-level data that remain constant during the    !
   ! integration (the "meteorological" forcing).                                           !
   !---------------------------------------------------------------------------------------!
   type lakemettype
      real(kind=8)  :: atm_rhos
      real(kind=8)  :: atm_tmp
      real(kind=8)  :: atm_theta
      real(kind=8)  :: atm_theiv
      real(kind=8)  :: atm_shv
      real(kind=8)  :: atm_rvap
      real(kind=8)  :: atm_rhv
      real(kind=8)  :: atm_co2
      real(kind=8)  :: atm_exner
      real(kind=8)  :: atm_prss
      real(kind=8)  :: atm_vels
      real(kind=8)  :: ucos
      real(kind=8)  :: usin
      real(kind=8)  :: geoht
      real(kind=8)  :: dsst_dt
      real(kind=8)  :: rshort
      real(kind=8)  :: rlong
      real(kind=8)  :: tanz
      real(kind=8)  :: lon
      real(kind=8)  :: lat
   end type lakemettype
   !---------------------------------------------------------------------------------------!


   !----- This is the actual buffer structure. --------------------------------------------!
   type(lakemettype)      :: lakemet
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    This variable structure has the site-level data that is updated during the         !
   ! integration (the "canopy" air space and fluxes).                                      !
   !---------------------------------------------------------------------------------------!
   type lakesitetype
      real(kind=8)  :: can_temp
      real(kind=8)  :: can_enthalpy
      real(kind=8)  :: can_theiv
      real(kind=8)  :: can_theta
      real(kind=8)  :: can_shv
      real(kind=8)  :: can_rhv
      real(kind=8)  :: can_ssh
      real(kind=8)  :: can_co2
      real(kind=8)  :: can_prss
      real(kind=8)  :: can_rvap
      real(kind=8)  :: can_exner
      real(kind=8)  :: can_rhos
      real(kind=8)  :: can_depth
      real(kind=8)  :: lake_temp
      real(kind=8)  :: lake_fliq
      real(kind=8)  :: lake_shv
      real(kind=8)  :: lake_ssh
      real(kind=8)  :: lake_rough
      real(kind=8)  :: ustar
      real(kind=8)  :: tstar
      real(kind=8)  :: estar
      real(kind=8)  :: qstar
      real(kind=8)  :: cstar
      real(kind=8)  :: zeta
      real(kind=8)  :: ribulk
      real(kind=8)  :: gglake
      real(kind=8)  :: avg_vapor_gc
      real(kind=8)  :: avg_vapor_ac
      real(kind=8)  :: avg_sensible_gc
      real(kind=8)  :: avg_sensible_ac
      real(kind=8)  :: avg_carbon_gc
      real(kind=8)  :: avg_carbon_ac
      real(kind=8)  :: avg_carbon_st
      real(kind=8)  :: avg_sflux_u
      real(kind=8)  :: avg_sflux_w
      real(kind=8)  :: avg_sflux_v
      real(kind=8)  :: avg_sflux_t
      real(kind=8)  :: avg_sflux_r
      real(kind=8)  :: avg_sflux_c
      real(kind=8)  :: avg_rshort_gnd
      real(kind=8)  :: avg_albedt
      real(kind=8)  :: avg_rlongup
      real(kind=8)  :: avg_ustar
      real(kind=8)  :: avg_tstar
      real(kind=8)  :: avg_qstar
      real(kind=8)  :: avg_cstar
   end type lakesitetype
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Structure with all the necessary buffers.                                             !
   !---------------------------------------------------------------------------------------!
   type integration_lake
      type(lakesitetype), pointer :: initp   ! The current state
      type(lakesitetype), pointer :: yscal   ! The scale for prognostic variables
      type(lakesitetype), pointer :: y       ! 
      type(lakesitetype), pointer :: dydx    ! 
      type(lakesitetype), pointer :: yerr    ! The error for the current guess
      type(lakesitetype), pointer :: ytemp   ! Temporary
      type(lakesitetype), pointer :: ak2     ! 
      type(lakesitetype), pointer :: ak3     !
      type(lakesitetype), pointer :: ak4     ! 
      type(lakesitetype), pointer :: ak5     ! 
      type(lakesitetype), pointer :: ak6     ! 
      type(lakesitetype), pointer :: ak7     ! 
   end type integration_lake
   !---------------------------------------------------------------------------------------!


   !----- This is the actual integration buffer structure. --------------------------------!
   type(integration_lake) :: lake_buff
   !---------------------------------------------------------------------------------------!


   !----- Local constants -----------------------------------------------------------------!
   real(kind=8), parameter :: z0fac_water = 1.6d-2 / dble(grav)
   !---------------------------------------------------------------------------------------!


   !----- Constants for albedo calculation. -----------------------------------------------!
   real(kind=8), parameter :: emiss_h2o  =  9.70d-1 ! Water emissivity
   real(kind=8), parameter :: albt_inter = -1.39d-2 ! Intercept for albedo
   real(kind=8), parameter :: albt_slope =  4.67d-2 ! Slope for albedo
   real(kind=8), parameter :: albt_min   =  3.00d-2 ! Minimum albedo
   real(kind=8), parameter :: albt_max   =  9.99d-1 ! Maximum albedo
   !---------------------------------------------------------------------------------------!


   !----- "Canopy" water and heat capacity variables. -------------------------------------!
   real(kind=8)    :: wcapcan
   real(kind=8)    :: wcapcani
   real(kind=8)    :: hcapcani
   real(kind=8)    :: ccapcani
   !---------------------------------------------------------------------------------------!


   !------ Time step stuff. ---------------------------------------------------------------!
   real(kind=8)    :: tlbeg
   real(kind=8)    :: tlend
   real(kind=8)    :: dtlake
   real(kind=8)    :: dtlakei
   !---------------------------------------------------------------------------------------!

   !----- Small number, to avoid singularities. -------------------------------------------!
   real(kind=8), parameter :: tiny_lakeoff = 1.d-20
   !----- Huge number, to bypass errmax checks. -------------------------------------------!
   real(kind=8), parameter :: huge_lakeoff = 1.d30
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This variable nullifies the integration buffer for a safe initialisation.         !
   !---------------------------------------------------------------------------------------!
   subroutine initial_lake_buff()
      implicit none
      nullify(lake_buff%initp)
      nullify(lake_buff%yscal)
      nullify(lake_buff%y    )
      nullify(lake_buff%dydx )
      nullify(lake_buff%yerr )
      nullify(lake_buff%ytemp)
      nullify(lake_buff%ak2  )
      nullify(lake_buff%ak3  )
      nullify(lake_buff%ak4  )
      nullify(lake_buff%ak5  )
      nullify(lake_buff%ak6  )
      nullify(lake_buff%ak7  )

     
      allocate(lake_buff%initp)
      allocate(lake_buff%yscal)
      allocate(lake_buff%y    )
      allocate(lake_buff%dydx )
      allocate(lake_buff%yerr )
      allocate(lake_buff%ytemp)
      allocate(lake_buff%ak2  )
      allocate(lake_buff%ak3  )
      allocate(lake_buff%ak4  )
      allocate(lake_buff%ak5  )
      allocate(lake_buff%ak6  )
      allocate(lake_buff%ak7  )

      return
   end subroutine initial_lake_buff
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This variable initialises the lakesite variable.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine zero_lakesite(lake)

      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(lakesitetype), target :: lake
      !------------------------------------------------------------------------------------!


      !----- Reset the variables. ---------------------------------------------------------!
      lake%can_temp        = 0.d0
      lake%can_enthalpy    = 0.d0
      lake%can_theiv       = 0.d0
      lake%can_theta       = 0.d0
      lake%can_shv         = 0.d0
      lake%can_rhv         = 0.d0
      lake%can_ssh         = 0.d0
      lake%can_co2         = 0.d0
      lake%can_prss        = 0.d0
      lake%can_rvap        = 0.d0
      lake%can_exner       = 0.d0
      lake%can_rhos        = 0.d0
      lake%can_depth       = 0.d0
      lake%lake_temp       = 0.d0
      lake%lake_fliq       = 0.d0
      lake%lake_shv        = 0.d0
      lake%lake_ssh        = 0.d0
      lake%lake_rough      = 0.d0
      lake%ustar           = 0.d0
      lake%tstar           = 0.d0
      lake%estar           = 0.d0
      lake%qstar           = 0.d0
      lake%cstar           = 0.d0
      lake%zeta            = 0.d0
      lake%ribulk          = 0.d0
      lake%gglake          = 0.d0
      lake%avg_vapor_gc    = 0.d0
      lake%avg_vapor_ac    = 0.d0
      lake%avg_sensible_gc = 0.d0
      lake%avg_sensible_ac = 0.d0
      lake%avg_carbon_gc   = 0.d0
      lake%avg_carbon_ac   = 0.d0
      lake%avg_carbon_st   = 0.d0
      lake%avg_sflux_u     = 0.d0
      lake%avg_sflux_w     = 0.d0
      lake%avg_sflux_v     = 0.d0
      lake%avg_sflux_t     = 0.d0
      lake%avg_sflux_r     = 0.d0
      lake%avg_sflux_c     = 0.d0
      lake%avg_rshort_gnd  = 0.d0
      lake%avg_albedt      = 0.d0
      lake%avg_rlongup     = 0.d0
      lake%avg_ustar       = 0.d0
      lake%avg_tstar       = 0.d0
      lake%avg_qstar       = 0.d0
      lake%avg_cstar       = 0.d0
      !------------------------------------------------------------------------------------!

      return
   end subroutine zero_lakesite
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This variable initialises the lakesite variable.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine clone_lakesite(lakein,lakeout)

      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(lakesitetype), target :: lakein
      type(lakesitetype), target :: lakeout
      !------------------------------------------------------------------------------------!


      !----- Reset the variables. ---------------------------------------------------------!
      lakeout%can_temp        = lakein%can_temp
      lakeout%can_enthalpy    = lakein%can_enthalpy
      lakeout%can_theiv       = lakein%can_theiv
      lakeout%can_theta       = lakein%can_theta
      lakeout%can_shv         = lakein%can_shv
      lakeout%can_rhv         = lakein%can_rhv
      lakeout%can_ssh         = lakein%can_ssh
      lakeout%can_co2         = lakein%can_co2
      lakeout%can_prss        = lakein%can_prss
      lakeout%can_rvap        = lakein%can_rvap
      lakeout%can_exner       = lakein%can_exner
      lakeout%can_rhos        = lakein%can_rhos
      lakeout%can_depth       = lakein%can_depth
      lakeout%lake_temp       = lakein%lake_temp  
      lakeout%lake_fliq       = lakein%lake_fliq  
      lakeout%lake_shv        = lakein%lake_shv   
      lakeout%lake_ssh        = lakein%lake_ssh   
      lakeout%lake_rough      = lakein%lake_rough 
      lakeout%ustar           = lakein%ustar
      lakeout%tstar           = lakein%tstar
      lakeout%estar           = lakein%estar
      lakeout%qstar           = lakein%qstar
      lakeout%cstar           = lakein%cstar
      lakeout%zeta            = lakein%zeta
      lakeout%ribulk          = lakein%ribulk
      lakeout%gglake          = lakein%gglake
      lakeout%avg_vapor_gc    = lakein%avg_vapor_gc
      lakeout%avg_vapor_ac    = lakein%avg_vapor_ac
      lakeout%avg_sensible_gc = lakein%avg_sensible_gc
      lakeout%avg_sensible_ac = lakein%avg_sensible_ac
      lakeout%avg_carbon_gc   = lakein%avg_carbon_gc
      lakeout%avg_carbon_ac   = lakein%avg_carbon_ac
      lakeout%avg_carbon_st   = lakein%avg_carbon_st
      lakeout%avg_sflux_u     = lakein%avg_sflux_u
      lakeout%avg_sflux_w     = lakein%avg_sflux_w
      lakeout%avg_sflux_v     = lakein%avg_sflux_v
      lakeout%avg_sflux_t     = lakein%avg_sflux_t
      lakeout%avg_sflux_r     = lakein%avg_sflux_r
      lakeout%avg_sflux_c     = lakein%avg_sflux_c
      lakeout%avg_rshort_gnd  = lakein%avg_rshort_gnd
      lakeout%avg_albedt      = lakein%avg_albedt
      lakeout%avg_rlongup     = lakein%avg_rlongup
      lakeout%avg_ustar       = lakein%avg_ustar
      lakeout%avg_tstar       = lakein%avg_tstar
      lakeout%avg_qstar       = lakein%avg_qstar
      lakeout%avg_cstar       = lakein%avg_cstar
      !------------------------------------------------------------------------------------!

      return
   end subroutine clone_lakesite
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This variable integrates the prognostic variables.                                !
   !---------------------------------------------------------------------------------------!
   subroutine integ_lakesite(lake,dlakedt,dtim)

      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(lakesitetype), target      :: lake
      type(lakesitetype), target      :: dlakedt
      real(kind=8)      , intent(in)  :: dtim
      !------------------------------------------------------------------------------------!


      !----- Integrate the variables. -----------------------------------------------------!
      lake%can_shv         = lake%can_shv          + dtim * dlakedt%can_shv
      lake%can_co2         = lake%can_co2          + dtim * dlakedt%can_co2
      lake%can_enthalpy    = lake%can_enthalpy     + dtim * dlakedt%can_enthalpy
      lake%lake_temp       = lake%lake_temp        + dtim * dlakedt%lake_temp
      lake%avg_vapor_gc    = lake%avg_vapor_gc     + dtim * dlakedt%avg_vapor_gc
      lake%avg_vapor_ac    = lake%avg_vapor_ac     + dtim * dlakedt%avg_vapor_ac
      lake%avg_sensible_gc = lake%avg_sensible_gc  + dtim * dlakedt%avg_sensible_gc
      lake%avg_sensible_ac = lake%avg_sensible_ac  + dtim * dlakedt%avg_sensible_ac
      lake%avg_carbon_gc   = lake%avg_carbon_gc    + dtim * dlakedt%avg_carbon_gc
      lake%avg_carbon_ac   = lake%avg_carbon_ac    + dtim * dlakedt%avg_carbon_ac
      lake%avg_carbon_st   = lake%avg_carbon_st    + dtim * dlakedt%avg_carbon_st
      lake%avg_sflux_u     = lake%avg_sflux_u      + dtim * dlakedt%avg_sflux_u
      lake%avg_sflux_w     = lake%avg_sflux_w      + dtim * dlakedt%avg_sflux_w
      lake%avg_sflux_v     = lake%avg_sflux_v      + dtim * dlakedt%avg_sflux_v
      lake%avg_sflux_t     = lake%avg_sflux_t      + dtim * dlakedt%avg_sflux_t
      lake%avg_sflux_r     = lake%avg_sflux_r      + dtim * dlakedt%avg_sflux_r
      lake%avg_sflux_c     = lake%avg_sflux_c      + dtim * dlakedt%avg_sflux_c
      lake%avg_rshort_gnd  = lake%avg_rshort_gnd   + dtim * dlakedt%avg_rshort_gnd
      lake%avg_albedt      = lake%avg_albedt       + dtim * dlakedt%avg_albedt
      lake%avg_rlongup     = lake%avg_rlongup      + dtim * dlakedt%avg_rlongup
      lake%avg_ustar       = lake%avg_ustar        + dtim * dlakedt%avg_ustar
      lake%avg_tstar       = lake%avg_tstar        + dtim * dlakedt%avg_tstar
      lake%avg_qstar       = lake%avg_qstar        + dtim * dlakedt%avg_qstar
      lake%avg_cstar       = lake%avg_cstar        + dtim * dlakedt%avg_cstar
      !------------------------------------------------------------------------------------!

      return
   end subroutine integ_lakesite
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This variable normalises some of the flux variables.                              !
   !---------------------------------------------------------------------------------------!
   subroutine normal_lakesite(lake,dttot)

      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(lakesitetype), target      :: lake
      real(kind=8)      , intent(in)  :: dttot
      !------ Local variables. ------------------------------------------------------------!
      real(kind=8)                      :: dttoti
      !------------------------------------------------------------------------------------!



      !----- Find the inverse of the full time step. --------------------------------------!
      dttoti = 1.d0 / dttot
      !------------------------------------------------------------------------------------!


      !----- Normalise the variables. -----------------------------------------------------!
      lake%avg_vapor_gc        = lake%avg_vapor_gc    * dttoti
      lake%avg_vapor_ac        = lake%avg_vapor_ac    * dttoti
      lake%avg_sensible_gc     = lake%avg_sensible_gc * dttoti
      lake%avg_sensible_ac     = lake%avg_sensible_ac * dttoti
      lake%avg_carbon_gc       = lake%avg_carbon_gc   * dttoti
      lake%avg_carbon_st       = lake%avg_carbon_st   * dttoti
      lake%avg_sflux_u         = lake%avg_sflux_u     * dttoti
      lake%avg_sflux_w         = lake%avg_sflux_w     * dttoti
      lake%avg_sflux_v         = lake%avg_sflux_v     * dttoti
      lake%avg_sflux_t         = lake%avg_sflux_t     * dttoti
      lake%avg_sflux_r         = lake%avg_sflux_r     * dttoti
      lake%avg_sflux_c         = lake%avg_sflux_c     * dttoti
      lake%avg_rshort_gnd      = lake%avg_rshort_gnd  * dttoti
      lake%avg_albedt          = lake%avg_albedt      * dttoti
      lake%avg_rlongup         = lake%avg_rlongup     * dttoti
      lake%avg_ustar           = lake%avg_ustar       * dttoti
      lake%avg_tstar           = lake%avg_tstar       * dttoti
      lake%avg_qstar           = lake%avg_qstar       * dttoti
      lake%avg_cstar           = lake%avg_cstar       * dttoti
      !------------------------------------------------------------------------------------!

      return
   end subroutine normal_lakesite
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine finds the error scale for the integrated variables, which will be  !
   ! later used to define the relative error.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine lake_yscal(lake,dlakedt,htry,lakescal)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(lakesitetype), target     :: lake
      type(lakesitetype), target     :: dlakedt
      type(lakesitetype), target     :: lakescal
      real(kind=8)      , intent(in) :: htry
      !------------------------------------------------------------------------------------!

      lakescal%can_enthalpy = abs(lake%can_enthalpy) + abs(dlakedt%can_enthalpy * htry)
      lakescal%can_shv      = abs(lake%can_shv     ) + abs(dlakedt%can_shv      * htry)
      lakescal%can_co2      = abs(lake%can_co2     ) + abs(dlakedt%can_co2      * htry)
      lakescal%lake_temp    = abs(lake%lake_temp   ) + abs(dlakedt%lake_temp    * htry)


      return
   end subroutine lake_yscal
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine loops through the integrating variables, seeking for the largest   !
   ! error.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine lake_errmax(errmax,lakescal,lakeerr)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8)      , intent(out) :: errmax
      type(lakesitetype), target      :: lakeerr
      type(lakesitetype), target      :: lakescal
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                    :: err
      !------------------------------------------------------------------------------------!



      !----- Initialise the error with an optimistic number... ----------------------------!
      errmax = 0.d0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Check each variable error, and keep track of the worst guess, which will be our !
      ! worst guess in the end.  We only check prognostic variables.                       !
      !------------------------------------------------------------------------------------!
      !----- 1. Log of potential temperature. ---------------------------------------------!
      err    = abs(lakeerr%can_enthalpy / lakescal%can_enthalpy)
      errmax = max(errmax,err)
      !----- 2. Specific humidity. --------------------------------------------------------!
      err    = abs(lakeerr%can_shv / lakescal%can_shv)
      errmax = max(errmax,err)
      !----- 3. CO2 mixing ratio. ---------------------------------------------------------!
      err    = abs(lakeerr%can_co2 / lakescal%can_co2)
      errmax = max(errmax,err)
      !----- 4. Lake temperature. ---------------------------------------------------------!
      err    = abs(lakeerr%lake_temp / lakescal%lake_temp)
      errmax = max(errmax,err)
      !------------------------------------------------------------------------------------!

      return
   end subroutine lake_errmax
   !=======================================================================================!
   !=======================================================================================!
end module lake_coms
!==========================================================================================!
!==========================================================================================!
