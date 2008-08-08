!==========================================================================================!
!    Module mem_scratch_grell - This module contains all Grell's scratch variables that do !
! not depend on ensemble dimensions. This module is shared between shallow and deep        !
! convection.                                                                              !
!==========================================================================================!
module mem_scratch_grell

   !SRF- feb-05-2002 : Variables for cumulus transport scheme
   !                   adapted in july-15-2002 for 5.x version
   !
   implicit none

   !------ Scalars, mostly grid definitions -----------------------------------------------!
   integer ::  mkx         & ! Number of cloud points this grid has;
              ,lpw         & ! Lower point in the vertical for this grid;
              ,kgoff       & ! Offset between BRAMS and Grell's grids
              ,kpbl        ! ! PBL when running Nakanishi/Niino
   real    ::  tscal_kf    ! ! Kain-Fritsch (1990) time scale, for ensemble.
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence (mgmzp), grid-related stuff --------------------------------------!
   real, allocatable, dimension(:) :: &
            z                         & ! Height above surface at model levels          [m]
           ,z_cup                     & ! Height above surface at cloud levels          [m]
           ,dzd_cld                   & ! Layer thickness for downdraft calculation     [m]
           ,dzu_cld                   ! ! Layer thickness for updraft calculation       [m]
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence, shallow/deep convection interface (mgmzp) -----------------------!
   real, allocatable, dimension(:) ::  &
            thetasta                   & ! Theta, with shallow cumulus effect
           ,rvsta                      ! ! Temporary rv, with shallow cumulus effect
   !---------------------------------------------------------------------------------------!


   !------ Scalars, important levels and flags --------------------------------------------!
   integer                            :: &
            ierr                         & ! Convection failure flag (0 means no failure).
           ,jmin                         & ! Level in which downdrafts originate
           ,k22                          & ! Level in which updrafts originate
           ,kbcon                        & ! cloud base
           ,kdet                         & ! Top of downdraft detrainemnt layer
           ,kstabi                       & ! cloud stable layer base
           ,kstabm                       & ! cloud stable layer top
           ,ktop                         ! ! cloud top
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence (mgmzp) ----------------------------------------------------------!
   real, allocatable, dimension(:) ::    &
            cdd                          & ! normalized downdraft detrainment rate  [  ---]
           ,cdu                          & ! normalized updraft detrainment rate    [  ---]
           ,mentrd_rate                  & ! normalized downdraft entrainment rate  [  ---]
           ,mentru_rate                  & ! normalized updraft entrainment rate    [  ---]
           ,etad_cld                     & ! normalized downdraft mass flux         [  ---]
           ,etau_cld                     & ! normalized updraft mass flux           [  ---]
           ,qd_cld                       & ! Water vapour mixing ratio at downdraft [kg/kg]
           ,qu_cld                       & ! Water vapour mixing ratio at updraft   [kg/kg]
           ,qld_cld                      & ! Liquid water mixing ratio at downdraft [kg/kg]
           ,qlu_cld                      & ! Liquid water mixing ratio at updraft   [kg/kg]
           ,td_cld                       & ! Downdraft temperature                  [    K]
           ,tu_cld                       ! ! updraft temperature                    [    K]
   !---------------------------------------------------------------------------------------!


   !------ Scalars, surface variables -----------------------------------------------------!
   real    ::  hesur     & ! Surface moist static energy                           [  J/kg]
              ,hessur    & ! Surface saturation moist static energy                [  J/kg]
              ,psur      & ! Surface pressure                                      [    Pa]
              ,qsur      & ! Surface mixing ratio                                  [ kg/kg]
              ,qessur    & ! Surface saturation mixing ratio                       [ kg/kg]
              ,tsur      ! ! Surface temperature                                   [     K]
   !---------------------------------------------------------------------------------------!

   !------ Scalars, for dynamic control ensemble calculation ------------------------------!
   real    ::  mconv     & ! Column integrated mass flux convergence              [kg/m²/s]
              ,dens_curr & ! Downdraft mass flux as density current               [kg/m²/s]
              ,prev_dnmf ! ! Downdraft mass flux at the previous call             [kg/m²/s]
   logical ::  upconv    ! ! I should assume that convection happened upstream (true/false)
   !---------------------------------------------------------------------------------------!


   !------ 1D dependence (mgmzp), variables with all forcings but convection --------------!
     real, allocatable, dimension(:) ::          &
            dexnerdt       & ! Temporary Exner function tendency                 [J/kg/K/s]
           ,drtdt          & ! Temporary mixing ratio tendency                   [ kg/kg/s]
           ,drtdt_shal     & ! Mixing ratio tendency due to shallower clouds     [ kg/kg/s]
           ,dthetadt       & ! Temporary temperature tendency                    [     K/s]
           ,dthetadt_shal  & ! temperature tendency due to shallower clouds      [     K/s]
           ,he             & ! Moist static energy                               [    J/kg]
           ,hes            & ! Saturation moist static energy                    [    J/kg]
           ,omeg           & ! Omega - Lagrangian pressure tendency              [    Pa/s]
           ,p              & ! Pressure                                          [      Pa]
           ,q              & ! Water vapour mixing ratio                         [   kg/kg]
           ,qes            & ! Saturation water vapour mixing ratio              [   kg/kg]
           ,rho            & ! Air density                                       [   kg/m³]
           ,rcpg           & ! Cloud droplet mixing ratio                        [   kg/kg]
           ,t              & ! Temperature                                       [       K]
           ,tkeg           & ! Turbulent kinetic energy                          [   m²/s²]
           ,uwind          & ! Zonal wind speed                                  [     m/s]
           ,vwind          & ! Meridional wind speed                             [     m/s]
           ,wwind          ! ! Vertical velocity                                 [     m/s]
   !------ Variables at cumulus level, used at the output ---------------------------------!
     real, allocatable, dimension(:) ::          &
            p_cup          & ! Pressure at Grell's grid                          [      Pa]
           ,q_cup          & ! Water vapour mixing ratio at Grell's grid         [   kg/kg]
           ,t_cup          ! ! Temperature at Grell's grid                       [       K]
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence (mgmzp), variables at current time step --------------------------!
   real, allocatable, dimension(:) :: &
            he0       & ! Moist static energy                                    [    J/kg]
           ,hes0      & ! Saturation moist static energy                         [    J/kg]
           ,p0        & ! Pressure with forcing                                   [    hPa]
           ,q0        & ! Water vapour mixing ratio                              [   kg/kg]
           ,qes0      & ! Saturation water vapour mixing ratio                   [   kg/kg]
           ,t0        ! ! Temperature                                             [      K]
   !---------------------------------------------------------------------------------------!


   !------ 1D dependence (mgmzp), forcing due to convectivion -----------------------------!
   real, allocatable, dimension(:) :: &
            outq      & ! Water vapour mixing ratio feedback                     [ kg/kg/s]
           ,outqc     & ! Condensated mixing ratio feedback                      [ kg/kg/s]
           ,outt      ! ! Potential temperature feedback                         [     K/s]  

   !----- Scalar, forcing due to convection -----------------------------------------------!
   real :: precip     ! ! Precipitation rate
   contains
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine alloc_scratch_grell(mgmzp)

     implicit none
     integer, intent(in)                    :: mgmzp

     allocate (z             (mgmzp))
     allocate (z_cup         (mgmzp))
     allocate (dzd_cld       (mgmzp))
     allocate (dzu_cld       (mgmzp))

     allocate (thetasta      (mgmzp))
     allocate (rvsta         (mgmzp))

     allocate (cdd           (mgmzp))
     allocate (cdu           (mgmzp))
     allocate (mentrd_rate   (mgmzp))
     allocate (mentru_rate   (mgmzp))
     allocate (etau_cld      (mgmzp))
     allocate (etad_cld      (mgmzp))
     allocate (qd_cld        (mgmzp))
     allocate (qu_cld        (mgmzp))
     allocate (qld_cld       (mgmzp))
     allocate (qlu_cld       (mgmzp))
     allocate (td_cld        (mgmzp))
     allocate (tu_cld        (mgmzp))


     allocate (dexnerdt      (mgmzp))
     allocate (drtdt         (mgmzp))
     allocate (drtdt_shal    (mgmzp))
     allocate (dthetadt      (mgmzp))
     allocate (dthetadt_shal (mgmzp))
     allocate (he            (mgmzp))
     allocate (hes           (mgmzp))
     allocate (omeg          (mgmzp))
     allocate (p             (mgmzp))
     allocate (q             (mgmzp))
     allocate (qes           (mgmzp))
     allocate (rho           (mgmzp))
     allocate (rcpg          (mgmzp))
     allocate (t             (mgmzp))
     allocate (tkeg          (mgmzp))
     allocate (uwind         (mgmzp))
     allocate (vwind         (mgmzp))
     allocate (wwind         (mgmzp))

     allocate (p_cup         (mgmzp))
     allocate (q_cup         (mgmzp))
     allocate (t_cup         (mgmzp))

     allocate (he0           (mgmzp))
     allocate (hes0          (mgmzp))
     allocate (p0            (mgmzp))
     allocate (q0            (mgmzp))
     allocate (qes0          (mgmzp))
     allocate (t0            (mgmzp))
     
     allocate (outq          (mgmzp))
     allocate (outqc         (mgmzp))
     allocate (outt          (mgmzp))

     return
  end subroutine alloc_scratch_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
  subroutine dealloc_scratch_grell()

     implicit none
    
     if(allocated(z             )) deallocate (z             )
     if(allocated(z_cup         )) deallocate (z_cup         )
     if(allocated(dzd_cld       )) deallocate (dzd_cld       )
     if(allocated(dzu_cld       )) deallocate (dzu_cld       )
    
     if(allocated(thetasta      )) deallocate (thetasta      )
     if(allocated(rvsta         )) deallocate (rvsta         )

     if(allocated(cdd           )) deallocate (cdd           )
     if(allocated(cdu           )) deallocate (cdu           )
     if(allocated(mentrd_rate   )) deallocate (mentrd_rate   )
     if(allocated(mentru_rate   )) deallocate (mentru_rate   )
     if(allocated(etad_cld      )) deallocate (etad_cld      )
     if(allocated(etau_cld      )) deallocate (etau_cld      )
     if(allocated(qd_cld        )) deallocate (qd_cld        )
     if(allocated(qu_cld        )) deallocate (qu_cld        )
     if(allocated(qld_cld       )) deallocate (qld_cld       )
     if(allocated(qlu_cld       )) deallocate (qlu_cld       )
     if(allocated(td_cld        )) deallocate (td_cld        )
     if(allocated(tu_cld        )) deallocate (tu_cld        )

     if(allocated(dexnerdt      )) deallocate (dexnerdt      )
     if(allocated(drtdt         )) deallocate (drtdt         )
     if(allocated(drtdt_shal    )) deallocate (drtdt_shal    )
     if(allocated(dthetadt      )) deallocate (dthetadt      )
     if(allocated(dthetadt_shal )) deallocate (dthetadt_shal )
     if(allocated(he            )) deallocate (he            )
     if(allocated(hes           )) deallocate (hes           )
     if(allocated(omeg          )) deallocate (omeg          )
     if(allocated(p             )) deallocate (p             )
     if(allocated(q             )) deallocate (q             )
     if(allocated(qes           )) deallocate (qes           )
     if(allocated(rcpg          )) deallocate (rcpg          )
     if(allocated(rho           )) deallocate (rho           )
     if(allocated(t             )) deallocate (t             )
     if(allocated(tkeg          )) deallocate (tkeg          )
     if(allocated(uwind         )) deallocate (uwind         )
     if(allocated(vwind         )) deallocate (vwind         )
     if(allocated(wwind         )) deallocate (wwind         )

     if(allocated(p_cup         )) deallocate (p_cup         )
     if(allocated(q_cup         )) deallocate (q_cup         )
     if(allocated(t_cup         )) deallocate (t_cup         )

     if(allocated(he0           )) deallocate (he0           )
     if(allocated(hes0          )) deallocate (hes0          )
     if(allocated(p0            )) deallocate (p0            )
     if(allocated(q0            )) deallocate (q0            )
     if(allocated(qes0          )) deallocate (qes0          )
     if(allocated(t0            )) deallocate (t0            )

     if(allocated(outq          )) deallocate (outq          )
     if(allocated(outqc         )) deallocate (outqc         )
     if(allocated(outt          )) deallocate (outt          )

     return
  end subroutine dealloc_scratch_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
 subroutine zero_scratch_grell()

   implicit none

     if(allocated(z             )) z             = 0.
     if(allocated(z_cup         )) z_cup         = 0.
     if(allocated(dzd_cld       )) dzd_cld       = 0.
     if(allocated(dzu_cld       )) dzu_cld       = 0.
    
     if(allocated(thetasta      )) thetasta      = 0.
     if(allocated(rvsta         )) rvsta         = 0.

     if(allocated(cdd           )) cdd           = 0.
     if(allocated(cdu           )) cdu           = 0.
     if(allocated(mentrd_rate   )) mentrd_rate   = 0.
     if(allocated(mentru_rate   )) mentru_rate   = 0.
     if(allocated(etad_cld      )) etad_cld      = 0.
     if(allocated(etau_cld      )) etau_cld      = 0.
     if(allocated(qd_cld        )) qd_cld        = 0.
     if(allocated(qu_cld        )) qu_cld        = 0.
     if(allocated(qld_cld       )) qld_cld       = 0.
     if(allocated(qlu_cld       )) qlu_cld       = 0.
     if(allocated(td_cld        )) td_cld        = 0.
     if(allocated(tu_cld        )) tu_cld        = 0.

     if(allocated(dexnerdt      )) dexnerdt      = 0.
     if(allocated(drtdt         )) drtdt         = 0.
     if(allocated(drtdt_shal    )) drtdt_shal    = 0.
     if(allocated(dthetadt      )) dthetadt      = 0.
     if(allocated(dthetadt_shal )) dthetadt_shal = 0.
     if(allocated(he            )) he            = 0.
     if(allocated(hes           )) hes           = 0.
     if(allocated(omeg          )) omeg          = 0.
     if(allocated(p             )) p             = 0.
     if(allocated(q             )) q             = 0.
     if(allocated(qes           )) qes           = 0.
     if(allocated(rho           )) rho           = 0.
     if(allocated(rcpg          )) rcpg          = 0.
     if(allocated(t             )) t             = 0.
     if(allocated(tkeg          )) tkeg          = 0.
     if(allocated(uwind         )) uwind         = 0.
     if(allocated(vwind         )) vwind         = 0.
     if(allocated(wwind         )) wwind         = 0.

     if(allocated(p_cup         )) p_cup         = 0.
     if(allocated(q_cup         )) q_cup         = 0.
     if(allocated(t_cup         )) t_cup         = 0.

     if(allocated(he0           )) he0           = 0.
     if(allocated(hes0          )) hes0          = 0.
     if(allocated(p0            )) p0            = 0.
     if(allocated(q0            )) q0            = 0.
     if(allocated(qes0          )) qes0          = 0.
     if(allocated(t0            )) t0            = 0.

     if(allocated(outq          )) outq          = 0.
     if(allocated(outqc         )) outqc         = 0.
     if(allocated(outt          )) outt          = 0.
     !-------------------------------------------------------------------------------------!
     ! Flushing scalars, we don't need to check for allocation here...                     !
     !-------------------------------------------------------------------------------------!
     !----- Integer variables -------------------------------------------------------------!
     mkx               = 0
     lpw               = 0
     kgoff             = 0
     kpbl              = 0

     ierr              = 0
     jmin              = 0
     k22               = 0
     kbcon             = 0
     kdet              = 0
     kstabi            = 0
     kstabm            = 0
     ktop              = 0

     !----- Real variables ----------------------------------------------------------------!
     hesur             = 0.
     hessur            = 0.
     psur              = 0.
     qsur              = 0.
     qessur            = 0.
     tsur              = 0.

     dens_curr         = 0.
     mconv             = 0.
     prev_dnmf         = 0.

     !----- Logical variables -------------------------------------------------------------!
     upconv            = .false.
     !-------------------------------------------------------------------------------------!

     return

  end subroutine zero_scratch_grell
!==========================================================================================!
!==========================================================================================!
end module mem_scratch_grell
