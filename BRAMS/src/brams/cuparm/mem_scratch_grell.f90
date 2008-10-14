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


   !------ Scalars, important levels and flags --------------------------------------------!
   integer                            :: &
            ierr                         & ! Convection failure flag (0 means no failure).
           ,jmin                         & ! Level in which downdrafts originate
           ,k22                          & ! Level in which updrafts originate
           ,kbcon                        & ! Level of free convection
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
           ,dbyd                         & ! Buoyancy associated with downdrafts    [ m/s²]
           ,dbyu                         & ! Buoyancy associated with updrafts      [ m/s²]
           ,etad_cld                     & ! normalized downdraft mass flux         [  ---]
           ,etau_cld                     & ! normalized updraft mass flux           [  ---]
           ,rhod_cld                     & ! Downdraft density                      [kg/m³]
           ,rhou_cld                     & ! Updraft density                        [kg/m³]
           ,qliqd_cld                    & ! Liquid water mixing ratio at downdraft [kg/kg]
           ,qliqu_cld                    & ! Liquid water mixing ratio at updraft   [kg/kg]
           ,qiced_cld                    & ! Ice mixing ratio at downdraft          [kg/kg]
           ,qiceu_cld                    ! ! Ice mixing ratio at updraft            [kg/kg]
   !---------------------------------------------------------------------------------------!

   !------ Scalar, static control variable ------------------------------------------------!
   real                            ::    &
            wbuoymin                     ! ! Minimum buoyant velocity               [  m/s]

   !------ Scalars, surface variables -----------------------------------------------------!
   real    ::  exnersur  & ! Surface: Exner function                               [J/kg/K]
              ,psur      & ! Surface: pressure                                     [    Pa]
              ,qtotsur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,qvapsur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,qliqsur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,qicesur   & ! Surface: mixing ratio                                 [ kg/kg]
              ,tsur      & ! Surface: temperature                                  [     K]
              ,theivsur  & ! Surface: ice-vapour equivalent potential temperature  [     K]
              ,thilsur   ! ! Surface: ice-liquid potential temperature             [     K]
   !---------------------------------------------------------------------------------------!

   !------ Scalars, for dynamic control ensemble calculation ------------------------------!
   real    ::  mconv     & ! Column integrated mass flux convergence              [kg/m²/s]
              ,dens_curr & ! Downdraft mass flux as density current               [kg/m²/s]
              ,prev_dnmf ! ! Downdraft mass flux at the previous call             [kg/m²/s]
   logical ::  upconv    ! ! I should assume that convection happened upstream (true/false)
   !---------------------------------------------------------------------------------------!


   !------ 1D dependence (mgmzp), variables with all forcings but convection --------------!
     real, allocatable, dimension(:) ::          &
            dqtotdt        & ! Temporary total mixing ratio tendency             [ kg/kg/s]
           ,dqtotdt_shal   & ! Mixing ratio tendency due to shallower clouds     [ kg/kg/s]
           ,dthildt        & ! Temporary temperature tendency                    [     K/s]
           ,dthildt_shal   & ! temperature tendency due to shallower clouds      [     K/s]
           ,dtkedt         & ! Temporary TKE tendency                            [  J/kg/s]
           ,exner          & ! Exner function                                    [  J/kg/K]
           ,omeg           & ! Omega - Lagrangian pressure tendency              [    Pa/s]
           ,p              & ! Pressure                                          [      Pa]
           ,rho            & ! Air density                                       [   kg/m³]
           ,qtot           & ! Total mixing ratio                                [   kg/kg]
           ,qice           & ! Ice mixing ratio                                  [   kg/kg]
           ,qliq           & ! Liquid water mixing ratio                         [   kg/kg]
           ,qvap           & ! Water vapour mixing ratio                         [   kg/kg]
           ,t              & ! Temperature                                       [       K]
           ,thil           & ! Ice-liquid potential temperature                  [       K]
           ,theiv          & ! Ice-vapour equivalent potential temperature       [       K]
           ,tke            & ! Turbulent kinetic energy                          [    J/kg]
           ,sigw           & ! Vertical velocity standard deviation              [     m/s]
           ,uwind          & ! Zonal wind speed                                  [     m/s]
           ,vwind          & ! Meridional wind speed                             [     m/s]
           ,wwind          ! ! Vertical velocity                                 [     m/s]
   !---------------------------------------------------------------------------------------!



   !------ 1D dependence (mgmzp), variables at current time step --------------------------!
   real, allocatable, dimension(:) :: &
            exner0    & ! Exner function                                         [  J/kg/K]
           ,p0        & ! Pressure with forcing                                  [     hPa]
           ,qtot0     & ! Total mixing ratio                                     [   kg/kg]
           ,qice0     & ! Ice mixing ratio                                       [   kg/kg]
           ,qliq0     & ! Liquid water mixing ratio                              [   kg/kg]
           ,qvap0     & ! Water vapour mixing ratio                              [   kg/kg]
           ,t0        & ! Temperature                                            [       K]
           ,thil0     & ! Ice-liquid potential temperature                       [       K]
           ,theiv0    & ! Ice-vapour equivalent potential temperature            [       K]
           ,tke0      ! ! Turbulent Kinetic Energy                               [    J/kg]
   !---------------------------------------------------------------------------------------!


   !------ 1D dependence (mgmzp), forcing due to convectivion -----------------------------!
   real, allocatable, dimension(:) :: &
            outqtot   & ! Total water mixing ratio tendency due to cumulus       [ kg/kg/s]
           ,outthil   ! ! Ice-liquid potential temperature tendency due to Cu    [     K/s]  

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

     allocate (cdd           (mgmzp))
     allocate (cdu           (mgmzp))
     allocate (mentrd_rate   (mgmzp))
     allocate (mentru_rate   (mgmzp))
     allocate (dbyd          (mgmzp))
     allocate (dbyu          (mgmzp))
     allocate (etad_cld      (mgmzp))
     allocate (etau_cld      (mgmzp))
     allocate (rhod_cld      (mgmzp))
     allocate (rhou_cld      (mgmzp))
     allocate (qliqd_cld     (mgmzp))
     allocate (qliqu_cld     (mgmzp))
     allocate (qiced_cld     (mgmzp))
     allocate (qiceu_cld     (mgmzp))


     allocate (dqtotdt       (mgmzp))
     allocate (dqtotdt_shal  (mgmzp))
     allocate (dthildt       (mgmzp))
     allocate (dthildt_shal  (mgmzp))
     allocate (dtkedt        (mgmzp))
     allocate (exner         (mgmzp))
     allocate (omeg          (mgmzp))
     allocate (p             (mgmzp))
     allocate (qtot          (mgmzp))
     allocate (qvap          (mgmzp))
     allocate (qliq          (mgmzp))
     allocate (qice          (mgmzp))
     allocate (rho           (mgmzp))
     allocate (t             (mgmzp))
     allocate (thil          (mgmzp))
     allocate (theiv         (mgmzp))
     allocate (tke           (mgmzp))
     allocate (sigw          (mgmzp))
     allocate (uwind         (mgmzp))
     allocate (vwind         (mgmzp))
     allocate (wwind         (mgmzp))

     allocate (exner0        (mgmzp))
     allocate (p0            (mgmzp))
     allocate (qtot0         (mgmzp))
     allocate (qvap0         (mgmzp))
     allocate (qliq0         (mgmzp))
     allocate (qice0         (mgmzp))
     allocate (t0            (mgmzp))
     allocate (thil0         (mgmzp))
     allocate (theiv0        (mgmzp))
     allocate (tke0          (mgmzp))
     
     allocate (outqtot       (mgmzp))
     allocate (outthil       (mgmzp))

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

     if(allocated(cdd           )) deallocate (cdd           )
     if(allocated(cdu           )) deallocate (cdu           )
     if(allocated(dbyd          )) deallocate (dbyd          )
     if(allocated(dbyu          )) deallocate (dbyu          )
     if(allocated(mentrd_rate   )) deallocate (mentrd_rate   )
     if(allocated(mentru_rate   )) deallocate (mentru_rate   )
     if(allocated(etad_cld      )) deallocate (etad_cld      )
     if(allocated(etau_cld      )) deallocate (etau_cld      )
     if(allocated(rhod_cld      )) deallocate (rhod_cld      )
     if(allocated(rhou_cld      )) deallocate (rhou_cld      )
     if(allocated(qiced_cld     )) deallocate (qiced_cld     )
     if(allocated(qiceu_cld     )) deallocate (qiceu_cld     )
     if(allocated(qliqd_cld     )) deallocate (qliqd_cld     )
     if(allocated(qliqu_cld     )) deallocate (qliqu_cld     )

     if(allocated(dqtotdt       )) deallocate (dqtotdt       )
     if(allocated(dqtotdt_shal  )) deallocate (dqtotdt_shal  )
     if(allocated(dthildt       )) deallocate (dthildt       )
     if(allocated(dthildt_shal  )) deallocate (dthildt_shal  )
     if(allocated(dtkedt        )) deallocate (dtkedt        )
     if(allocated(exner         )) deallocate (exner         )
     if(allocated(omeg          )) deallocate (omeg          )
     if(allocated(p             )) deallocate (p             )
     if(allocated(qtot          )) deallocate (qtot          )
     if(allocated(qvap          )) deallocate (qvap          )
     if(allocated(qliq          )) deallocate (qliq          )
     if(allocated(qice          )) deallocate (qice          )
     if(allocated(rho           )) deallocate (rho           )
     if(allocated(t             )) deallocate (t             )
     if(allocated(thil          )) deallocate (thil          )
     if(allocated(theiv         )) deallocate (theiv         )
     if(allocated(tke           )) deallocate (tke           )
     if(allocated(sigw          )) deallocate (sigw          )
     if(allocated(uwind         )) deallocate (uwind         )
     if(allocated(vwind         )) deallocate (vwind         )
     if(allocated(wwind         )) deallocate (wwind         )

     if(allocated(exner0        )) deallocate (exner0        )
     if(allocated(p0            )) deallocate (p0            )
     if(allocated(qtot0         )) deallocate (qtot0         )
     if(allocated(qvap0         )) deallocate (qvap0         )
     if(allocated(qliq0         )) deallocate (qliq0         )
     if(allocated(qice0         )) deallocate (qice0         )
     if(allocated(t0            )) deallocate (t0            )
     if(allocated(thil0         )) deallocate (thil0         )
     if(allocated(theiv0        )) deallocate (theiv0        )
     if(allocated(tke0          )) deallocate (tke0          )

     if(allocated(outqtot       )) deallocate (outqtot       )
     if(allocated(outthil       )) deallocate (outthil       )

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

     if(allocated(cdd           )) cdd           = 0.
     if(allocated(cdu           )) cdu           = 0.
     if(allocated(dbyd          )) dbyd          = 0.
     if(allocated(dbyu          )) dbyu          = 0.
     if(allocated(mentrd_rate   )) mentrd_rate   = 0.
     if(allocated(mentru_rate   )) mentru_rate   = 0.
     if(allocated(etad_cld      )) etad_cld      = 0.
     if(allocated(etau_cld      )) etau_cld      = 0.
     if(allocated(rhod_cld      )) rhod_cld      = 0.
     if(allocated(rhou_cld      )) rhou_cld      = 0.
     if(allocated(qiced_cld     )) qiced_cld     = 0.
     if(allocated(qiceu_cld     )) qiceu_cld     = 0.
     if(allocated(qliqd_cld     )) qliqd_cld     = 0.
     if(allocated(qliqu_cld     )) qliqu_cld     = 0.

     if(allocated(dqtotdt       )) dqtotdt       = 0.
     if(allocated(dqtotdt_shal  )) dqtotdt_shal  = 0.
     if(allocated(dthildt       )) dthildt       = 0.
     if(allocated(dthildt_shal  )) dthildt_shal  = 0.
     if(allocated(dtkedt        )) dtkedt        = 0.
     if(allocated(exner         )) exner         = 0.
     if(allocated(omeg          )) omeg          = 0.
     if(allocated(p             )) p             = 0.
     if(allocated(qtot          )) qtot          = 0.
     if(allocated(qvap          )) qvap          = 0.
     if(allocated(qliq          )) qliq          = 0.
     if(allocated(qice          )) qice          = 0.
     if(allocated(rho           )) rho           = 0.
     if(allocated(t             )) t             = 0.
     if(allocated(thil          )) thil          = 0.
     if(allocated(theiv         )) theiv         = 0.
     if(allocated(tke           )) tke           = 0.
     if(allocated(sigw          )) sigw          = 0.
     if(allocated(uwind         )) uwind         = 0.
     if(allocated(vwind         )) vwind         = 0.
     if(allocated(wwind         )) wwind         = 0.

     if(allocated(exner0        )) exner0        = 0.
     if(allocated(p0            )) p0            = 0.
     if(allocated(qtot0         )) qtot0         = 0.
     if(allocated(qvap0         )) qvap0         = 0.
     if(allocated(qliq0         )) qliq0         = 0.
     if(allocated(qice0         )) qice0         = 0.
     if(allocated(t0            )) t0            = 0.
     if(allocated(thil0         )) thil0         = 0.
     if(allocated(theiv0        )) theiv0        = 0.
     if(allocated(tke0          )) tke0          = 0.

     if(allocated(outqtot       )) outqtot       = 0.
     if(allocated(outthil       )) outthil       = 0.
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
     exnersur          = 0.
     psur              = 0.
     qtotsur           = 0.
     qvapsur           = 0.
     qliqsur           = 0.
     qicesur           = 0.
     tsur              = 0.
     theivsur          = 0.
     thilsur           = 0.

     dens_curr         = 0.
     mconv             = 0.
     prev_dnmf         = 0.

     precip            = 0.
     wbuoymin          = 0.

     !----- Logical variables -------------------------------------------------------------!
     upconv            = .false.
     !-------------------------------------------------------------------------------------!

     return

  end subroutine zero_scratch_grell
!==========================================================================================!
!==========================================================================================!
end module mem_scratch_grell
