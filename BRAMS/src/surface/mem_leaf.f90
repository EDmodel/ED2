!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
Module mem_leaf
   use grid_dims
  
   type leaf_vars
      
      !------------------------------------------------------------------------------------!
      !   Soil properties, dimensioned by (nzg,nxp,nyp,npatch), except for roughness and   !
      ! colour which are dimensioned by (nxp,nyp,npatch).                                  !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:,:), pointer :: soil_water  & ! Soil moisture        [    m³/m³]
                                         , soil_energy & ! Internal energy      [     J/m³]
                                         , soil_text   ! ! Soil texture class   [      ---]
      real, dimension(  :,:,:), pointer :: soil_rough  ! ! Soil roughness       [        m]
      real, dimension(  :,:,:), pointer :: soil_color  ! ! Soil colour          [      ---]

      !------------------------------------------------------------------------------------!
      !    Temporary surface water or snow layer properties, dimensioned by                !
      ! (nzs,nxp,nyp,npatch), except for the number of levels which is dimensioned by      !
      ! (nxp,nyp,npatch).                                                                  !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:,:), pointer :: sfcwater_mass   & ! Mass             [    kg/m²]
                                         , sfcwater_energy & ! Internal energy  [     J/kg]
                                         , sfcwater_depth  ! ! Layer depth      [        m]
      real, dimension(  :,:,:), pointer :: sfcwater_nlev   ! ! # of layers      [        m]


      !------------------------------------------------------------------------------------!
      !    Ground properties. It may refer to either the top soil layer or the snow/water  !
      ! top layer, depending on whether there is a water/snow layer or not.  Dimensioned   !
      ! by (nxp,nyp,npatch).                                                               !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: ground_rsat & ! Sat. mixing ratio      [    kg/kg]
                                       , ground_rvap & ! Vapour mixing ratio    [    kg/kg]
                                       , ground_temp & ! Temperature            [        K]
                                       , ground_fliq ! ! Liquid water fraction  [      ---]

      !------------------------------------------------------------------------------------!
      !     Vegetation properties, dimensioned by (nxp,nyp,npatch).                        !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: veg_fracarea & ! Fractional area       [      ---]
                                       , veg_lai      & ! Leaf area index       [    m²/m²]
                                       , veg_agb      & ! Above-ground biomass  [   kgC/m²]
                                       , veg_rough    & ! Roughness length      [        m]
                                       , veg_height   & ! Height                [        m]
                                       , veg_displace & ! Displacement height   [        m]
                                       , veg_albedo   & ! Albedo                [      ---]
                                       , veg_tai      & ! Tree area index       [    m²/m²]
                                       , veg_water    & ! Leaf surface water    [    kg/m²]
                                       , veg_hcap     & ! Heat capacity         [   J/m²/K]
                                       , veg_energy   & ! Internal energy       [     J/m²]
                                       , veg_ndvip    & ! Past NDVI             [      ---]
                                       , veg_ndvic    & ! Current NDVI          [      ---]
                                       , veg_ndvif    & ! Future NDVI           [      ---]
                                       , leaf_class   & ! Vegetation class      [      ---]
                                       , stom_condct  ! ! Stomatal resistance   [      ???]


      !------------------------------------------------------------------------------------!
      !     Canopy air properties, dimensioned by (nxp,nyp,npatch).                        !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: can_rvap     & ! Vapour Mixing ratio   [    kg/kg]
                                       , can_theta    & ! Potential temperature [        K]
                                       , can_theiv    & ! Theta_Eiv             [        K]
                                       , can_vpdef    & ! Vapour press. deficit [       Pa]
                                       , can_prss     & ! Pressure              [       Pa]
                                       , can_co2      ! ! CO2 mixing ratio      [ µmol/mol]


      !------------------------------------------------------------------------------------!
      !    "Star" variables. The stars are the characteristic friction scale, and are      !
      ! dimensioned by (nxp,nyp,npatch).                                                   !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: ustar & ! Friction velocity            [      m/s]
                                       , tstar & ! Potential temperature scale  [        K]
                                       , estar & ! Enthalpy scale               [     J/kg]
                                       , rstar & ! Mixing ratio scale           [    kg/kg]
                                       , cstar ! ! CO2 mixing ratio scale       [ µmol/mol]

      !------------------------------------------------------------------------------------!
      !    Auxilliary variables for diagnostics.                                           !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: zeta   & ! Dimensionless turbulent height scale
                                       , ribulk ! ! Bulk Richardson number

      !------------------------------------------------------------------------------------!
      !     Patch structural properties, dimensioned by (nxp,nyp,npatch).                  !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: patch_area   & ! Patch relative area   [      ---]
                                       , patch_rough  & ! Roughness length      [        m]
                                       , patch_wetind ! ! Wetness index         [      ???]

      !------------------------------------------------------------------------------------!
      !     Surface fluxes, dimensioned by (nxp,nyp,npatch).                               !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: gpp         & ! Gross primary prod.    [µmol/m²/s]
                                       , resphet     & ! Heterotrophic resp.    [µmol/m²/s]
                                       , plresp      & ! Plant respiration      [µmol/m²/s]
                                       , evap_gc     & ! Evaporation (Gnd->Can) [     W/m²]
                                       , evap_vc     & ! Evaporation (Veg->Can) [     W/m²]
                                       , transp      & ! Transpiration          [     W/m²]
                                       , sensible_gc & ! Sens. heat (Gnd->Can)  [     W/m²]
                                       , sensible_vc ! ! Sens. heat (Veg->Can)  [     W/m²]

      !------------------------------------------------------------------------------------!
      !     This is based on a 10-day running average of the relative soil moisture        !
      ! potential in the root zone, and it is used by all vegetation types that have       !
      ! drought phenology (phenology(nveg) = 4) and LAI is not to be read from standard    !
      ! files.                                                                             !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: psibar_10d

      !------------------------------------------------------------------------------------!
      !     Radiation variables, for diagnostics only (nxp,nyp,npatch).                    !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: rshort_gnd  & ! Absorbed SW radiation  [     W/m²]
                                       , rlong_gnd   ! ! Absorbed LW radiation  [     W/m²]

      !------------------------------------------------------------------------------------!
      !     Miscellaneous properties, dimensioned by (nxp,nyp,npatch).                     !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: R_aer   & ! Aerodynamic resistance     [      ???]
                                       , G_URBAN ! ! TEB_SPM variable...        [      ???]

      !------------------------------------------------------------------------------------!
      !     Miscellaneous variables to be dimensioned by (nxp,nyp).                                          !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:)  , pointer :: snow_mass  & ! Total snow mass         [    kg/m²]
                                       , snow_depth & ! Total snow depth        [        m]
                                       , seatp      & ! Past sea surface temp.  [        K]
                                       , seatf      ! ! Future sea sfc. temp.   [        K]
   end type leaf_vars
  
   type (leaf_vars), dimension(:), allocatable :: leaf_g  & ! Instantaneous variables.
                                                , leafm_g ! ! Averaged variables.

   !---------------------------------------------------------------------------------------!
   !     Other variables.                                                                  !
   !---------------------------------------------------------------------------------------!
   integer                 :: nslcon      ! Soil texture if constant for entire domain
   integer                 :: isoilcol    ! Number of vegetation types
   integer                 :: nvgcon      ! Vegetation class if constant for entire domain
   integer                 :: nvegpat     ! Number of vegetation types
   integer                 :: isfcl       ! Surface model (1. LEAF3, 2. LEAF-Hydro, 5. ED2)
   integer                 :: istar       ! Which surface layer model should I use?
                                          !    1. Louis (1979)
                                          !    2. Oncley and Dudhia (1995).
                                          !    3. Beljaars and Holtslag (1991)
                                          !    4. BH91, using OD95 to find zeta.
                                          !    5. OD95, using BH91 to find zeta
   integer                 :: igrndvap    ! Methods to find the ground -> canopy 
                                          !     conductance.  In all cases the beta term 
                                          !     is modified so it approaches zero as soil 
                                          !     moisture goes to dry air soil. 
                                          ! 0. Modified Lee Pielke (1992), adding field
                                          !    capacity, but using beta factor without the
                                          !    square, like in Noilhan and Planton (1989).
                                          !    This is the closest to the original ED-2.1
                                          ! 1. Test # 1 of Mahfouf and Noilhan (1991)
                                          ! 2. Test # 2 of Mahfouf and Noilhan (1991)
                                          ! 3. Test # 3 of Mahfouf and Noilhan (1991)
                                          ! 4. Test # 4 of Mahfouf and Noilhan (1991)
   integer                 :: isoilbc     ! Bottom soil boundary condition.
                                          !    0. Flat Bedrock (zero flow)
                                          !    1. Free drainage (gravity flow)
                                          !    2. Sink hole drainage (BC at -3.1MPa+z)
                                          !    3. Water table (BC at field capacity) 
                                          !    4. Aquifer (BC at porosity)
                                          !    5. Lateral drainage (gravity flow at slope)
                                          !       This requires sldrain, in future options
                                          !       0, 1, and 5 could be merged.
   real                    :: sldrain     ! Slope for lateral drainage in degrees.  Values
                                          !       can range between 0 (flat bedrock) and 90
                                          !       degrees (free drainage).
   integer                 :: ipercol     ! Percolation scheme:
                                          !    0. Original method, from LEAF-3.  Shed liquid 
                                          !       in excess of a 1:9 liquid-to-ice ratio 
                                          !       through percolation.
                                          !    1. Alternative "free" water calculation.  
                                          !       Anderson (1976), NOAA Tech Report NWS 19.
   real                    :: runoff_time ! Runoff time scale.
   real                    :: dtleaf      ! LEAF-3 target time step.  It will be either this
                                          !    this or the actual BRAMS time step (which-
                                          !    ever is the lowest).

   real, dimension(nzgmax) :: stgoff  ! Initial soil temperature offset
   real, dimension(nzgmax) :: slmstr  ! Initial soil moisture if constant for entire domain
   real, dimension(nzgmax) :: slz     ! Soil levels
   real                    :: zrough  ! Roughness if constant for entire domain
   real                    :: pctlcon ! Vegetation fraction if constant for entire domain
   real                    :: albedo  ! Albedo if constant for entire domain.
   real                    :: drtcon  ! Delta-Rt if constant for entire domain.
   real                    :: dthcon  ! Delta-Theta if constant for entire domain.
   real                    :: seatmp  ! Sea temperature if constant for entire domain.
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will allocate the pointers of the leaf structure, based on        !
   ! options if necessary.                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_leaf(leaf,nz,nx,ny,nzg,nzs,np,ng,teb_spm)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (leaf_vars), intent(inout) :: leaf
      integer         , intent(in)    :: nz
      integer         , intent(in)    :: nx
      integer         , intent(in)    :: ny
      integer         , intent(in)    :: nzg
      integer         , intent(in)    :: nzs
      integer         , intent(in)    :: np
      integer         , intent(in)    :: ng
      integer         , intent(in)    :: teb_spm
      !------------------------------------------------------------------------------------!


      allocate (leaf%soil_water       (nzg,nx,ny,np))
      allocate (leaf%soil_energy      (nzg,nx,ny,np))
      allocate (leaf%soil_text        (nzg,nx,ny,np))
      allocate (leaf%soil_rough       (    nx,ny,np))
      allocate (leaf%soil_color       (    nx,ny,np))

      allocate (leaf%sfcwater_mass    (nzs,nx,ny,np))
      allocate (leaf%sfcwater_energy  (nzs,nx,ny,np))
      allocate (leaf%sfcwater_depth   (nzs,nx,ny,np))
      allocate (leaf%sfcwater_nlev    (    nx,ny,np))

      allocate (leaf%ground_rsat      (    nx,ny,np))
      allocate (leaf%ground_rvap      (    nx,ny,np))
      allocate (leaf%ground_temp      (    nx,ny,np))
      allocate (leaf%ground_fliq      (    nx,ny,np))

      allocate (leaf%veg_fracarea     (    nx,ny,np))
      allocate (leaf%veg_lai          (    nx,ny,np))
      allocate (leaf%veg_agb          (    nx,ny,np))
      allocate (leaf%veg_rough        (    nx,ny,np))
      allocate (leaf%veg_height       (    nx,ny,np))
      allocate (leaf%veg_displace     (    nx,ny,np))
      allocate (leaf%veg_albedo       (    nx,ny,np))
      allocate (leaf%veg_tai          (    nx,ny,np))
      allocate (leaf%veg_water        (    nx,ny,np))
      allocate (leaf%veg_hcap         (    nx,ny,np))
      allocate (leaf%veg_energy       (    nx,ny,np))
      allocate (leaf%veg_ndvip        (    nx,ny,np))
      allocate (leaf%veg_ndvic        (    nx,ny,np))
      allocate (leaf%veg_ndvif        (    nx,ny,np))
      allocate (leaf%leaf_class       (    nx,ny,np))
      allocate (leaf%stom_condct      (    nx,ny,np))

      allocate (leaf%can_rvap         (    nx,ny,np))
      allocate (leaf%can_theta        (    nx,ny,np))
      allocate (leaf%can_theiv        (    nx,ny,np))
      allocate (leaf%can_vpdef        (    nx,ny,np))
      allocate (leaf%can_prss         (    nx,ny,np))
      allocate (leaf%can_co2          (    nx,ny,np))

      allocate (leaf%ustar            (    nx,ny,np))
      allocate (leaf%tstar            (    nx,ny,np))
      allocate (leaf%estar            (    nx,ny,np))
      allocate (leaf%rstar            (    nx,ny,np))
      allocate (leaf%cstar            (    nx,ny,np))

      allocate (leaf%zeta             (    nx,ny,np))
      allocate (leaf%ribulk           (    nx,ny,np))

      allocate (leaf%patch_area       (    nx,ny,np))
      allocate (leaf%patch_rough      (    nx,ny,np))
      allocate (leaf%patch_wetind     (    nx,ny,np))


      allocate (leaf%gpp              (    nx,ny,np))
      allocate (leaf%resphet          (    nx,ny,np))
      allocate (leaf%plresp           (    nx,ny,np))
      allocate (leaf%evap_gc          (    nx,ny,np))
      allocate (leaf%evap_vc          (    nx,ny,np))
      allocate (leaf%transp           (    nx,ny,np))
      allocate (leaf%sensible_gc      (    nx,ny,np))
      allocate (leaf%sensible_vc      (    nx,ny,np))
      allocate (leaf%psibar_10d       (    nx,ny,np))

      allocate (leaf%rshort_gnd       (    nx,ny,np))
      allocate (leaf%rlong_gnd        (    nx,ny,np))

      allocate (leaf%R_aer            (    nx,ny,np))

      if (teb_spm == 1) then
         allocate (leaf%G_URBAN       (    nx,ny,np))
      end if

      allocate (leaf%snow_mass        (    nx,ny   ))
      allocate (leaf%snow_depth       (    nx,ny   ))
      allocate (leaf%seatp            (    nx,ny   ))
      allocate (leaf%seatf            (    nx,ny   ))

      return
   end subroutine alloc_leaf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will nullify all pointers for a safe allocation.                   !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_leaf(leaf)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (leaf_vars), intent(inout) :: leaf
      !------------------------------------------------------------------------------------!


      nullify(leaf%soil_water       )
      nullify(leaf%soil_energy      )
      nullify(leaf%soil_text        )
      nullify(leaf%soil_rough       )
      nullify(leaf%soil_color       )

      nullify(leaf%sfcwater_mass    )
      nullify(leaf%sfcwater_energy  )
      nullify(leaf%sfcwater_depth   )
      nullify(leaf%sfcwater_nlev    )

      nullify(leaf%ground_rsat      )
      nullify(leaf%ground_rvap      )
      nullify(leaf%ground_temp      )
      nullify(leaf%ground_fliq      )

      nullify(leaf%veg_fracarea     )
      nullify(leaf%veg_lai          )
      nullify(leaf%veg_agb          )
      nullify(leaf%veg_rough        )
      nullify(leaf%veg_height       )
      nullify(leaf%veg_displace     )
      nullify(leaf%veg_albedo       )
      nullify(leaf%veg_tai          )
      nullify(leaf%veg_water        )
      nullify(leaf%veg_hcap         )
      nullify(leaf%veg_energy       )
      nullify(leaf%veg_ndvip        )
      nullify(leaf%veg_ndvic        )
      nullify(leaf%veg_ndvif        )
      nullify(leaf%leaf_class       )
      nullify(leaf%stom_condct      )

      nullify(leaf%can_rvap         )
      nullify(leaf%can_theta        )
      nullify(leaf%can_theiv        )
      nullify(leaf%can_vpdef        )
      nullify(leaf%can_prss         )
      nullify(leaf%can_co2          )

      nullify(leaf%ustar            )
      nullify(leaf%tstar            )
      nullify(leaf%estar            )
      nullify(leaf%rstar            )
      nullify(leaf%cstar            )

      nullify(leaf%zeta             )
      nullify(leaf%ribulk           )

      nullify(leaf%patch_area       )
      nullify(leaf%patch_rough      )
      nullify(leaf%patch_wetind     )


      nullify(leaf%gpp              )
      nullify(leaf%resphet          )
      nullify(leaf%plresp           )
      nullify(leaf%evap_gc          )
      nullify(leaf%evap_vc          )
      nullify(leaf%transp           )
      nullify(leaf%sensible_gc      )
      nullify(leaf%sensible_vc      )
      nullify(leaf%psibar_10d       )

      nullify(leaf%rshort_gnd       )
      nullify(leaf%rlong_gnd        )

      nullify(leaf%R_aer            )
      nullify(leaf%G_URBAN          )

      nullify(leaf%snow_mass        )
      nullify(leaf%snow_depth       )
      nullify(leaf%seatp            )
      nullify(leaf%seatf            )

      return
   end subroutine nullify_leaf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will assign zeroes to all variables.  This is to avoid some bogus  !
   ! values to variables that are never used.  The model only updates variables for        !
   ! patches that have a minimum area.                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine zero_leaf(leaf)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (leaf_vars), intent(inout) :: leaf
      !------------------------------------------------------------------------------------!


      if (associated(leaf%soil_water       ))  leaf%soil_water       = 0.0
      if (associated(leaf%soil_energy      ))  leaf%soil_energy      = 0.0
      if (associated(leaf%soil_text        ))  leaf%soil_text        = 0.0
      if (associated(leaf%soil_rough       ))  leaf%soil_rough       = 0.0
      if (associated(leaf%soil_color       ))  leaf%soil_color       = 0.0

      if (associated(leaf%sfcwater_mass    ))  leaf%sfcwater_mass    = 0.0
      if (associated(leaf%sfcwater_energy  ))  leaf%sfcwater_energy  = 0.0
      if (associated(leaf%sfcwater_depth   ))  leaf%sfcwater_depth   = 0.0
      if (associated(leaf%sfcwater_nlev    ))  leaf%sfcwater_nlev    = 0.0

      if (associated(leaf%ground_rsat      ))  leaf%ground_rsat      = 0.0
      if (associated(leaf%ground_rvap      ))  leaf%ground_rvap      = 0.0
      if (associated(leaf%ground_temp      ))  leaf%ground_temp      = 0.0
      if (associated(leaf%ground_fliq      ))  leaf%ground_fliq      = 0.0

      if (associated(leaf%veg_fracarea     ))  leaf%veg_fracarea     = 0.0
      if (associated(leaf%veg_lai          ))  leaf%veg_lai          = 0.0
      if (associated(leaf%veg_agb          ))  leaf%veg_agb          = 0.0
      if (associated(leaf%veg_rough        ))  leaf%veg_rough        = 0.0
      if (associated(leaf%veg_height       ))  leaf%veg_height       = 0.0
      if (associated(leaf%veg_displace     ))  leaf%veg_displace     = 0.0
      if (associated(leaf%veg_albedo       ))  leaf%veg_albedo       = 0.0
      if (associated(leaf%veg_tai          ))  leaf%veg_tai          = 0.0
      if (associated(leaf%veg_water        ))  leaf%veg_water        = 0.0
      if (associated(leaf%veg_hcap         ))  leaf%veg_hcap         = 0.0
      if (associated(leaf%veg_energy       ))  leaf%veg_energy       = 0.0
      if (associated(leaf%veg_ndvip        ))  leaf%veg_ndvip        = 0.0
      if (associated(leaf%veg_ndvic        ))  leaf%veg_ndvic        = 0.0
      if (associated(leaf%veg_ndvif        ))  leaf%veg_ndvif        = 0.0
      if (associated(leaf%leaf_class       ))  leaf%leaf_class       = 0.0
      if (associated(leaf%stom_condct      ))  leaf%stom_condct      = 0.0

      if (associated(leaf%can_rvap         ))  leaf%can_rvap         = 0.0
      if (associated(leaf%can_theta        ))  leaf%can_theta        = 0.0
      if (associated(leaf%can_theiv        ))  leaf%can_theiv        = 0.0
      if (associated(leaf%can_vpdef        ))  leaf%can_vpdef        = 0.0
      if (associated(leaf%can_prss         ))  leaf%can_prss         = 0.0
      if (associated(leaf%can_co2          ))  leaf%can_co2          = 0.0

      if (associated(leaf%ustar            ))  leaf%ustar            = 0.0
      if (associated(leaf%tstar            ))  leaf%tstar            = 0.0
      if (associated(leaf%estar            ))  leaf%estar            = 0.0
      if (associated(leaf%rstar            ))  leaf%rstar            = 0.0
      if (associated(leaf%cstar            ))  leaf%cstar            = 0.0

      if (associated(leaf%zeta             ))  leaf%zeta             = 0.0
      if (associated(leaf%ribulk           ))  leaf%ribulk           = 0.0

      if (associated(leaf%patch_area       ))  leaf%patch_area       = 0.0
      if (associated(leaf%patch_rough      ))  leaf%patch_rough      = 0.0
      if (associated(leaf%patch_wetind     ))  leaf%patch_wetind     = 0.0


      if (associated(leaf%gpp              ))  leaf%gpp              = 0.0
      if (associated(leaf%resphet          ))  leaf%resphet          = 0.0
      if (associated(leaf%plresp           ))  leaf%plresp           = 0.0
      if (associated(leaf%evap_gc          ))  leaf%evap_gc          = 0.0
      if (associated(leaf%evap_vc          ))  leaf%evap_vc          = 0.0
      if (associated(leaf%transp           ))  leaf%transp           = 0.0
      if (associated(leaf%sensible_gc      ))  leaf%sensible_gc      = 0.0
      if (associated(leaf%sensible_vc      ))  leaf%sensible_vc      = 0.0
      if (associated(leaf%psibar_10d       ))  leaf%psibar_10d       = 0.0

      if (associated(leaf%rshort_gnd       ))  leaf%rshort_gnd       = 0.0
      if (associated(leaf%rlong_gnd        ))  leaf%rlong_gnd        = 0.0

      if (associated(leaf%R_aer            ))  leaf%R_aer            = 0.0
      if (associated(leaf%G_URBAN          ))  leaf%G_URBAN          = 0.0

      if (associated(leaf%snow_mass        ))  leaf%snow_mass        = 0.0
      if (associated(leaf%snow_depth       ))  leaf%snow_depth       = 0.0
      if (associated(leaf%seatp            ))  leaf%seatp            = 0.0
      if (associated(leaf%seatf            ))  leaf%seatf            = 0.0

      return
   end subroutine zero_leaf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will deallocate all pointers before deallocating the structure.   !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_leaf(leaf)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (leaf_vars), intent(inout) :: leaf
      !------------------------------------------------------------------------------------!


      if (associated(leaf%soil_water       ))  deallocate(leaf%soil_water       )
      if (associated(leaf%soil_energy      ))  deallocate(leaf%soil_energy      )
      if (associated(leaf%soil_text        ))  deallocate(leaf%soil_text        )
      if (associated(leaf%soil_rough       ))  deallocate(leaf%soil_rough       )
      if (associated(leaf%soil_color       ))  deallocate(leaf%soil_color       )

      if (associated(leaf%sfcwater_mass    ))  deallocate(leaf%sfcwater_mass    )
      if (associated(leaf%sfcwater_energy  ))  deallocate(leaf%sfcwater_energy  )
      if (associated(leaf%sfcwater_depth   ))  deallocate(leaf%sfcwater_depth   )
      if (associated(leaf%sfcwater_nlev    ))  deallocate(leaf%sfcwater_nlev    )

      if (associated(leaf%ground_rsat      ))  deallocate(leaf%ground_rsat      )
      if (associated(leaf%ground_rvap      ))  deallocate(leaf%ground_rvap      )
      if (associated(leaf%ground_temp      ))  deallocate(leaf%ground_temp      )
      if (associated(leaf%ground_fliq      ))  deallocate(leaf%ground_fliq      )

      if (associated(leaf%veg_fracarea     ))  deallocate(leaf%veg_fracarea     )
      if (associated(leaf%veg_lai          ))  deallocate(leaf%veg_lai          )
      if (associated(leaf%veg_agb          ))  deallocate(leaf%veg_agb          )
      if (associated(leaf%veg_rough        ))  deallocate(leaf%veg_rough        )
      if (associated(leaf%veg_height       ))  deallocate(leaf%veg_height       )
      if (associated(leaf%veg_displace     ))  deallocate(leaf%veg_displace     )
      if (associated(leaf%veg_albedo       ))  deallocate(leaf%veg_albedo       )
      if (associated(leaf%veg_tai          ))  deallocate(leaf%veg_tai          )
      if (associated(leaf%veg_water        ))  deallocate(leaf%veg_water        )
      if (associated(leaf%veg_hcap         ))  deallocate(leaf%veg_hcap         )
      if (associated(leaf%veg_energy       ))  deallocate(leaf%veg_energy       )
      if (associated(leaf%veg_ndvip        ))  deallocate(leaf%veg_ndvip        )
      if (associated(leaf%veg_ndvic        ))  deallocate(leaf%veg_ndvic        )
      if (associated(leaf%veg_ndvif        ))  deallocate(leaf%veg_ndvif        )
      if (associated(leaf%leaf_class       ))  deallocate(leaf%leaf_class       )
      if (associated(leaf%stom_condct      ))  deallocate(leaf%stom_condct      )

      if (associated(leaf%can_rvap         ))  deallocate(leaf%can_rvap         )
      if (associated(leaf%can_theta        ))  deallocate(leaf%can_theta        )
      if (associated(leaf%can_theiv        ))  deallocate(leaf%can_theiv        )
      if (associated(leaf%can_vpdef        ))  deallocate(leaf%can_vpdef        )
      if (associated(leaf%can_prss         ))  deallocate(leaf%can_prss         )
      if (associated(leaf%can_co2          ))  deallocate(leaf%can_co2          )

      if (associated(leaf%ustar            ))  deallocate(leaf%ustar            )
      if (associated(leaf%tstar            ))  deallocate(leaf%tstar            )
      if (associated(leaf%estar            ))  deallocate(leaf%estar            )
      if (associated(leaf%rstar            ))  deallocate(leaf%rstar            )
      if (associated(leaf%cstar            ))  deallocate(leaf%cstar            )

      if (associated(leaf%zeta             ))  deallocate(leaf%zeta             )
      if (associated(leaf%ribulk           ))  deallocate(leaf%ribulk           )

      if (associated(leaf%patch_area       ))  deallocate(leaf%patch_area       )
      if (associated(leaf%patch_rough      ))  deallocate(leaf%patch_rough      )
      if (associated(leaf%patch_wetind     ))  deallocate(leaf%patch_wetind     )


      if (associated(leaf%gpp              ))  deallocate(leaf%gpp              )
      if (associated(leaf%resphet          ))  deallocate(leaf%resphet          )
      if (associated(leaf%plresp           ))  deallocate(leaf%plresp           )
      if (associated(leaf%evap_vc          ))  deallocate(leaf%evap_vc          )
      if (associated(leaf%evap_gc          ))  deallocate(leaf%evap_gc          )
      if (associated(leaf%transp           ))  deallocate(leaf%transp           )
      if (associated(leaf%sensible_gc      ))  deallocate(leaf%sensible_gc      )
      if (associated(leaf%sensible_vc      ))  deallocate(leaf%sensible_vc      )
      if (associated(leaf%psibar_10d       ))  deallocate(leaf%psibar_10d       )

      if (associated(leaf%rshort_gnd       ))  deallocate(leaf%rshort_gnd       )
      if (associated(leaf%rlong_gnd        ))  deallocate(leaf%rlong_gnd        )

      if (associated(leaf%R_aer            ))  deallocate(leaf%R_aer            )
      if (associated(leaf%G_URBAN          ))  deallocate(leaf%G_URBAN          )

      if (associated(leaf%snow_mass        ))  deallocate(leaf%snow_mass        )
      if (associated(leaf%snow_depth       ))  deallocate(leaf%snow_depth       )
      if (associated(leaf%seatp            ))  deallocate(leaf%seatp            )
      if (associated(leaf%seatf            ))  deallocate(leaf%seatf            )

      return
   end subroutine dealloc_leaf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will fill pointers to arrays into variable tables.                !
   !---------------------------------------------------------------------------------------!
   subroutine filltab_leaf(leaf,leafm,imean,nmz,nmx,nmy,nmzg,nmzs,nmpat,ng)
      use var_tables
      use io_params, only: ipastin ! INTENT(IN)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (leaf_vars), intent(inout) :: leaf
      type (leaf_vars), intent(inout) :: leafm
      integer         , intent(in)    :: imean
      integer         , intent(in)    :: nmz
      integer         , intent(in)    :: nmx
      integer         , intent(in)    :: nmy
      integer         , intent(in)    :: nmzg
      integer         , intent(in)    :: nmzs
      integer         , intent(in)    :: nmpat
      integer         , intent(in)    :: ng
      !----- Local variables. -------------------------------------------------------------!
      integer                         :: npts
      character(len=8)                :: str_recycle
      !------------------------------------------------------------------------------------!

      !----- Deciding whether we add the recycle argument in the keys. --------------------!
      select case (ipastin)
      case (1)
         str_recycle = 'recycle'
      case default
         str_recycle = ''
      end select

      !------------------------------------------------------------------------------------!
      !     4-D variables, dimensioned by nmzg*nmx*nmy*npat.                               !
      !------------------------------------------------------------------------------------!
      npts = nmzg * nmx * nmy * nmpat

      if (associated(leaf%soil_water))                                                     &
         call vtables2(leaf%soil_water,leafm%soil_water,ng,npts,imean                      &
                      ,'SOIL_WATER :4:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%soil_energy))                                                    &
         call vtables2(leaf%soil_energy,leafm%soil_energy,ng,npts,imean                    &
                      ,'SOIL_ENERGY :4:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%soil_text))                                                      &
         call vtables2(leaf%soil_text,leafm%soil_text,ng,npts,imean                        &
                      ,'SOIL_TEXT :4:hist:anal:mpti:mpt3'//trim(str_recycle))


      !------------------------------------------------------------------------------------!
      !     4-D variables, dimensioned by nmzs*nmx*nmy*nmpat.                              !
      !------------------------------------------------------------------------------------!
      npts = nmzs * nmx * nmy * nmpat

      if (associated(leaf%sfcwater_mass))                                                  &
         call vtables2(leaf%sfcwater_mass,leafm%sfcwater_mass                              &
                      ,ng,npts,imean                                                       &
                      ,'SFCWATER_MASS :5:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%sfcwater_energy))                                                &
         call vtables2(leaf%sfcwater_energy,leafm%sfcwater_energy                          &
                      ,ng,npts,imean                                                       &
                      ,'SFCWATER_ENERGY :5:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%sfcwater_depth))                                                 &
         call vtables2(leaf%sfcwater_depth,leafm%sfcwater_depth                            &
                      ,ng,npts,imean                                                       &
                      ,'SFCWATER_DEPTH :5:hist:anal:mpti:mpt3'//trim(str_recycle))

      !------------------------------------------------------------------------------------!
      !     3-D variables, dimensioned by nmx*nmy*nmpat.                                   !
      !------------------------------------------------------------------------------------!
      npts = nmx * nmy * nmpat

      if (associated(leaf%soil_rough))                                                     &
         call vtables2(leaf%soil_rough,leafm%soil_rough,ng,npts,imean                      &
                      ,'SOIL_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%soil_color))                                                     &
         call vtables2(leaf%soil_color,leafm%soil_color,ng,npts,imean                      &
                      ,'SOIL_COLOR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%sfcwater_nlev))                                                  &
         call vtables2(leaf%sfcwater_nlev,leafm%sfcwater_nlev,ng,npts,imean                &
                      ,'SFCWATER_NLEV :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ground_rsat))                                                    &
         call vtables2(leaf%ground_rsat,leafm%ground_rsat,ng,npts,imean                    &
                      ,'GROUND_RSAT :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ground_rvap))                                                    &
         call vtables2(leaf%ground_rvap,leafm%ground_rvap,ng,npts,imean                    &
                      ,'GROUND_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ground_temp))                                                    &
         call vtables2(leaf%ground_temp,leafm%ground_temp,ng,npts,imean                    &
                      ,'GROUND_TEMP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ground_fliq))                                                    &
         call vtables2(leaf%ground_fliq,leafm%ground_fliq,ng,npts,imean                    &
                      ,'GROUND_FLIQ :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_fracarea))                                                   &
         call vtables2(leaf%veg_fracarea,leafm%veg_fracarea,ng,npts,imean                  &
                      ,'VEG_FRACAREA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_lai))                                                        &
         call vtables2(leaf%veg_lai,leafm%veg_lai,ng,npts,imean                            &
                      ,'VEG_LAI :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_agb))                                                        &
         call vtables2(leaf%veg_agb,leafm%veg_agb,ng,npts,imean                            &
                      ,'VEG_AGB :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_rough))                                                      &
         call vtables2(leaf%veg_rough,leafm%veg_rough,ng,npts,imean                        &
                      ,'VEG_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_height))                                                     &
         call vtables2(leaf%veg_height,leafm%veg_height,ng,npts,imean                      &
                      ,'VEG_HEIGHT :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_displace))                                                   &
         call vtables2(leaf%veg_displace,leafm%veg_height,ng,npts,imean                    &
                      ,'VEG_DISPLACE :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_albedo))                                                     &
         call vtables2(leaf%veg_albedo,leafm%veg_albedo,ng,npts,imean                      &
                      ,'VEG_ALBEDO :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_tai))                                                        &
         call vtables2(leaf%veg_tai,leafm%veg_tai,ng,npts,imean                            &
                      ,'VEG_TAI :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_water))                                                      &
         call vtables2(leaf%veg_water,leafm%veg_water,ng,npts,imean                        &
                      ,'VEG_WATER :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_hcap))                                                       &
         call vtables2(leaf%veg_hcap,leafm%veg_hcap,ng,npts,imean                          &
                      ,'VEG_HCAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_energy))                                                     &
         call vtables2(leaf%veg_energy,leafm%veg_energy,ng,npts,imean                      &
                      ,'VEG_ENERGY :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_ndvip))                                                      &
         call vtables2(leaf%veg_ndvip,leafm%veg_ndvip,ng,npts,imean                        &
                      ,'VEG_NDVIP :6:hist:mpti:mpt3')

      if (associated(leaf%veg_ndvic))                                                      &
         call vtables2(leaf%veg_ndvic,leafm%veg_ndvic,ng,npts,imean                        &
                      ,'VEG_NDVIC :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_ndvif))                                                      &
         call vtables2(leaf%veg_ndvif,leafm%veg_ndvif,ng,npts,imean                        &
                      ,'VEG_NDVIF :6:hist:mpti:mpt3')

      if (associated(leaf%leaf_class))                                                     &
         call vtables2(leaf%leaf_class,leafm%leaf_class,ng,npts,imean                      &
                      ,'LEAF_CLASS :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%stom_condct))                                                    &
         call vtables2(leaf%stom_condct,leafm%stom_condct,ng,npts,imean                    &
                      ,'STOM_CONDCT :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_rvap))                                                       &
         call vtables2(leaf%can_rvap,leafm%can_rvap,ng,npts,imean                          &
                      ,'CAN_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_theta))                                                      &
         call vtables2(leaf%can_theta,leafm%can_theta,ng,npts,imean                        &
                      ,'CAN_THETA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_theiv))                                                      &
         call vtables2(leaf%can_theiv,leafm%can_theta,ng,npts,imean                        &
                      ,'CAN_THEIV :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_vpdef))                                                      &
         call vtables2(leaf%can_vpdef,leafm%can_vpdef,ng,npts,imean                        &
                      ,'CAN_VPDEF :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_prss))                                                       &
         call vtables2(leaf%can_prss,leafm%can_prss,ng,npts,imean                          &
                      ,'CAN_PRSS :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_co2))                                                        &
         call vtables2(leaf%can_co2,leafm%can_co2,ng,npts,imean                            &
                      ,'CAN_CO2 :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ustar))                                                          &
         call vtables2(leaf%ustar,leafm%ustar,ng,npts,imean                                &
                      ,'USTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%tstar))                                                          &
         call vtables2(leaf%tstar,leafm%tstar,ng,npts,imean                                &
                      ,'TSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%estar))                                                          &
         call vtables2(leaf%estar,leafm%estar,ng,npts,imean                                &
                      ,'ESTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%rstar))                                                          &
         call vtables2(leaf%rstar,leafm%rstar,ng,npts,imean                                &
                      ,'RSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%cstar))                                                          &
         call vtables2(leaf%cstar,leafm%cstar,ng,npts,imean                                &
                      ,'CSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%zeta))                                                           &
         call vtables2(leaf%zeta,leafm%zeta,ng,npts,imean                                  &
                      ,'ZETA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ribulk))                                                         &
         call vtables2(leaf%ribulk,leafm%ribulk,ng,npts,imean                              &
                      ,'RIBULK :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    
      if (associated(leaf%patch_area))                                                     &
         call vtables2(leaf%patch_area,leafm%patch_area,ng,npts,imean                      &
                      ,'PATCH_AREA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%patch_rough))                                                    &
         call vtables2(leaf%patch_rough,leafm%patch_rough,ng,npts,imean                    &
                      ,'PATCH_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%patch_wetind))                                                   &
         call vtables2(leaf%patch_wetind,leafm%patch_wetind,ng,npts,imean                  &
                      ,'PATCH_WETIND :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%gpp))                                                            &
         call vtables2(leaf%gpp,leafm%gpp,ng,npts,imean                                    &
                      ,'GPP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%resphet ))                                                       &
         call vtables2(leaf%resphet,leafm%resphet,ng,npts,imean                            &
                      ,'RESPHET :6:hist:anal:mpti:mpt3')

      if (associated(leaf%plresp))                                                         &
         call vtables2(leaf%plresp,leafm%plresp,ng,npts,imean                              &
                      ,'PLRESP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%evap_gc))                                                        &
         call vtables2(leaf%evap_gc,leafm%evap_gc,ng,npts,imean                            &
                      ,'EVAP_GC :6:hist:anal:mpti:mpt3')

      if (associated(leaf%evap_vc))                                                        &
         call vtables2(leaf%evap_vc,leafm%evap_vc,ng,npts,imean                            &
                      ,'EVAP_VC :6:hist:anal:mpti:mpt3')

      if (associated(leaf%transp))                                                         &
         call vtables2(leaf%transp,leafm%transp,ng,npts,imean                              &
                      ,'TRANSP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%sensible_gc))                                                    &
         call vtables2(leaf%sensible_gc,leafm%sensible_gc,ng,npts,imean                    &
                      ,'SENSIBLE_GC :6:hist:anal:mpti:mpt3')

      if (associated(leaf%sensible_vc))                                                    &
         call vtables2(leaf%sensible_vc,leafm%sensible_vc,ng,npts,imean                    &
                      ,'SENSIBLE_VC :6:hist:anal:mpti:mpt3')

      if (associated(leaf%psibar_10d))                                                     &
         call vtables2(leaf%psibar_10d,leafm%psibar_10d,ng,npts,imean                      &
                      ,'PSIBAR_10D :6:hist:anal:mpti:mpt3')

      if (associated(leaf%R_aer))                                                          &
         call vtables2(leaf%R_aer,leafm%R_aer,ng,npts,imean                                &
                      ,'R_AER :6:hist:mpti:mpt3')

      if (associated(leaf%G_URBAN))                                                        &
         call vtables2(leaf%G_URBAN,leafm%G_URBAN,ng,npts,imean                            &
                      ,'G_URBAN :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%rshort_gnd))                                                     &
         call vtables2(leaf%rshort_gnd,leafm%rshort_gnd,ng,npts,imean                      &
                      ,'RSHORT_GND :6:hist:anal:mpti:mpt3')

      if (associated(leaf%rlong_gnd))                                                      &
         call vtables2(leaf%rlong_gnd,leafm%rlong_gnd,ng,npts,imean                        &
                      ,'RLONG_GND :6:hist:anal:mpti:mpt3')


      !------------------------------------------------------------------------------------!
      !     2-D variables, dimensioned by nx*ny.                                           !
      !------------------------------------------------------------------------------------!
      npts=nmx*nmy

      if (associated(leaf%snow_mass))                                                      &
         call vtables2(leaf%snow_mass,leafm%snow_mass,ng,npts,imean                        &
                      ,'SNOW_MASS :2:mpti:mpt3')

      if (associated(leaf%snow_depth))                                                     &
         call vtables2(leaf%snow_depth,leafm%snow_depth,ng,npts,imean                      &
                      ,'SNOW_DEPTH :2:mpti:mpt3')

      if (associated(leaf%seatp))                                                          &
         call vtables2(leaf%seatp,leafm%seatp,ng,npts,imean                                &
                      ,'SEATP :2:mpti:mpt3')

      if (associated(leaf%seatf))                                                          &
         call vtables2(leaf%seatf,leafm%seatf,ng,npts,imean                                &
                      ,'SEATF :2:mpti:mpt3')

      return
   end subroutine filltab_leaf
   !=======================================================================================!
   !=======================================================================================!
end module mem_leaf
!==========================================================================================!
!==========================================================================================!
