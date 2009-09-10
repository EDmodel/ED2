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
      !   Soil properties, dimensioned by (nzg,nxp,nyp,npatch), except for roughness which !
      ! is dimensioned by (nxp,nyp,npatch).                                                !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:,:), pointer :: soil_water  & ! Soil moisture        [    m³/m³]
                                         , soil_energy & ! Internal energy      [     J/m³]
                                         , soil_text   ! ! Soil texture class   [      ---]
      real, dimension(  :,:,:), pointer :: soil_rough  ! ! Soil roughness       [        m]

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
                                       , ground_rvap ! ! Vapour mixing ratio    [    kg/kg]

      !------------------------------------------------------------------------------------!
      !     Vegetation properties, dimensioned by (nxp,nyp,npatch).                        !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: veg_fracarea & ! Fractional area       [      ---]
                                       , veg_lai      & ! Leaf area index       [    m²/m²]
                                       , veg_rough    & ! Roughness length      [        m]
                                       , veg_height   & ! Height                [        m]
                                       , veg_albedo   & ! Albedo                [      ---]
                                       , veg_tai      & ! Tree area index       [    m²/m²]
                                       , veg_water    & ! Leaf surface water    [    kg/m²]
                                       , veg_hcap     & ! Heat capacity         [   J/m²/K]
                                       , veg_energy   & ! Internal energy       [     J/m²]
                                       , veg_ndvip    & ! Past NDVI             [      ---]
                                       , veg_ndvic    & ! Current NDVI          [      ---]
                                       , veg_ndvif    & ! Future NDVI           [      ---]
                                       , leaf_class   & ! Vegetation class      [      ---]
                                       , stom_resist  ! ! Stomatal resistance   [      ???]


      !------------------------------------------------------------------------------------!
      !     Canopy air properties, dimensioned by (nxp,nyp,npatch).                        !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: can_rvap     & ! Vapour Mixing ratio   [    kg/kg]
                                       , can_theta    & ! Potential temperature [        K]
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
      !     Patch structural properties, dimensioned by (nxp,nyp,npatch).                  !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: patch_area   & ! Patch relative area   [      ---]
                                       , patch_rough  & ! Roughness length      [        m]
                                       , patch_wetind ! ! Wetness index         [      ???]

      !------------------------------------------------------------------------------------!
      !     Surface fluxes, dimensioned by (nxp,nyp,npatch).                               !
      !------------------------------------------------------------------------------------!
      real, dimension(:,:,:), pointer :: gpp      & ! Gross primary production  [µmol/m²/s]
                                       , resphet  & ! Heterotrophic respiration [µmol/m²/s]
                                       , plresp   & ! Plant respiration         [µmol/m²/s]
                                       , evap     & ! Evaporation               [     W/m²]
                                       , transp   & ! Transpiration             [     W/m²]
                                       , sensible ! ! Sensible heat flux        [     W/m²]

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
   integer                 :: nslcon  ! Soil texture if constant for entire domain
   integer                 :: nvgcon  ! Vegetation class if constant for entire domain
   integer                 :: nvegpat ! Number of vegetation types
   integer                 :: isfcl   ! Surface model
   real, dimension(nzgmax) :: stgoff  ! Initial soil temperature offset
   real, dimension(nzgmax) :: slmstr  ! Initial soil moisture if constant for entire domain
   real, dimension(nzgmax) :: slz     ! Soil levels
   real                    :: zrough  ! Roughness if constant for entire domain
   real                    :: pctlcon ! Vegetation fraction if constant for entire domain
   real                    :: albedo  ! Albedo if constant for entire domain
   real                    :: drtcon  ! ???
   real                    :: dthcon  ! ???
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

      allocate (leaf%sfcwater_mass    (nzs,nx,ny,np))
      allocate (leaf%sfcwater_energy  (nzs,nx,ny,np))
      allocate (leaf%sfcwater_depth   (nzs,nx,ny,np))
      allocate (leaf%sfcwater_nlev    (    nx,ny,np))

      allocate (leaf%ground_rsat      (    nx,ny,np))
      allocate (leaf%ground_rvap      (    nx,ny,np))

      allocate (leaf%veg_fracarea     (    nx,ny,np))
      allocate (leaf%veg_lai          (    nx,ny,np))
      allocate (leaf%veg_rough        (    nx,ny,np))
      allocate (leaf%veg_height       (    nx,ny,np))
      allocate (leaf%veg_albedo       (    nx,ny,np))
      allocate (leaf%veg_tai          (    nx,ny,np))
      allocate (leaf%veg_water        (    nx,ny,np))
      allocate (leaf%veg_hcap         (    nx,ny,np))
      allocate (leaf%veg_energy       (    nx,ny,np))
      allocate (leaf%veg_ndvip        (    nx,ny,np))
      allocate (leaf%veg_ndvic        (    nx,ny,np))
      allocate (leaf%veg_ndvif        (    nx,ny,np))
      allocate (leaf%leaf_class       (    nx,ny,np))
      allocate (leaf%stom_resist      (    nx,ny,np))

      allocate (leaf%can_rvap         (    nx,ny,np))
      allocate (leaf%can_theta        (    nx,ny,np))
      allocate (leaf%can_prss         (    nx,ny,np))
      allocate (leaf%can_co2          (    nx,ny,np))

      allocate (leaf%ustar            (    nx,ny,np))
      allocate (leaf%tstar            (    nx,ny,np))
      allocate (leaf%estar            (    nx,ny,np))
      allocate (leaf%rstar            (    nx,ny,np))
      allocate (leaf%cstar            (    nx,ny,np))

      allocate (leaf%patch_area       (    nx,ny,np))
      allocate (leaf%patch_rough      (    nx,ny,np))
      allocate (leaf%patch_wetind     (    nx,ny,np))


      allocate (leaf%gpp              (    nx,ny,np))
      allocate (leaf%resphet          (    nx,ny,np))
      allocate (leaf%plresp           (    nx,ny,np))
      allocate (leaf%evap             (    nx,ny,np))
      allocate (leaf%transp           (    nx,ny,np))
      allocate (leaf%sensible         (    nx,ny,np))

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


      if (associated(leaf%soil_water       ))  nullify(leaf%soil_water       )
      if (associated(leaf%soil_energy      ))  nullify(leaf%soil_energy      )
      if (associated(leaf%soil_text        ))  nullify(leaf%soil_text        )
      if (associated(leaf%soil_rough       ))  nullify(leaf%soil_rough       )

      if (associated(leaf%sfcwater_mass    ))  nullify(leaf%sfcwater_mass    )
      if (associated(leaf%sfcwater_energy  ))  nullify(leaf%sfcwater_energy  )
      if (associated(leaf%sfcwater_depth   ))  nullify(leaf%sfcwater_depth   )
      if (associated(leaf%sfcwater_nlev    ))  nullify(leaf%sfcwater_nlev    )

      if (associated(leaf%ground_rsat      ))  nullify(leaf%ground_rsat      )
      if (associated(leaf%ground_rvap      ))  nullify(leaf%ground_rvap      )

      if (associated(leaf%veg_fracarea     ))  nullify(leaf%veg_fracarea     )
      if (associated(leaf%veg_lai          ))  nullify(leaf%veg_lai          )
      if (associated(leaf%veg_rough        ))  nullify(leaf%veg_rough        )
      if (associated(leaf%veg_height       ))  nullify(leaf%veg_height       )
      if (associated(leaf%veg_albedo       ))  nullify(leaf%veg_albedo       )
      if (associated(leaf%veg_tai          ))  nullify(leaf%veg_tai          )
      if (associated(leaf%veg_water        ))  nullify(leaf%veg_water        )
      if (associated(leaf%veg_hcap         ))  nullify(leaf%veg_hcap         )
      if (associated(leaf%veg_energy       ))  nullify(leaf%veg_energy       )
      if (associated(leaf%veg_ndvip        ))  nullify(leaf%veg_ndvip        )
      if (associated(leaf%veg_ndvic        ))  nullify(leaf%veg_ndvic        )
      if (associated(leaf%veg_ndvif        ))  nullify(leaf%veg_ndvif        )
      if (associated(leaf%leaf_class       ))  nullify(leaf%leaf_class       )
      if (associated(leaf%stom_resist      ))  nullify(leaf%stom_resist      )

      if (associated(leaf%can_rvap         ))  nullify(leaf%can_rvap         )
      if (associated(leaf%can_theta        ))  nullify(leaf%can_theta        )
      if (associated(leaf%can_prss         ))  nullify(leaf%can_prss         )
      if (associated(leaf%can_co2          ))  nullify(leaf%can_co2          )

      if (associated(leaf%ustar            ))  nullify(leaf%ustar            )
      if (associated(leaf%tstar            ))  nullify(leaf%tstar            )
      if (associated(leaf%estar            ))  nullify(leaf%estar            )
      if (associated(leaf%rstar            ))  nullify(leaf%rstar            )
      if (associated(leaf%cstar            ))  nullify(leaf%cstar            )

      if (associated(leaf%patch_area       ))  nullify(leaf%patch_area       )
      if (associated(leaf%patch_rough      ))  nullify(leaf%patch_rough      )
      if (associated(leaf%patch_wetind     ))  nullify(leaf%patch_wetind     )


      if (associated(leaf%gpp              ))  nullify(leaf%gpp              )
      if (associated(leaf%resphet          ))  nullify(leaf%resphet          )
      if (associated(leaf%plresp           ))  nullify(leaf%plresp           )
      if (associated(leaf%evap             ))  nullify(leaf%evap             )
      if (associated(leaf%transp           ))  nullify(leaf%transp           )
      if (associated(leaf%sensible         ))  nullify(leaf%sensible         )

      if (associated(leaf%R_aer            ))  nullify(leaf%R_aer            )
      if (associated(leaf%G_URBAN          ))  nullify(leaf%G_URBAN          )

      if (associated(leaf%snow_mass        ))  nullify(leaf%snow_mass        )
      if (associated(leaf%snow_depth       ))  nullify(leaf%snow_depth       )
      if (associated(leaf%seatp            ))  nullify(leaf%seatp            )
      if (associated(leaf%seatf            ))  nullify(leaf%seatf            )

      return
   end subroutine nullify_leaf
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

      if (associated(leaf%sfcwater_mass    ))  deallocate(leaf%sfcwater_mass    )
      if (associated(leaf%sfcwater_energy  ))  deallocate(leaf%sfcwater_energy  )
      if (associated(leaf%sfcwater_depth   ))  deallocate(leaf%sfcwater_depth   )
      if (associated(leaf%sfcwater_nlev    ))  deallocate(leaf%sfcwater_nlev    )

      if (associated(leaf%ground_rsat      ))  deallocate(leaf%ground_rsat      )
      if (associated(leaf%ground_rvap      ))  deallocate(leaf%ground_rvap      )

      if (associated(leaf%veg_fracarea     ))  deallocate(leaf%veg_fracarea     )
      if (associated(leaf%veg_lai          ))  deallocate(leaf%veg_lai          )
      if (associated(leaf%veg_rough        ))  deallocate(leaf%veg_rough        )
      if (associated(leaf%veg_height       ))  deallocate(leaf%veg_height       )
      if (associated(leaf%veg_albedo       ))  deallocate(leaf%veg_albedo       )
      if (associated(leaf%veg_tai          ))  deallocate(leaf%veg_tai          )
      if (associated(leaf%veg_water        ))  deallocate(leaf%veg_water        )
      if (associated(leaf%veg_hcap         ))  deallocate(leaf%veg_hcap         )
      if (associated(leaf%veg_energy       ))  deallocate(leaf%veg_energy       )
      if (associated(leaf%veg_ndvip        ))  deallocate(leaf%veg_ndvip        )
      if (associated(leaf%veg_ndvic        ))  deallocate(leaf%veg_ndvic        )
      if (associated(leaf%veg_ndvif        ))  deallocate(leaf%veg_ndvif        )
      if (associated(leaf%leaf_class       ))  deallocate(leaf%leaf_class       )
      if (associated(leaf%stom_resist      ))  deallocate(leaf%stom_resist      )

      if (associated(leaf%can_rvap         ))  deallocate(leaf%can_rvap         )
      if (associated(leaf%can_theta        ))  deallocate(leaf%can_theta        )
      if (associated(leaf%can_prss         ))  deallocate(leaf%can_prss         )
      if (associated(leaf%can_co2          ))  deallocate(leaf%can_co2          )

      if (associated(leaf%ustar            ))  deallocate(leaf%ustar            )
      if (associated(leaf%tstar            ))  deallocate(leaf%tstar            )
      if (associated(leaf%estar            ))  deallocate(leaf%estar            )
      if (associated(leaf%rstar            ))  deallocate(leaf%rstar            )
      if (associated(leaf%cstar            ))  deallocate(leaf%cstar            )

      if (associated(leaf%patch_area       ))  deallocate(leaf%patch_area       )
      if (associated(leaf%patch_rough      ))  deallocate(leaf%patch_rough      )
      if (associated(leaf%patch_wetind     ))  deallocate(leaf%patch_wetind     )


      if (associated(leaf%gpp              ))  deallocate(leaf%gpp              )
      if (associated(leaf%resphet          ))  deallocate(leaf%resphet          )
      if (associated(leaf%plresp           ))  deallocate(leaf%plresp           )
      if (associated(leaf%evap             ))  deallocate(leaf%evap             )
      if (associated(leaf%transp           ))  deallocate(leaf%transp           )
      if (associated(leaf%sensible         ))  deallocate(leaf%sensible         )

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
   subroutine filltab_leaf(leaf,leafm,imean,nz,nx,ny,nzg,nzs,np,ng)
      use var_tables
      use io_params, only: ipastin ! INTENT(IN)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type (leaf_vars), intent(inout) :: leaf
      type (leaf_vars), intent(inout) :: leafm
      integer         , intent(in)    :: imean
      integer         , intent(in)    :: nz
      integer         , intent(in)    :: nx
      integer         , intent(in)    :: ny
      integer         , intent(in)    :: nzg
      integer         , intent(in)    :: nzs
      integer         , intent(in)    :: np
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
      !     4-D variables, dimensioned by nzg*nx*ny*np.                                    !
      !------------------------------------------------------------------------------------!
      npts = nzg * nx * ny * np

      if (associated(leaf%soil_water))                                                     &
         call vtables2(leaf%soil_water(1,1,1,1),leafm%soil_water(1,1,1,1),ng,npts,imean    &
                      ,'SOIL_WATER :4:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%soil_energy))                                                    &
         call vtables2(leaf%soil_energy(1,1,1,1),leafm%soil_energy(1,1,1,1),ng,npts,imean  &
                      ,'SOIL_ENERGY :4:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%soil_text))                                                      &
         call vtables2(leaf%soil_text(1,1,1,1),leafm%soil_text(1,1,1,1),ng,npts,imean      &
                      ,'SOIL_TEXT :4:hist:anal:mpti:mpt3'//trim(str_recycle))


      !------------------------------------------------------------------------------------!
      !     4-D variables, dimensioned by nzs*nx*ny*np.                                    !
      !------------------------------------------------------------------------------------!
      npts = nzs * nx * ny * np

      if (associated(leaf%sfcwater_mass))                                                  &
         call vtables2(leaf%sfcwater_mass(1,1,1,1),leafm%sfcwater_mass(1,1,1,1)            &
                      ,ng,npts,imean                                                       &
                      ,'SFCWATER_MASS :5:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%sfcwater_energy))                                                &
         call vtables2(leaf%sfcwater_energy(1,1,1,1),leafm%sfcwater_energy(1,1,1,1)        &
                      ,ng,npts,imean                                                       &
                      ,'SFCWATER_ENERGY :5:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%sfcwater_depth))                                                 &
         call vtables2(leaf%sfcwater_depth(1,1,1,1),leafm%sfcwater_depth(1,1,1,1)          &
                      ,ng,npts,imean                                                       &
                      ,'SFCWATER_DEPTH :5:hist:anal:mpti:mpt3'//trim(str_recycle))

      !------------------------------------------------------------------------------------!
      !     3-D variables, dimensioned by nx*ny*np.                                        !
      !------------------------------------------------------------------------------------!
      npts = nx * ny * np

      if (associated(leaf%soil_rough))                                                     &
         call vtables2(leaf%soil_rough(1,1,1),leafm%soil_rough(1,1,1),ng,npts,imean        &
                      ,'SOIL_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%sfcwater_nlev))                                                  &
         call vtables2(leaf%sfcwater_nlev(1,1,1),leafm%sfcwater_nlev(1,1,1),ng,npts,imean  &
                      ,'SFCWATER_NLEV :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ground_rsat))                                                    &
         call vtables2(leaf%ground_rsat(1,1,1),leafm%ground_rsat(1,1,1),ng,npts,imean      &
                      ,'GROUND_RSAT :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ground_rvap))                                                    &
         call vtables2(leaf%ground_rvap(1,1,1),leafm%ground_rvap(1,1,1),ng,npts,imean      &
                      ,'GROUND_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_fracarea))                                                   &
         call vtables2(leaf%veg_fracarea(1,1,1),leafm%veg_fracarea(1,1,1),ng,npts,imean    &
                      ,'VEG_FRACAREA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_lai))                                                        &
         call vtables2(leaf%veg_lai(1,1,1),leafm%veg_lai(1,1,1),ng,npts,imean              &
                      ,'VEG_LAI :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_rough))                                                      &
         call vtables2(leaf%veg_rough(1,1,1),leafm%veg_rough(1,1,1),ng,npts,imean          &
                      ,'VEG_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_height))                                                     &
         call vtables2(leaf%veg_height(1,1,1),leafm%veg_height(1,1,1),ng,npts,imean        &
                      ,'VEG_HEIGHT :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_albedo))                                                     &
         call vtables2(leaf%veg_albedo(1,1,1),leafm%veg_albedo(1,1,1),ng,npts,imean        &
                      ,'VEG_ALBEDO :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_tai))                                                        &
         call vtables2(leaf%veg_tai(1,1,1),leafm%veg_tai(1,1,1),ng,npts,imean              &
                      ,'VEG_TAI :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_water))                                                      &
         call vtables2(leaf%veg_water(1,1,1),leafm%veg_water(1,1,1),ng,npts,imean          &
                      ,'VEG_WATER :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_hcap))                                                       &
         call vtables2(leaf%veg_hcap(1,1,1),leafm%veg_hcap(1,1,1),ng,npts,imean            &
                      ,'VEG_HCAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_energy))                                                     &
         call vtables2(leaf%veg_energy(1,1,1),leafm%veg_energy(1,1,1),ng,npts,imean        &
                      ,'VEG_ENERGY :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_ndvip))                                                      &
         call vtables2(leaf%veg_ndvip(1,1,1),leafm%veg_ndvip(1,1,1),ng,npts,imean          &
                      ,'VEG_NDVIP :6:hist:mpti')

      if (associated(leaf%veg_ndvic))                                                      &
         call vtables2(leaf%veg_ndvic(1,1,1),leafm%veg_ndvic(1,1,1),ng,npts,imean          &
                      ,'VEG_NDVIC :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%veg_ndvif))                                                      &
         call vtables2(leaf%veg_ndvif(1,1,1),leafm%veg_ndvif(1,1,1),ng,npts,imean          &
                      ,'VEG_NDVIF :6:hist:mpti')

      if (associated(leaf%leaf_class))                                                     &
         call vtables2(leaf%leaf_class(1,1,1),leafm%leaf_class(1,1,1),ng,npts,imean        &
                      ,'LEAF_CLASS :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%stom_resist))                                                    &
         call vtables2(leaf%stom_resist(1,1,1),leafm%stom_resist(1,1,1),ng,npts,imean      &
                      ,'STOM_RESIST :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_rvap))                                                       &
         call vtables2(leaf%can_rvap(1,1,1),leafm%can_rvap(1,1,1),ng,npts,imean            &
                      ,'CAN_RVAP :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_theta))                                                      &
         call vtables2(leaf%can_theta(1,1,1),leafm%can_theta(1,1,1),ng,npts,imean          &
                      ,'CAN_THETA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_prss))                                                       &
         call vtables2(leaf%can_prss(1,1,1),leafm%can_prss(1,1,1),ng,npts,imean            &
                      ,'CAN_PRSS :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%can_co2))                                                        &
         call vtables2(leaf%can_co2(1,1,1),leafm%can_co2(1,1,1),ng,npts,imean              &
                      ,'CAN_CO2 :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%ustar))                                                          &
         call vtables2(leaf%ustar(1,1,1),leafm%ustar(1,1,1),ng,npts,imean                  &
                      ,'USTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%tstar))                                                          &
         call vtables2(leaf%tstar(1,1,1),leafm%tstar(1,1,1),ng,npts,imean                  &
                      ,'TSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%estar))                                                          &
         call vtables2(leaf%estar(1,1,1),leafm%estar(1,1,1),ng,npts,imean                  &
                      ,'ESTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%rstar))                                                          &
         call vtables2(leaf%rstar(1,1,1),leafm%rstar(1,1,1),ng,npts,imean                  &
                      ,'RSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%cstar))                                                          &
         call vtables2(leaf%cstar(1,1,1),leafm%cstar(1,1,1),ng,npts,imean                  &
                      ,'CSTAR :6:hist:anal:mpti:mpt3'//trim(str_recycle))
    
      if (associated(leaf%patch_area))                                                     &
         call vtables2(leaf%patch_area(1,1,1),leafm%patch_area(1,1,1),ng,npts,imean        &
                      ,'PATCH_AREA :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%patch_rough))                                                    &
         call vtables2(leaf%patch_rough(1,1,1),leafm%patch_rough(1,1,1),ng,npts,imean      &
                      ,'PATCH_ROUGH :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%patch_wetind))                                                   &
         call vtables2(leaf%patch_wetind(1,1,1),leafm%patch_wetind(1,1,1),ng,npts,imean    &
                      ,'PATCH_WETIND :6:hist:anal:mpti:mpt3'//trim(str_recycle))

      if (associated(leaf%gpp))                                                            &
         call vtables2(leaf%gpp(1,1,1),leafm%gpp(1,1,1),ng,npts,imean                      &
                      ,'GPP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%resphet ))                                                       &
         call vtables2(leaf%resphet(1,1,1),leafm%resphet(1,1,1),ng,npts,imean              &
                      ,'RESPHET :6:hist:anal:mpti:mpt3')

      if (associated(leaf%plresp))                                                         &
         call vtables2(leaf%plresp(1,1,1),leafm%plresp(1,1,1),ng,npts,imean                &
                      ,'PLRESP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%evap))                                                           &
         call vtables2(leaf%evap(1,1,1),leafm%evap(1,1,1),ng,npts,imean                    &
                      ,'EVAP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%transp))                                                         &
         call vtables2(leaf%transp(1,1,1),leafm%transp(1,1,1),ng,npts,imean                &
                      ,'TRANSP :6:hist:anal:mpti:mpt3')

      if (associated(leaf%sensible))                                                       &
         call vtables2(leaf%sensible(1,1,1),leafm%sensible(1,1,1),ng,npts,imean            &
                      ,'SENSIBLE :6:hist:anal:mpti:mpt3')

      if (associated(leaf%R_aer))                                                          &
         call vtables2(leaf%R_aer(1,1,1),leafm%R_aer(1,1,1),ng,npts,imean                  &
                      ,'R_AER :6:hist:mpti')

      if (associated(leaf%G_URBAN))                                                        &
         call vtables2(leaf%G_URBAN(1,1,1),leafm%G_URBAN(1,1,1),ng,npts,imean              &
                      ,'G_URBAN :6:hist:anal:mpti:mpt3'//trim(str_recycle))


      !------------------------------------------------------------------------------------!
      !     2-D variables, dimensioned by nx*ny.                                           !
      !------------------------------------------------------------------------------------!
      npts=nx*ny

      if (associated(leaf%snow_mass))                                                      &
         call vtables2(leaf%snow_mass(1,1),leafm%snow_mass(1,1),ng,npts,imean              &
                      ,'SNOW_MASS :2:mpti')

      if (associated(leaf%snow_depth))                                                     &
         call vtables2(leaf%snow_depth(1,1),leafm%snow_depth(1,1),ng,npts,imean            &
                      ,'SNOW_DEPTH :2:mpti')

      if (associated(leaf%seatp))                                                          &
         call vtables2(leaf%seatp(1,1),leafm%seatp(1,1),ng,npts,imean                      &
                      ,'SEATP :2:mpti')

      if (associated(leaf%seatf))                                                          &
         call vtables2(leaf%seatf(1,1),leafm%seatf(1,1),ng,npts,imean                      &
                      ,'SEATF :2:mpti')

      return
   end subroutine filltab_leaf
   !=======================================================================================!
   !=======================================================================================!
end module mem_leaf
!==========================================================================================!
!==========================================================================================!
