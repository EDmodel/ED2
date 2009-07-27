!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This module contains some variables used in LEAF-3.                                   !
!------------------------------------------------------------------------------------------!
module leaf_coms

   use grid_dims
   use rconstants, only: g

   integer :: niter_leaf   & ! number of leaf timesteps in model long timestep
            , niter_can    ! ! number of canopy timesteps in leaf timestep

   real    :: dtll         & ! leaf timestep
            , dtll_factor  & ! leaf timestep factor (leaf timestep / model timestep)
            , dtlc         & ! canopy timestep
            , dtlc_factor  & ! canopy timestep factor (canopy timestep / leaf timestep)
            , dtllohcc     & ! leaf timestep / canopy heat capacity
            , dtllowcc     & ! leaf timestep / canopy vapor capacity
            , dtlloccc     & ! leaf timestep / canopy CO2 capacity
            , dtlcohcc     & ! canopy timestep / canopy heat capacity
            , dtlcowcc     & ! canopy timestep / canopy vapor capacity
            , dtlcoccc     & ! canopy timestep / canopy CO2 capacity
            , dtlcohcv     & ! caonpy timestep / vegetation heat capacity
  
            , ups          & ! U velocity at top of surface layer [up(2,i,j)]    
            , vps          & ! V velocity at top of surface layer [vp(2,i,j)]
            , ths          & ! potential temperature at top of surface layer [theta(2,i,j)]
            , rvs          & ! vapor mixing ratio at top of surface layer [rv(2,i,j)]
            , rco2s        & ! CO2 mixing ratio at top of surface layer [rv(2,i,j)]
            , zts          & ! height at top of surface layer [zt(2)*rtgt(i,j)]
            , pis          & ! Exner function at surface
            , dens         & ! density at surface
            , prss         & ! pressure at surface
            , vels         & ! wind speed at top of surface layer
            , vels_pat     & ! vels with patch-dependent ubmin imposed as minimum
            , gzotheta     & ! (g*z/theta) at top of surface layer [for Richardson number]
            , pcpgl        & ! precip mass from cuparm and/or micphys in leaf timestep
            , qpcpgl       & ! precip energy from cuparm and/or micphys in leaf timestep
            , dpcpgl       & ! precip depth from cuparm and/or micphys in leaf timestep
            , pcpgc        & ! precip mass from cuparm and/or micphys in canopy timestep
            , qpcpgc       & ! precip energy from cuparm and/or micphys in canopy timestep
            
            , snowfac      & ! fraction of vegetation height covered by sfcwater
            , vf           & ! product of veg_fracarea and (1-snowfac)
            , thetacan     & ! canopy air potential temperature
            , transp       & ! transpiration flux [kg/m2/s]
            , timefac_sst  & ! time interpolation factor for SST
            
            , rb           & ! vegetation aerodynamic resistance
            , rd           & ! canopy to ground aerodynamic resistance
            , rdi          & ! canopy to ground aerodynamic conductance
            
            , rshort_g     & ! net SW radiation absorbed by grnd
            , rshort_v     & ! net SW radiation absorbed by veg
            , rshort_a     & ! net SW radiation reflected to atm by veg plus grnd
            
            , rlonga_v     & ! net atm LW radiation absorbed by veg
            , rlonga_gs    & ! net atm LW radiation absorbed by grnd OR snow
            , rlongv_gs    & ! net veg LW radiation absorbed by grnd OR snow
            , rlongv_a     & ! net veg LW radiation to atm
            , rlonggs_v    & ! net grnd OR snow LW radiation absorbed by veg
            , rlonggs_a    & ! net grnd OR snow LW radiation to atm
            , rlonga_a     & ! net atm LW radiation refl. to atm by veg plus grnd OR snow

            , hflxgc       & ! sensible heat from ground to canopy (J/m2)
            , wflxgc       & ! water vapor from ground to canopy (kg/m2)
            , cflxgc       & ! carbon from ground to canopy (umol/m2)
            , hflxvc       & ! sensible heat from vegetation to canopy (J/m2)
            , wflxvc       & ! water vapor from vegetation to canopy (kg/m2)
            , cflxvc       & ! carbon from vegetation to canopy (umol/m2)
            
            , wshed        & ! water shed from vegetation to ground (kg/m2)
            , qwshed       & ! energy from shed water (J/m2)
            , dewgnd       ! ! dew formation on ground (kg/m2)

   real, save, allocatable, dimension(:) ::  &
              dslz            & ! soil layer thickness at T point
            , dslzi           & ! (1. / soil layer thickness at T point)
            , dslzidt         & ! (dtll / soil layer thickness at T point)
            , slzt            & ! soil depth at T point
            , dslzt           & ! soil layer thickness at M point
            , dslzti          & ! (1. / soil layer thickness at M point)
            , dslztidt        & ! (dtll / soil layer thickness at M point)

            , rshort_s        & ! net SW radiation absorbed by snow layers (mzs)
            , sfcw_energy_ext & ! Extensive sfc. water internal energy (mzs)
            , tempk           & ! diagnosed temp (K) of soil and sfcwater
            , fracliq         & ! diagnosed liquid fraction (soil_water & sfcwater_mass)
      
            , hfluxgsc        & ! sensible heat flux between soil, sfcwater, canopy
            , psiplusz        & ! soil water potential plus geopotential [m]
            , half_soilair    & ! half of available airspace in soil [m]
            , rfactor         & ! soil, sfcwater thermal resistance    
            , wflux           & ! soil water flux [m]
            , soil_liq        & ! soil liquid water content [m]
            , qwflux          ! ! soil energy flux from water flux [J/m2] 
   !---------------------------------------------------------------------------------------!



   !----- Some dimensions -----------------------------------------------------------------!
   integer, parameter :: nstyp     = 12 ! # of soil types
   integer, parameter :: nvtyp     = 20 ! # of land use types
   integer, parameter :: nvtyp_teb =  1 ! # of TEB extra land use types (21 - Very urban).
   !---------------------------------------------------------------------------------------!



   !----- Soil properties -----------------------------------------------------------------!
   real, dimension(nstyp)           :: slden,slcpd,slbs,slcond,sfldcap,slcons,slmsts,slpots
   real, dimension(nstyp)           :: ssand,sclay,sorgan,sporo,soilcp,slfc,emisg,slcons00
   real, dimension(nstyp)           :: slcons0,fhydraul
   real, dimension(nzgmax,nstyp)    :: slcons1
   !---------------------------------------------------------------------------------------!



   !----- Plant properties ----------------------------------------------------------------!
   real   , dimension(nvtyp+nvtyp_teb)        :: albv_green,albv_brown,emisv,sr_max,tai_max
   real   , dimension(nvtyp+nvtyp_teb)        :: sai,veg_clump,veg_frac,veg_ht,glai_max
   real   , dimension(nvtyp+nvtyp_teb)        :: rcmin,dead_frac
   integer, dimension(nvtyp+nvtyp_teb)        :: kroot
   real   , dimension(nzgmax,nvtyp+nvtyp_teb) :: root
   !---------------------------------------------------------------------------------------!


   !----- Other variables -----------------------------------------------------------------!
   real                             :: cmin,corg,cwat,cair,cka,ckw
   !---------------------------------------------------------------------------------------!



   !----- Heat and vapour capacities ------------------------------------------------------!
   real, parameter :: hcapcan = 2.0e4   ! Canopy heat capacity
   real, parameter :: wcapcan = 2.0e1   ! Canopy vapour capacity
   real, parameter :: ccapcan = 6.904e2 ! Canopy CO2 capacity [wcapcan/molar mass of air]
   real, parameter :: hcapveg = 3.e4    ! Leaf heat capacity
   !---------------------------------------------------------------------------------------!



   !----- Roughness -----------------------------------------------------------------------!
   real, parameter :: z0fac_water    = .016 / g   ! Coefficient before ustar²
   real, parameter :: min_waterrough = .0001      ! Minimum water roughness height
   real, parameter :: snowrough      = .001       ! Snow roughness height
   !---------------------------------------------------------------------------------------!



   !----- Some constants to ensure the model good behaviour -------------------------------!
   real, parameter :: min_sfcwater_mass  = 1.e-6
   real, parameter :: min_sfcwater_depth = 1.e-9
   !---------------------------------------------------------------------------------------!



   !----- Parameters that used to be in LEAF-3 (leaftw) -----------------------------------!

   !-----------------------------------------------------------------------------
   ! parameters for new soil heat conductivity (8/17/00):  Move to sfcdata later
   real, dimension(nstyp), parameter :: soilcond0 = (/ 0.30, 0.30, 0.29, 0.27, 0.28, 0.28  &
                                                     , 0.26, 0.27, 0.27, 0.25, 0.25, 0.06 /)
   real, dimension(nstyp), parameter :: soilcond1 = (/ 4.80, 4.66, 4.27, 3.47, 3.63, 3.78  &
                                                     , 2.73, 3.23, 3.32, 2.58, 2.40, 0.46 /)
   real, dimension(nstyp), parameter :: soilcond2 = (/-2.70,-2.60,-2.31,-1.74,-1.85,-1.96  &
                                                     ,-1.20,-1.56,-1.63,-1.09,-0.96, 0.00 /)
   !-----------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------!


   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine allocates some specific scratch arrays used by LEAF-3.              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_leafcol(nzg,nzs)

      implicit none
      integer, intent(in) :: nzg,nzs

      !----- Allocate leaf column arrays. -------------------------------------------------!
      allocate (dslz            (nzg)      )
      allocate (dslzi           (nzg)      )
      allocate (dslzidt         (nzg)      )
      allocate (slzt            (nzg)      )
      allocate (dslzt           (nzg)      )
      allocate (dslzti          (nzg)      )
      allocate (dslztidt        (nzg)      )

      allocate (rshort_s        (nzs)      )
      allocate (sfcw_energy_ext (nzs)      )
      allocate (tempk           (nzg+nzs)  )
      allocate (fracliq         (nzg+nzs)  )

      allocate (hfluxgsc        (nzg+nzs+1))
      allocate (psiplusz        (nzg)      )
      allocate (half_soilair    (nzg)      )
      allocate (rfactor         (nzg+nzs)  )
      allocate (wflux           (nzg+1)    )
      allocate (qwflux          (nzg+1)    )
      allocate (soil_liq        (nzg)      )
               
      return
   end subroutine alloc_leafcol   
   !=======================================================================================!
   !=======================================================================================!
end module leaf_coms
!==========================================================================================!
!==========================================================================================!



