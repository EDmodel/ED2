!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module leaf_coms

use grid_dims

integer ::       &
    niter_leaf   & ! number of leaf timesteps in model long timestep
   ,niter_can      ! number of canopy timesteps in leaf timestep

real ::          &
    dtll         & ! leaf timestep
   ,dtll_factor  & ! leaf timestep factor (leaf timestep / model timestep)
   ,dtlc         & ! canopy timestep
   ,dtlc_factor  & ! canopy timestep factor (canopy timestep / leaf timestep)
   ,hcapcan      & ! canopy heat capacity
   ,wcapcan      & ! canopy vapor capacity
   ,hcapveg      & ! vegetation heat capacity
   ,dtllohcc     & ! leaf timestep / canopy heat capacity
   ,dtllowcc     & ! leaf timestep / canopy vapor capacity
   ,dtlcohcc     & ! canopy timestep / canopy heat capacity
   ,dtlcowcc     & ! canopy timestep / canopy vapor capacity
   ,dtlcohcv     & ! caonpy timestep / vegetation heat capacity
   
   ,ups          & ! U velocity at top of surface layer [up(2,i,j)]    
   ,vps          & ! V velocity at top of surface layer [vp(2,i,j)]
   ,ths          & ! potential temperature at top of surface layer [theta(2,i,j)]
   ,rvs          & ! vapor mixing ratio at top of surface layer [rv(2,i,j)]
   ,zts          & ! height at top of surface layer [zt(2)*rtgt(i,j)]
   ,pis          & ! Exner function at surface
   ,dens         & ! density at surface
   ,prss         & ! pressure at surface
   ,vels         & ! wind speed at top of surface layer
   ,vels_pat     & ! vels with patch-dependent ubmin imposed as minimum
   ,gzotheta     & ! (g*z/theta) at top of surface layer [for Richardson number]
   ,pcpgl        & ! precip mass from cuparm and/or micphys in leaf timestep
   ,qpcpgl       & ! precip energy from cuparm and/or micphys in leaf timestep
   ,dpcpgl       & ! precip depth from cuparm and/or micphys in leaf timestep
   ,pcpgc        & ! precip mass from cuparm and/or micphys in canopy timestep
   ,qpcpgc       & ! precip energy from cuparm and/or micphys in canopy timestep
   ,z0fac_water  & ! (.016 / g) factor of ustar^2 for z0 over water
   
   ,snowfac      & ! fraction of vegetation height covered by sfcwater
   ,vf           & ! product of veg_fracarea and (1-snowfac)
   ,thetacan     & ! canopy air potential temperature
   ,transp       & ! transpiration flux [kg/m2/s]
   ,snowrough    & ! snowcover roughness height
   ,timefac_sst  & ! time interpolation factor for SST
   
   ,rb           & ! vegetation aerodynamic resistance
   ,rd           & ! canopy to ground aerodynamic resistance
   ,rdi          & ! canopy to ground aerodynamic conductance
   
   ,rshort_g     & ! net SW radiation absorbed by grnd
   ,rshort_v     & ! net SW radiation absorbed by veg
   ,rshort_a     & ! net SW radiation reflected to atm by veg plus grnd
   
   ,rlonga_v     & ! net atm LW radiation absorbed by veg
   ,rlonga_gs    & ! net atm LW radiation absorbed by grnd OR snow
   ,rlongv_gs    & ! net veg LW radiation absorbed by grnd OR snow
   ,rlongv_a     & ! net veg LW radiation to atm
   ,rlonggs_v    & ! net grnd OR snow LW radiation absorbed by veg
   ,rlonggs_a    & ! net grnd OR snow LW radiation to atm
   ,rlonga_a     & ! net atm LW radiation reflected to atm by veg plus grnd OR snow

   ,hflxgc       & ! sensible heat from ground to canopy (J/m2)
   ,wflxgc       & ! water vapor from ground to canopy (kg/m2)
   ,hflxvc       & ! sensible heat from vegetation to canopy (J/m2)
   ,wflxvc       & ! water vapor from vegetation to canopy (kg/m2)
   
   ,wshed        & ! water shed from vegetation to ground (kg/m2)
   ,qwshed       & ! energy from shed water (J/m2)
   ,dewgnd         ! dew formation on ground (kg/m2)

real, save, allocatable, dimension(:) ::  &
    dslz         & ! soil layer thickness at T point
   ,dslzi        & ! (1. / soil layer thickness at T point)
   ,dslzidt      & ! (dtll / soil layer thickness at T point)
   ,slzt         & ! soil depth at T point
   ,dslzt        & ! soil layer thickness at M point
   ,dslzti       & ! (1. / soil layer thickness at M point)
   ,dslztidt     & ! (dtll / soil layer thickness at M point)

   ,rshort_s     & ! net SW radiation absorbed by snow layers (mzs)
   ,tempk        & ! diagnosed temp (K) of soil and sfcwater
   ,fracliq      & ! diagnosed liquid fraction of soil_water and sfcwater_mass
   
   ,hfluxgsc      & ! sensible heat flux between soil, sfcwater, canopy
   ,psiplusz      & ! soil water potential plus geopotential [m]
   ,half_soilair  & ! half of available airspace in soil [m]
   ,rfactor       & ! soil, sfcwater thermal resistance    
   ,wflux         & ! soil water flux [m]
   ,soil_liq      & ! soil liquid water content [m]
   ,qwflux          ! soil energy flux from water flux [J/m2] 


integer, parameter :: nstyp=12,nvtyp=20
! TEB
! Considering an extra urban type. Equivalent to code 21 - Very urban.
integer, parameter :: nvtyp_teb = 1

!-------------------------------------------------------------------------------

real, dimension(nstyp)           :: slden,slcpd,slbs,slcond,sfldcap  &
     ,slcons,slmsts,slpots,ssand,sclay  &
     ,sorgan,sporo,soilcp,slfc,emisg,slcons00

real, dimension(nvtyp+nvtyp_teb) :: albv_green,albv_brown,emisv,sr_max,tai_max &
     ,sai,veg_clump,veg_frac,veg_ht,glai_max  &
     ,dead_frac,rcmin

real                             :: cmin,corg,cwat,cair,cka,ckw

integer, dimension(nvtyp+nvtyp_teb) :: kroot

real, dimension(nzgmax,nvtyp+nvtyp_teb) :: root



End Module leaf_coms


subroutine alloc_leafcol(nzg,nzs)

use leaf_coms

implicit none
integer, intent(in) :: nzg,nzs

! Allocate leaf column arrays
allocate (dslz         (nzg)        &
         ,dslzi        (nzg)        &
         ,dslzidt      (nzg)        &
         ,slzt         (nzg)        &
         ,dslzt        (nzg)        &
         ,dslzti       (nzg)        &
         ,dslztidt     (nzg)        &

         ,rshort_s     (nzs)        &
         ,tempk        (nzg+nzs)    &
         ,fracliq      (nzg+nzs)    &

         ,hfluxgsc     (nzg+nzs+1)  &
         ,psiplusz     (nzg)        &
         ,half_soilair (nzg)        &
         ,rfactor      (nzg+nzs)    &
         ,wflux        (nzg+1)      &
         ,qwflux       (nzg+1)      &
         ,soil_liq     (nzg)        )
         
return
end subroutine     

