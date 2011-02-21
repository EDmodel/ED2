!==========================================================================================!
!==========================================================================================!
!      This sub-routine will solve the time step for leaf3 over the ocean and water        !
! regions in general.  Although we are not solving algae or water lilies at this time, we  !
! must update the "vegetation" variables so they will have something in case they are ever !
! used; we simply copy the canopy air space temperature and force water and internal       !
! energy to be zero.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine leaf3_ocean(mzg,ustar,tstar,rstar,cstar,can_prss,can_rvap,can_co2,sensible_gc   &
                      ,evap_gc,plresp,can_theta,can_theiv,veg_energy,veg_water,veg_hcap    &
                      ,ground_temp,ground_rsat,ground_rvap,ground_fliq)
   use leaf_coms , only : dtll           & ! intent(in)
                        , dtll_factor    & ! intent(in)
                        , ggbare         & ! intent(in)
                        , can_rhos       & ! intent(in)
                        , can_exner      & ! intent(in)
                        , tempk          & ! intent(in)
                        , dtllowcc       & ! intent(out)
                        , dtllohcc       & ! intent(out)
                        , dtlloccc       & ! intent(out)
                        , can_depth      & ! intent(out)
                        , can_temp       & ! intent(out)
                        , can_lntheta    & ! intent(out)
                        , can_shv        & ! intent(out)
                        , veg_temp       & ! intent(out)
                        , veg_fliq       & ! intent(out)
                        , ggveg          & ! intent(out)
                        , ggfact         & ! intent(out)
                        , ggnet          & ! intent(out)
                        , hflxgc         & ! intent(out)
                        , wflxgc         & ! intent(out)
                        , cflxgc         & ! intent(out)
                        , qwflxgc        & ! intent(out)
                        , rho_ustar      & ! intent(out)
                        , estar          & ! intent(out)
                        , eflxac         & ! intent(out)
                        , hflxac         & ! intent(out)
                        , wflxac         & ! intent(out)
                        , cflxac         ! ! intent(out)
   use rconstants, only : mmdry          & ! intent(in)
                        , mmdryi         & ! intent(in)
                        , cp             & ! intent(in)
                        , cpi            & ! intent(in)
                        , t3ple          & ! intent(in)
                        , alvl           ! ! intent(in)
   use therm_lib , only : rslif          & ! function
                        , thetaeiv       ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)    :: mzg
   real   , intent(in)    :: ustar
   real   , intent(in)    :: tstar
   real   , intent(in)    :: rstar
   real   , intent(in)    :: cstar
   real   , intent(in)    :: can_prss
   real   , intent(inout) :: can_rvap
   real   , intent(inout) :: can_co2
   real   , intent(inout) :: sensible_gc
   real   , intent(inout) :: evap_gc
   real   , intent(inout) :: plresp
   real   , intent(out)   :: can_theta
   real   , intent(out)   :: can_theiv
   real   , intent(out)   :: veg_energy
   real   , intent(out)   :: veg_water
   real   , intent(out)   :: veg_hcap
   real   , intent(out)   :: ground_temp
   real   , intent(out)   :: ground_rsat
   real   , intent(out)   :: ground_rvap
   real   , intent(out)   :: ground_fliq
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     For water patches, update temperature and moisture of "canopy" from divergence of !
   ! fluxes with water surface and atmosphere.                                             !
   !---------------------------------------------------------------------------------------!
   dtllowcc  = dtll / (can_depth * can_rhos)
   dtllohcc  = dtll / (can_depth * can_rhos * cp * can_temp)
   dtlloccc  = mmdry * dtllowcc
   

   !---------------------------------------------------------------------------------------!
   !    Find the ground properties.                                                        !
   !---------------------------------------------------------------------------------------!
   ground_temp = tempk(mzg)
   ground_rsat = rslif(can_prss,tempk(mzg))
   ground_rvap = ground_rsat
   if (tempk(mzg) >= t3ple) then
      ground_fliq = 1.
   else
      ground_fliq = 0.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    For water patches, the net conductance is the same as the bare ground.             !
   !---------------------------------------------------------------------------------------!
   ggveg   = 0.0
   ggnet   = ggbare * ggfact
   !---------------------------------------------------------------------------------------!

   !----- Compute the fluxes from water body to canopy. -----------------------------------!
   hflxgc  = ggnet * cp * can_rhos * (ground_temp - can_temp)
   wflxgc  = ggnet      * can_rhos * (ground_rsat - can_rvap)
   qwflxgc = wflxgc * alvl
   cflxgc  = 0. !----- No water carbon emission model available...

   !----- Compute the fluxes from atmosphere to canopy air space. -------------------------!
   rho_ustar = can_rhos * ustar
   eflxac    = rho_ustar * estar * cp * can_temp
   hflxac    = rho_ustar * tstar * can_exner
   wflxac    = rho_ustar * rstar
   cflxac    = rho_ustar * cstar * mmdryi

   !----- Integrate the state variables. --------------------------------------------------!
   can_lntheta  = can_lntheta + dtllohcc * (hflxgc + hflxac)
   can_rvap     = can_rvap    + dtllowcc * (wflxgc + wflxac)
   can_co2      = can_co2     + dtlloccc * (cflxgc + cflxac)

   !----- Integrate the fluxes. -----------------------------------------------------------!
   sensible_gc = sensible_gc + hflxgc*dtll_factor
   evap_gc     = evap_gc     + wflxgc*dtll_factor
   plresp      = plresp      + cflxgc*dtll_factor
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Update canopy air properties.  This is done inside the loop to ensure that we      !
   ! leave here with the right canopy air potential temperature.                           !
   !---------------------------------------------------------------------------------------!
   can_theta = exp(can_lntheta)
   can_shv   = can_rvap / (can_rvap + 1.)
   can_temp  = cpi * can_theta * can_exner
   can_theiv = thetaeiv(can_theta,can_prss,can_temp,can_rvap,can_rvap,-8)
   !---------------------------------------------------------------------------------------!


   !----- Update the vegetation temperature. ----------------------------------------------!
   veg_temp = can_temp
   if (veg_temp > t3ple) then
      veg_fliq = 1.0
   elseif (veg_temp < t3ple) then
      veg_fliq = 0.0
   else
      veg_fliq = 0.5
   end if
   veg_energy = 0.0
   veg_water  = 0.0
   veg_hcap   = 0.0
   !---------------------------------------------------------------------------------------!


   return
end subroutine leaf3_ocean
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine will find the internal energy of the ocean based on the previous and   !
! future sea surface temperature.                                                          !
!------------------------------------------------------------------------------------------!
subroutine ocean_energy(mzg,m2,m3,mpat,i,j,timefac,pastsst,futuresst,soil_energy)
   use rconstants, only : cliq       & ! intent(in)
                        , tsupercool ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                , intent(in)    :: mzg
   integer                                , intent(in)    :: m2
   integer                                , intent(in)    :: m3
   integer                                , intent(in)    :: mpat
   integer                                , intent(in)    :: i
   integer                                , intent(in)    :: j
   real                                   , intent(in)    :: timefac
   real        , dimension(m2,m3)         , intent(in)    :: pastsst
   real        , dimension(m2,m3)         , intent(in)    :: futuresst
   real        , dimension(mzg,m2,m3,mpat), intent(inout) :: soil_energy
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                :: izg
   integer                                                :: sst
   integer                                                :: ssq
   !---------------------------------------------------------------------------------------!


   !----- Find the sea surface temperature. -----------------------------------------------!
   sst = pastsst(i,j) + timefac * (futuresst(i,j) -pastsst(i,j))
   ssq = cliq * (sst - tsupercool)

   !----- Find the sea surface internal energy, assuming that it is always liquid. --------!
   do izg=1,mzg
      soil_energy(izg,i,j,1) = ssq
   end do

   return
end subroutine ocean_energy
!==========================================================================================!
!==========================================================================================!
