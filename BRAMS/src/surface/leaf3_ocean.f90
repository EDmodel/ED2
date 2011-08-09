!==========================================================================================!
!==========================================================================================!
!      This sub-routine will solve the time step for leaf3 over the ocean and water        !
! regions in general.  Although we are not solving algae or water lilies at this time, we  !
! must update the "vegetation" variables so they will have something in case they are ever !
! used; we simply copy the canopy air space temperature and force water and internal       !
! energy to be zero.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine leaf3_ocean(mzg,ustar,tstar,rstar,cstar,patch_rough,can_prss,can_rvap,can_co2   &
                      ,sensible_gc,evap_gc,plresp,ground_temp,ground_rsat,ground_rvap      &
                      ,ground_fliq)
   use leaf_coms , only : dtll           & ! intent(in)
                        , dtll_factor    & ! intent(in)
                        , ggbare         & ! intent(in)
                        , can_rhos       & ! intent(in)
                        , can_exner      & ! intent(in)
                        , soil_tempk     & ! intent(in)
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
                        , cflxac         & ! intent(out)
                        , ustmin         ! ! intent(in)
   use rconstants, only : mmdry          & ! intent(in)
                        , mmdryi         & ! intent(in)
                        , cp             & ! intent(in)
                        , cpi            & ! intent(in)
                        , t3ple          & ! intent(in)
                        , alvl           & ! intent(in)
                        , huge_num       ! ! intent(in)
   use therm_lib , only : rslif          & ! function
                        , thetaeiv       ! ! function
   use node_mod  , only : mynum          ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)    :: mzg
   real   , intent(in)    :: ustar
   real   , intent(in)    :: tstar
   real   , intent(in)    :: rstar
   real   , intent(in)    :: cstar
   real   , intent(in)    :: patch_rough
   real   , intent(in)    :: can_prss
   real   , intent(inout) :: can_rvap
   real   , intent(inout) :: can_co2
   real   , intent(inout) :: sensible_gc
   real   , intent(inout) :: evap_gc
   real   , intent(inout) :: plresp
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
   ground_temp = soil_tempk(mzg)
   ground_rsat = rslif(can_prss,soil_tempk(mzg))
   ground_rvap = ground_rsat
   if (soil_tempk(mzg) >= t3ple) then
      ground_fliq = 1.
   else
      ground_fliq = 0.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    For water patches, the net conductance is the same as the bare ground.             !
   !---------------------------------------------------------------------------------------!
   ggveg   = huge_num
   ggnet   = ggbare
   !---------------------------------------------------------------------------------------!

   !----- Compute the fluxes from water body to canopy. -----------------------------------!
   hflxgc  = ggnet * cp * can_rhos * (ground_temp - can_temp)
   wflxgc  = ggnet      * can_rhos * (ground_rsat - can_rvap)
   qwflxgc = wflxgc * alvl
   cflxgc  = 0. !----- No water carbon emission model available...

   !----- Compute the fluxes from atmosphere to canopy air space. -------------------------!
   rho_ustar = can_rhos  * ustar
   eflxac    = rho_ustar * estar * cp * can_temp
   hflxac    = rho_ustar * tstar * can_exner
   wflxac    = rho_ustar * rstar
   cflxac    = rho_ustar * cstar * mmdryi

   !----- Integrate the state variables. --------------------------------------------------!
   can_lntheta  = can_lntheta + dtllohcc * (hflxgc + hflxac)
   can_rvap     = can_rvap    + dtllowcc * (wflxgc + wflxac)
   can_co2      = can_co2     + dtlloccc * (cflxgc + cflxac)

   !----- Integrate the fluxes. -----------------------------------------------------------!
   sensible_gc = sensible_gc +  hflxgc * dtll_factor
   evap_gc     = evap_gc     + qwflxgc * dtll_factor
   plresp      = plresp      +  cflxgc * dtll_factor
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
subroutine leaf3_ocean_diag(ifm,mzg,pastsst,futuresst,soil_energy)
   use rconstants, only : cliq         & ! intent(in)
                        , tsupercool   ! ! intent(in)
   use mem_grid  , only : time         ! ! intent(in)
   use io_params , only : iupdsst      & ! intent(in)
                        , ssttime1     & ! intent(in)
                        , ssttime2     ! ! intent(in)
   use leaf_coms , only : timefac_sst  & ! intent(out)
                        , soil_tempk   & ! intent(in)
                        , soil_fracliq ! ! intent(in)
   use therm_lib , only : qtk          ! ! sub-routine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: ifm
   integer                     , intent(in)    :: mzg
   real                        , intent(in)    :: pastsst
   real                        , intent(in)    :: futuresst
   real        , dimension(mzg), intent(inout) :: soil_energy
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: izg
   integer                                     :: sst
   integer                                     :: ssq
   !---------------------------------------------------------------------------------------!




   !----- Find the time interpolation factor for updating SST. ----------------------------!
   if (iupdsst == 0) then
      timefac_sst = 0.
   else
      timefac_sst = sngl((time - ssttime1(ifm)) / (ssttime2(ifm) - ssttime1(ifm)))
   endif
   !---------------------------------------------------------------------------------------!


   !----- Find the sea surface temperature. -----------------------------------------------!
   sst = pastsst + timefac_sst * (futuresst -pastsst)
   ssq = cliq * (sst - tsupercool)

   !----- Find the sea surface internal energy, assuming that it is always liquid. --------!
   do izg=1,mzg
      soil_energy(izg) = ssq
      call qtk(soil_energy(izg),soil_tempk(izg),soil_fracliq(izg))
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_ocean_diag
!==========================================================================================!
!==========================================================================================!
