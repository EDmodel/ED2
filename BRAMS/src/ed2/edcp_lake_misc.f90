!==========================================================================================!
!==========================================================================================!
subroutine copy_lake_init(i,j,ifm,initp)
   use mem_leaf              , only : leaf_g           ! ! intent(in)
   use lake_coms             , only : lakesitetype     & ! structure
                                    , lakemet          & ! intent(in)
                                    , wcapcan          & ! intent(in)
                                    , wcapcani         & ! intent(in)
                                    , hcapcani         & ! intent(in)
                                    , ccapcani         ! ! intent(in)
   use consts_coms           , only : cp8              & ! intent(in)
                                    , cpi8             & ! intent(in)
                                    , p00i8            & ! intent(in)
                                    , rocp8            ! ! intent(in)
   use therm_lib8            , only : reducedpress8    & ! function
                                    , thetaeiv8        & ! function
                                    , idealdenssh8     & ! function
                                    , rehuil8          & ! function
                                    , rslif8           ! ! function
   use leaf_coms             , only : min_waterrough8  & ! intent(in)
                                    , z0fac_water8     & ! intent(in)
                                    , can_depth_min    ! ! intent(in)
   use canopy_air_coms       , only : ustmin8          ! ! intent(in)
   use canopy_struct_dynamics, only : ed_stars8        & ! intent(in)
                                    , can_whcap8       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer           , intent(in) :: i
   integer           , intent(in) :: j
   integer           , intent(in) :: ifm
   type(lakesitetype), target     :: initp
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)                   :: ustar8
   !---------------------------------------------------------------------------------------!


   !------ First we copy those variables that do not change with pressure. ----------------!
   initp%can_co2      = dble(leaf_g(ifm)%can_co2  (i,j,1))
   initp%can_theta    = dble(leaf_g(ifm)%can_theta(i,j,1))
   initp%can_theiv    = dble(leaf_g(ifm)%can_theiv(i,j,1))
   initp%can_rvap     = dble(leaf_g(ifm)%can_rvap (i,j,1))
   initp%can_depth    = dble(can_depth_min)

   !------ Convert mixing ratio to specific humidity. -------------------------------------!
   initp%can_shv      = initp%can_rvap / (1.d0 + initp%can_rvap)

   !------ Find the derived properties. ---------------------------------------------------!
   initp%can_lntheta  = log(initp%can_theta)

   initp%can_prss     = reducedpress8(lakemet%atm_prss,lakemet%atm_theta,lakemet%atm_shv   &
                                     ,lakemet%geoht,initp%can_theta,initp%can_shv          &
                                     ,initp%can_depth)

   initp%can_exner    = cp8 * (initp%can_prss * p00i8) ** rocp8

   initp%can_temp     = cpi8 * initp%can_theta * initp%can_exner

   initp%can_rhos     = idealdenssh8(initp%can_prss,initp%can_temp,initp%can_shv)
   initp%can_rhv      = rehuil8(initp%can_prss,initp%can_temp,initp%can_rvap)
   initp%can_ssh      = rslif8(initp%can_prss,initp%can_temp)
   initp%can_ssh      = initp%can_ssh / (initp%can_ssh + 1.d0)
   !---------------------------------------------------------------------------------------!

   !------ Copy in the sea surface temperature and find associated values. ----------------!
   initp%lake_temp    = dble(leaf_g(ifm)%ground_temp(i,j,1))
   initp%lake_fliq    = 1.d0 ! No sea ice for the time being...
   !------ Find the mixing ratio, then convert to specific humidity. ----------------------!
   initp%lake_ssh     = rslif8(initp%can_prss,initp%lake_temp)
   initp%lake_ssh     = initp%lake_ssh / (initp%lake_ssh + 1.d0)
   initp%lake_shv     = initp%lake_ssh
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Assign an initial value for roughness, which is based on the previous ustar.      !
   !---------------------------------------------------------------------------------------!
   ustar8             = max(ustmin8,dble(leaf_g(ifm)%ustar(i,j,1)))
   initp%lake_rough   = max(z0fac_water8 * ustar8 * ustar8,min_waterrough8)
   !---------------------------------------------------------------------------------------!

        
   !----- Find the characteristic scales (a.k.a. stars). ----------------------------------!
   call ed_stars8(lakemet%atm_theta,lakemet%atm_theiv,lakemet%atm_shv,lakemet%atm_co2      &
                 ,initp%can_theta  ,initp%can_theiv  ,initp%can_shv  ,initp%can_co2        &
                 ,lakemet%geoht,0.d0,lakemet%atm_vels,initp%lake_rough                     &
                 ,initp%ustar,initp%tstar,initp%estar,initp%qstar,initp%cstar              &
                 ,initp%zeta,initp%ribulk,initp%gglake)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Apply the conductance factor (should be removed soon).  Also, update the rough-  !
   ! ness so next time we use we have the most up to date value.                           !
   !---------------------------------------------------------------------------------------!
   initp%lake_rough = max(z0fac_water8 * initp%ustar * initp%ustar, min_waterrough8)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Calculate the heat and mass storage capacity of the canopy.                       !
   !---------------------------------------------------------------------------------------!
   call can_whcap8(initp%can_rhos,initp%can_temp,initp%can_depth                           &
                  ,wcapcan,wcapcani,hcapcani,ccapcani)
   !---------------------------------------------------------------------------------------!


   !----- Find the boundaries for the sanity check. ---------------------------------------!
   call lake_derived_thbounds(initp%can_rhos,initp%can_theta,initp%can_temp,initp%can_shv  &
                             ,initp%can_rvap,initp%can_prss,initp%can_depth)

   return
end subroutine copy_lake_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine updates the diagnostic variables.                                    !
!------------------------------------------------------------------------------------------!
subroutine lake_diagnostics(initp)
   use rk4_coms              , only : rk4min_can_shv        & ! intent(in)
                                    , rk4min_can_theta      & ! intent(in)
                                    , rk4max_can_theta      & ! intent(in)
                                    , rk4min_can_lntheta    & ! intent(in)
                                    , rk4max_can_lntheta    & ! intent(in)
                                    , rk4min_can_temp       & ! intent(in)
                                    , rk4max_can_shv        & ! intent(in)
                                    , rk4min_sfcw_temp      & ! intent(in)
                                    , rk4max_sfcw_temp      & ! intent(in)
                                    , tiny_offset           & ! intent(in)
                                    , force_idealgas        ! ! intent(in)
   use lake_coms             , only : lakesitetype          & ! structure
                                    , lakemet               & ! intent(in)
                                    , wcapcan               & ! intent(in)
                                    , wcapcani              & ! intent(in)
                                    , hcapcani              & ! intent(in)
                                    , ccapcani              ! ! intent(in)
   use consts_coms           , only : cp8                   & ! intent(in)
                                    , cpi8                  & ! intent(in)
                                    , p00i8                 & ! intent(in)
                                    , rdry8                 & ! intent(in)
                                    , epim18                & ! intent(in)
                                    , rocp8                 ! ! intent(in)
   use therm_lib8            , only : reducedpress8         & ! function
                                    , thetaeiv8             & ! function
                                    , idealdenssh8          & ! function
                                    , rehuil8               & ! function
                                    , rslif8                & ! function
                                    , thrhsh2temp8          ! ! function
   use leaf_coms             , only : min_waterrough8       & ! intent(in)
                                    , z0fac_water8          & ! intent(in)
                                    , can_depth_min         ! ! intent(in)
   use canopy_air_coms       , only : ustmin8               ! ! intent(in)
   use canopy_struct_dynamics, only : ed_stars8             & ! intent(in)
                                    , can_whcap8            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(lakesitetype), target       :: initp
   !----- Local variables. ----------------------------------------------------------------!
   logical                          :: ok_shv
   logical                          :: ok_theta
   logical                          :: ok_ground
   logical                          :: ok_flpoint
   !----- External functions. -------------------------------------------------------------!
   logical           , external     :: is_finite8
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Sanity check.                                                                     !
   !---------------------------------------------------------------------------------------!
   ok_flpoint = is_finite8(initp%can_shv)  .and. is_finite8(initp%can_lntheta) .and.       &
                is_finite8(initp%can_prss) .and. is_finite8(initp%can_co2)     .and.       &
                is_finite8(initp%lake_temp)
   if (.not. ok_flpoint) then
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) '  Something went wrong...                        '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' Meteorological forcing.'
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rhos           :',lakemet%atm_rhos
      write(unit=*,fmt='(a,1x,es12.5)') ' - Temp           :',lakemet%atm_tmp
      write(unit=*,fmt='(a,1x,es12.5)') ' - Theta          :',lakemet%atm_theta
      write(unit=*,fmt='(a,1x,es12.5)') ' - Theiv          :',lakemet%atm_theiv
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lntheta        :',lakemet%atm_lntheta
      write(unit=*,fmt='(a,1x,es12.5)') ' - Shv            :',lakemet%atm_shv
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rvap           :',lakemet%atm_rvap
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rel. hum.      :',lakemet%atm_rhv
      write(unit=*,fmt='(a,1x,es12.5)') ' - CO2            :',lakemet%atm_co2
      write(unit=*,fmt='(a,1x,es12.5)') ' - Exner          :',lakemet%atm_exner
      write(unit=*,fmt='(a,1x,es12.5)') ' - Press          :',lakemet%atm_prss
      write(unit=*,fmt='(a,1x,es12.5)') ' - Vels           :',lakemet%atm_vels
      write(unit=*,fmt='(a,1x,es12.5)') ' - ucos           :',lakemet%ucos
      write(unit=*,fmt='(a,1x,es12.5)') ' - usin           :',lakemet%usin
      write(unit=*,fmt='(a,1x,es12.5)') ' - Geoht          :',lakemet%geoht
      write(unit=*,fmt='(a,1x,es12.5)') ' - d(SST)/dt      :',lakemet%dsst_dt
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rshort         :',lakemet%rshort
      write(unit=*,fmt='(a,1x,es12.5)') ' - Rlong          :',lakemet%rlong
      write(unit=*,fmt='(a,1x,es12.5)') ' - Tanz           :',lakemet%tanz
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lon            :',lakemet%lon
      write(unit=*,fmt='(a,1x,es12.5)') ' - Lat            :',lakemet%lat
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      write(unit=*,fmt='(a)'          ) ' Runge-Kutta structure.'
      write(unit=*,fmt='(a)'          ) '  - Canopy air space.'
      write(unit=*,fmt='(a,1x,es12.5)') '    * Rhos        :',initp%can_rhos
      write(unit=*,fmt='(a,1x,es12.5)') '    * Temp        :',initp%can_temp
      write(unit=*,fmt='(a,1x,es12.5)') '    * Theta       :',initp%can_theta
      write(unit=*,fmt='(a,1x,es12.5)') '    * Theiv       :',initp%can_theiv
      write(unit=*,fmt='(a,1x,es12.5)') '    * Lntheta     :',initp%can_lntheta
      write(unit=*,fmt='(a,1x,es12.5)') '    * Shv         :',initp%can_shv
      write(unit=*,fmt='(a,1x,es12.5)') '    * Rvap        :',initp%can_rvap
      write(unit=*,fmt='(a,1x,es12.5)') '    * Rel. hum.   :',initp%can_rhv
      write(unit=*,fmt='(a,1x,es12.5)') '    * CO2         :',initp%can_co2
      write(unit=*,fmt='(a,1x,es12.5)') '    * Exner       :',initp%can_exner
      write(unit=*,fmt='(a,1x,es12.5)') '    * Press       :',initp%can_prss
      write(unit=*,fmt='(a,1x,es12.5)') '    * Depth       :',initp%can_depth
      write(unit=*,fmt='(a)'          ) '  - Lake.'
      write(unit=*,fmt='(a,1x,es12.5)') '    * Temp        :',initp%lake_temp
      write(unit=*,fmt='(a,1x,es12.5)') '    * Fliq        :',initp%lake_fliq
      write(unit=*,fmt='(a,1x,es12.5)') '    * Shv         :',initp%lake_shv
      write(unit=*,fmt='(a,1x,es12.5)') '    * Ssh         :',initp%lake_ssh
      write(unit=*,fmt='(a,1x,es12.5)') '    * Roughness   :',initp%lake_rough
      write(unit=*,fmt='(a)'          ) '  - Stars.'
      write(unit=*,fmt='(a,1x,es12.5)') '    * Ustar       :',initp%ustar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Tstar       :',initp%tstar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Estar       :',initp%estar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Qstar       :',initp%qstar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Cstar       :',initp%cstar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Zeta        :',initp%zeta
      write(unit=*,fmt='(a,1x,es12.5)') '    * Ribulk      :',initp%ribulk
      write(unit=*,fmt='(a,1x,es12.5)') '    * GGlake      :',initp%gglake
      write(unit=*,fmt='(a)'          ) '  - Partially integrated fluxes.'
      write(unit=*,fmt='(a,1x,es12.5)') '    * Vapour_gc   :',initp%avg_vapor_gc
      write(unit=*,fmt='(a,1x,es12.5)') '    * Vapour_ac   :',initp%avg_vapor_ac
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sensible_gc :',initp%avg_sensible_gc
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sensible_ac :',initp%avg_sensible_ac
      write(unit=*,fmt='(a,1x,es12.5)') '    * Carbon_gc   :',initp%avg_carbon_gc
      write(unit=*,fmt='(a,1x,es12.5)') '    * Carbon_ac   :',initp%avg_carbon_ac
      write(unit=*,fmt='(a,1x,es12.5)') '    * Carbon_st   :',initp%avg_carbon_st
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sflux_u     :',initp%avg_sflux_u
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sflux_v     :',initp%avg_sflux_u
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sflux_w     :',initp%avg_sflux_w
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sflux_t     :',initp%avg_sflux_t
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sflux_r     :',initp%avg_sflux_r
      write(unit=*,fmt='(a,1x,es12.5)') '    * Sflux_c     :',initp%avg_sflux_c
      write(unit=*,fmt='(a,1x,es12.5)') '    * Albedt      :',initp%avg_albedt
      write(unit=*,fmt='(a,1x,es12.5)') '    * Rlongup     :',initp%avg_rlongup
      write(unit=*,fmt='(a,1x,es12.5)') '    * Rshort_gnd  :',initp%avg_rshort_gnd
      write(unit=*,fmt='(a,1x,es12.5)') '    * Ustar       :',initp%avg_ustar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Tstar       :',initp%avg_tstar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Qstar       :',initp%avg_qstar
      write(unit=*,fmt='(a,1x,es12.5)') '    * Cstar       :',initp%avg_cstar
      write(unit=*,fmt='(a)'          ) ' '
      write(unit=*,fmt='(a)'          ) '-------------------------------------------------'
      call abort_run('Non-resolvable values','lake_diagnostics','edcp_lake_misc.f90')
   end if


   !----- Then we define some logicals to make the code cleaner. --------------------------!
   ok_shv    = initp%can_shv     >= rk4min_can_shv       .and.                             &
               initp%can_shv     <= rk4max_can_shv
   ok_theta  = initp%can_lntheta >= rk4min_can_lntheta   .and.                             &
               initp%can_lntheta <= rk4max_can_lntheta
   ok_ground = initp%lake_temp   >= rk4min_sfcw_temp     .and.                             &
               initp%lake_temp   <= rk4max_sfcw_temp
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Here we convert theta into temperature, potential temperature, and density, and   !
   ! ice-vapour equivalent potential temperature.  The latter variable (or its natural     !
   ! log) should eventually become the prognostic variable for canopy air space entropy    !
   ! when we add condensed/frozen water in the canopy air space.                           !
   !---------------------------------------------------------------------------------------!
   if (ok_shv .and. ok_theta) then

      !----- First, we update the canopy air potential temperature. -----------------------!
      initp%can_theta = exp(initp%can_lntheta)

      !----- Update the canop air 
      initp%can_rvap  = initp%can_shv / (1.d0 - initp%can_shv)

      !------------------------------------------------------------------------------------!
      !    Here we find the temperature in different ways depending on whether we are      !
      ! keeping pressure constant during one full time step or not.  If we are forcing     !
      ! ideal gas to be always respected, then we don't know the pressure until we have    !
      ! the temperature, so we compute the temperature based on potential temperature,     !
      ! density, and specific humidity, then update pressure.  Otherwise, we compute the   !
      ! temperature using the known pressure, even though this causes the ideal gas law to !
      ! not be always satisfied.                                                           !
      !------------------------------------------------------------------------------------!
      if (force_idealgas) then
         initp%can_temp  = thrhsh2temp8(initp%can_theta,initp%can_rhos,initp%can_shv)
         initp%can_prss  = initp%can_rhos * rdry8 * initp%can_temp                         &
                         * (1.d0 + epim18 * initp%can_shv)
         initp%can_exner = cp8 * (initp%can_prss * p00i8) ** rocp8
      else
         initp%can_temp  = cpi8 * initp%can_theta * initp%can_exner
      end if
      !------------------------------------------------------------------------------------!


      !----- Find derived properties. -----------------------------------------------------!
      initp%can_rhv   = rehuil8(initp%can_prss,initp%can_temp,initp%can_rvap)
      initp%can_ssh   = rslif8(initp%can_prss,initp%can_temp)
      initp%can_ssh   = initp%can_ssh / (initp%can_ssh + 1.d0)
      initp%can_theiv = thetaeiv8(initp%can_theta,initp%can_prss,initp%can_temp            &
                                 ,initp%can_rvap,initp%can_rvap)
      !------------------------------------------------------------------------------------!



      !----- Find the characteristic scales (a.k.a. stars). -------------------------------!
      call ed_stars8(lakemet%atm_theta,lakemet%atm_theiv,lakemet%atm_shv,lakemet%atm_co2   &
                    ,initp%can_theta  ,initp%can_theiv  ,initp%can_shv  ,initp%can_co2     &
                    ,lakemet%geoht,0.d0,lakemet%atm_vels,initp%lake_rough                  &
                    ,initp%ustar,initp%tstar,initp%estar,initp%qstar,initp%cstar           &
                    ,initp%zeta,initp%ribulk,initp%gglake)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Apply the conductance factor (should be removed soon).  Also, update the      !
      ! roughness so next time we use we have the most up to date value.                   !
      !------------------------------------------------------------------------------------!
      initp%lake_rough = max(z0fac_water8 * initp%ustar * initp%ustar,min_waterrough8)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Calculate the heat and mass storage capacity of the canopy.                    !
      !------------------------------------------------------------------------------------!
      call can_whcap8(initp%can_rhos,initp%can_temp,initp%can_depth                        &
                     ,wcapcan,wcapcani,hcapcani,ccapcani)
      !------------------------------------------------------------------------------------!

   elseif (initp%can_lntheta >= rk4max_can_lntheta) then
      !----- CAS is too hot, put a non-sense temperature so the sanity check fails. -------!
      initp%can_theta = rk4max_can_theta + 1.d0

   elseif (initp%can_lntheta <= rk4min_can_lntheta) then
      !----- CAS is too cold, put a non-sense temperature so the sanity check fails. ------!
      initp%can_theta = rk4min_can_theta - 1.d0
   end if
   !---------------------------------------------------------------------------------------!




   !------ Find the mixing ratio, then convert to specific humidity. ----------------------!
   if (ok_ground) then
      initp%lake_ssh     = rslif8(initp%can_prss,initp%lake_temp)
      initp%lake_ssh     = initp%lake_ssh / (1.d0 + initp%lake_ssh)
      initp%lake_shv     = initp%lake_ssh
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine lake_diagnostics
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine finds the derivatives for the prognostic variables.                  !
!------------------------------------------------------------------------------------------!
subroutine lake_derivs(initp,dinitp)
   use lake_coms             , only : lakesitetype        & ! structure
                                    , lakemet             & ! intent(in)
                                    , hcapcani            & ! intent(in)
                                    , wcapcani            & ! intent(in)
                                    , ccapcani            & ! intent(in)
                                    , emiss_h2o           & ! intent(in)
                                    , albt_inter          & ! intent(in)
                                    , albt_slope          & ! intent(in)
                                    , albt_min            & ! intent(in)
                                    , albt_max            ! ! intent(in)
   use consts_coms           , only : cp8                 & ! intent(in)
                                    , mmdryi8             & ! intent(in)
                                    , stefan8             ! ! intent(in)
   use canopy_struct_dynamics, only : vertical_vel_flux8  ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(lakesitetype) :: initp     ! ! Current guess
   type(lakesitetype) :: dinitp    ! ! Current derivative
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)       :: wflxac
   real(kind=8)       :: hflxac
   real(kind=8)       :: cflxac
   real(kind=8)       :: wflxgc
   real(kind=8)       :: hflxgc
   real(kind=8)       :: cflxgc
   real(kind=8)       :: rho_ustar
   !---------------------------------------------------------------------------------------!



   !----- Find the canopy => free atmosphere fluxes. --------------------------------------!
   rho_ustar = initp%can_rhos * initp%ustar
   wflxac    = rho_ustar * initp%qstar
   hflxac    = rho_ustar * initp%tstar * initp%can_exner
   cflxac    = rho_ustar * initp%cstar * mmdryi8
   !---------------------------------------------------------------------------------------!



   !----- Find the ground => canopy air space fluxes. -------------------------------------!
   wflxgc    =       initp%can_rhos * initp%gglake * (initp%lake_shv  - initp%can_shv )
   hflxgc    = cp8 * initp%can_rhos * initp%gglake * (initp%lake_temp - initp%can_temp)
   cflxgc    = 0.d0 ! We should add a simple ocean flux model at some point in the future.
   !---------------------------------------------------------------------------------------!



   !----- Find the derivatives of the canopy air space. -----------------------------------!
   dinitp%can_shv     = (wflxgc + wflxac) * wcapcani
   dinitp%can_lntheta = (hflxgc + hflxac) * hcapcani
   dinitp%can_co2     = (cflxgc + cflxac) * ccapcani
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      The time derivative of SST is already stored in lakemet, simply copy it to the   !
   ! buffer.                                                                               !
   !---------------------------------------------------------------------------------------!
   dinitp%lake_temp   = lakemet%dsst_dt
   !---------------------------------------------------------------------------------------!



   !----- Copy the fluxes to the integrals for average. -----------------------------------!
   dinitp%avg_vapor_gc     = wflxgc
   dinitp%avg_vapor_ac     = wflxac
   dinitp%avg_sensible_gc  = hflxgc
   dinitp%avg_sensible_ac  = hflxac
   dinitp%avg_carbon_gc    = cflxgc
   dinitp%avg_carbon_ac    = cflxac
   dinitp%avg_carbon_st    = cflxgc + cflxac
   !---------------------------------------------------------------------------------------!



   !----- Copy the fluxes to the integrals for average. -----------------------------------!
   dinitp%avg_sflux_u  = - initp%ustar * initp%ustar * lakemet%ucos
   dinitp%avg_sflux_v  = - initp%ustar * initp%ustar * lakemet%usin
   dinitp%avg_sflux_w  = vertical_vel_flux8(initp%zeta,initp%tstar,initp%ustar)
   dinitp%avg_sflux_t  = - initp%ustar * initp%tstar
   dinitp%avg_sflux_c  = - initp%ustar * initp%cstar
   !----- BRAMS needs the flux in terms of mixing ratio, not specific humidity. -----------!
   dinitp%avg_sflux_r  = - initp%ustar * initp%qstar                                       &
                       / ((1.d0-lakemet%atm_shv)*(1.d0 - initp%can_shv))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the contribution of this time step to the longwave fluxes.                   !
   !---------------------------------------------------------------------------------------!
   dinitp%avg_albedt     = min(max(albt_inter + albt_slope*lakemet%tanz,albt_min),albt_max)
   dinitp%avg_rlongup    = emiss_h2o * stefan8 * initp%lake_temp ** 4
   dinitp%avg_rshort_gnd = dinitp%avg_albedt * lakemet%rshort
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the contribution of this time step to the stars.                             !
   !---------------------------------------------------------------------------------------!
   dinitp%avg_ustar     = initp%ustar
   dinitp%avg_tstar     = initp%tstar
   dinitp%avg_qstar     = initp%qstar
   dinitp%avg_cstar     = initp%cstar
   !---------------------------------------------------------------------------------------!

   return
end subroutine lake_derivs
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine finds the derived properties for the Runge-Kutta/Euler tests on     !
! water regions.                                                                           !
!------------------------------------------------------------------------------------------!
subroutine lake_derived_thbounds(can_rhos,can_theta,can_temp,can_shv,can_rvap,can_prss     &
                                ,can_depth)
   use consts_coms , only : p008                & ! intent(in)
                          , rocp8               & ! intent(in)
                          , cp8                 & ! intent(in)
                          , rdry8               & ! intent(in)
                          , epim18              & ! intent(in)
                          , ep8                 & ! intent(in)
                          , mmdryi8             ! ! intent(in)
   use therm_lib8  , only : thetaeiv8           & ! function
                          , thetaeivs8          & ! function
                          , idealdenssh8        & ! function
                          , reducedpress8       & ! function
                          , eslif8              & ! function
                          , rslif8              ! ! function
   use lake_coms   , only : lakemet             ! ! intent(in)
   use rk4_coms    , only : rk4min_can_prss     & ! intent(in)
                          , rk4max_can_prss     & ! intent(in)
                          , rk4min_can_theta    & ! intent(in)
                          , rk4max_can_theta    & ! intent(in)
                          , rk4min_can_theiv    & ! intent(in)
                          , rk4max_can_theiv    & ! intent(in)
                          , rk4min_can_temp     & ! intent(in)
                          , rk4max_can_temp     & ! intent(in)
                          , rk4min_can_shv      & ! intent(in)
                          , rk4max_can_shv      & ! intent(in)
                          , rk4max_can_rvap      & ! intent(in)
                          , rk4min_can_lntheta  & ! intent(in)
                          , rk4max_can_lntheta  ! ! intent(in)
   implicit none
   !----- Arguments. -------------------------------------------------------------------!
   real(kind=8)                , intent(in) :: can_rhos
   real(kind=8)                , intent(in) :: can_theta
   real(kind=8)                , intent(in) :: can_temp
   real(kind=8)                , intent(in) :: can_shv
   real(kind=8)                , intent(in) :: can_rvap
   real(kind=8)                , intent(in) :: can_prss
   real(kind=8)                , intent(in) :: can_depth
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)                             :: can_prss_try
   real(kind=8)                             :: can_theta_try
   real(kind=8)                             :: can_theiv_try
   integer                                  :: k
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the bounds for pressure.  To avoid the pressure range to be too relaxed,     !
   ! switch one of the dependent variable a time, and use the current values for the       !
   ! other.  In addition, we force pressure to be bounded between the reduced pressure in  !
   ! case the reference height was different by the order of 10%.                          !
   !---------------------------------------------------------------------------------------!
   !----- 1. Initial value, the most extreme one. -----------------------------------------!
   rk4min_can_prss = reducedpress8(lakemet%atm_prss,lakemet%atm_theta,lakemet%atm_shv      &
                                  ,5.d-1*lakemet%geoht,can_theta,can_shv,can_depth)
   rk4max_can_prss = reducedpress8(lakemet%atm_prss,lakemet%atm_theta,lakemet%atm_shv      &
                                  ,1.1d0*lakemet%geoht,can_theta,can_shv,can_depth)
   !----- 2. Minimum temperature. ---------------------------------------------------------!
   can_prss_try    = can_rhos * rdry8 * rk4min_can_temp * (1.d0 + epim18 * can_shv)
   rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
   rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
   !----- 3. Maximum temperature. ---------------------------------------------------------!
   can_prss_try    = can_rhos * rdry8 * rk4max_can_temp * (1.d0 + epim18 * can_shv)
   rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
   rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
   !----- 4. Minimum specific humidity. ---------------------------------------------------!
   can_prss_try    = can_rhos * rdry8 * can_temp * (1.d0 + epim18 * rk4min_can_shv)
   rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
   rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
   !----- 5. Maximum specific humidity. ---------------------------------------------------!
   can_prss_try    = can_rhos * rdry8 * can_temp * (1.d0 + epim18 * rk4max_can_shv)
   rk4min_can_prss = min(rk4min_can_prss,can_prss_try)
   rk4max_can_prss = max(rk4max_can_prss,can_prss_try)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the bounds for potential temperature.  To avoid the pressure range to be too !
   ! relaxed, switch one of the dependent variable a time, and use the current values for  !
   ! the other.                                                                            !
   !---------------------------------------------------------------------------------------!
   !----- 1. Initial value, the most extreme one. -----------------------------------------!
   rk4min_can_theta   =  huge(1.d0)
   rk4max_can_theta   = -huge(1.d0)
   !----- 2. Minimum temperature. ---------------------------------------------------------!
   can_theta_try      = rk4min_can_temp * (p008 / can_prss) ** rocp8
   rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
   rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
   !----- 3. Maximum temperature. ---------------------------------------------------------!
   can_theta_try      = rk4max_can_temp * (p008 / can_prss) ** rocp8
   rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
   rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
   !----- 4. Minimum pressure. ------------------------------------------------------------!
   can_theta_try      = can_temp * (p008 / rk4min_can_prss) ** rocp8
   rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
   rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
   !----- 5. Maximum pressure. ------------------------------------------------------------!
   can_theta_try      = can_temp * (p008 / rk4max_can_prss) ** rocp8
   rk4min_can_theta   = min(rk4min_can_theta,can_theta_try)
   rk4max_can_theta   = max(rk4max_can_theta,can_theta_try)
   !----- 6. Find the logarithms. ---------------------------------------------------------!
   rk4min_can_lntheta = log(rk4min_can_theta)
   rk4max_can_lntheta = log(rk4max_can_theta)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum and maximum ice-vapour equivalent potential temperature.                 !
   !---------------------------------------------------------------------------------------!
   !----- 1. Initial value, the most extreme one. -----------------------------------------!
   rk4min_can_theiv = rk4min_can_theta
   rk4max_can_theiv = -huge(1.d0)
   !----- 2. Maximum temperature. ---------------------------------------------------------!
   can_theta_try    = rk4max_can_temp * (p008 / can_prss) ** rocp8
   can_theiv_try    = thetaeivs8(can_theta_try,rk4max_can_temp,can_rvap,0.d0,0.d0)
   rk4max_can_theiv = max(rk4max_can_theiv,can_theiv_try)
   !----- 3. Minimum pressure. ------------------------------------------------------------!
   can_theta_try    = can_temp * (p008 / rk4min_can_prss) ** rocp8
   can_theiv_try    = thetaeivs8(can_theta_try,can_temp,can_rvap,0.d0,0.d0)
   rk4max_can_theiv = max(rk4max_can_theiv,can_theiv_try)
   !----- 4. Maximum vapour mixing ratio. -------------------------------------------------!
   can_theta_try    = can_temp * (p008 / can_prss) ** rocp8
   can_theiv_try    = thetaeivs8(can_theta_try,can_temp,rk4max_can_rvap,0.d0,0.d0)
   rk4max_can_theiv = max(rk4max_can_theiv,can_theiv_try)
   !---------------------------------------------------------------------------------------!

   return
end subroutine lake_derived_thbounds
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will check for potentially serious problems.  Notice that the upper   !
! and lower bound are defined in rk4_coms.f90, so if you need to change any limit for      !
! some reason, you can adjust there.                                                       !
!------------------------------------------------------------------------------------------!
subroutine lake_sanity_check(y,reject_step,dydx,h,print_problems)
   use rk4_coms , only : rk4eps                & ! intent(in)
                       , rk4min_can_theiv      & ! intent(in)
                       , rk4max_can_theiv      & ! intent(in)
                       , rk4min_can_theta      & ! intent(in)
                       , rk4max_can_theta      & ! intent(in)
                       , rk4max_can_shv        & ! intent(in)
                       , rk4min_can_shv        & ! intent(in)
                       , rk4min_can_rhv        & ! intent(in)
                       , rk4max_can_rhv        & ! intent(in)
                       , rk4min_can_temp       & ! intent(in)
                       , rk4max_can_temp       & ! intent(in)
                       , rk4min_can_theiv      & ! intent(in)
                       , rk4max_can_theiv      & ! intent(in)
                       , rk4min_can_prss       & ! intent(in)
                       , rk4max_can_prss       & ! intent(in)
                       , rk4min_can_co2        & ! intent(in)
                       , rk4max_can_co2        & ! intent(in)
                       , rk4min_sfcw_temp      & ! intent(in)
                       , rk4max_sfcw_temp      ! ! intent(in)
   use lake_coms, only : lakesitetype          ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(lakesitetype), target      :: y
   type(lakesitetype), target      :: dydx
   logical           , intent(in)  :: print_problems
   logical           , intent(out) :: reject_step
   real(kind=8)      , intent(in)  :: h
   !----- Local variables -----------------------------------------------------------------!
   integer                         :: k
   !---------------------------------------------------------------------------------------!


   !----- Be optimistic and start assuming that things are fine. --------------------------!
   reject_step = .false.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air equivalent potential temperature is off.               !
   !---------------------------------------------------------------------------------------!
   if (y%can_theiv > rk4max_can_theiv .or. y%can_theiv < rk4min_can_theiv ) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Canopy air theta_Eiv is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air potential temperature is off.                          !
   !---------------------------------------------------------------------------------------!
   if (y%can_theta > rk4max_can_theta .or. y%can_theta < rk4min_can_theta ) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Canopy air pot. temp. is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air equivalent potential temperature is off.               !
   !---------------------------------------------------------------------------------------!
   if ( y%can_shv > rk4max_can_shv .or. y%can_shv < rk4min_can_shv  ) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Canopy air sp. humidity is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air equivalent potential temperature is off.               !
   !---------------------------------------------------------------------------------------!
   if (y%can_temp > rk4max_can_temp .or. y%can_temp < rk4min_can_temp) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Canopy air temperature is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air equivalent potential temperature is off.               !
   !---------------------------------------------------------------------------------------!
   if (y%can_prss > rk4max_can_prss .or. y%can_prss < rk4min_can_prss) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Canopy air pressure is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air equivalent potential temperature is off.               !
   !---------------------------------------------------------------------------------------!
   if (y%can_co2 > rk4max_can_co2 .or. y%can_co2 < rk4min_can_co2) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Canopy air CO2 is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !   Check whether the canopy air equivalent potential temperature is off.               !
   !---------------------------------------------------------------------------------------!
   if (y%lake_temp > rk4max_sfcw_temp .or. y%lake_temp < rk4min_sfcw_temp) then
      reject_step = .true.
      if (print_problems) then
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' + Lake(SST) temperature is off-track...'
         write(unit=*,fmt='(a)')           '-------------------------------------------'
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THEIV:         ',y%can_theiv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_THETA:         ',y%can_theta
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_SHV:           ',y%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHV:           ',y%can_rhv
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_TEMP:          ',y%can_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_RHOS:          ',y%can_rhos
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_CO2:           ',y%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_DEPTH:         ',y%can_depth
         write(unit=*,fmt='(a,1x,es12.4)') ' CAN_PRSS:          ',y%can_prss
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_TEMP:         ',y%lake_temp
         write(unit=*,fmt='(a,1x,es12.4)') ' LAKE_SHV:          ',y%lake_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_LNTHETA )/Dt:',dydx%can_lntheta
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_SHV     )/Dt:',dydx%can_shv
         write(unit=*,fmt='(a,1x,es12.4)') ' D(CAN_CO2     )/Dt:',dydx%can_co2
         write(unit=*,fmt='(a,1x,es12.4)') ' D(LAKE_TEMP   )/Dt:',dydx%lake_temp
         write(unit=*,fmt='(a)')           '==========================================='
         write(unit=*,fmt='(a)')           ' '
      else
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!


   if (reject_step .and. print_problems) then
      write(unit=*,fmt='(a)')           ' '
      write(unit=*,fmt='(78a)')         ('=',k=1,78)
      write(unit=*,fmt='(a,1x,es12.4)') ' TIMESTEP:          ',h
      write(unit=*,fmt='(a)')           ' '
      write(unit=*,fmt='(a)')           '         ---- SANITY CHECK BOUNDS ----'
      write(unit=*,fmt='(a)')           ' '
      write(unit=*,fmt='(a)')           ' 1. CANOPY AIR SPACE: '
      write(unit=*,fmt='(a)')           ' '
      write(unit=*,fmt='(6(a,1x))')     '   MIN_THEIV','   MAX_THEIV','     MIN_SHV'    &
                                       ,'     MAX_SHV','     MIN_RHV','     MAX_RHV'
      write(unit=*,fmt='(6(es12.5,1x))')  rk4min_can_theiv,rk4max_can_theiv             &
                                         ,rk4min_can_shv  ,rk4max_can_shv                &
                                         ,rk4min_can_rhv  ,rk4max_can_rhv
      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(4(a,1x))')     '    MIN_TEMP','    MAX_TEMP','   MIN_THETA'    &
                                       ,'   MAX_THETA'
      write(unit=*,fmt='(4(es12.5,1x))') rk4min_can_temp ,rk4max_can_temp               &
                                        ,rk4min_can_theta,rk4max_can_theta
      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(4(a,1x))')     '    MIN_PRSS','    MAX_PRSS','     MIN_CO2'    &
                                       ,'     MAX_CO2'
      write(unit=*,fmt='(4(es12.5,1x))') rk4min_can_prss ,rk4max_can_prss               &
                                        ,rk4min_can_co2  ,rk4max_can_co2
      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(78a)')         ('-',k=1,78)
      write(unit=*,fmt='(a)')           ' '
      write(unit=*,fmt='(a)')           ' 2. LAKE PROPERTIES: '
      write(unit=*,fmt='(2(a,1x))')     '    MIN_TEMP','    MAX_TEMP'
      write(unit=*,fmt='(2(es12.5,1x))') rk4min_sfcw_temp ,rk4max_sfcw_temp
      write(unit=*,fmt='(a)')           ' '
      write(unit=*,fmt='(78a)')         ('=',k=1,78)
      write(unit=*,fmt='(a)')           ' '
   end if

   return
end subroutine lake_sanity_check
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine prints the ocean version of everyone's favourite output message in    !
! ED.  After this message is printed, the model will stop.                                 !
!------------------------------------------------------------------------------------------!
subroutine print_lakesite(y,initp,htry)
   use lake_coms   , only : lakesitetype          & ! structure
                          , lakemet               ! ! intent(in)
   use ed_misc_coms, only : current_time          ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(lakesitetype) , target     :: y
   type(lakesitetype) , target     :: initp
   real(kind=8)       , intent(in) :: htry
   !----- Local variables. ----------------------------------------------------------------!
   integer                         :: k
   !---------------------------------------------------------------------------------------!


   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)') ' '
   write (unit=*,fmt='(a,1x,2(i2.2,a),i4.4,1x,f12.0,1x,a)')                                &
         'Time:',current_time%month,'/',current_time%date,'/',current_time%year            &
                ,current_time%time,'s'
   write(unit=*,fmt='(a,1x,es12.4)') 'Attempted step size:',htry
   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(80a)')         ('-',k=1,80)
   write (unit=*,fmt='(a)')           ' ATMOSPHERIC CONDITIONS: '
   write (unit=*,fmt='(a,1x,es12.4)') ' Air temperature       : ',lakemet%atm_tmp
   write (unit=*,fmt='(a,1x,es12.4)') ' Air potential temp.   : ',lakemet%atm_theta
   write (unit=*,fmt='(a,1x,es12.4)') ' Air theta_Eiv         : ',lakemet%atm_theiv
   write (unit=*,fmt='(a,1x,es12.4)') ' H2Ov mixing ratio     : ',lakemet%atm_shv
   write (unit=*,fmt='(a,1x,es12.4)') ' H2Ov rel. humidity    : ',lakemet%atm_rhv
   write (unit=*,fmt='(a,1x,es12.4)') ' CO2  mixing ratio     : ',lakemet%atm_co2
   write (unit=*,fmt='(a,1x,es12.4)') ' Density               : ',lakemet%atm_rhos
   write (unit=*,fmt='(a,1x,es12.4)') ' Pressure              : ',lakemet%atm_prss
   write (unit=*,fmt='(a,1x,es12.4)') ' Exner function        : ',lakemet%atm_exner
   write (unit=*,fmt='(a,1x,es12.4)') ' Wind speed            : ',lakemet%atm_vels
   write (unit=*,fmt='(a,1x,es12.4)') ' Height                : ',lakemet%geoht
   write (unit=*,fmt='(a,1x,es12.4)') ' d(SST)/dt             : ',lakemet%dsst_dt
   write (unit=*,fmt='(a,1x,es12.4)') ' Downward SW radiation : ',lakemet%rshort
   write (unit=*,fmt='(a,1x,es12.4)') ' Downward LW radiation : ',lakemet%rlong
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)') ' '


   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (lake_buff%y) ||||'

   write (unit=*,fmt='(7(a12,1x))')   '   LAKE_TEMP','    LAKE_SHV','  LAKE_ROUGH'         &
                                     ,'   CAN_DEPTH','     CAN_CO2','    CAN_PRSS'         &
                                     ,'      GGLAKE'
                                     
   write (unit=*,fmt='(7(es12.4,1x))') y%lake_temp,y%lake_shv,y%lake_rough,y%can_depth     &
                                      ,y%can_co2,y%can_prss,y%gglake
   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(8(a12,1x))')  '    CAN_RHOS','   CAN_THEIV','   CAN_THETA'          &
                                    ,'    CAN_TEMP','     CAN_SHV','     CAN_SSH'          &
                                    ,'    CAN_RVAP','     CAN_RHV'
                                     
                                     
   write (unit=*,fmt='(8(es12.4,1x))')   y%can_rhos , y%can_theiv, y%can_theta             &
                                       , y%can_temp , y%can_shv  , y%can_ssh               &
                                       , y%can_rvap , y%can_rhv
                                       

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                    ,'       TSTAR','       ESTAR','        ZETA'          &
                                    ,'     RI_BULK'
   write (unit=*,fmt='(7(es12.4,1x))') y%ustar,y%qstar,y%cstar,y%tstar,y%estar,y%zeta      &
                                      ,y%ribulk

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '


   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)

   write(unit=*,fmt='(a)')  ' |||| Printing PATCH information (lake_buff%initp) ||||'

   write (unit=*,fmt='(7(a12,1x))')   '   LAKE_TEMP','    LAKE_SHV','  LAKE_ROUGH'         &
                                     ,'   CAN_DEPTH','     CAN_CO2','    CAN_PRSS'         &
                                     ,'      GGLAKE'

   write (unit=*,fmt='(7(es12.4,1x))') initp%lake_temp,initp%lake_shv,initp%lake_rough     &
                                      ,initp%can_depth,initp%can_co2,initp%can_prss        &
                                      ,initp%gglake

   write (unit=*,fmt='(80a)') ('-',k=1,80)
   write (unit=*,fmt='(8(a12,1x))')  '    CAN_RHOS','   CAN_THEIV','   CAN_THETA'          &
                                    ,'    CAN_TEMP','     CAN_SHV','     CAN_SSH'          &
                                    ,'    CAN_RVAP','     CAN_RHV'

   write (unit=*,fmt='(8(es12.4,1x))')   initp%can_rhos , initp%can_theiv, initp%can_theta &
                                       , initp%can_temp , initp%can_shv  , initp%can_ssh   &
                                       , initp%can_rvap , initp%can_rhv

   write (unit=*,fmt='(80a)') ('-',k=1,80)

   write (unit=*,fmt='(7(a12,1x))')  '       USTAR','       QSTAR','       CSTAR'          &
                                    ,'       TSTAR','       ESTAR','        ZETA'          &
                                    ,'     RI_BULK'
   write (unit=*,fmt='(7(es12.4,1x))') initp%ustar,initp%qstar,initp%cstar,initp%tstar     &
                                      ,initp%estar,initp%zeta,initp%ribulk

   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) ' '

   call abort_run('IFLAG1 problem. The model didn''t converge!','print_lakesite'           &
                 ,'edcp_lake_misc.f90')
   return
end subroutine print_lakesite
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine prints a table with relative and absolute error when the integrator !
! fails.  Useful for debugging and upsetting the user.                                     !
!------------------------------------------------------------------------------------------!
subroutine print_lake_errmax(errmax,yerr,yscal,y,ytemp)
   use rk4_coms , only : rk4eps       ! ! intent(in)
   use lake_coms, only : lakesitetype & ! structure
                       , lakemet      ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(lakesitetype) , target       :: yerr
   type(lakesitetype) , target       :: yscal
   type(lakesitetype) , target       :: y
   type(lakesitetype) , target       :: ytemp
   real(kind=8)       , intent(out)  :: errmax
   !----- Local variables -----------------------------------------------------------------!
   integer                           :: k
   logical                           :: troublemaker
   real(kind=8)                      :: thiserr
   !----- Constants -----------------------------------------------------------------------!
   character(len=28)  , parameter    :: onefmt = '(a16,1x,5(es12.4,1x),11x,l1)'
   !----- Functions -----------------------------------------------------------------------!
   logical            , external     :: large_error
   !---------------------------------------------------------------------------------------!

   errmax = 0.d0

   write(unit=*,fmt='(80a)'    ) ('=',k=1,80)
   write(unit=*,fmt='(a)'      ) '  ..... PRINTING MAXIMUM ERROR INFORMATION: .....'
   write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
   write(unit=*,fmt='(a)'      ) 
   write(unit=*,fmt='(a)'      ) ' Patch level variables, single layer:'
   write(unit=*,fmt='(80a)'    ) ('-',k=1,80)
   write(unit=*,fmt='(7(a,1x))')  'Name            ','       ytemp','       yprev'         &
                                     ,'   Rel.Error','   Abs.Error','       Scale'         &
                                     ,'Problem(T|F)'

   thiserr      = abs(yerr%can_lntheta/yscal%can_lntheta)
   errmax       = max(errmax,thiserr)
   troublemaker = large_error(yerr%can_theiv,yscal%can_theiv)
   write(unit=*,fmt=onefmt) 'CAN_LNTHETA:',thiserr,ytemp%can_lntheta,y%can_lntheta         &
                                          ,yerr%can_lntheta,yscal%can_lntheta,troublemaker

   thiserr      = abs(yerr%can_shv/yscal%can_shv)
   errmax       = max(errmax,thiserr)
   troublemaker = large_error(yerr%can_shv,yscal%can_shv)
   write(unit=*,fmt=onefmt) 'CAN_SHV:',thiserr,ytemp%can_shv,y%can_shv,yerr%can_shv        &
                                      ,yscal%can_shv,troublemaker

   thiserr      = abs(yerr%can_co2/yscal%can_co2)
   errmax       = max(errmax,thiserr)
   troublemaker = large_error(yerr%can_co2,yscal%can_co2)
   write(unit=*,fmt=onefmt) 'CAN_CO2:',thiserr,ytemp%can_co2,y%can_co2,yerr%can_co2        &
                                      ,yscal%can_co2,troublemaker

   thiserr      = abs(yerr%lake_temp/yscal%lake_temp)
   errmax       = max(errmax,thiserr)
   troublemaker = large_error(yerr%lake_temp,yscal%lake_temp)
   write(unit=*,fmt=onefmt) 'LAKE_TEMP:',thiserr,ytemp%lake_temp,y%lake_temp               &
                                        ,yerr%lake_temp,yscal%lake_temp,troublemaker

   write(unit=*,fmt='(a)'  ) 
   write(unit=*,fmt='(80a)') ('=',k=1,80)
   write(unit=*,fmt='(a)'  ) 

   return
end subroutine print_lake_errmax
!==========================================================================================!
!==========================================================================================!
