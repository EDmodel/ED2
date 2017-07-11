!==========================================================================================!
!==========================================================================================!
!    This sub-routine estimates the soil respiration, based on ED-2.2 for most of it.  A   !
! few differences:                                                                         !
!                                                                                          !
!   Root biomass is defined as average LAI/SLA, where SLA follows the CLM-4 profile.       !
!                                                                                          !
!   Reference heterotrophic respiration is defined following:                              !
!                                                                                          !
!   Reichstein, M. et al. 2003: Modeling temporal and large-scale spatial variability of   !
!      soil respiration from soil water availability, temperature, and vegetation          !
!      productivity indices.  Glob. Biogeochem. Cycles, 117(4), 1104                       !
!      doi:10.1029/2003GB002035                                                            !
!                                                                                          !
!   exception is for bare ground, lakes, glaciers, where heterotrophic respiration is      !
!   set to zero.                                                                           !
!------------------------------------------------------------------------------------------!
subroutine leaf3_soil_resp(mzg,resolvable,soil_energy,soil_water,soil_text,leaf_class      &
                          ,veg_lai,veg_tai)
   use mem_leaf     , only : slz                      & ! intent(in)
                           , isfcl                    ! ! intent(in)
   use rconstants   , only : wdns                     & ! intent(in)
                           , umols_2_kgCyr            & ! intent(in)
                           , t00                      & ! intent(in)
                           , lnexp_min                & ! intent(in)
                           , lnexp_max                & ! intent(in)
                           , tiny_num8                ! ! intent(in)
   use leaf_coms    , only : kroot                    & ! intent(in)
                           , dslz                     & ! intent(in)
                           , sla_0                    & ! intent(in)
                           , sla_m                    & ! intent(in)
                           , tai_max                  & ! intent(in)
                           , sai                      & ! intent(in)
                           , rhmax_0                  & ! intent(in)
                           , rhmax_m                  & ! intent(in)
                           , slcpd                    & ! intent(in)
                           , slmsts                   & ! intent(in)
                           , soilcp                   & ! intent(in)
                           , k_hetresp                & ! intent(in)
                           , rr0_0                    & ! intent(in)
                           , rr0_qten                 & ! intent(in)
                           , rr0_dec                  & ! intent(in)
                           , decay_low_rh             & ! intent(in)
                           , decay_high_rh            & ! intent(in)
                           , low_temp_rh              & ! intent(in)
                           , high_temp_rh             & ! intent(in)
                           , decay_dry_rh             & ! intent(in)
                           , decay_wet_rh             & ! intent(in)
                           , dry_smoist_rh            & ! intent(in)
                           , wet_smoist_rh            & ! intent(in)
                           , phys_low_temp            & ! intent(in)
                           , phys_high_temp           & ! intent(in)
                           , tai_max                  & ! intent(in)
                           , sai                      & ! intent(in)
                           , tai_min                  & ! intent(in)
                           , soil_tempk               & ! intent(in)
                           , soil_fracliq             & ! intent(in)
                           , root_resp_o              & ! intent(out)
                           , het_resp_o               ! ! intent(out)
   use leaf3_physiol, only : leaf3_qten_fun           ! ! intent(in)
   use therm_lib    , only : uextcm2tl                ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)  :: mzg
   logical                     , intent(in)  :: resolvable
   real(kind=4), dimension(mzg), intent(in)  :: soil_energy
   real(kind=4), dimension(mzg), intent(in)  :: soil_water
   real(kind=4), dimension(mzg), intent(in)  :: soil_text
   real(kind=4)                , intent(in)  :: leaf_class
   real(kind=4)                , intent(in)  :: veg_lai
   real(kind=4)                , intent(in)  :: veg_tai
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: k
   integer                                   :: lrl
   integer                                   :: nveg
   integer                                   :: nsoil
   real(kind=4)                              :: Lc
   real(kind=4)                              :: rel_soil_moist
   real(kind=4)                              :: sum_soil_energy
   real(kind=4)                              :: sum_soil_water
   real(kind=4)                              :: sum_soil_hcap
   real(kind=4)                              :: sum_soil_slmsts
   real(kind=4)                              :: sum_soil_soilcp
   real(kind=4)                              :: avg_soil_temp
   real(kind=4)                              :: avg_soil_fliq
   real(kind=4)                              :: broot
   real(kind=4)                              :: rhmax
   real(kind=4)                              :: lnexplow
   real(kind=4)                              :: lnexphigh
   real(kind=4)                              :: tlow
   real(kind=4)                              :: thigh
   real(kind=4)                              :: tlow_fun
   real(kind=4)                              :: thigh_fun
   real(kind=4)                              :: temperature_limitation
   real(kind=4)                              :: water_limitation
   real(kind=4)                              :: lnexpdry
   real(kind=4)                              :: lnexpwet
   real(kind=4)                              :: smdry_fun
   real(kind=4)                              :: smwet_fun
   real(kind=8)                              :: rr0_08
   real(kind=8)                              :: rr0_qten8
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)                , external    :: sngloff
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Return with dummy output if this is not a "LEAF-4" run.                           !
   !---------------------------------------------------------------------------------------!
   if (isfcl /= 4) then
      root_resp_o = 0.
      het_resp_o  = 0.
      return
   end if
   !---------------------------------------------------------------------------------------!



   !----- Aliases to some useful variables. -----------------------------------------------!
   nveg      = nint(leaf_class)
   lrl       = kroot(nveg)
   tlow      = phys_low_temp (nveg) + t00
   thigh     = phys_high_temp(nveg) + t00
   rr0_08    = dble(rr0_0    (nveg))
   rr0_qten8 = dble(rr0_qten (nveg))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Solve respiration only if the vegetation is being resolved.                       !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then


      !----- Find the root biomass and the reference heterotrophic respiration. -----------!
      if (sla_m(nveg) == 0.) then
         broot = veg_lai / sla_0(nveg)
      else
         broot = log(1.0 + sla_m(nveg) * veg_lai / sla_0(nveg) ) / sla_m(nveg)
      end if
      rhmax = rhmax_0 + rhmax_m * ( tai_max(nveg) - sai(nveg) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Add "intensive" contribution of each layer, assuming that the roots are        !
      ! equally spread throughout the entire depth.                                        !
      !------------------------------------------------------------------------------------!
      root_resp_o = 0.0
      do k = lrl,mzg
         !---------------------------------------------------------------------------------!
         !    Compute the functions that will control the Rrf function for low and high    !
         ! temperature.  In order to avoid floating point exceptions, we check whether the !
         ! temperature will make the exponential too large or too small.                   !
         !---------------------------------------------------------------------------------!
         !----- Low temperature. ----------------------------------------------------------!
         lnexplow  = rr0_dec(nveg) * (tlow  - soil_tempk(k))
         lnexplow  = max(lnexp_min,min(lnexp_max,lnexplow))
         tlow_fun  = 1.0 +  exp(lnexplow)
         !----- High temperature. ---------------------------------------------------------!
         lnexphigh = rr0_dec(nveg) * (soil_tempk(k) - thigh)
         lnexphigh = max(lnexp_min,min(lnexp_max,lnexphigh))
         thigh_fun = 1.0 + exp(lnexphigh)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Decide which functional form to use based on the physiology.  This is just   !
         ! to make it look similar to the leaf respiration respiration.                    !
         !---------------------------------------------------------------------------------!
         root_resp_o = root_resp_o                                                         &
                     + sngloff(leaf3_qten_fun(dble(soil_tempk(k)),rr0_08,rr0_qten8)        &
                              , tiny_num8) / (tlow_fun * thigh_fun) * dslz(k)
         !---------------------------------------------------------------------------------!
      end do
      root_resp_o = root_resp_o * broot / abs(slz(lrl))
      !------------------------------------------------------------------------------------!
   else
      !----- Force respiration to be zero if the vegetation is not resolved. --------------!
      root_resp_o = 0.0
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     For heterotrophic respiration, check whether this is a place where things grow,   !
   ! otherwise, set it to zero.                                                            !
   !---------------------------------------------------------------------------------------!
   if ((tai_max(nveg)-sai(nveg)) > tai_min) then

      !------------------------------------------------------------------------------------!
      !     Integrate the soil extensive properties, plus the minimum and maximum possible !
      ! soil water content of the active layer.                                            !
      !------------------------------------------------------------------------------------!
      sum_soil_energy = 0.0
      sum_soil_hcap   = 0.0
      sum_soil_water  = 0.0
      sum_soil_slmsts = 0.0
      sum_soil_soilcp = 0.0
      do k = k_hetresp,mzg
         nsoil = nint(soil_text(k))

         !---------------------------------------------------------------------------------!
         !    Convert the units so energy is in J/m2, heat capacity in J/m2/K, and water   !
         ! in kg/m2.                                                                       !
         !---------------------------------------------------------------------------------!
         sum_soil_energy = sum_soil_energy + soil_energy(k)           * dslz(k)
         sum_soil_hcap   = sum_soil_hcap   + slcpd     (nsoil)        * dslz(k)
         sum_soil_water  = sum_soil_water  + soil_water(k)     * wdns * dslz(k)
         sum_soil_slmsts = sum_soil_slmsts + slmsts    (nsoil) * wdns * dslz(k)
         sum_soil_soilcp = sum_soil_soilcp + soilcp    (nsoil) * wdns * dslz(k)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !----- Find the average temperature and the relative soil moisture. -----------------!
      call uextcm2tl(sum_soil_energy,sum_soil_water,sum_soil_hcap,avg_soil_temp            &
                    ,avg_soil_fliq)
      rel_soil_moist = min( 1.0, max(0.0, ( sum_soil_water  - sum_soil_soilcp )            &
                                        / ( sum_soil_slmsts - sum_soil_soilcp ) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Similar to the original ED-1.0 formulation, which is based on the CENTURY     !
      ! model.  The change in the functional form is to avoid power of negative numbers,   !
      ! but the coefficients were tuned to give a similar curve.                           !
      !------------------------------------------------------------------------------------!
      !----- Low temperature limitation. --------------------------------------------------!
      lnexplow               = decay_low_rh * (low_temp_rh - avg_soil_temp)
      lnexplow               = max(lnexp_min,min(lnexp_max,lnexplow))
      tlow_fun               = 1.0 + exp(lnexplow)
      !----- High temperature limitation. -------------------------------------------------!
      lnexphigh              = decay_high_rh * (avg_soil_temp - high_temp_rh)
      lnexphigh              = max(lnexp_min,min(lnexp_max,lnexphigh))
      thigh_fun              = 1.0 + exp(lnexphigh)
      !----- Temperature limitation is a combination of both. -----------------------------!
      temperature_limitation = 1.0 / (tlow_fun * thigh_fun )
      !----- Dry soil limitation. ---------------------------------------------------------!
      lnexpdry               = decay_dry_rh * (dry_smoist_rh - rel_soil_moist)
      lnexpdry               = max(lnexp_min,min(lnexp_max,lnexpdry))
      smdry_fun              = 1.0 + exp(lnexpdry)
      !----- Wet soil limitation. ---------------------------------------------------------!
      lnexpwet               = decay_wet_rh * (rel_soil_moist - wet_smoist_rh)
      lnexpwet               = max(lnexp_min,min(lnexp_max,lnexpwet))
      smwet_fun              = 1.0 + exp(lnexpwet)
      !----- Soil moisture limitation is a combination of both. ---------------------------!
      water_limitation       = 1.0 / (smdry_fun * smwet_fun)
      !------------------------------------------------------------------------------------!



      !----- Find the average temperature and the relative soil moisture. -----------------!
      het_resp_o             = rhmax * temperature_limitation * water_limitation
      !------------------------------------------------------------------------------------!
   else
      !----- Set heterotrophic respiration to zero. ---------------------------------------!
      het_resp_o             = 0.0
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!


   return
end subroutine leaf3_soil_resp
!==========================================================================================!
!==========================================================================================!
