!==========================================================================================!
!==========================================================================================!
!     This routine is called by the radiation scheme and by LEAF.  It computes net surface !
! albedo plus radiative exchange between the atmosphere, vegetation, and the snow/ground   !
! given previously computed downward longwave and shortwave fluxes from the atmosphere.    !
! Also it finds functions of snowcover that are required for the above radiation           !
! calculations as well as other calculations in LEAF.                                      !
!     The shortwave parametrisation is only valid if the cosine of the zenith angle is     !
! greater than .03.  If running LEAF-4, it uses a two-stream radiation scheme largely      !
! based on CLM-4 and ED-2.2.                                                               !
!------------------------------------------------------------------------------------------!
subroutine leaf3_sfcrad(mzg,mzs,ip,soil_water,soil_color,soil_text,sfcwater_depth          &
                       ,patch_area,veg_fracarea,leaf_class,veg_lai,veg_tai                 &
                       ,veg_albedo_default,sfcwater_nlev,par_beam,par_diffuse,nir_beam     &
                       ,nir_diffuse,rshort,rshort_diffuse,rlong,cosz,albedt,rlongup        &
                       ,rshort_gnd,rlong_gnd,initial)
   use mem_leaf     , only : isfcl              & ! intent(in)
                           , albedo             ! ! intent(in)
   use leaf_coms    , only : tai_min            & ! intent(in)
                           , tai_min8           & ! intent(in)
                           , veg_frac           & ! intent(in)
                           , alb_oc_inter       & ! intent(in)
                           , alb_oc_slope       & ! intent(in)
                           , alb_oc_min         & ! intent(in)
                           , alb_oc_max         & ! intent(in)
                           , emiss_oc           & ! intent(in)
                           , slmsts             & ! intent(in)
                           , alb_nir_dry        & ! intent(in)
                           , alb_nir_wet        & ! intent(in)
                           , alb_vis_dry        & ! intent(in)
                           , alb_vis_wet        & ! intent(in)
                           , emisg              & ! intent(in)
                           , alb_damp           & ! intent(in)
                           , alb_snow           & ! intent(in)
                           , alb_snow_par       & ! intent(in)
                           , alb_snow_nir       & ! intent(in)
                           , emiss_snow         & ! intent(in)
                           , clumping_factor    & ! intent(in)
                           , leaf_scatter_vis   & ! intent(in)
                           , leaf_scatter_nir   & ! intent(in)
                           , wood_scatter_vis   & ! intent(in)
                           , wood_scatter_nir   & ! intent(in)
                           , phi1               & ! intent(in)
                           , phi2               & ! intent(in)
                           , emisv              & ! intent(in)
                           , sla_0              & ! intent(in)
                           , sla_m              & ! intent(in)
                           , vm0_0              & ! intent(in)
                           , dr_gamma           & ! intent(in)
                           , g_urban            & ! intent(in)
                           , alb_town           & ! intent(in)
                           , emis_town          & ! intent(in)
                           , snowfac            & ! intent(in)
                           , sfcwater_tempk     & ! intent(in)
                           , sfcwater_fracliq   & ! intent(in)
                           , soil_tempk         & ! intent(in)
                           , can_temp           & ! intent(in)
                           , veg_temp           & ! intent(in)
                           , ts_town            & ! intent(in)
                           , rshort_s           & ! intent(out)
                           , rshort_g           & ! intent(out)
                           , rshort_v           & ! intent(out)
                           , rlong_v            & ! intent(out)
                           , rlong_g            & ! intent(out)
                           , rlong_s            & ! intent(out)
                           , lai_ss             & ! intent(out)
                           , par_l_ss           & ! intent(out)
                           , sla_ss             & ! intent(out)
                           , vm0_ss             & ! intent(out)
                           , rd0_ss             ! ! intent(out)
   use rconstants   , only : stefan             & ! intent(in)
                           , lnexp_max          & ! intent(in)
                           , tiny_num8          ! ! intent(in)
   use mem_radiate  , only : rad_cosz_min       & ! intent(in)
                           , rad_rlong_min      ! ! intent(in)
   use therm_lib    , only : uextcm2tl          & ! subroutine
                           , idealdenssh        ! ! function
   use catt_start   , only : CATT               ! ! intent(in)
   use teb_spm_start, only : TEB_SPM            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in)     :: mzg
   integer                , intent(in)     :: mzs
   integer                , intent(in)     :: ip
   real   , dimension(mzg), intent(in)     :: soil_water
   real                   , intent(in)     :: soil_color
   real   , dimension(mzg), intent(in)     :: soil_text
   real   , dimension(mzs), intent(in)     :: sfcwater_depth 
   real                   , intent(in)     :: patch_area
   real                   , intent(in)     :: veg_fracarea
   real                   , intent(in)     :: leaf_class
   real                   , intent(in)     :: veg_lai
   real                   , intent(in)     :: veg_tai
   real                   , intent(in)     :: veg_albedo_default
   real                   , intent(in)     :: sfcwater_nlev
   real                   , intent(in)     :: par_beam
   real                   , intent(in)     :: par_diffuse
   real                   , intent(in)     :: nir_beam
   real                   , intent(in)     :: nir_diffuse
   real                   , intent(in)     :: rshort
   real                   , intent(in)     :: rshort_diffuse
   real                   , intent(in)     :: rlong
   real                   , intent(in)     :: cosz
   real                   , intent(inout)  :: albedt
   real                   , intent(inout)  :: rlongup
   real                   , intent(inout)  :: rshort_gnd
   real                   , intent(inout)  :: rlong_gnd
   logical                , intent(in)     :: initial
   !----- Local variables. ----------------------------------------------------------------!
   integer                                 :: nsoil
   integer                                 :: ncolour
   integer                                 :: k
   integer                                 :: ksn
   integer                                 :: nveg
   real(kind=4)                            :: fcpct
   real(kind=4)                            :: albedo_soil_par
   real(kind=4)                            :: albedo_soil_nir
   real(kind=4)                            :: albedo_soil
   real(kind=4)                            :: emiss_soil
   real(kind=4)                            :: albedo_damp_par
   real(kind=4)                            :: albedo_damp_nir
   real(kind=4)                            :: albedo_damp
   real(kind=4)                            :: albedo_sfcw_par
   real(kind=4)                            :: albedo_sfcw_nir
   real(kind=4)                            :: albedo_sfcw
   real(kind=4)                            :: rad_sfcw_par
   real(kind=4)                            :: rad_sfcw_nir
   real(kind=4)                            :: rad_sfcw
   real(kind=4)                            :: sfcw_odepth
   real(kind=4)                            :: fractrans_par
   real(kind=4)                            :: fractrans_nir
   real(kind=4)                            :: fractrans
   real(kind=4)           , dimension(mzs) :: abs_sfcw_par
   real(kind=4)           , dimension(mzs) :: abs_sfcw_nir
   real(kind=4)           , dimension(mzs) :: abs_sfcw
   real(kind=4)           , dimension(mzs) :: abs_sfcw_par_diffuse
   real(kind=4)           , dimension(mzs) :: abs_sfcw_par_beam
   real(kind=4)           , dimension(mzs) :: abs_sfcw_nir_diffuse
   real(kind=4)           , dimension(mzs) :: abs_sfcw_nir_beam
   real(kind=4)                            :: abs_ground_par
   real(kind=4)                            :: abs_ground_nir
   real(kind=4)                            :: abs_ground
   real(kind=4)                            :: albedo_ground_par
   real(kind=4)                            :: albedo_ground_nir
   real(kind=4)                            :: albedo_ground
   real(kind=4)                            :: albedo_veg_par
   real(kind=4)                            :: albedo_veg_nir
   real(kind=4)                            :: albedo_veg
   real(kind=4)                            :: emiss_veg
   real(kind=4)                            :: f_exposed
   real(kind=4)                            :: f_covered
   real(kind=4)                            :: rlost_veg
   real(kind=4)                            :: rlost_gs
   real(kind=4)                            :: downward_par_below_beam
   real(kind=4)                            :: downward_nir_below_beam
   real(kind=8)                            :: leaf_temp8
   real(kind=8)                            :: wood_temp8
   real(kind=8)                            :: lai8
   real(kind=8)                            :: wai8
   real(kind=8)                            :: tai8
   real(kind=4)                            :: par_v_beam
   real(kind=4)                            :: rshort_v_beam
   real(kind=4)                            :: par_v_diffuse
   real(kind=4)                            :: rshort_v_diffuse
   real(kind=4)                            :: lw_v
   real(kind=4)                            :: downward_par_below_diffuse
   real(kind=4)                            :: upward_par_above_diffuse
   real(kind=4)                            :: downward_nir_below_diffuse
   real(kind=4)                            :: upward_nir_above_diffuse 
   real(kind=4)                            :: T_surface
   real(kind=4)                            :: emissivity
   real(kind=4)                            :: par_beam_norm
   real(kind=4)                            :: par_diff_norm
   real(kind=4)                            :: nir_beam_norm
   real(kind=4)                            :: nir_diff_norm
   real(kind=4)                            :: sum_norm
   real(kind=4)                            :: downward_lw_below
   real(kind=4)                            :: upward_lw_below
   real(kind=4)                            :: upward_lw_above
   real(kind=4)                            :: downward_rshort_below_beam
   real(kind=4)                            :: downward_rshort_below_diffuse
   real(kind=4)                            :: upward_rshort_above_diffuse
   real(kind=4)                            :: surface_netabs_longwave
   real(kind=4)                            :: wleaf_vis
   real(kind=4)                            :: kl
   real(kind=4)                            :: ekl
   real(kind=4)                            :: kay
   real(kind=4)                            :: slai
   real(kind=4)                            :: stai
   real(kind=4)                            :: par_l_beam
   real(kind=4)                            :: par_l_diffuse
   real(kind=4)                            :: par_l
   real(kind=4)                            :: f_sun
   real(kind=4)                            :: f_shade
   real(kind=4)                            :: rlonga_v
   real(kind=4)                            :: rlonga_gs
   real(kind=4)                            :: rlongv_gs
   real(kind=4)                            :: rlongv_a
   real(kind=4)                            :: rlonggs_v
   real(kind=4)                            :: rlonggs_a
   real(kind=4)                            :: rlonga_a
   !----- External function. --------------------------------------------------------------!
   real(kind=4)           , external       :: sngloff
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Aliases to snow fraction corrected of veg_lai and veg_tai.                        !
   !---------------------------------------------------------------------------------------!
   slai = veg_lai * ( 1.0 - snowfac )
   stai = veg_tai * ( 1.0 - snowfac )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Decide which solution to use.                                                     !
   !---------------------------------------------------------------------------------------!
   select case (ip)
   case (1)
      !----- Albedo must be calculated during daytime. ------------------------------------!
      if (cosz > rad_cosz_min) then
         albedo_ground = alb_oc_inter + alb_oc_slope * tan(acos(cosz))
         albedo_ground = min( max( albedo_ground, alb_oc_min), alb_oc_max)
         albedt        = albedt + patch_area * albedo_ground
      else
         albedo_ground = 0.0
      end if
      rshort_gnd       = (1.0 - albedo_ground) * rshort
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Longwave radiation.  Assume black body unless this is a LEAF-4 run.            !
      !------------------------------------------------------------------------------------!
      if ( rlong >= rad_rlong_min ) then
         select case (isfcl)
         case (4,5)
            !------ Use a similar approach as CLM-4. --------------------------------------!
            upward_lw_above = emiss_oc * stefan * soil_tempk(mzg) ** 4
            rlong_gnd       = (1.0 - emiss_oc     ) * rlong
            !------------------------------------------------------------------------------!
         case default
            !------ Oceans are black bodies for longwave. ---------------------------------!
            upward_lw_above = stefan * soil_tempk(mzg) ** 4
            rlong_gnd       = 0.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      else
         !------ First time step, set emissions to zero. ----------------------------------!
         upward_lw_above = stefan * soil_tempk(mzg) ** 4
         rlong_gnd       = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !----- Update longwave radiation. ---------------------------------------------------!
      rlongup         = rlongup + patch_area * upward_lw_above
      !------------------------------------------------------------------------------------!


   case default
      !------------------------------------------------------------------------------------!
      !     Likely to be a land patch, decide how to solve based on the model used.        !
      !------------------------------------------------------------------------------------!
      select case (isfcl)
      case (0)
         !----- No vegetation model, use prescribed values and wish the user good luck. ---!
         albedt     = albedt  + patch_area * albedo
         rlongup    = rlongup + patch_area * stefan * can_temp ** 4
         rshort_gnd = (1.0 - albedo) * rshort
         rlong_gnd  = 0.0
         !---------------------------------------------------------------------------------!
      case (1,2)
         !---------------------------------------------------------------------------------!
         !      LEAF 3, simple radiation model.                                            !
         !---------------------------------------------------------------------------------!


         !------ Diagnose snow temperature and the influence of snow covering veg. --------!
         nveg  = nint(leaf_class)
         nsoil = nint(soil_text(mzg))
         ksn   = nint(sfcwater_nlev)
         !---------------------------------------------------------------------------------!


         !------ Defining the exposed area. -----------------------------------------------!
         f_covered = veg_fracarea * (1. - snowfac)
         f_exposed = 1. - f_covered
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Ground albedo.  Soil colour is based on CLM-4.  For now we cannot initial-  !
         ! ise the colour with maps, though.                                               !
         !---------------------------------------------------------------------------------!
         select case (nsoil)
         case (13)
            !----- Bedrock, use constants soil value for granite. -------------------------!
            albedo_soil = 0.32
            !------------------------------------------------------------------------------!
         case (12)
            !----- Peat, follow McCumber and Pielke (1981). -------------------------------!
            fcpct       = soil_water(mzg) / slmsts(nsoil)
            albedo_soil = max (0.07, 0.14 * (1.0 - fcpct))
            !------------------------------------------------------------------------------!
         case default
            !----- Other soils. -----------------------------------------------------------!
            ncolour = nint(soil_color)
            select case (ncolour)
            case (21)
               fcpct       = soil_water(mzg) / slmsts(nsoil)
               albedo_soil = max(0.14, 0.31 - 0.34 * fcpct)
            case default
               fcpct       = max (0.00, 0.11 - 0.40 * soil_water(mzg))
               albedo_soil = min ( 0.5 * (alb_nir_dry(ncolour) + alb_vis_dry(ncolour))     &
                                 , 0.5 * (alb_nir_wet(ncolour) + alb_vis_wet(ncolour))     &
                                 + fcpct )
            end select
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Vegetation albedo.                                                         !
         !---------------------------------------------------------------------------------!
         albedo_veg = veg_albedo_default
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !       Snow/surface water albedo.                                                !
         !---------------------------------------------------------------------------------!
         rad_sfcw      = 1.0
         albedo_ground = 1.0
         abs_ground    = 0.0
         emiss_veg     = emisv(nveg )
         emiss_soil    = emisg(nsoil)
         if (ksn == 0) then
            emissivity = emiss_soil
            T_surface  = soil_tempk(mzg)
         else
            !------------------------------------------------------------------------------!
            !      Sfcwater albedo ranges from wet-soil value for all-liquid to typical    !
            ! snow albedo for ice.                                                         !
            !------------------------------------------------------------------------------!
            albedo_sfcw = alb_damp + sfcwater_fracliq(ksn) * ( alb_snow - alb_damp )
            !------------------------------------------------------------------------------!



            !----- Fraction shortwave absorbed into sfcwater + soil. ----------------------!
            rad_sfcw = 1.0 - albedo_sfcw
            do k = ksn,1,-1
               !---------------------------------------------------------------------------!
               !      Fractrans is fraction of shortwave entering each sfcwater layer that !
               ! gets transmitted through that layer.                                      !
               !---------------------------------------------------------------------------!
               sfcw_odepth = min( lnexp_max,   20.0 * sfcwater_depth(k) )
               fractrans   = exp( - sfcw_odepth )
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      abs_sfcw(k) is fraction of total incident shortwave (at top of top   !
               ! sfcwater layer) that is absorbed in each sfcwater layer.                  !
               !---------------------------------------------------------------------------!
               abs_sfcw(k) = rad_sfcw * (1.0 - fractrans)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Rad is fraction of total incident shortwave (at top of top sfcwater  !
               ! layer) that remains at bottom of current sfcwater layer.                  !
               !---------------------------------------------------------------------------!
               rad_sfcw = rad_sfcw * fractrans
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Albedo_ground will ultimately be the albedo of the soil+sfcwater.    !
               ! So subtract out whatever is being absorbed by sfcwater.                   !
               !---------------------------------------------------------------------------!
               albedo_ground = albedo_ground - abs_sfcw(k)
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!

            !----- Long wave parameter if sfcwater exists. --------------------------------!
            emissivity = emiss_snow * snowfac + emiss_soil * (1.0 - snowfac)
            T_surface  = sqrt(sqrt( ( sfcwater_tempk(ksn)** 4 * emiss_snow*snowfac         &
                                    + soil_tempk    (mzg)** 4 * emiss_soil*(1.0-snowfac) ) &
                                  / emissivity ) )
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Get the total radiation absorbed by ground.                                 !
         !---------------------------------------------------------------------------------!
         abs_ground = (1.0 - albedo_soil) * (1.0 - snowfac + snowfac * rad_sfcw )
         !---------------------------------------------------------------------------------!




         !----- Subtract off ground absorption to obtain the soil+sfcwater albedo. --------!
         albedo_ground = albedo_ground - abs_ground
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Scale absorption by each layer by the fraction of landscape covered or     !
         ! exposed.                                                                        !
         !---------------------------------------------------------------------------------!
         do k=1,ksn
            rshort_s(k) = rshort * f_exposed * abs_sfcw(k)
         end do
         do k=ksn+1,mzs
            rshort_s(k) = 0.0
         end do
         rshort_g = rshort * f_exposed * abs_ground
         rshort_v = rshort * f_covered * (1. - albedo_veg + f_exposed * albedo_ground)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the net albedo for this patch and update total albedo.                 !
         !---------------------------------------------------------------------------------!
         albedo_ground = f_covered * albedo_veg + f_exposed * f_exposed * albedo_ground
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Add urban contribution if running TEB.                                     !
         !---------------------------------------------------------------------------------!
         if (teb_spm==1) then
            if (nint(g_urban) == 0) then
               albedt = albedt + patch_area * albedo_ground
            else
               albedt = albedt + patch_area * alb_town
            endif
         else
            albedt = albedt + patch_area * albedo_ground
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Longwave radiation calculations.                                            !
         !---------------------------------------------------------------------------------!
         rlost_gs  = emissivity * stefan * T_surface ** 4
         rlost_veg = emiss_veg  * stefan * veg_temp ** 4
         !----- Flux components. ----------------------------------------------------------!
         rlonga_v  = rlong     * f_covered * (emiss_veg + f_exposed * (1. - emissivity))
         rlonga_gs = rlong     * f_exposed * emissivity
         rlongv_gs = rlost_veg * f_covered * emissivity
         rlongv_a  = rlost_veg * f_covered                                                 &
                   * (2. - emissivity - f_covered + emissivity * f_covered)
         rlonggs_v = rlost_gs  * f_covered * emiss_veg
         rlonggs_a = rlost_gs  * f_exposed
         rlonga_a  = rlong     * ( f_covered * (1. - emiss_veg)                            &
                                 + f_exposed * f_exposed * (1. - emissivity) )
         !---------------------------------------------------------------------------------!



         !----- Add urban contribution if running TEB. ------------------------------------!
         if ( teb_spm == 1 ) then
            if (nint(g_urban) == 0) then
               upward_lw_above = rlongv_a + rlonggs_a + rlonga_a
            else
               upward_lw_above = emis_town * stefan * ts_town**4
            end if
         else
            upward_lw_above = rlongv_a + rlonggs_a + rlonga_a
         end if
         rlongup = rlongup + patch_area * upward_lw_above
         !---------------------------------------------------------------------------------!



         !----- Integrate the total absorbed light by ground (top soil plus TSW layers). --!
         rshort_gnd = rshort_g
         do k=1,ksn
            rshort_gnd = rshort_gnd + rshort_s(k)
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      In case rlong is not computed, zero out all longwave fluxes other than     !
         ! rlongup.  [On the first timestep, radiative fluxes may not be available until   !
         ! microphysics is called, and zeroing these fluxes avoids the imbalance of having !
         ! upward longwave without downward longwave in LEAF-3.  Also, this allows LEAF-3  !
         ! to run without radiation for all timesteps, if desired for testing.].           !
         !---------------------------------------------------------------------------------!
         if (rlong >= rad_rlong_min) then
            rlong_v   = rlonga_v  + rlonggs_v - rlongv_a  - rlongv_gs
            if (ksn > 0) then
               rlong_s  = rlonga_gs + rlongv_gs - rlonggs_a - rlonggs_a
               rlong_g  = 0.0
            else
               rlong_s  = 0.0
               rlong_g  = rlonga_gs + rlongv_gs - rlonggs_a - rlonggs_a
            end if
            rlong_gnd = rlonga_gs + rlongv_gs - rlonggs_a - rlonggs_a
         else
            rlong_v   = 0.0
            rlong_g   = 0.0
            rlong_s   = 0.0
            rlong_gnd = 0.0
         end if
         !---------------------------------------------------------------------------------!

      case (4,5)
         !---------------------------------------------------------------------------------!
         !      LEAF 4, use a two stream model.  ED-2.2 calls this only once, at the       !
         ! initialisation.                                                                 !
         !---------------------------------------------------------------------------------!


         !------ Aliases to some useful indices. ------------------------------------------!
         nveg    = nint(leaf_class)
         nsoil   = nint(soil_text(mzg))
         ksn     = nint(sfcwater_nlev)
         ncolour = nint(soil_color)
         !---------------------------------------------------------------------------------!


         !------ Covered and exposed areas. -----------------------------------------------!
         f_covered = veg_frac(nveg)
         f_exposed = 1. - f_covered
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Initialise the variables.                                                   !
         !---------------------------------------------------------------------------------!
         leaf_temp8           = dble(veg_temp)
         wood_temp8           = dble(veg_temp)
         if ( slai >= tai_min) then
            lai8              = dble(slai) / f_covered
         else
            lai8              = 0.d0
         end if
         if ( ( stai - slai ) >= tai_min) then
            wai8              = ( dble(stai) - dble(slai) ) / f_covered
         else
            wai8              = 0.d0
         end if
         tai8                 = lai8 + wai8
         par_v_beam           = 0.
         rshort_v_beam        = 0.
         par_v_diffuse        = 0.
         rshort_v_diffuse     = 0.
         lw_v                 = 0.
         abs_sfcw_par_diffuse = 0.
         abs_sfcw_par_beam    = 0.
         abs_sfcw_nir_diffuse = 0.
         abs_sfcw_nir_beam    = 0.
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Find the normalisation factors, but only when it is daytime (rshort > 0).  !
         !---------------------------------------------------------------------------------!
         if (cosz > rad_cosz_min .and. rshort > 0.) then
            par_beam_norm = max( 1.e-5 , par_beam    / rshort )
            par_diff_norm = max( 1.e-5 , par_diffuse / rshort )
            nir_beam_norm = max( 1.e-5 , nir_beam    / rshort )
            nir_diff_norm = max( 1.e-5 , nir_diffuse / rshort )
            sum_norm      = par_beam_norm + par_diff_norm + nir_beam_norm + nir_diff_norm
         else 
            !------------------------------------------------------------------------------!
            !     Night-time, nothing will happen, split equally between the 4 components. !
            !------------------------------------------------------------------------------!
            par_beam_norm = 0.25
            par_diff_norm = 0.25
            nir_beam_norm = 0.25
            nir_diff_norm = 0.25
            sum_norm      = 1.00
         end if
         !---------------------------------------------------------------------------------!
         !     Because we must tweak the radiation so none of the terms are zero, we       !
         ! must correct the normalised radiation variables so they add up to one.          !
         !---------------------------------------------------------------------------------!
         par_beam_norm = par_beam_norm / sum_norm
         par_diff_norm = par_diff_norm / sum_norm
         nir_beam_norm = nir_beam_norm / sum_norm
         nir_diff_norm = nir_diff_norm / sum_norm
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the ground albedo as a function of soil water relative moisture of the !
         ! top layer.                                                                      !
         !---------------------------------------------------------------------------------!
         select case (nsoil)
         case (13)
            !----- Bedrock, use constants soil value for granite. -------------------------!
            albedo_soil_par = alb_vis_dry(ncolour)
            albedo_soil_nir = alb_nir_dry(ncolour)
            !----- Damp soil, for temporary surface water albedo. -------------------------!
            albedo_damp_par = albedo_soil_par
            albedo_damp_nir = albedo_damp_nir
         case (12)
            !----- Peat, follow McCumber and Pielke (1981). -------------------------------!
            fcpct           = soil_water(mzg) / slmsts(nsoil)
            albedo_soil_par = max (0.07, 0.14 * (1.0 - fcpct))
            albedo_soil_nir = albedo_soil_par
            !----- Damp soil, for temporary surface water albedo. -------------------------!
            albedo_damp_par = 0.14
            albedo_damp_nir = 0.14
         case default
            select case (ncolour)
            case (21)
               !---------------------------------------------------------------------------!
               !     ED-2.1 soil colour.  Also, we use the ED-2.1 default method to        !
               ! determine the albedo.                                                     !
               !---------------------------------------------------------------------------!
               fcpct           = soil_water(mzg) / slmsts(nsoil)
               albedo_soil_par = max(0.14,0.31-0.34*fcpct)
               albedo_soil_nir = albedo_soil_par
               !----- Damp soil, for temporary surface water albedo. ----------------------!
               albedo_damp_par = 0.14
               albedo_damp_nir = 0.14
            case default
               !---------------------------------------------------------------------------!
               !      Other soils, we use the soil numbers from CLM-4.  The colour class   !
               ! must be given at RAMSIN.  At this point the value is the same for all     !
               ! points, but in the future we may read their files if the results are      !
               ! promising.                                                                !
               !---------------------------------------------------------------------------!
               fcpct           = max(0., 0.11 - 0.40 * soil_water(mzg))
               albedo_soil_par = min(alb_vis_dry(ncolour), alb_vis_wet(ncolour)  + fcpct)
               albedo_soil_nir = min(alb_nir_dry(ncolour), alb_nir_wet(ncolour)  + fcpct)
               !----- Damp soil, for temporary surface water albedo. ----------------------!
               fcpct           = max(0., 0.11 - 0.40 * slmsts(ncolour))
               albedo_damp_par = min(alb_vis_dry(ncolour), alb_vis_wet(ncolour)  + fcpct)
               albedo_damp_nir = min(alb_nir_dry(ncolour), alb_nir_wet(ncolour)  + fcpct)
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Decide what is our surface temperature.  When the soil is exposed, then     !
         ! that is the surface temperature.  Otherwise, we pick the temporary surface      !
         ! water or snow layer.                                                            !
         !---------------------------------------------------------------------------------!
         rad_sfcw_par      = 1.0
         rad_sfcw_nir      = 1.0
         albedo_ground_par = 1.0
         albedo_ground_nir = 1.0
         abs_sfcw_par      = 0.0
         abs_sfcw_nir      = 0.0
         emiss_soil        = emisg(nsoil)
         if (ksn == 0) then
            emissivity = emiss_soil
            T_surface  = soil_tempk(mzg)
         else
            !------------------------------------------------------------------------------!
            !      Sfcwater albedo ALS ranges from wet-soil value for all-liquid to        !
            ! typical snow albedo for ice.  In the future, we should consider a more       !
            ! realistic snow albedo model that takes snow age and snow melt into account.  !
            !                                                                              !
            !  Potential starting points:                                                  !
            !                                                                              !
            !  Roesch, A., et al., 2002: Comparison of spectral surface albedos and their  !
            !      impact on the general circulation model simulated surface climate.  J.  !
            !      Geophys. Res.-Atmosph., 107(D14), 4221, 10.1029/2001JD000809.           !
            !                                                                              !
            !  Oleson, K.W., et al., 2010: Technical description of version 4.0 of the     !
            !      Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR.        !
            !                                                                              !
            !------------------------------------------------------------------------------!
            albedo_sfcw_par = albedo_damp_par + sfcwater_fracliq(ksn)                      &
                                              * ( alb_snow_par - albedo_damp_par )
            albedo_sfcw_nir = albedo_damp_nir + sfcwater_fracliq(ksn)                      &
                                              * ( alb_snow_nir - albedo_damp_nir )
            !------------------------------------------------------------------------------!



            !----- Fraction shortwave absorbed into sfcwater + soil. ----------------------!
            rad_sfcw_par = 1.0 - albedo_sfcw_par
            rad_sfcw_nir = 1.0 - albedo_sfcw_nir
            do k = ksn,1,-1
               !---------------------------------------------------------------------------!
               !      Fractrans is fraction of shortwave entering each sfcwater layer that !
               ! gets transmitted through that layer.                                      !
               !---------------------------------------------------------------------------!
               sfcw_odepth   = min( lnexp_max,   20.0 * sfcwater_depth(k))
               fractrans_par = exp( - sfcw_odepth )
               fractrans_nir = fractrans_par
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      abs_sfcw_???(k) is fraction of total incident shortwave (atop of     !
               ! top sfcwater layer) that is absorbed in each sfcwater layer.              !
               !---------------------------------------------------------------------------!
               abs_sfcw_par(k) = rad_sfcw_par * (1.0 - fractrans_par)
               abs_sfcw_nir(k) = rad_sfcw_nir * (1.0 - fractrans_nir)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Rad is fraction of total incident shortwave (at top of top sfcwater  !
               ! layer) that remains at bottom of current sfcwater layer.                  !
               !---------------------------------------------------------------------------!
               rad_sfcw_par = rad_sfcw_par * fractrans_par
               rad_sfcw_nir = rad_sfcw_nir * fractrans_nir
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Albedo_ground will ultimately be the albedo of the soil+sfcwater.    !
               ! So subtract out whatever is being absorbed by sfcwater.                   !
               !---------------------------------------------------------------------------!
               albedo_ground_par = albedo_ground_par - abs_sfcw_par(k)
               albedo_ground_nir = albedo_ground_nir - abs_sfcw_nir(k)
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!

            !----- Long wave parameter if sfcwater exists. --------------------------------!
            emissivity = emiss_snow * snowfac + emiss_soil * (1.0 - snowfac)
            T_surface  = sqrt(sqrt( ( sfcwater_tempk(ksn)** 4 * emiss_snow*snowfac         &
                                    + soil_tempk    (mzg)** 4 * emiss_soil*(1.0-snowfac) ) &
                                  / emissivity ) )
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     This is the fraction of below-canopy radiation that is absorbed by the      !
         ! ground.                                                                         !
         !---------------------------------------------------------------------------------!
         abs_ground_par = (1.0 - albedo_soil_par)                                          &
                        * (1.0 - snowfac + snowfac * rad_sfcw_par )
         abs_ground_nir = (1.0 - albedo_soil_nir)                                          &
                        * (1.0 - snowfac + snowfac * rad_sfcw_nir )
         !---------------------------------------------------------------------------------!




         !----- Subtract off ground absorption to obtain the soil+sfcwater albedo. --------!
         albedo_ground_par = albedo_ground_par - abs_ground_par
         albedo_ground_nir = albedo_ground_nir - abs_ground_nir
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Decide what to do based on plant area index (TAI).  If TAI is too small,   !
         ! then we may bypass it entirely.                                                 !
         !---------------------------------------------------------------------------------!
         if (tai8 >= tai_min8) then
            !------------------------------------------------------------------------------!
            !    Call the two-stream model for the longwave radiation.                     !
            !------------------------------------------------------------------------------!
            call leaf3_2stream_lw(emissivity,T_surface,rlong,nveg,lai8,wai8 ,leaf_temp8    &
                                 ,wood_temp8,lw_v,downward_lw_below,upward_lw_below        &
                                 ,upward_lw_above)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Correct flux by vegetation fraction area.                                !
            !------------------------------------------------------------------------------!
            lw_v              = lw_v * f_covered
            downward_lw_below = downward_lw_below * f_covered + rlong * f_exposed
            upward_lw_below   = upward_lw_below * f_covered                                &
                              + ( (1.0 - emissivity) * rlong                               &
                                + emissivity * stefan * T_surface**4 ) * f_exposed
            upward_lw_above   = upward_lw_above * f_covered                                &
                              + ( (1.0 - emissivity) * rlong                               &
                                + emissivity * stefan * T_surface**4 ) * f_exposed
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Compute short wave only if it is daytime.                                !
            !------------------------------------------------------------------------------!
            if (cosz > rad_cosz_min) then
               !---------------------------------------------------------------------------!
               !    Call the two-stream model for PAR and NIR.                             !
               !---------------------------------------------------------------------------!
               call leaf3_2stream_sw(albedo_ground_par,albedo_ground_nir,cosz              &
                                    ,par_beam_norm,par_diff_norm,nir_beam_norm             &
                                    ,nir_diff_norm,nveg,lai8,wai8,par_v_beam               &
                                    ,par_v_diffuse,rshort_v_beam,rshort_v_diffuse          &
                                    ,downward_par_below_beam,downward_par_below_diffuse    &
                                    ,upward_par_above_diffuse,downward_nir_below_beam      &
                                    ,downward_nir_below_diffuse,upward_nir_above_diffuse)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Scale by the covered/exposed area.                                   !
               !---------------------------------------------------------------------------!
               par_v_beam                    = par_v_beam                 * f_covered
               par_v_diffuse                 = par_v_diffuse              * f_covered
               rshort_v_beam                 = rshort_v_beam              * f_covered
               rshort_v_diffuse              = rshort_v_diffuse           * f_covered
               downward_par_below_beam       = downward_par_below_beam    * f_covered      &
                                             + par_beam_norm              * f_exposed
               downward_par_below_diffuse    = downward_par_below_diffuse * f_covered      &
                                             + par_diff_norm              * f_exposed
               downward_nir_below_beam       = downward_nir_below_beam    * f_covered      &
                                             + nir_beam_norm              * f_exposed
               downward_nir_below_diffuse    = downward_nir_below_diffuse * f_covered      &
                                             + nir_diff_norm              * f_exposed
               upward_par_above_diffuse      = upward_par_above_diffuse   * f_covered      &
                                             + albedo_ground_par          * f_exposed      &
                                             * ( par_diff_norm + par_beam_norm )
               upward_nir_above_diffuse      = upward_nir_above_diffuse   * f_covered      &
                                             + albedo_ground_nir          * f_exposed      &
                                             * ( nir_diff_norm + nir_beam_norm )
               !---------------------------------------------------------------------------!



               !----- Above-canopy upwelling radiation. -----------------------------------!
               upward_rshort_above_diffuse = upward_par_above_diffuse                      &
                                           + upward_nir_above_diffuse
               !---------------------------------------------------------------------------!



               !----- Below-canopy downwelling radiation. ---------------------------------!
               downward_rshort_below_beam    = downward_par_below_beam                     &
                                             + downward_nir_below_beam
               downward_rshort_below_diffuse = downward_par_below_diffuse                  &
                                             + downward_nir_below_diffuse
               !---------------------------------------------------------------------------!


            else
               !----- The code expects values for these, even when it is not daytime. -----!
               par_v_beam                      = 0.
               par_v_diffuse                   = 0.
               rshort_v_beam                   = 0.
               rshort_v_diffuse                = 0.
               downward_par_below_beam         = par_beam_norm
               downward_par_below_diffuse      = par_diff_norm
               downward_nir_below_beam         = nir_beam_norm
               downward_nir_below_diffuse      = nir_diff_norm
               downward_rshort_below_beam      = par_beam_norm + nir_beam_norm
               downward_rshort_below_diffuse   = par_diff_norm + nir_diff_norm
               upward_par_above_diffuse        = albedo_ground_par                         &
                                               * ( par_diff_norm + par_beam_norm )
               upward_nir_above_diffuse        = albedo_ground_nir                         &
                                               * ( nir_diff_norm + nir_beam_norm )
               upward_rshort_above_diffuse     = upward_par_above_diffuse                  &
                                               + upward_nir_above_diffuse
               !---------------------------------------------------------------------------!

            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !    Fraction of visible light that is absorbed by leaves.                     !
            !------------------------------------------------------------------------------!
            wleaf_vis = sngloff ( ( clumping_factor(nveg)                                  &
                                  * (1.d0 - leaf_scatter_vis(nveg)) * lai8 )               &
                                / ( clumping_factor(nveg)                                  &
                                  * (1.d0 - leaf_scatter_vis(nveg)) * lai8                 &
                                  + (1.d0 - wood_scatter_vis(nveg)) * wai8 ), tiny_num8 )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find the fraction of PAR radiation that is absorbed by leaves.            !
            !------------------------------------------------------------------------------!
            par_l_beam    = par_v_beam    * wleaf_vis
            par_l_diffuse = par_v_diffuse * wleaf_vis
            !------------------------------------------------------------------------------!

         else
            !----- This is the case where there is no vegetation. -------------------------!
            downward_par_below_beam       = par_beam_norm
            downward_par_below_diffuse    = par_diff_norm
            downward_nir_below_beam       = nir_beam_norm
            downward_nir_below_diffuse    = nir_diff_norm
            downward_rshort_below_beam    = par_beam_norm + nir_beam_norm
            downward_rshort_below_diffuse = par_diff_norm + nir_diff_norm

            upward_par_above_diffuse      = albedo_ground_par                              &
                                          * ( par_diff_norm + par_beam_norm )
            upward_nir_above_diffuse      = albedo_ground_nir                              &
                                          * ( nir_diff_norm + nir_beam_norm )
            upward_rshort_above_diffuse   = upward_par_above_diffuse                       &
                                          + upward_nir_above_diffuse


            !----- Make sure net absorption by vegetation is zero. ------------------------!
            lw_v              = 0.0
            downward_lw_below = rlong
            upward_lw_below   = ( 1.0 - emissivity ) * downward_lw_below                   &
                              + emissivity * stefan * T_surface**4
            upward_lw_above   = upward_lw_below
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Set leaf absorption to zero.                                              !
            !------------------------------------------------------------------------------!
            par_l_beam    = 0.0
            par_l_diffuse = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Use the CLM-4 equation to find the fraction of sun leaves, based on         !
         ! Thornston and Zimmermann (2007) method.  Like in CLM-4, ensure the fractions    !
         ! are not too small to avoid singularities.                                       !
         !---------------------------------------------------------------------------------!
         if (cosz > rad_cosz_min .and. slai >= tai_min) then

            !------------------------------------------------------------------------------!
            !    Find the optical depth, but make sure it isn't too thick.  In case it is, !
            ! assume a maximum optical depth that won't give trouble due to floating point !
            ! errors.                                                                      !
            !------------------------------------------------------------------------------!
            kl  = sngloff( ( phi1(nveg) + phi2(nveg) * cosz ) * slai / cosz, tiny_num8 )
            kl  = min(kl,lnexp_max)
            ekl = exp(-kl)
            kay = - log(ekl) / slai
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Use the CLM-4 equation to find the fraction of sun leaves.                !
            !------------------------------------------------------------------------------!
            f_sun   = (1.0 - ekl) / kl
            f_shade = 1.0 - f_sun
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    To be safe, make sure the sun fraction isn't too small.                   !
            !------------------------------------------------------------------------------!
            if ( (f_sun * slai) < tai_min) then
               f_sun   = tai_min / slai
               f_shade = 1.0 - f_sun
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     If the shade fraction is too small, make 100% sun.                       !
            !------------------------------------------------------------------------------!
            if ( (f_shade * slai) < tai_min) then
               f_sun   = 1.0
               f_shade = 0.0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         else
            !----- Night time, or tiny cohort, we won't solve them. -----------------------!
            f_sun     = 1.0
            f_shade   = 0.0
            kl        = 0.0
            ekl       = 0.0
            kay       = 0.5
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Split LAI and absorbed PAR accordingly.                                     !
         !---------------------------------------------------------------------------------!
         !----- Sun leaves get all direct radiation plus a fraction of diffuse. -----------!
         lai_ss  (1) = f_sun   * slai
         par_l_ss(1) = ( f_sun   * par_l_diffuse + par_l_beam ) * rshort
         !----- Shade leaves get the remaining diffuse radiation. -------------------------!
         lai_ss  (2) = f_shade * slai
         par_l_ss(2) = f_shade * par_l_diffuse * rshort
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the specific leaf areas of sun and shade leaves using CLM-4 approach.  !
         !---------------------------------------------------------------------------------!
         !----- Sun leaves. ---------------------------------------------------------------!
         if (lai_ss(1) >= tai_min) then
            sla_ss(1) = - ( ekl * ( sla_m(nveg) * kl + sla_m(nveg) + sla_0(nveg) * kay )   &
                          - sla_m(nveg) - sla_0(nveg) * kay ) / ( kay * kay * lai_ss(1) )
         else
            sla_ss(1) = sla_0(nveg)
         end if
         !----- Shade leaves. -------------------------------------------------------------!
         if (lai_ss(2) >= tai_min) then
            sla_ss(2) = ( slai * ( sla_0(nveg) + 0.5 * sla_m(nveg) * slai )                &
                        - sla_ss(1) * lai_ss(1) ) / lai_ss(2)
         else
            sla_ss(2) = sla_ss(1)
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Find the photosynthetic capacity and dark respiration.                     !
         !---------------------------------------------------------------------------------!
         where (sla_ss(:) /= 0.)
            vm0_ss(:) = vm0_0(nveg) / (0.001 * sla_ss(:))
         elsewhere
            vm0_ss(:) = 0.0
         end where
         rd0_ss(:) = dr_gamma(nveg) * vm0_ss(:)
         !---------------------------------------------------------------------------------!


         !----- Net long wave absorption by either soil or sfcwater. ----------------------!
         surface_netabs_longwave          = downward_lw_below - upward_lw_below
         !---------------------------------------------------------------------------------!



         !----- Soil+sfcwater+veg albedo (PAR, NIR, and Total Shortwave). -----------------!
         albedo_ground_par = upward_par_above_diffuse / (par_beam_norm + par_diff_norm)
         albedo_ground_nir = upward_nir_above_diffuse / (nir_beam_norm + nir_diff_norm)
         albedo_ground     = upward_rshort_above_diffuse
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Add urban contribution if running TEB.                                     !
         !---------------------------------------------------------------------------------!
         if (teb_spm == 1) then
            if (nint(g_urban) /= 0) then
               albedo_ground   = alb_town
               upward_lw_above = emis_town * stefan * ts_town**4
            end if
         end if
         albedt  = albedt  + patch_area * albedo_ground
         rlongup = rlongup + patch_area * upward_lw_above
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Set the total absorption by each pool.                                      !
         !---------------------------------------------------------------------------------!
         !----- Short wave radiation. -----------------------------------------------------!
         !..... Vegetation. ...............................................................!
         rshort_v = rshort * ( rshort_v_beam + rshort_v_diffuse)
         !..... Ground. ...................................................................!
         rshort_g = rshort * ( downward_par_below_beam    * abs_ground_par                 &
                             + downward_nir_below_beam    * abs_ground_nir                 &
                             + downward_par_below_diffuse * abs_ground_par                 &
                             + downward_nir_below_diffuse * abs_ground_nir )
         !..... Water pounding or snow. ...................................................!
         do k=1,ksn
            rshort_s(k) =  rshort * ( downward_par_below_beam    * abs_sfcw_par(k)         &
                                    + downward_nir_below_beam    * abs_sfcw_nir(k)         &
                                    + downward_par_below_diffuse * abs_sfcw_par(k)         &
                                    + downward_nir_below_diffuse * abs_sfcw_nir(k) )
         end do
         do k=ksn+1,mzs
            rshort_s(k) = 0.
         end do
         !---------------------------------------------------------------------------------!




         !----- Integrate the total absorbed light by ground (top soil plus TSW layers). --!
         rshort_gnd = rshort_g
         do k=1,ksn
            rshort_gnd = rshort_gnd + rshort_s(k)
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      In case rlong is not computed, zero out all longwave fluxes other than     !
         ! rlongup.  [On the first timestep, radiative fluxes may not be available until   !
         ! microphysics is called, and zeroing these fluxes avoids the imbalance of having !
         ! upward longwave without downward longwave in LEAF-3.  Also, this allows LEAF-3  !
         ! to run without radiation for all timesteps, if desired for testing.].           !
         !---------------------------------------------------------------------------------!
         if (rlong >= rad_rlong_min) then
            !----- Net absorption. --------------------------------------------------------!
            rlong_v   = lw_v
            if (ksn > 0) then
               rlong_g = surface_netabs_longwave
               rlong_s = 0.0
            else
               rlong_g = 0.0
               rlong_s = surface_netabs_longwave
            end if
            rlong_gnd = surface_netabs_longwave
            !------------------------------------------------------------------------------!
         else
            !----- Nothing happens. -------------------------------------------------------!
            rlong_v   = 0.0
            rlong_g   = 0.0
            rlong_s   = 0.0
            rlong_gnd = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_sfcrad
!==========================================================================================!
!==========================================================================================!








!==========================================================================================!
!==========================================================================================!
!     This sub-routine solves the long-wave radiation method using the two-stream model,   !
! adapted from ED-2.  Contrary to short-wave radiation, we don't normalise long-wave       !
! radiation.                                                                               !
!                                                                                          !
! Longo, M., 2013: Amazon forest response to changes in rainfall regime: results from an   !
!    individual-based dynamic vegetation model.   Ph.D. dissertation, Harvard University,  !
!    Cambridge, MA, Dec 2013.                                                              !
!                                                                                          !
! Additional useful references:                                                            !
!                                                                                          !
! Liou, K. N., 2002: An introduction to atmospheric radiation. Academic Press, Intl.       !
!    Geophys. Series, vol. 84.  Chapter 6., p. 257-347. (L02)                              !
!                                                                                          !
! Oleson, K. W., and co-authors, 2004: Technical description of version 4.0 of the         !
!    community land model (CLM). NCAR Technical note NCAR/TN-478+STR. 257pp. (CLM10)       !
!                                                                                          !
! Sellers, P. J., 1985: Canopy reflectance, photosynthesis and transpiration. Intl. J.     !
!    of remote sensing, 6, 1335-1372. (S85)                                                !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine leaf3_2stream_lw(grnd_emiss4,grnd_temp4,rlong_top4,nveg,lai,wai,leaf_temp       &
                           ,wood_temp,tir,dw_tirlo,uw_tirlo,uw_tirhi)
   use leaf_coms , only : clumping_factor         & ! intent(in)
                        , orient_factor           & ! intent(in)
                        , leaf_emiss_tir          & ! intent(in)
                        , wood_emiss_tir          & ! intent(in)
                        , phi1                    & ! intent(in)
                        , phi2                    & ! intent(in)
                        , mu_bar                  & ! intent(in)
                        , leaf_backscatter_tir    & ! intent(in)
                        , wood_backscatter_tir    ! ! intent(in)
   use rconstants, only : tiny_num8               & ! intent(in)
                        , stefan8                 ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4)                            , intent(in)   :: grnd_emiss4
   real(kind=4)                            , intent(in)   :: grnd_temp4
   real(kind=4)                            , intent(in)   :: rlong_top4
   integer                                 , intent(in)   :: nveg
   real(kind=8)                            , intent(in)   :: lai
   real(kind=8)                            , intent(in)   :: wai
   real(kind=8)                            , intent(in)   :: leaf_temp
   real(kind=8)                            , intent(in)   :: wood_temp
   real(kind=4)                            , intent(out)  :: tir
   real(kind=4)                            , intent(out)  :: dw_tirlo
   real(kind=4)                            , intent(out)  :: uw_tirlo
   real(kind=4)                            , intent(out)  :: uw_tirhi
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                :: ncoh
   integer                                                :: ncohp1
   integer                                                :: ncoh2p2
   real(kind=8)                                           :: ncohi
   integer                                                :: nsiz
   integer                                                :: i
   integer                                                :: im1
   integer                                                :: ip1
   integer                                                :: i2
   integer                                                :: i2p1
   integer                                                :: i2m1
   integer                                                :: i2p2
   logical                                                :: sing
   real(kind=8), dimension(:)              , allocatable  :: black
   real(kind=8), dimension(:)              , allocatable  :: expl_plus
   real(kind=8), dimension(:)              , allocatable  :: expl_minus
   real(kind=8), dimension(:)              , allocatable  :: gamm_plus
   real(kind=8), dimension(:)              , allocatable  :: gamm_minus
   real(kind=8), dimension(:)              , allocatable  :: down
   real(kind=8), dimension(:)              , allocatable  :: up
   real(kind=8), dimension(:,:)            , allocatable  :: mmat
   real(kind=8), dimension(:)              , allocatable  :: yvec
   real(kind=8), dimension(:)              , allocatable  :: xvec
   real(kind=8), dimension(:)              , allocatable  :: elai
   real(kind=8), dimension(:)              , allocatable  :: etai
   real(kind=8), dimension(:)              , allocatable  :: mu
   real(kind=8), dimension(:)              , allocatable  :: leaf_weight
   real(kind=8), dimension(:)              , allocatable  :: wood_weight
   real(kind=8), dimension(:)              , allocatable  :: beta
   real(kind=8), dimension(:)              , allocatable  :: epsil
   real(kind=8), dimension(:)              , allocatable  :: iota
   real(kind=8), dimension(:)              , allocatable  :: lambda
   real(kind=8)                                           :: temiss_four
   real(kind=8)                                           :: iota_g
   real(kind=8)                                           :: black_g
   real(kind=8)                                           :: down_sky
   !----- External functions. -------------------------------------------------------------!
   real(kind=8)                            , external     :: eifun8
   real(kind=4)                            , external     :: sngloff
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the size of the solution matrix and arrays and allocate arrays.              !
   !---------------------------------------------------------------------------------------!
   ncoh    = ceiling(clumping_factor(nveg) * lai + wai)
   ncohp1  = ncoh + 1
   ncoh2p2 = 2*ncoh + 2
   ncohi   = dble(1./ncoh)
   nsiz    = ncoh2p2
   allocate(black      (ncohp1)          )
   allocate(expl_plus  (ncohp1)          )
   allocate(expl_minus (ncohp1)          )
   allocate(gamm_plus  (ncohp1)          )
   allocate(gamm_minus (ncohp1)          )
   allocate(down       (ncohp1)          )
   allocate(up         (ncohp1)          )
   allocate(mmat       (ncoh2p2,ncoh2p2) )
   allocate(yvec       (ncoh2p2)         )
   allocate(xvec       (ncoh2p2)         )
   allocate(elai       (ncohp1)          )
   allocate(etai       (ncohp1)          )
   allocate(mu         (ncohp1)          )
   allocate(leaf_weight(ncohp1)          )
   allocate(wood_weight(ncohp1)          )
   allocate(beta       (ncohp1)          )
   allocate(epsil      (ncohp1)          )
   allocate(iota       (ncohp1)          )
   allocate(lambda     (ncohp1)          )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert some properties to double precision.                                      !
   !---------------------------------------------------------------------------------------!
   iota_g   = 1.d0 - dble(grnd_emiss4)
   black_g  = stefan8 * dble(grnd_temp4) ** 4
   down_sky = dble(rlong_top4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert some leaf variables to double precision, then find some general           !
   ! properties of leaves/branches (i.e., properties that do not depend on which band we   !
   ! are solving):  the extinction coefficient lambda, the layer transmittance             !
   ! coefficients tau_beam and tau_diff, and the weight of LAI and WAI relative to the     !
   ! total area.                                                                           !
   !---------------------------------------------------------------------------------------!
   do i= 1, ncoh
      !----- Area indices, correct LAI by clumping factor. --------------------------------!
      elai(i) = ncohi * clumping_factor(nveg) * lai
      etai(i) = elai(i) + ncohi * wai
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    The weighting factors are defined by the relative LAI and WAI areas.            !
      !------------------------------------------------------------------------------------!
      leaf_weight(i)  =  elai(i) / etai(i)
      wood_weight(i)  = 1.d0 -  leaf_weight(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     mu is the corrected mu_bar to account for finite crown area, if it is          !
      ! provided, otherwise, it is just the same as mu_bar.                                !
      !------------------------------------------------------------------------------------!
      mu(i) =  mu_bar(nveg)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Black is the emission of a black body at the same temperature of the layer.   !
      !  The temperature is the weighted average where the weights are the product between !
      !  the area and the emissivity.                                                      !
      !------------------------------------------------------------------------------------!
      temiss_four = ( leaf_weight(i) * leaf_emiss_tir(nveg) * leaf_temp ** 4               &
                    + wood_weight(i) * wood_emiss_tir(nveg) * wood_temp ** 4 )             &
                  / ( leaf_weight(i) * leaf_emiss_tir(nveg)                                &
                    + wood_weight(i) * wood_emiss_tir(nveg) )
      black(i) = stefan8 * temiss_four
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      iota is the single scattering albedo.  Because we assume absorptivity to be   !
      ! the same as emissivity, we use the latter to define it.                            !
      !------------------------------------------------------------------------------------!
      iota(i) = 1.d0 - ( leaf_weight(i) * leaf_emiss_tir(nveg)                             &
                       + wood_weight(i) * wood_emiss_tir(nveg) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Beta is the backward scattering of diffuse radiation. Epsil = 1 - 2*beta.      !
      !------------------------------------------------------------------------------------!
      beta (i) = leaf_weight(i) * leaf_backscatter_tir(nveg)                               &
               + wood_weight(i) * wood_backscatter_tir(nveg)
      epsil(i) = 1.d0 - 2.d0 * beta(i)
      !------------------------------------------------------------------------------------!



      !----- lambda is the coefficient associated with the optical depth. -----------------!
      lambda(i) = sqrt( ( 1.d0 - epsil(i) * iota(i) ) * ( 1.d0 - iota(i) ) ) / mu(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    gamm_plus and gamm_minus are the coefficients that relate upwelling and down-   !
      ! welling radiation.                                                                 !
      !------------------------------------------------------------------------------------!
      gamm_plus (i) = 5.d-1 * ( 1.d0 + sqrt( ( 1.d0 -            iota(i) )                 &
                                           / ( 1.d0 - epsil(i) * iota(i) ) ) )
      gamm_minus(i) = 5.d-1 * ( 1.d0 - sqrt( ( 1.d0 -            iota(i) )                 &
                                           / ( 1.d0 - epsil(i) * iota(i) ) ) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    expl_plus and expl_minus are the transmitivity of diffuse light.                !
      !------------------------------------------------------------------------------------!
      expl_plus (i) = exp(   lambda(i) * etai(i) )
      expl_minus(i) = exp( - lambda(i) * etai(i) )
      !------------------------------------------------------------------------------------!
   end do


   !---------------------------------------------------------------------------------------!
   !     Define some boundary conditions for the vector properties above.                  !
   !---------------------------------------------------------------------------------------!
   elai       (ncohp1) = 0.d0
   etai       (ncohp1) = 0.d0
   leaf_weight(ncohp1) = 5.d-1
   wood_weight(ncohp1) = 5.d-1
   mu         (ncohp1) = 1.d0
   black      (ncohp1) = 0.d0
   iota       (ncohp1) = 1.d0
   beta       (ncohp1) = 0.d0
   epsil      (ncohp1) = 1.d0 - 2.d0 * beta(ncohp1)
   lambda     (ncohp1) = 0.d0
   gamm_plus  (ncohp1) = 1.d0
   gamm_minus (ncohp1) = 0.d0
   expl_plus  (ncohp1) = 1.d0
   expl_minus (ncohp1) = 1.d0
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Fill in the right hand side vector (Y) and the matrix (X)                         !
   !---------------------------------------------------------------------------------------!
   !------ Initialise the vector and the matrix.  The matrix is a sparse one... -----------!
   mmat(:,:) = 0.d0
   yvec(:)   = 0.d0
   !------ Add the bottom and top boundary conditions. ------------------------------------!
   mmat(1,1)         = ( gamm_minus(1) - iota_g * gamm_plus (1) ) * expl_minus(1)
   mmat(1,2)         = ( gamm_plus (1) - iota_g * gamm_minus(1) ) * expl_plus (1)
   mmat(nsiz,nsiz-1) = gamm_plus (ncohp1)
   mmat(nsiz,nsiz  ) = gamm_minus(ncohp1)
   yvec(1)           = (1.d0 - iota_g) * black_g - (1.d0 - iota_g) * black(1)
   yvec(nsiz)        = down_sky - black(ncoh+1)
   !---------------------------------------------------------------------------------------!
   do i=1,ncoh
      !----- Find auxiliary indices. ------------------------------------------------------!
      ip1  =  i + 1
      i2   =  2 * i
      i2m1 = i2 - 1
      i2p1 = i2 + 1
      i2p2 = i2 + 2
      !------------------------------------------------------------------------------------!

      !----- Make elements of C vector. ---------------------------------------------------!
      yvec(i2  ) = black(ip1) - black(i)
      yvec(i2p1) = black(ip1) - black(i)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Make elements of the tri-diagonal A array.                                     !
      !------------------------------------------------------------------------------------!
      mmat(i2  ,i2m1) =   gamm_plus   (i)
      mmat(i2  ,i2  ) =   gamm_minus  (i)
      mmat(i2  ,i2p1) = - gamm_plus (ip1) * expl_minus(ip1)
      mmat(i2  ,i2p2) = - gamm_minus(ip1) * expl_plus (ip1)
      mmat(i2p1,i2m1) =   gamm_minus  (i)
      mmat(i2p1,i2  ) =   gamm_plus   (i)
      mmat(i2p1,i2p1) = - gamm_minus(ip1) * expl_minus(ip1)
      mmat(i2p1,i2p2) = - gamm_plus (ip1) * expl_plus (ip1)
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !       Solve the linear system.  In the future we could use a band diagonal solver,    !
   ! which is a lot cheaper than the regular Gauss elimination, but for the time being, we !
   ! go with a tested method.                                                              !
   !---------------------------------------------------------------------------------------!
   call lisys_solver8(nsiz,mmat,yvec,xvec,sing)
   if (sing) then
      call abort_run('LW radiation failed... The matrix is singular!'                      &
                    ,'leaf3_2stream_lw','leaf3_radiate.f90')
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Copy the solution to the vectors, using the properties:                           !
   !  Un = Un(P)                                                                           !
   !  Dn = Dn(P)                                                                           !
   !---------------------------------------------------------------------------------------!
   do i=1,ncohp1

      !----- Auxiliary indices. -----------------------------------------------------------!
      i2   =  2 * i
      i2m1 = i2 - 1
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Retrieve the downward and upward diffuse (hemispheric) radiation, by using    !
      ! the solutions for full TAI                                                         !
      !------------------------------------------------------------------------------------!
      down(i) = xvec(i2m1) * gamm_plus (i) * expl_minus(i)                                 &
              + xvec  (i2) * gamm_minus(i) * expl_plus (i)                                 &
              + black  (i)
      up  (i) = xvec(i2m1) * gamm_minus(i) * expl_minus(i)                                 &
              + xvec  (i2) * gamm_plus (i) * expl_plus (i)                                 &
              + black  (i)
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!





   !------ Save the fluxes reaching the surface and leaving the top. ----------------------!
   dw_tirlo = sngloff(down(1)     , tiny_num8)
   uw_tirlo = sngloff(up  (1)     , tiny_num8)
   uw_tirhi = sngloff(up  (ncohp1), tiny_num8)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Save the radiation fluxes to the output variable.                                 !
   !---------------------------------------------------------------------------------------!
   tir = sngloff(down(ncohp1) - down(1) + up(1) - up(ncohp1), tiny_num8)
   !---------------------------------------------------------------------------------------!



   !----- Free memory. --------------------------------------------------------------------!
   deallocate(black      )
   deallocate(expl_plus  )
   deallocate(expl_minus )
   deallocate(gamm_plus  )
   deallocate(gamm_minus )
   deallocate(down       )
   deallocate(up         )
   deallocate(mmat       )
   deallocate(yvec       )
   deallocate(xvec       )
   deallocate(elai       )
   deallocate(etai       )
   deallocate(mu         )
   deallocate(leaf_weight)
   deallocate(wood_weight)
   deallocate(beta       )
   deallocate(epsil      )
   deallocate(iota       )
   deallocate(lambda     )
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_2stream_lw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine solves the short-wave radiation method using the two-stream model.  !
! Radiation fluxes are standardised here, and scaled in the main sub-routine.              !
!                                                                                          !
! Longo, M., 2013: Amazon forest response to changes in rainfall regime: results from an   !
!    individual-based dynamic vegetation model.   Ph.D. dissertation, Harvard University,  !
!    Cambridge, MA, Dec 2013.                                                              !
!                                                                                          !
! Additional useful references:                                                            !
!                                                                                          !
! Liou, K. N., 2002: An introduction to atmospheric radiation. Academic Press, Intl.       !
!    Geophys. Series, vol. 84.  Chapter 6., p. 257-347. (L02)                              !
!                                                                                          !
! Oleson, K. W., and co-authors, 2004: Technical description of the community land model   !
!   (CLM). NCAR Technical note NCAR/TN-461+STR. 186pp. (CLM04)                             !
!                                                                                          !
! Sellers, P. J., 1985: Canopy reflectance, photosynthesis and transpiration. Intl. J.     !
!    of remote sensing, 6, 1335-1372. (S85)                                                !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine leaf3_2stream_sw(grnd_alb_par4,grnd_alb_nir4,cosz4,par_beam_norm4               &
                           ,par_diff_norm4,nir_beam_norm4,nir_diff_norm4,nveg,lai,wai      &
                           ,par_beam,par_diff,sw_abs_beam,sw_abs_diff,dw_parlo_beam        &
                           ,dw_parlo_diff,uw_parhi_diff,dw_nirlo_beam,dw_nirlo_diff        &
                           ,uw_nirhi_diff)
   use rconstants , only : tiny_num8               ! ! intent(in)
   use leaf_coms  , only : clumping_factor         & ! intent(in)
                         , orient_factor           & ! intent(in)
                         , phi1                    & ! intent(in)
                         , phi2                    & ! intent(in)
                         , mu_bar                  & ! intent(in)
                         , leaf_scatter_vis        & ! intent(in)
                         , leaf_scatter_nir        & ! intent(in)
                         , leaf_backscatter_nir    & ! intent(in)
                         , leaf_backscatter_vis    & ! intent(in)
                         , wood_scatter_nir        & ! intent(in)
                         , wood_scatter_vis        & ! intent(in)
                         , wood_backscatter_nir    & ! intent(in)
                         , wood_backscatter_vis    ! ! intent(in)
   use mem_radiate, only : rad_cosz_min8           ! ! intent(in)
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   integer                                                :: ncoh
   integer                                                :: ncohp1
   integer                                                :: ncoh2p2
   real(kind=8)                                           :: ncohi
   real(kind=4)                            , intent(in)   :: grnd_alb_par4
   real(kind=4)                            , intent(in)   :: grnd_alb_nir4
   real(kind=4)                            , intent(in)   :: cosz4
   real(kind=4)                            , intent(in)   :: par_beam_norm4
   real(kind=4)                            , intent(in)   :: par_diff_norm4
   real(kind=4)                            , intent(in)   :: nir_beam_norm4
   real(kind=4)                            , intent(in)   :: nir_diff_norm4
   integer                                 , intent(in)   :: nveg
   real(kind=8)                            , intent(in)   :: lai
   real(kind=8)                            , intent(in)   :: wai
   real(kind=4)                            , intent(out)  :: par_beam
   real(kind=4)                            , intent(out)  :: par_diff
   real(kind=4)                            , intent(out)  :: sw_abs_beam
   real(kind=4)                            , intent(out)  :: sw_abs_diff
   real(kind=4)                            , intent(out)  :: uw_parhi_diff
   real(kind=4)                            , intent(out)  :: uw_nirhi_diff
   real(kind=4)                            , intent(out)  :: dw_parlo_beam
   real(kind=4)                            , intent(out)  :: dw_parlo_diff
   real(kind=4)                            , intent(out)  :: dw_nirlo_beam
   real(kind=4)                            , intent(out)  :: dw_nirlo_diff
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                :: iband
   integer                                                :: nsiz
   integer                                                :: i
   integer                                                :: ip1
   integer                                                :: i2
   integer                                                :: i2p1
   integer                                                :: i2m1
   integer                                                :: i2p2
   logical                                                :: sing
   real(kind=4)                                           :: nir_beam
   real(kind=4)                                           :: nir_diff
   real(kind=8)                                           :: leaf_scatter
   real(kind=8)                                           :: wood_scatter
   real(kind=8)                                           :: leaf_backscatter
   real(kind=8)                                           :: wood_backscatter
   real(kind=8), dimension(:)              , allocatable  :: expl_plus
   real(kind=8), dimension(:)              , allocatable  :: expl_minus
   real(kind=8), dimension(:)              , allocatable  :: expm0_minus
   real(kind=8), dimension(:)              , allocatable  :: gamm_plus
   real(kind=8), dimension(:)              , allocatable  :: gamm_minus
   real(kind=8), dimension(:)              , allocatable  :: delta
   real(kind=8), dimension(:)              , allocatable  :: upsilon
   real(kind=8), dimension(:)              , allocatable  :: down0
   real(kind=8), dimension(:)              , allocatable  :: down
   real(kind=8), dimension(:)              , allocatable  :: up
   real(kind=8), dimension(:,:)            , allocatable  :: mmat
   real(kind=8), dimension(:)              , allocatable  :: yvec
   real(kind=8), dimension(:)              , allocatable  :: xvec
   real(kind=8), dimension(:)              , allocatable  :: elai
   real(kind=8), dimension(:)              , allocatable  :: etai
   real(kind=8), dimension(:)              , allocatable  :: mu
   real(kind=8), dimension(:)              , allocatable  :: mu0
   real(kind=8), dimension(:)              , allocatable  :: leaf_weight
   real(kind=8), dimension(:)              , allocatable  :: wood_weight
   real(kind=8), dimension(:)              , allocatable  :: beta0
   real(kind=8), dimension(:)              , allocatable  :: beta
   real(kind=8), dimension(:)              , allocatable  :: epsil
   real(kind=8), dimension(:)              , allocatable  :: epsil0
   real(kind=8), dimension(:)              , allocatable  :: iota
   real(kind=8), dimension(:)              , allocatable  :: lambda
   real(kind=8), dimension(:)              , allocatable  :: iota_ratio
   real(kind=8), dimension(:)              , allocatable  :: proj_area
   real(kind=8), dimension(:)              , allocatable  :: a_aux
   real(kind=8), dimension(:)              , allocatable  :: s_aux
   real(kind=8)                                           :: iota_g
   real(kind=8)                                           :: iota_g_par
   real(kind=8)                                           :: iota_g_nir
   real(kind=8)                                           :: down_sky
   real(kind=8)                                           :: down0_sky
   real(kind=8)                                           :: czen
   real(kind=8)                                           :: par_beam_norm
   real(kind=8)                                           :: par_diff_norm
   real(kind=8)                                           :: nir_beam_norm
   real(kind=8)                                           :: nir_diff_norm
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)                            , external     :: sngloff
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the size of the solution matrix and array and allocate arrays.               !
   !---------------------------------------------------------------------------------------!
   ncoh    = ceiling(clumping_factor(nveg) * lai + wai)
   ncohp1  = ncoh + 1
   ncoh2p2 = 2*ncoh + 2
   ncohi   = dble(1./ncoh)
   nsiz    = ncoh2p2
   allocate(expl_plus  (ncohp1)         )
   allocate(expl_minus (ncohp1)         )
   allocate(expm0_minus(ncohp1)         )
   allocate(gamm_plus  (ncohp1)         )
   allocate(gamm_minus (ncohp1)         )
   allocate(delta      (ncohp1)         )
   allocate(upsilon    (ncohp1)         )
   allocate(down0      (ncohp1)         )
   allocate(down       (ncohp1)         )
   allocate(up         (ncohp1)         )
   allocate(mmat       (ncoh2p2,ncoh2p2))
   allocate(yvec       (ncoh2p2)        )
   allocate(xvec       (ncoh2p2)        )
   allocate(elai       (ncohp1)         )
   allocate(etai       (ncohp1)         )
   allocate(mu         (ncohp1)         )
   allocate(mu0        (ncohp1)         )
   allocate(leaf_weight(ncohp1)         )
   allocate(wood_weight(ncohp1)         )
   allocate(beta0      (ncohp1)         )
   allocate(beta       (ncohp1)         )
   allocate(epsil      (ncohp1)         )
   allocate(epsil0     (ncohp1)         )
   allocate(iota       (ncohp1)         )
   allocate(lambda     (ncohp1)         )
   allocate(iota_ratio (ncohp1)         )
   allocate(proj_area  (ncohp1)         )
   allocate(a_aux      (ncohp1)         )
   allocate(s_aux      (ncohp1)         )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert cosine of angle of incidence to double precision.                         !
   !---------------------------------------------------------------------------------------!
   czen          = max(rad_cosz_min8,dble(cosz4))
   iota_g_par    = dble(grnd_alb_par4)
   iota_g_nir    = dble(grnd_alb_nir4)
   par_beam_norm = dble(par_beam_norm4)
   par_diff_norm = dble(par_diff_norm4)
   nir_beam_norm = dble(nir_beam_norm4)
   nir_diff_norm = dble(nir_diff_norm4)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Here we loop over the two bands (band 1 is PAR and band 2 is NIR), and solve the  !
   ! same set of equations.                                                                !
   !---------------------------------------------------------------------------------------!
   bandloop: do iband = 1,2
      select case (iband)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Visible (PAR).                                                              !
         !---------------------------------------------------------------------------------!
         leaf_scatter     = leaf_scatter_vis    (nveg)
         wood_scatter     = wood_scatter_vis    (nveg)
         leaf_backscatter = leaf_backscatter_vis(nveg)
         wood_backscatter = wood_backscatter_vis(nveg)
         down_sky         = par_diff_norm
         down0_sky        = par_beam_norm
         iota_g           = iota_g_par
         !---------------------------------------------------------------------------------!
      case (2)
         !---------------------------------------------------------------------------------!
         !     Near-Infrared (NIR).                                                        !
         !---------------------------------------------------------------------------------!
         leaf_scatter     = leaf_scatter_nir    (nveg)
         wood_scatter     = wood_scatter_nir    (nveg)
         leaf_backscatter = leaf_backscatter_nir(nveg)
         wood_backscatter = wood_backscatter_nir(nveg)
         down_sky         = nir_diff_norm
         down0_sky        = nir_beam_norm
         iota_g           = iota_g_nir
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find some general properties of leaves/branches and also direct radiation      !
      ! properties.                                                                        !
      !------------------------------------------------------------------------------------!
      directloop: do i=1,ncoh
         !----- Area indices, correct LAI by clumping factor. -----------------------------!
         elai(i) = ncohi * clumping_factor(nveg) * lai
         etai(i) = elai(i) + ncohi * wai
         !---------------------------------------------------------------------------------!



         !----- The weighting factors are defined by the relative LAI and WAI areas. ------!
         leaf_weight(i) = elai(i) / etai(i)
         wood_weight(i) = 1.d0 - leaf_weight(i)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     We find the inverse optical depth of the direct radiation (mu0), following  !
         ! CLM10 (equation 3.3 and text after equation 3.3).                               !
         !---------------------------------------------------------------------------------!
         proj_area(i) = phi1(nveg) + phi2(nveg) * czen
         mu0      (i) = proj_area(i) / czen
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     The inverse optical depth of the diffuse radiation (mu), same as in long-   !
         ! wave radiation.                                                                 !
         !---------------------------------------------------------------------------------!
         mu (i) = mu_bar(nveg)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Find the backscatter for direct beam radiation, following CLM10, equations !
         ! (3.14) and (3.15).  The forward scattering was omitted on both equations        !
         ! because they are different for PAR and NIR, but in both cases they would cancel !
         ! out.  epsil0 = 1 - 2*beta0, no special definition...                            !
         !---------------------------------------------------------------------------------!
         iota_ratio(i) = 1.d0 / ( 2.d0 * ( 1.d0 + phi2(nveg) * mu0(i) ) )                  &
                       * ( 1.d0 - phi1(nveg) * mu0(i) / ( 1.d0 + phi2(nveg) * mu0(i) )     &
                                * log ( ( 1.d0  + ( phi1(nveg) + phi2(nveg) ) * mu0(i) )   &
                                      / ( phi1(nveg) * mu0(i) ) ) )
         beta0     (i) = iota_ratio(i) * ( 1.d0 + mu0(i) / mu(i) )
         epsil0    (i) = 1.d0 - 2.d0 * beta0(i)
         !---------------------------------------------------------------------------------!


         !----- Expm0_minus is the transmissivity of direct radiation. --------------------!
         expm0_minus(i) = exp( - etai  (i) / mu0 (i) )
         !---------------------------------------------------------------------------------!
      end do directloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Define some boundary conditions for the vector properties above.               !
      !------------------------------------------------------------------------------------!
      elai       (ncohp1) = 0.d0
      etai       (ncohp1) = 0.d0
      leaf_weight(ncohp1) = 5.d-1
      wood_weight(ncohp1) = 5.d-1
      proj_area  (ncohp1) = 5.d-1
      mu0        (ncohp1) = czen / proj_area(ncohp1)
      mu         (ncohp1) = 1.d0
      iota_ratio (ncohp1) = 5.d-1 * ( 1.d0 - 5.d-1 * mu0(ncohp1)                           &
                                      * log ( 1.d0  / ( 5.d-1 * mu0(ncohp1) ) + 1.d0 ) )
      beta0      (ncohp1) = iota_ratio(ncohp1) * ( mu0(ncohp1) + mu(ncohp1) ) / mu(ncohp1)
      epsil0     (ncohp1) = 1.d0 - 2.d0 * beta0(ncohp1)
      expm0_minus(ncohp1) = 1.d0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the direct radiation profile using the exponential attenuation curve.     !
      !------------------------------------------------------------------------------------!
      down0(ncohp1) = down0_sky
      do i = ncoh,1,-1
         ip1      = i + 1
         down0(i) = down0(ip1) * expm0_minus(i) 
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the diffuse radiation properties.                                         !
      !------------------------------------------------------------------------------------!
      diffuseloop: do i=1,ncoh
         !----- Scattering coefficient. ---------------------------------------------------!
         iota      (i) = leaf_weight(i) * leaf_scatter + wood_weight(i) * wood_scatter
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Beta is the backward scattering of diffuse radiation. Epsil = 1 - 2*beta.   !
         !---------------------------------------------------------------------------------!
         beta  (i) = leaf_weight(i) * leaf_backscatter + wood_weight(i) * wood_backscatter
         epsil (i) = 1.d0 - 2.d0 * beta(i)
         !---------------------------------------------------------------------------------!



         !----- lambda is the coefficient associated with the optical depth. --------------!
         lambda(i) = sqrt( ( 1.d0 - epsil(i) * iota(i) ) * ( 1.d0 - iota(i) ) ) / mu(i)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find some auxiliary variables to determine the right hand side.             !
         !---------------------------------------------------------------------------------!
         a_aux  (i) = - ( ( 1.d0 - epsil(i) * iota(i) ) * iota(i) / mu(i)                  &
                        + epsil0(i) * iota(i) / mu0(i) ) * down0(2) / mu0(i)
         s_aux  (i) = - ( ( 1.d0 - iota (i) ) * epsil0(i) * iota(i) / mu(i)                &
                        + iota(i) / mu0(i) ) * down0(2) / mu0(i)
         delta  (i) = ( a_aux(i) + s_aux(i) ) * mu0(i) * mu0(i)                            &
                    / ( 2.d0 * ( 1.d0 - lambda(i) * lambda(i) * mu0(i) * mu0(i) ) )
         upsilon(i) = ( a_aux(i) - s_aux(i) ) * mu0(i) * mu0(i)                            &
                    / ( 2.d0 * ( 1.d0 - lambda(i) * lambda(i) * mu0(i) * mu0(i) ) )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    gamm_plus and gamm_minus are the coefficients that relate upwelling and      !
         ! downwelling radiation.                                                          !
         !---------------------------------------------------------------------------------!
         gamm_plus (i) = 5.d-1 * ( 1.d0 + sqrt( ( 1.d0 -            iota(i) )              &
                                              / ( 1.d0 - epsil(i) * iota(i) ) ) )
         gamm_minus(i) = 5.d-1 * ( 1.d0 - sqrt( ( 1.d0 -            iota(i) )              &
                                              / ( 1.d0 - epsil(i) * iota(i) ) ) )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    expl_plus and expl_minus are the transmitivity of diffuse light.             !
         !---------------------------------------------------------------------------------!
         expl_plus  (i) = exp(   lambda(i) * etai(i) )
         expl_minus (i) = exp( - lambda(i) * etai(i) )
         !---------------------------------------------------------------------------------!
      end do diffuseloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Define some boundary conditions for the vector properties above.               !
      !------------------------------------------------------------------------------------!
      iota       (ncohp1) = 1.d0
      beta       (ncohp1) = 0.d0
      epsil      (ncohp1) = 1.d0 - 2.d0 * beta(ncohp1)
      lambda     (ncohp1) = 0.d0
      a_aux      (ncohp1) = 0.d0
      s_aux      (ncohp1) = 0.d0
      delta      (ncohp1) = 0.d0
      upsilon    (ncohp1) = 0.d0
      a_aux      (ncohp1) = -  epsil0(ncohp1) * down0_sky / ( mu0(ncohp1) * mu0(ncohp1) )
      s_aux      (ncohp1) = -  iota  (ncohp1) * down0_sky / ( mu0(ncohp1) * mu0(ncohp1) )
      delta      (ncohp1) = 5.d-1 * ( a_aux(ncohp1) + s_aux(ncohp1) )                      &
                                  * mu0(ncohp1) * mu0(ncohp1)
      upsilon    (ncohp1) = 5.d-1 * ( a_aux(ncohp1) - s_aux(ncohp1) )                      &
                                  * mu0(ncohp1) * mu0(ncohp1)
      gamm_plus  (ncohp1) = 1.d0
      gamm_minus (ncohp1) = 0.d0
      expl_plus  (ncohp1) = 1.d0
      expl_minus (ncohp1) = 1.d0
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Fill in the right hand side vector (Y) and the matrix (X)                      !
      !------------------------------------------------------------------------------------!
      !------ Initialise the vector and the matrix.  The matrix is a sparse one... --------!
      mmat(:,:)         = 0.d0
      yvec(:)           = 0.d0
      !------ Add the bottom and top boundary conditions. ---------------------------------!
      mmat(1,1)         = (gamm_minus(1) - iota_g * gamm_plus (1)) * expl_minus(1)
      mmat(1,2)         = (gamm_plus (1) - iota_g * gamm_minus(1)) * expl_plus (1)
      mmat(nsiz,nsiz-1) = gamm_plus (ncoh+1)
      mmat(nsiz,nsiz  ) = gamm_minus(ncoh+1)
      yvec(1)           = iota_g * down0(1)                                                &
                        - ( upsilon(1) - iota_g * delta(1) ) * expm0_minus(1)
      yvec(nsiz)        = down_sky - delta(ncoh+1)
      do i=1,ncoh
         !----- Find auxiliary indices. ---------------------------------------------------!
         ip1  =  i + 1
         i2   =  2 * i
         i2m1 = i2 - 1
         i2p1 = i2 + 1
         i2p2 = i2 + 2
         !---------------------------------------------------------------------------------!

         !----- Make elements of C vector. ------------------------------------------------!
         yvec(i2  ) = delta  (ip1) * expm0_minus(ip1) - delta  (i)
         yvec(i2p1) = upsilon(ip1) * expm0_minus(ip1) - upsilon(i)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Make elements of the tri-diagonal A array.                                  !
         !---------------------------------------------------------------------------------!
         mmat(i2  ,i2m1) =   gamm_plus   (i)
         mmat(i2  ,i2  ) =   gamm_minus  (i)
         mmat(i2  ,i2p1) = - gamm_plus (ip1) * expl_minus(ip1)
         mmat(i2  ,i2p2) = - gamm_minus(ip1) * expl_plus (ip1)
         mmat(i2p1,i2m1) =   gamm_minus  (i)
         mmat(i2p1,i2  ) =   gamm_plus   (i)
         mmat(i2p1,i2p1) = - gamm_minus(ip1) * expl_minus(ip1)
         mmat(i2p1,i2p2) = - gamm_plus (ip1) * expl_plus (ip1)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Solve the linear system.  In the future we could use a tridiagonal solver,   !
      ! which is a lot cheaper than the regular Gauss elimination, but for the time being, !
      ! we go with a tested method.                                                        !
      !------------------------------------------------------------------------------------!
      call lisys_solver8(nsiz,mmat,yvec,xvec,sing)
      if (sing) then
         call abort_run('SW radiation failed... The matrix is singular!'                   &
                       ,'leaf3_2stream_sw','leaf3_radiate.f90')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Copy the solution to the vectors, using the properties:                        !
      !  Un = Un(P)                                                                        !
      !  Dn = Dn(P)                                                                        !
      !------------------------------------------------------------------------------------!
      do i=1,ncoh+1

         !----- Auxiliary indices. --------------------------------------------------------!
         i2   =  2 * i
         i2m1 = i2 - 1
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Retrieve the downward and upward diffuse (hemispheric) radiation, by using !
         ! the solutions for full TAI                                                      !
         !---------------------------------------------------------------------------------!
         down(i) = xvec(i2m1) * gamm_plus (i) * expl_minus (i)                             &
                 + xvec  (i2) * gamm_minus(i) * expl_plus  (i)                             &
                              + delta     (i) * expm0_minus(i)
         up  (i) = xvec(i2m1) * gamm_minus(i) * expl_minus (i)                             &
                 + xvec  (i2) * gamm_plus (i) * expl_plus  (i)                             &
                              + upsilon   (i) * expm0_minus(i)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Save the fluxes that we will use outside this sub-routine.                     !
      !------------------------------------------------------------------------------------!
      select case (iband)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Visible (PAR).                                                              !
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Save the radiation fluxes to the output variable.                           !
         !---------------------------------------------------------------------------------!
         par_beam = sngloff( down0(ncohp1) - down0(1)                     , tiny_num8)
         par_diff = sngloff( down (ncohp1) - down (1) + up(1) - up(ncohp1), tiny_num8)
         !---------------------------------------------------------------------------------!



         !------ Save the fluxes reaching the surface and leaving the top. ----------------!
         dw_parlo_beam = sngloff(down0(1     ), tiny_num8)
         dw_parlo_diff = sngloff(down (1     ), tiny_num8)
         uw_parhi_diff = sngloff(up   (ncohp1), tiny_num8)
         !---------------------------------------------------------------------------------!

      case (2)
         !---------------------------------------------------------------------------------!
         !     Near-Infrared (NIR).                                                        !
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Save the radiation fluxes to the output variable.                           !
         !---------------------------------------------------------------------------------!
         nir_beam = sngloff( down0(ncohp1) - down0(1)                     , tiny_num8)
         nir_diff = sngloff( down (ncohp1) - down (1) + up(1) - up(ncohp1), tiny_num8)
         !---------------------------------------------------------------------------------!



         !------ Save the fluxes reaching the surface and leaving the top. ----------------!
         dw_nirlo_beam = sngloff(down0(1     ), tiny_num8)
         dw_nirlo_diff = sngloff(down (1     ), tiny_num8)
         uw_nirhi_diff = sngloff(up   (ncohp1), tiny_num8)
         !---------------------------------------------------------------------------------!

      end select
   end do bandloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Total normalised radiation.                                                       !
   !---------------------------------------------------------------------------------------!
   sw_abs_beam = par_beam + nir_beam
   sw_abs_diff = par_diff + nir_diff
   !---------------------------------------------------------------------------------------!



   !----- Free memory. --------------------------------------------------------------------!
   deallocate(expl_plus  )
   deallocate(expl_minus )
   deallocate(expm0_minus)
   deallocate(gamm_plus  )
   deallocate(gamm_minus )
   deallocate(delta      )
   deallocate(upsilon    )
   deallocate(down0      )
   deallocate(down       )
   deallocate(up         )
   deallocate(mmat       )
   deallocate(yvec       )
   deallocate(xvec       )
   deallocate(elai       )
   deallocate(etai       )
   deallocate(mu         )
   deallocate(mu0        )
   deallocate(leaf_weight)
   deallocate(wood_weight)
   deallocate(beta0      )
   deallocate(beta       )
   deallocate(epsil      )
   deallocate(epsil0     )
   deallocate(iota       )
   deallocate(lambda     )
   deallocate(iota_ratio )
   deallocate(proj_area  )
   deallocate(a_aux      )
   deallocate(s_aux      )
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_2stream_sw
!==========================================================================================!
!==========================================================================================!
