!==========================================================================================!
!==========================================================================================!
! MODULE RADIATE_UTILS
!------------------------------------------------------------------------------------------!
module radiate_utils
   contains
   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine allows the user to change the solar radiation splitting that    !
   ! comes from the meteorological forcing.  This is useful for when we don't really know  !how the   !
   ! radiation terms were found, or if we want to play around with the numbers.            !
   !---------------------------------------------------------------------------------------!
   subroutine solar_radiation_breakdown(cgrid,ipy)
      use ed_state_vars        , only : edtype        ! ! intent(in)
      use met_driver_coms      , only : imetrad       ! ! intent(in)
      use canopy_radiation_coms, only : fvis_beam_def & ! intent(in)
                                      , fvis_diff_def & ! intent(in)
                                      , fnir_beam_def & ! intent(in)
                                      , fnir_diff_def & ! intent(in)
                                      , cosz_min      ! ! intent(in)


      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype), target     :: cgrid
      integer     , intent(in) :: ipy
      !----- Local variables. -------------------------------------------------------------!
      real                     :: rshort_beam
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Add the radiation components, regardless of the method.                        !
      !------------------------------------------------------------------------------------!
      cgrid%met(ipy)%rshort_diffuse = cgrid%met(ipy)%par_diffuse                           &
                                    + cgrid%met(ipy)%nir_diffuse
      cgrid%met(ipy)%rshort         = cgrid%met(ipy)%rshort_diffuse                        &
                                    + cgrid%met(ipy)%par_beam + cgrid%met(ipy)%nir_beam
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Choose what to do with the radiation depending on the user choice.             !
      !------------------------------------------------------------------------------------!
      select case (imetrad)
      case (0)
         !---------------------------------------------------------------------------------!
         !    We just make sure there is no beam radiation in case the Sun is too low,     !
         ! otherwise, we leave radiation as is.                                            !
         !---------------------------------------------------------------------------------!
         if (cgrid%cosz(ipy) <= cosz_min) then
            cgrid%met(ipy)%rshort_diffuse = cgrid%met(ipy)%rshort
            cgrid%met(ipy)%nir_diffuse    = cgrid%met(ipy)%nir_diffuse                     &
                                          + cgrid%met(ipy)%nir_beam
            cgrid%met(ipy)%par_diffuse    = cgrid%met(ipy)%par_diffuse                     &
                                          + cgrid%met(ipy)%par_beam
            cgrid%met(ipy)%nir_beam       = 0.
            cgrid%met(ipy)%par_beam       = 0.
         end if
         !---------------------------------------------------------------------------------!

      case (1)
         !---------------------------------------------------------------------------------!
         !     Re-calculate the diffuse shortwave radiation using the SiB method.          !
         !---------------------------------------------------------------------------------!
         call short2diff_sib(cgrid%met(ipy)%rshort,cgrid%cosz(ipy)                         &
                            ,cgrid%met(ipy)%rshort_diffuse)
         rshort_beam = cgrid%met(ipy)%rshort - cgrid%met(ipy)%rshort_diffuse
         cgrid%met(ipy)%nir_diffuse = fnir_diff_def * cgrid%met(ipy)%rshort_diffuse
         cgrid%met(ipy)%par_diffuse = fvis_diff_def * cgrid%met(ipy)%rshort_diffuse
         cgrid%met(ipy)%nir_beam    = fnir_beam_def * rshort_beam
         cgrid%met(ipy)%par_beam    = fvis_beam_def * rshort_beam
         !---------------------------------------------------------------------------------!

      case (2)
         !---------------------------------------------------------------------------------!
         !     Re-calculate all components using the Weiss and Norman (1985) method.       !
         !---------------------------------------------------------------------------------!
         call short_bdown_weissnorman(cgrid%met(ipy)%rshort                                &
                                     ,cgrid%met(ipy)%prss                                  &
                                     ,cgrid%cosz(ipy)                                      &
                                     ,cgrid%met(ipy)%par_beam                              &
                                     ,cgrid%met(ipy)%par_diffuse                           &
                                     ,cgrid%met(ipy)%nir_beam                              &
                                     ,cgrid%met(ipy)%nir_diffuse                           &
                                     ,cgrid%met(ipy)%rshort_diffuse )
         !----- Make sure the total radiation is preserved. -------------------------------!
         cgrid%met(ipy)%rshort = cgrid%met(ipy)%par_beam + cgrid%met(ipy)%par_diffuse      &
                               + cgrid%met(ipy)%nir_beam + cgrid%met(ipy)%nir_diffuse
         !---------------------------------------------------------------------------------!

      case (3)
         !---------------------------------------------------------------------------------!
         !      Make all radiation diffuse.                                                !
         !---------------------------------------------------------------------------------!
         cgrid%met(ipy)%rshort_diffuse = cgrid%met(ipy)%rshort
         cgrid%met(ipy)%nir_diffuse    = cgrid%met(ipy)%nir_diffuse                        &
                                       + cgrid%met(ipy)%nir_beam
         cgrid%met(ipy)%par_diffuse    = cgrid%met(ipy)%par_diffuse                        &
                                       + cgrid%met(ipy)%par_beam
         cgrid%met(ipy)%nir_beam       = 0.
         cgrid%met(ipy)%par_beam       = 0.
         !---------------------------------------------------------------------------------!

      case (4)
         !---------------------------------------------------------------------------------!
         !     Make all radiation direct, except when the sun is too low, in which case we !
         ! must set to diffuse to avoid singularities.                                     !
         !---------------------------------------------------------------------------------!
         if (cgrid%cosz(ipy) > cosz_min) then
            cgrid%met(ipy)%rshort_diffuse = 0.
            cgrid%met(ipy)%nir_beam       = cgrid%met(ipy)%nir_diffuse                     &
                                          + cgrid%met(ipy)%nir_beam
            cgrid%met(ipy)%par_beam       = cgrid%met(ipy)%par_diffuse                     &
                                          + cgrid%met(ipy)%par_beam
            cgrid%met(ipy)%nir_diffuse    = 0.
            cgrid%met(ipy)%par_diffuse    = 0.
         else
            cgrid%met(ipy)%rshort_diffuse = cgrid%met(ipy)%rshort
            cgrid%met(ipy)%nir_diffuse    = cgrid%met(ipy)%nir_diffuse                     &
                                          + cgrid%met(ipy)%nir_beam
            cgrid%met(ipy)%par_diffuse    = cgrid%met(ipy)%par_diffuse                     &
                                          + cgrid%met(ipy)%par_beam
            cgrid%met(ipy)%nir_beam       = 0.
            cgrid%met(ipy)%par_beam       = 0.
         end if
         !---------------------------------------------------------------------------------!

      case (5)
         !---------------------------------------------------------------------------------!
         !     Re-calculate all components using the clearness index method.               !
         !---------------------------------------------------------------------------------!
         call short_bdown_clearidx(cgrid%met(ipy)%rshort                                   &
                                  ,cgrid%cosz(ipy)                                         &
                                  ,cgrid%met(ipy)%par_beam                                 &
                                  ,cgrid%met(ipy)%par_diffuse                              &
                                  ,cgrid%met(ipy)%nir_beam                                 &
                                  ,cgrid%met(ipy)%nir_diffuse                              &
                                  ,cgrid%met(ipy)%rshort_diffuse )
         !----- Make sure the total radiation is preserved. -------------------------------!
         cgrid%met(ipy)%rshort = cgrid%met(ipy)%par_beam + cgrid%met(ipy)%par_diffuse      &
                               + cgrid%met(ipy)%nir_beam + cgrid%met(ipy)%nir_diffuse
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine solar_radiation_breakdown
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine is an adaptation from the SiB2 model (Sellers et al. 1986), and   !
   ! it splits the radiation into diffuse and direct radiation.  This is currently used to !
   ! convert the radiation from BRAMS.  In the future, we should use the values that come  !
   ! from the radiation schemes.                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine short2diff_sib(rshort_tot,cosz,rshort_diff)
      use canopy_radiation_coms, only : cosz_min  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real, intent(in)    :: rshort_tot  ! Surface incident shortwave radiation
      real, intent(in)    :: cosz        ! Cosine of the zenith distance
      real, intent(out)   :: rshort_diff ! Surface incident diffuse shortwave radiation.
      !----- Local variables. -------------------------------------------------------------!
      real                :: cloud
      real                :: difrat
      real                :: vnrat
      !----- Local constants. -------------------------------------------------------------!
      real, parameter     :: c1 = 580.
      real, parameter     :: c2 = 464.
      real, parameter     :: c3 = 499.
      real, parameter     :: c4 = 963.
      real, parameter     :: c5 = 1160.
      !------------------------------------------------------------------------------------!

      if (cosz > cosz_min) then

         cloud       = min(1.,max(0.,(c5 * cosz - rshort_tot) / (c4 * cosz)))
         difrat      = min(1.,max(0.,0.0604 / ( cosz -0.0223 ) + 0.0683))

         difrat      = difrat + ( 1. - difrat ) * cloud
         vnrat       = ( c1 - cloud*c2 ) / ( ( c1 - cloud*c3 ) + ( c1 - cloud*c2 ))

         rshort_diff = difrat*vnrat*rshort_tot
      else
         rshort_diff = rshort_tot
      end if

      return
   end subroutine short2diff_sib
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine computes the split between direct and diffuse radiation, and     !
   ! between visible and near-infrared radiation using the method suggested by:            !
   !                                                                                       !
   ! Weiss, A., J. M. Norman, 1985: Partitioning solar radiation into direct and diffuse,  !
   !     visible and near-infrared components.  Agric. For. Meteorol., 34, 205-213. (WN85) !
   !---------------------------------------------------------------------------------------!
   subroutine short_bdown_weissnorman(rshort_full,atm_prss,cosz,par_beam,par_diff,nir_beam &
                                     ,nir_diff,rshort_diff)
      use consts_coms          , only : solar         & ! intent(in)
                                      , prefsea       & ! intent(in)
                                      , twothirds     ! ! intent(in)
      use canopy_radiation_coms, only : cosz_min      & ! intent(in)
                                      , fvis_beam_def & ! intent(in)
                                      , fvis_diff_def & ! intent(in)
                                      , fnir_beam_def & ! intent(in)
                                      , fnir_diff_def ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real, intent(in)    :: rshort_full       ! Incident SW radiation   (total)   [  W/m2]
      real, intent(in)    :: atm_prss          ! Atmospheric pressure              [    Pa]
      real, intent(in)    :: cosz              ! cos(zenith distance)              [   ---]
      real, intent(out)   :: par_beam          ! Incident PAR            (direct ) [  W/m2]
      real, intent(out)   :: par_diff          ! Incident PAR            (diffuse) [  W/m2]
      real, intent(out)   :: nir_beam          ! Incident near-infrared  (direct ) [  W/m2]
      real, intent(out)   :: nir_diff          ! Incident near-infrared  (diffuse) [  W/m2]
      real, intent(out)   :: rshort_diff       ! Incident SW radiation   (diffuse) [  W/m2]
      !----- Local variables. ----------------------------------------------------------------!
      real                :: par_beam_top      ! PAR at the top of atmosphere      [  W/m2]
      real                :: nir_beam_top      ! NIR at the top of atmosphere      [  W/m2]
      real                :: par_beam_pot      ! Potential incident PAR  (direct)  [  W/m2]
      real                :: par_diff_pot      ! Potential incident PAR  (diffuse) [  W/m2]
      real                :: par_full_pot      ! Potential incident PAR  (total)   [  W/m2]
      real                :: nir_beam_pot      ! Potential  PAR          (direct)  [  W/m2]
      real                :: nir_diff_pot      ! Potential  PAR          (diffuse) [  W/m2]
      real                :: nir_full_pot      ! Potential  PAR          (total)   [  W/m2]
      real                :: par_full          ! Actual incident PAR     (total)   [  W/m2]
      real                :: nir_full          ! Actual incident NIR     (total)   [  W/m2]
      real                :: fvis_beam_act     ! Actual beam fraction - PAR        [   ---]
      real                :: fnir_beam_act     ! Actual beam fraction - NIR        [   ---]
      real                :: fvis_diff_act     ! Actual diffuse fraction - PAR     [   ---]
      real                :: fnir_diff_act     ! Actual diffuse fraction - NIR     [   ---]
      real                :: ratio             ! Ratio between obs. and expected   [   ---]
      real                :: aux_par           ! Auxiliary variable                [   ---]
      real                :: aux_nir           ! Auxiliary variable                [   ---]
      real                :: secz              ! sec(zenith distance)              [   ---]
      real                :: log10secz         ! log10[sec(zenith distance)]       [   ---]
      real                :: w10               ! Minimum Water absorption          [  W/m2]
      !------------------------------------------------------------------------------------!
      !    Local constants.                                                                !
      !------------------------------------------------------------------------------------!
      !----- Extinction coefficient. (equations 1 and 4 of WN85) --------------------------!
      real              , parameter :: par_beam_expext  = -0.185
      real              , parameter :: nir_beam_expext  = -0.060
      !----- This is the typical conversion of diffuse radiation in sunny days. -----------!
      real              , parameter :: par2diff_sun = 0.400
      real              , parameter :: nir2diff_sun = 0.600
      !----- Coefficients for various equations in WN85. ----------------------------------!
      real, dimension(3), parameter :: wn85_06 = (/ -1.1950, 0.4459, -0.0345 /)
      real, dimension(2), parameter :: wn85_11 = (/  0.90, 0.70 /)
      real, dimension(2), parameter :: wn85_12 = (/  0.88, 0.68 /)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     First thing to check is whether this is daytime or "night-time".  If the       !
      ! zenith angle is too close to horizon, we assume it's dawn/dusk and all radiation   !
      ! goes to diffuse.                                                                   !
      !------------------------------------------------------------------------------------!
      if (cosz <= cosz_min) then
         rshort_diff  = rshort_full
         par_beam     = 0.0
         nir_beam     = 0.0
         par_diff     = fvis_diff_def * rshort_diff
         nir_diff     = fnir_diff_def * rshort_diff
      else
         !----- Save 1/cos(zen), which is the secant.  We will use this several times. ----!
         secz      = 1. / cosz
         log10secz = log10(secz)
         !---------------------------------------------------------------------------------!


        !----- Total radiation at the top [  W/m2], using ED defaults. -------------------!
         par_beam_top = fvis_beam_def * solar
         nir_beam_top = fnir_beam_def * solar
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Find the potential PAR components (beam, diffuse, total), using equations 1, !
         ! 3, and 9 of WN85.                                                               !
         !---------------------------------------------------------------------------------!
         par_beam_pot = par_beam_top                                                       &
                      * exp ( par_beam_expext * (atm_prss / prefsea) * secz) * cosz
         par_diff_pot = par2diff_sun * (par_beam_top - par_beam_pot) * cosz
         par_full_pot = par_beam_pot + par_diff_pot
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the NIR absorption of 10 mm of precipitable water, using WN85 eqn. 6.  !
         !---------------------------------------------------------------------------------!
         w10 = solar * 10 ** (wn85_06(1) + log10secz * (wn85_06(2) + wn85_06(3)*log10secz))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the potential direct and diffuse near-infrared radiation, using eqn.   !
         ! 4, 5, and 10 of WN85.                                                           !
         !---------------------------------------------------------------------------------!
         nir_beam_pot = ( nir_beam_top                                                     &
                        * exp ( nir_beam_expext * (atm_prss/prefsea) * secz) - w10 ) * cosz
         nir_diff_pot = nir2diff_sun * ( nir_beam_top - nir_beam_pot - w10 ) * cosz
         nir_full_pot = nir_beam_pot + nir_diff_pot
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the actual total for PAR and NIR, using equations 7 and 8.             !
         !---------------------------------------------------------------------------------!
         ratio    = rshort_full / (par_full_pot + nir_full_pot)
         par_full = ratio * par_full_pot
         nir_full = ratio * nir_full_pot
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the fraction of PAR and NIR that stays as beam, using equations 11 and !
         ! 12 of WN85.                                                                     !
         !---------------------------------------------------------------------------------!
         !----- Make sure that the ratio is bounded. --------------------------------------!
         aux_par       = min(wn85_11(1), max(0., ratio))
         aux_nir       = min(wn85_12(1), max(0., ratio))
         fvis_beam_act = min(1., max(0., par_beam_pot                                      &
                                    * (1. - ((wn85_11(1)-aux_par)/wn85_11(2))**twothirds)  &
                                    / par_full_pot ) )
         fnir_beam_act = min(1., max(0., nir_beam_pot                                      &
                                    * (1. - ((wn85_12(1)-aux_nir)/wn85_12(2))**twothirds)  &
                                    / nir_full_pot ) )
         fvis_diff_act = 1. - fvis_beam_act
         fnir_diff_act = 1. - fnir_beam_act
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the radiation components.                                              !
         !---------------------------------------------------------------------------------!
         par_beam    = fvis_beam_act * par_full
         par_diff    = fvis_diff_act * par_full
         nir_beam    = fnir_beam_act * nir_full
         nir_diff    = fnir_diff_act * nir_full
         rshort_diff = par_diff + nir_diff
         !---------------------------------------------------------------------------------!
      end if

      return
   end subroutine short_bdown_weissnorman
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine computes the split between direct and diffuse radiation, and     !
   ! between visible and near-infrared radiation using a combination of methods suggested  !
   ! by:                                                                                   !
   !                                                                                       !
   ! Boland, J., B. Ridley, B. Brown, 2008: Models of diffuse solar radiation.  Renew.     !
   !     Energy, 33, 575-584. doi:10.1016/j.renene.2007.04.012 (BD08).                     !
   !                                                                                       !
   ! Tsubo, M., S. Walker, 2005: Relationships between photosynthetically active radiation !
   !     and clearness index at Bloemfontein, South Africa.  Theor. Appl. Climatol., 80,   !
   !     17-25.  doi:10.1007/s00704-004-0080-5 (TW05).                                     !
   !                                                                                       !
   ! Bendix, J., B. Silva, K. Roos, D. O. Gottlicher, R. Rollenbeck, T. Nauß, E. Beck,     !
   !     2010: Model parameterization to simulate and compare the PAR absorption potential !
   !     of two competing plant species.  Int. J. Biometeorol. 54, 283-295.                !
   !     doi:10.1007/s00484-009-0279-3 (BX10).                                             !
   !---------------------------------------------------------------------------------------!
   subroutine short_bdown_clearidx(rshort_full,cosz,par_beam,par_diff,nir_beam,nir_diff    &
                                  ,rshort_diff)
      use consts_coms          , only : solar         & ! intent(in)
                                      , lnexp_min     & ! intent(in)
                                      , lnexp_max     ! ! intent(in)
      use canopy_radiation_coms, only : cosz_min      ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real, intent(in)    :: rshort_full       ! Incident SW radiation   (total)   [  W/m2]
      real, intent(in)    :: cosz              ! cos(zenith distance)              [   ---]
      real, intent(out)   :: par_beam          ! Incident PAR            (direct ) [  W/m2]
      real, intent(out)   :: par_diff          ! Incident PAR            (diffuse) [  W/m2]
      real, intent(out)   :: nir_beam          ! Incident near-infrared  (direct ) [  W/m2]
      real, intent(out)   :: nir_diff          ! Incident near-infrared  (diffuse) [  W/m2]
      real, intent(out)   :: rshort_diff       ! Incident SW radiation   (diffuse) [  W/m2]
      !----- Local variables. -------------------------------------------------------------!
      real                :: rshort_toa        ! Incident SW at the top of atmos.  [  W/m2]
      real                :: par_full          ! Incident PAR (Direct+Diffuse)     [  W/m2]
      real                :: nir_full          ! Incident NIR (Direct+Diffuse)     [  W/m2]
      real                :: fkt               ! Clearness index                   [   ---]
      real                :: fpar              ! PAR/SW ratio                      [   ---]
      real                :: fdiff_aux         ! Dummy      [   ---]
      real                :: fdiff_1st         ! Uncorrected diffuse fraction      [   ---]
      real                :: fdiff             ! Corrected diffuse fraction        [   ---]
      !------------------------------------------------------------------------------------!
      !    Local constants.                                                                !
      !------------------------------------------------------------------------------------!
      !----- Coefficients for Diffuse/total fraction. (equations 24,26, 32 of BD10) -------!
      real, dimension(2), parameter :: bd08_eqn32 = (/ -5.0033,  8.6025 /)
      !----- PAR/SW fraction coefficients (equation 2 of TW05). ---------------------------!
      real, dimension(3), parameter :: tw05_eqn02 = (/  0.613, -0.334,  0.121 /)
      !----- Flag to activate BX10's correction for sun angle. ----------------------------!
      logical                       :: apply_bx10_corr = .false.
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     First thing to check is whether this is daytime or "night-time".  If the       !
      ! zenith angle is too close to horizon, we assume it's dawn/dusk and all radiation   !
      ! goes to diffuse.                                                                   !
      !------------------------------------------------------------------------------------!
      if (cosz <= cosz_min) then
         rshort_diff  = rshort_full
         par_beam     = 0.0
         nir_beam     = 0.0
         par_diff     = tw05_eqn02(1) * rshort_diff
         nir_diff     = rshort_diff - par_diff
      else
         !----- Total radiation at the top of the atmosphere [  W/m2], using ED defaults. -!
         rshort_toa = solar * cosz
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the clearness index based on total shortwave radiation.                 !
         !---------------------------------------------------------------------------------!
         fkt = max(0.,min(1.,rshort_full / rshort_toa))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the fraction of PAR radiation as a function of total radiation and     !
         ! clearness index, following TW05, eqn. 2.                                        !
         !---------------------------------------------------------------------------------!
         fpar     = tw05_eqn02(1) + fkt * ( tw05_eqn02(2) + fkt * tw05_eqn02(3) )
         par_full = fpar * rshort_full
         nir_full = rshort_full - par_full
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the uncorrected diffuse fraction based on BD08.  We use this equation  !
         ! instead of the method proposed by BX10 because it is a continuous function and  !
         ! results are very similar to BX10 curve except at high fkt, where BD08 predicts  !
         ! higher fraction of direct -- a reasonable assumption because total radiation    !
         ! should be entirely direct if the it is the same as the TOA radiation.           !
         !---------------------------------------------------------------------------------!
         fdiff_aux = max(lnexp_min, min(lnexp_max, bd08_eqn32(1) + bd08_eqn32(2) * fkt ) )
         fdiff_1st = 1.0 / (1.0 + exp(fdiff_aux) )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     One option is to apply the correction term by BX10.  I'm still not sure if  !
         ! this term is really necessary, because the higher proportion of diffuse         !
         ! radiation at low sun angles should be due to thicker optical depth, and a lower !
         ! ratio between incident radiation at the ground and extraterrestrial radiation   !
         ! should be naturally lower.  If direct radiation is too high, then we may        !
         ! include this term.                                                              !
         !---------------------------------------------------------------------------------!
         if (apply_bx10_corr) then
            fdiff = fdiff_1st / ( (1.0-fdiff_1st) * cosz + fdiff_1st)
         else
            fdiff = fdiff_1st
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the radiation components.                                              !
         !---------------------------------------------------------------------------------!
         par_diff    = fdiff * par_full
         par_beam    = par_full - par_diff
         nir_diff    = fdiff * nir_full
         nir_beam    = nir_full - nir_diff
         rshort_diff = par_diff + nir_diff
         !---------------------------------------------------------------------------------!
      end if

      return
   end subroutine short_bdown_clearidx
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine updates the 10-day running average of radiation, which is used for !
   ! phenology.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine update_rad_avg(cgrid)
      use ed_state_vars , only : edtype      & ! structure
                               , polygontype & ! structure
                               , sitetype    ! ! structure
      use phenology_coms, only : radavg_wgt  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target    :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer   :: cpoly
      integer                      :: ipy
      integer                      :: isi
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the running averages.                                                   !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         siteloop: do isi = 1,cpoly%nsites
            cpoly%rad_avg(isi) = (1.0 - radavg_wgt) * cpoly%rad_avg(isi)                   &
                               +        radavg_wgt  * cpoly%met(isi)%rshort
         end do siteloop
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine update_rad_avg
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine prints the radiation information if the user wants detailed      !
   ! output.                                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine dump_radinfo(cgrid,ipy,mprev,mnext,prevmet_timea,nextmet_timea,met_frq       &
                          ,wprev,wnext,secz_prev,secz_next,fperp_prev,fperp_next           &
                          ,night_prev,night_next,variable)
      use ed_state_vars        , only : edtype        ! ! structure
      use ed_max_dims          , only : str_len       ! ! intent(in)
      use ed_misc_coms         , only : simtime       & ! intent(in)
                                      , current_time  ! ! intent(in)
      use met_driver_coms      , only : nbdsf_file    & ! intent(in)
                                      , nddsf_file    & ! intent(in)
                                      , vbdsf_file    & ! intent(in)
                                      , vddsf_file    ! ! intent(in)
      use canopy_radiation_coms, only : cosz_min      ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)          , target     :: cgrid
      integer               , intent(in) :: ipy
      type(simtime)         , intent(in) :: prevmet_timea
      type(simtime)         , intent(in) :: nextmet_timea
      integer               , intent(in) :: mprev
      integer               , intent(in) :: mnext
      real                  , intent(in) :: met_frq
      real                  , intent(in) :: wprev
      real                  , intent(in) :: wnext
      real                  , intent(in) :: secz_prev
      real                  , intent(in) :: secz_next
      real                  , intent(in) :: fperp_prev
      real                  , intent(in) :: fperp_next
      logical               , intent(in) :: night_prev
      logical               , intent(in) :: night_next
      character(len=5)      , intent(in) :: variable
      !----- Local variables. -------------------------------------------------------------!
      character(len=str_len)             :: thisfile
      character(len=22)                  :: when_prev
      character(len=22)                  :: when_next
      character(len=22)                  :: when_now
      real                               :: flux_prev
      real                               :: flux_next
      real                               :: flux_now
      real                               :: fperp_now
      logical                            :: night_now
      logical                            :: itsthere
      !----- Local constants. -------------------------------------------------------------!
      character(len=18)     , parameter  :: fmth = '(3(a,3x),15(a,1x))'
      character(len=35)     , parameter  :: fmtd = '(3(a,3x),3(11x,l1,1x),12(f12.3,1x))'
      character(len=29)     , parameter  :: fmtt = '(2(i2.2,a),i4.4,1x,3(i2.2,a))'
      !----- Locally saved variables, to delete and re-create scratch files. --------------!
      logical               , save       :: first_nbdsf = .true.
      logical               , save       :: first_nddsf = .true.
      logical               , save       :: first_vbdsf = .true.
      logical               , save       :: first_vddsf = .true.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    This copies the information to some common variables.                           !
      !------------------------------------------------------------------------------------!
      select case (variable)
      case ('nbdsf')

         flux_prev = cgrid%metinput(ipy)%nbdsf(mprev)
         flux_next = cgrid%metinput(ipy)%nbdsf(mnext)
         flux_now  = cgrid%met(ipy)%nir_beam
         thisfile  = nbdsf_file

         if (first_nbdsf) then

            !------------------------------------------------------------------------------!
            !     Check whether there was a file in there.  In case there was, delete and  !
            ! create a new header.                                                         !
            !------------------------------------------------------------------------------!
            inquire (file=trim(thisfile),exist=itsthere)
            if (itsthere) then
               open (unit=87,file=trim(thisfile),status='old',action='write')
               close(unit=87,status='delete')
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Make the header.                                                         !
            !------------------------------------------------------------------------------!
            open (unit=87,file=trim(thisfile),status='new',action='write')
            write (unit=87,fmt=fmth)                                                       &
                                  '             WHEN_PREV','              WHEN_NOW'        &
                                 ,'             WHEN_NEXT','  NIGHT_PREV','   NIGHT_NOW'   &
                                          ,'  NIGHT_NEXT','     MET_FRQ','       WPREV'    &
                                          ,'       WNEXT','   FLUX_PREV','   SECZ_PREV'    &
                                          ,'  FPERP_PREV','   FLUX_NEXT','   SECZ_NEXT'    &
                                          ,'  FPERP_NEXT','    FLUX_NOW','    COSZ_NOW'    &
                                          ,'   FPERP_NOW'
           close (unit=87,status='keep')
            !------------------------------------------------------------------------------!

            first_nbdsf = .false.
         end if

      case ('nddsf')

         flux_prev = cgrid%metinput(ipy)%nddsf(mprev)
         flux_next = cgrid%metinput(ipy)%nddsf(mnext)
         flux_now  = cgrid%met(ipy)%nir_diffuse
         thisfile  = nddsf_file

         if (first_nddsf) then

            !------------------------------------------------------------------------------!
            !     Check whether there was a file in there.  In case there was, delete and  !
            ! create a new header.                                                         !
            !------------------------------------------------------------------------------!
            inquire (file=trim(thisfile),exist=itsthere)
            if (itsthere) then
               open (unit=87,file=trim(thisfile),status='old',action='write')
               close(unit=87,status='delete')
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Make the header.                                                         !
            !------------------------------------------------------------------------------!
            open (unit=87,file=trim(thisfile),status='new',action='write')
            write (unit=87,fmt=fmth)                                                       &
                                  '             WHEN_PREV','              WHEN_NOW'        &
                                 ,'             WHEN_NEXT','  NIGHT_PREV','   NIGHT_NOW'   &
                                          ,'  NIGHT_NEXT','     MET_FRQ','       WPREV'    &
                                          ,'       WNEXT','   FLUX_PREV','   SECZ_PREV'    &
                                          ,'  FPERP_PREV','   FLUX_NEXT','   SECZ_NEXT'    &
                                          ,'  FPERP_NEXT','    FLUX_NOW','    COSZ_NOW'    &
                                          ,'   FPERP_NOW'
            close (unit=87,status='keep')
            !------------------------------------------------------------------------------!

            first_nddsf = .false.
         end if
      case ('vbdsf')

         flux_prev = cgrid%metinput(ipy)%vbdsf(mprev)
         flux_next = cgrid%metinput(ipy)%vbdsf(mnext)
         flux_now  = cgrid%met(ipy)%par_beam
         thisfile  = vbdsf_file

         if (first_vbdsf) then

            !------------------------------------------------------------------------------!
            !     Check whether there was a file in there.  In case there was, delete and  !
            ! create a new header.                                                         !
            !------------------------------------------------------------------------------!
            inquire (file=trim(thisfile),exist=itsthere)
            if (itsthere) then
               open (unit=87,file=trim(thisfile),status='old',action='write')
               close(unit=87,status='delete')
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Make the header.                                                         !
            !------------------------------------------------------------------------------!
            open (unit=87,file=trim(thisfile),status='new',action='write')
            write (unit=87,fmt=fmth)                                                       &
                                  '             WHEN_PREV','              WHEN_NOW'        &
                                 ,'             WHEN_NEXT','  NIGHT_PREV','   NIGHT_NOW'   &
                                          ,'  NIGHT_NEXT','     MET_FRQ','       WPREV'    &
                                          ,'       WNEXT','   FLUX_PREV','   SECZ_PREV'    &
                                          ,'  FPERP_PREV','   FLUX_NEXT','   SECZ_NEXT'    &
                                          ,'  FPERP_NEXT','    FLUX_NOW','    COSZ_NOW'    &
                                          ,'   FPERP_NOW'
            close (unit=87,status='keep')
            !------------------------------------------------------------------------------!

            first_vbdsf = .false.
         end if

      case ('vddsf')

         flux_prev = cgrid%metinput(ipy)%vddsf(mprev)
         flux_next = cgrid%metinput(ipy)%vddsf(mnext)
         flux_now  = cgrid%met(ipy)%par_diffuse
         thisfile  = vddsf_file

         if (first_vddsf) then

            !------------------------------------------------------------------------------!
            !     Check whether there was a file in there.  In case there was, delete and  !
            ! create a new header.                                                         !
            !------------------------------------------------------------------------------!
            inquire (file=trim(thisfile),exist=itsthere)
            if (itsthere) then
               open (unit=87,file=trim(thisfile),status='old',action='write')
               close(unit=87,status='delete')
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Make the header.                                                         !
            !------------------------------------------------------------------------------!
            open (unit=87,file=trim(thisfile),status='new',action='write')
            write (unit=87,fmt=fmth)                                                       &
                                  '             WHEN_PREV','              WHEN_NOW'        &
                                 ,'             WHEN_NEXT','  NIGHT_PREV','   NIGHT_NOW'   &
                                          ,'  NIGHT_NEXT','     MET_FRQ','       WPREV'    &
                                          ,'       WNEXT','   FLUX_PREV','   SECZ_PREV'    &
                                          ,'  FPERP_PREV','   FLUX_NEXT','   SECZ_NEXT'    &
                                          ,'  FPERP_NEXT','    FLUX_NOW','    COSZ_NOW'    &
                                          ,'   FPERP_NOW'
            close (unit=87,status='keep')
            !------------------------------------------------------------------------------!

            first_vddsf = .false.
         end if
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find some common variables.                                                   !
      !------------------------------------------------------------------------------------!
      write (when_prev,fmt=fmtt) prevmet_timea%month,'/',prevmet_timea%date,'/'            &
                             ,prevmet_timea%year,prevmet_timea%hour,'.',prevmet_timea%min  &
                             ,'.',prevmet_timea%sec,'UTC'
      write (when_now ,fmt=fmtt) current_time%month,'/',current_time%date,'/'              &
                             ,current_time%year,current_time%hour,'.',current_time%min,'.' &
                             ,current_time%sec,'UTC'
      write (when_next,fmt=fmtt) nextmet_timea%month,'/',nextmet_timea%date,'/'            &
                             ,nextmet_timea%year,nextmet_timea%hour,'.',nextmet_timea%min  &
                             ,'.',nextmet_timea%sec,'UTC'
      fperp_now = fperp_prev * wprev + fperp_next * wnext
      night_now = cgrid%cosz(ipy) <= cosz_min
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Write the information on the file.                                             !
      !------------------------------------------------------------------------------------!
      open(unit=87,file=trim(thisfile),status='old',action='write',position='append')
      write(unit=87,fmt=fmtd) when_prev,when_now,when_next,night_prev,night_now,night_next &
                             ,met_frq,wprev,wnext,flux_prev,secz_prev,fperp_prev,flux_next &
                             ,secz_next,fperp_next,flux_now,cgrid%cosz(ipy),fperp_now
      close(unit=87,status='keep')
      !------------------------------------------------------------------------------------!
      return
   end subroutine dump_radinfo
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates angle of incidence based on local slope and aspect.    !
   !---------------------------------------------------------------------------------------!
   subroutine angle_of_incid(aoi,cosz,solar_hour_aspect,slope,terrain_aspect)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real, intent(in)  :: cosz              ! Cosine of zenithal angle
      real, intent(in)  :: slope             ! Terrain slope
      real, intent(in)  :: solar_hour_aspect ! horizontal location of the sun defined with
                                             !    the same reference as terrain aspect.
      real, intent(in)  :: terrain_aspect    ! Terrain aspect
      real, intent(out) :: aoi               ! Angle of incidence
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)      :: cosz8     ! Double prec. counterpart of cosz
      real(kind=8)      :: sinz8     ! Sine of zenithal angle
      real(kind=8)      :: slope8    ! Double prec. counterpart of slope
      real(kind=8)      :: sh_asp8   ! Double prec. counterpart of solar_hour_aspect
      real(kind=8)      :: terr_asp8 ! Double prec. counterpart of terrain_aspect
      real(kind=8)      :: aoi8      ! Double prec. counterpart of aoi
      !----- Local parameters. ------------------------------------------------------------!
      real(kind=8), parameter :: tiny_offset=1.d-20
      !----- External functions. ----------------------------------------------------------!
      real        , external  :: sngloff
      !------------------------------------------------------------------------------------!

      cosz8     = dble(cosz)
      sinz8     = sqrt(1.d0-cosz8*cosz8)
      slope8    = dble(slope)
      sh_asp8   = dble(solar_hour_aspect)
      terr_asp8 = dble(terrain_aspect)
      if (cosz8 < 0.d0) then
         aoi8 = 0.d0
      else
         aoi8 = max(0.d0, cosz8*dcos(slope8) + sinz8*dsin(slope8)*dcos(sh_asp8-terr_asp8))
      end if

      aoi = sngloff(aoi8,tiny_offset)

      return
   end subroutine angle_of_incid
   !=======================================================================================!
   !=======================================================================================!









   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the cosine of the zenith angle.                            !
   !---------------------------------------------------------------------------------------!
   real function ed_zen(plon,plat,when)
      use ed_misc_coms , only : simtime      ! ! structure
      use consts_coms  , only : pio1808      & ! intent(in)
                              , twopi8       & ! intent(in)
                              , hr_sec8      & ! intent(in)
                              , tiny_num8    ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4) , intent(in) :: plon
      real(kind=4) , intent(in) :: plat
      type(simtime), intent(in) :: when
      !------ Local variables. ------------------------------------------------------------!
      integer                   :: doy     ! Day of year ("Julian" day)
      real(kind=8)              :: declin  ! Declination
      real(kind=8)              :: sdec    ! Sine of declination
      real(kind=8)              :: cdec    ! Cosine of declination
      real(kind=8)              :: dayhr   ! Hour of day 
      real(kind=8)              :: radlat  ! Latitude in radians
      real(kind=8)              :: clat    ! Cosine of latitude 
      real(kind=8)              :: slat    ! Sine of latitude
      real(kind=8)              :: dayhrr  ! Hour of day in radians
      real(kind=8)              :: hrangl  ! Hour angle
      !----- Local constants. -------------------------------------------------------------!
      real(kind=8), parameter :: capri    = -2.344d1 ! Tropic of Capricornium latitude
      real(kind=8), parameter :: ndaysnl  =  3.65d2  ! # of days of year (no leap years)
      real(kind=8), parameter :: ndayslp  =  3.66d2  ! # of days of year (leap years)
      integer     , parameter :: shsummer = -10      ! DoY of Southern Hemisphere summer
      !----- External functions. ----------------------------------------------------------!
      integer     , external  :: julday  ! Function to find day of year ("Julian" day)
      logical     , external  :: isleap  ! Function to determine whether the year is leap
      real        , external  :: sngloff ! Function to safely convert double to single prec.
      !------------------------------------------------------------------------------------!



      !----- Find the day of the year. ----------------------------------------------------!
      doy    = julday(when%month, when%date, when%year)
      !------------------------------------------------------------------------------------!



      !----- Find the hour angle, then get cosine of zenith angle. ------------------------!
      dayhr = dble(when%time) / hr_sec8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   declin is the solar latitude in degrees (also known as declination).             !
      !   sdec - sine of declination                                                       !
      !   cdec - cosine of declination.                                                    !
      !------------------------------------------------------------------------------------!
      if (isleap(when%year)) then
         declin = capri * cos(twopi8 * dble(doy - shsummer) / ndaysnl) * pio1808
      else
         declin = capri * cos(twopi8 * dble(doy - shsummer) / ndayslp) * pio1808
      end if
      sdec   = dsin(declin)
      cdec   = dcos(declin)
      !------------------------------------------------------------------------------------!

      !----- Find the latitude in radians. ------------------------------------------------!
      radlat = dble(plat) * pio1808
      clat   = dcos(radlat)
      slat   = dsin(radlat)
      !------------------------------------------------------------------------------------!

      !------ Find the hour angle. --------------------------------------------------------!
      dayhrr = dmod(dayhr+dble(plon)/1.5d1+2.4d1,2.4d1)
      hrangl = 1.5d1 * (dayhrr - 1.2d1) * pio1808

      ed_zen = sngloff(slat * sdec + clat * cdec * dcos(hrangl),tiny_num8)

      return
   end function ed_zen
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the average secant of the daytime zenith angle.  In case   !
   ! the period of integration is that accounts for the zenith angle is 0. or less than    !
   ! one time step, then the average is the actual value.  Night-time periods are ignored  !
   ! and if there is no daytime value, then we set it to 0.                                !
   !---------------------------------------------------------------------------------------!
   real function mean_daysecz(plon,plat,whena,dt,tmax)
      use update_derived_utils , only : update_model_time_dm ! ! sub-routine
      use ed_misc_coms         , only : simtime              ! ! structure
      use canopy_radiation_coms, only : cosz_min             ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4) , intent(in) :: plon
      real(kind=4) , intent(in) :: plat
      type(simtime), intent(in) :: whena
      real(kind=4) , intent(in) :: dt
      real(kind=4) , intent(in) :: tmax
      !------ Local variables. ------------------------------------------------------------!
      type(simtime)             :: now          ! Current time
      integer                   :: is           ! Step counter
      integer                   :: nsteps       ! Number of steps to perform the average
      real                      :: dtfit        ! Delta-t that nicely fits within tmax
      real                      :: dtnow        ! Delta-t for this time
      real                      :: cosz         ! Declination
      real                      :: daytot       ! Total time that was daytime
      real                      :: mean_daycosz ! Average cosine of zenith angle
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether tmax is less than the time step.  In case it is, we only have    !
      ! one time, so we don't need to do the average.                                      !
      !------------------------------------------------------------------------------------!
      if (dt >= tmax) then
         !----- Less than one time, no average necessary. ---------------------------------!
         cosz = ed_zen(plon,plat,whena)
         if (cosz > cosz_min) then
            mean_daysecz = 1.0 / cosz
         else
            !----- Night-time, set the mean to zero. --------------------------------------!
            mean_daysecz = 0.0
         end if

      else
         !---------------------------------------------------------------------------------!
         !     Several times, first find the number of steps, then the delta-t that fits   !
         ! nicely within the time span.                                                    !
         !---------------------------------------------------------------------------------!
         nsteps = ceiling(tmax / dt)
         dtfit  = tmax / real(nsteps)
         !---------------------------------------------------------------------------------!

         mean_daycosz = 0.0
         daytot       = 0.0
         do is=1,nsteps
            !----- Get the current time. --------------------------------------------------!
            now   = whena
            dtnow = dtfit * (real(is) - 0.5)
            call update_model_time_dm(now,dtnow)

            !----- Get the cosine of the zenith angle. ------------------------------------!
            cosz = ed_zen(plon,plat,now)

            !----- Add to the integral only if it this value is valid. --------------------!
            if (cosz > cosz_min) then
               mean_daycosz = mean_daycosz + dtfit * cosz 
               daytot       = daytot       + dtfit
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the normalisation factor.                                              !
         !---------------------------------------------------------------------------------!
         if (daytot > 0.0 .and. mean_daycosz > 0.0) then
            mean_daycosz = mean_daycosz / daytot
            mean_daysecz = 1.0 / mean_daycosz
         else
            mean_daysecz = 0.0
         end if
         !---------------------------------------------------------------------------------!
      end if

      return
   end function mean_daysecz
   !=======================================================================================!
   !=======================================================================================!
end module radiate_utils
!==========================================================================================!
!==========================================================================================!
