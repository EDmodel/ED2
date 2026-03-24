module stem_resp_driv
  contains

!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the stem respiration scheme                             !
! This is called every step, but not every sub-step.                                       !
!------------------------------------------------------------------------------------------!
subroutine stem_respiration(csite,ipa)
   use ed_state_vars  , only : sitetype           & ! structure
                             , patchtype          ! ! structure
   use pft_coms       , only : agf_bs             & ! intent(in)
                             , is_grass           ! ! intent(in)
   use consts_coms    , only : pi1                & ! intent(in)
                             , umols_2_kgCyr      ! ! intent(in)
   use ed_misc_coms   , only : dtlsm              & ! intent(in)
                             , frqsum             ! ! intent(in)
   use physiology_coms, only : istem_respiration_scheme ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite             ! Current site
   integer                   , intent(in)  :: ipa               ! Current patch #
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch             ! Current site
   integer                                 :: ico                ! Current cohort #
   integer                                 :: ipft
   real                                    :: stem_area        ! m2
   !----- Locally saved variables. --------------------------------------------------------!
   real                          , save    :: dtlsm_o_frqsum
   logical                       , save    :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- Assign the constant scaling factor. ---------------------------------------------!
   if (first_time) then
      first_time     = .false.
      dtlsm_o_frqsum = dtlsm / frqsum
   end if
   !---------------------------------------------------------------------------------------!


   !----- Point to the cohort structures --------------------------------------------------!
   cpatch => csite%patch(ipa)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Loop over all cohorts, from tallest to shortest.                                   !
   !---------------------------------------------------------------------------------------!
   cohortloop: do ico = 1,cpatch%ncohorts

        !----- Alias for PFT and root layer. ------------------------------------------!
        ipft  = cpatch%pft(ico)

        if ( is_grass(ipft) .or. istem_respiration_scheme == 0) then
            !----- Grass has no stem, set stem respiration to zero   ----------------------!
            cpatch%stem_respiration(ico) = 0.
        else
            !------------------------------------------------------------------------------!
            !  Woody PFTs, assume stem surface area is the cylindral stem plus WAI * pi    !
            !  The resulting stem surface area is comparable to reported values in the     !
            !  tropics, such as Chambers et al. 2004 Ecological Applicaitons               !
            !------------------------------------------------------------------------------!
            stem_area = ( cpatch%dbh(ico) * 1e-2 * cpatch%hite(ico) * cpatch%nplant(ico)   &
                        + cpatch%wai(ico) ) * pi1 / agf_bs(ipft)
            cpatch%stem_respiration(ico) = stem_resp_norm(ipft,cpatch%dbh(ico),cpatch%wood_temp(ico)) &
                                         * stem_area
                                         ! umol/m2 ground/s
        endif

        cpatch%today_stem_resp(ico) = cpatch%today_stem_resp(ico)                      &
                                    + cpatch%stem_respiration(ico) * dtlsm
                                    ! add this time step to the daily mean values
        
        !----- The output variable must be in [kgC/plant/yr]. -------------------------!
        cpatch%fmean_stem_resp(ico)  = cpatch%fmean_stem_resp (ico)                    &
                                        + cpatch%stem_respiration(ico)                    &
                                        * dtlsm_o_frqsum * umols_2_kgCyr                  &
                                        / cpatch%nplant          (ico)

   end do cohortloop
   !---------------------------------------------------------------------------------------!


   return
end subroutine stem_respiration
!==========================================================================================!
!==========================================================================================!

!==========================================================================================!
!==========================================================================================!
!     This function determines the normalised stem respiration (umol/m2 stem surface/s)    !
!------------------------------------------------------------------------------------------!
real function stem_resp_norm(ipft,dbh,wood_temp)
   use pft_coms       , only : stem_respiration_factor  & ! intent(in)
                             , stem_resp_size_scaler    & ! intent(in)
                             , srf_low_temp             & ! intent(in)
                             , srf_high_temp            & ! intent(in)
                             , srf_decay_elow           & ! intent(in)
                             , srf_decay_ehigh          & ! intent(in)
                             , srf_hor                  & ! intent(in)
                             , srf_q10                  ! ! intent(in)
   use farq_leuning   , only : arrhenius                & ! function
                             , collatz                  ! ! function
   use rk4_coms       , only : tiny_offset              ! ! intent(in)
   use physiology_coms, only : iphysiol                 ! ! intent(in)
   use consts_coms    , only : lnexp_min8               & ! intent(in)
                             , lnexp_max8               & ! intent(in)
                             , t008                     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in) :: ipft
   real(kind=4), intent(in) :: dbh
   real(kind=4), intent(in) :: wood_temp
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: wood_temp8
   real(kind=8)             :: srf08
   real(kind=8)             :: srf_low_temp8
   real(kind=8)             :: srf_high_temp8
   real(kind=8)             :: srf_decay_ehigh8
   real(kind=8)             :: srf_decay_elow8
   real(kind=8)             :: srf_hor8
   real(kind=8)             :: srf_q108
   real(kind=8)             :: lnexplow
   real(kind=8)             :: lnexphigh
   real(kind=8)             :: tlow_fun
   real(kind=8)             :: thigh_fun
   real(kind=8)             :: srf8
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)             :: sngloff
   !---------------------------------------------------------------------------------------!

   !----- Copy some variables to double precision temporaries. ----------------------------!
   wood_temp8      = dble(wood_temp                    )
   ! functional form from Chambers et al. 2004

   srf08           = 10. ** (                                                       &
                        dble(                                                       &
                         log10(stem_respiration_factor(ipft))                       & ! baseline
                        +stem_resp_size_scaler(ipft) * dbh                          & ! DBH effect
                        ))
   srf_low_temp8   = dble(srf_low_temp           (ipft)) + t008
   srf_high_temp8  = dble(srf_high_temp          (ipft)) + t008
   srf_decay_elow8  = dble(srf_decay_elow         (ipft))
   srf_decay_ehigh8 = dble(srf_decay_ehigh        (ipft))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Compute the functions that will control the Rrf function for low and high temper-  !
   ! ature.  In order to avoid floating point exceptions, we check whether the temperature !
   ! will make the exponential too large or too small.                                     !
   !---------------------------------------------------------------------------------------!
   !----- Low temperature. ----------------------------------------------------------------!
   lnexplow  = srf_decay_elow8 * (srf_low_temp8  - wood_temp8)
   lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
   tlow_fun  = 1.d0 +  exp(lnexplow)
   !----- High temperature. ---------------------------------------------------------------!
   lnexphigh = srf_decay_ehigh8 * (wood_temp8 - srf_high_temp8)
   lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
   thigh_fun = 1.d0 + exp(lnexphigh)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Decide which functional form to use based on the physiology.  This is just to make !
   ! it look similar to the leaf respiration respiration.                                  !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,1)
      srf_hor8 = dble(srf_hor(ipft))
      srf8     = arrhenius(wood_temp8,srf08,srf_hor8) / (tlow_fun * thigh_fun)
   case (2,3)
      srf_q108 = dble(srf_q10(ipft))
      srf8     = collatz(wood_temp8,srf08,srf_q108)   / (tlow_fun * thigh_fun)
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert result to single precision.                                               !
   !---------------------------------------------------------------------------------------!
   stem_resp_norm = sngloff(srf8,tiny_offset)
   !---------------------------------------------------------------------------------------!


   return
end function stem_resp_norm
!==========================================================================================!
!==========================================================================================!





end module stem_resp_driv
