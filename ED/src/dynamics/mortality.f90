!===================================================================
!=======================!
!===================================================================
!=======================!
module mortality
   !================================================================
   !=======================!
   !================================================================
   !=======================!


   contains



   !================================================================
      !=======================!
   !================================================================
      !=======================!
   !    This subroutine computes the total PFT-dependent mortality
      !     rate:                   !
   !-----------------------------------------------------------------
      !----------------------!
   subroutine mortality_rates(cpatch,ipa,ico,avg_daily_temp,&
         patch_age)
      use ed_state_vars , only : patchtype                  ! !
      !  Structure
      use pft_coms      , only : mort1                      & !
            !  intent(in)
                               , mort2                      & !
                               !  intent(in)
                               , mort3                      & !
                               !  intent(in)
                               , plant_min_temp             & !
                               !  intent(in)
                               , frost_mort                 ! !
      !  intent(in)
      use disturb_coms  , only : treefall_disturbance_rate  & !
            !  intent(in)
                               , treefall_hite_threshold    & !
                               !  intent(in)
                               , time2canopy                ! !
      !  intent(in)
      use ed_misc_coms  , only : current_time               ! !
      !  intent(in)
      use ed_max_dims   , only : n_pft                      ! !
      !  intent(in)
      use consts_coms   , only : lnexp_min                  & !
            !  intent(in)
                               , lnexp_max                  ! !
      !  intent(in)
      implicit none
      !----- Arguments ----------------------------------------------
      !----------------------!
      type(patchtype), target     :: cpatch          ! Current patch
      integer        , intent(in) :: ipa             ! Current patch
      !  ID
      integer        , intent(in) :: ico             ! Current cohort
      !  ID
      real           , intent(in) :: avg_daily_temp  ! Mean
      !  temperature yesterday
      real           , intent(in) :: patch_age
      !----- Local variables ----------------------------------------
      !----------------------!
      integer                     :: ipft            ! PFT 
      real                        :: temp_dep        ! Temp. function
      !  (frost mortality)
      real                        :: expmort         ! Carbon-balance
      !  term
      !--------------------------------------------------------------
      !----------------------!


      !----- Assume happy end, all plants survive... ----------------
      !----------------------!
      cpatch%mort_rate(1:4,ico) = 0.0
      ipft = cpatch%pft(ico)

      !--------------------------------------------------------------
      !----------------------!
      ! 1.  Ageing, PFT-dependent but otherwise constant.
      !                       !
      !--------------------------------------------------------------
      !----------------------!
      cpatch%mort_rate(1,ico) = mort3(ipft)
      !--------------------------------------------------------------
      !----------------------!


      !--------------------------------------------------------------
      !----------------------!
      ! 2.  Mortality rates due to negative carbon balance.
      !                       !
      !--------------------------------------------------------------
      !----------------------!
      expmort                 = max( lnexp_min, min( lnexp_max ,&
            mort2(ipft) * cpatch%cbr_bar(ico)))
      cpatch%mort_rate(2,ico) = mort1(ipft) / (1. + exp(expmort))
      
      !--------------------------------------------------------------
      !----------------------!


      !--------------------------------------------------------------
      !----------------------!
      ! 3.  Mortality due to treefall.
      !                       !
      !--------------------------------------------------------------
      !----------------------!
      if (cpatch%hite(ico) <= treefall_hite_threshold .and. patch_age&
            > time2canopy) then
         cpatch%mort_rate(3,ico) = treefall_disturbance_rate
      else
         cpatch%mort_rate(3,ico) = 0.
      end if
      !--------------------------------------------------------------
      !----------------------!



      !--------------------------------------------------------------
      !----------------------!
      ! 4.   Mortality due to cold, after:
      !                       !
      !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft,
      !       2006: The contributions !
      !           of land-use change, CO2 fertilization, and climate
      !            variability to the    !
      !           Eastern US carbon sink.  Glob. Change Biol., 12,
      !            2370-2390,              !
      !           doi: 10.1111/j.1365-2486.2006.01254.x
      !                                 !
      !--------------------------------------------------------------
      !----------------------!
      temp_dep = max( 0.0, min( 1.0, 1.0 - (avg_daily_temp -&
            plant_min_temp(ipft)) / 5.0) )
      cpatch%mort_rate(4,ico) = frost_mort(ipft) * temp_dep
      !--------------------------------------------------------------
      !----------------------!



      !--------------------------------------------------------------
      !----------------------!
      ! 5. Disturbance rate mortality.  This is not used by the
      !  cohort dynamics, instead   !
      !    this is just to account for the lost density due to the
      !     patch creation.  This   !
      !    mortality will be determined by the disturbance_mortality
      !     subroutine, not here. !
      !--------------------------------------------------------------
      !----------------------!
      !cpatch%mort_rate(5,ico) = TBD
      !--------------------------------------------------------------
      !----------------------!

      return
   end subroutine mortality_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine determines the mortality rates associated with the current        !
   ! disturbance.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine disturbance_mortality(csite,ipa,disturbance_rate,dest_type,poly_dest_type    &
                                   ,mindbh_harvest)
      use ed_state_vars, only : sitetype  & ! structure
                              , patchtype ! ! structure
      use ed_max_dims  , only : n_pft     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                   , target     :: csite
      integer                          , intent(in) :: ipa
      real                             , intent(in) :: disturbance_rate
      integer                          , intent(in) :: dest_type
      integer                          , intent(in) :: poly_dest_type
      real           , dimension(n_pft), intent(in) :: mindbh_harvest
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                  , pointer    :: cpatch
      integer                                       :: ico
      real                                          :: f_survival
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      do ico=1,cpatch%ncohorts
         f_survival = survivorship(dest_type,poly_dest_type,mindbh_harvest,csite,ipa,ico)
         cpatch%mort_rate(5,ico) = cpatch%mort_rate(5,ico)                                 &
                                 - log( f_survival                                         &
                                      + (1.0 - f_survival) * exp(- disturbance_rate) )
      end do
      return
   end subroutine disturbance_mortality
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the survivorship rate after some disturbance happens.      !
   !---------------------------------------------------------------------------------------!
   real function survivorship(dest_type,poly_dest_type,mindbh_harvest,csite,ipa,ico)
      use ed_state_vars, only : patchtype                & ! structure
                              , sitetype                 ! ! structure
      use disturb_coms , only : treefall_hite_threshold  & ! intent(in)
                              , fire_hite_threshold      ! ! intent(in)
      use pft_coms     , only : treefall_s_ltht          & ! intent(in)
                              , treefall_s_gtht          & ! intent(in)
                              , fire_s_ltht              & ! intent(in)
                              , fire_s_gtht              ! ! intent(in)
      use ed_max_dims  , only : n_pft                    ! ! intent(in)
      
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target     :: csite
      real          , dimension(n_pft), intent(in) :: mindbh_harvest
      integer                         , intent(in) :: ico
      integer                         , intent(in) :: ipa
      integer                         , intent(in) :: dest_type
      integer                         , intent(in) :: poly_dest_type
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                 , pointer    :: cpatch
      integer                                      :: ipft
      !------------------------------------------------------------------------------------!


      cpatch => csite%patch(ipa)
      ipft = cpatch%pft(ico)

      !----- Base the survivorship rates on the destination type. -------------------------!
      select case(dest_type)
      case (1) !----- Agriculture/cropland. -----------------------------------------------!
         survivorship = 0.0

      case (2) !----- Secondary land or forest plantation. --------------------------------!

         !----- Decide the fate based on the type of secondary disturbance. ---------------!
         select case (poly_dest_type)
         case (0) !----- Land abandonment, assume this is the last harvest. ---------------!
            survivorship = 0.0
         case (1) !----- Biomass logging, assume that nothing stays. ----------------------!
            survivorship = 0.0
         case (2) !----- Selective logging. -----------------------------------------------!
            !------------------------------------------------------------------------------!
            !     If the PFT DBH exceeds the minimum PFT for harvesting, the survivorship  !
            ! should be zero, otherwise, we assume survivorship similar to the treefall    !
            ! disturbance rate for short trees.                                            ! 
            !------------------------------------------------------------------------------!
            if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
               survivorship = 0.0
            else
               survivorship = treefall_s_ltht(ipft)
            end if
         end select

      case (3) !----- Primary land. -------------------------------------------------------!

         !----- Decide the fate based on the type of natural disturbance. -----------------!
         select case (poly_dest_type)
         case (0) !----- Treefall, we must check the cohort height. -----------------------!
            if (cpatch%hite(ico) < treefall_hite_threshold) then
               survivorship = treefall_s_ltht(ipft)
            else
               survivorship = treefall_s_gtht(ipft)
            end if

         case (1) !----- Fire, also check the cohort height. ------------------------------!
            if (cpatch%hite(ico) < fire_hite_threshold) then
               survivorship = fire_s_ltht(ipft)
            else
               survivorship = fire_s_gtht(ipft)
            end if
         end select
      end select

      return
   end function survivorship
   !=======================================================================================!
   !=======================================================================================!
end module mortality
!==========================================================================================!
!==========================================================================================!
