!==========================================================================================!
!==========================================================================================!
!    This subroutine will evaluate whether fire conditions exist, and if that is the case, !
! it will calculate the disturbance rate due to fire.                                      !
!------------------------------------------------------------------------------------------!
subroutine fire_frequency(cgrid)
   use ed_state_vars , only : edtype                 & ! structure
                            , polygontype            & ! structure
                            , sitetype               & ! structure
                            , patchtype              ! ! structure
   use ed_misc_coms  , only : simtime                & ! intent(in)
                            , current_time           & ! intent(in)
                            , dtlsm                  ! ! intent(in)
   use grid_coms     , only : nzg                    ! ! intent(in)
   use soil_coms     , only : slz                    & ! intent(in)
                            , soil                   & ! intent(in)
                            , dslz                   ! ! intent(in)
   use disturb_coms  , only : include_fire           & ! intent(in)
                            , fire_dryness_threshold & ! intent(in)
                            , fire_smoist_depth      & ! intent(in)
                            , k_fire_first           & ! intent(in)
                            , fire_parameter         ! ! intent(in)
   use allometry     , only : ed_biomass             ! ! function
   use consts_coms   , only : wdns                   & ! intent(in)
                            , day_sec                ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target     :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer    :: cpoly
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   type(simtime)                  :: lastmonth
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   integer                        :: k
   integer                        :: nsoil
   real                           :: ndaysi
   real                           :: normfac
   real                           :: fire_wmass_threshold
   real                           :: fuel
   real                           :: ignition_rate
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Find the number of days of last month so we can normalise the integrated ground   !
   ! water.                                                                                !
   !---------------------------------------------------------------------------------------!
   call lastmonthdate(current_time,lastmonth,ndaysi)
   normfac = dtlsm * ndaysi / (day_sec)
   !---------------------------------------------------------------------------------------!



   !----- Loop over polygons and sites. ---------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !----- Initialize ignition rate (a site variable). -------------------------------!
         ignition_rate = 0.0
         !---------------------------------------------------------------------------------!

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            !----- Normalise the monthly mean ground water. -------------------------------!
            csite%avg_monthly_gndwater(ipa) = csite%avg_monthly_gndwater(ipa) * normfac
            !------------------------------------------------------------------------------!

            !----- Initialize patch fuel. -------------------------------------------------!
            fuel = 0.0
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Loop through all cohorts in this patch, and compute the fuel.  Fuel will  !
            ! be defined as the above-ground biomass per unit area.                        !
            !------------------------------------------------------------------------------!
            cohortloop: do ico = 1,cpatch%ncohorts
               fuel = fuel + cpatch%nplant(ico) * cpatch%agb(ico)
            end do cohortloop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Determine the correct threshold to ignite fires according to the fire    !
            ! method.                                                                      !
            !------------------------------------------------------------------------------!
            select case (include_fire)
            case (0)
               !------ Set water mass threshold to infinity, so fires will never happen. --!
               fire_wmass_threshold = huge(1.)
               !---------------------------------------------------------------------------!

            case (1)
               !---------------------------------------------------------------------------!
               !     The fire threshold is equivalent to the dryness factor, converted to  !
               ! kg/m2.  This will be compared to the full column, so if the soil is too   !
               ! deep fires will be nearly impossible.                                     !
               !---------------------------------------------------------------------------!
               fire_wmass_threshold = fire_dryness_threshold * wdns
               !---------------------------------------------------------------------------!

            case (2)
               !---------------------------------------------------------------------------!
               !     We now compute the minimum amount of water in kg/m2 that the soil     !
               ! must have to avoid fires, using the soil properties and the soil moisture !
               ! fraction threshold.                                                       !
               !---------------------------------------------------------------------------!
               fire_wmass_threshold = 0
               do k = k_fire_first, nzg
                  nsoil                = cpoly%ntext_soil(k,isi)
                  fire_wmass_threshold = fire_wmass_threshold                              &
                                       + soil(nsoil)%soilfr * dslz(k) * wdns
               end do
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    If the soil is dry, then calculate patch contribution to the ignition     !
            ! rate.                                                                        !
            !------------------------------------------------------------------------------!
            if (csite%avg_monthly_gndwater(ipa) < fire_wmass_threshold) then
               ignition_rate = ignition_rate + fuel * csite%area(ipa)
            end if
            !------------------------------------------------------------------------------!


            !----- Reset the ground water for next month. ---------------------------------!
            csite%avg_monthly_gndwater(ipa) = 0.
            !------------------------------------------------------------------------------!

         end do patchloop
         !---------------------------------------------------------------------------------!



         !----- Calculate fire disturbance rate [1/month]. --------------------------------!
         cpoly%lambda_fire  (current_time%month,isi) = fire_parameter * ignition_rate
         cpoly%ignition_rate                   (isi) = ignition_rate
         !---------------------------------------------------------------------------------!
      end do siteloop
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine fire_frequency
!==========================================================================================!
!==========================================================================================!
