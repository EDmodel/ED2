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
                            , dslz                   & ! intent(in)
                            , dslzi                  ! ! intent(in)
   use disturb_coms  , only : include_fire           & ! intent(in)
                            , sm_fire                & ! intent(in)
                            , fire_dryness_threshold & ! intent(in)
                            , fire_smoist_depth      & ! intent(in)
                            , k_fire_first           & ! intent(in)
                            , fire_parameter         ! ! intent(in)
   use consts_coms   , only : wdns                   & ! intent(in)
                            , wdnsi                  & ! intent(in)
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
   integer                        :: imon
   integer                        :: k
   integer                        :: nsoil
   real                           :: ndaysi
   real                           :: normfac
   real                           :: fire_wmass_threshold
   real                           :: fire_intensity
   real                           :: fuel
   real                           :: avg_slmst
   real                           :: avg_slpot
   real                           :: ignition_rate
   real                           :: fire_scale
   real                           :: mean_fire_intensity
   real                           :: sum_pcpg
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the number of days of last month so we can normalise the integrated ground   !
   ! water.                                                                                !
   !---------------------------------------------------------------------------------------!
   call lastmonthdate(current_time,lastmonth,ndaysi)
   normfac = dtlsm * ndaysi / (day_sec)
   !---------------------------------------------------------------------------------------!


   !----- Current month. ------------------------------------------------------------------!
   imon = current_time%month
   !---------------------------------------------------------------------------------------!


   !----- Loop over polygons and sites. ---------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)


      !------------------------------------------------------------------------------------!
      !     Loop over all sites.                                                           !
      !------------------------------------------------------------------------------------!
      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !     Find the total rainfall of the past year and reset the counter for this     !
         ! month.                                                                          !
         !---------------------------------------------------------------------------------!
         sum_pcpg                         = sum(cpoly%avg_monthly_pcpg(:,isi))
         cpoly%avg_monthly_pcpg(imon,isi) = 0.
         !---------------------------------------------------------------------------------!



         !----- Initialize ignition rate and the mean fire intensity (site variables). ----!
         ignition_rate       = 0.0
         mean_fire_intensity = 0.0
         !---------------------------------------------------------------------------------!

         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            !----- Normalise the monthly mean ground water. -------------------------------!
            csite%avg_monthly_gndwater(ipa) = csite%avg_monthly_gndwater(ipa) * normfac
            !------------------------------------------------------------------------------!


            !----- Normalise the monthly mean ground water. -------------------------------!
            csite%avg_monthly_waterdef(ipa) = max(0.0,csite%avg_monthly_waterdef(ipa))
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
               fire_intensity       = 0.
               !---------------------------------------------------------------------------!

            case (1)
               !---------------------------------------------------------------------------!
               !     The fire threshold is equivalent to the dryness factor, converted to  !
               ! kg/m2.  This will be compared to the full column, so if the soil is too   !
               ! deep fires will be nearly impossible.                                     !
               !---------------------------------------------------------------------------!
               fire_wmass_threshold = fire_dryness_threshold * wdns
               if (csite%avg_monthly_gndwater(ipa) < fire_wmass_threshold) then
                  fire_intensity      = fire_parameter
                  mean_fire_intensity = mean_fire_intensity                                &
                                      + fire_intensity * csite%area(ipa)
               else
                  fire_intensity      = 0.0
               end if
               !---------------------------------------------------------------------------!

            case (2)
               !---------------------------------------------------------------------------!
               !     We now compute the minimum amount of water in kg/m2 that the soil     !
               ! must have to avoid fires, using the soil properties and the soil moisture !
               ! fraction threshold.                                                       !
               !---------------------------------------------------------------------------!
               fire_wmass_threshold = 0.
               do k = k_fire_first, nzg
                  nsoil                = cpoly%ntext_soil(k,isi)
                  fire_wmass_threshold = fire_wmass_threshold                              &
                                       + soil(nsoil)%soilfr * dslz(k) * wdns
               end do
               if (csite%avg_monthly_gndwater(ipa) < fire_wmass_threshold) then
                  fire_intensity      = fire_parameter
                  mean_fire_intensity = mean_fire_intensity                                &
                                      + fire_intensity * csite%area(ipa)
               else
                  fire_intensity      = 0.0
               end if
               !---------------------------------------------------------------------------!

            case (3)
               !---------------------------------------------------------------------------!
               !     The threshold is independent on soil moisture.  We use climatological !
               ! water deficit instead.                                                    !
               !---------------------------------------------------------------------------!
               if (csite%avg_monthly_waterdef(ipa) >= sm_fire * sum_pcpg) then
                  fire_intensity      = fire_parameter
                  mean_fire_intensity = mean_fire_intensity                                &
                                      + fire_intensity * csite%area(ipa)
               else
                  fire_intensity      = 0.0
               end if
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    If the soil is dry, then calculate patch contribution to the ignition     !
            ! rate.                                                                        !
            !------------------------------------------------------------------------------!
            ignition_rate = ignition_rate + fire_intensity * fuel * csite%area(ipa)
            !------------------------------------------------------------------------------!


            !----- Reset the ground water for next month. ---------------------------------!
            csite%avg_monthly_gndwater(ipa) = 0.
            !------------------------------------------------------------------------------!

         end do patchloop
         !---------------------------------------------------------------------------------!



         !----- Calculate fire disturbance rate [1/month]. --------------------------------!
         cpoly%lambda_fire  (imon,isi) = ignition_rate
         if (mean_fire_intensity > 0.) then
            cpoly%ignition_rate (isi) = ignition_rate / mean_fire_intensity
         else
            cpoly%ignition_rate (isi) = 0.0
         end if
         !---------------------------------------------------------------------------------!
      end do siteloop
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine fire_frequency
!==========================================================================================!
!==========================================================================================!
