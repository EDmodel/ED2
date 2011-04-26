!==========================================================================================!
!==========================================================================================!
!    This subroutine will evaluate whether fire conditions exist, and if that is the case, !
! it will calculate the disturbance rate due to fire.                                      !
!------------------------------------------------------------------------------------------!
subroutine fire_frequency(month, cgrid)
   use ed_state_vars , only : edtype                 & ! structure
                            , polygontype            & ! structure
                            , sitetype               & ! structure
                            , patchtype              ! ! structure
   use grid_coms     , only : nzg                    ! ! intent(in)
   use soil_coms     , only : slz                    & ! intent(in)
                            , soil                   & ! intent(in)
                            , dslz                   ! ! intent(in)
   use disturb_coms  , only : include_fire           & ! intent(in)
                            , fire_dryness_threshold & ! intent(in)
                            , fire_smoist_threshold  & ! intent(in)
                            , fire_smoist_depth      & ! intent(in)
                            , k_fire_first           & ! intent(in)
                            , fire_parameter         ! ! intent(in)
   use allometry     , only : ed_biomass             ! ! function
   use consts_coms   , only : wdnsi                  & ! intent(in)
                            , wdns                   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target     :: cgrid
   integer           , intent(in) :: month
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer    :: cpoly
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   integer                        :: k
   integer                        :: ka
   integer                        :: nsoil
   real                           :: babove
   real                           :: fire_wmass_threshold
   real                           :: fuel
   real                           :: ignition_rate
   real                           :: patch_water_depth
   real                           :: patch_water_mass
   !----- Locally saved variables. --------------------------------------------------------!
   logical           , save       :: first_time=.true.
   !---------------------------------------------------------------------------------------!

   if (first_time) then
      kfireloop: do k_fire_first=nzg-1,1,-1
         if (slz(k_fire_first) < fire_smoist_depth) exit kfireloop
      end do kfireloop
      k_fire_first = k_fire_first + 1
      first_time = .false.
   end if


   !----- Loop over polygons and sites. ---------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !----- Initialize ignition rate (a site variable). -------------------------------!
         ignition_rate = 0.0
         
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            !----- Initialize patch fuel. -------------------------------------------------!
            fuel = 0.0
            
            !------------------------------------------------------------------------------!
            !    Loop through all cohorts in this patch, and compute the fuel.  Fuel will  !
            ! be defined as the above-ground biomass per unit area.                        !
            !------------------------------------------------------------------------------!
            cohortloop: do ico = 1,cpatch%ncohorts
               babove = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico),cpatch%bleaf(ico)  &
                                  ,cpatch%pft(ico),cpatch%hite(ico),cpatch%bstorage(ico)   &
                                  ,cpatch%bsapwood(ico))                                   &
                      * cpatch%nplant(ico)
               fuel   = fuel + babove
            end do cohortloop

            select case (include_fire)
            case (1)
               !---------------------------------------------------------------------------!
               !     Calculate the patch-level equivalent depth of ground liquid water.    !
               ! This is done by combining the total water mass from both the temporary    !
               ! surface water/snow layers plus the soil moisture.                         !
               !     SFCWATER_DEPTH is the depth of the temporary surface water/snow       !
               ! layer, and the depth calculated here will be the same as SFCWATER_DEPTH   !
               ! is the layer is liquid, but this will not be true if the ground is        !
               ! covered with snow.  Therefore we compute the depth based on the liquid    !
               ! water density [kg/m3] and surface water mass [kg/m2].                     !
               !     Similarly, the underground water is converted to meters using the     !
               ! soil moisture [m3/m3] and the depth of the soil layer [m].                !
               !---------------------------------------------------------------------------!
               patch_water_depth = 0.0
               do k = 1, csite%nlev_sfcwater(ipa)
                  patch_water_depth = patch_water_depth                                    &
                                    + csite%sfcwater_mass(k,ipa) * wdnsi
               end do
               do k = cpoly%lsl(isi), nzg
                  patch_water_depth = patch_water_depth                                    &
                                    + csite%soil_water(k,ipa)    * dslz(k)
               end do

               !---------------------------------------------------------------------------!
               !    If the soil is dry, then calculate patch contribution to the ignition  !
               ! rate.                                                                     !
               !---------------------------------------------------------------------------!
               if (patch_water_depth < fire_dryness_threshold) then
                  ignition_rate = ignition_rate + fuel * csite%area(ipa)
               end if

            case (2)
               !---------------------------------------------------------------------------!
               !    Compute the total (ground + underground) water in kg/m2.               !
               !---------------------------------------------------------------------------!
               patch_water_mass = 0.0
               ka = max(cpoly%lsl(isi),k_fire_first)

               do k = 1, csite%nlev_sfcwater(ipa)
                  patch_water_mass = patch_water_mass + csite%sfcwater_mass(k,ipa)
               end do
               do k = ka, nzg
                  patch_water_mass = patch_water_mass                                      &
                                   + csite%soil_water(k,ipa) * dslz(k) * wdns
               end do

               !---------------------------------------------------------------------------!
               !     We now compute the minimum amount of water in kg/m2 that the soil     !
               ! must have to avoid fires, using the soil properties and the soil moisture !
               ! fraction threshold.                                                       !
               !---------------------------------------------------------------------------!
               fire_wmass_threshold = 0
               do k = ka, nzg
                  nsoil                = cpoly%ntext_soil(k,isi)
                  fire_wmass_threshold = fire_wmass_threshold                              &
                                       + ( fire_smoist_threshold                           &
                                         * (soil(nsoil)%slmsts - soil(nsoil)%soilcp)       &
                                         + soil(nsoil)%soilcp ) * dslz(k) * wdns
               end do

               !---------------------------------------------------------------------------!
               !    If the soil is dry, then calculate patch contribution to the ignition  !
               ! rate.                                                                     !
               !---------------------------------------------------------------------------!
               if (patch_water_mass < fire_wmass_threshold) then
                  ignition_rate = ignition_rate + fuel * csite%area(ipa)
               end if

            end select
         end do patchloop
         
         !----- Calculate fire dist rate [1/month]. ---------------------------------------!
         cpoly%lambda_fire(month,isi) = fire_parameter * ignition_rate
         cpoly%ignition_rate(isi)     = ignition_rate

      end do siteloop
   end do polyloop

   return
end subroutine fire_frequency
!==========================================================================================!
!==========================================================================================!
