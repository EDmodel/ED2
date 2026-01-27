!==========================================================================================!
!==========================================================================================!
! MODULE FIRE
!
!> \brief This module contains routines to obtain fire disturbance rates
!> \details These subroutines are intended to calculate fire intensity and burned area
!!          which are or can be used to obtain fire disturbance rate and survivorship
!> \author  Paul Moorcroft, converted to fortran by David Medvigy
!> \author  10 Jan 2018.  MLO converted it into module so the code compiles with ifort 17.
!!          Also implementing process-based model, step by step.
!------------------------------------------------------------------------------------------!
module fire

   contains 

   !=======================================================================================!
   !=======================================================================================!
   ! SUB-ROUTINE FIRE_FREQUENCY
   !> This subroutine will evaluate whether fire conditions exist, and if that is the
   !! case, it will calculate the disturbance rate due to fire.
   !---------------------------------------------------------------------------------------!
   subroutine fire_frequency(cgrid)
      use ed_state_vars , only : edtype                 & ! structure
                               , polygontype            & ! structure
                               , sitetype               & ! structure
                               , patchtype              ! ! structure
      use ed_misc_coms  , only : simtime                & ! intent(in)
                               , current_time           & ! intent(in)
                               , dtlsm                  ! ! intent(in)
      use grid_coms     , only : nzg                    ! ! intent(in)
      use soil_coms     , only : soil                   & ! intent(in)
                               , dslz                   ! ! intent(in)
      use pft_coms      , only : is_grass               ! ! intent(in)
      use disturb_coms  , only : include_fire           & ! intent(in)
                               , fire_dryness_threshold & ! intent(in)
                               , k_fire_first           & ! intent(in)
                               , fire_parameter         & ! intent(in)
                               , fuel_height_max        ! ! intent(in)
      use consts_coms   , only : wdns                   & ! intent(in)
                               , wdnsi                  & ! intent(in)
                               , day_sec                ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      type(simtime)                  :: lastmonth
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      integer                        :: imon
      integer                        :: ipft
      integer                        :: k
      integer                        :: nsoil
      real                           :: ndaysi
      real                           :: normfac
      real                           :: fire_wmass_threshold
      real                           :: fire_intensity
      real                           :: fuel
      real                           :: ignition_rate
      real                           :: mean_fire_intensity
      real                           :: sum_accp
      logical                        :: people_around
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the number of days of last month so we can normalise the integrated       !
      ! ground water.                                                                      !
      !------------------------------------------------------------------------------------!
      call lastmonthdate(current_time,lastmonth,ndaysi)
      normfac = dtlsm * ndaysi / (day_sec)
      !------------------------------------------------------------------------------------!


      !----- Current month. ---------------------------------------------------------------!
      imon = current_time%month
      !------------------------------------------------------------------------------------!


      !----- Loop over polygons and sites. ------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !---------------------------------------------------------------------------------!
         !     Loop over all sites.                                                        !
         !---------------------------------------------------------------------------------!
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !     Find the total rainfall of the past year and reset the counter for this  !
            ! month.                                                                       !
            !------------------------------------------------------------------------------!
            sum_accp                         = sum(cpoly%avg_monthly_accp(:,isi))
            cpoly%avg_monthly_accp(imon,isi) = 0.
            !------------------------------------------------------------------------------!



            !----- Initialize ignition rate and the mean fire intensity (site variables). -!
            ignition_rate       = 0.0
            mean_fire_intensity = 0.0
            !------------------------------------------------------------------------------!



            !----- Temporary check for human activities. ----------------------------------!
            people_around = .false.
            humanloop: do ipa=1,csite%npatches
               select case (csite%dist_type(ipa))
               case (3)
                  continue
               case default
                  people_around = .true.
                  exit humanloop
               end select
            end do humanloop
            !----- Allow fires to ignite in intact forests. -------------------------------!
            ! people_around = .true.
            !------------------------------------------------------------------------------!


            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Normalise the monthly mean ground water. ----------------------------!
               csite%avg_monthly_gndwater(ipa) = csite%avg_monthly_gndwater(ipa) * normfac
               !---------------------------------------------------------------------------!


               !----- Normalise the monthly mean ground water. ----------------------------!
               csite%avg_monthly_waterdef(ipa) = max(0.0,csite%avg_monthly_waterdef(ipa))
               !------------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Obtain fuel stocks.  The original fire model would consider all       !
               ! above-ground biomass and no litter.  When include_fire is set to 3, fuel  !
               ! is the sum of all above-ground fast soil C, and biomass from grasses and  !
               ! small individuals (up to fuel_max_height).                                !
               !---------------------------------------------------------------------------!
               select case (include_fire)
               case (3)
                  fuel = csite%fast_grnd_C(ipa) + csite%structural_grnd_C(ipa)
                  fuelcohloop_3: do ico = 1,cpatch%ncohorts
                     ipft = cpatch%pft(ico)
                     if (is_grass(ipft) .or. cpatch%height(ico) <= fuel_height_max) then
                        fuel = fuel + cpatch%nplant(ico) * cpatch%agb(ico)
                     end if
                  end do fuelcohloop_3
               case default
                  fuel = 0.0
                  fuelcohloop_d: do ico = 1,cpatch%ncohorts
                     fuel = fuel + cpatch%nplant(ico) * cpatch%agb(ico)
                  end do fuelcohloop_d
               end select
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Determine the correct threshold to ignite fires according to the fire !
               ! method.                                                                   !
               !---------------------------------------------------------------------------!
               select case (include_fire)
               case (0)
                  !------ Set water mass threshold to infinity, so fires can't happen. ----!
                  fire_wmass_threshold = huge(1.)
                  fire_intensity       = 0.
                  !------------------------------------------------------------------------!

               case (1)
                  !------------------------------------------------------------------------!
                  !     The fire threshold is equivalent to the dryness factor, converted  !
                  ! to kg/m2.  This will be compared to the full column, so if the soil is !
                  ! too deep then fires would be nearly impossible.                        !
                  !------------------------------------------------------------------------!
                  fire_wmass_threshold = fire_dryness_threshold * wdns
                  if (csite%avg_monthly_gndwater(ipa) < fire_wmass_threshold) then
                     fire_intensity      = fire_parameter
                     mean_fire_intensity = mean_fire_intensity                             &
                                         + fire_intensity * csite%area(ipa)
                  else
                     fire_intensity      = 0.0
                  end if
                  !------------------------------------------------------------------------!

               case (2)
                  !------------------------------------------------------------------------!
                  !     We now compute the minimum amount of water in kg/m2 that the soil  !
                  ! must have to avoid fires, using the soil properties and the soil       !
                  ! moisture fraction threshold.                                           !
                  !------------------------------------------------------------------------!
                  fire_wmass_threshold = 0.
                  do k = k_fire_first, nzg
                     nsoil                = cpoly%ntext_soil(k,isi)
                     fire_wmass_threshold = fire_wmass_threshold                           &
                                          + soil(nsoil)%soilfr * dslz(k) * wdns
                  end do
                  if (csite%avg_monthly_gndwater(ipa) < fire_wmass_threshold) then
                     fire_intensity      = fire_parameter
                     mean_fire_intensity = mean_fire_intensity                             &
                                         + fire_intensity * csite%area(ipa)
                  else
                     fire_intensity      = 0.0
                  end if
                  !------------------------------------------------------------------------!

               case (3)
                  !------------------------------------------------------------------------!
                  !     Set fire intensity the same as method 2, except that we prevent    !
                  ! fires until the site has any anthropogenic disturbance.  This will     !
                  ! change in the future to allow natural fires.                           !
                  !------------------------------------------------------------------------!
                  if (people_around) then
                     fire_wmass_threshold = 0.
                     do k = k_fire_first, nzg
                        nsoil                = cpoly%ntext_soil(k,isi)
                        fire_wmass_threshold = fire_wmass_threshold                        &
                                             + soil(nsoil)%soilfr * dslz(k) * wdns
                     end do
                     if (csite%avg_monthly_gndwater(ipa) < fire_wmass_threshold) then
                        fire_intensity      = fire_parameter
                        mean_fire_intensity = mean_fire_intensity                          &
                                            + fire_intensity * csite%area(ipa)
                     else
                        fire_intensity      = 0.0
                     end if
                  else
                     !------ Set water mass threshold to infinity, so fires won't happen. -!
                     fire_wmass_threshold = huge(1.)
                     fire_intensity       = 0.
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end select
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    If the soil is dry, then calculate patch contribution to the ignition  !
               ! rate.                                                                     !
               !---------------------------------------------------------------------------!
               ignition_rate = ignition_rate + fire_intensity * fuel * csite%area(ipa)
               !---------------------------------------------------------------------------!


               !----- Reset the ground water for next month. ------------------------------!
               csite%avg_monthly_gndwater(ipa) = 0.
               !---------------------------------------------------------------------------!

            end do patchloop
            !------------------------------------------------------------------------------!



            !----- Calculate fire disturbance rate [1/month]. -----------------------------!
            cpoly%lambda_fire  (imon,isi) = ignition_rate
            if (mean_fire_intensity > 0.) then
               cpoly%ignition_rate (isi) = ignition_rate / mean_fire_intensity
            else
               cpoly%ignition_rate (isi) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine fire_frequency
   !=======================================================================================!
   !=======================================================================================!
end module fire
!==========================================================================================!
!==========================================================================================!

