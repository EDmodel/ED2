!==========================================================================================!
!==========================================================================================!
!    This module contains some subroutines used to derive turbulence-related properties of !
! the canopy air space.                                                                    !
!------------------------------------------------------------------------------------------!
module canopy_struct_dynamics

   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     The following subroutine calculates the turbulent transfer of scalar quantities   !
   ! using Monin Obukhov Similarity Theory under a few different flavors.  The basic       !
   ! assumptions are that the ground surface and canopy space is modeled as a rough wall   !
   ! boundary above which exists a dynamic sublayer, above which exists the free atmo-     !
   ! sphere.  We assume a logarithmic decay of wind speed over the dynamic sublayer, which !
   ! is due to the assumption of constant shear (ustar) in the dynamic sublayer.           !
   !                                                                                       !
   ! 0.  We apply the wind profile based on the similarity theory for the tallest cohort,  !
   !     with further reduction due to the presence of obstacles (trees), similar to       !
   !     Leuning et al. (1995) but considering a finite crown area.  The aerodynamic       !
   !     resistance from ground to canopy air is found in a similar way as in LEAF-3.      !
   !     Displacement height and roughness are linear functions of the vegetation height,  !
   !     where vegetation height is the weighted average of cohort heights, the weights    !
   !     being the basal area.                                                             !
   !     References:                                                                       !
   !                                                                                       !
   !     Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf          !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   !                                                                                       !
   ! 1. This is similar to option 0, where displacement height and rougness are linear     !
   !    functions of mean canopy height, where mean canopy height is defined by the        !
   !    dead-biomass weighted average of canopy height.  Wind speed and diffusivity decay  !
   !    exponentially in the canopy on parameter alpha.  HOWEVER, if the windspeed at      !
   !    reference height is less than the canopy height (WHICH IS REALLY INNAPROPRIATE TO  !
   !    BEGIN WITH, BUT OH WELL), then revert to a calculation of ustar that uses no       !
   !    displacement height.  After ustar is calculated, then use displacement height in   !
   !    the calculation of the eddy length scale at the canopy top to find diffusivity at  !
   !    canopy top. THIS IS THE METHOD USED IN LEAF3 AND ED2.0.                            !
   !                                                                                       !
   ! 2. Follow the methods of Massman 1997.  The canopy top is taken as the height of the  !
   !    tallest cohort in the patch. The displacement height is determined by a center of  !
   !    pressure method. The wind speed inside the canopy decays exponentially on          !
   !    cumulative foliar drag elements, where the equation is scaled such that it solves  !
   !    for the balance of integrated canopy shear with dissipation.  To close the         !
   !    equations M97 assumes that the surface drag (U(h)/Ustar) is a function of the      !
   !    cumulative foliar drag surfaces.  With this, we can explicitly solve for roughness !
   !    length. Assumptions are made about the profile of canopy drag as a function of     !
   !    cohort attributes, please see the technical manual by RGK2009 for more info.       !
   !    Please also refer to the original paper by Massman for information regarding the   !
   !    basic principles of the closure scheme.                                            !
   !                                                                                       !
   !    Massman, WJ. An Analytical One-Dimensional Model of Momentum Transfer By           !
   !      Vegetation Of Arbitrary Structure. Boundary Layer Meteorology, 83,               !
   !      p. 407-421, 1997.                                                                !
   !                                                                                       !
   ! 3.  Same as 2, but the drag and sheltering parameters are functions of height.        !
   !                                                                                       !
   !     Ultimately, this routine solves for the resistance of water vapor and sensible    !
   ! heat from the soil surface to canopy air space and from the leaf surfaces to canopy   !
   ! air space.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_turbulence(cpoly,isi,ipa)
      use ed_state_vars  , only : polygontype          & ! structure
                                , sitetype             & ! structure
                                , patchtype            ! ! structure
      use met_driver_coms, only : met_driv_state       ! ! structure
      use rk4_coms       , only : ibranch_thermo       ! ! intent(in)
      use canopy_air_coms, only : icanturb             & ! intent(in), can. turb. scheme
                                , ustmin               & ! intent(in)
                                , ugbmin               & ! intent(in)
                                , ubmin                & ! intent(in)
                                , gamh                 & ! intent(in)
                                , exar                 & ! intent(in)
                                , dz_m97               & ! intent(in)
                                , cdrag0               & ! intent(in)
                                , pm0                  & ! intent(in)
                                , c1_m97               & ! intent(in)
                                , c2_m97               & ! intent(in)
                                , c3_m97               & ! intent(in)
                                , kvwake               & ! intent(in)
                                , alpha1_m97           & ! intent(in)
                                , alpha2_m97           & ! intent(in)
                                , psi_m97              & ! intent(in)
                                , zztop                & ! intent(in)
                                , zzmid                & ! intent(in)
                                , lad                  & ! intent(in)
                                , dladdz               & ! intent(in)
                                , cdrag                & ! intent(in)
                                , pshelter             & ! intent(in)
                                , cumldrag             & ! intent(in)
                                , windm97              & ! intent(in)
                                , rb_inter             & ! intent(in)
                                , rb_slope             & ! intent(in)
                                , ggfact               & ! intent(in)
                                , zoobukhov            ! ! function
      use consts_coms    , only : vonk                 & ! intent(in)
                                , cp                   & ! intent(in)
                                , cpi                  & ! intent(in)
                                , grav                 & ! intent(in)
                                , epim1                & ! intent(in)
                                , sqrt2o2              ! ! intent(in)
      use soil_coms      , only : snow_rough           & ! intent(in)
                                , soil_rough           ! ! intent(in)
      use allometry      , only : h2trunkh             & ! function
                                , dbh2bl               ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)   , target      :: cpoly         ! Current polygon
      integer             , intent(in)  :: isi           ! Site index
      integer             , intent(in)  :: ipa           ! Patch index
      !----- Pointers ---------------------------------------------------------------------!
      type(sitetype)      , pointer     :: csite         ! Current site
      type(patchtype)     , pointer     :: cpatch        ! Current patch
      type(met_driv_state), pointer     :: cmet          ! Current meteorology
      !----- Local variables --------------------------------------------------------------!
      integer        :: ico          ! Cohort loop
      integer        :: ipft         ! PFT alias
      integer        :: k            ! Elevation index
      integer        :: zcan         ! Index of canopy top elevation
      integer        :: kafull       ! First layer fully occupied by crown
      integer        :: kzfull       ! Last layer fully occupied by crown
      integer        :: kapartial    ! First layer partially occupied by crown
      integer        :: kzpartial    ! Last layer partially occupied by crown
      logical        :: stable       ! Stable canopy air space
      real           :: rasveg       ! Resistance of vegetated ground           [      s/m]
      real           :: atm_thetav   ! Free atmosphere virtual potential temp.  [        K]
      real           :: can_thetav   ! Free atmosphere virtual potential temp.  [        K]
      real           :: ldga_bk      ! Cumulative leaf drag area                [      ---]
      real           :: lyrhalf      ! Half the contrib. of this layer to zeta  [      1/m]
      real           :: sigmakm      ! Km coefficient at z=h                    [        m]
      real           :: K_top        ! Diffusivity at canopy top z=h            [     m2/s]
      real           :: Kdiff        ! Diffusivity                              [     m2/s]
      real           :: surf_rough   ! Roughness length of the bare ground 
                                     !     at canopy bottom                     [        m]
      real           :: uh           ! Wind speed at the canopy top (z=h)       [      m/s]
      real           :: factv        ! Wind-dependent term for old rasveg
      real           :: aux          ! Aux. variable
      real           :: estar        ! Equivalent potential temperature         [        K]
      real           :: gbhmos_min   ! Minimum boundary layer heat conductance. [      m/s]
      real           :: wcapcan      ! Canopy air space water capacity          [    kg/m2]
      real           :: wcapcani     ! Inverse of the guy above                 [    m2/kg]
      real           :: hcapcani     ! Inverse of canopy air space heat cap.    [   m2.K/J]
      real           :: ccapcani     ! Inverse of canopy air space CO2 capacity [   m2/mol]
      real           :: ustarouh     ! The ratio of ustar over u(h)             [      ---]
      real           :: nn           ! In-canopy wind attenuation scal. param.  [      ---]
      real           :: waiuse       ! Wood area index                          [    m2/m2]
      real           :: htopcrown    ! height at the top of the crown           [        m]
      real           :: hmidcrown    ! Height at the middle of the crown        [        m]
      real           :: hbotcrown    ! Height at the bottom of the crown        [        m]
      real           :: htop         ! Height of the topmost layer              [        m]
      real           :: dzcrown      ! Depth that contains leaves/branches      [        m]
      real           :: d0ohgt       ! d0/height                                [      ---]
      real           :: z0ohgt       ! z0/height                                [      ---]
      real           :: phih         ! Correction term for unstable cases       [      ---]
      real           :: ladcohort    ! Leaf Area Density of this cohort         [    m2/m3]
      real           :: ribcan       ! Ground-to-canopy bulk Richardson number  [      ---]
      real           :: hgtoz0       ! height/z0                                [      ---]
      real           :: lnhgtoz0     ! log(height/z0)                           [      ---]
      real           :: zetacan      ! Estimate of z/L within the canopy        [      ---]
      !------------------------------------------------------------------------------------!

      !----- Assign some pointers. --------------------------------------------------------!
      csite  => cpoly%site(isi)
      cmet   => cpoly%met(isi)
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the virtual potential temperatures and decide whether the canopy air is   !
      ! stable or not.                                                                     !
      !------------------------------------------------------------------------------------!
      atm_thetav = cmet%atm_theta       * (1. + epim1 * cmet%atm_shv      )
      can_thetav = csite%can_theta(ipa) * (1. + epim1 * csite%can_shv(ipa))
      stable     = atm_thetav >= can_thetav



      !------------------------------------------------------------------------------------!
      !     If there is no vegetation in this patch, then we apply turbulence to bare      !
      ! soil, no d0 and exit.                                                              !
      !------------------------------------------------------------------------------------!
      if (cpatch%ncohorts == 0) then
         
         !----- Get the appropriate characteristic wind speed. ----------------------------!
         if (stable) then
            cmet%vels = cmet%vels_stab
         else
            cmet%vels = cmet%vels_unstab
         end if

         !----- Calculate the surface roughness inside the canopy. ------------------------!
         csite%rough(ipa) = soil_rough * (1.0 - csite%snowfac(ipa))                        &
                          + snow_rough * csite%snowfac(ipa)
         
         !----- Finding the characteristic scales (a.k.a. stars). -------------------------!
         call ed_stars(cmet%atm_theta,cmet%atm_theiv,cmet%atm_shv,cmet%atm_co2             &
                      ,csite%can_theta(ipa),csite%can_theiv(ipa),csite%can_shv(ipa)        &
                      ,csite%can_co2(ipa),cmet%geoht,csite%veg_displace(ipa),cmet%vels     &
                      ,csite%rough(ipa),csite%ustar(ipa),csite%tstar(ipa),estar            &
                      ,csite%qstar(ipa),csite%cstar(ipa),csite%zeta(ipa),csite%ribulk(ipa) &
                      ,csite%ggbare(ipa))

         !---------------------------------------------------------------------------------!
         !      This is a bare ground cohort, so there is no vegetated ground.  Assign     !
         ! zero to the conductance, but it shouldn't be used at all.                       !
         !---------------------------------------------------------------------------------!
         csite%ggveg(ipa) = 0.0
         csite%ggnet(ipa) = ggfact * csite%ggbare(ipa)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap(csite%can_rhos(ipa),csite%can_temp(ipa),csite%can_depth(ipa)       &
                       ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!
         return
      end if
      !------------------------------------------------------------------------------------!



      !----- Check whether the reference height is high enough . --------------------------!
      if (cpatch%hite(1) >= cmet%geoht) then
         write (unit=*,fmt='(a)') ' Your reference height is too low...'
         write (unit=*,fmt='(a,1x,es12.5)') ' Ref. height:           ',cmet%geoht
         write (unit=*,fmt='(a,1x,es12.5)') ' Tallest cohort height: ',cpatch%hite(1)
         call fatal_error('Bad reference height','canopy_turbulence'                       &
                         ,'canopy_struct_dynamics.f90')
      end if
      !------------------------------------------------------------------------------------!




      !---- Find the minimum leaf boundary layer heat conductance. ------------------------!
      if (any(cpatch%resolvable)) then
         gbhmos_min = 1. / (rb_inter + rb_slope * (csite%lai(ipa) + csite%wai(ipa)))
      else
         gbhmos_min = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we do have cohorts, choose which method we use to compute the          !
      ! resistance.                                                                        !
      !------------------------------------------------------------------------------------!
      select case (icanturb)

      !------------------------------------------------------------------------------------!
      ! LEAF-3 Case: This approach is very similar to older implementations of ED-2, and   !
      !              it is very similar to LEAF-3, and to option 1, except that option 1   !
      !              computes the vegetated ground conductance and wind profile different- !
      !              ly.                                                                   !
      !------------------------------------------------------------------------------------!
      case (0)

         !---------------------------------------------------------------------------------!
         !     Find the roughness as the average between the bare ground and vegetated     !
         ! ground, and apply the snow cover to further scale it.  The weighting factors    !
         ! are the fraction of open canopy and the fraction of the canopy buried in snow.  !
         !---------------------------------------------------------------------------------!
         csite%rough(ipa) = snow_rough * csite%snowfac(ipa)                                &
                          + ( soil_rough           * csite%opencan_frac(ipa)               &
                            + csite%veg_rough(ipa) * (1.0 - csite%opencan_frac(ipa)) )     &
                          * (1.0 - csite%snowfac(ipa))
         !---------------------------------------------------------------------------------!



         !----- Get the appropriate characteristic wind speed. ----------------------------!
         if (stable) then
            cmet%vels = cmet%vels_stab
         else
            cmet%vels = cmet%vels_unstab
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!
         call ed_stars(cmet%atm_theta,cmet%atm_theiv,cmet%atm_shv,cmet%atm_co2             &
                      ,csite%can_theta(ipa),csite%can_theiv(ipa),csite%can_shv(ipa)        &
                      ,csite%can_co2(ipa),cmet%geoht,csite%veg_displace(ipa),cmet%vels     &
                      ,csite%rough(ipa),csite%ustar(ipa),csite%tstar(ipa),estar            &
                      ,csite%qstar(ipa),csite%cstar(ipa),csite%zeta(ipa),csite%ribulk(ipa) &
                      ,csite%ggbare(ipa))

         if (csite%snowfac(ipa) < 0.9) then
            factv  = log((cmet%geoht - csite%veg_displace(ipa)) / csite%rough(ipa))        &
                   / (vonk * vonk * cmet%vels)
            aux    = exp(exar * (1. - (csite%veg_displace(ipa) + csite%rough(ipa))         &
                                    / csite%veg_height(ipa)) )
            csite%ggveg(ipa) = (exar * (csite%veg_height(ipa) - csite%veg_displace(ipa) )) &
                             / (factv * csite%veg_height(ipa) * (exp(exar) - aux))
         else 
            csite%ggveg(ipa) = 0.
         end if


         !---------------------------------------------------------------------------------!
         !     This part of the code initializes the geometry of the canopy air space, the !
         ! structure of the vegetation and its attenuation effects and the heat and water  !
         ! capacities.                                                                     !
         !---------------------------------------------------------------------------------!
         uh = reduced_wind(csite%ustar(ipa),csite%zeta(ipa),csite%ribulk(ipa),cmet%geoht   &
                          ,csite%veg_displace(ipa),cpatch%hite(1),csite%rough(ipa))

         do ico=1,cpatch%ncohorts

            !----- Calculate the wind speed at height z. ----------------------------------!
            ipft       = cpatch%pft(ico)

            cpatch%veg_wind(ico) = max(uh,ugbmin)
            uh = uh * ( cpatch%crown_area(ico)                                             &
                      * exp(- 0.5 * cpatch%lai(ico) / cpatch%crown_area(ico))              &
                      + 1.0 - cpatch%crown_area(ico)) 

            if (cpatch%resolvable(ico)) then

               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call aerodynamic_conductances(ipft,cpatch%veg_wind(ico)                     &
                                            ,cpatch%veg_temp(ico),csite%can_temp(ipa)      &
                                            ,csite%can_shv(ipa),csite%can_rhos(ipa)        &
                                            ,gbhmos_min,cpatch%gbh(ico),cpatch%gbw(ico))
               !---------------------------------------------------------------------------!
            else
               cpatch%gbh(ico)      = 0.0
               cpatch%gbw(ico)      = 0.0
            end if
         end do

         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap(csite%can_rhos(ipa),csite%can_temp(ipa),csite%can_depth(ipa)       &
                       ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (csite%opencan_frac(ipa) > 0.999 .or. csite%snowfac(ipa) >= 0.9) then
            csite%ggnet(ipa) = ggfact * csite%ggbare(ipa)
         else
            csite%ggnet(ipa) = ggfact * csite%ggbare(ipa) * csite%ggveg(ipa)               &
                             / ( csite%ggveg(ipa) + (1. - csite%opencan_frac(ipa))         &
                                                  * csite%ggbare(ipa))
         end if
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      ! ED-2.1 Case: There uses a slightly different approach for wind profile and ground  !
      !              conductance.                                                          !
      !------------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     The roughness is found by combining two weighted averages.  The first one   !
         ! checks the fraction of the patch that has closed canopy, and averages between   !
         ! soil and vegetation roughness.  The other is the fraction of the vegetation     !
         ! that is covered in snow.                                                        !
         !---------------------------------------------------------------------------------!
         csite%rough(ipa) = snow_rough * csite%snowfac(ipa)                                &
                          + ( soil_rough           * csite%opencan_frac(ipa)               &
                            + csite%veg_rough(ipa) * (1.0 - csite%opencan_frac(ipa)) )     &
                          * (1.0 - csite%snowfac(ipa))
         !---------------------------------------------------------------------------------!



         !----- Calculate the soil surface roughness inside the canopy. -------------------!
         surf_rough = soil_rough * (1.0 - csite%snowfac(ipa))                              &
                    + snow_rough * csite%snowfac(ipa)
         !---------------------------------------------------------------------------------!



         !----- Get the appropriate characteristic wind speed. ----------------------------!
         if (stable) then
            cmet%vels = cmet%vels_stab
         else
            cmet%vels = cmet%vels_unstab
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!
         call ed_stars(cmet%atm_theta,cmet%atm_theiv,cmet%atm_shv,cmet%atm_co2             &
                      ,csite%can_theta(ipa),csite%can_theiv(ipa),csite%can_shv(ipa)        &
                      ,csite%can_co2(ipa),cmet%geoht,csite%veg_displace(ipa),cmet%vels     &
                      ,csite%rough(ipa),csite%ustar(ipa),csite%tstar(ipa),estar            &
                      ,csite%qstar(ipa),csite%cstar(ipa),csite%zeta(ipa),csite%ribulk(ipa) &
                      ,csite%ggbare(ipa))
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     The surface resistance in the sub-canopy layer is the integration of        !
         ! inverse K, from the rough soil surface, to the reference point in the canopy    !
         ! where the state variable is integrated (canopy top ~ h).                        !
         !---------------------------------------------------------------------------------!
         K_top            = vonk * csite%ustar(ipa)                                        &
                          * (csite%veg_height(ipa) - csite%veg_displace(ipa))
         csite%ggveg(ipa) = ( exar * K_top )                                               &
                          / ( exp(exar) * csite%veg_height(ipa)                            &
                            * (exp(-exar * surf_rough/csite%veg_height(ipa)) - exp(-exar)))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     This part of the code finds the geometry of the canopy air space, the       !
         ! structure of the vegetation and its attenuation effects and the heat and water  !
         ! capacities.                                                                     !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh   = reduced_wind(csite%ustar(ipa),csite%zeta(ipa),csite%ribulk(ipa),cmet%geoht &
                            ,csite%veg_displace(ipa),cpatch%hite(1),csite%rough(ipa))
         htop = cpatch%hite(1)
         do ico=1,cpatch%ncohorts

            !----- Alias for PFT type. ----------------------------------------------------!
            ipft  = cpatch%pft(ico)
            !------------------------------------------------------------------------------!



            !----- Estimate the height at the crown midpoint. -----------------------------!
            htopcrown = cpatch%hite(ico)
            hbotcrown = h2trunkh(cpatch%hite(ico),ipft)
            hmidcrown = 0.5 * (htopcrown + hbotcrown)
            !------------------------------------------------------------------------------!



            !----- Calculate the wind speed at height z. ----------------------------------!
            cpatch%veg_wind(ico) = max( ugbmin                                             &
                                      , uh * exp(-exar * (1.0 - hmidcrown/htop)))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find the aerodynamic conductances for heat and water at the leaf boundary !
            ! layer.                                                                       !
            !------------------------------------------------------------------------------!
            if (cpatch%resolvable(ico)) then
               call aerodynamic_conductances(ipft,cpatch%veg_wind(ico)                     &
                                            ,cpatch%veg_temp(ico)                          &
                                            ,csite%can_temp(ipa),csite%can_shv(ipa)        &
                                            ,csite%can_rhos(ipa),gbhmos_min                &
                                            ,cpatch%gbh(ico),cpatch%gbw(ico))
            else
               cpatch%gbh(ico) = 0.0
               cpatch%gbw(ico) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap(csite%can_rhos(ipa),csite%can_temp(ipa),csite%can_depth(ipa)       &
                       ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (csite%opencan_frac(ipa) > 0.999 .or. csite%snowfac(ipa) >= 0.9) then
            csite%ggnet(ipa) = ggfact * csite%ggbare(ipa)
         else
            csite%ggnet(ipa) = ggfact * csite%ggbare(ipa) * csite%ggveg(ipa)               &
                             / ( csite%ggveg(ipa)  + (1. - csite%opencan_frac(ipa))        &
                                                   * csite%ggbare(ipa) )
         end if
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Use the methods of Massman (1997).  Using option 2 will make the within-      !
      ! canopy drag and sheltering factors to be constant, whereas option 3 will allow the !
      ! values to vary.                                                                    !
      !------------------------------------------------------------------------------------!
      case (2,3)
         !----- Get the appropriate characteristic wind speed. ----------------------------!
         if (stable) then
            cmet%vels = cmet%vels_stab
         else
            cmet%vels = cmet%vels_unstab
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the top layer and the top height.                                       !
         !---------------------------------------------------------------------------------!
         zcan = ceiling(cpatch%hite(1) / dz_m97)
         htop = zztop(zcan)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Loop through cohorts, and integrate the leaf area density.  Notice that we  !
         ! will be solving all cohorts here, because even those with no leaves have        !
         ! effects on the canopy air turbulence.  We only need to decide whether we have   !
         ! branches or not, and we do it outside the loop to make it more efficient.       !
         !---------------------------------------------------------------------------------!
         !----- Reset the leaf area density array. ----------------------------------------!
         lad(:) = 0.0
         select case (ibranch_thermo)
         case (0)
            !------------------------------------------------------------------------------!
            !     Although we are not solving branch thermodynamics, we have some drag     !
            ! effect due to the branches.  Assume branch area index to be 15% of on-       !
            ! -allometry LAI.  Even for grasses we may want to keep some residual due to   !
            ! dead material that remains around.                                           !
            !------------------------------------------------------------------------------!
            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)


               !------ Estimate WAI. ------------------------------------------------------!
               waiuse = 0.10 * cpatch%nplant(ico) * cpatch%sla(ico)                        &
                      * dbh2bl(cpatch%dbh(ico),ipft)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the heights, and compute the LAD of this cohort.                 !
               !---------------------------------------------------------------------------!
               htopcrown = cpatch%hite(ico)
               hbotcrown = h2trunkh(cpatch%hite(ico),ipft)
               ladcohort = (cpatch%lai(ico) + waiuse) / (htopcrown - hbotcrown)
               kapartial = floor  (hbotcrown/dz_m97) + 1
               kafull    = ceiling(hbotcrown/dz_m97) + 1
               kzpartial = ceiling(htopcrown/dz_m97)
               kzfull    = floor  (htopcrown/dz_m97)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Add the LAD for the full layers.                                      !
               !---------------------------------------------------------------------------!
               do k = kafull,kzfull
                  lad(k) = lad(k) + ladcohort
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the LAD for the partial layers.  The only special case is when   !
               ! they are both the same layer, which must be done separately.              !
               !---------------------------------------------------------------------------!
               if (kapartial == kzpartial) then
                  lad(kapartial) = lad(kapartial)                                          &
                                 + ladcohort * (htopcrown        - hbotcrown         )     &
                                             / (zztop(kapartial) - zztop(kapartial-1))
               else
                  lad(kapartial) = lad(kapartial)                                          &
                                 + ladcohort * (zztop(kapartial) - hbotcrown         )     &
                                             / (zztop(kapartial) - zztop(kapartial-1))
                  lad(kzpartial) = lad(kzpartial)                                          &
                                 + ladcohort * (htopcrown        - zztop(kzpartial-1))     &
                                             / (zztop(kzpartial) - zztop(kzpartial-1))
               end if
               !---------------------------------------------------------------------------!
            end do

         case default
            !----- Use the default wood area index. ---------------------------------------!
            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)

               !---------------------------------------------------------------------------!
               !     Find the heights, and compute the LAD of this cohort.                 !
               !---------------------------------------------------------------------------!
               htopcrown = cpatch%hite(ico)
               hbotcrown = h2trunkh(cpatch%hite(ico),ipft)
               ladcohort = (cpatch%lai(ico) + cpatch%wai(ico)) / (htopcrown - hbotcrown)
               kapartial = floor  (hbotcrown/dz_m97) + 1
               kafull    = ceiling(hbotcrown/dz_m97) + 1
               kzpartial = ceiling(htopcrown/dz_m97)
               kzfull    = floor  (htopcrown/dz_m97)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Add the LAD for the full layers.                                      !
               !---------------------------------------------------------------------------!
               do k = kafull,kzfull
                  lad(k) = lad(k) + ladcohort
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the LAD for the partial layers.  The only special case is when   !
               ! they are both the same layer, which must be done separately.              !
               !---------------------------------------------------------------------------!
               if (kapartial == kzpartial) then
                  lad(kapartial) = lad(kapartial)                                          &
                                 + ladcohort * (htopcrown        - hbotcrown         )     &
                                             / (zztop(kapartial) - zztop(kapartial-1))
               else
                  lad(kapartial) = lad(kapartial)                                          &
                                 + ladcohort * (zztop(kapartial) - hbotcrown         )     &
                                             / (zztop(kapartial) - zztop(kapartial-1))
                  lad(kzpartial) = lad(kzpartial)                                          &
                                 + ladcohort * (htopcrown        - zztop(kzpartial-1))     &
                                             / (zztop(kzpartial) - zztop(kzpartial-1))
               end if
               !---------------------------------------------------------------------------!
            end do
         end select
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     In this loop we find the drag coefficient, the sheltering factor, and the   !
         ! cumulative leaf drag area.  This must be done in a separate loop because we     !
         ! need the full integral of the leaf area density before we determine these       !
         ! variables.                                                                      !
         !---------------------------------------------------------------------------------!
         select case (icanturb)
         case (2)
            !----- Constant drag and sheltering factor. -----------------------------------!
            cdrag   (:) = cdrag0
            pshelter(:) = pm0
            ldga_bk     = 0.0
            do k = 1,zcan
               !---------------------------------------------------------------------------!
               !     Add the contribution of this layer to Massman's zeta (which we call   !
               ! cumldrag here to not confuse with the other zeta from the similarity      !
               ! theory.  We integrate in three steps so we save the value in the middle   !
               ! of the layer.                                                             !
               !---------------------------------------------------------------------------!
               lyrhalf      = 0.5 * lad(k) * cdrag(k) / pshelter(k) * dz_m97
               cumldrag(k)  = ldga_bk + lyrhalf
               ldga_bk      = ldga_bk + 2.0 * lyrhalf
               !---------------------------------------------------------------------------!
            end do
         case (3)
            ldga_bk      = 0.0
            !----- Drag and sheltering factor are height-dependent. -----------------------!
            do k = 1,zcan
               !------ Functional form as in Massman (1997). ------------------------------!
               cdrag   (k) = cdrag0 * exp( - alpha2_m97 * (1.0 - zzmid(k) / htop))
               pshelter(k) = pm0 / (1.0 + alpha1_m97 * htop * lad(k))

               !---------------------------------------------------------------------------!
               !     Add the contribution of this layer to Massman's zeta (which we call   !
               ! cumldrag here to not confuse with the other zeta from the similarity      !
               ! theory.  We integrate in three steps so we save the value in the middle   !
               ! of the layer.                                                             !
               !---------------------------------------------------------------------------!
               lyrhalf     = 0.5 * lad(k) * cdrag(k) / pshelter(k) * dz_m97
               cumldrag(k) = ldga_bk + lyrhalf
               ldga_bk     = ldga_bk + 2.0 * lyrhalf
               !---------------------------------------------------------------------------!
            end do
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the ratio between u* and u at the top cohort, using Massman's equation  !
         ! (6).                                                                            !
         !---------------------------------------------------------------------------------!
         ustarouh = (c1_m97 - c2_m97 * exp(-c3_m97 * cumldrag(zcan)))
         !---------------------------------------------------------------------------------!



         !----- NN is Massman's n, the coefficient of attenuation. ------------------------!
         nn = 0.5 * cumldrag(zcan) / (ustarouh * ustarouh)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the ratio between zero-plane displacement height and height, and the   !
         ! ratio between roughness and height.  The actual values will be defined based    !
         ! on the patch-level averaged height, not the height of the tallest cohort,       !
         ! because it may be an extremely sparse cohort that has little impact on the      !
         ! drag.                                                                           !
         !---------------------------------------------------------------------------------!
         d0ohgt = 1.0
         do k=1,zcan
            d0ohgt = d0ohgt - dz_m97 / htop                                                &
                            * exp(-2.0 * nn *(1.0 - cumldrag(k) / cumldrag(zcan)))
         end do
         z0ohgt = (1.0 - d0ohgt) * exp(- vonk / ustarouh + psi_m97)
         !---------------------------------------------------------------------------------!




         !----- Find the actual displacement height and roughness. ------------------------!
         csite%veg_displace(ipa) = max(0.,d0ohgt) * csite%veg_height(ipa)
         csite%rough(ipa)        = max(soil_rough, z0ohgt * csite%veg_height(ipa))
         !---------------------------------------------------------------------------------!



         !----- Find the characteristic scales (a.k.a. stars). ----------------------------!
         call ed_stars(cmet%atm_theta,cmet%atm_theiv,cmet%atm_shv,cmet%atm_co2             &
                      ,csite%can_theta(ipa),csite%can_theiv(ipa),csite%can_shv(ipa)        &
                      ,csite%can_co2(ipa),cmet%geoht,csite%veg_displace(ipa),cmet%vels     &
                      ,csite%rough(ipa),csite%ustar(ipa),csite%tstar(ipa),estar            &
                      ,csite%qstar(ipa),csite%cstar(ipa),csite%zeta(ipa),csite%ribulk(ipa) &
                      ,csite%ggbare(ipa))
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate the leaf level aerodynamic resistance.                            !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh = reduced_wind(csite%ustar(ipa),csite%zeta(ipa),csite%ribulk(ipa),cmet%geoht   &
                          ,csite%veg_displace(ipa),htop,csite%rough(ipa))
         !---------------------------------------------------------------------------------!
         do ico=1,cpatch%ncohorts
            ipft      = cpatch%pft(ico)

            !----- Find the crown relevant heights. ---------------------------------------!
            htopcrown = cpatch%hite(ico)
            hbotcrown = h2trunkh(cpatch%hite(ico),cpatch%pft(ico))
            hmidcrown = 0.5 * (hbotcrown + htopcrown)
            !------------------------------------------------------------------------------!



            !----- Determine which layer we should use for wind reduction. ----------------!
            k = max(1,ceiling(hmidcrown / dz_m97))
            !------------------------------------------------------------------------------!



            !----- Calculate the wind speed at height z. ----------------------------------!
            cpatch%veg_wind(ico) = max( ugbmin                                             &
                                      , uh * exp(-nn * (1. - cumldrag(k)/cumldrag(zcan))))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find the aerodynamic conductances for heat and water at the leaf boundary !
            ! layer.                                                                       !
            !------------------------------------------------------------------------------!
            if (cpatch%resolvable(ico)) then
               call aerodynamic_conductances(ipft,cpatch%veg_wind(ico)                     &
                                            ,cpatch%veg_temp(ico),csite%can_temp(ipa)      &
                                            ,csite%can_shv(ipa),csite%can_rhos(ipa)        &
                                            ,gbhmos_min,cpatch%gbh(ico),cpatch%gbw(ico))
            else
               cpatch%gbh(ico) = 0.0
               cpatch%gbw(ico) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the bulk resistance by integrating the inverse of the diffusivity for  !
         ! each layer (e.g. Sellers et al. (1986)).  We assumed Km = sigma u, which is the !
         ! preferred approach by Sellers et al. (1986).  To find sigma we used that        !
         ! tau = rho * ustar^2 at z=h, and that tau = rho * Km * dU/dz.                    !
         !---------------------------------------------------------------------------------!
         rasveg  = 0.0
         sigmakm = ustarouh * ustarouh * cumldrag(zcan) * pshelter(zcan)                   &
                 / ( nn * cdrag(zcan) * lad(zcan) )
         do k=1,zcan
            windm97(k) = max(ugbmin, uh * exp(- nn * (1.0 - cumldrag(k)/cumldrag(zcan)) ))
            Kdiff      = sigmakm * windm97(k) + kvwake
            rasveg     = rasveg + dz_m97 / Kdiff
         end do
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         !    According to Sellers et al (1986), equation A15, we must correct the         !
         ! resistance for the unstable case.  There is a typo on their equation (a gravity !
         ! term is missing, and g appears in their code in SiB-2.0 that is coupled to the  !
         ! original BRAMS-4.0).  Instead of their approximation that is ultimately based   !
         ! on Businger et al. (1971), we find an estimate of zeta using the same method    !
         ! that is applied to get zeta for the Beljaars-Holtslag method.  According to     !
         ! Sellers et al. (1986), we should apply this correction for the unstable case    !
         ! only.                                                                           !
         !---------------------------------------------------------------------------------!
         if (csite%ground_temp(ipa) > csite%can_temp(ipa)) then
            ribcan   = 2. * grav * (csite%can_temp(ipa) - csite%ground_temp(ipa))          &
                     * htop * (1. - d0ohgt - z0ohgt)                                       &
                     / ( (csite%can_temp(ipa) + csite%ground_temp(ipa)) * uh * uh)
            hgtoz0   = (1. - d0ohgt) / z0ohgt
            lnhgtoz0 = log(hgtoz0)
            zetacan  = zoobukhov(ribcan,htop * (1.0 - d0ohgt),htop * z0ohgt                &
                                ,hgtoz0,lnhgtoz0,hgtoz0,lnhgtoz0,.false.)
            phih     = sqrt(1.0 - gamh * zetacan)
         else
            !----- Stable case, no correction needed. -------------------------------------!
            phih     = 1.0
         end if
         csite%ggveg(ipa) = phih / rasveg
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (csite%opencan_frac(ipa) > 0.999 .or. csite%snowfac(ipa) >= 0.9) then
            csite%ggnet(ipa) = ggfact * csite%ggbare(ipa)
         else
            csite%ggnet(ipa) = ggfact * csite%ggbare(ipa) * csite%ggveg(ipa)               &
                             / ( csite%ggveg(ipa) + (1. - csite%opencan_frac(ipa))         &
                                                  * csite%ggbare(ipa))
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap(csite%can_rhos(ipa),csite%can_temp(ipa),csite%can_depth(ipa)       &
                       ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!


      return
   end subroutine canopy_turbulence
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This is the double precision version of the above sub-routine.                    !
   !                                                                                       !
   !     The following subroutine calculates the turbulent transfer of scalar quantities   !
   ! using Monin Obukhov Similarity Theory under a few different flavors.  The basic       !
   ! assumptions are that the ground surface and canopy space is modeled as a rough wall   !
   ! boundary above which exists a dynamic sublayer, above which exists the free atmo-     !
   ! sphere.  We assume a logarithmic decay of wind speed over the dynamic sublayer, which !
   ! is due to the assumption of constant shear (ustar) in the dynamic sublayer.           !
   !                                                                                       !
   ! 0.  We apply the wind profile based on the similarity theory for the tallest cohort,  !
   !     with further reduction due to the presence of obstacles (trees), similar to       !
   !     Leuning et al. (1995) but considering a finite crown area.  The aerodynamic       !
   !     resistance from ground to canopy air is found in a similar way as in LEAF-3.      !
   !     Displacement height and roughness are linear functions of the vegetation height,  !
   !     where vegetation height is the weighted average of cohort heights, the weights    !
   !     being the basal area.                                                             !
   !     References:                                                                       !
   !                                                                                       !
   !     Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf          !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   !                                                                                       !
   ! 1. This is similar to option 0, where displacement height and rougness are linear     !
   !    functions of mean canopy height, where mean canopy height is defined by the        !
   !    dead-biomass weighted average of canopy height.  Wind speed and diffusivity decay  !
   !    exponentially in the canopy on parameter alpha.  HOWEVER, if the windspeed at      !
   !    reference height is less than the canopy height (WHICH IS REALLY INNAPROPRIATE TO  !
   !    BEGIN WITH, BUT OH WELL), then revert to a calculation of ustar that uses no       !
   !    displacement height.  After ustar is calculated, then use displacement height in   !
   !    the calculation of the eddy length scale at the canopy top to find diffusivity at  !
   !    canopy top. THIS IS THE METHOD USED IN LEAF3 AND ED2.0.                            !
   !                                                                                       !
   ! 2. Follow the methods of Massman 1997.  The canopy top is taken as the height of the  !
   !    tallest cohort in the patch. The displacement height is determined by a center of  !
   !    pressure method. The wind speed inside the canopy decays exponentially on          !
   !    cumulative foliar drag elements, where the equation is scaled such that it solves  !
   !    for the balance of integrated canopy shear with dissipation.  To close the         !
   !    equations M97 assumes that the surface drag (U(h)/Ustar) is a function of the      !
   !    cumulative foliar drag surfaces.  With this, we can explicitly solve for roughness !
   !    length. Assumptions are made about the profile of canopy drag as a function of     !
   !    cohort attributes, please see the technical manual by RGK2009 for more info.       !
   !    Please also refer to the original paper by Massman for information regarding the   !
   !    basic principles of the closure scheme.                                            !
   !                                                                                       !
   !    Massman, WJ. An Analytical One-Dimensional Model of Momentum Transfer By           !
   !      Vegetation Of Arbitrary Structure. Boundary Layer Meteorology, 83,               !
   !      p. 407-421, 1997.                                                                !
   !                                                                                       !
   ! 3.  Same as 2, but the drag and sheltering parameters are functions of height.        !
   !                                                                                       !
   !     Ultimately, this routine solves for the resistance of water vapor and sensible    !
   ! heat from the soil surface to canopy air space and from the leaf surfaces to canopy   !
   ! air space.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_turbulence8(csite,initp,ipa)
      use ed_state_vars  , only : polygontype          & ! structure
                                , sitetype             & ! structure
                                , patchtype            ! ! structure
      use rk4_coms       , only : ibranch_thermo       & ! intent(in)
                                , rk4patchtype         & ! structure
                                , rk4site              & ! intent(in)
                                , tiny_offset          & ! intent(in)
                                , ibranch_thermo       & ! intent(in)
                                , wcapcan              & ! intent(out)
                                , wcapcani             & ! intent(out)
                                , hcapcani             & ! intent(out)
                                , ccapcani             ! ! intent(out)
      use canopy_air_coms, only : icanturb             & ! intent(in), can. turb. scheme
                                , ustmin8              & ! intent(in)
                                , ugbmin8              & ! intent(in)
                                , ubmin8               & ! intent(in)
                                , exar8                & ! intent(in)
                                , gamh8                & ! intent(in)
                                , dz_m978              & ! intent(in)
                                , cdrag08              & ! intent(in)
                                , pm08                 & ! intent(in)
                                , c1_m978              & ! intent(in)
                                , c2_m978              & ! intent(in)
                                , c3_m978              & ! intent(in)
                                , kvwake8              & ! intent(in)
                                , alpha1_m978          & ! intent(in)
                                , alpha2_m978          & ! intent(in)
                                , psi_m978             & ! intent(in)
                                , zztop8               & ! intent(in)
                                , zzmid8               & ! intent(in)
                                , lad8                 & ! intent(in)
                                , dladdz8              & ! intent(in)
                                , cdrag8               & ! intent(in)
                                , pshelter8            & ! intent(in)
                                , cumldrag8            & ! intent(in)
                                , windm978             & ! intent(in)
                                , rb_inter             & ! intent(in)
                                , rb_slope             & ! intent(in)
                                , ggfact8              & ! intent(in)
                                , zoobukhov8           ! ! intent(in)
      use consts_coms    , only : vonk8                & ! intent(in)
                                , cpi8                 & ! intent(in)
                                , grav8                & ! intent(in)
                                , epim18               & ! intent(in)
                                , sqrt2o28             ! ! intent(in)
      use soil_coms      , only : snow_rough8          & ! intent(in)
                                , soil_rough8          ! ! intent(in)
      use allometry      , only : h2trunkh             & ! function
                                , dbh2bl               ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)     , target     :: csite         ! Current site
      type(rk4patchtype) , target     :: initp         ! Current integrator
      integer            , intent(in) :: ipa           ! Patch loop
      !----- Pointers ---------------------------------------------------------------------!
      type(patchtype)    , pointer    :: cpatch        ! Current patch
      !----- Local variables --------------------------------------------------------------!
      integer        :: ico        ! Cohort loop
      integer        :: ipft       ! PFT alias
      integer        :: k          ! Elevation index
      integer        :: zcan       ! Index of canopy top elevation
      integer        :: kafull     ! First layer fully occupied by crown
      integer        :: kzfull     ! Last layer fully occupied by crown
      integer        :: kapartial  ! First layer partially occupied by crown
      integer        :: kzpartial  ! Last layer partially occupied by crown
      logical        :: stable     ! Stable canopy air space
      real(kind=8)   :: rasveg     ! Resistance of vegetated ground             [      s/m]
      real(kind=8)   :: atm_thetav ! Free atmosphere virtual potential temp.    [        K]
      real(kind=8)   :: can_thetav ! Free atmosphere virtual potential temp.    [        K]
      real(kind=8)   :: ldga_bk    ! Cumulative zeta function                   [      ---]
      real(kind=8)   :: lyrhalf    ! Half the contrib. of this layer to zeta    [      1/m]
      real(kind=8)   :: sigmakm    ! Km coefficient at z=h                      [        m]
      real(kind=8)   :: K_top      ! Diffusivity at canopy top z=h              [     m2/s]
      real(kind=8)   :: kdiff      ! Diffusivity                                [     m2/s]
      real(kind=8)   :: surf_rough ! Roughness length of the bare ground 
                                   !     at canopy bottom                       [        m]
      real(kind=8)   :: uh         ! Wind speed at the canopy top (z=h)         [      m/s]
      real(kind=8)   :: factv      ! Wind-dependent term for old rasveg
      real(kind=8)   :: aux        ! Aux. variable
      real(kind=8)   :: gbhmos_min ! Minimum leaf boundary layer heat condct.   [      m/s]
      real(kind=8)   :: ustarouh   ! The ratio of ustar over u(h)               [      ---]
      real(kind=8)   :: nn         ! In-canopy wind attenuation scal. param.    [      ---]
      real(kind=8)   :: waiuse     ! Wood area index                            [    m2/m2]
      real(kind=8)   :: htopcrown  ! height at the top of the crown             [        m]
      real(kind=8)   :: hmidcrown  ! Height at the middle of the crown          [        m]
      real(kind=8)   :: hbotcrown  ! Height at the bottom of the crown          [        m]
      real(kind=8)   :: htop       ! Height of the topmost layer                [        m]
      real(kind=8)   :: dzcrown    ! Depth that contains leaves/branches        [        m]
      real(kind=8)   :: d0ohgt     ! d0/height                                  [      ---]
      real(kind=8)   :: z0ohgt     ! z0/height                                  [      ---]
      real(kind=8)   :: phih       ! Correction term for unstable cases         [      ---]
      real(kind=8)   :: ribcan     ! Ground-to-canopy bulk Richardson number    [      ---]
      real(kind=8)   :: hgtoz0     ! height/z0                                  [      ---]
      real(kind=8)   :: lnhgtoz0   ! log(height/z0)                             [      ---]
      real(kind=8)   :: zetacan    ! Estimate of z/L within the canopy          [      ---]
      real(kind=8)   :: ladcohort  ! Leaf Area Density of this cohort           [    m2/m3]
      !------ External procedures ---------------------------------------------------------!
      real        , external             :: sngloff  ! Safe double -> simple precision.
      !------------------------------------------------------------------------------------!



      !----- Assign some pointers. --------------------------------------------------------!
      cpatch=>csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the virtual potential temperatures and decide whether the canopy air is   !
      ! stable or not.                                                                     !
      !------------------------------------------------------------------------------------!
      atm_thetav = rk4site%atm_theta * (1.d0 + epim18 * rk4site%atm_shv)
      can_thetav = initp%can_theta   * (1.d0 + epim18 * initp%can_shv  )
      stable     = atm_thetav >= can_thetav

      !------------------------------------------------------------------------------------!
      !     If there is no vegetation in this patch, then we apply turbulence to bare      !
      ! soil, no d0 and exit.                                                              !
      !------------------------------------------------------------------------------------!
      if (cpatch%ncohorts == 0) then

         !----- Calculate the surface roughness inside the canopy. ------------------------!
         initp%rough = soil_rough8 *(1.d0 - initp%snowfac) + snow_rough8 * initp%snowfac
         
         !----- Finding the characteristic scales (a.k.a. stars). -------------------------!
         call ed_stars8(rk4site%atm_theta,rk4site%atm_theiv,rk4site%atm_shv                &
                       ,rk4site%atm_co2,initp%can_theta ,initp%can_theiv ,initp%can_shv    &
                       ,initp%can_co2,rk4site%geoht,initp%veg_displace,rk4site%vels        &
                       ,initp%rough,initp%ustar,initp%tstar,initp%estar,initp%qstar        &
                       ,initp%cstar,initp%zeta,initp%ribulk,initp%ggbare)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      This patch is empty, so we can't define a conductance for vegetated        !
         ! grounds.  Assign it to be zero, and set the net conductance to be the bare      !
         ! ground.                                                                         !
         !---------------------------------------------------------------------------------!
         initp%ggveg = 0.d0
         initp%ggnet = ggfact8 * initp%ggbare
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap8(initp%can_rhos,initp%can_temp,initp%can_depth                     &
                        ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!
         
         return
      end if
      !------------------------------------------------------------------------------------!




      !----- Check whether the reference height is high enough . --------------------------!
      if (dble(cpatch%hite(1)) >= rk4site%geoht) then
         write (unit=*,fmt='(a)') ' Your reference height is too low...'
         write (unit=*,fmt='(a,1x,es12.5)') ' Ref. height:           ',rk4site%geoht
         write (unit=*,fmt='(a,1x,es12.5)') ' Tallest cohort height: ',cpatch%hite(1)
         call fatal_error('Bad reference height','canopy_turbulence8'                      &
                         ,'canopy_struct_dynamics.f90')
      end if
      !------------------------------------------------------------------------------------!


      !---- Find the minimum boundary layer heat conductance. -----------------------------!
      if (any(initp%resolvable)) then
         gbhmos_min = 1.d0 / ( dble(rb_inter) + dble(rb_slope)                             &
                             * (dble(csite%lai(ipa)) + dble(csite%wai(ipa))))
      else
         gbhmos_min = 0.d0
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we do have cohorts, choose which method we use to compute the          !
      ! resistance.                                                                        !
      !------------------------------------------------------------------------------------!
      select case (icanturb)



      !------------------------------------------------------------------------------------!
      ! LEAF-3 Case: This approach is very similar to older implementations of ED-2, and   !
      !              it is very similar to LEAF-3, and to option 1, except that option 1   !
      !              computes the vegetated ground conductance and wind profile different- !
      !              ly.                                                                   !
      !------------------------------------------------------------------------------------!
      case (0) 

         !---------------------------------------------------------------------------------!
         !     The roughness is found by combining two weighted averages.  The first one   !
         ! checks the fraction of the patch that has closed canopy, and averages between   !
         ! soil and vegetation roughness.  The other is the fraction of the vegetation     !
         ! that is covered in snow.                                                        !
         !---------------------------------------------------------------------------------!
         initp%rough = snow_rough8 * initp%snowfac                                         &
                     + ( soil_rough8     * initp%opencan_frac                              &
                       + initp%veg_rough * (1.d0 - initp%opencan_frac) )                   &
                     * (1.d0 - initp%snowfac)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!
         call ed_stars8(rk4site%atm_theta,rk4site%atm_theiv,rk4site%atm_shv                &
                       ,rk4site%atm_co2,initp%can_theta ,initp%can_theiv,initp%can_shv     &
                       ,initp%can_co2,rk4site%geoht,initp%veg_displace,rk4site%vels        &
                       ,initp%rough,initp%ustar,initp%tstar,initp%estar,initp%qstar        &
                       ,initp%cstar,initp%zeta,initp%ribulk,initp%ggbare)
         !---------------------------------------------------------------------------------!

         if (initp%snowfac < 9.d-1) then
            factv        = log((rk4site%geoht-initp%veg_displace) / initp%rough)           &
                         / (vonk8 * vonk8 * rk4site%vels)
            aux          = exp(exar8 * (1.d0 - (initp%veg_displace + initp%rough)          &
                                             / initp%veg_height) )
            initp%ggveg  = (exar8 * (initp%veg_height - initp%veg_displace))               &
                         / (factv * initp%veg_displace * (exp(exar8) - aux))
         else 
            initp%ggveg = 0.d0
         end if


         !---------------------------------------------------------------------------------!
         !     This part of the code initializes the geometry of the canopy air space, the !
         ! structure of the vegetation and its attenuation effects and the heat and water  !
         ! capacities.                                                                     !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh = reduced_wind8(initp%ustar,initp%zeta,initp%ribulk,rk4site%geoht              &
                           ,initp%veg_displace,dble(cpatch%hite(1)),initp%rough)

         do ico=1,cpatch%ncohorts

            ipft  = cpatch%pft(ico)

            !----- Calculate the wind speed at height z. ----------------------------------!
            initp%veg_wind(ico) = max(uh,ugbmin8)
            uh = uh * ( initp%crown_area(ico)                                              &
                      * exp(- 5.d-1 * initp%lai(ico) / initp%crown_area(ico))              &
                      + 1.d0 - initp%crown_area(ico) )
            !------------------------------------------------------------------------------!

            if (initp%resolvable(ico)) then


               !---------------------------------------------------------------------------!
               !    Find the aerodynamic conductances for heat and water at the leaf       !
               ! boundary layer.                                                           !
               !---------------------------------------------------------------------------!
               call aerodynamic_conductances8(ipft,initp%veg_wind(ico),initp%veg_temp(ico) &
                                             ,initp%can_temp,initp%can_shv,initp%can_rhos  &
                                             ,gbhmos_min,initp%gbh(ico),initp%gbw(ico)     &
                                             ,initp%veg_reynolds(ico)                      &
                                             ,initp%veg_grashof(ico)                       &
                                             ,initp%veg_nussfree(ico)                      &
                                             ,initp%veg_nussforc(ico) )
               cpatch%gbh(ico) = sngloff(initp%gbh(ico),tiny_offset)
               cpatch%gbw(ico) = sngloff(initp%gbw(ico),tiny_offset)
               !---------------------------------------------------------------------------!
            else
               initp%veg_reynolds(ico) = 0.d0
               initp%veg_grashof (ico) = 0.d0
               initp%veg_nussfree(ico) = 0.d0
               initp%veg_nussforc(ico) = 0.d0
               initp%gbh(ico)          = 0.d0
               initp%gbw(ico)          = 0.d0
               cpatch%gbh(ico)         = 0.0
               cpatch%gbw(ico)         = 0.0
            end if
         end do


         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap8(initp%can_rhos,initp%can_temp,initp%can_depth                     &
                        ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (initp%opencan_frac > 9.99d-1 .or. initp%snowfac >= 9.d-1) then
            initp%ggnet = ggfact8 * initp%ggbare
         else
            initp%ggnet = ggfact8 * initp%ggbare * initp%ggveg                             &
                        / ( initp%ggveg + (1.d0 - initp%opencan_frac) * initp%ggbare )
         end if
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      ! ED-2.1 Case: There uses a slightly different approach for wind profile and ground  !
      !              conductance.                                                          !
      !------------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     The roughness is found by combining two weighted averages.  The first one   !
         ! checks the fraction of the patch that has closed canopy, and averages between   !
         ! soil and vegetation roughness.  The other is the fraction of the vegetation     !
         ! that is covered in snow.                                                        !
         !---------------------------------------------------------------------------------!
         initp%rough = snow_rough8 * initp%snowfac                                         &
                     + ( soil_rough8     * initp%opencan_frac                              &
                       + initp%veg_rough * (1.d0 - initp%opencan_frac) )                   &
                     * (1.d0 - initp%snowfac)
         !---------------------------------------------------------------------------------!
         
         !----- Calculate the soil surface roughness inside the canopy. -------------------!
         surf_rough = soil_rough8 * (1.d0 - initp%snowfac) + snow_rough8 * initp%snowfac
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!
         call ed_stars8(rk4site%atm_theta,rk4site%atm_theiv,rk4site%atm_shv                &
                       ,rk4site%atm_co2,initp%can_theta,initp%can_theiv,initp%can_shv      &
                       ,initp%can_co2,rk4site%geoht,initp%veg_displace,rk4site%vels        &
                       ,initp%rough,initp%ustar,initp%tstar,initp%estar,initp%qstar        &
                       ,initp%cstar,initp%zeta,initp%ribulk,initp%ggbare)


         !---------------------------------------------------------------------------------!
         !     The surface conductance in the sub-canopy layer is the integration of       !
         ! K, from the rough soil surface, to the reference point in the canopy            !
         ! where the state variable is integrated (canopy top ~ h).                        !
         !---------------------------------------------------------------------------------!
         K_top = vonk8 * initp%ustar * (initp%veg_height - initp%veg_displace)
         initp%ggveg = (exar8 * K_top)                                                     &
                     / ( exp(exar8) * initp%veg_height                                     &
                       * (exp(-exar8 * surf_rough/initp%veg_height) - exp(-exar8)))
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     This part of the code finds the geometry of the canopy air space, the       !
         ! structure of the vegetation and its attenuation effects and the heat and water  !
         ! capacities.                                                                     !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh   = reduced_wind8(initp%ustar,initp%zeta,initp%ribulk,rk4site%geoht            &
                             ,initp%veg_displace,dble(cpatch%hite(1)),initp%rough)
         htop = dble(cpatch%hite(1))
         do ico=1,cpatch%ncohorts

            !----- Alias for PFT type. ----------------------------------------------------!
            ipft  = cpatch%pft(ico)
            !------------------------------------------------------------------------------!



            !----- Estimate the height center of the crown. -------------------------------!
            htopcrown = dble(cpatch%hite(ico))
            hbotcrown = dble(h2trunkh(cpatch%hite(ico),ipft))
            hmidcrown = 5.d-1 * (htopcrown + hbotcrown)
            !------------------------------------------------------------------------------!



            !----- Calculate the wind speed at height z. ----------------------------------!
            initp%veg_wind(ico) = max(ugbmin8                                              &
                                     , uh * exp(-exar8 * (1.d0 - hmidcrown/ htop)))

            !------------------------------------------------------------------------------!
            !    Find the aerodynamic conductances for heat and water at the leaf boundary !
            ! layer.                                                                       !
            !------------------------------------------------------------------------------!
            if (initp%resolvable(ico)) then
               call aerodynamic_conductances8(ipft,initp%veg_wind(ico),initp%veg_temp(ico) &
                                             ,initp%can_temp,initp%can_shv,initp%can_rhos  &
                                             ,gbhmos_min,initp%gbh(ico),initp%gbw(ico)     &
                                             ,initp%veg_reynolds(ico)                      &
                                             ,initp%veg_grashof(ico)                       &
                                             ,initp%veg_nussfree(ico)                      &
                                             ,initp%veg_nussforc(ico))
               cpatch%gbh(ico) = sngloff(initp%gbh(ico),tiny_offset)
               cpatch%gbw(ico) = sngloff(initp%gbw(ico),tiny_offset)
            else
               initp%veg_reynolds(ico) = 0.d0
               initp%veg_grashof (ico) = 0.d0
               initp%veg_nussfree(ico) = 0.d0
               initp%veg_nussforc(ico) = 0.d0
               initp%gbh         (ico) = 0.d0
               initp%gbw         (ico) = 0.d0
               cpatch%gbh        (ico) = 0.0
               cpatch%gbw        (ico) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap8(initp%can_rhos,initp%can_temp,initp%can_depth                     &
                        ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (initp%opencan_frac > 9.99d-1 .or. initp%snowfac >= 9.d-1) then
            initp%ggnet = ggfact8 * initp%ggbare
         else
            initp%ggnet = ggfact8 * initp%ggbare * initp%ggveg                             &
                        / (initp%ggveg + (1.d0 - initp%opencan_frac) * initp%ggbare )
         end if
         !---------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Use the methods of Massman (1997).  Using option 2 will make the within-      !
      ! canopy drag and sheltering factors to be constant, whereas option 3 will allow the !
      ! values to vary.                                                                    !
      !------------------------------------------------------------------------------------!
      case (2,3)
         !---------------------------------------------------------------------------------!
         !    Find the top layer and the top height.                                       !
         !---------------------------------------------------------------------------------!
         zcan = ceiling(dble(cpatch%hite(1)) / dz_m978)
         htop = zztop8(zcan)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Loop through cohorts, and integrate the leaf area density.  Notice that we  !
         ! will be solving all cohorts here, because even those with no leaves have        !
         ! effects on the canopy air turbulence.  We only need to decide whether we have   !
         ! branches or not, and we do it outside the loop to make it more efficient.       !
         !---------------------------------------------------------------------------------!
         !----- Reset the leaf area density array. ----------------------------------------!
         lad8(:) = 0.d0
         select case (ibranch_thermo)
         case (0)
            !------------------------------------------------------------------------------!
            !     Although we are not solving branch thermodynamics, we have some drag     !
            ! effect due to the branches.  Assume branch area index to be 15% of on-       !
            ! -allometry LAI.  Even for grasses we may want to keep some residual due to   !
            ! dead material that remains around.                                           !
            !------------------------------------------------------------------------------!
            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)


               !------ Estimate WAI. ------------------------------------------------------!
               waiuse = 1.d-1 * initp%nplant(ico) * dble(cpatch%sla(ico))                  &
                      * dble(dbh2bl(cpatch%dbh(ico),ipft))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the heights, and compute the LAD of this cohort.                 !
               !---------------------------------------------------------------------------!
               htopcrown = dble(cpatch%hite(ico))
               hbotcrown = dble(h2trunkh(cpatch%hite(ico),cpatch%pft(ico)))
               ladcohort = (initp%lai(ico) + waiuse) / (htopcrown - hbotcrown)
               kapartial = floor  (hbotcrown/dz_m978) + 1
               kafull    = ceiling(hbotcrown/dz_m978) + 1
               kzpartial = ceiling(htopcrown/dz_m978)
               kzfull    = floor  (htopcrown/dz_m978)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Add the LAD for the full layers.                                      !
               !---------------------------------------------------------------------------!
               do k = kafull,kzfull
                  lad8(k) = lad8(k) + ladcohort
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the LAD for the partial layers.  The only special case is when   !
               ! they are both the same layer, which must be done separately.              !
               !---------------------------------------------------------------------------!
               if (kapartial == kzpartial) then
                  lad8(kapartial) = lad8(kapartial)                                        &
                                  + ladcohort * (htopcrown         - hbotcrown          )  &
                                              / (zztop8(kapartial) - zztop8(kapartial-1))
               else
                  lad8(kapartial) = lad8(kapartial)                                        &
                                  + ladcohort * (zztop8(kapartial) - hbotcrown          )  &
                                              / (zztop8(kapartial) - zztop8(kapartial-1))
                  lad8(kzpartial) = lad8(kzpartial)                                        &
                                  + ladcohort * (htopcrown         - zztop8(kzpartial-1))  &
                                              / (zztop8(kzpartial) - zztop8(kzpartial-1))
               end if
               !---------------------------------------------------------------------------!
            end do

         case default
            !----- Use the default wood area index. ---------------------------------------!
            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)

               !---------------------------------------------------------------------------!
               !     Find the heights, and compute the LAD of this cohort.                 !
               !---------------------------------------------------------------------------!
               htopcrown = dble(cpatch%hite(ico))
               hbotcrown = dble(h2trunkh(cpatch%hite(ico),ipft))
               ladcohort = (initp%lai(ico) + initp%wai(ico)) / (htopcrown - hbotcrown)
               kapartial = floor  (hbotcrown/dz_m978) + 1
               kafull    = ceiling(hbotcrown/dz_m978) + 1
               kzpartial = ceiling(htopcrown/dz_m978)
               kzfull    = floor  (htopcrown/dz_m978)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Add the LAD for the full layers.                                      !
               !---------------------------------------------------------------------------!
               do k = kafull,kzfull
                  lad8(k) = lad8(k) + ladcohort
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Add the LAD for the partial layers.  The only special case is when   !
               ! they are both the same layer, which must be done separately.              !
               !---------------------------------------------------------------------------!
               if (kapartial == kzpartial) then
                  lad8(kapartial) = lad8(kapartial)                                        &
                                  + ladcohort * (htopcrown         - hbotcrown          )  &
                                              / (zztop8(kapartial) - zztop8(kapartial-1))
               else
                  lad8(kapartial) = lad8(kapartial)                                        &
                                  + ladcohort * (zztop8(kapartial) - hbotcrown          )  &
                                              / (zztop8(kapartial) - zztop8(kapartial-1))
                  lad8(kzpartial) = lad8(kzpartial)                                        &
                                  + ladcohort * (htopcrown         - zztop8(kzpartial-1))  &
                                              / (zztop8(kzpartial) - zztop8(kzpartial-1))
               end if
               !---------------------------------------------------------------------------!
            end do
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     In this loop we find the drag coefficient, the sheltering factor, and the   !
         ! cumulative leaf drag area.  This must be done in a separate loop because we     !
         ! need the full integral of the leaf area density before we determine these       !
         ! variables.                                                                      !
         !---------------------------------------------------------------------------------!
         select case (icanturb)
         case (2)
            !----- Constant drag and sheltering factor. -----------------------------------!
            cdrag8   (:) = cdrag08
            pshelter8(:) = pm08
            ldga_bk      = 0.0
            do k = 1,zcan
               !---------------------------------------------------------------------------!
               !     Add the contribution of this layer to Massman's zeta (which we call   !
               ! cumldrag here to not confuse with the other zeta from the similarity      !
               ! theory.  We integrate in three steps so we save the value in the middle   !
               ! of the layer.                                                             !
               !---------------------------------------------------------------------------!
               lyrhalf      = 5.d-1 * lad8(k) * cdrag8(k) / pshelter8(k) * dz_m978
               cumldrag8(k) = ldga_bk + lyrhalf
               ldga_bk      = ldga_bk + 2.d0 * lyrhalf
               !---------------------------------------------------------------------------!
            end do
         case (3)
            ldga_bk      = 0.0
            !----- Drag and sheltering factor are height-dependent. -----------------------!
            do k = 1,zcan
               !------ Functional form as in Massman (1997). ------------------------------!
               cdrag8   (k) = cdrag08 * exp( - alpha2_m978 * (1.d0 - zzmid8(k) / htop))
               pshelter8(k) = pm08 / (1.d0 + alpha1_m978 * htop * lad8(k))

               !---------------------------------------------------------------------------!
               !     Add the contribution of this layer to Massman's zeta (which we call   !
               ! cumldrag here to not confuse with the other zeta from the similarity      !
               ! theory.  We integrate in three steps so we save the value in the middle   !
               ! of the layer.                                                             !
               !---------------------------------------------------------------------------!
               lyrhalf      = 5.d-1 * lad8(k) * cdrag8(k) / pshelter8(k) * dz_m978
               cumldrag8(k) = ldga_bk + lyrhalf
               ldga_bk      = ldga_bk + 2.d0 * lyrhalf
               !---------------------------------------------------------------------------!
            end do
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find the ratio between u* and u at the top cohort, using Massman's equation  !
         ! (6).                                                                            !
         !---------------------------------------------------------------------------------!
         ustarouh = (c1_m978 - c2_m978 * exp(-c3_m978 * cumldrag8(zcan)))
         !---------------------------------------------------------------------------------!



         !----- NN is Massman's n, the coefficient of attenuation. ------------------------!
         nn = 5.d-1 * cumldrag8(zcan) / (ustarouh * ustarouh)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the ratio between zero-plane displacement height and height, and the   !
         ! ratio between roughness and height.  The actual values will be defined based    !
         ! on the patch-level averaged height, not the height of the tallest cohort,       !
         ! because it may be an extremely sparse cohort that has little impact on the      !
         ! drag.                                                                           !
         !---------------------------------------------------------------------------------!
         d0ohgt = 1.d0
         do k=1,zcan
            d0ohgt = d0ohgt - dz_m978 / htop                                               &
                            * exp(-2.d0 * nn *(1.d0 - cumldrag8(k) / cumldrag8(zcan)))
         end do
         z0ohgt = (1.d0 - d0ohgt) * exp(- vonk8 / ustarouh + psi_m978)
         !---------------------------------------------------------------------------------!




         !----- Find the actual displacement height and roughness. ------------------------!
         initp%veg_displace = max(0.d0,d0ohgt) * initp%veg_height
         initp%rough        = max(soil_rough8, z0ohgt * initp%veg_height)
         !---------------------------------------------------------------------------------!


         !----- Calculate ustar, tstar, qstar, and cstar. ---------------------------------!
         call ed_stars8(rk4site%atm_theta,rk4site%atm_theiv,rk4site%atm_shv                &
                       ,rk4site%atm_co2,initp%can_theta,initp%can_theiv,initp%can_shv      &
                       ,initp%can_co2,rk4site%geoht,initp%veg_displace,rk4site%vels        &
                       ,initp%rough,initp%ustar,initp%tstar,initp%estar,initp%qstar        &
                       ,initp%cstar,initp%zeta,initp%ribulk,initp%ggbare)
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Calculate the leaf level aerodynamic resistance.                            !
         !---------------------------------------------------------------------------------!
         !----- Top of canopy wind speed. -------------------------------------------------!
         uh = reduced_wind8(initp%ustar,initp%zeta,initp%ribulk,rk4site%geoht              &
                           ,initp%veg_displace,htop,initp%rough)
         do ico=1,cpatch%ncohorts
            ipft = cpatch%pft(ico)

            !----- Find the crown relevant heights. ---------------------------------------!
            htopcrown = dble(cpatch%hite(ico))
            hbotcrown = dble(h2trunkh(cpatch%hite(ico),cpatch%pft(ico)))
            hmidcrown = 5.d-1 * (hbotcrown + htopcrown)
            !------------------------------------------------------------------------------!



            !----- Determine which layer we should use for wind reduction. ----------------!
            k = max(1,ceiling(hmidcrown / dz_m978))
            !------------------------------------------------------------------------------!



            !----- Calculate the wind speed at height z. ----------------------------------!
            initp%veg_wind(ico) = max(ugbmin8                                              &
                                     ,uh * exp(-nn*(1.d0 - cumldrag8(k)/cumldrag8(zcan))))
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find the aerodynamic conductances for heat and water at the leaf boundary !
            ! layer.                                                                       !
            !------------------------------------------------------------------------------!
            if (initp%resolvable(ico)) then
               call aerodynamic_conductances8(ipft,initp%veg_wind(ico),initp%veg_temp(ico) &
                                             ,initp%can_temp,initp%can_shv,initp%can_rhos  &
                                             ,gbhmos_min,initp%gbh(ico),initp%gbw(ico)     &
                                             ,initp%veg_reynolds(ico)                      &
                                             ,initp%veg_grashof(ico)                       &
                                             ,initp%veg_nussfree(ico)                      &
                                             ,initp%veg_nussforc(ico) )
               cpatch%gbh(ico) = sngloff(initp%gbh(ico),tiny_offset)
               cpatch%gbw(ico) = sngloff(initp%gbw(ico),tiny_offset)
            else
               initp%veg_reynolds(ico) = 0.d0
               initp%veg_grashof (ico) = 0.d0
               initp%veg_nussfree(ico) = 0.d0
               initp%veg_nussforc(ico) = 0.d0
               initp%gbh(ico)          = 0.d0
               initp%gbw(ico)          = 0.d0
               cpatch%gbh(ico)         = 0.0
               cpatch%gbw(ico)         = 0.0
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the bulk resistance by integrating the inverse of the diffusivity for  !
         ! each layer (e.g. Sellers et al. (1986)).  We assumed Km = sigma u, which is the !
         ! preferred approach by Sellers et al. (1986).  To find sigma we used that        !
         ! tau = rho * ustar^2 at z=h, and that tau = rho * Km * dU/dz.                    !
         !---------------------------------------------------------------------------------!
         rasveg  = 0.d0
         sigmakm = ustarouh * ustarouh * cumldrag8(zcan) * pshelter8(zcan)                 &
                 / ( nn * cdrag8(zcan) * lad8(zcan) )
         do k=1,zcan
            windm978(k) = max( ugbmin8                                                     &
                             , uh * exp(- nn * (1.d0 - cumldrag8(k)/cumldrag8(zcan))) )
            Kdiff       = sigmakm * windm978(k) + kvwake8
            rasveg      = rasveg + dz_m978 / Kdiff
         end do
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !    According to Sellers et al (1986), equation A15, we must correct the         !
         ! resistance for the unstable case.  There is a typo on their equation (a gravity !
         ! term is missing, and g appears in their code in SiB-2.0 that is coupled to the  !
         ! original BRAMS-4.0).  Instead of their approximation that is ultimately based   !
         ! on Businger et al. (1971), we find an estimate of zeta using the same method    !
         ! that is applied to get zeta for the Beljaars-Holtslag method.  According to     !
         ! Sellers et al. (1986), we should apply this correction for the unstable case    !
         ! only.                                                                           !
         !---------------------------------------------------------------------------------!
         if (initp%ground_temp > initp%can_temp) then
            ribcan   = 2.d0 * grav8 * (initp%can_temp - initp%ground_temp)                 &
                     * htop * (1.d0 - d0ohgt - z0ohgt)                                     &
                     / ( (initp%can_temp + initp%ground_temp) * uh * uh)
            hgtoz0   = (1.d0 - d0ohgt) / z0ohgt
            lnhgtoz0 = log(hgtoz0)
            zetacan  = zoobukhov8(ribcan,htop * (1.d0 - d0ohgt),htop * z0ohgt              &
                                 ,hgtoz0,lnhgtoz0,hgtoz0,lnhgtoz0,.false.)
            phih     = sqrt(1.d0 - gamh8 * zetacan)
         else
            !----- Stable case, no correction needed. -------------------------------------!
            phih     = 1.d0
         end if
         initp%ggveg = phih / rasveg
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !     Find the net ground conductance.  The net conductance is derived from the   !
         ! net resistance, which is, in turn, the weighted average of the resistances in   !
         ! bare and vegetated grounds.                                                     !
         !---------------------------------------------------------------------------------!
         if (initp%opencan_frac > 9.99d-1 .or. initp%snowfac >= 9.d-1) then
            initp%ggnet = ggfact8 * initp%ggbare
         else
            initp%ggnet = ggfact8 * initp%ggbare * initp%ggveg                             &
                        / (initp%ggveg + (1.d0 - initp%opencan_frac) * initp%ggbare )
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap8(initp%can_rhos,initp%can_temp,initp%can_depth                     &
                        ,wcapcan,wcapcani,hcapcani,ccapcani)
         !---------------------------------------------------------------------------------!

      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine canopy_turbulence8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the characteristic scales based on surface layer          !
   ! parameterization.  Assume that stability is calculated on the potential  gradient     !
   ! from the surface to the atmosphere at reference height.  Apply the instability        !
   ! parameters to calculate the friction velocity.  Use the instability parameters for    !
   ! momentum and scalars, that were calculated over the distance from surface to refer-   !
   ! ence height to determine the heat, moisture and carbon flux rates at the canopy to    !
   ! atmosphere at reference height.                                                       !
   !    Two models are available, and the user can choose between them by setting the      !
   ! variable ISFCLYRM in ED2IN (if running the coupled model, this is done in ISTAR).     ! 
   !                                                                                       !
   ! 1. Based on L79;                                                                      !
   ! 2. Based on OD95, but not assuming that z>>z0, so the zeta0 terms shall be computed   !
   !    as presented in L79 and B71 to avoid singularities.                                !
   ! 3. Based on BH91, using an iterative method to find zeta, and using the modified      !
   !    equation for stable layers.                                                        !
   !                                                                                       !
   ! References:                                                                           !
   ! B71.  BUSINGER, J.A, et. al; Flux-Profile relationships in the atmospheric surface    !
   !           layer. J. Atmos. Sci., 28, 181-189, 1971.                                   !
   ! L79.  LOUIS, J.F.; Parametric Model of vertical eddy fluxes in the atmosphere.        !
   !           Boundary-Layer Meteor., 17, 187-202, 1979.                                  !
   ! BH91. BELJAARS, A.C.M.; HOLTSLAG, A.A.M.; Flux parameterization over land surfaces    !
   !           for atmospheric models. J. Appl. Meteor., 30, 327-341, 1991.                !
   ! OD95. ONCLEY, S.P.; DUDHIA, J.; Evaluation of surface fluxes from MM5 using observa-  !
   !           tions.  Mon. Wea. Rev., 123, 3344-3357, 1995.                               !
   !---------------------------------------------------------------------------------------!
   subroutine ed_stars(theta_atm,theiv_atm,shv_atm,co2_atm,theta_can,theiv_can             &
                      ,shv_can,co2_can,zref,dheight,uref,rough,ustar,tstar,estar,qstar     &
                      ,cstar,zeta,rib,ggbare)
      use consts_coms     , only : grav          & ! intent(in)
                                 , vonk          & ! intent(in)
                                 , epim1         & ! intent(in)
                                 , halfpi        ! ! intent(in)
      use canopy_air_coms , only : isfclyrm      & ! intent(in)
                                 , ustmin        & ! intent(in)
                                 , bl79          & ! intent(in)
                                 , csm           & ! intent(in)
                                 , csh           & ! intent(in)
                                 , dl79          & ! intent(in)
                                 , ribmaxod95    & ! intent(in)
                                 , ribmaxbh91    & ! intent(in)
                                 , tprandtl      & ! intent(in)
                                 , z0moz0h       & ! intent(in)
                                 , z0hoz0m       & ! intent(in)
                                 , psim          & ! function
                                 , psih          & ! function
                                 , zoobukhov     ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: theta_atm    ! Above canopy air pot. temperature     [        K]
      real, intent(in)  :: theiv_atm    ! Above canopy air eq. pot. temperature [        K]
      real, intent(in)  :: shv_atm      ! Above canopy vapour spec. hum.        [kg/kg_air]
      real, intent(in)  :: co2_atm      ! CO2 specific volume                   [  µmol/m³]
      real, intent(in)  :: theta_can    ! Canopy air potential temperature      [        K]
      real, intent(in)  :: theiv_can    ! Canopy air eq. pot. temperature       [        K]
      real, intent(in)  :: shv_can      ! Canopy air vapour spec. humidity      [kg/kg_air]
      real, intent(in)  :: co2_can      ! Canopy air CO2 specific volume        [  µmol/m³]
      real, intent(in)  :: zref         ! Height at reference point             [        m]
      real, intent(in)  :: dheight      ! Zero-plane displacement height        [        m]
      real, intent(in)  :: uref         ! Wind speed at reference height        [      m/s]
      real, intent(in)  :: rough        ! Roughness                             [        m]
      real, intent(out) :: ustar        ! U*, friction velocity                 [      m/s]
      real, intent(out) :: qstar        ! Specific humidity friction scale      [kg/kg_air]
      real, intent(out) :: tstar        ! Temperature friction scale            [        K]
      real, intent(out) :: estar        ! Equivalent pot. temp. scale           [        K]
      real, intent(out) :: cstar        ! CO2 spec. volume friction scale       [  µmol/m³]
      real, intent(out) :: zeta         ! z/(Obukhov length).                   [    -----]
      real, intent(out) :: rib          ! Bulk richardson number.               [    -----]
      real, intent(out) :: ggbare       ! Ground conductance                    [      m/s]
      !----- Local variables, used by L79. ------------------------------------------------!
      logical           :: stable       ! Stable state
      real              :: zoz0m        ! zref/rough(momentum)
      real              :: lnzoz0m      ! ln[zref/rough(momentum)]
      real              :: zoz0h        ! zref/rough(heat)
      real              :: lnzoz0h      ! ln[zref/rough(heat)]
      real              :: c3           ! coefficient to find the other stars
      !----- Local variables --------------------------------------------------------------!
      real              :: a2           ! Drag coefficient in neutral conditions
      real              :: c1           ! a2 * vels
      real              :: fm           ! Stability parameter for momentum
      real              :: fh           ! Stability parameter for heat
      real              :: c2           ! Part of the c coeff. common to momentum & heat.
      real              :: cm           ! c coefficient times |Rib|^1/2 for momentum.
      real              :: ch           ! c coefficient times |Rib|^1/2 for heat.
      real              :: ee           ! (z/z0)^1/3 -1. for eqn. 20 w/o assuming z/z0 >> 1.
      !----- Local variables, used by OD95. -----------------------------------------------!
      real              :: zeta0m       ! roughness(momentum)/(Obukhov length).
      real              :: zeta0h       ! roughness(heat)/(Obukhov length).
      !----- Aux. environment conditions. -------------------------------------------------!
      real              :: thetav_atm   ! Atmos. virtual potential temperature  [        K]
      real              :: thetav_can   ! Canopy air virtual pot. temperature   [        K]
      !----- External functions. ----------------------------------------------------------!
      real, external    :: cbrt         ! Cubic root
      !------------------------------------------------------------------------------------!

      !----- Finding the variables common to both methods. --------------------------------!
      thetav_atm = theta_atm * (1. + epim1 * shv_atm)
      thetav_can = theta_can * (1. + epim1 * shv_can)
      zoz0m      = (zref-dheight)/rough
      lnzoz0m    = log(zoz0m)
      zoz0h      = z0moz0h * zoz0m
      lnzoz0h    = log(zoz0h)
      rib        = 2.0 * grav * (zref-dheight-rough) * (thetav_atm-thetav_can)             &
                 / ( (thetav_atm+thetav_can) * uref * uref)
      stable     = thetav_atm >= thetav_can

      !------------------------------------------------------------------------------------!
      !     Here we find u* and the coefficient to find the other stars based on the       !
      ! chosen surface model.                                                              !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (1)

         !----- Compute the a-square factor and the coefficient to find theta*. -----------!
         a2   = vonk * vonk / (lnzoz0m * lnzoz0m)
         c1   = a2 * uref

         if (stable) then
            !----- Stable case ------------------------------------------------------------!
     
            fm = 1.0 / (1.0 + (2.0 * bl79 * rib / sqrt(1.0 + dl79 * rib)))
            fh = 1.0 / (1.0 + (3.0 * bl79 * rib * sqrt(1.0 + dl79 * rib)))

         else

            !------------------------------------------------------------------------------!
            !     Unstable case.  The only difference from the original method is that we  !
            ! no longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.    !
            !------------------------------------------------------------------------------!
            ee = cbrt(zoz0m) - 1.
            c2 = bl79 * a2 * ee * sqrt(ee * abs(rib))
            cm = csm * c2
            ch = csh * c2
            fm = (1.0 - 2.0 * bl79 * rib / (1.0 + 2.0 * cm))
            fh = (1.0 - 3.0 * bl79 * rib / (1.0 + 3.0 * ch))
         end if

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max(ustmin,sqrt(c1 * uref * fm))
         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3 = c1 * fh / ustar

         !---------------------------------------------------------------------------------!

         !----- Compute zeta from u* and T* -----------------------------------------------!
         zeta = grav * vonk * c3 * (thetav_atm - thetav_can) / (thetav_atm * ustar * ustar)


      case (2,4)
         !---------------------------------------------------------------------------------!
         ! 2. Here we use the model proposed by OD95, the standard for MM5, but with some  !
         !    terms that were computed in B71 (namely, the "0" terms). which prevent sin-  !
         !    gularities.  We use OD95 to estimate zeta, which avoids the computation      !
         !    of the Obukhov length L , we can't compute zeta0 by its definition(z0/L).    !
         !    However we know zeta, so zeta0 can be written as z0/z * zeta.                !
         ! 4. We use the model proposed by BH91, but we find zeta using the approximation  !
         !    given by OD95.                                                               !
         !---------------------------------------------------------------------------------!

         !----- Make sure that the bulk Richardson number is not above ribmax. ------------!
         rib = min(rib,ribmaxod95)
         
         !----- We now compute the stability correction functions. ------------------------!
         if (stable) then
            !----- Stable case. -----------------------------------------------------------!
            zeta  = rib * lnzoz0m / (1.1 - 5.0 * rib)
         else
            !----- Unstable case. ---------------------------------------------------------!
            zeta = rib * lnzoz0m
         end if
         zeta0m = rough * zeta / (zref - dheight)

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max (ustmin, vonk * uref                                                  &
                            / (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable)))

         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3    = vonk / (tprandtl * (lnzoz0m - psih(zeta,stable) + psih(zeta0m,stable)))

         !---------------------------------------------------------------------------------!


      case (3,5)
         !---------------------------------------------------------------------------------!
         ! 3. Here we use the model proposed by BH91, which is almost the same as the OD95 !
         !    method, with the two following (important) differences.                      !
         !    a. Zeta (z/L) is actually found using the iterative method.                  !
         !    b. Stable functions are computed in a more generic way.  BH91 claim that the !
         !       oft-used approximation (-beta*zeta) can cause poor ventilation of the     !
         !       stable layer, leading to decoupling between the atmosphere and the canopy !
         !       air space and excessive cooling.                                          !
         ! 5. Similar as 3, but we compute the stable functions the same way as OD95.      !
         !---------------------------------------------------------------------------------!

         !----- Make sure that the bulk Richardson number is not above ribmax. ------------!
         rib = min(rib,ribmaxbh91)
         
         !----- We now compute the stability correction functions. ------------------------!
         zeta   = zoobukhov(rib,zref-dheight,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
         zeta0m = rough * zeta / (zref-dheight)
         zeta0h = z0hoz0m * zeta0m

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max (ustmin, vonk * uref                                                  &
                            / (lnzoz0m - psim(zeta,stable) + psim(zeta0m,stable)))

         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3    = vonk / (tprandtl * (lnzoz0h - psih(zeta,stable) + psih(zeta0h,stable)))

         !---------------------------------------------------------------------------------!
      end select

      !----- Compute the other scales. ----------------------------------------------------!
      qstar = c3 *    (shv_atm   - shv_can   )
      tstar = c3 *    (theta_atm - theta_can )
      estar = c3 * log(theiv_atm / theiv_can )
      cstar = c3 *    (co2_atm   - co2_can   )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute the bare ground conductance.  This equation is similar to the original, !
      ! except that we don't assume the ratio between the gradient and the characteristic  !
      ! scale to be 0.2; instead we use the actual ratio that is computed here.            !
      !------------------------------------------------------------------------------------!
      ggbare  = c3 * ustar
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_stars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Double precision version of the above-described subroutine.  This subroutine       !
   ! computes the characteristic scales based on surface layer parameterization.  Assume   !
   ! that stability is calculated on the potential gradient from the surface to the        !
   ! atmosphere at reference height.  Apply the instability parameters to calculate the    !
   ! friction velocity.  Use the instability parameters for momentum and scalars, that     !
   ! were calculated over the distance from surface to reference height to determine the   !
   ! heat, moisture and carbon flux rates at the canopy to atmosphere at reference height. !
   !    Two models are available, and the user can choose between them by setting the      !
   ! variable ISFCLYRM in ED2IN (if running the coupled model, this is done in ISTAR).     ! 
   !                                                                                       !
   ! 1. Based on L79;                                                                      !
   ! 2. Based on: OD95, but with some terms computed as in L79 and B71 to avoid singular-  !
   !    ities.                                                                             !
   ! 3. Based on BH91, using an iterative method to find zeta, and using the modified      !
   !    equation for stable layers.                                                        !
   !                                                                                       !
   ! References:                                                                           !
   ! B71.  BUSINGER, J.A, et. al; Flux-Profile relationships in the atmospheric surface    !
   !           layer. J. Atmos. Sci., 28, 181-189, 1971.                                   !
   ! L79.  LOUIS, J.F.; Parametric Model of vertical eddy fluxes in the atmosphere.        !
   !           Boundary-Layer Meteor., 17, 187-202, 1979.                                  !
   ! BH91. BELJAARS, A.C.M.; HOLTSLAG, A.A.M.; Flux parameterization over land surfaces    !
   !           for atmospheric models. J. Appl. Meteor., 30, 327-341, 1991.                !
   ! OD95. ONCLEY, S.P.; DUDHIA, J.; Evaluation of surface fluxes from MM5 using observa-  !
   !           tions.  Mon. Wea. Rev., 123, 3344-3357, 1995.                               !
   !---------------------------------------------------------------------------------------!
   subroutine ed_stars8(theta_atm,theiv_atm,shv_atm,co2_atm                                &
                       ,theta_can,theiv_can,shv_can,co2_can                                &
                       ,zref,dheight,uref,rough,ustar,tstar,estar,qstar,cstar,zeta,rib     &
                       ,ggbare)
      use consts_coms     , only : grav8         & ! intent(in)
                                 , vonk8         & ! intent(in)
                                 , epim18        & ! intent(in)
                                 , halfpi8       ! ! intent(in)
      use canopy_air_coms , only : isfclyrm      & ! intent(in)
                                 , ustmin8       & ! intent(in)
                                 , bl798         & ! intent(in)
                                 , csm8          & ! intent(in)
                                 , csh8          & ! intent(in)
                                 , dl798         & ! intent(in)
                                 , ribmaxod958   & ! intent(in)
                                 , ribmaxbh918   & ! intent(in)
                                 , tprandtl8     & ! intent(in)
                                 , z0moz0h8      & ! intent(in)
                                 , z0hoz0m8      & ! intent(in)
                                 , psim8         & ! function
                                 , psih8         & ! function
                                 , zoobukhov8    ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: theta_atm    ! Above canopy air pot. temp.    [        K]
      real(kind=8), intent(in)  :: theiv_atm    ! Above canopy air eq. pot. T    [        K]
      real(kind=8), intent(in)  :: shv_atm      ! Above canopy vapour spec. hum. [kg/kg_air]
      real(kind=8), intent(in)  :: co2_atm      ! CO2 specific volume            [  µmol/m³]
      real(kind=8), intent(in)  :: theta_can    ! Canopy air potential temp.     [        K]
      real(kind=8), intent(in)  :: theiv_can    ! Canopy air eq. pot. temp.      [        K]
      real(kind=8), intent(in)  :: shv_can      ! Canopy air vapour spec. hum.    [kg/kg_air]
      real(kind=8), intent(in)  :: co2_can      ! Canopy air CO2 specific volume [  µmol/m³]
      real(kind=8), intent(in)  :: zref         ! Height at reference point      [        m]
      real(kind=8), intent(in)  :: dheight      ! Zero-plane displacement height [        m]
      real(kind=8), intent(in)  :: uref         ! Wind speed at reference height [      m/s]
      real(kind=8), intent(in)  :: rough        ! Roughness                      [        m]
      real(kind=8), intent(out) :: ustar        ! U*, friction velocity          [      m/s]
      real(kind=8), intent(out) :: qstar        ! Specific hum. friction scale   [kg/kg_air]
      real(kind=8), intent(out) :: tstar        ! Temperature friction scale     [        K]
      real(kind=8), intent(out) :: estar        ! Theta_E friction scale         [        K]
      real(kind=8), intent(out) :: cstar        ! CO2 spec. volume friction scale[  µmol/m³]
      real(kind=8), intent(out) :: zeta         ! z/(Obukhov length)             [      ---]
      real(kind=8), intent(out) :: rib          ! Bulk richardson number.        [      ---]
      real(kind=8), intent(out) :: ggbare       ! Ground conductance.            [      m/s]
      !----- Local variables, used by L79. ------------------------------------------------!
      logical                   :: stable       ! Stable state
      real(kind=8)              :: zoz0m        ! zref/rough(momentum)
      real(kind=8)              :: lnzoz0m      ! ln[zref/rough(momentum)]
      real(kind=8)              :: zoz0h        ! zref/rough(heat)
      real(kind=8)              :: lnzoz0h      ! ln[zref/rough(heat)]
      real(kind=8)              :: c3           ! coefficient to find the other stars
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)      :: a2           ! Drag coefficient in neutral conditions
      real(kind=8)      :: c1           ! a2 * vels
      real(kind=8)      :: fm           ! Stability parameter for momentum
      real(kind=8)      :: fh           ! Stability parameter for heat
      real(kind=8)      :: c2           ! Part of the c coeff. common to momentum & heat.
      real(kind=8)      :: cm           ! c coefficient times |Rib|^1/2 for momentum.
      real(kind=8)      :: ch           ! c coefficient times |Rib|^1/2 for heat.
      real(kind=8)      :: ee           ! (z/z0)^1/3 -1. for eqn. 20 w/o assuming z/z0 >> 1.
      !----- Local variables, used by OD95. -----------------------------------------------!
      real(kind=8)      :: zeta0m       ! roughness(momentum)/(Obukhov length).
      real(kind=8)      :: zeta0h       ! roughness(heat)/(Obukhov length).
      !----- Aux. environment conditions. -------------------------------------------------!
      real(kind=8)      :: thetav_atm   ! Atmos. virtual potential temperature  [        K]
      real(kind=8)      :: thetav_can   ! Canopy air virtual pot. temperature   [        K]
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external :: cbrt8   ! Cubic root
      !------------------------------------------------------------------------------------!

      !----- Finding the variables common to both methods. --------------------------------!
      thetav_atm = theta_atm * (1.d0 + epim18 * shv_atm)
      thetav_can = theta_can * (1.d0 + epim18 * shv_can)
      zoz0m      = (zref-dheight)/rough
      lnzoz0m    = log(zoz0m)
      zoz0h      = z0moz0h8 * zoz0m
      lnzoz0h    = log(zoz0h)
      rib        = 2.d0 * grav8 * (zref-dheight-rough) * (thetav_atm-thetav_can)           &
                 / ( (thetav_atm+thetav_can) * uref * uref)
      stable     = thetav_atm >= thetav_can

      !------------------------------------------------------------------------------------!
      !     Here we find u* and the coefficient to find the other stars based on the       !
      ! chosen surface model.                                                              !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (1)

         !----- Compute the a-square factor and the coefficient to find theta*. -----------!
         a2   = vonk8 * vonk8 / (lnzoz0m * lnzoz0m)
         c1   = a2 * uref

         if (stable) then
            !----- Stable case ------------------------------------------------------------!
     
            fm = 1.d0 / (1.d0 + (2.d0 * bl798 * rib / sqrt(1.d0 + dl798 * rib)))
            fh = 1.d0 / (1.d0 + (3.d0 * bl798 * rib * sqrt(1.d0 + dl798 * rib)))

         else

            !------------------------------------------------------------------------------!
            !     Unstable case.  The only difference from the original method is that we  !
            ! no longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.    !
            !------------------------------------------------------------------------------!
            ee = cbrt8(zoz0m) - 1.d0
            c2 = bl798 * a2 * ee * sqrt(ee * abs(rib))
            cm = csm8 * c2
            ch = csh8 * c2
            fm = (1.d0 - 2.d0 * bl798 * rib / (1.d0 + 2.d0 * cm))
            fh = (1.d0 - 3.d0 * bl798 * rib / (1.d0 + 3.d0 * ch))
         end if

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max(ustmin8,sqrt(c1 * uref * fm))
         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3 = c1 * fh / ustar
         !---------------------------------------------------------------------------------!

         !----- Compute zeta from u* and T* -----------------------------------------------!
         zeta = grav8 * vonk8 * c3 * (thetav_atm - thetav_can)                             &
              / (thetav_atm * ustar * ustar)


      case (2,4)
         !---------------------------------------------------------------------------------!
         ! 2. Here we use the model proposed by OD95, the standard for MM5, but with some  !
         !    terms that were computed in B71 (namely, the "0" terms). which prevent sin-  !
         !    gularities.  We use OD95 to estimate zeta, which avoids the computation      !
         !    of the Obukhov length L , we can't compute zeta0 by its definition(z0/L).    !
         !    However we know zeta, so zeta0 can be written as z0/z * zeta.                !
         ! 4. We use the model proposed by BH91, but we find zeta using the approximation  !
         !    given by OD95.                                                               !
         !---------------------------------------------------------------------------------!

         !----- Make sure that the bulk Richardson number is not above ribmax. ------------!
         rib = min(rib,ribmaxod958)
         
         !----- We now compute the stability correction functions. ------------------------!
         if (stable) then
            !----- Stable case. -----------------------------------------------------------!
            zeta  = rib * lnzoz0m / (1.1d0 - 5.0d0 * rib)
         else
            !----- Unstable case. ---------------------------------------------------------!
            zeta  = rib * lnzoz0m
         end if

         zeta0m = rough * zeta / (zref - dheight)

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max (ustmin8, vonk8 * uref                                                &
                             / (lnzoz0m - psim8(zeta,stable) + psim8(zeta0m,stable)))

         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3    = vonk8 / (tprandtl8 * (lnzoz0m - psih8(zeta,stable) + psih8(zeta0m,stable)))

         !---------------------------------------------------------------------------------!


      case (3,5)
         !---------------------------------------------------------------------------------!
         ! 3. Here we use the model proposed by BH91, which is almost the same as the OD95 !
         !    method, with the two following (important) differences.                      !
         !    a. Zeta (z/L) is actually found using the iterative method.                  !
         !    b. Stable functions are computed in a more generic way.  BH91 claim that the !
         !       oft-used approximation (-beta*zeta) can cause poor ventilation of the     !
         !       stable layer, leading to decoupling between the atmosphere and the canopy !
         !       air space and excessive cooling.                                          !
         ! 5. Similar as 3, but we compute the stable functions the same way as OD95.      !
         !---------------------------------------------------------------------------------!

         !----- Make sure that the bulk Richardson number is not above ribmax. ------------!
         rib = min(rib,ribmaxbh918)

         !----- We now compute the stability correction functions. ------------------------!
         zeta   = zoobukhov8(rib,zref-dheight,rough,zoz0m,lnzoz0m,zoz0h,lnzoz0h,stable)
         zeta0m = rough * zeta / (zref - dheight)
         zeta0h = z0hoz0m8 * zeta0m

         !----- Finding ustar, making sure it is not too small. ---------------------------!
         ustar = max (ustmin8, vonk8 * uref                                                &
                             / (lnzoz0m - psim8(zeta,stable) + psim8(zeta0m,stable)))

         !----- Finding the coefficient to scale the other stars. -------------------------!
         c3    = vonk8                                                                     &
               / (tprandtl8 * (lnzoz0h - psih8(zeta,stable) + psih8(zeta0h,stable)))

         !---------------------------------------------------------------------------------!
      end select

      !----- Computing the other scales. --------------------------------------------------!
      qstar = c3 *    (shv_atm   - shv_can   )
      tstar = c3 *    (theta_atm - theta_can )
      estar = c3 * log(theiv_atm / theiv_can )
      cstar = c3 *    (co2_atm   - co2_can   )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute the bare ground conductance.  This equation is similar to the original, !
      ! except that we don't assume the ratio between the gradient and the characteristic  !
      ! scale to be 0.2; instead we use the actual ratio that is computed here.            !
      !------------------------------------------------------------------------------------!
      ggbare = c3 * ustar
      !------------------------------------------------------------------------------------!

      return
   end subroutine ed_stars8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function determines the wind at a given height, given that the stars are al-  !
   ! ready known, as well as the Richardson number and the zetas.                          !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function reduced_wind(ustar,zeta,rib,zref,dheight,height,rough)
      use consts_coms    , only : vonk     ! ! intent(in)
      use canopy_air_coms, only : isfclyrm & ! intent(in)
                                , bl79     & ! intent(in)
                                , csm      & ! intent(in)
                                , csh      & ! intent(in)
                                , dl79     & ! intent(in)
                                , ugbmin   & ! intent(in)
                                , psim     ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: ustar     ! Friction velocity                   [    m/s]
      real(kind=4), intent(in) :: zeta      ! Normalised height                   [    ---]
      real(kind=4), intent(in) :: rib       ! Bulk Richardson number              [    ---]
      real(kind=4), intent(in) :: zref      ! Reference height                    [      m]
      real(kind=4), intent(in) :: dheight   ! Displacement height                 [      m]
      real(kind=4), intent(in) :: height    ! Height to determine the red. wind   [      m]
      real(kind=4), intent(in) :: rough     ! Roughness scale                     [      m]
      !----- Local variables. -------------------------------------------------------------!
      logical                  :: stable    ! Canopy air space is stable          [    T|F]
      real(kind=4)             :: zetah     ! Zeta for h=height                   [    ---]
      real(kind=4)             :: zeta0     ! Zeta for h=rough                    [    ---]
      real(kind=4)             :: hoz0      ! ((h-dheight)/z0)                    [    ---]
      real(kind=4)             :: lnhoz0    ! ln ((h-dheight)/z0)                 [    ---]
      real(kind=4)             :: a2        ! Drag coeff. in neutral conditions
      real(kind=4)             :: fm        ! Stability parameter for momentum
      real(kind=4)             :: c2        ! Part of the c coefficient.
      real(kind=4)             :: cm        ! c coefficient times |Rib|^1/2
      real(kind=4)             :: ee        ! (z/z0)^1/3 -1. for eqn. 20 (L79)
      !----- External functions. ----------------------------------------------------------!
      real(kind=4), external   :: cbrt      ! Cubic root
      !------------------------------------------------------------------------------------!



      !----- Define whether the layer is stable or not. -----------------------------------!
      stable    = rib >= 0.
      !------------------------------------------------------------------------------------!



      !----- Find the log for the log-height interpolation of wind. -----------------------!
      hoz0      = (height - dheight)/rough
      lnhoz0    = log(hoz0)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The wind at a given height is found by using the same definition of wind speed !
      ! at a given height.                                                                 !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (1) !---- Louis (1979) method. -------------------------------------------------!

         !----- Compute the a-square factor and the coefficient to find theta*. -----------!
         a2   = vonk * vonk / (lnhoz0 * lnhoz0)

         if (stable) then
            !----- Stable case ------------------------------------------------------------!
            fm = 1.0 / (1.0 + (2.0 * bl79 * rib / sqrt(1.0 + dl79 * rib)))

         else
            !------------------------------------------------------------------------------!
            !     Unstable case.  The only difference from the original method is that we  !
            ! no longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.    !
            !------------------------------------------------------------------------------!
            ee = cbrt(hoz0) - 1.
            c2 = bl79 * a2 * ee * sqrt(ee * abs(rib))
            cm = csm * c2
            fm = (1.0 - 2.0 * bl79 * rib / (1.0 + 2.0 * cm))
         end if
         
         !----- Find the wind. ------------------------------------------------------------!
         reduced_wind = (ustar/vonk) * (lnhoz0/sqrt(fm))

      case default  !----- Other methods. -------------------------------------------------!

         !----- Determine zeta for the sought height and for roughness height. ------------!
         zetah = zeta * (height-dheight) / (zref-dheight)
         zeta0 = zeta * rough            / (zref-dheight)
         !---------------------------------------------------------------------------------!

         reduced_wind = (ustar/vonk) * (lnhoz0 - psim(zetah,stable) + psim(zeta0,stable))

      end select
      !------------------------------------------------------------------------------------!



      !----- Impose the wind to be more than the minimum. ---------------------------------!
      reduced_wind = max(reduced_wind, ugbmin)
      !------------------------------------------------------------------------------------!


      return
   end function reduced_wind
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function determines the wind at a given height, given that the stars are al-  !
   ! ready known, as well as the Richardson number and the zetas.                          !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function reduced_wind8(ustar,zeta,rib,zref,dheight,height,rough)
      use consts_coms    , only : vonk8     ! ! intent(in)
      use canopy_air_coms, only : isfclyrm  & ! intent(in)
                                , bl798     & ! intent(in)
                                , csm8      & ! intent(in)
                                , csh8      & ! intent(in)
                                , dl798     & ! intent(in)
                                , ugbmin8   & ! intent(in)
                                , psim8     ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: ustar     ! Friction velocity                   [    m/s]
      real(kind=8), intent(in) :: zeta      ! Normalised height                   [    ---]
      real(kind=8), intent(in) :: rib       ! Bulk Richardson number              [    ---]
      real(kind=8), intent(in) :: zref      ! Reference height                    [      m]
      real(kind=8), intent(in) :: dheight   ! Displacement height                 [      m]
      real(kind=8), intent(in) :: height    ! Height to determine the red. wind   [      m]
      real(kind=8), intent(in) :: rough     ! Roughness scale                     [      m]
      !----- Local variables. -------------------------------------------------------------!
      logical                  :: stable    ! Canopy air space is stable          [    T|F]
      real(kind=8)             :: zetah     ! Zeta for h=height                   [    ---]
      real(kind=8)             :: zeta0     ! Zeta for h=rough                    [    ---]
      real(kind=8)             :: hoz0      ! ((h-d0)/z0)                         [    ---]
      real(kind=8)             :: lnhoz0    ! ln ((h-d0)/z0)                      [    ---]
      real(kind=8)             :: a2        ! Drag coeff. in neutral conditions
      real(kind=8)             :: fm        ! Stability parameter for momentum
      real(kind=8)             :: c2        ! Part of the c coefficient.
      real(kind=8)             :: cm        ! c coefficient times |Rib|^1/2
      real(kind=8)             :: ee        ! (z/z0)^1/3 -1. for eqn. 20 (L79)
      !----- External functions. ----------------------------------------------------------!
      real(kind=8), external   :: cbrt8     ! Cubic root
      !------------------------------------------------------------------------------------!



      !----- Define whether the layer is stable or not. -----------------------------------!
      stable    = rib >= 0.d0
      !------------------------------------------------------------------------------------!



      !----- Find the log for the log-height interpolation of wind. -----------------------!
      hoz0      = (height-dheight)/rough
      lnhoz0    = log(hoz0)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The wind at a given height is found by using the same definition of wind speed !
      ! at a given height.                                                                 !
      !------------------------------------------------------------------------------------!
      select case (isfclyrm)
      case (1) !---- Louis (1979) method. -------------------------------------------------!

         !----- Compute the a-square factor and the coefficient to find theta*. -----------!
         a2   = vonk8 * vonk8 / (lnhoz0 * lnhoz0)

         if (stable) then
            !----- Stable case ------------------------------------------------------------!
            fm = 1.d0 / (1.d0 + (2.d0 * bl798 * rib / sqrt(1.d0 + dl798 * rib)))

         else
            !------------------------------------------------------------------------------!
            !     Unstable case.  The only difference from the original method is that we  !
            ! no longer assume z >> z0, so the "c" coefficient uses the full z/z0 term.    !
            !------------------------------------------------------------------------------!
            ee = cbrt8(hoz0) - 1.d0
            c2 = bl798 * a2 * ee * sqrt(ee * abs(rib))
            cm = csm8 * c2
            fm = (1.d0 - 2.d0 * bl798 * rib / (1.d0 + 2.d0 * cm))
         end if
         
         !----- Find the wind. ------------------------------------------------------------!
         reduced_wind8 = (ustar/vonk8) * (lnhoz0/sqrt(fm))

      case default  !----- Other methods. -------------------------------------------------!

         !----- Determine zeta for the sought height and for roughness height. ------------!
         zetah = zeta * (height-dheight) / (zref-dheight)
         zeta0 = zeta * rough            / (zref-dheight)
         !---------------------------------------------------------------------------------!

         reduced_wind8 = (ustar/vonk8)                                                     &
                       * (lnhoz0 - psim8(zetah,stable) + psim8(zeta0,stable))

      end select
      !------------------------------------------------------------------------------------!



      !----- Impose the wind to be more than the minimum. ---------------------------------!
      reduced_wind8 = max(reduced_wind8,ugbmin8)
      !------------------------------------------------------------------------------------!

      return
   end function reduced_wind8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function vertical_vel_flux(zeta,tstar,ustar)
      use consts_coms , only : vonk ! intent(in)
     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)    :: zeta
      real, intent(in)    :: ustar
      real, intent(in)    :: tstar
      !----- Local variables --------------------------------------------------------------!
      real                :: cx
      real                :: psin
      !----- Constants --------------------------------------------------------------------!
      real, parameter     :: wtol = 1.e-20
      !------------------------------------------------------------------------------------!

      if (zeta < 0.0)then
         cx = zeta * sqrt(sqrt(1.0 - 15.0 * zeta))
      else
         cx = zeta / (1.0 + 4.7 * zeta)
      endif
     
      psin = sqrt((1.0-2.86 * cx) / (1.0 + cx * (-5.390 + cx * 6.9980 )))
      vertical_vel_flux = ( 0.27 * max(6.25 * (1.0 - cx) * psin,wtol)                      &
                          - 1.180 * cx * psin) * ustar * ustar
     
      return
   end function vertical_vel_flux
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real(kind=8) function vertical_vel_flux8(zeta,tstar,ustar)
      use consts_coms , only : vonk8 ! intent(in)
     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)    :: zeta
      real(kind=8), intent(in)    :: ustar
      real(kind=8), intent(in)    :: tstar
      !----- Local variables --------------------------------------------------------------!
      real(kind=8) :: cx
      real(kind=8) :: psin
      !----- Constants --------------------------------------------------------------------!
      real(kind=8), parameter     :: wtol = 1.d-20
      !------------------------------------------------------------------------------------!
     
     
      if (zeta < 0.d0)then
         cx = zeta * sqrt(sqrt(1.d0 - 1.5d1 * zeta))
      else
         cx = zeta / (1.0d0 + 4.7d0 * zeta)
      endif
     
      psin = sqrt((1.d0-2.86d0 * cx) / (1.d0 + cx * (-5.39d0 + cx * 6.998d0 )))
      vertical_vel_flux8 = ( 2.7d-1 * max(6.25d0 * (1.d0 - cx) * psin,wtol)                &
                           - 1.18d0 * cx * psin) * ustar * ustar
     
      return
   end function vertical_vel_flux8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Calculate some canopy air space properties, such as the total mass and depth, and !
   ! also the total capacities (carbon, water, and heat).                                  !
   !---------------------------------------------------------------------------------------!
   subroutine can_whcap(can_rhos,can_temp,can_depth,wcapcan,wcapcani,hcapcani,ccapcani)
      use consts_coms, only : cpi    & ! intent(in)
                            , mmdry  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)  :: can_rhos
      real(kind=4), intent(in)  :: can_temp
      real(kind=4), intent(in)  :: can_depth
      real(kind=4), intent(out) :: wcapcan
      real(kind=4), intent(out) :: wcapcani
      real(kind=4), intent(out) :: hcapcani
      real(kind=4), intent(out) :: ccapcani
      !------------------------------------------------------------------------------------!

      wcapcan  = can_rhos * can_depth
      wcapcani = 1.0 / wcapcan
      hcapcani = cpi * wcapcani / can_temp
      ccapcani = mmdry * wcapcani

      return
   end subroutine can_whcap
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     Calculate some canopy air space properties, such as the total mass and depth, and !
   ! also the total capacities (carbon, water, and heat).                                  !
   !     Canopy air space depth must change because we are imposing that pressure and      !
   ! density are constants. Since mass and temperature can freely change, in order to      !
   ! guarantee that the ideal gas law and the first law of thermodynamic apply, the volume !
   ! (or the canopy depth) must be allowed to change over time, so work can be done by the !
   ! canopy or into the canopy.                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine can_whcap8(can_rhos,can_temp,can_depth,wcapcan,wcapcani,hcapcani,ccapcani)
      use consts_coms, only : cpi8    & ! intent(in)
                            , rdry8   & ! intent(in)
                            , ep8     & ! intent(in)
                            , mmdry8  ! ! intent(in)
      
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: can_rhos
      real(kind=8), intent(in)  :: can_temp
      real(kind=8), intent(in)  :: can_depth
      real(kind=8), intent(out) :: wcapcan
      real(kind=8), intent(out) :: wcapcani
      real(kind=8), intent(out) :: hcapcani
      real(kind=8), intent(out) :: ccapcani
      !------------------------------------------------------------------------------------!

      wcapcan  = can_rhos * can_depth
      wcapcani = 1.d0 / wcapcan
      hcapcani = cpi8 * wcapcani / can_temp
      ccapcani = mmdry8 * wcapcani
      return
   end subroutine can_whcap8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine computes the aerodynamic conductance between leaf and canopy     !
   ! air space for both heat and water vapour, based on:                                   !
   !                                                                                       !
   ! L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf        !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !                                                                                       !
   ! Notice that the units are somewhat different from L95.                                !
   ! - gbh is in J/(K m2 s), and                                                           !
   ! - gbw is in kg_H2O/m2/s.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine aerodynamic_conductances(ipft,veg_wind,veg_temp,can_temp,can_shv,can_rhos    &
                                      ,gbhmos_min,gbh,gbw)
      use pft_coms       , only : leaf_width ! ! intent(in)
      use canopy_air_coms, only : aflat_turb & ! intent(in)
                                , aflat_lami & ! intent(in)
                                , nflat_turb & ! intent(in)
                                , nflat_lami & ! intent(in)
                                , bflat_turb & ! intent(in)
                                , bflat_lami & ! intent(in)
                                , mflat_turb & ! intent(in)
                                , mflat_lami & ! intent(in)
                                , beta_r1    & ! intent(in)
                                , beta_r2    & ! intent(in)
                                , beta_re0   & ! intent(in)
                                , beta_g1    & ! intent(in)
                                , beta_g2    & ! intent(in)
                                , beta_gr0   ! ! intent(in)
      use consts_coms    , only : gr_coeff   & ! intent(in)
                                , th_diffi   & ! intent(in)
                                , th_diff    & ! intent(in)
                                , cp         ! ! intent(in)
      use physiology_coms, only : gbh_2_gbw  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                      :: ipft            ! Plant functional type   [      ---]
      real(kind=4)   , intent(in)  :: veg_wind        ! Wind at cohort height   [      m/s]
      real(kind=4)   , intent(in)  :: veg_temp        ! Leaf temperature        [        K]
      real(kind=4)   , intent(in)  :: can_temp        ! Canopy air temperature  [        K]
      real(kind=4)   , intent(in)  :: can_shv         ! Canopy air spec. hum.   [    kg/kg]
      real(kind=4)   , intent(in)  :: can_rhos        ! Canopy air density      [    kg/m³]
      real(kind=4)   , intent(in)  :: gbhmos_min      ! Min. Heat  conductance  [      m/s]
      real(kind=4)   , intent(out) :: gbh             ! Heat  conductance       [ J/K/m²/s]
      real(kind=4)   , intent(out) :: gbw             ! Water conductance       [  kg/m²/s]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                 :: lwidth          ! Leaf width              [        m]
      real(kind=4)                 :: grashof         ! Grashof number          [      ---]
      real(kind=4)                 :: reynolds        ! Reynolds number         [      ---]
      real(kind=4)                 :: nusselt_lami    ! Nusselt number (laminar)[      ---]
      real(kind=4)                 :: nusselt_turb    ! Nusselt number (turb.)  [      ---]
      real(kind=4)                 :: nusselt         ! Nusselt number          [      ---]
      real(kind=4)                 :: beta_forced     ! Correct.  term (forced) [      ---]
      real(kind=4)                 :: beta_free       ! Correct.  term (free)   [      ---]
      real(kind=4)                 :: forced_gbh_mos  ! Forced convection cond. [      m/s]
      real(kind=4)                 :: free_gbh_mos    ! Free convection cond.   [      m/s]
      real(kind=4)                 :: gbh_mos         ! Total convection cond.  [      m/s]
      !------------------------------------------------------------------------------------!


      !----- Save the leaf width of this PFT. ---------------------------------------------!
      lwidth = leaf_width(ipft)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s, associated with forced convection.               !
      !------------------------------------------------------------------------------------!
      !----- 1. Compute the Reynolds number. ----------------------------------------------!
      reynolds        = veg_wind * lwidth * th_diffi
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = aflat_lami * reynolds ** nflat_lami
      nusselt_turb    = aflat_turb * reynolds ** nflat_turb
      !----- 3. Compute the correction term for the theoretical Nusselt numbers. ----------!
      beta_forced     = beta_r1 + beta_r2 * tanh(log(reynolds/beta_re0))
      !----- 4. The right Nusselt number is the largest. ----------------------------------!
      nusselt         = beta_forced * max(nusselt_lami,nusselt_turb)
      !----- 5. The conductance is given by MU08 - equation 10.4 --------------------------!
      forced_gbh_mos  = th_diff * nusselt / lwidth
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s,  associated with free convection.                !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the Grashof number. --------------------------------------------------!
      grashof         = gr_coeff  * abs(veg_temp - can_temp) * lwidth * lwidth * lwidth
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = bflat_lami * grashof ** mflat_lami
      nusselt_turb    = bflat_turb * grashof ** mflat_turb
      !----- 3. Compute the correction term for the theoretical Nusselt numbers. ----------!
      if (grashof == 0.0) then
         beta_free    = beta_g1 - beta_g2
      else
         beta_free    = beta_g1 + beta_g2 * tanh(log(grashof/beta_gr0))
      end if
      !----- 4. The right Nusselt number is the largest. ----------------------------------!
      nusselt         = beta_free * max(nusselt_lami,nusselt_turb)
      !----- 5. The conductance is given by MU08 - equation 10.4 --------------------------!
      free_gbh_mos    = th_diff * nusselt / lwidth
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     The heat conductance for the thermodynamic budget is the sum of conductances,  !
      ! because we assume both forms of convection happen parallelly.  The conversion from !
      ! heat to water conductance (in m/s) can be found in L95, page 1198, after equation  !
      ! E5.  For the ED purposes, the output variables are converted to the units of       !
      ! entropy and water fluxes [J/K/m²/s and kg/m²/s, respectively].                     !
      !------------------------------------------------------------------------------------!
      gbh_mos = max(gbhmos_min, free_gbh_mos + forced_gbh_mos)
      gbh     =             gbh_mos * can_rhos * cp
      gbw     = gbh_2_gbw * gbh_mos * can_rhos
      !------------------------------------------------------------------------------------!

      return
   end subroutine aerodynamic_conductances
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine computes the aerodynamic conductance between leaf and canopy     !
   ! air space for both heat and water vapour, based on:                                   !
   !                                                                                       !
   ! L95 - Leuning, R., F. M. Kelliher, D. G. G. de Pury, E. D. Schulze, 1995: Leaf        !
   !       nitrogen, photosynthesis, conductance and transpiration: scaling from leaves to !
   !       canopies.  Plant, Cell and Environ., 18, 1183-1200.                             !
   ! M08 - Monteith, J. L., M. H. Unsworth, 2008. Principles of Environmental Physics,     !
   !       3rd. edition, Academic Press, Amsterdam, 418pp.  (Mostly Chapter 10).           !
   !                                                                                       !
   ! Notice that the units are somewhat different from L95.                                !
   ! - gbh is in J/(K m2 s), and                                                           !
   ! - gbw is in kg_H2O/m2/s.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine aerodynamic_conductances8(ipft,veg_wind,veg_temp,can_temp,can_shv,can_rhos   &
                                       ,gbhmos_min,gbh,gbw,reynolds,grashof,nusselt_free   &
                                       ,nusselt_forced)
      use pft_coms       , only : leaf_width  ! ! intent(in)
      use canopy_air_coms, only : aflat_turb8 & ! intent(in)
                                , aflat_lami8 & ! intent(in)
                                , nflat_turb8 & ! intent(in)
                                , nflat_lami8 & ! intent(in)
                                , bflat_turb8 & ! intent(in)
                                , bflat_lami8 & ! intent(in)
                                , mflat_turb8 & ! intent(in)
                                , mflat_lami8 & ! intent(in)
                                , beta_lami8  & ! intent(in)
                                , beta_turb8  & ! intent(in)
                                , beta_r18    & ! intent(in)
                                , beta_r28    & ! intent(in)
                                , beta_re08   & ! intent(in)
                                , beta_g18    & ! intent(in)
                                , beta_g28    & ! intent(in)
                                , beta_gr08   ! ! intent(in)
      use consts_coms    , only : gr_coeff8   & ! intent(in)
                                , th_diffi8   & ! intent(in)
                                , th_diff8    & ! intent(in)
                                , cp8         ! ! intent(in)
      use physiology_coms, only : gbh_2_gbw8  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                      :: ipft            ! Plant functional type   [      ---]
      real(kind=8)   , intent(in)  :: veg_wind        ! Wind at cohort height   [      m/s]
      real(kind=8)   , intent(in)  :: veg_temp        ! Leaf temperature        [        K]
      real(kind=8)   , intent(in)  :: can_temp        ! Canopy air temperature  [        K]
      real(kind=8)   , intent(in)  :: can_shv         ! Canopy air spec. hum.   [    kg/kg]
      real(kind=8)   , intent(in)  :: can_rhos        ! Canopy air density      [    kg/m³]
      real(kind=8)   , intent(in)  :: gbhmos_min      ! Min. heat  conductance  [      m/s]
      real(kind=8)   , intent(out) :: gbh             ! Heat  conductance       [ J/K/m²/s]
      real(kind=8)   , intent(out) :: gbw             ! Water conductance       [  kg/m²/s]
      real(kind=8)   , intent(out) :: grashof         ! Grashof number          [      ---]
      real(kind=8)   , intent(out) :: reynolds        ! Reynolds number         [      ---]
      real(kind=8)   , intent(out) :: nusselt_free    ! Nusselt number (free)   [      ---]
      real(kind=8)   , intent(out) :: nusselt_forced  ! Nusselt number (forced) [      ---]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                 :: lwidth          ! Leaf width              [        m]
      real(kind=8)                 :: nusselt_lami    ! Nusselt number (laminar)[      ---]
      real(kind=8)                 :: nusselt_turb    ! Nusselt number (turb.)  [      ---]
      real(kind=8)                 :: beta_forced     ! Correct.  term (forced) [      ---]
      real(kind=8)                 :: beta_free       ! Correct.  term (free)   [      ---]
      real(kind=8)                 :: forced_gbh_mos  ! Forced convection cond. [      m/s]
      real(kind=8)                 :: free_gbh_mos    ! Free convection cond.   [      m/s]
      real(kind=8)                 :: gbh_mos         ! Total convection cond.  [      m/s]
      !------------------------------------------------------------------------------------!


      !----- Save the leaf width of this PFT. ---------------------------------------------!
      lwidth = dble(leaf_width(ipft))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s, associated with forced convection.               !
      !------------------------------------------------------------------------------------!
      !----- 1. Compute the Reynolds number. ----------------------------------------------!
      reynolds        = veg_wind * lwidth * th_diffi8
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = aflat_lami8 * reynolds ** nflat_lami8
      nusselt_turb    = aflat_turb8 * reynolds ** nflat_turb8
      !----- 3. Compute the correction term for the theoretical Nusselt numbers. ----------!
      beta_forced     = beta_r18 + beta_r28 * tanh(log(reynolds/beta_re08))
      !----- 4. The right Nusselt number is the largest of the both. ----------------------!
      nusselt_forced  = beta_forced * max(nusselt_lami,nusselt_turb)
      !----- 5. The conductance is given by MU08 - equation 10.4 --------------------------!
      forced_gbh_mos  = th_diff8 * nusselt_forced / lwidth
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the conductance, in m/s,  associated with free convection.                !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the Grashof number. --------------------------------------------------!
      grashof         = gr_coeff8 * abs(veg_temp - can_temp) * lwidth * lwidth * lwidth
      !----- 2. Compute the Nusselt number for both the laminar and turbulent case. -------!
      nusselt_lami    = bflat_lami8 * grashof ** mflat_lami8
      nusselt_turb    = bflat_turb8 * grashof ** mflat_turb8
      !----- 3. Compute the correction term for the theoretical Nusselt numbers. ----------!
      if (grashof == 0.d0) then
         beta_free    = beta_g18 - beta_g28
      else
         beta_free    = beta_g18 + beta_g28 * tanh(log(grashof/beta_gr08))
      end if
      !----- 4. The right Nusselt number is the largest of the both. ----------------------!
      nusselt_free    = beta_free * max(nusselt_lami,nusselt_turb)
      !----- 5. The conductance is given by MU08 - equation 10.4 --------------------------!
      free_gbh_mos    = th_diff8 * nusselt_free / lwidth
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The heat conductance for the thermodynamic budget is the sum of conductances,  !
      ! because we assume both forms of convection happen parallelly.  The conversion from !
      ! heat to water conductance (in m/s) can be found in L95, page 1198, after equation  !
      ! E5.  For the ED purposes, the output variables are converted to the units of       !
      ! entropy and water fluxes [J/K/m²/s and kg/m²/s, respectively].                     !
      !------------------------------------------------------------------------------------!
      gbh_mos = max(gbhmos_min, free_gbh_mos + forced_gbh_mos)
      gbh     =              gbh_mos * can_rhos * cp8
      gbw     = gbh_2_gbw8 * gbh_mos * can_rhos
      !------------------------------------------------------------------------------------!

      return
   end subroutine aerodynamic_conductances8
   !=======================================================================================!
   !=======================================================================================!
end module canopy_struct_dynamics
!==========================================================================================!
!==========================================================================================!
