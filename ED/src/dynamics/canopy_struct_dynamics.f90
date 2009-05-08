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
   ! is due to the assumption of constant shear (ustar) in the dynamic sublayer.  An       !
   ! appropriate closure the MOST equations uses a reference wind velocity well above the  !
   ! top of the rough wall boundary.  However sometimes we do not have wind speed at the   !
   ! correct reference height, which brings us to the first two variants of this canopy    !
   ! turbulence scheme.                                                                    !
   !                                                                                       !
   ! 0. Uses the default MOST method, where displacement height and rougness are linear    !
   !    functions of mean canopy height, where mean canopy height is defined by the        !
   !    dead-biomass weighted average of canopy height.  Wind speed and diffusivity decay  !
   !    exponentially in the canopy on parameter alpha.  HOWEVER, if the windspeed at      !
   !    reference height is less than the canopy height (WHICH IS REALLY INNAPROPRIATE TO  !
   !    BEGIN WITH, BUT OH WELL), then revert to a calculation of ustar that uses no       !
   !    displacement height.  After ustar is calculated, then use displacement height in   !
   !    the calculation of the eddy length scale at the canopy top to find diffusivity at  !
   !    canopy top. THIS IS THE METHOD USED IN LEAF3 AND ED2.0.                            !
   !                                                                                       !
   ! 1. This is the same as method 0, with the following difference.  If the wind-speed at !
   !    reference height is less than the canopy height, first, calculate wind speed at    !
   !    canopy top using the exponential decay function (reversely).  Then solve for ustar !
   !    with this wind speed, the canopy top as reference height, and the default zero     !
   !    plane displacement height and roughness length.                                    !
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
   !     Ultimately, this routine solves for the resistance of water vapor and sensible    !
   ! heat from the soil surface to canopy air space and from the leaf surfaces to canopy   !
   ! air space.                                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_turbulence(csite,initp,isi,ipa,get_flow_geom)
      use ed_misc_coms   , only : icanturb             ! ! intent(in), can. turb. scheme
      use ed_state_vars  , only : polygontype          & ! structure
                                , sitetype             & ! structure
                                , patchtype            ! ! structure
      use rk4_coms       , only : rk4patchtype         & ! structure
                                , rk4met               & ! intent(in)
                                , tiny_offset          ! ! intent(in)
      use pft_coms       , only : crown_depth_fraction & ! intent(in)
                                , leaf_width           ! ! intent(in)
      use canopy_air_coms, only : exar8                & ! intent(in)
                                , dz                   & ! intent(in)
                                , mu0                  & ! intent(in)
                                , mu0_c                & ! intent(in)
                                , Cd0                  & ! intent(in)
                                , Pm                   & ! intent(in)
                                , c1_m97               & ! intent(in)
                                , c2_m97               & ! intent(in)
                                , c3_m97               & ! intent(in)
                                , kvwake               & ! intent(in)
                                , rb_inter             & ! intent(in)
                                , rb_slope             ! ! intent(in)
      use consts_coms    , only : vonk8                & ! intent(in)
                                , cpi8                 & ! intent(in)
                                , sqrt2o28             ! ! intent(in)
      use soil_coms      , only : snow_rough           & ! intent(in)
                                , soil_rough           ! ! intent(in)
      use allometry      , only : h2trunkh             ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)     , target     :: csite
      type(patchtype)    , pointer    :: cpatch
      type(rk4patchtype) , target     :: initp
      integer            , intent(in) :: isi           ! Site loop
      integer            , intent(in) :: ipa           ! Patch loop
      logical            , intent(in) :: get_flow_geom
      !----- Local variables --------------------------------------------------------------!
      integer        :: ico        ! Cohort loop
      integer        :: ipft       ! PFT alias
      integer        :: k          ! Elevation index
      integer        :: zcan       ! Index of canopy top elevation
      real(kind=8)   :: sigma_e    ! Vortex penetration depth                   [        m]
      real(kind=8)   :: crown_d    ! Diameter of a plant's crown                [        m]
      real(kind=8)   :: Re_c       ! Canopy Reynolds Number                     [      ---]
      real(kind=8)   :: Cd         ! Canopy Drag Coefficient                    [      ---]
      real(kind=8)   :: a_front    ! Frontal area / vol.of flow exposed to drag [    m2/m3]
      real(kind=8)   :: zetac      ! Cumulative zeta function                   [      ---]
      real(kind=8)   :: K_top      ! Diffusivity at canopy top z=h              [     m2/s]
      real(kind=8)   :: crowndepth ! Depth of vegetation leafy crown            [        m]
      real(kind=8)   :: h          ! Canopy height                              [        m]
      real(kind=8)   :: layertai   ! The total leaf area of that discrete layer [         ]
      real(kind=8)   :: idz        ! Height of the lowest discrete canopy level [        m]
      real(kind=8)   :: z          ! Elevation in the canopy                    [        m]
      real(kind=8)   :: Kdiff      ! Diffusivity                                [     m2/s]
      real(kind=8)   :: z_om       ! Roughness length that imposes drag 
                                   !     on the ABL winds                       [        m]
      real(kind=8)   :: zref       ! Reference height                           [        m]
      real(kind=8)   :: vels_ref   ! Wind speed at the refernce height
      real(kind=8)   :: surf_rough ! Roughness length of the bare ground 
                                   !     at canopy bottom                       [        m]
      real(kind=8)   :: uh         ! Wind speed at the canopy top (z=h)         [      m/s]
      real(kind=8)   :: uz         ! Wind speed at some height z below h        [      m/s]
      real(kind=8)   :: fm         ! Stability parameter (multiplicative) using RIb
      real(kind=8)   :: can_theta  ! Canopy air potential temperature           [        K]
      real(kind=8)   :: rb_max     ! Maximum aerodynamic resistance.            [      s/m]
      real(kind=8)   :: hite8      ! Double precision version of height.        [        m]
      !----- Saved variables --------------------------------------------------------------!
      real(kind=8), dimension(200), save :: zeta     ! Attenuation factor for sub-canopy K 
                                                     !    and u.  A vector size of 200, 
                                                     !    allows for trees 100 meters tall
      real(kind=8),save                  :: d0       ! Zero-plane displacement height (m)
      real(kind=8),save                  :: ustarouh ! The ratio of ustar over u(h)
      real(kind=8),save                  :: eta      ! The in-canopy wind attenuation scal-
                                                     !    ing parameter
      !------ External procedures ---------------------------------------------------------!
      real        , external             :: sngloff  ! Safe double -> simple precision.
      !------------------------------------------------------------------------------------!

      cpatch=>csite%patch(ipa)

      !---- Finding the maximum aerodynamic resistance. -----------------------------------!
      rb_max = dble(rb_inter)                                                              &
             + dble(rb_slope) * (dble(csite%lai(ipa)) + dble(csite%wpa(ipa)))

      !------------------------------------------------------------------------------------!
      !     If there is no vegetation in this patch, then we apply turbulence to bare      !
      ! soil, no d0 and exit.                                                              !
      !------------------------------------------------------------------------------------!
      if (cpatch%ncohorts == 0) then
         
         vels_ref = rk4met%vels
         zref     = rk4met%geoht
         h        = 0.d0
         d0       = 0.d0

         !----- Calculate the surface roughness inside the canopy. ------------------------!
         initp%rough = dble(soil_rough)*(1.d0 - dble(csite%snowfac(ipa)))                  &
                     + dble(snow_rough)*dble(csite%snowfac(ipa))
         
         !----- Finding the characteristic scales (a.k.a. stars). -------------------------!
         can_theta = cpi8 * rk4met%exner * initp%can_temp
         call ed_stars8(rk4met%atm_theta,rk4met%atm_shv,rk4met%atm_co2, can_theta          &
                       ,initp%can_shv,initp%can_co2,zref,d0,vels_ref,initp%rough           &
                       ,initp%ustar,initp%tstar,initp%qstar,initp%cstar,fm)

         !---------------------------------------------------------------------------------!
         !      The surface resistance inside vegetated canopies is inconsequential, so    !
         ! just give it a nominal zero value.                                              !
         !---------------------------------------------------------------------------------!
         initp%rasveg = 0.0  


         !---------------------------------------------------------------------------------!
         !     Calculate the heat and mass storage capacity of the canopy.                 !
         !---------------------------------------------------------------------------------!
         call can_whcap(csite,ipa,h,zref)
         return
      end if


      !------------------------------------------------------------------------------------!
      !     In case we do have cohorts, choose which method we use to compute the          !
      ! resistance.                                                                        !
      !------------------------------------------------------------------------------------!
      select case (icanturb)

      !------------------------------------------------------------------------------------!
      ! Default Case: There is no use of zero-plane displacement in calculating ustar from !
      !               log-scaling the constant shear layer.  So even if the reference      !
      !               height is inside the canopy, the mathematics will be well behaved,   !
      !               even though it is nonsense.                                          !
      !------------------------------------------------------------------------------------!
      case (0)
         h        = dble(csite%veg_height(ipa)) ! Canopy height
         d0       = 6.3d-1 * h                  ! 0-plane displacement
         vels_ref = rk4met%vels
         zref     = rk4met%geoht

         !---------------------------------------------------------------------------------!
         !      Calculate a surface roughness that is visible to the ABL.  Roughness of    !
         ! the rough wall.                                                                 !
         !---------------------------------------------------------------------------------!
         initp%rough = max(dble(soil_rough),dble(csite%veg_rough(ipa)))                    &
                     * (1.d0 - dble(csite%snowfac(ipa))) + dble(snow_rough)
         
         !----- Calculate the soil surface roughness inside the canopy. -------------------!
         surf_rough = dble(soil_rough) * (1.d0 - dble(csite%snowfac(ipa)))                 &
                    + dble(snow_rough)*dble(csite%snowfac(ipa))


         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !                                                                                 !
         ! IMPORTANT NOTE: the displacement height here (d0), is used only for scaling the !
         !                 vortices when determining diffusivity.  The default method      !
         !                 taken from LEAF-3, as applied here, assumes that the zero plane !
         !                 is at the ground surface when computing the log wind profile,   !
         !                 hence the 0.0 as the argument to ed_stars8.                     !
         !---------------------------------------------------------------------------------!
         can_theta = cpi8 * rk4met%exner * initp%can_temp
         call ed_stars8(rk4met%atm_theta,rk4met%atm_shv,rk4met%atm_co2, can_theta          &
                       ,initp%can_shv,initp%can_co2,zref,0.d0,vels_ref,initp%rough         &
                       ,initp%ustar,initp%tstar,initp%qstar,initp%cstar,fm)

         K_top = vonk8 * initp%ustar*(h-d0)

         !---------------------------------------------------------------------------------!
         !     The surface resistance in the sub-canopy layer is the integration of        !
         ! inverse K, from the rough soil surface, to the reference point in the canopy    !
         ! where the state variable is integrated (canopy top ~ h).                        !
         !---------------------------------------------------------------------------------!
         initp%rasveg = exp(exar8) * h / (exar8 * K_top)                                   &
                      * (exp(-exar8 * surf_rough/h)-exp(-exar8))


         !---------------------------------------------------------------------------------!
         !     This part of the code initializes the geometry of the canopy air space, the !
         ! structure of the vegetation and its attenuation effects and the heat and water  !
         ! capacities.                                                                     !
         !---------------------------------------------------------------------------------!
         if(get_flow_geom) then

            !----- Top of canopy wind speed. ----------------------------------------------!
            uh = (initp%ustar/vonk8)*(log(h/initp%rough)/sqrt(fm))
            
            do ico=1,cpatch%ncohorts

               ipft  = cpatch%pft(ico)
               hite8 = dble(cpatch%hite(ico))

               !----- Estimate the height center of the crown. ----------------------------!
               !z = hite8 * (1.d0 - 5.d-1 * dble(crown_depth_fraction(ipft)))
               z = 5.d-1 * (hite8 + dble(h2trunkh(cpatch%hite(ico))))

               !----- Calculate the wind speed at height z. -------------------------------!
               uz = uh * exp(-exar8 * (1.d0 - z/h))
               
               !----- Find the aerodynamic resistance. [s/m] ------------------------------!
               initp%rb(ico) = min(rb_max                                                  &
                                  ,1.d0 / ( 3.d-3 * sqrt(uz/dble(leaf_width(ipft)))        &
                                          + 5.d-1 * 2.06d-5                                &
                                          * (1.6d8*abs(initp%veg_temp(ico)-initp%can_temp) &
                                            * dble(leaf_width(ipft))**3 ) **2.5d-1         &
                                          / dble(leaf_width(ipft)) ) )
               cpatch%rb(ico) = sngloff(initp%rb(ico),tiny_offset)
            end do

            !------------------------------------------------------------------------------!
            !    Calculate the heat and mass storage capacity of the canopy and inter-     !
            ! facial air spaces.  This is a tough call, because the reference height is    !
            ! allowed to be abnormally low in this case, and it is possible that it is     !
            ! even lower than the top of the canopy.  So... we will set the top of the     !
            ! interfacial layer as the "reference elevation plus the top of the canopy".   !
            ! An alternative could be to make a conditional like in case(1).               !
            !------------------------------------------------------------------------------!
            call can_whcap(csite,ipa,h,zref+h)
         end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If the reference height is less than the height of the canopy, then we have to !
      ! realize that we are in the exponential extinction layer.  If we are above the      !
      ! canpopy, we are in a log extinction zone (interfacial layer).  So if zref<h, use   !
      ! the exponential decay function to calculate u at the canopy top, and then use that !
      ! velocity to estimate ustar.                                                        !
      !------------------------------------------------------------------------------------!
      case(1)
         h    = dble(csite%veg_height(ipa)) ! Canopy height
         d0   = 6.3d-1 * h                  ! 0-plane displacement
         
         !----- Calculate a surface roughness that is visible to the ABL. -----------------!
         initp%rough = max(dble(soil_rough),dble(csite%veg_rough(ipa)))                    &
                     * (1.d0 - dble(csite%snowfac(ipa))) + dble(snow_rough)
         
         !----- Calculate the soil surface roughness inside the canopy. -------------------!
         surf_rough = dble(soil_rough) * (1.d0 - dble(csite%snowfac(ipa)))                 &
                    + dble(snow_rough)*dble(csite%snowfac(ipa))
         
         !----- Check what is the relative position of our reference data. ----------------!
         if (rk4met%geoht < h) then
            !----- First, find the wind speed at the canopy top. --------------------------!
            vels_ref = rk4met%vels / exp(-exar8 *(1.d0 - zref/h))
            !----- Assume a new reference elevation at the canopy top. --------------------!
            zref = h
            if (get_flow_geom) call can_whcap(csite,ipa,h,zref+h)

         else
            vels_ref = rk4met%vels
            zref     = rk4met%geoht
            if(get_flow_geom) call can_whcap(csite,ipa,h,zref)
         end if
         
         !---------------------------------------------------------------------------------!
         !      Get ustar for the ABL, assume it is a dynamic shear layer that generates a !
         ! logarithmic profile of velocity.                                                !
         !---------------------------------------------------------------------------------!

         can_theta = cpi8 * rk4met%exner * initp%can_temp
         call ed_stars8(rk4met%atm_theta,rk4met%atm_shv,rk4met%atm_co2, can_theta          &
                       ,initp%can_shv,initp%can_co2,zref,d0,vels_ref,initp%rough           &
                       ,initp%ustar,initp%tstar,initp%qstar,initp%cstar,fm)

         K_top = vonk8 * initp%ustar * (h-d0)

         !---------------------------------------------------------------------------------!
         !     The surface resistance in the sub-canopy layer is the integration of        !
         ! inverse K, from the rough soil surface, to the reference point in the canopy    !
         ! where the state variable is integrated (canopy top ~ h).                        !
         !---------------------------------------------------------------------------------!
         initp%rasveg = exp(exar8) * h / (exar8 * K_top)                                   &
                      * (exp(-exar8 * surf_rough/h)-exp(-exar8))

         !---------------------------------------------------------------------------------!
         !     Calculate the leaf level aerodynamic resistance.                            !
         !---------------------------------------------------------------------------------!
         if(get_flow_geom) then

            !----- Top of canopy wind speed. ----------------------------------------------!
            uh = (initp%ustar/vonk8)*(log(h/initp%rough)/sqrt(fm))

            do ico=1,cpatch%ncohorts
               ipft  = cpatch%pft(ico)
               hite8 = dble(cpatch%hite(ico))
               !----- Estimate the height center of the crown. ----------------------------!
               z = hite8 * (1.d0-5.d-1*dble(crown_depth_fraction(ipft)))

               !----- Calculate the wind speed at height z. -------------------------------!
               uz = uh * exp(-exar8 * (1.d0 - z/h))
               
               !----- Find the aerodynamic resistance. [s/m] ------------------------------!
               initp%rb(ico) = min(rb_max                                                  &
                                  ,1.d0 / ( 3.d-3 * sqrt(uz/dble(leaf_width(ipft)))        &
                                          + 5.d-1 * 2.06d-5                                &
                                          * (1.6d8*abs(initp%veg_temp(ico)-initp%can_temp) &
                                            * dble(leaf_width(ipft))**3 ) **2.5d-1         &
                                          / dble(leaf_width(ipft)) ) )
               cpatch%rb(ico) = sngloff(initp%rb(ico),tiny_offset)
            end do
         end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Use the methods of Massman 97, to close the equations of estimating turbulent !
      ! transfer.  If the velocity reference height is not greater than the canopy height, !
      ! we can do a few things:                                                            !
      ! 1. Stop the run.                                                                   !
      ! 2. Say that zref is really h+zref (sketchy...)                                     !
      ! 3. Say that zref is 2*h           (sketchy...)                                     !
      !------------------------------------------------------------------------------------!
      case(2)
         
         !----- Calculate the soil surface roughness inside the canopy. -------------------!
         surf_rough = dble(soil_rough) * (1.d0 - dble(csite%snowfac(ipa)))                 &
                    + dble(snow_rough)*dble(csite%snowfac(ipa))


         !---------------------------------------------------------------------------------!
         !    LAI base drag calculation for center of pressure, d0.  Cumulative LAI based  !
         ! velocity attenuation in the canopy.                                             !
         !---------------------------------------------------------------------------------!
         h = dble(cpatch%hite(1))

         if (rk4met%geoht < h) then

            !----- 1. Stop the run. -------------------------------------------------------!
            write (unit=*,fmt='(a)') ' Your reference height is too low...'
            write (unit=*,fmt='(a,1x,es12.5)') ' Ref. height:           ',rk4met%geoht
            write (unit=*,fmt='(a,1x,es12.5)') ' Tallest cohort height: ',h
            call fatal_error('Bad reference height for M97','canopy_turbulence'            &
                            ,'canopy_struct_dynamics.f90')
            !----- 2. Say that zref is really h+zref (sketchy...). ------------------------!
            !zref     = rk4met%geoht + h
            !vels_ref = rk4met%vels
            !----- 3. Say that zref is 2*h (sketchy...). ----------------------------------!
            !zref     = 2.d0 * rk4met%geoht
            !vels_ref = rk4met%vels
         else
            vels_ref = rk4met%vels
            zref     = rk4met%geoht
         end if

         zcan = ceiling(h/dz)
         idz  = 5.d-1 * dz     ! lowest element center

         if (get_flow_geom) then
            
            !------------------------------------------------------------------------------!
            !    Only go through the double loop once per step, do not go through it       !
            ! during all of the derivative steps.  This loop is only to get the displace-  !
            ! ment height and attenuation which has not changed.                           !
            !------------------------------------------------------------------------------!
            zetac = 0.d0  ! Cumulative zeta

            !------------------------------------------------------------------------------!
            !     Loop through the canopy at equal increments.  At each increment,         !
            ! determine the frontal area drag surface, and the drag force zeta.            !
            !------------------------------------------------------------------------------!
            do k = 1, zcan
               
               !---------------------------------------------------------------------------!
               !    Determine the combined leaf area from all cohort crowns residing in    !
               ! this layer.                                                               !
               !---------------------------------------------------------------------------!
               layertai = 0.d0
               z = dble(k-1) * dz + idz  ! elevation of discrete layer
               do ico=1,cpatch%ncohorts
                  !----- Some aliases. ----------------------------------------------------!
                  ipft  = cpatch%pft(ico)
                  hite8 = dble(cpatch%hite(ico))

                  crowndepth = max(dz,dble(crown_depth_fraction(ipft))*hite8)

                  if( z < hite8 .and. z >= (hite8-crowndepth)) then

                     !---------------------------------------------------------------------!
                     !     Assume that at full leaf-out, there is sheltering of branches.  !
                     ! When leaves are not at full out, then the stems and branches start  !
                     ! to become visible to the fluid flow.  Assume that when leaves are   !
                     ! gone, then the branches contribute about 50% of the drag surface,   !
                     ! everything in between is a linear combination.                      !
                     !---------------------------------------------------------------------!
                     ! layertai=layertai + cpatch%lai(ico)*(dz/crowndepth) 
                     ! layertai=layertai + cpatch%nplant(ico)*0.5*(dz/crowndepth)

                     !---------------------------------------------------------------------!
                     !    Use LAI and WPA to define the frontal area of drag surface.  If  !
                     ! the user decided to ignore branches, ignore them here too.          !
                     !---------------------------------------------------------------------!
                     layertai = layertai + (initp%lai(ico) + initp%wpa(ico))               &
                                           * (dz /crowndepth)
                  end if
               end do

               a_front  = layertai/dz ! Frontal area of drag surface
               zetac    = zetac + 5.d-1 * a_front * (Cd0 / Pm) * dz
               zeta(k)  = zetac       ! Use a centered drag
               zetac    = zetac + 5.d-1 * a_front * (Cd0 / Pm) *dz
            end do
            
            !----- The following constains the ratio of ustar over u. ---------------------!
            ustarouh = (c1_m97 - c2_m97 * exp(-c3_m97*zeta(zcan)))
            
            !----- Eta, coefficient of attenuation. ---------------------------------------!
            eta = 5.d-1 * zeta(zcan) / (ustarouh*ustarouh)
            
            !------------------------------------------------------------------------------!
            !     Displacement height.  Notice that if Pm and Cd0 are uniform, we can      !
            ! estimate zeta(h) without a loop, if this routine is taking to long, we can   !
            ! merge the loops.                                                             !
            !------------------------------------------------------------------------------!
            d0 = h
            do k=1,zcan
               d0 = d0 - dz*exp(-2.d0*eta*(1-zeta(k)/zeta(zcan)))
            end do

            !----- Calculate the roughness lengths zo,zt,zr. ------------------------------!
            initp%rough = max((h-d0)*exp(-vonk8/ustarouh),dble(soil_rough))
         end if

         !----- Calculate ustar, tstar, qstar, and cstar. ---------------------------------!
         can_theta = cpi8 * rk4met%exner * initp%can_temp
         call ed_stars8(rk4met%atm_theta,rk4met%atm_shv,rk4met%atm_co2, can_theta          &
                       ,initp%can_shv,initp%can_co2,zref,d0,vels_ref,initp%rough           &
                       ,initp%ustar,initp%tstar,initp%qstar,initp%cstar,fm)

         if(get_flow_geom) then
            
            !----- Calculate the diffusivity at the canopy top. ---------------------------!
            K_top = vonk8 * initp%ustar * (h-d0)

            initp%rasveg=0.d0

            !----- Numerically integrate the inverse diffusivity. -------------------------!
            do k=1,zcan
               Kdiff        = K_top * exp(-eta * (1.d0-zeta(k)/zeta(zcan))) + kvwake
               initp%rasveg = initp%rasveg + dz / Kdiff
            end do

            !------------------------------------------------------------------------------!
            !     Calculate the leaf level aerodynamic resistance.                         !
            !------------------------------------------------------------------------------!
            
            !----- Top of canopy wind speed. ----------------------------------------------!
            uh = (initp%ustar/vonk8)*(log((h-d0)/initp%rough)/sqrt(fm))

            do ico=1,cpatch%ncohorts
               ipft = cpatch%pft(ico)
               hite8 = dble(cpatch%hite(ico))

               !----- Estimate the height center of the crown. ----------------------------!
               z = hite8 * (1.d0 - 5.d-1 * dble(crown_depth_fraction(ipft)))

               !----- Determine the zeta index. -------------------------------------------!
               k = ceiling(z/dz)

               !----- Calculate the wind speed at height z. -------------------------------!
               uz = uh * exp(-eta * (1.d0 - zeta(k)/zeta(zcan) ))
               
               !----- Find the aerodynamic resistance. [s/m] ------------------------------!
               initp%rb(ico) = min(rb_max                                                  &
                                  ,1.d0 / ( 3.d-3 * sqrt(uz/dble(leaf_width(ipft)))        &
                                          + 5.d-1 * 2.06d-5                                &
                                          * (1.6d8*abs(initp%veg_temp(ico)-initp%can_temp) &
                                            * dble(leaf_width(ipft))**3 ) **2.5d-1         &
                                          / dble(leaf_width(ipft)) ) )
               cpatch%rb(ico) = sngloff(initp%rb(ico),tiny_offset)
            end do

            !------------------------------------------------------------------------------!
            ! Calculate the heat and mass storage capacity of the canopy and interfacial   !
            ! air spaces.                                                                  !
            !------------------------------------------------------------------------------!
            call can_whcap(csite,ipa,h,zref)
         end if

      end select

      return
   end subroutine canopy_turbulence
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine computes the characteristic scales based on  Louis (1981) surface   !
   ! layer parameterization.  Assume that stability is calculated on the potential         !
   ! gradient from the surface to the atmosphere at reference height.  Apply the instabi-  !
   ! lity parameters to calculate the friction velocity.  Use the instability parameters   !
   ! for momentum and scalars, that were calculated over the distance from surface to      !
   ! refernece height to determine the heat, moisture and carbon flux rates at the canopy  !
   ! to atmosphere at reference height.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine ed_stars(theta_atm,shv_atm,co2_atm,theta_can,shv_can,co2_can,zref,d0,uref    &
                      ,rough,ustar,tstar,qstar,cstar,fm)
      use consts_coms     , only : grav   & ! intent(in)
                                 , vonk   ! ! intent(in)
      use canopy_air_coms , only : ustmin & ! intent(in)
                                 , ubmin  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)  :: theta_atm ! Above canopy air pot. temperature        [        K]
      real, intent(in)  :: shv_atm   ! Above canopy vapour spec. hum.           [kg/kg_air]
      real, intent(in)  :: co2_atm   ! CO2 specific volume                      [  痠ol/m設
      real, intent(in)  :: theta_can ! Canopy air potential temperature         [        K]
      real, intent(in)  :: shv_can   ! Canopy air vapour spec. humidity         [kg/kg_air]
      real, intent(in)  :: co2_can   ! Canopy air CO2 specific volume           [  痠ol/m設
      real, intent(in)  :: zref      ! Height at reference point                [        m]
      real, intent(in)  :: d0        ! Zero-plane displacement height           [        m]
      real, intent(in)  :: uref      ! Wind speed at reference height           [      m/s]
      real, intent(in)  :: rough     ! Roughness                                [        m]
      real, intent(out) :: ustar     ! U*, friction velocity                    [      m/s]
      real, intent(out) :: qstar     ! Specific humidity friction scale         [kg/kg_air]
      real, intent(out) :: tstar     ! Temperature friction scale               [        K]
      real, intent(out) :: cstar     ! CO2 spec. volume friction scale          [  痠ol/m設
      real, intent(out) :: fm        ! Stability parameter for momentum
      !----- Local variables --------------------------------------------------------------!
      real              :: a2        ! Drag coefficient in neutral conditions, 
                                     !     here same for h/m
      real              :: c1        !
      real              :: ri        ! Bulk richardson numer, eq. 3.45 in Garratt
      real              :: fh        ! Stability parameter for heat
      real              :: c2        !
      real              :: cm        !
      real              :: ch        !
      real              :: c3        !
      real              :: vels_pat  !
      !----- Constants --------------------------------------------------------------------!
      real, parameter   :: b   = 5.0 !
      real, parameter   :: csm = 7.5 !
      real, parameter   :: csh = 5.0 !
      real, parameter   :: d   = 5.0 !
      !------------------------------------------------------------------------------------!

      vels_pat = max(uref,ubmin)

      !----- Make the log profile (constant shear assumption). ----------------------------!
      a2 = (vonk / log((zref-d0)/rough)) ** 2.
      c1 = a2 * vels_pat

      ri = grav * (zref-d0) * (theta_atm - theta_can)                                      &
         / (0.5 * (theta_atm + theta_can) * vels_pat * vels_pat )
     
      if (theta_atm - theta_can > 0.0) then
         !----- Stable case ---------------------------------------------------------------!
         fm = 1.0 / (1.0 + (2.0 * b * ri / sqrt(1.0 + d * ri)))
         fh = 1.0 / (1.0 + (3.0 * b * ri * sqrt(1.0 + d * ri)))
      else
         !----- Unstable case -------------------------------------------------------------!
         c2 = b * a2 * sqrt( (zref-d0)/rough * (abs(ri)))
         cm = csm * c2
         ch = csh * c2
         fm = (1.0 - 2.0 * b * ri / (1.0 + 2.0 * cm))
         fh = (1.0 - 3.0 * b * ri / (1.0 + 3.0 * ch))
      end if

      ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
      c3 = c1 * fh / ustar

      qstar = c3 * (shv_atm - shv_can)
      tstar = c3 * (theta_atm - theta_can)
      cstar = c3 * (co2_atm - co2_can)

      return
   end subroutine ed_stars
   !==========================================================================================!
   !==========================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine computes the characteristic scales based on  Louis (1981) surface   !
   ! layer parameterization.  Assume that stability is calculated on the potential         !
   ! gradient from the surface to the atmosphere at reference height.  Apply the instabi-  !
   ! lity parameters to calculate the friction velocity.  Use the instability parameters   !
   ! for momentum and scalars, that were calculated over the distance from surface to      !
   ! refernece height to determine the heat, moisture and carbon flux rates at the canopy  !
   ! to atmosphere at reference height.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine ed_stars8(theta_atm,shv_atm,co2_atm,theta_can,shv_can,co2_can,zref,d0,uref   &
                       ,rough,ustar,tstar,qstar,cstar,fm)
      use consts_coms     , only : grav8   & ! intent(in)
                                 , vonk8   ! ! intent(in)
      use canopy_air_coms , only : ustmin8 & ! intent(in)
                                 , ubmin8  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: theta_atm ! Above canopy air pot. temperature [        K]
      real(kind=8), intent(in)  :: shv_atm   ! Above canopy vapour spec. hum.    [kg/kg_air]
      real(kind=8), intent(in)  :: co2_atm   ! CO2 specific volume               [  痠ol/m設
      real(kind=8), intent(in)  :: theta_can ! Canopy air potential temperature  [        K]
      real(kind=8), intent(in)  :: shv_can   ! Canopy air vapour spec. humidity  [kg/kg_air]
      real(kind=8), intent(in)  :: co2_can   ! Canopy air CO2 specific volume    [  痠ol/m設
      real(kind=8), intent(in)  :: zref      ! Height at reference point         [        m]
      real(kind=8), intent(in)  :: d0        ! Zero-plane displacement height    [        m]
      real(kind=8), intent(in)  :: uref      ! Wind speed at reference height    [      m/s]
      real(kind=8), intent(in)  :: rough     ! Roughness                         [        m]
      real(kind=8), intent(out) :: ustar     ! U*, friction velocity             [      m/s]
      real(kind=8), intent(out) :: qstar     ! Specific humidity friction scale  [kg/kg_air]
      real(kind=8), intent(out) :: tstar     ! Temperature friction scale        [        K]
      real(kind=8), intent(out) :: cstar     ! CO2 spec. volume friction scale   [  痠ol/m設
      real(kind=8), intent(out) :: fm        ! Stability parameter for momentum
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)              :: a2        ! Drag coefficient in neutral conditions, 
                                             !     here same for h/m
      real(kind=8)              :: c1        !
      real(kind=8)              :: ri        ! Bulk richardson numer, eq. 3.45 in Garratt
      real(kind=8)              :: fh        ! Stability parameter for heat
      real(kind=8)              :: c2        !
      real(kind=8)              :: cm        !
      real(kind=8)              :: ch        !
      real(kind=8)              :: c3        !
      real(kind=8)              :: vels_pat  !
      !----- Constants --------------------------------------------------------------------!
      real(kind=8), parameter   :: b   = 5.0d0 !
      real(kind=8), parameter   :: csm = 7.5d0 !
      real(kind=8), parameter   :: csh = 5.0d0 !
      real(kind=8), parameter   :: d   = 5.0d0 !
      !------------------------------------------------------------------------------------!

      vels_pat = max(uref,ubmin8)

      !----- Make the log profile (constant shear assumption). ----------------------------!
      a2 = (vonk8 / log((zref-d0)/rough)) ** 2.
      c1 = a2 * vels_pat

      ri = grav8 * (zref-d0) * (theta_atm - theta_can)                                     &
         / (5.d-1 * (theta_atm + theta_can) * vels_pat * vels_pat )
     
      if (theta_atm - theta_can > 0.d0) then
         !----- Stable case ---------------------------------------------------------------!
         fm = 1.d0 / (1.d0 + (2.d0 * b * ri / sqrt(1.d0 + d * ri)))
         fh = 1.d0 / (1.d0 + (3.d0 * b * ri * sqrt(1.d0 + d * ri)))
      else
         !----- Unstable case -------------------------------------------------------------!
         c2 = b * a2 * sqrt( (zref-d0)/rough * (abs(ri)))
         cm = csm * c2
         ch = csh * c2
         fm = (1.d0 - 2.d0 * b * ri / (1.d0 + 2.d0 * cm))
         fh = (1.d0 - 3.d0 * b * ri / (1.d0 + 3.d0 * ch))
      end if

      ustar = max(ustmin8,sqrt(c1 * vels_pat * fm))
      c3 = c1 * fh / ustar

      qstar = c3 * (shv_atm   - shv_can   )
      tstar = c3 * (theta_atm - theta_can )
      cstar = c3 * (co2_atm   - co2_can   )

      return
   end subroutine ed_stars8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function vertical_vel_flux(gzotheta,tstar,ustar)
      use consts_coms , only : vonk ! intent(in)
     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real, intent(in)    :: ustar
      real, intent(in)    :: tstar
      real, intent(in)    :: gzotheta
      !----- Local variables --------------------------------------------------------------!
      real                :: zoverl
      real                :: cx
      real                :: psin
      !----- Constants --------------------------------------------------------------------!
      real, parameter     :: wtol = 1.e-20
      !------------------------------------------------------------------------------------!
     
     
      zoverl = gzotheta * vonk * tstar / (ustar * ustar)
     
      if (zoverl < 0.0)then
         cx = zoverl * sqrt(sqrt(1.0 - 15.0 * zoverl))
      else
         cx = zoverl / (1.0 + 4.7 * zoverl)
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
   real(kind=8) function vertical_vel_flux8(gzotheta,tstar,ustar)
      use consts_coms , only : vonk8 ! intent(in)
     
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)    :: ustar
      real(kind=8), intent(in)    :: tstar
      real(kind=8), intent(in)    :: gzotheta
      !----- Local variables --------------------------------------------------------------!
      real(kind=8) :: zoverl
      real(kind=8) :: cx
      real(kind=8) :: psin
      !----- Constants --------------------------------------------------------------------!
      real(kind=8), parameter     :: wtol = 1.d-20
      !------------------------------------------------------------------------------------!
     
     
      zoverl = gzotheta * vonk8 * tstar / (ustar * ustar)
     
      if (zoverl < 0.d0)then
         cx = zoverl * sqrt(sqrt(1.d0 - 1.5d1 * zoverl))
      else
         cx = zoverl / (1.0d0 + 4.7d0 * zoverl)
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
   !     Calculate the heat can general scalar storage capacity of the canopy air and      !
   ! interfacial air spaces.  The capacites are essentiall the amount of air mass stored   !
   ! in each space.                                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine can_whcap(csite,ipa,veg_height,zref)

      use rk4_coms             , only : rk4met                & ! intent(in)
                                      , wcapcan               & ! intent(out)
                                      , wcapcani              & ! intent(out)
                                      , hcapcani              ! ! intent(out)
      use canopy_air_coms      , only : minimum_canopy_depth  ! ! intent(in)
      use ed_state_vars        , only : sitetype              ! ! structure
      use consts_coms          , only : cpi8                  ! ! intent(in)
      
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      real(kind=8)   , intent(in) :: zref
      real(kind=8)   , intent(in) :: veg_height
      !----- Local variables --------------------------------------------------------------!
      real                    :: canopy_depth
      !------------------------------------------------------------------------------------!
      
      canopy_depth  = max(veg_height,minimum_canopy_depth)
      
      wcapcan  = rk4met%rhos * canopy_depth
      wcapcani = 1.d0 / wcapcan
      hcapcani = cpi8 * wcapcani

      return
   end subroutine can_whcap
   !=======================================================================================!
   !=======================================================================================!
end module canopy_struct_dynamics
!==========================================================================================!
!==========================================================================================!
