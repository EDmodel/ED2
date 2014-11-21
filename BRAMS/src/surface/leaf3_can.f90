!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the variables that depend on heat, water, and (eventually)  !
! carbon exchanges with the canopy air space and vegetation.  This solves only the         !
! prognostic variables, the diagnostic ones will be updated in a separate sub-routine.     !
!------------------------------------------------------------------------------------------!
subroutine leaf3_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass         &
                       ,ustar,tstar,rstar,cstar,zeta,ribulk,soil_rough,veg_rough           &
                       ,patch_rough,veg_height,veg_displace,veg_lai,veg_tai,veg_water      &
                       ,veg_hcap,veg_energy,leaf_class,veg_fracarea,stom_condct,can_prss   &
                       ,can_rvap,can_co2,hflxac_out,wflxac_out,qwflxac_out,eflxac_out      &
                       ,cflxac_out,hflxgc_out,wflxgc_out,qwflxgc_out,hflxvc_out,wflxvc_out &
                       ,qwflxvc_out,transp_out,qtransp_out,intercepted_out                 &
                       ,qintercepted_out,wshed_out,qwshed_out,throughfall_out              &
                       ,qthroughfall_out,gpp_out,plresp_out,resphet_out,growresp           &
                       ,ground_rsat,ground_rvap,ground_temp,ground_fliq,available_water    &
                       ,rshort)
   use leaf_coms , only : dtl3                 & ! intent(in)
                        , dtl3_factor          & ! intent(in)
                        , dtvg                 & ! intent(in)
                        , ndtvegi              & ! intent(in)
                        , supersat_ok          & ! intent(in)
                        , resolvable           & ! intent(in)
                        , gr_factor            & ! intent(in)
                        , leaf_maxwhc          & ! intent(in)
                        , gbh                  & ! intent(in)
                        , gbw                  & ! intent(in)
                        , gsw                  & ! intent(in)
                        , snowfac              & ! intent(in)
                        , exar                 & ! intent(in)
                        , geoht                & ! intent(in)
                        , atm_vels             & ! intent(in)
                        , atm_temp_zcan        & ! intent(in)
                        , pcpgl                & ! intent(in)
                        , qpcpgl               & ! intent(in)
                        , dpcpgl               & ! intent(in)
                        , can_exner            & ! intent(in)
                        , can_depth            & ! intent(in)
                        , can_temp             & ! intent(in)
                        , can_rhv              & ! intent(in)
                        , can_shv              & ! intent(in)
                        , can_rhos             & ! intent(in)
                        , can_cp               & ! intent(in)
                        , veg_temp             & ! intent(in)
                        , veg_fliq             & ! intent(in)
                        , veg_wind             & ! intent(in)
                        , sfcwater_tempk       & ! intent(in)
                        , sfcwater_fracliq     & ! intent(in)
                        , soil_tempk           & ! intent(in)
                        , soil_fracliq         & ! intent(in)
                        , rshort_v             & ! intent(in)
                        , rlong_v              & ! intent(in)
                        , ggsoil               & ! intent(in)
                        , transp_o             & ! intent(in)
                        , gpp_o                & ! intent(in)
                        , leaf_resp_o          & ! intent(in)
                        , root_resp_o          & ! intent(in)
                        , het_resp_o           & ! intent(in)
                        , drysoil              & ! intent(in)
                        , satsoil              & ! intent(in)
                        , can_enthalpy         & ! intent(inout)
                        , ggbare               & ! intent(out)
                        , dtl3owcc             & ! intent(out)
                        , dtl3ohcc             & ! intent(out)
                        , dtl3occc             & ! intent(out)
                        , rho_ustar            & ! intent(out)
                        , eflxac               & ! intent(out)
                        , hflxac               & ! intent(out)
                        , hflxvc               & ! intent(out)
                        , hflxgc               & ! intent(out)
                        , hflxsc               & ! intent(out)
                        , wflxac               & ! intent(out)
                        , qwflxac              & ! intent(out)
                        , wflxgc               & ! intent(out)
                        , qwflxgc              & ! intent(out)
                        , wflxsc               & ! intent(out)
                        , qwflxsc              & ! intent(out)
                        , wflxvc               & ! intent(out)
                        , qwflxvc              & ! intent(out)
                        , transp               & ! intent(out)
                        , qtransp              & ! intent(out)
                        , hflxvc_tot           & ! intent(out)
                        , wflxvc_tot           & ! intent(out)
                        , qwflxvc_tot          & ! intent(out)
                        , cflxvc_tot           & ! intent(out)
                        , cflxgc_tot           & ! intent(out)
                        , transp_tot           & ! intent(out)
                        , qtransp_tot          & ! intent(out)
                        , cflxac               & ! intent(out)
                        , cflxgc               & ! intent(out)
                        , cflxvc               & ! intent(out)
                        , estar                & ! intent(out)
                        , qstar                & ! intent(out)
                        , ggveg                & ! intent(out)
                        , ggnet                & ! intent(out)
                        , wshed                & ! intent(out)
                        , qwshed               & ! intent(out)
                        , dwshed               & ! intent(out)
                        , intercepted_tot      & ! intent(out)
                        , qintercepted_tot     & ! intent(out)
                        , dintercepted_tot     & ! intent(out)
                        , throughfall_tot      & ! intent(out)
                        , qthroughfall_tot     & ! intent(out)
                        , dthroughfall_tot     & ! intent(out)
                        , wshed_tot            & ! intent(out)
                        , qwshed_tot           & ! intent(out)
                        , dwshed_tot           & ! intent(out)
                        , dewgnd_tot           & ! intent(out)
                        , qdewgnd_tot          & ! intent(out)
                        , ddewgnd_tot          ! ! intent(out)
   use rconstants, only : wdns                 & ! intent(in)
                        , wdnsi                & ! intent(in)
                        , fdnsi                & ! intent(in)
                        , mmdry                & ! intent(in)
                        , mmdryi               & ! intent(in)
                        , day_sec              & ! intent(in)
                        , vonk                 & ! intent(in)
                        , huge_num             & ! intent(in)
                        , toodry               & ! intent(in)
                        , twothirds            & ! intent(in)
                        , pi1                  & ! intent(in)
                        , tiny_num             ! ! intent(in)
   use therm_lib , only : qslif                & ! function
                        , rslif                & ! function
                        , thetaeiv             & ! function
                        , tslif                & ! function
                        , uextcm2tl            & ! subroutine
                        , tl2uint              & ! function
                        , tq2enthalpy          ! ! function
   use mem_leaf  , only : ndtveg               & ! intent(in)
                        , isfcl                ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer             , intent(in)    :: mzg
   integer             , intent(in)    :: mzs
   integer             , intent(in)    :: ksn
   real, dimension(mzg), intent(in)    :: soil_energy
   real, dimension(mzg), intent(in)    :: soil_water
   real, dimension(mzg), intent(in)    :: soil_text
   real, dimension(mzs), intent(in)    :: sfcwater_mass
   real                , intent(in)    :: ustar
   real                , intent(in)    :: tstar
   real                , intent(in)    :: rstar
   real                , intent(in)    :: cstar
   real                , intent(in)    :: zeta
   real                , intent(in)    :: ribulk
   real                , intent(in)    :: soil_rough
   real                , intent(in)    :: veg_rough
   real                , intent(in)    :: patch_rough
   real                , intent(in)    :: veg_height
   real                , intent(in)    :: veg_displace
   real                , intent(in)    :: veg_lai
   real                , intent(in)    :: veg_tai
   real                , intent(in)    :: leaf_class
   real                , intent(in)    :: veg_fracarea
   real                , intent(in)    :: ground_rsat
   real                , intent(in)    :: ground_rvap
   real                , intent(in)    :: ground_temp
   real                , intent(in)    :: ground_fliq
   real                , intent(in)    :: available_water
   real                , intent(in)    :: rshort
   real                , intent(inout) :: veg_water
   real                , intent(inout) :: veg_hcap
   real                , intent(inout) :: veg_energy
   real                , intent(inout) :: stom_condct
   real                , intent(inout) :: can_prss
   real                , intent(inout) :: can_rvap
   real                , intent(inout) :: can_co2
   real                , intent(inout) :: hflxac_out
   real                , intent(inout) :: wflxac_out
   real                , intent(inout) :: qwflxac_out
   real                , intent(inout) :: eflxac_out
   real                , intent(inout) :: cflxac_out
   real                , intent(inout) :: hflxgc_out
   real                , intent(inout) :: wflxgc_out
   real                , intent(inout) :: qwflxgc_out
   real                , intent(inout) :: hflxvc_out
   real                , intent(inout) :: wflxvc_out
   real                , intent(inout) :: qwflxvc_out
   real                , intent(inout) :: transp_out
   real                , intent(inout) :: qtransp_out
   real                , intent(inout) :: intercepted_out
   real                , intent(inout) :: qintercepted_out
   real                , intent(inout) :: wshed_out
   real                , intent(inout) :: qwshed_out
   real                , intent(inout) :: throughfall_out
   real                , intent(inout) :: qthroughfall_out
   real                , intent(inout) :: gpp_out
   real                , intent(inout) :: plresp_out
   real                , intent(inout) :: resphet_out
   real                , intent(inout) :: growresp
   !----- Local arguments. ----------------------------------------------------------------!
   integer                             :: k
   integer                             :: kk
   integer                             :: nveg
   integer                             :: nsoil
   integer                             :: nit
   logical                             :: is_dew_cp
   logical                             :: is_dew_cs
   real                                :: aux
   real                                :: dsm
   real                                :: fac
   real                                :: factv
   real                                :: qwtot
   real                                :: sigmaw 
   real                                :: slai
   real                                :: stai
   real                                :: slpotv
   real                                :: swp
   real                                :: vpd
   real                                :: wtemp
   real                                :: wtroot
   real                                :: x
   real                                :: zognd
   real                                :: zoveg
   real                                :: zdisp
   real                                :: zveg
   real                                :: wflx
   real                                :: old_veg_water
   real                                :: old_veg_energy
   real                                :: shv_gradient
   real                                :: sfcwater_qsat
   real                                :: rfac
   real                                :: gr_weight
   real                                :: ground_shv
   real                                :: ground_qsat
   real                                :: vfarea_eff
   real                                :: dew_now
   !----- Local constants. ----------------------------------------------------------------!
   character(len=9)      , parameter   :: fmti='(a,1x,i6)'
   character(len=13)     , parameter   :: fmtf='(a,1x,es12.5)'
   character(len=3)      , parameter   :: fmtc='(a)'
   character(len=9)      , parameter   :: fmtl='(a,1x,l1)'
   !----- External functions. -------------------------------------------------------------!
   real                  , external    :: leaf3_reduced_wind
   !---------------------------------------------------------------------------------------!



   !----- Find the time step auxiliary variables. -----------------------------------------!
   dtl3owcc  = dtl3 / (can_depth * can_rhos)
   dtl3ohcc  = dtl3 / (can_depth * can_rhos)
   dtl3occc  = dtl3owcc * mmdry
   gr_weight = dtl3 / day_sec
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the atmosphere -> canopy fluxes.                                             !
   !---------------------------------------------------------------------------------------!
   rho_ustar  = can_rhos  * ustar
   eflxac     = rho_ustar * estar                 ! Enthalpy exchange
   hflxac     = rho_ustar * tstar * can_exner     ! Sensible heat exchange
   wflxac     = rho_ustar * qstar                 ! Water vapour exchange
   !----- "Latent" heat flux. -------------------------------------------------------------!
   if (wflxac >= 0.0) then
      qwflxac = wflxac * tq2enthalpy(atm_temp_zcan,1.0,.true.)
   else
      qwflxac = wflxac * tq2enthalpy(can_temp,1.0,.true.)
   end if
   cflxac     = rho_ustar * cstar * mmdryi        ! CO2 exchange
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Compute ground-canopy resistance rgnd.  Assume zognd not affected by snow.        !
   ! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining the values          !
   ! (zognd,0) when vegetation is fully buried in snow.                                    !
   !---------------------------------------------------------------------------------------!
   zognd = soil_rough
   zoveg = veg_rough    * (1. - snowfac) + zognd * snowfac
   zdisp = veg_displace * (1. - snowfac)
   zveg  = veg_height   * (1. - snowfac)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the ground conductance due to vegetation and find the total conductance      !
   ! (the inverse of the sum due to bare ground plus vegetation).                          !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then

      !------------------------------------------------------------------------------------!
      !    If vegetation is sufficiently abundant and not covered by snow, compute  heat   !
      ! and moisture fluxes from vegetation to canopy, and flux resistance from soil or    !
      ! snow to canopy.                                                                    !
      !------------------------------------------------------------------------------------!
      factv = log((geoht-zdisp) / zoveg) / (vonk * vonk * atm_vels)
      aux   = exp(exar * (1. - (zdisp + zoveg) / zveg))
      ggveg = (exar * (zveg - zdisp)) / (factv * zveg  * (exp(exar) - aux))
      ggnet = ggbare * ggveg / (ggbare + ggveg)
   else
      ggveg = huge_num
      ggnet = ggbare
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Split the total precipitation into intercepted water and through-fall, according  !
   ! to the fraction of open canopy.                                                       !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then
      vfarea_eff       = veg_fracarea * (1. - snowfac)
      intercepted_tot  = pcpgl  * vfarea_eff
      qintercepted_tot = qpcpgl * vfarea_eff
      dintercepted_tot = dpcpgl * vfarea_eff
   else
      vfarea_eff       = 0.0
      intercepted_tot  = 0.0
      qintercepted_tot = 0.0
      dintercepted_tot = 0.0
   end if
   !----- The throughfall is what is not intercepted.  Shedding is zero at this point. ----!
   throughfall_tot  = pcpgl  - intercepted_tot
   qthroughfall_tot = qpcpgl - qintercepted_tot
   dthroughfall_tot = dpcpgl - dintercepted_tot
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between the ground and the canopy air   !
   ! space.  The ground may be either the top soil layer or the temporary surface water    !
   ! snow surface.  First we check whether the vapour fluxes to ground are positive or     !
   ! negative.  Later, three variables will define the flux between ground and canopy air  !
   ! space.  They are either positive or zero:                                             !
   !                                                                                       !
   ! - wflxsc [kg/m2/s] is going to be the evaporation flux from pounding water or         !
   !                    snowpack to canopy.                                                !
   ! - wflxgc [kg/m2/s] is going to be the evaporation flux from top soil layer to canopy. !
   ! - dewgnd [kg/m2/s] is going to be the dew/frost flux from canopy to ground (this      !
   !                    encompasses dew formation on both because dew always goes to the   !
   !                    virtual layer).                                                    !
   !                                                                                       !
   !     Sensible heat is defined by two variables, hflxgc (soil) and hflxsc (TSW)         !
   ! [J/m2/s], which can be either positive or negative.                                   !
   !---------------------------------------------------------------------------------------!
   if (ksn > 0) then
      hflxsc        = snowfac * ggnet * can_rhos * can_cp * (sfcwater_tempk(ksn) - can_temp)
      sfcwater_qsat = qslif(can_prss,sfcwater_tempk(ksn))
      is_dew_cp     = sfcwater_qsat <= can_shv
   else
      hflxsc        = 0.0
      sfcwater_qsat = can_shv
      is_dew_cp     = .false.
   end if
   hflxgc      = (1.0 - snowfac) * ggnet * can_rhos * can_cp * (ground_temp - can_temp)
   ground_qsat = ground_rsat / ( 1.0 + ground_rsat )
   is_dew_cs   = ground_qsat <= can_shv
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Here we will decide how to compute the evaporation and condensation fluxes        !
   ! between temporary surface water and canopy air space, based on the sign of water      !
   ! flux.                                                                                 !
   !---------------------------------------------------------------------------------------!
   if (is_dew_cp) then
      !----- Dew should be formed due to contact between top TSW and canopy air space. ----!
      dew_now     = snowfac * ggnet * can_rhos * ( can_shv - sfcwater_qsat )
      dewgnd_tot  = dew_now
      qdewgnd_tot = dew_now * tq2enthalpy(sfcwater_tempk(ksn),1.00,.true.)
      ddewgnd_tot = dew_now * (         sfcwater_fracliq(ksn)   * wdnsi                    &
                              + ( 1.0 - sfcwater_fracliq(ksn) ) * fdnsi )
      wflxsc      = 0.0
      qwflxsc     = 0.0
   elseif (ksn > 0 .and. ( can_rhv < 1.0 .or. supersat_ok )) then
      !------------------------------------------------------------------------------------!
      !     Evaporation should occur, and canopy air space is not yet super-saturated (or  !
      ! the user is fine with super-saturation).                                           !
      !------------------------------------------------------------------------------------!
      dewgnd_tot  = 0.0
      qdewgnd_tot = 0.0
      ddewgnd_tot = 0.0
      wflxsc      = snowfac * ggnet * can_rhos * ( sfcwater_qsat - can_shv )
      qwflxsc     = wflxsc * tq2enthalpy(sfcwater_tempk(ksn),1.0,.true.)
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Similarly, we will decide how to compute the evaporation and condensation fluxes  !
   ! between top soil layer and canopy air space, based on the sign of water flux.         !
   !---------------------------------------------------------------------------------------!
   if (is_dew_cs) then
      !----- Dew should be formed due to contact between top soil and canopy air space. ---!
      dew_now     = ( 1.d0 - snowfac ) * ggnet * can_rhos * ( can_shv - ground_qsat )
      dewgnd_tot  = dewgnd_tot  + dew_now
      qdewgnd_tot = qdewgnd_tot + dew_now * tq2enthalpy(ground_temp,1.0,.true.)
      ddewgnd_tot = ddewgnd_tot + dew_now * (         ground_fliq   * wdnsi                &
                                            + ( 1.0 - ground_fliq ) * fdnsi )
      wflxgc      = 0.0
      qwflxgc     = 0.0
   else if (drysoil(mzg)) then
      !------------------------------------------------------------------------------------!
      !    There should be evaporation, except that there is no water left to be extracted !
      ! from the ground... Set both evaporation and condensation fluxes to zero.           !
      !------------------------------------------------------------------------------------!
      wflxgc      = 0.0
      qwflxgc     = 0.0
      !------------------------------------------------------------------------------------!

   elseif (can_rhv < 1.0 .or. supersat_ok) then
      !------------------------------------------------------------------------------------!
      !     Evaporation should occur, and canopy air space is not yet super-saturated (or  !
      ! the user is fine with super-saturation).                                           !
      !------------------------------------------------------------------------------------!
      wflxgc      = ( 1.0 - snowfac ) * ggnet * can_rhos * ( ground_qsat - can_shv )       &
                  * ( 1.d0 / (1.d0 + ggnet / ggsoil) )
      qwflxgc     = wflxgc * tq2enthalpy(ground_temp,1.0,.true.)
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Respiration driver.                                                               !
   !---------------------------------------------------------------------------------------!
   call leaf3_soil_resp(mzg,resolvable,soil_energy,soil_water,soil_text,leaf_class,veg_lai &
                       ,veg_tai)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     If vegetation is sufficiently abundant and not covered by snow, compute heat and  !
   ! moisture fluxes from vegetation to canopy, and flux resistance from soil or snow to   !
   ! canopy.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then
      !----- Save the vegetation class in an alias. ---------------------------------------!
      nveg = nint(leaf_class)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the "perceived" leaf area and tree area indices, dependending on snow     !
      ! cover.                                                                             !
      !------------------------------------------------------------------------------------!
      stai = veg_tai * ( 1.0 - snowfac )
      slai = veg_lai * ( 1.0 - snowfac )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the wind speed at the top of the canopy, using the similarity theory.     !
      !------------------------------------------------------------------------------------!
      veg_wind = leaf3_reduced_wind(ustar,zeta,ribulk,geoht,veg_displace,veg_height        &
                                  ,patch_rough,vfarea_eff,stai)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Vegetation time step.  The fluxes have already been flushed to zero before    !
      ! entering here, so in principle there is no need to flush them again.               !
      !------------------------------------------------------------------------------------!
      vegloop: do nit = 1,ndtveg
         !---------------------------------------------------------------------------------!
         !    Find the aerodynamic conductances for heat and water at the leaf boundary    !
         ! layer.                                                                          !
         !---------------------------------------------------------------------------------!
         call leaf3_aerodynamic_conductances(leaf_class,veg_wind,veg_temp,can_temp         &
                                            ,can_rvap,can_rhos,can_cp)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the intercellular vapour pressure and mixing ratio; they are both      !
         ! assumed to be the saturation value at the leaf temperature.                     !
         !---------------------------------------------------------------------------------!
         shv_gradient  = qslif(can_prss,veg_temp) - can_shv
         rfac          = (1. + rslif(can_prss,veg_temp)) * (1. + can_rvap)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Photosynthesis driver.                                                      !
         !---------------------------------------------------------------------------------!
         call leaf3_photo_driver(mzg,soil_water,soil_text,leaf_class,rshort,veg_lai        &
                                ,veg_tai,can_prss,can_rvap,can_co2,stom_condct)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Sensible heat flux from leaf to canopy.  The 2*lai means that heat         !
         ! exchange occurs on both sides of the leaf, whereas pi1*(stai-slai) means that   !
         ! heat exchange occurs through the entire perimeter of branches and twigs.        !
         !---------------------------------------------------------------------------------!
         hflxvc = ( 2. * slai + pi1 * ( stai - slai )) * gbh * (veg_temp - can_temp)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check whether evaporation or condensation should happen.                    !
         !---------------------------------------------------------------------------------!
         if (shv_gradient > 0.0 .and. can_rhv >= 1.0 .and. (.not. supersat_ok)) then
            !------------------------------------------------------------------------------!
            !     Stop flux to canopy air space if relative humidity is 100% and the user  !
            ! doesn't want super-saturation to happen.                                     !
            !------------------------------------------------------------------------------!
            wflxvc = 0.0
            !------------------------------------------------------------------------------!
         elseif (shv_gradient >= 0.0) then
            !------------------------------------------------------------------------------!
            !     Saturation at leaf temperature is greater than canopy air space mixing   !
            ! ratio.  If some water is sitting on the leaf surface, let it evaporate.      !
            ! Evaporation should happen from only one side of the leaf.                    !
            !------------------------------------------------------------------------------!
            !----- Compute the fraction of leaves covered in water. -----------------------!
            sigmaw = min(1.0, (veg_water / (.2 * stai)) ** twothirds)
            !----- Find the evaporation rate. ---------------------------------------------!
            wflxvc = min(stai * gbw * shv_gradient * rfac * sigmaw,veg_water/dtvg)
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     Flow is towards the leaf, dew/frost formation on one side of the leaf.   !
            ! We don't let the dew formation completely eliminate water from the canopy    !
            ! air space, and impose 0. transpiration.                                      !
            !------------------------------------------------------------------------------!
            wflxvc = max(stai * gbw * shv_gradient * rfac, (10.*toodry-can_shv) / dtl3owcc)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Update vegetation moisture from precipitation, vapor flux with canopy,      !
         ! and shedding of excess moisture.                                                !
         !---------------------------------------------------------------------------------!
         old_veg_water = veg_water
         wtemp         = veg_water + ( intercepted_tot - wflxvc ) * dtvg
         !----- If this result will lead to negative solution, then reduce evaporation. ---!
         if (wtemp < 0.005 * leaf_maxwhc * stai) then
            wflxvc     = veg_water / dtvg + intercepted_tot
            veg_water  = 0.
            wshed      = 0.
            qwshed     = 0.
            dwshed     = 0.
         elseif (wtemp > leaf_maxwhc * stai) then
            !----- Shed water in excess of the leaf holding water capacity. ---------------!
            wshed      = ( wtemp - leaf_maxwhc * stai ) / dtvg
            qwshed     = wshed * tl2uint(veg_temp,veg_fliq)
            dwshed     = wshed * (veg_fliq * wdnsi + (1.0 - veg_fliq) * fdnsi)
            veg_water  = leaf_maxwhc * stai
         else
            wshed      = 0.
            qwshed     = 0.
            dwshed     = 0.
            veg_water  = wtemp
         end if
         !---------------------------------------------------------------------------------!



         !------ Find the associated latent heat flux from vegetation to canopy. ----------!
         if ( wflxvc >= 0.0 ) then
            qwflxvc = wflxvc * tq2enthalpy(veg_temp,1.0,.true.)
         else
            qwflxvc = wflxvc * tq2enthalpy(can_temp,1.0,.true.)
         end if
         transp     = transp_o
         qtransp    = transp * tq2enthalpy(veg_temp,1.0,.true.)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Update vegetation internal energy from radiative fluxes and sensible and   !
         ! latent heat transfer with canopy air.                                           !
         !---------------------------------------------------------------------------------!
         old_veg_energy = veg_energy
         veg_energy     = veg_energy                                                       &
                        + dtvg * ( rshort_v + rlong_v - hflxvc - qwflxvc - qtransp         &
                                 + qintercepted_tot - qwshed )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Update the exponential smoothing of growth respiration.  In case the value  !
         ! is becoming too small, flush it to zero to avoid floating point exceptions.     !
         !---------------------------------------------------------------------------------!
         growresp = growresp * (1. - gr_weight)                                            &
                  + gr_factor(nveg) * (gpp_o - leaf_resp_o - root_resp_o ) * gr_weight
         if (growresp < tiny_num) growresp = 0.0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Update the fluxes.                                                           !
         !---------------------------------------------------------------------------------!
         hflxvc_tot  = hflxvc_tot  + hflxvc  * ndtvegi
         wflxvc_tot  = wflxvc_tot  + wflxvc  * ndtvegi
         qwflxvc_tot = qwflxvc_tot + qwflxvc * ndtvegi
         transp_tot  = transp_tot  + transp  * ndtvegi
         qtransp_tot = qtransp_tot + qtransp * ndtvegi
         wshed_tot   = wshed_tot   + wshed   * ndtvegi
         qwshed_tot  = qwshed_tot  + qwshed  * ndtvegi
         dwshed_tot  = dwshed_tot  + dwshed  * ndtvegi
         cflxgc_tot  = cflxgc_tot  + (root_resp_o + 0.3 * growresp + het_resp_o) * ndtvegi
         cflxvc_tot  = cflxvc_tot  + (leaf_resp_o + 0.7 * growresp - gpp_o     ) * ndtvegi
         !----- Output GPP and plant respiration for output. ------------------------------!
         gpp_out     = gpp_out     + gpp_o   * dtl3_factor * ndtvegi
         plresp_out  = plresp_out  + (leaf_resp_o + root_resp_o + growresp)                &
                                   * dtl3_factor * ndtvegi
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Update vegetation temperature and liquid water fraction.                     !
         !---------------------------------------------------------------------------------!
         call uextcm2tl(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)
         !---------------------------------------------------------------------------------!
      end do vegloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Update enthalpy, CO2, and canopy mixing ratio.                                !
      !------------------------------------------------------------------------------------!
      can_enthalpy = can_enthalpy                                                          &
                   + dtl3ohcc * ( hflxsc + hflxgc + hflxvc_tot + qwflxsc + qwflxgc         &
                                - qdewgnd_tot + qwflxvc_tot + qtransp_tot + eflxac)
      can_shv      = can_shv                                                               &
                   + dtl3owcc * ( wflxsc + wflxgc - dewgnd_tot + wflxvc_tot + transp_tot   &
                                + wflxac)
      can_rvap     = can_shv / ( 1.0 - can_shv )
      can_co2    = can_co2  + dtl3occc * (cflxgc_tot + cflxvc_tot + cflxac)
      !------------------------------------------------------------------------------------!



      !----- Update the fluxes. -----------------------------------------------------------!
      hflxac_out       = hflxac_out       + hflxac                     * dtl3_factor
      wflxac_out       = wflxac_out       + wflxac                     * dtl3_factor
      qwflxac_out      = qwflxac_out      + qwflxac                    * dtl3_factor
      eflxac_out       = eflxac_out       + eflxac                     * dtl3_factor
      cflxac_out       = cflxac_out       + cflxac                     * dtl3_factor
      hflxgc_out       = hflxgc_out       + ( hflxgc  + hflxsc  )      * dtl3_factor
      wflxgc_out       = wflxgc_out       + ( wflxgc  + wflxsc  - dewgnd_tot  )            &
                                                                       * dtl3_factor
      qwflxgc_out      = qwflxgc_out      + ( qwflxgc + qwflxsc - qdewgnd_tot )            &
                                                                       * dtl3_factor
      hflxvc_out       = hflxvc_out       + hflxvc_tot                 * dtl3_factor
      wflxvc_out       = wflxvc_out       + wflxvc_tot                 * dtl3_factor
      qwflxvc_out      = qwflxvc_out      + qwflxvc_tot                * dtl3_factor
      transp_out       = transp_out       + transp_tot                 * dtl3_factor
      qtransp_out      = qtransp_out      + qtransp_tot                * dtl3_factor
      intercepted_out  = intercepted_out  + intercepted_tot            * dtl3_factor
      qintercepted_out = qintercepted_out + qintercepted_tot           * dtl3_factor
      wshed_out        = wshed_out        + wshed_tot                  * dtl3_factor
      qwshed_out       = qwshed_out       + qwshed_tot                 * dtl3_factor
      throughfall_out  = throughfall_out  + throughfall_tot            * dtl3_factor
      qthroughfall_out = qthroughfall_out + qthroughfall_tot           * dtl3_factor
      resphet_out      = resphet_out      + het_resp_o                 * dtl3_factor
      !------------------------------------------------------------------------------------!

   else
      !------------------------------------------------------------------------------------!
      !    No vegetation, or vegetation buried in snow...                                  !
      !------------------------------------------------------------------------------------!

      !----- Decay the growth respiration. ------------------------------------------------!
      growresp   = growresp * (1. - gr_weight)
      if (growresp < tiny_num) growresp = 0.0
      !------------------------------------------------------------------------------------!


      !----- Update the canopy prognostic variables. --------------------------------------!
      can_enthalpy = can_enthalpy                                                          &
                   + dtl3ohcc * (hflxsc + hflxgc + qwflxsc + qwflxgc - qdewgnd_tot + eflxac)
      can_shv      = can_shv + dtl3owcc * (wflxsc + wflxgc - dewgnd_tot + wflxac)
      can_rvap     = can_shv / ( 1.0 - can_shv )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel       !
      ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.      !
      !------------------------------------------------------------------------------------!
      cflxgc       = het_resp_o
      cflxvc_tot   = 0.
      can_co2      = can_co2      + dtl3occc * (cflxgc + cflxvc_tot + cflxac)
      !------------------------------------------------------------------------------------!



      !----- Update the fluxes. -----------------------------------------------------------!
      hflxac_out       = hflxac_out       + hflxac                     * dtl3_factor
      wflxac_out       = wflxac_out       + wflxac                     * dtl3_factor
      qwflxac_out      = qwflxac_out      + qwflxac                    * dtl3_factor
      eflxac_out       = eflxac_out       + eflxac                     * dtl3_factor
      cflxac_out       = cflxac_out       + cflxac                     * dtl3_factor
      hflxgc_out       = hflxgc_out       + ( hflxgc  + hflxsc      )  * dtl3_factor
      wflxgc_out       = wflxgc_out       + ( wflxgc  + wflxsc  - dewgnd_tot  )            &
                                          * dtl3_factor
      qwflxgc_out      = qwflxgc_out      + ( qwflxgc + qwflxsc - qdewgnd_tot )            &
                                          * dtl3_factor
      throughfall_out  = throughfall_out  + throughfall_tot            * dtl3_factor
      qthroughfall_out = qthroughfall_out + qthroughfall_tot           * dtl3_factor
      resphet_out      = resphet_out      + het_resp_o                 * dtl3_factor
      !------------------------------------------------------------------------------------!

   end if

   return
end subroutine leaf3_canopy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine initialises a couple of canopy air space variables for a given     !
! patch, before the time step iteration loop.                                              !
!------------------------------------------------------------------------------------------!
subroutine leaf3_can_diag(ip,can_theta,can_theiv,can_vpdef,can_rvap,leaf_class,can_prss    &
                         ,initial)
   use leaf_coms , only : atm_prss      & ! intent(in)
                        , atm_theta     & ! intent(in)
                        , atm_shv       & ! intent(in)
                        , atm_temp_zcan & ! intent(in)
                        , atm_enthalpy  & ! intent(in)
                        , geoht         & ! intent(in)
                        , veg_ht        & ! intent(in)
                        , can_shv       & ! intent(out)
                        , can_rsat      & ! intent(out)
                        , can_rhv       & ! intent(out)
                        , can_depth     & ! intent(inout)
                        , can_depth_min & ! intent(in)
                        , can_temp      & ! intent(out)
                        , can_exner     & ! intent(inout)
                        , can_enthalpy  & ! intent(inout)
                        , can_rhos      & ! intent(inout)
                        , can_cp        ! ! intent(inout)
   use rconstants, only : cpdry         & ! intent(in)
                        , cph2o         & ! intent(in)
                        , ep            & ! intent(in)
                        , p00i          & ! intent(in)
                        , rocp          ! ! intent(in)
   use therm_lib , only : reducedpress  & ! function
                        , rslif         & ! function
                        , idealdenssh   & ! function
                        , thetaeiv      & ! function
                        , vpdefil       & ! function
                        , press2exner   & ! function
                        , extheta2temp  & ! function
                        , extemp2theta  & ! function
                        , tq2enthalpy   & ! function
                        , hq2temp       ! ! function

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)      :: ip
   real   , intent(inout)   :: can_theta
   real   , intent(out)     :: can_theiv
   real   , intent(out)     :: can_vpdef
   real   , intent(in)      :: can_rvap
   real   , intent(in)      :: leaf_class
   real   , intent(inout)   :: can_prss
   logical, intent(in)      :: initial
   !----- Local variables. ----------------------------------------------------------------!
   integer                  :: nveg
   !---------------------------------------------------------------------------------------!




   !----- Define the canopy depth. --------------------------------------------------------!
   if (initial) then
      nveg      = nint(leaf_class)
      if (nveg == 0 .or. ip == 1) then
         can_depth = can_depth_min
      else
         can_depth = max(can_depth_min,veg_ht(nveg))
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !----- Set initial value of specific humidity ------------------------------------------!
   can_shv     = can_rvap / (can_rvap + 1.)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update canopy air pressure and the Exner function here.  Canopy air pressure is    !
   ! assumed to remain constant during one LEAF full timestep, which means that heat       !
   ! equals to change in enthalpy.  Between time steps, atmospheric pressure changes so    !
   ! enthalpy is no longer conserved.  Potential temperature and equivalent potential      !
   ! temperature, on the other hand, are conserved if no heat happens, which is the case   !
   ! between successive calls.  Therefore, we track theta and theta_eiv.                   !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      can_prss  = reducedpress(atm_prss,atm_theta,atm_shv,geoht,can_theta,can_shv          &
                              ,can_depth)
      can_exner = press2exner(can_prss)

      !----- Also, find the specific enthalpy of the air aloft. ---------------------------!
      atm_temp_zcan = extheta2temp(can_exner,atm_theta)
      atm_enthalpy  = tq2enthalpy(atm_temp_zcan,atm_shv,.true.)
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     If this is the initial step, we initialise the canopy air specific enthalpy.      !
   ! Otherwise, we update temperature and potential temperature from enthalpy.             !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      can_temp     = extheta2temp(can_exner,can_theta)
      can_enthalpy = tq2enthalpy(can_temp,can_shv,.true.)
   else
      can_temp    = hq2temp(can_enthalpy,can_shv,.true.)
      can_theta   = extemp2theta(can_exner,can_temp)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Even though this leads to small violations of the ideal gas law, we are also      !
   ! going to keep the canopy air density constant over one LEAF time step.  This is to    !
   ! avoid that non-linearities sprout and cause problems to the code.                     !
   !---------------------------------------------------------------------------------------!
   if (initial) then
      can_rhos = idealdenssh (can_prss,can_temp,can_shv)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the canopy air saturation mixing ratio and the relative humidity.  These are !
   ! used to decide whether the canopy air space can hold more water or not.               !
   !---------------------------------------------------------------------------------------!
   can_rsat = rslif(can_prss,can_temp)
   can_rhv  = can_rvap * (ep + can_rsat) / (can_rsat * (ep + can_rvap))
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the canopy air space ice-vapour equivalent potential temperature.  This is   !
   ! currently a diagnostic variable only, but it should become the main variable if we    !
   ! ever switch to foggy canopy air space.                                                !
   !---------------------------------------------------------------------------------------!
   can_theiv = thetaeiv(can_theta,can_prss,can_temp,can_rvap,can_rvap)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the vapour pressure deficit.  This is a diagnostic variable.                 !
   !---------------------------------------------------------------------------------------!
   can_vpdef = vpdefil(can_prss,can_temp,can_shv,.true.)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the canopy air space specific heat at constant pressure.                     !
   !---------------------------------------------------------------------------------------!
   can_cp    = (1.0 - can_shv) * cpdry + can_shv * cph2o
   !---------------------------------------------------------------------------------------!



   return
end subroutine leaf3_can_diag
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine initialises a couple of leaf-level variables for a given patch,    !
! before the time step iteration loop.                                                     !
!------------------------------------------------------------------------------------------!
subroutine leaf3_veg_diag(veg_energy,veg_water,veg_hcap)
   use leaf_coms , only : can_temp         & ! intent(in)
                        , min_patch_area   & ! intent(in)
                        , veg_temp         & ! intent(out)
                        , veg_fliq         & ! intent(out)
                        , resolvable       ! ! intent(in)
   use rconstants, only : t3ple            ! ! intent(in)
   use therm_lib , only : uextcm2tl        ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(inout) :: veg_energy
   real   , intent(inout) :: veg_water
   real   , intent(in)    :: veg_hcap
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find leaf temperature and liquid fraction of water, but only if we are solving    !
   ! leaves.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (resolvable) then
      call uextcm2tl(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)
   else
      veg_temp = can_temp
      if (veg_temp > t3ple) then
         veg_fliq = 1.0
      elseif (veg_temp < t3ple) then
         veg_fliq = 0.0
      else
         veg_fliq = 0.5
      end if
      veg_energy = veg_hcap * veg_temp
      veg_water  = 0.0
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_veg_diag
!==========================================================================================!
!==========================================================================================!
