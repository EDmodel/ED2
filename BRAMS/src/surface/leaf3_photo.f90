!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the photosynthesis scheme (Farquar and Leuning).  The   !
! implementation in LEAF-3 is very similar to ED-2.2 (Longo 2013) for most of it, but a    !
! few things such as the water stress and the division between sun and shade leaves is     !
! modelled as in CLM-4 (Oleson et al. 2010).                                               !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine leaf3_photo_driver(mzg,soil_water,soil_text,leaf_class,rshort,veg_lai,veg_tai   &
                             ,can_prss,can_rvap,can_co2,stom_condct)
   use mem_leaf     , only : slz                    & ! intent(in)
                           , isfcl                  ! ! intent(in)
   use leaf_coms    , only : dtvg                   & ! intent(in)
                           , supersat_ok            & ! intent(in)
                           , slzt                   & ! intent(in)
                           , dslz                   & ! intent(in)
                           , tai_min                & ! intent(in)
                           , gsw_max                & ! intent(in)
                           , transp_max             & ! intent(in)
                           , kroot                  & ! intent(in)
                           , veg_ht                 & ! intent(in)
                           , brad                   & ! intent(in)
                           , srad                   & ! intent(in)
                           , btlo                   & ! intent(in)
                           , stlo                   & ! intent(in)
                           , bthi                   & ! intent(in)
                           , sthi                   & ! intent(in)
                           , bvpd                   & ! intent(in)
                           , svpd                   & ! intent(in)
                           , bsmp                   & ! intent(in)
                           , ssmp                   & ! intent(in)
                           , gsw_0                  & ! intent(in)
                           , lai_ss                 & ! intent(in)
                           , par_l_ss               & ! intent(in)
                           , vm0_ss                 & ! intent(in)
                           , rd0_ss                 & ! intent(in)
                           , stom_side              & ! intent(in)
                           , slmsts                 & ! intent(in)
                           , slpots                 & ! intent(in)
                           , slbs                   & ! intent(in)
                           , soilwp                 & ! intent(in)
                           , snowfac                & ! intent(in)
                           , psifc                  & ! intent(in)
                           , psiwp                  & ! intent(in)
                           , gbw                    & ! intent(in)
                           , can_rhos               & ! intent(in)
                           , can_rhv                & ! intent(in)
                           , can_shv                & ! intent(in)
                           , veg_temp               & ! intent(in)
                           , soil_tempk             & ! intent(in)
                           , soil_fracliq           & ! intent(in)
                           , gpp_ss                 & ! intent(out)
                           , leaf_resp_ss           & ! intent(out)
                           , transp_ss              & ! intent(out)
                           , gpp_o                  & ! intent(out)
                           , leaf_resp_o            & ! intent(out)
                           , transp_o               ! ! intent(out)
   use rconstants   , only : twothirds              & ! intent(in)
                           , ep                     & ! intent(in)
                           , huge_num               & ! intent(in)
                           , wdns                   ! ! intent(in)
   use therm_lib    , only : qslif                  & ! function
                           , eslif                  ! ! function
   use leaf3_physiol, only : leaf3_farquhar_leuning ! ! subroutine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: mzg
   real(kind=4), dimension(mzg), intent(in)    :: soil_water
   real(kind=4), dimension(mzg), intent(in)    :: soil_text
   real(kind=4)                , intent(in)    :: rshort
   real(kind=4)                , intent(in)    :: leaf_class
   real(kind=4)                , intent(in)    :: veg_lai
   real(kind=4)                , intent(in)    :: veg_tai
   real(kind=4)                , intent(in)    :: can_prss
   real(kind=4)                , intent(in)    :: can_rvap
   real(kind=4)                , intent(in)    :: can_co2
   real(kind=4)                , intent(inout) :: stom_condct
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: nveg
   integer                                     :: nsoil
   integer                                     :: k
   integer                                     :: lrl
   real(kind=4)                                :: mcheight
   real(kind=4)                                :: wgpfrac
   real(kind=4)                                :: psiplusz
   real(kind=4)                                :: fs_open
   real(kind=4)                                :: fs_closed
   real(kind=4)                                :: leaf_par
   real(kind=4)                                :: lint_shv
   real(kind=4)                                :: lint_pvap
   real(kind=4)                                :: lint_rvap
   real(kind=4)                                :: lsfc_pvap
   real(kind=4)                                :: lsfc_rvap
   real(kind=4)                                :: a_open
   real(kind=4)                                :: a_closed
   real(kind=4)                                :: gsw_open
   real(kind=4)                                :: gsw_closed
   real(kind=4)                                :: gsw_inf
   real(kind=4)                                :: gsw
   real(kind=4)                                :: gleaf_open
   real(kind=4)                                :: gleaf_closed
   real(kind=4)                                :: rvap_gradient
   real(kind=4)                                :: shv_gradient
   real(kind=4)                                :: smpot
   real(kind=4)                                :: smpot_wet
   real(kind=4)                                :: vpdef
   real(kind=4)                                :: ftlo
   real(kind=4)                                :: fthi
   real(kind=4)                                :: frad
   real(kind=4)                                :: fsmp
   real(kind=4)                                :: fvpd
   real(kind=4)                                :: slai
   real(kind=4)                                :: stai
   real(kind=4)                                :: transp_wilt
   real(kind=4)                                :: transp_test
   real(kind=4)                                :: gpp_test
   real(kind=4)                                :: rfac
   real(kind=4)                                :: available_water
   !---------------------------------------------------------------------------------------!


   !----- Aliases to some useful variables. -----------------------------------------------$
   nveg     = nint(leaf_class)
   lrl      = kroot(nveg)
   mcheight = twothirds * veg_ht(nveg)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the "perceived" leaf area and tree area indices, dependending on snow cover. !
   !---------------------------------------------------------------------------------------!
   slai = veg_lai * ( 1.0 - snowfac )
   stai = veg_tai * ( 1.0 - snowfac )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find the available water for this patch.                                         !
   !---------------------------------------------------------------------------------------!
   fs_open         =    0.0
   smpot_wet       = -200.0
   available_water =    0.0
   do k=mzg,lrl,-1
      nsoil = nint(soil_text(k))


      !----- Find the potential for this layer. -------------------------------------------!
      wgpfrac  = min(soil_water(k)/slmsts(nsoil), 1.0)
      smpot    = slpots(nsoil) / wgpfrac ** slbs(nsoil)
      psiplusz = slzt(k) - mcheight + smpot
      !------------------------------------------------------------------------------------!


      !----- Update the wettest matric potential. -----------------------------------------!
      if (smpot > smpot_wet) smpot_wet = smpot
      !------------------------------------------------------------------------------------!



      !----- Integrate available water. ---------------------------------------------------!
      available_water = available_water                                                    &
                      + dslz(k) * wdns                                                     &
                      * max(0.0, soil_fracliq(k) * (soil_water(k)-soilwp(nsoil)))
      !------------------------------------------------------------------------------------!



      !----- Find the available water factor for this layer. ------------------------------!
      fs_open = fs_open                                                    &
              + max(0.,(psiplusz - psiwp(nsoil)) / (psifc(nsoil) - psiwp(nsoil))) &
              * soil_fracliq(k) * dslz(k)
      !------------------------------------------------------------------------------------!
   end do
   fs_open   = max( 0., min(1., fs_open  / abs(slz(lrl)) ) )
   fs_closed = 1.0 - fs_open
   !----- Find the maximum transpiration allowed. -----------------------------------------!
   if (can_rhv >= 1.0 .and. (.not. supersat_ok)) then
      transp_wilt = 0.0
   else
      select case (isfcl)
      case (1,2)
         transp_wilt  = min(available_water / dtvg, transp_max)
      case (4,5)
         transp_wilt  = available_water / dtvg
      end select
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!





   !---- Find specific humidity at the intercellular space. -------------------------------!
   lint_shv      = qslif(can_prss,veg_temp)
   lint_rvap     = lint_shv / (1.0 - lint_shv)
   lint_pvap     = eslif(veg_temp)
   lsfc_rvap     = ( gbw * can_rvap + stom_condct * lint_rvap ) / ( gbw + stom_condct )
   lsfc_pvap     = lsfc_rvap * can_prss / (ep + lsfc_rvap)
   rvap_gradient = max( lint_rvap - can_rvap , 0.0 )
   shv_gradient  = max( lint_shv  - can_shv  , 0.0 )
   vpdef         = max( lint_pvap - lsfc_pvap, 0.0 )
   rfac          = ( 1.0 + can_rvap ) * ( 1.0 + lint_rvap )
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Decide whether to use the old (transpiration only) or the new scheme.             !
   !---------------------------------------------------------------------------------------!
   select case (isfcl)
   case (1,2)
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     Original scheme, no photosynthesis is provided.                                !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !      Use a combination of 5 environmental factors that may control stomatal        !
      ! conductance to estimate the target value for this new time step.                   !
      !------------------------------------------------------------------------------------!
      ftlo    = 1. + exp(-stlo * (veg_temp  - btlo))
      fthi    = 1. + exp(-sthi * (veg_temp  - bthi))
      frad    = 1. + exp(-srad * (rshort    - brad))
      fsmp    = 1. + exp(-ssmp * (smpot_wet - bsmp))
      fvpd    = 1. + exp(-svpd * (vpdef     - bvpd))
      if (slai > tai_min) then
         gsw_inf = can_rhos * gsw_max(nveg) / (ftlo * fthi * frad * fvpd * fsmp) / slai
      else
         gsw_inf = 0.0
      end if
      !------------------------------------------------------------------------------------!

      !----- 15-minute response time for stomatal conductance (must have dtvg <= 900.). ---!
      gsw = max(gsw_0(nveg), stom_condct + dtvg * (gsw_inf - stom_condct) / max(dtvg,900.))
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Limit maximum transpiration to be <= transp_max and less than the maximum water !
      ! available by decreasing stomatal conductance if necessary.                         !
      !------------------------------------------------------------------------------------!
      transp_test = stom_side(nveg) * gbw * gsw * slai * shv_gradient * rfac / ( gbw + gsw )
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In the originial LEAF scheme, transpiration is not solved as a direct function  !
      ! of photosynthesis.  Since we don't have intercellular CO[2], we just assign zero.  !
      !------------------------------------------------------------------------------------!
      gpp_test    = 0.0
      leaf_resp_o = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
   case (4)
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !     New scheme, using the Farquhar-Leuning method, based on ED.                    !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Loop over sun and shade, and find the photosynthetic rates per unit of leaf.   !
      ! Make sure that incredibly sparse vegetation is not called.                         !
      !------------------------------------------------------------------------------------!
      do k=1,2
         if (lai_ss(k) > tai_min) then
            !---- Find energy per unit leaf. ----------------------------------------------!
            leaf_par = par_l_ss(k) / lai_ss(k)
            !------------------------------------------------------------------------------!



            !----- Call the photosynthesis scheme. ----------------------------------------!
            call leaf3_farquhar_leuning(can_prss,can_rhos,can_shv,can_co2,leaf_class       &
                                       ,leaf_par,veg_temp,lint_shv,vm0_ss(k),rd0_ss(k)     &
                                       ,a_open,a_closed,gsw_open,gsw_closed                &
                                       ,leaf_resp_ss(k))
            !------------------------------------------------------------------------------!


            !----- Scale leaf respiration by LAI. -----------------------------------------!
            leaf_resp_ss(k) = lai_ss(k) * leaf_resp_ss(k)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Find the stomatal conductance scaled by the wilting factor.             !
            !------------------------------------------------------------------------------!
            gpp_ss(k) = lai_ss(k) * ( fs_open * a_open + fs_closed * a_closed )            &
                      + leaf_resp_ss(k)
            gpp_ss(k) = max(0.,gpp_ss(k))
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Find the stomatal conductance and transpiration.                        !
            !------------------------------------------------------------------------------!
            gleaf_open   = gbw * gsw_open   / ( gbw + gsw_open   )
            gleaf_closed = gbw * gsw_closed / ( gbw + gsw_closed )
            transp_ss(k) = stom_side(nveg) * lai_ss(k) * max(0.,shv_gradient) * rfac       &
                         * ( gleaf_open * fs_open + gleaf_closed * fs_closed )
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     Nothing happens.                                                         !
            !------------------------------------------------------------------------------!
            leaf_resp_ss(k) = 0.0
            gpp_ss      (k) = 0.0
            transp_ss   (k) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Total transpiration is the sum of both layers.  This value may be modified      !
      ! later if transpiration  would extract too much water.                              !
      !------------------------------------------------------------------------------------!
      transp_test = sum(transp_ss)
      gpp_test    = sum(gpp_ss   )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Leaf respiration doesn't depend on whether the maximum allowed transpiration.  !
      !------------------------------------------------------------------------------------!
      leaf_resp_o = sum(leaf_resp_ss)
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Scale both the GPP and transpiration.                                             !
   !---------------------------------------------------------------------------------------!
   if (transp_test <= 0. .or. transp_wilt <= 0.) then
      !----- Transpiration would be negative.  Set it to zero. ----------------------------!
      gpp_o    = 0.0
      transp_o = 0.0
      !------------------------------------------------------------------------------------!

   elseif (transp_test > transp_wilt) then
      !----- Transpiration would require more water than available.  Reduce both. ---------!
      gpp_o    = (transp_wilt / transp_test) * gpp_test
      transp_o = transp_wilt
      !------------------------------------------------------------------------------------!

   else
      !----- Transpiration is fine. -------------------------------------------------------!
      gpp_o    = gpp_test
      transp_o = transp_test
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Find the net stomatal conductance.  Make sure that the singularity case is taken !
   ! care of as well.                                                                      !
   !---------------------------------------------------------------------------------------!
   if (transp_o == gbw * stom_side(nveg) * slai * shv_gradient * rfac) then
      stom_condct = huge_num
      continue
   else
      stom_condct = gbw * transp_o                                                         &
                  / ( stom_side(nveg) * gbw * slai * shv_gradient * rfac - transp_o )
   end if
   gsw = stom_condct
   !---------------------------------------------------------------------------------------!

   return
end subroutine leaf3_photo_driver
!==========================================================================================!
!==========================================================================================!
