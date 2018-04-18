!==========================================================================================!
!==========================================================================================!
!    This module is a library with several allometric relationships.                       !
!------------------------------------------------------------------------------------------!
module allometry

   implicit none

   contains
   !=======================================================================================!
   !=======================================================================================!
   real function h2dbh(h,ipft)

      use pft_coms    , only : is_tropical & ! intent(in)
                             , is_liana    & ! intent(in)
                             , b1Ht        & ! intent(in), lookup table
                             , b2Ht        & ! intent(in), lookup table
                             , hgt_ref     ! ! intent(in)
      use ed_misc_coms, only : iallom      ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Select which type of model we are running.                                    !
      !------------------------------------------------------------------------------------!
      !----- Size- and age-structure (typical ED model). ----------------------------------!
      if (is_liana(ipft)) then
         !----- Inverse Weibull function, similar to Poorter et al. (2006). ---------------!
         h2dbh =  ( log(hgt_ref(ipft) / ( hgt_ref(ipft) - h ) ) / b1Ht(ipft) )             &
               ** ( 1.0 / b2Ht(ipft) )
         !---------------------------------------------------------------------------------!
      elseif (is_tropical(ipft)) then
         select case (iallom)
            case (0,1)
               !----- Default ED-2.1 allometry. -------------------------------------------!
               h2dbh = exp((log(h)-b1Ht(ipft))/b2Ht(ipft))
            case default
               !----- Poorter et al. (2006) allometry. ------------------------------------!
               h2dbh =  ( log(hgt_ref(ipft) / ( hgt_ref(ipft) - h ) ) / b1Ht(ipft) )       &
                  ** ( 1.0 / b2Ht(ipft) )
         end select
      else ! Temperate
         h2dbh = log(1.0-(h-hgt_ref(ipft))/b1Ht(ipft))/b2Ht(ipft)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function h2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2h(ipft, dbh)
      use pft_coms     , only : is_tropical & ! intent(in)
                              , is_liana    & ! intent(in)
                              , dbh_crit    & ! intent(in)
                              , b1Ht        & ! intent(in)
                              , b2Ht        & ! intent(in)
                              , hgt_ref     & ! intent(in)
                              , hgt_max     ! ! intent(in)
      use ed_misc_coms , only : iallom      & ! intent(in)
                              , ibigleaf    ! ! intent(in)
      use ed_state_vars, only : patchtype   ! ! structure
      use consts_coms  , only : lnexp_max   & ! intent(in)
                              , lnexp_min   ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      integer       , intent(in) :: ipft
      real          , intent(in) :: dbh
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: lnexp
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Select which type of model we are running.                                    !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
         case (0)
            !----- Size- and age-structure (typical ED model). ----------------------------!
            if (is_liana(ipft)) then
               !----- Weibull function, similar to Poorter et al. (2006). -----------------!
               mdbh = min(dbh,dbh_crit(ipft))
               dbh2h = hgt_ref(ipft) * (1. - exp(-b1Ht(ipft) * mdbh ** b2Ht(ipft)))
            elseif (is_tropical(ipft)) then
               mdbh = min(dbh,dbh_crit(ipft))
               select case (iallom)
               case (0,1)
                  !----- Default ED-2.1 allometry. ----------------------------------------!
                  dbh2h = exp (b1Ht(ipft) + b2Ht(ipft) * log(mdbh) )
               case default
                  !----- Poorter et al. (2006) allometry. ---------------------------------!
                  lnexp = max(lnexp_min,min(lnexp_max,b1Ht(ipft) * mdbh ** b2Ht(ipft)))
                  dbh2h = hgt_ref(ipft) * (1. - exp(-lnexp))
               end select
            else !----- Temperate PFT allometry. ------------------------------------------!
               lnexp = max(lnexp_min,min(lnexp_max,b2Ht(ipft) * dbh))
               dbh2h = hgt_ref(ipft) + b1Ht(ipft) * (1.0 - exp(lnexp))
            end if

         case (1)
            !------------------------------------------------------------------------------!
            !     Big-leaf version of ED. DBH is not really meaningful, but in the big-leaf!
            ! model the typical allometry doesn't really make sense so we impose maximum   !
            ! height.                                                                      !
            !------------------------------------------------------------------------------!
            dbh2h = hgt_max(ipft)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function dbh2h
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     Function that finds Bdead from DBH.                                               !
   !---------------------------------------------------------------------------------------!
   real function size2bd(dbh,hite,ipft)

      use pft_coms    , only : C2B         & ! intent(in)
                             , dbh_adult   & ! intent(in)
                             , dbh_crit    & ! intent(in)
                             , b1Bs_small  & ! intent(in)
                             , b2Bs_small  & ! intent(in)
                             , b1Bs_large  & ! intent(in)
                             , b2Bs_large  & ! intent(in)
                             , is_grass    & ! intent(in)
                             , is_tropical & ! intent(in)
                             , is_liana    ! ! intent(in)
      use ed_misc_coms, only : igrass      & ! intent(in)
                             , iallom      ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: hite
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Structural biomass depends on the allometry and the size.                      !
      !------------------------------------------------------------------------------------!
      if (igrass == 1 .and. is_grass(ipft)   ) then
         size2bd = 0.0
      else if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
         !----- Decide parameters based on seedling/adult size. ---------------------------!
         if (dbh <= dbh_adult(ipft)) then
            size2bd = b1Bs_small(ipft) / C2B * (dbh*dbh*hite) ** b2Bs_small(ipft)
         else
            size2bd = b1Bs_large(ipft) / C2B * (dbh*dbh*hite) ** b2Bs_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      else 
         !----- Decide parameters baded on the maximum height. ----------------------------!
         if (dbh <= dbh_crit(ipft)) then
            size2bd = b1Bs_small(ipft) / C2B * dbh ** b2Bs_small(ipft)
         else
            size2bd = b1Bs_large(ipft) / C2B * dbh ** b2Bs_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2bd
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the diameter at the breast height given the structural     !
   ! (dead) biomass and the PFT type.  In addition to these inputs, it is also required to !
   ! give a first guess for DBH, which can be the previous value, for example.   This is   !
   ! only used for the new allometry, because we must solve an iterative root-finding      !
   ! method.                                                                               !
   !---------------------------------------------------------------------------------------!
   real function bd2dbh(ipft, bdead)
      use pft_coms    , only : b1Bs_small  & ! intent(in)
                             , b2Bs_small  & ! intent(in)
                             , b1Bs_large  & ! intent(in)
                             , b2Bs_large  & ! intent(in)
                             , bdead_adult & ! intent(in)
                             , bdead_crit  & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , bdead_lut   & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut & ! intent(out)
                             , is_tropical & ! intent(in)
                             , is_liana    & ! intent(in)
                             , C2B         ! ! intent(in)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in) :: ipft      ! PFT type                            [         ---]
      real   , intent(in) :: bdead     ! Structural (dead) biomass           [   kgC/plant]
      !----- Local variables. -------------------------------------------------------------!
      integer             :: ilwr        ! Lower index of the lookup table
      integer             :: iupr        ! Upper index of the lookup table
      real                :: finterp     ! Interpolation factor
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide which coefficients to use based on the critical bdead.                   !
      !------------------------------------------------------------------------------------!
      if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft)) ) then
         if (bdead <= bdead_lut(1,ipft)) then
            !----- Use the look-up table to find the best dbh. ----------------------------!
            bd2dbh = dbh_lut(1,ipft) * bdead / bdead_lut(1,ipft)
            !------------------------------------------------------------------------------!
         else if (bdead >= bdead_crit(ipft)) then
            !------------------------------------------------------------------------------!
            !     Bdead is above critical value, height is known.  It is also very likely  !
            ! that bdead is above minimum adult value, but just to be safe, we also check  !
            ! it.                                                                          !
            !------------------------------------------------------------------------------!
            if (bdead < bdead_adult(ipft)) then
               bd2dbh =  ( C2B * bdead                                                     &
                         / ( b1Bs_small(ipft) * hgt_max(ipft)**b2Bs_small(ipft) ))         &
                      ** ( 1. / (2. * b2Bs_small(ipft) ) )
            else
               bd2dbh =  ( C2B * bdead                                                     &
                         / ( b1Bs_large(ipft) * hgt_max(ipft)**b2Bs_large(ipft) ))         &
                      ** ( 1. / (2. * b2Bs_large(ipft) ) )
            end if
            !------------------------------------------------------------------------------!
         else
            !----- Use the look-up table to find the best dbh. ----------------------------!
            le_mask_lut(:) = bdead <= bdead_lut(:,ipft)
            ge_mask_lut(:) = bdead >= bdead_lut(:,ipft)
            ilwr  = maxloc (bdead_lut(:,ipft),dim=1,mask=le_mask_lut)
            iupr  = minloc (bdead_lut(:,ipft),dim=1,mask=ge_mask_lut)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      In case ilwr and iupr are the same, we have an exact estimate.  Other-  !
            ! wise, use log-linear interpolation.                                          !
            !------------------------------------------------------------------------------!
            if (ilwr == iupr) then
               bd2dbh = dbh_lut(ilwr,ipft)
            else
               finterp = log( dbh_lut  (iupr,ipft) / dbh_lut  (ilwr,ipft))                 &
                       * log( bdead                / bdead_lut(ilwr,ipft))                 &
                       / log( bdead_lut(iupr,ipft) / bdead_lut(ilwr,ipft))
               finterp = max(lnexp_min,min(lnexp_max,finterp))
               bd2dbh  = dbh_lut(ilwr,ipft) * exp(finterp)
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      else
         !----- Bdead is not a direct function of height, simple inversion is sufficient. -!
         if (bdead <= bdead_crit(ipft)) then
            bd2dbh = (bdead / b1Bs_small(ipft) * C2B)**(1.0/b2Bs_small(ipft))
         else
            bd2dbh = (bdead / b1Bs_large(ipft) * C2B)**(1.0/b2Bs_large(ipft))
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function bd2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function computes the maximum leaf biomass [kgC/m2] given the size (either  !
   ! height or DBH).  Trees would always use DBH as the starting point, but grasses may    !
   ! use DBH (old style) or height (new style).  This replaces dbh2bl and h2bl with a      !
   ! single generic function that should be used by all plants.                            !
   !---------------------------------------------------------------------------------------!
   real function size2bl(dbh,hite,ipft)
      use pft_coms     , only : dbh_crit       & ! intent(in)
                              , dbh_adult      & ! intent(in)
                              , C2B            & ! intent(in)
                              , b1Bl_small     & ! intent(in)
                              , b2Bl_small     & ! intent(in)
                              , b1Bl_large     & ! intent(in)
                              , b2Bl_large     & ! intent(in)
                              , is_liana       & ! intent(in)
                              , is_grass       & ! intent(in)
                              , is_tropical    & ! intent(in)
                              , liana_dbh_crit ! ! intent(in)
      use ed_misc_coms , only : igrass         & ! intent(in)
                              , iallom         ! ! intent(in)
      use ed_state_vars, only : patchtype      ! ! structure

      !----- Arguments --------------------------------------------------------------------!
      real          , intent(in) :: dbh
      real          , intent(in) :: hite
      integer       , intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                       :: mdbh
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Get the DBH, or potential DBH in case of grasses.                              !
      !------------------------------------------------------------------------------------!
      if (igrass == 1 .and. is_grass(ipft)) then 
         !---- Use height for new grasses. ------------------------------------------------!
         mdbh   = min(dbh,h2dbh(hite,ipft))
      elseif (is_liana(ipft)) then
         mdbh   = min(dbh,liana_dbh_crit)
      else
         !---- Use dbh for trees. ---------------------------------------------------------!
         mdbh   = min(dbh,dbh_crit(ipft))
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find leaf biomass depending on the tree size.                                  !
      !------------------------------------------------------------------------------------!
      if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
         !----- Decide parameters based on seedling/adult size. ---------------------------!
         if (mdbh < dbh_adult(ipft)) then
            size2bl = b1Bl_small(ipft) / C2B * (mdbh*mdbh*hite) ** b2Bl_small(ipft)
         else
            size2bl = b1Bl_large(ipft) / C2B * (mdbh*mdbh*hite) ** b2Bl_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      else 
         !----- Decide parameters based on seedling/adult size. ---------------------------!
         if (mdbh < dbh_adult(ipft)) then
            size2bl = b1Bl_small(ipft) / C2B * mdbh ** b2Bl_small(ipft)
         else
            size2bl = b1Bl_large(ipft) / C2B * mdbh ** b2Bl_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      return
   end function size2bl
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine finds height given the biomass of living tissues, assuming       !
   ! minimum sapwood biomass given DBH.                                                    !
   !---------------------------------------------------------------------------------------!
   real function ba2h(balive,ipft)
      use pft_coms    , only : hgt_min     & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , balive_crit & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , balive_lut  & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut ! ! intent(out)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: balive      ! Live biomass
      integer, intent(in) :: ipft        ! PFT
      !----- Local variables. -------------------------------------------------------------!
      integer             :: ilwr        ! Lower index of the lookup table
      integer             :: iupr        ! Upper index of the lookup table
      real                :: dbh         ! Best guess for DBH.
      real                :: finterp     ! Interpolation factor
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    In case there is more biomass than balive_crit, assume that the plant is at     !
      ! maximum height.                                                                    !
      !------------------------------------------------------------------------------------!
      if (balive >= balive_crit(ipft)) then
         ba2h        = hgt_max(ipft)
      elseif (balive <= balive_lut(1,ipft)) then
         !----- Do not let height be less than minimum. -----------------------------------!
         ba2h        = hgt_min(ipft)
         !---------------------------------------------------------------------------------!
      else
         !----- Use the look-up table to find the best dbh. -------------------------------!
         le_mask_lut(:) = balive <= balive_lut(:,ipft)
         ge_mask_lut(:) = balive >= balive_lut(:,ipft)
         ilwr  = maxloc (balive_lut(:,ipft),dim=1,mask=le_mask_lut)
         iupr  = minloc (balive_lut(:,ipft),dim=1,mask=ge_mask_lut)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      In case ilwr and iupr are the same, we have an exact estimate.  Otherwise, !
         ! use log-linear interpolation.                                                   !
         !---------------------------------------------------------------------------------!
         if (ilwr == iupr) then
            dbh     = dbh_lut(ilwr,ipft)
         else
            finterp = log( dbh_lut   (iupr,ipft) / dbh_lut   (ilwr,ipft))                  &
                    * log( balive                / balive_lut(ilwr,ipft))                  &
                    / log( balive_lut(iupr,ipft) / balive_lut(ilwr,ipft))
            finterp = max(lnexp_min,min(lnexp_max,finterp))
            dbh     = dbh_lut(ilwr,ipft) * exp(finterp)
         end if
         ba2h = dbh2h(ipft,dbh)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
    end function ba2h
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !      This function determines an effective DBH for grasses given their leaf biomass.  !
   ! DBH has no real meaning for grasses with the new allometry.                           !
   !---------------------------------------------------------------------------------------!
   real function bl2dbh(bleaf,ipft)
      use pft_coms,     only : dbh_crit    & ! intent(in)
                             , hgt_max     & ! intent(in)
                             , is_grass    & ! intent(in)
                             , C2B         & ! intent(in)
                             , bleaf_adult & ! intent(in)
                             , bleaf_crit  & ! intent(in)
                             , b1Bl_small  & ! intent(in)
                             , b2Bl_small  & ! intent(in)
                             , b1Bl_large  & ! intent(in)
                             , b2Bl_large  & ! intent(in)
                             , dbh_lut     & ! intent(in)
                             , bleaf_lut   & ! intent(in)
                             , le_mask_lut & ! intent(out)
                             , ge_mask_lut & ! intent(out)
                             , is_tropical & ! intent(in)
                             , is_liana    ! ! intent(in)
      use ed_misc_coms, only : igrass      & ! intent(in)
                             , iallom      ! ! intent(in)
      use consts_coms , only : lnexp_min   & ! intent(in)
                             , lnexp_max   ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: bleaf
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      integer             :: ilwr        ! Lower index of the lookup table
      integer             :: iupr        ! Upper index of the lookup table
      real                :: finterp     ! Interpolation factor
      !------------------------------------------------------------------------------------!


      if ( iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft)) ) then
         if (bleaf <= bleaf_lut(1,ipft)) then
            !----- Use the look-up table to find the best dbh. ----------------------------!
            bl2dbh = dbh_lut(1,ipft) * bleaf / bleaf_lut(1,ipft)
            !------------------------------------------------------------------------------!
         else if (bleaf >= bleaf_crit(ipft)) then
            !----- Bleaf should not exceed bleaf_crit.  Set dbh to maximum. ---------------!
            bl2dbh = dbh_crit(ipft)
            !------------------------------------------------------------------------------!
         else
            !----- Use the look-up table to find the best dbh. ----------------------------!
            le_mask_lut(:) = bleaf <= bleaf_lut(:,ipft)
            ge_mask_lut(:) = bleaf >= bleaf_lut(:,ipft)
            ilwr  = maxloc (bleaf_lut(:,ipft),dim=1,mask=le_mask_lut)
            iupr  = minloc (bleaf_lut(:,ipft),dim=1,mask=ge_mask_lut)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      In case ilwr and iupr are the same, we have an exact estimate.  Other-  !
            ! wise, use log-linear interpolation.                                          !
            !------------------------------------------------------------------------------!
            if (ilwr == iupr) then
               bl2dbh = dbh_lut(ilwr,ipft)
            else
               finterp = log( dbh_lut  (iupr,ipft) / dbh_lut  (ilwr,ipft))                 &
                       * log( bleaf                / bleaf_lut(ilwr,ipft))                 &
                       / log( bleaf_lut(iupr,ipft) / bleaf_lut(ilwr,ipft))
               finterp = max(lnexp_min,min(lnexp_max,finterp))
               bl2dbh  = dbh_lut(ilwr,ipft) * exp(finterp)
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!

      else
         !----- Find out whether this is an adult tree or a sapling/grass. ----------------!
         if (bleaf < bleaf_adult(ipft)) then
            mdbh = (bleaf * C2B / b1Bl_small(ipft) ) ** (1./b2Bl_small(ipft))
         else
            mdbh = (bleaf * C2B / b1Bl_large(ipft) ) ** (1./b2Bl_large(ipft))
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     For grasses, limit maximum effective dbh by maximum height.                 !
         !---------------------------------------------------------------------------------!
         if (is_grass(ipft) .and. igrass == 1) then
            bl2dbh = min(mdbh, h2dbh(hgt_max(ipft),ipft))
         else
            bl2dbh = min(mdbh, dbh_crit(ipft))
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function bl2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This function determines the height for grasses given their leaf biomass.        !
   !---------------------------------------------------------------------------------------!
   real function bl2h(bleaf,ipft)
      use pft_coms,      only:  hgt_max    ! ! intent(in), lookup table
      
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in)      :: bleaf
      integer, intent(in)      :: ipft
      !------------------------------------------------------------------------------------!

      !----- Use existing allometric equations to convert leaves to height. ---------------!
      bl2h = min(hgt_max(ipft),dbh2h(ipft,bl2dbh(bleaf,ipft)))

      return
   end function bl2h
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !    Canopy Area allometry from Dietze and Clark (2008).                                !
   !---------------------------------------------------------------------------------------!
   real function size2ca(dbh,hite,sla,ipft)
      use ed_misc_coms , only : iallom         ! ! intent(in)
      use pft_coms     , only : dbh_adult      & ! intent(in)
                              , dbh_crit       & ! intent(in)
                              , hgt_max        & ! intent(in)
                              , is_grass       & ! intent(in)
                              , is_liana       & ! intent(in)
                              , is_tropical    & ! intent(in)
                              , b1Ca           & ! intent(in)
                              , b2Ca           & ! intent(in)
                              , liana_dbh_crit ! ! intent(in)
      use ed_misc_coms , only : igrass         & ! intent(in)
                              , iallom         ! ! intent(in)
      use ed_state_vars, only : patchtype      ! ! structure
      use consts_coms  , only : tiny_num       ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real          , intent(in) :: dbh
      real          , intent(in) :: hite
      real          , intent(in) :: sla
      integer       , intent(in) :: ipft
      !----- Internal variables -----------------------------------------------------------!
      real                       :: loclai ! The maximum local LAI for a given DBH
      real                       :: mdbh   ! The maximum DBH
      !------------------------------------------------------------------------------------!

      if (dbh < tiny_num) then
         loclai  = 0.0
         size2ca = 0.0
      else

         !----- make this function generic to size, not just dbh. -------------------------!
         loclai = sla * size2bl(dbh,hite,ipft) 
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check whether to impose a maximum height.                                   !
         !---------------------------------------------------------------------------------!
         select case (iallom)
         case (0)
            !----- No upper bound in the allometry. ---------------------------------------!
            mdbh = dbh
            !------------------------------------------------------------------------------!

         case default
            !----- Check whether this is grass or tree, impose maximum DBH. ---------------!
            if (is_grass(ipft) .and. igrass == 1) then
               mdbh = min(dbh,h2dbh(hgt_max(ipft),ipft) )
            elseif (is_liana(ipft)) then
               mdbh = min(dbh,liana_dbh_crit)
            else
               mdbh = min(dbh,dbh_crit(ipft))
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !----- Find the nominal crown area. ----------------------------------------------!
         if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
            if (mdbh < dbh_adult(ipft)) then
               size2ca = loclai
            else
               size2ca = b1Ca(ipft) * ( mdbh * mdbh * hite ) ** b2Ca(ipft)
            end if
         else
            size2ca = b1Ca(ipft) * mdbh ** b2Ca(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !----- Local LAI / Crown area should never be less than one. ------------------------!
      size2ca = min (loclai, size2ca)
      !------------------------------------------------------------------------------------!

      return
   end function size2ca
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the commercial volume of a tree (stem/bole volume, in m3).  !
   ! The current equation is a re-fit from Nogueira et al. (2008) so a single set of       !
   ! parameters can be used.  Their equation is for tropical trees only, so temperate and  !
   ! boreal forests may need a different set of parameters.  Grasses are assumed to have   !
   ! no commercial volume.                                                                 !
   !                                                                                       !
   ! Nogueira, E. M., et al.  Estimates of forest biomass in the Brazilian Amazon: new     !
   !    allometric equations and adjustments to biomass from wood-volume inventories.      !
   !    Forest Ecol. Manag., 256(11), 1853-1867, Nov. 2008,                                !
   !    doi:10.1016/j.foreco.2008.07.022.                                                  !
   !---------------------------------------------------------------------------------------!
   real function size2vol(dbh,hgt,ipft)
      use pft_coms    , only : b1Vol    & ! intent(in)
                             , b2Vol    ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: hgt
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      size2vol = b1Vol(ipft) * (dbh * dbh * hgt) ** b2Vol(ipft)

      return
   end function size2vol
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the commercial timber biomass for different PFTs.           !
   !---------------------------------------------------------------------------------------!
   real function size2bt(dbh,hgt,bdead,bsapwooda,bbark,ipft)
      use pft_coms    , only : is_grass    & ! intent(in)
                             , is_tropical & ! intent(in)
                             , rho         & ! intent(in)
                             , agf_bs      & ! intent(in)
                             , C2B         ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: dbh
      real(kind=4), intent(in) :: hgt
      real(kind=4), intent(in) :: bdead
      real(kind=4), intent(in) :: bsapwooda
      real(kind=4), intent(in) :: bbark
      integer     , intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: btimber
      real(kind=4)             :: bagwood
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Calculate commercial timber biomass based on the life form and habitat.       !
      !------------------------------------------------------------------------------------!
      if (is_grass(ipft)) then
         !---------------------------------------------------------------------------------!
         !     Grasses don't have commercial timber.                                       !
         !---------------------------------------------------------------------------------!
         size2bt = 0.0
         !---------------------------------------------------------------------------------!
      else if (is_tropical(ipft)) then
         !---------------------------------------------------------------------------------!
         !     Tropical trees.  We use the tree volume scaled with typical wood density.   !
         !---------------------------------------------------------------------------------!
         btimber  = 1000. * rho(ipft) * size2vol(dbh,hgt,ipft) / C2B
         bagwood  = agf_bs(ipft) * (bdead + bbark) + bsapwooda
         size2bt  = min(btimber,bagwood)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !     Temperate trees.  For the time being we use total aboveground wood biomass, !
         ! so results are consistent with former implementation.                           !
         !---------------------------------------------------------------------------------!
         size2bt  = agf_bs(ipft) * (bdead + bbark) + bsapwooda
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2bt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes bark thickness.  To obtain the actual thickness, we compare !
   ! the actual bark biomass with on-allometry bark biomass.                               !
   !---------------------------------------------------------------------------------------!
   real function size2xb(dbh,hgt,bbark,ipft)
      use pft_coms   , only : qbark       & ! intent(in)
                            , b1Xb        ! ! intent(in)
      use consts_coms, only : tiny_num    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: dbh
      real(kind=4), intent(in) :: hgt
      real(kind=4), intent(in) :: bbark
      integer     , intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: bleaf_max
      real(kind=4)             :: bbark_max
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Calculate commercial timber biomass based on the life form and habitat.       !
      !------------------------------------------------------------------------------------!
      if (qbark(ipft) < tiny_num) then
         size2xb = 0.0
      else
         !----- Find on-allometry bark biomass. -------------------------------------------!
         bleaf_max = size2bl(dbh,hgt,ipft) 
         bbark_max = qbark(ipft) * hgt * bleaf_max
         !---------------------------------------------------------------------------------!

         !----- Scale bark thickness with the ratio between actual and on-allometry. ------!
         size2xb   = b1Xb(ipft) * dbh * bbark / bbark_max
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function size2xb
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   integer function size2krdepth(hite,dbh,ipft,lsl)
      use ed_misc_coms, only : iallom   ! ! intent(in)
      use grid_coms   , only : nzg      ! ! intent(in)
      use soil_coms   , only : slz      ! ! intent(in)
      use pft_coms    , only : b1Rd     & ! intent(in)
                             , b2Rd     ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hite
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      integer, intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      real                :: root_depth
      integer             :: k
      !------------------------------------------------------------------------------------!


      select case (iallom)
      case (0)
         !---------------------------------------------------------------------------------!
         !    Original ED-2.1 (I don't know the source for this equation, though).         !
         ! Grasses get a fixed rooting depth of 70cm.                                      !
         !---------------------------------------------------------------------------------!
         root_depth = b1Rd(ipft)  * (hite * dbh * dbh) ** b2Rd(ipft)
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !    This is just a test allometry, that imposes root depth to be 0.5 m for       !
         ! plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.           !
         !---------------------------------------------------------------------------------!
         root_depth = b1Rd(ipft) * hite ** b2Rd(ipft)
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Root depth is the maximum root depth if the soil is that deep.  Find what is   !
      ! the deepest soil layer this root can go.                                           !
      !------------------------------------------------------------------------------------!
      size2krdepth = nzg
      do k=nzg,lsl+1,-1
         if (root_depth < slz(k)) size2krdepth = k-1
      end do
      !------------------------------------------------------------------------------------!

      return
   end function size2krdepth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function finds the trunk height.  Currently this is based on the following    !
   ! reference, which is for a site in Bolivia:                                            !
   !                                                                                       !
   ! Poorter L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree        !
   !     species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.       !
   !                                                                                       !
   ! Crown length is defined as the height of the tree minus the height of the first leaf  !
   !---------------------------------------------------------------------------------------!
   real function h2crownbh(height,ipft)
      use pft_coms, only : b1Cl & ! intent(in)
                         , b2Cl ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: height
      integer, intent(in) :: ipft
      !----- Local variables. -------------------------------------------------------------!
      real                :: crown_length
      !------------------------------------------------------------------------------------!
      
      crown_length = b1Cl(ipft) * height ** b2Cl(ipft)
      h2crownbh    = max(0.05,height - crown_length)

      return
   end function h2crownbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass (wood + leaves)              !
   !---------------------------------------------------------------------------------------!
   real function ed_biomass(cpatch,ico)
      use pft_coms,      only : agf_bs     ! ! intent(in)
      use ed_state_vars, only : patchtype  ! ! Structure

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target :: cpatch
      integer, intent(in)     :: ico
      !------------------------------------------------------------------------------------!

      ed_biomass = cpatch%bleaf(ico) + cpatch%bsapwooda(ico)                               &
                 + ( cpatch%bdead(ico) + cpatch%bbark(ico) ) * agf_bs(cpatch%pft(ico))

      return
   end function ed_biomass
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine estimates the tree area indices, namely leaf, branch(plus twigs), !
   ! and stem.  For the leaf area index (LAI), we use the specific leaf area (SLA), a      !
   ! constant.  The wood area index WAI is found either from Hormann et al. (2003) or      !
   ! based on the average WAI/LAI ratio from Olivas et al. (2013).                         !
   !                                                                                       !
   ! G. Hormann, S. Irrgan, H. Jochheim, M. Lukes, H. Meesenburg, J. Muller, B. Scheler,   !
   !    J. Scherzer, G. Schuler, B. Schultze, B. Strohbach, F. Suckow, M. Wegehenkel, and  !
   !    G. Wessolek.   Wasserhaushalt von waldokosystemen: methodenleitfaden zur           !
   !    bestimmung der wasserhaushaltskomponenten auf level II-flachen. Technical note,    !
   !    Bundesministerium fur Verbraucherschutz, Ernahrung und Landwirtschaft (BMVEL),     !
   !    Bonn, Germany, 2003. URL http://www.wasklim.de/download/Methodenband.pdf.          !
   !                                                                                       !
   ! P. C. Olivas, S. F. Oberbauer, D. B. Clark, D. A. Clark, M. G. Ryan, J. J. O'Brien,   !
   !    and H. Ordonez.  Comparison of direct and indirect methods for assessing leaf area !
   !    index across a tropical rain forest landscape.  Agric. For. Meteorol., 177,        !
   !    110--116, 2013.                                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine area_indices(cpatch, ico)
      use ed_state_vars, only : patchtype       ! ! Structure
      use pft_coms     , only : dbh_crit        & ! intent(in)
                              , dbh_adult       & ! intent(in)
                              , is_liana        & ! intent(in)
                              , is_tropical     & ! intent(in)
                              , is_grass        & ! intent(in)
                              , b1WAI_small     & ! intent(in)
                              , b2WAI_small     & ! intent(in)
                              , b1WAI_large     & ! intent(in)
                              , b2WAI_large     & ! intent(in)
                              , liana_dbh_crit  ! ! intent(in)
      use rk4_coms     , only : ibranch_thermo  ! ! intent(in)
      use ed_misc_coms , only : igrass          & ! intent(in)
                              , iallom          ! ! intent(in)

      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target :: cpatch
      integer, intent(in)     :: ico
      !----- Local variables --------------------------------------------------------------!
      real                    :: loccai
      real                    :: mdbh
      integer                 :: ipft
      !------------------------------------------------------------------------------------!


      ipft   = cpatch%pft(ico)


      !------------------------------------------------------------------------------------!
      !     First, we compute the LAI.  Because SLA may not be constant, we must use the   !
      ! cohort-level variable as opposed to the pft look-up table.  Likewise, we use       !
      ! actual bleaf because leaves may not be fully flushed.                              !
      !------------------------------------------------------------------------------------!
      cpatch%lai(ico) = cpatch%bleaf(ico) * cpatch%nplant(ico) * cpatch%sla(ico)
      !------------------------------------------------------------------------------------!

      !----- Find the crown area. ---------------------------------------------------------!
      loccai = size2ca(cpatch%dbh(ico),cpatch%hite(ico),cpatch%sla(ico),ipft)
      cpatch%crown_area(ico) = min(1.0, cpatch%nplant(ico) * loccai)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Here we check whether we need to compute the branch, stem, and effective       !
      ! branch area indices.  These are only needed when branch thermodynamics is used,    !
      ! otherwise, simply assign zeroes to them.                                           !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (0) 
         !----- Ignore branches and trunk. ------------------------------------------------!
         cpatch%wai(ico)  = 0.
         !---------------------------------------------------------------------------------!

      case (1,2)
         !---------------------------------------------------------------------------------!
         !     Get the DBH, or potential DBH in case of grasses.                           !
         !---------------------------------------------------------------------------------!
         if (is_grass(ipft) .and. igrass == 1) then 
            !---- Use height for new grasses. ---------------------------------------------!
            mdbh   = min(cpatch%dbh(ico),h2dbh(cpatch%hite(ico),ipft))
         elseif (is_liana(ipft)) then
            mdbh   = min(cpatch%dbh(ico),liana_dbh_crit)
         else
            !---- Use dbh for trees. ------------------------------------------------------!
            mdbh   = min(cpatch%dbh(ico),dbh_crit(ipft))
         end if
         !---------------------------------------------------------------------------------!


         !-----Find WAI. ------------------------------------------------------------------!
         if (iallom == 3 .and. is_tropical(ipft) .and. (.not. is_liana(ipft))) then
            cpatch%wai(ico) = cpatch%nplant(ico) * b1WAI_small(ipft)                       &
                            * ( mdbh * mdbh * cpatch%hite(ico) ) ** b2WAI_small(ipft)
         else if (mdbh < dbh_adult(ipft)) then
            cpatch%wai(ico) = cpatch%nplant(ico)                                           &
                            * b1WAI_small(ipft) * mdbh ** b2WAI_small(ipft)
         else
            cpatch%wai(ico) = cpatch%nplant(ico)                                           &
                            * b1WAI_large(ipft) * mdbh ** b2WAI_large(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine area_indices
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


