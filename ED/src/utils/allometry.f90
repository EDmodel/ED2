!==========================================================================================!
!==========================================================================================!
!    This module is a library with several allometric relationships.                       !
!------------------------------------------------------------------------------------------!
module allometry

   contains
   !=======================================================================================!
   !=======================================================================================!
   real function h2dbh(h,ipft)

      use pft_coms, only:  is_tropical & ! intent(in)
                         , b1Ht        & ! intent(in), lookup table
                         , b2Ht        & ! intent(in), lookup table
                         , hgt_ref     ! ! intent(in), lookup table

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      if (is_tropical(ipft)) then
         h2dbh = exp((log(h)-b1Ht(ipft))/b2Ht(ipft))
      else ! Temperate
         h2dbh = log(1.0-(h-hgt_ref(ipft))/b1Ht(ipft))/b2Ht(ipft)
      end if

      return
   end function h2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2h(ipft, dbh)
      use pft_coms, only:  is_tropical & ! intent(in)
                         , rho         & ! intent(in)
                         , max_dbh     & ! intent(in)
                         , b1Ht        & ! intent(in)
                         , b2Ht        & ! intent(in)
                         , hgt_ref     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer , intent(in) :: ipft
      real    , intent(in) :: dbh
      !------------------------------------------------------------------------------------!

      if (is_tropical(ipft)) then
         dbh2h = exp (b1Ht(ipft) + b2Ht(ipft) * log(min(dbh,max_dbh(ipft))) )
      else !----- Temperate PFT allometry. ------------------------------------------------!
         dbh2h = hgt_ref(ipft) + b1Ht(ipft) * (1.0 - exp(b2Ht(ipft) * dbh))
      end if

      return
   end function dbh2h
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2bd(dbh,ipft)

      use pft_coms    , only : C2B         & ! intent(in)
                             , b1Bs_small  & ! intent(in), lookup table
                             , b2Bs_small  & ! intent(in), lookup table
                             , b1Bs_big    & ! intent(in), lookup table
                             , b2Bs_big    & ! intent(in), lookup table
                             , max_dbh     ! ! intent(in), lookup table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: agb
      real                :: qd
      !------------------------------------------------------------------------------------!

      if (dbh <= max_dbh(ipft)) then
         dbh2bd = b1Bs_small(ipft) / C2B * dbh ** b2Bs_small(ipft)
      else
         dbh2bd = b1Bs_big(ipft) / C2B * dbh ** b2Bs_big(ipft)
      end if

      return
   end function dbh2bd
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the diameter at the breast height given the structural   !
   ! (dead) biomass and the PFT type.  In addition to these inputs, it is also required to !
   ! give a first guess for DBH, which can be the previous value, for example.   This is   !
   ! only used for the new allometry, because we must solve an iterative root-finding      !
   ! method.                                                                               !
   !---------------------------------------------------------------------------------------!
   real function bd2dbh(ipft, bdead)
      use pft_coms    , only : b1Bs_small  & ! intent(in), lookup table
                             , b2Bs_small  & ! intent(in), lookup table
                             , b1Bs_big    & ! intent(in), lookup table
                             , b2Bs_big    & ! intent(in), lookup table
                             , bdead_crit  & ! intent(in), lookup table
                             , C2B         ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in) :: ipft      ! PFT type                            [         ---]
      real   , intent(in) :: bdead     ! Structural (dead) biomass           [   kgC/plant]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Decide which coefficients to use based on the critical bdead.                   !
      !------------------------------------------------------------------------------------!
      if (bdead <= bdead_crit(ipft)) then
         bd2dbh = (bdead / b1Bs_small(ipft) * C2B)**(1.0/b2Bs_small(ipft))
      else
         bd2dbh = (bdead / b1Bs_big(ipft) * C2B)**(1.0/b2Bs_big(ipft))
      end if
      !------------------------------------------------------------------------------------!

      return
   end function bd2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2bl(dbh,ipft)
      use pft_coms    , only : max_dbh     & ! intent(in), lookup table
                             , C2B         & ! intent(in)
                             , b1Bl        & ! intent(in), lookup table
                             , b2Bl        ! ! intent(in), lookup table

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: agb
      real                :: qd
      !------------------------------------------------------------------------------------!


      mdbh   = min(dbh,max_dbh(ipft))
      dbh2bl = b1Bl(ipft) / C2B * mdbh ** b2Bl(ipft)

      return
   end function dbh2bl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Canopy Area allometry from Dietze and Clark (2008).                                !
   !---------------------------------------------------------------------------------------!
   real function dbh2ca(dbh,sla,ipft)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use pft_coms    , only : max_dbh     & ! intent(in)
                             , is_tropical & ! intent(in)
                             , is_grass    & ! intent(in)
                             , b1Ca        & ! intent(in)
                             , b2Ca        ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: sla
      integer, intent(in) :: ipft
      !----- Internal variables -----------------------------------------------------------!
      real                :: loclai ! The maximum local LAI for a given DBH
      real                :: dmbh   ! The maximum 
      !------------------------------------------------------------------------------------!
      if (dbh < tiny(1.0)) then
         loclai = 0.0
         dbh2ca = 0.0
      !----- Based on Poorter et al. (2006) -----------------------------------------------!
      !elseif(is_tropical(ipft) .or. is_grass(ipft)) then
      !   hite   = dbh2h(ipft,dbh)
      !   dbh2ca = 0.156766*hite**1.888
      !----- Based on Dietze and Clark (2008). --------------------------------------------!
      else
         loclai = sla * dbh2bl(dbh,ipft)

         select case (iallom)
         case (0)
            !----- No upper bound in the allometry. ---------------------------------------!
            dbh2ca = b1Ca(ipft) * dbh ** b2Ca(ipft)

         case default
            !----- Impose a maximum crown area. -------------------------------------------!
            dbh2ca = b1Ca(ipft) * min(dbh,max_dbh(ipft)) ** b2Ca(ipft)

         end select
      end if

      !----- Local LAI / Crown area should never be less than one. ------------------------!
      dbh2ca = min (loclai, dbh2ca)

      return
   end function dbh2ca
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the standing volume of a tree.                              !
   !---------------------------------------------------------------------------------------!
   real function dbh2vol(hgt,dbh,ipft)
      use pft_coms    , only : b1Vol    & ! intent(in)
                             , b2Vol    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hgt
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      dbh2vol = b1Vol(ipft) * hgt * dbh ** b2Vol(ipft)

      return
   end function dbh2vol
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   integer function dbh2krdepth(hgt,dbh,ipft,lsl)
      use ed_misc_coms, only : iallom   ! ! intent(in)
      use grid_coms   , only : nzg      ! ! intent(in)
      use soil_coms   , only : slz      ! ! intent(in)
      use pft_coms    , only : b1Rd     & ! intent(in)
                             , b2Rd     ! ! intent(in)
      implicit none 
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: hgt
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      integer, intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      real                :: volume
      real                :: root_depth
      integer             :: k
      !------------------------------------------------------------------------------------!


      !----- Grasses get a fixed rooting depth of 70 cm. ----------------------------------!
      select case (iallom)
      case (0:2)
         !---------------------------------------------------------------------------------!
         !    Original ED-2.1 (I don't know the source for this equation, though).         !
         !---------------------------------------------------------------------------------!
         volume     = dbh2vol(hgt,dbh,ipft) 
         root_depth = b1Rd(ipft)  * volume ** b2Rd(ipft)

      case (3)
         !---------------------------------------------------------------------------------!
         !     This equation is based on Kenzo et al (2009), figure 4e, although a         !
         ! correction was needed to make root depth a function of DBH rather than basal    !
         ! diameter.                                                                       !
         ! Source: Kenzo, T., and co-authors: 2009.  Development of allometric relation-   !
         !            ships for accurate estimation of above- and below-ground biomass in  !
         !            tropical secondary forests in Sarawak, Malaysia. J. Trop. Ecology,   !
         !            25, 371-386.                                                         !
         !---------------------------------------------------------------------------------!
         root_depth = b1Rd(ipft)  * dbh ** b2Rd(ipft)

      case (4)
         !---------------------------------------------------------------------------------!
         !    This is just a test allometry, that imposes root depth to be 0.5 m for       !
         ! plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.           !
         !---------------------------------------------------------------------------------!
         root_depth = b1Rd(ipft) * hgt ** b2Rd(ipft)
      end select


      !------------------------------------------------------------------------------------!
      !     Root depth is the maximum root depth if the soil is that deep.  Find what is   !
      ! the deepest soil layer this root can go.                                           !
      !------------------------------------------------------------------------------------!
      dbh2krdepth = nzg
      do k=nzg,lsl+1,-1
         if (root_depth < slz(k)) dbh2krdepth = k-1
      end do
      !------------------------------------------------------------------------------------!

      return
   end function dbh2krdepth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function finds the trunk height.  Currently this is based on the following    !
   ! reference, which is for a site in Bolivia:                                            !
   !                                                                                       !
   ! Poorter L., L. Bongers, F. Bongers, 2006: Architecture of 54 moist-forest tree        !
   !     species: traits, trade-offs, and functional groups. Ecology, 87, 1289-1301.       !
   !---------------------------------------------------------------------------------------!
   real function h2crownbh(height,ipft)
      use pft_coms, only : b1Cl & ! intent(in)
                         , b2Cl ! ! intent(in)
      implicit none
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
   !     This subroutine finds the total above ground biomass corresponding to stems.      !
   !---------------------------------------------------------------------------------------!
   real function wood_biomass(bdead, bsapwood, pft)
      use pft_coms, only:  agf_bs          ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: bdead
      real    , intent(in) :: bsapwood
      integer , intent(in) :: pft
      !----- Local variables --------------------------------------------------------------!
      real                 :: bstem
      real                 :: absapwood
      !------------------------------------------------------------------------------------!

      bstem        = agf_bs * bdead
      absapwood    = agf_bs * bsapwood
      wood_biomass = bstem + absapwood
      return
   end function wood_biomass
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass (wood + leaves)              !
   !---------------------------------------------------------------------------------------!
   real function ed_biomass(bdead, balive, bleaf, pft, hite, bstorage, bsapwood)
      use pft_coms, only:  agf_bs ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: bdead
      real    , intent(in) :: balive
      real    , intent(in) :: bleaf
      real    , intent(in) :: hite
      real    , intent(in) :: bstorage
      real    , intent(in) :: bsapwood
      integer , intent(in) :: pft
      !----- Local variables --------------------------------------------------------------!
      real                 :: bwood
      !------------------------------------------------------------------------------------!

      bwood      = wood_biomass(bdead, bsapwood, pft)
      ed_biomass = bleaf + bwood

      return
   end function ed_biomass
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine estimates the tree area indices, namely leaf, branch(plus twigs), !
   ! and stem.  For the leaf area index (LAI), we use the specific leaf area (SLA), a      !
   ! constant.  The wood area index WAI is found using the model proposed by Järvelä       !
   ! (2004) to find the specific projected area.                                           !
   !                                                                                       !
   ! Järvelä, J., 2004: Determination of flow resistance caused by non-submerged woody     !
   !                    vegetation. Intl. J. River Basin Management, 2, 61-70.             !
   !                                                                                       !
   !     There is also a very simplified estimation of branch area index, which is just a  !
   ! simple curve adjusted with the information I found in Conijn (1995), which is actual- !
   ! ly for the Sahel...                                                                   !
   !                                                                                       !
   ! Conijn, J.G., 1995: RECAFS: a model for resource competition and cycling in agro-     !
   !                     forestry systems. Rapports Production Soudano-Sahélienne.         !
   !                     Wageningen, 1995.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine area_indices(nplant,bleaf,bdead,balive,dbh,hite,pft,sla,lai,wpa,wai          &
                          ,crown_area,bsapwood)
      use pft_coms    , only : is_tropical     & ! intent(in)
                             , is_grass        & ! intent(in)
                             , rho             & ! intent(in)
                             , C2B             & ! intent(in)
                             , horiz_branch    & ! intent(in)
                             , rbranch         & ! intent(in)
                             , rdiamet         & ! intent(in)
                             , rlength         & ! intent(in)
                             , diammin         & ! intent(in)
                             , ntrunk          & ! intent(in)
                             , conijn_a        & ! intent(in)
                             , conijn_b        & ! intent(in)
                             , conijn_c        & ! intent(in)
                             , conijn_d        ! ! intent(in)
      use consts_coms , only : onethird        & ! intent(in)
                             , pi1             ! ! intent(in)
      use rk4_coms    , only : ibranch_thermo  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer , intent(in)  :: pft        ! Plant functional type            [         ---]
      real    , intent(in)  :: nplant     ! Number of plants                 [    plant/m²]
      real    , intent(in)  :: bleaf      ! Specific leaf biomass            [   kgC/plant]
      real    , intent(in)  :: bdead      ! Specific structural              [   kgC/plant]
      real    , intent(in)  :: balive     ! Specific live tissue biomass     [   kgC/plant]
      real    , intent(in)  :: bsapwood   ! Specific sapwood biomass         [   kgC/plant]
      real    , intent(in)  :: dbh        ! Diameter at breast height        [          cm]
      real    , intent(in)  :: hite       ! Plant height                     [           m]
      real    , intent(in)  :: sla        ! Specific leaf area               [m²leaf/plant]
      real    , intent(out) :: lai        ! Leaf area index                  [   m²leaf/m²]
      real    , intent(out) :: wpa        ! Wood projected area              [   m²wood/m²]
      real    , intent(out) :: wai        ! Wood area index                  [   m²wood/m²]
      real    , intent(out) :: crown_area ! Crown area                       [  m²crown/m²]
      !----- Local variables --------------------------------------------------------------!
      real                  :: bwood      ! Wood biomass                     [   kgC/plant]
      real                  :: swa        ! Specific wood area               [    m²/plant]
      real                  :: bdiamet    ! Diameter of current branch       [           m]
      real                  :: blength    ! Length of each branch            [           m]
      real                  :: nbranch    ! Number of branches               [        ----]
      real                  :: bdmin      ! Minimum diameter                 [           m]
      !----- External functions -----------------------------------------------------------!
      real    , external    :: errorfun ! Error function.
      !------------------------------------------------------------------------------------!
      
      !----- First, we compute the LAI ----------------------------------------------------!
      lai = bleaf * nplant * sla


      !----- Find the crown area. ---------------------------------------------------------!
      crown_area = min(1.0, nplant * dbh2ca(dbh,sla,pft))

      !------------------------------------------------------------------------------------!
      !     Here we check whether we need to compute the branch, stem, and effective       !
      ! branch area indices.  These are only needed when branch thermodynamics is used,    !
      ! otherwise, simply assign zeroes to them.                                           !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      case (0) 
         !----- Ignore branches and trunk. ------------------------------------------------!
         wpa  = 0.
         wai  = 0.
      !------------------------------------------------------------------------------------!


      case (1) 
         !---------------------------------------------------------------------------------!
         !     Use curve fit based on Conijn (1995) model.   Find the total wood biomass   !
         ! and the fraction corresponding to branches.                                     !
         !---------------------------------------------------------------------------------!
         bwood   = wood_biomass(bdead, bsapwood, pft)
         if (is_grass(pft)) then
            swa = conijn_a(pft)
         else
            swa = conijn_a(pft)                                                            &
                + conijn_b(pft) * errorfun(conijn_c(pft)*C2B*bwood + conijn_d(pft))
         end if
         wai = nplant * bwood * swa
         wpa = wai * dbh2ca(dbh,sla,pft)
         !---------------------------------------------------------------------------------!


      case (2) 
         !---------------------------------------------------------------------------------!
         !     Use  Järvelä (2004) method.                                                 !
         !---------------------------------------------------------------------------------!
         !----- Now we check the first branching height, that will be the trunk height. ---!
         blength = h2crownbh(hite,pft)
         !----- Main branch diameter is DBH (in meters) -----------------------------------!
         bdiamet = dbh * 0.01
         !----- Minimum branch diameter (in meters) ---------------------------------------!
         bdmin   = diammin(pft) * 0.01
         !----- Number of main "branches" (trunk), this is usually 1. ---------------------!
         nbranch = ntrunk(pft)

         swa = nbranch * blength * bdiamet
         !---------------------------------------------------------------------------------!
         !     Initialize branch values with trunk.                                        !
         !---------------------------------------------------------------------------------!
         branchloop: do
            if (bdiamet < bdmin) exit branchloop
            !----- Updating branch habits. ------------------------------------------------!
            bdiamet = bdiamet / rdiamet(pft)
            blength = blength / rlength(pft)
            nbranch = nbranch * rbranch(pft)
            swa     = swa + nbranch * blength * bdiamet
         end do branchloop
         !----- The wood projected area and the wood area index. --------------------------!
         wpa = nplant       * swa
         wai = horiz_branch(pft) * wpa
      case (3) 
         !---------------------------------------------------------------------------------!
         !    Use the equation by:                                                         !
         !    Ahrends, B., C. Penne, O. Panferov, 2010: Impact of target diameter          !
         !        harvesting on spatial and temporal pattern of drought risk in forest     !
         !        ecosystems under climate change conditions.  The Open Geography Journal, !
         !        3, 91-102  (they didn't develop the allometry, but the original          !
         !        reference is in German...)                                               !
         !---------------------------------------------------------------------------------!
         select case (pft)
         case (6:8,17)
             !----- Conifers. -------------------------------------------------------------!
             wai = nplant * 0.0553 * 0.5 * dbh ** 1.9769
         case default
             !----- Broadleaf trees and grasses. ------------------------------------------!
             wai = nplant * 0.0192 * 0.5 * dbh ** 2.0947
         end select
         !---------------------------------------------------------------------------------!
         wpa = wai * dbh2ca(dbh,sla,pft)

      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!

      return
   end subroutine area_indices
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


