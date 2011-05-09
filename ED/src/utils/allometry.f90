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

         case (1)
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
   real function calc_root_depth(h,dbh,ipft)
      use consts_coms, only: pi1
      implicit none 
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: volume
      !------------------------------------------------------------------------------------!

      select case (ipft)
      case(1,5) !----- Grasses get a fixed rooting depth of 70 cm. ------------------------!
         calc_root_depth = -0.70
      case default
         volume           = h * 0.65 * pi1 * (dbh*0.11)*(dbh*0.11)
         calc_root_depth = -10.0**(0.545 + 0.277*log10(volume))
      end select

      return
   end function calc_root_depth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   integer function assign_root_depth(rd, lsl)
      use grid_coms, only: nzg
      use soil_coms, only: slz
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: rd
      integer , intent(in) :: lsl
      !----- Local variables --------------------------------------------------------------!
      integer :: k
      !------------------------------------------------------------------------------------!
      
      assign_root_depth = nzg
      do k=nzg,lsl+1,-1
         if (rd < slz(k)) assign_root_depth = k-1
      end do
      return
   end function assign_root_depth
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This function finds the trunk height, based on Antonarakis (2008) estimation for   !
   ! secondary forests.                                                                    !
   !---------------------------------------------------------------------------------------!
   real function h2trunkh(h,ipft)
      use pft_coms, only : b1Tht & ! intent(in)
                         , b2Tht ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      
      h2trunkh = b1Tht(ipft) * h ** b2THt(ipft)

      return
   end function h2trunkh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass corresponding to wood (stem  !
   ! + branches).                                                                          !
   !---------------------------------------------------------------------------------------!
   real function wood_biomass(bdead, balive, pft, hite, bsapwood)
      use pft_coms, only:  agf_bs          ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: bdead
      real    , intent(in) :: balive
      real    , intent(in) :: bsapwood
      real    , intent(in) :: hite
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

      bwood      = wood_biomass(bdead, balive, pft, hite, bsapwood)
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
      !----- Ignore branches and trunk. ---------------------------------------------------!
      case (0) 
         wpa  = 0.
         wai  = 0.
      !------------------------------------------------------------------------------------!


      !----- Use curve fit based on Conijn (1995) model. ----------------------------------!
      case (1) 
         !---------------------------------------------------------------------------------!
         !     Find the total wood biomass and the fraction corresponding to branches.     !
         !---------------------------------------------------------------------------------!
         bwood   = wood_biomass(bdead, balive, pft, hite, bsapwood)
         if (is_grass(pft)) then
            swa = conijn_a(pft)
         else
            swa = conijn_a(pft)                                                            &
                + conijn_b(pft) * errorfun(conijn_c(pft)*C2B*bwood + conijn_d(pft))
         end if
         wai = nplant * bwood * swa
         wpa = wai * dbh2ca(dbh,sla,pft)
         !---------------------------------------------------------------------------------!


      !----- Use  Järvelä (2004) method. --------------------------------------------------!
      case (2) 
         !----- Now we check the first branching height, that will be the trunk height. ---!
         blength = h2trunkh(hite,pft)
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
      !------------------------------------------------------------------------------------!
      end select

      return
   end subroutine area_indices
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


