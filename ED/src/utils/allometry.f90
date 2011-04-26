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
   real function dbh2bd(dbh,height,ipft)

      use ed_misc_coms, only : iallom      ! ! intent(in)
      use pft_coms    , only : is_tropical & ! intent(in), lookup table
                             , C2B         & ! intent(in)
                             , agf_bs      & ! intent(in)
                             , b1Bs_small  & ! intent(in), lookup table
                             , b2Bs_small  & ! intent(in), lookup table
                             , b1Bs_big    & ! intent(in), lookup table
                             , b2Bs_big    & ! intent(in), lookup table
                             , max_dbh     & ! intent(in), lookup table
                             , q           & ! intent(in), lookup table
                             , qsw         ! ! intent(in), lookup table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: height
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: agb
      real                :: qd
      !------------------------------------------------------------------------------------!

      if (is_tropical(ipft) .and. iallom == 1) then
         agb    = dbh2agb(dbh,ipft)
         qd     = dbh2qd(dbh,ipft)
         dbh2bd = qd * agb / (agf_bs * (1. + q(ipft) + qsw(ipft) * height + qd))
      elseif (dbh <= max_dbh(ipft)) then
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
   real function dbh2bl(dbh,height,ipft)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use pft_coms    , only : is_tropical & ! intent(in), lookup table
                             , max_dbh     & ! intent(in), lookup table
                             , C2B         & ! intent(in)
                             , agf_bs      & ! intent(in)
                             , b1Bl        & ! intent(in), lookup table
                             , b2Bl        & ! intent(in), lookup table
                             , q           & ! intent(in), lookup table
                             , qsw         ! ! intent(in), lookup table
   
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: height
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: agb
      real                :: qd
      !------------------------------------------------------------------------------------!


      if (is_tropical(ipft) .and. iallom == 1) then
         agb    = dbh2agb(dbh,ipft)
         qd     = dbh2qd(dbh,ipft)
         dbh2bl = agb / (agf_bs * (1. + q(ipft) + qsw(ipft) * height + qd))
      else
         mdbh   = min(dbh,max_dbh(ipft))
         dbh2bl = b1Bl(ipft) / C2B * mdbh ** b2Bl(ipft)
      end if

      return
   end function dbh2bl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Canopy Area allometry from Dietze and Clark (2008).                                !
   !---------------------------------------------------------------------------------------!
   real function dbh2ca(dbh,height,sla,ipft)
      use pft_coms, only : is_tropical & ! intent(in)
                         , is_grass    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: height
      real   , intent(in) :: sla
      integer, intent(in) :: ipft
      !----- Internal variables -----------------------------------------------------------!
      real                :: loclai ! The maximum local LAI for a given DBH
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
         loclai = sla * dbh2bl(dbh,height,ipft)
         dbh2ca = 2.490154*dbh**0.8068806
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
   real function h2trunkh(h)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real , intent(in) :: h
      !------------------------------------------------------------------------------------!
      
      h2trunkh = 0.4359 * h ** 0.878

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
   !     This function finds the ratio between dead and leaf biomass for a given DBH.      !
   !---------------------------------------------------------------------------------------!
   real function dbh2qd(dbh,ipft)
      use pft_coms, only : max_dbh    & ! intent(in), look-up table
                         , b1Bs_small & ! intent(in), look-up table
                         , b1Bs_big   & ! intent(in), look-up table
                         , b2Bs_small & ! intent(in), look-up table
                         , b2Bs_big   & ! intent(in), look-up table
                         , b1Bl       & ! intent(in), look-up table
                         , b2Bl       ! ! intent(in), look-up table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      if (dbh <= max_dbh(ipft)) then
         dbh2qd = (b1Bs_small(ipft)/b1Bl(ipft)) * dbh ** (b2Bs_small(ipft) - b2Bl(ipft))
      else
         dbh2qd = (b1Bs_big(ipft)/b1Bl(ipft)) * dbh ** b2Bs_big(ipft)                      &
                / max_dbh(ipft) ** b2Bl (ipft)
      end if

      return
   end function dbh2qd
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the above-ground biomass based on the DBH.  This function  !
   ! comes from equation 2 of:                                                             !
   !                                                                                       !
   ! Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial      !
   !    patterns in Amazonian forest biomass.  Glob. Change Biol., 10, 545-562.            !
   !---------------------------------------------------------------------------------------!
   real function dbh2agb(dbh,ipft)

      use pft_coms, only : is_tropical & ! intent(in), lookup table
                         , rho         & ! intent(in), lookup table
                         , C2B         & ! intent(in)
                         , b1agb       & ! intent(in), lookup table
                         , b2agb       & ! intent(in), lookup table
                         , b3agb       ! ! intent(in), lookup table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!

      if (is_tropical(ipft)) then
         dbh2agb = rho(ipft) * exp(b1agb(ipft) + b2agb(ipft) * log(dbh))                   &
                 / (b3agb(ipft) * C2B)
      else
         write(unit=*,fmt='(a,1x,es12.5)') ' DBH      =',dbh
         write(unit=*,fmt='(a,1x,i12)')    ' PFT      =',ipft
         write(unit=*,fmt='(a,11x,l1)')    ' TROPICAL =',is_tropical(ipft)
         
         call fatal_error('This function cannot be called for temperate PFTs!'             &
                         ,'dbh2agb','allometry.f90')
      end if

      return
   end function dbh2agb
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the DBH based on the above-ground biomass.  This function  !
   ! comes from equation 2 of:                                                             !
   !                                                                                       !
   ! Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial      !
   !    patterns in Amazonian forest biomass.  Glob. Change Biol., 10, 545-562.            !
   !---------------------------------------------------------------------------------------!
   real function agb2dbh(agb,ipft)

      use pft_coms, only : is_tropical & ! intent(in), lookup table
                         , rho         & ! intent(in), lookup table
                         , C2B         & ! intent(in)
                         , b1agb       & ! intent(in), lookup table
                         , b2agb       & ! intent(in), lookup table
                         , b3agb       ! ! intent(in), lookup table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: agb
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      if (is_tropical(ipft)) then
         agb2dbh = exp( ( -b1agb(ipft) + log(agb * b3agb(ipft) * C2B / rho(ipft)))         &
                      / b2agb(ipft)) 
      else
         write(unit=*,fmt='(a,1x,es12.5)') ' AGB      =',agb
         write(unit=*,fmt='(a,1x,i12)')    ' PFT      =',ipft
         write(unit=*,fmt='(a,11x,l1)')    ' TROPICAL =',is_tropical(ipft)
         
         call fatal_error('This function cannot be called for temperate PFTs!'             &
                         ,'agb2dbh','allometry.f90')
      end if

      return
   end function agb2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the derivative of the ratio between dead and leaf biomass for !
   ! a given DBH.                                                                          !
   !---------------------------------------------------------------------------------------!
   real function dqd_ddbh(dbh,qd,ipft)
      use pft_coms, only : max_dbh    & ! intent(in), look-up table
                         , b1Bs_small & ! intent(in), look-up table
                         , b1Bs_big   & ! intent(in), look-up table
                         , b2Bs_small & ! intent(in), look-up table
                         , b2Bs_big   & ! intent(in), look-up table
                         , b1Bl       & ! intent(in), look-up table
                         , b2Bl       ! ! intent(in), look-up table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: qd
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      if (dbh <= max_dbh(ipft)) then
         dqd_ddbh = (b2Bs_small(ipft) - b2Bl(ipft)) * qd / dbh
      else
         dqd_ddbh = b2Bs_big(ipft) * qd / dbh
      end if

      return
   end function dqd_ddbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the derivative of the height for a given DBH.                 !
   !---------------------------------------------------------------------------------------!
   real function dh_ddbh(dbh,height,ipft)
      use pft_coms, only : max_dbh    & ! intent(in), look-up table
                         , b2Ht       ! ! intent(in), look-up table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: height
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      if (dbh <= max_dbh(ipft)) then
         dh_ddbh = b2Ht(ipft) * height / dbh
      else
         dh_ddbh = 0.
      end if

      return
   end function dh_ddbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the derivative of the above-ground biomass for a given DBH.   !
   ! This function is derived from equation 2 of:                                          !
   !                                                                                       !
   ! Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial      !
   !    patterns in Amazonian forest biomass.  Glob. Change Biol., 10, 545-562.            !
   !---------------------------------------------------------------------------------------!
   real function dagb_ddbh(dbh,agb,ipft)
      use pft_coms, only : is_tropical & ! intent(in), look-up table
                         , max_dbh     & ! intent(in), look-up table
                         , b2agb       ! ! intent(in), look-up table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: agb
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!


      if (is_tropical(ipft)) then
         dagb_ddbh = b2agb(ipft) * agb / dbh
      else
         write(unit=*,fmt='(a,1x,es12.5)') ' DBH      =',dbh
         write(unit=*,fmt='(a,1x,es12.5)') ' AGB      =',agb
         write(unit=*,fmt='(a,1x,i12)')    ' PFT      =',ipft
         write(unit=*,fmt='(a,11x,l1)')    ' TROPICAL =',is_tropical(ipft)
         
         call fatal_error('This function cannot be called for temperate PFTs!'             &
                         ,'agb2dbh','allometry.f90')
      end if

      return
   end function dagb_ddbh
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
      crown_area = min(1.0, nplant * dbh2ca(dbh,hite,sla,pft))

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
         !     Finding the total wood biomass and the fraction corresponding to branches.  !
         !---------------------------------------------------------------------------------!
         bwood   = wood_biomass(bdead, balive, pft, hite, bsapwood)
         if (is_grass(pft)) then
            swa = conijn_a(pft)
         else
            swa = conijn_a(pft)                                                            &
                + conijn_b(pft) * errorfun(conijn_c(pft)*C2B*bwood + conijn_d(pft))
         end if
         wai = nplant * bwood * swa
         wpa = wai * dbh2ca(dbh,hite,sla,pft)
         !---------------------------------------------------------------------------------!


      !----- Use  Järvelä (2004) method. --------------------------------------------------!
      case (2) 
         !----- Now we check the first branching height, that will be the trunk height. ---!
         blength = h2trunkh(hite)
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






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine computes the diameter at the breast height given the structural   !
   ! (dead) biomass and the PFT type.  In addition to these inputs, it is also required to !
   ! give a first guess for DBH, which can be the previous value, for example.   This is   !
   ! only used for the new allometry, because we must solve an iterative root-finding      !
   ! method.                                                                               !
   !---------------------------------------------------------------------------------------!
   real function bd2dbh(ipft, bdead,dbh1st)
      use ed_misc_coms, only : iallom      ! ! intent(in)
      use pft_coms    , only : is_tropical & ! intent(in), lookup table
                             , b1Bs_small  & ! intent(in), lookup table
                             , b2Bs_small  & ! intent(in), lookup table
                             , b1Bs_big    & ! intent(in), lookup table
                             , b2Bs_big    & ! intent(in), lookup table
                             , bdead_crit  & ! intent(in), lookup table
                             , C2B         & ! intent(in)
                             , agf_bsi     & ! intent(in)
                             , hgt_max     & ! intent(in), lookup table
                             , min_dbh     & ! intent(in), lookup table
                             , max_dbh     & ! intent(in), lookup table
                             , q           & ! intent(in), lookup table
                             , qsw         ! ! intent(in), lookup table
      use therm_lib   , only : maxfpo      & ! intent(in)
                             , toler       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in) :: ipft      ! PFT type                            [         ---]
      real   , intent(in) :: bdead     ! Structural (dead) biomass           [   kgC/plant]
      real   , intent(in) :: dbh1st    ! 1st guess for DBH                   [          cm]
      !----- Local variables. -------------------------------------------------------------!
      real                :: deriv     ! Function derivative                 [kgC/plant/cm]
      real                :: fun       ! Function for which we seek a root.  [   kgC/plant]
      real                :: qd        ! Bdead:Bleaf ratio                   [         ---]
      real                :: agb       ! On-allometry above-ground biomass   [   kgC/plant]
      real                :: h         ! Height                              [         ---]
      real                :: qd_prime  ! Bdead:Bleaf ratio derivative        [        1/cm]
      real                :: agb_prime ! Above-ground biomass derivative     [kgC/plant/cm]
      real                :: h_prime   ! Height derivative                   [         ---]
      real                :: salloci   ! Inverse of sum of allocation ratios [         ---]
      real                :: funa      ! Smallest  guess function            [   kgC/plant]
      real                :: funz      ! Largest   guess function            [   kgC/plant]
      real                :: dbha      ! Smallest guess (or previous guess)  [          cm]
      real                :: dbhz      ! Largest   guess (or new guess)      [          cm]
      real                :: delta     ! Aux. var (2nd guess for bisection)  [         ---]
      integer             :: itn       ! Iteration counter                   [         ---]
      integer             :: itb       ! Iteration counter                   [         ---]
      logical             :: converged ! Method converged                    [         T|F]
      logical             :: zside     ! Flag to check for 1-sided approach  [         T|F]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    If we don't need iterative methods (old allometry or temperate PFTs), we solve  !
      ! here and quit the function.                                                        !
      !------------------------------------------------------------------------------------!
      if ( (.not. is_tropical(ipft)) .or. iallom /= 1) then
         if (bdead <= bdead_crit(ipft)) then
            bd2dbh = (bdead / b1Bs_small(ipft) * C2B)**(1.0/b2Bs_small(ipft))
         else
            bd2dbh = (bdead / b1Bs_big(ipft) * C2B)**(1.0/b2Bs_big(ipft))
         end if
         return
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     New tropical allometry.  We start with the first guess and evaluate the        !
      ! function and the derivative.                                                       !
      !------------------------------------------------------------------------------------!
      dbhz      = dbh1st
      !----- Auxiliary variables, and function. -------------------------------------------!
      h         = dbh2h(ipft,dbhz)
      qd        = dbh2qd(dbhz,ipft)
      agb       = dbh2agb(dbhz,ipft)
      salloci   = 1. / (1. + q(ipft) + qsw(ipft) * h + qd)
      funz      = qd * agb * salloci * agf_bsi / bdead - 1.
      !----- Auxiliary variables, and function derivative. --------------------------------!
      qd_prime  = dqd_ddbh(dbhz,qd,ipft)
      agb_prime = dagb_ddbh(dbhz,agb,ipft)
      h_prime   = dh_ddbh(dbhz,h,ipft)
      deriv     = (funz + 1.) * ( qd_prime / qd + agb_prime / agb                          &
                                - (qsw(ipft) * h_prime + qd_prime) * salloci )
      !----- Copy dbha to dbhz just in case it fails right at the 1st Newton's iteration. -!
      dbha      = dbhz
      funa      = funz
      fun       = funz

      !------------------------------------------------------------------------------------!
      !    Enter Newton's method loop.                                                     !
      !------------------------------------------------------------------------------------!
      converged = .false.
      newloop: do itn = 1,maxfpo/6
         !---------------------------------------------------------------------------------!
         !    Check whether the derivative is not too flat, or whether the next step could !
         ! make the next DBH guess negative.                                               !
         !---------------------------------------------------------------------------------!
         if (abs(deriv) < toler) then 
            exit newloop
         elseif (dbhz - fun / deriv < toler) then
            exit newloop
         end if
         !---------------------------------------------------------------------------------!


         !----- Copy current guess to previous. -------------------------------------------!
         dbha = dbhz
         funa = fun
         !---------------------------------------------------------------------------------!

         !----- New guess, its function, and derivative evaluation. -----------------------!
         dbhz = dbha - fun / deriv
         !----- Auxiliary variables, and function. ----------------------------------------!
         h         = dbh2h(ipft,dbhz)
         qd        = dbh2qd(dbhz,ipft)
         agb       = dbh2agb(dbhz,ipft)
         salloci   = 1. / (1. + q(ipft) + qsw(ipft) * h + qd)
         fun       = qd * agb * salloci * agf_bsi / bdead - 1.
         !----- Auxiliary variables, and function derivative. -----------------------------!
         qd_prime  = dqd_ddbh(dbhz,qd,ipft)
         agb_prime = dagb_ddbh(dbhz,agb,ipft)
         h_prime   = dh_ddbh(dbhz,h,ipft)
         deriv     = (funz + 1.) * ( qd_prime / qd + agb_prime / agb                       &
                                   - (qsw(ipft) * h_prime + qd_prime) * salloci )
         !---------------------------------------------------------------------------------!



         !----- Check whether it converged to a solution. ---------------------------------!
         converged = 2.0 * abs(dbhz - dbha) < toler * (abs(dbhz)+abs(dbha))
         if (converged) then
            bd2dbh = 0.5 * (dbha + dbhz)
            return
         elseif (fun == 0.) then !Converged by luck!
            bd2dbh = dbhz
            return
         end if
         !---------------------------------------------------------------------------------!
      end do newloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If we've reached this point then Newton's method failed.  We will try a bisec- !
      ! tion approach instead, but first we must find two guesses, with opposite signs. We !
      ! may have it already, if not, we try guesses increasingly farther from the last     !
      ! guess.                                                                             !
      !------------------------------------------------------------------------------------!
      if (funa * fun < 0.) then
         funz  = fun
         zside = .true.
      else
         if (abs(fun-funa) < 100.*toler*abs(dbha)) then
            delta = 100.*toler*abs(dbha)
         else
            delta = max(abs(funa * (dbhz-dbha)/(fun-funa)),100.*toler*abs(dbha))
         end if
         dbhz  = dbha + delta
         zside = .false.
         zgssloop: do itb=1,maxfpo
            dbhz    = max(toler, dbha + real((-1)**itb * (itb+3)/2) * delta)
            qd      = dbh2qd(dbhz,ipft)
            agb     = dbh2agb(dbhz,ipft)
            salloci = 1. / (1. + q(ipft) + qsw(ipft) * h + qd)
            funz    = qd * agb * salloci * agf_bsi / bdead - 1.
            zside   = funa*funz < 0.
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)'          ) '------------------------------------------'
            write (unit=*,fmt='(a)'          ) ' Failed finding second guess...'
            write (unit=*,fmt='(a)'          ) ' Input data'
            write (unit=*,fmt='(a,1x,i5)'    ) ' - PFT     =',ipft
            write (unit=*,fmt='(a,1x,es12.5)') ' - BDEAD   =',bdead
            write (unit=*,fmt='(a,1x,es12.5)') ' - DBH1ST  =',dbh1st
            write (unit=*,fmt='(a)'          ) ' '
            write (unit=*,fmt='(a)'          ) ' Previous iteration'
            write (unit=*,fmt='(a,1x,i5)'    ) ' - ITN     =',itn
            write (unit=*,fmt='(a,1x,es12.5)') ' - DBHA    =',dbha
            write (unit=*,fmt='(a,1x,es12.5)') ' - DBHZ    =',dbhz
            write (unit=*,fmt='(a,1x,es12.5)') ' - FUNA    =',funa
            write (unit=*,fmt='(a,1x,es12.5)') ' - FUNZ    =',funz
            write (unit=*,fmt='(a,1x,es12.5)') ' - DELTA   =',delta
            write (unit=*,fmt='(a)'          ) '------------------------------------------'
            call fatal_error('Failed finding the second guess for regula falsi'            &
                            ,'bd2dbh','allometry.f90')
         end if
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    The mofied Regula Falsi (Illinois) method.  This is similar to the normal       !
      ! bisection method, except that it tries to avoid one-sided (i.e. slow) approach.    !
      !------------------------------------------------------------------------------------!
      bisloop: do itb=itn,maxfpo
         bd2dbh =  (funz*dbha - funa*dbhz)/(funz-funa)

         !---------------------------------------------------------------------------------!
         !     Now that we updated the guess, check whether they are really close. If so,  !
         ! it converged, I can use this as my guess.                                       !
         !---------------------------------------------------------------------------------!
         converged = 2.0 * abs(bd2dbh-dbha) < toler * (abs(dbha)+abs(dbhz))
         if (converged) exit bisloop

         !------ Find the new function evaluation. ----------------------------------------!
         h         = dbh2h(ipft,bd2dbh)
         qd        = dbh2qd(bd2dbh,ipft)
         agb       = dbh2agb(bd2dbh,ipft)
         salloci   = 1. / (1. + q(ipft) + qsw(ipft) * h + qd)
         fun       = qd * agb * salloci * agf_bsi / bdead - 1.

         !------ Defining my new interval based on the intermediate value theorem. --------!
         if (fun*funa < 0. ) then
            dbhz  = bd2dbh
            funz  = fun
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa = funa * 0.5
            !----- We've just updated zside, set zside to true. ---------------------------!
            zside = .true.
         else
            dbha   = bd2dbh
            funa   = fun
            !----- If we are updating aside again, modify aside (Illinois method) ---------!
            if (.not. zside) funz = funz * 0.5
            !----- We've just updated aside, set zside to true. ---------------------------!
            zside = .false.
         end if
      end do bisloop

      if (.not.converged) then
         write (unit=*,fmt='(a)'          ) '------------------------------------------'
         write (unit=*,fmt='(a)'          ) ' Failed finding DBH...'
         write (unit=*,fmt='(a)'          ) '------------------------------------------'
         write (unit=*,fmt='(a)'          ) ' Input data'
         write (unit=*,fmt='(a,1x,i5)'    ) ' - PFT     =',ipft
         write (unit=*,fmt='(a,1x,es12.5)') ' - BDEAD   =',bdead
         write (unit=*,fmt='(a,1x,es12.5)') ' - DBH1ST  =',dbh1st
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a)'          ) ' Previous iteration'
         write (unit=*,fmt='(a,1x,es12.5)') ' - DBHA    =',dbha
         write (unit=*,fmt='(a,1x,es12.5)') ' - DBHZ    =',dbhz
         write (unit=*,fmt='(a,1x,es12.5)') ' - FUNA    =',funa
         write (unit=*,fmt='(a,1x,es12.5)') ' - FUNZ    =',funz
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(a)'          ) ' Current iteration'
         write (unit=*,fmt='(a,1x,es12.5)') ' - HEIGHT  =',h
         write (unit=*,fmt='(a,1x,es12.5)') ' - QD      =',qd
         write (unit=*,fmt='(a,1x,es12.5)') ' - AGB     =',agb
         write (unit=*,fmt='(a,1x,es12.5)') ' - SALLOCI =',agb
         write (unit=*,fmt='(a,1x,es12.5)') ' - FUN     =',fun
         write (unit=*,fmt='(a)'          ) '------------------------------------------'
         call fatal_error('DBH didn''t converge, I gave up, sorry!!!'                      &
                         ,'bd2dbh','allometry.f90')
      end if

      return
   end function bd2dbh
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


