!==========================================================================================!
!==========================================================================================!
!    This module is a library with several allometric relationships.                       !
!------------------------------------------------------------------------------------------!
module allometry

   !----- Constants shared by both bdead and bleaf ----------------------------------------!
   real, parameter :: a1    =  -1.981
   real, parameter :: b1    =   1.047
   real, parameter :: dcrit = 100.0
   real, parameter :: ff    =   0.640
   real, parameter :: gg    =   0.370
   !----- Constants used by bdead only ----------------------------------------------------!
   real, parameter :: c1d   =   0.572
   real, parameter :: d1d   =   0.931
   real, parameter :: a2d   =  -1.086
   real, parameter :: b2d   =   0.876
   real, parameter :: c2d   =   0.604
   real, parameter :: d2d   =   0.871
   !----- Constants used by bleaf only ----------------------------------------------------!
   real, parameter :: c1l   =  -0.584
   real, parameter :: d1l   =   0.550
   real, parameter :: a2l   =  -4.111
   real, parameter :: b2l   =   0.605
   real, parameter :: c2l   =   0.848
   real, parameter :: d2l   =   0.438
   !---------------------------------------------------------------------------------------!

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
         h2dbh = 10.0**((log10(h)-0.37)/0.64)
      else ! Temperate
         h2dbh = log(1.0-(h-hgt_ref(ipft))/b1Ht(ipft))/b2Ht(ipft)
      end if

      return
   end function h2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2bd(dbh,h,ipft)

      use pft_coms, only:  is_tropical & ! intent(in), lookup table
                         , rho         & ! intent(in), lookup table
                         , C2B         & ! intent(in)
                         , b1Bs        & ! intent(in), lookup table
                         , b2Bs        & ! intent(in), lookup table
                         , max_dbh     ! ! intent(in), lookup table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: p,r,qq
      !------------------------------------------------------------------------------------!

      if (is_tropical(ipft)) then

         if (dbh > max_dbh(ipft)) then
            p  = a1 + c1d * log(h) + d1d * log(rho(ipft))
            r  = ( (a2d - a1) + (c2d - c1d)*log(h) + log(rho(ipft))                        &
                 * (d2d - d1d)) * (1.0/log(dcrit))
            qq = 2.0 * b2d + r
         else
            p  = a1 + c1d * gg * log(10.0) + d1d * log(rho(ipft))
            r  = ( (a2d - a1) + gg * log(10.0) * (c2d - c1d) + log(rho(ipft))              &
                 * (d2d - d1d)) * (1.0/log(dcrit))
            qq = 2.0 * b2d + c2d * ff + r
         end if
         dbh2bd = exp(p) / C2B * dbh**qq
      else !----- Temperate PFT -----------------------------------------------------------!
         dbh2bd = b1Bs(ipft) / C2B * dbh**b2Bs(ipft)
      end if

      return
   end function dbh2bd
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function dbh2bl(dbh,ipft)
      use pft_coms, only:  is_tropical & ! intent(in), lookup table
                         , rho         & ! intent(in), lookup table
                         , max_dbh     & ! intent(in), lookup table
                         , C2B         & ! intent(in)
                         , b1Bl        & ! intent(in), lookup table
                         , b2Bl        ! ! intent(in), lookup table
   
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: p,r,qq
      !------------------------------------------------------------------------------------!

      if (is_tropical(ipft)) then
         p  = a1 + c1l * gg * log(10.0) + d1l * log(rho(ipft))
         r  = ( (a2l - a1) + gg * log(10.0) * (c2l - c1l) + log(rho(ipft))                 &
              * (d2l - d1l)) * (1.0/log(dcrit))
         qq = 2.0 * b2l + c2l * ff + r  

         if(dbh <= max_dbh(ipft))then
            dbh2bl = exp(p) / C2B * dbh**qq
         else 
            dbh2bl = exp(p) / C2B * max_dbh(ipft)**qq
         end if
      else !----- Temperate ---------------------------------------------------------------!
         mdbh   = min(dbh,max_dbh(ipft))
         dbh2bl = b1Bl(ipft) /C2B * mdbh**b2Bl(ipft)
      end if

      return
   end function dbh2bl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Canopy Area allometry from Dietze and Clark (2008).                                !
   !---------------------------------------------------------------------------------------!
   real function dbh2ca(dbh,ipft)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      if(dbh < tiny(1.0)) then
         dbh2ca = 0.0
      else
         dbh2ca = 2.490154*dbh**0.8068806
      end if
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
         calc_root_depth = -0.7
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
         if (dbh <= max_dbh(ipft)) then
            !----- This means that height is below its maximum. ---------------------------!
            dbh2h = 10.0 ** (log10(dbh) * 0.64 + 0.37)
         else
            !----- Height is at maximum. --------------------------------------------------!
            dbh2h = 10.0 ** (log10(max_dbh(ipft)) * 0.64 + 0.37)
         end if
      else !----- Temperate PFT allometry. ------------------------------------------------!
         dbh2h = hgt_ref(ipft) + b1Ht(ipft) * (1.0 - exp(b2Ht(ipft) * dbh))
      end if

      return
   end function dbh2h
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass corresponding to wood (stem  !
   ! + branches).                                                                          !
   !---------------------------------------------------------------------------------------!
   real function wood_biomass(bdead, balive, pft, hite)
      use pft_coms, only:  agf_bs          & ! intent(in)
                         , q               & ! intent(in)
                         , qsw             ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: bdead
      real    , intent(in) :: balive
      real    , intent(in) :: hite
      integer , intent(in) :: pft
      !----- Local variables --------------------------------------------------------------!
      real                 :: bstem
      real                 :: bsapwood
      !------------------------------------------------------------------------------------!

      bstem        = agf_bs * bdead
      bsapwood     = agf_bs * balive * qsw(pft) * hite / (1.0 + q(pft) + qsw(pft) * hite)
      wood_biomass = bstem + bsapwood
      return
   end function wood_biomass
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine finds the total above ground biomass (wood + leaves)              !
   !---------------------------------------------------------------------------------------!
   real function ed_biomass(bdead, balive, bleaf, pft, hite, bstorage)
      use pft_coms, only:  agf_bs & ! intent(in)
                         , q      & ! intent(in)
                         , qsw    ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real    , intent(in) :: bdead
      real    , intent(in) :: balive
      real    , intent(in) :: bleaf
      real    , intent(in) :: hite
      real    , intent(in) :: bstorage
      integer , intent(in) :: pft
      !----- Local variables --------------------------------------------------------------!
      real                 :: bwood
      !------------------------------------------------------------------------------------!

      bwood      = wood_biomass(bdead, balive, pft, hite)
      ed_biomass = bleaf + bwood

      return
   end function ed_biomass
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function bd2dbh(ipft, bdead)

      use pft_coms, only:  is_tropical & ! intent(in), lookup table
                         , rho         & ! intent(in), lookup table
                         , b1Bs        & ! intent(in), lookup table
                         , b2Bs        & ! intent(in), lookup table
                         , C2B         & ! intent(in)
                         , max_dbh     ! ! intent(in), lookup table

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in) :: ipft
      real   , intent(in) :: bdead
      !----- Local variables --------------------------------------------------------------!
      real :: p
      real :: q
      real :: r
      real :: dbh_pot
      real :: h_max
      !------------------------------------------------------------------------------------!

      if (is_tropical(ipft)) then

         p = a1 + c1d * gg * log(10.0) + d1d * log(rho(ipft))
         r = ( (a2d - a1 ) + gg * log(10.0) * (c2d - c1d) + log(rho(ipft))                 &
             * (d2d - d1d)) /log(dcrit)
         q = 2.0 * b2d + c2d * ff + r
         dbh_pot = (bdead * 2.0 * exp(-p))**(1.0/q)     
         
         if (dbh_pot <= max_dbh(ipft)) then
            bd2dbh = dbh_pot
         else
            h_max = dbh2h(ipft, max_dbh(ipft))
            p = a1 + c1d * log(h_max) + d1d * log(rho(ipft))
            r = ( (a2d - a1 ) + (c2d - c1d) * log(h_max) + log(rho(ipft))                  &
                * (d2d - d1d))/log(dcrit)
            q = 2.0 * b2d + r
            bd2dbh = (bdead * 2.0 * exp(-p))**(1.0/q)
         endif
      else !----- Temperate PFT -----------------------------------------------------------!
         bd2dbh = (bdead / b1Bs(ipft) * C2B)**(1.0/b2Bs(ipft))
      end if

      return
   end function bd2dbh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine estimates the tree area indices, namely leaf, branch(plus twigs), !
   ! and stem.  For the leaf area index (LAI), we use the specific leaf area (SLA), a      !
   ! constant.  The branch/twig and stem (trunk) area indices (BAI and SAI) are found      !
   ! using the model proposed by Järvelä (2004) to find the specific projected area.       !
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
   subroutine area_indices(nplant,bleaf,bdead,balive,dbh,hite,pft,lai,bai,sai)
      use pft_coms    , only : is_tropical     & ! intent(in)
                             , is_grass        & ! intent(in)
                             , rho             & ! intent(in)
                             , C2B             & ! intent(in)
                             , sla             & ! intent(in)
                             , branch_fraction & ! intent(in)
                             , rbranch         & ! intent(in)
                             , rdiamet         & ! intent(in)
                             , rlength         & ! intent(in)
                             , h1stbr          & ! intent(in)
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
      integer , intent(in)  :: pft     ! Plant functional type                [        ---]
      real    , intent(in)  :: nplant  ! Number of plants                     [   plant/m²]
      real    , intent(in)  :: bleaf   ! Specific leaf biomass                [  kgC/plant]
      real    , intent(in)  :: bdead   ! Specific structural                  [  kgC/plant]
      real    , intent(in)  :: balive  ! Specific live tissue biomass         [  kgC/plant]
      real    , intent(in)  :: dbh     ! Diameter at breast height            [         cm]
      real    , intent(in)  :: hite    ! Plant height                         [          m]
      real    , intent(out) :: lai     ! Leaf area index                      [  m²leaf/m²]
      real    , intent(out) :: sai     ! Stem(trunk) projected area index     [  m²wood/m²]
      real    , intent(out) :: bai     ! Branch (+twigs) projected area index [  m²wood/m²]
      !----- Local variables --------------------------------------------------------------!
      real                  :: bwood   ! Wood biomass                         [  kgC/plant]
      real                  :: btrunk  ! Trunk biomass                        [  kgC/plant]
      real                  :: bbranch ! Branch biomass                       [  kgC/plant]
      real                  :: htrunk  ! Trunk height (until first branching) [          m]
      real                  :: sta     ! Specific trunk area                  [   m²/plant]
      real                  :: sba     ! Specific branch area                 [   m²/plant]
      real                  :: bdiamet ! Diameter of current branch           [          m]
      real                  :: blength ! Length of each branch of this order  [          m]
      real                  :: nbranch ! Number of branches of this order     [       ----]
      !----- External functions -----------------------------------------------------------!
      real    , external    :: errorfun ! Error function.
      !------------------------------------------------------------------------------------!
      
      !----- First, we compute the LAI ----------------------------------------------------!
      lai = bleaf * nplant * sla(pft)
      
      !------------------------------------------------------------------------------------!
      !     Here we check whether we need to compute the branch, stem, and effective       !
      ! branch area indices.  These are only needed when branch thermodynamics is used,    !
      ! otherwise, simply assign zeroes to them.                                           !
      !------------------------------------------------------------------------------------!
      select case (ibranch_thermo)
      !----- Ignore branches and trunk. ---------------------------------------------------!
      case (0) 
         bai  = 0.
         sai  = 0.
      !------------------------------------------------------------------------------------!


      !----- Use curve fit based on Conijn (1995) model. ----------------------------------!
      case (1) 
         !---------------------------------------------------------------------------------!
         !     Finding the total wood biomass and the fraction corresponding to branches.  !
         !---------------------------------------------------------------------------------!
         bwood   = wood_biomass(bdead, balive, pft, hite)
         bbranch = branch_fraction(pft) * bwood
         if (is_grass(pft)) then
            sba = conijn_a(pft)
         else
            sba = conijn_a(pft)                                                            &
                + conijn_b(pft) * errorfun(conijn_c(pft)*C2B*bbranch + conijn_d(pft))
         end if
         bai = nplant * bbranch * sba
         sai = 0. !----- No stem information available. -----------------------------------!
      !------------------------------------------------------------------------------------!


      !----- Use . Järvelä (2004) method. -------------------------------------------------!
      case (2) 
         !---------------------------------------------------------------------------------!
         !     Finding the total wood biomass and its partition between bole and branches. !
         !---------------------------------------------------------------------------------!
         bwood   = wood_biomass(bdead, balive, pft, hite)
         bbranch = branch_fraction(pft) * bwood
         btrunk  = (1.-branch_fraction(pft)) * bwood

         !----- Now we check the first branching height -----------------------------------!
         if (is_tropical(pft) .or. is_grass(pft)) then
             !-----------------------------------------------------------------------------!
             !     Htrunk is found assuming the trunk is a cilinder of diameter dbh.  The  !
             ! 40. term comes from the implicit unit convertion:                           !
             ! 1. dbh is diameter, so the cross-section area is pi*dbh²/4;                 !
             ! 2. dbh is given in cm, so dbh²=dbh_SI²/10000.;                              !
             ! 3. wood density is given in g/cm³, and 1g/cm³ = 1000. kg/m³;                !
             !-----------------------------------------------------------------------------!
             htrunk = 40. * btrunk * C2B / (rho(pft) * pi1 * dbh * dbh)
         else
             htrunk = h1stbr(pft) * hite
         end if
         !----- STA is the specific trunk area.  SAI is the stem(trunk) area index. -------!
         sta = ntrunk(pft) * htrunk * (0.01 * dbh)
         sai = nplant * sta

         if (branch_fraction(pft) == 0.) then
            !------------------------------------------------------------------------------!
            !    If there is no branch biomass, skip it and make bai = 0.                  !
            !------------------------------------------------------------------------------!
            sba  = 0.
            bai  = 0.
         else
            !------------------------------------------------------------------------------!
            !     Now we use Järvelä (2004) method to find the lower-order branches        !
            ! contribution to the specific branch area.  Since we don't know how many      !
            ! orders we have beforehand, we use a loop until the diameter gets too small.  !
            !------------------------------------------------------------------------------!
            sba = 0.
            !----- Initialize branch values with trunk. -----------------------------------!
            nbranch = ntrunk(pft)
            bdiamet = dbh * 0.01
            blength = htrunk
            branchloop: do
               if (bdiamet < diammin(pft)) exit branchloop
               !----- Updating branch habits. ---------------------------------------------!
               bdiamet = bdiamet / rdiamet(pft)
               blength = blength / rlength(pft)
               nbranch = nbranch * rbranch(pft)
               sba     = sba + nbranch * blength * bdiamet
            end do branchloop
            !----- The branch area index, similarly to SAI. -------------------------------!
            bai = nplant * sba
         end if
      !------------------------------------------------------------------------------------!

      end select

      return
   end subroutine area_indices
   !=======================================================================================!
   !=======================================================================================!
end module allometry
!==========================================================================================!
!==========================================================================================!


