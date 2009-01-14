!==========================================================================================!
!==========================================================================================!
!    This module is a library with several allometric relationships.                       !
!------------------------------------------------------------------------------------------!
module allometry

   !----- Constants shared by both bdead and bleaf ------------------------------------------!
   real, parameter :: a1    =  -1.981
   real, parameter :: b1    =   1.047
   real, parameter :: dcrit = 100.0
   real, parameter :: ff    =   0.640
   real, parameter :: gg    =   0.370
   !----- Constants used by bdead only ------------------------------------------------------!
   real, parameter :: c1d   =   0.572
   real, parameter :: d1d   =   0.931
   real, parameter :: a2d   =  -1.086
   real, parameter :: b2d   =   0.876
   real, parameter :: c2d   =   0.604
   real, parameter :: d2d   =   0.871
   !----- Constants used by bleaf only ------------------------------------------------------!
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

      use pft_coms, only:  rho     & ! intent(in)
                         , b1Ht    & ! intent(in), lookup table
                         , b2Ht    & ! intent(in), lookup table
                         , hgt_ref ! ! intent(in), lookup table

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      if (rho(ipft) /= 0.0)then  ! Tropical 
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

      use pft_coms, only:  rho     & ! intent(in), lookup table
                         , C2B     & ! intent(in)
                         , b1Bs    & ! intent(in), lookup table
                         , b2Bs    & ! intent(in), lookup table
                         , max_dbh ! ! intent(in), lookup table
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      real   , intent(in) :: h
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: p,r,qq
      !------------------------------------------------------------------------------------!

      if (rho(ipft) /= 0.0) then !----- Tropical PFT --------------------------------------!

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
      use pft_coms, only:  rho     & ! intent(in), lookup table
                         , max_dbh & ! intent(in), lookup table
                         , C2B     & ! intent(in)
                         , b1Bl    & ! intent(in), lookup table
                         , b2Bl    ! ! intent(in), lookup table
   
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real   , intent(in) :: dbh
      integer, intent(in) :: ipft
      !----- Local variables --------------------------------------------------------------!
      real                :: mdbh
      real                :: p,r,qq
      !------------------------------------------------------------------------------------!

      if (rho(ipft) /= 0.0) then !----- Tropical PFT --------------------------------------!
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
         volume          = h * 0.65 * pi1 * (dbh*0.11)**2
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
      use pft_coms, only:  rho      & ! intent(in)
                         , max_dbh  & ! intent(in)
                         , b1Ht     & ! intent(in)
                         , b2Ht     & ! intent(in)
                         , hgt_ref  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer , intent(in) :: ipft
      real    , intent(in) :: dbh
      !------------------------------------------------------------------------------------!

      if (rho(ipft) /= 0.0) then !----- Tropical PFT allometry ----------------------------!
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
      real                 :: bstem
      real                 :: bsw
      !------------------------------------------------------------------------------------!

      bstem      = agf_bs * bdead
      bsw        = agf_bs * balive * qsw(pft) * hite / (1.0 + q(pft) + qsw(pft) * hite)
      ed_biomass = bstem + bleaf + bsw
      return
   end function ed_biomass
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   real function bd2dbh(ipft, bdead)

      use pft_coms, only:  rho     & ! intent(in), lookup table
                         , b1Bs    & ! intent(in), lookup table
                         , b2Bs    & ! intent(in), lookup table
                         , C2B     & ! intent(in)
                         , max_dbh ! ! intent(in), lookup table

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

      if (rho(ipft) /= 0.0) then !----- Tropical PFT --------------------------------------!

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
end module allometry
!==========================================================================================!
!==========================================================================================!


