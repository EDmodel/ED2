real function h2dbh(h,ipft)

  use pft_coms, only: rho, b1Ht, b2Ht

  implicit none
  
  real :: h
  integer :: ipft

  if(rho(ipft).ne.0.0)then  ! Tropical 
     h2dbh = 10.0**((log10(h)-0.37)/0.64)
  else ! Temperate
     h2dbh = log(1.0-(h-1.3)/b1Ht(ipft))/b2Ht(ipft)
  endif

  return
end function h2dbh
!----------------------------------------------------------------
real function dbh2bd(dbh,h,ipft)

  use pft_coms, only: rho, C2B, b1Bs, b2Bs, max_dbh

  implicit none

  real :: dbh,h,p,r,qq
  integer ipft
  real, parameter :: a1 = -1.981
  real, parameter :: b1 = 1.047
  real, parameter :: c1 = 0.572
  real, parameter :: d1 = 0.931
  real, parameter :: dcrit = 100.0
  real, parameter :: a2 = -1.086
  real, parameter :: b2 = 0.876
  real, parameter :: c2 = 0.604
  real, parameter :: d2 = 0.871
  real, parameter :: f = 0.64
  real, parameter :: g = 0.37

  if(rho(ipft).ne.0.0)then

     if(dbh .gt. max_dbh(ipft))then
        p  = a1 + c1 * log(h) + d1 * log(rho(ipft))
        r  = ((a2 - a1) + (c2 - c1)*log(h) + log(rho(ipft))  &
             * (d2 - d1)) * (1/log(dcrit))
        qq = 2.0 * b2 + r
     else
        p  = a1 + c1 * g * log(10.0) + d1 * log(rho(ipft))
        r  = ((a2 - a1) + g * log(10.0) * (c2 - c1) + log(rho(ipft))  &
             * (d2 - d1)) * (1/log(dcrit))
        qq = 2.0 * b2 + c2 * f + r
     endif
  
     dbh2bd = exp(p) / C2B * dbh**qq

  else ! Temperate

     dbh2bd = b1Bs(ipft) / C2B * dbh**b2Bs(ipft)

  endif

  return
end function dbh2bd

!==========================================================================

real function dbh2bl(dbh,ipft)
  use pft_coms, only: rho, max_dbh, C2B, b1Bl, b2Bl

  implicit none

  real :: dbh,mdbh
  integer :: ipft
  real :: p,r,qq
  real, parameter :: a1 = -1.981
  real, parameter :: b1 = 1.047
  real, parameter :: c1 = -0.584
  real, parameter :: d1 = 0.55
  real, parameter :: dcrit = 100.0
  real, parameter :: a2 = -4.111
  real, parameter :: b2 = 0.605
  real, parameter :: c2 = 0.848
  real, parameter :: d2 = 0.438
  real, parameter :: f = 0.64
  real, parameter :: g = 0.37

  if(rho(ipft).ne.0.0)then 
     ! Tropics
     p  = a1 + c1 * g * log(10.0) + d1 * log(rho(ipft))
     r  = ((a2 - a1) + g * log(10.0) * (c2 - c1) + log(rho(ipft))   &
          * (d2 - d1)) * (1/log(dcrit))
     qq = 2.0 * b2 + c2 * f + r  

     if(dbh .le. max_dbh(ipft))then
        dbh2bl = exp(p) / C2B * dbh**qq
     else 
        dbh2bl = exp(p) / C2B * max_dbh(ipft)**qq
     endif
  else
     ! Temperate
     mdbh = min(dbh,max_dbh(ipft))
     dbh2bl = b1Bl(ipft) /C2B * mdbh**b2Bl(ipft)
  endif


  return
end function dbh2bl

!========================================================================

real function calc_root_depth(h,dbh,ipft)
  implicit none 
  real :: h,volume,dbh
  integer :: ipft

  if(ipft.ne.1 .and. ipft.ne.5)then
     volume = h * 0.65 * 3.1415 * (dbh*0.11)**2
     calc_root_depth = -10.0**(0.545 + 0.277*log10(volume))
  else
     ! grasses get a fixed rooting depth of 70 cm.
     calc_root_depth = -0.7
  endif

  return
end function calc_root_depth


!================================================================

integer function assign_root_depth(rd, lsl)

  use grid_coms, only: nzg
  use soil_coms, only: slz

  implicit none

  real, intent(in) :: rd
  integer :: k
  integer, intent(in) :: lsl

  assign_root_depth = nzg
  do k=nzg,lsl+1,-1
     if(rd < slz(k))assign_root_depth = k-1
  enddo

  return
end function assign_root_depth

!=====================================================================

real function dbh2h(ipft, dbh)

  use pft_coms, only: rho, max_dbh, b1Ht, b2Ht

  implicit none

  integer, intent(in) :: ipft
  real, intent(in) :: dbh

  if(rho(ipft) /= 0.0)then

     ! Amazon-type allometry
     if(dbh <= max_dbh(ipft))then

        ! This means that height is below its maximum.
        dbh2h = 10.0 ** (log10(dbh) * 0.64 + 0.37)

     else

        ! Height is at maximum.
        dbh2h = 10.0 ** (log10(max_dbh(ipft)) * 0.64 + 0.37)

     endif

  else

     ! North America-type allometry
     dbh2h = 1.3 + b1Ht(ipft) * (1.0 - exp(b2Ht(ipft) * dbh))

  endif

  return
end function dbh2h

!===================================================================

real function ed_biomass(bdead, balive, bleaf, pft, hite, bstorage)

  use pft_coms, only: agf_bs, q, qsw

  implicit none

  real, intent(in) :: bdead
  real, intent(in) :: balive
  real, intent(in) :: bleaf
  real, intent(in) :: hite
  real, intent(in) :: bstorage
  integer, intent(in) :: pft
  real :: bstem
  real :: bsw

  bstem = agf_bs * bdead
  bsw = agf_bs * balive * qsw(pft) * hite / (1.0 + q(pft) + qsw(pft) * hite)

  ed_biomass = bstem + bleaf + bsw

  return
end function ed_biomass

!=====================================================================

real function bd2dbh(ipft, bdead)

  use pft_coms, only: rho, b1Bs, b2Bs, C2B, max_dbh

  implicit none

  integer, intent(in) :: ipft
  real, intent(in) :: bdead
  real, parameter :: a1 = -1.981
  real, parameter :: b1 = 1.047
  real, parameter :: c1 = 0.572
  real, parameter :: d1 = 0.931
  real, parameter :: dcrit = 100.0
  real, parameter :: a2 = -1.086
  real, parameter :: b2 = 0.876
  real, parameter :: c2 = 0.604
  real, parameter :: d2 = 0.871
  real, parameter :: f = 0.64
  real, parameter :: g = 0.37
  real :: p
  real :: q
  real :: r
  real :: dbh_pot
  real :: h_max
  real, external :: dbh2h

  if(rho(ipft) /= 0.0)then

     p = a1 + c1 * g * log(10.0) + d1 * log(rho(ipft))
     r = ((a2 - a1) + g * log(10.0) * (c2 - c1) + log(rho(ipft)) *   &
          (d2 - d1)) /log(dcrit)
     q = 2.0 * b2 + c2 * f + r
     dbh_pot = (bdead * 2.0 * exp(-p))**(1.0/q)     
     if(dbh_pot <= max_dbh(ipft))then
        bd2dbh = dbh_pot
     else
        h_max = dbh2h(ipft, max_dbh(ipft))
        p = a1 + c1 * log(h_max) + d1 * log(rho(ipft))
        r = ((a2 - a1) + (c2 - c1) * log(h_max) + log(rho(ipft)) *   &
             (d2 - d1))/log(dcrit)
        q = 2 * b2 + r
        bd2dbh = (bdead * 2.0 * exp(-p))**(1.0/q)
     endif

  else

     bd2dbh = (bdead / b1Bs(ipft) * C2B)**(1.0/b2Bs(ipft))
     
  endif

  return
end function bd2dbh


