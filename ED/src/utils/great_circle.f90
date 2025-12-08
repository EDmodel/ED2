!==========================================================================================!
!==========================================================================================!
!   Function to compute the great circle distance between two points: the s suffix denotes !
! source point, and f denotes the destination - "forepoint"). The results are given in     !
! metres. The formula is intended to be accurate for both small and large distances and    !
! uses double precision to avoid ill-conditioned behaviour of sin and cos for numbers      !
! close to the n*pi/2.                                                                     !
!------------------------------------------------------------------------------------------!
real function dist_gc(slons,slonf,slats,slatf)
   use consts_coms, only : erad    & ! intent(in)
                         , pio1808 ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real, intent(in) :: slons
   real, intent(in) :: slonf
   real, intent(in) :: slats
   real, intent(in) :: slatf
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)     :: lons
   real(kind=8)     :: lonf
   real(kind=8)     :: lats
   real(kind=8)     :: latf
   real(kind=8)     :: dlon
   real(kind=8)     :: dlat
   real(kind=8)     :: x
   real(kind=8)     :: y
   !---------------------------------------------------------------------------------------!

   !----- Convert the co-ordinates to double precision and to radians. --------------------!
   lons = dble(slons) * pio1808
   lonf = dble(slonf) * pio1808
   lats = dble(slats) * pio1808
   latf = dble(slatf) * pio1808
   dlon = lonf - lons
   dlat = latf - lats

   !----- Find the arcs. ------------------------------------------------------------------!
   x    = dsin(lats) * dsin(latf) + dcos(lats) * dcos(latf) * dcos(dlon)
   y    = dsqrt( (dcos(latf)*dsin(dlon)) * (dcos(latf)*dsin(dlon))                         &
               + (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon))                  &
               * (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon)) )

   !----- Convert the arcs to actual distance. --------------------------------------------!
   dist_gc = erad*sngl(datan2(y,x))

   return
end function dist_gc
!==========================================================================================!
!==========================================================================================!



!==========================================================================================!
!==========================================================================================!
!   Function to compute the area associated with a solid angle defined by a pair of co-    !
! ordinates corresponding to the SW and NE corners of a "rectangle" in longitude/latitude  !
! coordinates.                                                                             !
!------------------------------------------------------------------------------------------!
real function solid_area(wlon,slat,elon,nlat)
   use consts_coms, only : erad8     & ! intent(in)
                         , pio1808   & ! intent(in)
                         , tiny_num8 ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4), intent(in) :: wlon
   real(kind=4), intent(in) :: slat
   real(kind=4), intent(in) :: elon
   real(kind=4), intent(in) :: nlat
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)            :: wlon8
   real(kind=8)            :: slat8
   real(kind=8)            :: elon8
   real(kind=8)            :: nlat8
   real(kind=8)            :: area8
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external  :: sngloff
   !---------------------------------------------------------------------------------------!

   !----- Convert the co-ordinates to double precision and to radians. --------------------!
   wlon8 = dble(wlon) * pio1808
   elon8 = dble(elon) * pio1808
   slat8 = dble(slat) * pio1808
   nlat8 = dble(nlat) * pio1808
   !---------------------------------------------------------------------------------------!


   !----- Find the double-precision area, and convert to single precision for output. -----!
   area8      = (sin(nlat8) - sin(slat8)) * (elon8 - wlon8) * erad8 * erad8
   solid_area = sngloff(area8,tiny_num8)
   !---------------------------------------------------------------------------------------!

   return
end function solid_area
!==========================================================================================!
!==========================================================================================!
