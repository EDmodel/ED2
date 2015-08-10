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
