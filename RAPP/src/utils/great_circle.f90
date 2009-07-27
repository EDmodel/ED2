!==========================================================================================!
!==========================================================================================!
!   Function to compute the great circle distance between two points: the s suffix denotes !
! source point, and f denotes the destination - "forepoint"). The results are given in     !
! metres. The formula is intended to be accurate for both small and large distances and    !
! uses double precision to avoid ill-conditioned behaviour of sin and cos for numbers      !
! close to the n*pi/2.                                                                     !
!------------------------------------------------------------------------------------------!
real function dist_gc(slons,slonf,slats,slatf)
   use rconstants , only : erad    & ! intent(in)
                         , pio180  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
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


   !----- Converting the coordinates to radians, in double precision. ---------------------!
   lons=dble(slons * pio180)
   lonf=dble(slonf * pio180)
   lats=dble(slats * pio180)
   latf=dble(slatf * pio180)
   !---------------------------------------------------------------------------------------!

   !----- Finding the delta-lon and delta-lat. --------------------------------------------!
   dlon=lonf-lons
   dlat=latf-lats

   x=dsin(lats)*dsin(latf)+dcos(lats)*dcos(latf)*dcos(dlon)
   y=dsqrt((dcos(latf)*dsin(dlon))*(dcos(latf)*dsin(dlon))+           &
           (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon))*  &
           (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon)))
   dist_gc = erad*sngl(datan2(y,x))

   return
end function dist_gc
!==========================================================================================!
!==========================================================================================!
