real function dist_gc(slons,slonf,slats,slatf)
!------------------------------------------------------------------------------------------!
!   Function to compute the great circle distance between two points: the s suffix denotes !
! source point, and f denotes the destination - "forepoint"). The results are given in     !
! metres. The formula is intended to be accurate for both small and large distances and    !
! uses double precision to avoid ill-conditioned behaviour of sin and cos for numbers      !
! close to the n*pi/2.                                                                     !
!------------------------------------------------------------------------------------------!
   use consts_coms, only: erad,pio180
   implicit none
   real, intent(in) :: slons,slonf,slats,slatf
   real(kind=8) :: lons,lonf,lats,latf,dlon,dlat
   real(kind=8) :: x,y
   lons=dble(slons*pio180)
   lonf=dble(slonf*pio180)
   lats=dble(slats*pio180)
   latf=dble(slatf*pio180)
   dlon=lonf-lons
   dlat=latf-lats
   x=dsin(lats)*dsin(latf)+dcos(lats)*dcos(latf)*dcos(dlon)
   y=dsqrt((dcos(latf)*dsin(dlon))*(dcos(latf)*dsin(dlon))+           &
           (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon))*  &
           (dcos(lats)*dsin(latf)-dsin(lats)*dcos(latf)*dcos(dlon)))
   dist_gc = erad*sngl(datan2(y,x))
   return
end function dist_gc

