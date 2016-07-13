module init_hydro_sites
  contains

!!! translate from C to fortran
! ===============================================================


!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calc_flow_routing(cgrid,ipy)
   use ed_state_vars, only: polygontype, edtype
   type(edtype), target :: cgrid ! Alias for current grid
   integer, intent(in) :: ipy    ! Current polygon Polygon ID
   integer :: i,ihys,ines ! ihys -> current hydro site ID. ines-> next hydro site ID
   type(polygontype),pointer :: cpoly ! Alias for current polygon
   real :: side = 10.0  !cell side length (m) used for calculating adjacency
   real(kind=8) :: row_sum
   integer, external :: find_rank
   integer, allocatable, dimension(:) :: hyrank
   integer :: nsites

   !if we loaded an adjacency matrix from file, we don't need to calculate flow routing
   if(cgrid%load_adjacency(ipy) == 1) return   

   cpoly => cgrid%polygon(ipy)
   nsites=cpoly%nsites
   
   allocate (hyrank(nsites))

   ! set flow routing -> sort by TCI, it will give the sites the hydro order
   call rank_up(nsites,cpoly%TCI,hyrank)


   !init structures and recalculate sitenums list based on hydro order
   do ihys=1,nsites
      do ines=1,nsites
         cgrid%site_adjacency(ihys,ines,ipy) = 0.0
      end do
   end do
   !routing established, now calc approximate adjacency matrix
   do i=1,nsites-1
      ihys=find_rank(i,nsites,hyrank)
      ines=find_rank(i+1,nsites,hyrank)
      cgrid%site_adjacency(ihys,ihys,ipy) = 2.0*side*(side-1.0)*cpoly%area(ihys) &
               + side/2.0*cpoly%area(ihys)/(cpoly%area(ihys)+cpoly%area(ines))
      cgrid%site_adjacency(ihys,ines,ipy) = side*cpoly%area(ines)/(cpoly%area(ihys)+cpoly%area(ihys))
   end do
   !Now find the adjacency for the last rank site:
   ihys=find_rank(nsites,nsites,hyrank)
   ines=find_rank(1,nsites,hyrank)
   cgrid%site_adjacency(ihys,ihys,ipy) = 2.0*side*(side-1.0)*cpoly%area(ihys) &
            + side/2.0*cpoly%area(ihys)/(cpoly%area(ihys)+cpoly%area(ines))
   cgrid%site_adjacency(ihys,ines,ipy) = side*cpoly%area(ines)/(cpoly%area(ihys)+cpoly%area(ihys))

   !normalize routing (rows sum to 1)  
   do ihys=1,nsites
      row_sum=sum(dble(cgrid%site_adjacency(ihys,:,ipy)))
      cgrid%site_adjacency(ihys,:,ipy)=real(dble(cgrid%site_adjacency(ihys,:,ipy))/row_sum)
   end do

   !check routing
!   print*,"CHECK ROUTING" 
!   do i=1,myPolygon%nsites
!      print*,cgrid%site_adjacency(i,:,ipy)
!   end do

   deallocate(hyrank)
   return
end subroutine calc_flow_routing


end module init_hydro_sites