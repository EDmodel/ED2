!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine hemintrp_cof(n2,n3,qlatt,qlatu,qlatv,qlont,qlonu,qlonv)

use mem_grid
use grid_dims, only: maxhp

implicit none
integer :: n2,n3
real :: qlatt(n2,n3),qlatu(n2,n3),qlatv(n2,n3),qlont(n2,n3)  &
   ,qlonu(n2,n3),qlonv(n2,n3)

integer :: i,j

! Initialize interpolation coefficients for hemispheric grid communication

!cccccccccccccccccccccccccc
!c          RETURN
!ccccccccccccccccccccccccccc

if (nhemgrd2 .le. 1) return

call newgrid(1)
do j = 1,n3
   do i = 1,n2

! Get coordinates (qlat,qlon), the "rotated" latitude and longitude relative
! to the hemispheric grid (#1) pole point, for t, u, and v points.  If the
! pole point is not located at a geographic pole, qlat and qlon are not
! the geographic latitude and longitude.

      call xy_ll(qlatt(i,j),qlont(i,j),90.,0.,xt(i),yt(j))
      call xy_ll(qlatu(i,j),qlonu(i,j),90.,0.,xm(i),yt(j))
      call xy_ll(qlatv(i,j),qlonv(i,j),90.,0.,xt(i),ym(j))

   enddo
enddo
! Get set of points from opposite hemisphere that need to be interpolated and
! the interpolation points and weights from the present hemisphere for
! t, u, and v points.

! Interpolation of T points to T points.

call hemintrp_table (n2,n3,maxhp,qlatt,qlont,qlatt,qlont  &
   ,xt(1),yt(1),deltax,deltay  &
   ,nhemt,ihem1tt,jhem1tt,whem1tt,ihem2tt,jhem2tt,hlatt,hlont)

! Interpolation of U points to U points.

call hemintrp_table (n2,n3,maxhp,qlatu,qlonu,qlatu,qlonu  &
   ,xm(1),yt(1),deltax,deltay  &
   ,nhemu,ihem1uu,jhem1uu,whem1uu,ihem2uu,jhem2uu,hlatu,hlonu)

! Interpolation of V points to U points.

call hemintrp_table (n2,n3,maxhp,qlatv,qlonv,qlatu,qlonu  &
   ,xt(1),ym(1),deltax,deltay  &
   ,nhemu,ihem1vu,jhem1vu,whem1vu,ihem2vu,jhem2vu,hlatu,hlonu)

! Interpolation of U points to V points.

call hemintrp_table (n2,n3,maxhp,qlatu,qlonu,qlatv,qlonv  &
   ,xm(1),yt(1),deltax,deltay  &
   ,nhemv,ihem1uv,jhem1uv,whem1uv,ihem2uv,jhem2uv,hlatv,hlonv)

! Interpolation of V points to V points.

call hemintrp_table (n2,n3,maxhp,qlatv,qlonv,qlatv,qlonv  &
   ,xt(1),ym(1),deltax,deltay  &
   ,nhemv,ihem1vv,jhem1vv,whem1vv,ihem2vv,jhem2vv,hlatv,hlonv)

return
end subroutine hemintrp_cof

!     ****************************************************************

subroutine hemintrp_table (n2,n3,maxhp,qlat1,qlon1,qlat2,qlon2  &
                          ,xmt1,ymt1,deltax,deltay  &
                          ,ind,ihem1,jhem1,whem1,ihem2,jhem2,hlat,hlon)
  implicit none
  integer, intent(in)         :: n2,n3,maxhp
  real, dimension(n2,n3)      :: qlat1,qlon1,qlat2,qlon2
  real, dimension(4,maxhp)    :: whem1
  real, dimension(maxhp)      :: hlat,hlon
  integer, dimension(4,maxhp) :: ihem1,jhem1
  integer, dimension(maxhp)   :: ihem2,jhem2
  real,intent(in)             :: xmt1,ymt1,deltax,deltay
  integer, intent(out)        :: ind
  ! Local variables
  integer :: i,j,ntri,ii,jj,i1,i2,i3,i4,j1,j2,j3,j4,iyes
  real    :: x,x0,x1,x2,x3,x4,x5,y,y0,y1,y2,y3,y4,y5,w34

  ! Loop through (i,j) points of hemisphere 1 that need to be interpolated
  ! (those with rotated latitude less than 0 degrees).  Compute (lat,lon)
  ! coordinates of corresponding (i,j) points of hemisphere 2, as reflections
  ! hemisphere 1 points, and store these in (hlat,hlon).  Compute the (x,y)
  ! coordinates of those points in hemisphere 1 for rapid searching for a
  ! local set of triangles for interpolation.

  ind = 0
  do j = 1,n3
     do i = 1,n2

  ! T points: interpolation from 4 T points

        if (qlat2(i,j) .lt. 0.) then
           ind = ind + 1

           hlat(ind) =  - qlat2(i,j)
           hlon(ind) =  - qlon2(i,j)

           if (hlon(ind) .le. -180.) hlon(ind) = hlon(ind) + 360.

           call ll_xy(hlat(ind),hlon(ind),90.,0.,x,y)

  ! Estimate (ii,jj) index of "ct" point in which (x,y) coordinates lie.

           ii = int ((x - xmt1) / deltax) + 1
           jj = int ((y - ymt1) / deltay) + 1

  ! search over 8 triangles surrounding "ct" point.

           triloop: do ntri = 1,8
              call get_triinds (ntri,ii,jj,i1,i2,i3,i4,j1,j2,j3,j4)

              x0 = hlat(ind)
              x1 = qlat1(i1,j1)
              x2 = qlat1(i2,j2)
              x3 = qlat1(i3,j3)
              x4 = qlat1(i4,j4)
              x5 = .25 * (x1 + x2 + x3 + x4)

              y0 = hlon(ind)
              y1 = qlon1(i1,j1)
              y2 = qlon1(i2,j2)
              y3 = qlon1(i3,j3)
              y4 = qlon1(i4,j4)

  ! If {y1,y2,y3,y4} crosses the 180 degree meridian, add 360 to all negative
  ! values.

              if (max(y1,y2,y3,y4) .gt. 90. .and.  &
                  min(y1,y2,y3,y4) .lt. -90.) then

                 if (y0 .lt. 0.) y0 = y0 + 360.
                 if (y1 .lt. 0.) y1 = y1 + 360.
                 if (y2 .lt. 0.) y2 = y2 + 360.
                 if (y3 .lt. 0.) y3 = y3 + 360.
                 if (y4 .lt. 0.) y4 = y4 + 360.

              endif
              y5 = .25 * (y1 + y2 + y3 + y4)

              call tricheck (x0,y0,x1,x2,x5,y1,y2,y5,iyes)
              if (iyes .eq. 1) exit triloop
           end do triloop
           if (iyes /= 1) then
             print*, 'no interpolation triangles found for T point'
             print*, 'i,j',i,j
             stop 'hemintrp_table'
           end if

  ! Get interpolation indices and coefficients.

           ihem2(ind) = i
           jhem2(ind) = j

           ihem1(1,ind) = i1
           ihem1(2,ind) = i2
           ihem1(3,ind) = i3
           ihem1(4,ind) = i4

           jhem1(1,ind) = j1
           jhem1(2,ind) = j2
           jhem1(3,ind) = j3
           jhem1(4,ind) = j4

           call tri_intp_cofs (x0,y0,x1,x2,x5,y1,y2,y5  &
              ,whem1(1,ind),whem1(2,ind),whem1(3,ind))

           w34 = .25 * whem1(3,ind)

           whem1(1,ind) = whem1(1,ind) + w34
           whem1(2,ind) = whem1(2,ind) + w34
           whem1(3,ind) = w34
           whem1(4,ind) = w34

        endif
     enddo
  enddo

  return
end

!     ****************************************************************

subroutine get_triinds (ntri,ii,jj,i1,i2,i3,i4,j1,j2,j3,j4)
implicit none
integer :: ntri,ii,jj,i1,i2,i3,i4,j1,j2,j3,j4

i1 = ii
i2 = ii
i3 = ii
i4 = ii

j1 = jj
j2 = jj
j3 = jj
j4 = jj

if (ntri .eq. 1) then
   i1 = ii+1
   i2 = ii+1
   j2 = jj+1
   j3 = jj+1
elseif (ntri .eq. 2) then
   i1 = ii+1
   i4 = ii+1
   j1 = jj+1
   j2 = jj+1
elseif (ntri .eq. 3) then
   i3 = ii+1
   i4 = ii+1
   j1 = jj+1
   j4 = jj+1
elseif (ntri .eq. 4) then
   i2 = ii+1
   i3 = ii+1
   j3 = jj+1
   j4 = jj+1
elseif (ntri .eq. 5) then
   i1 = ii+1
   i2 = ii+1
   i3 = ii+2
   i4 = ii+2
   j1 = jj+1
   j4 = jj+1
elseif (ntri .eq. 6) then
   i2 = ii+1
   i3 = ii+1
   j1 = jj+1
   j2 = jj+1
   j3 = jj+2
   j4 = jj+2
elseif (ntri .eq. 7) then
   i3 = ii-1
   i4 = ii-1
   j2 = jj+1
   j3 = jj+1
elseif (ntri .eq. 8) then
   i1 = ii+1
   i4 = ii+1
   j3 = jj-1
   j4 = jj-1
endif

return
end

!     ****************************************************************

subroutine tricheck (x0,y0,x1,x2,x5,y1,y2,y5,iyes)
implicit none
integer :: iyes
real :: x0,y0,x1,x2,x5,y1,y2,y5

! check to see if point (x0,y0) is located within planar triangle defined
! by vertices (x1,y1), (x2,y2), (x5,y5) which are listed in counterclockwise
! order.
! 6/9/99:  Allow point (x0,y0) to be up to 1.e-4 degrees (10 meters) outside
! triangle to accomodate roundoff error in subroutine xy_ll.

iyes = 0

if ((x2-x1) * (y0-y1) - (x0-x1) * (y2-y1) .le. 1.e-4 .and.  &
    (x5-x2) * (y0-y2) - (x0-x2) * (y5-y2) .le. 1.e-4 .and.  &
    (x1-x5) * (y0-y5) - (x0-x5) * (y1-y5) .le. 1.e-4) iyes = 1

!       print*, 'one', (x2-x1) * (y0-y1) - (x0-x1) * (y2-y1)
!       print*, 'two', (x5-x2) * (y0-y2) - (x0-x2) * (y5-y2)
!       print*, 'tre', (x1-x5) * (y0-y5) - (x0-x5) * (y1-y5)

return
end

!     ****************************************************************

subroutine tri_intp_cofs (x,y,x1,x2,x5,y1,y2,y5,tic1,tic2,tic5)
implicit none
real :: x,y,x1,x2,x5,y1,y2,y5,tic1,tic2,tic5,area

area = (x5 - x2) * (y1 - y5) - (x1 - x5) * (y5 - y2)

tic1 = ((x5-x2) * (y-y2) - (y5-y2) * (x-x2)) / area
tic2 = ((x1-x5) * (y-y1) - (y1-y5) * (x-x1)) / area
tic5 = ((x2-x1) * (y-y1) - (y2-y1) * (x-x1)) / area

return
end

!     *******************************************************************

subroutine hemintrp

use mem_grid
use mem_basic
use var_tables

implicit none

integer :: nscal

! Interpolate prognostic variables between hemispheric grids.
! This subroutine is for non-parallel runs only.

if (nhemgrd2 .le. 1) return

! Interpolate grid 1 interior values to grid nhemgrd2 boundaries.

call hemintuv(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(1)%up(1,1,1)         &
   ,basic_g(1)%vp(1,1,1)         &
   ,basic_g(nhemgrd2)%up(1,1,1)  &
   ,basic_g(nhemgrd2)%vp(1,1,1)  )

call hemintuv(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(1)%uc(1,1,1)         &
   ,basic_g(1)%vc(1,1,1)         &
   ,basic_g(nhemgrd2)%uc(1,1,1)  &
   ,basic_g(nhemgrd2)%vc(1,1,1)  )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(1)%wp(1,1,1)         &
   ,basic_g(nhemgrd2)%wp(1,1,1)  )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(1)%wc(1,1,1)         &
   ,basic_g(nhemgrd2)%wc(1,1,1)  )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(1)%pp(1,1,1)         &
   ,basic_g(nhemgrd2)%pp(1,1,1)  )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(1)%pc(1,1,1)         &
   ,basic_g(nhemgrd2)%pc(1,1,1)  )

do nscal = 1,num_scalar(1)
   call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
      ,scalar_tab(nscal,1)%var_p         &
      ,scalar_tab(nscal,nhemgrd2)%var_p  )
enddo

! Interpolate grid nhemgrd2 interior values to grid 1 boundaries.

call hemintuv(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(nhemgrd2)%up(1,1,1)  &
   ,basic_g(nhemgrd2)%vp(1,1,1)  &
   ,basic_g(1)%up(1,1,1)         &
   ,basic_g(1)%vp(1,1,1)         )

call hemintuv(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(nhemgrd2)%uc(1,1,1)  &
   ,basic_g(nhemgrd2)%vc(1,1,1)  &
   ,basic_g(1)%uc(1,1,1)         &
   ,basic_g(1)%vc(1,1,1)         )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(nhemgrd2)%wp(1,1,1)  &
   ,basic_g(1)%wp(1,1,1)         )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(nhemgrd2)%wc(1,1,1)  &
   ,basic_g(1)%wc(1,1,1)         )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(nhemgrd2)%pp(1,1,1)  &
   ,basic_g(1)%pp(1,1,1)         )

call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
   ,basic_g(nhemgrd2)%pc(1,1,1)  &
   ,basic_g(1)%pc(1,1,1)         )

do nscal = 1,num_scalar(1)
   call hemintt(nnzp(1),nnxp(1),nnyp(1)  &
      ,scalar_tab(nscal,nhemgrd2)%var_p  &
      ,scalar_tab(nscal,1)%var_p )
enddo

return
end

!     ****************************************************************

subroutine hemintuv(n1,n2,n3,u1,v1,u2,v2)

use mem_grid

implicit none
integer :: n1,n2,n3
real :: u1(n1,n2,n3),v1(n1,n2,n3),u2(n1,n2,n3),v2(n1,n2,n3)

integer :: ij,k
real :: u2r,v2r,ue,ve,u,v

do ij = 1,nhemu
   do k = 1,n1

      u2r = whem1uu(1,ij) * u1(k,ihem1uu(1,ij),jhem1uu(1,ij))  &
          + whem1uu(2,ij) * u1(k,ihem1uu(2,ij),jhem1uu(2,ij))  &
          + whem1uu(3,ij) * u1(k,ihem1uu(3,ij),jhem1uu(3,ij))  &
          + whem1uu(4,ij) * u1(k,ihem1uu(4,ij),jhem1uu(4,ij))

      v2r = whem1vu(1,ij) * v1(k,ihem1vu(1,ij),jhem1vu(1,ij))  &
          + whem1vu(2,ij) * v1(k,ihem1vu(2,ij),jhem1vu(2,ij))  &
          + whem1vu(3,ij) * v1(k,ihem1vu(3,ij),jhem1vu(3,ij))  &
          + whem1vu(4,ij) * v1(k,ihem1vu(4,ij),jhem1vu(4,ij))

      call uvtoueve(u2r,v2r,ue,ve,hlatu(ij),hlonu(ij)  &
         ,90.,0.)
      call uevetouv(u,v,ue,ve,hlatu(ij),hlonu(ij)  &
         ,-90.,180.)

      u2(k,ihem2uu(ij),jhem2uu(ij)) = u

   enddo
enddo

do ij = 1,nhemv
   do k = 1,n1

      u2r = whem1uv(1,ij) * u1(k,ihem1uv(1,ij),jhem1uv(1,ij))  &
          + whem1uv(2,ij) * u1(k,ihem1uv(2,ij),jhem1uv(2,ij))  &
          + whem1uv(3,ij) * u1(k,ihem1uv(3,ij),jhem1uv(3,ij))  &
          + whem1uv(4,ij) * u1(k,ihem1uv(4,ij),jhem1uv(4,ij))

      v2r = whem1vv(1,ij) * v1(k,ihem1vv(1,ij),jhem1vv(1,ij))  &
          + whem1vv(2,ij) * v1(k,ihem1vv(2,ij),jhem1vv(2,ij))  &
          + whem1vv(3,ij) * v1(k,ihem1vv(3,ij),jhem1vv(3,ij))  &
          + whem1vv(4,ij) * v1(k,ihem1vv(4,ij),jhem1vv(4,ij))

      call uvtoueve(u2r,v2r,ue,ve,hlatv(ij),hlonv(ij),90.,0.)
      call uevetouv(u,v,ue,ve,hlatv(ij),hlonv(ij),-90.,180.)

      v2(k,ihem2vv(ij),jhem2vv(ij)) = v

   enddo
enddo
return
end

!     ****************************************************************

subroutine hemintt(n1,n2,n3,t1,t2)

use mem_grid

implicit none
integer :: n1,n2,n3
real :: t1(n1,n2,n3),t2(n1,n2,n3)

integer :: ij,k

do ij = 1,nhemt
   do k = 1,n1

      t2(k,ihem2tt(ij),jhem2tt(ij))  &
         = whem1tt(1,ij) * t1(k,ihem1tt(1,ij),jhem1tt(1,ij))  &
         + whem1tt(2,ij) * t1(k,ihem1tt(2,ij),jhem1tt(2,ij))  &
         + whem1tt(3,ij) * t1(k,ihem1tt(3,ij),jhem1tt(3,ij))  &
         + whem1tt(4,ij) * t1(k,ihem1tt(4,ij),jhem1tt(4,ij))

   enddo
enddo
return
end

