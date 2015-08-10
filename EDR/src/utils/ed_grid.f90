!==========================================================================================!
!==========================================================================================!
subroutine ed_gridset(ngra)


  ! ed_gridset: grid initialization or moving a nested grid.
  !   initialization iff ngra = 1
  !   moving grid "ngra" ow

  use grid_coms, only: &
       jdim,           &
       ngrids,         &
       polelat,        &
       polelon,        &
       platn,          &
       plonn,          &
       centlat,        &
       centlon,        &
       deltax,         &
       deltaxn,        &
       deltay,         &
       deltayn,        &
       nstratx,        &
       nstraty,        &
       ninest,         &
       xmn,            &
       xtn,            &
       njnest,         &
       ymn,            &
       ytn,            &
       nxtnest,        &
       nnxp,           &
       nnyp


  use ed_node_coms,only : mxp,myp,mmxp,mmyp

  use ed_max_dims, only: &
       maxgrds,        &
       nzpmax

  use consts_coms, only: erad

  implicit none

  integer, intent(in) :: ngra

  integer :: ifm
  integer :: icm
  integer :: i
  integer :: j
  integer :: ngrb
  real :: centx1
  real :: centy1
  real :: centx
  real :: centy

  ! Grid initialization:

  if (ngra == 1) then

     ! Set platn plonn arrays, the values of polelat and polelon for each grid

     do ifm = 1,ngrids
        icm = nxtnest(ifm)
        if (ifm == 1) then
           platn(ifm) = polelat
           plonn(ifm) = polelon
        else
           platn(ifm) = platn(icm)
           plonn(ifm) = plonn(icm)
        endif
     enddo

     ! just for the coarser grid:
     ! set deltax and deltay (deltaxn and deltayn)
     ! find grid center point coordinates at polar stereographic projection
     ! from these, find coordinates of the lower left coarser grid cell
     ! (xmn, ymn) are coordinates of the "higher" cell boundary
     ! these coordinates are required for computing ninest, njnest of nested grids

     deltaxn(1) = deltax
     deltayn(1) = deltay
     call ed_ll_xy(centlat(1),centlon(1),platn(1),plonn(1),centx1,centy1)
     xmn(1,1) = centx1 - 0.5 * float(nnxp(1)-2) * deltaxn(1)
     ymn(1,1) = centy1 - 0.5 * float(nnyp(1)-2) * deltayn(1)


     if (nnyp(1) == 1) ymn(1,1) = centy1

     ! for all nested grids:
     ! set deltaxn and deltayn of the nested grid as multiples of deltaxn and
     ! deltayn of the coarser grid.
     ! convert (centlat, centlon) into (ninest, njnest) or use given (ninest, njnest);
     ! with (ninest, njnest), place all interior cells of the nested grid 
     ! fully covering inner cells of the coarser grid. This is guaranteed by
     ! placing the first cell of the nested grid just outside one cell of the
     ! coarser grid and by the fact that the number of inner cells of the nested
     ! grid is a multiple of nstartx and nstraty (as enforced by opspec2)

     do ifm = 1,ngrids       ! fine grid number
        icm = nxtnest(ifm)   ! coarse grid number
        if (icm >= 1) then
           deltaxn(ifm) = deltaxn(icm) / float(nstratx(ifm))
           deltayn(ifm) = deltayn(icm) / float(nstraty(ifm))
           if (ninest(ifm) <= 0 .or. njnest(ifm) <= 0)  &
                call ed_ll_xy(centlat(ifm),centlon(ifm),platn(ifm)  &
                ,plonn(ifm),centx,centy)
           if (ninest(ifm) <= 0) then
              xmn(1,ifm) = centx - 0.5 * float(nnxp(ifm)-2) * deltaxn(ifm)
              ninest(ifm) = int((xmn(1,ifm) - xmn(1,icm)) / deltaxn(icm) + 1.5)
           endif
           if (njnest(ifm) <= 0) then
              ymn(1,ifm) = centy - 0.5 * float(nnyp(ifm)-2) * deltayn(ifm)
              njnest(ifm) = int((ymn(1,ifm) - ymn(1,icm)) / deltayn(icm) + 1.5)
           endif
           xmn(1,ifm) = xmn(1,icm) + deltaxn(icm) * float(ninest(ifm)-1)
           ymn(1,ifm) = ymn(1,icm) + deltayn(icm) * float(njnest(ifm)-1)
        endif
     enddo
  endif

  ! Given the first point of xmn and ymn for the outer
  ! grid, compute remaining outer grid cell locations


  if (ngra == 1) then
     call ed_newgrid(1)

     do i = 2,mxp
        xmn(i,1) = xmn(i-1,1) + deltaxn(1)
     enddo

     do j = 2,myp
        ymn(j,1) = ymn(j-1,1) + deltayn(1)
     end do

  endif

  ! compute xmn and ymn for any required nested grids, and xtn and ytn for
  ! any required grids.

  ngrb = ngrids
  if (ngra > 1) ngrb = ngra

  do ifm = ngra,ngrb
     call ed_newgrid(ifm)
     icm = nxtnest(ifm)
     if (icm >= 1) then
        xmn(1,ifm) = xmn(ninest(ifm),icm)
        ymn(1,ifm) = ymn(njnest(ifm),icm)
     endif

     do i = 2,mxp
        xmn(i,ifm) = xmn(i-1,ifm) + deltaxn(ifm)
        xtn(i,ifm) = .5 * (xmn(i,ifm) + xmn(i-1,ifm))
     enddo
     xtn(1,ifm) = 1.5 * xmn(1,ifm) - .5 * xmn(2,ifm)

     if (jdim == 1) then
        do j = 2,myp
           ymn(j,ifm) = ymn(j-1,ifm) + deltayn(ifm)
           ytn(j,ifm) = .5 * (ymn(j,ifm) + ymn(j-1,ifm))
        enddo
        ytn(1,ifm) = 1.5 * ymn(1,ifm) - .5 * ymn(2,ifm)
     else
        ytn(1,ifm) = ymn(1,ifm)
     endif
  enddo

  return
end subroutine ed_gridset
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_polarst(n2,n3,glat,glon)

   use grid_coms, only: xt,yt,platn,plonn,ngrid

   implicit none

   integer , intent(in)                    :: n2,n3
   real    , intent(out), dimension(n2,n3) :: glat,glon

   integer :: i,j

   !  Calculates map factors and inverse map factors at u,v,t,m-points and
   !  geographical lat/lon at t-points for a given polar stereographic grid

   do j = 1,n3
      do i = 1,n2
         call ed_xy_ll(glat(i,j),glon(i,j),platn(ngrid),plonn(ngrid),xt(i),yt(j))
      enddo
   enddo

end subroutine ed_polarst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_newgrid(ngr)

   use grid_coms, only : ngrid,nxp,nnxp,nyp,nnyp,deltax,deltaxn,xt,xtn,xm,xmn              &
                        ,deltay,deltayn,yt,ytn,ym,ymn
   use ed_node_coms , only : mxp,myp,mmxp,mmyp,ia,mia,iz,miz,ja,mja,jz,mjz,i0,mi0,j0,mj0   &
                            ,ibcon,mibcon

   implicit none
   integer :: ngr

   integer :: i,j

   !     +----------------------------------------------------------------
   !     !    Fill the single and 1D variables that the rest of the model
   !     !      uses from the nest arrays and change grid level in the I/O.
   !     +----------------------------------------------------------------

   ngrid = ngr

   !         grid point references
   !         x - direction

   nxp=nnxp(ngr)
   nyp=nnyp(ngr)
   !
   !          grid spacings

   deltax=deltaxn(ngr)
   do i=1,nxp
      xt(i)=xtn(i,ngr)
      xm(i)=xmn(i,ngr)
   enddo
   !
   deltay=deltayn(ngr)
   do j=1,nyp
      yt(j)=ytn(j,ngr)
      ym(j)=ymn(j,ngr)
   enddo
   !
   !        node gridpoint info

   mxp=mmxp(ngr)
   myp=mmyp(ngr)

   ia=mia(ngr)
   iz=miz(ngr)

   ja=mja(ngr)
   jz=mjz(ngr)
  
   i0=mi0(ngr)
   j0=mj0(ngr)
   ibcon=mibcon(ngr)

   return
end subroutine ed_newgrid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_xy_ll (qlat,qlon,polelat,polelon,x,y)
  use consts_coms, only: erad,erad2,pio180
  implicit none
  real, intent(out) :: qlat
  real, intent(out) :: qlon
  real, intent(in ) :: polelat
  real, intent(in ) :: polelon
  real, intent(in ) :: x
  real, intent(in ) :: y

  real :: sinplat
  real :: cosplat
  real :: sinplon
  real :: cosplon
  real :: x3p
  real :: y3p
  real :: z3p
  real :: z3q
  real :: x3q
  real :: y3q
  real :: xq
  real :: yq
  real :: zq
  real :: t
  real :: d
  real :: alpha
  real :: r3q


  ! Evaluate sine and cosine of latitude and longitude of pole point p.

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (x3,y3,z3) coordinates of the pole point where the origin is the
  ! center of the earth, the z axis is the north pole, the x axis is the
  ! equator and prime meridian, and the y axis is the equator and 90 E.

  x3p = erad * cosplat * cosplon
  y3p = erad * cosplat * sinplon
  z3p = erad * sinplat

  ! Compute distance d from given point R on the polar stereographic plane
  ! to the pole point P:

  d = sqrt (x ** 2 + y ** 2)

  ! Compute angle QCP where C is the center of the Earth.  This is twice
  ! angle QAP where A is the antipodal point.  Angle QAP is the same as
  ! angle RAP:

  alpha = 2. * atan2(d,erad2)

  ! Compute zq, the height of Q relative to the polar stereographic plane:

  zq = erad * (cos(alpha) - 1.)

  ! Compute the parameter t which is the the distance ratio AQ:AR

  t = (erad2 + zq) / erad2

  ! Compute xq and yq, the x and y coordinates of Q in polar stereographic space:

  xq = t * x
  yq = t * y

  ! Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):

  x3q = x3p - xq * sinplon - yq * cosplon * sinplat  &
       + zq * cosplat * cosplon
  y3q = y3p + xq * cosplon - yq * sinplon * sinplat  &
       + zq * cosplat * sinplon
  z3q = z3p + yq * cosplat + zq * sinplat

  ! Compute the latitude and longitude of Q:

  qlon = atan2(y3q,x3q) / pio180
  r3q = sqrt(x3q ** 2 + y3q ** 2)
  qlat = atan2(z3q,r3q) / pio180

end subroutine ed_xy_ll
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_ll_xy (qlat,qlon,polelat,polelon,x,y)
  use consts_coms, only : erad,erad2,pio180
  ! ll_xy: projects earth surface point (given by latitude,longitude in degrees)
  !        into a plane tangent to earths surface at a pole point (given by
  !        latitude, longitude in degrees) using a polar stereographic projection;
  !        returns point position in the projected cartesian system (in meters).

  implicit none
  real, intent(in ) :: qlat
  real, intent(in ) :: qlon
  real, intent(in ) :: polelat
  real, intent(in ) :: polelon
  real, intent(out) :: x
  real, intent(out) :: y

  real :: sinplat
  real :: cosplat
  real :: sinplon
  real :: cosplon
  real :: sinqlat
  real :: cosqlat
  real :: sinqlon
  real :: cosqlon
  real :: x3p
  real :: y3p
  real :: z3p
  real :: z3q
  real :: x3q
  real :: y3q
  real :: xq
  real :: yq
  real :: zq
  real :: t


  ! Evaluate sine and cosine of latitude and longitude of pole point p and
  ! input point q.

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  sinqlat = sin(qlat * pio180)
  cosqlat = cos(qlat * pio180)
  sinqlon = sin(qlon * pio180)
  cosqlon = cos(qlon * pio180)

  ! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
  ! the z axis is the north pole, the x axis is the equator and prime
  ! meridian, and the y axis is the equator and 90 E.

  ! For the pole point, these are:

  x3p = erad * cosplat * cosplon
  y3p = erad * cosplat * sinplon
  z3p = erad * sinplat

  ! For the given lat,lon point, these are:

  z3q = erad * sinqlat
  x3q = erad * cosqlat * cosqlon
  y3q = erad * cosqlat * sinqlon

  ! Transform q point from (x3,y3,z3) coordinates in the above system to
  ! polar stereographic coordinates (x,y,z):

  xq = - sinplon * (x3q-x3p) + cosplon * (y3q-y3p)
  yq =   cosplat * (z3q-z3p)  &
       - sinplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )
  zq =   sinplat * (z3q-z3p)  &
       + cosplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )

  ! Parametric equation for line from antipodal point at (0,0,-2 erad) to
  ! point q has the following parameter (t) value on the polar stereographic
  ! plane:

  t = erad2 / (erad2 + zq)

  ! This gives the following x and y coordinates for the projection of point q
  ! onto the polar stereographic plane:

  x = xq * t
  y = yq * t
end subroutine ed_ll_xy
!==========================================================================================!
!==========================================================================================!
