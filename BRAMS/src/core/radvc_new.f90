subroutine advect_vec(iMax, ibDim,   &  
     ibStart, ibEnd, ibEndU, ibEndV,     &
     kDim, kStart, kEnd, kEndW,          &
     uc,vc,wc,ut,vt,wt,dn0,dn0u,dn0v,    &
     dxt,dxu,dxv,dyt,dyu,dyv,            &
     rtgt,rtgu,rtgv,f13t,f23t,           &
     fmapt,fmapu,fmapv,fmapui,fmapvi,    &
     itopo, hw4, dzt, dzm,               &
     realData, noWorkU, noWorkV,         &
     ! For DEBUG if necessary
     iPerIb, jPerIb)
     !

  ! FIELDS DATA STRUCTURE
  !
  ! 3D fields (i,j,k) are represented by 3D arrays 
  ! addressed by (ib,k,jb). First dimension is a set 
  ! of horizontal points. Second dimension is vertical levels.
  ! Third dimension is how many first dimensions are 
  ! required to represent an entire horizontal plane.
  !
  ! FIRST DIMENSION:
  ! 
  ! Mapping of (i,j) into (ib,jb) is such that, for
  ! every (i,j) mapped into (ib,jb) for a fixed jb, then:
  ! (i+1,j) is mapped into (ib+1,jb)
  ! (i-1,j) is mapped into (ib-1,jb)
  ! (i,j+1) is mapped into (ij+iMax,jb)
  ! (i,j-1) is mapped into (ij-iMax,jb)
  ! for (i,j) declared (iMax,jMax).
  !
  ! Given (i,j) declared as (iMax,jMax), with a ghost
  ! zone of lenght 1 at each direction, meaning that
  ! real data is stored at (2:iMax-1,2:jMax-1),
  ! assume that the section of (i,j) to be mapped into a single jb
  ! takes all data from (iStart,jStart) to (iEnd,jEnd)
  ! in array element order, ghost zone requires to store
  ! all data from (iStart-1,jStart-1) to (iEnd+1,jEnd+1)
  ! in array element order.
  !
  ! Since array element order of (i,j) for an array
  ! declared (iMax,jMax) is (j-1)*iMax+i, then the mapping
  ! (i,j) -> (j-1)*iMax+i
  ! satisfy the above mapping properties for
  ! all pairs (i,j) from (iStart,jStart) into (iEnd,jEnd)
  ! in array element order.
  !
  !
  ! 2D array will have size  ijMax>=iMax*jMax
  ! First data point will be ijStart>=iMax+2
  ! Last data point will be  ijEnd<=ijMax-iMax-1
  ! Ghost zone starts at     ijGhostStart=ijStart-iMax-1
  ! Ghost zone ends at       ijGhostEnd=ijEnd+iMax+1
  !
  ! Take from an array dimensioned (5,6) the section (3,3) to (2,5),
  ! including the ghost zone:
  !     123456
  !   1 --oxxo
  !   2 -ooxxo
  !   3 -oxxoo
  !   4 -oxxo-
  !   5 -oxxo-
  ! This section is mapped into an array of size 30,
  ! first data point at 13, last data point at 22.
  !
  ! SECOND DIMENSION
  !
  ! Second Dimension is vertical, dimensioned kDim.
  ! Real Data is at verticals kStart to kEnd, with
  ! one vertical of ghost zone. That implies
  ! kStart >= 2 and kEnd <= kDim-1
  !
  ! GHOST ZONE IDENTIFICATION.
  !
  ! The ij ghost zone must be identified. This is done by
  ! the RealData array, a 2D real array dimensioned
  ! (ijMax)
  ! with entries 0 at the ghost zone and 1 at real points.
  !
  !
  implicit none

  ! Arguments:
  integer, intent(in   ) :: iMax
  integer, intent(in   ) :: ibDim
  integer, intent(in   ) :: ibStart
  integer, intent(in   ) :: ibEnd
  integer, intent(in   ) :: ibEndU
  integer, intent(in   ) :: ibEndV
  integer, intent(in   ) :: kDim
  integer, intent(in   ) :: kStart
  integer, intent(in   ) :: kEnd
  integer, intent(in   ) :: kEndW
  real,    intent(in   ) :: uc(ibDim,kDim)
  real,    intent(in   ) :: vc(ibDim,kDim)
  real,    intent(in   ) :: wc(ibDim,kDim)
  real,    intent(inout) :: ut(ibDim,kDim)
  real,    intent(inout) :: vt(ibDim,kDim)
  real,    intent(inout) :: wt(ibDim,kDim)
  real,    intent(in   ) :: dn0(ibDim,kDim)
  real,    intent(in   ) :: dn0u(ibDim,kDim)
  real,    intent(in   ) :: dn0v(ibDim,kDim)
  real,    intent(in   ) :: dxt(ibDim)
  real,    intent(in   ) :: dxu(ibDim)
  real,    intent(in   ) :: dxv(ibDim)
  real,    intent(in   ) :: dyt(ibDim)
  real,    intent(in   ) :: dyu(ibDim)
  real,    intent(in   ) :: dyv(ibDim)
  real,    intent(in   ) :: rtgt(ibDim)
  real,    intent(in   ) :: rtgu(ibDim)
  real,    intent(in   ) :: rtgv(ibDim)
  real,    intent(in   ) :: f13t(ibDim)
  real,    intent(in   ) :: f23t(ibDim)
  real,    intent(in   ) :: fmapt(ibDim)
  real,    intent(in   ) :: fmapu(ibDim)
  real,    intent(in   ) :: fmapv(ibDim)
  real,    intent(in   ) :: fmapui(ibDim)
  real,    intent(in   ) :: fmapvi(ibDim)
  integer, intent(in   ) :: itopo
  real,    intent(in   ) :: hw4(kDim)
  real,    intent(in   ) :: dzt(kDim)
  real,    intent(in   ) :: dzm(kDim)
  real,    intent(in   ) :: realData(ibDim)
  real,    intent(in   ) :: noWorkU(ibDim)
  real,    intent(in   ) :: noWorkV(ibDim)
  ! For DEBUG if necessary
  integer, intent(in   ) :: iPerIB(ibDim)
  integer, intent(in   ) :: jPerIB(ibDim)
  !

  ! Local Variables:
  real :: c1x_u(ibDim)
  real :: c1y_u(ibDim)
  real :: c1z_u(ibDim)
  real :: c1x_v(ibDim)
  real :: c1y_v(ibDim)
  real :: c1z_v(ibDim)
  real :: c1x_w(ibDim)
  real :: c1y_w(ibDim)
  real :: c1z_w(ibDim)
  real :: flxu(ibDim,4)
  real :: flxv(ibDim,4)
  real :: flxw(ibDim,4)
  integer :: i
  integer :: ib
  integer :: k
  integer :: ibEndFlx
  integer :: kAuxi
  integer :: kPrev
  integer :: kThis
  integer :: kNext
  integer :: kFutu
  real    :: dninv
  real    :: vel2
  character(len=10) :: c0, c1, c2
  character(len=*), parameter :: h="**(vel_advect_2D)**"
  logical, parameter :: dumpLocal=.FALSE.

  ibEndFlx = max(ibEndU+1, ibEndV+iMax, ibEnd)

  if (dumpLocal) then
     write(c0,"(i10)") ibDim
     write(c1,"(i10)") kDim
     write(c2,"(i10)") iMax
!!$     write(c3,"(i2.2)") mynum
     write(*,"(a)") &
          h//&
!!$          h//trim(c3)//&
          "** fields dimensioned ("//trim(adjustl(c0))//&
          ","//trim(adjustl(c1))//") with iMax="//trim(adjustl(c2))
     write(c0,"(i10)") ibStart
     write(c1,"(i10)") ibEnd
!!$     write(c4,"(i2.2,1x,i2.2)") iPerIb(ibStart),jPerIb(ibStart)
!!$     write(c5,"(i2.2,1x,i2.2)") iPerIb(ibEnd),jPerIb(ibEnd)
!!$     write(*,"(a)") h//trim(c3)//" ibStart="//trim(adjustl(c0))//&
!!$          " i,j="//trim(c4)//"; ibEnd="//trim(adjustl(c1))//&
!!$          " i,j="//trim(c5)
     write(c0,"(i10)") ibEndU 
!!$     write(c4,"(i2.2,1x,i2.2)") iPerIb(ibEndU),jPerIb(ibEndV)
     write(c1,"(i10)") ibEndV
     write(c2,"(i10)") ibEndFlx
!!$     write(*,"(a)") h//trim(c3)//" ibEndU="//trim(adjustl(c0))//&
!!$          " i,j="//trim(c4)//"; ibEndV="//trim(adjustl(c1))//&
!!$          "; ibEndFlx="//trim(adjustl(c2))
     write(c0,"(i10)") KStart
     write(c1,"(i10)") KEnd
     write(c2,"(i10)") KEndW
     write(*,"(a)") &
          h//" kStart="//trim(adjustl(c0))//&
!!$          h//trim(c3)//" kStart="//trim(adjustl(c0))//&
          "; kEnd="//trim(adjustl(c1))//"; kEndW="//trim(adjustl(c2))
!!$     write(*,"(a)") h//" RealData"
!!$     do ib = 1, ibDim, 10
!!$        write(*,"(a,10f4.1)") h, (RealData(i), i=ib, min(ib+9,ibDim))
!!$     end do
!!$     write(*,"(a)") h//" noWorkV"
!!$     do ib = 1, ibDim, 10
!!$        write(*,"(a,10f4.1)") h, (noWorkV(i), i=ib, min(ib+9,ibDim))
!!$     end do
!!$     write(*,"(a)") h//" noWorkU"
!!$     do ib = 1, ibDim, 10
!!$        write(*,"(a,10f4.1)") h, (noWorkU(i), i=ib, min(ib+9,ibDim))
!!$     end do
  end if

  ! vertically independent constants for ut
  ! include noWorkU to keep old values of ut at border points

  do ib = ibStart, ibEndU
     c1z_u(ib) = noWorkU(ib)*.25 / rtgu(ib)
     c1x_u(ib) = noWorkU(ib)*c1z_u(ib) * fmapu(ib) * dxu(ib)
     c1y_u(ib) = noWorkU(ib)*c1z_u(ib) * fmapu(ib) * dyu(ib)
  end do

  ! vertically independent constants for vt
  ! include noWorkV to keep old values of vt at border points

  do ib = ibStart, ibEndV
     c1z_v(ib) = noWorkV(ib)*.25 / rtgv(ib)
     c1x_v(ib) = noWorkV(ib)*c1z_v(ib) * fmapv(ib) * dxv(ib)
     c1y_v(ib) = noWorkV(ib)*c1z_v(ib) * fmapv(ib) * dyv(ib)
  end do

  ! vertically independent constants for wt
  ! include realData to keep old values of wt at border points

  do ib = ibStart, ibEnd
     c1z_w(ib) = realData(ib)*.5 / rtgt(ib)
     c1x_w(ib) = realData(ib)*c1z_w(ib) * fmapt(ib) * dxt(ib)
     c1y_w(ib) = realData(ib)*c1z_w(ib) * fmapt(ib) * dyt(ib)
  end do

  kPrev = 1
  kThis = 2
  kNext = 3
  kFutu = 4

  ! Momentum fluxes flxu and flxv
  !    Compute at k since flxu and flxv at (k) and (k+1)
  !    is required for flxw at (k)

  do ib = 1, ibEndFlx
     flxu(ib,kThis) = uc(ib,kStart-1)*dn0u(ib,kStart-1)*rtgu(ib)*fmapui(ib)
     flxu(ib,kNext) = uc(ib,kStart  )*dn0u(ib,kStart  )*rtgu(ib)*fmapui(ib)
     flxv(ib,kThis) = vc(ib,kStart-1)*dn0v(ib,kStart-1)*rtgv(ib)*fmapvi(ib)
     flxv(ib,kNext) = vc(ib,kStart  )*dn0v(ib,kStart  )*rtgv(ib)*fmapvi(ib)
  enddo

  ! Momentum fluxes flxw
  !    Compute at kStart-1 since flxw(k-1) 
  !    is required for ut, vt, wt at (k)

  if(itopo == 0) then
     do ib = ibStart, ibEndFlx
        flxw(ib,kThis) = wc(ib,kStart-1)*.5*(dn0(ib,kStart-1)+dn0(ib,kStart))
     end do
  else
     do ib = ibStart, ibEndFlx
        flxw(ib,kThis) = wc(ib,kStart-1)  &
             * .5 * (dn0(ib,kStart-1) + dn0(ib,kStart))  &
             + hw4(kStart-1) * ((flxu(ib,kThis) + flxu(ib,kNext)  &
             + flxu(ib-1,kThis) + flxu(ib-1,kNext)) * f13t(ib)  &
             + (flxv(ib,kThis) + flxv(ib,kNext)  &
             + flxv(ib-iMax,kThis) + flxv(ib-iMax,kNext)) * f23t(ib))
     end do
  end if

  ! Main loop for K
  do k = kStart-1, kEnd


     if (k <= kEnd-1) then

        ! Momentum fluxes flxu and flxv
        !    compute flxu(ib,k), for 
        !            ib=ibStart-1, max(ibEndU+1, ibEndV+iMax, ibEnd)
        !            k =kStart-1, max(kEnd,kEndW+1)+1
        !    compute flxv(ib,k), for 
        !            ib=ibStart-iMax, max(ibEndU+1, ibEndV+iMax, ibEnd)
        !            k =kStart-1, max(kEnd,kEndW+1)+1
        !    for a single loop, compute both fluxes at
        !            ib=ibStart-iMax, max(ibEndU+1, ibEndV+iMax, ibEnd)
        !            k =kStart-1, max(kEnd,kEndW+1)+1
        !
        !    Compute at (k+2) since already computed at (k) and (k+1) and
        !    flxu and flxv are required at (k+1) and (k+2) to compute
        !    flxw at (k+1)

        do ib = 1, ibEndFlx
           flxu(ib,kFutu) = uc(ib,k+2)*dn0u(ib,k+2)*rtgu(ib)*fmapui(ib)
           flxv(ib,kFutu) = vc(ib,k+2)*dn0v(ib,k+2)*rtgv(ib)*fmapvi(ib)
        end do

        ! Momentum flux flxw
        !    compute flxw(ib,k), for 
        !            ib=ibStart , max(ibEndU+1, ibEndV+iMax, ibEnd)
        !            k =kStart-1, max(kEnd,kEndW+1)
        !    requires flxu(ib-   1:ib,    k:k+1   )
        !    requires flxv(ib-iMax:ib,    k:k+1    )
        !
        !    Compute at (k+1) since already computed at (k) and
        !    flxw(k+1) is required to compute wt(k)

        if(itopo == 0) then
           do ib = ibStart, ibEndFlx
              flxw(ib,kNext) = wc(ib,k+1)*.5*(dn0(ib,k+1)+dn0(ib,k+2))
           end do
        else
           do ib = ibStart, ibEndFlx
              flxw(ib,kNext) = wc(ib,k+1)  &
                   * .5 * (dn0(ib,k+1) + dn0(ib,k+2))  &
                   + hw4(k+1) * ((flxu(ib,kNext) + flxu(ib,kFutu)  &
                   + flxu(ib-1,kNext) + flxu(ib-1,kFutu)) * f13t(ib)  &
                   + (flxv(ib,kNext) + flxv(ib,kFutu)  &
                   + flxv(ib-iMax,kNext) + flxv(ib-iMax,kFutu)) * f23t(ib))
           end do
        end if
     end if

     if (k>=kStart) then

        ! Advection contribution to U tendency
        !    compute ut(ib,k), for ib=ibStart:ibEndU, k=kStart:kEnd
        !    requires flxu(ib-   1:ib+   1,    k   )
        !    requires flxv(ib-iMax:ib+   1,    k    )
        !    requires flxw(ib     :ib+   1,  k-1:k  )

        do ib = ibStart, ibEndU
!!$           dninv = 1.0/dn0u(ib,k)
!!$           vel2 = 2.*uc(ib,k)
!!$           ut(ib,k) = ut(ib,k) + c1x_u(ib) *dninv * (  &
!!$                (flxu(ib,kThis) + flxu(ib-1,kThis))  &
!!$                * (uc(ib,k) + uc(ib-1,k))  &
!!$                - (flxu(ib,kThis) + flxu(ib+1,kThis))  &
!!$                * (uc(ib,k) + uc(ib+1,k))  &
!!$                + (flxu(ib+1,kThis) - flxu(ib-1,kThis)) * vel2)
!!$           ut(ib,k) = ut(ib,k) + c1y_u(ib) * dninv * (  &
!!$                (flxv(ib-iMax,kThis) + flxv(ib+1-iMax,kThis))  &
!!$                * (uc(ib,k) + uc(ib-iMax,k))  &
!!$                - (flxv(ib,kThis) + flxv(ib+1,kThis))  &
!!$                * (uc(ib,k) + uc(ib+iMax,k))&
!!$                + (flxv(ib,kThis) + flxv(ib+1,kThis) - flxv(ib-iMax,kThis)  &
!!$                - flxv(ib+1-iMax,kThis)) * vel2 )
!!$           ut(ib,k) = ut(ib,k) + c1z_u(ib) * dzt(k) *dninv * (  &
!!$                (flxw(ib,kPrev) + flxw(ib+1,kPrev))  &
!!$                * (uc(ib,k) + uc(ib,k-1))  &
!!$                - (flxw(ib,kThis) + flxw(ib+1,kThis))  &
!!$                * (uc(ib,k) + uc(ib,k+1))   &
!!$                + (flxw(ib,kThis) + flxw(ib+1,kThis) - flxw(ib,kPrev)  &
!!$                - flxw(ib+1,kPrev)) * vel2 )
           ! sem modificacao na aritmetica
           ut(ib,k) = ut(ib,k) + c1x_u(ib)/dn0u(ib,k) * (  &
                (flxu(ib,kThis) + flxu(ib-1,kThis))  &
                * (uc(ib,k) + uc(ib-1,k))  &
                - (flxu(ib,kThis) + flxu(ib+1,kThis))  &
                * (uc(ib,k) + uc(ib+1,k))  &
                + (flxu(ib+1,kThis) - flxu(ib-1,kThis)) * 2.*uc(ib,k))
           ut(ib,k) = ut(ib,k) + c1y_u(ib)/dn0u(ib,k) * (  &
                (flxv(ib-iMax,kThis) + flxv(ib+1-iMax,kThis))  &
                * (uc(ib,k) + uc(ib-iMax,k))  &
                - (flxv(ib,kThis) + flxv(ib+1,kThis))  &
                * (uc(ib,k) + uc(ib+iMax,k))&
                + (flxv(ib,kThis) + flxv(ib+1,kThis) - flxv(ib-iMax,kThis)  &
                - flxv(ib+1-iMax,kThis)) * 2.*uc(ib,k))
           ut(ib,k) = ut(ib,k) + c1z_u(ib) * dzt(k)/dn0u(ib,k) * (  &
                (flxw(ib,kPrev) + flxw(ib+1,kPrev))  &
                * (uc(ib,k) + uc(ib,k-1))  &
                - (flxw(ib,kThis) + flxw(ib+1,kThis))  &
                * (uc(ib,k) + uc(ib,k+1))   &
                + (flxw(ib,kThis) + flxw(ib+1,kThis) - flxw(ib,kPrev)  &
                - flxw(ib+1,kPrev)) * 2.*uc(ib,k))
        end do


        ! Advection contribution to V tendency
        !    compute vt(ib,k), for ib=ibStart:ibEndV, k=kStart:kEnd
        !    requires flxu(ib-   1:ib+iMax,    k   )
        !    requires flxv(ib-iMax:ib+iMax,    k   )
        !    requires flxw(ib     :ib+iMax, k-1:k  )

        do ib = ibStart, ibEndV
           dninv = 1.0/dn0v(ib,k)
           vel2 = 2.*vc(ib,k)
!!$           vt(ib,k) = vt(ib,k) + c1x_v(ib) *dninv * (  &
!!$                (flxu(ib-1,kThis) + flxu(ib-1+iMax,kThis))  &
!!$                * (vc(ib,k) + vc(ib-1,k))  &
!!$                - (flxu(ib,kThis) + flxu(ib+iMax,kThis))  &
!!$                * (vc(ib,k) + vc(ib+1,k))  &
!!$                + (flxu(ib,kThis) + flxu(ib+iMax,kThis) - flxu(ib-1,kThis)  &
!!$                - flxu(ib-1+iMax,kThis)) * vel2 )
!!$           vt(ib,k) = vt(ib,k) + c1y_v(ib) *dninv * (  &
!!$                (flxv(ib,kThis) + flxv(ib-iMax,kThis))  &
!!$                * (vc(ib,k) + vc(ib-iMax,k))  &
!!$                - (flxv(ib,kThis) + flxv(ib+iMax,kThis))  &
!!$                * (vc(ib,k) + vc(ib+iMax,k))  &
!!$                + (flxv(ib+iMax,kThis) - flxv(ib-iMax,kThis))  &
!!$                * vel2 )
!!$           vt(ib,k) = vt(ib,k) + c1z_v(ib) * dzt(k) *dninv * (  &
!!$                (flxw(ib,kPrev) + flxw(ib+iMax,kPrev))  &
!!$                * (vc(ib,k) + vc(ib,k-1))  &
!!$                - (flxw(ib,kThis) + flxw(ib+iMax,kThis))  &
!!$                * (vc(ib,k) + vc(ib,k+1))  &
!!$                + (flxw(ib,kThis) + flxw(ib+iMax,kThis) - flxw(ib,kPrev)  &
!!$                - flxw(ib+iMax,kPrev)) * vel2 )
           ! sem modificacao na aritmetica
           vt(ib,k) = vt(ib,k) + c1x_v(ib)/dn0v(ib,k) * (  &
                (flxu(ib-1,kThis) + flxu(ib-1+iMax,kThis))  &
                * (vc(ib,k) + vc(ib-1,k))  &
                - (flxu(ib,kThis) + flxu(ib+iMax,kThis))  &
                * (vc(ib,k) + vc(ib+1,k))  &
                + (flxu(ib,kThis) + flxu(ib+iMax,kThis) - flxu(ib-1,kThis)  &
                - flxu(ib-1+iMax,kThis)) * 2.*vc(ib,k) )
           vt(ib,k) = vt(ib,k) + c1y_v(ib)/dn0v(ib,k) * (  &
                (flxv(ib,kThis) + flxv(ib-iMax,kThis))  &
                * (vc(ib,k) + vc(ib-iMax,k))  &
                - (flxv(ib,kThis) + flxv(ib+iMax,kThis))  &
                * (vc(ib,k) + vc(ib+iMax,k))  &
                + (flxv(ib+iMax,kThis) - flxv(ib-iMax,kThis))  &
                * 2.*vc(ib,k) )
           vt(ib,k) = vt(ib,k) + c1z_v(ib) * dzt(k)/dn0v(ib,k) * (  &
                (flxw(ib,kPrev) + flxw(ib+iMax,kPrev))  &
                * (vc(ib,k) + vc(ib,k-1))  &
                - (flxw(ib,kThis) + flxw(ib+iMax,kThis))  &
                * (vc(ib,k) + vc(ib,k+1))  &
                + (flxw(ib,kThis) + flxw(ib+iMax,kThis) - flxw(ib,kPrev)  &
                - flxw(ib+iMax,kPrev)) * 2.*vc(ib,k) )
        end do


        if (k <= kEndW) then

           ! Advection contribution to W tendency
           !    compute wt(ib,k), for ib=ibStart:ibEnd, k=kStart:kEndW
           !    requires flxu(ib-   1:ib     , k  :k+1)
           !    requires flxv(ib-iMax:ib     , k  :k+1)
           !    requires flxw(     ib:ib     , k-1:k+1)

           do ib = ibStart, ibEnd
!!$              dninv = 1.0/(dn0(ib,k)+dn0(ib,k+1))
!!$              vel2 = 2.*wc(ib,k)
!!$              wt(ib,k) = wt(ib,k)  &
!!$                   + c1x_w(ib) * dninv * (  &
!!$                   (flxu(ib-1,kThis) + flxu(ib-1,kNext))  &
!!$                   * (wc(ib,k) + wc(ib-1,k))  &
!!$                   - (flxu(ib,kThis) + flxu(ib,kNext))  &
!!$                   * (wc(ib,k) + wc(ib+1,k))  &
!!$                   + (flxu(ib,kThis) + flxu(ib,kNext) - flxu(ib-1,kThis)  &
!!$                   - flxu(ib-1,kNext)) * vel2 )
!!$              wt(ib,k) = wt(ib,k)  &
!!$                   + c1y_w(ib) * dninv * (  &
!!$                   (flxv(ib-iMax,kThis) + flxv(ib-iMax,kNext))  &
!!$                   * (wc(ib,k) + wc(ib-iMax,k))  &
!!$                   - (flxv(ib,kThis) + flxv(ib,kNext))  &
!!$                   * (wc(ib,k) + wc(ib+iMax,k))  &
!!$                   + (flxv(ib,kThis) + flxv(ib,kNext) - flxv(ib-iMax,kThis)  &
!!$                   - flxv(ib-iMax,kNext)) * vel2 )
!!$              wt(ib,k) = wt(ib,k)  &
!!$                   + c1z_w(ib) * dzm(k) * dninv * (  &
!!$                   (flxw(ib,kThis) + flxw(ib,kPrev))  &
!!$                   * (wc(ib,k) + wc(ib,k-1))  &
!!$                   - (flxw(ib,kThis) + flxw(ib,kNext))  &
!!$                   * (wc(ib,k) + wc(ib,k+1))   &
!!$                   + (flxw(ib,kNext) - flxw(ib,kPrev)) * vel2 )
              ! sem modificacao na aritmetica
              wt(ib,k) = wt(ib,k)  &
                   + c1x_w(ib)/(dn0(ib,k)+dn0(ib,k+1)) * (  &
                   (flxu(ib-1,kThis) + flxu(ib-1,kNext))  &
                   * (wc(ib,k) + wc(ib-1,k))  &
                   - (flxu(ib,kThis) + flxu(ib,kNext))  &
                   * (wc(ib,k) + wc(ib+1,k))  &
                   + (flxu(ib,kThis) + flxu(ib,kNext) - flxu(ib-1,kThis)  &
                   - flxu(ib-1,kNext)) * 2.*wc(ib,k) )
              wt(ib,k) = wt(ib,k)  &
                   + c1y_w(ib)/(dn0(ib,k)+dn0(ib,k+1)) * (  &
                   (flxv(ib-iMax,kThis) + flxv(ib-iMax,kNext))  &
                   * (wc(ib,k) + wc(ib-iMax,k))  &
                   - (flxv(ib,kThis) + flxv(ib,kNext))  &
                   * (wc(ib,k) + wc(ib+iMax,k))  &
                   + (flxv(ib,kThis) + flxv(ib,kNext) - flxv(ib-iMax,kThis)  &
                   - flxv(ib-iMax,kNext)) * 2.*wc(ib,k) )
              wt(ib,k) = wt(ib,k)  &
                   + c1z_w(ib) * dzm(k)/(dn0(ib,k)+dn0(ib,k+1)) * (  &
                   (flxw(ib,kThis) + flxw(ib,kPrev))  &
                   * (wc(ib,k) + wc(ib,k-1))  &
                   - (flxw(ib,kThis) + flxw(ib,kNext))  &
                   * (wc(ib,k) + wc(ib,k+1))   &
                   + (flxw(ib,kNext) - flxw(ib,kPrev)) * 2.*wc(ib,k) )
           end do
        end if
     end if
     kAuxi = kPrev
     kPrev = kThis
     kThis = kNext
     kNext = kFutu
     kFutu = kAuxi
  end do ! K Main Loop
end subroutine advect_vec





subroutine advect_sca(ibDim, kDim, iDim, &
     ibStart, ibEnd, kStart, kEnd, nScalar, &
     dtlt, rtgt, rtgu, rtgv, dxu, dyv, dzm, zt, zm, hw4, &
     dn0, dn0u, dn0v, f13t, f23t, dxt, dyt, dzt, &
     fmapt, fmapui, fmapvi, &
     up, vp, wp, uc, vc, wc, RealData, noGhostI, &
     noBottonGhostI, scalarp, scalart,         &
     ! For DEBUG if necessary
     iPerIb, jPerIb)
     !

  implicit none

  ! Arguments:
  integer, intent(in) :: ibDim
  integer, intent(in) :: kDim
  integer, intent(in) :: iDim
  integer, intent(in) :: ibStart
  integer, intent(in) :: ibEnd
  integer, intent(in) :: kStart
  integer, intent(in) :: kEnd
  integer, intent(in) :: nScalar
  real, intent(in) :: dtlt
  real, intent(in) :: rtgt(ibDim)
  real, intent(in) :: rtgu(ibDim)
  real, intent(in) :: rtgv(ibDim)
  real, intent(in) :: dxu(ibDim)
  real, intent(in) :: dyv(ibDim)
  real, intent(in) :: dzm(kDim)
  real, intent(in) :: zt(kDim)
  real, intent(in) :: zm(kDim)
  real, intent(in) :: hw4(kDim)
  real, intent(in) :: fmapt(ibDim)
  real, intent(in) :: fmapui(ibDim)
  real, intent(in) :: fmapvi(ibDim)
  real, intent(in) :: dn0(ibDim,kDim)
  real, intent(in) :: dn0u(ibDim,kDim)
  real, intent(in) :: dn0v(ibDim,kDim)
  real, intent(in) :: f13t(ibDim)
  real, intent(in) :: f23t(ibDim)
  real, intent(in) :: dxt(ibDim)
  real, intent(in) :: dyt(ibDim)
  real, intent(in) :: dzt(kDim)
  real, intent(in) :: up(ibDim,kDim)
  real, intent(in) :: vp(ibDim,kDim)
  real, intent(in) :: wp(ibDim,kDim)
  real, intent(in) :: uc(ibDim,kDim)
  real, intent(in) :: vc(ibDim,kDim)
  real, intent(in) :: wc(ibDim,kDim)
  real, intent(in) :: RealData(ibDim)
  real, intent(in) :: noGhostI(ibDim)
  real, intent(in) :: noBottonGhostI(ibDim)
  real, intent(in) :: scalarp(ibDim,kDim,nScalar)
  real, intent(inout) :: scalart(ibDim,kDim,nScalar)
  ! For DEBUG if necessary
  integer, intent(in   ) :: iPerIB(ibDim)
  integer, intent(in   ) :: jPerIB(ibDim)
  !
  
  ! Local Variables:
  real :: dtlto2
  real :: dtli(ibDim)
  real :: dfact
  real :: vctr3
  real :: rtgti(ibDim)
  integer :: ib, ibGhostStart, ibGhostEnd
  integer :: k
  integer :: nSca
  real :: vt3da(ibDim,kdim)
  real :: vt3db(ibDim,kdim)
  real :: vt3dc(ibDim,kdim)
  real :: vt3dd(ibDim,kdim)
  real :: vt3de(ibDim,kdim)
  real :: vt3df(ibDim,kdim)
  real :: vt3dh(ibDim,kdim)
  real :: vt3di(ibDim,kdim)
  real :: vt3dj(ibDim,kdim)
  real :: vt3dk(ibDim,kdim)
  real :: vt3dl(ibDim,kdim)
  real :: scr1(ibDim,2)
  real :: vt3dg_1D(ibDim)
  real :: vt3dg_2D(ibDim,2)
  real :: c1(ibDim)
  real :: c2(ibDim)
  real :: c3(ibDim)
  real :: c4(ibDim)
  real :: vctr1(kDim)
  real :: vctr2(kDim)
  real, parameter :: jdim_original=1
  integer :: kPrev, kThis, kAux

  ! Advect  scalars

  dtlto2 = .5 * dtlt
  dfact = .5

  ibGhostStart = ibStart-iDim-1
  ibGhostEnd   = ibEnd+iDim+1

  do ib = ibStart, ibEnd
     dtli(ib) = RealData(ib)/dtlt
  end do

  do ib = ibGhostStart, ibGhostEnd
     rtgti(ib) = 1.0 / rtgt(ib)
     c1(ib)    = .5 * dxu(ib)
     c2(ib)    = .5 * dyv(ib)
     c3(ib)    = dxt(ib) * fmapt(ib) * rtgti(ib)
     c4(ib)    = dyt(ib) * fmapt(ib) * rtgti(ib)
  end do

  do ib = ibGhostStart, ibGhostEnd
     vt3di(ib,1) = noBottonGhostI(ib) * .5
     vt3di(ib,2) = noBottonGhostI(ib) * .5
     vt3di(ib,3) = noGhostI(ib)*.5
     vt3di(ib,4) = noGhostI(ib)*.5
  end do

  do k = kStart-1, kEnd
     vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
     vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
  end do

  do k = kStart-1, kEnd+1
     do ib = ibGhostStart, ibGhostEnd
        vt3da(ib,k) = (up(ib,k)+uc(ib,k)) * dtlto2
        vt3db(ib,k) = (vp(ib,k)+vc(ib,k)) * dtlto2
        vt3dc(ib,k) = (wp(ib,k)+wc(ib,k)) * dtlto2
     end do
  end do

  do k = kStart-1, kEnd
     do ib = ibStart, ibEnd
        vt3dc(ib,k) = (&
             (vt3da(ib,k) + vt3da(ib,k+1) + vt3da(ib-1,k) + vt3da(ib-1,k+1)) * f13t(ib)  &
             + (vt3db(ib,k) + vt3db(ib,k+1) + vt3db(ib-iDim,k) + vt3db(ib-iDim,k+1)) * f23t(ib)) * hw4(k)  &
             + vt3dc(ib,k) * rtgti(ib)
        vt3df(ib,k) = .5 * vt3dc(ib,k) * dzm(k)
     end do

     do ib = ibGhostStart, ibGhostEnd
        vt3dd(ib,k) = noBottonGhostI(ib) * c1(ib) * vt3da(ib,k)
        vt3de(ib,k) = noGhostI(ib)*c2(ib) * vt3db(ib,k)
        vctr3 = 1.0/dn0(ib,k)
        !vt3dh(ib,k) = noGhostI(ib) * c3(ib) * vctr3
        vt3dh(ib,k) = c3(ib) * vctr3

        vt3dj(ib,k) = c4(ib) * vctr3
        vt3dl(ib,k) = RealData(ib) * vt3dj(ib,k)
        vt3dk(ib,k) = RealData(ib) * dzt(k) * vctr3
     end do
  end do

  ! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
  ! into mass fluxes times dtlt.

  do ib = ibGhostStart, ibGhostEnd
     c1(ib) = fmapui(ib) * rtgu(ib)
     c2(ib) = fmapvi(ib) * rtgv(ib)
  end do

  do k = kStart-1, kEnd
     do ib = ibGhostStart, ibGhostEnd
        vt3da(ib,k) = vt3da(ib,k) * c1(ib) * dn0u(ib,k)
        vt3db(ib,k) = vt3db(ib,k) * c2(ib) * dn0v(ib,k)
        vt3dc(ib,k) = vt3dc(ib,k) * .5  &
             * (dn0(ib,k) + dn0(ib,k+1))
     end do
  end do

!!$call ftrace_region_end("init")
!!$call ftrace_region_begin("sca")

  do nSca = 1, nScalar

     kThis = 1
     kPrev = 2
     do k = kStart-2, kEnd
        do ib = ibGhostStart, ibGhostEnd

           scr1(ib,kThis) = scalarp(ib,k+1,nSca)

        end do


        ! Compute scalar flux times dtlt [VT3DG]

        if (k >= kStart-1 .and. k <= kEnd-1) then
           do ib = ibGhostStart, ibGhostEnd-1

              vt3dg_1D(ib) = vt3da(ib,k+1)  &
                   * (vt3di(ib,1) * scr1(ib,kThis)  &
                   +  vt3di(ib,2) * scr1(ib+1,kThis)  &
                   +  vt3dd(ib,k+1) * (scr1(ib,kThis) - scr1(ib+1,kThis)))

           end do

           
           ! Modify fluxes to retain positive-definiteness on scalar quantities.
           !    If a flux will remove 1/2 quantity during a timestep,
           !    reduce to first order flux. This will remain positive-definite
           !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
           !    both fluxes are evacuating the box.
           
           do ib = ibGhostStart, ibGhostEnd-1

              if (vt3da(ib,k+1) .gt. 0.) then
                 if (vt3dg_1D(ib) * vt3dh(ib,k+1) .gt.  &
                      dfact * scr1(ib,kThis)) then

                    vt3dg_1D(ib) = vt3da(ib,k+1) * scr1(ib,kThis)

                 end if
              else
                 if (vt3da(ib,k+1) .lt. 0.) then
                    if (-vt3dg_1D(ib) * vt3dh(ib+1,k+1) .gt.  &
                         dfact * scr1(ib+1,kThis)) then

                       vt3dg_1D(ib) = vt3da(ib,k+1) * scr1(ib+1,kThis)

                    end if
                 end if
              end if

           end do
           
           ! Compute flux divergence
           
           do ib = ibGhostStart+1, ibGhostEnd-1

              scr1(ib,kThis) = scr1(ib,kThis)  &
                   + vt3dh(ib,k+1) * (vt3dg_1D(ib-1) - vt3dg_1D(ib)  &
                   + scalarp(ib,k+1,nSca) * (vt3da(ib,k+1) - vt3da(ib-1,k+1)))

           end do
           
           if (jdim_original == 1)  then
              
              ! y component of flux divergence
              !    compute  vt3dg(ib,k), for ib=ibStart:ibEnd, k=kStart:kEnd
              !    requires vt3db(ib,k)
              !    requires scr1 (ib:ib+iDim,k)
              !    requires vt3de(ib,k)
              !    requires vt3di(ib,3:4)
              !    Avoid changing Ghost Zone points by inserting noGhostI
              !    in vt3di(:,3:4) and vt3de(:,k)
        
              do ib = ibGhostStart, ibEnd
                 vt3dg_1D(ib) = vt3db(ib,k+1)  &
                      * (vt3di(ib,3) * scr1(ib,kThis)  &
                      + vt3di(ib,4) * scr1(ib+iDim,kThis)  &
                      + vt3de(ib,k+1) * (scr1(ib,kThis) - scr1(ib+iDim,kThis)))
              end do
              
              ! Modify fluxes to retain positive-definiteness on scalar 
              ! quantities.
              ! If a flux will remove 1/2 quantity during a timestep,
              ! reduce to first order flux. This will remain positive-definite
              ! under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
              ! both fluxes are evacuating the box.
              
              do ib = ibGhostStart, ibEnd
                 if (vt3db(ib,k+1) .gt. 0.) then
                    if (vt3dg_1D(ib) * vt3dj(ib,k+1) .gt.  &
                         dfact * scr1(ib,kThis)) then
                       vt3dg_1D(ib) = vt3db(ib,k+1) * scr1(ib,kThis)
                    end if
                 else
                    if (vt3db(ib,k+1) .lt. 0.) then
                       if (-vt3dg_1D(ib) * vt3dj(ib+iDim,k+1) .gt.  &
                            dfact * scr1(ib+iDim,kThis)) then
                          vt3dg_1D(ib) = vt3db(ib,k+1) * scr1(ib+iDim,kThis)
                       end if
                    end if
                 end if
              end do
              
              ! y component of flux divergence
              !    compute  scr1   (ib,k), for ib=ibStart:ibEnd, k=kStart:kEnd
              !    requires scr1   (ib,k)
              !    requires scalarp(ib,k)
              !    requires vt3dl  (ib,k)
              !    requires vt3dg  (ib-iDim:ib)
              !    requires vt3db  (ib-iDim:ib,k)
              !    Results change if inserting RealData in vt3dj;build vt3dl as
              !    a copy of vt3dj with RealData
        
              do ib = ibStart, ibEnd
                 scr1(ib,kThis) = scr1(ib,kThis)  &
                      + vt3dl(ib,k+1) * (vt3dg_1D(ib-iDim) - vt3dg_1D(ib)  &
                      + scalarp(ib,k+1,nSca) * &
                      (vt3db(ib,k+1) - vt3db(ib-iDim,k+1)))
              end do
           end if
        end if

        if (k>=kStart-1 .and. k<=kEnd) then

           ! z component of flux divergence
           !    compute  vt3dg(ib,k), for ib=ibStart:ibEnd, k=kStart-1:kEnd
           !    requires vt3dc(ib,k)
           !    requires scr1 (ib,k:k+1)
           !    requires vt3df(ib,k)
        
           do ib = ibStart, ibEnd
              vt3dg_2D(ib,kThis) = vt3dc(ib,k)  &
                   * RealData(ib)*(vctr1(k) * scr1(ib,kPrev)  &
                   +  vctr2(k) * scr1(ib,kThis)  &
                   +  vt3df(ib,k) * (scr1(ib,kPrev) - scr1(ib,kThis)))
           end do
        
        
           ! Modify fluxes to retain positive-definiteness on scalar quantities
           !    If a flux will remove 1/2 quantity during a timestep,
           !    reduce to first order flux. This will remain positive-definite
           !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
           !    both fluxes are evacuating the box.
           
           do ib = ibStart, ibEnd
              if (vt3dc(ib,k) .gt. 0.) then
                 if (vt3dg_2D(ib,kThis) * vt3dk(ib,k) .gt.  &
                      dfact * scr1(ib,kPrev)) then
                    vt3dg_2D(ib,kThis) = vt3dc(ib,k) * scr1(ib,kPrev)
                 end if
              else
                 if (vt3dc(ib,k) .lt. 0.) then
                    if (-vt3dg_2D(ib,kThis) * vt3dk(ib,k+1) .gt.  &
                         dfact * scr1(ib,kThis)) then
                       vt3dg_2D(ib,kThis) = vt3dc(ib,k) * scr1(ib,kThis)
                    end if
                 end if
              end if
           end do
        end if
        
        if (k >= kStart .and. k<=kEnd) then
        
           ! z component of flux divergence
           !    compute  scr1   (ib,k), for ib=ibStart:ibEnd, k=kStart:kEnd
           !    requires scr1   (ib,k)
           !    requires scalarp(ib,k)
           !    requires vt3dk  (ib,k)
           !    requires vt3dg  (ib,k-1:k)
           !    requires vt3dc  (ib,k-1:k)
           !    Avoid changing Ghost Zone points by inserting RealData
           !    in vt3dk
        
           do ib = ibStart, ibEnd
              scr1(ib,kPrev) = scr1(ib,kPrev)  &
                   + vt3dk(ib,k) * (vt3dg_2D(ib,kPrev) - vt3dg_2D(ib,kThis)  &
                   + scalarp(ib,k,nSca) * (vt3dc(ib,k) - vt3dc(ib,k-1)))
           end do
           
           ! Update tendency with advection
           !    compute scalart(ib,k), for ib=ibStart:ibEnd, k=kStart:kEnd
           !    requires scalart, scr1, scalarp, all (ib,k)
           !    Avoid changing Ghost Zone points by inserting RealData
           !    in dtli

           do ib = ibStart, ibEnd
              scalart(ib,k,nSca) = scalart(ib,k,nSca) + &
                   (scr1(ib,kPrev)-scalarp(ib,k,nSca)) * dtli(ib)
           end do

        end if

        kAux  = kThis
        kThis = kPrev
        kPrev = kAux

     end do
  end do
end subroutine advect_sca
