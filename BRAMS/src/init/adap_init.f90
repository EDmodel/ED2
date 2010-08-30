!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine lpuvw_init(n2,n3,flpu,flpv,flpw)
  implicit none
  integer :: n2,n3
  real, dimension(n2,n3) :: flpu,flpv,flpw

  integer :: i,j

  do j = 1,n3
     do i = 1,n2
        flpu(i,j) = 2.
        flpv(i,j) = 2.
        flpw(i,j) = 2.
     enddo
  enddo
end subroutine lpuvw_init



subroutine ctrlvols (n1,n2,n3,nstbot,dzt,xm,ym,zm,platn,plonn  &
     ,aru,arv,arw,volt,volu,volv,volw  &
     ,flpu,flpv,flpw,dxu,dxv,dxt,dyu,dyv,dyt,topma,topm,ifm, nzg, npatch)

  implicit none

  integer :: n1,n2,n3,nstbot,ifm, nzg, npatch
  real,dimension(n1,n2,n3) :: aru,arv,arw,volt,volu,volv,volw
  real, dimension(n2,n3) :: dxu,dxv,dxt,dyu,dyv,dyt,topma,topm
  real, dimension(*) :: dzt,xm,ym,zm
  real :: platn,plonn
  real, dimension(n2,n3) :: flpu,flpv,flpw

  integer :: i,j,im,jm,iter,k,km,imin,jmin,ip,jp,kp
  real :: hmin,hmax

  ! Define face areas and volumes for all cells
  ! This algorithm applies to regular input topography defined at M points,
  !    and must be replaced or modified for more complex boundary topologies
  ! Small corrections to topography are applied to prevent very small apertures

  do j = 1,n3
     do i = 1,n2
        topm(i,j) = topma(i,j)      
     enddo
  enddo

  do iter = 1,2

     do j = 1,n3
        jm = max(j-1,1)
        do i = 1,n2
           im = max(i-1,1)
           ! aru:
           if (topma(i,j) > topma(i,jm)) then
              hmin = topma(i,jm)
              hmax = topma(i,j)
              jmin = jm
           elseif (topma(i,j) < topma(i,jm)) then
              hmin = topma(i,j)
              hmax = topma(i,jm)
              jmin = j
           else
              hmin = topma(i,j)
              hmax = topma(i,j)
           endif

           do k = 1,n1
              km = max(k-1,1)

              if (zm(k) <= hmin) then
                 aru(k,i,j) = 0.
              elseif (zm(km) >= hmax) then
                 aru(k,i,j) = 1. / (dyu(i,j) * dzt(k))
              elseif(zm(k) <  hmax .and. zm(km) < hmin) then
                 aru(k,i,j) = .5 * (zm(k) - hmin) ** 2  &
                      / (dyu(i,j) * (hmax - hmin))

                 if (aru(k,i,j) < .05 / (dyu(i,j) * dzt(k))) then
                    topm(i,jmin) = zm(k)
                    if (iter == 1) write(6,301)i,jmin,topma(i,jmin),topm(i,jmin)
                    if (iter == 2) write(6,311)i,jmin,topma(i,jmin),topm(i,jmin)
                 endif

              elseif(zm(k) <  hmax .and. zm(km) >=  hmin) then
                 aru(k,i,j) = (.5 * (zm(k) + zm(km)) - hmin)  &
                      / ((hmax - hmin) * dyu(i,j) * dzt(k))
              elseif(zm(k) >= hmax .and. zm(km) < hmin) then
                 aru(k,i,j) = (zm(k) - .5 * (hmax + hmin)) / dyu(i,j)

                 if (aru(k,i,j) < .05 / (dyu(i,j) * dzt(k))) then
                    topm(i,jm) = zm(k)
                    topm(i,j) = zm(k)
                    if (iter == 1) write(6,302)i,jm,topma(i,jm),topm(i,jm)
                    if (iter == 1) write(6,303)i,j,topma(i,j),topm(i,j)
                    if (iter == 2) write(6,312)i,jm,topma(i,jm),topm(i,jm)
                    if (iter == 2) write(6,313)i,j,topma(i,j),topm(i,j)
                 endif

              elseif(zm(k) >= hmax .and. zm(km) >=  hmin) then
                 aru(k,i,j) = 1. / (dyu(i,j) * dzt(k))  &
                      - .5 * (hmax - zm(km)) ** 2 / ((hmax - hmin) * dyu(i,j))
              else
                 print*, 'aru option not reached ',k,i,j,zm(k),zm(km),hmax,hmin
              endif
           enddo
           ! arv:
           if (topma(i,j) > topma(im,j)) then
              hmin = topma(im,j)
              hmax = topma(i,j)
              imin = im
           elseif (topma(i,j) < topma(im,j)) then
              hmin = topma(i,j)
              hmax = topma(im,j)
              imin = i
           else
              hmin = topma(i,j)
              hmax = topma(i,j)
           endif


           do k = 1,n1
              km = max(k-1,1)

              if (zm(k) <= hmin) then
                 arv(k,i,j) = 0.
              elseif (zm(km) >= hmax) then
                 arv(k,i,j) = 1. / (dxv(i,j) * dzt(k))
              elseif(zm(k) <  hmax .and. zm(km) < hmin) then
                 arv(k,i,j) = .5 * (zm(k) - hmin) ** 2  &
                      / (dxv(i,j) * (hmax - hmin))

                 if (arv(k,i,j) < .05 / (dxv(i,j) * dzt(k))) then
                    topm(imin,j) = zm(k)
                    if (iter == 1) write(6,304)imin,j,topma(imin,j),topm(imin,j)
                    if (iter == 2) write(6,314)imin,j,topma(imin,j),topm(imin,j)
                 endif

              elseif(zm(k) <  hmax .and. zm(km) >=  hmin) then
                 arv(k,i,j) = (.5 * (zm(k) + zm(km)) - hmin)  &
                      / ((hmax - hmin) * dxv(i,j) * dzt(k))
              elseif(zm(k) >= hmax .and. zm(km) < hmin) then
                 arv(k,i,j) = (zm(k) - .5 * (hmax + hmin)) / dxv(i,j)

                 if (arv(k,i,j) < .05 / (dxv(i,j) * dzt(k))) then
                    topm(im,j) = zm(k)
                    topm(i,j) = zm(k)
                    if (iter == 1) write(6,305)im,j,topma(im,j),topm(im,j)
                    if (iter == 1) write(6,306)i,j,topma(i,j),topm(i,j)
                    if (iter == 2) write(6,315)im,j,topma(im,j),topm(im,j)
                    if (iter == 2) write(6,316)i,j,topma(i,j),topm(i,j)
                 endif

              elseif(zm(k) >= hmax .and. zm(km) >=  hmin) then
                 arv(k,i,j) = 1. / (dxv(i,j) * dzt(k))  &
                      - .5 * (hmax - zm(km)) ** 2 / ((hmax - hmin) * dxv(i,j))
              else
                 print*, 'arv option not reached ',k,i,j,zm(k),zm(km),hmax,hmin
              endif
           enddo
        enddo
     enddo

     !   do j = 1,n3
     !      do i = 1,n2
     !         if (abs(topma(i,j) - topm(i,j)) > .001) then
     !            if (iter == 1) write(6,303)i,j,topma(i,j),topm(i,j)
     !            if (iter == 2) write(6,304)i,j,topma(i,j),topm(i,j)
     !         endif
     !     enddo
     !   enddo   

301  format('aru1 - first topm adjustment performed ',2i4,2f8.1)
311  format('aru1 - second topm adjustment performed ',2i4,2f8.1)
302  format('aru2 - first topm adjustment performed ',2i4,2f8.1)
312  format('aru2 - second topm adjustment performed ',2i4,2f8.1)
303  format('aru3 - first topm adjustment performed ',2i4,2f8.1)
313  format('aru3 - second topm adjustment performed ',2i4,2f8.1)
304  format('arv1 - first topm adjustment performed ',2i4,2f8.1)
314  format('arv1 - second topm adjustment performed ',2i4,2f8.1)
305  format('arv2 - first topm adjustment performed ',2i4,2f8.1)
315  format('arv2 - second topm adjustment performed ',2i4,2f8.1)
306  format('arv3 - first topm adjustment performed ',2i4,2f8.1)
316  format('arv3 - second topm adjustment performed ',2i4,2f8.1)

     if (iter == 1) then 
        do j = 1,n3
           do i = 1,n2
              topma(i,j) = topm(i,j)
           enddo
        enddo
     endif

  enddo

  do j = 1,n3
     do i = 1,n2
        topm(i,j) = 0.
     enddo
  enddo

  ! arw:

  do j = 1,n3
     jm = max(j-1,1)
     do i = 1,n2
        im = max(i-1,1)

        hmin = min(topma(i,j),topma(im,j),topma(i,jm),topma(im,jm))
        hmax = max(topma(i,j),topma(im,j),topma(i,jm),topma(im,jm))

        do k = 1,n1
           if (zm(k) >= hmax) then
              arw(k,i,j) = 1. / (dxt(i,j) * dyt(i,j))
           elseif (zm(k) <= hmin) then
              arw(k,i,j) = 0.
           else
              arw(k,i,j) = 0.
              if (zm(k) >= topma(i,j)) arw(k,i,j) = arw(k,i,j)  &
                   + .25 / (dxt(i,j) * dyt(i,j))
              if (zm(k) >= topma(im,j)) arw(k,i,j) = arw(k,i,j)  &
                   + .25 / (dxt(i,j) * dyt(i,j))
              if (zm(k) >= topma(i,jm)) arw(k,i,j) = arw(k,i,j)  &
                   + .25 / (dxt(i,j) * dyt(i,j))
              if (zm(k) >= topma(im,jm)) arw(k,i,j) = arw(k,i,j)  &
                   + .25 / (dxt(i,j) * dyt(i,j))
           endif

        enddo
     enddo
  enddo

  ! Special blocked cells:

  if (nstbot == 1) then

     do j = 1,n3
        do i = 1,n2
           aru(1,i,j) = 0.
           arv(1,i,j) = 0.
           arw(1,i,j) = 0.
        enddo
     enddo

  endif

  ! Read and apply building file

  call adap_bldg(n1,n2,n3,xm(1),ym(1),zm(1),platn,plonn  &
       ,aru(1,1,1),arv(1,1,1),arw(1,1,1),topma(1,1),ifm, nzg, npatch)

  ! Define lowest prognostic level for U, V, and W.

  do j = 1,n3
     do i = 1,n2
        do k = n1-1,2,-1
           if (aru(k,i,j) > .001 / (dyu(i,j) * dzt(k)))   flpu(i,j) = real(k)
           if (arv(k,i,j) > .001 / (dxv(i,j) * dzt(k)))   flpv(i,j) = real(k)
           if (arw(k,i,j) > .001 / (dxt(i,j) * dyt(i,j))) flpw(i,j) = real(k)
        enddo
     enddo
  enddo

  ! Define volt

  do j = 1,n3
     jm = max(j-1,1)
     do i = 1,n2
        im = max(i-1,1)
        do k = 1,n1
           volt(k,i,j) = max(1.e-6                          &
                ,aru(k,i,j) * dyu(i,j) * dzt(k)     &
                ,aru(k,im,j) * dyu(im,j) * dzt(k)   &
                ,arv(k,i,j) * dxv(i,j) * dzt(k)     &
                ,arv(k,i,jm) * dxv(i,jm) * dzt(k)   &
                ,arw(k,i,j) * dxt(i,j) * dyt(i,j))  &
                / (dxt(i,j) * dyt(i,j) * dzt(k))
        enddo
     enddo
  enddo

  ! Define volu, volv, and volw

  do j = 1,n3
     jp = min(j+1,n3)
     do i = 1,n2
        ip = min(i+1,n2)

        do k = 1,n1
           kp = min(k+1,n1)
           volu(k,i,j) = .5 * (volt(k,i,j) + volt(k,ip,j))
           volv(k,i,j) = .5 * (volt(k,i,j) + volt(k,i,jp)) 
           volw(k,i,j) = .5 * (volt(k,i,j) + volt(kp,i,j)) 
        enddo
     enddo
  enddo
end subroutine ctrlvols







subroutine adap_bldg(n1,n2,n3,xm,ym,zm,platn,plonn,aru,arv,arw,topma  &
     ,ifm, nzg, npatch)

  use mem_leaf

  ! subroutine adap_bldg 

  ! Purpose - read list of building locations and sizes from file and
  !   modify aru, arv, and arw arrays accordingly.  

  ! Algorithm - Buildings will completely block out one or more grid cells
  !    and will not partially block out any.  Thus, buildings are assumed
  !    rectangular with vertical walls and a flat top, and are oriented 
  !    with the model grid coordinate axes.  Any building, no matter how 
  !    small, will completely block out at least one grid cell, and once
  !    the first grid cell is blocked, any additional cell at least 10%
  !    covered by the building in the x and/or y direction will also be 
  !    blocked.

  ! Input variables

  ! n1,n2,n3 - grid dimensions
  ! xm,ym,zm - grid coordinates (m)
  ! aru,arv,arw - grid cell apertures (m^2)
  ! topma - topography height (m)
  ! platn,plonn - grid pole point latitude and longitude (deg)

  ! Output variables

  ! aru,arv,arw - grid cell apertures (m^2)

  ! Local variables

  ! blat,blon           - latitude, longitude of each building center(deg)
  ! bdx,bdy,bdz         - x,y,z dimensions of each building (m)
  ! bx,by               - x,y location of each building center
  ! bi,bj               - real building coordinates in i,j space
  ! ib,jb               - integer part of bi,bj
  ! nbdi,nbdj           - x,y building dimension in number of grid points
  ! ib1,ib2,jb1,jb2,kb2 - min,max i,j,k index blocked out by a building
  ! topma_t             - average of 4 topma values at t point
  ! topma_avg           - average of topma_t values in 1 building
  ! btop                - absolute elevation of building top

  implicit none

  integer :: n1,n2,n3,ifm, nzg, npatch
  real :: platn,plonn
  real, dimension(n1,n2,n3) :: aru,arv,arw
  real, dimension(n2,n3) :: topma
  real, dimension(*) :: xm,ym,zm

  integer, parameter :: maxbldgs = 1000
  integer :: lb,numbldgs,nbdi,nbdj,ib,jb,ib1,ib2,jb1,jb2,kb2,i,j,k
  real :: bx,by,bi,bj,topma_avg,topma_t,btop
  real, dimension(maxbldgs) :: blat,blon,bdx,bdy,bdz
  integer, dimension(maxbldgs) :: itype

  integer :: marker,iver
  logical :: there


  ! Open and read list of building sizes and locations
  inquire (file='buildings5',exist=there)

  if (.not. there) return

  open(18,file='buildings5',status='old',form='formatted')

  read(18,*) marker,iver
  read(18,*) numbldgs

  do lb = 1,numbldgs
     read(18,*) itype(lb),blat(lb),blon(lb),bdx(lb),bdy(lb),bdz(lb)
  enddo

  close(18)

  ! Loop through building list

  do lb = 1,numbldgs

     ! If building type entered as positive, location is lat-lon
     ! If building type entered as negative, location is x/y

     if(itype(lb) > 0) then
        call ll_xy(blat(lb),blon(lb),platn,plonn,bx,by)
     else
        bx=blat(lb)
        by=blon(lb)
     endif

     nbdi = max(1,int(bdx(lb) / (xm(2) - xm(1)) + .9))
     nbdj = max(1,int(bdy(lb) / (ym(2) - ym(1)) + .9))

     print*, 'nbdi,nbdj,bdx(lb),bdy(lb) ',nbdi,nbdj,bdx(lb),bdy(lb),xm(2),xm(1)

     bi = (bx - xm(1)) / (xm(2) - xm(1)) + 1.
     bj = (by - ym(1)) / (ym(2) - ym(1)) + 1. 

     ib = int(bi)
     jb = int(bj)

     if (bi - float(ib) > .5) then
        ib1 = ib - (nbdi - 1) / 2
     else
        ib1 = ib - nbdi / 2
     endif

     if (bj - float(jb) > .5) then
        jb1 = jb - (nbdj - 1) / 2
     else
        jb1 = jb - nbdj / 2
     endif

     ib2 = ib1 + nbdi - 1
     jb2 = jb1 + nbdj - 1

     if (ib1 < 3 .or. ib2 > n2 - 2 .or. jb1 < 3 .or. jb2 > n3 - 2) then

        print*, 'Building location out of bounds for this grid.'
        print*, 'ib1,ib2,jb1,jb2,nnxp,nnyp = ',ib1,ib2,jb1,jb2,n2,n3
        stop 'bldg loc'
     endif

     topma_avg = 0.
     do j = jb1,jb2
        do i = ib1,ib2

           topma_t = .25 * (topma(i,j) + topma(i-1,j)  &
                + topma(i,j-1) + topma(i-1,j-1))

           topma_avg = topma_avg + topma_t
        enddo
     enddo
     topma_avg = topma_avg / float(nbdi * nbdj)

     btop = topma_avg + bdz(lb)

     do k = 1,n1-1
        if (zm(k) < btop) kb2 = k + 1
     enddo

     ! Close all 6 apertures for any grid cell in (or under) a building 

     print*, 'jb1,jb2,ib1,ib2,kb2,topma_avg ',jb1,jb2,ib1,ib2,kb2,topma_avg

     do j = jb1,jb2
        do i = ib1,ib2
           do k = 1,kb2
              aru(k,i-1,j) = 0.
              aru(k,i,j)   = 0.
              arv(k,i,j-1) = 0.
              arv(k,i,j)   = 0.
              arw(k,i,j)   = 0.
           enddo

           ! Set the vegetation/soil characteristics for the building roof
           leaf_g(ifm)%patch_area(i,j,1:npatch)=0.
           leaf_g(ifm)%patch_area(i,j,2)=1.
           leaf_g(ifm)%leaf_class(i,j,2)=3.
           leaf_g(ifm)%soil_text(1:nzg,i,j,1:npatch)=11.

        enddo
     enddo

  enddo
end subroutine adap_bldg






subroutine setgv(m1,m2,m3,flpu,flpv,flpw,up,vp,wp,uc,vc,wc)

  implicit none

  integer :: m1,m2,m3
  real, dimension(m1,m2,m3) :: up,vp,wp,uc,vc,wc
  real, dimension(m2,m3) :: flpu,flpv,flpw
  integer :: lpu,lpv,lpw
  integer :: i,j,k

  do j = 1,m3
     do i = 1,m2
        lpu = nint(flpu(i,j))
        do k = 1,lpu-1
           up(k,i,j) = 0.
           uc(k,i,j) = 0.
        enddo
        lpv = nint(flpv(i,j))
        do k = 1,lpv-1
           vp(k,i,j) = 0.
           vc(k,i,j) = 0.
        enddo
        lpw = nint(flpw(i,j))
        do k = 1,lpw-1
           wp(k,i,j) = 0.
           wc(k,i,j) = 0.
        enddo

     enddo
  enddo
end subroutine setgv





subroutine advected(m1,m2,m3,uc,vc,wc,ucb,vcb,wcb,flpu,flpv,flpw,i,j,jdim)

  implicit none

  integer :: m1,m2,m3,i,j,jdim
  real, dimension(m1,m2,m3) :: uc,vc,wc,ucb,vcb,wcb
  real, dimension(m2,m3) :: flpu,flpv,flpw
  integer :: lpw,lpv,lpu
  integer :: lpw1,lpv1,lpu1
  integer :: k

  ! Version 1: zero gradient condition

  ! Set boundary conditions for neighboring U points below and to sides

  k = nint(flpu(i,j))
  ucb(k-1,i,j) = uc(k,i,j)
  ! ucb(k-1,i,j) = 0.

  lpu  = nint(flpu(i,j))
  lpu1 = nint(flpu(i+1,j))

  do k = lpu,lpu1-1
     ucb(k,i+1,j) = uc(k,i,j)
     !   ucb(k,i+1,j) = 0.
  enddo

  lpu  = nint(flpu(i,j))
  lpu1 = nint(flpu(i-1,j))

  do k = lpu,lpu1-1
     ucb(k,i-1,j) = uc(k,i,j)
     !   ucb(k,i-1,j) = 0.
  enddo

  if (jdim == 1) then

     lpu  = nint(flpu(i,j))
     lpu1 = nint(flpu(i,j+1))

     do k = lpu,lpu1-1
        ucb(k,i,j+1) = uc(k,i,j)
        !   ucb(k,i,j+1) = 0.
     enddo

     lpu  = nint(flpu(i,j))
     lpu1 = nint(flpu(i,j-1))

     do k = lpu,lpu1-1
        ucb(k,i,j-1) = uc(k,i,j)
        !   ucb(k,i,j-1) = 0.
     enddo
  endif

  ! Set boundary conditions for neighboring V points below and to sides

  k = nint(flpv(i,j))
  vcb(k-1,i,j) = vc(k,i,j)

  lpv  = nint(flpv(i,j))
  lpv1 = nint(flpv(i+1,j))

  do k =lpv,lpv-1
     vcb(k,i+1,j) = vc(k,i,j)
  enddo

  lpv  = nint(flpv(i,j))
  lpv1 = nint(flpv(i-1,j))

  do k = lpv,lpv1-1
     vcb(k,i-1,j) = vc(k,i,j)
  enddo

  if (jdim == 1) then
     lpv  = nint(flpv(i,j))
     lpv1 = nint(flpv(i,j+1))
     do k = lpv,lpv1-1
        vcb(k,i,j+1) = vc(k,i,j)
     enddo

     lpv  = nint(flpv(i,j))
     lpv1 = nint(flpv(i,j-1))
     do k = lpv,lpv1-1
        vcb(k,i,j-1) = vc(k,i,j)
     enddo
  endif

  ! Set boundary conditions for neighboring W points below and to sides

  k = nint(flpw(i,j))
  wcb(k-1,i,j) = wc(k,i,j)
  !  wcb(k-1,i,j) = 0.

  lpw  = nint(flpw(i,j))
  lpw1 = nint(flpw(i+1,j))
  do k = lpw,lpw1-1
     wcb(k,i+1,j) = wc(k,i,j)
     !   wcb(k,i+1,j) = 0.
  enddo

  lpw =  nint(flpw(i,j))
  lpw1 = nint(flpw(i-1,j))
  do k = lpw,lpw1-1
     wcb(k,i-1,j) = wc(k,i,j)
     !     wcb(k,i-1,j) = 0.
  enddo

  if (jdim == 1) then

     lpw  = nint(flpw(i,j))
     lpw1 = nint(flpw(i,j+1))

     do k = lpw,lpw1-1
        wcb(k,i,j+1) = wc(k,i,j)
        !   wcb(k,i,j+1) = 0.
     enddo

     lpw  = nint(flpw(i,j))
     lpw1 = nint(flpw(i,j-1))

     do k = lpw,lpw1-1
        wcb(k,i,j-1) = wc(k,i,j)
        !   wcb(k,i,j-1) = 0.
     enddo
  endif
end subroutine advected
