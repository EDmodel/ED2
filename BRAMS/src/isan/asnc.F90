!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine nc_pressure_stage(n1,n2,nhem,glat,glon,glat2,glon2)

  use isan_coms
  use rconstants
#if USENC
  use netcdf
#endif
  implicit none

  integer :: n1,n2,nhem
  real :: glat(n1,n2),glon(n1,n2),glat2(n1,n2),glon2(n1,n2)

  real :: fnprx,fnpry,grx,gry,gglat,gglon,thmax,thmin  &
       ,xswlon_east,cntlon_east,rr
  integer :: i,j,k,lv,ifm,loop,n,iunit
  real, external :: rs

  !------ Netcdf specific variables -------------------!
  integer,parameter :: mxdims=10
  integer :: nf_stat                  ! error flag
  integer :: var_id                   ! variable number
  character(len=30) :: dummy_name     ! dummy variable name
  integer :: xtype                    ! data type of variable (real,double,int,etc)
  integer :: ndims                    ! variable's number of dimensions
  integer :: natts                    ! number of attributes
  integer,dimension(mxdims) :: dim_id      ! the dimension's number
  character(len=30) :: dim_name(mxdims)    ! the dimension's name
  integer,dimension(mxdims) :: dim_len     ! the dimension's length
  real :: last
  real(kind=8),allocatable :: lat(:)
  real(kind=8),allocatable :: lon(:)
  real(kind=8),allocatable :: lev(:)
  !----------------------------------------------------!

  ! Read the header of input pressure file.
#if USENC
    print 91,innpr(1:len_trim(innpr))
  91 format(//,' Reading netcdf pressure gridded data',/,a,//)


    nf_stat = nf90_open(innpr,0,iunit)
    if (nf_stat.ne.0) then
       print*,"netcdf input failed at checkpoint 0"
       stop
    endif

    inproj = 1


    ! Things we need...
    ! nprz: Number of pressure levels  integer
    ! nprx: Number of longitude points integer
    ! npry: Number of latitude points  integer
    ! xswlon: The smallest longitude   real
    ! xswlat: The smallest latitude    real
    ! gdatdx: The grid resolution in the x (longitude) real
    ! gdatdy: The grid resolution in the y (latitude)  real
    ! -----------------------------------------------------

    ! --- Open latitudinal iformation ----


    nf_stat = nf90_inq_varid(iunit,'lat',var_id)
    if (nf_stat.ne.0) then
       print*,"netcdf input failed at checkpoint 1"
       stop
    endif
    nf_stat = nf90_Inquire_Variable(iunit,var_id,dummy_name, &
         xtype,ndims,dim_id,natts)
    nf_stat = nf90_inquire_dimension(iunit,dim_id(1),dim_name(1),dim_len(1))
    allocate(lat(dim_len(1)))
    nf_stat = nf90_get_var(iunit,var_id,lat)
    npry = dim_len(1)
    xswlat = 90
    gdatdy = 0.0
    last = real(lat(1)-(lat(2)-lat(1)))
    do i=1,npry
       xswlat = min(xswlat,real(lat(i)))
       gdatdy = gdatdy + (real(lat(i))-last)
       last=real(lat(i))
    enddo
    gdatdy = gdatdy/real(npry)

    ! --- Open longitudinal information ----

    nf_stat = nf90_inq_varid(iunit,'lon',var_id)
    if (nf_stat.ne.0) then
       print*,"netcdf input failed at checkpoint 1"
       stop
    endif
    nf_stat = nf90_Inquire_Variable(iunit,var_id,dummy_name, &
         xtype,ndims,dim_id,natts)
    nf_stat = nf90_inquire_dimension(iunit,dim_id(1),dim_name(1),dim_len(1))
    allocate(lon(dim_len(1)))
    nf_stat = nf90_get_var(iunit,var_id,lon)
    nprx = dim_len(1)
    xswlon = 360
    gdatdx = 0.0
    last = real(lon(1)-(lon(2)-lon(1)))
    do i=1,nprx
       gdatdx = gdatdx + (real(lon(i))-last)
       
       last=real(lon(i))
       if(real(lon(i)).gt.180) lon(i)=lon(i)-360
       xswlon = min(xswlon,real(lon(i)))
    enddo
    gdatdx = gdatdx/real(nprx)


    ! --- Open pressure level information ----

    nf_stat = nf90_inq_varid(iunit,'lev',var_id)
    if (nf_stat.ne.0) then
       print*,"netcdf input failed at checkpoint 1"
       stop
    endif
    nf_stat = nf90_Inquire_Variable(iunit,var_id,dummy_name, &
         xtype,ndims,dim_id,natts)
    nf_stat = nf90_inquire_dimension(iunit,dim_id(1),dim_name(1),dim_len(1))
    allocate(lev(dim_len(1)))
    nf_stat = nf90_get_var(iunit,var_id,lev)

    nprz=dim_len(1)
    
    do i=1,dim_len(1)
       levpr(i)=real(lev(nprz-i+1))
    enddo
    

    deallocate(lat,lon,lev)


    ! Check pressure data domain size and location

    if (inproj/=1.and.nhem>0) then
       print*,'You must input a lat-lon pressure grid '  &
            ,'to run a global simulation !!!!' 
       stop 'glob-no-press'
    endif


    if (inproj.eq.1) then

       ! If necessary, convert longitude specification to values in the range
       ! [-180.,180.].

       xswlon=mod(xswlon+900.,360.)-180.
       xnelon=mod(xnelon-900.,-360.)+180.

       if(xswlon.lt.-180..or.xswlon.ge.180.  .or.  &
            xnelon.lt.-180..or.xnelon.gt.180.01.or.  &
            xswlat.lt.-90. .or.xnelat.gt.90.01) then
          print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          print*,'!!! BAD DOMAIN SPECIFICATION   !!'
          print*,'!!! xswlat,xswlon - ',xswlat,xswlon
          print*,'!!! xnelat,xnelon-xswlon - ',xnelat,xnelon-xswlon
          print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          stop 'bad_domain'
       endif

       fnprx=float(nprx)
       fnpry=float(npry)

       ! Set "global domain" flags to determine whether 2 extra rows will be
       ! added to borders of the input gridded data.  iglobew = 1 if doing
       ! full 360 degrees of longitude.  iglobs = 1 if including the south
       ! pole and if iglobew = 1, and iglobn = 1 if including
       ! the north pole and if iglobew = 1.

       iglobew=0
       iglobn=0
       iglobs=0
       if(fnprx*gdatdx.gt.359.9) then
          if(abs(xswlon+180.).gt.1.e-3) then
             print*,'When extracting all longitudes of pressure data,'
             print*,'xswlon must be set to -180.'
             !stop 'xswlon'
          endif
          iglobew=1
       endif
       if(iglobew.eq.1.and.xswlat.lt.-89.9) iglobs=1
       if(iglobew.eq.1.and.xnelat.gt. 89.9) iglobn=1
       idatelin=0
       if(xnelon.lt.xswlon+.01) idatelin=1

       ! Search coarse grid points for any that are outside the bounds of
       ! the pressure data.  If any are found, stop.

       loop=1
       if(nhem>0) loop=2
       do ifm = 1,loop
          do j = 1,n2
             do i = 1,n1
                if(ifm==1) then
                   gglat=glat(i,j)
                   gglon=glon(i,j)
                elseif(ifm==2) then
                   gglat=glat2(i,j)
                   gglon=glon2(i,j)
                endif
                gry = (gglat - xswlat) / gdatdy + 1.
                if(gry.lt.2-iglobs.or.gry.gt.fnpry+iglobn-1) then
                   print*,'Model grid point latitude must be'
                   print*,'at least 1 gridpoint inside'
                   print*,'pressure data area'
                   print*,'ifm,i,j,glat,xswlat,xnelat'  &
                        ,ifm,i,j,gglat,xswlat,xnelat
                   stop 'isan: outside p-data 1'
                endif
                if(idatelin.eq.1.and.gglon.lt.xswlon) gglon=gglon+360.
                grx=(gglon-xswlon)/gdatdx+1.
                if(grx.lt.2-iglobew.or.grx.gt.fnprx+iglobew-1) then
                   print*,'Model grid point longitude must be'
                   print*,'at least 1 gridpoint inside'
                   print*,'pressure data area'
                   print*,'ifm,i,j,glon,xswlon,xnelon'  &
                        ,ifm,i,j,gglon,xswlon,xnelon,iglobew,iglobs,iglobn
                   stop 'isan: outside p-data 2'
                endif
             enddo
          enddo
       enddo
    endif

    ! Deallocate memory for the pressure data 

    if(allocated(p_u))   deallocate(p_u)
    if(allocated(p_v))   deallocate(p_v)
    if(allocated(p_t))   deallocate(p_t)
    if(allocated(p_z))   deallocate(p_z)
    if(allocated(p_r))   deallocate(p_r)
    if(allocated(p_ur))  deallocate(p_ur)
    if(allocated(p_vr))  deallocate(p_vr)
    if(allocated(p_lat)) deallocate(p_lat)
    if(allocated(p_lon)) deallocate(p_lon)
    if(allocated(p_slp)) deallocate(p_slp)
    if(allocated(p_sfp)) deallocate(p_sfp)
    if(allocated(p_sft)) deallocate(p_sft)
    if(allocated(p_snow)) deallocate(p_snow)
    if(allocated(p_sst)) deallocate(p_sst)

    ! Allocate memory for the pressure data

    allocate(p_u(nprx,npry,nprz))
    allocate(p_v(nprx,npry,nprz))
    allocate(p_t(nprx,npry,nprz))
    allocate(p_z(nprx,npry,nprz))
    allocate(p_r(nprx,npry,nprz))

    ! p_ur,p_vr arrays for rotated winds used later

    allocate(p_ur(nprx,npry,nprz))
    allocate(p_vr(nprx,npry,nprz))

    allocate(p_lat(nprx,npry))
    allocate(p_lon(nprx,npry))

    allocate(p_slp(nprx,npry))
    allocate(p_sfp(nprx,npry))
    allocate(p_sft(nprx,npry))
    allocate(p_snow(nprx,npry))
    allocate(p_sst(nprx,npry))


    ! Fill with missing in case they are not present
    p_slp (1:nprx,1:npry)=1.e30
    p_sfp (1:nprx,1:npry)=1.e30
    p_sft (1:nprx,1:npry)=1.e30
    p_snow(1:nprx,1:npry)=1.e30
    p_sst (1:nprx,1:npry)=1.e30


    print 92
  92 format(////,1x,70('*')/  &
         ,'  Access coarse resolution pressure data',/,1X,70('*')///)

    do k=1,nprz
       pnpr(k)=levpr(k)*100
    enddo

    ! Call routine to fill pressure arrays from the chosen dataset.

    call nc_get_press (iunit)


  !!!!!!!! Be careful !!!!!!!!!
    !  Check input humidity variable. Assume that if the max of the field is greater
    !  than 1.1 (allow for some machine roundoff), 
    !  it is specific humidity (in g/kg) which needs to be converted to rh

    !  Netcdf data is in kg/kg

    if(maxval(p_r(1:nprx,1:npry,1:nprz)) > 1.1 ) then
      

       print*,'------------------------------------------------'
       print*,' Converting specific humidity to rh',nprx,npry,nprz
       print*,'------------------------------------------------'
       do k=1,nprz
          do j=1,npry
             do i=1,nprx

                rr=p_r(i,j,k)*0.001 / (1.-p_r(i,j,k)*0.001 )

                p_r(i,j,k)=rr/rs(pnpr(k),p_t(i,j,k))


  !              print*,rr,rs(pnpr(k),p_t(i,j,k)),p_r(i,j,k)



                if(p_r(i,j,k) <= 0.05) p_r(i,j,k)=.05
                if(i.eq.25.and.j.eq.25)print*,k,pnpr(k),p_r(i,j,k),1000*rr,1000*rs(pnpr(k),p_t(i,j,k)),p_t(i,j,k)-t00
             enddo
          enddo
       enddo
    endif

    if( maxval(p_r(1:nprx,1:npry,1:nprz)) < 0.05 )then
       print*,'------------------------------------------------'
       print*,' Converting specific humidity to rh',nprx,npry,nprz
       print*,'------------------------------------------------'
       do k=1,nprz
          do j=1,npry
             do i=1,nprx
                rr=p_r(i,j,k) / (1.-p_r(i,j,k) )
                p_r(i,j,k)=rr/rs(pnpr(k),p_t(i,j,k))
                if(p_r(i,j,k) <= 0.05) p_r(i,j,k)=.05
                if(i.eq.25.and.j.eq.25)print*,k, pnpr(k),p_r(i,j,k), 1000*rr, 1000*rs(pnpr(k),p_t(i,j,k)), p_t(i,j,k)-t00
             enddo
          enddo
       enddo
    endif



    ! Find max-min theta at bottom and top levels
    thmax=1.
    thmin=1000.
    do j=1,npry
       do i=1,nprx
          if(p_t(i,j,nprz).lt.1.e20)  &
               thmax=max(thmax,p_t(i,j,nprz)*(p00/pnpr(nprz))**rocp)
          if(p_t(i,j,1).lt.1.e20)  &
               thmin=min(thmin,p_t(i,j,1)*(p00/pnpr(1))**rocp)
       enddo
    enddo

    print 300,levpr(1),thmin,levpr(nprz),thmax
  300 format(//,' Minimum THETA at ',I4,' mb - ',F8.2/  &
         ,' Maximum THETA at ',I4,' mb - ',F8.2)


    !  Compute lat-lon at input pressure data points

    if (inproj==1) then
       do j=1,npry
          do i=1,nprx
             p_lat(i,j)=xswlat+(j-1)*gdatdy
             p_lon(i,j)=xswlon+(i-1)*gdatdx

          enddo
       enddo
    elseif (inproj==2) then
       print*,'lamb-con:',xswlat,xswlon,cntlat,cntlon,gdatdx,gdatdy
       do j=1,npry
          do i=1,nprx
             call lc2ll(cntlat,cntlon,p_lat(i,j),p_lon(i,j)  &
                  ,float(i),float(j),xswlat,xswlon,gdatdx)
             !         print*,i,j,p_lat(i,j),p_lon(i,j)
          enddo
       enddo
    elseif (inproj==3) then
       print*,'polar:',cntlat,cntlon,gdatdx,gdatdy,xswlat,xswlon
       xswlon_east=xswlon
       if(xswlon<0.) xswlon_east=xswlon+360.
       cntlon_east=cntlon
       if(cntlon<0.) cntlon_east=cntlon+360.
       do j=1,npry
          do i=1,nprx
             call w3fb07(float(i),float(j),xswlat,xswlon_east,gdatdx,cntlon_east  &
                  ,p_lat(i,j),p_lon(i,j))
             !         print*,i,j,p_lat(i,j),p_lon(i,j)
          enddo
       enddo
    endif

    close(iunit)
#else
  write (*,fmt='(a)') 'You haven''t provided the variable NC_LIBS in your include.mk.??? file, '
  write (*,fmt='(a)') ' so you cannot use netCDF as input data, I''ll have to stop here...'
  stop 'nc_pressure_stage'
#endif
  return
end subroutine nc_pressure_stage

! ===================================================================================

subroutine nc_get_press (iunit)

  use isan_coms
#if USENC
  use netcdf
#endif
  implicit none

  integer :: iunit

  real,allocatable::as(:,:)
  integer :: ithere(maxpr,5),isfthere(5)
  character*4 field,idat(5)
  data idat/'T','R','U','V','H'/

  integer :: i,j,k,nv,nvar,misstot,lv,n

  ! ----- Netcdf type variables ----
  integer,parameter :: mxdims = 10
  integer :: nf_stat
  integer :: ndims
  integer :: nvars
  integer :: natts
  integer :: var_id
  integer :: xtype
  integer :: nvdims,latd,lond,levd,timd
  integer,dimension(mxdims) :: var_dim_id,var_dim_len
  character(len=30) :: var_dim_name(mxdims)
  character(len=30) :: dummy_name,var_name
  real,        allocatable :: in_var_ar4(:,:,:,:)


#if USENC
  !  Read upper air fields

  ! --- Determine the number of variables

  nf_stat = nf90_inquire(iunit,ndims,nvars,natts)
  if (nf_stat.ne.0) print*,"A"


  ! --- Determine the scratch space needed to allocate data arrays

  ! Use T_L
  nf_stat = nf90_inq_varid(iunit,'T_L',var_id)
  if(nf_stat.ne.0) print*,"B"


  nf_stat = nf90_Inquire_Variable(iunit,var_id,dummy_name, &
       xtype,nvdims,var_dim_id,natts)
  if(nf_stat.ne.0) print*,"C"

  do i=1,nvdims
     nf_stat = nf90_inquire_dimension(iunit,var_dim_id(i),var_dim_name(i),var_dim_len(i))

     if(var_dim_name(i).eq.'lon') lond = i
     if(var_dim_name(i).eq.'lat') latd = i
     if(var_dim_name(i).eq.'lev') levd = i
     if(var_dim_name(i).eq.'time') timd = i
 
     if(nf_stat.ne.0) print*,"D"
  enddo

  ! --- Allocate the input array
  
  allocate(in_var_ar4(var_dim_len(1),var_dim_len(2),var_dim_len(3),var_dim_len(4)))

  ! Figure out which is which...,lat,lon,lev,time
  ! We want to loop through the slabs
  
  
  do nvar = 1,nvars

     nf_stat = nf90_Inquire_Variable(iunit,nvar,var_name, &
          xtype,ndims,var_dim_id,natts)
     if (nf_stat.ne.0) print*,"E"

      
     if (var_name.eq.'U_L'  ) then
        nf_stat = nf90_get_var(iunit,nvar,in_var_ar4)
        call nc_prfill(nprx,npry,nprz,in_var_ar4,p_u)

     elseif (var_name.eq.'V_L') then
        nf_stat = nf90_get_var(iunit,nvar,in_var_ar4)
        call nc_prfill(nprx,npry,nprz,in_var_ar4,p_v)

     elseif (var_name.eq.'T_L') then
        nf_stat = nf90_get_var(iunit,nvar,in_var_ar4)
        call nc_prfill(nprx,npry,nprz,in_var_ar4,p_t)

     elseif (var_name.eq.'Z3_L') then
        nf_stat = nf90_get_var(iunit,nvar,in_var_ar4)
        call nc_prfill(nprx,npry,nprz,in_var_ar4,p_z)

     elseif (var_name.eq.'Q_L') then
        nf_stat = nf90_get_var(iunit,nvar,in_var_ar4)
        call nc_prfill(nprx,npry,nprz,in_var_ar4,p_r)

     endif

  enddo
        
  
  deallocate(in_var_ar4)

  nf_stat = nf90_close(iunit)

  ! Check for levels that may be all missing

  ithere(1:nprz,1:5)=0
  isfthere(1:5)=0

  do k=1,nprz
     do j=1,npry
        do i=1,nprx
           if(p_t(i,j,k) > 1e20) ithere(k,1)=ithere(k,1)+1
           if(p_r(i,j,k) > 1e20) ithere(k,2)=ithere(k,2)+1
           if(p_u(i,j,k) > 1e20) ithere(k,3)=ithere(k,3)+1
           if(p_v(i,j,k) > 1e20) ithere(k,4)=ithere(k,4)+1
           if(p_z(i,j,k) > 1e20) ithere(k,5)=ithere(k,5)+1
        enddo
     enddo
  enddo
  do nv=1,5
     do k=1,nprz
        if(ithere(k,nv) < nprx*npry) then
           ithere(k,nv)=1
        else
           ithere(k,nv)=0
        endif
     enddo
  enddo

  misstot=0
  do nv=1,5
     do k=1,nprz
        if(ithere(k,nv) == 0) misstot=misstot+1
     enddo
  enddo

  do j=1,npry
     do i=1,nprx
        if(p_slp (i,j) > 1e20) isfthere(1)=isfthere(1)+1
        if(p_sfp (i,j) > 1e20) isfthere(2)=isfthere(2)+1
        if(p_sft (i,j) > 1e20) isfthere(3)=isfthere(3)+1
        if(p_snow(i,j) > 1e20) isfthere(4)=isfthere(4)+1
        if(p_sst (i,j) > 1e20) isfthere(5)=isfthere(5)+1
     enddo
  enddo
  do nv=1,5
     if(isfthere(nv) < nprx*npry) then
        isfthere(nv)=1
     else
        isfthere(nv)=0
     endif
  enddo

  print*,'------------------------------------------------'
  print*,' Missing parameters: 0 = all missing'
  print*,'------------------------------------------------'
  print '(t20,5(a1,6x))',(idat(n),n=1,5)

  do k=1,nprz
     print '(f10.1,t14,5(i7))',pnpr(k),(ithere(k,n),n=1,5)
!     print '(f10.1,t14,5(f7.2))',pnpr(k),(p_r(nprx,n,k),n=npry-5,npry)

  enddo

  print*,'------------------------------------------------'
  print*,' Missing surface parameters: 0 = all missing'
  print*,'------------------------------------------------'
  print*,'              SLP    SFP    SFT    SNOW   SST'

  print '(t14,5(i7))',(isfthere(n),n=1,5)


  if(misstot.gt.0) then
     ! Let's see if we can get creative and make up data for the missing fields

     call press_miss(nprx,npry,nprz,p_u,p_v,p_t,p_z,p_r  &
          ,ithere,maxpr,levpr)

     ! Check again for missing fields

     misstot=0
     do nv=1,5
        do k=1,nprz
           if(ithere(k,nv).eq.0) misstot=misstot+1
        enddo
     enddo

     print*,'------------------------------------------------'
     print*,' After missing parameters check: 0 = all missing'
     print*,'------------------------------------------------'
     print '(t20,5(a1,6x))',(idat(n),n=1,5)

     do k=1,nprz
        print '(f10.1,t14,5(i7))',pnpr(k),(ithere(k,n),n=1,5)
     enddo

  endif
#else
  write (*,fmt='(a)') 'You haven''t provided the variable NC_LIBS in your include.mk.??? file, '
  write (*,fmt='(a)') ' so you cannot use netCDF as input data, I''ll have to stop here...'
  stop 'nc_get_press'
#endif

  return
end subroutine nc_get_press

!***************************************************************************

subroutine nc_prfill (nprx,npry,nprz,xx,dn)
  
  implicit none
  
  integer :: nprx,npry,nprz
  real :: xx(nprx,npry,nprz,1),dn(nprx,npry,nprz)
  
  integer :: i,j,k
  
  do j=1,npry
     do i=1,nprx
        do k=1,nprz
           dn(i,j,k)=xx(i,j,nprz-k+1,1)
           if(dn(i,j,k).lt.-998.) dn(i,j,k)=1.e30
        enddo
     enddo
  enddo
  
  return
end subroutine nc_prfill
