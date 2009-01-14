!==========================================================================================!
!==========================================================================================!
subroutine leaf_database(ofn, nlandsea, iaction, lat, lon, idatp)

  use hdf5_utils,  only: shdf5_open_f, shdf5_close_f, shdf5_irec_f

  implicit none

  character(len=*), intent(in) :: ofn
  integer, intent(in) :: nlandsea
  integer :: nio
  integer :: njo
  integer :: nperdeg
  real :: offpix
  integer :: niosh
  integer :: njosh
  integer :: ifiles
  integer :: jfiles
  integer, allocatable, dimension(:,:) :: nump
  integer, allocatable, dimension(:,:) :: numpind1
  integer, allocatable, dimension(:,:) :: numpind2
  integer(kind=4), allocatable, dimension(:,:) :: idato
  integer :: ifile
  integer :: jfile
  real :: wlat1
  real :: wlon1
  real :: rio_full
  real :: rjo_full
  integer :: jo_full
  integer :: io_full
  integer :: ind
  integer, dimension(4*nlandsea) :: ptable
  integer, dimension(4*nlandsea) :: ictable,jctable
  integer :: ind1
  integer :: ind2
  integer :: iwoc
  integer :: isoc
  integer :: isocpt
  integer :: isocpo
  integer :: iwocpt
  integer :: iwocpo
  integer :: iwocph
  character(len=3) :: title1
  character(len=4) :: title2
  character(len=128) :: title3
  logical :: l1
  logical :: l2
  integer :: ndims
  integer, dimension(2) :: idims
  character(len=*), intent(in) :: iaction
  integer :: i,j,i0,i1,j0,j1,ic,jc,dj,di
  integer :: nt
  integer :: ilandsea
  integer :: dq
  integer :: stext,best
  real, dimension(3,nlandsea), intent(in) :: lat
  real, dimension(3,nlandsea), intent(in) :: lon
  
  integer, dimension(0:20,nlandsea) :: hgramtypes
  
  integer, intent(out), dimension(nlandsea) :: idatp

  ! Initialize idatp
  hgramtypes = 0

  ! Read header file

  print*,trim(ofn)//'HEADER'
  

  open(29,file=trim(ofn)//'HEADER',form='formatted',status='old')
  read(29,*) nio, njo, nperdeg
  close(29)

  ! Compute number of pixels in a shift to adjacent file (niosh and njosh).
  ! Compute number of files in database that span all latitudes and
  ! longitudes on earth. [This algorithm will change when multiple resolutions
  ! of the SRTM data become available.]
  
  if (mod(nio,nperdeg) == 2) then
     offpix = .5
     niosh = nio - 2
     njosh = njo - 2
  else
     offpix = 0.
     niosh = nio - 1
     njosh = njo - 1
  endif

  ! Compute number of files in input dataset that span all latitudes and
  ! longitudes on earth.  
  
  ifiles = 360 * nperdeg / niosh
  jfiles = 180 * nperdeg / njosh
  
  ! Allocate 5 arrays.
  
  allocate (nump    (ifiles,jfiles))
  allocate (numpind1(ifiles,jfiles))
  allocate (numpind2(ifiles,jfiles))
  allocate (idato   (nio,njo))
  
  do jfile = 1,jfiles
     do ifile = 1,ifiles
        nump(ifile,jfile) = 0
     enddo
  enddo

  ! Get file index (ifile,jfile) within full dataset and count number of p 
  ! points (nump) that occur in each file

  do ilandsea = 1, nlandsea

     ! Add all four corners to the search process
     do ic=2,3
        do jc=2,3

           ! Set the lat and lon for file lookup
           wlat1 = max(-89.9999,min(89.9999,lat(jc,ilandsea)))
           wlon1 = lon(ic,ilandsea)
           if(wlon1 >= 180.) wlon1 = wlon1 - 360.
           if(wlon1 < -180.) wlon1 = wlon1 + 360.
           wlon1 = max(-179.9999,min(179.9999,wlon1))
           
           rio_full = (wlon1 + 180.) * real(nperdeg) ! must ignore pixel offset here
           rjo_full = (wlat1 +  90.) * real(nperdeg) ! must ignore pixel offset here
           
           io_full = int(rio_full)
           jo_full = int(rjo_full)
           
           ifile = io_full / niosh + 1
           jfile = jo_full / njosh + 1
           
           ! Summation to # pts filled by file (ifile,jfile) in dataset 
           nump(ifile,jfile) = nump(ifile,jfile) + 1  
           
        enddo
     enddo
     
  enddo
  
  ! Set up array index values for ptable array
  
  ind = 1
  do jfile = 1,jfiles
     do ifile = 1,ifiles
        numpind1(ifile,jfile) = ind
        numpind2(ifile,jfile) = ind
        ind = ind + nump(ifile,jfile)
     enddo
  enddo
  
  ! Fill ptable array

  do ilandsea = 1, nlandsea
     
     ! Add all four corners
     do ic=2,3
        do jc=2,3
           
           wlat1 = max(-89.9999,min(89.9999,lat(jc,ilandsea)))
           wlon1 = lon(ic,ilandsea)
           
           if(wlon1 >= 180.) wlon1 = wlon1 - 360.
           if(wlon1 < -180.) wlon1 = wlon1 + 360.
           wlon1 = max(-179.9999,min(179.9999,wlon1))
           
           rio_full = (wlon1 + 180.) * real(nperdeg) ! must ignore pixel offset here
           rjo_full = (wlat1 +  90.) * real(nperdeg) ! must ignore pixel offset here
           
           io_full = int(rio_full)
           jo_full = int(rjo_full)
           
           ifile = io_full / niosh + 1
           jfile = jo_full / njosh + 1
           
           ind = numpind2(ifile,jfile)
           ptable(ind) = ilandsea
           ictable(ind)  = ic
           jctable(ind)  = jc
           numpind2(ifile,jfile) = numpind2(ifile,jfile) + 1 ! Sums to ending index
           
        enddo
     enddo
  enddo

  ! Read files and extract data
  
  do jfile = 1,jfiles
     do ifile = 1,ifiles

        ind1 = numpind1(ifile,jfile)
        ind2 = numpind2(ifile,jfile)
        
        if (ind2 > ind1) then
           ! SW longitude of current file
           iwoc = (ifile - 1) * niosh/nperdeg - 180
           ! SW latitude of current file
           isoc = (jfile - 1) * njosh/nperdeg -  90

           ! Construct filename
           isocpt = abs(isoc) / 10
           isocpo = abs(isoc) - isocpt*10
           iwocph = abs(iwoc) / 100
           iwocpt = (abs(iwoc) - iwocph * 100) / 10
           iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10
           
           if (isoc >= 0) then
              write(title1,'(2i1,a1)') isocpt,isocpo,'N'
           else
              write(title1,'(2i1,a1)') isocpt,isocpo,'S'
           endif
           
           if (iwoc >= 0) then
              write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
           else
              write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
           endif
           
           title3 = trim(ofn)//title1//title2//'.h5'
           inquire(file=trim(title3),exist=l1,opened=l2)

           ! Read file

           if (l1) then
              print*, 'getting file ',trim(title3)    
              
              call shdf5_open_f(trim(title3),'R')

              ndims = 2
              idims(1) = nio
              idims(2) = njo
              
              if (trim(iaction) == 'leaf_class') then
                 call shdf5_irec_f(ndims,idims,'oge2',ivara=idato)
              elseif (trim(iaction) == 'soil_text') then
                 call shdf5_irec_f(ndims,idims,'fao',ivara=idato)
              else
                 print*, 'incorrect action specified in leaf_database'
                 print*, 'stopping run'
                 stop 'stop landuse_input1'
              endif

              call shdf5_close_f()
           else
              print*, 'In landuse_input, ',iaction,' file is missing'
              print*, 'Filename = ',trim(title3)
              print*, 'Stopping model run'
              stop 'stop_landuse_input2'
           endif

           do ind = ind1,ind2-1

              ilandsea = ptable(ind)
              ic       = ictable(ind)
              jc       = jctable(ind)
              
              ! ic=2,jc=2 Top Left
              if (ic==2 .and. jc==2) then
                 nt=0
                 call lat2index(j0,nperdeg,njosh,offpix,nt,0     ,lat(2,ilandsea))
                 call lat2index(j1,nperdeg,njosh,offpix,nt,1     ,lat(3,ilandsea))
                 nt=0
                 call lon2index(i0,nperdeg,niosh,offpix,nt,0     ,lon(2,ilandsea))
                 call lon2index(i1,nperdeg,niosh,offpix,nt,niosh ,lon(3,ilandsea))
                 
              ! ic=2,jc=3 Bot Left  - Go up    (j then i)
              elseif (ic==2 .and. jc==3) then
                 nt=0;
                 call lat2index(j0,nperdeg,njosh,offpix,nt,0    ,lat(3,ilandsea))
                 call lat2index(j1,nperdeg,njosh,offpix,nt,njosh,lat(2,ilandsea))
                 nt=0;
                 call lon2index(i0,nperdeg,niosh,offpix,nt,0    ,lon(2,ilandsea))
                 call lon2index(i1,nperdeg,niosh,offpix,nt,niosh,lon(3,ilandsea))
          
              ! ic=3,jc=2 Top Right - Go left  (i then j)
              elseif (ic==3 .and. jc==2) then
                 nt=0;
                 call lat2index(j0,nperdeg,njosh,offpix,nt,0,lat(2,ilandsea))
                 call lat2index(j1,nperdeg,njosh,offpix,nt,1,lat(3,ilandsea))
                 nt=0;
                 call lon2index(i0,nperdeg,niosh,offpix,nt,0,lon(3,ilandsea))
                 call lon2index(i1,nperdeg,niosh,offpix,nt,1,lon(2,ilandsea))
          
              ! ic=3,jc=3 Bottom Right  - Go down  (j then i)
              elseif (ic==3 .and. jc==3) then
                 nt=0;
                 call lat2index(j0,nperdeg,njosh,offpix,nt,0,lat(3,ilandsea))
                 call lat2index(j1,nperdeg,njosh,offpix,nt,njosh,lat(2,ilandsea))
                 nt=0;
                 call lon2index(i0,nperdeg,niosh,offpix,nt,0    ,lon(3,ilandsea))
                 call lon2index(i1,nperdeg,niosh,offpix,nt,1    ,lon(2,ilandsea))
                 

              else
                 print*,"INDEXING IS FLAWED"
                 print*,"ABORTING"
                 stop

              endif

              if (trim(iaction) == 'leaf_class' ) then

                 dj=1
                 di=1
                 if(j0>j1)dj=-1
                 if(i0>i1)di=-1

                 do j=j0,j1,dj
                    do i=i0,i1,di
                       call datp2datq(idato(i,j),dq)
                       hgramtypes(dq,ilandsea) = hgramtypes(dq,ilandsea) + 1
                    enddo
                 enddo
                 
              elseif(trim(iaction) == 'soil_text') then
                 dj=1
                 di=1
                 if(j0>j1)dj=-1
                 if(i0>i1)di=-1
                 do j=j0,j1,dj
                    do i=i0,i1,di
                       call datp2datsoil(idato(i,j),dq)
                       hgramtypes(dq,ilandsea) = hgramtypes(dq,ilandsea) + 1
                    enddo
                 enddo
              endif
              
!              print*,"POLY: ",ilandsea," IC: ",ic," JC: ",jc
!              print*,"I:  ",i0,"-",i1,niosh," J ",j0,"-",j1,njosh
!              print*,"LON:",lon(2,ilandsea),lon(3,ilandsea)," LAT: ",lat(2,ilandsea),lat(3,ilandsea)

           enddo

        endif
     enddo
  enddo

  ! If this is a land percentage calculation then compare the first two classes to all others
  if (trim(iaction) == 'leaf_class' ) then
     
     do i=1,nlandsea
        idatp(i) = int(100.* real(sum(hgramtypes(2:20,i))) /real(sum(hgramtypes(:,i))))
     enddo
     
  elseif(trim(iaction) == 'soil_text') then
     
     do i=1,nlandsea
        best=0
        stext=0
        do j=0,20
           if(hgramtypes(j,i)>best)then
              stext=j
              best=hgramtypes(j,i)
           endif
        enddo
        idatp(i) = stext
     enddo

  endif


  deallocate(nump,numpind1,numpind2,idato)

  return
end subroutine leaf_database

!==========================================================================================!
!==========================================================================================!

subroutine lat2index(jo,nperdeg,njosh,offpix,nthtile,alt,lat)
  
  implicit none

  real :: lat,wlat1
  real :: rjo_full
  real :: offpix
  integer :: alt
  integer :: nperdeg,njosh
  integer :: jo_full
  integer :: jo1,jo2
  real    :: wjo1,wjo2
  integer,intent(out) :: jo
  integer :: nthtile,thistile
  

  wlat1 = max(-89.9999,min(89.9999,lat))
  
  rjo_full = (wlat1 +  90.) * real(nperdeg) ! must ignore pixel offset here

  jo_full = int(rjo_full + 0.001)

  thistile = ceiling(rjo_full/real(njosh))

  if(nthtile==0) then
     
     nthtile = ceiling(rjo_full/real(njosh))
     
     jo1 = mod(jo_full,njosh) + 1
     
  else

     if(thistile/=nthtile) then
        jo = alt
        return
     else
        jo1 = mod(jo_full,njosh) + 1
     end if
     
  endif

  
  wjo2 = rjo_full - real(jo_full) + offpix
  
  ! At this point, io1, jo1, wio2, and wjo2 are correct if offpix = 0,
  ! but need correction if offpix = .5
              
  if (wjo2 > 1.) then
     wjo2 = wjo2 - 1.
     jo1 = jo1 + 1
  endif
  
  ! This is end of correction
  jo2 = jo1 + 1

  wjo1 = 1. - wjo2
  
  jo = jo2
  
  if (wjo2 < .5) jo = jo1

  if (jo==njosh+1) jo=njosh

  if (jo .lt. 1 .or. jo .gt. njosh) then
     print*,"JO IS MESSED UP",jo,njosh,lat,wlat1,rjo_full,alt 
     print*,"TRY SHIFTING YOUR GRID CENTER EVER SO SLIGHTLY, IE +-0.01"
     print*,"THIS ALGORITHM DOES NOT LIKE IT WHEN GRID POINTS"
     print*,"LIE EXACTLY ON TILE BOUNDARIES"
     stop
  endif


  
  return 
end subroutine lat2index

!=========================================================================================!

subroutine lon2index(io,nperdeg,niosh,offpix,nthtile,alt,lon)
  
  implicit none

  real    :: lon
  real    :: rio_full
  real    :: offpix
  integer :: nperdeg
  integer :: io_full
  integer :: io1,io2
  real    :: wio1,wio2
  integer,intent(out) :: io
  integer :: nthtile,thistile
  integer :: alt,niosh

  if (lon >=  180.) lon = lon - 360.
  if (lon <= -180.) lon = lon + 360.
  
  rio_full = (lon + 180.) * real(nperdeg) ! must ignore pixel offset here
  
  io_full = int(rio_full + 0.001)

  thistile = ceiling(rio_full/real(niosh))

  if(nthtile==0) then
     
     nthtile = ceiling(rio_full/real(niosh))
     
     io1 = mod(io_full,niosh) + 1
     
  else

     if(thistile/=nthtile) then
        io = alt
        return
     else
        io1 = mod(io_full,niosh) + 1
     end if
     
  endif

  wio2 = rio_full - real(io_full) + offpix


  ! At this point, io1, jo1, wio2, and wjo2 are correct if offpix = 0,
  ! but need correction if offpix = .5
  
  if (wio2 > 1.) then
     wio2 = wio2 - 1.
     io1 = io1 + 1
  endif
  
  ! This is end of correction
  
  io2 = io1 + 1
  
  wio1 = 1. - wio2
  
  ! Use nearest data point - do not interpolate
  
  io = io2
  
  if (wio2 < .5) io = io1

  if (io==niosh+1) io=niosh

  if (io .lt. 1 .or. io .gt. niosh) then
     print*,"IO IS MESSED UP",io,niosh,lon,rio_full,io_full,nthtile,alt
     print*,"TRY SHIFTING YOUR GRID CENTER EVER SO SLIGHTLY, IE +-0.01"
     print*,"THIS ALGORITHM DOES NOT LIKE IT WHEN GRID POINTS"
     print*,"LIE EXACTLY ON TILE BOUNDARIES"
     stop
  endif
  return
end subroutine lon2index


!==========================================================================================!
!==========================================================================================!
subroutine datp2datq(datp,datq)

! This subroutine maps the input datp classes to a smaller set datq
! which represents the full set of LEAF-2 or LEAF-3 classes for which 
! LSP values are
! defined.

   implicit none

   integer, intent(in) :: datp
   integer, intent(out) :: datq


!  Olson Global Ecosystems dataset OGE_2 (96 classes) mapped to LEAF-3 classes
!  (see leaf3_document).
! 97 & 98  - not used
! 99       - is Goode Homolosine empty space
! 100      - is missing data
! Map all of these to ocean (datq=0)

! THE FOLLOWING ARRAY WAS FOUND IN DATP_DATQ FOM THE LEAF3_INIT.F90
! SUBROUTINE. IT DIFFERS FROM THE ONE IMPORTED BY DMM
! ANY IDEAS ON WHICH ONE WE SHOULD USE???
!
! It shouldn't matter, the only "vegetation" that ED is looking for is ocean,
! which has the same index in LEAF-2 and LEAF-3. 
   integer, parameter, dimension(0:100) :: catb=(/ &
   !-------------------------------------------!
   !          1  2  3  4  5  6  7  8  9  0
                                         0,  & !   0
             19, 8, 4, 5, 6, 7, 9, 3,11,16,  & !  10
             10, 2,17, 1, 0,12,13,14,18, 4,  & !  20
              4, 4,14,14, 6, 6, 4, 7, 7,15,  & !  30
             15, 6, 7, 7,15,16,16,16,16, 8,  & !  40
              8, 8,18,17,17,12,12, 7,10, 3,  & !  50
             10,10,11,14,18,18,18,18,13, 6,  & !  60
              5, 4,11,12, 0, 0, 0, 0, 3, 2,  & !  70
              3,20, 0,17,17,17, 4,14, 7, 3,  & !  80
              3, 3, 3, 3, 3, 3, 8,12, 7, 6,  & !  90
             18,15,15,15, 4, 5, 0, 0, 0, 0   /)! 100
   !-------------------------------------------!     
   !          1  2  3  4  5  6  7  8  9  0           

   datq = catb(datp)

   return
 end subroutine datp2datq
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine datp2datsoil(datp,datsoil)

! This subroutine maps the input datp soil classes to a smaller set datsoil
! which represents the full set of LEAF-2 classes for which soil parameter
! values are defined.

implicit none
integer, intent(in)  :: datp
integer, intent(out) :: datsoil
! (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our
! own missing value designated 133, to the USDA soil textural classes.  FAO
! classes [0] (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not
! used in original FAO dataset), [132] (water), and [133] (our missing value)
! are all mapped to a default class of sandy clay loam in case they happen to
! correspond to a land surface area in the landuse dataset that RAMS uses to
! define land area.  We wrote missing value class 133 to the RAMS FAO files
! whenever a negative value (which is out of range of defined values) was
! found in the original FAO dataset, which occurred in about 2.6% of the
! pixels.  For the remaining FAO classes, a cross reference table to Zobler
! soil texture classes that was provided, plus our own cross referencing table
! from Zobler to USDA classes listed below, provides the mapping from FAO to
! USDA.  In this mapping, we use only organic USDA classes and omit nonorganic
! classes since most soils do contain organic matter, and organic content
! information is not provided in the Zobler classes.

!  Zobler Class              USDA Class

!  1  Coarse                 2  Loamy sand
!  2  Medium                 4  Silt loam
!  3  Fine                   8  Clay loam
!  4  Coarse-medium          3  Sandy loam
!  5  Coarse-fine            6  Sandy clay loam
!  6  Medium-fine            7  Silty clay loam
!  7  Coarse-medium-fine     6  Sandy clay loam
!  8  Organic matter         5  Loam

!                            1  Sand (not used)
!                            9  Sandy clay (not used)
!                           10  Silty clay (not used)
!                           11  Clay (not used)
!                           12  Peat (not used)
integer, parameter, dimension(0:133) :: catb=(/ &
                                      6  & !0
         , 6, 4, 4, 7, 7, 8, 6, 4, 4, 4  & !10
         , 7, 4, 4, 4, 8, 4, 8, 4, 4, 8  & !20
         , 4, 2, 4, 4, 4, 4, 6, 8, 8, 8  & !30
         , 4, 8, 8, 2, 6, 4, 7, 4, 4, 3  & !40
         , 4, 6, 7, 4, 4, 4, 4, 4, 4, 4  & !50
         , 4, 4, 4, 4, 4, 4, 2, 4, 4, 2  & !60
         , 4, 3, 4, 2, 7, 6, 4, 4, 6, 8  & !70
         , 8, 7, 2, 5, 4, 5, 6, 6, 4, 2  & !80
         , 2, 2, 4, 6, 2, 2, 2, 2, 2, 4  & !90
         , 2, 2, 2, 4, 2, 4, 3, 6, 2, 7  & !100
         , 4, 4, 4, 8, 8, 8, 3, 7, 4, 4  & !110
         , 4, 3, 6, 4, 2, 4, 4, 4, 2, 2  & !120
         , 2, 4, 6, 4, 4, 7, 7, 6, 3, 2  & !130
         , 2, 6, 6 /)                      !133

datsoil = catb(datp)

return
end subroutine datp2datsoil
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_soil_depth()
   use soil_coms, only: soildepth_db,isoildepthflg,nlon_lyr,nlat_lyr,layer_index,slz
   use grid_coms, only: nzg
   implicit none
   real, dimension(nlat_lyr,nlon_lyr) :: soil_depth
   logical                            :: found
   real                               :: rdum
   integer                            :: ilat,k
!------------------------------------------------------------------------------------------!


!----- First thing, allocate the scratch soildepth index layer and initialize it. ---------!
   if (allocated(layer_index)) deallocate(layer_index)
   allocate(layer_index(nlat_lyr,nlon_lyr))
   
!----- If I don't want to use the soil depth database, I set layer index to 1 and leave. --!
   if (isoildepthflg /= 1) then 
      layer_index=1
      return
   else !Simply initialize layer_index to 2 layers, the minimum it is allowed.
      layer_index=nzg-1
   end if

!----- Opening the file, provided it exists. If it doesn't, abort the run.-----------------!
   inquire (file=trim(soildepth_db),exist=found)
   if(.not.found) then
      write(unit=*,fmt='(a)') 'Your ED2IN has isoildepthflg=1, however the file'
      write(unit=*,fmt='(a)') trim(soildepth_db)//' doesn''t exist.'
      call fatal_error('File specified in SOILDEPTH_DB not found!'                         &
                      ,'read_soil_depth','leaf_database.f90')
   end if
   open(unit=12,file=trim(soildepth_db),form='formatted',action='read',status='old')

!----- Reading the database. First latitude bin: 90N-89N; first longitude bin: 180W-179W --!
   do ilat=1,nlat_lyr
      read(unit=12,fmt=*) rdum,soil_depth(ilat,1:nlon_lyr)
   end do
   close(unit=12,status='keep')

!----- Converting it to metres and using the same notation as ED (negative values) --------!
   soil_depth=-soil_depth*0.01

!----- Do something in case the atmospheric land cell happens to fall in an FAO water cell !
   where (soil_depth == 0.0) soil_depth = 1.0

!----- Loop through the layers, and find points in which the soil is shallower. -----------!
   do k=nzg-2,1,-1
      where (soil_depth < slz(k+1)) layer_index = k
   end do
   
   return
end subroutine read_soil_depth
