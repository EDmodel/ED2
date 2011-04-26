!==========================================================================================!
!==========================================================================================!
subroutine leaf_database(ofn,nlandsea,iaction,lat,lon,idatp)

   use hdf5_utils , only : shdf5_open_f  & ! subroutine
                         , shdf5_close_f & ! subroutine
                         , shdf5_irec_f  ! ! subroutine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)                         , intent(in)  :: ofn
   character(len=*)                         , intent(in)  :: iaction
   integer                                  , intent(in)  :: nlandsea
   real           , dimension(3,nlandsea)   , intent(in)  :: lat
   real           , dimension(3,nlandsea)   , intent(in)  :: lon
   integer        , dimension(nlandsea)     , intent(out) :: idatp
   !----- Variables describing the files. -------------------------------------------------!
   integer                                                :: nio
   integer                                                :: njo
   integer                                                :: nperdeg
   integer                                                :: niosh
   integer                                                :: njosh
   integer                                                :: ifiles
   integer                                                :: jfiles
   real                                                   :: offpix
   integer                                                :: west
   integer                                                :: south
   !----- Variables for file handling. ----------------------------------------------------!
   integer                                                :: iwoc
   integer                                                :: isoc
   integer                                                :: isocpt
   integer                                                :: isocpo
   integer                                                :: iwocpt
   integer                                                :: iwocpo
   integer                                                :: iwocph
   character(len=3)                                       :: title1
   character(len=4)                                       :: title2
   character(len=128)                                     :: title3
   logical                                                :: l1
   logical                                                :: l2
   integer                                                :: ndims
   integer        , dimension(2)                          :: idims
   integer                                                :: ierr
   !----- Variables for file indexing. ----------------------------------------------------!
   integer        , dimension(:,:)          , allocatable :: nump
   integer        , dimension(:,:,:)        , allocatable :: cnr_ind
   integer(kind=4), dimension(:,:)          , allocatable :: idato
   integer        , dimension(4*nlandsea)                 :: cnr_i1
   integer        , dimension(4*nlandsea)                 :: cnr_i2
   integer        , dimension(4*nlandsea)                 :: cnr_j1
   integer        , dimension(4*nlandsea)                 :: cnr_j2
   integer        , dimension(4*nlandsea)                 :: cnr_ipoly
   integer        , dimension(0:20,nlandsea)              :: hgramtypes
   integer                                                :: ifile,ifile2
   integer                                                :: jfile,jfile2
   integer                                                :: jo_full
   integer                                                :: io_full
   integer                                                :: ind
   integer                                                :: io_loc,jo_loc
   integer                                                :: max_per_file
   integer                                                :: cid,i,j,i1,i2,j1,j2,ic,jc
   integer                                                :: nt
   integer                                                :: ilandsea
   integer                                                :: dq
   integer                                                :: stext,best 
   real                                                   :: rio_full
   real                                                   :: rjo_full
   !---------------------------------------------------------------------------------------!




   !----- Initialize idatp. ---------------------------------------------------------------!
   hgramtypes = 0

   !---------------------------------------------------------------------------------------!
   !     Read header file.                                                                 !
   !---------------------------------------------------------------------------------------!
   write(unit=*,fmt='(2a)') trim(ofn),'HEADER'
   open (unit=29,file=trim(ofn)//'HEADER',form='formatted',status='old')
   read (unit=29,fmt=*,iostat=ierr) nio, njo, nperdeg,west,south
   !----- Check whether this was read sucessfully.  If not, rewind it and assume global. --!
   if (ierr /= 0) then
      rewind(unit=29)
      read  (unit=29,fmt=*) nio, njo, nperdeg
      west  = -180
      south =  -90
   end if
   close(unit=29)
  
   !---------------------------------------------------------------------------------------!
   !     Compute number of pixels in a shift to adjacent file (niosh and njosh).  Compute  !
   ! number of files in database that span all latitudes and longitudes on Earth.  [This   !
   ! algorithm will change when multiple resolutions of the SRTM data become available.]   !
   !---------------------------------------------------------------------------------------!
   if (mod(nio,nperdeg) == 2) then
      offpix = .5
      niosh = nio - 2
      njosh = njo - 2
   else
      offpix = 0.
      niosh = nio - 1
      njosh = njo - 1
   endif

   !---------------------------------------------------------------------------------------!
   !     Compute number of files in input dataset that span all latitudes and longitudes   !
   ! on Earth.  This could be done more efficiently, but for the time being we will only   !
   ! check the southernmost and westernmost points, and assume that the data set could go  !
   ! to the North Pole and the International Date Line.                                    !
   !---------------------------------------------------------------------------------------!
   ifiles = (180 -  west) * nperdeg / niosh
   jfiles = ( 90 - south) * nperdeg / njosh
  
  ! Allocate 5 arrays.
  
  allocate (nump    (ifiles,jfiles))
  allocate (idato   (nio,njo))

  do jfile = 1,jfiles
     do ifile = 1,ifiles
        nump(ifile,jfile) = 0
     enddo
  enddo
  
  ! Get file index (ifile,jfile) within full dataset and count
  ! the number of corners that appear in every file
  max_per_file = 0
  do ilandsea = 1, nlandsea
     do ic=2,3
        do jc=2,3
           call get_file_indices(lat(jc,ilandsea),lon(ic,ilandsea),west,south, &
                nperdeg,niosh,njosh, &
                io_full,jo_full,     &
                io_loc,jo_loc,       &
                ifile,jfile)

           nump(ifile,jfile) = nump(ifile,jfile) + 1
           if(nump(ifile,jfile)>max_per_file) max_per_file=nump(ifile,jfile)

        enddo
     enddo
  enddo
  
  ! Allocate the mapping of file indices to global indices
  allocate (cnr_ind(ifiles,jfiles,max_per_file))
  
  ! Re-zero nump
  nump = 0
  ind  = 0
  do ilandsea = 1, nlandsea
     
     ! Add all four corners to the search process
     do ic=2,3
        do jc=2,3

           ind = ind+1
           i2 = 5-ic
           j2 = 5-jc
           
           ! Get the file indices and the local indices within the file
           
           call get_file_indices(lat(jc,ilandsea),lon(ic,ilandsea),west,south, &
                nperdeg,niosh,njosh, &
                io_full,jo_full,     &
                io_loc,jo_loc,       &
                ifile,jfile)

           nump(ifile,jfile) = nump(ifile,jfile) + 1

           cnr_ind(ifile,jfile,nump(ifile,jfile)) = ind

           cnr_ipoly(ind) = ilandsea
           
           ! This is the local index where you find the starting point
           
           cnr_i1(ind) = io_loc
           cnr_j1(ind) = jo_loc
           
           ! Get the local indices of the other corner

           call get_file_indices(lat(j2,ilandsea),lon(i2,ilandsea),west,south, &
                nperdeg,niosh,njosh, &
                io_full,jo_full,     &
                io_loc,jo_loc,       &
                ifile2,jfile2)

           ! This is the local index where you find the second i point
           
           if (ifile2==ifile) then
              cnr_i2(ind) = io_loc
           else if (ifile2>ifile) then
              cnr_i2(ind) = niosh+1
           else
              cnr_i2(ind) = 2
           endif

           ! This is the local index where you find the second j point
           
           if (jfile2==jfile) then
              cnr_j2(ind) = jo_loc
           else if (jfile2>jfile) then
              cnr_j2(ind) = njosh+1
           else
              cnr_j2(ind) = 2
           endif
           

           
        enddo
     enddo
     
  enddo


  ! Read files and extract data
  
  do jfile = 1,jfiles
     do ifile = 1,ifiles
        
        if(nump(ifile,jfile)>0) then
           

        ! SW longitude of current file
        iwoc = (ifile - 1) * niosh/nperdeg + west
        ! SW latitude of current file
        isoc = (jfile - 1) * njosh/nperdeg + south
        
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
              call fatal_error('Incorrect action specified in leaf_database'               &
                              ,'leaf_database','leaf_database.f90')
           endif
           
           call shdf5_close_f()
        else
           call fatal_error( 'In landuse_input, '//trim(iaction)//' file is missing'       &
                           ,'leaf_database','leaf_database.f90')
        endif


        ! Loop through all of the corner points that are inside this
        ! file

        do cid = 1,nump(ifile,jfile)

           ! get the global index of this point
           
           ind = cnr_ind(ifile,jfile,cid)
           
           ! get the polygon index of this corner

           ilandsea = cnr_ipoly(ind)

           if (cnr_i1(ind)<cnr_i2(ind)) then
              i1 = cnr_i1(ind)
              i2 = cnr_i2(ind)
           else
              i1 = cnr_i2(ind)
              i2 = cnr_i1(ind)
           endif
           
           if (cnr_j1(ind)<cnr_j2(ind)) then
              j1 = cnr_j1(ind)
              j2 = cnr_j2(ind)
           else
              j1 = cnr_j2(ind)
              j2 = cnr_j1(ind)
           endif

           ! now walk throught he point spaces for soil and vegetation class
           
           if (trim(iaction) == 'leaf_class' ) then

              do i=i1,i2
                 do j=j1,j2
                    call datp2datq(idato(i,j),dq)
                    hgramtypes(dq,ilandsea) = hgramtypes(dq,ilandsea) + 1
                 enddo
              enddo
                 
           elseif(trim(iaction) == 'soil_text') then

              do j=j1,j2
                 do i=i1,i2
                    ! call datp2datsoil(idato(i,j),dq)
                    dq = idato(i,j)
                    hgramtypes(dq,ilandsea) = hgramtypes(dq,ilandsea) + 1
                 enddo
              enddo
           endif

        end do
     endif
  end do
end do


  
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
  
  deallocate(nump,idato,cnr_ind)

  return
end subroutine leaf_database

!==========================================================================================!
!==========================================================================================!

subroutine get_file_indices(lat,lon,west,south,nper,ni,nj,io_full,jo_full,io_loc,jo_loc    &
                           ,ifile,jfile)
  
  implicit none
  
  real,intent(in)     :: lat,lon
  integer,intent(in)  :: west,south
  integer,intent(in)  :: nper,ni,nj
  integer,intent(out) :: io_full,jo_full
  integer,intent(out) :: io_loc,jo_loc
  integer,intent(out) :: ifile,jfile
  
  real :: wlat,wlon
  real :: rio_full,rjo_full
  
  ! Set the lat and lon for file lookup
  wlat = max(-89.9999,min(89.9999,lat))
  wlon = lon
  if(wlon >= 180.) wlon = wlon - 360.
  if(wlon < -180.) wlon = wlon + 360.
  wlon = max(-179.9999,min(179.9999,wlon))
  
  rio_full = (wlon -  west) * real(nper) ! must ignore pixel offset here
  rjo_full = (wlat - south) * real(nper) ! must ignore pixel offset here
  
  io_full = int(rio_full)
  jo_full = int(rjo_full)
  
  ifile = io_full / ni + 1
  jfile = jo_full / nj + 1
  
  io_loc = mod(io_full,ni) + 1
  jo_loc = mod(jo_full,nj) + 1
  
  return
end subroutine get_file_indices

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
