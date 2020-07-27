!==========================================================================================!
!==========================================================================================!
subroutine leaf_database(ofn,nsite,nlandsea,iaction,lat,lon,classout,pctout)

   use hdf5_utils , only : shdf5_open_f  & ! subroutine
                         , shdf5_close_f & ! subroutine
                         , shdf5_irec_f  ! ! subroutine
   use soil_coms  , only : nslcon        & ! intent(in)
                         , isoilcol      & ! intent(in)
                         , ed_nstyp      & ! intent(in)
                         , ed_nvtyp      ! ! intent(in)
   use grid_coms  , only : nzg        ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)                          , intent(in)  :: ofn
   character(len=*)                          , intent(in)  :: iaction
   integer                                   , intent(in)  :: nsite
   integer                                   , intent(in)  :: nlandsea
   real           , dimension(    3,nlandsea), intent(in)  :: lat
   real           , dimension(    3,nlandsea), intent(in)  :: lon
   integer , dimension(nsite,nlandsea,nzg), intent(out) :: classout
   real           , dimension(nsite,nlandsea), intent(out) :: pctout
   !----- Variables describing the files. -------------------------------------------------!
   integer                                                 :: nio
   integer                                                 :: njo
   integer                                                 :: nperdeg
   integer                                                 :: niosh
   integer                                                 :: njosh
   integer                                                 :: ifiles
   integer                                                 :: jfiles
   real                                                    :: offpix
   integer                                                 :: west
   integer                                                 :: south
   !----- Variables for file handling. ----------------------------------------------------!
   integer                                                 :: iwoc
   integer                                                 :: isoc
   integer                                                 :: isocpt
   integer                                                 :: isocpo
   integer                                                 :: iwocpt
   integer                                                 :: iwocpo
   integer                                                 :: iwocph
   character(len=3)                                        :: title1
   character(len=4)                                        :: title2
   character(len=128)                                      :: title3
   logical                                                 :: l1
   logical                                                 :: l2
   integer                                                 :: ndims
   integer        , dimension(3)                           :: idims
   integer                                                 :: ierr
   !----- Variables for file indexing. ----------------------------------------------------!
   integer        , dimension(:,:)           , allocatable :: nump
   integer        , dimension(:,:,:)         , allocatable :: cnr_ind
   integer(kind=4), dimension(:,:)           , allocatable :: idato
   integer(kind=4), dimension(:,:,:)         , allocatable :: idato2
   integer        , dimension(4*nlandsea)                  :: cnr_i1
   integer        , dimension(4*nlandsea)                  :: cnr_i2
   integer        , dimension(4*nlandsea)                  :: cnr_j1
   integer        , dimension(4*nlandsea)                  :: cnr_j2
   integer        , dimension(4*nlandsea)                  :: cnr_ipoly
   integer        , dimension(:,:)           , allocatable :: class_count
   integer        , dimension(:)             , allocatable :: rankpct
   real           , dimension(:)             , allocatable :: fraction
   real                                                    :: total_count
   integer                                                 :: ifile
   integer                                                 :: ifile2
   integer                                                 :: jfile
   integer                                                 :: jfile2
   integer                                                 :: jo_full
   integer                                                 :: io_full
   integer                                                 :: ind
   integer                                                 :: io_loc
   integer                                                 :: jo_loc
   integer                                                 :: max_per_file
   integer                                                 :: cid
   integer                                                 :: i
   integer                                                 :: irank
   integer                                                 :: itext
   integer                                                 :: j
   integer                                                 :: i1
   integer                                                 :: i2
   integer                                                 :: j1
   integer                                                 :: j2
   integer                                                 :: ic
   integer                                                 :: jc
   integer                                                 :: ilandsea
   integer                                                 :: dq
   integer                                                 :: ed_nctyp
   integer                                                 :: z
   !------ External functions. ------------------------------------------------------------!
   integer, external                                       :: find_rank
   !---------------------------------------------------------------------------------------!




   !----- Allocate some auxiliary variables and initialise them. --------------------------!
   ed_nctyp = max(ed_nvtyp,ed_nstyp)
   allocate(class_count(0:ed_nctyp,nlandsea))
   allocate(rankpct(ed_nstyp))
   allocate(fraction(ed_nstyp))
   class_count(:,:) = 0
   rankpct    (:)   = 0
   fraction   (:)   = 0.
   !---------------------------------------------------------------------------------------!



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
   close(unit=29,status='keep')
   !---------------------------------------------------------------------------------------!



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
   end if
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Compute number of files in input dataset that span all latitudes and longitudes   !
   ! on Earth.  This could be done more efficiently, but for the time being we will only   !
   ! check the southernmost and westernmost points, and assume that the data set could go  !
   ! to the North Pole and the International Date Line.                                    !
   !---------------------------------------------------------------------------------------!
   ifiles = (180 -  west) * nperdeg / niosh
   jfiles = ( 90 - south) * nperdeg / njosh
   !---------------------------------------------------------------------------------------!



   !----- Allocate nump and idato. --------------------------------------------------------!
   allocate (nump    (ifiles,jfiles))
   allocate (idato   (nio,njo))
   allocate ( idato2 (nio,njo,nzg))

   do jfile = 1,jfiles
      do ifile = 1,ifiles
         nump(ifile,jfile) = 0
      end do
   end do


   !---------------------------------------------------------------------------------------!
   !      Get file index (ifile,jfile) within full dataset and count the number of corners !
   ! that appear in every file.                                                            !
   !---------------------------------------------------------------------------------------!
   max_per_file = 0
   do ilandsea = 1, nlandsea
      do ic=2,3
         do jc=2,3
            call get_file_indices(lat(jc,ilandsea),lon(ic,ilandsea),west,south,nperdeg     &
                                 ,niosh,njosh,io_full,jo_full,io_loc,jo_loc,ifile,jfile)

            nump(ifile,jfile) = nump(ifile,jfile) + 1
            if(nump(ifile,jfile)>max_per_file) max_per_file = nump(ifile,jfile)
         end do
      end do
   end do

   !----- Allocate the mapping of file indices to global indices. -------------------------!
   allocate (cnr_ind(ifiles,jfiles,max_per_file))
   
   !----- Re-zero nump. -------------------------------------------------------------------!
   nump(:,:) = 0
   ind       = 0
   do ilandsea = 1, nlandsea
      
      !------ Add all four corners to the search process. ---------------------------------!
      do ic=2,3
         do jc=2,3

            ind = ind+1
            i2  = 5-ic
            j2  = 5-jc

            !----- Get the file indices and the local indices within the file. ------------!
            call get_file_indices(lat(jc,ilandsea),lon(ic,ilandsea),west,south,nperdeg     &
                                 ,niosh,njosh,io_full,jo_full,io_loc,jo_loc,ifile,jfile)

            nump(ifile,jfile) = nump(ifile,jfile) + 1

            cnr_ind(ifile,jfile,nump(ifile,jfile)) = ind

            cnr_ipoly(ind) = ilandsea
            
            !----- This is the local index where you find the starting point. -------------!
            cnr_i1(ind) = io_loc
            cnr_j1(ind) = jo_loc

            !----- Get the local indices of the other corner. -----------------------------!
            call get_file_indices(lat(j2,ilandsea),lon(i2,ilandsea),west,south,nperdeg      &
                                 ,niosh,njosh,io_full,jo_full,io_loc,jo_loc,ifile2,jfile2)

            !----- This is the local index where you find the second i point. -------------!
            if (ifile2 == ifile) then
               cnr_i2(ind) = io_loc
            else if (ifile2 > ifile) then
               cnr_i2(ind) = niosh+1
            else
               cnr_i2(ind) = 2
            end if

            !----- This is the local index where you find the second j point. -------------!
            if (jfile2 == jfile) then
               cnr_j2(ind) = jo_loc
            else if (jfile2 > jfile) then
               cnr_j2(ind) = njosh+1
            else
               cnr_j2(ind) = 2
            end if
         end do
      end do
   end do


   !----- Read files and extract data. ----------------------------------------------------!
   do jfile = 1,jfiles
      do ifile = 1,ifiles
         
         if  (nump(ifile,jfile) > 0) then
            !----- SW longitude of current file. ------------------------------------------!
            iwoc = (ifile - 1) * niosh/nperdeg + west
            !----- SW latitude of current file. -------------------------------------------!
            isoc = (jfile - 1) * njosh/nperdeg + south
            
            !----- Construct filename. ----------------------------------------------------!
            isocpt = abs(isoc) / 10
            isocpo = abs(isoc) - isocpt*10
            iwocph = abs(iwoc) / 100
            iwocpt = (abs(iwoc) - iwocph * 100) / 10
            iwocpo = abs(iwoc) - iwocph * 100 - iwocpt * 10
            
            if (isoc >= 0) then
               write(title1,'(2i1,a1)') isocpt,isocpo,'N'
            else
               write(title1,'(2i1,a1)') isocpt,isocpo,'S'
            end if
            
            if (iwoc >= 0) then
               write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'E'
            else
               write(title2,'(3i1,a1)') iwocph,iwocpt,iwocpo,'W'
            end if
            
            title3 = trim(ofn)//title1//title2//'.h5'
            inquire(file=trim(title3),exist=l1,opened=l2)
            
            !----- Read file. -------------------------------------------------------------!
            if (l1) then
               write (unit=*,fmt='(a,1x,a,a)') ' -> Getting file:',trim(title3),'...'
               
               call shdf5_open_f(trim(title3),'R')
               
               ndims = 2
               idims(1) = nio
               idims(2) = njo
               
               select case (trim(iaction))
               case ('leaf_class')
                  call shdf5_irec_f(ndims,idims,'oge2',ivara=idato)
               case ('soil_text')
                  ndims = 3
                  idims(1) = nio
                  idims(2) = njo
                  idims(3) = nzg
                  call shdf5_irec_f(ndims,idims,'fao',ivara=idato2)
                  idato(:,:)=idato2(:,:,1)
               case ('soil_col')
                  call shdf5_irec_f(ndims,idims,'colour',ivara=idato)
               case default
                  call fatal_error('Incorrect action specified in leaf_database'           &
                                  ,'leaf_database','leaf_database.f90')
               end select

               call shdf5_close_f()
            else
               write (unit=*,fmt='(a,1x,a,a)') ' -> File:',trim(title3)                    &
                                              ,'doesn''t exist.  Using default values...'
               select case (trim(iaction))
               case ('leaf_class')
                  idato(:,:) = 0
               case ('soil_text')
                 do z = 1,nzg 
                    idato2(:,:,z) = nslcon(z)
                 end do
                    idato(:,:) = nslcon(nzg-1)
               case ('soil_col')
                  idato(:,:) = isoilcol
               case default
                  call fatal_error('Incorrect action specified in leaf_database'           &
                                  ,'leaf_database','leaf_database.f90')
               end select
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Loop through all of the corner points that are inside this file.         !
            !------------------------------------------------------------------------------!
            do cid = 1,nump(ifile,jfile)

               !------ Get the global index of this point. --------------------------------!
               ind = cnr_ind(ifile,jfile,cid)
               !---------------------------------------------------------------------------!


               !------ Get the polygon index of this corner. ------------------------------!
               ilandsea = cnr_ipoly(ind)
               !---------------------------------------------------------------------------!


               if (cnr_i1(ind) < cnr_i2(ind)) then
                  i1 = cnr_i1(ind)
                  i2 = cnr_i2(ind)
               else
                  i1 = cnr_i2(ind)
                  i2 = cnr_i1(ind)
               end if
               
               if (cnr_j1(ind) < cnr_j2(ind)) then
                  j1 = cnr_j1(ind)
                  j2 = cnr_j2(ind)
               else
                  j1 = cnr_j2(ind)
                  j2 = cnr_j1(ind)
               endif

               !----- Now walk through the point spaces for soil or vegetation class. -----!
               select case (trim(iaction))
               case ('leaf_class')
                  do j=j1,j2
                     do i=i1,i2
                        call ed_datp_datq(idato(i,j),dq)
                        class_count(dq,ilandsea) = class_count(dq,ilandsea) + 1
                     end do
                  end do
                     
               case ('soil_text','soil_col')
                  do j=j1,j2
                     do i=i1,i2
                        dq = idato(i,j)
                        class_count(dq,ilandsea) = class_count(dq,ilandsea) + 1
                     end do
                  end do

               end select
               !---------------------------------------------------------------------------!
            end do
         end if
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Decide what to do based on whether we are loading land classes or soil classes.  !
   !---------------------------------------------------------------------------------------!
   select case (trim(iaction))
   case ('leaf_class')
      !------------------------------------------------------------------------------------!
      !     Leaf class.  We really don't care about the leaf class, we just want to know   !
      ! the percentage that is not water.                                                  !
      !------------------------------------------------------------------------------------!
      do ilandsea = 1,nlandsea

         total_count          = real(sum(class_count(0:ed_nvtyp,ilandsea)))
         pctout  (:,ilandsea) = 0.
         pctout  (1,ilandsea) = real(sum(class_count(2:ed_nvtyp,ilandsea))) / total_count
         !----- This won't be used, assign any number. ------------------------------------!
         classout(:,ilandsea,:) = 6
      end do
      !------------------------------------------------------------------------------------!

   case ('soil_text')
      !------------------------------------------------------------------------------------!
      !     We must determine the abundance of each soil classes so that the output will   !
      ! be have the commonest soil type as the first, and the rarest as the last.          !
      !------------------------------------------------------------------------------------!
      do ilandsea = 1,nlandsea
         !---------------------------------------------------------------------------------!
         !     We don't consider soil type 0 as water is determined by other means.  If    !
         ! that causes the total_count variable to be zero, we assume that only one site   !
         ! can exist, with the user-defined default.                                       !
         !---------------------------------------------------------------------------------!
         total_count         = real(sum(class_count(1:ed_nstyp,ilandsea)))
         if (total_count == 0.) then
            pctout(:,ilandsea) = 0.
            pctout(1,ilandsea) = 1.
            !------------------------------------------------------------------------------!
            !    Assign the default class as this may be used if soil and land use maps    !
            ! disagree about whether a place is land or water.                             !
            !------------------------------------------------------------------------------!
            do z = 1,nzg 
               classout(:,ilandsea,z) = nslcon(z)
            end do
            !------------------------------------------------------------------------------!

         else
            !------ Normalise the areas. --------------------------------------------------!
            fraction(1:ed_nstyp) = real(class_count(1:ed_nstyp,ilandsea)) / total_count
            !------------------------------------------------------------------------------!

            !------ Rank the sites (1 = commonest, ed_nstyp = rarest). --------------------!
            call rank_down(ed_nstyp,fraction,rankpct)
            !------------------------------------------------------------------------------!


            !------ Fill in sites with the right order. -----------------------------------!
            do irank = 1,nsite
               !----- Itext is the texture type that has rank irank. ----------------------!
               itext = find_rank(irank,ed_nstyp,rankpct)
               
               classout(irank,ilandsea,:) = itext
               pctout  (irank,ilandsea) = fraction(itext)
            end do
            !------------------------------------------------------------------------------!
         end if
      end do
   end select
   !---------------------------------------------------------------------------------------!



   !----- Free memory so it won't leak it. ------------------------------------------------!
   deallocate(nump       )
   deallocate(idato      )
   deallocate(idato2     )
   deallocate(cnr_ind    )
   deallocate(class_count)
   deallocate(rankpct    )
   deallocate(fraction   )
   !---------------------------------------------------------------------------------------!
  return
end subroutine leaf_database
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine looks for the indices of a given longitude and latitude within an    !
! input file.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine get_file_indices(lat,lon,west,south,nper,ni,nj,io_full,jo_full,io_loc,jo_loc    &
                           ,ifile,jfile)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in)  :: lat
   real   , intent(in)  :: lon
   integer, intent(in)  :: west
   integer, intent(in)  :: south
   integer, intent(in)  :: nper
   integer, intent(in)  :: ni
   integer, intent(in)  :: nj
   integer, intent(out) :: io_full
   integer, intent(out) :: jo_full
   integer, intent(out) :: io_loc
   integer, intent(out) :: jo_loc
   integer, intent(out) :: ifile
   integer, intent(out) :: jfile
   !----- Local variables. ----------------------------------------------------------------!
   real                 :: wlat
   real                 :: wlon
   real                 :: rio_full
   real                 :: rjo_full
   !---------------------------------------------------------------------------------------!



   !----- Set the lat and lon for file lookup. --------------------------------------------!
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






!==========================================================================================!
!==========================================================================================!
!      This subroutine maps the input datp classes to a smaller set datq which represents  !
! the full set of LEAF-2 or LEAF-3 classes for which LSP values are defined.  We don't use !
! the vegetation information for anything regarding ED plant communities: all we want to   !
! know is whether the place is water or land.                                              !
!------------------------------------------------------------------------------------------!
subroutine ed_datp_datq(datp,datq)


   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)       :: datp      ! Input data, the Olson Global Ecosystems dataset
   integer, intent(out)      :: datq      ! Output data, the LEAF-3 classification. 
   !----- Local variable. -----------------------------------------------------------------!
   integer, dimension(0:100) :: oge2leaf3 ! Conversion table
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Olson Global Ecosystems dataset OGE_2 (96 classes) mapped to LEAF-3 classes (see  !
   ! leaf3_document). 97 and 98 are not used, 99 is Goode Homolosine empty data, and 100   !
   ! is missing data. They are all mapped to ocean.                                        !
   !---------------------------------------------------------------------------------------!
   oge2leaf3 = (/                                      0, &  !   0 -   0
                  19,  8,  4,  5,  6,  7,  9,  3, 11, 16, &  !   1 -  10
                  10,  2, 20,  0,  0, 12, 13, 14, 18,  4, &  !  11 -  20
                   4,  4, 14, 14,  6,  6,  4,  7,  7, 15, &  !  21 -  30
                  15,  6,  7,  7, 15, 16, 16, 16, 16,  8, &  !  31 -  40
                   8,  8, 18, 17, 17, 12, 12,  7, 10,  3, &  !  41 -  50
                  10, 10, 11, 14, 18, 18, 18, 18, 13,  6, &  !  51 -  60
                   5,  4, 11, 12,  0,  0,  0,  0,  3,  2, &  !  61 -  70
                   3, 20,  0, 17, 17, 17,  4, 14,  7,  3, &  !  71 -  80
                   3,  3,  3,  3,  3,  3,  8, 12,  7,  6, &  !  81 -  90
                  18, 15, 15, 15, 19,  0,  0,  0,  0,  0  /) !  91 - 100
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert the data.  In the unlikely case that the dataset lies outside the range,  !
   ! re-assign it to ocean.                                                                !
   !---------------------------------------------------------------------------------------!
   select case (datp)
   case (0:95)
      datq = oge2leaf3(datp)
   case default
      datq = 0
   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_datp_datq
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine maps the input datp soil classes to a smaller set datsoil which     !
! represents the full set of LEAF-2 classes for which soil parameter values are defined.   !
! This sub-routine is now deprecated as we converted the soil types in the input files to  !
! have the same classification.                                                            !
!------------------------------------------------------------------------------------------!
subroutine ed_datp_datsoil(datp,datsoil)


   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)       :: datp      ! Input data, the FAO dataset
   integer, intent(out)      :: datsoil   ! Output data, the LEAF-3 classification. 
   !----- Local variable. -----------------------------------------------------------------!
   integer, dimension(0:133) :: fao2leaf3 ! Conversion table
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      (Bob 9/14/2000) This table maps FAO soil units numbered 0 to 132, plus our own   !
   ! missing value designated 133, to the USDA soil textural classes.  FAO classes [0]     !
   ! (ocean), [1, 7, 27, 42, 66, 69, 77, 78, 84, 98, 113] (soils not used in original FAO  !
   ! dataset), [132] (water), and [133] (our missing value) are all mapped to a default    !
   ! class of sandy clay loam in case they happen to correspond to a land surface area in  !
   ! the landuse dataset that RAMS uses to define land area.  We wrote missing value class !
   ! 133 to the RAMS FAO files  whenever a negative value (which is out of range of        !
   ! defined values) was found in the original FAO dataset, which occurred in about 2.6%   !
   ! of the pixels.  For the remaining FAO classes, a cross reference table to Zobler soil !
   ! texture classes that was provided, plus our own cross referencing table from Zobler   !
   ! to USDA classes listed below, provides the mapping from FAO to USDA.  In this         !
   ! mapping, we use only organic USDA classes and omit nonorganic classes since most      !
   ! soils do contain organic matter, and organic content information is not provided in   !
   ! the Zobler classes.                                                                   !
   !                                                                                       !
   !          |--------------------------------+--------------------------------|          !
   !          |          Zobler Class          |           USDA Class           |          !
   !          |--------------------------------+--------------------------------|          !
   !          | 1  Coarse                      |  2  Loamy sand                 |          !
   !          | 2  Medium                      |  4  Silt loam                  |          !
   !          | 3  Fine                        |  8  Clay loam                  |          !
   !          | 4  Coarse-medium               |  3  Sandy loam                 |          !
   !          | 5  Coarse-fine                 |  6  Sandy clay loam            |          !
   !          | 6  Medium-fine                 |  7  Silty clay loam            |          !
   !          | 7  Coarse-medium-fine          |  6  Sandy clay loam            |          !
   !          | 8  Organic matter              |  5  Loam                       |          !
   !          |--------------------------------+--------------------------------|          !
   !          |                                | ... not used:                  |          !
   !          |                                |  1  Sand                       |          !
   !          |                                |  9  Sandy clay                 |          !
   !          |                                | 10  Silty clay                 |          !
   !          |                                | 11  Clay                       |          !
   !          |                                | 12  Peat                       |          !
   !          |--------------------------------+--------------------------------|          !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!

   fao2leaf3 = (/                                      6  &  !   0 -   0
                ,  6,  4,  4,  7,  7,  8,  6,  4,  4,  4  &  !   1 -  10
                ,  7,  4,  4,  4,  8,  4,  8,  4,  4,  8  &  !  11 -  20
                ,  4,  2,  4,  4,  4,  4,  6,  8,  8,  8  &  !  21 -  30
                ,  4,  8,  8,  2,  6,  4,  7,  4,  4,  3  &  !  31 -  40
                ,  4,  6,  7,  4,  4,  4,  4,  4,  4,  4  &  !  41 -  50
                ,  4,  4,  4,  4,  4,  4,  2,  4,  4,  2  &  !  51 -  60
                ,  4,  3,  4,  2,  7,  6,  4,  4,  6,  8  &  !  61 -  70
                ,  8,  7,  2,  5,  4,  5,  6,  6,  4,  2  &  !  71 -  80
                ,  2,  2,  4,  6,  2,  2,  2,  2,  2,  4  &  !  81 -  90
                ,  2,  2,  2,  4,  2,  4,  3,  6,  2,  7  &  !  91 - 100
                ,  4,  4,  4,  8,  8,  8,  3,  7,  4,  4  &  ! 101 - 110
                ,  4,  3,  6,  4,  2,  4,  4,  4,  2,  2  &  ! 111 - 120
                ,  2,  4,  6,  4,  4,  7,  7,  6,  3,  2  &  ! 121 - 130
                ,  2,  6,  6                              /) ! 131 - 133

   datsoil = fao2leaf3(datp)

   return
end subroutine ed_datp_datsoil
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine reads in the soil depth files.                                      !
!------------------------------------------------------------------------------------------!
subroutine read_soil_depth()
   use soil_coms, only : soildepth_db   & ! intent(in)
                       , isoildepthflg  & ! intent(in)
                       , nlon_lyr       & ! intent(in)
                       , nlat_lyr       & ! intent(in)
                       , layer_index    & ! intent(in)
                       , slz            ! ! intent(in)
   use grid_coms, only : nzg            ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real, dimension(nlat_lyr,nlon_lyr) :: soil_depth
   logical                            :: found
   real                               :: rdum
   integer                            :: ilat
   integer                            :: k
   !---------------------------------------------------------------------------------------!


   !----- Allocate the scratch array containing the soil depth index layer. ---------------!
   if (allocated(layer_index)) deallocate(layer_index)
   allocate(layer_index(nlat_lyr,nlon_lyr))
   
   !---------------------------------------------------------------------------------------! 
   !     If the user doesn't want to use the soil depth database, set all layer indices to !
   ! 1 and quit.  This will make all soil layers active across all polygons.               !
   !---------------------------------------------------------------------------------------! 
   if (isoildepthflg /= 1) then 
      layer_index(:,:) = 1
      return
   end if
   !---------------------------------------------------------------------------------------! 


   !---------------------------------------------------------------------------------------! 
   !    If we reach this point, the user provided soil depth data.  Initialise the array   !
   ! by assigning just the top 2 layers, which is the minimum required. 
   !---------------------------------------------------------------------------------------! 
   layer_index=nzg-1
   !---------------------------------------------------------------------------------------! 



   !----- Open the file, provided it exists. If it doesn't, abort the run.-----------------!
   inquire (file=trim(soildepth_db),exist=found)
   if(.not.found) then
      write(unit=*,fmt='(a)') 'Your ED2IN has isoildepthflg=1, however the file'
      write(unit=*,fmt='(a)') trim(soildepth_db)//' doesn''t exist.'
      call fatal_error('File specified in SOILDEPTH_DB not found!'                         &
                      ,'read_soil_depth','leaf_database.f90')
   end if
   open (unit=12,file=trim(soildepth_db),form='formatted',action='read',status='old')
   !---------------------------------------------------------------------------------------! 

   !---------------------------------------------------------------------------------------! 
   !      Read in the database.  At this point we are assuming the grid resolution of the  !
   ! input data.  We may need to re-visit this if we ever find another soil depth dataset. !
   ! For the current data, the first latitude bin is 90N-89N, and the first longitude bin  !
   ! is 180W-179W.                                                                         !
   !---------------------------------------------------------------------------------------! 
   do ilat=1,nlat_lyr
      read(unit=12,fmt=*) rdum,soil_depth(ilat,1:nlon_lyr)
   end do
   !---------------------------------------------------------------------------------------! 


   !------ Close the file, we got all that we need from it. -------------------------------!
   close(unit=12,status='keep')
   !---------------------------------------------------------------------------------------! 



   !----- Convert it to metres and use the same notation as ED (negative values) ----------!
   soil_depth(:,:) = - soil_depth(:,:) * 0.01
   !---------------------------------------------------------------------------------------! 


   !---------------------------------------------------------------------------------------!
   !     Do something in case the atmospheric land cell happens to fall in an FAO water    !
   ! cell.                                                                                 !
   !---------------------------------------------------------------------------------------!
   where (soil_depth(:,:) == 0.0) soil_depth(:,:) = slz(1)
   !---------------------------------------------------------------------------------------!


   !----- Loop through the layers, and find points in which the soil is shallower. --------!
   do k = nzg-2,1,-1
      where (soil_depth(:,:) < slz(k+1)) layer_index(:,:) = k
   end do
   !---------------------------------------------------------------------------------------!
   
   return
end subroutine read_soil_depth
!==========================================================================================!
!==========================================================================================!
