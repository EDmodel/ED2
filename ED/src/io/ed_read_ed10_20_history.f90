!==========================================================================================!
!==========================================================================================!
!     This subroutine reads ED-1.0 and ED-2.0 history files (aka restart), for a polygon   !
! level run with previous ecological state, but with standard initial thermodynamic        !
! conditions.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine read_ed10_ed20_history_file


   use ed_max_dims    , only : n_pft               & ! intent(in)
                             , huge_patch          & ! intent(in)
                             , huge_cohort         & ! intent(in)
                             , max_water           & ! intent(in)
                             , str_len             & ! intent(in)
                             , maxfiles            & ! intent(in)
                             , maxlist             ! ! intent(in)
   use pft_coms       , only : SLA                 & ! intent(in)
                             , q                   & ! intent(in)
                             , qsw                 & ! intent(in)
                             , hgt_min             & ! intent(in)
                             , min_dbh             & ! intent(in)
                             , min_bdead           & ! intent(in)
                             , is_grass            & ! intent(in)
                             , include_pft         & ! intent(in)
                             , include_pft_ag      & ! intent(in)
                             , pft_1st_check       & ! intent(in)
                             , agf_bs              & ! intent(in)
                             , include_these_pft   ! ! intent(in)
   use ed_misc_coms   , only : sfilin              & ! intent(in)
                             , ied_init_mode       ! ! intent(in)
   use mem_polygons   , only : grid_res            & ! intent(in)
                             , edres               ! ! intent(in)
   use consts_coms    , only : pio180              & ! intent(in)
                             , pio4                ! ! intent(in)
   use ed_misc_coms   , only : use_target_year     & ! intent(in)
                             , restart_target_year ! ! intent(in)
   use ed_state_vars  , only : polygontype         & ! variable type
                             , sitetype            & ! variable type
                             , patchtype           & ! variable type
                             , edtype              & ! variable type
                             , edgrid_g            & ! variable type
                             , allocate_sitetype   & ! subroutine
                             , allocate_patchtype  ! ! subroutine
   use grid_coms      , only : ngrids              ! ! intent(in)
   use allometry      , only : bd2dbh              & ! function
                             , dbh2h               & ! function
                             , dbh2bd              & ! function
                             , size2bl             & ! function
                             , ed_biomass          & ! function
                             , area_indices        ! ! subroutine
   use fuse_fiss_utils, only : sort_cohorts        & ! subroutine
                             , sort_patches        ! ! subroutine
   use disturb_coms   , only : ianth_disturb       ! ! intent(in)
   implicit none

   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter :: min_area = 1.d-7           ! Minimum acceptable area.
   real(kind=8), parameter :: min_ok   = 1.d-20          ! Minimum acceptable value for 
                                                         !    any restart variable.
   logical     , parameter :: harvard_override = .false. ! Overwrite some initial values
                                                         !    for a specific Harvard run?
                                                         !    (this should never be used
                                                         !     unless you are 110% sure of
                                                         !     what you are doing.)
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)          , pointer                        :: cgrid
   type(polygontype)     , pointer                        :: cpoly
   type(sitetype)        , pointer                        :: csite
   type(patchtype)       , pointer                        :: cpatch
   character(len=str_len), dimension(maxlist)             :: full_list
   character(len=str_len), dimension(maxfiles)            :: site_list
   character(len=str_len), dimension(maxfiles)            :: pss_list
   character(len=str_len), dimension(maxfiles)            :: css_list
   character(len=str_len), dimension(huge_patch)          :: pname
   character(len=str_len)                                 :: pss_name
   character(len=str_len)                                 :: css_name
   character(len=str_len)                                 :: cdum
   character(len=str_len), dimension(huge_cohort)         :: cname
   character(len=str_len), dimension(huge_cohort)         :: cpname
   integer               , dimension(n_pft)               :: include_pft_ep
   integer               , dimension(huge_patch)          :: trk
   integer               , dimension(huge_patch)          :: sitenum
   integer               , dimension(huge_cohort)         :: leaves_on
   integer               , dimension(huge_cohort)         :: ipft
   integer                                                :: year
   integer                                                :: pft
   integer                                                :: igr
   integer                                                :: ipy
   integer                                                :: isi
   integer                                                :: ipa
   integer                                                :: ico
   integer                                                :: ip
   integer                                                :: ip2
   integer                                                :: ic
   integer                                                :: ic2
   integer                                                :: nwater
   integer                                                :: ierr
   integer                                                :: nf
   integer                                                :: nflist
   integer                                                :: nflsite
   integer                                                :: nflpss
   integer                                                :: nflcss
   integer                                                :: nclosest
   integer                                                :: ncohorts
   integer                                                :: npatchco
   integer                                                :: npatches
   integer                                                :: nsitepat
   integer                                                :: npatch2
   integer                                                :: nw
   logical              , dimension(huge_cohort)          :: add_this_cohort
   logical                                                :: renumber_pfts
   logical                                                :: site_match
   real                 , dimension(max_water)            :: depth
   real                 , dimension(huge_patch)           :: time
   real                 , dimension(huge_patch)           :: age
   real                 , dimension(huge_patch)           :: area
   real                 , dimension(huge_patch)           :: fsc
   real                 , dimension(huge_patch)           :: stsc
   real                 , dimension(huge_patch)           :: stsl
   real                 , dimension(huge_patch)           :: ssc
   real                 , dimension(huge_patch)           :: psc
   real                 , dimension(huge_patch)           :: msn
   real                 , dimension(huge_patch)           :: fsn
   real                 , dimension(max_water,huge_patch) :: water
   real                 , dimension(12,huge_cohort)       :: cb
   real                 , dimension(12,huge_cohort)       :: cb_max
   real                 , dimension(huge_cohort)          :: balive
   real                 , dimension(huge_cohort)          :: avgRg
   real                 , dimension(huge_cohort)          :: bdead
   real                 , dimension(huge_cohort)          :: nplant
   real                 , dimension(huge_cohort)          :: hite
   real                 , dimension(huge_cohort)          :: dbh
   real                 , dimension(huge_cohort)          :: ctime
   real                 , dimension(maxfiles)             :: slon_list,slat_list
   real                 , dimension(maxfiles)             :: plon_list,plat_list
   real                 , dimension(maxfiles)             :: clon_list,clat_list
   real                 , dimension(maxfiles)             :: file_pdist,file_cdist
   real                                                   :: dummy
   real                                                   :: area_tot
   real                                                   :: area_sum
   real(kind=8)         , dimension(max_water)            :: dwater
   real(kind=8)                                           :: dage
   real(kind=8)                                           :: darea
   real(kind=8)                                           :: dfsc
   real(kind=8)                                           :: dstsc
   real(kind=8)                                           :: dstsl
   real(kind=8)                                           :: dssc
   real(kind=8)                                           :: dpsc
   real(kind=8)                                           :: dmsn
   real(kind=8)                                           :: dfsn
   !----- External function. --------------------------------------------------------------!
   real                 , external                        :: sngloff
   real                 , external                        :: dist_gc
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !      Now we loop over all all grids and polygons, and fill them with patches and      !
   ! cohorts from the closest polygon.                                                     !
   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids

      !----- Retrieve all files with the specified prefix. --------------------------------!
      call ed_filelist(full_list,sfilin(igr),nflist)

      !----- Retrieve LON/LAT information for sites ---------------------------------------!
      select case (ied_init_mode)
      case (3)
         renumber_pfts = .false.
         call ed1_fileinfo('.site',nflist,full_list,nflsite,site_list,slon_list,slat_list)
      case (6)
         renumber_pfts = .false.
      case default
         renumber_pfts = .true.
      end select

      !----- Retrieve LON/LAT information for patches and cohorts -------------------------!
      call ed1_fileinfo('.pss',nflist,full_list,nflpss,pss_list,plon_list,plat_list)
      call ed1_fileinfo('.css',nflist,full_list,nflcss,css_list,clon_list,clat_list)

      cgrid => edgrid_g(igr)

      polyloop: do ipy = 1,cgrid%npolygons
         
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !     Especial set up of soil types for Harvard forest when running some specific !
         ! test.                                                                           !
         !---------------------------------------------------------------------------------!
         if (harvard_override) then
            cpoly%lsl(1)            = 1
            cpoly%ntext_soil(1,isi) = 2
            cpoly%ntext_soil(2,isi) = 2
            cpoly%ntext_soil(3,isi) = 3
            cpoly%ntext_soil(4,isi) = 3
            cpoly%ncol_soil(isi)    = 10
         end if


         !---------------------------------------------------------------------------------!
         !    Initialise the distances as very large numbers, so if we don't fill all the  !
         ! patches and cohorts, we will not going to take non-sense as a valid polygon.    !
         !---------------------------------------------------------------------------------!
         file_pdist(:) = 1.e20
         file_cdist(:) = 1.e20

         !---------------------------------------------------------------------------------!
         !    Compute the distances between every polygon in the restart files and the     !
         ! current polygon.                                                                !
         !---------------------------------------------------------------------------------!
         do nf=1,nflpss
            file_pdist(nf) = dist_gc(cgrid%lon(ipy),plon_list(nf)                          &
                                    ,cgrid%lat(ipy),plat_list(nf) )
         end do
         do nf=1,nflcss
            file_cdist(nf) = dist_gc(cgrid%lon(ipy),clon_list(nf)                          &
                                    ,cgrid%lat(ipy),clat_list(nf) )
         end do
         !---------------------------------------------------------------------------------!


         find_nonwater: do nf=1,nflpss
            !------------------------------------------------------------------------------!
            !     Find the file that is the closest to the current polygon, based on the   !
            ! distance vector.                                                             !
            !------------------------------------------------------------------------------!
            nclosest = minloc(file_pdist,dim=1)
            pss_name = trim(pss_list(nclosest))
            write (unit=*,fmt='(2a)') 'Using patch file: ',trim(pss_name)

            !------------------------------------------------------------------------------!
            !    Open the patch file and read in all patches.                              !
            !------------------------------------------------------------------------------!
            open(unit=12,file=trim(pss_name),form='formatted',status='old',action='read')
            read(unit=12,fmt='(a4)')  cdum
            
            !----- Read the other information from the header (if there is any...). -------!
            nwater = 1
            if (ied_init_mode == 1 ) then
               read (unit=12,fmt=*) cdum,nwater
               read (unit=12,fmt=*) cdum,depth(1:nwater)
               read (unit=12,fmt=*)
            end if
            
            !------------------------------------------------------------------------------!
            !     Now we loop over all patches and decide whether they should be included  !
            ! or not.                                                                      !
            !------------------------------------------------------------------------------!
            ip      = 1
            sitenum = 0
            count_patches: do
               
               !---------------------------------------------------------------------------!
               !     We must check whether we are not exceeding the maximum number of      !
               ! patches that we can read.                                                 !
               !---------------------------------------------------------------------------!
               if (ip > huge_patch) then
                  write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(pss_name)
                  write (unit=*,fmt='(a)')       ' Number of patches is > HUGE_PATCH...'
                  write (unit=*,fmt='(a,1x,i7)') ' HUGE_PATCH:',huge_patch
                  write (unit=*,fmt='(a)')       ' Increase HUGE_PATCH to read this...'
                  call fatal_error('Too many patches to be read...'                        &
                                  ,'read_ed10_ed20_history_file','ed_history_io.f90')
               end if

               !---------------------------------------------------------------------------!
               !    If we can still add new patches, read the next one, according to the   !
               ! input format type.                                                        !
               !---------------------------------------------------------------------------!
               select case (ied_init_mode)
               case (1)
                  !----- ED-1.0 file. -----------------------------------------------------!
                  read(unit=12,fmt=*,iostat=ierr) time(ip),pname(ip),trk(ip),dage,darea    &
                                                 ,dfsc,dstsc,dstsl,dssc,dpsc,dmsn,dfsn     &
                                                 ,dwater(1:nwater)

                  !------------------------------------------------------------------------!
                  !     Check whether the file has hit the end, and if so, leave the loop. !
                  !------------------------------------------------------------------------!
                  if (ierr /= 0) exit count_patches

                  !----- Copy the double-precision scratch variables to the arrays. -------!
                  area  (ip) = sngloff(darea      ,min_area)
                  age   (ip) = sngloff(dage       ,min_ok  )
                  fsc   (ip) = sngloff(dfsc       ,min_ok  )
                  stsc  (ip) = sngloff(dstsc      ,min_ok  )
                  stsl  (ip) = sngloff(dstsl      ,min_ok  )
                  ssc   (ip) = sngloff(dssc       ,min_ok  )
                  psc   (ip) = sngloff(dpsc       ,min_ok  )
                  msn   (ip) = sngloff(dmsn       ,min_ok  )
                  fsn   (ip) = sngloff(dfsn       ,min_ok  )
                  do nw=1,nwater
                     water(nw,ip) = sngloff(dwater(nw) ,min_ok  )
                  end do

               case (2,6)
                  !----- Standard ED-2.0 file. --------------------------------------------!
                  read(unit=12,fmt=*,iostat=ierr) time(ip),pname(ip),trk(ip),dage,darea    &
                                                 ,dwater(1),dfsc,dstsc,dstsl,dssc,dummy    &
                                                 ,dmsn,dfsn

                  !------------------------------------------------------------------------!
                  !     Check whether the file has hit the end, and if so, leave the loop. !
                  !------------------------------------------------------------------------!
                  if(ierr /= 0)exit count_patches

                  !----- Copy the double-precision scratch variables to the arrays. -------!
                  area   (ip) = sngloff(darea    ,min_area)
                  age    (ip) = sngloff(dage     ,min_ok  )
                  fsc    (ip) = sngloff(dfsc     ,min_ok  )
                  stsc   (ip) = sngloff(dstsc    ,min_ok  )
                  stsl   (ip) = sngloff(dstsl    ,min_ok  )
                  ssc    (ip) = sngloff(dssc     ,min_ok  )
                  msn    (ip) = sngloff(dmsn     ,min_ok  )
                  fsn    (ip) = sngloff(dfsn     ,min_ok  )
                  water(1,ip) = sngloff(dwater(1),min_ok  )
                  
               case (3)
                  !----- ED-2.0 file, with site information. ------------------------------!
                  read(unit=12,fmt=*,iostat=ierr) sitenum(ip),time(ip),pname(ip),trk(ip)   &
                                                 ,age(ip),darea,water(1,ip),fsc(ip)        &
                                                 ,stsc(ip),stsl(ip),ssc(ip),psc(ip)        &
                                                 ,msn(ip),fsn(ip)

                  !------------------------------------------------------------------------!
                  !     Check whether the file has hit the end, and if so, leave the loop. !
                  !------------------------------------------------------------------------!
                  if (ierr /= 0) exit count_patches


                  area(ip)=sngloff(darea, min_area)

                  site_match = .false.
                  do isi = 1,cpoly%nsites
                     if (sitenum(ip).eq.cpoly%sitenum(isi)) then
                        site_match = .true.
                     end if
                  end do
                  if (.not.site_match) then
                     write (unit=*,fmt='(a,1x,i6)') ' ISI           = ',isi
                     write (unit=*,fmt='(a,1x,i6)') ' CPOLY%NSITES  = ',cpoly%nsites
                     write (unit=*,fmt='(a,1x,i6)') ' IP            = ',ip
                     write (unit=*,fmt='(a,1x,i6)') ' CPOLY%SITENUM = ',cpoly%sitenum
                     write (unit=*,fmt='(a,1x,a)')  ' Patch file    = ',trim(pss_name)
                     write (unit=*,fmt='(a,1x,i5,1x,a)')                                   &
                                                    ' Site number', sitenum,'not found!'
                     call fatal_error('Error reading from patch file'                      &
                                     ,'read_ed10_ed20_history_file','ed_history_io.F90')
                  end if
               case default !Nearly bare ground
                  exit count_patches
               end select
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      By adding ip only if the area is above minimum, we avoid including   !
               ! patches that are tiny since they will be overwritten by the next patch.   !
               !---------------------------------------------------------------------------!
               if (area(ip) > min_area) ip = ip + 1
               !---------------------------------------------------------------------------!

            end do count_patches
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Here we determine the number of patches.  Unless there are no patches,  !
            ! (this may happen if the file used was all water.  In that case, then we      !
            ! should use the next closest file.                                            !
            !------------------------------------------------------------------------------!
            npatches = max(ip-1,0)

            close(unit=12,status='keep')

            if (npatches > 0) then

               !---------------------------------------------------------------------------!
               !      We have found a suitable file, we will leave the loop, after we find !
               ! the name of the corresponding cohort file.                                !
               !---------------------------------------------------------------------------!
               nclosest = minloc(abs(file_pdist(nclosest)-file_cdist),dim=1)
               css_name = trim(css_list(nclosest))
               write (unit=*,fmt='(2a)') 'Using cohort file: ',trim(css_name)
               exit find_nonwater
            else
               !----- The closest file was no good, so we make it far away for now. -------!
               file_pdist(nclosest) = 1.e20 
            end if
         end do find_nonwater
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Allocate the patches.                                                       !
         !---------------------------------------------------------------------------------!
         select case (ied_init_mode)
         case (3) !---- Multiple-site run, find the number of patches in each site. -------!

            !----- Loop over all sites. ---------------------------------------------------!
            do isi = 1,cpoly%nsites

               npatch2 = 0

               do ip=1,npatches
                  if(sitenum(ip) == cpoly%sitenum(isi)) npatch2 = npatch2 + 1
               end do

               csite => cpoly%site(isi)
               
               !----- Allocate the patches in this site. ----------------------------------!
               call allocate_sitetype(csite,npatch2)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Go through the patch list and fill the sites' patches with data.      !
               !---------------------------------------------------------------------------!
               ip2 = 0
               do ip=1,npatches

                  if (sitenum(ip) == cpoly%sitenum(isi)) then
                     ip2 = ip2 + 1

                     !----- Translate the land use type. ----------------------------------!
                     select case (trk(ip))
                     case (0)
                        !----- Agriculture. -----------------------------------------------!
                        csite%dist_type (ip2) = 1
                        !------------------------------------------------------------------!
                     case (1)
                        !------------------------------------------------------------------!
                        !     Secondary forest.  This was ambiguous in previous versions,  !
                        ! so we decide which type to point this patch based on whether     !
                        ! anthropogenic disturbance is turned on or not.  Using flags 5 or !
                        ! 6 is unambiguous.                                                !
                        !------------------------------------------------------------------!
                        select case (ianth_disturb)
                        case (0)
                           !----- No anthropogenic, assume abandoned lands. ---------------!
                           csite%dist_type (ip2) = 5
                           !---------------------------------------------------------------!
                        case (1)
                           !----- Anthropogenic, assume logging. --------------------------!
                           csite%dist_type (ip2) = 6
                           !---------------------------------------------------------------!
                        end select
                        !------------------------------------------------------------------!
                     case (2)
                        !----- Primary forest. --------------------------------------------!
                        csite%dist_type (ip2) = 3
                        !------------------------------------------------------------------!
                     case (3)
                        !----- Forest plantation. -----------------------------------------!
                        csite%dist_type (ip2) = 2
                        !------------------------------------------------------------------!
                     case (4)
                        !----- Burnt patch. -----------------------------------------------!
                        csite%dist_type (ip2) = 4
                        !------------------------------------------------------------------!
                     case (5)
                        !----- Abandoned land (secondary growth). -------------------------!
                        csite%dist_type (ip2) = 5
                        !------------------------------------------------------------------!
                     case (6)
                        !----- Logged forest. ---------------------------------------------!
                        csite%dist_type (ip2) = 6
                        !------------------------------------------------------------------!
                     end select
                     !---------------------------------------------------------------------!

                     csite%age               (ip2) = age (ip)
                     csite%area              (ip2) = area(ip)
                     csite%fast_soil_C       (ip2) = fsc (ip)
                     csite%slow_soil_C       (ip2) = ssc (ip)
                     csite%structural_soil_C (ip2) = stsc(ip)
                     csite%structural_soil_L (ip2) = stsl(ip)
                     csite%mineralized_soil_N(ip2) = msn (ip)
                     csite%fast_soil_N       (ip2) = fsn (ip)
                     csite%pname             (ip2) = trim(pname(ip))
                     csite%sum_dgd           (ip2) = 0.0
                     csite%sum_chd           (ip2) = 0.0
                     csite%cohort_count      (ip2) = 0
                  end if
               end do
               !---- Initialize the cohort counts per patch. ------------------------------!
               csite%cohort_count(:) = 0
            end do

         case default
            !------------------------------------------------------------------------------!
            !      We don't have information on sites, assume all sites to be the same.    !
            !------------------------------------------------------------------------------!
            do isi = 1, cpoly%nsites
               csite => cpoly%site(isi)

               !----- Allocate the patches in this site. ----------------------------------!
               call allocate_sitetype(csite,npatches)

               do ip=1,npatches

                  !----- Translate the land use type. -------------------------------------!
                  select case (trk(ip))
                  case (0)
                     !----- Agriculture. --------------------------------------------------!
                     csite%dist_type (ip) = 1
                     !---------------------------------------------------------------------!
                  case (1)
                     !---------------------------------------------------------------------!
                     !     Secondary forest.  This was ambiguous in previous versions,     !
                     ! so we decide which type to point this patch based on whether        !
                     ! anthropogenic disturbance is turned on or not.  Using flags 5 or    !
                     ! 6 is unambiguous.                                                   !
                     !---------------------------------------------------------------------!
                     select case (ianth_disturb)
                     case (0)
                        !----- No anthropogenic, assume abandoned lands. ------------------!
                        csite%dist_type (ip) = 5
                        !------------------------------------------------------------------!
                     case (1)
                        !----- Anthropogenic, assume logging. -----------------------------!
                        csite%dist_type (ip) = 6
                        !------------------------------------------------------------------!
                     end select
                     !---------------------------------------------------------------------!
                  case (2)
                     !----- Primary forest. -----------------------------------------------!
                     csite%dist_type (ip) = 3
                     !---------------------------------------------------------------------!
                  case (3)
                     !----- Forest plantation. --------------------------------------------!
                     csite%dist_type (ip) = 2
                     !---------------------------------------------------------------------!
                  case (4)
                     !----- Burnt patch. --------------------------------------------------!
                     csite%dist_type (ip) = 4
                     !---------------------------------------------------------------------!
                  case (5)
                     !----- Abandoned land (secondary growth). ----------------------------!
                     csite%dist_type (ip) = 5
                     !---------------------------------------------------------------------!
                  case (6)
                     !----- Logged forest. ------------------------------------------------!
                     csite%dist_type (ip) = 6
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!

                  csite%age               (ip) = age (ip)
                  csite%area              (ip) = area(ip)
                  csite%fast_soil_C       (ip) = fsc (ip)
                  csite%slow_soil_C       (ip) = ssc (ip)
                  csite%structural_soil_C (ip) = stsc(ip)
                  csite%structural_soil_L (ip) = stsl(ip)
                  csite%mineralized_soil_N(ip) = msn (ip)
                  csite%fast_soil_N       (ip) = fsn (ip)
                  csite%pname             (ip) = trim(pname(ip))
                  csite%sum_dgd           (ip) = 0.0
                  csite%sum_chd           (ip) = 0.0
                  csite%cohort_count      (ip) = 0
               end do
               !---------------------------------------------------------------------------!

               !---- Initialize the cohort counts per patch. ------------------------------!
               csite%cohort_count(:) = 0
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

         close(unit=12,status='keep')

         !---------------------------------------------------------------------------------!
         !    Open the cohort file and read in all cohorts.                                !
         !---------------------------------------------------------------------------------!
         open(unit=12,file=trim(css_name),form='formatted',status='old')
         read(unit=12,fmt='(a4)')  cdum

         if (ied_init_mode == 1) then
            read(unit=12,fmt=*) !---- Skip second line. -----------------------------------!
         end if
         ic = 0

         !---------------------------------------------------------------------------------!
         !     Now we loop over all patches and decide whether they should be included     !
         ! or not.                                                                         !
         !---------------------------------------------------------------------------------!
         read_cohorts: do

            ic = ic + 1
            add_this_cohort(ic) = .true.
            
            !------------------------------------------------------------------------------!
            !     We must check whether we are not exceeding the maximum number of patches !
            ! that we can read.                                                            !
            !------------------------------------------------------------------------------!
            if (ic > huge_cohort) then
               write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(css_name)
               write (unit=*,fmt='(a)')       ' Number of cohorts is > HUGE_COHORT...'
               write (unit=*,fmt='(a,1x,i7)') ' HUGE_COHORT:',huge_cohort
               write (unit=*,fmt='(a)')       ' Increase HUGE_COHORT to read this...'
               call fatal_error('Too many cohorts to be read...'                           &
                               ,'read_ed10_ed20_history_file','ed_history_io.f90')
            end if

            !------------------------------------------------------------------------------!
            !    If we can still add new cohorts, read the next one, according to the      !
            ! input format type.                                                           !
            !------------------------------------------------------------------------------!
            select case (ied_init_mode)
            case (1)
               !----- ED-1.0 file. --------------------------------------------------------!
               read(unit=12,fmt=*,iostat=ierr)  ctime(ic),cpname(ic),cname(ic),dbh(ic)     &
                                               ,hite(ic),ipft(ic),nplant(ic),bdead(ic)     &
                                               ,balive(ic),avgRg(ic),leaves_on(ic)         &
                                               ,cb(1:12,ic),cb_max(1:12,ic)

               !---------------------------------------------------------------------------!
               !     Check whether the file has hit the end, and if so, leave the loop.    !
               !---------------------------------------------------------------------------!
               if(ierr /= 0) exit read_cohorts  

            case (2,3,6)
               !----- ED-2.0 file. --------------------------------------------------------!
               read(unit=12,fmt=*,iostat=ierr) ctime(ic),cpname(ic),cname(ic),dbh(ic)      &
                                              ,hite(ic),ipft(ic),nplant(ic),bdead(ic)      &
                                              ,balive(ic),avgRg(ic)
               !---------------------------------------------------------------------------!
               !     Check whether the file has hit the end, and if so, leave the loop.    !
               !---------------------------------------------------------------------------!
               if(ierr /= 0)exit read_cohorts


               !----- No carbon balance information.  Assign 1. ---------------------------!
               cb(1:12,ic)     = 1.0
               cb_max(1:12,ic) = 1.0
            end select

            !------------------------------------------------------------------------------!
            !     The PFT classes has changed between different ED versions, here we       !
            ! standardise this.                                                            !
            !------------------------------------------------------------------------------!
            if(renumber_pfts) then
               if (ipft(ic) < 100) then
                  ipft(ic) = ipft(ic) + 1
                  if(ipft(ic) >= 5) ipft(ic) = ipft(ic) - 3
               else
                  ipft(ic) = ipft(ic) - 100
               end if
            end if

            !----- Check if the year matches.  If not, we will ignore this cohort. --------!
            year = int(ctime(ic))
            if(use_target_year == 1 .and. year /= restart_target_year) then
               add_this_cohort(ic) = .false. 
            end if

            !----- Remove cohort in case nplant > 0. --------------------------------------!
            if(nplant(ic) < tiny(1.0)) add_this_cohort(ic) = .false.


            !------------------------------------------------------------------------------!
            !     Find site and patch to which this cohort belong, and start counting how  !
            ! many we need to allocate.                                                    !
            !------------------------------------------------------------------------------!
            put_cohort:do isi=1,cpoly%nsites
               csite => cpoly%site(isi)
               do ipa=1,csite%npatches

                  !------------------------------------------------------------------------!
                  !    Here we test whether the PFT of this cohort is expected to be       !
                  ! included.                                                              !
                  !------------------------------------------------------------------------!
                  if (.not. include_pft(ipft(ic))) then
                     !----- This PFT wasn't expected... -----------------------------------!
                     select case (pft_1st_check)
                     case (0)
                        !----- Stop the run. ----------------------------------------------!
                        write (unit=*,fmt='(a,1x,i5,1x,a)')                                &
                             'I found a cohort with PFT=',ipft(ic)                         &
                            ,' and it is not in your include_these_pft...'
                        call fatal_error('Invalid PFT in history file'                     &
                                        ,'read_ed10_ed20_history_file','ed_history_io.f90')

                     case (1)
                        !----- Add the unexpected PFT to the list of possible PFTs. -------!
                        write (unit=*,fmt='(a,1x,i5,1x,a)')                                &
                             'I found a cohort with PFT=',ipft(ic)                         &
                            ,'... Including this PFT in your include_these_pft...'
                        include_pft(ipft(ic))                 = .true.
                        include_these_pft(count(include_pft)) = ipft(ic)
                        call sort_up(include_these_pft,n_pft)
                        if (is_grass(ipft(ic))) include_pft_ag(ipft(ic)) = .true.

                     case (2)
                        !----- Ignore the cohort. -----------------------------------------!
                        write (unit=*,fmt='(a,1x,i5,1x,a)')                                &
                             'I found a cohort with PFT=',ipft(ic),'... Ignoring it...'
                        add_this_cohort(ic) = .false.
                     end select
                  end if
                  
                  if (trim(csite%pname(ipa)) == trim(cpname(ic) ) .and.                    &
                      add_this_cohort(ic)                              ) then
                     csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
                     exit put_cohort
                  end if
               end do
            end do put_cohort
         end do read_cohorts

         !----- Find the total number of cohorts. -----------------------------------------!
         ncohorts = max(ic-1,0)
         
         
         close (unit=12,status='keep')

         loop_sites: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            loop_patches: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               if (csite%cohort_count(ipa) /= 0) then

                  call allocate_patchtype(cpatch,csite%cohort_count(ipa))
                  csite%plant_ag_biomass(ipa) = 0.
                  ic2 = 0
                  do ic = 1,ncohorts

                     if (trim(csite%pname(ipa)) == trim(cpname(ic) ) .and.                 &
                         add_this_cohort(ic)                              ) then

                        ic2 = ic2 + 1
                        cpatch%pft(ic2) = ipft(ic)


                        !------------------------------------------------------------------!
                        !      There were some differences on how a few variables were     !
                        ! defined between ED-1.0 and ED-2.0.  Here we will put everything  !
                        ! in the ED-2.1 standard.                                          !
                        !------------------------------------------------------------------!
                        !----- Plant density, it should be in [plants/m2]. ----------------!
                        if(ied_init_mode == 1) then
                           cpatch%nplant(ic2) = nplant(ic) / (csite%area(ipa) )
                        else
                           cpatch%nplant(ic2) = nplant(ic)
                        end if


                        !------------------------------------------------------------------!
                        !     Select which history file we are using.  We use DBH only     !
                        ! when ied_init_mode is 6 (inventory initialisation, when DBH is   !
                        ! most likely to be the actual measured variable), otherwise we    !
                        ! use BDEAD instead.                                               !
                        !------------------------------------------------------------------!
                        select case(ied_init_mode)
                        case (6)
                           !----- Inventory.  Read DBH and find the other stuff. ----------!
                           cpatch%dbh(ic2)   = max(dbh(ic),min_dbh(ipft(ic)))
                           cpatch%hite(ic2)  = dbh2h(ipft(ic),dbh(ic))
                           cpatch%bdead(ic2) = dbh2bd(dbh(ic),ipft(ic))

                        case default
                           !---------------------------------------------------------------!
                           !    Old ED files.  Check whether bdead is valid.  If it is, we !
                           ! initialise DBH and height from BDEAD, othewise we use DBH     !
                           ! instead.  This is because allometry may be different, and ED  !
                           ! biomass varies less than DBH under different allometric       !
                           ! equations.                                                    !
                           !---------------------------------------------------------------!
                           if (bdead(ic) > 0.0) then
                              cpatch%bdead(ic2) = max(bdead(ic),min_bdead(ipft(ic)))
                              cpatch%dbh(ic2)   = bd2dbh(ipft(ic),bdead(ic))
                              cpatch%hite(ic2)  = dbh2h(ipft(ic),dbh(ic))
                           else
                              cpatch%dbh(ic2)   = dbh(ic)
                              cpatch%hite(ic2)  = dbh2h(ipft(ic),dbh(ic))
                              cpatch%bdead(ic2) = dbh2bd(dbh(ic),ipft(ic))
                           end if
                        end select
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Use allometry to define leaf and the other live biomass      !
                        ! pools.                                                           !
                        !------------------------------------------------------------------!
                        cpatch%bleaf(ic2)     = size2bl(dbh(ic), hite(ic),ipft(ic))
                        cpatch%balive(ic2)    = cpatch%bleaf(ic2) * (1.0 + q(ipft(ic))     &
                                              + qsw(ipft(ic)) * cpatch%hite(ic2))
                        cpatch%broot(ic2)     = cpatch%balive(ic2) * q(ipft(ic))           &
                                              / ( 1.0 + q(ipft(ic)) + qsw(ipft(ic))        &
                                               * cpatch%hite(ic2))
                        cpatch%bsapwooda(ic2) = agf_bs(ipft(ic)) * cpatch%balive(ic2)                &
                                             * qsw(ipft(ic))* cpatch%hite(ic2)             &
                                             / ( 1.0 + q(ipft(ic)) + qsw(ipft(ic))         &
                                               * cpatch%hite(ic2))
                        cpatch%bsapwoodb(ic2) = (1.-agf_bs(ipft(ic))) * cpatch%balive(ic2)           &
                                             * qsw(ipft(ic))* cpatch%hite(ic2)             &
                                             / ( 1.0 + q(ipft(ic)) + qsw(ipft(ic))         &
                                               * cpatch%hite(ic2))


                        !------------------------------------------------------------------!
                        !     Start plants with full phenology, we will take care of       !
                        ! phenology after this sub-routine.                                !
                        !------------------------------------------------------------------!
                        cpatch%phenology_status(ic2) = 0
                        cpatch%bstorage        (ic2) = 0.0
                        !------------------------------------------------------------------!



                        !----- Assign LAI, WAI, and CAI -----------------------------------!
                        call area_indices(cpatch%nplant(ic2),cpatch%bleaf(ic2)             &
                                         ,cpatch%bdead(ic2),cpatch%balive(ic2)             &
                                         ,cpatch%dbh(ic2), cpatch%hite(ic2)                &
                                         ,cpatch%pft(ic2), SLA(cpatch%pft(ic2))            &
                                         ,cpatch%lai(ic2), cpatch%wai(ic2)                 &
                                         ,cpatch%crown_area(ic2),cpatch%bsapwooda(ic2))

                        !------------------------------------------------------------------!
                        !     Initialise the carbon balance.  We ignore the carbon balance !
                        ! even for ED-1.0, the models are so different that there is no    !
                        ! reason to use the stored value.                                  !
                        !------------------------------------------------------------------!
                        cpatch%cb         (1:12,ic2) = 1.0
                        cpatch%cb_lightmax(1:12,ic2) = 1.0
                        cpatch%cb_moistmax(1:12,ic2) = 1.0
                        cpatch%cb_mlmax   (1:12,ic2) = 1.0
                        cpatch%cb         (  13,ic2) = 0.0
                        cpatch%cb_lightmax(  13,ic2) = 0.0
                        cpatch%cb_moistmax(  13,ic2) = 0.0
                        cpatch%cb_mlmax   (  13,ic2) = 0.0
                        !------------------------------------------------------------------!

                        !----- Above ground biomass, use the allometry. -------------------!
                        cpatch%agb(ic2) = ed_biomass(cpatch%bdead(ic2),cpatch%bleaf(ic2)   &
                                                    ,cpatch%bsapwooda(ic2),cpatch%pft(ic2))
                        cpatch%basarea(ic2)  = pio4 * cpatch%dbh(ic2) * cpatch%dbh(ic2)

                        !----- Growth rates, start with zero. -----------------------------!
                        cpatch%dagb_dt  (ic2)  = 0.
                        cpatch%dlnagb_dt(ic2)  = 0.
                        cpatch%dba_dt   (ic2)  = 0.
                        cpatch%dlnba_dt (ic2)  = 0.
                        cpatch%ddbh_dt  (ic2)  = 0.
                        cpatch%dlndbh_dt(ic2)  = 0.
                        
                        !------------------------------------------------------------------!
                        !      Initialise other cohort variables.  Some of them won't be   !
                        ! updated unless the lai exceeds lai_min.                          !
                        !------------------------------------------------------------------!
                        cpatch%fsw(ic2)   = 1.0
                        cpatch%gpp(ic2)   = 0.0
                        cpatch%par_l(ic2) = 0.0

                        !----- Update the patch level above-ground biomass. ---------------!
                        csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)          &
                                                    + cpatch%agb(ic2) * cpatch%nplant(ic2)
                     end if
                  end do
               end if
            end do loop_patches
         end do loop_sites

         close (unit=12,status='keep')

         !----- Initialise all the other site-, patch-, and cohort-level variables. -------!
         do isi = 1,cpoly%nsites
            
            area_sum = 0.0
            ncohorts = 0


            !----- Make sure that the total patch area is 1. ------------------------------!
            csite => cpoly%site(isi)
            area_tot      = sum(csite%area(1:csite%npatches))
            csite%area(:) = csite%area(:)/area_tot

            !----- Find the patch-level LAI, WAI, and CAI. --------------------------------!
            do ipa=1,csite%npatches
               area_sum        = area_sum + csite%area(ipa)

               cpatch => csite%patch(ipa)
               do ico = 1,cpatch%ncohorts
                  ncohorts        = ncohorts + 1
               end do
            end do

            !----- Initialise the cohort variables, then sort them by size. ---------------!
            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)
               do ico = 1,cpatch%ncohorts                 
                  call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
               end do

               !----- Make sure that cohorts are organised from tallest to shortest. ------!
               call sort_cohorts(cpatch)
            end do

            !----- Initialise the patch-level variables. ----------------------------------!
            call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))

            !----- Make sure that patches are organised from oldest to youngest. ----------!
            call sort_patches(csite)
         end do

         !----- Initialise site-level variables. ------------------------------------------!
         call init_ed_site_vars(cpoly,cgrid%lat(ipy))


         !----- Get a diagnostic on the polygon's vegetation. -----------------------------!
         ncohorts = 0

         do isi = 1,cpoly%nsites
            nsitepat = 0
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               npatchco        = 0

               cpatch => csite%patch(ipa)
               do ico = 1,cpatch%ncohorts
                  ncohorts        = ncohorts+1
                  npatchco        = npatchco+1
               end do
               csite%cohort_count(ipa) = npatchco
               nsitepat                = nsitepat + 1
            end do

           ! cpoly%patch_count(isi) = nsitepat
         end do
      end do polyloop

      !----- Initialise the polygon-level variables. --------------------------------------!
      call init_ed_poly_vars(cgrid)
   end do gridloop

   return
end subroutine read_ed10_ed20_history_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will create the ED1.0/ED2.0 file name based on the longitude and     !
! latitude.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine create_ed10_ed20_fname(lat,ed_res,lon,sfilin,pss_name,css_name,site_name)
   use ed_max_dims, only : str_len
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real                  , intent(in)  :: lat
   real                  , intent(in)  :: lon
   real                  , intent(in)  :: ed_res
   character(len=*)      , intent(in)  :: sfilin
   character(len=str_len), intent(out) :: pss_name
   character(len=str_len), intent(out) :: css_name
   character(len=str_len), intent(out) :: site_name
   !----- Local variables. ----------------------------------------------------------------!
   real                                :: flon
   real                                :: flat
   character(len=str_len)              :: ed_fname
   !---------------------------------------------------------------------------------------!



   !----- Find the latitude and longitude corresponding to the dataset. -------------------!
   if (lat >= 0.0) then
      flat =   ed_res * real(int( lat / ed_res)) + 0.5 * ed_res 
   else
      flat = - ed_res * real(int(-lat / ed_res)) - 0.5 * ed_res
   end if
   if (lon >= 0.0) then
      flon =   ed_res * real(int( lon / ed_res)) + 0.5 * ed_res 
   else
      flon = - ed_res * real(int(-lon / ed_res)) - 0.5 * ed_res
   endif

   if (ed_res > 0.999 .and. ed_res < 1.001) then

      if (flat <= -10.0) then
         write(ed_fname,'(a,f5.1)')trim(sfilin)//'lat',flat
      elseif(flat < 0.0 .or. flat >= 10.0)then
         write(ed_fname,'(a,f4.1)')trim(sfilin)//'lat',flat
      else
         write(ed_fname,'(a,f3.1)')trim(sfilin)//'lat',flat
      end if
      if (flon <= -100.0) then
         write(ed_fname,'(a,f6.1)')trim(ed_fname)//'lon',flon
      elseif(flon <= -10.0 .or. flon >= 100.0)then
         write(ed_fname,'(a,f5.1)')trim(ed_fname)//'lon',flon
      elseif(flon < 0.0)then
         write(ed_fname,'(a,f4.1)')trim(ed_fname)//'lon',flon
      elseif(flon < 10.0)then
         write(ed_fname,'(a,f3.1)')trim(ed_fname)//'lon',flon
      else
         write(ed_fname,'(a,f4.1)')trim(ed_fname)//'lon',flon
      end if
   elseif (ed_res > 0.0999 .and. ed_res < 0.1001) then

      if (flat <= -10.0) then
         write(ed_fname,'(a,f6.2)')trim(sfilin)//'lat',flat
      elseif(flat < 0.0 .or. flat >= 10.0)then
         write(ed_fname,'(a,f5.2)')trim(sfilin)//'lat',flat
      else
         write(ed_fname,'(a,f4.2)')trim(sfilin)//'lat',flat
      end if
      if (flon <= -100.0) then
         write(ed_fname,'(a,f7.2)')trim(ed_fname)//'lon',flon
      elseif(flon <= -10.0 .or. flon >= 100.0)then
         write(ed_fname,'(a,f6.2)')trim(ed_fname)//'lon',flon
      elseif(flon < 0.0)then
         write(ed_fname,'(a,f5.2)')trim(ed_fname)//'lon',flon
      elseif(flon < 10.0)then
         write(ed_fname,'(a,f4.2)')trim(ed_fname)//'lon',flon
      else
         write(ed_fname,'(a,f5.2)')trim(ed_fname)//'lon',flon
      end if
   else
      call fatal_error('ed_res must be 0.1 or 1.0 and yours isn''t... '                    &
                      ,'create_ed10_ed20_fname','ed_read_ed10_ed20_history.f90')
   end if

   pss_name  = trim(ed_fname)//'.pss'
   css_name  = trim(ed_fname)//'.css'
   site_name = trim(ed_fname)//'.site'

   return
end subroutine create_ed10_ed20_fname
!==========================================================================================!
!==========================================================================================!
