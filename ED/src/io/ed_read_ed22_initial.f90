!==========================================================================================!
!==========================================================================================!
!     This subroutine reads ED-2.2 initial files.  The format is very similar to the       !
! ED-1.0/ED-2.0 "history" format, except that it has additional fields to initialise the   !
! necromass, and a site file.                                                              !
!------------------------------------------------------------------------------------------!
subroutine read_ed22_initial_file

   use ed_max_dims         , only : n_pft                       & ! intent(in)
                                  , huge_site                   & ! intent(in)
                                  , huge_patch                  & ! intent(in)
                                  , huge_cohort                 & ! intent(in)
                                  , max_water                   & ! intent(in)
                                  , str_len                     & ! intent(in)
                                  , maxfiles                    & ! intent(in)
                                  , maxlist                     & ! intent(in)
                                  , undef_character             & ! intent(in)
                                  , undef_integer               & ! intent(in)
                                  , undef_real                  ! ! intent(in)
   use pft_coms            , only : q                           & ! intent(in)
                                  , qsw                         & ! intent(in)
                                  , qbark                       & ! intent(in)
                                  , SLA                         & ! intent(in)
                                  , min_dbh                     & ! intent(in)
                                  , is_grass                    & ! intent(in)
                                  , include_pft                 & ! intent(in)
                                  , include_pft_ag              & ! intent(in)
                                  , pft_1st_check               & ! intent(in)
                                  , agf_bs                      & ! intent(in)
                                  , f_bstorage_init             & ! intent(in)
                                  , include_these_pft           & ! intent(in)
                                  , leaf_turnover_rate          & ! intent(in)
                                  , vm0                         & ! intent(in)
                                  , rd0                         & ! intent(in)
                                  , negligible_nplant           ! ! intent(in)
   use ed_misc_coms        , only : sfilin                      ! ! intent(in)
   use consts_coms         , only : pio180                      & ! intent(in)
                                  , pio4                        & ! intent(in)
                                  , almost_zero                 & ! intent(in)
                                  , tiny_num                    ! ! intent(in)
   use ed_misc_coms        , only : use_target_year             & ! intent(in)
                                  , restart_target_year         ! ! intent(in)
   use ed_state_vars       , only : polygontype                 & ! variable type
                                  , sitetype                    & ! variable type
                                  , patchtype                   & ! variable type
                                  , edtype                      & ! variable type
                                  , edgrid_g                    & ! variable type
                                  , allocate_polygontype        & ! subroutine
                                  , allocate_sitetype           & ! subroutine
                                  , allocate_patchtype          ! ! subroutine
   use grid_coms           , only : ngrids                      & ! intent(in)
                                  , nzg                         ! ! intent(in)
   use soil_coms           , only : soil_hydro_scheme           & ! intent(in)
                                  , slz                         & ! intent(in)
                                  , slxkey_ref                  & ! intent(inout)
                                  , slxsand_ref                 & ! intent(inout)
                                  , slxsilt_ref                 & ! intent(inout)
                                  , slxclay_ref                 & ! intent(inout)
                                  , slhydro_ref                 & ! intent(inout)
                                  , slsoc_ref                   & ! intent(inout)
                                  , slph_ref                    & ! intent(inout)
                                  , slcec_ref                   & ! intent(inout)
                                  , sldbd_ref                   & ! intent(inout)
                                  , ed_gen_soil_table           ! ! subroutine
   use allometry           , only : bd2dbh                      & ! function
                                  , dbh2h                       & ! function
                                  , size2bd                     & ! function
                                  , size2bl                     & ! function
                                  , size2bt                     & ! function
                                  , size2xb                     & ! function
                                  , ed_balive                   & ! function
                                  , ed_biomass                  & ! function
                                  , area_indices                ! ! subroutine
   use fuse_fiss_utils     , only : sort_cohorts                & ! subroutine
                                  , sort_patches                ! ! subroutine
   use decomp_coms         , only : decomp_scheme               & ! intent(in)
                                  , c2n_structural              ! ! intent(in)
   use phenology_coms      , only : llspan_inf                  ! ! intent(in)
   use physiology_coms     , only : iddmort_scheme              & ! intent(in)
                                  , trait_plasticity_scheme     ! ! intent(in)
   use update_derived_utils, only : update_cohort_plastic_trait ! ! subroutine
   use ed_init             , only : soil_default_fill           & ! sub-routine
                                  , sfcdata_ed                  ! ! sub-routine
   use ed_type_init        , only : init_ed_cohort_vars         & ! subroutine
                                  , init_ed_patch_vars          & ! subroutine
                                  , init_ed_site_vars           & ! subroutine
                                  , init_ed_poly_vars           ! ! subroutine
   use ed_init             , only : calc_flow_routing           ! ! subroutine
   implicit none

   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter :: min_area = 1.d-7           ! Minimum acceptable area.
   real(kind=8), parameter :: min_ok   = 1.d-20          ! Minimum acceptable value for
                                                         !    any restart variable.
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)          , pointer                             :: cgrid
   type(polygontype)     , pointer                             :: cpoly
   type(sitetype)        , pointer                             :: csite
   type(patchtype)       , pointer                             :: cpatch
   character(len=str_len), dimension(maxlist)                  :: full_list
   character(len=str_len), dimension(maxfiles)                 :: sss_list
   character(len=str_len), dimension(maxfiles)                 :: pss_list
   character(len=str_len), dimension(maxfiles)                 :: css_list
   character(len=str_len), dimension(huge_site)                :: sname
   character(len=str_len), dimension(huge_patch)               :: psname
   character(len=str_len), dimension(huge_patch)               :: pname
   character(len=str_len), dimension(huge_cohort)              :: csname
   character(len=str_len), dimension(huge_cohort)              :: cpname
   character(len=str_len), dimension(huge_cohort)              :: cname
   character(len=str_len)                                      :: sss_name
   character(len=str_len)                                      :: pss_name
   character(len=str_len)                                      :: css_name
   character(len=str_len)                                      :: cdum
   integer               , dimension(huge_site)                :: nscol
   integer               , dimension(huge_site)                :: ntext
   integer               , dimension(huge_site)                :: lsl
   integer               , dimension(huge_site)                :: patch_count
   integer               , dimension(huge_site)                :: last_ipa
   integer               , dimension(huge_patch)               :: dtype
   integer               , dimension(huge_patch)               :: psite_id
   integer               , dimension(huge_patch)               :: ppatch_id
   integer               , dimension(huge_patch)               :: cohort_count
   integer               , dimension(huge_patch)               :: last_ico
   integer               , dimension(huge_cohort)              :: ipft
   integer               , dimension(huge_cohort)              :: cpatch_id
   integer               , dimension(huge_cohort)              :: csite_id
   integer                                                     :: year
   integer                                                     :: igr
   integer                                                     :: ipy
   integer                                                     :: gsi
   integer                                                     :: gpa
   integer                                                     :: gco
   integer                                                     :: isi
   integer                                                     :: ipa
   integer                                                     :: ico
   integer                                                     :: apa
   integer                                                     :: aco
   integer                                                     :: ierr
   integer                                                     :: nf
   integer                                                     :: nflist
   integer                                                     :: nflsss
   integer                                                     :: nflpss
   integer                                                     :: nflcss
   integer                                                     :: nclosest
   integer                                                     :: ncohorts
   integer                                                     :: npatches
   integer                                                     :: nsites
   logical              , dimension(n_pft)                     :: discarded_pft
   logical              , dimension(:)           , allocatable :: shmask
   logical                                                     :: single_poi
   real(kind=8)                                                :: darea
   real                 , dimension(huge_site)                 :: s_area
   real                 , dimension(huge_site)                 :: depth
   real                 , dimension(huge_site)                 :: sand
   real                 , dimension(huge_site)                 :: clay
   real                 , dimension(huge_site)                 :: slsoc
   real                 , dimension(huge_site)                 :: slph
   real                 , dimension(huge_site)                 :: slcec
   real                 , dimension(huge_site)                 :: sldbd
   real                 , dimension(huge_site)                 :: elevation
   real                 , dimension(huge_site)                 :: slope
   real                 , dimension(huge_site)                 :: aspect
   real                 , dimension(huge_site)                 :: tci
   real                 , dimension(huge_site)                 :: moist_f
   real                 , dimension(huge_site)                 :: moist_w
   real                 , dimension(huge_patch)                :: time
   real                 , dimension(huge_patch)                :: age
   real                 , dimension(huge_patch)                :: p_area
   real                 , dimension(huge_patch)                :: fgc
   real                 , dimension(huge_patch)                :: fsc
   real                 , dimension(huge_patch)                :: stgc
   real                 , dimension(huge_patch)                :: stgl
   real                 , dimension(huge_patch)                :: stsc
   real                 , dimension(huge_patch)                :: stsl
   real                 , dimension(huge_patch)                :: msc
   real                 , dimension(huge_patch)                :: ssc
   real                 , dimension(huge_patch)                :: psc
   real                 , dimension(huge_patch)                :: fsn
   real                 , dimension(huge_patch)                :: msn
   real                 , dimension(huge_cohort)               :: balive
   real                 , dimension(huge_cohort)               :: bdead
   real                 , dimension(huge_cohort)               :: nplant
   real                 , dimension(huge_cohort)               :: height
   real                 , dimension(huge_cohort)               :: dbh
   real                 , dimension(huge_cohort)               :: ctime
   real                 , dimension(maxfiles)                  :: slon_list
   real                 , dimension(maxfiles)                  :: slat_list
   real                 , dimension(maxfiles)                  :: plon_list
   real                 , dimension(maxfiles)                  :: plat_list
   real                 , dimension(maxfiles)                  :: clon_list
   real                 , dimension(maxfiles)                  :: clat_list
   real                 , dimension(maxfiles)                  :: file_sdist
   real                 , dimension(maxfiles)                  :: file_pdist
   real                 , dimension(maxfiles)                  :: file_cdist
   real                 , dimension(n_pft)                     :: leaf_lifespan
   real                 , dimension(:)           , allocatable :: ed_slz
   real                                                        :: dummy
   real                                                        :: area_sum
   !----- External function. --------------------------------------------------------------!
   real                 , external                             :: sngloff
   real                 , external                             :: dist_gc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Define PFT-dependent leaf life span, used for initialisation.                     !
   !---------------------------------------------------------------------------------------!
   leaf_lifespan(:) = merge( 12.0 / leaf_turnover_rate(:)                                  &
                           , llspan_inf                                                    &
                           , leaf_turnover_rate(:) > 0.0  )
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !        Allocate temporary variables used during initialisation.                       !
   !---------------------------------------------------------------------------------------!
   allocate(shmask(nzg))
   allocate(ed_slz(nzg))
   shmask(:) = .false.
   ed_slz(:) = slz(1:nzg)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Now we loop over all all grids and polygons, and fill them with patches and      !
   ! cohorts from the closest polygon.                                                     !
   !---------------------------------------------------------------------------------------!
   main_gridloop: do igr = 1,ngrids
      cgrid => edgrid_g(igr)

      !----- Retrieve all files with the specified prefix. --------------------------------!
      call ed_filelist(full_list,sfilin(igr),nflist)
      !------------------------------------------------------------------------------------!

      !----- Retrieve LON/LAT information for sites, patches and cohorts ------------------!
      call ed1_fileinfo('.sss',nflist,full_list,nflsss,sss_list,slon_list,slat_list)
      call ed1_fileinfo('.pss',nflist,full_list,nflpss,pss_list,plon_list,plat_list)
      call ed1_fileinfo('.css',nflist,full_list,nflcss,css_list,clon_list,clat_list)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Save logical flag to decide whether or not this is a single-grid, single-     !
      ! -polygon simulation.  This information allows us to change the default soil        !
      ! properties so they are site-specific.                                              !
      !------------------------------------------------------------------------------------!
      single_poi = (ngrids == 1) .and. (cgrid%npolygons == 1)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop through every polygon.                                                    !
      !------------------------------------------------------------------------------------!
      main_polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         
         !----- Initialise load adjacency with dummy value. -------------------------------!
         cgrid%load_adjacency(ipy) = 0
         cgrid%wbar          (ipy) = 0.0
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Reset patch and cohort count.                                               !
         !---------------------------------------------------------------------------------!
         patch_count (:) = 0
         cohort_count(:) = 0
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Initialise the distances as very large numbers, so if we don't fill all the  !
         ! patches and cohorts, we will not going to take non-sense as a valid polygon.    !
         !---------------------------------------------------------------------------------!
         file_sdist(:) = 1.e20
         file_pdist(:) = 1.e20
         file_cdist(:) = 1.e20
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !       Reset the other placeholder variables.                                    !
         !---------------------------------------------------------------------------------!
         !------ Site-level variables. ----------------------------------------------------!
         sname        (:) = undef_character
         nscol        (:) = undef_integer
         ntext        (:) = undef_integer
         lsl          (:) = undef_integer
         patch_count  (:) = undef_integer
         dtype        (:) = undef_integer
         s_area       (:) = undef_real
         depth        (:) = undef_real
         sand         (:) = undef_real
         clay         (:) = undef_real
         slsoc        (:) = undef_real
         slph         (:) = undef_real
         slcec        (:) = undef_real
         sldbd        (:) = undef_real
         elevation    (:) = undef_real
         slope        (:) = undef_real
         aspect       (:) = undef_real
         tci          (:) = undef_real
         moist_f      (:) = undef_real
         moist_w      (:) = undef_real
         !------ Patch-level variables. ---------------------------------------------------!
         psname       (:) = undef_character
         pname        (:) = undef_character
         psite_id     (:) = undef_integer
         ppatch_id    (:) = undef_integer
         cohort_count (:) = undef_integer
         time         (:) = undef_real
         age          (:) = undef_real
         p_area       (:) = undef_real
         fgc          (:) = undef_real
         fsc          (:) = undef_real
         stgc         (:) = undef_real
         stgl         (:) = undef_real
         stsc         (:) = undef_real
         stsl         (:) = undef_real
         msc          (:) = undef_real
         ssc          (:) = undef_real
         psc          (:) = undef_real
         fsn          (:) = undef_real
         msn          (:) = undef_real
         !------ Cohort-level variables. --------------------------------------------------!
         csname       (:) = undef_character
         cpname       (:) = undef_character
         cname        (:) = undef_character
         csite_id     (:) = undef_integer
         cpatch_id    (:) = undef_integer
         ipft         (:) = undef_integer
         balive       (:) = undef_real
         bdead        (:) = undef_real
         nplant       (:) = undef_real
         height       (:) = undef_real
         dbh          (:) = undef_real
         ctime        (:) = undef_real
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !    Compute the distances between every polygon in the initial files and the     !
         ! current polygon.                                                                !
         !---------------------------------------------------------------------------------!
         do nf=1,nflsss
            file_sdist(nf) = dist_gc(cgrid%lon(ipy),slon_list(nf)                          &
                                    ,cgrid%lat(ipy),slat_list(nf) )
         end do
         do nf=1,nflpss
            file_pdist(nf) = dist_gc(cgrid%lon(ipy),plon_list(nf)                          &
                                    ,cgrid%lat(ipy),plat_list(nf) )
         end do
         do nf=1,nflcss
            file_cdist(nf) = dist_gc(cgrid%lon(ipy),clon_list(nf)                          &
                                    ,cgrid%lat(ipy),clat_list(nf) )
         end do
         !---------------------------------------------------------------------------------!





         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !      Read site file.                                                            !
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Find the file that is the closest to the current polygon, based on the      !
         ! distance vector.                                                                !
         !---------------------------------------------------------------------------------!
         nclosest = minloc(file_sdist,dim=1)
         sss_name = trim(sss_list(nclosest))
         write (unit=*,fmt='(2a)') '+ Using site file: ',trim(sss_name)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Open the patch file and read skip the header.                                !
         !---------------------------------------------------------------------------------!
         open(unit=12,file=trim(sss_name),form='formatted',status='old',action='read')
         read(unit=12,fmt='(a4)')  cdum
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Now we loop over all patches and decide whether they should be included or  !
         ! not.                                                                            !
         !---------------------------------------------------------------------------------!
         gsi     = 1
         read_sites: do

            !------------------------------------------------------------------------------!
            !     We must check whether we are not exceeding the maximum number of sites   !
            ! that we can read.                                                            !
            !------------------------------------------------------------------------------!
            if (gsi > huge_site) then
               write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(sss_name)
               write (unit=*,fmt='(a)')       ' Number of sites is > HUGE_SITE...'
               write (unit=*,fmt='(a,1x,i7)') ' HUGE_SITE:',huge_site
               write (unit=*,fmt='(a)')       ' Increase HUGE_SITE to read this...'
               call fatal_error('Too many patches to be read...'                           &
                               ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
            end if
            !------------------------------------------------------------------------------!


            !----- Read line.  Exit loop when finished reading sites. ---------------------!
            read(unit=12,fmt=*,iostat=ierr) sname(gsi),darea,depth(gsi),nscol(gsi)         &
                                           ,ntext(gsi),sand(gsi),clay(gsi),slsoc(gsi)      &
                                           ,slph(gsi),slcec(gsi),sldbd(gsi),elevation(gsi) &
                                           ,slope(gsi),aspect(gsi),tci(gsi),moist_f(gsi)   &
                                           ,moist_w(gsi)
            if (ierr /= 0) exit read_sites
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Add site only if the area is above minimum, we avoid including          !
            ! sites that are tiny since they will be overwritten by the next site.         !
            !------------------------------------------------------------------------------!
            s_area(gsi) = sngloff(darea, min_area)
            if (s_area(gsi) > min_area) gsi = gsi + 1
            !------------------------------------------------------------------------------!
         end do read_sites
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!



         !------ Close the file. ----------------------------------------------------------!
         close(unit=12,status='keep')
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Here we determine the number of sites.  We also make sure that there is at !
         ! least one valid site, otherwise we issue an error.                              !
         !---------------------------------------------------------------------------------!
         nsites = gsi - 1
         if (nsites <= 0) then
            write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(sss_name)
            write (unit=*,fmt='(a)')       ' Invalid number of sites: ',nsites
            write (unit=*,fmt='(a)')       ' File is probably corrupted...'
         else
            !------ Make sure the sum of site areas is 1. ---------------------------------!
            area_sum         = sum(s_area(1:nsites))
            s_area(1:nsites) = s_area(1:nsites) / area_sum
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!






         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !      Read patch file.                                                           !
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Find the file that is the closest to the current polygon, based on the      !
         ! distance vector.                                                                !
         !---------------------------------------------------------------------------------!
         nclosest = minloc(file_pdist,dim=1)
         pss_name = trim(pss_list(nclosest))
         write (unit=*,fmt='(2a)') 'Using patch file: ',trim(pss_name)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Open the patch file and skip the header.                                     !
         !---------------------------------------------------------------------------------!
         open(unit=12,file=trim(pss_name),form='formatted',status='old',action='read')
         read(unit=12,fmt='(a4)')  cdum
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Loop over all patches and decide whether they should be included or not.     !
         !---------------------------------------------------------------------------------!
         gpa = 1
         read_patches: do
            !------------------------------------------------------------------------------!
            !     We must check whether we are not exceeding the maximum number of patches !
            ! that we can read.                                                            !
            !------------------------------------------------------------------------------!
            if (gpa > huge_patch) then
               write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(pss_name)
               write (unit=*,fmt='(a)')       ' Number of patches is > HUGE_PATCH...'
               write (unit=*,fmt='(a,1x,i7)') ' HUGE_PATCH:',huge_patch
               write (unit=*,fmt='(a)')       ' Increase HUGE_PATCH to read this...'
               call fatal_error('Too many patches to be read...'                           &
                               ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
            end if
            !------------------------------------------------------------------------------!



            !----- Read line.  Exit loop when finished reading patches. -------------------!
            read(unit=12,fmt=*,iostat=ierr) time(gpa),psname(gpa),pname(gpa),dtype(gpa)    &
                                           ,age(gpa),darea,fgc(gpa),fsc(gpa),stgc(gpa)     &
                                           ,stgl(gpa),stsc(gpa),stsl(gpa),msc(gpa)         &
                                           ,ssc(gpa),psc(gpa),fsn(gpa),msn(gpa)            &
                                           ,dummy,dummy,dummy,dummy
            if (ierr /= 0) exit read_patches
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Add site only if the area is above minimum, we avoid including          !
            ! sites that are tiny since they will be overwritten by the next site.         !
            !------------------------------------------------------------------------------!
            p_area(gpa) = sngloff(darea, min_area)
            if (p_area(gpa) > min_area) gpa = gpa + 1
            !------------------------------------------------------------------------------!
         end do read_patches
         !---------------------------------------------------------------------------------!


         !------ Close the file. ----------------------------------------------------------!
         close(unit=12,status='keep')
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Here we determine the number of patches.  We also make sure that there is  !
         ! at least one valid patch, otherwise we issue an error.                          !
         !---------------------------------------------------------------------------------!
         npatches = gpa - 1
         if (npatches <= 0) then
            write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(pss_name)
            write (unit=*,fmt='(a)')       ' Invalid number of patches: ',npatches
            write (unit=*,fmt='(a)')       ' File is probably corrupted...'
         end if
         !---------------------------------------------------------------------------------!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!





         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !      Read cohort file.                                                          !
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the file that is the closest to the current polygon, based on the      !
         ! distance vector.                                                                !
         !---------------------------------------------------------------------------------!
         nclosest = minloc(abs(file_pdist(nclosest)-file_cdist),dim=1)
         css_name = trim(css_list(nclosest))
         write (unit=*,fmt='(2a)') 'Using cohort file: ',trim(css_name)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Open the cohort file and read in all cohorts.                                !
         !---------------------------------------------------------------------------------!
         open(unit=12,file=trim(css_name),form='formatted',status='old')
         read(unit=12,fmt='(a4)')  cdum
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Loop over all cohorts and decide whether they should be included or not.     !
         !---------------------------------------------------------------------------------!
         discarded_pft(:) = .false.
         gco              = 1
         read_cohorts: do

            !------------------------------------------------------------------------------!
            !     We must check whether we are not exceeding the maximum number of patches !
            ! that we can read.                                                            !
            !------------------------------------------------------------------------------!
            if (gco > huge_cohort) then
               write (unit=*,fmt='(a,1x,a)' ) ' In file:',trim(css_name)
               write (unit=*,fmt='(a)')       ' Number of cohorts is > HUGE_COHORT...'
               write (unit=*,fmt='(a,1x,i7)') ' HUGE_COHORT:',huge_cohort
               write (unit=*,fmt='(a)')       ' Increase HUGE_COHORT to read this...'
               call fatal_error('Too many cohorts to be read...'                           &
                               ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
            end if
            !------------------------------------------------------------------------------!





            !----- Read line.  Exit loop when finished reading cohorts. -------------------!
            read(unit=12,fmt=*,iostat=ierr) ctime(gco),csname(gco),cpname(gco),cname(gco)  &
                                           ,dbh(gco),height(gco),ipft(gco),nplant(gco)     &
                                           ,bdead(gco),balive(gco),dummy,dummy
            if (ierr /= 0) exit read_cohorts
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Make sure this cohort qualifies to be included (if not qualified, we     !
            ! cycle the read_cohorts loop without updating gco.                            !
            !------------------------------------------------------------------------------!
            year = int(ctime(gco))
            if (use_target_year == 1 .and. year /= restart_target_year) then
               !----- User wants css year to match the restart target year. ---------------!
               cycle read_cohorts
               !---------------------------------------------------------------------------!
            else if (nplant(gco) < negligible_nplant(ipft(gco))) then
               !----- Cohort population is negligible. ------------------------------------!
               cycle read_cohorts
               !---------------------------------------------------------------------------!
            else if (.not. include_pft(ipft(gco))) then
               !---------------------------------------------------------------------------!
               !     This PFT is not in the list of PFTs to include.  Decide what to do    !
               ! based on the PFT_1ST_CHECK settings.                                      !
               !---------------------------------------------------------------------------!
               select case (pft_1st_check)
               case (0)
                  !----- Stop the run. ----------------------------------------------------!
                  write (unit=*,fmt='(a,1x,a)' ) ' In file:',trim(css_name)
                  write (unit=*,fmt='(a,1x,i5,1x,a)')                                      &
                             'There are cohorts of PFT =',ipft(gco)                        &
                           ,', which are not defined in NL%INCLUDE_THESE_PFT.'
                  write (unit=*,fmt='(a,1x,a)') ' Either edit NL%INCLUDE_THESE_PFT or'     &
                                               ,'set NL%PFT_1ST_CHECK to 1 or 2.'
                  call fatal_error('Invalid PFT in initial file'                           &
                                  ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
                  !------------------------------------------------------------------------!
               case (1)
                  !----- Warn user about the unexpected PFT. ------------------------------!
                  write (unit=*,fmt='(a,1x,a)' ) ' In file:',trim(css_name)
                  write (unit=*,fmt='(a,1x,i5,1x,a)')                                      &
                             'There are cohorts of PFT =',ipft(gco)                        &
                           ,', which are not defined in NL%INCLUDE_THESE_PFT.'
                  write (unit=*,fmt='(a,1x,a)') ' Changing include_these_pft to'           &
                                               ,'incorporate the PFT.'
                  call warning('Unexpected PFT in initial file'                            &
                              ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
                  !------------------------------------------------------------------------!


                  !----- Add the unexpected PFT to the list of possible PFTs. -------------!
                  include_pft(ipft(gco))                 = .true.
                  include_these_pft(count(include_pft)) = ipft(gco)
                  call sort_up(include_these_pft,n_pft)
                  if (is_grass(ipft(gco))) include_pft_ag(ipft(gco)) = .true.
                  !------------------------------------------------------------------------!
               case (2)
                  if (.not. discarded_pft(ipft(gco))) then
                     !---------------------------------------------------------------------!
                     !     In case this is the first time finding this unexpected PFT,     !
                     ! warn user.                                                          !
                     !---------------------------------------------------------------------!
                     write (unit=*,fmt='(a,1x,a)' ) ' In file:',trim(css_name)
                     write (unit=*,fmt='(a,1x,i5,1x,a)')                                   &
                                'There are cohorts of PFT =',ipft(gco)                     &
                              ,', which are not defined in NL%INCLUDE_THESE_PFT.'
                     write (unit=*,fmt='(a,1x,a)') ' Discarding PFT.'
                     call warning('Unexpected PFT in initial file'                         &
                                 ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
                     !---------------------------------------------------------------------!


                     !------ Switch flag so we don't overwhelm the output with warnings. --!
                     discarded_pft(ipft(gco)) = .true.
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

                  !----- Skip the cohort. -------------------------------------------------!
                  cycle read_cohorts
                  !------------------------------------------------------------------------!
               end select
               !---------------------------------------------------------------------------!
            else
               !----- Keep the cohort. ----------------------------------------------------!
               gco = gco + 1
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do read_cohorts
         !---------------------------------------------------------------------------------!


         !------ Close the file. ----------------------------------------------------------!
         close(unit=12,status='keep')
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Here we determine the number of sites.  We also make sure that there is at !
         ! least one valid site, otherwise we issue an error.                              !
         !---------------------------------------------------------------------------------!
         ncohorts = gco - 1
         if (npatches < 0) then
            write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(pss_name)
            write (unit=*,fmt='(a)')       ' Invalid number of patches: ',npatches
            write (unit=*,fmt='(a)')       ' File is probably corrupted...'
         end if
         !---------------------------------------------------------------------------------!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!




         !---------------------------------------------------------------------------------!
         !    Link patches and cohorts to sites.  To reduce the computational burden of    !
         ! long patch/cohort files, we skip the beginning of the loop when all the initial !
         ! patches and cohorts have already been assigned.  This will only be effective    !
         ! if the files are organised (most of the time they are, but the code doesn't     !
         ! assume they are).                                                               !
         !---------------------------------------------------------------------------------!
         apa = 1
         aco = 1
         do gsi=1,nsites
            !----- Flag all patches associated with this site. ----------------------------!
            do gpa=apa,npatches
               !----- Link patch to site if the names match. ------------------------------!
               if (trim(psname(gpa)) == trim(sname(gsi))) then
                  psite_id(gpa) = gsi


                  !------------------------------------------------------------------------!
                  !       Update apa in case all patch elements up to this point have been !
                  ! assigned to a site.                                                    !
                  !------------------------------------------------------------------------!
                  if (gpa == apa) then
                     apa = apa + 1
                  end if
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Link cohorts to patches.                                               !
               !---------------------------------------------------------------------------!
               do gco=aco,ncohorts
                  !----- Flag all cohorts associated with this patch and site. ------------!
                  if ( (trim(cpname(gco)) == trim(pname (gpa))) .and.                      &
                       (trim(csname(gco)) == trim(psname(gpa))) .and.                      &
                       (trim(csname(gco)) == trim(sname (gsi)))       ) then
                     cpatch_id(gco) = gpa
                     csite_id (gco) = gsi

                     !---------------------------------------------------------------------!
                     !       Update apa in case all cohort elements up to this point have  !
                     ! been assigned to a patch/site.                                      !
                     !---------------------------------------------------------------------!
                     if (gco == aco) then
                        aco = aco + 1
                     end if
                     !---------------------------------------------------------------------!

                  end if
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!


               !----- Count cohorts belonging to this patch. ------------------------------!
               cohort_count(gpa) = count(cpatch_id(:) == gpa)
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!


            !----- Count patches belonging to this site. ----------------------------------!
            patch_count(gsi) = count(psite_id(:) == gsi)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     In case this site does not have any patch, issue an error.  Sites should !
            ! always have at least one patch.                                              !
            !------------------------------------------------------------------------------!
            if (patch_count(gsi) == 0) then
               write (unit=*,fmt='(a)'      )  '=========================================='
               write (unit=*,fmt='(a)'      )  '  Site without patches found!'
               write (unit=*,fmt='(a)'      )  '=========================================='
               write (unit=*,fmt='(a,1x,a)' )  ' Site file: ',trim(sss_name)
               write (unit=*,fmt='(a,1x,a)' )  ' Patch file:',trim(pss_name)
               write (unit=*,fmt='(a,1x,i5)')  ' Site ID:   ',gsi
               write (unit=*,fmt='(a,1x,a)' )  ' Site name: ',sname(gsi)
               write (unit=*,fmt='(a)'      )  '=========================================='
               call fatal_error('A site without corresponding patches was found!'          &
                               ,'read_ed22_initial_file','ed_read_ed22_initial.f90')
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !      Allocate sites.                                                            !
         !---------------------------------------------------------------------------------!
         call allocate_polygontype(cpoly,nsites)
         call soil_default_fill(cgrid,igr,ipy)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !       Loop through sites and allocate patches and cohorts.                      !
         !---------------------------------------------------------------------------------!
         gsi = 0
         init_sites: do isi=1,cpoly%nsites
            !------ Update pointers and counters. -----------------------------------------!
            csite => cpoly%site(isi)
            gsi = gsi + 1
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     In case this is a single-grid, single-polygon simulation, we rewrite     !
            ! the site properties of the first soil classes with the site-level            !
            ! information.                                                                 !
            !------------------------------------------------------------------------------!
            if (single_poi) then
               !----- Replace texture with the site ID. -----------------------------------!
               ntext(gsi) = isi
               !---------------------------------------------------------------------------!


               !----- Overwrite reference properties of the isi-th site. ------------------!
               slxkey_ref (isi) = 'Site'
               slhydro_ref(isi) = soil_hydro_scheme
               slxsand_ref(isi) = sand (gsi)
               slxclay_ref(isi) = clay (gsi)
               slxsilt_ref(isi) = 1. - slxsand_ref(isi) - slxclay_ref(isi)
               slsoc_ref  (isi) = slsoc(gsi)
               slph_ref   (isi) = slph (gsi)
               slcec_ref  (isi) = slcec(gsi)
               sldbd_ref  (isi) = sldbd(gsi)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Initialise site variables.  If TOPMODEL becomes again functional, this    !
            ! list may need to be expanded.                                                !
            !------------------------------------------------------------------------------!
            cpoly%area        (isi) = s_area   (gsi)
            cpoly%ntext_soil(:,isi) = ntext    (gsi)
            cpoly%ncol_soil   (isi) = nscol    (gsi)
            cpoly%elevation   (isi) = elevation(gsi)
            cpoly%slope       (isi) = slope    (gsi)
            cpoly%aspect      (isi) = aspect   (gsi)
            cpoly%tci         (isi) = tci      (gsi)
            cpoly%moist_f     (isi) = moist_f  (gsi)
            cpoly%moist_w     (isi) = moist_w  (gsi)
            !------------------------------------------------------------------------------!


            !------ Dummy variables. ------------------------------------------------------!
            cpoly%sitenum    (isi) = isi
            !------------------------------------------------------------------------------!


            !------ Find the lowest soil level to simulate. -------------------------------!
            shmask(:)      = ed_slz(:) <= - abs(depth(gsi))
            cpoly%lsl(isi) = min(max(1,maxloc(ed_slz(:),dim=1,mask=shmask)),nzg-1)
            !------------------------------------------------------------------------------!


            !------ Copy patch count. -----------------------------------------------------!
            cpoly%patch_count(isi) = patch_count(gsi)
            !------------------------------------------------------------------------------!


            !------ Update the polygon average topographic moisture index. ----------------!
            cgrid%wbar(ipy) = cgrid%wbar(ipy) + cpoly%moist_w(isi) * cpoly%area(isi)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Allocate patches for this site.                                          !
            !------------------------------------------------------------------------------!
            call allocate_sitetype(csite,patch_count(gsi))
            !------------------------------------------------------------------------------!
         end do init_sites
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !       Loop through all patches, and initialise them at the appropriate site.    !
         !---------------------------------------------------------------------------------!
         last_ipa(:) = 0
         init_patches: do gpa=1,npatches
            !------------------------------------------------------------------------------!
            !    Try setting the site.  If the site had a tiny area, this patch may be     !
            ! orphaned, in which case we skip it.                                          !
            !------------------------------------------------------------------------------!
            if (psite_id(gpa) == undef_integer) then
               !----- Invalid patch, skip it. ---------------------------------------------!
               cycle init_patches
               !---------------------------------------------------------------------------!
            else
               !----- Valid patch, assign site. -------------------------------------------!
               isi =  psite_id(gpa)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !----- Set other counters and pointers. ---------------------------------------!
            ipa            =  last_ipa(isi) + 1
            ppatch_id(gpa) =  ipa
            csite          => cpoly%site (isi)
            cpatch         => csite%patch(ipa)
            !------------------------------------------------------------------------------!



            !------ Set patch-level variables. --------------------------------------------!
            csite%pname             (ipa) = trim(pname  (gpa))
            csite%dist_type         (ipa) = dtype       (gpa)
            csite%age               (ipa) = age         (gpa)
            csite%area              (ipa) = p_area      (gpa)
            csite%fast_grnd_C       (ipa) = fgc         (gpa)
            csite%fast_soil_C       (ipa) = fsc         (gpa)
            csite%structural_grnd_C (ipa) = stgc        (gpa)
            csite%structural_grnd_L (ipa) = stgl        (gpa)
            csite%structural_soil_C (ipa) = stsc        (gpa)
            csite%structural_soil_L (ipa) = stsl        (gpa)
            csite%mineralized_soil_N(ipa) = msn         (gpa)
            csite%cohort_count      (ipa) = cohort_count(gpa)
            !------------------------------------------------------------------------------!



            !------ Set fast nitrogen pools so they are proportional to the carbon pools. -!
            if ( (fgc(gpa)+fsc(gpa)) > tiny_num) then
               csite%fast_grnd_N      (ipa) = fgc(gpa) * fsn(gpa) / (fgc(gpa)+fsc(gpa))
               csite%fast_soil_N      (ipa) = fsc(gpa) * fsn(gpa) / (fgc(gpa)+fsc(gpa))
            else
               csite%fast_grnd_C      (ipa) = 0.0
               csite%fast_soil_C      (ipa) = 0.0
               csite%fast_grnd_N      (ipa) = 0.0
               csite%fast_soil_N      (ipa) = 0.0
            end if
            !------------------------------------------------------------------------------!



            !------ Use stoichiometry to derive structural N pools. -----------------------!
            csite%structural_grnd_N (ipa) = csite%structural_grnd_C (ipa) / c2n_structural
            csite%structural_soil_N (ipa) = csite%structural_soil_C (ipa) / c2n_structural
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Check decomposition scheme before assigning microbial carbon.            !
            !------------------------------------------------------------------------------!
            select case (decomp_scheme)
            case (5)
               csite%microbial_soil_C(ipa) = msc(gpa)
               csite%slow_soil_C     (ipa) = ssc(gpa)
               csite%passive_soil_C  (ipa) = psc(gpa)
            case default
               csite%microbial_soil_C(ipa) = 0.0
               csite%slow_soil_C     (ipa) = msc(gpa) + ssc(gpa) + psc(gpa)
               csite%passive_soil_C  (ipa) = 0.0
            end select
            !------------------------------------------------------------------------------!


            !----- Initialise other properties. -------------------------------------------!
            csite%fbeam             (ipa) = 1.0
            csite%light_type        (ipa) = 1
            csite%sum_dgd           (ipa) = 0.0
            csite%sum_chd           (ipa) = 0.0
            csite%plant_ag_biomass  (ipa) = 0.0
            !------------------------------------------------------------------------------!


            !------ Allocate cohorts for this patch. --------------------------------------!
            if ( cohort_count(gpa) /= 0) then
               call allocate_patchtype(cpatch,cohort_count(gpa))
            end if
            !------------------------------------------------------------------------------!


            !------ Update last_ipa for this site. ----------------------------------------!
            last_ipa(isi) = ipa
            !------------------------------------------------------------------------------!
         end do init_patches
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !       Loop through all cohort, and initialise them at the appropriate patch and !
         ! site.                                                                           !
         !---------------------------------------------------------------------------------!
         last_ico(:) = 0
         init_cohorts: do gco=1,ncohorts
            !------------------------------------------------------------------------------!
            !    Try setting the site and patch.  If the site or the patch had a tiny      !
            ! area, this cohort may be orphaned, in which case we skip it.                 !
            !------------------------------------------------------------------------------!
            if (csite_id(gco) == undef_integer .or. cpatch_id(gco) == undef_integer) then
               !----- Invalid cohort, skip it. --------------------------------------------!
               cycle init_cohorts
               !---------------------------------------------------------------------------!
            else
               !----- Valid cohort, assign site and global patch. -------------------------!
               isi = csite_id (gco)
               gpa = cpatch_id(gco)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !----- Set other counters and pointers. ---------------------------------------!
            ipa    =  ppatch_id(gpa)
            ico    =  last_ico (gpa) + 1
            csite  => cpoly%site(isi)
            cpatch => csite%patch(ipa)
            !------------------------------------------------------------------------------!




            !------ Copy data from files to cohort. ---------------------------------------!
            cpatch%nplant(ico) = nplant(gco)
            cpatch%pft   (ico) = ipft  (gco)
            cpatch%dbh   (ico) = max(dbh(gco),min_dbh(ipft(gco)))
            !------------------------------------------------------------------------------!




            !------ Update allometry to define height and heartwood. ----------------------!
            cpatch%height(ico) = dbh2h(cpatch%pft(ico),cpatch%dbh(ico))
            bdead        (gco) = size2bd(cpatch%dbh(ico),cpatch%height(ico)                &
                                        ,cpatch%pft(ico))
            cpatch%bdeada(ico) =        agf_bs(cpatch%pft(ico))  * bdead(gco)
            cpatch%bdeadb(ico) = (1.0 - agf_bs(cpatch%pft(ico))) * bdead(gco)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Initialise SLA, Vm0, Rd0, and  with the look-up table value.  These      !
            ! variables may be updated during phenology initialisation, or trait           !
            ! plasticity, but they must have an initial assignment so we can even          !
            ! calculate the initial area indices and inicial trait values needed for the   !
            ! trait update.                                                                !
            !------------------------------------------------------------------------------!
            cpatch%sla   (ico) = SLA               (ipft(gco))
            cpatch%vm_bar(ico) = Vm0               (ipft(gco))
            cpatch%rd_bar(ico) = Rd0               (ipft(gco))
            cpatch%llspan(ico) = leaf_turnover_rate(ipft(gco))
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Use allometry to define leaf and the other live biomass pools.           !
            !------------------------------------------------------------------------------!
            cpatch%bleaf    (ico) = size2bl(cpatch%dbh(ico),cpatch%height(ico)             &
                                           ,cpatch%sla(ico),ipft(gco))
            cpatch%broot    (ico) = cpatch%bleaf(ico) * q(ipft(gco))
            cpatch%bsapwooda(ico) = agf_bs(ipft(gco))                                      &
                                  * cpatch%bleaf(ico) * qsw(ipft(gco))                     &
                                  * cpatch%height(ico)
            cpatch%bsapwoodb(ico) = (1.-agf_bs(ipft(gco)))                                 &
                                  * cpatch%bleaf(ico) * qsw(ipft(gco))                     &
                                  * cpatch%height(ico)
            cpatch%bbarka(ico)    = agf_bs(ipft(gco))                                      &
                                  * cpatch%bleaf(ico) * qbark(ipft(gco))                   &
                                  * cpatch%height(ico)
            cpatch%bbarkb(ico)    = (1.-agf_bs(ipft(gco)))                                 &
                                  * cpatch%bleaf(ico) * qbark(ipft(gco))                   &
                                  * cpatch%height(ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Start plants with full phenology, we will take care of phenology after   !
            ! this sub-routine.                                                            !
            !------------------------------------------------------------------------------!
            cpatch%phenology_status(ico) = 0
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     In case we are representing trait plasticity, update traits (SLA, Vm0).  !
            ! This must be done before calculating LAI and before ed_balive.               !
            !------------------------------------------------------------------------------!
            select case (trait_plasticity_scheme)
            case (0)
               continue
            case default
               call update_cohort_plastic_trait(cpatch,ico,.true.                          &
                                               ,leaf_lifespan(ipft(gco))                   &
                                               ,vm0          (ipft(gco))                   &
                                               ,rd0          (ipft(gco))                   &
                                               ,sla          (ipft(gco)) )
            end select
            !------------------------------------------------------------------------------!




            !----- Assign biomass of living tissues. --------------------------------------!
            cpatch%balive(ico) = ed_balive(cpatch, ico)
            !------------------------------------------------------------------------------!


            !----- Initialise storage biomass (after setting balive). ---------------------!
            cpatch%bstorage(ico)  = max( almost_zero, f_bstorage_init(ipft(gco)))          &
                                  * cpatch%balive(ico)
            !------------------------------------------------------------------------------!



            !----- Assign LAI, WAI, and CAI -----------------------------------------------!
            call area_indices(cpatch, ico)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Initialise the carbon balance.  For initial conditions, we always assume !
            ! storage biomass for the previous months so the scale is correct (carbon      !
            ! balance is given in kgC/pl).  The current month carbon balance must be       !
            ! initialised consistently with the iddmort_scheme we are using.               !
            !------------------------------------------------------------------------------!
            cpatch%cb         (1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_lightmax(1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_moistmax(1:12,ico) = cpatch%bstorage(ico)
            cpatch%cb_mlmax   (1:12,ico) = cpatch%bstorage(ico)
            select case (iddmort_scheme)
            case (0)
               !------ Storage is not accounted. ------------------------------------------!
               cpatch%cb         (13,ico) = 0.0
               cpatch%cb_lightmax(13,ico) = 0.0
               cpatch%cb_moistmax(13,ico) = 0.0
               cpatch%cb_mlmax   (13,ico) = 0.0
               !---------------------------------------------------------------------------!
            case (1)
               !------ Storage is accounted. ----------------------------------------------!
               cpatch%cb         (13,ico) = cpatch%bstorage(ico)
               cpatch%cb_lightmax(13,ico) = cpatch%bstorage(ico)
               cpatch%cb_moistmax(13,ico) = cpatch%bstorage(ico)
               cpatch%cb_mlmax   (13,ico) = cpatch%bstorage(ico)
               !---------------------------------------------------------------------------!
            end select
            cpatch%cbr_bar          (ico) = 1.0
            !------------------------------------------------------------------------------!



            !----- Above ground biomass, use the allometry. -------------------------------!
            cpatch%agb    (ico) = ed_biomass(cpatch, ico)
            cpatch%basarea(ico) = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)
            cpatch%btimber(ico) = size2bt( cpatch%dbh       (ico)                          &
                                         , cpatch%height    (ico)                          &
                                         , cpatch%bdeada    (ico)                          &
                                         , cpatch%bsapwooda (ico)                          &
                                         , cpatch%bbarka    (ico)                          &
                                         , cpatch%pft       (ico) )
            cpatch%thbark (ico) = size2xb( cpatch%dbh       (ico)                          &
                                         , cpatch%height    (ico)                          &
                                         , cpatch%bbarka    (ico)                          &
                                         , cpatch%bbarkb    (ico)                          &
                                         , cpatch%sla       (ico)                          &
                                         , cpatch%pft       (ico) )
            !------------------------------------------------------------------------------!



            !----- Growth rates, start with zero. -----------------------------------------!
            cpatch%dagb_dt  (ico)  = 0.
            cpatch%dlnagb_dt(ico)  = 0.
            cpatch%dba_dt   (ico)  = 0.
            cpatch%dlnba_dt (ico)  = 0.
            cpatch%ddbh_dt  (ico)  = 0.
            cpatch%dlndbh_dt(ico)  = 0.
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Initialise other cohort variables.  Some of them won't be updated       !
            ! unless the lai exceeds lai_min.                                              !
            !------------------------------------------------------------------------------!
            cpatch%fsw(ico)   = 1.0
            cpatch%gpp(ico)   = 0.0
            cpatch%par_l(ico) = 0.0
            !------------------------------------------------------------------------------!


            !----- Update the patch level above-ground biomass. ---------------------------!
            csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                      &
                                        + cpatch%agb(ico) * cpatch%nplant(ico)
            !------------------------------------------------------------------------------!



            !------ Update last_ico for this patch. ---------------------------------------!
            last_ico(gpa) = ico
            !------------------------------------------------------------------------------!
         end do init_cohorts
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find site adjacency (unlikely to be needed here).
         !---------------------------------------------------------------------------------!
         if (cgrid%load_adjacency(ipy) /= 0) then
            call calc_flow_routing(cgrid,ipy)
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Initialise additional site-, patch-, and cohort-level variables.           !
         !---------------------------------------------------------------------------------!
         init2_sites: do isi = 1,cpoly%nsites

            !----- Make sure that the total patch area is 1. ------------------------------!
            csite                        => cpoly%site(isi)
            area_sum                     =  sum(csite%area(1:csite%npatches))
            csite%area(1:csite%npatches) =  csite%area(1:csite%npatches) / area_sum
            !------------------------------------------------------------------------------!




            !----- Initialise the cohort variables, then sort them by size. ---------------!
            init2_patches: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Initialise additional cohort variables. -----------------------------!
               init2_cohorts: do ico = 1,cpatch%ncohorts
                  call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi),nzg                   &
                                          ,cpoly%ntext_soil(:,isi))
               end do init2_cohorts
               !---------------------------------------------------------------------------!

               !----- Make sure that cohorts are organised from tallest to shortest. ------!
               call sort_cohorts(cpatch)
               !---------------------------------------------------------------------------!
            end do init2_patches
            !------------------------------------------------------------------------------!



            !----- Initialise the remaining patch-level variables. ------------------------!
            call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
            !------------------------------------------------------------------------------!

            !----- Make sure that patches are organised from oldest to youngest. ----------!
            call sort_patches(csite)
            !------------------------------------------------------------------------------!
         end do init2_sites
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      In case this is a single site simulation, we must update the soil          !
         ! properties before initialising the site-level variables, because soil           !
         ! characteristics are read from the site file.                                    !
         !---------------------------------------------------------------------------------!
         if (single_poi) then
            call sfcdata_ed()
         end if
         !---------------------------------------------------------------------------------!



         !----- Initialise the remaining site-level variables. ----------------------------!
         call init_ed_site_vars(cpoly)
         !---------------------------------------------------------------------------------!
      end do main_polyloop
      !------------------------------------------------------------------------------------!



      !----- Initialise the polygon-level variables. --------------------------------------!
      call init_ed_poly_vars(cgrid)
      !------------------------------------------------------------------------------------!
   end do main_gridloop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Update the soil parameter table as soil properties have been overwritten for     !
   ! a few sites.                                                                          !
   !---------------------------------------------------------------------------------------!
   call ed_gen_soil_table()
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !        Free memory.                                                                   !
   !---------------------------------------------------------------------------------------!
   deallocate(shmask)
   deallocate(ed_slz)
   !---------------------------------------------------------------------------------------!


   return
end subroutine read_ed22_initial_file
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
!==========================================================================================!
