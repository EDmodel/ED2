module landuse_init
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This file will read the landuse files and assign the anthropogenic disturbance    !
   ! matrices.                                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine read_landuse_matrix

      use ed_state_vars , only : edtype              & ! structure
                               , polygontype         & ! structure
                               , sitetype            & ! structure
                               , edgrid_g            ! ! structure
      use pft_coms      , only : is_grass            & ! intent(in)
                               , is_liana            ! ! intent(in)
      use consts_coms   , only : erad                & ! intent(in)
                               , pio180              & ! intent(in)
                               , lnexp_max           ! ! intent(in)
      use disturb_coms  , only : lutime              & ! intent(in)
                               , max_lu_years        & ! intent(in)
                               , num_lu_trans        & ! intent(in)
                               , ianth_disturb       & ! intent(in)
                               , lu_database         & ! intent(in)
                               , sl_pft              & ! intent(in)
                               , sl_scale            & ! intent(in)
                               , sl_nyrs             & ! intent(in)
                               , sl_yr_first         & ! intent(in)
                               , sl_skid_rel_area    & ! intent(in)
                               , sl_skid_dbh_thresh  & ! intent(in)
                               , sl_skid_s_gtharv    & ! intent(in)
                               , sl_skid_s_ltharv    & ! intent(in)
                               , sl_felling_s_ltharv & ! intent(in)
                               , sl_mindbh_harvest   & ! intent(in)
                               , sl_prob_harvest     ! ! intent(in)
      use ed_misc_coms  , only : iyeara              & ! intent(in)
                               , iyearz              ! ! intent(in)
      use grid_coms     , only : ngrids              ! ! intent(in)
      use ed_max_dims   , only : str_len             & ! intent(in)
                               , huge_lu             & ! intent(in)
                               , n_pft               & ! intent(in)
                               , maxlist             & ! intent(in)
                               , undef_real          & ! intent(in)
                               , undef_integer       ! ! intent(in)
      use detailed_coms , only : idetailed           ! ! intent(in)

      implicit none
      !----- Local variables --------------------------------------------------------------!
      type(edtype)          , pointer                 :: cgrid
      type(polygontype)     , pointer                 :: cpoly
      type(sitetype)        , pointer                 :: csite
      type(lutime)          , pointer                 :: clutime
      type(lutime)          , pointer                 :: onelutime
      character(len=str_len), dimension(maxlist)      :: full_list
      character(len=str_len), dimension(maxlist)      :: lu_list
      real                  , dimension(maxlist)      :: llon_list
      real                  , dimension(maxlist)      :: llat_list
      real                  , dimension(maxlist)      :: file_ldist
      character(len=6)                                :: hform
      character(len=str_len)                          :: lu_name
      character(len=str_len)                          :: cdum
      character(len=str_len)                          :: vkey
      character(len=13)                               :: hifmt
      character(len=15)                               :: hffmt
      integer                                         :: nharvest
      integer               , dimension(n_pft)        :: harvest_pft
      integer                                         :: nf
      integer                                         :: nflist
      integer                                         :: nfllu
      integer                                         :: ncl
      integer                                         :: iyear
      integer                                         :: igr
      integer                                         :: ipy
      integer                                         :: isi
      integer                                         :: ipft
      integer                                         :: h
      integer                                         :: sim_years
      integer                                         :: yd_1st
      integer                                         :: yd_this
      integer                                         :: yd_last
      integer                                         :: poseq
      logical                                         :: inside
      logical                                         :: write_lu_settings
      real                                            :: skid_rel_area
      real                  , dimension(n_pft)        :: mindbh_slog
      real                  , dimension(n_pft)        :: harvprob_slog
      real                  , dimension(n_pft)        :: mindbh_fplt
      real                  , dimension(n_pft)        :: harvprob_fplt
      real                  , dimension(n_pft)        :: skid_dbh_thresh
      real                  , dimension(n_pft)        :: felling_s_gtharv
      real                  , dimension(n_pft)        :: felling_s_ltharv
      real                  , dimension(n_pft)        :: skid_s_gtharv
      real                  , dimension(n_pft)        :: skid_s_ltharv
      real                  , dimension(num_lu_trans) :: landuse_now
      real                                            :: lu_area
      real                                            :: lu_area_i
      real                                            :: wlon
      real                                            :: elon
      real                                            :: slat
      real                                            :: nlat
      !----- Local constants. -------------------------------------------------------------!
      character(len=12)     , parameter :: fffmt    = '(a,1x,f12.5)'
      character(len=13)     , parameter :: esfmt    = '(a,1x,es12.5)'
      real                  , parameter :: huge_dbh = huge(1.)
      character(len=str_len), parameter :: lu_table = "anth_disturb_table.txt"
      !----- External function. -----------------------------------------------------------!
      real                  , external  :: dist_gc
      real                  , external  :: solid_area
      !------------------------------------------------------------------------------------!



      !----- Find number of simulation years ----------------------------------------------!
      sim_years = iyearz-iyeara+1
      !------------------------------------------------------------------------------------!

      !------ Decide whether to write lu settings. ----------------------------------------!
      write_lu_settings = btest(idetailed,6) .and. ianth_disturb /= 0
      !------------------------------------------------------------------------------------!


      !----- Crashing the run if the user set up a very long run... -----------------------!
      if (ianth_disturb /= 0 .and. sim_years > max_lu_years) then
         write (unit=*,fmt='(a,1x,i5)') 'IYEARA       (From namelist)        :',iyeara
         write (unit=*,fmt='(a,1x,i5)') 'IYEARZ       (From namelist)        : ',iyearz
         write (unit=*,fmt='(a,1x,i5)') 'MAX_LU_YEARS (From disturb_coms.f90): '           &
                                                                              ,max_lu_years
         write (unit=*,fmt='(a)') ' Your run is too long.  Try increasing max_lu_years,'
         write (unit=*,fmt='(a)') ' so MAX_LU_YEARS >= IYEARZ-IYEARA+1, then recompile'
         write (unit=*,fmt='(a)') ' your model.'
         call fatal_error ('Simulation is too long for anthropogenic disturbance.'         &
                          ,'landuse_init','landuse_init.f90')
      end if
      !------------------------------------------------------------------------------------!


      gridloop: do igr = 1,ngrids

         !---------------------------------------------------------------------------------!
         !     Find the list of disturbance rate files.                                    !
         !---------------------------------------------------------------------------------!
         if (ianth_disturb == 1) then
            call ed_filelist(full_list,lu_database(igr),nflist)
            call ed1_fileinfo('.lu',nflist,full_list,nfllu,lu_list,llon_list,llat_list)
         end if
         !---------------------------------------------------------------------------------!

         cgrid=>edgrid_g(igr)

         polyloop: do ipy = 1,cgrid%npolygons
            cpoly => cgrid%polygon(ipy)

            select case (ianth_disturb)
            case (0)
               !----- No plantations. -----------------------------------------------------!
               cpoly%plantation(:) = 0
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Anthropogenic disturbance is not used this time, allocate only a     !
               ! single landuse year.                                                      !
               !---------------------------------------------------------------------------!
               allocate(cpoly%clutimes(1,cpoly%nsites))
               !---------------------------------------------------------------------------!

               !----- Set the parameters in a way that no logging/ploughing will happen. --!
               do isi = 1,cpoly%nsites
                  cpoly%num_landuse_years       (isi) = 1
                  cpoly%mindbh_harvest  (:,isi) = huge_dbh
                  cpoly%prob_harvest    (:,isi) = 0.
                  cpoly%skid_dbh_thresh (:,isi) = huge_dbh
                  cpoly%skid_s_gtharv   (:,isi) = 1.
                  cpoly%skid_s_ltharv   (:,isi) = 1.
                  cpoly%felling_s_gtharv(:,isi) = 1.
                  cpoly%felling_s_ltharv(:,isi) = 1.

                  cpoly%clutimes(1,isi)%landuse_year            = iyeara
                  cpoly%clutimes(1,isi)%landuse(1:num_lu_trans) = 0.0
               end do
               !---------------------------------------------------------------------------!

            case (1)

               !---------------------------------------------------------------------------!
               !    Initialise plantation patches if plantation information is available.  !
               !---------------------------------------------------------------------------!
               cpoly%plantation(:) = 0
               call read_plantation_fractions(cpoly,cgrid%lon(ipy),cgrid%lat(ipy),igr)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Initialise the PFT-dependent arrays with namelist data, but keep      !
               ! harvest probability zero and minimum DBH for harvesting to infinity, to   !
               ! prevent felling of these PFTs.                                            !
               !---------------------------------------------------------------------------!
               do isi = 1,cpoly%nsites
                  cpoly%mindbh_harvest  (:,isi) = huge_dbh
                  cpoly%prob_harvest    (:,isi) = 0.
                  cpoly%skid_dbh_thresh (:,isi) = merge(huge_dbh,sl_skid_dbh_thresh        &
                                                                              ,is_grass(:))
                  cpoly%skid_s_gtharv   (:,isi) = merge( 1.00,sl_skid_s_gtharv,is_grass(:))
                  cpoly%skid_s_ltharv   (:,isi) = merge( 1.00,sl_skid_s_ltharv,is_grass(:))
                  cpoly%felling_s_gtharv(:,isi) = merge( 0.70,            0.00,is_grass(:))
                  cpoly%felling_s_ltharv(:,isi) = merge( 0.70,sl_felling_s_ltharv          &
                                                                              ,is_grass(:))
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Compute the distance between the current polygon and all the files.   !
               !---------------------------------------------------------------------------!
               do nf=1,nfllu
                  file_ldist(nf) = dist_gc(cgrid%lon(ipy),llon_list(nf)                    &
                                          ,cgrid%lat(ipy),llat_list(nf) )
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Pick the closest file.  This is not a guarantee that it will be used   !
               ! because the closest polygon area must contain the point associated to the !
               ! current polygon.                                                          !
               !---------------------------------------------------------------------------!
               ncl     = minloc(file_ldist(1:nfllu),dim=1)
               lu_name = lu_list(ncl)
               write (unit=*,fmt='(2a)') 'Using land use file: ',trim(lu_name)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Open the patch file and read in all patches.                           !
               !---------------------------------------------------------------------------!
               open(unit=12,file=trim(lu_name),form='formatted',status='old',action='read')
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Top header variables. Initialise them with dummy values.  Except for   !
               ! skid_area, they must be all assigned.  Skid area by default the value     !
               ! from the namelist.                                                        !
               !---------------------------------------------------------------------------!
               wlon          = undef_real
               elon          = undef_real
               slat          = undef_real
               nlat          = undef_real
               lu_area       = undef_real
               yd_1st        = undef_integer
               yd_last       = undef_integer
               skid_rel_area = sl_skid_rel_area
               nharvest      = undef_integer
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Initialise the temporary arrays with dummy values. These are all      !
               ! optional variables, but they are useful to define logging strategies in   !
               ! more detail.                                                              !
               !---------------------------------------------------------------------------!
               harvest_pft     (1:n_pft) = undef_integer
               mindbh_slog     (1:n_pft) = undef_real
               harvprob_slog   (1:n_pft) = undef_real
               mindbh_fplt     (1:n_pft) = undef_real
               harvprob_fplt   (1:n_pft) = undef_real
               skid_dbh_thresh (1:n_pft) = undef_real
               skid_s_gtharv   (1:n_pft) = undef_real
               skid_s_ltharv   (1:n_pft) = undef_real
               felling_s_gtharv(1:n_pft) = undef_real
               felling_s_ltharv(1:n_pft) = undef_real
               !---------------------------------------------------------------------------!


               !----- Define the format for the header. -----------------------------------!
               write(hform,fmt='(a,i3.3,a)') '(a',str_len,')'
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Read the top header.                                                  !
               !---------------------------------------------------------------------------!
               read_tophead: do
                  !---- Read line. --------------------------------------------------------!
                  read (unit=12,fmt=hform) cdum
                  poseq = index(cdum,'=')
                  !------------------------------------------------------------------------!

                  !---- If we reached the header, exit the loop. --------------------------!
                  if (poseq == 0) exit read_tophead
                  !------------------------------------------------------------------------!


                  !---- If not, read the variable key for now and remove it from line. ----!
                  vkey = cdum(:poseq-1)
                  cdum = cdum(poseq+1:)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Identify which variable to read (multiple options are for          !
                  ! back compatibility).                                                   !
                  !------------------------------------------------------------------------!
                  select case (trim(vkey))
                  case ('WEST_LONGITUDE','WEST.LONGITUDE')
                     !---- Longitude of the west edge of the domain. ----------------------!
                     read(cdum,fmt=*) wlon
                     !---------------------------------------------------------------------!
                  case ('EAST_LONGITUDE','EAST.LONGITUDE')
                     !---- Longitude of the east edge of the domain. ----------------------!
                     read(cdum,fmt=*) elon
                     !---------------------------------------------------------------------!
                  case ('SOUTH_LATITUDE','SOUTH.LATITUDE')
                     !---- Latitude of the south edge of the domain. ----------------------!
                     read(cdum,fmt=*) slat
                     !---------------------------------------------------------------------!
                  case ('NORTH_LATITUDE','NORTH.LATITUDE')
                     !---- Latitude of the north edge of the domain. ----------------------!
                     read(cdum,fmt=*) nlat
                     !---------------------------------------------------------------------!
                  case ('BLOCK_AREA'    ,'BLOCK.AREA'    )
                     !---- Domain area. ---------------------------------------------------!
                     read(cdum,fmt=*) lu_area
                     !---------------------------------------------------------------------!
                  case ('FIRST_LUYEAR'  ,'FIRST.LUYEAR'  )
                     !---- First year with disturbance data. ------------------------------!
                     read(cdum,fmt=*) yd_1st
                     !---------------------------------------------------------------------!
                  case ('LAST_LUYEAR'   ,'LAST.LUYEAR'   )
                     !---- Last year with disturbance data. -------------------------------!
                     read(cdum,fmt=*) yd_last
                     !---------------------------------------------------------------------!
                  case ('SKID_AREA'     ,'SKID.AREA'     )
                     !---- Damage area due to skid trails and roads. ----------------------!
                     read(cdum,fmt=*) skid_rel_area
                     !---------------------------------------------------------------------!
                  case ('N_PFT_HARVEST' ,'N.PFT.HARVEST' )
                     !---- Number of PFTs to harvest. -------------------------------------!
                     read(cdum,fmt=*) nharvest
                     !---------------------------------------------------------------------!

                     !----- This should precede the PFT-specific instructions. ------------!
                     exit read_tophead
                     !---------------------------------------------------------------------!
                  case default
                     !----- Key is not recognised. Stop the model. ------------------------!
                     write (unit=*,fmt='(a)'        ) '-----------------------------------'
                     write (unit=*,fmt='(a,1x,a)'   ) 'File: ',trim(lu_name)
                     write (unit=*,fmt='(a,2(1x,a))') 'Key:  ',trim(vkey),'is invalid!'
                     write (unit=*,fmt='(a)'        ) '-----------------------------------'
                     call fatal_error(' Problems reading input file (not an ED2 bug)!'     &
                                     ,'read_landuse_matrix','landuse_init.f90')
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!
               end do read_tophead
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !    Before proceeding, ensure that all top header variables have been      !
               ! properly read.                                                            !
               !---------------------------------------------------------------------------!
               if ( wlon     == undef_real    .or. elon     == undef_real    .or.          &
                    slat     == undef_real    .or. nlat     == undef_real    .or.          &
                    lu_area  == undef_real    .or. yd_1st   == undef_integer .or.          &
                    yd_last  == undef_integer .or. nharvest == undef_integer      ) then
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(2(a,1x))')     ' - File: ',trim(lu_name)
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(a,1x,es12.5)') ' - WEST_LONGITUDE = ',wlon
                  write (unit=*,fmt='(a,1x,es12.5)') ' - EAST_LONGITUDE = ',elon
                  write (unit=*,fmt='(a,1x,es12.5)') ' - SOUTH_LATITUDE = ',slat
                  write (unit=*,fmt='(a,1x,es12.5)') ' - NORTH_LATITUDE = ',nlat
                  write (unit=*,fmt='(a,1x,es12.5)') ' - BLOCK_AREA     = ',lu_area
                  write (unit=*,fmt='(a,1x,i6)')     ' - FIRST_LUYEAR   = ',yd_1st
                  write (unit=*,fmt='(a,1x,i6)')     ' - LAST_LUYEAR    = ',yd_last
                  write (unit=*,fmt='(a,1x,es12.5)') ' - SKID_AREA      = ',skid_rel_area
                  write (unit=*,fmt='(a,1x,i6)')     ' - N_PFT_HARVEST  = ',nharvest
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  call fatal_error('Missing variables in input file (not an ED2 bug)!'     &
                                  ,'read_landuse_matrix','landuse_init.f90')
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Decide whether or not to read PFT specific harvest instructions.      !
               !---------------------------------------------------------------------------!
               if (nharvest > 0 ) then
                  !------------------------------------------------------------------------!
                  !     If nharvest is not 0, then the harvesting is not PFT-agnostic,     !
                  ! read the PFT information in the file.  To allow for back-compatibility !
                  ! we must check every variable until we reach the header.                !
                  !------------------------------------------------------------------------!
                  read_harvest: do
                     !---- Read line. -----------------------------------------------------!
                     read (unit=12,fmt=hform) cdum
                     poseq = index(cdum,'=')
                     !---------------------------------------------------------------------!

                     !---- If we reached the header, exit the loop. -----------------------!
                     if (poseq == 0) exit read_harvest
                     !---------------------------------------------------------------------!


                     !---- If not, read the variable key for now and remove it from line. -!
                     vkey = cdum(:poseq-1)
                     cdum = cdum(poseq+1:)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Identify which variable to read (multiple options are for       !
                     ! back compatibility).                                                !
                     !---------------------------------------------------------------------!
                     select case (trim(vkey))
                     case ('HARVEST_PFT','HARVEST.PFT')
                        !---- PFTs to read. -----------------------------------------------!
                        read (cdum, fmt=*) (harvest_pft(h)  ,h=1,nharvest)
                        !------------------------------------------------------------------!

                        !------------------------------------------------------------------!
                        !     Assign defaults to variables in case they have not been      !
                        ! read yet.                                                        !
                        !------------------------------------------------------------------!
                        mindbh_slog      = merge( mindbh_slog     , 0.                     &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  mindbh_slog      /= undef_real         )
                        harvprob_slog    = merge( harvprob_slog   , 1.                     &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  harvprob_slog    /= undef_real         )
                        mindbh_fplt      = merge( mindbh_fplt     , 0.                     &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  mindbh_fplt      /= undef_real         )
                        harvprob_fplt    = merge( harvprob_fplt   , 1.                     &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  harvprob_fplt    /= undef_real         )
                        skid_dbh_thresh  = merge( skid_dbh_thresh , sl_skid_dbh_thresh     &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  skid_dbh_thresh  /= undef_real         )
                        skid_s_gtharv    = merge( skid_s_gtharv   , sl_skid_s_gtharv       &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  skid_s_gtharv    /= undef_real         )
                        skid_s_ltharv    = merge( skid_s_ltharv   , sl_skid_s_ltharv       &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  skid_s_ltharv    /= undef_real         )
                        felling_s_ltharv = merge( felling_s_ltharv, sl_felling_s_ltharv    &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  felling_s_ltharv /= undef_real          )
                        felling_s_gtharv = merge( felling_s_gtharv, 0.                     &
                                                , harvest_pft      /= undef_integer .and.  &
                                                  felling_s_gtharv /= undef_real          )
                        !------------------------------------------------------------------!
                        
                     case ('MINDBH_SLOG','MINDBH.SLOG','MINDBH_1ARY','MINDBH.1ARY')
                        !---- Minimum DBH for harvesting: (selective) logging. ------------!
                        read (cdum, fmt=*) (mindbh_slog     (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('HARVPROB_SLOG','HARVPROB.SLOG','HARVPROB_1ARY','HARVPROB.1ARY')
                        !---- Harvest probability: (selective) logging. -------------------!
                        read (cdum, fmt=*) (harvprob_slog   (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('MINDBH_FPLT','MINDBH.FPLT','MINDBH_2ARY','MINDBH.2ARY')
                        !---- Minimum DBH for harvesting: forest plantation. --------------!
                        read (cdum, fmt=*) (mindbh_fplt     (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('HARVPROB_FPLT','HARVPROB.FPLT','HARVPROB_2ARY','HARVPROB.2ARY')
                        !---- Harvest probability: forest plantation. ---------------------!
                        read (cdum, fmt=*) (harvprob_fplt   (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('SKID_DBH_THRESH','SKID.DBH.THRESH')
                        !---- Skid damage: DBH threshold for small/large tree. ------------!
                        read (cdum, fmt=*) (skid_dbh_thresh (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('SKID_S_GTHARV','SKID.S.GTHARV')
                        !---- Skid damage: Survivorship of large trees. -------------------!
                        read (cdum, fmt=*) (skid_s_gtharv   (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('SKID_S_LTHARV','SKID.S.LTHARV')
                        !---- Skid damage: Survivorship of small trees. -------------------!
                        read (cdum, fmt=*) (skid_s_ltharv   (h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('FELLING_S_LTHARV','FELLING.S.LTHARV')
                        !---- Tree felling: Survivorship of small trees. ------------------!
                        read (cdum, fmt=*) (felling_s_ltharv(h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case ('FELLING_S_GTHARV','FELLING.S.GTHARV')
                        !------------------------------------------------------------------!
                        !    Tree felling: Survivorship of large trees.  Leaving this as   !
                        ! a future option, though there is no strong reason to make this   !
                        ! survivorship anything other than 0.                              !
                        !------------------------------------------------------------------!
                        read (cdum, fmt=*) (felling_s_ltharv(h),h=1,nharvest)
                        !------------------------------------------------------------------!
                     case default
                        !----- Key is not recognised. Stop the model. ---------------------!
                        write (unit=*,fmt='(a)'        ) '--------------------------------'
                        write (unit=*,fmt='(a,1x,a)'   ) 'File: ',trim(lu_name)
                        write (unit=*,fmt='(a,2(1x,a))') 'Key:  ',trim(vkey),'is invalid!'
                        write (unit=*,fmt='(a)'        ) '--------------------------------'
                        call fatal_error(' Problems reading input file (not an ED2 bug)!'  &
                                        ,'read_landuse_matrix','landuse_init.f90')
                        !------------------------------------------------------------------!
                     end select
                     !---------------------------------------------------------------------!
                  end do read_harvest
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !     No specific PFT information was given, this is likely to be a case !
                  ! in which the logging is based on absolute target biomass (PFT- and     !
                  ! DBH-blind).                                                            !
                  !------------------------------------------------------------------------!
                  h = 0
                  nopftloop: do ipft=1,n_pft
                     !--- Skip grasses and lianas. ----------------------------------------!
                     if ( is_grass(ipft) .or. is_liana(ipft) ) cycle nopftloop
                     !---------------------------------------------------------------------!


                     !----- Fill in with default. -----------------------------------------!
                     h                   = h + 1
                     mindbh_slog     (h) = 0.
                     harvprob_slog   (h) = 1.
                     mindbh_fplt     (h) = 0.
                     harvprob_fplt   (h) = 1.
                     skid_dbh_thresh (h) = sl_skid_dbh_thresh
                     skid_s_gtharv   (h) = sl_skid_s_gtharv
                     skid_s_ltharv   (h) = sl_skid_s_ltharv
                     felling_s_ltharv(h) = sl_felling_s_ltharv
                     felling_s_gtharv(h) = 0.
                     !---------------------------------------------------------------------!
                  end do nopftloop
                  !------------------------------------------------------------------------!


                  !--- Skip header and be ready for reading transitions. ------------------!
                  read (unit=12,fmt=*) 
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

               !----- Use file_lat to compute the physical area sampled by the file. ------!
               if (lu_area == 0.) then
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  write (unit=*,fmt='(2(a,1x))')     ' - File:    ',trim(lu_name)
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Wlon:    ',wlon
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Elon:    ',elon
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Slat:    ',slat
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Nlat:    ',nlat
                  write (unit=*,fmt='(a,1x,es12.5)') ' - Lu_area: ',lu_area
                  write (unit=*,fmt='(a,1x,i6)')     ' - Yd_1st:  ',yd_1st
                  write (unit=*,fmt='(a,1x,i6)')     ' - Yd_last: ',yd_last
                  write (unit=*,fmt='(a,1x,i6)')     ' - Nharvest:',nharvest
                  write (unit=*,fmt='(a)')           '------------------------------------'
                  call fatal_error('Land use area is zero, it doesn''t make any sense!'    &
                                  ,'landuse_init','landuse_init.f90')
               else
                  lu_area_i = 1. / lu_area
               end if

               !----- Determine whether this block contains the current polygon. ----------!
               inside = cgrid%lon(ipy) >= wlon .and. cgrid%lon(ipy) <= elon .and.          &
                        cgrid%lat(ipy) >= slat .and. cgrid%lat(ipy) <= nlat


               !---------------------------------------------------------------------------!
               !     Here we will only use the land information if the polygon centre is   !
               ! inside the block of land use disturbance we are about to read.            !
               !---------------------------------------------------------------------------!
               if (inside) then
                  !----- File exists, allocate the maximum number of years. ---------------!
                  allocate(cpoly%clutimes(max_lu_years,cpoly%nsites))


                  !----- Copy the file information to the first site. ---------------------!
                  isi = 1
                  csite => cpoly%site(isi)
                  !------------------------------------------------------------------------!


                  !----- Determine the number of disturbance years. -----------------------!
                  cpoly%num_landuse_years(isi) = max(yd_last,iyearz)-min(yd_1st,iyeara) + 1
                  !------------------------------------------------------------------------!


                  !----- Define the degree of damage relative to felling. -----------------!
                  cpoly%skid_rel_area(isi) = skid_rel_area
                  !------------------------------------------------------------------------!


                  !----- Fill the arrays with the appropriate PFT. ------------------------!
                  select case(cpoly%plantation(isi))
                  case (0)
                     harvloop_slog: do h=1,nharvest
                        ipft = harvest_pft(h)
                        if (ipft >= 1 .and. ipft <= n_pft) then
                           cpoly%mindbh_harvest(ipft,isi) = mindbh_slog  (h)
                           cpoly%prob_harvest  (ipft,isi) = harvprob_slog(h)
                        end if
                     end do harvloop_slog
                  case (1)
                     harvloop_fplt: do h=1,nharvest
                        ipft = harvest_pft(h)
                        if (ipft >= 1 .and. ipft <= n_pft) then
                           cpoly%mindbh_harvest(ipft,isi) = mindbh_fplt  (h)
                           cpoly%prob_harvest  (ipft,isi) = harvprob_fplt(h)
                        end if
                     end do harvloop_fplt
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     These variables are filled the same way for selective logging      !
                  ! and forest plantations.                                                !
                  !------------------------------------------------------------------------!
                  skidloop_site: do h=1,nharvest
                     ipft = harvest_pft(h)
                     if (ipft >= 1 .and. ipft <= n_pft) then
                        cpoly%skid_dbh_thresh (ipft,isi) = skid_dbh_thresh (h)
                        cpoly%skid_s_gtharv   (ipft,isi) = skid_s_gtharv   (h)
                        cpoly%skid_s_ltharv   (ipft,isi) = skid_s_ltharv   (h)
                        cpoly%felling_s_ltharv(ipft,isi) = felling_s_ltharv(h)
                        cpoly%felling_s_gtharv(ipft,isi) = felling_s_gtharv(h)
                     end if
                  end do skidloop_site
                  !------------------------------------------------------------------------!


                  !----- Padding disturbances with zero before first available lu year. ---!
                  iyear = 0
                  do yd_this = iyeara,(yd_1st-1)
                     iyear = iyear + 1
                     clutime => cpoly%clutimes(iyear,isi)
                     clutime%landuse_year            = yd_this
                     clutime%landuse(1:num_lu_trans) = 0.0
                  end do

                  !---- Reading the years that have data ----------------------------------!
                  do yd_this = yd_1st,yd_last
                     iyear = iyear + 1
                     clutime => cpoly%clutimes(iyear,isi)

                     read(unit=12,fmt=*)  clutime%landuse_year                             &
                                        , clutime%landuse(1:num_lu_trans)
                     
                     !---------------------------------------------------------------------!
                     !    Here we normalise by the area, except when landuse(12) and/or    !
                     ! landuse(14) are negative (a special flag to use the selective       !
                     ! logging.                                                            !
                     !---------------------------------------------------------------------!
                     if ( clutime%landuse(12) > 0. ) then
                        clutime%landuse(12) = lu_area_i * clutime%landuse(12)
                     end if
                     if ( clutime%landuse(14) > 0. ) then
                        clutime%landuse(14) = lu_area_i * clutime%landuse(14) 
                     end if
                     clutime%landuse(16) = lu_area_i * clutime%landuse(16)
                     clutime%landuse(18) = lu_area_i * clutime%landuse(18)
                  end do

                  !----- Padding disturbances with zero after last available lu year. -----!
                  do yd_this = (yd_last+1),iyearz
                     iyear = iyear + 1
                     clutime => cpoly%clutimes(iyear,isi)
                     clutime%landuse_year            = yd_this
                     clutime%landuse(1:num_lu_trans) = 0.0
                  end do
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Copy the information from the first site to the other, if they    !
                  ! exist.                                                                 !
                  !------------------------------------------------------------------------!
                  siteloop_one: do isi = 2,cpoly%nsites
                     csite => cpoly%site(isi)

                     !----- Determine the number of disturbance years. --------------------!
                     cpoly%num_landuse_years(isi) = cpoly%num_landuse_years(1)
                     !---------------------------------------------------------------------!


                     !----- Define the degree of damage relative to felling. --------------!
                     cpoly%skid_rel_area(isi) = cpoly%skid_rel_area(1)
                     !---------------------------------------------------------------------!


                     !----- PFT-dependent harvest characteristics. ------------------------!
                     cpoly%mindbh_harvest  (:,isi) = cpoly%mindbh_harvest  (:,1)
                     cpoly%prob_harvest    (:,isi) = cpoly%prob_harvest    (:,1)
                     cpoly%skid_dbh_thresh (:,isi) = cpoly%skid_dbh_thresh (:,1)
                     cpoly%skid_s_gtharv   (:,isi) = cpoly%skid_s_gtharv   (:,1)
                     cpoly%skid_s_ltharv   (:,isi) = cpoly%skid_s_ltharv   (:,1)
                     cpoly%felling_s_ltharv(:,isi) = cpoly%felling_s_ltharv(:,1)
                     cpoly%felling_s_gtharv(:,isi) = cpoly%felling_s_gtharv(:,1)
                     !---------------------------------------------------------------------!



                     !----- Disturbances. -------------------------------------------------!
                     do iyear = 1,cpoly%num_landuse_years(isi)
                        clutime   => cpoly%clutimes(iyear,isi)
                        onelutime => cpoly%clutimes(iyear,1)
                        clutime%landuse_year            = onelutime%landuse_year 
                        clutime%landuse(1:num_lu_trans) = onelutime%landuse(1:num_lu_trans)
                     end do
                     !---------------------------------------------------------------------!

                  end do siteloop_one
               else
                  !------------------------------------------------------------------------!
                  !      No GLU data for this site.  Probably water.                       !
                  !------------------------------------------------------------------------!
                  write (unit=*,fmt='(a)') '----------------------------------------------'
                  write (unit=*,fmt='(a)') ' The closest land use point is too far away...'
                  write (unit=*,fmt='(a)') ' - File:                     ',trim(lu_name)
                  write (unit=*,fmt=fffmt) ' - Polygon longitude:        ',cgrid%lon(ipy)
                  write (unit=*,fmt=fffmt) ' - Polygon latitude:         ',cgrid%lat(ipy)
                  write (unit=*,fmt=fffmt) ' - Closest LU central long.: ',llon_list(ncl)
                  write (unit=*,fmt=fffmt) ' - Closest LU central lat.:  ',llat_list(ncl)
                  write (unit=*,fmt=fffmt) ' - Closest LU western edge:  ',wlon
                  write (unit=*,fmt=fffmt) ' - Closest LU eastern edge:  ',elon
                  write (unit=*,fmt=fffmt) ' - Closest LU southern edge: ',slat
                  write (unit=*,fmt=fffmt) ' - Closest LU northern edge: ',nlat
                  write (unit=*,fmt=esfmt) ' - Distance:                 ',file_ldist(ncl)
                  write (unit=*,fmt=*)     ' '
                  write (unit=*,fmt='(a)') ' We will assign no land use disturbance rate.'
                  write (unit=*,fmt='(a)') '----------------------------------------------'

                  !----- Allocate just 1 landuse year. ------------------------------------!
                  allocate(cpoly%clutimes(1,cpoly%nsites))

                  !------------------------------------------------------------------------!
                  !    Set the parameters in a way that no logging/ploughing will happen.  !
                  !------------------------------------------------------------------------!
                  do isi = 1,cpoly%nsites
                     cpoly%clutimes(1,isi)%landuse_year            = iyeara
                     cpoly%clutimes(1,isi)%landuse(1:num_lu_trans) = 0.0
                     cpoly%skid_rel_area     (isi)                 = 0.
                     cpoly%mindbh_harvest  (:,isi)                 = huge_dbh
                     cpoly%prob_harvest    (:,isi)                 = 0.
                     cpoly%skid_dbh_thresh (:,isi)                 = huge_dbh
                     cpoly%skid_s_gtharv   (:,isi)                 = 1.
                     cpoly%skid_s_ltharv   (:,isi)                 = 1.
                     cpoly%felling_s_ltharv(:,isi)                 = 1.
                     cpoly%felling_s_gtharv(:,isi)                 = 1.
                  end do
               end if
               !---------------------------------------------------------------------------!


               !----- Close the land use file, outside the if statement. ------------------!
               close(unit=12,status='keep')
               !---------------------------------------------------------------------------!

            case (2)
               !---------------------------------------------------------------------------!
               !      Make the land use data based on ED2IN.                               !
               !      Work with the first site, then copy the data to the others.          !
               !---------------------------------------------------------------------------!
               isi   = 1
               csite => cpoly%site(isi)
               !---------------------------------------------------------------------------!


               !----- No plantations. -----------------------------------------------------!
               cpoly%plantation(:) = 0
               !---------------------------------------------------------------------------!



               !----- Determine the number of disturbance years. --------------------------!
               cpoly%num_landuse_years(isi) = sim_years
               !---------------------------------------------------------------------------!



               !----- File exists, allocate the maximum number of years. ------------------!
               allocate(cpoly%clutimes(sim_years,cpoly%nsites))
               !---------------------------------------------------------------------------!


               !----- Define the degree of damage relative to felling. --------------------!
               cpoly%skid_rel_area(isi) = cpoly%skid_rel_area(1)
               !---------------------------------------------------------------------------!


               !----- Initialise the PFT-dependent arrays. --------------------------------!
               cpoly%mindbh_harvest  (:,isi) = huge_dbh
               cpoly%prob_harvest    (:,isi) = 0.
               cpoly%skid_dbh_thresh (:,isi) = merge(huge_dbh,sl_skid_dbh_thresh           &
                                                                              ,is_grass(:))
               cpoly%skid_s_gtharv   (:,isi) = merge( 1.00,sl_skid_s_gtharv   ,is_grass(:))
               cpoly%skid_s_ltharv   (:,isi) = merge( 1.00,sl_skid_s_ltharv   ,is_grass(:))
               cpoly%felling_s_gtharv(:,isi) = merge( 0.70,            0.00   ,is_grass(:))
               cpoly%felling_s_ltharv(:,isi) = merge( 0.70,sl_felling_s_ltharv,is_grass(:))
               !---------------------------------------------------------------------------!



               !------ Find the number of PFT that can be harvested. ----------------------!
               nharvest = count(sl_pft >= 1 .and. sl_pft <= n_pft)
               !---------------------------------------------------------------------------!

               !----- Fill the arrays with the appropriate PFT. ---------------------------!
               harvloop_two: do h=1,nharvest
                  ipft = sl_pft(h)
                  if (ipft >= 1 .and. ipft <= n_pft) then
                     cpoly%mindbh_harvest(ipft,isi) = sl_mindbh_harvest(h)
                     cpoly%prob_harvest  (ipft,isi) = sl_prob_harvest  (h)
                  end if
               end do harvloop_two
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Fill in the disturbance matrices and biomass target.                 !
               !---------------------------------------------------------------------------!
               iyear = 0
               do yd_this = iyeara,iyearz
                  iyear = iyear + 1
                  clutime => cpoly%clutimes(iyear,isi)

                  clutime%landuse_year            = yd_this
                  clutime%landuse(1:num_lu_trans) = 0.

                  !------------------------------------------------------------------------!
                  !     Decide whether to include logging disturbance in this year.        !
                  !------------------------------------------------------------------------!
                  if (yd_this >= sl_yr_first) then
                     if ( (sl_scale == 1) .or. (mod(yd_this-sl_yr_first,sl_nyrs) == 0) )   &
                     then
                        clutime%landuse(11) = lnexp_max
                        clutime%landuse(12) = -1.0
                        clutime%landuse(14) = -1.0
                     end if
                  end if
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Copy the information from the first site to the others.              !
               !---------------------------------------------------------------------------!
               siteloop_two: do isi = 2,cpoly%nsites
                  csite => cpoly%site(isi)

                  !----- Determine the number of disturbance years. -----------------------!
                  cpoly%num_landuse_years(isi) = cpoly%num_landuse_years(1)
                  !------------------------------------------------------------------------!


                  !----- Define the degree of damage relative to felling. -----------------!
                  cpoly%skid_rel_area(isi) = cpoly%skid_rel_area(1)
                  !------------------------------------------------------------------------!


                  !----- PFT-dependent harvest characteristics. ---------------------------!
                  cpoly%mindbh_harvest  (:,isi) = cpoly%mindbh_harvest  (:,1)
                  cpoly%prob_harvest    (:,isi) = cpoly%prob_harvest    (:,1)
                  cpoly%skid_dbh_thresh (:,isi) = cpoly%skid_dbh_thresh (:,1)
                  cpoly%skid_s_gtharv   (:,isi) = cpoly%skid_s_gtharv   (:,1)
                  cpoly%skid_s_ltharv   (:,isi) = cpoly%skid_s_ltharv   (:,1)
                  cpoly%felling_s_gtharv(:,isi) = cpoly%felling_s_gtharv(:,1)
                  cpoly%felling_s_ltharv(:,isi) = cpoly%felling_s_ltharv(:,1)
                  !------------------------------------------------------------------------!



                  !----- Disturbances. ----------------------------------------------------!
                  do iyear = 1,cpoly%num_landuse_years(isi)
                     clutime   => cpoly%clutimes(iyear,isi)
                     onelutime => cpoly%clutimes(iyear,1)
                     clutime%landuse_year            = onelutime%landuse_year 
                     clutime%landuse(1:num_lu_trans) = onelutime%landuse(1:num_lu_trans)
                  end do
                  !------------------------------------------------------------------------!

               end do siteloop_two
            end select
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !       Write a table somewhat similar to lu file.                             !
            !------------------------------------------------------------------------------!
            if (write_lu_settings) then
               cpoly   => cgrid%polygon(ipy)
               isi     = 1
            
               wlon    = 0.1 * real(nint(10.*cgrid%lon(ipy))) - 0.5
               elon    = 0.1 * real(nint(10.*cgrid%lon(ipy))) + 0.5
               slat    = 0.1 * real(nint(10.*cgrid%lat(ipy))) - 0.5
               nlat    = 0.1 * real(nint(10.*cgrid%lat(ipy))) + 0.5
               lu_area = solid_area(wlon,slat,elon,nlat)
               
               !---------------------------------------------------------------------------!
               !      Find which PFTs to harvest.  If none is provided, we skip this part  !
               ! of the header.                                                            !
               !---------------------------------------------------------------------------!
               nharvest = count(cpoly%prob_harvest(:,isi) > 0.0)
               write(hifmt,fmt='(a,i2.2,a)') '(a,',nharvest,'(i2,1x))'
               write(hffmt,fmt='(a,i2.2,a)') '(a,',nharvest,'(f8.3,1x))'
               h = 0
               do ipft=1,n_pft
                  if (cpoly%prob_harvest(ipft,isi) > 0.0) then
                     h = h + 1
                     harvest_pft     (h) = ipft
                     mindbh_slog     (h) = cpoly%mindbh_harvest  (ipft,isi)
                     harvprob_slog   (h) = cpoly%prob_harvest    (ipft,isi)
                     mindbh_fplt     (h) = cpoly%mindbh_harvest  (ipft,isi)
                     harvprob_fplt   (h) = cpoly%prob_harvest    (ipft,isi)
                     skid_dbh_thresh (h) = cpoly%skid_dbh_thresh (ipft,isi)
                     skid_s_gtharv   (h) = cpoly%skid_s_gtharv   (ipft,isi)
                     skid_s_ltharv   (h) = cpoly%skid_s_ltharv   (ipft,isi)
                     felling_s_ltharv(h) = cpoly%felling_s_ltharv(ipft,isi)
                     felling_s_gtharv(h) = cpoly%felling_s_gtharv(ipft,isi)
                  end if
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Write the LU-like file.                                              !
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Write the LU-like file.                                              !
               !---------------------------------------------------------------------------!
               open(unit=16,file=trim(lu_table),status='replace',action='write')
               write(unit=16,fmt='(a,1x,f8.3)' ) 'WEST_LONGITUDE = ',wlon
               write(unit=16,fmt='(a,1x,f8.3)' ) 'EAST_LONGITUDE = ',elon
               write(unit=16,fmt='(a,1x,f8.3)' ) 'SOUTH_LATITUDE = ',slat
               write(unit=16,fmt='(a,1x,f8.3)' ) 'NORTH_LATITUDE = ',nlat
               write(unit=16,fmt='(a,1x,f20.5)') 'BLOCK_AREA     = ',lu_area
               write(unit=16,fmt='(a,1x,i4.4)' ) 'FIRST_LUYEAR   = ',iyeara
               write(unit=16,fmt='(a,1x,i4.4)' ) 'LAST_LUYEAR    = ',iyearz
               write(unit=16,fmt='(a,1x,f8.3)' ) 'SKID_AREA      = ',skid_rel_area
               write(unit=16,fmt='(a,1x,i2)'   ) 'N_PFT_HARVEST  = ',nharvest
               if (nharvest > 0) then
                  write(unit=16,fmt=hifmt)                                                 &
                                   'HARVEST_PFT      = ',(harvest_pft     (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'MINDBH_SLOG      = ',(mindbh_slog     (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'HARVPROB_SLOG    = ',(harvprob_slog   (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'MINDBH_FPLT      = ',(mindbh_fplt     (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'HARVPROB_FPLT    = ',(harvprob_fplt   (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'SKID_DBH_THRESH  = ',(skid_dbh_thresh (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'SKID_S_GTHARV    = ',(skid_s_gtharv   (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'SKID_S_LTHARV    = ',(skid_s_ltharv   (h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'FELLING_S_LTHARV = ',(felling_s_ltharv(h),h=1,nharvest)
                  write(unit=16,fmt=hffmt)                                                 &
                                   'FELLING_S_GTHARV = ',(felling_s_gtharv(h),h=1,nharvest)
               end if
               write(unit=16,fmt='(a,19(1x,a))')  'YEAR','     CPL_PST','     PST_CPL'     &
                          ,'     PST_VEG','     VEG_PST','     VEG_CLP','     CPL_VEG'     &
                          ,'     SEC_CPL','     CPL_SEC','     SEC_PST','     PST_SEC'     &
                          ,'     VEG_SEC','  BT_MAT_SEC','  FL_MAT_SEC','  BT_MAT_VEG'     &
                          ,'  FL_MAT_VEG','  BT_YNG_SEC','  FL_YNG_SEC','  BT_YNG_VEG'     &
                          ,'  FL_YNG_VEG'
               !----- Disturbances. -------------------------------------------------------!
               do iyear = 1,cpoly%num_landuse_years(isi)
                  clutime     => cpoly%clutimes(iyear,isi)
                  landuse_now =  clutime%landuse
                  if (landuse_now(12) > 0.) landuse_now(12) = lu_area * landuse_now(12)
                  if (landuse_now(14) > 0.) landuse_now(14) = lu_area * landuse_now(14)
                  landuse_now(16) = lu_area * landuse_now(16)
                  landuse_now(18) = lu_area * landuse_now(18)

                  write(unit=16,fmt='(i4.4,19(1x,es12.5))')                                &
                                     clutime%landuse_year,(landuse_now(h),h=1,num_lu_trans)
                  !------------------------------------------------------------------------!
               end do
               close(unit=16,status='keep')
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do polyloop
         !---------------------------------------------------------------------------------!
      end do gridloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine read_landuse_matrix
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine reads the plantation fraction.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine read_plantation_fractions(cpoly,polylon,polylat,igr)
      use ed_state_vars , only : polygontype         ! ! structure
      use ed_max_dims   , only : str_len             ! ! intent(in)
      use disturb_coms  , only : plantation_file     & ! intent(in)
                               , max_plantation_dist & ! intent(in)
                               , min_plantation_frac ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(polygontype)     , target      :: cpoly
      real                  , intent(in)  :: polylon
      real                  , intent(in)  :: polylat
      integer               , intent(in)  :: igr
      !----- Local variables --------------------------------------------------------------!
      logical                             :: exans
      integer                             :: ierr
      integer                             :: isi
      real, dimension(:)    , allocatable :: plantlon
      real, dimension(:)    , allocatable :: plantlat
      real, dimension(:)    , allocatable :: plantdist
      real, dimension(:)    , allocatable :: fracplant
      real                                :: rdum
      integer                             :: ndat
      integer                             :: n
      !----- Local constants. -------------------------------------------------------------!
      character(len=12)    , parameter    :: fffmt='(a,1x,f12.5)'
      character(len=13)    , parameter    :: esfmt='(a,1x,es12.5)'
      !----- External functions. ----------------------------------------------------------!
      real                  , external    :: dist_gc
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     If the user left the plantation file empty, it means that they don't want to   !
      ! use plantation files, skip the subroutine and don't print any warnings.            !
      !------------------------------------------------------------------------------------!
      if (len_trim(plantation_file(igr)) == 0) then
         return
      else
         !----- Check whether plantation file exists. -------------------------------------!
         inquire(file=trim(plantation_file(igr)),exist=exans)
         if (.not.exans)then
            write (unit=*,fmt='(a)') '----------------------------------------------------'
            write (unit=*,fmt='(a)') 'File :'//trim(plantation_file(igr))//' not found...'
            write (unit=*,fmt='(a)') 'Assuming that there are no plantations'
            write (unit=*,fmt='(a)') '----------------------------------------------------'
            return
         end if
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If we reach this point, then there is a plantation file. Open it.              !
      !------------------------------------------------------------------------------------!
      open (unit=12, file=trim(plantation_file(igr)), form='formatted', status='old')
      !----- The first loop will determine how many points are available. -----------------!
      ndat = 0
      countloop: do
         read (unit=12,fmt=*,iostat=ierr) rdum
         if (ierr /= 0) exit countloop
         
         ndat = ndat + 1
      end do countloop
      rewind(unit=12)

      !----- Now that we know the number of lines, allocate the temporary vectors. --------!
      allocate (plantlon(ndat),plantlat(ndat),plantdist(ndat),fracplant(ndat))

      !----- Now we read the data. --------------------------------------------------------!
      read_plantation: do n=1,ndat
         read (unit=12,fmt=*) plantlat(n), plantlon(n), fracplant(n)
         plantdist(n) = dist_gc(polylon,plantlon(n),polylat,plantlat(n))
      end do read_plantation
      close(unit=12,status='keep')

      !----- Determine which point was the closest one. -----------------------------------!
      n = minloc(plantdist,dim=1)

      !----- Check whether the distance is not too large.  If it is, don't use it. --------!
      if (plantdist(n) > max_plantation_dist) then
         write (unit=*,fmt='(a)') '--------------------------------------------------------'
         write (unit=*,fmt=*)     ' The closest plantation point is too far away...'
         write (unit=*,fmt=fffmt) ' - Polygon longitude:                   ',polylon
         write (unit=*,fmt=fffmt) ' - Polygon latitude:                    ',polylat
         write (unit=*,fmt=fffmt) ' - Plantation closest longitude:        ',plantlon(n)
         write (unit=*,fmt=fffmt) ' - Plantation closest latitude:         ',plantlat(n)
         write (unit=*,fmt=esfmt) ' - Distance from polygon to plantation: ',plantdist(n)
         write (unit=*,fmt=esfmt) ' - Maximum accepted distance:               '           &
                                     ,max_plantation_dist
         write (unit=*,fmt=*)     ' '
         write (unit=*,fmt='(a)') ' We will assume no plantations in this polygon.'
         write (unit=*,fmt='(a)') '-------------------------------------------------------'

      elseif (fracplant(n) > min_plantation_frac) then
         !----- Assume that all sites are plantation sites. -------------------------------!
         do isi = 1,cpoly%nsites
            cpoly%plantation(isi) = 1
         end do
      end if
      !------------------------------------------------------------------------------------!

      deallocate(plantlon,plantlat,fracplant,plantdist)

      return
   end subroutine read_plantation_fractions
   !=======================================================================================!
   !=======================================================================================!

end module landuse_init
