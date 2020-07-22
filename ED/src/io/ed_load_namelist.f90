!==========================================================================================!
!==========================================================================================!
subroutine read_nl(namelist_name)
   use ename_coms, only : nl              & ! intent(inout)
                        , init_ename_vars ! ! subroutine
  
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*), intent(in) :: namelist_name
   !----- Local variables. ----------------------------------------------------------------!
   logical                      :: fexists
   !----- Name lists. ---------------------------------------------------------------------!
   namelist /ED_NL/ nl
   !---------------------------------------------------------------------------------------!

   !----- Open the namelist file. ---------------------------------------------------------!
   inquire (file=trim(namelist_name),exist=fexists)
   if (.not. fexists) then
      call fatal_error('The namelist file '//trim(namelist_name)//' is missing.'           &
                      ,'read_nl','ed_load_namelist.f90')
   end if

   !----- Initialise the name list with absurd, undefined values. -------------------------!
   call init_ename_vars(nl) 
  
   !----- Read grid point and options information from the namelist. ----------------------!
   open (unit=10, status='OLD', file=namelist_name)
   read (unit=10, nml=ED_NL)
   close(unit=10)

   return
end subroutine read_nl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine copy_nl(copy_type)

   use ed_max_dims          , only : n_pft                     & ! intent(in)
                                   , nzgmax                    & ! intent(in)
                                   , undef_integer             & ! intent(in)
                                   , skip_integer              & ! intent(in)
                                   , skip_real                 & ! intent(in)
                                   , maxgrds                   ! ! intent(in)
   use ename_coms           , only : nl                        ! ! intent(in)
   use soil_coms            , only : find_soil_class           & ! function
                                   , isoilflg                  & ! intent(out)
                                   , islcolflg                 & ! intent(out)
                                   , nslcon                    & ! intent(out)
                                   , isoilcol                  & ! intent(out)
                                   , slxclay                   & ! intent(out)
                                   , slxsand                   & ! intent(out)
                                   , slmstr                    & ! intent(out)
                                   , stgoff                    & ! intent(out)
                                   , zrough                    & ! intent(out)
                                   , soil_database             & ! intent(out)
                                   , isoilstateinit            & ! intent(out)
                                   , isoildepthflg             & ! intent(out)
                                   , isoilbc                   & ! intent(out)
                                   , sldrain                   & ! intent(out)
                                   , soilstate_db              & ! intent(out)
                                   , soildepth_db              & ! intent(out)
                                   , runoff_time               & ! intent(out)
                                   , slz                       & ! intent(out)
                                   , veg_database              ! ! intent(out)
   use met_driver_coms      , only : ed_met_driver_db          & ! intent(out)
                                   , ishuffle                  & ! intent(out)
                                   , metcyc1                   & ! intent(out)
                                   , metcycf                   & ! intent(out)
                                   , imettype                  & ! intent(out)
                                   , imetavg                   & ! intent(out)
                                   , imetrad                   & ! intent(out)
                                   , initial_co2               & ! intent(out)
                                   , lapse_scheme              ! ! intent(out)
   use mem_polygons         , only : n_poi                     & ! intent(out)
                                   , poi_lat                   & ! intent(out)
                                   , poi_lon                   & ! intent(out)
                                   , poi_res                   & ! intent(out)
                                   , n_ed_region               & ! intent(out)
                                   , ed_reg_latmin             & ! intent(out)
                                   , ed_reg_latmax             & ! intent(out)
                                   , ed_reg_lonmin             & ! intent(out)
                                   , ed_reg_lonmax             & ! intent(out)
                                   , grid_res                  & ! intent(out)
                                   , grid_type                 & ! intent(out)
                                   , edres                     & ! intent(out)
                                   , maxsite                   & ! intent(out)
                                   , maxpatch                  & ! intent(out)
                                   , maxcohort                 ! ! intent(out)
   use physiology_coms      , only : iphysiol                  & ! intent(out)
                                   , h2o_plant_lim             & ! intent(out)
                                   , plant_hydro_scheme        & ! intent(out)
                                   , istomata_scheme           & ! intent(out)
                                   , istruct_growth_scheme     & ! intent(out)
                                   , istem_respiration_scheme  & ! intent(out)
                                   , trait_plasticity_scheme   & ! intent(out)
                                   , iddmort_scheme            & ! intent(out)
                                   , cbr_scheme                & ! intent(out)
                                   , ddmort_const              & ! intent(out)
                                   , carbon_mortality_scheme   & ! intent(out)
                                   , hydraulic_mortality_scheme& ! intent(out)
                                   , n_plant_lim               & ! intent(out)
                                   , vmfact_c3                 & ! intent(out)
                                   , vmfact_c4                 & ! intent(out)
                                   , mphoto_trc3               & ! intent(out)
                                   , mphoto_tec3               & ! intent(out)
                                   , mphoto_c4                 & ! intent(out)
                                   , bphoto_blc3               & ! intent(out)
                                   , bphoto_nlc3               & ! intent(out)
                                   , bphoto_c4                 & ! intent(out)
                                   , kw_grass                  & ! intent(out)
                                   , kw_tree                   & ! intent(out)
                                   , gamma_c3                  & ! intent(out)
                                   , gamma_c4                  & ! intent(out)
                                   , d0_grass                  & ! intent(out)
                                   , d0_tree                   & ! intent(out)
                                   , alpha_c3                  & ! intent(out)
                                   , alpha_c4                  & ! intent(out)
                                   , klowco2in                 & ! intent(out)
                                   , rrffact                   & ! intent(out)
                                   , growthresp                & ! intent(out)
                                   , q10_c3                    & ! intent(out)
                                   , q10_c4                    & ! intent(out)
                                   , quantum_efficiency_T      ! ! intent(out)
   use phenology_coms       , only : iphen_scheme              & ! intent(out)
                                   , iphenys1                  & ! intent(out)
                                   , iphenysf                  & ! intent(out)
                                   , iphenyf1                  & ! intent(out)
                                   , iphenyff                  & ! intent(out)
                                   , phenpath                  & ! intent(out)
                                   , repro_scheme              & ! intent(out)
                                   , radint                    & ! intent(out)
                                   , radslp                    & ! intent(out)
                                   , thetacrit                 ! ! intent(out)
   use decomp_coms          , only : n_decomp_lim              & ! intent(out)
                                   , decomp_scheme             ! ! intent(out)
   use disturb_coms         , only : include_fire              & ! intent(out)
                                   , fire_parameter            & ! intent(out)
                                   , ianth_disturb             & ! intent(out)
                                   , treefall_disturbance_rate & ! intent(out)
                                   , lu_database               & ! intent(out)
                                   , plantation_file           & ! intent(out)
                                   , lu_rescale_file           & ! intent(out)
                                   , sm_fire                   & ! intent(out)
                                   , time2canopy               & ! intent(out)
                                   , min_patch_area            & ! intent(out)
                                   , sl_scale                  & ! intent(out)
                                   , sl_yr_first               & ! intent(out)
                                   , sl_nyrs                   & ! intent(out)
                                   , sl_pft                    & ! intent(out)
                                   , sl_prob_harvest           & ! intent(out)
                                   , sl_mindbh_harvest         & ! intent(out)
                                   , sl_biomass_harvest        & ! intent(out)
                                   , sl_skid_rel_area          & ! intent(out)
                                   , sl_skid_s_gtharv          & ! intent(out)
                                   , sl_skid_s_ltharv          & ! intent(out)
                                   , sl_felling_s_ltharv       & ! intent(out)
                                   , cl_fseeds_harvest         & ! intent(out)
                                   , cl_fstorage_harvest       & ! intent(out)
                                   , cl_fleaf_harvest          ! ! intent(out)
   use pft_coms             , only : include_these_pft         & ! intent(out)
                                   , pasture_stock             & ! intent(out)
                                   , agri_stock                & ! intent(out)
                                   , plantation_stock          & ! intent(out)
                                   , pft_1st_check             ! ! intent(out)
   use ed_misc_coms         , only : expnme                    & ! intent(out)
                                   , runtype                   & ! intent(out)
                                   , itimez                    & ! intent(out)
                                   , idatez                    & ! intent(out)
                                   , imonthz                   & ! intent(out)
                                   , iyearz                    & ! intent(out)
                                   , itimea                    & ! intent(out)
                                   , idatea                    & ! intent(out)
                                   , imontha                   & ! intent(out)
                                   , iyeara                    & ! intent(out)
                                   , itimeh                    & ! intent(out)
                                   , idateh                    & ! intent(out)
                                   , imonthh                   & ! intent(out)
                                   , iyearh                    & ! intent(out)
                                   , ifoutput                  & ! intent(out)
                                   , iclobber                  & ! intent(out)
                                   , frqfast                   & ! intent(out)
                                   , ndcycle                   & ! intent(out)
                                   , sfilin                    & ! intent(out)
                                   , ied_init_mode             & ! intent(out)
                                   , current_time              & ! intent(out)
                                   , thsums_database           & ! intent(out)
                                   , end_time                  & ! intent(out)
                                   , radfrq                    & ! intent(out)
                                   , ivegt_dynamics            & ! intent(out)
                                   , ibigleaf                  & ! intent(out)
                                   , integration_scheme        & ! intent(out)
                                   , nsub_euler                & ! intent(out)
                                   , ffilout                   & ! intent(out)
                                   , idoutput                  & ! intent(out)
                                   , imoutput                  & ! intent(out)
                                   , iyoutput                  & ! intent(out)
                                   , iqoutput                  & ! intent(out)
                                   , itoutput                  & ! intent(out)
                                   , iooutput                  & ! intent(out)
                                   , igoutput                  & ! intent(out)
                                   , obstime_db                & ! intent(out)
                                   , dtlsm                     & ! intent(out)
                                   , month_yrstep              & ! intent(out)
                                   , frqstate                  & ! intent(out)
                                   , sfilout                   & ! intent(out)
                                   , isoutput                  & ! intent(out)
                                   , gfilout                   & ! intent(out)
                                   , iadd_site_means           & ! intent(out)
                                   , iadd_patch_means          & ! intent(out)
                                   , iadd_cohort_means         & ! intent(out)
                                   , iprintpolys               & ! intent(out)
                                   , printvars                 & ! intent(out)
                                   , npvars                    & ! intent(out)
                                   , pfmtstr                   & ! intent(out)
                                   , ipmax                     & ! intent(out)
                                   , ipmin                     & ! intent(out)
                                   , iedcnfgf                  & ! intent(out)
                                   , outfast                   & ! intent(out)
                                   , outstate                  & ! intent(out)
                                   , unitfast                  & ! intent(out)
                                   , unitstate                 & ! intent(out)
                                   , event_file                & ! intent(out)
                                   , iallom                    & ! intent(out)
                                   , economics_scheme          & ! intent(out)
                                   , igrass                    & ! intent(out)
                                   , min_site_area             & ! intent(out)
                                   , fast_diagnostics          & ! intent(out)
                                   , writing_dail              & ! intent(out)
                                   , writing_mont              & ! intent(out)
                                   , writing_dcyc              & ! intent(out)
                                   , writing_year              & ! intent(out)
                                   , writing_long              & ! intent(out) 
                                   , writing_eorq              & ! intent(out)
                                   , history_fast              & ! intent(out) 
                                   , history_dail              & ! intent(out) 
                                   , history_eorq              & ! intent(out)
                                   , growth_resp_scheme        & ! intent(out)
                                   , storage_resp_scheme       & ! intent(out)
                                   , attach_metadata           ! ! intent(out)
   use grid_coms            , only : time                      & ! intent(out)
                                   , centlon                   & ! intent(out)
                                   , centlat                   & ! intent(out)
                                   , deltax                    & ! intent(out)
                                   , deltay                    & ! intent(out)
                                   , nnxp                      & ! intent(out)
                                   , nnyp                      & ! intent(out)
                                   , nstratx                   & ! intent(out)
                                   , nstraty                   & ! intent(out)
                                   , polelat                   & ! intent(out)
                                   , polelon                   & ! intent(out)
                                   , ngrids                    & ! intent(out)
                                   , timmax                    & ! intent(out)
                                   , time                      & ! intent(out)
                                   , nzg                       & ! intent(out)
                                   , nzs                       ! ! intent(out)
   use canopy_air_coms      , only : icanturb                  & ! intent(out)
                                   , isfclyrm                  & ! intent(out)
                                   , ied_grndvap               & ! intent(out)
                                   , ubmin                     & ! intent(out)
                                   , ugbmin                    & ! intent(out)
                                   , ustmin                    & ! intent(out)
                                   , gamm                      & ! intent(out)
                                   , gamh                      & ! intent(out)
                                   , tprandtl                  & ! intent(out)
                                   , ribmax                    & ! intent(out)
                                   , lwidth_grass              & ! intent(out)
                                   , lwidth_bltree             & ! intent(out)
                                   , lwidth_nltree             & ! intent(out)
                                   , leaf_maxwhc               ! ! intent(out)
   use canopy_layer_coms    , only : crown_mod                 ! ! intent(out)
   use canopy_radiation_coms, only : icanrad                   & ! intent(out)
                                   , ihrzrad                   & ! intent(out)
                                   , ltrans_vis                & ! intent(out)
                                   , ltrans_nir                & ! intent(out)
                                   , lreflect_vis              & ! intent(out)
                                   , lreflect_nir              & ! intent(out)
                                   , orient_tree               & ! intent(out)
                                   , orient_grass              & ! intent(out)
                                   , clump_tree                & ! intent(out)
                                   , clump_grass               ! ! intent(out)
   use rk4_coms             , only : ibranch_thermo            & ! intent(out)
                                   , ipercol                   & ! intent(out)
                                   , rk4_tolerance             ! ! intent(out)
   use ed_para_coms         , only : loadmeth                  ! ! intent(out)
   use detailed_coms        , only : dt_census                 & ! intent(out)
                                   , yr1st_census              & ! intent(out)
                                   , mon1st_census             & ! intent(out)
                                   , min_recruit_dbh           & ! intent(out)
                                   , idetailed                 & ! intent(out)
                                   , patch_keep                ! ! intent(out)
   use consts_coms          , only : vonk                      & ! intent(in)
                                   , day_sec                   & ! intent(in)
                                   , hr_sec                    & ! intent(in)
                                   , min_sec                   ! ! intent(in)
   use fusion_fission_coms  , only : ifusion                   ! ! intent(out)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*), intent(in) :: copy_type
   !----- Internal variables. -------------------------------------------------------------!
   integer                      :: ifm
   integer, dimension(n_pft)    :: idx
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Here we decide which variables we should copy based on the input variable.        !
   !---------------------------------------------------------------------------------------!
   select case (trim(copy_type))
   case ('ALL_CASES')
      !------------------------------------------------------------------------------------!
      !     The namelist variables in this section will always be used from the namelist,  !
      ! regardless of which type of run this is (a history start or not).  This allows     !
      ! model options to be changed if it is a history start.                              !
      !------------------------------------------------------------------------------------!
      expnme                    = nl%expnme
      runtype                   = nl%runtype
      loadmeth                  = nl%loadmeth

      itimez                    = nl%itimez
      idatez                    = nl%idatez
      imonthz                   = nl%imonthz
      iyearz                    = nl%iyearz
      dtlsm                     = nl%dtlsm
      radfrq                    = nl%radfrq
      month_yrstep              = nl%month_yrstep

      ifoutput                  = nl%ifoutput
      idoutput                  = nl%idoutput
      imoutput                  = nl%imoutput
      iqoutput                  = nl%iqoutput
      iyoutput                  = nl%iyoutput
      itoutput                  = nl%itoutput
      iooutput                  = nl%iooutput
      isoutput                  = nl%isoutput

      iadd_site_means           = nl%iadd_site_means
      iadd_patch_means          = nl%iadd_patch_means
      iadd_cohort_means         = nl%iadd_cohort_means

      attach_metadata           = nl%attach_metadata

      iclobber                  = nl%iclobber
      unitfast                  = nl%unitfast
      unitstate                 = nl%unitstate
      frqfast                   = nl%frqfast
      frqstate                  = nl%frqstate
      outfast                   = nl%outfast
      outstate                  = nl%outstate
      
      sfilin                    = nl%sfilin

      itimeh                    = nl%itimeh
      idateh                    = nl%idateh
      imonthh                   = nl%imonthh
      iyearh                    = nl%iyearh
      
      ffilout                   = nl%ffilout
      sfilout                   = nl%sfilout
      ied_init_mode             = nl%ied_init_mode

      isoilflg                  = nl%isoilflg
      islcolflg                 = nl%islcolflg
      nslcon                    = nl%nslcon
      isoilcol                  = nl%isoilcol
      slxclay                   = nl%slxclay
      slxsand                   = nl%slxsand
      slmstr(1:nzgmax)          = nl%slmstr(1:nzgmax)
      stgoff(1:nzgmax)          = nl%stgoff(1:nzgmax)

      soil_database             = nl%soil_database
      veg_database              = nl%veg_database
      lu_database               = nl%lu_database
      plantation_file           = nl%plantation_file
      lu_rescale_file           = nl%lu_rescale_file
      thsums_database           = nl%thsums_database

      ed_met_driver_db          = nl%ed_met_driver_db
      obstime_db                = nl%obstime_db
      soilstate_db              = nl%soilstate_db
      soildepth_db              = nl%soildepth_db
      isoilstateinit            = nl%isoilstateinit
      isoildepthflg             = nl%isoildepthflg
      isoilbc                   = nl%isoilbc
      sldrain                   = nl%sldrain

      n_poi                     = nl%n_poi
      n_ed_region               = nl%n_ed_region
      grid_res                  = nl%grid_res
      grid_type                 = nl%grid_type
      poi_lat                   = nl%poi_lat
      poi_lon                   = nl%poi_lon
      poi_res                   = nl%poi_res
      ed_reg_latmin             = nl%ed_reg_latmin
      ed_reg_latmax             = nl%ed_reg_latmax
      ed_reg_lonmin             = nl%ed_reg_lonmin
      ed_reg_lonmax             = nl%ed_reg_lonmax

      ivegt_dynamics            = nl%ivegt_dynamics
      ibigleaf                  = nl%ibigleaf
      integration_scheme        = nl%integration_scheme
      nsub_euler                = nl%nsub_euler
      rk4_tolerance             = nl%rk4_tolerance
      ibranch_thermo            = nl%ibranch_thermo
      iphysiol                  = nl%iphysiol
      iallom                    = nl%iallom
      economics_scheme          = nl%economics_scheme
      igrass                    = nl%igrass
      iphen_scheme              = nl%iphen_scheme
      repro_scheme              = nl%repro_scheme
      lapse_scheme              = nl%lapse_scheme
      crown_mod                 = nl%crown_mod
      icanrad                   = nl%icanrad
      ihrzrad                   = nl%ihrzrad
      ltrans_vis                = nl%ltrans_vis
      ltrans_nir                = nl%ltrans_nir
      lreflect_vis              = nl%lreflect_vis
      lreflect_nir              = nl%lreflect_nir
      orient_tree               = nl%orient_tree
      orient_grass              = nl%orient_grass
      clump_tree                = nl%clump_tree
      clump_grass               = nl%clump_grass
      igoutput                  = nl%igoutput
      gfilout                   = nl%gfilout
      h2o_plant_lim             = nl%h2o_plant_lim
      plant_hydro_scheme        = nl%plant_hydro_scheme
      istomata_scheme           = nl%istomata_scheme
      istruct_growth_scheme     = nl%istruct_growth_scheme
      istem_respiration_scheme  = nl%istem_respiration_scheme
      trait_plasticity_scheme   = nl%trait_plasticity_scheme
      iddmort_scheme            = nl%iddmort_scheme
      cbr_scheme                = nl%cbr_scheme
      ddmort_const              = nl%ddmort_const
      carbon_mortality_scheme   = nl%carbon_mortality_scheme
      hydraulic_mortality_scheme= nl%hydraulic_mortality_scheme
      vmfact_c3                 = nl%vmfact_c3
      vmfact_c4                 = nl%vmfact_c4
      mphoto_trc3               = nl%mphoto_trc3
      mphoto_tec3               = nl%mphoto_tec3
      mphoto_c4                 = nl%mphoto_c4
      bphoto_blc3               = nl%bphoto_blc3
      bphoto_nlc3               = nl%bphoto_nlc3
      bphoto_c4                 = nl%bphoto_c4
      kw_grass                  = nl%kw_grass
      kw_tree                   = nl%kw_tree
      gamma_c3                  = nl%gamma_c3
      gamma_c4                  = nl%gamma_c4
      d0_grass                  = nl%d0_grass
      d0_tree                   = nl%d0_tree
      alpha_c3                  = nl%alpha_c3
      alpha_c4                  = nl%alpha_c4
      klowco2in                 = nl%klowco2in
      rrffact                   = nl%rrffact
      growthresp                = nl%growthresp
      lwidth_grass              = nl%lwidth_grass
      lwidth_bltree             = nl%lwidth_bltree
      lwidth_nltree             = nl%lwidth_nltree
      q10_c3                    = nl%q10_c3
      q10_c4                    = nl%q10_c4
      thetacrit                 = nl%thetacrit
      quantum_efficiency_T      = nl%quantum_efficiency_T
      radint                    = nl%radint
      radslp                    = nl%radslp
      n_plant_lim               = nl%n_plant_lim
      n_decomp_lim              = nl%n_decomp_lim
      include_fire              = nl%include_fire
      fire_parameter            = nl%fire_parameter
      sm_fire                   = nl%sm_fire
      ianth_disturb             = nl%ianth_disturb
      sl_scale                  = nl%sl_scale
      sl_yr_first               = nl%sl_yr_first
      sl_nyrs                   = nl%sl_nyrs
      sl_pft                    = nl%sl_pft
      sl_prob_harvest           = nl%sl_prob_harvest
      sl_mindbh_harvest         = nl%sl_mindbh_harvest
      sl_biomass_harvest        = nl%sl_biomass_harvest
      sl_skid_rel_area          = nl%sl_skid_rel_area
      sl_skid_s_gtharv          = nl%sl_skid_s_gtharv
      sl_skid_s_ltharv          = nl%sl_skid_s_ltharv
      sl_felling_s_ltharv       = nl%sl_felling_s_ltharv
      cl_fseeds_harvest         = nl%cl_fseeds_harvest
      cl_fstorage_harvest       = nl%cl_fstorage_harvest
      cl_fleaf_harvest          = nl%cl_fleaf_harvest

      decomp_scheme             = nl%decomp_scheme

      icanturb                  = nl%icanturb
      isfclyrm                  = nl%isfclyrm
      ied_grndvap               = nl%ied_grndvap
      gamm                      = nl%gamm
      gamh                      = nl%gamh
      tprandtl                  = nl%tprandtl
      ribmax                    = nl%ribmax
      leaf_maxwhc               = nl%leaf_maxwhc
      ipercol                   = nl%ipercol

      include_these_pft         = nl%include_these_pft
      pasture_stock             = nl%pasture_stock
      agri_stock                = nl%agri_stock
      plantation_stock          = nl%plantation_stock
      pft_1st_check             = nl%pft_1st_check
      
      treefall_disturbance_rate = nl%treefall_disturbance_rate
      time2canopy               = nl%time2canopy
      runoff_time               = nl%runoff_time
      ubmin                     = nl%ubmin
      ugbmin                    = nl%ugbmin
      ustmin                    = nl%ustmin
      
      growth_resp_scheme        = nl%growth_resp_scheme
      storage_resp_scheme       = nl%storage_resp_scheme

      !----- Print control parameters. ----------------------------------------------------!
      iprintpolys               = nl%iprintpolys
      npvars                    = nl%npvars
      printvars                 = nl%printvars
      pfmtstr                   = nl%pfmtstr
      ipmin                     = nl%ipmin
      ipmax                     = nl%ipmax

      imettype                  = nl%imettype
      ishuffle                  = nl%ishuffle
      metcyc1                   = nl%metcyc1
      metcycf                   = nl%metcycf
      imetavg                   = nl%imetavg
      imetrad                   = nl%imetrad
      initial_co2               = nl%initial_co2
      
      iphenys1                  = nl%iphenys1
      iphenysf                  = nl%iphenysf
      iphenyf1                  = nl%iphenyf1
      iphenyff                  = nl%iphenyff
      
      iedcnfgf                  = nl%iedcnfgf
      event_file                = nl%event_file
      phenpath                  = nl%phenpath
      ifusion                   = nl%ifusion
      maxsite                   = nl%maxsite
      maxpatch                  = nl%maxpatch
      maxcohort                 = nl%maxcohort
      min_site_area             = nl%min_site_area
      min_patch_area            = nl%min_patch_area
      zrough                    = nl%zrough

      dt_census                 = nl%dt_census
      yr1st_census              = nl%yr1st_census
      mon1st_census             = nl%mon1st_census
      min_recruit_dbh           = nl%min_recruit_dbh
      idetailed                 = nl%idetailed
      patch_keep                = nl%patch_keep

      nnxp                      = nl%nnxp
      nnyp                      = nl%nnyp

      deltax                    = nl%deltax
      deltay                    = nl%deltay
      
      polelat                   = nl%polelat
      polelon                   = nl%polelon
      
      centlat                   = nl%centlat
      centlon                   = nl%centlon
      
      nstratx                   = nl%nstratx
      nstraty                   = nl%nstraty
      
      edres                     = nl%edres

      !------------------------------------------------------------------------------------!
      !     If the grid type is lat/lon, then we reset nnxp and nnyp to fit this new grid. !
      ! This is going to be useful to distribute the polygons across the nodes.            !
      !------------------------------------------------------------------------------------!
      ngrids = n_ed_region + n_poi
      do ifm=1,n_ed_region
         if (grid_type == 0) then
            nnxp(ifm)=floor( real(nstratx(ifm)) * (ed_reg_lonmax(ifm)-ed_reg_lonmin(ifm))  &
                           / grid_res)
            nnyp(ifm)=floor( real(nstratx(ifm)) * (ed_reg_latmax(ifm)-ed_reg_latmin(ifm))  &
                           / grid_res)
         end if
      end do

      do ifm=n_ed_region+1,ngrids
         nnxp(ifm) = 1
         nnyp(ifm) = 1
         nstratx(ifm)=1
         nstraty(ifm)=1
      end do

      !------------------------------------------------------------------------------------!
      !      Set current time to initial time here.  If this is a history run, reset       !
      ! current time in subroutine history_start.                                          !
      !------------------------------------------------------------------------------------!
      end_time%year  = iyearz
      end_time%month = imonthz
      end_time%date  = idatez
      end_time%time  = real(int(real(itimez) * 0.01)) * hr_sec                             &
                     + (real(itimez) * 0.01 - real(int(real(itimez)*0.01)))                &
                     * 100.0 * min_sec

   case ('NOT_HISTORY')
      !------------------------------------------------------------------------------------!
      !      The namelist variables in this section either must not be changed on a        !
      ! history restart or changing them would be irrelevant.  Thus, they are only copied  !
      ! to main model memory if this is not a history restart.                             !
      !------------------------------------------------------------------------------------!

      itimea        = nl%itimea
      idatea        = nl%idatea
      imontha       = nl%imontha
      iyeara        = nl%iyeara

      nzg           = nl%nzg
      nzs           = nl%nzs

      slz(1:nzgmax) = nl%slz(1:nzgmax)
      
      !------------------------------------------------------------------------------------!
      !      Set current time to initial time here.  If this is a history run, then reset  !
      ! current time in subroutine history_start.                                          !
      !------------------------------------------------------------------------------------!
      current_time%year  = iyeara
      current_time%month = imontha
      current_time%date  = idatea
      current_time%time  = real(int(real(itimea) * 0.01)) * hr_sec                         &
                         + (real(itimea) * 0.01 - real(int(real(itimea)*0.01)))            &
                         * 100.0 * min_sec
      time               = 0.0d0

   case ('HISTORY')
      !------------------------------------------------------------------------------------!
      !      The namelist variables in this section either must not be changed on a        !
      ! history restart or changing them would be irrelevant.  Thus, they are only copied  !
      ! to main model memory if this is not a history restart.                             !
      !------------------------------------------------------------------------------------!
      itimea        = nl%itimea
      idatea        = nl%idatea
      imontha       = nl%imontha
      iyeara        = nl%iyeara

      nzg           = nl%nzg
      nzs           = nl%nzs

      slz(1:nzgmax) = nl%slz(1:nzgmax)
      
      !------------------------------------------------------------------------------------!
      !      Set current time to initial time here.  If this is a history run, reset       !
      ! current time in subroutine history_start.                             !
      !------------------------------------------------------------------------------------!
      current_time%year  = nl%iyearh
      current_time%month = nl%imonthh
      current_time%date  = nl%idateh
      current_time%time  = real(int(real(nl%itimeh) * 0.01)) * hr_sec                      &
                         + (real(nl%itimeh) * 0.01 - real(int(real(nl%itimeh)*0.01)))      &
                         * 100.0 * min_sec
      !----- Calculate the current time. --------------------------------------------------!
      call date_2_seconds (nl%iyearh,nl%imonthh,nl%idateh,nl%itimeh*100,iyeara,imontha     &
                          ,idatea,itimea*100,time)

   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The following variable will be used to allocate the mean diurnal cycle.  It will  !
   ! be set to 1 in case the user doesn't want the mean diurnal cycle, or if frqanl is     !
   ! invalid.                                                                              !
   !---------------------------------------------------------------------------------------!
   if (iqoutput == 0 .or. frqfast <= 0 .or. unitfast > 0) then
      ndcycle = 1 
   else
      ndcycle = max(1,int(day_sec / frqfast))
   end if
   !---------------------------------------------------------------------------------------!




   !----- Sort up the chosen PFTs. --------------------------------------------------------!
   where (include_these_pft < 1 .or. include_these_pft > n_pft) 
      include_these_pft = skip_integer
   end where
   call sort_up(include_these_pft,n_pft)
   !---------------------------------------------------------------------------------------!




   !----- Sort up the PFTs that may be logged. --------------------------------------------!
   where (sl_pft < 1 .or. sl_pft > n_pft)
      sl_pft            = skip_integer
      sl_mindbh_harvest = skip_real
      sl_prob_harvest   = skip_real
   end where
   call rank_up_i(n_pft,sl_pft,idx)
   sl_pft           (idx) = sl_pft(:)
   sl_mindbh_harvest(idx) = sl_mindbh_harvest(:)
   sl_prob_harvest  (idx) = sl_prob_harvest(:)
   !---------------------------------------------------------------------------------------!

      
   !----- Determine the length of simuation. ----------------------------------------------!
   call date_2_seconds (iyearz,imonthz,idatez,itimez*100,iyeara,imontha,idatea,itimea*100  &
                       ,timmax)

   !---------------------------------------------------------------------------------------!
   !     For the following databases, we must check whether only the first grid was given. !
   ! In this case, we copy the values of the first grid to the other grids.                !
   !---------------------------------------------------------------------------------------!
   call copy_path_from_grid_1(ngrids,'soil_database'  ,soil_database  )
   call copy_path_from_grid_1(ngrids,'veg_database'   ,veg_database   )
   call copy_path_from_grid_1(ngrids,'lu_database'    ,lu_database    )
   call copy_path_from_grid_1(ngrids,'plantation_file',plantation_file)
   call copy_path_from_grid_1(ngrids,'lu_rescale_file',lu_rescale_file)

   call copy_path_from_grid_1(ngrids,'sfilin'         ,sfilin         )
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !       Check whether we must re-define the default soil type (nslcon) based on the     !
   ! fractions of sand and clay.                                                           !
   !---------------------------------------------------------------------------------------!
   if ( any(isoilflg == 2) .and. slxclay > 0. .and. slxsand > 0. .and.                     &
        (slxclay + slxsand) <= 1. ) then
      nslcon = find_soil_class(slxsand,slxclay)
   end if
   !---------------------------------------------------------------------------------------!



   !----- Define some useful variables that control the output. ---------------------------!
   writing_dail     = idoutput > 0
   writing_mont     = imoutput > 0
   writing_dcyc     = iqoutput > 0
   writing_year     = iyoutput > 0
   writing_long     = writing_dail .or. writing_mont .or. writing_dcyc
   writing_eorq     = writing_mont .or. writing_dcyc
   history_fast     = ifoutput == 0  .and.                                                 &
                      ( ( unitstate == 0 .and. frqstate <= day_sec ) .or.                  &
                        ( unitstate == 1 .and. frqstate == 1.      )      )
   history_dail     = unitstate == 0 .and. frqstate < day_sec
   history_eorq     = unitstate <= 1
   fast_diagnostics = ifoutput /= 0 .or. idoutput /= 0 .or. iooutput /= 0 .or.             &
                      imoutput /= 0 .or. iqoutput /= 0 .or. itoutput /= 0
   !---------------------------------------------------------------------------------------!


   return
end subroutine copy_nl
!==========================================================================================!
!==========================================================================================!
