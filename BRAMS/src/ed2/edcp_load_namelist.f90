!==========================================================================================!
!==========================================================================================!
!     This subroutine will read the ED parameters from the ED2_INFO namelist (part of the  !
! RAMSIN file).                                                                            !
!------------------------------------------------------------------------------------------!
subroutine read_ednl(iunit,filename)
   !----- ED2 modules. --------------------------------------------------------------------!
   use ed_max_dims          , only : n_pft                                 & ! intent(in)
                                   , undef_integer                         & ! intent(in)
                                   , undef_path                            & ! intent(in)
                                   , maxgrds                               ! ! intent(in)
   use soil_coms            , only : ed_zrough => zrough                   & ! intent(out)
                                   , soil_database                         & ! intent(out)
                                   , isoilstateinit                        & ! intent(out)
                                   , isoildepthflg                         & ! intent(out)
                                   , isoilbc                               & ! intent(out)
                                   , soilstate_db                          & ! intent(out)
                                   , soildepth_db                          & ! intent(out)
                                   , runoff_time                           & ! intent(out)
                                   , veg_database                          & ! intent(out)
                                   , slxclay                               & ! intent(out)
                                   , slxsand                               ! ! intent(out)
   use met_driver_coms      , only : ed_met_driver_db                      & ! intent(out)
                                   , imettype                              & ! intent(out)
                                   , metcyc1                               & ! intent(out)
                                   , metcycf                               & ! intent(out)
                                   , imetavg                               & ! intent(out)
                                   , imetrad                               & ! intent(out)
                                   , lapse_scheme                          ! ! intent(out)
   use mem_polygons         , only : n_poi                                 & ! intent(out)
                                   , poi_lat                               & ! intent(out)
                                   , poi_lon                               & ! intent(out)
                                   , poi_res                               & ! intent(out)
                                   , n_ed_region                           & ! intent(out)
                                   , ed_reg_latmin                         & ! intent(out)
                                   , ed_reg_latmax                         & ! intent(out)
                                   , ed_reg_lonmin                         & ! intent(out)
                                   , ed_reg_lonmax                         & ! intent(out)
                                   , grid_res                              & ! intent(out)
                                   , grid_type                             & ! intent(out)
                                   , edres                                 & ! intent(out)
                                   , maxsite                               & ! intent(out)
                                   , maxpatch                              & ! intent(out)
                                   , maxcohort                             ! ! intent(out)
   use physiology_coms      , only : iphysiol                              & ! intent(out)
                                   , quantum_efficiency_t                  & ! intent(out)
                                   , h2o_plant_lim                         & ! intent(out)
                                   , n_plant_lim                           & ! intent(out)
                                   , vmfact_c3                             & ! intent(out)
                                   , vmfact_c4                             & ! intent(out)
                                   , mphoto_trc3                           & ! intent(out)
                                   , mphoto_tec3                           & ! intent(out)
                                   , mphoto_c4                             & ! intent(out)
                                   , bphoto_blc3                           & ! intent(out)
                                   , bphoto_nlc3                           & ! intent(out)
                                   , bphoto_c4                             & ! intent(out)
                                   , kw_grass                              & ! intent(out)
                                   , kw_tree                               & ! intent(out)
                                   , gamma_c3                              & ! intent(out)
                                   , gamma_c4                              & ! intent(out)
                                   , d0_grass                              & ! intent(out)
                                   , d0_tree                               & ! intent(out)
                                   , alpha_c3                              & ! intent(out)
                                   , alpha_c4                              & ! intent(out)
                                   , klowco2in                             & ! intent(out)
                                   , rrffact                               & ! intent(out)
                                   , growthresp                            & ! intent(out)
                                   , lwidth_grass                          & ! intent(out)
                                   , lwidth_bltree                         & ! intent(out)
                                   , lwidth_nltree                         & ! intent(out)
                                   , q10_c3                                & ! intent(out)
                                   , q10_c4                                ! ! intent(out)
   use phenology_coms       , only : iphen_scheme                          & ! intent(out)
                                   , repro_scheme                          & ! intent(out)
                                   , iphenys1                              & ! intent(out)
                                   , iphenysf                              & ! intent(out)
                                   , iphenyf1                              & ! intent(out)
                                   , iphenyff                              & ! intent(out)
                                   , phenpath                              & ! intent(out)
                                   , radint                                & ! intent(out)
                                   , radslp                                & ! intent(out)
                                   , thetacrit                             ! ! intent(out)
   use decomp_coms          , only : n_decomp_lim                          & ! intent(out)
                                   , LloydTaylor                           ! ! intent(out)
   use disturb_coms         , only : include_fire                          & ! intent(out)
                                   , ianth_disturb                         & ! intent(out)
                                   , lu_database                           & ! intent(out)
                                   , plantation_file                       & ! intent(out)
                                   , lu_rescale_file                       & ! intent(out)
                                   , treefall_disturbance_rate             & ! intent(out)
                                   , time2canopy                           & ! intent(out)
                                   , sm_fire                               & ! intent(out)
                                   , min_patch_area                        ! ! intent(out)
   use pft_coms             , only : include_these_pft                     & ! intent(out)
                                   , agri_stock                            & ! intent(out)
                                   , plantation_stock                      & ! intent(out)
                                   , pft_1st_check                         ! ! intent(out)
   use ed_misc_coms         , only : ifoutput                              & ! intent(out)
                                   , idoutput                              & ! intent(out)
                                   , imoutput                              & ! intent(out)
                                   , iqoutput                              & ! intent(out)
                                   , iyoutput                              & ! intent(out)
                                   , itoutput                              & ! intent(out)
                                   , isoutput                              & ! intent(out)
                                   , frqfast                               & ! intent(out)
                                   , frqstate                              & ! intent(out)
                                   , outfast                               & ! intent(out)
                                   , outstate                              & ! intent(out)
                                   , unitfast                              & ! intent(out)
                                   , unitstate                             & ! intent(out)
                                   , ndcycle                               & ! intent(out)
                                   , ied_init_mode                         & ! intent(out)
                                   , current_time                          & ! intent(out)
                                   , thsums_database                       & ! intent(out)
                                   , end_time                              & ! intent(out)
                                   , ivegt_dynamics                        & ! intent(out)
                                   , integration_scheme                    & ! intent(out)
                                   , ffilout                               & ! intent(out)
                                   , dtlsm                                 & ! intent(out)
                                   , iprintpolys                           & ! intent(out)
                                   , printvars                             & ! intent(out)
                                   , npvars                                & ! intent(out)
                                   , pfmtstr                               & ! intent(out)
                                   , ipmax                                 & ! intent(out)
                                   , ipmin                                 & ! intent(out)
                                   , iedcnfgf                              & ! intent(out)
                                   , ffilout                               & ! intent(out)
                                   , sfilout                               & ! intent(out)
                                   , sfilin                                & ! intent(out)
                                   , event_file                            & ! intent(out)
                                   , attach_metadata                       & ! intent(out)
                                   , iallom                                & ! intent(out)
                                   , min_site_area                         ! ! intent(out)
   use canopy_air_coms      , only : icanturb                              & ! intent(out)
                                   , isfclyrm                              & ! intent(out)
                                   , ied_grndvap                           ! ! intent(out)
   use grid_coms            , only : timmax                                & ! intent(out)
                                   , time                                  ! ! intent(out)
   use optimiz_coms         , only : ioptinpt                              ! ! intent(out)
   use rk4_coms             , only : rk4_tolerance                         & ! intent(out)
                                   , ibranch_thermo                        ! ! intent(out)
   use canopy_layer_coms    , only : crown_mod                             ! ! intent(out)
   use canopy_radiation_coms, only : icanrad                               & ! intent(out)
                                   , ltrans_vis                            & ! intent(out)
                                   , ltrans_nir                            & ! intent(out)
                                   , lreflect_vis                          & ! intent(out)
                                   , lreflect_nir                          & ! intent(out)
                                   , orient_tree                           & ! intent(out)
                                   , orient_grass                          & ! intent(out)
                                   , clump_tree                            & ! intent(out)
                                   , clump_grass                           ! ! intent(out)
   !----- Coupled ED-BRAMS modules. -------------------------------------------------------!
   use mem_edcp             , only : co2_offset                            ! ! intent(out)
   !----- BRAMS modules. ------------------------------------------------------------------!
   use mem_grid             , only : expnme                                & ! intent(in)
                                   , runtype                               & ! intent(in)
                                   , itimez                                & ! intent(in)
                                   , idatez                                & ! intent(in)
                                   , imonthz                               & ! intent(in)
                                   , iyearz                                & ! intent(in)
                                   , itimea                                & ! intent(in)
                                   , idatea                                & ! intent(in)
                                   , imontha                               & ! intent(in)
                                   , iyeara                                & ! intent(in)
                                   , centlon                               & ! intent(in)
                                   , centlat                               & ! intent(in)
                                   , deltax                                & ! intent(in)
                                   , deltay                                & ! intent(in)
                                   , nnxp                                  & ! intent(in)
                                   , nnyp                                  & ! intent(in)
                                   , nstratx                               & ! intent(in)
                                   , nstraty                               & ! intent(in)
                                   , polelat                               & ! intent(in)
                                   , polelon                               & ! intent(in)
                                   , ngrids                                & ! intent(in)
                                   , nzg                                   & ! intent(in)
                                   , nzs                                   & ! intent(in)
                                   , npatch                                ! ! intent(in)
   use io_params            , only : ioutput                               & ! intent(in)
                                   , iclobber                              & ! intent(in)
                                   , frqanl                                & ! intent(in)
                                   , frqhis                                & ! intent(in)
                                   , iyearh                                & ! intent(in)
                                   , idateh                                & ! intent(in)
                                   , itimeh                                & ! intent(in)
                                   , imonthh                               & ! intent(in)
                                   , isoilflg                              ! ! intent(in)
   use mem_leaf             , only : nslcon                                & ! intent(in)
                                   , isoilcol                              & ! intent(in)
                                   , slz                                   & ! intent(in)
                                   , stgoff                                & ! intent(in)
                                   , slmstr                                & ! intent(in)
                                   , isfcl                                 & ! intent(in)
                                   , nvegpat                               & ! intent(in)
                                   , istar                                 & ! intent(in)
                                   , igrndvap                              & ! intent(in)
                                   , leaf_zrough         => zrough         & ! intent(in)
                                   , leaf_isoilbc        => isoilbc        & ! intent(in)
                                   , leaf_ipercol        => ipercol        & ! intent(in)
                                   , leaf_runoff_time    => runoff_time    ! ! intent(in)
   use leaf_coms            , only : leaf_ubmin          => ubmin          & ! intent(in)
                                   , leaf_ugbmin         => ugbmin         & ! intent(in)
                                   , leaf_ustmin         => ustmin         & ! intent(in)
                                   , leaf_gamm           => gamm           & ! intent(in)
                                   , leaf_gamh           => gamh           & ! intent(in)
                                   , leaf_tprandtl       => tprandtl       & ! intent(in)
                                   , leaf_ribmax         => ribmax         & ! intent(in)
                                   , leaf_leaf_maxwhc    => leaf_maxwhc    & ! intent(in)
                                   , leaf_min_patch_area => min_patch_area ! ! intent(in)
   use mem_radiate          , only : radfrq                                ! ! intent(in)
   use consts_coms          , only : day_sec                               ! ! intent(in)
   use detailed_coms        , only : idetailed                             & ! intent(in)
                                   , patch_keep                            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer         , intent(in) :: iunit    ! Namelist unit number
   character(len=*), intent(in) :: filename ! file name with namelists
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: i
   integer                      :: err
   integer                      :: decomp_scheme
   logical                      :: fexists
   logical                      :: op
   !----- Namelist. -----------------------------------------------------------------------!
   namelist /ED2_INFO/  dtlsm,co2_offset,ifoutput,idoutput,imoutput,iqoutput,iyoutput      &
                       ,itoutput,isoutput,attach_metadata,outfast,outstate,ffilout,sfilout &
                       ,ied_init_mode,edres,sfilin,veg_database,soil_database,lu_database  &
                       ,plantation_file,lu_rescale_file,thsums_database,soilstate_db       &
                       ,soildepth_db,isoilstateinit,isoildepthflg,ivegt_dynamics           &
                       ,integration_scheme,rk4_tolerance,ibranch_thermo,iphysiol,iallom    &
                       ,iphen_scheme,radint,radslp,repro_scheme,lapse_scheme,crown_mod     &
                       ,icanrad,ltrans_vis,ltrans_nir,lreflect_vis,lreflect_nir            &
                       ,orient_tree,orient_grass,clump_tree,clump_grass,decomp_scheme      &
                       ,h2o_plant_lim,vmfact_c3,vmfact_c4,mphoto_trc3,mphoto_tec3          &
                       ,mphoto_c4,bphoto_blc3,bphoto_nlc3,bphoto_c4,kw_grass               &
                       ,kw_tree,gamma_c3,gamma_c4,d0_grass,d0_tree,alpha_c3                &
                       ,alpha_c4,klowco2in,rrffact,growthresp,lwidth_grass,lwidth_bltree   &
                       ,lwidth_nltree,q10_c3,q10_c4,thetacrit,quantum_efficiency_t         &
                       ,n_plant_lim,n_decomp_lim,include_fire,sm_fire,ianth_disturb        &
                       ,icanturb,include_these_pft,agri_stock,plantation_stock             &
                       ,pft_1st_check,maxpatch,maxcohort,min_patch_area                    &
                       ,treefall_disturbance_rate,time2canopy,iprintpolys,npvars,printvars &
                       ,pfmtstr,ipmin,ipmax,imetrad,iphenys1,iphenysf,iphenyf1,iphenyff    &
                       ,iedcnfgf,event_file,phenpath

   !----- Initialise some database variables with a non-sense path. -----------------------!
   soil_database   (:) = undef_path
   veg_database    (:) = undef_path
   lu_database     (:) = undef_path
   plantation_file (:) = undef_path
   lu_rescale_file (:) = undef_path

   sfilin          (:) = undef_path
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Read in MODEL_FILE_INFO namelist.                                                !
   !---------------------------------------------------------------------------------------!
   read (unit=iunit, iostat=err, NML=ED2_INFO)
   if (err /= 0) then
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ' ERROR reading the namelist! '
      write (unit=*,fmt='(a,1x,a)')   ' Namelist file   : ',trim(filename)
      write (unit=*,fmt='(a,1x,a)')   ' Namelist section: ','MODEL_FILE_INFO'
      write (unit=*,fmt='(a)')        ' Compare values read with file contents...'
      write (unit=*,fmt='(a)')        '--------------------------------------------------'
      write (unit=*,fmt='(a)')        ''
      write (unit=*,fmt=*) ' dtlsm                     =',dtlsm
      write (unit=*,fmt=*) ' co2_offset                =',co2_offset
      write (unit=*,fmt=*) ' ifoutput                  =',ifoutput
      write (unit=*,fmt=*) ' idoutput                  =',idoutput
      write (unit=*,fmt=*) ' imoutput                  =',imoutput
      write (unit=*,fmt=*) ' iqoutput                  =',iqoutput
      write (unit=*,fmt=*) ' iyoutput                  =',iyoutput
      write (unit=*,fmt=*) ' itoutput                  =',itoutput
      write (unit=*,fmt=*) ' isoutput                  =',isoutput
      write (unit=*,fmt=*) ' attach_metadata           =',attach_metadata
      write (unit=*,fmt=*) ' outfast                   =',outfast
      write (unit=*,fmt=*) ' outstate                  =',outstate
      write (unit=*,fmt=*) ' ffilout                   =',trim(ffilout)
      write (unit=*,fmt=*) ' sfilout                   =',trim(sfilout)
      write (unit=*,fmt=*) ' ied_init_mode             =',ied_init_mode
      write (unit=*,fmt=*) ' edres                     =',edres
      write (unit=*,fmt=*) ' sfilin                    =',(trim(sfilin(i))//';'            &
                                                          ,i=1,size(sfilin))
      write (unit=*,fmt=*) ' veg_database              =',(trim(veg_database(i))//';'      &
                                                          ,i=1,size(veg_database))
      write (unit=*,fmt=*) ' soil_database             =',(trim(soil_database(i))//';'     &
                                                          ,i=1,size(soil_database))
      write (unit=*,fmt=*) ' lu_database               =',(trim(lu_database(i))//';'       &
                                                          ,i=1,size(lu_database))
      write (unit=*,fmt=*) ' plantation_file           =',(trim(plantation_file(i))//';'   &
                                                          ,i=1,size(plantation_file))
      write (unit=*,fmt=*) ' lu_rescale_file           =',(trim(lu_rescale_file(i))//';'   &
                                                          ,i=1,size(lu_rescale_file))
      write (unit=*,fmt=*) ' thsums_database           =',trim(thsums_database)
      write (unit=*,fmt=*) ' soilstate_db              =',trim(soilstate_db)
      write (unit=*,fmt=*) ' soildepth_db              =',trim(soildepth_db)
      write (unit=*,fmt=*) ' isoilstateinit            =',isoilstateinit
      write (unit=*,fmt=*) ' isoildepthflg             =',isoildepthflg
      write (unit=*,fmt=*) ' ivegt_dynamics            =',ivegt_dynamics
      write (unit=*,fmt=*) ' integration_scheme        =',integration_scheme
      write (unit=*,fmt=*) ' rk4_tolerance             =',rk4_tolerance
      write (unit=*,fmt=*) ' ibranch_thermo            =',ibranch_thermo
      write (unit=*,fmt=*) ' iphysiol                  =',iphysiol
      write (unit=*,fmt=*) ' iallom                    =',iallom
      write (unit=*,fmt=*) ' iphen_scheme              =',iphen_scheme
      write (unit=*,fmt=*) ' radint                    =',radint
      write (unit=*,fmt=*) ' radslp                    =',radslp
      write (unit=*,fmt=*) ' repro_scheme              =',repro_scheme
      write (unit=*,fmt=*) ' lapse_scheme              =',lapse_scheme
      write (unit=*,fmt=*) ' crown_mod                 =',crown_mod
      write (unit=*,fmt=*) ' icanrad                   =',icanrad
      write (unit=*,fmt=*) ' ltrans_vis                =',ltrans_vis
      write (unit=*,fmt=*) ' ltrans_nir                =',ltrans_nir
      write (unit=*,fmt=*) ' lreflect_vis              =',lreflect_vis
      write (unit=*,fmt=*) ' lreflect_nir              =',lreflect_nir
      write (unit=*,fmt=*) ' orient_tree               =',orient_tree
      write (unit=*,fmt=*) ' orient_grass              =',orient_grass
      write (unit=*,fmt=*) ' clump_tree                =',clump_tree
      write (unit=*,fmt=*) ' clump_grass               =',clump_grass
      write (unit=*,fmt=*) ' decomp_scheme             =',decomp_scheme
      write (unit=*,fmt=*) ' h2o_plant_lim             =',h2o_plant_lim
      write (unit=*,fmt=*) ' vmfact_c3                 =',vmfact_c3
      write (unit=*,fmt=*) ' vmfact_c4                 =',vmfact_c4
      write (unit=*,fmt=*) ' mphoto_trc3               =',mphoto_trc3
      write (unit=*,fmt=*) ' mphoto_tec3               =',mphoto_tec3
      write (unit=*,fmt=*) ' mphoto_c4                 =',mphoto_c4
      write (unit=*,fmt=*) ' bphoto_blc3               =',bphoto_blc3
      write (unit=*,fmt=*) ' bphoto_nlc3               =',bphoto_nlc3
      write (unit=*,fmt=*) ' bphoto_c4                 =',bphoto_c4
      write (unit=*,fmt=*) ' kw_grass                  =',kw_grass
      write (unit=*,fmt=*) ' kw_tree                   =',kw_tree
      write (unit=*,fmt=*) ' gamma_c3                  =',gamma_c3
      write (unit=*,fmt=*) ' gamma_c4                  =',gamma_c4
      write (unit=*,fmt=*) ' d0_grass                  =',d0_grass
      write (unit=*,fmt=*) ' d0_tree                   =',d0_tree
      write (unit=*,fmt=*) ' alpha_c3                  =',alpha_c3
      write (unit=*,fmt=*) ' alpha_c4                  =',alpha_c4
      write (unit=*,fmt=*) ' klowco2in                 =',klowco2in
      write (unit=*,fmt=*) ' rrffact                   =',rrffact
      write (unit=*,fmt=*) ' growthresp                =',growthresp
      write (unit=*,fmt=*) ' lwidth_grass              =',lwidth_grass
      write (unit=*,fmt=*) ' lwidth_bltree             =',lwidth_bltree
      write (unit=*,fmt=*) ' lwidth_nltree             =',lwidth_nltree
      write (unit=*,fmt=*) ' q10_c3                    =',q10_c3
      write (unit=*,fmt=*) ' q10_c4                    =',q10_c4
      write (unit=*,fmt=*) ' thetacrit                 =',thetacrit
      write (unit=*,fmt=*) ' quantum_efficiency_t      =',quantum_efficiency_t
      write (unit=*,fmt=*) ' n_plant_lim               =',n_plant_lim
      write (unit=*,fmt=*) ' n_decomp_lim              =',n_decomp_lim
      write (unit=*,fmt=*) ' include_fire              =',include_fire
      write (unit=*,fmt=*) ' sm_fire                   =',sm_fire
      write (unit=*,fmt=*) ' ianth_disturb             =',ianth_disturb
      write (unit=*,fmt=*) ' icanturb                  =',icanturb
      write (unit=*,fmt=*) ' include_these_pft         =',include_these_pft
      write (unit=*,fmt=*) ' agri_stock                =',agri_stock
      write (unit=*,fmt=*) ' plantation_stock          =',plantation_stock
      write (unit=*,fmt=*) ' pft_1st_check             =',pft_1st_check
      write (unit=*,fmt=*) ' maxsite                   =',maxsite
      write (unit=*,fmt=*) ' maxpatch                  =',maxpatch
      write (unit=*,fmt=*) ' maxcohort                 =',maxcohort
      write (unit=*,fmt=*) ' min_patch_area            =',min_patch_area
      write (unit=*,fmt=*) ' treefall_disturbance_rate =',treefall_disturbance_rate
      write (unit=*,fmt=*) ' time2canopy               =',time2canopy
      write (unit=*,fmt=*) ' iprintpolys               =',iprintpolys
      write (unit=*,fmt=*) ' npvars                    =',npvars
      write (unit=*,fmt=*) ' printvars                 =',(trim(printvars(i))//';'         &
                                                          ,i=1,size(printvars))
      write (unit=*,fmt=*) ' pfmtstr                   =',(trim(pfmtstr(i))//';'           &
                                                          ,i=1,size(pfmtstr))
      write (unit=*,fmt=*) ' ipmin                     =',ipmin
      write (unit=*,fmt=*) ' ipmax                     =',ipmax
      write (unit=*,fmt=*) ' imetrad                   =',imetrad
      write (unit=*,fmt=*) ' iphenys1                  =',iphenys1
      write (unit=*,fmt=*) ' iphenysf                  =',iphenysf
      write (unit=*,fmt=*) ' iphenyf1                  =',iphenyf1
      write (unit=*,fmt=*) ' iphenyff                  =',iphenyff
      write (unit=*,fmt=*) ' iedcnfgf                  =',trim(iedcnfgf)
      write (unit=*,fmt=*) ' event_file                =',trim(event_file)
      write (unit=*,fmt=*) ' phenpath                  =',trim(phenpath)
      call abort_run('Error reading namelist, ED2_INFO block.','read_ednl'                 &
                    ,'edcp_load_namelist.f90')
   end if
   
   !---------------------------------------------------------------------------------------!
   !    Decomposition scheme is not a real variable in the model, internally we use        !
   ! Lloyd_Taylor instead.                                                                 !
   !---------------------------------------------------------------------------------------!
   LloydTaylor = decomp_scheme == 1

   !---------------------------------------------------------------------------------------!
   !     Some variables that ED needs are also defined and used by other BRAMS modules.    !
   ! To avoid defining these variabes twice, or even worse, to avoid mismatches, we copy   !
   ! the variables from BRAMS modules to the ED modules.                                   !
   !---------------------------------------------------------------------------------------!
   call copy_in_bramsnl(expnme,runtype,itimez,idatez,imonthz,iyearz,itimea,idatea,imontha  &
                       ,iyeara,itimeh,idateh,imonthh,iyearh,radfrq,nnxp,nnyp,deltax        &
                       ,deltay,polelat,polelon,centlat,centlon,nstratx,nstraty,iclobber    &
                       ,nzg,nzs,isoilflg,nslcon,isoilcol,slz,slmstr,stgoff,leaf_zrough     &
                       ,ngrids,leaf_ubmin,leaf_ugbmin,leaf_ustmin,leaf_isoilbc             &
                       ,leaf_ipercol,leaf_runoff_time,leaf_gamm,leaf_gamh,leaf_tprandtl    &
                       ,leaf_ribmax,leaf_leaf_maxwhc)

   !---------------------------------------------------------------------------------------!
   !      The following variables can be defined in the regular ED2IN file for stand-alone !
   ! runs, but they cannot be changed in the coupled simulation (or they are never used    !
   ! in the coupled run).  We assign some standard values to these variables.              !
   !---------------------------------------------------------------------------------------!
   n_ed_region   = ngrids   ! No POI is allowed in coupled runs.
   n_poi = 0                ! No POI is allowed in coupled runs.
   grid_res      = 0        ! Not used: this is for lat-lon style.
   grid_type     = 1        ! We reinforce polar-stereo for now, it actually grabs the
                            !   BRAMS coordinates.
   ed_reg_latmin = 0        ! Not used in coupled runs.
   ed_reg_latmax = 0        ! Not used in coupled runs.
   ed_reg_lonmin = 0        ! Not used in coupled runs.
   ed_reg_lonmax = 0        ! Not used in coupled runs.
   poi_lat = 0.             ! POI is not an option in coupled runs.
   poi_lon = 0.             ! POI is not an option in coupled runs.
   poi_res = 0.             ! POI is not an option in coupled runs.
   ed_met_driver_db = ''    ! BRAMS is the meteorology driver... 
   imettype = 1             ! BRAMS is the meteorology driver...
   metcyc1  = 0000          ! BRAMS is the meteorology driver...
   metcycf  = 0000          ! BRAMS is the meteorology driver...
   imetavg  = 0             ! BRAMS is the meteorology driver...
   ioptinpt = ''            ! It will be used once optimization is 
                            !    implemented in ED-2.1.
   unitfast  = 0            ! Since BRAMS uses frqanl and frqhist in seconds, there is no
   unitstate = 0            !     reason to ask the user for units for outfast and 
                            !     outstate, the special flags cover all possibilities.
   slxclay   = -1.          ! This is not going to be used in coupled runs because the 
   slxsand   = -1.          !     soil should come from lon/lat maps.
   idetailed =  0           ! No detailed output in coupled runs (it is already too slow 
                            !     with the normal output...)
   patch_keep = 0           ! Keep all patches.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      We make sure that the maximum number of sites per polygon in ED2 is equivalent   !
   ! to the number of land patches (NPATCH-1) as defined in RAMSIN.  We fill some LEAF-3   !
   ! variables that are used by other BRAMS routines with ED values.  As in LEAF-3, the    !
   ! first leaf "patch" is reserved for water.  Also, we copy the minimum patch area to ED !
   ! minimum site area, so both models are consistent.                                     !
   !---------------------------------------------------------------------------------------!
   nvegpat       = npatch - 1
   maxsite       = nvegpat
   min_site_area = leaf_min_patch_area
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      These are ED2 variables that have the same function as certain BRAMS namelist    !
   ! variables under a different name.                                                     !
   !---------------------------------------------------------------------------------------!
   frqfast     = frqanl
   frqstate    = frqhis
   isfclyrm    = istar
   ied_grndvap = igrndvap
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     The following variable will be used to allocate the mean diurnal cycle.  It will  !
   ! be set to 1 in case the user doesn't want the mean diurnal cycle, or if frqanl is     !
   ! invalid.                                                                              !
   !---------------------------------------------------------------------------------------!
   if (iqoutput == 0 .or. frqfast <= 0.) then
      ndcycle = 1 
   else
      ndcycle = max(1,int(day_sec / frqfast))
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Set current time to initial time here.  If this is a history run, reset current  !
   ! time in subroutine history_start.                                                     !
   !---------------------------------------------------------------------------------------!
   end_time%year  = iyearz
   end_time%month = imonthz
   end_time%date  = idatez
   end_time%time  = int(itimez * 0.01) * 3600.0                                            &
                  + (itimez * 0.01 - int(itimez*0.01))*100.0*60.0

   !---------------------------------------------------------------------------------------!
   !     Here we set up the current time.  The time depends on whether this is a history   !
   ! or initial run.                                                                       !
   !---------------------------------------------------------------------------------------!
   select case(trim(runtype))
   case('INITIAL')
      current_time%year  = iyeara
      current_time%month = imontha
      current_time%date  = idatea
      current_time%time  = int(itimea * 0.01) * 3600.0                                     &
                         + (itimea * 0.01 - int(itimea*0.01))*100.0*60.0
      time               = 0.0

   case ('HISTORY')
      current_time%year  = iyearh
      current_time%month = imonthh
      current_time%date  = idateh
      current_time%time  = int(itimeh * 0.01) * 3600.0                                     &
                         + (itimeh * 0.01 - int(itimeh*0.01))*100.0*60.0

      !----- Calculate the current time. --------------------------------------------------!
      call date_2_seconds (iyearh,imonthh,idateh,itimeh*100,iyeara,imontha,idatea          &
                          ,itimea*100,time)

   end select
      
   !----- Sort up the chosen PFTs. --------------------------------------------------------!
   where (include_these_pft < 1 .or. include_these_pft == undef_integer) 
      include_these_pft = huge(1)
   end where
   call sort_up(include_these_pft,n_pft)

   !----- Determine the length of simuation. ----------------------------------------------!
   call date_2_seconds(iyearz,imonthz,idatez,itimez*100,iyeara,imontha,idatea,itimea*100   &
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

   return
end subroutine read_ednl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies variables that have the same name in both BRAMS and ED2       !
! modules, from BRAMS to ED.                                                               !
!------------------------------------------------------------------------------------------!
subroutine copy_in_bramsnl(expnme_b,runtype_b,itimez_b,idatez_b,imonthz_b,iyearz_b         &
                          ,itimea_b,idatea_b,imontha_b,iyeara_b,itimeh_b,idateh_b          &
                          ,imonthh_b,iyearh_b,radfrq_b,nnxp_b,nnyp_b,deltax_b,deltay_b     &
                          ,polelat_b,polelon_b,centlat_b,centlon_b,nstratx_b,nstraty_b     &
                          ,iclobber_b,nzg_b,nzs_b,isoilflg_b,nslcon_b,isoilcol_b,slz_b     &
                          ,slmstr_b,stgoff_b,zrough_b,ngrids_b,ubmin_b,ugbmin_b,ustmin_b   &
                          ,isoilbc_b,ipercol_b,runoff_time_b,gamm_b,gamh_b,tprandtl_b      &
                          ,ribmax_b,leaf_maxwhc_b)
   use ed_misc_coms   , only : expnme            & ! intent(out)
                             , runtype           & ! intent(out)
                             , itimez            & ! intent(out)
                             , idatez            & ! intent(out)
                             , imonthz           & ! intent(out)
                             , iyearz            & ! intent(out)
                             , itimea            & ! intent(out)
                             , idatea            & ! intent(out)
                             , imontha           & ! intent(out)
                             , iyeara            & ! intent(out)
                             , itimeh            & ! intent(out)
                             , idateh            & ! intent(out)
                             , imonthh           & ! intent(out)
                             , iyearh            & ! intent(out)
                             , iclobber          & ! intent(out)
                             , radfrq            ! ! intent(out)
   use grid_coms      , only : centlon           & ! intent(out)
                             , centlat           & ! intent(out)
                             , deltax            & ! intent(out)
                             , deltay            & ! intent(out)
                             , nnxp              & ! intent(out)
                             , nnyp              & ! intent(out)
                             , nstratx           & ! intent(out)
                             , nstraty           & ! intent(out)
                             , polelat           & ! intent(out)
                             , polelon           & ! intent(out)
                             , ngrids            & ! intent(out)
                             , timmax            & ! intent(out)
                             , time              & ! intent(out)
                             , nzg               & ! intent(out)
                             , nzs               ! ! intent(out)
   use soil_coms      , only : isoilflg          & ! intent(out)
                             , nslcon            & ! intent(out)
                             , isoilcol          & ! intent(out)
                             , slmstr            & ! intent(out)
                             , zrough            & ! intent(out)
                             , slz               & ! intent(out)
                             , stgoff            & ! intent(out)
                             , isoilbc           & ! intent(in)
                             , runoff_time       ! ! intent(in)
   use grid_dims      , only : maxgrds           & ! intent(out)
                             , nzgmax            ! ! intent(out)
   use canopy_air_coms, only : ubmin             & ! intent(out)
                             , ugbmin            & ! intent(out)
                             , ustmin            & ! intent(out)
                             , gamm              & ! intent(out)
                             , gamh              & ! intent(out)
                             , tprandtl          & ! intent(out)
                             , ribmax            & ! intent(out)
                             , leaf_maxwhc       ! ! intent(out)
   use rk4_coms       , only : ipercol           ! ! intent(out)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=64)         , intent(in) :: expnme_b      ! Experiment name
   integer                   , intent(in) :: ngrids_b      ! Number of grids
   character(len=16)         , intent(in) :: runtype_b     ! Simulation type
   integer                   , intent(in) :: iyeara_b      ! Initial year
   integer                   , intent(in) :: imontha_b     ! Initial month
   integer                   , intent(in) :: idatea_b      ! Initial day
   integer                   , intent(in) :: itimea_b      ! Initial hour
   integer                   , intent(in) :: iyearz_b      ! Final year
   integer                   , intent(in) :: imonthz_b     ! Final month
   integer                   , intent(in) :: idatez_b      ! Final day
   integer                   , intent(in) :: itimez_b      ! Final hour
   integer                   , intent(in) :: iyearh_b      ! History year
   integer                   , intent(in) :: imonthh_b     ! History month
   integer                   , intent(in) :: idateh_b      ! History day
   integer                   , intent(in) :: itimeh_b      ! History hour
   real                      , intent(in) :: radfrq_b      ! Radiation time step
   integer                   , intent(in) :: iclobber_b    ! Flag to decide whether files
                                                           !   can be overwritten or not
   integer,dimension(maxgrds), intent(in) :: nnxp_b        ! Number of points in X
   integer,dimension(maxgrds), intent(in) :: nnyp_b        ! Number of points in Y
   real                      , intent(in) :: deltax_b      ! Grid res. close to the pole
   real                      , intent(in) :: deltay_b      ! Grid res. close to the pole
   real                      , intent(in) :: polelat_b     ! Pole latitude (deg)
   real                      , intent(in) :: polelon_b     ! Pole longitude (deg)
   real   ,dimension(maxgrds), intent(in) :: centlat_b     ! Grid center latitude (deg)
   real   ,dimension(maxgrds), intent(in) :: centlon_b     ! Grid center longitude (deg)
   integer,dimension(maxgrds), intent(in) :: nstratx_b     ! Nest ratio for next 
                                                           !    coarser grid
   integer,dimension(maxgrds), intent(in) :: nstraty_b     ! Nest ratio for next
                                                           !    coarser grid
   integer                   , intent(in) :: nzg_b         ! Number of soil layers
   integer                   , intent(in) :: nzs_b         ! Number of snow layers
   integer,dimension(maxgrds), intent(in) :: isoilflg_b    ! Method to initialise soil type
   integer                   , intent(in) :: nslcon_b      ! Soil type if constant for 
                                                           !    all grids
   integer                   , intent(in) :: isoilcol_b    ! Soil colour if constant for 
                                                           !    all grids
   real                      , intent(in) :: zrough_b      ! Soil roughness if constant...
   real, dimension(nzgmax)   , intent(in) :: slmstr_b      ! Initial soil moist. if const.
   real, dimension(nzgmax)   , intent(in) :: stgoff_b      ! Initial soil temp. offset
   real, dimension(nzgmax)   , intent(in) :: slz_b         ! Soil layers
   real                      , intent(in) :: ubmin_b       ! Minimum u
   real                      , intent(in) :: ugbmin_b      ! Minimum u at leaf level
   real                      , intent(in) :: ustmin_b      ! Minimum u*
   integer                   , intent(in) :: isoilbc_b     ! Bottom soil boundary condition
   integer                   , intent(in) :: ipercol_b     ! Percolation scheme.
   real                      , intent(in) :: runoff_time_b ! Runoff time scale.
   real                      , intent(in) :: ribmax_b      ! Maximum bulk Richardson number
   real                      , intent(in) :: gamm_b        ! Sfc. lyr. Gamma for momentum
   real                      , intent(in) :: gamh_b        ! Sfc. lyr. Gamma for heat
   real                      , intent(in) :: tprandtl_b    ! Turbulent Prandtl number
   real                      , intent(in) :: leaf_maxwhc_b ! Leaf max. water holding cap.
   !---------------------------------------------------------------------------------------!



   !----- Copy the variables. -------------------------------------------------------------!
   expnme      = expnme_b
   ngrids      = ngrids_b
   runtype     = runtype_b
   itimez      = itimez_b
   idatez      = idatez_b
   imonthz     = imonthz_b
   iyearz      = iyearz_b
   itimeh      = itimeh_b
   idateh      = idateh_b
   imonthh     = imonthh_b
   iyearh      = iyearh_b
   itimea      = itimea_b
   idatea      = idatea_b
   imontha     = imontha_b
   iyeara      = iyeara_b

   radfrq      = radfrq_b

   iclobber    = iclobber_b

   centlon     = centlon_b
   centlat     = centlat_b
   deltax      = deltax_b
   deltay      = deltay_b
   nnxp        = nnxp_b
   nnyp        = nnyp_b
   nstratx     = nstratx_b
   nstraty     = nstraty_b
   polelat     = polelat_b
   polelon     = polelon_b
   ngrids      = ngrids_b

   nzg         = nzg_b
   nzs         = nzs_b

   slz         = slz_b
   slmstr      = slmstr_b
   stgoff      = stgoff_b
   zrough      = zrough_b

   isoilflg    = isoilflg_b


   ubmin       = ubmin_b
   ugbmin      = ugbmin_b
   ustmin      = ustmin_b
   isoilbc     = isoilbc_b
   ipercol     = ipercol_b
   runoff_time = runoff_time_b

   gamm        = gamm_b
   gamh        = gamh_b
   tprandtl    = tprandtl_b

   ribmax      = ribmax_b
   leaf_maxwhc = leaf_maxwhc_b
   isoilcol    = isoilcol_b
   !---------------------------------------------------------------------------------------!

   return
end subroutine copy_in_bramsnl
!==========================================================================================!
!==========================================================================================!
