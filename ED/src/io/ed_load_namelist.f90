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
                                   , maxgrds                   ! ! intent(in)
   use ename_coms           , only : nl                        ! ! intent(in)
   use soil_coms            , only : isoilflg                  & ! intent(out)
                                   , nslcon                    & ! intent(out)
                                   , slxclay                   & ! intent(out)
                                   , slxsand                   & ! intent(out)
                                   , slmstr                    & ! intent(out)
                                   , stgoff                    & ! intent(out)
                                   , zrough                    & ! intent(out)
                                   , soil_database             & ! intent(out)
                                   , isoilstateinit            & ! intent(out)
                                   , isoildepthflg             & ! intent(out)
                                   , isoilbc                   & ! intent(out)
                                   , soilstate_db              & ! intent(out)
                                   , soildepth_db              & ! intent(out)
                                   , runoff_time               & ! intent(out)
                                   , betapower                 & ! intent(out)
                                   , slz                       & ! intent(out)
                                   , veg_database              ! ! intent(out)
   use met_driver_coms      , only : ed_met_driver_db          & ! intent(out)
                                   , ishuffle                  & ! intent(out)
                                   , metcyc1                   & ! intent(out)
                                   , metcycf                   & ! intent(out)
                                   , imettype                  & ! intent(out)
                                   , initial_co2               & ! intent(out)
                                   , lapse_scheme              ! ! intent(out)
   use mem_polygons         , only : n_poi                     & ! intent(out)
                                   , poi_lat                   & ! intent(out)
                                   , poi_lon                   & ! intent(out)
                                   , n_ed_region               & ! intent(out)
                                   , ed_reg_latmin             & ! intent(out)
                                   , ed_reg_latmax             & ! intent(out)
                                   , ed_reg_lonmin             & ! intent(out)
                                   , ed_reg_lonmax             & ! intent(out)
                                   , grid_res                  & ! intent(out)
                                   , grid_type                 & ! intent(out)
                                   , edres                     & ! intent(out)
                                   , maxpatch                  & ! intent(out)
                                   , maxcohort                 ! ! intent(out)
   use physiology_coms      , only : istoma_scheme             & ! intent(out)
                                   , h2o_plant_lim             & ! intent(out)
                                   , n_plant_lim               & ! intent(out)
                                   , vmfact                    & ! intent(out)
                                   , mfact                     & ! intent(out)
                                   , kfact                     & ! intent(out)
                                   , gamfact                   & ! intent(out)
                                   , lwfact                    & ! intent(out)
                                   , thioff                    & ! intent(out)
                                   , icomppt                   & ! intent(out)
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
                                   , LloydTaylor               ! ! intent(out)
   use disturb_coms         , only : include_fire              & ! intent(out)
                                   , ianth_disturb             & ! intent(out)
                                   , treefall_disturbance_rate & ! intent(out)
                                   , lu_database               & ! intent(out)
                                   , plantation_file           & ! intent(out)
                                   , lu_rescale_file           & ! intent(out)
                                   , sm_fire                   ! ! intent(out)
   use pft_coms             , only : include_these_pft         & ! intent(out)
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
                                   , integration_scheme        & ! intent(out)
                                   , ffilout                   & ! intent(out)
                                   , idoutput                  & ! intent(out)
                                   , imoutput                  & ! intent(out)
                                   , iyoutput                  & ! intent(out)
                                   , iqoutput                  & ! intent(out)
                                   , itoutput                  & ! intent(out)
                                   , dtlsm                     & ! intent(out)
                                   , frqstate                  & ! intent(out)
                                   , sfilout                   & ! intent(out)
                                   , isoutput                  & ! intent(out)
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
                                   , event_file                ! ! intent(out)
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
   use ed_misc_coms         , only : attach_metadata           ! ! intent(out)
   use canopy_air_coms      , only : icanturb                  & ! intent(out)
                                   , isfclyrm                  & ! intent(out)
                                   , i_blyr_condct             & ! intent(out)
                                   , ustmin                    & ! intent(out)
                                   , ggfact                    & ! intent(out)
                                   , gamm                      & ! intent(out)
                                   , gamh                      & ! intent(out)
                                   , tprandtl                  & ! intent(out)
                                   , vkopr                     & ! intent(out)
                                   , vh2vr                     & ! intent(out)
                                   , vh2dh                     ! ! intent(out)
   use optimiz_coms         , only : ioptinpt                  ! ! intent(out)
   use canopy_radiation_coms, only : crown_mod                 ! ! intent(out)
   use rk4_coms             , only : ibranch_thermo            & ! intent(out)
                                   , ipercol                   & ! intent(out)
                                   , rk4_tolerance             ! ! intent(out)
   use ed_para_coms         , only : loadmeth                  ! ! intent(out)
   use consts_coms          , only : vonk                      & ! intent(in)
                                   , day_sec                   & ! intent(in)
                                   , hr_sec                    & ! intent(in)
                                   , min_sec                   ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*), intent(in) :: copy_type
   !----- Internal variables. -------------------------------------------------------------!
   integer                      :: ifm
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

      ifoutput                  = nl%ifoutput
      idoutput                  = nl%idoutput
      imoutput                  = nl%imoutput
      iqoutput                  = nl%iqoutput
      iyoutput                  = nl%iyoutput
      itoutput                  = nl%itoutput
      isoutput                  = nl%isoutput

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
      nslcon                    = nl%nslcon
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
      soilstate_db              = nl%soilstate_db
      soildepth_db              = nl%soildepth_db
      isoilstateinit            = nl%isoilstateinit
      isoildepthflg             = nl%isoildepthflg
      isoilbc                   = nl%isoilbc

      n_poi                     = nl%n_poi
      n_ed_region               = nl%n_ed_region
      grid_res                  = nl%grid_res
      grid_type                 = nl%grid_type
      poi_lat                   = nl%poi_lat
      poi_lon                   = nl%poi_lon
      ed_reg_latmin             = nl%ed_reg_latmin
      ed_reg_latmax             = nl%ed_reg_latmax
      ed_reg_lonmin             = nl%ed_reg_lonmin
      ed_reg_lonmax             = nl%ed_reg_lonmax

      integration_scheme        = nl%integration_scheme
      rk4_tolerance             = nl%rk4_tolerance
      ibranch_thermo            = nl%ibranch_thermo
      istoma_scheme             = nl%istoma_scheme
      iphen_scheme              = nl%iphen_scheme
      repro_scheme              = nl%repro_scheme
      lapse_scheme              = nl%lapse_scheme
      crown_mod                 = nl%crown_mod
      h2o_plant_lim             = nl%h2o_plant_lim
      vmfact                    = nl%vmfact
      mfact                     = nl%mfact
      kfact                     = nl%kfact
      gamfact                   = nl%gamfact
      thetacrit                 = nl%thetacrit
      lwfact                    = nl%lwfact
      thioff                    = nl%thioff
      icomppt                   = nl%icomppt
      quantum_efficiency_T      = nl%quantum_efficiency_T
      radint                    = nl%radint
      radslp                    = nl%radslp
      n_plant_lim               = nl%n_plant_lim
      n_decomp_lim              = nl%n_decomp_lim
      include_fire              = nl%include_fire
      sm_fire                   = nl%sm_fire
      ianth_disturb             = nl%ianth_disturb

      !----- Decomp_scheme is not a true ED variable, we save it in LloydTaylor instead. --!
      LloydTaylor               = nl%decomp_scheme == 1
      
      icanturb                  = nl%icanturb
      i_blyr_condct             = nl%i_blyr_condct
      isfclyrm                  = nl%isfclyrm
      gamm                      = nl%gamm
      gamh                      = nl%gamh
      tprandtl                  = nl%tprandtl
      vh2vr                     = nl%vh2vr
      vh2dh                     = nl%vh2dh
      ipercol                   = nl%ipercol

      include_these_pft         = nl%include_these_pft
      agri_stock                = nl%agri_stock
      plantation_stock          = nl%plantation_stock
      pft_1st_check             = nl%pft_1st_check
      
      treefall_disturbance_rate = nl%treefall_disturbance_rate
      runoff_time               = nl%runoff_time
      betapower                 = nl%betapower
      ustmin                    = nl%ustmin
      ggfact                    = nl%ggfact

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
      initial_co2               = nl%initial_co2
      
      iphenys1                  = nl%iphenys1
      iphenysf                  = nl%iphenysf
      iphenyf1                  = nl%iphenyf1
      iphenyff                  = nl%iphenyff
      
      iedcnfgf                  = nl%iedcnfgf
      event_file                = nl%event_file
      phenpath                  = nl%phenpath
      maxpatch                  = nl%maxpatch
      maxcohort                 = nl%maxcohort
      ioptinpt                  = nl%ioptinpt
      zrough                    = nl%zrough
      
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
   where (include_these_pft < 1 .or. include_these_pft == undef_integer) 
      include_these_pft=huge(1)
   end where
   call sort_up(include_these_pft,n_pft)
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



   !----- Find von-Karman/Prandtl number ratio. -------------------------------------------!
   if (tprandtl /= 0.0) then
      vkopr = vonk / tprandtl
   else
      !------------------------------------------------------------------------------------!
      !     It doesn't make sense, but tprandtl is wrong and the run will crash at         !
      ! ed_opspec.                                                                         !
      !------------------------------------------------------------------------------------!
      vkopr = 0.0
   end if 
   !---------------------------------------------------------------------------------------!

   return
end subroutine copy_nl
!==========================================================================================!
!==========================================================================================!
