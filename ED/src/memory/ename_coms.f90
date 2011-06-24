!==========================================================================================!
!==========================================================================================!
!     Module ename_coms.  This module contains the namelist structure, and it will be the  !
! first place where the namelist variables will be stored before the variables are loaded. !
!------------------------------------------------------------------------------------------!
module ename_coms

   use ed_max_dims, only : max_poi        & ! intent(in)
                         , max_ed_regions & ! intent(in)
                         , str_len        & ! intent(in)
                         , n_pft          & ! intent(in)
                         , maxgrds        & ! intent(in)
                         , nzgmax         & ! intent(in)
                         , maxpvars       ! ! intent(in)
   implicit none

   !---------------------------------------------------------------------------------------!
   !      This is the namelist structure.  Please assign the appropriate non-sense initial !
   ! value at the subroutine below in case you add a new variable.  This is to make the    !
   ! user aware that he or she needs to define a new variable.                             !
   !---------------------------------------------------------------------------------------!
   type ename_vars
      !----- This is just to give some name to this simulation, otherwise not used. -------!
      character(len=str_len)                            :: expnme 

      !----- Type of simulation (INITIAL or HISTORY). -------------------------------------!
      character(len=str_len)                            :: runtype

      !------ Node balance method if this is a regional, parallel run. --------------------!
      integer                                           :: loadmeth

      !------ Simulation starting time. ---------------------------------------------------!
      integer                                           :: itimea
      integer                                           :: idatea
      integer                                           :: imontha
      integer                                           :: iyeara

      !----- Simulation ending time. ------------------------------------------------------!
      integer                                           :: itimez
      integer                                           :: idatez
      integer                                           :: imonthz
      integer                                           :: iyearz

      !----- Timestep specification. ------------------------------------------------------!
      real                                              :: dtlsm
      real                                              :: radfrq

      !----- Analysis/history files. ------------------------------------------------------!
      integer                                           :: ifoutput
      integer                                           :: idoutput
      integer                                           :: imoutput
      integer                                           :: iqoutput
      integer                                           :: iyoutput
      integer                                           :: itoutput
      integer                                           :: isoutput
      integer                                           :: attach_metadata
      integer                                           :: iclobber
      integer                                           :: unitfast
      integer                                           :: unitstate
      real                                              :: frqfast
      real                                              :: frqstate
      real                                              :: outfast
      real                                              :: outstate

      character(len=str_len)                            :: ffilout
      integer                                           :: ied_init_mode

      character(len=str_len), dimension(maxgrds)        :: sfilin
      character(len=str_len)                            :: sfilout

      integer                                           :: itimeh
      integer                                           :: idateh
      integer                                           :: imonthh
      integer                                           :: iyearh

      !----- Polar stereographic projection parameters. -----------------------------------!
      integer                                           :: ngrids
      integer               , dimension(maxgrds)        :: nnxp
      integer               , dimension(maxgrds)        :: nnyp
      real                                              :: deltay
      real                                              :: deltax
      integer               , dimension(maxgrds)        :: nstratx
      integer               , dimension(maxgrds)        :: nstraty
      real                                              :: polelat
      real                                              :: polelon
      real                  , dimension(maxgrds)        :: centlat
      real                  , dimension(maxgrds)        :: centlon
      integer               , dimension(maxgrds)        :: ninest
      integer               , dimension(maxgrds)        :: njnest

      !----- Soil/surface water variables. ------------------------------------------------!
      integer                                           :: nzg
      integer                                           :: nzs
      integer               , dimension(maxgrds)        :: isoilflg
      integer                                           :: nslcon
      real                                              :: slxclay
      real                                              :: slxsand
      integer                                           :: isoilstateinit
      integer                                           :: isoildepthflg
      integer                                           :: isoilbc

      real                  , dimension(nzgmax)         :: slz
      real                  , dimension(nzgmax)         :: slmstr
      real                  , dimension(nzgmax)         :: stgoff

      !----- Input databases. -------------------------------------------------------------!
      character(len=str_len), dimension(maxgrds)        :: soil_database
      character(len=str_len), dimension(maxgrds)        :: veg_database
      character(len=str_len), dimension(maxgrds)        :: lu_database
      character(len=str_len), dimension(maxgrds)        :: plantation_file
      character(len=str_len), dimension(maxgrds)        :: lu_rescale_file
      character(len=str_len)                            :: thsums_database 
      character(len=str_len)                            :: soilstate_db
      character(len=str_len)                            :: soildepth_db
      character(len=str_len)                            :: ed_met_driver_db

      !----- ED polygon specification. ----------------------------------------------------!
      integer                                           :: n_poi
      integer                                           :: n_ed_region
      integer                                           :: grid_type
      real                                              :: grid_res
      real                  , dimension(max_poi)        :: poi_lat
      real                  , dimension(max_poi)        :: poi_lon
      real                  , dimension(max_poi)        :: poi_res
      real                  , dimension(max_ed_regions) :: ed_reg_latmin
      real                  , dimension(max_ed_regions) :: ed_reg_latmax
      real                  , dimension(max_ed_regions) :: ed_reg_lonmin
      real                  , dimension(max_ed_regions) :: ed_reg_lonmax
  

      !----- Options for model dynamics. --------------------------------------------------!
      integer                                           :: integration_scheme
      real                                              :: rk4_tolerance
      integer                                           :: ibranch_thermo
      integer                                           :: iphysiol
      integer                                           :: istoma_scheme
      integer                                           :: iallom
      integer                                           :: iphen_scheme
      integer                                           :: repro_scheme
      integer                                           :: lapse_scheme
      integer                                           :: crown_mod
      integer                                           :: ican_swrad
      integer                                           :: h2o_plant_lim
      real                                              :: vmfact
      real                                              :: mfact
      real                                              :: kfact
      real                                              :: gamfact
      real                                              :: d0fact
      real                                              :: alphafact
      real                                              :: thetacrit
      real                                              :: lwfact
      real                                              :: thioff
      integer                                           :: quantum_efficiency_T
      integer                                           :: n_plant_lim
      integer                                           :: n_decomp_lim
      integer                                           :: decomp_scheme
      integer                                           :: include_fire
      real                                              :: sm_fire
      integer                                           :: ianth_disturb
      integer                                           :: icanturb
      integer                                           :: i_blyr_condct
      integer                                           :: isfclyrm
      integer                                           :: ied_grndvap
      integer                                           :: ipercol
      integer               , dimension(n_pft)          :: include_these_pft
      integer                                           :: agri_stock
      integer                                           :: plantation_stock
      integer                                           :: pft_1st_check
      real                                              :: treefall_disturbance_rate
      real                                              :: Time2Canopy
      real                                              :: runoff_time
      real                                              :: betapower
      real                                              :: ustmin
      real                                              :: gamm
      real                                              :: gamh
      real                                              :: tprandtl
      real                                              :: vh2vr
      real                                              :: vh2dh
      real                                              :: ribmax
      real                                              :: leaf_maxwhc
      real                                              :: ggfact

      !----- Options for printing polygon vectors/arrays to standard output. --------------!
      integer                                           :: iprintpolys
      integer                                           :: npvars
      character(len=str_len), dimension(maxpvars)       :: printvars
      character(len=str_len), dimension(maxpvars)       :: pfmtstr
      integer                                           :: ipmin
      integer                                           :: ipmax

      !----- Options controlling meteorological forcing. ----------------------------------!
      integer                                           :: imettype
      integer                                           :: ishuffle
      integer                                           :: metcyc1
      integer                                           :: metcycf
      integer                                           :: imetavg
      integer                                           :: imetrad
      real                                              :: initial_co2

      !------ Options controlling prescribed phenology forcing. ---------------------------!
      integer                                           :: iphenys1
      integer                                           :: iphenysf
      integer                                           :: iphenyf1
      integer                                           :: iphenyff
      real                                              :: radint
      real                                              :: radslp

      !------ XML, phenology, and prescribed event files. ---------------------------------!
      character(len=str_len)                            :: iedcnfgf
      character(len=str_len)                            :: phenpath
      character(len=str_len)                            :: event_file

      !----- Variables that control the sought number of patches and cohorts. -------------!
      integer                                           :: maxsite
      integer                                           :: maxpatch
      integer                                           :: maxcohort
      real                                              :: min_site_area

      !----- Directory for optimizer inputs. ----------------------------------------------!
      character(len=str_len)                            :: ioptinpt

      !----- Roughness length. ------------------------------------------------------------!
      real                                              :: zrough

      !----- ED restart grid resolution. --------------------------------------------------!
      real                                              :: edres
   end Type ename_vars

   !----- This is the name of the structure containing the namelist. ----------------------!
   type(ename_vars), save :: nl

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !      The namelist structure will be read by the namelist, but we first assign         !
   ! default values, which don't make any sense, so the user will know if he or she is not !
   ! including all variables in his or her ED2IN.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine init_ename_vars(enl)
      use ed_max_dims, only : undef_real      & ! intent(in)
                            , undef_real      & ! intent(in)
                            , undef_integer   & ! intent(in)
                            , undef_character & ! intent(in)
                            , undef_path      & ! intent(in)
                            , undef_logical   ! ! intent(in)

      !----- Arguments. -------------------------------------------------------------------!
      type(ename_vars), intent(out) :: enl
      !----- Local variables. -------------------------------------------------------------!
      integer                       :: i
      !------------------------------------------------------------------------------------!

      enl%expnme                    = undef_character
      enl%runtype                   = undef_character
      enl%loadmeth                  = undef_integer
      enl%itimea                    = undef_integer
      enl%idatea                    = undef_integer
      enl%imontha                   = undef_integer
      enl%iyeara                    = undef_integer
      enl%itimez                    = undef_integer
      enl%idatez                    = undef_integer
      enl%imonthz                   = undef_integer
      enl%iyearz                    = undef_integer

      enl%dtlsm                     = undef_real
      enl%radfrq                    = undef_real

      enl%ifoutput                  = undef_integer
      enl%idoutput                  = undef_integer
      enl%imoutput                  = undef_integer
      enl%iqoutput                  = undef_integer
      enl%iyoutput                  = undef_integer
      enl%itoutput                  = undef_integer
      enl%isoutput                  = undef_integer

      enl%attach_metadata           = undef_integer

      enl%iclobber                  = undef_integer
      enl%unitfast                  = undef_integer
      enl%unitstate                 = undef_integer
      enl%frqfast                   = undef_real
      enl%frqstate                  = undef_real
      enl%outfast                   = undef_real
      enl%outstate                  = undef_real

      enl%ffilout                   = undef_path
      enl%ied_init_mode             = undef_integer

      enl%sfilin                    = (/ (undef_path, i=1,maxgrds) /)
      enl%sfilout                   = undef_path

      enl%itimeh                    = undef_integer
      enl%idateh                    = undef_integer
      enl%imonthh                   = undef_integer
      enl%iyearh                    = undef_integer

      enl%ngrids                    = undef_integer
      enl%nnxp                      = (/ (undef_integer, i=1,maxgrds) /)
      enl%nnyp                      = (/ (undef_integer, i=1,maxgrds) /)
      enl%deltay                    = undef_real
      enl%deltax                    = undef_real
      enl%nstratx                   = (/ (undef_integer, i=1,maxgrds) /)
      enl%nstraty                   = (/ (undef_integer, i=1,maxgrds) /)
      enl%polelat                   = undef_real
      enl%polelon                   = undef_real
      enl%centlat                   = (/ (undef_real, i=1,maxgrds) /)
      enl%centlon                   = (/ (undef_real, i=1,maxgrds) /)
      enl%ninest                    = (/ (undef_integer, i=1,maxgrds) /)
      enl%njnest                    = (/ (undef_integer, i=1,maxgrds) /)

      enl%nzg                       = undef_integer
      enl%nzs                       = undef_integer
      enl%isoilflg                  = (/ (undef_integer,i=1,maxgrds) /)
      enl%nslcon                    = undef_integer
      enl%slxclay                   = undef_real
      enl%slxsand                   = undef_real
      enl%isoilstateinit            = undef_integer
      enl%isoildepthflg             = undef_integer
      enl%isoilbc                   = undef_integer

      enl%slz                       = (/ (-undef_real, i=1,nzgmax) /)
      enl%slmstr                    = (/ ( undef_real, i=1,nzgmax) /)
      enl%stgoff                    = (/ ( undef_real, i=1,nzgmax) /)


      enl%soil_database             = (/ (undef_path, i=1,maxgrds) /)
      enl%veg_database              = (/ (undef_path, i=1,maxgrds) /)
      enl%lu_database               = (/ (undef_path, i=1,maxgrds) /)
      enl%plantation_file           = (/ (undef_path, i=1,maxgrds) /)
      enl%lu_rescale_file           = (/ (undef_path, i=1,maxgrds) /)

      enl%thsums_database           =     undef_path
      enl%soilstate_db              =     undef_path
      enl%soildepth_db              =     undef_path
      enl%ed_met_driver_db          =     undef_path

      enl%n_poi                     = undef_integer
      enl%n_ed_region               = undef_integer
      enl%grid_type                 = undef_integer
      enl%grid_res                  = undef_real
      enl%poi_lat                   = (/ (undef_real, i=1,max_poi) /)
      enl%poi_lon                   = (/ (undef_real, i=1,max_poi) /)
      enl%poi_res                   = (/ (undef_real, i=1,max_poi) /)
      enl%ed_reg_latmin             = (/ (undef_real,i=1,max_ed_regions) /)
      enl%ed_reg_latmax             = (/ (undef_real,i=1,max_ed_regions) /)
      enl%ed_reg_lonmin             = (/ (undef_real,i=1,max_ed_regions) /)
      enl%ed_reg_lonmax             = (/ (undef_real,i=1,max_ed_regions) /)
 

      enl%integration_scheme        = undef_integer
      enl%rk4_tolerance             = undef_real
      enl%ibranch_thermo            = undef_integer
      enl%iphysiol                  = undef_integer
      enl%istoma_scheme             = undef_integer
      enl%iallom                    = undef_integer
      enl%iphen_scheme              = undef_integer
      enl%repro_scheme              = undef_integer
      enl%lapse_scheme              = undef_integer
      enl%crown_mod                 = undef_integer
      enl%ican_swrad                = undef_integer
      enl%h2o_plant_lim             = undef_integer
      enl%vmfact                    = undef_real
      enl%mfact                     = undef_real
      enl%kfact                     = undef_real
      enl%gamfact                   = undef_real
      enl%d0fact                    = undef_real
      enl%alphafact                 = undef_real
      enl%thetacrit                 = undef_real
      enl%lwfact                    = undef_real
      enl%thioff                    = undef_real
      enl%quantum_efficiency_T      = undef_integer
      enl%n_plant_lim               = undef_integer
      enl%n_decomp_lim              = undef_integer
      enl%decomp_scheme             = undef_integer
      enl%include_fire              = undef_integer
      enl%sm_fire                   = undef_real
      enl%ianth_disturb             = undef_integer
      enl%icanturb                  = undef_integer
      enl%i_blyr_condct             = undef_integer
      enl%isfclyrm                  = undef_integer
      enl%ipercol                   = undef_integer

      enl%include_these_pft         = (/(undef_integer,i=1,n_pft)/) 
      enl%agri_stock                = undef_integer
      enl%plantation_stock          = undef_integer
      enl%pft_1st_check             = undef_integer

      enl%treefall_disturbance_rate = undef_real
      enl%Time2Canopy               = undef_real
      enl%runoff_time               = undef_real
      enl%betapower                 = undef_real
      enl%ustmin                    = undef_real
      enl%gamm                      = undef_real
      enl%gamh                      = undef_real
      enl%tprandtl                  = undef_real
      enl%vh2vr                     = undef_real
      enl%vh2dh                     = undef_real
      enl%ribmax                    = undef_real
      enl%leaf_maxwhc               = undef_real
      enl%ggfact                    = undef_real

      enl%iprintpolys               = undef_integer
      enl%npvars                    = undef_integer
      enl%printvars                 = (/ (undef_character, i=1,maxpvars) /)
      enl%pfmtstr                   = (/ (undef_character, i=1,maxpvars) /)
      enl%ipmin                     = undef_integer
      enl%ipmax                     = undef_integer

      enl%imettype                  = undef_integer
      enl%ishuffle                  = undef_integer
      enl%metcyc1                   = undef_integer
      enl%metcycf                   = undef_integer
      enl%imetavg                   = undef_integer
      enl%imetrad                   = undef_integer
      enl%initial_co2               = undef_real

      enl%iphenys1                  = undef_integer
      enl%iphenysf                  = undef_integer
      enl%iphenyf1                  = undef_integer
      enl%iphenyff                  = undef_integer
      enl%radint                    = undef_real
      enl%radslp                    = undef_real

      enl%iedcnfgf                  = undef_path
      enl%phenpath                  = undef_path
      enl%event_file                = undef_path

      enl%maxsite                   = undef_integer
      enl%maxpatch                  = undef_integer
      enl%maxcohort                 = undef_integer
      enl%min_site_area             = undef_real

      enl%ioptinpt                  = undef_path
      enl%zrough                    = undef_real
      enl%edres                     = undef_real 

      return
   end subroutine init_ename_vars
   !=======================================================================================!
   !=======================================================================================!
end module ename_coms
!==========================================================================================!
!==========================================================================================!
