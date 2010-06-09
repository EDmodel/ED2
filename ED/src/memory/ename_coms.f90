Module ename_coms

  use ed_max_dims, only: max_poi, max_ed_regions, str_len,n_pft,maxgrds, nzgmax,maxpvars

  implicit none

  integer :: i
  
  ! The following variables will be read by the namelist.  The default is to use non-sense
  ! numbers so the user will know if he or she is not including all variables in his or 
  ! her ED2IN

  Type ename_vars

!!    RUNTYPE/NAME

      character(len=str_len) :: expnme   = 'ED2 test'
      character(len=str_len) :: runtype  = ''

!!   LOAD METHOD (default is to equally split poly's on nodes)

      integer :: loadmeth = -1
      

!!    START OF SIMULATION

      integer           :: itimea  = -999
      integer           :: idatea  = -999
      integer           :: imontha = -999
      integer           :: iyeara  = -999

!!    SIMULATION ENDING TIME

      integer           :: itimez  = -999
      integer           :: idatez  = -999
      integer           :: imonthz = -999
      integer           :: iyearz  = -999
      
!!    TIMESTEP SPECIFICATION

      real              :: dtlsm  = -999.9
      real              :: radfrq = -999.9

!!    ANALYSIS/HISTORY FILES

      integer :: ifoutput  = -999
      integer :: idoutput  = -999
      integer :: imoutput  = -999
      integer :: iyoutput  = -999
      integer :: itoutput  = -999
      integer :: isoutput  = -999
      
      integer :: attach_metadata = -999

      integer :: iclobber    = -999
      integer :: unitfast    = -999
      integer :: unitstate   = -999
      real    :: frqfast     = -999.9
      real    :: frqstate    = -999.9
      real    :: outfast     = -999.9
      real    :: outstate    = -999.9

      character(len=str_len) :: ffilout   = '/nowhere/X'
      integer                :: ied_init_mode = -999

      character(len=str_len) :: sfilin    = '/nowhere/X'
      character(len=str_len) :: sfilout   = '/nowhere/X'

      integer           :: itimeh  = -999
      integer           :: idateh  = -999
      integer           :: imonthh = -999
      integer           :: iyearh  = -999


!!    POLAR STEREOGRAPHIC PROJECTION PARAMETERS
      integer                     :: ngrids  = -999
      integer, dimension(maxgrds) :: nnxp    = (/ (-999, i=1,maxgrds) /)
      integer, dimension(maxgrds) :: nnyp    = (/ (-999, i=1,maxgrds) /)
      real                        :: deltay  = -999.9
      real                        :: deltax  = -999.9
      integer, dimension(maxgrds) :: nstratx = (/ (-999, i=1,maxgrds) /)
      integer, dimension(maxgrds) :: nstraty = (/ (-999, i=1,maxgrds) /)
      real                        :: polelat = -999.9
      real                        :: polelon = -999.9
      real   , dimension(maxgrds) :: centlat = (/ (-999.9, i=1,maxgrds) /)
      real   , dimension(maxgrds) :: centlon = (/ (-999.9, i=1,maxgrds) /)
      integer, dimension(maxgrds) :: ninest  = (/ (-999, i=1,maxgrds) /)
      integer, dimension(maxgrds) :: njnest  = (/ (-999, i=1,maxgrds) /)

!!    SOIL/SURFACE WATER VARIABLES
      integer :: nzg      = -999
      integer :: nzs      = -999
      integer, dimension(maxgrds) :: isoilflg = (/ (-999,i=1,maxgrds) /)
      integer :: nslcon   = -999
      real    :: slxclay  = -999.           ! Use default clay fraction (NML)
      real    :: slxsand  = -999.           ! Use default sand fraction (NML)
      integer :: isoilstateinit = -999
      integer :: isoildepthflg  = -999
      integer :: isoilbc        = -999

      real, dimension(nzgmax) :: slz    = (/ ( 999.9, i=1,nzgmax) /)
      real, dimension(nzgmax) :: slmstr = (/ (-999.9, i=1,nzgmax) /)
      real, dimension(nzgmax) :: stgoff = (/ (-999.9, i=1,nzgmax) /)

!!    INPUT DATABASES

      character(len=str_len) :: soil_database    = '/nowhere/X'
      character(len=str_len) :: veg_database     = '/nowhere/X'
      character(len=str_len) :: ed_inputs_dir    = '/nowhere/X'
      character(len=str_len) :: ed_met_driver_db = '/nowhere/X'
      character(len=str_len) :: soilstate_db     = '/nowhere/X'
      character(len=str_len) :: soildepth_db     = '/nowhere/'

!!    ED SITE SPECIFICATION

      integer           :: n_poi         = -999
      integer           :: n_ed_region   = -999
      integer           :: grid_type     = -999
      real              :: grid_res      = -999.9
      real, dimension(max_poi) :: poi_lat = (/ (-999.9, i=1,max_poi) /)
      real, dimension(max_poi) :: poi_lon = (/ (-999.9, i=1,max_poi) /)
      real, dimension(max_ed_regions) :: ed_reg_latmin = (/ (-999.9,i=1,max_ed_regions) /)
      real, dimension(max_ed_regions) :: ed_reg_latmax = (/ (-999.9,i=1,max_ed_regions) /)
      real, dimension(max_ed_regions) :: ed_reg_lonmin = (/ (-999.9,i=1,max_ed_regions) /)
      real, dimension(max_ed_regions) :: ed_reg_lonmax = (/ (-999.9,i=1,max_ed_regions) /)
 

!!    OPTIONS FOR MODEL DYNAMICS
      integer                   :: integration_scheme = -999
      real                      :: rk4_tolerance      = -999.9
      integer                   :: ibranch_thermo     = -999
      integer                   :: istoma_scheme      = -999
      integer                   :: iphen_scheme       = -999
      integer                   :: repro_scheme       = -999
      integer                   :: lapse_scheme       = -999
      integer                   :: crown_mod          = -999
      integer                   :: n_plant_lim        = -999
      integer                   :: n_decomp_lim       = -999
      integer                   :: decomp_scheme      = -999
      integer                   :: include_fire       = -999
      integer                   :: ianth_disturb      = -999
      integer                   :: icanturb           = -999
      integer                   :: isfclyrm           = -999
      
      ! Huge(1) will initialize with the maximum representable number, which 
      !   will be ignored by ED, which include pfts that are <= n_pft only.
      integer, dimension(n_pft) :: include_these_pft = (/(huge(1),i=1,n_pft)/) 
      integer                   :: agri_stock        = -999
      integer                   :: plantation_stock  = -999
      integer                   :: pft_1st_check     = -999
      
      real                      :: treefall_disturbance_rate = -999.9
      real                      :: runoff_time               = -999.9

!!    OPTIONS FOR PRINTING POLYGON VECTORS/ARRAYS TO STANDARD OUTPUT
      integer :: iprintpolys = -999
      integer :: npvars = -999
      character(len=str_len),dimension(maxpvars) :: printvars = (/ ('', i=1,maxpvars) /)
      character(len=str_len),dimension(maxpvars) :: pfmtstr = (/ ('', i=1,maxpvars) /)
      integer :: ipmin = -999
      integer :: ipmax = -999

!!    OPTIONS CONTROLLING METEOROLOGICAL FORCING
      integer           :: imettype      = -999
      integer           :: ishuffle      = -999
      integer           :: metcyc1       = -999
      integer           :: metcycf       = -999
      real              :: initial_co2   = -999.

!!    OPTIONS CONTROLLING PRESCRIBED PHENOLOGY FORCING
      integer           :: iphenys1       = -999
      integer           :: iphenysf       = -999
      integer           :: iphenyf1       = -999
      integer           :: iphenyff       = -999

!!    XML CONFIGURATION FILE
      character(len=str_len) :: iedcnfgf='/nowhere'

!!    phenology file
      character(len=str_len) :: phenpath='/nowhere'

!!    XML EVENT FILE
      character(len=str_len) :: event_file='/nowhere'

!!    Variables that control the sought number of patches and cohorts
      integer :: maxpatch  = 999 ! Maximum # of patches
      integer :: maxcohort = 999 ! Maximum # of cohorts
      
  
!! Directory for optimizer inputs
      character(len=str_len) :: ioptinpt  ='/nowhere'

!! Roughness length
      real :: zrough=-999.9

!! ED restart grid resolution
      real :: edres=-999.9
      
   End Type ename_vars

   type(ename_vars),save :: nl

end Module ename_coms
