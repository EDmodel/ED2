Module ename_coms

  use max_dims, only: max_soi, max_ed_regions, str_len,n_pft,maxgrds, nzgmax,maxpvars

  implicit none

  integer :: i
  
  ! DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - NOT ALL COMPILERS WILL ACTUALLY INITIALIZE THEM
  ! THE FOLLOWING ARE INITIALIZATIONS OF NON-PARAMETERS IN THEIR MODULE DEFINING THEM
  ! IT IS LIKELY THEY WILL BE OVER WRITTEN IMMEDIATELY DURING THE NAMELIST READING
  ! BUT THESE ARE EARMARKED TO BE INITIALIZED OUTSIDE THIS ROUTINE RGK 7-18-08

  Type ename_vars

!!    RUNTYPE/NAME

      character(len=str_len) :: expnme   = 'ED2 test'
      character(len=str_len) :: runtype  = ''

!!    START OF SIMULATION

      integer           :: itimea = 0
      integer           :: idatea = 0
      integer           :: imontha = 0
      integer           :: iyeara = 0

!!    SIMULATION ENDING TIME

      integer           :: itimez = 0
      integer           :: idatez = 0
      integer           :: imonthz = 0
      integer           :: iyearz = 0
      
!!    TIMESTEP SPECIFICATION

      real              :: dtlsm = 0.0
      real              :: radfrq = 1800.0

!!    ANALYSIS/HISTORY FILES

      integer :: ifoutput  = 1      
      integer :: idoutput  = 1
      integer :: imoutput  = 1
      integer :: iyoutput  = 0
      integer :: isoutput  = 1
      
      integer :: attach_metadata=0

      integer :: iclobber    = 0
      integer :: unitfast    = 0
      integer :: unitstate   = 0
      real    :: frqfast     = 3600.0
      real    :: frqstate    = 21600.0
      real    :: outfast     = -1.
      real    :: outstate    = -1.

      character(len=str_len) :: ffilout   = 'F/'
      integer           :: ied_init_mode = 0

      character(len=str_len) :: sfilin    = ''
      character(len=str_len) :: sfilout   = 'S/'

      integer           :: itimeh = 0
      integer           :: idateh = 0
      integer           :: imonthh = 0
      integer           :: iyearh = 0


!!    POLAR STEREOGRAPHIC PROJECTION PARAMETERS
      integer :: ngrids
      integer, dimension(maxgrds) :: nnxp =(/ (0, i=1,maxgrds) /)
      integer, dimension(maxgrds) :: nnyp =(/ (0, i=1,maxgrds) /)
      real    :: deltay = 0.0
      real    :: deltax = 0.0
      integer, dimension(maxgrds) :: nstratx = (/ (0, i=1,maxgrds) /)
      integer, dimension(maxgrds) :: nstraty = (/ (0, i=1,maxgrds) /)
      real    :: polelat = 0.
      real    :: polelon = 0.
      real, dimension(maxgrds)    :: centlat = (/ (0.0, i=1,maxgrds) /)
      real, dimension(maxgrds)    :: centlon = (/ (0.0, i=1,maxgrds) /)
      integer, dimension(maxgrds)    :: ninest = (/ (0, i=1,maxgrds) /)
      integer, dimension(maxgrds)    :: njnest = (/ (0, i=1,maxgrds) /)

!!    SOIL/SURFACE WATER VARIABLES
      integer :: nzg      = 11
      integer :: nzs      = 1
      integer, dimension(maxgrds) :: isoilflg = (/ (2,i=1,maxgrds) /)
      integer :: nslcon   = 6
      integer :: isoilstateinit = 0
      integer :: isoildepthflg  = 0
      integer :: isoilbc        = 0

      real, dimension(nzgmax) :: slz = (/ -1.00, -.85, -.70, -.60, -.50, -.40, &
           -.30, -.20, -.15, -.10, -.05, (0.0, i=12,nzgmax) /)
      
      real, dimension(nzgmax) :: slmstr = (/ .35, .35, .35, .35, .35, .35, .35, &
           .35, .35, .35, .35, (0.0, i=12,nzgmax) /)

      real, dimension(nzgmax) :: stgoff = (/ (0.0, i=1,nzgmax) /)

!!    INPUT DATABASES

      character(len=str_len) :: soil_database = ''
      character(len=str_len) :: veg_database = ''
      character(len=str_len) :: ed_inputs_dir = ''
      character(len=str_len) :: ed_met_driver_db = ''
      character(len=str_len) :: soilstate_db  = ''
      character(len=str_len) :: soildepth_db  = ''

!!    ED SITE SPECIFICATION

      integer           :: n_soi         = 1
      integer           :: n_ed_region   = 0
      integer           :: grid_type     = 0
      real              :: grid_res      = 1.0
      real, dimension(max_soi) :: soi_lat = (/ (0.0, i=1,max_soi) /)
      real, dimension(max_soi) :: soi_lon = (/ (0.0, i=1,max_soi) /)
      real, dimension(max_ed_regions) :: ed_reg_latmin =   &
           (/ (0.0, i=1,max_ed_regions) /)
      real, dimension(max_ed_regions) :: ed_reg_latmax = &
           (/ (0.0, i=1,max_ed_regions) /)
      real, dimension(max_ed_regions) :: ed_reg_lonmin = &
           (/ (0.0, i=1,max_ed_regions) /)
      real, dimension(max_ed_regions) :: ed_reg_lonmax = &
           (/ (0.0, i=1,max_ed_regions) /)
 

!!    OPTIONS FOR MODEL DYNAMICS
      integer                   :: integration_scheme = 0
      integer                   :: ibranch_thermo     = 0
      integer                   :: istoma_scheme      = 0
      integer                   :: iphen_scheme       = 0
      integer                   :: repro_scheme       = 1
      integer                   :: lapse_scheme       = 0
      integer                   :: crown_mod          = 0
      integer                   :: n_plant_lim        = 0
      integer                   :: n_decomp_lim       = 0
      integer                   :: include_fire       = 0
      integer                   :: ianth_disturb      = 0
      integer                   :: icanturb           = 0
      
      ! Huge(1) will initialize with the maximum representable number, which 
      !   will be ignored by ED, which include pfts that are <= n_pft only.
      integer, dimension(n_pft) :: include_these_pft = (/(huge(1),i=1,n_pft)/) 
      integer                   :: agri_stock = 0
      integer                   :: plantation_stock = 0
      integer                   :: pft_1st_check = 0
      
      real              :: treefall_disturbance_rate = 0.0
      real              :: runoff_time   = 0.0

!!    OPTIONS FOR PRINTING POLYGON VECTORS/ARRAYS TO STANDARD OUTPUT
      integer :: iprintpolys = 0
      integer :: npvars = 0
      character(len=32),dimension(maxpvars) :: printvars = (/ ('', i=1,maxpvars) /)
      character(len=32),dimension(maxpvars) :: pfmtstr = (/ ('', i=1,maxpvars) /)
      integer :: ipmin = 0
      integer :: ipmax = 0

!!    OPTIONS CONTROLLING METEOROLOGICAL FORCING
      integer           :: imettype      = 0
      integer           :: metcyc1       = 0
      integer           :: metcycf       = 0
      real              :: initial_co2   = 370.

!!    OPTIONS CONTROLLING PRESCRIBED PHENOLOGY FORCING
      integer           :: iphenys1       = 0
      integer           :: iphenysf       = 0
      integer           :: iphenyf1       = 0
      integer           :: iphenyff       = 0

!!    XML CONFIGURATION FILE
      character(len=str_len) :: iedcnfgf=''

!!    phenology file
      character(len=str_len) :: phenpath=''

!!    XML EVENT FILE
      character(len=str_len) :: event_file=''

!!    VARIABLES THAT WILL EVENTUALLY DISAPPEAR
      integer :: maxpatch  = 12 ! Maximum # of patches
      integer :: maxcohort = 12 ! Maximum # of cohorts
      
  
!! Directory for optimizer inputs
      character(len=str_len) :: ioptinpt  =''

!! Roughness length
      real :: zrough=0.1

!! ED restart grid resolution
      real :: edres=1.0
      
   End Type ename_vars

   type(ename_vars),save :: nl

end Module ename_coms
